"""
tools/peaks.py
============
MDC and EDC peak fitting, with no Qt: a band model that carries its own
identity and the energy range it exists over, hand-placed seeds
interpolated into a starting guess for every line, a variable-projection
solver, and -- for EDCs -- the Fermi cut-off in the model rather than
divided out of the data first.

What this takes from the lab's MATLAB tools, and what it does differently
--------------------------------------------------------------------------
**Seeding, from Han Peng's** ``fit_MDC_demo``. Mark a band at four or five
energies and the rest is interpolated. Automatic peak finding fails exactly
where it matters -- band crossings, weak features, two peaks merging -- and
a person can point at the band in seconds. This is the single most useful
idea in those files. The interpolation here is **PCHIP, not a spline**: a
cubic spline overshoots between widely spaced seeds and can put a starting
centre outside the data, and a monotone interpolant cannot.

**Bands persist, and know where they live.** A band is one object across
the whole series, with a range taken from its seeds, so a band that only
exists below E_F is fitted only there. A fixed peak count over a whole cut
is wrong for real data.

**Variable projection, from O'Haver's** ``peakfit``. Only centres and widths
go to the optimiser; the heights and the background are solved exactly by
linear least squares at every step. That halves the nonlinear dimension and
takes the strongly correlated height parameters out of the search.

**Bounds that actually apply.** Both MATLAB fitters compute bounds and then
call ``fminsearch``, which cannot accept them -- in ``fit_MDC_demo`` the
arrays ``x_lb`` and ``x_ub`` are built and then simply dropped. Here they go
to ``least_squares``, which enforces them.

**Errors on everything**, from the full Jacobian at the solution, with
Poisson weights. Neither MATLAB tool reports an uncertainty on anything.
``peakfit`` has a "bootstrap", but it randomly swaps *adjacent* points
rather than resampling, which on smooth oversampled data barely perturbs
the fit and so reports a standard deviation far smaller than the truth.

**The Fermi cut-off belongs in the model.** An EDC near E_F is
``[A(E) f(E,T)] * R``, not ``A(E)`` -- the occupation multiplies the
spectrum *before* the resolution smears it, so a peak sitting at E_F comes
out shifted and clipped if the cut-off is ignored. Dividing it out of the
data first is worse still: it divides the noise up as well, and above E_F
it divides by a number known only to within that noise.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import optimize
from scipy.interpolate import PchipInterpolator
from scipy.special import wofz

from tools.fermi import K_B, SIGMA_TO_FWHM, fermi_function, _convolution_offsets
from tools.process import shirley_background as _shirley_background

#: Line shapes a band can take.
PEAK_SHAPES = ("lorentzian", "gaussian", "voigt")
#: Backgrounds that can sit under the peaks of one line.
BACKGROUNDS = ("none", "constant", "linear", "shirley")


# ==========================================================================
# Line shapes
# ==========================================================================
def peak_profile(shape: str, x, centre: float, width: float,
                 resolution: float = 0.0):
    """One unit-height peak. ``width`` is the FWHM throughout.

    FWHM, never the half width: ``fit_MDC_demo`` uses the half width in its
    Lorentzian and draws its width handle at +/- that, so a width carried
    from there to anywhere else is wrong by a factor of two. One convention,
    stated once.
    """
    x = np.asarray(x, dtype=float)
    width = max(abs(float(width)), 1e-12)
    resolution = max(float(resolution), 0.0)

    if shape == "gaussian":
        full = float(np.hypot(width, resolution))
        sigma = full / SIGMA_TO_FWHM
        return np.exp(-0.5 * ((x - centre) / sigma) ** 2)
    if shape == "lorentzian" and resolution <= 0:
        half = width / 2.0
        return half ** 2 / ((x - centre) ** 2 + half ** 2)

    # Voigt: the true profile from the Faddeeva function, not the
    # pseudo-Voigt of peakfit.m -- and parameterised the way an ARPES
    # measurement actually constrains it. The Gaussian width is the
    # instrument resolution, known in meV; the Lorentzian width is the
    # unknown being fitted. peakfit.m takes their *ratio* as a fixed input,
    # which is the one thing nobody knows in advance.
    sigma = max(resolution, 1e-12) / SIGMA_TO_FWHM
    gamma = width / 2.0
    scale = sigma * np.sqrt(2.0)
    profile = np.real(wofz(((x - centre) + 1j * gamma) / scale))
    peak = float(np.real(wofz(1j * gamma / scale)))
    return profile / max(peak, 1e-30)


#: The Shirley inelastic background of one spectrum. It is the same curve the
#: 2-D background tab subtracts, so a fit that removes it and a picture that
#: removes it agree; it lives in :mod:`tools.process` because that is where the
#: other background models are.
shirley_background = _shirley_background


# ==========================================================================
# Bands and their seeds
# ==========================================================================
@dataclass
class Seed:
    """One hand-placed point: "this band is here, this wide, at this line"."""
    position: float           # where on the axis the lines are spaced along
    centre: float             # where the peak sits along the line
    width: float              # FWHM
    height: float = 0.0


@dataclass
class Band:
    """One band, across the whole series.

    ``share_width`` names another band whose width this one copies, which is
    how the two branches of a symmetric dispersion are held to the same
    width -- the constrained-shape idea from ``peakfit``'s equal-width
    variants, but per band rather than all-or-nothing.
    """
    name: str
    seeds: list = field(default_factory=list)
    shape: str = "lorentzian"
    share_width: str = ""
    fix_width: bool = False
    colour: int = 0

    def sorted_seeds(self):
        return sorted(self.seeds, key=lambda s: s.position)

    def span(self):
        """The range of lines this band exists over, from its own seeds."""
        if not self.seeds:
            return None
        positions = [s.position for s in self.seeds]
        return (min(positions), max(positions))

    def covers(self, position: float, tol: float = 0.0) -> bool:
        span = self.span()
        if span is None:
            return False
        return (span[0] - tol) <= float(position) <= (span[1] + tol)

    def guess(self, positions):
        """Starting centre, width and height at each of ``positions``.

        PCHIP rather than a cubic spline. Between seeds placed a few tens of
        meV apart a spline overshoots, and an overshooting starting centre
        lands outside the data where the fit has nothing to hold on to;
        a monotone interpolant cannot overshoot by construction. Outside the
        seeded range the value is held at the nearest seed rather than
        extrapolated, for the same reason.
        """
        positions = np.atleast_1d(np.asarray(positions, dtype=float))
        seeds = self.sorted_seeds()
        if not seeds:
            raise ValueError(f"band {self.name!r} has no seeds")
        known = np.array([s.position for s in seeds], dtype=float)
        out = []
        for values in (np.array([s.centre for s in seeds], dtype=float),
                       np.array([s.width for s in seeds], dtype=float),
                       np.array([s.height for s in seeds], dtype=float)):
            if known.size == 1:
                out.append(np.full(positions.shape, float(values[0])))
                continue
            curve = PchipInterpolator(known, values, extrapolate=False)
            sampled = curve(np.clip(positions, known[0], known[-1]))
            out.append(np.asarray(sampled, dtype=float))
        return out[0], out[1], out[2]


@dataclass
class FitSettings:
    """Everything about the fit that is not a band."""
    background: str = "linear"
    resolution: float = 0.0             # Gaussian FWHM, instrument
    weighting: str = "poisson"
    #: EDC only: multiply the model by the Fermi-Dirac occupation.
    fermi: bool = False
    ef: float = 0.0
    temperature: float = 30.0
    max_iterations: int = 400

    def background_columns(self) -> int:
        return {"none": 0, "shirley": 0, "constant": 1, "linear": 2}[
            self.background]


# ==========================================================================
# One line
# ==========================================================================
@dataclass
class LineFit:
    """One fitted MDC or EDC."""
    x: np.ndarray
    y: np.ndarray
    model: np.ndarray
    background: np.ndarray
    names: list
    centres: np.ndarray
    centre_errors: np.ndarray
    widths: np.ndarray
    width_errors: np.ndarray
    heights: np.ndarray
    height_errors: np.ndarray
    areas: np.ndarray
    position: float = float("nan")
    chi2: float = float("nan")
    r_squared: float = float("nan")
    shapes: list = field(default_factory=list)
    resolution: float = 0.0
    success: bool = True
    message: str = ""

    @property
    def residual(self) -> np.ndarray:
        return self.y - self.model

    def component(self, index: int) -> np.ndarray:
        """One band's own peak, on top of the background."""
        return self.background + self.heights[index] * self._shape(index)

    def _shape(self, index: int) -> np.ndarray:
        return peak_profile(self.shapes[index], self.x, self.centres[index],
                            self.widths[index], self.resolution)


def _weights(y, weighting: str):
    if weighting == "poisson":
        # sqrt(N) is the counting error; the floor of 1 keeps an empty
        # channel from being given infinite weight.
        return np.sqrt(np.maximum(np.abs(y), 1.0))
    if weighting == "none":
        return np.ones_like(y)
    raise ValueError(f"weighting {weighting!r} is not poisson or none")


def _occupation(x, settings: FitSettings):
    """``(offsets, weights, f)`` for the resolution/Fermi convolution.

    Returned as a grid so the same machinery serves every column: a column
    is ``sum_o w_o * profile(x - o) * f(x - o)``, which is the occupation
    multiplying the spectrum *before* the resolution smears it. Doing it the
    other way round -- smearing each separately -- puts a peak at E_F in the
    wrong place.
    """
    x = np.asarray(x, dtype=float)
    sigma = max(float(settings.resolution), 0.0) / SIGMA_TO_FWHM
    if not settings.fermi:
        return None, None, None
    if sigma <= 0:
        return (np.zeros(1), np.ones(1),
                fermi_function(x, settings.ef, settings.temperature)[:, None])
    kt = K_B * max(float(settings.temperature), 1e-6)
    offsets, weights = _convolution_offsets(sigma, kt)
    grid = x[:, None] - offsets[None, :]
    return offsets, weights, fermi_function(grid, settings.ef,
                                            settings.temperature)


def _design(x, centres, widths, shapes, settings: FitSettings, fermi_cache):
    """The linear model's columns: one per band, then the background."""
    x = np.asarray(x, dtype=float)
    columns = []
    offsets, weights, occupation = fermi_cache
    for centre, width, shape in zip(centres, widths, shapes):
        if occupation is None:
            columns.append(peak_profile(shape, x, centre, width,
                                        settings.resolution))
        else:
            grid = x[:, None] - offsets[None, :]
            profile = peak_profile(shape, grid, centre, width, 0.0)
            columns.append((profile * occupation) @ weights)
    n_background = settings.background_columns()
    if n_background >= 1:
        # The background is detector dark counts and secondaries. It is not
        # modulated by the occupation, so it does not get the Fermi factor.
        columns.append(np.ones_like(x))
    if n_background >= 2:
        columns.append(x - float(np.mean(x)))
    return np.column_stack(columns) if columns else np.zeros((x.size, 0))


def _solve_linear(design, y, sigma, n_peaks):
    """Heights and background by weighted linear least squares.

    Heights are held at or above zero; the background is free. ``peakfit``
    does the non-negativity with ``abs(A\\y)``, which flips the sign of a
    negative height instead of forbidding it and so hides the fit that
    wanted one. A bounded solve refuses it honestly.
    """
    if design.shape[1] == 0:
        return np.zeros(0)
    scaled = design / sigma[:, None]
    target = y / sigma
    lower = np.concatenate([np.zeros(n_peaks),
                            np.full(design.shape[1] - n_peaks, -np.inf)])
    upper = np.full(design.shape[1], np.inf)
    try:
        return optimize.lsq_linear(scaled, target, bounds=(lower, upper),
                                   method="bvls").x
    except Exception:
        return np.linalg.lstsq(scaled, target, rcond=None)[0]


def fit_line(x, y, bands, settings: FitSettings, *, guesses=None,
             position: float = float("nan")) -> LineFit:
    """Fit one MDC or EDC with the given bands.

    ``guesses`` is an optional ``{name: (centre, width, height)}`` from the
    seeds or from the previous line; without it each band's own seeds are
    interpolated to ``position``.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    good = np.isfinite(x) & np.isfinite(y)
    if good.sum() < 3 * len(bands) + 2:
        raise ValueError("not enough valid points for that many bands")
    order = np.argsort(x[good])
    xf, yf = x[good][order], y[good][order]

    bands = list(bands)
    shapes = [b.shape for b in bands]
    names = [b.name for b in bands]
    span = float(abs(xf[-1] - xf[0])) or 1.0
    step = span / max(xf.size - 1, 1)

    fixed = np.zeros_like(yf)
    if settings.background == "shirley":
        fixed = shirley_background(xf, yf)

    # -- starting values and bounds -------------------------------------
    centres0, widths0 = [], []
    for band in bands:
        if guesses and band.name in guesses:
            centre, width, _ = guesses[band.name]
        else:
            centre, width, _ = (float(v[0]) for v in band.guess([position]))
        centres0.append(float(centre))
        widths0.append(max(float(width), step))

    # Which parameters are actually free. A band that shares its width with
    # another contributes no width parameter of its own, which is what makes
    # the constraint a constraint rather than a penalty.
    width_owner = {}
    for i, band in enumerate(bands):
        if band.share_width and band.share_width in names:
            width_owner[i] = names.index(band.share_width)
        else:
            width_owner[i] = i
    free_widths = [i for i, band in enumerate(bands)
                   if width_owner[i] == i and not band.fix_width]

    params0, lower, upper = [], [], []
    for i in range(len(bands)):
        params0.append(centres0[i])
        # The bound is taken from the *seed*, not from the whole axis: a
        # peak is not allowed to wander more than its own width away from
        # where it was said to be, which is what stops two bands swapping
        # identity at a crossing.
        lower.append(centres0[i] - 1.5 * widths0[i])
        upper.append(centres0[i] + 1.5 * widths0[i])
    for i in free_widths:
        params0.append(widths0[i])
        lower.append(step / 2.0)
        upper.append(5.0 * span)
    params0 = np.asarray(params0, dtype=float)
    lower = np.asarray(lower, dtype=float)
    upper = np.asarray(upper, dtype=float)
    params0 = np.clip(params0, lower + 1e-12, upper - 1e-12)

    sigma = _weights(yf, settings.weighting)
    fermi_cache = _occupation(xf, settings)
    n_bands = len(bands)

    def unpack(params):
        centres = np.asarray(params[:n_bands], dtype=float)
        widths = np.array(widths0, dtype=float)
        for slot, i in enumerate(free_widths):
            widths[i] = params[n_bands + slot]
        for i in range(n_bands):
            widths[i] = widths[width_owner[i]]
        return centres, widths

    def residual(params):
        centres, widths = unpack(params)
        design = _design(xf, centres, widths, shapes, settings, fermi_cache)
        coeffs = _solve_linear(design, yf - fixed, sigma, n_bands)
        model = design @ coeffs + fixed
        return (model - yf) / sigma

    result = optimize.least_squares(residual, params0, bounds=(lower, upper),
                                    max_nfev=settings.max_iterations * 10)
    centres, widths = unpack(result.x)
    design = _design(xf, centres, widths, shapes, settings, fermi_cache)
    coeffs = _solve_linear(design, yf - fixed, sigma, n_bands)
    model = design @ coeffs + fixed
    heights = np.asarray(coeffs[:n_bands], dtype=float)
    background = fixed + (design[:, n_bands:] @ coeffs[n_bands:]
                          if design.shape[1] > n_bands else 0.0)

    # -- errors -----------------------------------------------------------
    # Variable projection is for the *search*; the covariance is taken from
    # the full parameter set -- centres, widths and heights together --
    # because those are what the errors are wanted on, and the heights are
    # correlated with the widths in a way a projected Jacobian hides.
    full = np.concatenate([centres, [widths[i] for i in free_widths], heights,
                           coeffs[n_bands:]])

    def full_model(params):
        c = np.asarray(params[:n_bands], dtype=float)
        w = np.array(widths0, dtype=float)
        for slot, i in enumerate(free_widths):
            w[i] = params[n_bands + slot]
        for i in range(n_bands):
            w[i] = w[width_owner[i]]
        rest = params[n_bands + len(free_widths):]
        d = _design(xf, c, w, shapes, settings, fermi_cache)
        return (d @ rest + fixed - yf) / sigma

    dof = max(int(yf.size - full.size), 1)
    chi2 = float(np.sum(((model - yf) / sigma) ** 2) / dof)
    errors = np.full(full.size, np.nan)
    try:
        jacobian = _numeric_jacobian(full_model, full)
        _, singular, vt = np.linalg.svd(jacobian, full_matrices=False)
        keep = singular > (np.finfo(float).eps * max(jacobian.shape)
                           * singular[0])
        covariance = (vt[keep].T / singular[keep] ** 2) @ vt[keep]
        errors = np.sqrt(np.maximum(np.diag(covariance) * chi2, 0.0))
    except (np.linalg.LinAlgError, IndexError, ValueError, FloatingPointError):
        pass

    centre_errors = errors[:n_bands]
    width_errors = np.full(n_bands, np.nan)
    for slot, i in enumerate(free_widths):
        width_errors[i] = errors[n_bands + slot]
    for i in range(n_bands):
        width_errors[i] = width_errors[width_owner[i]]
    height_errors = errors[n_bands + len(free_widths):
                           n_bands + len(free_widths) + n_bands]

    spread = float(np.sum((yf - np.mean(yf)) ** 2))
    r_squared = (1.0 - float(np.sum((yf - model) ** 2)) / spread
                 if spread > 0 else float("nan"))
    areas = np.array([heights[i] * _area_factor(shapes[i], widths[i],
                                                settings.resolution)
                      for i in range(n_bands)])

    return LineFit(xf, yf, model, np.asarray(background) * np.ones_like(xf),
                   names, centres, centre_errors, widths, width_errors,
                   heights, height_errors, areas, float(position), chi2,
                   r_squared, shapes, settings.resolution,
                   bool(result.success), str(result.message))


def _area_factor(shape: str, width: float, resolution: float) -> float:
    """Integral of a unit-height peak of this shape and width."""
    if shape == "gaussian":
        full = float(np.hypot(width, resolution))
        return full * np.sqrt(np.pi / (4.0 * np.log(2.0)))
    if shape == "lorentzian" and resolution <= 0:
        return float(np.pi * width / 2.0)
    # Voigt: integrate the profile once, which is exact enough and avoids a
    # closed form that is only valid for the normalised shape.
    grid = np.linspace(-20.0 * (width + resolution), 20.0 * (width + resolution),
                       4001)
    return float(np.trapezoid(peak_profile(shape, grid, 0.0, width, resolution),
                              grid))


def _numeric_jacobian(function, params, relative: float = 1e-6):
    """Central-difference Jacobian, stepped relative to each parameter."""
    params = np.asarray(params, dtype=float)
    base = function(params)
    jacobian = np.empty((base.size, params.size))
    for i in range(params.size):
        step = relative * max(abs(params[i]), 1.0)
        forward, backward = params.copy(), params.copy()
        forward[i] += step
        backward[i] -= step
        jacobian[:, i] = (function(forward) - function(backward)) / (2 * step)
    return jacobian


# ==========================================================================
# A whole series
# ==========================================================================
@dataclass
class BandSeries:
    """Every line's fit, and the dispersions that come out of them."""
    positions: np.ndarray
    names: list
    centres: np.ndarray                  # (n_lines, n_bands), NaN where absent
    centre_errors: np.ndarray
    widths: np.ndarray
    width_errors: np.ndarray
    heights: np.ndarray
    areas: np.ndarray
    chi2: np.ndarray
    fits: list = field(default_factory=list)
    direction: str = "mdc"
    x_label: str = ""
    y_label: str = ""
    settings: FitSettings = None

    def band(self, index: int = 0):
        """``(k, E, k_error)`` for one band, in plotting order.

        For an MDC series the fitted quantity is the momentum and the
        position is the energy; for an EDC series it is the other way round.
        Which one carries the error decides how the dispersion must be
        fitted, so it is returned explicitly rather than left to the caller
        to remember.
        """
        good = np.isfinite(self.centres[:, index])
        if self.direction == "mdc":
            return (self.centres[good, index], self.positions[good],
                    self.centre_errors[good, index])
        return (self.positions[good], self.centres[good, index],
                self.centre_errors[good, index])

    def error_axis(self) -> str:
        """Which of the two axes the fitted error belongs to."""
        return "x" if self.direction == "mdc" else "y"


def fit_series(values, axes, bands, settings: FitSettings, *,
               direction: str = "mdc", step: int = 1, combine: int = 1,
               bounds=None, refine: bool = True, progress=None) -> BandSeries:
    """Fit every line of a cut, seeding each from the bands' own seeds.

    ``refine`` runs a second pass using the first pass's answer as the
    starting guess. The first pass follows the seeds, which are sparse; the
    second follows the fitted dispersion, which is dense. Doing only the
    second -- starting each line from its neighbour, as the previous
    implementation here did -- loses the band at a crossing and never
    recovers, because nothing afterwards pulls it back.
    """
    values = np.asarray(values, dtype=float)
    x = np.asarray(axes[0], dtype=float)
    y = np.asarray(axes[1], dtype=float)
    bands = [b for b in bands if b.seeds]
    if not bands:
        raise ValueError("place at least one seed on one band first")

    if direction == "mdc":
        along, across = x, y
        line_of = lambda i: values[:, i]
    elif direction == "edc":
        along, across = y, x
        line_of = lambda i: values[i, :]
    else:
        raise ValueError(f"direction {direction!r} is not mdc or edc")

    lo, hi = 0, across.size - 1
    if bounds is not None and bounds != (None, None):
        want_lo = -np.inf if bounds[0] is None else float(bounds[0])
        want_hi = np.inf if bounds[1] is None else float(bounds[1])
        keep = np.flatnonzero((across >= min(want_lo, want_hi)) &
                              (across <= max(want_lo, want_hi)))
        if keep.size == 0:
            raise ValueError("that range contains no lines")
        lo, hi = int(keep[0]), int(keep[-1])

    half = max(0, (int(combine) - 1) // 2)
    indices = list(range(lo, hi + 1, max(1, int(step))))
    tolerance = abs(across[1] - across[0]) if across.size > 1 else 0.0

    def one_pass(previous=None):
        fits = []
        for count, i in enumerate(indices):
            position = float(across[i])
            active = [b for b in bands if b.covers(position, tolerance)]
            if not active:
                fits.append(None)
                continue
            first = max(0, i - half)
            last = min(across.size, i + half + 1)
            line = np.nansum(np.array([line_of(j) for j in range(first, last)]),
                             axis=0)
            guesses = None
            if previous is not None and previous[count] is not None:
                old = previous[count]
                guesses = {name: (old.centres[j], old.widths[j], old.heights[j])
                           for j, name in enumerate(old.names)}
            try:
                fits.append(fit_line(along, line, active, settings,
                                     guesses=guesses, position=position))
            except Exception:
                fits.append(None)
            if progress is not None and not progress(count + 1, len(indices)):
                fits.extend([None] * (len(indices) - len(fits)))
                return fits, False
        return fits, True

    fits, finished = one_pass()
    if refine and finished:
        refined, _ = one_pass(fits)
        fits = [new if new is not None else old
                for new, old in zip(refined, fits)]

    names = [b.name for b in bands]
    blank = np.full((len(fits), len(names)), np.nan)
    centres, centre_errors = blank.copy(), blank.copy()
    widths, width_errors = blank.copy(), blank.copy()
    heights, areas = blank.copy(), blank.copy()
    chi2 = np.full(len(fits), np.nan)
    for row, fit in enumerate(fits):
        if fit is None:
            continue
        for j, name in enumerate(fit.names):
            column = names.index(name)
            centres[row, column] = fit.centres[j]
            centre_errors[row, column] = fit.centre_errors[j]
            widths[row, column] = fit.widths[j]
            width_errors[row, column] = fit.width_errors[j]
            heights[row, column] = fit.heights[j]
            areas[row, column] = fit.areas[j]
        chi2[row] = fit.chi2

    return BandSeries(np.array([float(across[i]) for i in indices]), names,
                      centres, centre_errors, widths, width_errors, heights,
                      areas, chi2, fits, direction, settings=settings)


# ==========================================================================
# Seeding helpers
# ==========================================================================
def suggest_seed(x, y, centre: float, *, window: float = None):
    """Height and width for a seed dropped at ``centre``.

    The user points at the band; the width and height are read off the data
    rather than asked for, since a click carries only a position.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    span = float(abs(x[-1] - x[0])) or 1.0
    window = float(window) if window else span / 8.0
    near = np.abs(x - float(centre)) <= window
    if near.sum() < 3:
        near = np.ones(x.shape, dtype=bool)
    local = y[near]
    floor = float(np.nanmin(local))
    peak = float(np.nanmax(local))
    height = max(peak - floor, 1e-9)
    above = x[near][y[near] - floor >= height / 2.0]
    width = float(above.max() - above.min()) if above.size > 1 else span / 20.0
    return height, max(width, span / 200.0)
