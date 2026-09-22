"""
tools/fermi.py
=============
Fermi-edge fitting, as plain numpy/scipy with no Qt -- the same split as
``tools.kspace`` and ``tools.analysis``.

Model (all in eV, energies as the axis carries them):

    I(E) = [ (a0 + a1*(E-EF)) * f(E; EF, T) + (b0 + b1*(E-EF)) ] (*) G(FWHM)

    f(E) = 1 / (exp((E-EF)/kT) + 1)

so the occupied density of states is a straight line through EF, the
background (dark counts, secondaries, higher-order light) is a straight line
that is *not* multiplied by the occupation, and the whole thing is smeared
by the analyser's Gaussian resolution.

What is different from the lab's MATLAB ``FitFermiSurface``, and why:

- **The convolution is a quadrature, not a padded grid.** The MATLAB code
  builds an extended uniform grid, convolves and slices the middle out,
  which needs the axis to be ascending and uniformly spaced, gets its own
  step wrong by one sample (``range/N`` instead of ``range/(N-1)``) and can
  slip a bin when the colon operator rounds. Integrating
  ``m(E - t) G(t) dt`` at each requested energy has none of those failure
  modes and works on any axis, ascending or not.
  (Gauss-Hermite quadrature would be the elegant version and is what this
  module first tried, but it under-resolves the edge whenever kT is small
  compared with the resolution -- exactly the nano-ARPES case -- so the
  integral is taken on a fine grid whose step follows the *narrower* of the
  two widths.)
- **Poisson weighting.** Counts have variance ~ counts, so an unweighted
  least squares over-weights the occupied side. This biases the width, and
  the width is the whole point.
- **Uncertainties.** ``least_squares`` gives a Jacobian, so every parameter
  comes back with an error bar and the fit with a reduced chi-squared. The
  MATLAB tool reports no uncertainty at all.
- **T and the resolution are degenerate** on a single edge: the measured
  width is roughly sqrt((3.53 kT)^2 + FWHM^2), so fitting both freely gives
  two numbers neither of which means anything on its own. The temperature
  is therefore *held* by default (the cryostat reading is in the file), and
  when both are free the fit reports their correlation so the degeneracy is
  visible rather than silent.
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass, field

import numpy as np

try:
    from scipy.optimize import least_squares
    from scipy.special import expit
except ImportError as exc:      # pragma: no cover - scipy is a hard dependency
    raise ImportError("tools.fermi needs scipy") from exc


#: Boltzmann's constant in eV/K.
K_B = 8.617333262e-5
#: FWHM = SIGMA_TO_FWHM * sigma for a Gaussian.
SIGMA_TO_FWHM = 2.0 * np.sqrt(2.0 * np.log(2.0))
#: The parameters of the model, in the order the fitter uses.
PARAMETERS = ("ef", "temperature", "resolution", "dos0", "dos1", "bkg0", "bkg1")
#: How each one is labelled in the interface, with its unit.
PARAMETER_LABELS = {
    "ef": "E_F (eV)",
    "temperature": "Temperature (K)",
    "resolution": "Resolution FWHM (eV)",
    "dos0": "DOS at E_F",
    "dos1": "DOS slope (/eV)",
    "bkg0": "Background",
    "bkg1": "Background slope (/eV)",
}


# --------------------------------------------------------------------------
# The model
# --------------------------------------------------------------------------
def fermi_function(energy, ef: float, temperature: float) -> np.ndarray:
    """Fermi-Dirac occupation. ``expit(-x)`` rather than ``1/(exp(x)+1)``:
    same function, but it does not overflow when the exponent is large."""
    energy = np.asarray(energy, dtype=float)
    kt = K_B * max(float(temperature), 1e-6)
    return expit(-(energy - float(ef)) / kt)


def _unbroadened(energy, ef, temperature, dos0, dos1, bkg0, bkg1):
    shifted = np.asarray(energy, dtype=float) - float(ef)
    return ((dos0 + dos1 * shifted) * fermi_function(energy, ef, temperature)
            + bkg0 + bkg1 * shifted)


def _convolution_offsets(sigma: float, kt: float, n_sigma: float = 6.0,
                         max_points: int = 1201):
    """Offsets and weights for integrating against a Gaussian of width
    ``sigma``.

    The step follows the *narrower* of the two scales in the problem -- the
    resolution and the thermal width -- because a step tuned to the
    resolution alone walks straight over a cold Fermi edge.
    """
    scale = min(sigma, kt) if kt > 0 else sigma
    step = max(scale / 6.0, 2.0 * n_sigma * sigma / (max_points - 1))
    half = int(np.ceil(n_sigma * sigma / step))
    offsets = np.arange(-half, half + 1) * step
    weights = np.exp(-0.5 * (offsets / sigma) ** 2)
    return offsets, weights / weights.sum()


def fermi_edge_model(energy, ef: float, temperature: float, resolution: float,
                     dos0: float, dos1: float = 0.0, bkg0: float = 0.0,
                     bkg1: float = 0.0) -> np.ndarray:
    """The full model above, evaluated at ``energy`` (any order, any
    spacing). ``resolution`` is the Gaussian **FWHM**."""
    energy = np.asarray(energy, dtype=float)
    sigma = float(resolution) / SIGMA_TO_FWHM
    if not np.isfinite(sigma) or sigma <= 0:
        return _unbroadened(energy, ef, temperature, dos0, dos1, bkg0, bkg1)
    offsets, weights = _convolution_offsets(sigma, K_B * max(temperature, 1e-6))
    grid = energy[:, None] - offsets[None, :]
    values = _unbroadened(grid, ef, temperature, dos0, dos1, bkg0, bkg1)
    return values @ weights


def broadened_fermi(energy, ef: float, temperature: float,
                    resolution: float) -> np.ndarray:
    """Just the occupation, resolution-broadened and normalised to 1 well
    below E_F. This -- not the whole model -- is what dividing a spectrum by
    the Fermi cut-off should divide by."""
    return fermi_edge_model(energy, ef, temperature, resolution,
                            dos0=1.0, dos1=0.0, bkg0=0.0, bkg1=0.0)


# --------------------------------------------------------------------------
# Starting values
# --------------------------------------------------------------------------
def _smooth(values: np.ndarray, width: int) -> np.ndarray:
    width = max(3, int(width) | 1)
    if values.size < width:
        return values
    kernel = np.ones(width) / width
    return np.convolve(values, kernel, mode="same")


def initial_guess(energy, intensity, temperature: float = 30.0) -> dict:
    """Starting values read off the data itself, so nothing has to be typed.

    E_F is the steepest fall of the smoothed spectrum; the background and
    the density of states are the levels well above and well below it; the
    resolution comes from the 16-84 % width of the edge with the thermal
    part taken out in quadrature (an edge is 3.32 kT wide between those two
    levels, a Gaussian 2 sigma, and the two add in quadrature).
    """
    energy = np.asarray(energy, dtype=float)
    intensity = np.asarray(intensity, dtype=float)
    order = np.argsort(energy)
    e, i = energy[order], intensity[order]
    finite = np.isfinite(i)
    if finite.sum() < 5:
        raise ValueError("not enough finite points to guess from")
    e, i = e[finite], i[finite]

    smoothed = _smooth(i, max(3, e.size // 50))
    slope = np.gradient(smoothed, e)
    ef = float(e[int(np.argmin(slope))])

    kt = K_B * max(float(temperature), 1e-6)
    above = e > ef + 10 * kt
    below = e < ef - 10 * kt
    bkg0 = float(np.median(i[above])) if above.sum() > 2 else float(np.min(i))
    level = float(np.median(i[below])) if below.sum() > 2 else float(np.max(i))
    dos0 = max(level - bkg0, 1e-12)

    # 16-84 % width of the edge, minus the thermal part
    lo_target, hi_target = bkg0 + 0.16 * dos0, bkg0 + 0.84 * dos0
    crossings = []
    for target in (hi_target, lo_target):
        idx = np.argmin(np.abs(smoothed - target))
        crossings.append(float(e[idx]))
    width = abs(crossings[1] - crossings[0])
    thermal = 3.32 * kt
    sigma = np.sqrt(max(width ** 2 - thermal ** 2, 0.0)) / 2.0
    resolution = max(sigma * SIGMA_TO_FWHM, float(np.mean(np.diff(e))) if e.size > 1 else 1e-3)

    dos1 = 0.0
    if below.sum() > 3:                       # slope of the occupied states
        dos1 = float(np.polyfit(e[below] - ef, i[below] - bkg0, 1)[0])
    bkg1 = 0.0
    if above.sum() > 3:
        bkg1 = float(np.polyfit(e[above] - ef, i[above], 1)[0])

    return {"ef": ef, "temperature": float(temperature), "resolution": float(resolution),
            "dos0": dos0, "dos1": dos1, "bkg0": bkg0, "bkg1": bkg1}


# --------------------------------------------------------------------------
# The fit
# --------------------------------------------------------------------------
@dataclass
class FermiFit:
    """What a fit produced: values, their uncertainties, and enough about
    the fit itself to judge whether to believe them."""
    values: dict
    errors: dict
    fixed: tuple = ()
    reduced_chi2: float = float("nan")
    correlation: dict = field(default_factory=dict)
    window: tuple = (None, None)
    n_points: int = 0
    success: bool = False
    message: str = ""

    def model(self, energy) -> np.ndarray:
        return fermi_edge_model(energy, **{k: self.values[k] for k in PARAMETERS})

    @property
    def thermal_width(self) -> float:
        """The 10-90 % width the temperature alone would give (eV)."""
        return 3.53 * K_B * self.values["temperature"]

    @property
    def combined_width(self) -> float:
        """Thermal and instrumental widths added in quadrature -- the width
        actually measured, and the only one determined when both are free."""
        return float(np.hypot(self.thermal_width, self.values["resolution"]))

    def summary(self) -> str:
        ef, def_ = self.values["ef"], self.errors.get("ef", float("nan"))
        res, dres = self.values["resolution"], self.errors.get("resolution", float("nan"))
        text = (f"E_F = {ef:.5g} ± {def_:.2g} eV\n"
                f"resolution (FWHM) = {res * 1000:.4g} ± {dres * 1000:.2g} meV\n"
                f"T = {self.values['temperature']:.4g} K"
                f"{' (held)' if 'temperature' in self.fixed else ''}, "
                f"thermal width {self.thermal_width * 1000:.3g} meV\n"
                f"reduced chi2 = {self.reduced_chi2:.4g} over {self.n_points} points")
        pair = self.correlation.get(("temperature", "resolution"))
        if pair is not None:
            text += (f"\ncorrelation(T, resolution) = {pair:+.2f}"
                     + ("  -- they are trading off; hold one of them"
                        if abs(pair) > 0.9 else ""))
        return text


def _weights(intensity: np.ndarray, weighting: str) -> np.ndarray:
    if weighting == "uniform":
        return np.ones_like(intensity)
    # Poisson: sigma = sqrt(counts), with a floor so empty channels do not
    # get infinite weight.
    floor = max(1.0, float(np.nanmax(intensity)) * 1e-6)
    return 1.0 / np.sqrt(np.maximum(intensity, floor))


def fit_fermi_edge(energy, intensity, *, start: dict = None,
                   fixed=("temperature",), temperature: float = 30.0,
                   window=(None, None), weighting: str = "poisson",
                   max_nfev: int = 2000) -> FermiFit:
    """Fit the model to one EDC.

    ``fixed`` names the parameters to hold at their starting values; the
    temperature is held by default, because on a single edge it cannot be
    separated from the resolution. ``window`` restricts the fit to an energy
    range, which keeps band structure far from E_F out of it.
    """
    energy = np.asarray(energy, dtype=float)
    intensity = np.asarray(intensity, dtype=float)
    if energy.shape != intensity.shape:
        raise ValueError("energy and intensity must have the same shape")

    keep = np.isfinite(energy) & np.isfinite(intensity)
    lo, hi = window
    if lo is not None:
        keep &= energy >= float(lo)
    if hi is not None:
        keep &= energy <= float(hi)
    e, i = energy[keep], intensity[keep]
    if e.size < 8:
        raise ValueError(f"only {e.size} points inside the fit window")

    start = dict(start or initial_guess(e, i, temperature=temperature))
    start.setdefault("temperature", temperature)
    fixed = tuple(name for name in fixed if name in PARAMETERS)
    free = [name for name in PARAMETERS if name not in fixed]
    if not free:
        raise ValueError("every parameter is held -- nothing to fit")

    span = float(np.ptp(e)) or 1.0
    step = float(np.mean(np.abs(np.diff(np.sort(e))))) if e.size > 1 else span / 100
    scale = float(np.nanmax(np.abs(i))) or 1.0
    limits = {
        "ef": (float(e.min()), float(e.max())),
        "temperature": (0.1, 2000.0),
        # A resolution finer than a fraction of the sampling cannot be seen,
        # and zero would divide by zero in the kernel.
        "resolution": (step / 10.0, span),
        "dos0": (-10 * scale, 10 * scale),
        "dos1": (-1e4 * scale / span, 1e4 * scale / span),
        "bkg0": (-10 * scale, 10 * scale),
        "bkg1": (-1e4 * scale / span, 1e4 * scale / span),
    }
    x0 = np.array([np.clip(start[name], *limits[name]) for name in free], dtype=float)
    lower = np.array([limits[name][0] for name in free])
    upper = np.array([limits[name][1] for name in free])
    weights = _weights(i, weighting)

    def residuals(x):
        values = dict(start)
        values.update(dict(zip(free, x)))
        model = fermi_edge_model(e, **{k: values[k] for k in PARAMETERS})
        return (model - i) * weights

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = least_squares(residuals, x0, bounds=(lower, upper),
                               x_scale="jac", max_nfev=max_nfev)

    values = dict(start)
    values.update(dict(zip(free, result.x)))

    # Uncertainties from the Jacobian, scaled by the fit's own scatter (the
    # convention curve_fit uses by default): trustworthy relative to each
    # other whether or not the absolute weights are right.
    dof = max(e.size - len(free), 1)
    chi2 = float(2 * result.cost)
    reduced = chi2 / dof
    errors = {name: float("nan") for name in PARAMETERS}
    correlation = {}
    try:
        _, s, vt = np.linalg.svd(result.jac, full_matrices=False)
        good = s > max(result.jac.shape) * np.finfo(float).eps * s[0]
        cov = (vt[good].T / s[good] ** 2) @ vt[good] * reduced
        for index, name in enumerate(free):
            errors[name] = float(np.sqrt(abs(cov[index, index])))
        for a in range(len(free)):
            for b in range(a + 1, len(free)):
                denom = np.sqrt(abs(cov[a, a] * cov[b, b]))
                if denom > 0:
                    correlation[(free[a], free[b])] = float(cov[a, b] / denom)
    except (np.linalg.LinAlgError, IndexError, ValueError):
        pass
    for name in fixed:
        errors[name] = 0.0

    return FermiFit(values=values, errors=errors, fixed=fixed,
                    reduced_chi2=reduced, correlation=correlation,
                    window=(float(e.min()), float(e.max())), n_points=int(e.size),
                    success=bool(result.success), message=str(result.message))


# --------------------------------------------------------------------------
# Channel by channel, for a Fermi-surface correction from a reference
# --------------------------------------------------------------------------
def fit_channels(angles, energy, frame, *, half_width: int = 0, step: int = 1,
                 temperature: float = 30.0, fixed=("temperature",),
                 window=(None, None), weighting: str = "poisson",
                 start: dict = None, progress=None):
    """Fit the edge in every detector channel of ``frame`` (angle, energy).

    ``half_width`` adds that many neighbouring channels on each side to the
    EDC before fitting -- the usual answer to a reference spectrum that is
    too thin to fit channel by channel. The channel's own position is
    unchanged; only its statistics improve.

    ``step`` fits only every Nth channel. The curvature is smooth and a
    thousand fits take minutes, so sampling it and letting the polynomial
    through the points do the interpolation costs nothing real.

    Returns ``(ef, ef_error, ok)``, all the length of ``angles``: the edge
    position per channel, its uncertainty, and whether that channel's fit is
    worth using (converged, E_F inside the window, and an uncertainty
    smaller than the resolution). Channels that were skipped or failed are
    NaN with ``ok`` false.
    """
    angles = np.asarray(angles, dtype=float)
    frame = np.asarray(frame, dtype=float)
    energy = np.asarray(energy, dtype=float)
    if frame.shape != (angles.size, energy.size):
        raise ValueError(f"frame {frame.shape} does not match "
                         f"({angles.size}, {energy.size})")

    half_width = max(0, int(half_width))
    ef = np.full(angles.size, np.nan)
    err = np.full(angles.size, np.nan)
    ok = np.zeros(angles.size, dtype=bool)
    guess = dict(start) if start else None

    wanted = list(range(0, angles.size, max(1, int(step))))
    if wanted and wanted[-1] != angles.size - 1:
        wanted.append(angles.size - 1)      # always anchor both ends
    for count, index in enumerate(wanted):
        lo = max(0, index - half_width)
        hi = min(angles.size, index + half_width + 1)
        edc = np.nansum(frame[lo:hi], axis=0)
        try:
            fit = fit_fermi_edge(energy, edc, start=guess, fixed=fixed,
                                 temperature=temperature, window=window,
                                 weighting=weighting)
        except Exception:
            if progress is not None and not progress(count + 1, len(wanted)):
                break
            continue
        ef[index] = fit.values["ef"]
        err[index] = fit.errors.get("ef", np.nan)
        inside = float(np.min(energy)) < ef[index] < float(np.max(energy))
        ok[index] = bool(fit.success and inside
                         and (not np.isfinite(err[index])
                              or err[index] < max(fit.values["resolution"], 1e-9)))
        if ok[index]:
            # Walk along the detector: the next channel starts from this
            # one's answer, which is both faster and steadier than starting
            # every channel from scratch.
            guess = dict(fit.values)
        if progress is not None and not progress(count + 1, len(wanted)):
            break                           # the caller asked us to stop
    return ef, err, ok


# --------------------------------------------------------------------------
# Dividing the Fermi cut-off out
# --------------------------------------------------------------------------
def divide_fermi(energy, values, ef: float, temperature: float,
                 resolution: float, *, energy_dim: int = -1,
                 cutoff_kt: float = 4.0, background=None):
    """Divide a spectrum by the resolution-broadened Fermi function.

    Only the **occupation** is divided out -- not the fitted density of
    states, which the MATLAB version includes and so removes real spectral
    weight along with the cut-off.

    ``background`` is the fitted ``(b0, b1)``, and if it is given it is taken
    out before the division and put back afterwards. That matters: dark
    counts and secondaries are not modulated by the occupation, so dividing
    them by it inflates them near E_F by exactly the factor the division is
    meant to remove from the signal.

    Everything above ``ef + cutoff_kt * kT`` is set to NaN: the divisor there
    is a small number known only to within the noise, and dividing by it
    produces spikes, not spectra. The viewer draws NaN as missing, which is
    what it is.
    """
    energy = np.asarray(energy, dtype=float)
    values = np.asarray(values, dtype=float)
    divisor = broadened_fermi(energy, ef, temperature, resolution)
    cut = float(ef) + float(cutoff_kt) * K_B * max(float(temperature), 1e-6)
    divisor = np.where(energy <= cut, divisor, np.nan)
    shape = [1] * values.ndim
    shape[energy_dim] = energy.size
    with np.errstate(invalid="ignore", divide="ignore"):
        divisor = np.where(np.abs(divisor) > 1e-6, divisor, np.nan).reshape(shape)
        if background is None:
            return values / divisor
        b0, b1 = (float(v) for v in background)
        bkg = (b0 + b1 * (energy - float(ef))).reshape(shape)
        return (values - bkg) / divisor + bkg
