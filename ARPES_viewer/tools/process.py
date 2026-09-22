"""
tools/process.py
==============
The 2-D data operations of the lab's MATLAB tools, rewritten as plain numpy
functions with no Qt: smoothing, derivatives, curvature, symmetrisation,
background subtraction, normalisation and despiking.

MDC/EDC peak fitting, band dispersion and self-energy extraction used to
live here too; they now have their own modules, :mod:`tools.peaks` and
:mod:`tools.dispersion`, because the fitting grew a Fermi-Dirac occupation
model and a variable-projection solver that have nothing to do with the
image operations below.

Everything here takes ``(values, axes)`` and returns ``(values, axes)`` or a
result object, so the same function serves the interactive panel, a script
and the tests. What the axes *mean* is the caller's business -- these
functions only need to know their steps.

Three deliberate departures from the MATLAB originals
-----------------------------------------------------
**Curvature is dimensionless.** ``Curvature.m`` feeds ``I_x`` (counts per
A^-1) and ``I_y`` (counts per eV) into the same ``(a0 + I_x^2 + I_y^2)``.
Those differ by one to two orders of magnitude, so its "isotropic" curvature
is nothing of the sort and the best ``a0`` moves with every change of axis
range or intensity scale. Here both axes are mapped to [0, 1] and the
intensity to [0, 1] before the derivatives are taken, so ``a0`` is a pure
number that means the same thing on the next dataset. :func:`curvature`
still takes an explicit ``ratio`` for the cases where the anisotropy is
wanted on purpose.

**Derivatives are Savitzky-Golay.** ``smooth_derivation_v1.m`` takes nested
one-sided ``diff``s -- which offsets each order by half a pixel and offsets
the cross terms by more -- and then forces the borders to zero, putting a
false zero frame on the picture and poisoning the auto-levels. A local
polynomial fit smooths and differentiates in the same step, keeps every
stencil centred, and handles the borders by fitting the truncated window
instead of inventing a value.

**No histogram clipping.** ``Curvature.m`` ends by multiplying its result
by a mask built from a 100-bin histogram, which silently sets every negative
curvature to zero -- half the structure, and the half that says where the
band is *not*. Clipping is a display decision; it belongs in the levels, and
the panel offers it there.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import ndimage, signal

#: Boltzmann's constant in eV/K -- the same value tools.fermi uses.
K_B = 8.617333262e-5


# ==========================================================================
# Shared helpers
# ==========================================================================
def axis_step(axis) -> float:
    """The mean step of an axis, positive, never zero.

    Axes here are regular by construction (the file parser builds them from
    start/stop/n), so the mean is the step; taking it from the ends rather
    than from ``diff`` keeps it exact for a descending axis too.
    """
    axis = np.asarray(axis, dtype=float).ravel()
    if axis.size < 2:
        return 1.0
    step = abs(float(axis[-1] - axis[0])) / (axis.size - 1)
    return step if step > 0 else 1.0


def _as_2d(values):
    values = np.asarray(values, dtype=float)
    if values.ndim != 2:
        raise ValueError(f"this operation works on 2-D data, not {values.ndim}-D")
    return values


def _fill_nans(values: np.ndarray) -> np.ndarray:
    """Replace every NaN with the value of its nearest finite neighbour.

    Filters need a value everywhere. Filling with zero would drag a band's
    edge down towards zero and put a step where the data merely stops;
    filling from the nearest finite point keeps the local level, and the
    caller puts the NaNs back afterwards, so nothing is invented in the
    output -- only in the window that produced it.
    """
    finite = np.isfinite(values)
    if finite.all():
        return values
    if not finite.any():
        return np.zeros_like(values)
    index = ndimage.distance_transform_edt(~finite, return_distances=False,
                                           return_indices=True)
    return values[tuple(index)]


def _restore_nans(result: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """Put the original missing points back into a filtered result."""
    if mask.all():
        return result
    out = np.asarray(result, dtype=float).copy()
    out[~mask] = np.nan
    return out


def _odd(n: int, lo: int = 3) -> int:
    """The nearest odd integer at least ``lo`` -- window lengths must be odd
    for the fit to be centred on the point it replaces."""
    n = int(round(n))
    if n % 2 == 0:
        n += 1
    return max(lo, n)


def normalise_unit(values: np.ndarray):
    """Map finite values to [0, 1]; returns ``(scaled, lo, span)``.

    Curvature and the gradient enhancement compare derivatives against a
    constant, so the intensity has to have a known scale before that
    constant means anything.
    """
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return values, 0.0, 1.0
    lo = float(finite.min())
    span = float(finite.max()) - lo
    if span <= 0:
        return np.zeros_like(values), lo, 1.0
    return (values - lo) / span, lo, span


# ==========================================================================
# 1. Smoothing
# ==========================================================================
def gaussian_smooth(values, axes, sigmas, *, in_units: bool = True):
    """NaN-aware Gaussian blur, with the widths given in the axes' own units.

    ``sigmas`` is one width per dimension -- in eV, A^-1, degrees, whatever
    the axis is in, when ``in_units``; in pixels otherwise. A width in
    physical units is the point: a resolution of 15 meV means the same blur
    whether the cut was taken on a 200-point or a 1000-point energy axis,
    which ``csaps``'s ``p = 1 - 10^-k`` never could.

    Missing points do not leak: the data and a 0/1 mask are blurred with the
    same kernel and divided, so a point beside a gap is the weighted mean of
    the neighbours that exist rather than being dragged towards zero. The
    gap itself stays NaN.
    """
    values = np.asarray(values, dtype=float)
    if len(axes) != values.ndim or len(sigmas) != values.ndim:
        raise ValueError(f"{values.ndim}-D data needs {values.ndim} axes and widths")
    pixels = []
    for dim, (axis, sigma) in enumerate(zip(axes, sigmas)):
        sigma = float(sigma)
        pixels.append(sigma / axis_step(axis) if in_units else sigma)
    if all(p <= 0 for p in pixels):
        return values.copy()

    mask = np.isfinite(values)
    filled = np.where(mask, values, 0.0)
    num = ndimage.gaussian_filter(filled, pixels, mode="nearest")
    den = ndimage.gaussian_filter(mask.astype(float), pixels, mode="nearest")
    with np.errstate(invalid="ignore", divide="ignore"):
        out = np.where(den > 1e-9, num / den, np.nan)
    return _restore_nans(out, mask)


def savgol_smooth(values, axes, windows, order: int = 2, *, in_units: bool = True):
    """Savitzky-Golay smoothing, applied along each axis in turn.

    Unlike a Gaussian this keeps peak height and width: it fits a local
    polynomial rather than averaging, so an MDC that is about to be fitted
    is not made artificially wider by the smoothing that preceded it.
    ``windows`` is the full width of the fitting window per axis, in the
    axis's units when ``in_units``.
    """
    values = np.asarray(values, dtype=float)
    if len(axes) != values.ndim or len(windows) != values.ndim:
        raise ValueError(f"{values.ndim}-D data needs {values.ndim} axes and windows")
    mask = np.isfinite(values)
    out = _fill_nans(values)
    for dim, (axis, window) in enumerate(zip(axes, windows)):
        width = float(window)
        if width <= 0:
            continue
        n = _odd(width / axis_step(axis) if in_units else width, lo=order + 2)
        n = min(n, _odd(out.shape[dim], lo=3))
        if n <= order + 1:
            continue
        out = signal.savgol_filter(out, n, order, axis=dim, mode="interp")
    return _restore_nans(out, mask)


def box_smooth(values, axes, windows, *, in_units: bool = True):
    """Plain moving average -- the honest version of a detector rebin that
    keeps the sampling. NaN-aware in the same way as :func:`gaussian_smooth`."""
    values = np.asarray(values, dtype=float)
    sizes = []
    for axis, window in zip(axes, windows):
        width = float(window)
        n = _odd(width / axis_step(axis) if in_units else width, lo=1)
        sizes.append(max(1, n))
    if all(s <= 1 for s in sizes):
        return values.copy()
    mask = np.isfinite(values)
    filled = np.where(mask, values, 0.0)
    num = ndimage.uniform_filter(filled, sizes, mode="nearest")
    den = ndimage.uniform_filter(mask.astype(float), sizes, mode="nearest")
    with np.errstate(invalid="ignore", divide="ignore"):
        out = np.where(den > 1e-9, num / den, np.nan)
    return _restore_nans(out, mask)


SMOOTHERS = {"gaussian": gaussian_smooth, "savgol": savgol_smooth, "box": box_smooth}


# ==========================================================================
# 2. Derivatives
# ==========================================================================
def derivative(values, axes, dim: int, order: int = 2, *, window=None,
               poly: int = None, in_units: bool = True, negate: bool = True):
    """The ``order``-th derivative along one axis, by local polynomial fit.

    ``window`` is the fitting width in the axis's units (its default is
    seven points, the narrowest that supports a stable second derivative).
    ``negate`` returns ``-d^n/dx^n``, which is the convention every ARPES
    paper prints: a peak then shows as a peak.

    ``poly`` defaults to ``order + 1`` rather than to ``order``, and that
    matters only at the border -- but it matters a lot there. The border is
    handled by fitting the truncated window and evaluating it, so with
    ``poly == order`` the derivative is a single constant across the whole
    edge window and cannot follow the curve: on a test function the error
    is 9% at the edge against 0.07% inside. One more degree makes it 0.07%
    everywhere. Either way it beats forcing the border to zero the way
    ``smooth_derivation_v1.m`` does -- a zero frame is a feature that is not
    in the data, and it drags the auto-levels with it.
    """
    values = np.asarray(values, dtype=float)
    dim = int(dim) % values.ndim
    step = axis_step(axes[dim])
    if poly is None:
        poly = order + 1
    if window is None:
        n = 7
    else:
        n = _odd(float(window) / step if in_units else float(window), lo=order + 2)
    n = max(n, order + 2 if (order + 2) % 2 else order + 3)
    n = min(_odd(n), _odd(values.shape[dim]))
    if n <= poly:
        poly = max(1, n - 1)
    poly = max(poly, order)

    mask = np.isfinite(values)
    filled = _fill_nans(values)
    out = signal.savgol_filter(filled, n, poly, deriv=order, delta=step,
                               axis=dim, mode="interp")
    if negate:
        out = -out
    return _restore_nans(out, mask)


def laplacian(values, axes, *, window=None, poly: int = None, in_units: bool = True,
              weights=None, negate: bool = True):
    """``-(d^2/dx^2 + d^2/dy^2)``, each term normalised before it is added.

    ``del2`` in MATLAB silently carries a factor of 1/4 and adds the two
    second derivatives in their own units, so whichever axis happens to be
    in the smaller unit dominates. Here each term is divided by its own
    RMS first (or by ``weights``, if the caller wants a specific balance),
    which is what makes the sum mean "curvature in both directions" rather
    than "curvature along whichever axis has the coarser step".
    """
    values = _as_2d(values)
    dxx = derivative(values, axes, 0, 2, window=window, poly=poly,
                     in_units=in_units, negate=False)
    dyy = derivative(values, axes, 1, 2, window=window, poly=poly,
                     in_units=in_units, negate=False)
    if weights is None:
        wx = _rms(dxx)
        wy = _rms(dyy)
        dxx = dxx / wx
        dyy = dyy / wy
    else:
        dxx = dxx * float(weights[0])
        dyy = dyy * float(weights[1])
    out = dxx + dyy
    return -out if negate else out


def _rms(values: np.ndarray) -> float:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return 1.0
    value = float(np.sqrt(np.mean(finite ** 2)))
    return value if value > 0 else 1.0


# ==========================================================================
# 3. Curvature  (Zhang et al., Rev. Sci. Instrum. 82, 043712 (2011))
# ==========================================================================
def curvature(values, axes, *, mode: str = "2d", a0: float = 1.0,
              ratio: float = 1.0, window=None, poly: int = None,
              normalise: bool = True, negate: bool = True):
    """Curvature of the intensity surface, with a dimensionless ``a0``.

    ``mode`` is ``"2d"`` for the isotropic two-dimensional form (Fermi
    surfaces, constant-energy maps), ``"x"`` or ``"y"`` for the
    one-dimensional form along that axis (dispersions -- ``"y"`` along
    energy sharpens EDC peaks, ``"x"`` along momentum sharpens MDC peaks).

    With ``normalise`` both axes are mapped to [0, 1] and the intensity to
    [0, 1] first. That is the whole difference from the MATLAB version:
    ``a0`` becomes a pure number, so the value that worked yesterday works
    today, and the 2-D form really is isotropic. ``ratio`` reintroduces a
    deliberate anisotropy -- ``C_y = C_x / ratio^2`` -- for the cases where
    one direction should be weighted more.

    ``a0`` is Zhang et al.'s ``C0``, and it interpolates between the two
    things people actually mean by "sharpening". Where ``a0`` is much larger
    than the squared slope the denominator is effectively constant and the
    result is the plain second derivative, scaled: sharp, and as noisy as a
    second derivative always is. As ``a0`` falls towards the squared slope
    the denominator starts to bite and the result approaches the true
    geometric curvature of the intensity surface, which is damped wherever
    the surface is steep. :func:`suggest_a0` returns the scale at which that
    crossover happens, which is the useful place to start.

    Returned so that a **peak is positive**, which is how it is printed.
    """
    values = _as_2d(values)
    x_axis = np.asarray(axes[0], dtype=float)
    y_axis = np.asarray(axes[1], dtype=float)

    if normalise:
        work, _, _ = normalise_unit(values)
        span_x = abs(float(x_axis[-1] - x_axis[0])) or 1.0
        span_y = abs(float(y_axis[-1] - y_axis[0])) or 1.0
        unit_x = np.linspace(0.0, 1.0, x_axis.size)
        unit_y = np.linspace(0.0, 1.0, y_axis.size)
        work_axes = [unit_x, unit_y]
        # A window given in real units has to follow the axes into [0, 1].
        if window is not None:
            window = (float(window[0]) / span_x, float(window[1]) / span_y)
    else:
        work = values
        work_axes = [x_axis, y_axis]

    wx = None if window is None else window[0]
    wy = None if window is None else window[1]

    a0 = max(float(a0), 1e-12)
    cx = 1.0 / a0
    cy = cx / max(float(ratio), 1e-12) ** 2

    ix = derivative(work, work_axes, 0, 1, window=wx, poly=poly, negate=False)
    iy = derivative(work, work_axes, 1, 1, window=wy, poly=poly, negate=False)
    ixx = derivative(work, work_axes, 0, 2, window=wx, poly=poly, negate=False)
    iyy = derivative(work, work_axes, 1, 2, window=wy, poly=poly, negate=False)

    with np.errstate(invalid="ignore", divide="ignore"):
        if mode == "x":
            out = cx * ixx / np.power(1.0 + cx * ix ** 2, 1.5)
        elif mode == "y":
            out = cy * iyy / np.power(1.0 + cy * iy ** 2, 1.5)
        elif mode == "2d":
            ixy = derivative(ix, work_axes, 1, 1, window=wy, poly=poly, negate=False)
            num = ((1.0 + cx * ix ** 2) * cy * iyy
                   - 2.0 * cx * cy * ix * iy * ixy
                   + (1.0 + cy * iy ** 2) * cx * ixx)
            out = num / np.power(1.0 + cx * ix ** 2 + cy * iy ** 2, 1.5)
        else:
            raise ValueError(f"curvature mode {mode!r} is not one of 2d, x, y")
    return -out if negate else out


def suggest_a0(values, axes, *, mode: str = "2d", factor: float = 1.0) -> float:
    """A starting ``a0``: the mean square slope of the normalised surface.

    Zhang et al. recommend choosing ``a0`` near the scale of ``I'^2`` -- the
    point where the denominator starts to matter at all. Because the axes
    and the intensity are already dimensionless here, that recommendation
    turns into a single number the panel can offer as a default, which is
    what ``Curvature.m``'s log slider was really searching for by hand.
    """
    values = _as_2d(values)
    work, _, _ = normalise_unit(values)
    unit = [np.linspace(0.0, 1.0, values.shape[0]),
            np.linspace(0.0, 1.0, values.shape[1])]
    ix = derivative(work, unit, 0, 1, negate=False)
    iy = derivative(work, unit, 1, 1, negate=False)
    if mode == "x":
        scale = _rms(ix) ** 2
    elif mode == "y":
        scale = _rms(iy) ** 2
    else:
        scale = _rms(ix) ** 2 + _rms(iy) ** 2
    return float(max(scale, 1e-9) * float(factor))


def gradient_enhance(values, axes, *, floor: float = 0.05, window=None,
                     normalise: bool = True):
    """``I / |grad I|``, the contrast trick from ``Gradient.m``, made safe.

    Two repairs. ``Gradient.m`` sums the squares of **eight** neighbour
    differences, mixing the diagonals in at the same weight as the axial
    ones, which is not the modulus of a gradient and is not isotropic; this
    uses ``sqrt(I_x^2 + I_y^2)`` on the normalised axes. And it divides by a
    quantity that goes to zero wherever the image is flat, so the flattest
    -- usually emptiest -- parts of the picture come out brightest; here the
    divisor is floored at ``floor`` times its own RMS.

    It stays a display trick with no absolute meaning, which is why the
    panel offers curvature first.
    """
    values = _as_2d(values)
    if normalise:
        work, _, _ = normalise_unit(values)
        work_axes = [np.linspace(0.0, 1.0, values.shape[0]),
                     np.linspace(0.0, 1.0, values.shape[1])]
    else:
        work, work_axes = values, [np.asarray(a, float) for a in axes]
    wx = None if window is None else window[0]
    wy = None if window is None else window[1]
    ix = derivative(work, work_axes, 0, 1, window=wx, negate=False)
    iy = derivative(work, work_axes, 1, 1, window=wy, negate=False)
    grad = np.sqrt(ix ** 2 + iy ** 2)
    guard = max(float(floor), 1e-6) * _rms(grad)
    with np.errstate(invalid="ignore", divide="ignore"):
        return work / np.maximum(grad, guard)


# ==========================================================================
# 4. Symmetrisation
# ==========================================================================
@dataclass
class SymmetryResult:
    """A symmetrised map and the evidence for judging it."""
    values: np.ndarray
    x: np.ndarray
    y: np.ndarray
    coverage: np.ndarray            # how many copies contributed per pixel
    centre: tuple                   # the origin actually used
    residual: float = 0.0           # RMS(sym - original) / RMS(original)

    def difference(self, original) -> np.ndarray:
        """``symmetrised - original``, for the third panel of the three that
        should always be shown together."""
        return self.values - np.asarray(original, dtype=float)


def _bilinear(values, x_axis, y_axis, xq, yq):
    """Sample a regular grid at arbitrary points, NaN-aware.

    A regular grid means a map_coordinates call rather than a Delaunay
    triangulation: ``TriScatteredInterp`` in ``data_symmetrization_pro.m``
    rebuilds a triangulation of every pixel to interpolate a grid it already
    has, which is both slower and blurrier than sampling it directly.
    """
    values = np.asarray(values, dtype=float)
    x0 = float(x_axis[0])
    y0 = float(y_axis[0])
    dx = (float(x_axis[-1]) - x0) / max(x_axis.size - 1, 1)
    dy = (float(y_axis[-1]) - y0) / max(y_axis.size - 1, 1)
    if dx == 0:
        dx = 1.0
    if dy == 0:
        dy = 1.0
    ix = (np.asarray(xq, float) - x0) / dx
    iy = (np.asarray(yq, float) - y0) / dy

    mask = np.isfinite(values)
    filled = np.where(mask, values, 0.0)
    coords = np.array([ix.ravel(), iy.ravel()])
    num = ndimage.map_coordinates(filled, coords, order=1, mode="constant", cval=0.0)
    den = ndimage.map_coordinates(mask.astype(float), coords, order=1,
                                  mode="constant", cval=0.0)
    inside = ((ix >= 0) & (ix <= values.shape[0] - 1) &
              (iy >= 0) & (iy <= values.shape[1] - 1)).ravel()
    with np.errstate(invalid="ignore", divide="ignore"):
        out = np.where((den > 0.5) & inside, num / np.maximum(den, 1e-12), np.nan)
    return out.reshape(np.shape(xq))


def symmetry_operations(fold: int = 1, *, mirror_angles=(), inversion: bool = False):
    """The list of (angle, mirror) operations a symmetrisation will apply.

    ``fold`` is the rotation order (``6`` for a hexagonal surface),
    ``mirror_angles`` are the directions of any mirror lines in degrees, and
    ``inversion`` adds ``k -> -k``. ``data_symmetrization_pro.m`` only ever
    rotates; mirrors and inversion are the other half of a surface's point
    group and are what let a half-measured map be completed.
    """
    fold = max(1, int(fold))
    rotations = [360.0 * i / fold for i in range(fold)]
    ops = [(angle, None) for angle in rotations]
    for line in mirror_angles:
        for angle in rotations:
            ops.append((angle, float(line)))
    if inversion:
        ops.extend([(angle + 180.0, None) for angle in rotations])
    return ops


def symmetrise(values, x_axis, y_axis, *, fold: int = 1, centre=(0.0, 0.0),
               mirror_angles=(), inversion: bool = False,
               sector=None) -> SymmetryResult:
    """Average a map over its symmetry operations.

    The origin is ``centre``, and it is always the caller's. An automatic
    search for it was written and then removed: it produced a beautifully
    symmetric picture whenever it was slightly wrong, and nothing about the
    result showed that it had happened. The panels let the origin be clicked
    off the image instead, which is as quick and is visible.

    ``sector`` is an optional ``(phi1, phi2)`` wedge in degrees: only source
    points inside it are used, which is how a map with one good quadrant is
    unfolded into a whole one.

    Each output pixel is the mean of the operations that had data there, and
    ``coverage`` says how many those were. That map is the honest companion
    to the picture: a region built from one copy is a measurement, a region
    built from six is an average, and they should not be read the same way.
    """
    values = _as_2d(values)
    x_axis = np.asarray(x_axis, dtype=float)
    y_axis = np.asarray(y_axis, dtype=float)

    cx, cy = float(centre[0]), float(centre[1])

    gx, gy = np.meshgrid(x_axis - cx, y_axis - cy, indexing="ij")
    radius = np.hypot(gx, gy)
    phi = np.degrees(np.arctan2(gy, gx))

    total = np.zeros_like(values)
    count = np.zeros(values.shape, dtype=float)
    for angle, mirror in symmetry_operations(fold, mirror_angles=mirror_angles,
                                             inversion=inversion):
        if mirror is None:
            theta = phi - angle
        else:
            # Reflect about the mirror line, then undo the rotation.
            theta = 2.0 * mirror - phi - angle
        if sector is not None:
            lo, hi = (float(sector[0]), float(sector[1]))
            wrapped = (theta - lo) % 360.0
            inside = wrapped <= ((hi - lo) % 360.0 or 360.0)
        else:
            inside = np.ones(theta.shape, dtype=bool)
        rad = np.radians(theta)
        sample = _bilinear(values, x_axis, y_axis,
                           cx + radius * np.cos(rad), cy + radius * np.sin(rad))
        good = np.isfinite(sample) & inside
        total = np.where(good, total + np.nan_to_num(sample), total)
        count = count + good

    with np.errstate(invalid="ignore", divide="ignore"):
        out = np.where(count > 0, total / np.maximum(count, 1e-12), np.nan)

    both = np.isfinite(out) & np.isfinite(values)
    if both.any():
        scale = _rms(values[both])
        residual = float(np.sqrt(np.mean((out[both] - values[both]) ** 2)) / scale)
    else:
        residual = float("nan")
    return SymmetryResult(out, x_axis, y_axis, count, (cx, cy), residual)


# ==========================================================================
# 5. Background subtraction
# ==========================================================================
def shirley_background(energy, spectrum, *, iterations: int = 30,
                       tol: float = 1e-6):
    """The Shirley inelastic background of one spectrum.

    The background at each point is proportional to the total spectral
    weight above it, iterated to convergence; the standard choice for a core
    level or a valence-band EDC, and the one thing ``EDC_MDC_plot.m``'s
    straight line cannot represent.
    """
    energy = np.asarray(energy, dtype=float)
    spectrum = np.asarray(spectrum, dtype=float)
    order = np.argsort(energy)
    e = energy[order]
    y = np.nan_to_num(spectrum[order], nan=0.0)
    lo, hi = float(y[0]), float(y[-1])
    # The background has to *end* at each endpoint's own intensity: B(e0) =
    # y[0] and B(e_last) = y[-1], rising towards whichever side carries the
    # inelastic tail. Writing it the other way round produces a curve that
    # crosses the data and is not even monotonic.
    background = np.full(y.shape, hi)
    for _ in range(max(1, int(iterations))):
        above = y - background
        integral = np.concatenate(([0.0], np.cumsum(
            (above[1:] + above[:-1]) / 2.0 * np.diff(e))))
        total = integral[-1]
        if abs(total) < 1e-30:
            break
        new = hi + (lo - hi) * (1.0 - integral / total)
        if np.max(np.abs(new - background)) < tol * max(abs(hi - lo), 1e-12):
            background = new
            break
        background = new
    out = np.empty_like(background)
    out[order] = background
    return out


def tougaard_background(energy, spectrum, *, b: float = 2866.0, c: float = 1643.0):
    """The two-parameter Tougaard universal cross-section background.

    Physically better founded than Shirley for a wide energy window, and
    worth having when the spectrum runs far enough below the peak for the
    loss tail to matter. ``b`` and ``c`` are in eV^2, with the usual metal
    values as defaults.
    """
    energy = np.asarray(energy, dtype=float)
    spectrum = np.asarray(spectrum, dtype=float)
    order = np.argsort(energy)
    e = energy[order]
    y = np.nan_to_num(spectrum[order], nan=0.0)
    step = np.gradient(e)
    background = np.zeros_like(y)
    for i in range(y.size):
        loss = e[i:] - e[i]
        kernel = b * loss / ((c + loss ** 2) ** 2)
        background[i] = np.sum(kernel * y[i:] * step[i:])
    out = np.empty_like(background)
    out[order] = background
    return out


def percentile_background(values, axes, dim: int, *, percentile: float = 5.0,
                          smooth: float = 0.0):
    """The angle-independent background: a low percentile across ``dim``.

    For each position along the other axis, take a low quantile of the
    intensity across ``dim`` and call that the background. For a cut this
    is "the darkest part of each energy", which removes the detector's dark
    level and the secondaries without touching the dispersion -- and unlike
    a fitted straight line it needs no assumption about the shape.
    """
    values = _as_2d(values)
    dim = int(dim) % 2
    level = np.nanpercentile(np.where(np.isfinite(values), values, np.nan),
                             float(percentile), axis=dim, keepdims=True)
    level = np.nan_to_num(level, nan=0.0)
    if smooth > 0:
        other = 1 - dim
        sigma = float(smooth) / axis_step(axes[other])
        level = ndimage.gaussian_filter1d(level, sigma, axis=other, mode="nearest")
    return np.broadcast_to(level, values.shape).copy()


def polynomial_background(values, axes, *, order: int = 2, mask=None):
    """A smooth two-dimensional polynomial fitted to the whole image.

    ``mask`` restricts the fit to the pixels that are background (the
    corners of a map, say); everything else is still subtracted. Useful for
    a slowly varying illumination or a detector gradient that no per-line
    method can see.
    """
    values = _as_2d(values)
    x = np.asarray(axes[0], dtype=float)
    y = np.asarray(axes[1], dtype=float)
    gx, gy = np.meshgrid((x - x.mean()) / (np.ptp(x) or 1.0),
                         (y - y.mean()) / (np.ptp(y) or 1.0), indexing="ij")
    terms = [gx ** i * gy ** j
             for i in range(order + 1) for j in range(order + 1 - i)]
    design = np.column_stack([t.ravel() for t in terms])
    good = np.isfinite(values).ravel()
    if mask is not None:
        good &= np.asarray(mask, dtype=bool).ravel()
    if good.sum() < design.shape[1]:
        raise ValueError("not enough valid points to fit that polynomial order")
    coeffs, *_ = np.linalg.lstsq(design[good], values.ravel()[good], rcond=None)
    return (design @ coeffs).reshape(values.shape)


BACKGROUNDS = ("shirley", "tougaard", "percentile", "polynomial")


def subtract_background(values, axes, method: str, *, dim: int = 1, clip: bool = False,
                        **kwargs):
    """Compute and subtract one of the backgrounds above.

    ``dim`` is the direction a one-dimensional background runs along (1, the
    energy axis, for a cut). ``clip`` holds the result at zero, which is
    what a counting experiment implies but which also hides an over-
    subtraction, so it is off by default.
    """
    values = _as_2d(values)
    if method == "percentile":
        background = percentile_background(values, axes, dim,
                                           percentile=kwargs.get("percentile", 5.0),
                                           smooth=kwargs.get("smooth", 0.0))
    elif method == "polynomial":
        background = polynomial_background(values, axes,
                                           order=kwargs.get("order", 2),
                                           mask=kwargs.get("mask"))
    elif method in ("shirley", "tougaard"):
        axis = np.asarray(axes[dim], dtype=float)
        background = np.empty_like(values)
        function = shirley_background if method == "shirley" else tougaard_background
        extra = {k: v for k, v in kwargs.items() if k in ("iterations", "tol", "b", "c")}
        lines = values.shape[1 - dim]
        for i in range(lines):
            index = (i, slice(None)) if dim == 1 else (slice(None), i)
            background[index] = function(axis, values[index], **extra)
    else:
        raise ValueError(f"background method {method!r} is not one of {BACKGROUNDS}")
    out = values - background
    return (np.maximum(out, 0.0) if clip else out), background


# ==========================================================================
# 6. Fermi-Dirac division  (thin front for tools.fermi.divide_fermi)
# ==========================================================================
def divide_fermi_edge(values, axes, *, ef: float, temperature: float,
                      resolution: float, dim: int = 1, cutoff_kt: float = 4.0):
    """Divide out the resolution-broadened Fermi cut-off along ``dim``.

    This is how the states just above E_F are shown -- they are there, at a
    few percent of the occupied weight, and dividing by the occupation is
    what brings them up without inventing them. Above ``cutoff_kt`` thermal
    widths the divisor is smaller than its own noise, so that region comes
    back as NaN rather than as spikes.
    """
    from tools.fermi import divide_fermi

    values = _as_2d(values)
    return divide_fermi(np.asarray(axes[dim], dtype=float), values, ef,
                        temperature, resolution, energy_dim=int(dim),
                        cutoff_kt=cutoff_kt)


# ==========================================================================
# 7. Normalisation
# ==========================================================================
def normalise(values, axes, dim: int, *, mode: str = "area", window=None):
    """Divide every line along ``dim`` by its own scale.

    ``mode`` is ``"area"`` (the line's integral), ``"max"`` (its peak) or
    ``"mean"``. Normalising every MDC to its own area is the classic way to
    show a dispersion whose intensity falls off with binding energy: it is
    the shape that carries the physics, and the matrix element that carries
    the brightness.

    ``window`` is an optional ``(lo, hi)`` on ``dim``'s own axis limiting
    which points set the scale -- so a normalisation can be tied to a
    feature (a core level, the flat part above E_F) rather than to whatever
    happens to be in the frame.
    """
    values = np.asarray(values, dtype=float)
    dim = int(dim) % values.ndim
    axis = np.asarray(axes[dim], dtype=float)

    region = values
    if window is not None and window != (None, None):
        lo = -np.inf if window[0] is None else float(window[0])
        hi = np.inf if window[1] is None else float(window[1])
        keep = (axis >= min(lo, hi)) & (axis <= max(lo, hi))
        if not keep.any():
            raise ValueError("that window contains no points")
        region = np.take(values, np.flatnonzero(keep), axis=dim)

    with np.errstate(invalid="ignore", divide="ignore"):
        if mode == "max":
            scale = np.nanmax(region, axis=dim, keepdims=True)
        elif mode == "mean":
            scale = np.nanmean(region, axis=dim, keepdims=True)
        elif mode == "area":
            step = axis_step(axis)
            scale = np.nansum(region, axis=dim, keepdims=True) * step
        else:
            raise ValueError(f"normalisation mode {mode!r} is not area, max or mean")
        scale = np.where(np.isfinite(scale) & (np.abs(scale) > 1e-30), scale, np.nan)
        return values / scale


# ==========================================================================
# 8. Despiking
# ==========================================================================
def despike(values, *, threshold: float = 6.0, size: int = 3,
            scale: str = "poisson", replace: bool = True):
    """Find and repair cosmic rays and dead pixels.

    A pixel is a spike when it differs from the median of its neighbourhood
    by more than ``threshold`` noise widths -- a median test, so a genuine
    sharp band is not flattened the way a mean-and-standard-deviation test
    would flatten it.

    How the noise width is estimated is the part that decides whether this
    is usable. Taking it from the same 3x3 neighbourhood, as is tempting,
    estimates it from nine numbers: on pure Gaussian noise that flagged 37
    pixels of a 3600-pixel test image where two were real, because wherever
    those nine happened to agree the threshold collapsed. So:

    ``"poisson"`` (the default) uses ``sqrt(local median)``, which is the
    actual noise of a counting detector and follows the image brightness;
    ``"global"`` uses one robust MAD over the whole frame, right for data
    that has already been normalised or divided; ``"local"`` is the 3x3
    estimate, kept for the rare case of a strongly structured background,
    and then ``size`` is worth raising to 5 or 7.

    Returns ``(values, spike_mask)``; with ``replace`` the spikes carry the
    local median, otherwise they come back as NaN so nothing is invented.
    """
    values = np.asarray(values, dtype=float)
    mask = np.isfinite(values)
    filled = _fill_nans(values)
    median = ndimage.median_filter(filled, size=int(size), mode="nearest")
    deviation = filled - median

    if scale == "poisson":
        width = np.sqrt(np.maximum(median, 1.0))
    elif scale == "global":
        finite = np.abs(deviation[mask])
        width = np.full(values.shape,
                        1.4826 * float(np.median(finite)) if finite.size else 1.0)
    elif scale == "local":
        width = 1.4826 * ndimage.median_filter(np.abs(deviation), size=int(size),
                                               mode="nearest")
    else:
        raise ValueError(f"despike scale {scale!r} is not poisson, global or local")

    spikes = (np.abs(deviation) > float(threshold) * np.maximum(width, 1e-12)) & mask
    out = values.copy()
    out[spikes] = median[spikes] if replace else np.nan
    return out, spikes


# ==========================================================================
# 9. Provenance
# ==========================================================================
@dataclass
class Step:
    """One operation, recorded so the dataset can say how it was made."""
    name: str
    parameters: dict = field(default_factory=dict)
    source: str = ""

    def describe(self) -> str:
        if not self.parameters:
            return self.name
        bits = ", ".join(f"{k}={_short(v)}" for k, v in sorted(self.parameters.items()))
        return f"{self.name}({bits})"


def _short(value) -> str:
    if isinstance(value, float):
        return f"{value:.6g}"
    if isinstance(value, (list, tuple)):
        return "[" + ", ".join(_short(v) for v in value) + "]"
    return str(value)


def history_of(data) -> list:
    """Every processing step recorded on a dataset, oldest first.

    Read back off the ``proc.step.N`` keys in its info, so a dataset saved
    to a file and read back still carries its history -- the thing the
    MATLAB tools never kept, and the reason a figure could not be reproduced
    six months later.
    """
    info = getattr(getattr(data, "scan", None), "info", None) or {}
    steps = []
    for key in sorted(k for k in info if k.startswith("proc.step.")):
        try:
            index = int(key.rsplit(".", 1)[1])
        except ValueError:
            continue
        steps.append((index, str(info[key])))
    return [text for _, text in sorted(steps)]


def record_step(source_info: dict, step: Step) -> dict:
    """Return ``source_info`` plus one more step, ready for MemoryData.

    Numbered rather than appended to one string so the chain survives a
    round trip through the file's flat key/value metadata.
    """
    info = dict(source_info or {})
    existing = [k for k in info if k.startswith("proc.step.")]
    index = len(existing) + 1
    info[f"proc.step.{index}"] = step.describe()
    if step.source:
        info[f"proc.step.{index}.from"] = step.source
    return info
