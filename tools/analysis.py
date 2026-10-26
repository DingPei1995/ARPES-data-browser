"""
tools/analysis.py
================
The two analysis operations ported from the lab's MATLAB tools, as plain
numpy/scipy functions with no Qt anywhere -- the same split as
``tools/kspace.py``, so both can be tested (and reused) without a display.

- **Arbitrary-direction cut** (``arbitrary_cut``), from
  ``arbi_cut_plot_demo.m``: sample a cube along a path of up to six points
  picked on the constant-energy contour and return intensity against
  distance along that path.
- **Fermi-surface correction** (``fs_correction``), from ``correction.m``:
  fit a polynomial through points picked along a feature that should be
  flat (a Fermi edge, a band bottom) and shift every angle column in energy
  so that it is.

Deviations from the MATLAB originals are deliberate and noted at each
function; they are all cases where the original's mechanics get in the way
of the result rather than defining it.
"""
from __future__ import annotations

import numpy as np

try:
    from scipy.interpolate import RegularGridInterpolator
except ImportError as exc:      # pragma: no cover - scipy is a hard dependency
    raise ImportError("tools.analysis needs scipy (RegularGridInterpolator)") from exc


# --------------------------------------------------------------------------
# Arbitrary-direction cut
# --------------------------------------------------------------------------
def _axis_step(axis) -> float:
    axis = np.asarray(axis, dtype=float)
    if axis.size < 2:
        return 1.0
    return float(abs(np.mean(np.diff(axis))))


def line_samples(x_axis, y_axis, p0, p1):
    """Points along the segment p0 -> p1, sampled at the grid's own spacing.

    Returns ``(xs, ys, r)`` with ``r`` the distance from ``p0``.

    The `.m` file samples with ``linspace`` over *x* and swaps the axes when
    the line is within 1e-5 of vertical, because a line parameterised by x
    degenerates as it steepens. Sampling along the segment itself has no
    such special case and gives evenly spaced points for every direction,
    so that is what this does -- the two agree wherever the MATLAB version
    is well behaved.
    """
    x0, y0 = (float(v) for v in p0)
    x1, y1 = (float(v) for v in p1)
    length = float(np.hypot(x1 - x0, y1 - y0))
    step = min(_axis_step(x_axis), _axis_step(y_axis))
    if length <= 0 or step <= 0:
        return np.array([x0]), np.array([y0]), np.array([0.0])
    n = int(np.ceil(length / step)) + 1
    t = np.linspace(0.0, 1.0, n)
    return x0 + t * (x1 - x0), y0 + t * (y1 - y0), t * length


def arbitrary_cut(x_axis, y_axis, z_axis, cube, points):
    """Intensity along a path of 2..6 points through an (x, y, z) cube.

    ``cube`` is indexed ``[x, y, z]`` -- for a deflector map, (deflector
    angle, slit angle, energy). ``points`` is a list of (x, y) corners; the
    path runs through them in order.

    Returns ``(distance, z_axis, values, joints)``:
    ``values`` has shape ``(len(distance), len(z_axis))``, ``distance`` is
    cumulative along the whole path (so the segments are laid end to end),
    and ``joints`` gives the distance at each corner between two segments,
    for drawing the boundary markers the MATLAB GUI offers.

    Segments are joined by dropping each one's first sample, which is the
    previous segment's last sample. The `.m` file instead writes a NaN at
    every joint and calls ``inpaint_nans`` to fill it back in; not
    duplicating the point in the first place needs no repair.
    """
    x_axis = np.asarray(x_axis, dtype=float)
    y_axis = np.asarray(y_axis, dtype=float)
    z_axis = np.asarray(z_axis, dtype=float)
    cube = np.asarray(cube, dtype=float)
    pts = [(float(px), float(py)) for px, py in points]
    if len(pts) < 2:
        raise ValueError("an arbitrary cut needs at least two points")
    if cube.shape[:2] != (x_axis.size, y_axis.size):
        raise ValueError(f"cube {cube.shape} does not match axes "
                         f"({x_axis.size}, {y_axis.size})")

    # NaNs would poison every interpolated sample around them, so they are
    # carried as weights instead: interpolate the values and the mask, then
    # divide, exactly as the display-side smoothing does.
    finite = np.isfinite(cube)
    interp = RegularGridInterpolator((x_axis, y_axis), np.where(finite, cube, 0.0),
                                     bounds_error=False, fill_value=np.nan)
    weight = RegularGridInterpolator((x_axis, y_axis), finite.astype(float),
                                     bounds_error=False, fill_value=np.nan)

    chunks, dists, joints = [], [], []
    offset = 0.0
    for i in range(len(pts) - 1):
        xs, ys, r = line_samples(x_axis, y_axis, pts[i], pts[i + 1])
        if i:                       # drop the duplicated corner sample
            xs, ys, r = xs[1:], ys[1:], r[1:]
            joints.append(offset)
        sample = np.column_stack([np.clip(xs, x_axis.min(), x_axis.max()),
                                  np.clip(ys, y_axis.min(), y_axis.max())])
        num = interp(sample)
        den = weight(sample)
        with np.errstate(invalid="ignore", divide="ignore"):
            values = num / den
        chunks.append(np.where(den > 1e-9, values, np.nan))
        dists.append(offset + r)
        offset += float(r[-1]) if r.size else 0.0

    return (np.concatenate(dists), z_axis, np.concatenate(chunks, axis=0),
            np.asarray(joints, dtype=float))


# --------------------------------------------------------------------------
# Fermi-surface correction
# --------------------------------------------------------------------------
def fit_feature(points, order: int = 2):
    """Least-squares polynomial through the picked (angle, energy) points,
    as ``correction.m``'s ``polyfit(m, n, 2)``. Returns the coefficients,
    highest power first, ready for ``np.polyval``."""
    pts = np.asarray(points, dtype=float)
    if pts.ndim != 2 or pts.shape[1] != 2:
        raise ValueError("points must be a list of (angle, energy) pairs")
    order = int(order)
    if pts.shape[0] < order + 1:
        raise ValueError(f"an order-{order} fit needs at least {order + 1} points")
    return np.polyfit(pts[:, 0], pts[:, 1], order)


def shift_bins(angle_axis, energy_axis, coeffs):
    """How many energy bins each angle column has to move so the fitted
    feature comes out flat.

    The column where the fit is highest keeps its energies (shift 0) and
    every other column moves down by the fit's deficit there, which is the
    MATLAB code's convention -- both of its branches (positive and negative
    quadratic term) reduce to exactly this, they only differ in which end
    the padding lands on, and it trims that away afterwards anyway.
    """
    angle_axis = np.asarray(angle_axis, dtype=float)
    step = _axis_step(energy_axis)
    if step <= 0:
        raise ValueError("the energy axis has no usable step")
    fitted = np.polyval(coeffs, angle_axis)
    return np.round((fitted.max() - fitted) / step).astype(int), fitted


def fs_correction(values, angle_axis, energy_axis, coeffs, *,
                  angle_dim: int = 0, energy_dim: int = -1):
    """Straighten a curved feature by shifting each angle column in energy.

    Works on a 2D frame (angle, energy) and on a 3D cube alike: ``angle_dim``
    names the axis the fit was made along and every other axis is carried
    through untouched, which is how a fit made on one slit cut corrects a
    whole (deflector, slit, energy) map.

    Returns ``(corrected, new_energy_axis)``. The energy axis grows by the
    span of the fit, so nothing is cropped; the bins that no column reaches
    are NaN, and are left NaN rather than filled -- there is no measurement
    there, and the display already treats NaN as missing.

    Deviation: the bin width is taken from the axis itself
    (``mean(diff(E))``) rather than the `.m` file's ``range/N``, which is off
    by one bin in N and skews the shift slightly across a long axis.
    """
    values = np.asarray(values, dtype=float)
    energy_axis = np.asarray(energy_axis, dtype=float)
    shifts, fitted = shift_bins(angle_axis, energy_axis, coeffs)

    moved = np.moveaxis(values, (angle_dim, energy_dim), (0, -1))
    n_angle, n_energy = moved.shape[0], moved.shape[-1]
    if shifts.size != n_angle:
        raise ValueError(f"{shifts.size} shifts for {n_angle} angle points")
    new_n = int(n_energy + shifts.max())
    out = np.full(moved.shape[:-1] + (new_n,), np.nan)
    for i, s in enumerate(shifts):
        out[i, ..., s:s + n_energy] = moved[i, ...]

    # Trim energy bins no column reached (belt and braces: with this shift
    # convention the first and last are always covered).
    used = np.any(np.isfinite(out), axis=tuple(range(out.ndim - 1)))
    if used.any():
        lo, hi = int(np.argmax(used)), int(len(used) - np.argmax(used[::-1]))
        out = out[..., lo:hi]
    else:
        lo, hi = 0, new_n

    step = _axis_step(energy_axis)
    new_axis = float(energy_axis[0]) + (np.arange(lo, hi) * step)
    return np.moveaxis(out, (0, -1), (angle_dim, energy_dim)), new_axis
