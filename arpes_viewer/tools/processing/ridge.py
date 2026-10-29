"""
tools/processing/ridge.py
=========================
Bands traced as ridges of the intensity.

A band in a cut (or a Fermi-surface contour in a constant-energy map) is a
ridge of the photoemission intensity: across it the intensity has a
maximum, along it it hardly changes. Ridge detection finds exactly those
points from the Hessian of the smoothed image -- the direction of strongest
downward curvature is "across the band", and a ridge point is where the
slope along that direction passes through zero. Unlike the second
derivative or the curvature it does not draw a band as a pair of lobes and
it leaves no false lines between two close bands; and every ridge point is
located to a fraction of a pixel, so the result can be fitted.

Ported from the lab's ``ridge_detection.m`` (and ``ridge_detection_3.m``),
with what made that one depend on things it should not:

* **Smoothing in each axis' own units.** The .m file used one radius,
  ``radius / sqrt(dx * dy)`` pixels for both axes -- a number that changes
  with the units the axes happen to be in (degrees or Å⁻¹, eV or meV) and
  smooths an axis with a coarse step as much as one with a fine one. Here
  the smoothing width is given per axis in that axis' unit.
* **A Hessian that does not depend on the units either.** The .m file took
  finite differences in physical units, so which direction counted as
  "across the band" turned on the ratio of eV to degrees. Here the
  derivatives are taken in units of the smoothing width along each axis
  (scale-normalised coordinates), which is the natural metric of the
  smoothed image and the same whatever the axes are measured in.
* **Exact derivatives of the smoothed image.** Gaussian-derivative filters
  instead of differences of a smoothed image: the same smoothing, one pass
  per derivative, no extra noise from the differencing.
* **A ridge test that does not depend on the intensity scale.** The .m file
  accepted a pixel when the slope was below an absolute tolerance, a value
  that had to be retuned for every image. Here a ridge is where the slope
  across the band changes sign: the sub-pixel offset ``-L_p / L_pp`` of the
  maximum along that direction must fall inside the pixel. The sharpness
  threshold is a fraction of the strongest ridge.
* **Vectorised.** The .m file looped over pixels with ``arrayfun``; this is
  a few whole-array operations (a 1000 x 1000 cut takes a fraction of a
  second, a map is done plane by plane).

Qt-free, so it can be tested on its own.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage

__all__ = ["RidgeResult", "ridges", "ridge_cube", "default_sigmas"]


class RidgeResult:
    """``strength``: the ridge image (0 off the ridges), on the input grid;
    ``points``: ``(N, 3)`` of sub-pixel ridge points ``x, y, strength``."""

    def __init__(self, strength, points):
        self.strength = strength
        self.points = points


def _step(axis) -> float:
    axis = np.asarray(axis, dtype=float)
    return float(np.median(np.diff(axis))) if axis.size > 1 else 1.0


def default_sigmas(x_axis, y_axis, pixels: float = 2.0):
    """A starting smoothing width per axis: ``pixels`` steps of each."""
    return abs(_step(x_axis)) * pixels, abs(_step(y_axis)) * pixels


def ridges(values, x_axis, y_axis, sigma_x: float, sigma_y: float, *,
           threshold: float = 0.05, gamma: float = 1.0,
           kind: str = "max") -> RidgeResult:
    """Ridge detection on one 2-D image ``values[x, y]``.

    ``sigma_x``/``sigma_y`` are the Gaussian smoothing widths in the units
    of ``x_axis``/``y_axis``. ``threshold`` drops ridges flatter than that
    fraction of the sharpest one. ``gamma`` shapes the returned image
    (``strength ** gamma``, normalised to 1) for display, not the points.
    ``kind`` "max" traces maxima (bands), "min" valleys (gaps).
    """
    values = np.asarray(values, dtype=float)
    if values.ndim != 2 or min(values.shape) < 3:
        raise ValueError("ridge detection needs a 2-D image of at least 3 x 3")
    x_axis = np.asarray(x_axis, dtype=float)
    y_axis = np.asarray(y_axis, dtype=float)
    dx, dy = _step(x_axis), _step(y_axis)
    sx = max(0.5, abs(sigma_x / dx)) if dx else 1.0
    sy = max(0.5, abs(sigma_y / dy)) if dy else 1.0

    finite = np.isfinite(values)
    image = np.where(finite, values, np.nanmean(values) if finite.any() else 0)
    lo, hi = float(np.min(image)), float(np.max(image))
    image = (image - lo) / (hi - lo) if hi > lo else np.zeros_like(image)
    if kind == "min":
        image = 1.0 - image

    sigma = (sx, sy)
    # Derivatives in units of the smoothing width along each axis:
    # d/du = sx d/di, so second derivatives pick up sx^2, sx*sy, sy^2.
    lu = ndimage.gaussian_filter(image, sigma, order=(1, 0)) * sx
    lv = ndimage.gaussian_filter(image, sigma, order=(0, 1)) * sy
    luu = ndimage.gaussian_filter(image, sigma, order=(2, 0)) * sx * sx
    lvv = ndimage.gaussian_filter(image, sigma, order=(0, 2)) * sy * sy
    luv = ndimage.gaussian_filter(image, sigma, order=(1, 1)) * sx * sy

    # Eigenvalues of [[luu, luv], [luv, lvv]], lam1 <= lam2, and the unit
    # eigenvector (pu, pv) of lam1: the direction across the ridge.
    half_trace = 0.5 * (luu + lvv)
    root = np.sqrt((0.5 * (luu - lvv)) ** 2 + luv ** 2)
    lam1, lam2 = half_trace - root, half_trace + root
    pu = luv
    pv = lam1 - luu
    # Where luv ~ 0 the formula degenerates; pick the axis directly.
    flat = (np.abs(pu) + np.abs(pv)) < 1e-12
    pu = np.where(flat, np.where(luu <= lvv, 1.0, 0.0), pu)
    pv = np.where(flat, np.where(luu <= lvv, 0.0, 1.0), pv)
    norm = np.hypot(pu, pv)
    pu, pv = pu / norm, pv / norm

    lp = pu * lu + pv * lv
    with np.errstate(divide="ignore", invalid="ignore"):
        t = -lp / lam1
    # The maximum along p, in pixels of each axis.
    di, dj = t * pu * sx, t * pv * sy
    candidate = ((lam1 < 0) & (np.abs(lam1) > np.abs(lam2))
                 & (np.abs(di) <= 0.5) & (np.abs(dj) <= 0.5) & finite)
    strength = np.where(candidate, -lam1, 0.0)
    peak = float(strength.max()) if strength.size else 0.0
    if peak > 0:
        strength = np.where(strength >= threshold * peak, strength, 0.0)
    keep = strength > 0
    ii, jj = np.nonzero(keep)
    xs = np.interp(ii + di[keep], np.arange(x_axis.size), x_axis)
    ys = np.interp(jj + dj[keep], np.arange(y_axis.size), y_axis)
    points = np.column_stack([xs, ys, strength[keep]]) if ii.size else np.zeros((0, 3))
    shown = strength / peak if peak > 0 else strength
    if gamma != 1.0:
        shown = shown ** float(gamma)
    return RidgeResult(shown, points)


def ridge_cube(values, axes, hold: int, sigmas, **options):
    """Ridge detection plane by plane through a cube.

    ``hold`` is the axis held constant (2 for constant-energy planes of an
    (x, y, E) cube); ``sigmas`` the smoothing widths of the other two axes,
    in order. Returns the ridge cube on the input grid, and the points as
    ``(N, 4)``: the held axis' value first, then the plane's x, y and
    strength.
    """
    values = np.asarray(values, dtype=float)
    plane_axes = [d for d in range(3) if d != hold]
    out = np.zeros_like(values)
    points = []
    for index, held in enumerate(np.asarray(axes[hold], dtype=float)):
        slicer = [slice(None)] * 3
        slicer[hold] = index
        plane = values[tuple(slicer)]
        result = ridges(plane, axes[plane_axes[0]], axes[plane_axes[1]],
                        sigmas[0], sigmas[1], **options)
        out[tuple(slicer)] = result.strength
        if result.points.size:
            points.append(np.column_stack([np.full(len(result.points), held),
                                           result.points]))
    return out, (np.vstack(points) if points else np.zeros((0, 4)))
