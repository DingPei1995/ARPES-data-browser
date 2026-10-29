"""
tools/momentum/rebin.py
=======================
Averaging for a conversion onto a grid coarser than the measurement.

A conversion reads its output points off the measured grid by
interpolation. That is right when the output is as fine as the data or
finer, and wrong when it is coarser: an output point then stands for
several measured samples, and interpolation *picks* one or two of them
instead of averaging them -- throwing counts away (a noisier result than
the same grid binned) and letting detail finer than the new step alias
into it. So a reducing conversion averages first:

* along an axis that is resampled directly (energy), each output point is
  the mean over its own bin, with partly covered samples weighted by the
  overlap (:func:`bin_average`);
* where the output is read off the measured plane at computed positions
  (the angle axes), each plane is box-averaged over the footprint of one
  output pixel before it is sampled (:func:`box_filter`).

Steps within :data:`REDUCTION_THRESHOLD` of the measured one are left
alone: the measured step in k itself varies by a few per cent across an
energy window (k scales as sqrt(E)) and across the angles, so "as many
points as measured" is that range, and converting at the measured
sampling must stay exactly what it was.

Both are NaN-aware: a missing sample is left out of its average rather
than spreading.
"""
from __future__ import annotations

import numpy as np

#: An output step this many times the measured one, or less, is the
#: measured sampling and is interpolated as before; beyond it, averaged.
REDUCTION_THRESHOLD = 1.1


def is_reduction(width: float) -> bool:
    """Whether an output step of ``width`` measured steps is a reduction."""
    return float(width) > REDUCTION_THRESHOLD


def box_kernel(width: float) -> np.ndarray:
    """Weights of a box ``width`` samples wide, centred on a sample.

    Each sample covers one unit cell, so a fractional width takes the
    overlapping part of the outermost samples: 2.5 gives
    ``[0.25, 1, 1, 0.25]``-ish weights over five taps, normalised.
    """
    half = 0.5 * float(width)
    reach = int(np.ceil(half - 0.5))
    offsets = np.arange(-reach, reach + 1, dtype=float)
    weights = np.clip(np.minimum(offsets + 0.5, half) - np.maximum(offsets - 0.5, -half),
                      0.0, None)
    return weights / weights.sum()


def box_filter(plane: np.ndarray, widths) -> np.ndarray:
    """Box-average a 2-D ``plane`` over ``widths`` samples per axis.

    An axis whose width is not a reduction is left untouched. Missing
    samples are excluded from each average (and a point with no measured
    neighbour at all stays NaN)."""
    from scipy.ndimage import convolve1d

    widths = [float(w) for w in widths]
    if not any(is_reduction(w) for w in widths):
        return plane
    finite = np.isfinite(plane)
    values = np.where(finite, plane, 0.0).astype(float)
    weight = finite.astype(float)
    for axis, width in enumerate(widths):
        if not is_reduction(width):
            continue
        kernel = box_kernel(width)
        values = convolve1d(values, kernel, axis=axis, mode="constant", cval=0.0)
        weight = convolve1d(weight, kernel, axis=axis, mode="constant", cval=0.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = values / weight
    return np.where(weight > 1e-12, out, np.nan)


def bin_average(cube: np.ndarray, axis_values: np.ndarray, target: np.ndarray,
                axis: int = -1, dtype=float) -> np.ndarray:
    """Average ``cube`` along ``axis`` into bins centred on ``target``.

    The bins tile the target range (each one target step wide, the end ones
    trimmed to the measured range), and each measured sample counts by the
    part of its own cell that falls inside the bin -- so a constant stays
    that constant, the mean is conserved, and nothing between two bins is
    lost or counted twice. ``axis_values`` may be descending or uneven.
    One output plane at a time: no full-size temporaries.
    """
    cube = np.moveaxis(np.asarray(cube), axis, -1)
    axis_values = np.asarray(axis_values, dtype=float)
    target = np.asarray(target, dtype=float)
    order = np.argsort(axis_values)
    src = axis_values[order]
    n = src.size
    index = np.arange(n, dtype=float)
    # bin edges halfway between targets, the ends half a step outside
    t = np.sort(target)
    step = np.diff(t)
    edges = np.concatenate([[t[0] - 0.5 * step[0]], t[:-1] + 0.5 * step,
                            [t[-1] + 0.5 * step[-1]]])
    # into fractional source indices (cells are [i - 0.5, i + 0.5])
    lo_edge = np.interp(edges, src, index, left=-0.5, right=n - 0.5)
    lo_edge = np.clip(lo_edge, -0.5, n - 0.5)
    rank = np.argsort(np.argsort(target))       # target -> its sorted slot
    out = np.full(cube.shape[:-1] + (target.size,), np.nan, dtype=dtype)
    for j in range(target.size):
        k = rank[j]
        a, b = lo_edge[k], lo_edge[k + 1]
        first = int(np.floor(a + 0.5))
        last = min(int(np.ceil(b - 0.5)), n - 1)
        total = np.zeros(cube.shape[:-1])
        weight = np.zeros(cube.shape[:-1])
        for i in range(max(first, 0), last + 1):
            w = min(i + 0.5, b) - max(i - 0.5, a)
            if w <= 0:
                continue
            plane = cube[..., order[i]]
            ok = np.isfinite(plane)
            total += np.where(ok, plane, 0.0) * w
            weight += ok * w
        with np.errstate(invalid="ignore", divide="ignore"):
            out[..., j] = np.where(weight > 0, total / np.where(weight > 0, weight, 1.0), np.nan)
    return np.moveaxis(out, -1, axis)
