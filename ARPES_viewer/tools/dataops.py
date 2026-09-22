"""
tools/dataops.py
===============
The three bulk operations from the lab's MATLAB data tools, as plain numpy
functions with no Qt: **truncate** (``td_demo.m``), **self-normalisation**
(``self_normalization_demo.m``) and **compress** (``data_comb_resamp_demo.m``,
its "combine" half -- resampling is deliberately left out, see below).

All three take an array plus its axes and return the same, so they work the
same way on a 2D cut, a 3D map and a 4D spatial scan: the axis *meanings*
are the caller's business.

Why compress rather than resample: resampling keeps every Nth point and
throws the rest away, which on a photon-starved ARPES scan throws away most
of the counts. Combining sums each block instead, so nothing is lost and the
noise actually improves -- which is what "compress" is for here.
"""
from __future__ import annotations

import numpy as np

from loader.nxs_file import AXIS_SLOTS


#: Which scan axes correspond to the array's dimensions, in order, for each
#: data kind. A spatial scan's cube is stored (y, x, k, E) -- the order the
#: file's own 4D dataset is put in -- so that is the order an operation sees
#: it in, whatever the information panel calls the columns.
#:
#: Derived from ``loader.nxs_file.AXIS_SLOTS`` rather than written out again:
#: these two tables and three more scattered copies used to say the same
#: thing separately, so a new kind could be added to one and missed in the
#: others. The names are kept because they are what the rest of the program
#: imports.
ARRAY_AXES = {kind: slots["array"] for kind, slots in AXIS_SLOTS.items()}

#: The order MemoryData's constructor takes axes in, per kind -- not always
#: the array order above (a spatial scan is constructed x, y, k, z).
CONSTRUCTOR_AXES = {kind: slots["constructor"]
                    for kind, slots in AXIS_SLOTS.items()}


# --------------------------------------------------------------------------
# Shared
# --------------------------------------------------------------------------
def axis_index(axis, value, default: int) -> int:
    """Index of ``value`` on ``axis``, the way the MATLAB tools resolve a
    "real scale" bound: round to the nearest step from the axis start.

    ``None`` means "leave it at the end you were heading for", which is how
    the blank/"Min"/"Max" boxes behave there.
    """
    if value is None:
        return int(default)
    axis = np.asarray(axis, dtype=float)
    if axis.size < 2:
        return 0
    step = (axis[-1] - axis[0]) / (axis.size - 1)
    if step == 0:
        return 0
    return int(round((float(value) - axis[0]) / step))


def resolve_bounds(axis, lo, hi, by_index: bool):
    """Turn one axis's (lo, hi) request into a half-open index slice.

    ``by_index`` picks whether the numbers are 1-based indices (as the
    MATLAB boxes are) or values in the axis's own units. Either bound may be
    None for "as far as the data goes". Out-of-range bounds are clamped
    rather than refused -- asking for more than exists means all of it.
    """
    axis = np.asarray(axis, dtype=float)
    n = axis.size
    if by_index:
        first = 0 if lo is None else int(round(lo)) - 1
        last = n - 1 if hi is None else int(round(hi)) - 1
    else:
        first = axis_index(axis, lo, 0)
        last = axis_index(axis, hi, n - 1)
    if first > last:
        first, last = last, first
    first = int(np.clip(first, 0, n - 1))
    last = int(np.clip(last, 0, n - 1))
    return slice(first, last + 1)


# --------------------------------------------------------------------------
# Truncate
# --------------------------------------------------------------------------
def truncate(values, axes, bounds, by_index: bool = False):
    """Cut a rectangular piece out of the data (``td_demo.m``).

    ``axes`` is one array per dimension of ``values``; ``bounds`` is one
    ``(lo, hi)`` pair per axis, either of which may be None. Returns
    ``(values, axes)`` cut to those bounds -- views' worth of data, copied so
    the result owns its memory.
    """
    values = np.asarray(values)
    if len(axes) != values.ndim or len(bounds) != values.ndim:
        raise ValueError(f"{values.ndim}D data needs {values.ndim} axes and bounds")
    slices = tuple(resolve_bounds(axis, lo, hi, by_index)
                   for axis, (lo, hi) in zip(axes, bounds))
    new_axes = [np.asarray(axis, dtype=float)[sl] for axis, sl in zip(axes, slices)]
    if any(axis.size == 0 for axis in new_axes):
        raise ValueError("those bounds leave nothing behind")
    return np.array(values[slices], copy=True), new_axes


# --------------------------------------------------------------------------
# Self-normalisation
# --------------------------------------------------------------------------
def self_normalize(values, axes, dims, *, window=None, by_index: bool = False,
                   to_peak: bool = False):
    """Divide the data by its own intensity, line by line or plane by plane.

    ``dims`` names the directions the normalisation runs *along*: one axis
    (the MATLAB "x"/"y"/"z" buttons) divides every line along it by that
    line's own total; two axes (its "xy"/"yz"/"xz" buttons) divide every
    plane by the plane's total. This is what takes a beam-current drift, or
    a detector whose sensitivity falls off along the slit, out of a scan.

    ``window`` is an optional ``(lo, hi)`` on the *first* named axis
    restricting which points are summed (energies inside a core level, say);
    ``to_peak`` divides by the maximum instead of the sum. NaNs are excluded
    from the sums and put back afterwards, so a missing point neither
    poisons its line nor is invented.
    """
    values = np.asarray(values, dtype=float)
    dims = tuple(int(d) % values.ndim for d in dims)
    if not dims:
        raise ValueError("pick at least one direction to normalise along")
    if len(set(dims)) != len(dims):
        raise ValueError("the same direction twice is not a direction")

    finite = np.isfinite(values)
    filled = np.where(finite, values, 0.0)

    if window is not None and window != (None, None):
        axis = axes[dims[0]]
        sl = resolve_bounds(axis, window[0], window[1], by_index)
        index = [slice(None)] * values.ndim
        index[dims[0]] = sl
        region = filled[tuple(index)]
    else:
        region = filled

    if to_peak:
        norm = np.nanmax(np.where(np.isfinite(region), region, -np.inf),
                         axis=dims, keepdims=True)
        norm = np.where(np.isfinite(norm), norm, np.nan)
    else:
        norm = region.sum(axis=dims, keepdims=True)

    with np.errstate(invalid="ignore", divide="ignore"):
        out = filled / norm
    out = np.where(np.isfinite(out), out, 0.0)
    # A line that was entirely zero (or entirely missing) has no scale of its
    # own; leaving it at zero is honest, dividing by zero is not.
    return np.where(finite, out, np.nan)


# --------------------------------------------------------------------------
# Compress
# --------------------------------------------------------------------------
def compress(values, axes, factors):
    """Bin the data by summing blocks of neighbouring points.

    One factor per axis; 1 leaves that axis alone. Values are **summed**
    over each block and the axis takes the block's **mean**, exactly as
    ``data_comb_resamp_demo.m``'s combine does, so total counts are
    preserved and each new point sits at the centre of what it came from.
    A remainder that does not fill a whole block is dropped (the MATLAB code
    floors the same way): a final point made of half a block would be half
    as bright as its neighbours for no physical reason.
    """
    values = np.asarray(values, dtype=float)
    factors = [max(1, int(f)) for f in factors]
    if len(axes) != values.ndim or len(factors) != values.ndim:
        raise ValueError(f"{values.ndim}D data needs {values.ndim} axes and factors")

    out = values
    new_axes = []
    for dim, (axis, factor) in enumerate(zip(axes, factors)):
        axis = np.asarray(axis, dtype=float)
        n_blocks = out.shape[dim] // factor
        if n_blocks < 1:
            raise ValueError(f"a factor of {factor} is larger than axis {dim} "
                             f"({out.shape[dim]} points)")
        kept = n_blocks * factor
        trimmed = np.take(out, np.arange(kept), axis=dim)
        shape = list(trimmed.shape)
        shape[dim: dim + 1] = [n_blocks, factor]
        out = np.nansum(trimmed.reshape(shape), axis=dim + 1)
        new_axes.append(axis[:kept].reshape(n_blocks, factor).mean(axis=1))
    return out, new_axes


# --------------------------------------------------------------------------
# What a dataset looks like, for the browser's "Data Information" panel and
# for deciding whether several datasets can be operated on together
# --------------------------------------------------------------------------
def axis_summary(axis):
    """``(min, max, num, step)`` for one axis, or None if it has none."""
    if axis is None:
        return None
    axis = np.asarray(axis, dtype=float).ravel()
    if axis.size == 0:
        return None
    step = ((axis[-1] - axis[0]) / (axis.size - 1)) if axis.size > 1 else 0.0
    return (float(axis.min()), float(axis.max()), int(axis.size), float(step))


def same_format(shapes) -> bool:
    """Whether several datasets can be put through the same operation.

    Identical means identical: same kind, same number of axes and the same
    number of points on each. Anything else and one set of bounds or factors
    would mean different things for different datasets, so the caller is
    told to choose again rather than being handed a silently wrong result.
    """
    shapes = list(shapes)
    return bool(shapes) and all(s == shapes[0] for s in shapes)
