"""
tools/processing/combine.py
===========================
Datasets measured in pieces, joined back into one.

A map whose deflector range was taken in two or three scans, a kz map
whose photon-energy range was split over several runs, a cut whose energy
window was measured in two halves: each piece is the same kind of data on
the same axes, except that one axis covers a different stretch in each.
:func:`join` puts them end to end along that axis.

What it does with the axes
--------------------------
*The axis joined along* becomes one regular axis from the lowest value of
any piece to the highest, with the first piece's step. When every piece
sits on that one lattice (the usual case: same step, ranges that abut or
overlap by whole points) each value is copied, not interpolated. Otherwise
-- a different step, a half-step offset -- each piece is linearly
interpolated onto the common axis, and the notes say so. Where pieces
overlap, the result is their average, or the first / last piece in the
chosen order. A stretch no piece covers is a gap: it is filled with zeros
and reported, since inventing data there would be worse.

*The other axes* must describe the same thing. When they agree point for
point (within :data:`TOLERANCE` of a step) they are kept as they are; when
they do not -- an energy window a few meV different between two runs, a
different number of slit channels -- every piece is resampled onto the
first piece's points inside the range all pieces cover, and the notes say
so.

Nothing here imports Qt; :func:`plan` describes what :func:`join` would do
without touching the data, which is what the dialog shows.
"""
from __future__ import annotations

import numpy as np

__all__ = ["TOLERANCE", "OVERLAP_MODES", "plan", "join", "interp_axis"]

#: How far two axis values may differ, as a fraction of a step, and still
#: be called the same point.
TOLERANCE = 0.02

#: What to do where two pieces cover the same points.
OVERLAP_MODES = ("average", "first", "last")


def _step(axis) -> float:
    axis = np.asarray(axis, dtype=float)
    if axis.size < 2:
        return 0.0
    return float(np.median(np.abs(np.diff(axis))))


def _same_axis(a, b) -> bool:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.size != b.size:
        return False
    if a.size < 2:
        return bool(np.allclose(a, b))
    step = _step(a) or 1.0
    off = float(np.nanmax(np.abs(a - b)))
    return np.isfinite(off) and off <= TOLERANCE * step


def interp_axis(values, axis_values, new_axis, axis: int) -> np.ndarray:
    """``values`` linearly interpolated along ``axis`` from the points
    ``axis_values`` onto ``new_axis`` (which must lie inside their range).
    Works for ascending and descending ``axis_values``."""
    values = np.asarray(values)
    old = np.asarray(axis_values, dtype=float)
    new = np.asarray(new_axis, dtype=float)
    if old.size == 1:
        return np.repeat(np.take(values, [0], axis=axis), new.size, axis=axis)
    if old[0] > old[-1]:
        old = old[::-1]
        values = np.flip(values, axis=axis)
    position = np.interp(new, old, np.arange(old.size, dtype=float))
    low = np.clip(np.floor(position).astype(int), 0, old.size - 1)
    high = np.clip(low + 1, 0, old.size - 1)
    weight = (position - low).astype(np.float32)
    shape = [1] * values.ndim
    shape[axis] = new.size
    weight = weight.reshape(shape)
    a = np.take(values, low, axis=axis).astype(np.float32, copy=False)
    b = np.take(values, high, axis=axis).astype(np.float32, copy=False)
    return a * (1 - weight) + b * weight


def _common_axis(axes, index: int):
    """The points of the first piece's axis ``index`` that every piece
    covers, and whether resampling is needed at all."""
    first = np.asarray(axes[0][index], dtype=float)
    if all(_same_axis(first, piece[index]) for piece in axes[1:]):
        return first, False
    lo = max(float(np.min(piece[index])) for piece in axes)
    hi = min(float(np.max(piece[index])) for piece in axes)
    slack = TOLERANCE * (_step(first) or 0.0)
    keep = (first >= lo - slack) & (first <= hi + slack)
    common = first[keep]
    if common.size < 2:
        raise ValueError(
            f"axis {index + 1}: the pieces have no common range "
            f"({lo:.4g} to {hi:.4g}) to put them on")
    return np.clip(common, lo, hi), True


def guess_join_axis(axes) -> int:
    """The axis pieces with these ``axes`` were most likely split along.

    The one whose range differs most between pieces, measured in steps: a
    deflector range scanned in two parts differs only along the deflector,
    a cut taken in two energy windows only along the energy. When every
    axis agrees -- nothing to tell them apart -- the first, which is the
    scanned axis of a map or kz map.
    """
    if not axes:
        return 0
    best, best_shift = 0, 0.0
    for index in range(len(axes[0])):
        first = np.asarray(axes[0][index], dtype=float)
        step = _step(first) or 1.0
        shift = 0.0
        for piece in axes[1:]:
            other = np.asarray(piece[index], dtype=float)
            if _same_axis(first, other) or not (first.size and other.size):
                continue
            shift = max(shift,
                        abs(float(np.min(other)) - float(np.min(first))) / step,
                        abs(float(np.max(other)) - float(np.max(first))) / step)
        if shift > best_shift:
            best, best_shift = index, shift
    return best


def plan(axes, axis: int, names=None) -> dict:
    """What joining pieces with these ``axes`` along ``axis`` gives.

    ``axes`` is ``[(axis0, axis1, ...), ...]``, one tuple per piece, in
    the chosen order. Returns ``{"axes": [...], "exact": bool,
    "resampled": [indices], "gap_points": int, "overlap_points": int,
    "notes": [str]}``; raises ``ValueError`` (with a message to show)
    when they cannot be joined.
    """
    axes = [tuple(np.asarray(a, dtype=float).ravel() for a in piece)
            for piece in axes]
    names = list(names or [f"piece {i + 1}" for i in range(len(axes))])
    if len(axes) < 2:
        raise ValueError("Select at least two datasets to combine.")
    ndim = len(axes[0])
    if any(len(piece) != ndim for piece in axes):
        raise ValueError("The datasets do not have the same number of axes.")
    if not 0 <= axis < ndim:
        raise ValueError(f"There is no axis {axis + 1} to join along.")

    out_axes, resampled, notes = [], [], []
    for index in range(ndim):
        if index == axis:
            out_axes.append(None)
            continue
        common, changed = _common_axis(axes, index)
        if changed:
            resampled.append(index)
            notes.append(
                f"Axis {index + 1} differs between the pieces: all are "
                f"resampled onto {common.size} points, "
                f"{common[0]:.4g} to {common[-1]:.4g} (the range they all "
                f"cover).")
        out_axes.append(common)

    # -- the joined axis ------------------------------------------------------
    pieces = [piece[axis] for piece in axes]
    step = next((s for s in (_step(p) for p in pieces) if s > 0), 0.0)
    if step <= 0:
        raise ValueError("Every piece has a single point along this axis; "
                         "there is no step to join them with. (To put "
                         "single cuts together, stack them instead.)")
    origin = float(pieces[0][0])
    every = np.concatenate(pieces)
    offsets = (every - origin) / step
    exact = bool(np.all(np.abs(offsets - np.round(offsets)) <= TOLERANCE * 4))
    for name, piece in zip(names[1:], pieces[1:]):
        other = _step(piece)
        if other > 0 and abs(other - step) > TOLERANCE * step:
            exact = False
            notes.append(f"“{name}” has a step of {other:.4g}, not "
                         f"{step:.4g}: it is interpolated onto the common "
                         f"axis.")
    if exact:
        first, last = int(np.round(offsets.min())), int(np.round(offsets.max()))
        joined = origin + step * np.arange(first, last + 1, dtype=float)
    else:
        lo, hi = float(every.min()), float(every.max())
        count = int(np.floor((hi - lo) / step + TOLERANCE)) + 1
        joined = lo + step * np.arange(count, dtype=float)
        if not any("interpolated" in n for n in notes):
            notes.append("The pieces do not sit on one grid along this "
                         "axis: each is interpolated onto the common axis.")
    out_axes[axis] = joined

    covered = np.zeros(joined.size, dtype=int)
    for piece in pieces:
        covered += _coverage(piece, joined, step)
    gaps = int(np.sum(covered == 0))
    overlaps = int(np.sum(covered > 1))
    if gaps:
        notes.append(f"{gaps} point(s) along the joined axis are in no "
                     f"piece: they are left at zero.")
    if overlaps:
        notes.append(f"{overlaps} point(s) along the joined axis are in "
                     f"more than one piece.")
    return {"axes": out_axes, "exact": exact, "resampled": resampled,
            "gap_points": gaps, "overlap_points": overlaps, "notes": notes,
            "step": step}


def _coverage(piece, joined, step) -> np.ndarray:
    slack = TOLERANCE * 4 * step
    lo, hi = float(np.min(piece)), float(np.max(piece))
    return ((joined >= lo - slack) & (joined <= hi + slack)).astype(int)


def join(pieces, axis: int, overlap: str = "average", names=None):
    """Join ``pieces`` end to end along ``axis``.

    ``pieces`` is ``[(axes, values), ...]`` in the chosen order, ``axes``
    a tuple of 1-D arrays and ``values`` the array they index. ``overlap``
    is one of :data:`OVERLAP_MODES`. Returns ``(axes, values, plan)`` with
    ``values`` as float32.
    """
    if overlap not in OVERLAP_MODES:
        raise ValueError(f"overlap must be one of {OVERLAP_MODES}")
    info = plan([a for a, _ in pieces], axis, names)
    out_axes = info["axes"]
    joined, step = out_axes[axis], info["step"]
    shape = tuple(a.size for a in out_axes)
    total = np.zeros(shape, dtype=np.float32)
    count = np.zeros(joined.size, dtype=np.float32)
    bshape = [1] * len(shape)
    bshape[axis] = joined.size

    order = list(range(len(pieces)))
    if overlap == "last":
        order = order[::-1]
    for index in order:
        axes, values = pieces[index]
        axes = [np.asarray(a, dtype=float).ravel() for a in axes]
        values = np.asarray(values)
        for other in info["resampled"]:
            values = interp_axis(values, axes[other], out_axes[other], other)
        mask = _coverage(axes[axis], joined, step).astype(bool)
        if overlap in ("first", "last"):
            mask &= count == 0
        if not mask.any():
            continue
        targets = joined[mask]
        if info["exact"]:
            # Whole points: pick each target's own sample.
            positions = np.rint((targets - axes[axis][0]) /
                                (axes[axis][1] - axes[axis][0])
                                if axes[axis].size > 1 else
                                np.zeros(targets.size)).astype(int)
            positions = np.clip(positions, 0, axes[axis].size - 1)
            part = np.take(values, positions, axis=axis).astype(np.float32)
        else:
            inside = np.clip(targets, np.min(axes[axis]), np.max(axes[axis]))
            part = interp_axis(values, axes[axis], inside, axis)
        where = np.nonzero(mask)[0]
        slicer = [slice(None)] * len(shape)
        slicer[axis] = where
        total[tuple(slicer)] += part
        count[where] += 1
    scale = np.where(count > 0, 1.0 / np.maximum(count, 1), 0.0)
    total *= scale.astype(np.float32).reshape(bshape)
    return out_axes, total, info
