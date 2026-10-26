"""
tools/processing/stackcuts.py
=============================
Several cuts stacked into one cube along an axis the user names.

A temperature series, a gate-voltage series, cuts taken one photon energy at
a time: each is a set of cuts measured with the same analyser settings, and
what changes between them is a quantity the files often do not record in any
usable way. Stacking them gives a (scanned, angle, energy) cube that opens
in the map viewer, so the series can be sliced at constant energy or along
the new axis like any other map.

What this module insists on
---------------------------
The cuts must share their grid: the same number of points on both axes,
the same axis values (within a small fraction of a step), and the same
units. Stacking is then exact -- no resampling, nothing interpolated -- and
the result's angle and energy axes are simply the first cut's. Cuts that
differ are refused with the first difference found, rather than resampled
onto each other: in a temperature series a shifted energy window is
usually a sign that the cuts are not what the user thinks they are.

The new axis is ``start + step * i`` in the chosen order, with a name and a
unit of the user's own. Its *role* (see ``loader.registry.AXIS_ROLES``)
decides what the cube is: a photon-energy series is a ``kz_map``, so the
k_z tools take it; anything else is a ``map`` whose first axis is marked as
not an angle, so nothing offers to convert it to momentum.

Nothing here imports Qt.
"""
from __future__ import annotations

import numpy as np

from tools.processing.cutops import axis_unit

__all__ = ["compatibility", "stack_axis", "stack", "axis_label",
           "kind_for_role"]

#: How far two axes may differ, as a fraction of their step, and still be
#: called the same grid. Readbacks jitter; a real offset does not.
TOLERANCE = 0.02


def _axis_problem(name: str, first, other, label: str) -> str:
    first = np.asarray(first, dtype=float)
    other = np.asarray(other, dtype=float)
    if first.size != other.size:
        return (f"its {label} axis has {other.size} points, "
                f"not {first.size}")
    if first.size < 2:
        return ""
    step = np.median(np.abs(np.diff(first))) or 1.0
    off = float(np.nanmax(np.abs(first - other)))
    if not np.isfinite(off) or off > TOLERANCE * step:
        return (f"its {label} axis runs {other[0]:.4g} to {other[-1]:.4g}, "
                f"not {first[0]:.4g} to {first[-1]:.4g}")
    return ""


def compatibility(cuts) -> str:
    """Why ``cuts`` cannot be stacked, or ``""`` if they can.

    ``cuts`` is ``[(name, (x, y), values, labels), ...]`` with ``labels`` a
    dict holding ``"x"`` and ``"y"``. The first cut is the reference; the
    message names the first one that differs from it and how.
    """
    if len(cuts) < 2:
        return "Select at least two cuts to combine."
    name0, (x0, y0), v0, labels0 = cuts[0]
    shape0 = np.shape(v0)
    if len(shape0) != 2:
        return f"“{name0}” is not a 2-D cut."
    for name, (x, y), values, labels in cuts[1:]:
        if np.shape(values) != shape0:
            return (f"“{name}” is {np.shape(values)[0]} × "
                    f"{np.shape(values)[1]} points, but “{name0}” is "
                    f"{shape0[0]} × {shape0[1]}.")
        for slot, what in (("x", "angle / momentum"), ("y", "energy")):
            a, b = axis_unit(labels0.get(slot, "")), axis_unit(
                labels.get(slot, ""))
            if a and b and a != b:
                return (f"“{name}”'s {what} axis is in {b}, but "
                        f"“{name0}”'s is in {a}.")
        for axis0, axis, what in ((x0, x, "angle / momentum"),
                                  (y0, y, "energy")):
            problem = _axis_problem(name, axis0, axis, what)
            if problem:
                return f"“{name}” does not match “{name0}”: {problem}."
    return ""


def stack_axis(count: int, start: float, step: float) -> np.ndarray:
    """``start, start + step, ...`` -- one value per cut, in stacking order."""
    if not np.isfinite(step) or step == 0:
        raise ValueError("the step of the new axis must not be zero")
    return float(start) + float(step) * np.arange(int(count), dtype=float)


def stack(values) -> np.ndarray:
    """The cuts' arrays as one (cut, angle, energy) cube, in the order
    given."""
    return np.stack([np.asarray(v, dtype=np.float32) for v in values], axis=0)


def axis_label(name: str, unit: str) -> str:
    """``"Temperature (K)"``; just the name when there is no unit."""
    name = (name or "").strip() or "Scanned variable"
    unit = (unit or "").strip()
    return f"{name} ({unit})" if unit else name


def kind_for_role(role: str, angle_label: str = "") -> str:
    """What the stacked cube is. A photon-energy series of *angle* cuts is a
    ``kz_map``, which is what the k_z conversion takes; anything else --
    including a photon-energy series of cuts already in momentum, which the
    k_z conversion would read as angles -- is a ``map``."""
    if role == "photon_energy" and axis_unit(angle_label) != "invA":
        return "kz_map"
    return "map"
