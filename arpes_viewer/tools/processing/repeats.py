"""
tools/processing/repeats.py
===========================
One measurement taken in several batches, added back together.

A fine cut or a fine map is routinely acquired as many short scans rather
than one long one, so that a sample that ages under the beam (adsorbates,
a drifting surface potential, radiation damage) shows up as a scan that no
longer looks like the others -- and can be left out -- instead of being
smeared into the whole. The MATLAB tool for this (``add_two_value_arrays.m``)
added exactly two arrays, element by element, with no check that they were
the same measurement and no way to see which batch had changed.

Here, for any number of repeats of one kind (cuts, maps, kz maps, k-maps):

1. **Are they repeats?** The grids are compared axis by axis. The same grid
   point for point is summed as it is; a grid that differs (a kinetic-
   energy window shifted by a photon-energy drift, one more point) is
   interpolated onto the first scan's grid, over the range every scan
   covers. The recorded metadata -- photon energy, polarisation, angles,
   position, temperature, pass energy, lens mode -- are compared against
   the first scan (the cut-arithmetic table, ``cutops.COMPARED``).
2. **Has the sample changed?** Each scan is compared with the sum of all
   the *others*, with the counting statistics of both:

       chi2/N = mean over pixels of (a - s b)^2 / (a + s^2 b),
       s = sum(a) / sum(b),

   which is 1 for a scan that differs from the rest only by counting noise
   and grows with every real change in the spectrum (a band that moved or
   broadened, a background that grew). Comparing against the *others*
   rather than against the total keeps a changed scan from hiding in a
   reference it is part of. Where the counts are not Poisson (an MCP's
   gain spreads them) every scan's value is raised alike, so a scan is
   flagged by how far it stands above the median of the set, not above 1.
   The total count rate per scan is shown alongside -- a falling rate
   with an unchanged shape is beam or flux, not the sample.
3. **Has the energy scale drifted?** The angle-integrated EDC of each scan
   is cross-correlated with that of the others through their derivative
   (so the Fermi edge and the band edges, not the background, set it),
   with a parabola through the peak for the sub-step part; the others are
   first lined up by their own shifts (a few rounds), so a shifted scan does
   not drag the reference with it, and the shifts are given relative to the
   set's median. The shifts can
   be undone before summing; the energy range is then cropped to what
   every shifted scan still covers.

The result keeps counts as counts (a plain sum) unless the mean is asked
for, and records which scans went in, which were left out and why, the
shifts applied and each scan's chi2.

Qt-free, so it can be tested on its own.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage

__all__ = ["grid_plan", "energy_shift", "diagnose", "sum_repeats",
           "flag_changed", "Diagnosis"]

#: How far two axes may differ, as a fraction of a step, and still be the
#: same grid.
SAME_GRID = 0.02


def _step(axis) -> float:
    axis = np.asarray(axis, dtype=float)
    return float(np.median(np.abs(np.diff(axis)))) if axis.size > 1 else 1.0


def grid_plan(axes_list):
    """How the scans' grids relate.

    ``axes_list`` holds each scan's axes (a tuple per scan, one array per
    dimension). Returns ``{"exact": bool, "target": axes, "notes": [...]}``
    -- ``target`` the grid the sum is built on: the first scan's own when
    all agree, otherwise the first scan's points inside the range every
    scan covers. Raises ``ValueError`` when the scans cannot be repeats of
    one another (a different number of axes, or no common range).
    """
    first = [np.asarray(a, dtype=float) for a in axes_list[0]]
    exact = True
    notes = []
    target = []
    for dim, axis in enumerate(first):
        lo, hi = float(np.min(axis)), float(np.max(axis))
        step = _step(axis)
        same = True
        for i, axes in enumerate(axes_list[1:], start=2):
            if len(axes) != len(first):
                raise ValueError(f"scan {i} has {len(axes)} axes, the first "
                                 f"{len(first)}: not the same measurement")
            other = np.asarray(axes[dim], dtype=float)
            lo, hi = max(lo, float(np.min(other))), min(hi, float(np.max(other)))
            if (other.size != axis.size
                    or np.max(np.abs(other - axis)) > SAME_GRID * step):
                same = False
        if hi < lo:
            raise ValueError(f"the scans have no range of axis {dim + 1} in "
                             f"common")
        if same:
            target.append(axis)
        else:
            exact = False
            keep = (axis >= lo - 1e-9 * step) & (axis <= hi + 1e-9 * step)
            target.append(axis[keep])
            notes.append(f"axis {dim + 1} differs between scans: "
                         f"interpolated onto the first scan's points from "
                         f"{lo:.6g} to {hi:.6g}")
    return {"exact": exact, "target": tuple(target), "notes": notes}


def _fractional_index(axis, values):
    axis = np.asarray(axis, dtype=float)
    order = np.argsort(axis)
    return np.interp(values, axis[order], order.astype(float))


def onto_grid(values, from_axes, to_axes, order: int = 1):
    """``values`` on ``from_axes`` resampled onto ``to_axes`` (each
    monotonic), linearly; points outside become NaN."""
    values = np.asarray(values, dtype=float)
    same = all(np.asarray(f).size == np.asarray(t).size
               and np.allclose(f, t) for f, t in zip(from_axes, to_axes))
    if same:
        return values
    coords = np.meshgrid(*[_fractional_index(f, np.asarray(t, dtype=float))
                           for f, t in zip(from_axes, to_axes)],
                         indexing="ij")
    out = ndimage.map_coordinates(values, coords, order=order, mode="nearest")
    for dim, (f, t) in enumerate(zip(from_axes, to_axes)):
        t = np.asarray(t, dtype=float)
        outside = (t < np.min(f) - 1e-9) | (t > np.max(f) + 1e-9)
        if outside.any():
            index = [slice(None)] * out.ndim
            index[dim] = outside
            out[tuple(index)] = np.nan
    return out


def edc_profile(values, energy_dim: int):
    """The angle-integrated EDC: everything summed but the energy axis."""
    values = np.asarray(values, dtype=float)
    other = tuple(d for d in range(values.ndim) if d != energy_dim)
    return np.nansum(values, axis=other)


def energy_shift(reference, profile, energy, max_fraction: float = 0.15,
                 min_correlation: float = 0.3):
    """How far ``profile`` sits from ``reference`` along ``energy``, in the
    axis' units (positive: ``profile`` is at higher energy).

    Cross-correlation of the smoothed derivatives, both normalised, over
    shifts up to ``max_fraction`` of the range; a parabola through the
    three points around the maximum gives the part below one step. When
    the two share no structure (the best correlation is below
    ``min_correlation`` -- a flat EDC, or noise) there is nothing to line
    up and the shift is 0.
    """
    energy = np.asarray(energy, dtype=float)
    step = _step(energy)
    n = energy.size
    if n < 8:
        return 0.0

    def prepared(curve):
        curve = ndimage.gaussian_filter1d(np.nan_to_num(np.asarray(curve, float)), 1.5)
        slope = np.gradient(curve)
        slope -= slope.mean()
        norm = np.sqrt(np.sum(slope ** 2))
        return slope / norm if norm > 0 else slope

    a, b = prepared(reference), prepared(profile)
    reach = max(2, int(max_fraction * n))
    lags = np.arange(-reach, reach + 1)
    scores = np.array([np.sum(a[max(0, -lag):n - max(0, lag)]
                              * b[max(0, lag):n - max(0, -lag)])
                       for lag in lags])
    best = int(np.argmax(scores))
    # Two unrelated noisy curves correlate by about 1/sqrt(independent
    # points) -- after the smoothing, about a third of the points -- so a
    # short curve needs a clearer peak than a long one.
    if scores[best] < max(min_correlation, 3.0 / np.sqrt(max(n / 3.0, 1.0))):
        return 0.0
    lag = float(lags[best])
    if 0 < best < lags.size - 1:
        y0, y1, y2 = scores[best - 1:best + 2]
        curve = y0 - 2 * y1 + y2
        if curve < 0:
            lag += 0.5 * (y0 - y2) / curve
    direction = 1.0 if energy[-1] >= energy[0] else -1.0
    return direction * lag * step


def _moved(profile, energy, shift):
    """``profile`` moved back by ``shift`` (energy units): sampled at
    E + shift."""
    order = np.argsort(energy)
    return np.interp(energy + shift, energy[order], np.asarray(profile)[order],
                     left=np.nan, right=np.nan)


def _aligned_shifts(profiles, energy, rounds: int = 3):
    """Each EDC's shift against the others, the others first brought into
    line with their own current estimates -- so a shifted scan does not
    drag the reference along with it, and the shifts are relative to the
    set's median scan."""
    energy = np.asarray(energy, dtype=float)
    count = len(profiles)
    shifts = np.zeros(count)
    for _round in range(rounds):
        aligned = [np.nan_to_num(_moved(p, energy, s))
                   for p, s in zip(profiles, shifts)]
        total = np.sum(aligned, axis=0)
        shifts = np.array([energy_shift(total - aligned[i], profiles[i], energy)
                           for i in range(count)])
        shifts -= np.median(shifts)
    return shifts


class Diagnosis:
    """What :func:`diagnose` found, one entry per scan."""

    def __init__(self, totals, chi2, shifts, flags, notes):
        self.totals = np.asarray(totals, dtype=float)
        self.chi2 = np.asarray(chi2, dtype=float)
        self.shifts = np.asarray(shifts, dtype=float)
        self.flags = list(flags)
        self.notes = list(notes)

    @property
    def relative_totals(self):
        mean = np.mean(self.totals) if self.totals.size else 1.0
        return self.totals / mean if mean else self.totals


def flag_changed(chi2, factor: float = 1.5, floor: float = 1.2):
    """Which scans stand out: chi2 above ``factor`` times the median of
    the set, and above ``floor`` (so a set of identical noise-only scans
    flags nothing)."""
    chi2 = np.asarray(chi2, dtype=float)
    finite = chi2[np.isfinite(chi2)]
    median = float(np.median(finite)) if finite.size else 1.0
    return [bool(np.isfinite(c) and c > factor * median and c > floor)
            for c in chi2]


def diagnose(get, count: int, axes_list, energy_dim: int, plan=None,
             progress=None) -> Diagnosis:
    """Total counts, chi2 against the others, and energy shift per scan.

    ``get(i)`` returns scan ``i``'s values (read when asked, so a set of
    large maps is never all in memory at once); ``axes_list`` its axes.
    """
    plan = plan or grid_plan(axes_list)
    target = plan["target"]

    def onto(i):
        return onto_grid(get(i), axes_list[i], target)

    total = None
    totals, profiles = [], []
    for i in range(count):
        values = onto(i)
        totals.append(float(np.nansum(values)))
        profiles.append(edc_profile(values, energy_dim))
        total = values.copy() if total is None else total + values
        if progress:
            progress(0.5 * (i + 1) / count, f"Reading scan {i + 1} of {count}")
    energy = target[energy_dim]
    shifts = _aligned_shifts(profiles, energy) if count > 1 else [0.0] * count
    chi2 = []
    for i in range(count):
        a = onto(i)
        b = total - a
        good = np.isfinite(a) & np.isfinite(b) & ((a + b) > 0)
        sa, sb = float(np.sum(a[good])), float(np.sum(b[good]))
        if count < 2 or sb <= 0 or not good.any():
            chi2.append(np.nan)
        else:
            s = sa / sb
            var = np.clip(a[good], 0, None) + s * s * np.clip(b[good], 0, None)
            ok = var > 0
            chi2.append(float(np.mean((a[good][ok] - s * b[good][ok]) ** 2
                                      / var[ok])))
        if progress:
            progress(0.5 + 0.5 * (i + 1) / count, f"Comparing scan {i + 1}")
    flags = flag_changed(chi2)
    notes = list(plan["notes"])
    if count >= 3 and totals[0] > 0:
        slope = np.polyfit(np.arange(count), np.asarray(totals) / np.mean(totals), 1)[0]
        if abs(slope) * (count - 1) > 0.1:
            notes.append(f"the count rate {'falls' if slope < 0 else 'rises'} "
                         f"by {abs(slope) * (count - 1) * 100:.0f} % from the "
                         f"first scan to the last")
    return Diagnosis(totals, chi2, shifts, flags, notes)


def sum_repeats(get, axes_list, include, energy_dim: int, *, shifts=None,
                mode: str = "sum", plan=None):
    """The included scans added up.

    Returns ``(axes, values, notes)``. With ``shifts`` (one per scan, in
    energy units, as :func:`energy_shift` measures them) every scan is moved
    back by its own shift first and the energy axis cropped to what all of
    them still cover. ``mode`` is "sum" (counts stay counts) or "mean".
    """
    plan = plan or grid_plan([axes_list[i] for i in include])
    target = list(plan["target"])
    notes = list(plan["notes"])
    include = list(include)
    if not include:
        raise ValueError("no scan is included")
    energy = np.asarray(target[energy_dim], dtype=float)
    step = _step(energy)
    sign = 1.0 if energy.size < 2 or energy[-1] >= energy[0] else -1.0
    used = np.zeros(len(axes_list)) if shifts is None else np.asarray(shifts, float)
    used = used - (np.mean(used[include]) if shifts is not None else 0.0)
    total = None
    for i in include:
        values = onto_grid(get(i), axes_list[i], target)
        if shifts is not None and abs(used[i]) > 1e-3 * step:
            # Moving a scan down by its shift: sample it at E + shift.
            offset = [0.0] * values.ndim
            offset[energy_dim] = sign * used[i] / step
            values = ndimage.shift(np.nan_to_num(values), [-o for o in offset],
                                   order=1, mode="constant", cval=np.nan)
        total = values if total is None else total + values
    if shifts is not None and np.any(np.abs(used[include]) > 1e-3 * step):
        reach = int(np.ceil(np.max(np.abs(used[include])) / step))
        if reach * 2 >= energy.size:
            raise ValueError("the energy shifts are larger than the range")
        keep = slice(reach, energy.size - reach)
        index = [slice(None)] * total.ndim
        index[energy_dim] = keep
        total = total[tuple(index)]
        target[energy_dim] = energy[keep]
        notes.append(f"energy shifts undone (up to {np.max(np.abs(used[include])) * 1000:.1f} "
                     f"meV); {reach} points cropped at each end of the energy axis")
    if mode == "mean":
        total = total / len(include)
    return tuple(target), np.nan_to_num(total), notes
