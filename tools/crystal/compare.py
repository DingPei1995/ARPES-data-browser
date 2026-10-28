"""
tools/crystal/compare.py
========================
Putting a calculation and a measurement on the same axes.

A band-structure figure shows the calculation along a path of
high-symmetry points (Γ-X-M-Γ-Z-R-A-Z), and the measurement under it
along the *same* path. A measurement is not taken along such a path: it is
a set of maps -- constant-k_z planes from photon energies chosen for the
zone centre and the zone top, and a photon-energy scan for the k_z
direction -- so the measured picture along the path is assembled from
them:

* :class:`PlaneMap` -- a converted map (k_x, k_y, E) measured on one k_z
  plane; it serves the segments lying in that plane (Γ-X-M-Γ on the Γ
  plane, Z-R-A-Z on the zone-top plane);
* :class:`KzMap` -- a converted photon-energy scan (k_z, k_parallel, E)
  along one in-plane direction; it serves the segments along k_z (Γ-Z).

:func:`experiment_along_path` samples each path point from whichever map
holds it, using the zone's symmetry when the map covers only part of the
plane (a point outside the map is looked up at an equivalent point inside
it). The calculation along the same path is
:func:`tools.crystal.bands.named_path` and the model's ``energies``; the
two then share their x axis, the path length.

:func:`band_points` traces a measured cut's bands (ridge detection) for
:func:`tools.crystal.bands.fit_alignment`, the Z and shift that put the
calculation on them.

No Qt here; the window is :mod:`ui.panels.bands`.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["PlaneMap", "KzMap", "in_plane_operations", "experiment_along_path",
           "band_points", "band_segments", "symmetry_fold"]


@dataclass(eq=False)
class PlaneMap:
    """A map on one k_z plane: ``values[x, y, E]`` on ``x`` (k_x), ``y``
    (k_y) in Å⁻¹ and ``energy`` (eV, E_F = 0). ``rotation_deg`` is the
    crystal's a axis measured from the map's k_x (the sample azimuth),
    ``kz`` the plane's k_z (Å⁻¹)."""
    values: np.ndarray
    x: np.ndarray
    y: np.ndarray
    energy: np.ndarray
    kz: float = 0.0
    rotation_deg: float = 0.0
    name: str = ""


@dataclass(eq=False)
class KzMap:
    """A photon-energy scan converted to ``values[kz, k_par, E]``, taken
    along ``direction_deg`` (from the crystal's a axis) in the plane."""
    values: np.ndarray
    kz: np.ndarray
    kpar: np.ndarray
    energy: np.ndarray
    direction_deg: float = 0.0
    name: str = ""


def in_plane_operations(fold: int, mirror: bool = True):
    """The 2 x 2 matrices of an n-fold rotation axis along k_z (and, with
    ``mirror``, the mirror lines through it): identity first."""
    fold = max(1, int(fold))
    ops = []
    for n in range(fold):
        c, s = np.cos(2 * np.pi * n / fold), np.sin(2 * np.pi * n / fold)
        ops.append(np.array([[c, -s], [s, c]]))
    if mirror:
        m = np.array([[1.0, 0.0], [0.0, -1.0]])
        ops += [op @ m for op in ops]
    return ops


def symmetry_fold(space_group: int) -> int:
    """The rotation order about the c axis of a space group's crystal
    family -- 6 hexagonal/trigonal (3 with the mirrors' 6 equivalents is
    enough to fill a plane), 4 tetragonal, 2 lower; cubic 4."""
    n = int(space_group or 0)
    if 143 <= n <= 194:
        return 6
    if 75 <= n <= 142 or n >= 195:
        return 4
    return 2


def _rotate(points, degrees):
    c, s = np.cos(np.deg2rad(degrees)), np.sin(np.deg2rad(degrees))
    return points @ np.array([[c, s], [-s, c]])


def _inside(axis, values):
    lo, hi = float(np.min(axis)), float(np.max(axis))
    return (values >= lo) & (values <= hi)


def _interpolator(values, axes):
    """A NaN-outside linear interpolator, on axes that may run downwards."""
    from scipy.interpolate import RegularGridInterpolator

    values = np.asarray(values, dtype=float)
    axes = [np.asarray(a, dtype=float) for a in axes]
    for dim, axis in enumerate(axes):
        if axis.size > 1 and axis[1] < axis[0]:
            axes[dim] = axis[::-1]
            values = np.flip(values, axis=dim)
    return RegularGridInterpolator(tuple(axes), values, bounds_error=False,
                                   fill_value=np.nan)


def experiment_along_path(k_path, energy, sources, *, fold: int = 1,
                          mirror: bool = True, kz_period: float = None,
                          kz_tolerance: float = None,
                          line_tolerance: float = 0.03):
    """The measured intensity ``(N, n_E)`` at the path's points ``k_path``
    (N, 3), on the energy axis ``energy``, taken from ``sources``
    (:class:`PlaneMap` / :class:`KzMap`). Also returns the index of the
    source each point came from (-1: none, the column is NaN).

    A point belongs to a plane map whose k_z equals its own within
    ``kz_tolerance`` (k_z counted modulo ``kz_period``, the reciprocal
    lattice's period along k_z, and up to sign); it is looked up at the
    first of its symmetry-equivalent points (``fold``-fold axis, mirrors)
    that the map covers. A point not in any plane comes from a k_z map
    whose line passes within ``line_tolerance`` Å⁻¹ of it.
    """
    k_path = np.atleast_2d(np.asarray(k_path, dtype=float))
    if k_path.shape[1] == 2:
        k_path = np.column_stack([k_path, np.zeros(len(k_path))])
    energy = np.asarray(energy, dtype=float)
    n = len(k_path)
    image = np.full((n, energy.size), np.nan)
    owner = np.full(n, -1, dtype=int)
    if kz_tolerance is None:
        kz_tolerance = 0.06 * kz_period if kz_period else 0.03
    ops = in_plane_operations(fold, mirror)

    def same_kz(a, b):
        if kz_period:
            d = np.abs(((a - b) + 0.5 * kz_period) % kz_period - 0.5 * kz_period)
            e = np.abs(((a + b) + 0.5 * kz_period) % kz_period - 0.5 * kz_period)
            return np.minimum(d, e) <= kz_tolerance
        return np.minimum(np.abs(a - b), np.abs(a + b)) <= kz_tolerance

    for index, source in enumerate(sources):
        if not isinstance(source, PlaneMap):
            continue
        free = (owner < 0) & same_kz(k_path[:, 2], float(source.kz))
        if not free.any():
            continue
        rows = np.flatnonzero(free)
        inplane = k_path[rows, :2]
        found = np.full((len(rows), 2), np.nan)
        for op in ops:
            moved = _rotate(inplane @ op.T, -source.rotation_deg)
            todo = np.isnan(found[:, 0])
            ok = todo & _inside(source.x, moved[:, 0]) & _inside(source.y, moved[:, 1])
            found[ok] = moved[ok]
        ok = ~np.isnan(found[:, 0])
        if not ok.any():
            continue
        rows, found = rows[ok], found[ok]
        look = _interpolator(source.values, (source.x, source.y, source.energy))
        grid = np.column_stack([np.repeat(found, energy.size, axis=0),
                                np.tile(energy, len(rows))])
        image[rows] = look(grid).reshape(len(rows), energy.size)
        owner[rows] = index

    for index, source in enumerate(sources):
        if not isinstance(source, KzMap):
            continue
        free = np.flatnonzero(owner < 0)
        if free.size == 0:
            break
        found = np.full((free.size, 2), np.nan)
        for op in ops:
            inplane = k_path[free, :2] @ op.T
            phi = np.deg2rad(source.direction_deg)
            along = inplane @ np.array([np.cos(phi), np.sin(phi)])
            across = inplane @ np.array([-np.sin(phi), np.cos(phi)])
            kz = k_path[free, 2]
            candidates = [kz, -kz]
            if kz_period:
                # the scan's k_z is absolute (several zones out): bring each
                # point to the periods the scan's range holds
                lo, hi = float(np.min(source.kz)), float(np.max(source.kz))
                first = int(np.floor(lo / kz_period)) - 1
                last = int(np.ceil(hi / kz_period)) + 1
                candidates = [sign * kz + n * kz_period for n in range(first, last + 1)
                              for sign in (1, -1)]
            for kz_try in candidates:
                todo = np.isnan(found[:, 0])
                ok = (todo & (np.abs(across) <= line_tolerance)
                      & _inside(source.kz, kz_try) & _inside(source.kpar, along))
                found[ok] = np.column_stack([kz_try, along])[ok]
        ok = ~np.isnan(found[:, 0])
        if not ok.any():
            continue
        rows, found = free[ok], found[ok]
        look = _interpolator(source.values, (source.kz, source.kpar, source.energy))
        grid = np.column_stack([np.repeat(found, energy.size, axis=0),
                                np.tile(energy, len(rows))])
        image[rows] = look(grid).reshape(len(rows), energy.size)
        owner[rows] = index
    return image, owner


def band_points(values, k_axis, energy, *, sigma_k: float = None,
                sigma_e: float = None, threshold: float = 0.15,
                energy_window=(-np.inf, np.inf), along: str = "both"):
    """Points ``(k, E)`` on a measured cut's bands, by ridge detection:
    the input of :func:`tools.crystal.bands.fit_alignment`.

    ``values[k, E]``. ``along`` "mdc" keeps ridges running across k
    (dispersive bands), "edc" those running across E (flat bands), "both"
    all. Points outside ``energy_window`` (eV) are dropped -- above E_F
    there is nothing to match.
    """
    from tools.processing import ridge as RD

    k_axis = np.asarray(k_axis, dtype=float)
    energy = np.asarray(energy, dtype=float)
    if sigma_k is None or sigma_e is None:
        dk, de = RD.default_sigmas(k_axis, energy, 2.0)
        sigma_k = dk if sigma_k is None else sigma_k
        sigma_e = de if sigma_e is None else sigma_e
    result = RD.ridges(values, k_axis, energy, sigma_k, sigma_e,
                       threshold=threshold)
    points = result.points
    if not len(points):
        return np.zeros(0), np.zeros(0)
    lo, hi = energy_window
    keep = (points[:, 1] >= lo) & (points[:, 1] <= hi)
    return points[keep, 0], points[keep, 1]


def band_segments(distance, energies, ticks=(), window=None):
    """The calculated bands as polylines ``[(band, xs, ys)]`` for an
    overlay: each band split where the path jumps (a tick label "X|U"),
    and, with ``window`` = (lo, hi) eV, clipped to it (NaN outside, so a
    line is not drawn across the frame's edge)."""
    distance = np.asarray(distance, dtype=float)
    energies = np.atleast_2d(np.asarray(energies, dtype=float))
    # a jump leaves its two points at the same distance: split between them
    edges = [0]
    for d, label in ticks:
        if "|" in str(label):
            same = np.flatnonzero((np.abs(distance[:-1] - d) < 1e-9)
                                  & (np.abs(distance[1:] - d) < 1e-9))
            if same.size:
                edges.append(int(same[0]) + 1)
    edges.append(distance.size)
    out = []
    for band, row in enumerate(energies):
        row = row.copy()
        if window is not None:
            lo, hi = window
            pad = 0.05 * (hi - lo)
            row[(row < lo - pad) | (row > hi + pad)] = np.nan
        for a, b in zip(edges[:-1], edges[1:]):
            ys = row[a:b]
            if np.isfinite(ys).sum() >= 2:
                out.append((band, distance[a:b], ys))
    return out
