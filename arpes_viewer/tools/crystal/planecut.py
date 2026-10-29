"""
tools/crystal/planecut.py
=========================
A calculation cut by any plane of k-space, the way the Brillouin-zone tool
cuts the zone: a surface normal given as Miller indices (hkl), an offset
along it (stepped through the planes that hold a high-symmetry point), and
a k_z window to integrate over -- what an ARPES map, with its finite k_z
resolution, really sees.

Everything is in the calculation's own Cartesian frame. The Miller indices
are read against the conventional reciprocal cell, as the Brillouin-zone
tool reads them; the named points are that tool's (``highsym``), brought
into the calculation's frame through the fractional coordinates the two
share. No Qt.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from arpes_viewer.tools.crystal import bands as B

__all__ = ["Plane", "to_model_frame", "plane_of", "named_points", "through_points",
           "sample", "contours", "spectral_map", "zone_outlines"]


@dataclass
class Plane:
    """The plane ``k . normal = offset``, with in-plane axes ``e1, e2``
    (orthonormal; ``e1 x e2 = normal``), drawn over ``-reach..reach``."""
    normal: np.ndarray
    e1: np.ndarray
    e2: np.ndarray
    offset: float = 0.0
    reach: float = 1.0

    def points(self, u, v, shift: float = 0.0) -> np.ndarray:
        """3-D k of in-plane coordinates ``u``, ``v`` (any shape, the same),
        moved ``shift`` along the normal."""
        u, v = np.asarray(u, float), np.asarray(v, float)
        return ((self.offset + shift) * self.normal + u[..., None] * self.e1
                + v[..., None] * self.e2)


def to_model_frame(params, reciprocal) -> np.ndarray:
    """The matrix taking a Cartesian k of the crystal's standard frame
    (``params``) into the calculation's (rows ``reciprocal``): through the
    fractional coordinates both share."""
    from arpes_viewer.tools.crystal.lattice import primitive_vectors, reciprocal_vectors

    b_std = reciprocal_vectors(primitive_vectors(params))
    return np.linalg.inv(b_std) @ np.asarray(reciprocal, float)


def plane_of(params, reciprocal, hkl, offset: float = 0.0,
             reach: float = None, rotation_deg: float = 0.0) -> Plane:
    """The plane with normal (hkl) -- conventional reciprocal axes, as in
    the Brillouin-zone tool -- in the calculation's frame. Its in-plane
    axes are the zone tool's (:func:`tools.crystal.bz3d.plane_basis`), so
    (001) gives the familiar k_x, k_y -- turned by ``rotation_deg`` about
    the normal (counter-clockwise as drawn), to match a map's azimuth."""
    from arpes_viewer.tools.crystal.bz3d import plane_basis
    from arpes_viewer.tools.crystal.lattice import conventional_vectors, reciprocal_vectors

    convert = to_model_frame(params, reciprocal)
    normal = np.asarray(hkl, float) @ reciprocal_vectors(conventional_vectors(params)) @ convert
    length = float(np.linalg.norm(normal))
    if length < 1e-9:
        raise ValueError("the plane's normal (h, k, l) cannot be (0, 0, 0)")
    normal /= length
    e1, e2 = plane_basis(normal)
    turn = np.deg2rad(float(rotation_deg))
    e1, e2 = (np.cos(turn) * e1 + np.sin(turn) * e2,
              -np.sin(turn) * e1 + np.cos(turn) * e2)
    if reach is None:
        reach = 0.8 * float(np.linalg.norm(reciprocal, axis=1).max())
    return Plane(normal, e1, e2, float(offset), float(reach))


def named_points(params, reciprocal):
    """Every copy of the zone's named points in the first zone, in the
    calculation's frame: ``[HighSymmetryPoint]``."""
    from arpes_viewer.tools.crystal.highsym import HighSymmetryPoint, equivalent_points

    convert = to_model_frame(params, reciprocal)
    return [HighSymmetryPoint(p.label, np.asarray(p.k) @ convert)
            for p in equivalent_points(params)]


def through_points(points, normal):
    """The offsets at which a plane of this ``normal`` holds a named point:
    ``[(offset, [labels])]``, from 0 up (the Brillouin-zone tool's "Through
    point" list)."""
    from arpes_viewer.tools.crystal.highsym import plane_offsets

    return plane_offsets(points, normal)


def sample(model, plane: Plane, n: int = 161, kz_window: float = 0.0,
           slices: int = 7, progress=None):
    """The bands on the plane: ``(axis, energies)`` with ``axis`` the ``n``
    in-plane coordinates (from -reach to reach) and ``energies`` of shape
    ``(slices, bands, n, n)`` relative to E_F -- one plane, or ``slices``
    planes spread evenly over ``kz_window`` (Å⁻¹) about it."""
    axis = np.linspace(-plane.reach, plane.reach, n)
    gu, gv = np.meshgrid(axis, axis, indexing="ij")
    shifts = (np.linspace(-kz_window / 2, kz_window / 2, max(int(slices), 2))
              if kz_window > 0 else np.zeros(1))
    out = np.empty((len(shifts), model.band_count(), n, n))
    step = 2000
    for s, shift in enumerate(shifts):
        k = plane.points(gu, gv, shift).reshape(-1, 3)
        flat = np.empty((model.band_count(), len(k)))
        for start in range(0, len(k), step):
            flat[:, start:start + step] = model.energies(k[start:start + step])
            if progress:
                progress((s + (start + step) / len(k)) / len(shifts),
                         f"The bands on the plane ({s + 1} of {len(shifts)})")
        out[s] = flat.reshape(-1, n, n)
    return axis, out


def contours(energies, axis, energy: float):
    """The contours at ``energy`` of the middle plane: ``[(band, xs, ys)]``."""
    middle = energies[len(energies) // 2]
    out = []
    for band, grid in enumerate(middle):
        if grid.min() <= energy <= grid.max():
            for xs, ys in B.contour_lines(grid, axis, axis, energy):
                out.append((band, xs, ys))
    return out


def spectral_map(energies, energy: float, width: float = 0.02):
    """The calculated intensity at ``energy`` -- Lorentzians of half-width
    ``width`` (eV) about every band -- averaged over the planes of the k_z
    window: ``(n, n)``, what a map integrating over k_z would see."""
    weight = width / ((energies - energy) ** 2 + width ** 2) / np.pi
    return weight.sum(axis=1).mean(axis=0)


def zone_outlines(reciprocal, plane: Plane):
    """Where the zone and its neighbours cut the plane: ``[(xs, ys)]`` in
    the plane's (e1, e2) coordinates (the zone tool's cut)."""
    from arpes_viewer.tools.crystal.bz3d import cuts_in_frame, wigner_seitz_cell

    zone = wigner_seitz_cell(np.asarray(reciprocal, float))
    window = (-plane.reach, plane.reach)
    return cuts_in_frame(zone, reciprocal, plane.offset * plane.normal,
                         plane.e1, plane.e2, window, window)
