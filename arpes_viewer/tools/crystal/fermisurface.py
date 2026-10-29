"""
tools/crystal/fermisurface.py
=============================
The 3-D Fermi surface (or any constant-energy surface) of a calculation
that knows its bands everywhere -- a band grid or a Wannier Hamiltonian --
as triangles inside the first Brillouin zone, ready to be drawn.

How: the bands are sampled on a Cartesian grid over the zone's bounding
box (a Hamiltonian first on a fractional grid it is then interpolated
from, since diagonalising at every point of a fine grid is slow), each band
that crosses the energy is meshed by marching tetrahedra, and the triangles
outside the zone are dropped -- what the MATLAB scripts of the theory side
do with ``isosurface`` and ``BZconfine``. Each triangle carries the normal
of the band's gradient at its centre, for smooth shading.

Marching *tetrahedra* rather than cubes: no 256-case table, every case is
one or two triangles, and it vectorises over the whole grid in a few
lines. It makes about twice the triangles of marching cubes, which the
drawing copes with.

No Qt.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import ndimage

from arpes_viewer.tools.crystal import bands as B

__all__ = ["Surface", "FermiSurface", "marching_tetrahedra", "fermi_surface",
           "inside_zone", "zone_distance", "clip_to_zone", "crossing_bands",
           "sampled"]

#: the six tetrahedra of a cube, all sharing its diagonal 0-6; corners
#: numbered (0,0,0) (1,0,0) (1,1,0) (0,1,0) (0,0,1) (1,0,1) (1,1,1) (0,1,1)
_CORNERS = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                     [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]])
_TETRAHEDRA = ((0, 1, 2, 6), (0, 2, 3, 6), (0, 3, 7, 6), (0, 7, 4, 6),
               (0, 4, 5, 6), (0, 5, 1, 6))


@dataclass
class Surface:
    """One band's surface: ``triangles`` (M, 3, 3) Cartesian (Å⁻¹) and
    ``normals`` (M, 3), unit."""
    band: int                       # index in the calculation (0-based)
    number: int                     # as the file counts it (1-based)
    triangles: np.ndarray
    normals: np.ndarray


@dataclass
class FermiSurface:
    energy: float
    surfaces: list = field(default_factory=list)
    zone: object = None             # tools.crystal.bz3d.BrillouinZone
    reciprocal: np.ndarray = None


def marching_tetrahedra(values, level: float):
    """The surface ``values == level`` of a 3-D array, as triangles in
    index coordinates: ``(M, 3, 3)``."""
    values = np.asarray(values, dtype=float)
    n1, n2, n3 = values.shape
    if min(n1, n2, n3) < 2:
        return np.zeros((0, 3, 3))
    idx = np.stack(np.meshgrid(np.arange(n1 - 1), np.arange(n2 - 1),
                               np.arange(n3 - 1), indexing="ij"), -1).reshape(-1, 3)
    # keep the cubes the level passes through
    corner_values = np.stack([values[idx[:, 0] + c[0], idx[:, 1] + c[1], idx[:, 2] + c[2]]
                              for c in _CORNERS], 1)                 # (N, 8)
    crossed = (corner_values.min(1) < level) & (corner_values.max(1) >= level)
    idx, corner_values = idx[crossed], corner_values[crossed]
    out = []

    def edge(points_a, points_b, va, vb):
        t = (level - va) / np.where(vb != va, vb - va, 1.0)
        return points_a + t[:, None] * (points_b - points_a)

    for tet in _TETRAHEDRA:
        v = corner_values[:, tet]                                    # (N, 4)
        p = idx[:, None, :] + _CORNERS[list(tet)][None]              # (N, 4, 3)
        below = v < level
        count = below.sum(1)
        # one corner apart from the other three: a triangle
        for lone in range(4):
            others = [i for i in range(4) if i != lone]
            mask = ((count == 1) & below[:, lone]) | ((count == 3) & ~below[:, lone])
            if not mask.any():
                continue
            pm, vm = p[mask], v[mask]
            out.append(np.stack([edge(pm[:, lone], pm[:, o], vm[:, lone], vm[:, o])
                                 for o in others], 1))
        # two and two: a quadrilateral, as two triangles
        for a, b in ((0, 1), (0, 2), (0, 3)):
            c, d = [i for i in range(4) if i not in (a, b)]
            mask = (count == 2) & (below[:, a] == below[:, b])
            if not mask.any():
                continue
            pm, vm = p[mask], v[mask]
            ac = edge(pm[:, a], pm[:, c], vm[:, a], vm[:, c])
            ad = edge(pm[:, a], pm[:, d], vm[:, a], vm[:, d])
            bd = edge(pm[:, b], pm[:, d], vm[:, b], vm[:, d])
            bc = edge(pm[:, b], pm[:, c], vm[:, b], vm[:, c])
            out.append(np.stack([ac, ad, bd], 1))
            out.append(np.stack([ac, bd, bc], 1))
    return np.concatenate(out) if out else np.zeros((0, 3, 3))


def _neighbours(reciprocal):
    b = np.asarray(reciprocal, float)
    g = [i * b[0] + j * b[1] + k * b[2] for i in (-1, 0, 1) for j in (-1, 0, 1)
         for k in (-1, 0, 1) if (i, j, k) != (0, 0, 0)]
    return np.array(g)


def inside_zone(points, reciprocal, tol: float = 1e-6):
    """Which of ``points`` (N, 3) lie in the first Brillouin zone (nearer Γ
    than any other reciprocal lattice point)."""
    g = _neighbours(reciprocal)
    points = np.asarray(points, float)
    return np.all(points @ g.T <= 0.5 * np.sum(g * g, 1) * (1 + tol) + tol, axis=1)


def zone_distance(points, reciprocal):
    """How far inside the first zone each of ``points`` lies (Å⁻¹;
    negative outside): the distance to the nearest of its bounding
    planes."""
    g = _neighbours(reciprocal)
    norms = np.linalg.norm(g, axis=1)
    points = np.asarray(points, float)
    return np.min((0.5 * norms ** 2 - points @ g.T) / norms, axis=-1)


def clip_to_zone(triangles, reciprocal):
    """``triangles`` (M, 3, 3) cut at the zone boundary: the parts inside,
    as triangles. Each is cut where the distance to the boundary (linear
    along an edge, short of the zone's corners) changes sign, so the
    surface ends on the zone's faces instead of in a ragged edge."""
    if not len(triangles):
        return triangles, np.zeros(0, dtype=int)
    d = zone_distance(triangles.reshape(-1, 3), reciprocal).reshape(-1, 3)
    inside = d >= 0
    count = inside.sum(1)
    pieces, owners = [triangles[count == 3]], [np.flatnonzero(count == 3)]

    def cut(p, q, dp, dq):
        t = dp / np.where(dp != dq, dp - dq, 1.0)
        return p + t[:, None] * (q - p)

    for lone in range(3):
        a, b = [i for i in range(3) if i != lone]
        # only the lone corner inside: the small triangle at it
        mask = (count == 1) & inside[:, lone]
        if mask.any():
            t, dd = triangles[mask], d[mask]
            pieces.append(np.stack([t[:, lone],
                                    cut(t[:, lone], t[:, a], dd[:, lone], dd[:, a]),
                                    cut(t[:, lone], t[:, b], dd[:, lone], dd[:, b])], 1))
            owners.append(np.flatnonzero(mask))
        # only the lone corner outside: the quadrilateral without it
        mask = (count == 2) & ~inside[:, lone]
        if mask.any():
            t, dd = triangles[mask], d[mask]
            pa = cut(t[:, a], t[:, lone], dd[:, a], dd[:, lone])
            pb = cut(t[:, b], t[:, lone], dd[:, b], dd[:, lone])
            pieces.append(np.stack([t[:, a], t[:, b], pb], 1))
            pieces.append(np.stack([t[:, a], pb, pa], 1))
            owners += [np.flatnonzero(mask)] * 2
    return np.concatenate(pieces), np.concatenate(owners)


def sampled(model, n: int = 28, progress=None):
    """A :class:`~tools.crystal.bands.BandGrid` for ``model``: itself, or a
    Hamiltonian sampled on an ``n``³ fractional grid."""
    if isinstance(model, B.BandGrid):
        return model
    frac = np.stack(np.meshgrid(*(np.arange(n) / n,) * 3, indexing="ij"), -1).reshape(-1, 3)
    k = frac @ model.reciprocal
    energies = np.empty((model.band_count(), len(k)))
    step = 400
    for start in range(0, len(k), step):
        energies[:, start:start + step] = model.energies(k[start:start + step]) + model.fermi
        if progress:
            progress(0.7 * min(1.0, (start + step) / len(k)),
                     f"Diagonalising the Hamiltonian on a {n}³ grid")
    grid = B.BandGrid(model.reciprocal, energies.reshape(-1, n, n, n),
                      fermi=model.fermi, name=getattr(model, "name", ""))
    return grid


def crossing_bands(model, energy: float, n: int = 24, progress=None):
    """The indices of the bands that reach ``energy`` (eV, relative to
    E_F) somewhere in the zone."""
    grid = sampled(model, n, progress)
    rel = grid.values - grid.fermi
    flat = rel.reshape(len(rel), -1)
    return [int(i) for i in np.flatnonzero((flat.min(1) <= energy) & (flat.max(1) >= energy))]


def fermi_surface(model, energy: float = 0.0, *, n: int = 48, bands=None,
                  sample: int = 28, progress=None) -> FermiSurface:
    """The surfaces at ``energy`` (eV, relative to E_F) of ``bands`` (all
    that cross it when None), inside the first zone, on a Cartesian grid
    of about ``n`` points along its longest side; a Hamiltonian is first
    sampled on a ``sample``³ fractional grid."""
    from arpes_viewer.tools.crystal.bz3d import wigner_seitz_cell

    grid = sampled(model, sample, progress)
    b = np.asarray(grid.reciprocal, float)
    zone = wigner_seitz_cell(b)
    lo = zone.vertices.min(0) * 1.02
    hi = zone.vertices.max(0) * 1.02
    span = hi - lo
    counts = np.maximum(8, np.round(n * span / span.max()).astype(int))
    axes = [np.linspace(lo[i], hi[i], counts[i]) for i in range(3)]
    points = np.stack(np.meshgrid(*axes, indexing="ij"), -1).reshape(-1, 3)
    rel = grid.values - grid.fermi
    if bands is None:
        flat = rel.reshape(len(rel), -1)
        bands = [int(i) for i in np.flatnonzero((flat.min(1) <= energy)
                                                & (flat.max(1) >= energy))]
    numbers = getattr(model, "band_numbers", None)
    result = FermiSurface(energy=float(energy), zone=zone, reciprocal=b)
    spacing = span / (counts - 1)
    for count, band in enumerate(bands):
        if progress:
            progress(0.7 + 0.3 * count / max(len(bands), 1),
                     f"Meshing band {band + 1} ({count + 1} of {len(bands)})")
        single = B.BandGrid(b, grid.values[band:band + 1], origin=grid.origin,
                            fermi=grid.fermi)
        cube = single.energies(points)[0].reshape(*counts)
        triangles = marching_tetrahedra(cube, energy)
        if not len(triangles):
            continue
        triangles, owner = clip_to_zone(lo + triangles * spacing, b)
        if not len(triangles):
            continue
        centres_index = ((triangles.mean(1) - lo) / spacing)
        gradient = [ndimage.map_coordinates(np.gradient(cube, axis=i) / spacing[i],
                                            centres_index.T, order=1, mode="nearest")
                    for i in range(3)]
        normals = np.stack(gradient, 1)
        flat = np.cross(triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0])
        weak = np.linalg.norm(normals, axis=1) < 1e-12
        normals[weak] = flat[weak]
        normals /= np.maximum(np.linalg.norm(normals, axis=1, keepdims=True), 1e-300)
        number = int(numbers[band]) if numbers is not None else band + 1
        result.surfaces.append(Surface(band, number, triangles, normals))
    if progress:
        progress(1.0, "Done")
    return result
