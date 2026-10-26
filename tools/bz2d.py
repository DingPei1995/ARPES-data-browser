"""
tools/bz2d.py
===========
The 2-D Wigner-Seitz cell of a plane lattice, and tiling copies of it across
a field of view. Used directly by the Brillouin-zone tool's "2D layer"
mode (a monolayer/flat lattice given as a, b, gamma rather than a 3-D space
group), by :func:`tools.bz3d.plane_cut` indirectly through nothing (that
module has its own, unrelated 3-D routine), and by :mod:`tools.moire`, which
builds a moire *reciprocal* lattice and hands it here for its own
Wigner-Seitz cell.

Relation to the lab's MATLAB tools
-----------------------------------
This generalizes ``plotBZ_demo_old`` (in ``MiniBZ_plotter.m``): same idea
-- take a small tiling of reciprocal-lattice points around the origin and
compute the origin point's Voronoi cell -- but that version hardcodes
``gamma = 120`` (hexagonal only). Here ``g1``/``g2`` are arbitrary, so a
square, rectangular, centred-rectangular or oblique plane lattice works
exactly the same way as a hexagonal one.
"""
from __future__ import annotations

import numpy as np
from scipy.spatial import Voronoi

__all__ = ["reciprocal_vectors_2d", "wigner_seitz_cell_2d", "tile_2d",
           "lattice_point_group_2d", "irreducible_cell_2d"]


def reciprocal_vectors_2d(a1: np.ndarray, a2: np.ndarray):
    """2-D reciprocal vectors (physicist's convention, 2*pi included) for
    the plane lattice spanned by real-space ``a1``, ``a2``.

    Uses the standard 3-D cross-product construction with both vectors
    embedded in the z = 0 plane, which is exactly equivalent to (and
    simpler than re-deriving) the 2x2 matrix-inverse formula.
    """
    a1_3d = np.array([a1[0], a1[1], 0.0])
    a2_3d = np.array([a2[0], a2[1], 0.0])
    z = np.array([0.0, 0.0, 1.0])
    area = float(np.cross(a1_3d, a2_3d)[2])
    if abs(area) < 1e-12:
        raise ValueError("a1 and a2 are parallel (degenerate 2-D lattice)")
    b1_3d = 2.0 * np.pi * np.cross(z, a2_3d) / area
    b2_3d = 2.0 * np.pi * np.cross(a1_3d, z) / area
    return b1_3d[:2], b2_3d[:2]


def wigner_seitz_cell_2d(g1: np.ndarray, g2: np.ndarray, shell: int = 2) -> np.ndarray:
    """The origin's Voronoi cell in the lattice spanned by ``g1``, ``g2``:
    the Wigner-Seitz cell of a 2-D reciprocal lattice, i.e. the 2-D
    Brillouin zone. Returns a closed ``(n+1, 2)`` array of vertices.

    ``shell`` neighbour rings (``-shell..shell`` in each lattice index) are
    given to :class:`scipy.spatial.Voronoi`; grown automatically if the
    origin's region turns out unbounded (touches a Voronoi ridge at
    infinity), which only happens for an unreasonably small ``shell``.
    """
    g1, g2 = np.asarray(g1, float), np.asarray(g2, float)
    for attempt_shell in range(shell, shell + 4):
        idx = np.arange(-attempt_shell, attempt_shell + 1)
        I, J = np.meshgrid(idx, idx, indexing="ij")
        I, J = I.ravel(), J.ravel()
        points = I[:, None] * g1 + J[:, None] * g2
        center = np.argmin(np.einsum("ij,ij->i", points, points))  # the (0,0) point

        vor = Voronoi(points)
        region = vor.regions[vor.point_region[center]]
        if -1 in region or len(region) < 3:
            continue  # unbounded region: not enough neighbours yet
        verts = vor.vertices[region]
        centroid = verts.mean(axis=0)
        angles = np.arctan2(verts[:, 1] - centroid[1], verts[:, 0] - centroid[0])
        verts = verts[np.argsort(angles)]
        return np.vstack([verts, verts[:1]])
    raise RuntimeError(
        f"2-D Wigner-Seitz cell did not converge within {shell + 3} "
        "neighbour shells -- check g1, g2 for degeneracy.")


def _cross2(u, v) -> float:
    """The scalar cross product of two plane vectors.

    Spelled out rather than calling :func:`numpy.cross`, which deprecated
    2-component input in NumPy 2.0 and will eventually refuse it.
    """
    return float(u[0] * v[1] - u[1] * v[0])


def lattice_point_group_2d(g1: np.ndarray, g2: np.ndarray,
                           max_index: int = 3) -> np.ndarray:
    """The plane lattice's own point group (its holohedry) as ``(N, 2, 2)``
    Cartesian matrices: every rotation and reflection that maps the lattice
    spanned by ``g1``, ``g2`` onto itself.

    Found rather than classified. A map sending the basis to two lattice
    vectors is a symmetry exactly when it preserves all three inner products
    ``g1.g1``, ``g2.g2``, ``g1.g2`` (the Gram matrix) -- then it is
    orthogonal, and maps the lattice onto a sublattice of equal cell area,
    i.e. onto the whole lattice. So: enumerate the lattice points within
    ``max_index`` shells, keep the pairs whose Gram matrix matches, and read
    off the matrix each pair defines. This gets the oblique (2), rectangular
    and centred-rectangular (2mm), square (4mm) and hexagonal (6mm) cases
    right without a table of special cases, including a centred-rectangular
    cell handed over as ``a = b`` with an odd angle, which a classifier
    keyed on "is gamma 90 or 120" would call oblique.

    Note this is the *lattice's* symmetry. A real layer's point group is
    this one or a subgroup of it, depending on what sits at each lattice
    site; the irreducible wedge below is therefore the smallest one the
    lattice allows, and a layer with a less symmetric basis has a larger
    one. The Brillouin-zone dialog says so where it offers it, since a 2-D
    layer is given as three numbers with no basis to work from.
    """
    g1, g2 = np.asarray(g1, float), np.asarray(g2, float)
    basis = np.column_stack([g1, g2])
    if abs(float(np.linalg.det(basis))) < 1e-12:
        raise ValueError("g1 and g2 are parallel (degenerate 2-D lattice)")
    gram = basis.T @ basis
    scale = float(np.max(np.abs(gram)))
    inverse = np.linalg.inv(basis)

    idx = np.arange(-max_index, max_index + 1)
    points = np.array([m * g1 + n * g2 for m in idx for n in idx])
    lengths = np.einsum("ij,ij->i", points, points)

    first = points[np.abs(lengths - gram[0, 0]) < 1e-9 * scale]
    second = points[np.abs(lengths - gram[1, 1]) < 1e-9 * scale]
    operations = []
    for v1 in first:
        for v2 in second:
            if abs(float(v1 @ v2) - gram[0, 1]) > 1e-9 * scale:
                continue
            matrix = np.column_stack([v1, v2]) @ inverse
            if np.allclose(matrix @ matrix.T, np.eye(2), atol=1e-9):
                operations.append(matrix)
    if not operations:
        raise RuntimeError("no lattice symmetry found, not even the identity")
    return np.array(operations)


def _clip_half_plane(polygon: np.ndarray, normal: np.ndarray, tol: float):
    """The part of a convex ``polygon`` (open, no repeated last point) with
    ``normal . k <= 0``, by Sutherland-Hodgman clipping."""
    kept = []
    n = len(polygon)
    for i in range(n):
        current, following = polygon[i], polygon[(i + 1) % n]
        here, there = float(normal @ current), float(normal @ following)
        if here <= tol:
            kept.append(current)
        if (here > tol) != (there > tol):
            kept.append(current + here / (here - there) * (following - current))
    return np.array(kept) if kept else np.zeros((0, 2))


def _simplify_polygon(polygon: np.ndarray, tol: float) -> np.ndarray:
    """Drop repeated and collinear vertices from a convex polygon (open).

    Clipping a polygon with a line that runs exactly through one of its
    vertices reports that vertex twice -- once as a kept point, once as the
    crossing -- and clipping a wedge out of a many-sided cell leaves several
    points strung along what is really one straight edge. Neither changes
    the shape or its area, but both make a Gamma-M-K wedge come back as a
    thirteen-sided polygon, so they are removed rather than drawn.
    """
    kept = []
    for point in polygon:
        if not kept or np.linalg.norm(point - kept[-1]) > tol:
            kept.append(point)
    while len(kept) > 1 and np.linalg.norm(kept[0] - kept[-1]) <= tol:
        kept.pop()
    if len(kept) < 3:
        return np.array(kept)
    straight = []
    n = len(kept)
    for i in range(n):
        before, here, after = kept[i - 1], kept[i], kept[(i + 1) % n]
        edge_in, edge_out = here - before, after - here
        lengths = np.linalg.norm(edge_in) * np.linalg.norm(edge_out)
        if lengths <= 0:
            continue
        if abs(_cross2(edge_in, edge_out)) > tol * np.sqrt(lengths):
            straight.append(here)
    return np.array(straight) if len(straight) >= 3 else np.array(kept)


def irreducible_cell_2d(polygon: np.ndarray, operations) -> np.ndarray:
    """The irreducible 2-D Brillouin zone: the wedge of ``polygon`` that the
    ``operations`` (from :func:`lattice_point_group_2d`) repeat into the
    whole of it. Returns a closed ``(n+1, 2)`` array, like
    :func:`wigner_seitz_cell_2d`.

    Same construction as :func:`tools.bz3d.irreducible_wedge` one dimension
    down -- the Dirichlet cell of a generic point's orbit, which here is a
    wedge cut by half-planes ``k . (g.p - p) <= 0`` through the origin -- so
    the same check applies: the wedge's area is the zone's area divided by
    the number of operations, which ``test_bz2d.py`` asserts.

    Time reversal needs no special treatment in 2-D the way it does in 3-D:
    every plane lattice's holohedry already contains the two-fold rotation,
    which *is* inversion in the plane.
    """
    ops = np.asarray(operations, dtype=float)
    if ops.ndim != 3 or ops.shape[1:] != (2, 2):
        raise ValueError("operations must be an (N, 2, 2) array")
    polygon = np.asarray(polygon, dtype=float)
    open_polygon = polygon[:-1] if np.allclose(polygon[0], polygon[-1]) else polygon
    if len(ops) < 2:
        return np.vstack([open_polygon, open_polygon[:1]])

    # distance from the origin to the nearest edge: |edge x (origin - start)|
    # over |edge|, the usual point-to-line distance
    inradius = float("inf")
    for i, start in enumerate(open_polygon):
        edge = open_polygon[(i + 1) % len(open_polygon)] - start
        length = float(np.linalg.norm(edge))
        if length < 1e-12:
            continue
        inradius = min(inradius, abs(_cross2(edge, -start)) / length)
    if not np.isfinite(inradius) or inradius <= 0:
        raise ValueError("the zone polygon is degenerate")

    seed = None
    for direction in ((0.31286, 0.21763), (0.27182, 0.16180), (0.41421, 0.17320)):
        candidate = np.array(direction, dtype=float)
        candidate = 0.4 * inradius * candidate / np.linalg.norm(candidate)
        orbit = ops @ candidate
        spread = np.linalg.norm(orbit[:, None, :] - orbit[None, :, :], axis=-1)
        np.fill_diagonal(spread, np.inf)
        if spread.min() > 1e-6 * inradius:
            seed = candidate
            break
    if seed is None:
        raise RuntimeError("no point of the cell is moved to distinct places "
                           "by every operation -- check the operation set")

    wedge = open_polygon
    tol = 1e-12 * max(inradius, 1e-12)
    for image in ops @ seed:
        normal = image - seed
        if float(np.linalg.norm(normal)) < 1e-9 * inradius:
            continue
        wedge = _clip_half_plane(wedge, normal / np.linalg.norm(normal), tol)
        if len(wedge) < 3:
            raise RuntimeError("the irreducible wedge collapsed -- the "
                               "operations may not form a group")
    wedge = _simplify_polygon(wedge, 1e-9 * inradius)
    if len(wedge) < 3:
        raise RuntimeError("the irreducible wedge collapsed to a line")
    centroid = wedge.mean(axis=0)
    angles = np.arctan2(wedge[:, 1] - centroid[1], wedge[:, 0] - centroid[0])
    wedge = wedge[np.argsort(angles)]
    return np.vstack([wedge, wedge[:1]])


def tile_2d(polygon: np.ndarray, g1: np.ndarray, g2: np.ndarray,
           x_range, y_range, margin: int = 2):
    """Copies of ``polygon`` (as returned by :func:`wigner_seitz_cell_2d`,
    or any closed polygon centred on the origin) translated to every point
    of the lattice spanned by ``g1``, ``g2`` that falls near
    ``x_range``/``y_range`` -- "the same cell, repeated" across the field
    of view. Returns a list of closed ``(n+1, 2)`` arrays.
    """
    g1, g2 = np.asarray(g1, float), np.asarray(g2, float)
    span = max(x_range[1] - x_range[0], y_range[1] - y_range[0], 1e-9)
    n_i = min(int(np.ceil(span / max(np.linalg.norm(g1), 1e-9))) + margin, 40)
    n_j = min(int(np.ceil(span / max(np.linalg.norm(g2), 1e-9))) + margin, 40)

    poly_extent = float(np.max(np.linalg.norm(polygon, axis=1)))
    tiles = []
    for i in range(-n_i, n_i + 1):
        for j in range(-n_j, n_j + 1):
            shift = i * g1 + j * g2
            shifted = polygon + shift
            if (shifted[:, 0].max() < x_range[0] - poly_extent or
                    shifted[:, 0].min() > x_range[1] + poly_extent or
                    shifted[:, 1].max() < y_range[0] - poly_extent or
                    shifted[:, 1].min() > y_range[1] + poly_extent):
                continue
            tiles.append(shifted)
    return tiles
