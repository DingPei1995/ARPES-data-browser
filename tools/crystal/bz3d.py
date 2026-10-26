"""
tools/crystal/bz3d.py
=====================
The general 3-D Brillouin zone (the reciprocal lattice's Wigner-Seitz cell),
an arbitrary-orientation planar cut through it, and tiling copies of it
across a field of view -- the geometry engine behind the Brillouin-zone
tool's "3D crystal" mode.

Relation to the lab's MATLAB tools
-----------------------------------
``gen_brillouin.m`` builds the same cell by the same idea -- start from a
box, clip it face by face with the perpendicular-bisector plane to each
nearby reciprocal-lattice point -- but only ever considers the 26 nearest
neighbours (``-1:1`` in each direction). That is enough for most lattices,
but the true Wigner-Seitz cell of a sufficiently oblique cell can need a
face from the *second* shell of neighbours (a well-known pitfall of this
construction, not specific to this implementation of it); ``gen_brillouin.m``
has no check for that. :func:`wigner_seitz_cell` here grows the neighbour
shell (1, 2, 3, ...) and stops only once the resulting volume has converged,
so it does not silently under-count faces for an oblique lattice.

The clipping itself is also done differently: rather than intersecting each
candidate plane with the *previous* polyhedron in sequence (``gen_brillouin.m``'s
approach, and ``cut_brillouin.m``'s plane-vs-polygon-list routine for cutting
it afterwards), everything here goes through
:class:`scipy.spatial.HalfspaceIntersection`, which solves for all of the
half-spaces' common intersection at once via its interior point (the origin
-- trivially interior, since every half-space here is "closer to Gamma than
to G"). :func:`plane_cut`, likewise, does not chain intersection segments
face-by-face into a loop the way ``cut_brillouin.m`` does (matching endpoints
within a tolerance, which is exactly the kind of thing that gets fragile at
a cut that grazes an edge or a vertex): because the Wigner-Seitz cell is
convex, any planar cut through it is a single convex polygon, so its points
can just be angle-sorted around their own centroid -- simpler, and does not
depend on tuning an epsilon for "these two points are the same vertex".
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.spatial import ConvexHull, HalfspaceIntersection

__all__ = ["BrillouinZone", "wigner_seitz_cell", "plane_basis", "plane_cut",
           "tile_and_cut", "face_planes", "irreducible_wedge", "cut_points_3d",
           "cuts_in_frame", "kz_frame", "kz_zone_cuts"]


@dataclass
class BrillouinZone:
    """A convex polyhedron: :attr:`faces` for drawing (each an ordered loop
    of vertices, already closed -- the first point is repeated at the end),
    :attr:`vertices` for bounding-box/tiling-range estimates."""
    faces: list          # list of (n_i, 3) arrays, each face's vertices in order
    vertices: np.ndarray  # (N, 3), the polyhedron's corner points
    volume: float
    shells_used: int


def _face_loop(points: np.ndarray, normal: np.ndarray) -> np.ndarray:
    """Order a face's (unordered) vertices into a closed loop by angle
    around their centroid, using ``normal`` to fix a consistent winding.
    """
    centroid = points.mean(axis=0)
    # any two orthogonal in-plane axes will do; the loop order they give is
    # only used for drawing, not for any signed-area computation
    ref = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(ref, normal)) > 0.9:
        ref = np.array([0.0, 1.0, 0.0])
    u = ref - np.dot(ref, normal) * normal
    u /= np.linalg.norm(u)
    v = np.cross(normal, u)
    rel = points - centroid
    angles = np.arctan2(rel @ v, rel @ u)
    order = np.argsort(angles)
    loop = points[order]
    return np.vstack([loop, loop[:1]])


def _faces_from_halfspaces(verts: np.ndarray, normals: np.ndarray,
                           offsets: np.ndarray, scale: float):
    """The polyhedron's faces, as ordered vertex loops: for each bounding
    plane ``normal . x = offset``, the vertices lying on it.

    Planes that only touch the polytope at an edge or a vertex (fewer than
    three vertices on them) are redundant and dropped, which is what makes
    this usable with a deliberately over-complete set of half-spaces.
    """
    faces = []
    tol = 1e-7 * scale
    for normal, offset in zip(normals, offsets):
        on_plane = np.abs(verts @ normal - offset) < tol
        if on_plane.sum() < 3:
            continue
        faces.append(_face_loop(verts[on_plane], normal / np.linalg.norm(normal)))
    return faces


def face_planes(faces, tol: float = 1e-12):
    """``(normals, offsets)`` of the planes a polyhedron's ``faces`` lie in,
    as unit normals pointing away from the origin (``normal . x = offset``
    with ``offset > 0``).

    The inverse of :func:`_faces_from_halfspaces`: it recovers the bounding
    half-spaces from the vertex loops, so a zone built once can be fed back
    into another half-space intersection (which is how
    :func:`irreducible_wedge` carves it up) without keeping a separate copy
    of the planes around. Assumes the origin is strictly inside, which is
    true of every Wigner-Seitz cell by construction.
    """
    normals, offsets = [], []
    for loop in faces:
        verts = loop[:-1]
        centroid = verts.mean(axis=0)
        rel = verts - centroid
        # the largest cross product among the spokes to the centroid: the
        # most numerically robust normal this face offers, and immune to
        # three nearly-collinear vertices that a fixed choice would trip on
        best, normal = 0.0, None
        for i in range(len(rel)):
            for j in range(i + 1, len(rel)):
                candidate = np.cross(rel[i], rel[j])
                length = float(np.linalg.norm(candidate))
                if length > best:
                    best, normal = length, candidate
        if normal is None or best < tol:
            continue
        normal = normal / np.linalg.norm(normal)
        offset = float(np.dot(normal, centroid))
        if offset < 0:
            normal, offset = -normal, -offset
        normals.append(normal)
        offsets.append(offset)
    return np.array(normals), np.array(offsets)


def wigner_seitz_cell(b: np.ndarray, max_shell: int = 5,
                      rel_tol: float = 1e-9) -> BrillouinZone:
    """The first Brillouin zone of the reciprocal lattice spanned by ``b``
    (rows ``b1, b2, b3``): the Wigner-Seitz cell about the origin.

    Grows the neighbour shell used for the bounding half-spaces (``-n..n``
    in each of the three lattice indices, excluding the origin) from 1
    until the polytope's volume changes by less than ``rel_tol`` between
    successive shells, or until ``max_shell`` is reached (a warning-worthy
    situation for any physically normal cell; raises past that point rather
    than silently returning an under-converged shape).
    """
    b1, b2, b3 = b
    target = abs(float(np.linalg.det(np.asarray(b, dtype=float))))
    previous_volume = None
    result = None
    for shell in range(1, max_shell + 1):
        idx = np.arange(-shell, shell + 1)
        I, J, K = np.meshgrid(idx, idx, idx, indexing="ij")
        I, J, K = I.ravel(), J.ravel(), K.ravel()
        keep = ~((I == 0) & (J == 0) & (K == 0))
        I, J, K = I[keep], J[keep], K[keep]
        G = I[:, None] * b1 + J[:, None] * b2 + K[:, None] * b3
        g2 = np.einsum("ij,ij->i", G, G)
        # half-space "closer to Gamma than to G": G.x - |G|^2/2 <= 0
        halfspaces = np.hstack([G, (-0.5 * g2)[:, None]])

        hs = HalfspaceIntersection(halfspaces, np.zeros(3))
        verts = hs.intersections
        volume = ConvexHull(verts).volume

        faces = _faces_from_halfspaces(verts, G, 0.5 * g2, np.sqrt(g2.max()))

        result = BrillouinZone(faces=faces, vertices=verts, volume=volume,
                               shells_used=shell)
        # The intersection always *contains* the true cell, and the true
        # cell's volume is the reciprocal cell's, so equal volume means
        # done -- an exact test, where "the same as last shell" could in
        # principle stall on two equally incomplete shells.
        if abs(volume - target) <= 1e-7 * target:
            return result
        if previous_volume is not None and \
                abs(volume - previous_volume) <= rel_tol * previous_volume:
            return result
        previous_volume = volume
    raise RuntimeError(
        f"Wigner-Seitz cell did not converge within {max_shell} neighbour "
        "shells -- check the lattice vectors (an extremely oblique cell can "
        "need a larger max_shell).")


#: directions tried, in order, for the generic seed point
#: :func:`irreducible_wedge` builds its fundamental domain around. The first
#: is deliberately "irrational-looking" so that it misses every mirror plane
#: and rotation axis of every crystallographic point group; the rest are
#: fallbacks in case a low-symmetry lattice happens to align with it.
_SEED_DIRECTIONS = (
    (0.31286, 0.21763, 0.14159),
    (0.27182, 0.31831, 0.16180),
    (0.41421, 0.17320, 0.22360),
)


def irreducible_wedge(zone: BrillouinZone, operations) -> BrillouinZone:
    """The irreducible Brillouin zone: the part of ``zone`` from which the
    whole of it can be rebuilt by the crystal's own symmetry.

    ``operations`` are Cartesian 3x3 rotation matrices -- the Laue class, as
    :func:`tools.crystal.lattice.point_group_operations` returns it. Returns the same
    :class:`BrillouinZone` shape as the full zone, so it cuts, tiles and
    draws through exactly the same code.

    How it is built: pick a point ``p`` inside the zone that no operation
    leaves fixed, and keep the points of the zone that are at least as close
    to ``p`` as to every image ``g.p`` -- the Dirichlet (Voronoi) cell of the
    orbit of ``p``, which is a standard fundamental domain for a finite group
    of isometries fixing the origin. Because every operation preserves length
    (``|g.p| = |p|``), each of those bisector planes passes through the
    origin and reduces to the half-space ``k . (g.p - p) <= 0``, so the
    wedge is a cone from Gamma clipped by the zone -- convex, and obtainable
    from the very same :class:`scipy.spatial.HalfspaceIntersection` the zone
    itself came from.

    This is preferred here over the usual textbook route of naming the
    high-symmetry points of each lattice type and hard-coding the wedge
    between them (which is a table per Bravais lattice, and silently wrong
    for any setting the table did not anticipate): nothing below knows which
    crystal system it is looking at, and the result is checkable --
    ``volume == zone.volume / len(operations)`` exactly, which is what
    ``test_bz3d.py`` asserts for the cubic, hexagonal and triclinic cases.
    """
    ops = np.asarray(operations, dtype=float)
    if ops.ndim != 3 or ops.shape[1:] != (3, 3):
        raise ValueError("operations must be an (N, 3, 3) array of rotations")
    if len(ops) < 2:
        return zone                       # nothing to reduce: C1 (plus nothing)

    normals, offsets = face_planes(zone.faces)
    if len(normals) == 0:
        raise ValueError("the zone has no usable faces")
    inradius = float(np.min(offsets))

    seed = None
    for direction in _SEED_DIRECTIONS:
        candidate = np.array(direction, dtype=float)
        candidate = 0.4 * inradius * candidate / np.linalg.norm(candidate)
        orbit = ops @ candidate
        spread = np.linalg.norm(orbit[:, None, :] - orbit[None, :, :], axis=-1)
        np.fill_diagonal(spread, np.inf)
        if spread.min() > 1e-6 * inradius:
            seed = candidate
            break
    if seed is None:
        raise RuntimeError(
            "could not find a point of the zone that the symmetry operations "
            "move to distinct places -- the operations may not form a group")

    # k . (g.p - p) <= 0 for every g, de-duplicated by direction: several
    # operations routinely give the same bisector, and handing qhull the same
    # plane many times is asking for trouble for no gain.
    wedge_normals = []
    for image in ops @ seed:
        normal = image - seed
        length = float(np.linalg.norm(normal))
        if length < 1e-9 * inradius:
            continue                      # g fixes the seed: no constraint
        normal = normal / length
        if not any(np.allclose(normal, seen, atol=1e-9) for seen in wedge_normals):
            wedge_normals.append(normal)
    if not wedge_normals:
        return zone

    all_normals = np.vstack([normals, np.array(wedge_normals)])
    all_offsets = np.concatenate([offsets, np.zeros(len(wedge_normals))])
    halfspaces = np.hstack([all_normals, -all_offsets[:, None]])

    hs = HalfspaceIntersection(halfspaces, seed)
    verts = hs.intersections
    hull = ConvexHull(verts)
    faces = _faces_from_halfspaces(verts, all_normals, all_offsets,
                                   float(np.max(np.abs(verts))))
    return BrillouinZone(faces=faces, vertices=verts, volume=hull.volume,
                         shells_used=zone.shells_used)


def plane_basis(normal: np.ndarray):
    """Two orthonormal in-plane axes (u, v) for the plane with this unit
    normal, matching the convention ``cut_brillouin.m``/``update_cut`` use
    (u, v pick themselves up from z-hat unless the normal already is
    z-hat), so a cut along (0, 0, 1) reproduces the familiar kx/ky layout.
    """
    normal = normal / np.linalg.norm(normal)
    z = np.array([0.0, 0.0, 1.0])
    if np.allclose(normal, z) or np.allclose(normal, -z):
        u = np.array([1.0, 0.0, 0.0])
    else:
        u = np.cross(z, normal)
        u /= np.linalg.norm(u)
    v = np.cross(normal, u)
    return u, v


def _project_and_sort(points: np.ndarray, normal: np.ndarray, point: np.ndarray):
    """Project coplanar 3-D ``points`` into the plane's own (u, v) basis and
    order them into a closed convex loop by angle about their centroid."""
    u_axis, v_axis = plane_basis(normal)
    u = (points - point) @ u_axis
    v = (points - point) @ v_axis
    centroid_u, centroid_v = u.mean(), v.mean()
    order = np.argsort(np.arctan2(v - centroid_v, u - centroid_u))
    u, v = u[order], v[order]
    return np.append(u, u[0]), np.append(v, v[0])


def plane_cut(faces, normal: np.ndarray, point: np.ndarray, tol: float = 1e-9):
    """The polygon where the plane through ``point`` with unit ``normal``
    cuts the convex polyhedron given by ``faces`` (as returned in
    :attr:`BrillouinZone.faces`, already translated to wherever this
    particular copy of the zone sits).

    Returns ``(u, v)`` coordinates (each a 1-D array, the polygon already
    closed) in the plane's own basis (:func:`plane_basis`), or ``(None,
    None)`` if this copy of the polyhedron does not reach the plane at all.

    A plane that happens to *contain* one of the faces is handled first and
    separately: there the cut is that whole face, but no edge of it crosses
    from one side to the other, so the edge-walk below would find nothing
    and wrongly report an empty cut. Far from a corner case, this is the
    normal situation for an irreducible wedge, whose faces are the crystal's
    own mirror planes -- exactly the planes a user asks to cut along.
    """
    normal = normal / np.linalg.norm(normal)
    d = np.dot(normal, point)

    scale = max((float(np.max(np.abs(loop))) for loop in faces), default=1.0)
    flat = max(tol, 1e-9 * max(scale, 1.0))
    for loop in faces:
        verts = loop[:-1]
        if len(verts) >= 3 and np.all(np.abs(verts @ normal - d) < flat):
            return _project_and_sort(verts, normal, point)

    pts = []
    for loop in faces:
        verts = loop[:-1]           # drop the repeated closing point
        signed = verts @ normal - d
        inside = signed < -tol
        n = len(verts)
        for i in range(n):
            j = (i + 1) % n
            if inside[i] == inside[j]:
                continue
            s1, s2 = signed[i], signed[j]
            if abs(s2 - s1) < tol:
                continue
            t = -s1 / (s2 - s1)
            pts.append(verts[i] + t * (verts[j] - verts[i]))
    if len(pts) < 3:
        return None, None

    pts = np.asarray(pts)
    # de-duplicate points where the plane passes through a shared edge/vertex
    # of several faces (each reports the same crossing independently)
    uniq = []
    for p in pts:
        if not any(np.linalg.norm(p - q) < 1e-7 for q in uniq):
            uniq.append(p)
    if len(uniq) < 3:
        return None, None
    return _project_and_sort(np.asarray(uniq), normal, point)


def cut_points_3d(u, v, normal: np.ndarray, point: np.ndarray) -> np.ndarray:
    """Put a cut polygon back into 3-D: the inverse of the projection at the
    end of :func:`plane_cut`, for drawing the cut where it actually sits
    inside the zone (which is what the 3-D preview of the cut plane shows).
    """
    normal = normal / np.linalg.norm(normal)
    u_axis, v_axis = plane_basis(normal)
    u = np.asarray(u, dtype=float)[:, None]
    v = np.asarray(v, dtype=float)[:, None]
    return np.asarray(point, dtype=float) + u * u_axis + v * v_axis


def tile_and_cut(zone: BrillouinZone, b: np.ndarray, normal: np.ndarray,
                 point: np.ndarray, u_range, v_range, margin: int = 2):
    """The cut polygons of every translated copy of ``zone`` whose cut
    reaches within ``u_range``/``v_range`` (each a ``(lo, hi)`` pair, in the
    plane's own (u, v) coordinates from :func:`plane_basis`) of the plane
    -- i.e. "the same Brillouin zone, repeated" rather than higher-order
    zones, since every copy is the same shape.

    ``margin`` is no longer used and is kept so existing callers still
    work; see below.

    Returns a list of ``(u, v)`` polygons (each closed, as from
    :func:`plane_cut`).

    The copies are found by :func:`cuts_in_frame`, from the window and the
    plane's actual position. This used to search a box of ``span / |b_i|``
    lattice indices around Gamma, which missed copies -- or every copy --
    once the plane was slid further along its normal than the window is
    wide (a rhombohedral (001) cut at 6 Å⁻¹ came back empty and was reported
    as "the plane does not meet the zone"), and under-counted for an
    oblique cell, where a short distance can take large lattice indices.
    """
    u_axis, v_axis = plane_basis(normal)
    span = max(u_range[1] - u_range[0], v_range[1] - v_range[0], 1e-9)
    # a little beyond the window, so a copy whose corner just pokes in at
    # the edge of the field of view is not dropped
    pad = 0.05 * span
    return cuts_in_frame(zone, b, point, u_axis, v_axis,
                         (u_range[0] - pad, u_range[1] + pad),
                         (v_range[0] - pad, v_range[1] + pad))


# -- cuts containing the surface normal: the zone a k_z map sees -------------
def cuts_in_frame(zone: BrillouinZone, b: np.ndarray, origin, e1, e2,
                  x_range, y_range, limit: int = 20000):
    """Cut every translated copy of ``zone`` by the plane through ``origin``
    spanned by the orthonormal ``e1`` and ``e2``, and return the polygons
    that reach the rectangle ``x_range`` x ``y_range`` -- in that frame's
    own coordinates, ``x = (p - origin) . e1`` and ``y = (p - origin) . e2``.

    :func:`tile_and_cut` answers the same question in the plane's own
    :func:`plane_basis`, which is fine for a constant-energy map (the
    azimuth turns it afterwards) but useless when the two drawn axes are
    fixed by the measurement, as a k_z map's are: there one axis *is* the
    surface normal and the other the slit direction, and a cut has to come
    out in exactly those.

    It also does not search a box of indices around the origin: which
    copies can touch the window is worked out from the window itself, so a
    k_z map that starts at 3 A^-1 costs no more than one starting at 0.
    """
    b = np.asarray(b, dtype=float)
    e1 = np.asarray(e1, dtype=float)
    e2 = np.asarray(e2, dtype=float)
    origin = np.asarray(origin, dtype=float)
    normal = np.cross(e1, e2)
    normal /= np.linalg.norm(normal)
    reach = float(np.max(np.linalg.norm(zone.vertices, axis=1)))

    # every point within ``reach`` of the window, as a box of corners, in
    # lattice coordinates: any copy whose centre is outside its index range
    # cannot reach the window at all
    x0, x1 = float(min(x_range)) - reach, float(max(x_range)) + reach
    y0, y1 = float(min(y_range)) - reach, float(max(y_range)) + reach
    corners = np.array([origin + x * e1 + y * e2 + s * normal
                        for x in (x0, x1) for y in (y0, y1)
                        for s in (-reach, reach)])
    indices = corners @ np.linalg.inv(b)
    low = np.floor(indices.min(axis=0)).astype(int)
    high = np.ceil(indices.max(axis=0)).astype(int)
    if np.prod(high - low + 1) > 50 * limit:
        raise ValueError("the window is too large for this zone: it would "
                         "take too many copies of the zone to fill")
    grid = np.stack(np.meshgrid(*(np.arange(lo, hi + 1)
                                  for lo, hi in zip(low, high)),
                                indexing="ij"), axis=-1).reshape(-1, 3)
    centres = grid @ b
    relative = centres - origin
    along_x, along_y = relative @ e1, relative @ e2
    keep = ((np.abs(relative @ normal) <= reach + 1e-9)
            & (along_x >= x0) & (along_x <= x1)
            & (along_y >= y0) & (along_y <= y1))
    centres = centres[keep]
    if len(centres) > limit:
        raise ValueError(f"{len(centres)} copies of the zone reach this "
                         f"window; that is too many to draw")

    polygons = []
    for centre in centres:
        u, v = plane_cut([loop + centre for loop in zone.faces], normal, origin)
        if u is None:
            continue
        points = cut_points_3d(u, v, normal, origin) - origin
        xs, ys = points @ e1, points @ e2
        if (xs.max() < min(x_range) or xs.min() > max(x_range)
                or ys.max() < min(y_range) or ys.min() > max(y_range)):
            continue
        polygons.append((xs, ys))
    return polygons


def kz_frame(params, surface_hkl, inplane_hkl=None, rotation_deg: float = 0.0):
    """The two axes of a k_z map, as unit vectors in Cartesian k-space.

    Returns ``(along_kz, along_kpar, across)``: the surface normal, the
    in-plane direction the slit measures along, and the third axis that
    completes them (the one a deflector or manipulator tilt would move
    along).

    ``surface_hkl`` is read against the **conventional** reciprocal cell,
    like every Miller index the user types (see
    :meth:`ui.viewers.windows.BrillouinZoneDialog._zone_3d`). ``inplane_hkl`` names
    a reciprocal-lattice direction whose component perpendicular to the
    normal becomes k_par; ``None`` (or ``(0, 0, 0)``) picks the first
    conventional reciprocal axis that is not along the normal -- a* for a
    (001) surface, which is the usual Gamma-M of a hexagonal crystal.
    ``rotation_deg`` then turns k_par about the normal, for a slit that sits
    between two named directions.
    """
    from tools.crystal.lattice import conventional_vectors, reciprocal_vectors

    b_conventional = reciprocal_vectors(conventional_vectors(params))
    normal = np.asarray(surface_hkl, dtype=float) @ b_conventional
    length = float(np.linalg.norm(normal))
    if length < 1e-9:
        raise ValueError("the surface normal (h, k, l) cannot be (0, 0, 0)")
    normal = normal / length

    def perpendicular(vector):
        vector = np.asarray(vector, dtype=float)
        scale = float(np.linalg.norm(vector))
        if scale < 1e-12:
            return None
        flat = vector - np.dot(vector, normal) * normal
        if np.linalg.norm(flat) < 1e-6 * scale:
            return None
        return flat / np.linalg.norm(flat)

    if inplane_hkl is not None and any(int(v) for v in inplane_hkl):
        along = perpendicular(np.asarray(inplane_hkl, dtype=float) @ b_conventional)
        if along is None:
            raise ValueError("the k_par direction is parallel to the surface "
                             "normal, so it has no in-plane part")
    else:
        along = None
        for axis in b_conventional:
            along = perpendicular(axis)
            if along is not None:
                break
    across = np.cross(normal, along)
    radians = np.radians(float(rotation_deg))
    along, across = (np.cos(radians) * along + np.sin(radians) * across,
                     -np.sin(radians) * along + np.cos(radians) * across)
    return normal, along, across


def kz_zone_cuts(zone: BrillouinZone, b: np.ndarray, frame, kz_range,
                 kpar_range, offset: float = 0.0, tile: bool = True):
    """The zone as a k_z map shows it: polygons in ``(k_z, k_par)``.

    ``frame`` is :func:`kz_frame`'s three axes. ``offset`` slides the
    measured plane along the third, out-of-slit axis, from Gamma -- a map
    taken off normal emission in the other direction. With ``tile`` off
    only the one copy whose centre is nearest the middle of the window is
    returned, rather than the one at Gamma, which a k_z map several zones
    out from the surface would not even contain.
    """
    along_kz, along_kpar, across = (np.asarray(v, dtype=float) for v in frame)
    origin = float(offset) * across
    polygons = cuts_in_frame(zone, b, origin, along_kz, along_kpar,
                             kz_range, kpar_range)
    if tile or not polygons:
        return polygons
    middle = np.array([0.5 * (kz_range[0] + kz_range[1]),
                       0.5 * (kpar_range[0] + kpar_range[1])])
    distances = [np.hypot(xs[:-1].mean() - middle[0], ys[:-1].mean() - middle[1])
                 for xs, ys in polygons]
    return [polygons[int(np.argmin(distances))]]
