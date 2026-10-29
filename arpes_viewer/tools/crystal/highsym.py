"""
tools/crystal/highsym.py
========================
The named high-symmetry points of the first Brillouin zone (Gamma, M, K,
Z, L, ...), for the Brillouin-zone tool: labels in its 3-D view of the cut
plane, labels on the contour, and the list of offsets along a cut normal
that pass through one of them.

A table, not a search
---------------------
Which points are special, and what they are called, is fixed by the Bravais
lattice -- and so by the space-group number, which :mod:`tools.crystal.spacegroups`
already maps to one. So the points are a static table, :data:`POINT_TABLE`,
in fractional coordinates of the *primitive* reciprocal vectors
(:func:`tools.crystal.lattice.primitive_vectors` -> :func:`tools.crystal.lattice.reciprocal_vectors`),
following Setyawan & Curtarolo, Comput. Mater. Sci. 49, 299 (2010) -- the
labelling most band-structure codes (and papers) use.

For most lattices the table is plain numbers. For a few the positions move
with the cell's shape -- a rhombohedral zone's B, P, Q and X slide with its
angle, a body-centred tetragonal one's Z and Sigma with c/a -- so those rows
are small formulas of the cell, evaluated when asked. Which variant applies
(RHL1 or RHL2, BCT1 or BCT2) is the cell's shape too.

The primitive cells here are the same as Setyawan & Curtarolo's for every
lattice in the table (the centring transforms in :mod:`tools.crystal.lattice` are
the standard ones), with one difference that does not matter: the
rhombohedral cell is built from hexagonal axes, so its three vectors are
the paper's up to a relabelling and a rotation -- and since all three are
equivalent by symmetry, fractional coordinates read against either set give
the same points. ``test_highsym.py`` checks every tabulated point lands on
the zone's surface, which is what would catch a wrong table entry.

Lattices without an entry (face-, body- and base-centred orthorhombic,
monoclinic, triclinic, whose labels depend on axis orderings this program
does not impose) fall back to Gamma plus the centre of each zone face,
labelled by its coordinates in the *conventional* reciprocal axes, e.g.
"(½ 0 0)" -- unnamed, but still the points a cut most often wants to pass
through.

No Qt in here.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from arpes_viewer.tools.crystal.lattice import (LatticeParams, bravais_symbol,
                                   conventional_vectors, primitive_vectors,
                                   reciprocal_vectors, point_group_operations)

__all__ = ["HighSymmetryPoint", "POINT_TABLE", "lattice_variant",
           "high_symmetry_points", "equivalent_points", "points_on_plane",
           "plane_offsets", "points_in_frame"]


@dataclass(frozen=True)
class HighSymmetryPoint:
    label: str
    k: np.ndarray        # Cartesian, Å⁻¹ (2*pi included)

    def __hash__(self):
        return hash((self.label, tuple(np.round(self.k, 9))))

    def __eq__(self, other):
        return (isinstance(other, HighSymmetryPoint)
                and self.label == other.label
                and np.allclose(self.k, other.k, atol=1e-9))


def _rhl1(alpha):
    eta = (1 + 4 * np.cos(alpha)) / (2 + 4 * np.cos(alpha))
    nu = 0.75 - eta / 2
    return {"Γ": (0, 0, 0), "B": (eta, 0.5, 1 - eta),
            "B1": (0.5, 1 - eta, eta - 1), "F": (0.5, 0.5, 0),
            "L": (0.5, 0, 0), "L1": (0, 0, -0.5), "P": (eta, nu, nu),
            "P1": (1 - nu, 1 - nu, 1 - eta), "P2": (nu, nu, eta - 1),
            "Q": (1 - nu, nu, 0), "X": (nu, 0, -nu), "Z": (0.5, 0.5, 0.5)}


def _rhl2(alpha):
    eta = 1 / (2 * np.tan(alpha / 2) ** 2)
    nu = 0.75 - eta / 2
    return {"Γ": (0, 0, 0), "F": (0.5, -0.5, 0), "L": (0.5, 0, 0),
            "P": (1 - nu, -nu, 1 - nu), "P1": (nu, nu - 1, nu - 1),
            "Q": (eta, eta, eta), "Q1": (1 - eta, -eta, -eta),
            "Z": (0.5, -0.5, 0.5)}


def _bct1(a, c):
    eta = (1 + c ** 2 / a ** 2) / 4
    return {"Γ": (0, 0, 0), "M": (-0.5, 0.5, 0.5), "N": (0, 0.5, 0),
            "P": (0.25, 0.25, 0.25), "X": (0, 0, 0.5), "Z": (eta, eta, -eta),
            "Z1": (-eta, 1 - eta, eta)}


def _bct2(a, c):
    eta = (1 + a ** 2 / c ** 2) / 4
    zeta = a ** 2 / (2 * c ** 2)
    return {"Γ": (0, 0, 0), "N": (0, 0.5, 0), "P": (0.25, 0.25, 0.25),
            "Σ": (-eta, eta, eta), "Σ1": (eta, 1 - eta, -eta),
            "X": (0, 0, 0.5), "Y": (-zeta, zeta, 0.5),
            "Y1": (0.5, 0.5, -zeta), "Z": (0.5, 0.5, -0.5)}


#: Fractional coordinates (of the primitive reciprocal vectors) of each
#: lattice's named points, by Setyawan & Curtarolo variant. A dict is a
#: fixed table; a callable takes the cell and returns one.
POINT_TABLE = {
    "CUB": {"Γ": (0, 0, 0), "M": (0.5, 0.5, 0), "R": (0.5, 0.5, 0.5),
            "X": (0, 0.5, 0)},
    "FCC": {"Γ": (0, 0, 0), "K": (3 / 8, 3 / 8, 3 / 4), "L": (0.5, 0.5, 0.5),
            "U": (5 / 8, 1 / 4, 5 / 8), "W": (0.5, 0.25, 0.75),
            "X": (0.5, 0, 0.5)},
    "BCC": {"Γ": (0, 0, 0), "H": (0.5, -0.5, 0.5), "P": (0.25, 0.25, 0.25),
            "N": (0, 0, 0.5)},
    "TET": {"Γ": (0, 0, 0), "A": (0.5, 0.5, 0.5), "M": (0.5, 0.5, 0),
            "R": (0, 0.5, 0.5), "X": (0, 0.5, 0), "Z": (0, 0, 0.5)},
    "BCT1": _bct1,
    "BCT2": _bct2,
    # Named along the cell's own a*, b*, c* -- Setyawan & Curtarolo sort the
    # axes so that a < b < c first, which would rename X/Y/Z behind the
    # user's back for a cell typed in any other order.
    "ORC": {"Γ": (0, 0, 0), "R": (0.5, 0.5, 0.5), "S": (0.5, 0.5, 0),
            "T": (0, 0.5, 0.5), "U": (0.5, 0, 0.5), "X": (0.5, 0, 0),
            "Y": (0, 0.5, 0), "Z": (0, 0, 0.5)},
    "HEX": {"Γ": (0, 0, 0), "A": (0, 0, 0.5), "H": (1 / 3, 1 / 3, 0.5),
            "K": (1 / 3, 1 / 3, 0), "L": (0.5, 0, 0.5), "M": (0.5, 0, 0)},
    "RHL1": _rhl1,
    "RHL2": _rhl2,
}


def lattice_variant(params: LatticeParams) -> "str | None":
    """Which row of :data:`POINT_TABLE` this cell uses, or None if there is
    none (the generic fallback then applies)."""
    symbol = bravais_symbol(params)
    if symbol == "cP":
        return "CUB"
    if symbol == "cF":
        return "FCC"
    if symbol == "cI":
        return "BCC"
    if symbol == "tP":
        return "TET"
    if symbol == "tI":
        return "BCT1" if params.c < params.a else "BCT2"
    if symbol == "oP":
        return "ORC"
    if symbol == "hP":
        return "HEX"
    if symbol == "hR":
        return "RHL1" if _rhombohedral_angle(params) < np.pi / 2 else "RHL2"
    return None


def _rhombohedral_angle(params: LatticeParams) -> float:
    a1, a2, _a3 = primitive_vectors(params)
    cosine = np.dot(a1, a2) / (np.linalg.norm(a1) * np.linalg.norm(a2))
    return float(np.arccos(np.clip(cosine, -1.0, 1.0)))


def _table_for(variant: str, params: LatticeParams) -> dict:
    entry = POINT_TABLE[variant]
    if not callable(entry):
        return entry
    if variant.startswith("RHL"):
        return entry(_rhombohedral_angle(params))
    return entry(params.a, params.c)


def high_symmetry_points(params: LatticeParams):
    """The named points, one representative each, as
    :class:`HighSymmetryPoint` -- Gamma first.

    The representatives are the table's; every other copy of a point (the
    six M of a hexagonal zone, say) comes from :func:`equivalent_points`.
    """
    b = reciprocal_vectors(primitive_vectors(params))
    variant = lattice_variant(params)
    if variant is not None:
        table = _table_for(variant, params)
        return [HighSymmetryPoint(label, np.asarray(frac, dtype=float) @ b)
                for label, frac in table.items()]
    return _generic_points(params, b)


def _fraction(value: float) -> str:
    for text, number in (("0", 0.0), ("½", 0.5), ("-½", -0.5), ("1", 1.0),
                         ("-1", -1.0), ("¼", 0.25), ("-¼", -0.25),
                         ("¾", 0.75), ("-¾", -0.75)):
        if abs(value - number) < 1e-6:
            return text
    return f"{value:.3g}"


def _generic_points(params: LatticeParams, b: np.ndarray):
    """Gamma, and the centre G/2 of each face of the zone (a
    time-reversal-invariant point, whatever the lattice), labelled by its
    coordinates in the conventional reciprocal axes. One per symmetry
    orbit."""
    from arpes_viewer.tools.crystal.bz3d import wigner_seitz_cell, face_planes
    zone = wigner_seitz_cell(b)
    normals, offsets = face_planes(zone.faces)
    b_conventional = reciprocal_vectors(conventional_vectors(params))
    to_conventional = np.linalg.inv(b_conventional)
    try:
        ops = point_group_operations(params)
    except Exception:                                    # noqa: BLE001
        ops = np.array([np.eye(3), -np.eye(3)])
    points = [HighSymmetryPoint("Γ", np.zeros(3))]
    seen = []
    for normal, offset in sorted(zip(normals, offsets), key=lambda p: p[1]):
        centre = normal * offset
        if any(np.allclose(op @ centre, other, atol=1e-7)
               for op in ops for other in seen):
            continue
        seen.append(centre)
        frac = centre @ to_conventional
        label = "(" + " ".join(_fraction(v) for v in frac) + ")"
        points.append(HighSymmetryPoint(label, centre))
    return points


# -- every copy of each point -------------------------------------------------
def _inside(zone_planes, point, tol):
    normals, offsets = zone_planes
    return bool(np.all(normals @ point <= offsets + tol))


def equivalent_points(params: LatticeParams, points=None, zone=None):
    """Every copy of each named point that lies in (or on the surface of)
    the first zone: its images under the Laue class, and those images
    moved by a reciprocal-lattice vector back onto the zone where that
    lands them on its surface (K and K' at a hexagonal zone's corners).

    Returns a list of :class:`HighSymmetryPoint`, each label as many times
    as the zone holds that point. A table entry whose copies are exactly
    another entry's (a second name for the same point) is left out.
    """
    from arpes_viewer.tools.crystal.bz3d import wigner_seitz_cell, face_planes
    b = reciprocal_vectors(primitive_vectors(params))
    if zone is None:
        zone = wigner_seitz_cell(b)
    if points is None:
        points = high_symmetry_points(params)
    planes = face_planes(zone.faces)
    scale = float(np.max(np.linalg.norm(zone.vertices, axis=1)))
    tol = 1e-7 * max(scale, 1.0)
    try:
        ops = point_group_operations(params)
    except Exception:                                    # noqa: BLE001
        ops = np.array([np.eye(3), -np.eye(3)])
    shifts = np.array([i * b[0] + j * b[1] + k * b[2]
                       for i in (-1, 0, 1) for j in (-1, 0, 1)
                       for k in (-1, 0, 1)])

    out, groups = [], []
    for point in points:
        found = []
        for image in ops @ point.k:
            for shift in shifts:
                candidate = image + shift
                if not _inside(planes, candidate, tol):
                    continue
                if not any(np.linalg.norm(candidate - other) < 1e-6
                           for other in found):
                    found.append(candidate)
        # Setyawan & Curtarolo name some copies of one point twice over for
        # the sake of band paths (L and L1 of a rhombohedral zone, K and U of
        # an FCC one): the same places, so they are labelled once, by the
        # first name the table gives them.
        if any(len(found) == len(group) and all(
                any(np.linalg.norm(k - other) < 1e-6 for other in group)
                for k in found) for group in groups):
            continue
        groups.append(found)
        out.extend(HighSymmetryPoint(point.label, k) for k in found)
    return out


def points_on_plane(points, normal, offset: float, tol: float = 1e-6):
    """The points of ``points`` lying in the plane ``normal . k = offset``
    (``normal`` a unit vector)."""
    normal = np.asarray(normal, dtype=float)
    return [p for p in points if abs(float(np.dot(p.k, normal)) - offset) < tol]


def plane_offsets(points, normal, tol: float = 1e-4):
    """The offsets along ``normal`` (unit) at which a plane passes through
    a high-symmetry point of the first zone: ``[(offset, labels), ...]``,
    in increasing order from 0, each offset once with every point it
    meets. Only offsets >= 0 -- the Laue class includes inversion, so the
    negative ones are the same planes mirrored.
    """
    normal = np.asarray(normal, dtype=float)
    normal = normal / np.linalg.norm(normal)
    found = []
    for point in points:
        offset = float(np.dot(point.k, normal))
        if offset < -tol:
            continue
        offset = max(offset, 0.0)
        for entry in found:
            if abs(entry[0] - offset) < tol:
                if point.label not in entry[1]:
                    entry[1].append(point.label)
                break
        else:
            found.append([offset, [point.label]])
    found.sort(key=lambda entry: entry[0])
    order = {}
    for point in points:
        order.setdefault(point.label, len(order))
    return [(offset, sorted(labels, key=lambda l: order[l]))
            for offset, labels in found]


def points_in_frame(points, b, origin, e1, e2, x_range, y_range,
                    tol: float = 1e-6, limit: int = 400):
    """The copies of ``points`` (a first-zone set from
    :func:`equivalent_points`), translated by every reciprocal-lattice
    vector, that lie in the plane through ``origin`` spanned by the
    orthonormal ``e1``, ``e2`` and inside ``x_range`` x ``y_range`` of it.

    Returns ``[(label, x, y), ...]`` in that frame's coordinates, the same
    ones :func:`tools.crystal.bz3d.cuts_in_frame` draws the zone in, each place
    once (a point on a face shared by two copies of the zone is one
    point). At most ``limit`` are returned, nearest the window's centre
    first, so that a very large window does not bury the picture in text.
    """
    b = np.asarray(b, dtype=float)
    e1, e2 = np.asarray(e1, dtype=float), np.asarray(e2, dtype=float)
    origin = np.asarray(origin, dtype=float)
    normal = np.cross(e1, e2)
    normal = normal / np.linalg.norm(normal)
    x0, x1 = float(min(x_range)), float(max(x_range))
    y0, y1 = float(min(y_range)), float(max(y_range))
    reach = max((float(np.linalg.norm(p.k)) for p in points), default=0.0) + 1e-9
    corners = np.array([origin + x * e1 + y * e2 + s * normal
                        for x in (x0 - reach, x1 + reach)
                        for y in (y0 - reach, y1 + reach)
                        for s in (-reach, reach)])
    indices = corners @ np.linalg.inv(b)
    low = np.floor(indices.min(axis=0)).astype(int)
    high = np.ceil(indices.max(axis=0)).astype(int)
    if np.prod(high - low + 1) > 2_000_000:
        return []
    grid = np.stack(np.meshgrid(*(np.arange(lo, hi + 1)
                                  for lo, hi in zip(low, high)),
                                indexing="ij"), axis=-1).reshape(-1, 3)
    centres = grid @ b
    near = np.abs((centres - origin) @ normal) <= reach
    centres = centres[near]

    labelled = {}
    for point in points:
        positions = centres + point.k - origin
        distance = positions @ normal
        on = np.abs(distance) < tol
        xs, ys = positions[on] @ e1, positions[on] @ e2
        inside = (xs >= x0 - tol) & (xs <= x1 + tol) & (ys >= y0 - tol) & (ys <= y1 + tol)
        for x, y in zip(xs[inside], ys[inside]):
            key = (round(float(x), 5), round(float(y), 5))
            labelled.setdefault(key, point.label)
    middle = (0.5 * (x0 + x1), 0.5 * (y0 + y1))
    found = sorted(((label, x, y) for (x, y), label in labelled.items()),
                   key=lambda item: np.hypot(item[1] - middle[0], item[2] - middle[1]))
    return found[:limit]
