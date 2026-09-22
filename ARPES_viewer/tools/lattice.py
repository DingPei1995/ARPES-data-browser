"""
tools/lattice.py
==============
Real-space and reciprocal-space lattice vectors from unit-cell parameters
and (optionally) a space-group number, for the Brillouin-zone tool.

This has no counterpart in the lab's MATLAB tools. The closest thing,
``brillouin_zone_plot.m``'s ``gen_realspace_basevector.m``, only reads
``beta`` and ``gamma`` -- it has no ``alpha`` input at all (the line that
would read it is commented out) -- and reuses ``cos(beta)`` where the
general triclinic formula needs ``cos(alpha)``. That is invisible for the
crystal systems with alpha = beta = 90 deg (cubic, tetragonal, orthorhombic,
hexagonal), where both formulas happen to give zero, which is presumably why
it went unnoticed; it gives the wrong basis for monoclinic, triclinic and
trigonal/rhombohedral cells. Its bcc/fcc centering, likewise, is only the
conventional-cubic-cell transform (a fixed matrix that assumes a = b = c,
alpha = beta = gamma = 90 deg) applied unconditionally -- correct only for
cubic I/F, silently wrong (not flagged) for a body- or face-centered cell of
any other crystal system, and with no base-centered (A/B/C) option at all.

Here the two are separated properly:

* :func:`conventional_vectors` builds the general triclinic conventional
  cell from the six parameters, correct for any crystal system.
* :func:`primitive_vectors` applies the centering transform for whichever
  of the 14 Bravais lattices applies (from a space-group number, via
  :mod:`tools.spacegroups`, or given directly), not just I/F-cubic.

Space-group number, not a full symmetry search
-----------------------------------------------
The Brillouin-zone tool asks the user for a space-group number directly
(1-230) rather than an atomic basis, so there is no structure for a package
like ``spglib`` to search for symmetry in -- the space group is already
known. What is still useful from that number is purely reference data: which
of the 7 crystal systems and 14 Bravais lattices it belongs to, which is
exactly what :mod:`tools.spacegroups` provides. That module is a static table
generated once from spglib's own database (see ``devtools/gen_spacegroups.py``); the
running program does not import spglib at all.

Conventions
-----------
``a, b, c`` in angstrom; ``alpha, beta, gamma`` in degrees, ``alpha`` the
angle between b and c, ``beta`` between a and c, ``gamma`` between a and b
(the standard crystallographic convention). Real-space vectors are rows of a
(3, 3) array; so are reciprocal vectors, in units of 2*pi/angstrom (i.e. the
physicist's convention used throughout this program, matching
``tools.kspace``'s Å⁻¹ axes -- not the crystallographer's 1/angstrom without
the 2*pi).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from tools.spacegroups import spacegroup_info, point_group_rotations

__all__ = ["LatticeParams", "bravais_symbol", "validate_lattice_parameters",
           "conventional_vectors", "primitive_vectors", "reciprocal_vectors",
           "point_group_operations", "free_parameters",
           "CENTERING_TRANSFORMS", "LATTICE_POINTS_PER_CONVENTIONAL_CELL"]


@dataclass
class LatticeParams:
    """The six cell parameters plus how the user identified the Bravais
    lattice: either a space-group number, or an explicit centering letter
    for someone who knows the lattice type but not (or does not want to
    type) the exact space group."""
    a: float
    b: float
    c: float
    alpha: float = 90.0
    beta: float = 90.0
    gamma: float = 90.0
    space_group: "int | None" = None
    centering: "str | None" = None   # P/A/B/C/I/F/R, used only if space_group is None

    def resolved_centering(self) -> str:
        if self.space_group is not None:
            return spacegroup_info(self.space_group).centering
        return (self.centering or "P").upper()

    def crystal_system(self) -> "str | None":
        """None when only a bare centering letter was given -- a centering
        letter alone does not say which crystal system it belongs to (e.g.
        'P' is every one of them), so there is nothing to validate the
        angles against in that case."""
        if self.space_group is not None:
            return spacegroup_info(self.space_group).crystal_system
        return None


def bravais_symbol(params: LatticeParams) -> str:
    """e.g. 'cF', 'hR', 'mC' -- for display, and to key the centering
    transform table."""
    if params.space_group is not None:
        return spacegroup_info(params.space_group).bravais_symbol
    letter = params.resolved_centering()
    system = params.crystal_system() or "triclinic"
    family = {"triclinic": "a", "monoclinic": "m", "orthorhombic": "o",
             "tetragonal": "t", "trigonal": "h", "hexagonal": "h",
             "cubic": "c"}[system]
    return family + letter


# -- validating the six parameters against the crystal system --------------
# What each crystal system requires (constraint checks only -- nothing here
# alters the numbers, it just says whether they are self-consistent).
def _requirements(system: str):
    """(equal_length_groups, required_angles) for a crystal system.

    ``equal_length_groups`` is a list of tuples of parameter names that must
    be equal; ``required_angles`` maps an angle name to the value it must
    equal (in degrees), for angles that are not free.
    """
    if system == "triclinic":
        return [], {}
    if system == "monoclinic":
        return [], {"alpha": 90.0, "gamma": 90.0}   # beta free (unique axis b)
    if system == "orthorhombic":
        return [], {"alpha": 90.0, "beta": 90.0, "gamma": 90.0}
    if system == "tetragonal":
        return [("a", "b")], {"alpha": 90.0, "beta": 90.0, "gamma": 90.0}
    if system == "hexagonal":
        return [("a", "b")], {"alpha": 90.0, "beta": 90.0, "gamma": 120.0}
    if system == "trigonal":
        # Covers both settings: hexagonal axes (a=b, gamma=120) for hP, and
        # the rhombohedral axes (a=b=c, alpha=beta=gamma) some tables quote
        # for hR. Which applies is decided by the centering, not asked here.
        return [("a", "b")], {"alpha": 90.0, "beta": 90.0, "gamma": 120.0}
    if system == "cubic":
        return [("a", "b"), ("b", "c")], {"alpha": 90.0, "beta": 90.0, "gamma": 90.0}
    raise ValueError(f"unknown crystal system {system!r}")


def validate_lattice_parameters(params: LatticeParams, rel_tol: float = 1e-3):
    """Warnings (not exceptions) about ``params`` being inconsistent with
    its crystal system -- e.g. a cubic space group with a != b.

    Returns a list of human-readable strings, empty if everything is
    consistent (or if the crystal system is unknown, i.e. a bare centering
    letter with no space group). The Brillouin-zone dialog shows these
    rather than silently "correcting" the numbers: the user typed them, and
    which one is wrong (the lengths, or the choice of space group) is their
    call, not this module's.
    """
    system = params.crystal_system()
    if params.resolved_centering() == "R" and system == "trigonal":
        # The rhombohedral-axes alternative (a=b=c, alpha=beta=gamma=/=90)
        # is equally valid for hR and cannot be told apart from the six
        # numbers alone, so trigonal+R is not checked here at all rather
        # than flagging a perfectly good rhombohedral-axes cell as wrong.
        return []
    if system is None:
        return []
    equal_groups, angles = _requirements(system)
    warnings = []
    lengths = {"a": params.a, "b": params.b, "c": params.c}
    for x, y in equal_groups:
        lx, ly = lengths[x], lengths[y]
        if abs(lx - ly) > rel_tol * max(lx, ly, 1e-12):
            warnings.append(
                f"{system} requires {x} = {y}, got {x}={lx:.6g}, {y}={ly:.6g}")
    values = {"alpha": params.alpha, "beta": params.beta, "gamma": params.gamma}
    for name, required in angles.items():
        v = values[name]
        if abs(v - required) > rel_tol * max(required, 1.0) * 90.0 / 90.0 \
                and abs(v - required) > 0.05:
            warnings.append(
                f"{system} requires {name} = {required:g} deg, got {v:.6g} deg")
    return warnings


# -- conventional cell ------------------------------------------------------
def conventional_vectors(params: LatticeParams) -> np.ndarray:
    """The conventional-cell real-space basis vectors as rows of a (3, 3)
    array, for the general triclinic case (every other crystal system is a
    special case of these six parameters).

    Standard construction: a1 along x, a2 in the xy-plane, a3 completing the
    set from its length and its two angles to a1 and a2::

        a1 = (a, 0, 0)
        a2 = (b cos(gamma), b sin(gamma), 0)
        a3 = (c cos(beta), c (cos(alpha) - cos(beta) cos(gamma)) / sin(gamma), c*h)

    with ``h`` fixed by |a3| = c. The cell volume computed from these vectors
    is cross-checked against the closed-form
    ``V = a b c sqrt(1 - cos^2(a) - cos^2(b) - cos^2(g) + 2 cos(a)cos(b)cos(g))``
    (equal to machine precision when the angles describe a valid cell).
    """
    a, b, c = params.a, params.b, params.c
    if a <= 0 or b <= 0 or c <= 0:
        raise ValueError(f"lattice lengths must be positive, got a={a}, b={b}, c={c}")
    alpha, beta, gamma = (np.deg2rad(params.alpha), np.deg2rad(params.beta),
                         np.deg2rad(params.gamma))
    ca, cb, cg, sg = np.cos(alpha), np.cos(beta), np.cos(gamma), np.sin(gamma)
    if abs(sg) < 1e-10:
        raise ValueError(f"gamma = {params.gamma} deg gives a degenerate a-b plane")

    a1 = np.array([a, 0.0, 0.0])
    a2 = np.array([b * cg, b * sg, 0.0])
    a3x = c * cb
    a3y = c * (ca - cb * cg) / sg
    h2 = 1.0 - cb ** 2 - ((ca - cb * cg) / sg) ** 2
    if h2 <= -1e-9:
        raise ValueError(
            f"alpha={params.alpha}, beta={params.beta}, gamma={params.gamma} deg "
            "do not describe a valid unit cell (negative cell volume)")
    a3z = c * np.sqrt(max(h2, 0.0))
    a3 = np.array([a3x, a3y, a3z])

    vectors = np.array([a1, a2, a3])
    volume = float(np.dot(vectors[0], np.cross(vectors[1], vectors[2])))
    volume_closed_form = a * b * c * np.sqrt(max(
        1.0 - ca ** 2 - cb ** 2 - cg ** 2 + 2.0 * ca * cb * cg, 0.0))
    if not np.isclose(volume, volume_closed_form, rtol=1e-6, atol=1e-9 * a * b * c):
        raise AssertionError(
            "internal error: conventional_vectors' volume does not match the "
            f"closed-form triclinic volume ({volume:.6g} vs {volume_closed_form:.6g})")
    return vectors


# -- centering: conventional -> primitive -----------------------------------
# Each row gives one primitive vector as fractional coefficients of the
# conventional (a, b, c). Standard International Tables transformations
# (the same ones ASE, pymatgen and VESTA use); see e.g. ITA Vol. A, Table
# 5.1.3.1. 'R' is the hexagonal-axes "obverse" setting (a=b, gamma=120,
# independent c) -> rhombohedral primitive cell, 1/3 the conventional volume.
CENTERING_TRANSFORMS = {
    "P": np.eye(3),
    "C": np.array([[0.5, -0.5, 0.0], [0.5, 0.5, 0.0], [0.0, 0.0, 1.0]]),
    "A": np.array([[1.0, 0.0, 0.0], [0.0, 0.5, -0.5], [0.0, 0.5, 0.5]]),
    "B": np.array([[0.5, 0.0, -0.5], [0.0, 1.0, 0.0], [0.5, 0.0, 0.5]]),
    "I": np.array([[-0.5, 0.5, 0.5], [0.5, -0.5, 0.5], [0.5, 0.5, -0.5]]),
    "F": np.array([[0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]]),
    "R": np.array([[2.0, 1.0, 1.0], [-1.0, 1.0, 1.0], [-1.0, -2.0, 1.0]]) / 3.0,
}

#: how many conventional cells' worth of lattice points sit in one primitive
#: cell -- i.e. 1 / (conventional volume / primitive volume). Used as a
#: sanity check in tests and by anyone wanting a quick point count.
LATTICE_POINTS_PER_CONVENTIONAL_CELL = {"P": 1, "A": 2, "B": 2, "C": 2,
                                        "I": 2, "F": 4, "R": 3}


def primitive_vectors(params: LatticeParams) -> np.ndarray:
    """The primitive-cell real-space basis, as rows of a (3, 3) array:
    :func:`conventional_vectors` with the centering transform for this
    lattice's Bravais type applied.
    """
    conventional = conventional_vectors(params)
    letter = params.resolved_centering()
    try:
        transform = CENTERING_TRANSFORMS[letter]
    except KeyError:
        raise ValueError(
            f"unknown centering {letter!r}; expected one of "
            f"{sorted(CENTERING_TRANSFORMS)}") from None
    return transform @ conventional


def reciprocal_vectors(vectors: np.ndarray) -> np.ndarray:
    """Reciprocal-lattice vectors (rows) in the physicist's convention
    (2*pi included), matching the Å⁻¹ axes the rest of this program uses.

    ``vectors`` are real-space basis vectors as rows (either
    :func:`conventional_vectors` or :func:`primitive_vectors` -- the
    reciprocal of the primitive cell is what the Brillouin-zone construction
    needs).
    """
    a1, a2, a3 = vectors
    volume = float(np.dot(a1, np.cross(a2, a3)))
    if abs(volume) < 1e-12:
        raise ValueError("degenerate real-space cell (zero volume)")
    b1 = 2.0 * np.pi * np.cross(a2, a3) / volume
    b2 = 2.0 * np.pi * np.cross(a3, a1) / volume
    b3 = 2.0 * np.pi * np.cross(a1, a2) / volume
    return np.array([b1, b2, b3])


# -- which cell parameters the space group leaves free ----------------------
@dataclass(frozen=True)
class CellConstraints:
    """What a space group fixes about the six cell parameters.

    ``mirrors`` maps a length that is *not* independent to the one it must
    equal (``{"b": "a", "c": "a"}`` for cubic); ``fixed_angles`` maps an
    angle to the value symmetry forces on it. Anything not mentioned is
    free. The Brillouin-zone dialog uses this to disable the boxes that are
    not the user's to choose, rather than accepting a cell that contradicts
    the space group and warning about it afterwards.
    """
    mirrors: dict
    fixed_angles: dict
    note: str = ""


def free_parameters(space_group: int) -> CellConstraints:
    """The cell-parameter constraints implied by a space-group number.

    The conventional (not primitive) cell is what is constrained, which is
    what :func:`conventional_vectors` takes -- e.g. a rhombohedral (R) group
    is constrained to the *hexagonal* axes setting (a = b, gamma = 120, c
    free), because that is the setting ``CENTERING_TRANSFORMS["R"]``
    inverts. The equally valid rhombohedral-axes setting (a = b = c,
    alpha = beta = gamma) is a different description of the same lattice and
    is deliberately not offered, so that one set of six numbers always means
    one cell.
    """
    info = spacegroup_info(space_group)
    system = info.crystal_system
    right = {"alpha": 90.0, "beta": 90.0, "gamma": 90.0}
    if system == "triclinic":
        return CellConstraints({}, {}, "all six parameters free")
    if system == "monoclinic":
        return CellConstraints({}, {"alpha": 90.0, "gamma": 90.0},
                               "unique axis b: beta is the free angle")
    if system == "orthorhombic":
        return CellConstraints({}, dict(right), "a, b, c free; all angles 90")
    if system == "tetragonal":
        return CellConstraints({"b": "a"}, dict(right), "a = b, c free")
    if system in ("trigonal", "hexagonal"):
        note = "a = b, c free, gamma = 120"
        if info.centering == "R":
            note += " (hexagonal axes, obverse setting)"
        return CellConstraints({"b": "a"},
                               {"alpha": 90.0, "beta": 90.0, "gamma": 120.0},
                               note)
    if system == "cubic":
        return CellConstraints({"b": "a", "c": "a"}, dict(right),
                               "a = b = c, all angles 90")
    raise ValueError(f"unknown crystal system {system!r}")


# -- point group, in Cartesian reciprocal space -----------------------------
def point_group_operations(params: LatticeParams, laue: bool = True) -> np.ndarray:
    """The space group's point-group rotations as Cartesian 3x3 matrices,
    ``(N, 3, 3)`` -- what acts on a k-vector.

    :mod:`tools.spacegroups` stores them as integer matrices ``W`` in
    conventional-cell *fractional* coordinates, where they are exact. A
    Cartesian position is ``r = A^T x`` for fractional column ``x`` and
    conventional vectors ``A`` (as rows), so ``x -> W x`` is
    ``r -> A^T W (A^T)^-1 r``. Rotations are orthogonal in Cartesian space,
    and reciprocal space is rotated by the same matrices as real space (the
    reciprocal lattice is the same lattice's dual, and an orthogonal map
    commutes with dualising), so these are applied to k directly.

    With ``laue=True`` (the default) inversion is added if the group does not
    already contain it, giving the Laue class. For a non-magnetic crystal
    time reversal makes ``E(k) = E(-k)`` whether or not the crystal itself is
    centrosymmetric, so the Laue class -- not the bare point group -- is what
    the irreducible Brillouin zone is a fundamental domain of.

    Raises if the resulting matrices are not orthogonal, which means the six
    cell parameters contradict the space group (a cubic group with a != b,
    say); a non-orthogonal "rotation" would silently produce a meaningless
    fundamental domain rather than an obviously wrong one.
    """
    if params.space_group is None:
        raise ValueError("point-group operations need a space-group number")
    conventional = conventional_vectors(params).T          # columns = a1,a2,a3
    inverse = np.linalg.inv(conventional)
    cartesian = []
    for w in point_group_rotations(params.space_group):
        r = conventional @ np.asarray(w, dtype=float) @ inverse
        if not np.allclose(r @ r.T, np.eye(3), atol=1e-8):
            raise ValueError(
                "the cell parameters contradict space group "
                f"{params.space_group} ({spacegroup_info(params.space_group).symbol}): "
                "its symmetry operations are not rotations of this cell")
        cartesian.append(r)
    cartesian = np.array(cartesian)

    if laue and not any(np.allclose(r, -np.eye(3), atol=1e-8) for r in cartesian):
        cartesian = np.concatenate([cartesian, -cartesian])
    return cartesian
