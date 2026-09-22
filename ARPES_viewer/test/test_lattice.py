"""Tests for tools/lattice.py: the real/reciprocal lattice vectors behind the
Brillouin-zone tool, and the space-group-consistency checker.
"""
import numpy as np
import pytest

from tools.lattice import (LatticeParams, bravais_symbol, conventional_vectors,
                         primitive_vectors, reciprocal_vectors,
                         validate_lattice_parameters, free_parameters,
                         point_group_operations,
                         LATTICE_POINTS_PER_CONVENTIONAL_CELL)


def _volume(vectors):
    return abs(np.dot(vectors[0], np.cross(vectors[1], vectors[2])))


def _angles_deg(vectors):
    a1, a2, a3 = vectors

    def angle(u, v):
        return np.degrees(np.arccos(np.dot(u, v) / np.linalg.norm(u) / np.linalg.norm(v)))
    return angle(a2, a3), angle(a1, a3), angle(a1, a2)  # alpha, beta, gamma


# -- conventional_vectors: the fix for the .m file's missing-alpha bug -----
def test_conventional_vectors_recovers_all_three_angles_independently():
    """A general triclinic cell with alpha, beta, gamma all different must
    come back out unchanged. ``gen_realspace_basevector.m`` fails exactly
    this case (it has no alpha input and uses cos(beta) in alpha's place),
    though the failure is invisible for alpha = beta = 90 deg -- this test
    uses angles far from that to make sure the general case is right, not
    just the degenerate one.
    """
    params = LatticeParams(5.0, 6.0, 7.0, alpha=80.0, beta=100.0, gamma=110.0)
    vectors = conventional_vectors(params)
    alpha, beta, gamma = _angles_deg(vectors)
    assert alpha == pytest.approx(80.0)
    assert beta == pytest.approx(100.0)
    assert gamma == pytest.approx(110.0)
    lengths = [np.linalg.norm(v) for v in vectors]
    assert lengths == pytest.approx([5.0, 6.0, 7.0])


def test_conventional_vectors_matches_closed_form_volume():
    params = LatticeParams(5.0, 6.0, 7.0, alpha=80.0, beta=100.0, gamma=110.0)
    vectors = conventional_vectors(params)
    a, b, c = 5.0, 6.0, 7.0
    ca, cb, cg = (np.cos(np.deg2rad(x)) for x in (80.0, 100.0, 110.0))
    expected = a * b * c * np.sqrt(1 - ca ** 2 - cb ** 2 - cg ** 2 + 2 * ca * cb * cg)
    assert _volume(vectors) == pytest.approx(expected)


def test_conventional_vectors_rejects_impossible_angles():
    with pytest.raises(ValueError):
        conventional_vectors(LatticeParams(5, 5, 5, alpha=10, beta=10, gamma=170))


# -- centering transforms: volume ratios and orthogonality -----------------
@pytest.mark.parametrize("letter,space_group,params", [
    ("P", 1, LatticeParams(4.0, 5.0, 6.0, 80.0, 100.0, 110.0, space_group=1)),
    ("C", 15, LatticeParams(5.0, 6.0, 7.0, 90.0, 100.0, 90.0, space_group=15)),
    ("I", 229, LatticeParams(4.0, 4.0, 4.0, 90.0, 90.0, 90.0, space_group=229)),
    ("F", 225, LatticeParams(4.0, 4.0, 4.0, 90.0, 90.0, 90.0, space_group=225)),
])
def test_primitive_cell_volume_ratio(letter, space_group, params):
    conv = conventional_vectors(params)
    prim = primitive_vectors(params)
    ratio = _volume(conv) / _volume(prim)
    assert ratio == pytest.approx(LATTICE_POINTS_PER_CONVENTIONAL_CELL[letter])


def test_rhombohedral_centering_volume_ratio():
    params = LatticeParams(3.0, 3.0, 10.0, 90.0, 90.0, 120.0, space_group=146)  # R3
    assert bravais_symbol(params) == "hR"
    conv = conventional_vectors(params)
    prim = primitive_vectors(params)
    assert _volume(conv) / _volume(prim) == pytest.approx(3.0)


def test_reciprocal_vectors_are_dual_to_real_space():
    params = LatticeParams(4.0, 5.0, 6.0, 80.0, 100.0, 110.0, space_group=1)
    prim = primitive_vectors(params)
    recip = reciprocal_vectors(prim)
    assert prim @ recip.T == pytest.approx(2 * np.pi * np.eye(3))


# -- validation against the crystal system ----------------------------------
def test_validate_flags_inconsistent_cubic_lengths():
    params = LatticeParams(4.0, 4.5, 4.5, 90, 90, 90, space_group=225)  # Fm-3m
    warnings = validate_lattice_parameters(params)
    assert len(warnings) == 1
    assert "a = b" in warnings[0]


def test_validate_accepts_consistent_hexagonal():
    params = LatticeParams(2.46, 2.46, 10.0, 90, 90, 120, space_group=194)
    assert validate_lattice_parameters(params) == []


def test_validate_leaves_monoclinic_beta_free():
    params = LatticeParams(5, 6, 7, 90, 73.5, 90, space_group=15)  # C2/c
    assert validate_lattice_parameters(params) == []


def test_validate_skips_rhombohedral_axes_setting():
    """A trigonal R space group in rhombohedral axes (a=b=c, alpha=beta=
    gamma != 90) is equally valid to the hexagonal-axes setting and cannot
    be told apart from the six numbers alone, so it must not be flagged."""
    params = LatticeParams(5, 5, 5, 60, 60, 60, space_group=146)  # R3
    assert validate_lattice_parameters(params) == []


def test_bravais_symbol_from_manual_centering_without_space_group():
    """A bare centering letter with no space group cannot by itself say
    which crystal system it belongs to (F-centering exists for both cubic
    and orthorhombic, with different length/angle constraints) -- so
    ``bravais_symbol`` falls back to the family-less "a" (triclinic) prefix
    rather than guessing, and there is nothing for the validator to check.
    """
    params = LatticeParams(4, 4, 4, 90, 90, 90, centering="F")
    assert bravais_symbol(params) == "aF"
    assert params.crystal_system() is None
    assert validate_lattice_parameters(
        LatticeParams(4, 5, 6, 90, 90, 90, centering="F")) == []


# -- what the space group fixes about the cell ------------------------------
@pytest.mark.parametrize("space_group,mirrors,angles", [
    (1, {}, {}),                                              # P1
    (5, {}, {"alpha": 90.0, "gamma": 90.0}),                  # C2, unique axis b
    (62, {}, {"alpha": 90.0, "beta": 90.0, "gamma": 90.0}),   # Pnma
    (139, {"b": "a"}, {"alpha": 90.0, "beta": 90.0, "gamma": 90.0}),   # I4/mmm
    (194, {"b": "a"}, {"alpha": 90.0, "beta": 90.0, "gamma": 120.0}),  # P6_3/mmc
    (166, {"b": "a"}, {"alpha": 90.0, "beta": 90.0, "gamma": 120.0}),  # R-3m
    (225, {"b": "a", "c": "a"}, {"alpha": 90.0, "beta": 90.0, "gamma": 90.0}),
])
def test_free_parameters_match_the_crystal_system(space_group, mirrors, angles):
    constraints = free_parameters(space_group)
    assert constraints.mirrors == mirrors
    assert constraints.fixed_angles == angles


def test_free_parameters_constrain_rhombohedral_to_hexagonal_axes():
    """R-3m has two equally valid descriptions; only the hexagonal-axes one
    is offered, because that is the setting CENTERING_TRANSFORMS['R']
    inverts -- and the note says so rather than leaving it implicit."""
    constraints = free_parameters(166)
    assert constraints.fixed_angles["gamma"] == 120.0
    assert "hexagonal axes" in constraints.note


def test_a_cell_built_to_the_constraints_passes_validation():
    """The two halves agree: a cell that satisfies free_parameters has
    nothing for validate_lattice_parameters to complain about. This is what
    lets the dialog disable the constrained boxes instead of warning."""
    for space_group in (1, 5, 62, 139, 194, 225):
        constraints = free_parameters(space_group)
        values = {"a": 4.0, "b": 5.0, "c": 6.0,
                  "alpha": 88.0, "beta": 95.0, "gamma": 103.0}
        values.update(constraints.fixed_angles)
        for target, source in constraints.mirrors.items():
            values[target] = values[source]
        params = LatticeParams(space_group=space_group, **values)
        assert validate_lattice_parameters(params) == []


# -- the point group, in Cartesian reciprocal space --------------------------
@pytest.mark.parametrize("params,bare,laue", [
    # P1 has only the identity, but time reversal still gives E(k) = E(-k),
    # so its Laue class is -1: the one case where the two orders differ by
    # more than nothing.
    (LatticeParams(4.0, 5.0, 6.0, 88.0, 95.0, 103.0, space_group=1), 1, 2),
    (LatticeParams(4.0, 4.0, 4.0, space_group=195), 12, 24),   # P23 -> m-3
    (LatticeParams(4.0, 4.0, 4.0, space_group=221), 48, 48),   # Pm-3m
    (LatticeParams(4.0, 4.0, 4.0, space_group=225), 48, 48),   # Fm-3m
    (LatticeParams(2.46, 2.46, 10.0, 90.0, 90.0, 120.0, space_group=194), 24, 24),
    (LatticeParams(5.0, 6.0, 4.0, 90.0, 105.0, 90.0, space_group=12), 4, 4),
])
def test_point_group_operation_counts(params, bare, laue):
    assert len(point_group_operations(params, laue=False)) == bare
    assert len(point_group_operations(params)) == laue


def test_point_group_operations_are_orthogonal_and_close():
    """They must be genuine rotations in Cartesian space (not merely integer
    matrices in the lattice basis) and must form a group -- that is what the
    irreducible-wedge construction relies on."""
    params = LatticeParams(2.46, 2.46, 10.0, 90.0, 90.0, 120.0, space_group=194)
    ops = point_group_operations(params)
    for op in ops:
        assert np.allclose(op @ op.T, np.eye(3), atol=1e-9)
    known = {tuple(np.round(op.ravel(), 6)) for op in ops}
    for x in ops:
        for y in ops:
            assert tuple(np.round((x @ y).ravel(), 6)) in known


def test_point_group_operations_map_the_reciprocal_lattice_onto_itself():
    """The real check that the basis conversion is right: every operation
    must send each reciprocal-lattice vector to another one, i.e. to integer
    coefficients in the same basis."""
    params = LatticeParams(4.0, 4.0, 4.0, space_group=225)
    b = reciprocal_vectors(primitive_vectors(params))
    for op in point_group_operations(params):
        coefficients = np.linalg.solve(b.T, (op @ b.T))
        assert np.allclose(coefficients, np.round(coefficients), atol=1e-9)


def test_point_group_operations_reject_a_cell_that_contradicts_the_group():
    """A cubic space group with a != b is not a cubic cell; the operations
    would not be rotations of it, and saying so is better than returning
    matrices that quietly are not."""
    params = LatticeParams(4.0, 5.0, 4.0, space_group=225)
    with pytest.raises(ValueError, match="contradict"):
        point_group_operations(params)
