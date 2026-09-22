"""Tests for tools/bz3d.py against closed-form 3-D Brillouin zone shapes, and
the plane-cut / tiling routines built on top of it.

The four cubic-family cases below are textbook results (e.g. Ashcroft &
Mermin ch. 9): simple cubic's BZ is a cube; body-centred cubic's reciprocal
lattice is face-centred cubic, whose Wigner-Seitz cell is a truncated
octahedron (14 faces, 24 vertices); face-centred cubic's reciprocal lattice
is body-centred cubic, whose cell is a rhombic dodecahedron (12 faces, 14
vertices) -- this is also exactly why ``gen_brillouin.m``'s fixed +-1
neighbour shell is not a universal shortcut: it happens to be enough for
these high-symmetry lattices, but nothing in this test relies on that,
since :func:`wigner_seitz_cell` grows the shell itself.
"""
import numpy as np
import pytest

from tools.lattice import (LatticeParams, primitive_vectors, reciprocal_vectors,
                         point_group_operations)
from tools.bz3d import (wigner_seitz_cell, plane_cut, plane_basis, tile_and_cut,
                      face_planes, irreducible_wedge, cut_points_3d)


def _bz(space_group, a=4.0, c=None, gamma=90.0):
    c = a if c is None else c
    params = LatticeParams(a, a, c, 90.0, 90.0, gamma, space_group=space_group)
    prim = primitive_vectors(params)
    b = reciprocal_vectors(prim)
    return wigner_seitz_cell(b), b, prim


@pytest.mark.parametrize("space_group,n_faces,n_vertices", [
    (195, 6, 8),     # cubic P -> cube
    (229, 12, 14),   # cubic I (bcc) -> rhombic dodecahedron
    (225, 14, 24),   # cubic F (fcc) -> truncated octahedron
])
def test_cubic_family_shapes(space_group, n_faces, n_vertices):
    bz, _, _ = _bz(space_group)
    assert len(bz.faces) == n_faces
    assert len(bz.vertices) == n_vertices


def test_hexagonal_prism():
    bz, _, _ = _bz(194, a=2.46, c=10.0, gamma=120.0)
    assert len(bz.faces) == 8       # 2 hexagons + 6 rectangles
    assert len(bz.vertices) == 12


@pytest.mark.parametrize("space_group,a,c,gamma", [
    (195, 4.0, 4.0, 90.0), (229, 4.0, 4.0, 90.0), (225, 4.0, 4.0, 90.0),
    (194, 2.46, 10.0, 120.0),
])
def test_bz_volume_matches_reciprocal_cell_volume(space_group, a, c, gamma):
    """The Wigner-Seitz cell's volume must equal the primitive reciprocal
    cell's own volume |b1.(b2 x b3)| -- it is a different-shaped choice of
    primitive cell for the same lattice, not a different lattice."""
    bz, b, _ = _bz(space_group, a=a, c=c, gamma=gamma)
    expected = abs(np.dot(b[0], np.cross(b[1], b[2])))
    assert bz.volume == pytest.approx(expected, rel=1e-9)


def test_shell_growth_needed_for_oblique_cell():
    """A very oblique triclinic cell needs more than the single +-1
    neighbour shell ``gen_brillouin.m`` always uses -- this is the failure
    mode that motivates growing the shell automatically."""
    params = LatticeParams(4.0, 4.0, 4.0, alpha=75.0, beta=75.0, gamma=75.0,
                           space_group=1)
    prim = primitive_vectors(params)
    b = reciprocal_vectors(prim)
    bz = wigner_seitz_cell(b)
    assert bz.shells_used >= 2
    # and the volume is still exactly the reciprocal cell volume
    expected = abs(np.dot(b[0], np.cross(b[1], b[2])))
    assert bz.volume == pytest.approx(expected, rel=1e-9)


# -- plane_cut ---------------------------------------------------------------
def test_plane_cut_cubic_p_gives_square_half_width_pi_over_a():
    a = 4.0
    bz, b, _ = _bz(195, a=a)
    u, v = plane_cut(bz.faces, np.array([0.0, 0.0, 1.0]), np.array([0.0, 0.0, 0.0]))
    assert u is not None
    assert len(u) == 5  # square, closed (4 + repeat)
    assert np.max(np.abs(u)) == pytest.approx(np.pi / a)
    assert np.max(np.abs(v)) == pytest.approx(np.pi / a)


def test_plane_cut_fcc_001_gives_octagon_at_2pi_over_a():
    """The (001) cross-section of the truncated-octahedron FCC zone through
    Gamma is the well-known "square with the corners cut off" octagon,
    reaching the X point at 2*pi/a along <100>."""
    a = 4.0
    bz, b, _ = _bz(225, a=a)
    u, v = plane_cut(bz.faces, np.array([0.0, 0.0, 1.0]), np.array([0.0, 0.0, 0.0]))
    assert len(u) - 1 == 8
    assert np.max(np.abs(u)) == pytest.approx(2 * np.pi / a)


def test_plane_cut_far_outside_the_zone_is_empty():
    bz, b, _ = _bz(195, a=4.0)
    u, v = plane_cut(bz.faces, np.array([0.0, 0.0, 1.0]), np.array([0.0, 0.0, 100.0]))
    assert u is None and v is None


# -- tile_and_cut -------------------------------------------------------------
def test_tile_and_cut_reproduces_a_regular_grid_for_cubic_p():
    a = 4.0
    bz, b, _ = _bz(195, a=a)
    spacing = 2 * np.pi / a
    window = (-3 * spacing, 3 * spacing)
    polygons = tile_and_cut(bz, b, np.array([0.0, 0.0, 1.0]), np.array([0.0, 0.0, 0.0]),
                            window, window)
    centroids = sorted({(round(u.mean(), 6), round(v.mean(), 6)) for u, v in polygons})
    # every centroid should sit on the expected square grid, spacing apart
    xs = sorted({c[0] for c in centroids})
    diffs = np.diff(xs)
    assert np.allclose(diffs, spacing, atol=1e-6)
    # and every polygon is a unit square of the right size
    for u, v in polygons:
        assert (u.max() - u.min()) == pytest.approx(spacing)
        assert (v.max() - v.min()) == pytest.approx(spacing)


# -- irreducible wedge --------------------------------------------------------
def _zone_and_ops(params):
    prim = primitive_vectors(params)
    b = reciprocal_vectors(prim)
    return wigner_seitz_cell(b), point_group_operations(params), b


@pytest.mark.parametrize("params,order", [
    # The defining property: the wedge repeated by the |Laue| operations of
    # the crystal fills the zone exactly once, so its volume is the zone's
    # divided by that order. Pm-3m/Fm-3m/Im-3m are all Laue class m-3m (48);
    # P6_3/mmc is 6/mmm (24); P1 has only the identity, but time reversal
    # still gives inversion, so its Laue class is -1 (2).
    (LatticeParams(4.0, 4.0, 4.0, space_group=221), 48),
    (LatticeParams(4.0, 4.0, 4.0, space_group=225), 48),
    (LatticeParams(4.0, 4.0, 4.0, space_group=229), 48),
    (LatticeParams(2.46, 2.46, 10.0, 90.0, 90.0, 120.0, space_group=194), 24),
    (LatticeParams(4.0, 4.0, 6.0, space_group=139), 16),
    (LatticeParams(3.0, 3.0, 20.0, 90.0, 90.0, 120.0, space_group=166), 12),
    (LatticeParams(5.0, 6.0, 4.0, 90.0, 105.0, 90.0, space_group=12), 4),
    (LatticeParams(4.0, 5.0, 6.0, 88.0, 95.0, 103.0, space_group=1), 2),
])
def test_irreducible_wedge_volume_is_zone_volume_over_group_order(params, order):
    zone, ops, _ = _zone_and_ops(params)
    assert len(ops) == order
    wedge = irreducible_wedge(zone, ops)
    assert wedge.volume == pytest.approx(zone.volume / order, rel=1e-9)


def test_fcc_irreducible_wedge_is_the_textbook_six_vertex_shape():
    """The 1/48 wedge of the FCC zone is the familiar Gamma-X-W-K-U-L solid:
    six corners, five faces, and those corners are exactly the standard
    high-symmetry points, whose distances from Gamma are known in closed
    form (Ashcroft & Mermin ch. 9 / any band-structure path table):

        Gamma 0     L (pi/a)(1,1,1)        X (2pi/a)(1,0,0)
        K (2pi/a)(3/4,3/4,0)   U (2pi/a)(1,1/4,1/4)   W (2pi/a)(1,1/2,0)

    K and U are the same distance out (they are one reciprocal-lattice
    translation apart), which is why two of the six radii coincide.
    """
    a = 4.0
    params = LatticeParams(a, a, a, space_group=225)
    zone, ops, _ = _zone_and_ops(params)
    wedge = irreducible_wedge(zone, ops)
    assert len(wedge.vertices) == 6
    assert len(wedge.faces) == 5

    k = 2 * np.pi / a
    expected = sorted([
        0.0,                                  # Gamma
        np.sqrt(3) * np.pi / a,               # L
        k,                                    # X
        k * 0.75 * np.sqrt(2),                # K
        k * np.sqrt(1 + 2 * 0.25 ** 2),       # U
        k * np.sqrt(1 + 0.5 ** 2),            # W
    ])
    assert sorted(np.linalg.norm(wedge.vertices, axis=1)) == pytest.approx(
        expected, abs=1e-9)


def test_irreducible_wedge_of_a_zone_sits_inside_that_zone():
    params = LatticeParams(2.46, 2.46, 10.0, 90.0, 90.0, 120.0, space_group=194)
    zone, ops, _ = _zone_and_ops(params)
    wedge = irreducible_wedge(zone, ops)
    normals, offsets = face_planes(zone.faces)
    # every wedge vertex satisfies every one of the zone's own face planes
    assert np.all(wedge.vertices @ normals.T <= offsets[None, :] + 1e-9)


def test_identity_only_operations_leave_the_zone_alone():
    params = LatticeParams(4.0, 5.0, 6.0, 88.0, 95.0, 103.0, space_group=1)
    zone, _, _ = _zone_and_ops(params)
    same = irreducible_wedge(zone, np.eye(3)[None, :, :])
    assert same.volume == pytest.approx(zone.volume)


def test_face_planes_round_trips_the_zone_it_came_from():
    """Each face plane must touch the zone (offset equal to the largest
    projection of any vertex onto its normal) and contain no vertex beyond
    it -- i.e. they really are the zone's supporting half-spaces."""
    zone, _, _ = _bz(225, a=4.0)
    normals, offsets = face_planes(zone.faces)
    projections = zone.vertices @ normals.T
    assert np.all(projections <= offsets[None, :] + 1e-9)
    assert np.allclose(projections.max(axis=0), offsets, atol=1e-9)


# -- cut_points_3d ------------------------------------------------------------
def test_cut_points_3d_puts_the_cut_back_where_it_came_from():
    bz, _, _ = _bz(225, a=4.0)
    normal = np.array([1.0, 1.0, 1.0]) / np.sqrt(3)
    point = normal * 0.2
    u, v = plane_cut(bz.faces, normal, point)
    points = cut_points_3d(u, v, normal, point)
    # every reconstructed point lies in the cutting plane...
    assert np.allclose(points @ normal, np.dot(normal, point), atol=1e-9)
    # ...and projecting it back gives the (u, v) it started from
    u_axis, v_axis = plane_basis(normal)
    assert np.allclose((points - point) @ u_axis, u, atol=1e-9)
    assert np.allclose((points - point) @ v_axis, v, atol=1e-9)


def test_plane_cut_offset_along_the_normal_shrinks_a_cubic_zone_linearly():
    """Sliding the plane along its normal is the parameter the dialog offers
    for moving off Gamma; for a cubic zone cut along (0,0,1) the square
    cross-section keeps its size until the plane leaves the zone at pi/a."""
    a = 4.0
    bz, _, _ = _bz(195, a=a)
    normal = np.array([0.0, 0.0, 1.0])
    for offset in (0.0, 0.3, 0.6):
        u, v = plane_cut(bz.faces, normal, normal * offset)
        assert np.max(np.abs(u)) == pytest.approx(np.pi / a)
    beyond = plane_cut(bz.faces, normal, normal * (np.pi / a + 0.05))
    assert beyond == (None, None)


def test_plane_basis_is_orthonormal_and_perpendicular_to_normal():
    for normal in [np.array([0.0, 0.0, 1.0]), np.array([1.0, 1.0, 1.0]) / np.sqrt(3),
                  np.array([0.3, -0.7, 0.2])]:
        normal = normal / np.linalg.norm(normal)
        u, v = plane_basis(normal)
        assert np.linalg.norm(u) == pytest.approx(1.0)
        assert np.linalg.norm(v) == pytest.approx(1.0)
        assert np.dot(u, v) == pytest.approx(0.0, abs=1e-9)
        assert np.dot(u, normal) == pytest.approx(0.0, abs=1e-9)
        assert np.dot(v, normal) == pytest.approx(0.0, abs=1e-9)


def test_plane_cut_along_one_of_the_polyhedrons_own_faces_returns_that_face():
    """A plane that contains a face cuts the solid exactly in that face. No
    edge crosses the plane there, so the sign-change walk finds nothing and
    the coplanar case has to be caught first -- which is the *normal* case
    for an irreducible wedge, whose faces are the crystal's mirror planes.
    """
    a = 4.0
    bz, _, _ = _bz(195, a=a)                       # simple cubic: a cube
    normal = np.array([0.0, 0.0, 1.0])
    u, v = plane_cut(bz.faces, normal, normal * (np.pi / a))   # the top face
    assert u is not None
    assert len(u) - 1 == 4
    assert np.max(np.abs(u)) == pytest.approx(np.pi / a)
    assert np.max(np.abs(v)) == pytest.approx(np.pi / a)


def test_irreducible_wedge_cuts_along_its_own_mirror_plane():
    """The FCC wedge has a face in the kz = 0 plane -- cutting along it is
    the first thing anyone does with a wedge, and it must come back as that
    face rather than as an empty cut.

    That face is the quadrilateral Gamma-X-W-K (L and U are the two corners
    with kz =/= 0), so its corners sit at the known distances for those four
    points.
    """
    a = 4.0
    params = LatticeParams(a, a, a, space_group=225)
    zone, ops, _ = _zone_and_ops(params)
    wedge = irreducible_wedge(zone, ops)
    u, v = plane_cut(wedge.faces, np.array([0.0, 0.0, 1.0]), np.zeros(3))
    assert u is not None
    assert len(u) - 1 == 4

    k = 2 * np.pi / a
    assert sorted(np.hypot(u[:-1], v[:-1])) == pytest.approx(sorted([
        0.0,                          # Gamma
        k,                            # X   (1, 0, 0)
        k * 0.75 * np.sqrt(2),        # K   (3/4, 3/4, 0)
        k * np.sqrt(1 + 0.5 ** 2),    # W   (1, 1/2, 0)
    ]), abs=1e-9)
