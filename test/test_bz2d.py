"""Tests for tools/bz2d.py: the general 2-D Wigner-Seitz cell (generalizing
``plotBZ_demo_old``'s hard-coded-hexagonal version) and its tiling."""
import numpy as np
import pytest

from tools.bz2d import (reciprocal_vectors_2d, wigner_seitz_cell_2d, tile_2d,
                      lattice_point_group_2d, irreducible_cell_2d)


def _area(polygon):
    x, y = polygon[:, 0], polygon[:, 1]
    return 0.5 * abs(np.sum(x[:-1] * y[1:] - x[1:] * y[:-1]))


def _hex_vectors(a, rotation_deg=0.0):
    r = np.deg2rad(rotation_deg)
    rot = np.array([[np.cos(r), -np.sin(r)], [np.sin(r), np.cos(r)]])
    a1 = rot @ np.array([a, 0.0])
    a2 = rot @ np.array([a * np.cos(np.deg2rad(120.0)), a * np.sin(np.deg2rad(120.0))])
    return a1, a2


def test_hexagonal_cell_has_six_sides_and_known_high_symmetry_distances():
    a = 2.46
    a1, a2 = _hex_vectors(a)
    g1, g2 = reciprocal_vectors_2d(a1, a2)
    poly = wigner_seitz_cell_2d(g1, g2)
    assert len(poly) - 1 == 6
    # Gamma-K (a vertex) and Gamma-M (an edge midpoint), textbook results
    gamma_k = np.linalg.norm(poly[0])
    gamma_m = np.linalg.norm((poly[0] + poly[1]) / 2.0)
    assert gamma_k == pytest.approx(4 * np.pi / (3 * a))
    assert gamma_m == pytest.approx(2 * np.pi / (a * np.sqrt(3)))


def test_square_cell_is_a_square_of_known_half_width():
    a = 3.0
    a1, a2 = np.array([a, 0.0]), np.array([0.0, a])
    g1, g2 = reciprocal_vectors_2d(a1, a2)
    poly = wigner_seitz_cell_2d(g1, g2)
    assert len(poly) - 1 == 4
    assert np.max(np.abs(poly)) == pytest.approx(np.pi / a)


def test_rectangular_cell_has_two_different_half_widths():
    a, b = 3.0, 5.0
    g1, g2 = reciprocal_vectors_2d(np.array([a, 0.0]), np.array([0.0, b]))
    poly = wigner_seitz_cell_2d(g1, g2)
    assert np.max(np.abs(poly[:, 0])) == pytest.approx(np.pi / a)
    assert np.max(np.abs(poly[:, 1])) == pytest.approx(np.pi / b)


def test_oblique_cell_area_matches_reciprocal_unit_cell_area():
    """Whatever its shape, the Wigner-Seitz cell's area must equal
    |g1 x g2| -- a different-shaped primitive cell of the same lattice."""
    a1 = np.array([3.0, 0.0])
    a2 = np.array([1.2, 4.1])
    g1, g2 = reciprocal_vectors_2d(a1, a2)
    poly = wigner_seitz_cell_2d(g1, g2)
    # shoelace formula
    x, y = poly[:, 0], poly[:, 1]
    area = 0.5 * abs(np.sum(x[:-1] * y[1:] - x[1:] * y[:-1]))
    expected = abs(g1[0] * g2[1] - g1[1] * g2[0])
    assert area == pytest.approx(expected, rel=1e-9)


def test_tile_2d_covers_a_regular_grid_for_square_lattice():
    a = 2.0
    g1, g2 = reciprocal_vectors_2d(np.array([a, 0.0]), np.array([0.0, a]))
    poly = wigner_seitz_cell_2d(g1, g2)
    spacing = np.pi / a * 2  # full square width
    window = (-2 * spacing, 2 * spacing)
    tiles = tile_2d(poly, g1, g2, window, window)
    centroids = {(round(t[:-1, 0].mean(), 6), round(t[:-1, 1].mean(), 6)) for t in tiles}
    xs = sorted({c[0] for c in centroids})
    assert np.allclose(np.diff(xs), spacing, atol=1e-6)


# -- the plane lattice's own point group, and the irreducible wedge -----------
def _oblique(a, b, gamma_deg):
    return (np.array([a, 0.0]),
            np.array([b * np.cos(np.deg2rad(gamma_deg)),
                      b * np.sin(np.deg2rad(gamma_deg))]))


@pytest.mark.parametrize("vectors,order", [
    # The five plane Bravais lattices and their holohedries: oblique (2),
    # rectangular and centred-rectangular (2mm), square (4mm), hexagonal
    # (6mm). The centred-rectangular case is given here the way a user
    # would type it -- equal lengths and an odd angle -- which is exactly
    # the case a "is gamma 90 or 120?" classifier would call oblique.
    (_oblique(3.0, 4.1, 73.0), 2),
    (_oblique(3.0, 5.0, 90.0), 4),
    (_oblique(3.0, 3.0, 100.0), 4),
    (_oblique(3.0, 3.0, 90.0), 8),
    (_oblique(2.46, 2.46, 120.0), 12),
])
def test_plane_lattice_point_group_orders(vectors, order):
    g1, g2 = reciprocal_vectors_2d(*vectors)
    ops = lattice_point_group_2d(g1, g2)
    assert len(ops) == order
    # and it really is a group of orthogonal maps
    for op in ops:
        assert np.allclose(op @ op.T, np.eye(2), atol=1e-9)
    known = {tuple(np.round(op.ravel(), 9)) for op in ops}
    for x in ops:
        for y in ops:
            assert tuple(np.round((x @ y).ravel(), 9)) in known


@pytest.mark.parametrize("vectors,order", [
    (_oblique(3.0, 4.1, 73.0), 2),
    (_oblique(3.0, 5.0, 90.0), 4),
    (_oblique(3.0, 3.0, 90.0), 8),
    (_oblique(2.46, 2.46, 120.0), 12),
])
def test_irreducible_cell_area_is_the_zone_area_over_the_group_order(vectors, order):
    g1, g2 = reciprocal_vectors_2d(*vectors)
    zone = wigner_seitz_cell_2d(g1, g2)
    wedge = irreducible_cell_2d(zone, lattice_point_group_2d(g1, g2))
    assert _area(wedge) == pytest.approx(_area(zone) / order, rel=1e-9)


def test_hexagonal_irreducible_cell_is_the_gamma_m_k_triangle():
    """The 1/12 wedge of a hexagonal zone is the triangle Gamma-M-K, whose
    two non-zero corners sit at the textbook distances."""
    a = 2.46
    g1, g2 = reciprocal_vectors_2d(*_oblique(a, a, 120.0))
    wedge = irreducible_cell_2d(wigner_seitz_cell_2d(g1, g2),
                                lattice_point_group_2d(g1, g2))
    assert len(wedge) - 1 == 3
    radii = sorted(np.linalg.norm(wedge[:-1], axis=1))
    assert radii == pytest.approx([0.0,
                                   2 * np.pi / (a * np.sqrt(3)),   # Gamma-M
                                   4 * np.pi / (3 * a)],           # Gamma-K
                                  abs=1e-9)


def test_square_irreducible_cell_is_the_gamma_x_m_triangle():
    a = 3.0
    g1, g2 = reciprocal_vectors_2d(*_oblique(a, a, 90.0))
    wedge = irreducible_cell_2d(wigner_seitz_cell_2d(g1, g2),
                                lattice_point_group_2d(g1, g2))
    assert len(wedge) - 1 == 3
    radii = sorted(np.linalg.norm(wedge[:-1], axis=1))
    assert radii == pytest.approx([0.0,
                                   np.pi / a,                      # Gamma-X
                                   np.sqrt(2) * np.pi / a],        # Gamma-M
                                  abs=1e-9)


def test_irreducible_cell_lies_inside_the_zone_it_came_from():
    g1, g2 = reciprocal_vectors_2d(*_oblique(2.46, 2.46, 120.0))
    zone = wigner_seitz_cell_2d(g1, g2)
    wedge = irreducible_cell_2d(zone, lattice_point_group_2d(g1, g2))
    # a convex zone contains a point iff the point is behind all its edges
    for i in range(len(zone) - 1):
        edge = zone[i + 1] - zone[i]
        outward = np.array([edge[1], -edge[0]])
        outward = outward / np.linalg.norm(outward)
        if outward @ (zone[i] - zone[:-1].mean(axis=0)) < 0:
            outward = -outward
        assert np.all((wedge - zone[i]) @ outward <= 1e-9)
