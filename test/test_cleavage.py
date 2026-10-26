"""Reading a cleavage plane off a measured k_z period.

The tests are against lattices whose answers are known by hand. The one
that matters most is the centred case: in a body-centred cubic crystal the
period along [001] is ``4*pi/a``, not ``2*pi/a``, because the (001) planes
of the conventional cell are not all lattice planes. Getting that wrong
gives an answer that is exactly a factor of two out and looks entirely
plausible.
"""
import numpy as np
import pytest

from tools import cleavage
from tools.lattice import LatticeParams


def cubic(a, space_group):
    return LatticeParams(a=a, b=a, c=a, space_group=space_group)


def hexagonal(a, c, space_group=194):
    return LatticeParams(a=a, b=a, c=c, alpha=90.0, beta=90.0, gamma=120.0,
                         space_group=space_group)


def along(found, direction):
    """The entry whose (hkl) is ``direction``, or None."""
    for hkl, length, _unit in found:
        if tuple(int(v) for v in hkl) == tuple(direction):
            return length
    return None


# --------------------------------------------------------------------------
# The reciprocal lengths themselves
# --------------------------------------------------------------------------
def test_a_primitive_cube_repeats_every_two_pi_over_a():
    found = cleavage.reciprocal_lengths(cubic(4.0, 221), max_index=2)
    assert along(found, (0, 0, 1)) == pytest.approx(2 * np.pi / 4.0)
    assert along(found, (1, 0, 0)) == pytest.approx(2 * np.pi / 4.0)


def test_body_centring_halves_the_repeat_along_the_axes():
    """The real point of enumerating the primitive reciprocal lattice.

    (001) is a forbidden reflection for an I lattice; the first one along
    that normal is (002), so the k_z period is 4*pi/a.
    """
    found = cleavage.reciprocal_lengths(cubic(4.0, 229), max_index=2)
    assert along(found, (0, 0, 1)) is None
    assert along(found, (0, 0, 2)) == pytest.approx(4 * np.pi / 4.0)


def test_face_centring_does_the_same():
    found = cleavage.reciprocal_lengths(cubic(5.4, 225), max_index=2)
    assert along(found, (0, 0, 1)) is None
    assert along(found, (0, 0, 2)) == pytest.approx(4 * np.pi / 5.4)
    # ... but (111) is allowed for F, and is the close-packed direction
    assert along(found, (1, 1, 1)) == pytest.approx(
        2 * np.pi * np.sqrt(3) / 5.4)


def test_a_layered_hexagonal_crystal_repeats_with_c():
    found = cleavage.reciprocal_lengths(hexagonal(3.2, 6.0), max_index=2)
    assert along(found, (0, 0, 1)) == pytest.approx(2 * np.pi / 6.0)


def test_every_direction_appears_once_with_its_shortest_vector():
    found = cleavage.reciprocal_lengths(cubic(4.0, 221), max_index=3)
    directions = [tuple(np.round(unit, 6)) for _hkl, _length, unit in found]
    assert len(directions) == len(set(directions))
    # (002) is the second order of (001) and must not be listed separately
    assert along(found, (0, 0, 2)) is None


def test_the_indices_are_never_written_with_a_leading_minus():
    found = cleavage.reciprocal_lengths(hexagonal(3.2, 6.0), max_index=2)
    for hkl, _length, _unit in found:
        first = next((v for v in hkl if v), 0)
        assert first > 0, hkl


# --------------------------------------------------------------------------
# Matching a measurement
# --------------------------------------------------------------------------
def test_the_right_plane_is_found_from_its_own_period():
    params = hexagonal(3.2, 6.0)
    found = cleavage.candidates(2 * np.pi / 6.0, params)
    assert found
    assert found[0].hkl == (0, 0, 1)
    assert abs(found[0].error) < 1e-6
    assert found[0].spacing == pytest.approx(6.0)


def test_a_period_slightly_off_still_matches():
    found = cleavage.candidates(1.05, hexagonal(3.2, 6.0), tolerance=0.15)
    assert found[0].hkl == (0, 0, 1)
    assert abs(found[0].error) < 0.01


def test_a_period_well_outside_the_window_does_not():
    assert cleavage.candidates(0.2, hexagonal(3.2, 6.0), tolerance=0.15) == []


def test_picking_two_zones_apart_is_recognised():
    found = cleavage.candidates(2 * 2 * np.pi / 6.0, hexagonal(3.2, 6.0))
    best = found[0]
    assert best.hkl == (0, 0, 1)
    assert best.orders == 2
    assert "2 zones apart" in best.describe()


def test_the_tolerance_is_respected():
    params = hexagonal(3.2, 6.0)
    exact = 2 * np.pi / 6.0
    assert cleavage.candidates(exact * 1.10, params, tolerance=0.15)
    assert not any(c.hkl == (0, 0, 1)
                   for c in cleavage.candidates(exact * 1.10, params,
                                                tolerance=0.05))


def test_each_plane_appears_once_however_many_orders_fit():
    found = cleavage.candidates(2 * np.pi / 6.0, hexagonal(3.2, 6.0),
                                max_orders=3)
    names = [c.hkl for c in found]
    assert len(names) == len(set(names))


def test_the_results_are_ordered_by_how_well_they_fit():
    found = cleavage.candidates(1.1, hexagonal(3.2, 6.0), tolerance=0.5)
    errors = [abs(c.error) for c in found]
    assert errors == sorted(errors)


def test_a_centred_lattice_reports_the_family_it_belongs_to():
    found = cleavage.candidates(4 * np.pi / 4.0, cubic(4.0, 229))
    best = [c for c in found if c.hkl == (0, 0, 2)]
    assert best
    assert best[0].family == (0, 0, 1)
    assert "family" in best[0].describe()


def test_an_uncentred_plane_does_not_mention_a_family():
    found = cleavage.candidates(2 * np.pi / 6.0, hexagonal(3.2, 6.0))
    assert found[0].family == (0, 0, 1)
    assert "family" not in found[0].describe()


@pytest.mark.parametrize("bad", [0.0, -1.0, float("nan")])
def test_a_meaningless_distance_is_refused(bad):
    with pytest.raises(ValueError, match="not something to match"):
        cleavage.candidates(bad, hexagonal(3.2, 6.0))


def test_the_description_carries_the_numbers_someone_would_check():
    found = cleavage.candidates(1.05, hexagonal(3.2, 6.0))
    text = found[0].describe()
    assert "(0 0 1)" in text
    assert "1.047" in text          # the period
    assert "6.000" in text          # the real-space repeat
    assert "%" in text              # how far off
