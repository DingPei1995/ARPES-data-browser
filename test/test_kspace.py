"""The sample rotation read off a picked direction.

There was no test file here, which is why the bug these are written against
shipped: ``points_to_azimuth`` folded its answer to (-180, 180], so clicking
the same two points in the other order -- or picking a pair on the far side
of the origin -- returned a rotation 180 degrees away. That number looks
entirely reasonable and converts the map upside down.

Two points mark a *line*. A line has no direction, so nothing about how it
was picked may reach the answer, and the answer is the smaller of the two
turns that stand it vertical.
"""
import numpy as np
import pytest

from tools.kspace import points_to_azimuth


def direction(degrees, length=1.0):
    """Two points along a line at ``degrees`` to the kx axis."""
    radians = np.radians(degrees)
    return [(0.0, 0.0), (length * np.cos(radians), length * np.sin(radians))]


# --------------------------------------------------------------------------
# The property the feature is about
# --------------------------------------------------------------------------
@pytest.mark.parametrize("degrees", [0, 1, 17, 30, 45, 60, 89, 90, 91, 120,
                                     135, 150, 179, 200, 270, 315, 359])
def test_the_order_of_the_two_points_does_not_matter(degrees):
    first, second = direction(degrees)
    assert points_to_azimuth([first, second]) == pytest.approx(
        points_to_azimuth([second, first]), abs=1e-9)


@pytest.mark.parametrize("degrees", [0, 23, 45, 67, 90, 111, 135, 168])
def test_reversing_the_direction_does_not_matter(degrees):
    """Picking the pair on the other side of the origin traces the same
    line the other way, and must give the same rotation."""
    forward = points_to_azimuth(direction(degrees))
    backward = points_to_azimuth(direction(degrees + 180.0))
    assert forward == pytest.approx(backward, abs=1e-9)


@pytest.mark.parametrize("degrees", [5, 35, 95, 140, 200, 350])
def test_where_the_pair_sits_does_not_matter(degrees):
    """Only the direction counts, not the offset: the rotation is about the
    origin, so a pair translated anywhere gives the same answer."""
    (x0, y0), (x1, y1) = direction(degrees)
    moved = [(x0 + 3.7, y0 - 1.2), (x1 + 3.7, y1 - 1.2)]
    assert points_to_azimuth(direction(degrees)) == pytest.approx(
        points_to_azimuth(moved), abs=1e-9)


@pytest.mark.parametrize("degrees", [4, 33, 78, 96, 155, 233, 301])
def test_how_far_apart_they_are_does_not_matter(degrees):
    assert points_to_azimuth(direction(degrees, 0.05)) == pytest.approx(
        points_to_azimuth(direction(degrees, 40.0)), abs=1e-9)


# --------------------------------------------------------------------------
# The range
# --------------------------------------------------------------------------
def test_the_answer_is_always_the_smaller_turn():
    for degrees in np.linspace(0.0, 360.0, 1441):
        rotation = points_to_azimuth(direction(degrees))
        assert -90.0 < rotation <= 90.0, (degrees, rotation)


def test_a_vertical_line_needs_no_rotation():
    assert points_to_azimuth(direction(90)) == pytest.approx(0.0, abs=1e-9)
    assert points_to_azimuth(direction(270)) == pytest.approx(0.0, abs=1e-9)


def test_a_horizontal_line_always_gives_the_same_quarter_turn():
    """A quarter turn either way is equally small, so this is a real tie --
    and ``arctan2`` lands on opposite sides of it depending on the sign of a
    sine that should be zero. Every way of describing the same line has to
    give one number."""
    answers = {
        points_to_azimuth([(0.0, 0.0), (1.0, 0.0)]),
        points_to_azimuth([(0.0, 0.0), (-1.0, 0.0)]),
        points_to_azimuth([(1.0, 0.0), (0.0, 0.0)]),
        points_to_azimuth([(-2.0, 0.5), (3.0, 0.5)]),
        points_to_azimuth(direction(0.0)),
        points_to_azimuth(direction(180.0)),
        points_to_azimuth(direction(360.0)),
    }
    assert answers == {90.0}


def test_zero_is_written_without_a_minus_sign():
    """``-0.0`` formats as "-0.000", which reads like a real offset."""
    assert np.copysign(1.0, points_to_azimuth(direction(90))) > 0


# --------------------------------------------------------------------------
# What it actually computes
# --------------------------------------------------------------------------
@pytest.mark.parametrize("degrees,expected", [
    (90, 0.0), (60, 30.0), (45, 45.0), (30, 60.0), (0, 90.0),
    (120, -30.0), (135, -45.0), (150, -60.0), (179, -89.0),
])
def test_the_rotation_stands_the_line_upright(degrees, expected):
    assert points_to_azimuth(direction(degrees)) == pytest.approx(expected,
                                                                  abs=1e-9)


def test_one_point_is_the_line_from_the_origin():
    assert points_to_azimuth([(1.0, 1.0)]) == pytest.approx(45.0, abs=1e-9)
    # ... and the far side of the origin is the same line
    assert points_to_azimuth([(-1.0, -1.0)]) == pytest.approx(45.0, abs=1e-9)


def test_no_points_is_an_error():
    with pytest.raises(ValueError, match="at least one point"):
        points_to_azimuth([])


def test_two_identical_points_are_an_error():
    with pytest.raises(ValueError, match="define no direction"):
        points_to_azimuth([(0.3, -0.2), (0.3, -0.2)])


# --------------------------------------------------------------------------
# Does the rotation do what it claims, all the way through a conversion?
# --------------------------------------------------------------------------
@pytest.mark.parametrize("degrees", [20.0, 55.0, 130.0])
def test_the_rotation_really_does_stand_the_direction_up(degrees):
    """End to end against the rotation the conversion applies, rather than
    against the formula that produced it."""
    from tools.kspace import rotation_matrix

    rotation = points_to_azimuth(direction(degrees))
    radians = np.radians(degrees)
    vector = np.array([np.cos(radians), np.sin(radians), 0.0])
    turned = rotation_matrix(0.0, 0.0, rotation) @ vector
    # Up the ky axis, give or take which end -- the line has no direction.
    assert abs(turned[0]) < 1e-9
    assert abs(abs(turned[1]) - 1.0) < 1e-9
