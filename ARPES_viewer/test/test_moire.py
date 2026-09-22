"""Tests for tools/moire.py: the general reciprocal-vector-difference moire
construction, cross-checked against ``moire_lattice.m``'s closed-form
hexagonal-bilayer formula (:func:`hex_moire_lattice_fast`).
"""
import numpy as np
import pytest

from tools.bz2d import reciprocal_vectors_2d
from tools.moire import (hex_moire_lattice_fast, moire_reciprocal_vectors,
                       moire_bz)


def _hex_vectors(a, rotation_deg):
    r = np.deg2rad(rotation_deg)
    rot = np.array([[np.cos(r), -np.sin(r)], [np.sin(r), np.cos(r)]])
    a1 = rot @ np.array([a, 0.0])
    a2 = rot @ np.array([a * np.cos(np.deg2rad(120.0)), a * np.sin(np.deg2rad(120.0))])
    return a1, a2


def _moire_period_general(a_top, a_bot, twist_deg):
    g1t, g2t = reciprocal_vectors_2d(*_hex_vectors(a_top, 0.0))
    g1b, g2b = reciprocal_vectors_2d(*_hex_vectors(a_bot, twist_deg))
    gm1, gm2 = moire_reciprocal_vectors(g1t, g2t, g1b, g2b)
    am1, am2 = reciprocal_vectors_2d(gm1, gm2)
    return np.linalg.norm(am1)


@pytest.mark.parametrize("twist", [1.0, 3.0, 5.0, 10.0, 20.0, 29.0])
def test_general_method_matches_closed_form_equal_hexagonal_layers(twist):
    a = 2.46
    general = _moire_period_general(a, a, twist)
    analytic, _ = hex_moire_lattice_fast(a, a, twist)
    assert general == pytest.approx(analytic, rel=1e-9)


@pytest.mark.parametrize("twist,folded", [
    (31.0, 29.0), (35.0, 25.0), (40.0, 20.0), (45.0, 15.0), (55.0, 5.0),
])
def test_general_method_folds_past_30_degrees_by_hexagonal_symmetry(twist, folded):
    """Twisting a hexagonal bilayer by theta gives the same moire pattern as
    60 - theta (a real consequence of the hexagonal point group); the
    closed-form formula does not know this and only matches up to 30
    degrees (see the other test), but the general method finds the smaller
    of the two automatically."""
    a = 2.46
    general_at_twist = _moire_period_general(a, a, twist)
    general_at_fold = _moire_period_general(a, a, folded)
    assert general_at_twist == pytest.approx(general_at_fold, rel=1e-9)


def test_zero_twist_pure_mismatch_matches_the_two_comb_beat_length():
    """For zero twist, two lattices with slightly different spacing beat
    against each other exactly like two combs: period = a1*a2/|a1-a2|. This
    is the standard textbook zero-twist limit, independent of hexagonal
    geometry (it only uses the two nearest-neighbour reciprocal vectors,
    same as a 1-D lattice would); the general method must reduce to it.
    """
    a_top, a_bot = 2.46, 2.50
    general = _moire_period_general(a_top, a_bot, 0.0)
    expected = a_top * a_bot / abs(a_top - a_bot)
    assert general == pytest.approx(expected, rel=1e-6)


def test_moire_bz_is_hexagonal_for_two_hexagonal_layers():
    a = 2.46
    g1t, g2t = reciprocal_vectors_2d(*_hex_vectors(a, 0.0))
    g1b, g2b = reciprocal_vectors_2d(*_hex_vectors(a, 8.0))
    gm1, gm2, poly = moire_bz(g1t, g2t, g1b, g2b)
    assert len(poly) - 1 == 6


def test_identical_layers_have_no_moire_pattern():
    a = 2.46
    g1, g2 = reciprocal_vectors_2d(*_hex_vectors(a, 0.0))
    with pytest.raises(ValueError, match="no moire pattern"):
        moire_reciprocal_vectors(g1, g2, g1, g2)
