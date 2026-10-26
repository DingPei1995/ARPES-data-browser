"""Tests for tools/degrid.py on synthetic maps whose grid, contrast and
shifts are known: moving bands, a Fermi edge that does not move, a
hexagonal grid with a per-slice contrast and sub-pixel shift, Poisson."""
import numpy as np
import pytest

from tools import degrid as D

N0, N1, NS = 150, 190, 24          # angle px, energy px, slices
PERIOD = 5.2


def hex_grid(shape, shift=(0.0, 0.0), contrast=0.03):
    i0, i1 = np.meshgrid(np.arange(shape[0]) - shift[0],
                         np.arange(shape[1]) - shift[1], indexing="ij")
    g = np.zeros(shape)
    for angle in (10.0, 70.0, 130.0):
        t = np.radians(angle)
        g += np.cos(2 * np.pi / PERIOD * (np.cos(t) * i0 + np.sin(t) * i1))
    return 1.0 + contrast * g / 1.5


def photoemission(slice_index, level=400.0):
    k = np.linspace(-1, 1, N0)[:, None] + 0.04 * slice_index
    e = np.linspace(0, 1, N1)[None, :]
    band = np.exp(-((e - 0.35 - 0.3 * k ** 2) / 0.03) ** 2)
    band2 = np.exp(-((k - 0.2 * np.sin(6 * e)) / 0.05) ** 2) * (e < 0.8)
    fermi = 1 / (1 + np.exp((e - 0.8) / 0.01))          # never moves
    return level * (0.3 + band + 0.6 * band2) * fermi + 0.02 * level


BETAS = 1.0 + 0.15 * np.sin(np.arange(NS) / 3.0)
SHIFTS = np.column_stack([0.25 * np.sin(np.arange(NS) / 4.0),
                          -0.2 * np.cos(np.arange(NS) / 5.0)])


@pytest.fixture(scope="module")
def synthetic():
    rng = np.random.default_rng(0)
    truth = np.array([photoemission(s) for s in range(NS)])
    grids = np.array([1 + BETAS[s] * (hex_grid((N0, N1), SHIFTS[s]) - 1)
                      for s in range(NS)])
    measured = rng.poisson(truth * grids).astype(float)
    return truth, grids, measured


@pytest.fixture(scope="module")
def degridded(synthetic):
    return D.degrid_map(synthetic[2])


def test_the_grid_is_found(degridded):
    periods = [1 / np.hypot(a, b) for _s, a, b in degridded.model.peaks[:3]]
    np.testing.assert_allclose(periods, PERIOD, rtol=0.03)
    assert degridded.model.region.mean() < 0.05


def test_the_grid_is_removed(synthetic, degridded):
    truth, grids, measured = synthetic
    assert degridded.contrast_before > 3
    assert 0.8 < degridded.contrast_after < 1.3
    box = degridded.model.box
    # what is left, compared with the grid that was there and with counting
    # noise alone
    residual = (degridded.values[:, box[0], box[1]] / truth[:, box[0], box[1]] - 1)
    noise_only = (measured / grids)[:, box[0], box[1]] / truth[:, box[0], box[1]] - 1
    assert np.std(residual) < 1.05 * np.std(noise_only)


def test_contrast_and_shift_per_slice_recovered(degridded):
    betas = np.array([f.beta for f in degridded.fits])
    shifts = np.array([f.shift for f in degridded.fits])
    # the map grid is the average of the slices', so β and d come back
    # relative to that average
    assert np.corrcoef(betas, BETAS)[0, 1] > 0.9
    for axis in (0, 1):
        rel = SHIFTS[:, axis] - SHIFTS[:, axis].mean()
        got = shifts[:, axis] - shifts[:, axis].mean()
        assert np.corrcoef(rel, got)[0, 1] > 0.9
        assert np.median(np.abs(rel - got)) < 0.06


def test_what_does_not_move_is_left_alone(synthetic, degridded):
    """The Fermi edge is at the same pixel in every slice; averaging would
    take it for grid if the energy-only part were not removed."""
    truth, _grids, _measured = synthetic
    box = degridded.model.box
    edge_truth = truth[:, box[0], box[1]].sum(axis=(0, 1))
    edge_out = degridded.values[:, box[0], box[1]].sum(axis=(0, 1))
    np.testing.assert_allclose(edge_out / edge_truth, 1.0, atol=0.01)


def test_a_cut_with_the_map_grid(synthetic, degridded):
    rng = np.random.default_rng(5)
    cut_truth = photoemission(7.5)
    cut = rng.poisson(cut_truth * (1 + 1.1 * (hex_grid((N0, N1), (0.1, -0.1)) - 1))).astype(float)
    result = D.degrid_cut_with_grid(cut, degridded.model.full())
    assert result.contrast_before > 2 and 0.8 < result.contrast_after < 1.3
    with pytest.raises(ValueError, match="detector window"):
        D.degrid_cut_with_grid(cut[:-5], degridded.model.full())


def test_tiles_scale_with_the_image(synthetic, degridded):
    """150 × 190 px is too small for 5 × 5 tiles: one global fit is used."""
    model = degridded.model
    image = synthetic[2][3][model.box]
    sm = D._smooth(image, 8.0)
    G5, _ = D.fit_to(model, image, sm, tiles=5)
    G0, _ = D.fit_to(model, image, sm, tiles=0)
    np.testing.assert_allclose(G5, G0)


def test_the_notch_fallback():
    rng = np.random.default_rng(6)
    cut = rng.poisson(photoemission(3, level=3000) * hex_grid((N0, N1))).astype(float)
    result = D.degrid_cut_notch(cut)
    assert result.contrast_after < 1.0 < result.contrast_before
    assert any("photoemission" in note for note in result.notes)


def test_no_grid_no_invention():
    rng = np.random.default_rng(7)
    cube = np.array([rng.poisson(photoemission(s)) for s in range(10)]).astype(float)
    with pytest.raises(D.GridNotFound):
        D.degrid_map(cube)
    with pytest.raises(D.GridNotFound):
        D.degrid_cut_notch(cube[0])


def test_stored_pattern_round_trip(degridded):
    model = D.model_from_pattern(degridded.model.full())
    assert model.box == degridded.model.box
    np.testing.assert_allclose(model.pattern(), degridded.model.pattern(), atol=1e-9)
    assert model.region.mean() < 0.1


def test_illuminated_box():
    total = np.zeros((100, 120))
    total[20:80, 10:110] = 1.0
    box = D.illuminated_box(total, margin=5)
    assert box == (slice(25, 75), slice(15, 105))
    assert D.illuminated_box(np.zeros((40, 40))) == (slice(0, 40), slice(0, 40))


@pytest.mark.parametrize("info,blocked", [
    ({}, False),
    ({"proc.step.1": "truncate(x=1)"}, False),
    ({"fscorr.order": 2}, True),
    ({"kconv.work_function": 4.5}, True),
    ({"kz_align.window": 1}, True),
    ({"arbcut.joints": 1}, True),
    ({"degrid.method": "x"}, True),
    ({"proc.step.1": "curvature(a0=1)"}, True),
    ({"cassiopee.energy_reference": "per_member"}, True),
    ({"cassiopee.energy_reference": "common"}, False),
])
def test_pixel_lock_check(info, blocked):
    assert (D.not_pixel_locked(info) is not None) is blocked


def test_settings_match_and_default_source():
    from ui.degrid import default_source, settings_match
    a = {"lens_mode": "A30", "pass_energy_eV": 20.0}
    assert settings_match(a, dict(a))
    assert not settings_match(a, {"lens_mode": "A30", "pass_energy_eV": 50.0})
    assert not settings_match({}, {})             # nothing to compare: not chosen
    assert settings_match({"MBS.lens_mode": "L4", "MBS.passenergy": "PE100"},
                          {"MBS.lens_mode": "L4", "MBS.passenergy": "PE100"})
    assert default_source([("g1", None, None, False), ("g2", None, None, True)]) == 2
    assert default_source([("g1", None, None, False)]) == 0
