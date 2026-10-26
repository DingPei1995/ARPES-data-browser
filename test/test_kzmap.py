"""Putting a photon-energy scan onto one Fermi level.

These tests build a kz map whose Fermi level is *known*, spectrum by
spectrum, and then check that the calibration finds it back. That is the
only kind of test worth writing here: the algorithm's whole job is to
recover a number nothing in the data states outright, and the failure mode
is not a crash but a plausible-looking cube aligned to the wrong place.

The synthetic scan carries the two things that make the real one awkward --
a Fermi level that drifts with photon energy, and a beamline flux that
varies by nearly an order of magnitude across it.
"""
import numpy as np
import pytest

from tools import kzmap
from tools.fermi import fermi_edge_model, fit_fermi_edge


TEMPERATURE = 25.0


def make_scan(ef_per_spectrum, *, n_angle=24, n_energy=180, flux=None,
              resolution=0.025, noise=True, seed=0):
    """A kz map with a known Fermi level in every spectrum.

    ``flux`` scales each spectrum, standing in for the photon flux and
    analyser transmission, which between them vary by a large factor across
    a real scan and are the reason normalising is offered at all.
    """
    ef_per_spectrum = np.asarray(ef_per_spectrum, dtype=float)
    count = ef_per_spectrum.size
    energy = np.linspace(-1.2, 0.5, n_energy)
    angles = np.linspace(-15.0, 15.0, n_angle)
    if flux is None:
        flux = np.ones(count)
    band = 1.0 + 0.4 * np.cos(np.linspace(-np.pi, np.pi, n_angle))

    cube = np.empty((count, n_angle, n_energy))
    for index, ef in enumerate(ef_per_spectrum):
        edge = fermi_edge_model(energy, ef=float(ef), temperature=TEMPERATURE,
                                resolution=resolution, dos0=120.0, dos1=0.0,
                                bkg0=3.0, bkg1=0.0)
        cube[index] = float(flux[index]) * np.outer(band, edge)
    if noise:
        cube = np.random.default_rng(seed).poisson(np.clip(cube, 0, None)).astype(float)
    return cube, angles, energy


def whole(cube):
    """The box that covers everything, as inclusive index bounds."""
    return ((0, cube.shape[1] - 1), (0, cube.shape[2] - 1))


# --------------------------------------------------------------------------
# Fitting
# --------------------------------------------------------------------------
def test_the_fitted_levels_are_the_ones_that_were_put_in():
    true_ef = np.linspace(-0.3, 0.2, 9)
    cube, _angles, energy = make_scan(true_ef)
    ef, ok, _fits = kzmap.fit_levels(cube, energy, whole(cube),
                                     temperature=TEMPERATURE)
    assert ok.all()
    assert np.abs(ef - true_ef).max() < 0.01


def test_the_edc_sums_only_the_picked_channels():
    """The angle range is what separates a clean edge from a noisy one on a
    detector whose ends see nothing."""
    cube, _angles, _energy = make_scan([0.0], n_angle=10)
    edc = kzmap.region_edc(cube[0], ((2, 5), (0, 4)))
    assert np.allclose(edc, cube[0][2:6, :].sum(axis=0))


def test_the_energy_range_narrows_the_fit_not_the_data():
    """Picking a narrow energy window must change where the fit looks, not
    what the cube contains."""
    true_ef = np.linspace(-0.1, 0.1, 5)
    cube, _angles, energy = make_scan(true_ef)
    near = int(np.argmin(np.abs(energy + 0.4))), int(np.argmin(np.abs(energy - 0.4)))
    ef, ok, _ = kzmap.fit_levels(cube, energy, ((0, cube.shape[1] - 1), near),
                                 temperature=TEMPERATURE)
    assert ok.all()
    assert np.abs(ef - true_ef).max() < 0.01


def test_a_spectrum_that_cannot_be_fitted_is_interpolated_not_dropped():
    """One dead spectrum in the middle of a scan should not take the cube
    with it, and the fill-in has to be flagged rather than silent."""
    true_ef = np.linspace(-0.2, 0.2, 7)
    cube, _angles, energy = make_scan(true_ef)
    cube[3] = 0.0                       # nothing to fit: a dropped acquisition
    ef, ok, _fits = kzmap.fit_levels(cube, energy, whole(cube),
                                     temperature=TEMPERATURE)
    assert not ok[3]
    assert ok.sum() == 6
    assert np.isfinite(ef[3])
    # interpolated from its neighbours, so it lands between them
    assert min(ef[2], ef[4]) - 1e-9 <= ef[3] <= max(ef[2], ef[4]) + 1e-9


def test_a_scan_with_no_edge_anywhere_is_an_error():
    """Better than returning a cube aligned to noise.

    Flat data does not make the fitter raise -- it converges happily and
    slides the edge off the end of the window, reporting success. Catching
    that is why fit_levels checks the answer lands inside the window it
    fitted over.
    """
    cube = np.ones((4, 10, 60))
    energy = np.linspace(-1.0, 1.0, 60)
    with pytest.raises(ValueError, match="could not be fitted in any"):
        kzmap.fit_levels(cube, energy, whole(cube), temperature=TEMPERATURE)


def test_a_level_outside_the_fitted_window_is_not_believed():
    """The runaway case, seen on real data: a wide window lets the fit reach
    the band structure below the edge and settle somewhere off the end. The
    number still looks like an answer, and shifting a spectrum by it wrecks
    the stack."""
    true_ef = np.array([0.0, 0.0, 0.0])
    cube, _angles, energy = make_scan(true_ef, noise=False)
    cube[1] = 1.0                       # featureless: nothing to lock onto
    ef, ok, _fits = kzmap.fit_levels(cube, energy, whole(cube),
                                     temperature=TEMPERATURE)
    assert ok[0] and ok[2]
    assert not ok[1]
    lo, hi = float(energy.min()), float(energy.max())
    assert lo <= ef[1] <= hi


# --------------------------------------------------------------------------
# Aligning and cropping
# --------------------------------------------------------------------------
def test_aligning_puts_every_edge_at_zero():
    true_ef = np.array([-0.25, -0.1, 0.05, 0.2])
    cube, _angles, energy = make_scan(true_ef, noise=False)
    aligned, axis, _trimmed = kzmap.align(cube, energy, true_ef)
    for index in range(cube.shape[0]):
        fit = fit_fermi_edge(axis, aligned[index].sum(axis=0),
                             temperature=TEMPERATURE)
        assert abs(fit.values["ef"]) < 0.005


def test_the_crop_is_exactly_the_spread_of_the_levels():
    """Shifting each spectrum by its own amount leaves the stack ragged;
    what survives is the original width less the spread."""
    ef = np.array([-0.3, 0.0, 0.2])
    cube, _angles, energy = make_scan(ef, noise=False)
    _aligned, axis, trimmed = kzmap.align(cube, energy, ef)
    step = abs(energy[1] - energy[0])
    assert trimmed == pytest.approx(ef.max() - ef.min(), abs=2 * step)
    kept = axis[-1] - axis[0]
    assert kept == pytest.approx((energy[-1] - energy[0]) - trimmed, abs=2 * step)


def test_the_kept_range_is_covered_by_every_spectrum():
    """The point of cropping: no spectrum may be extrapolated into."""
    ef = np.array([-0.3, 0.0, 0.25])
    cube, _angles, energy = make_scan(ef, noise=False)
    _aligned, axis, _trimmed = kzmap.align(cube, energy, ef)
    for value in ef:
        assert axis[0] >= energy[0] - value - 1e-9
        assert axis[-1] <= energy[-1] - value + 1e-9


def test_nothing_is_blank_after_the_crop():
    ef = np.linspace(-0.3, 0.3, 6)
    cube, _angles, energy = make_scan(ef)
    aligned, _axis, _trimmed = kzmap.align(cube, energy, ef)
    assert np.isfinite(aligned).all()


def test_a_zero_shift_leaves_the_data_alone():
    cube, _angles, energy = make_scan([0.0, 0.0, 0.0], noise=False)
    aligned, axis, trimmed = kzmap.align(cube, energy, np.zeros(3))
    assert trimmed == pytest.approx(0.0, abs=1e-9)
    assert np.allclose(axis, energy)
    assert np.allclose(aligned, cube)


def test_a_descending_energy_axis_works_too():
    """A scale written high-to-low is legal and does happen."""
    ef = np.array([-0.2, 0.0, 0.15])
    cube, _angles, energy = make_scan(ef, noise=False)
    flipped = cube[:, :, ::-1]
    aligned, axis, _trimmed = kzmap.align(flipped, energy[::-1], ef)
    for index in range(3):
        fit = fit_fermi_edge(axis, aligned[index].sum(axis=0),
                             temperature=TEMPERATURE)
        assert abs(fit.values["ef"]) < 0.01


def test_levels_spread_wider_than_the_window_is_an_error():
    """Nothing would be left in common, and a cube of nothing is not an
    improvement on a misaligned one."""
    cube, _angles, energy = make_scan([0.0, 0.0], noise=False)
    with pytest.raises(ValueError, match="no energy range is left"):
        kzmap.align(cube, energy, np.array([-5.0, 5.0]))


def test_one_level_per_spectrum_is_required():
    cube, _angles, energy = make_scan([0.0, 0.1], noise=False)
    with pytest.raises(ValueError, match="Fermi levels for"):
        kzmap.align(cube, energy, np.array([0.0]))


# --------------------------------------------------------------------------
# Normalising
# --------------------------------------------------------------------------
def test_normalising_makes_every_spectrum_sum_to_one():
    cube, _angles, _energy = make_scan([0.0] * 5, flux=[1, 3, 9, 2, 5])
    out = kzmap.normalise_totals(cube)
    assert np.allclose(out.sum(axis=(1, 2)), 1.0)


def test_normalising_removes_the_flux_variation_and_keeps_the_shape():
    """A tenfold flux difference is the beamline, not the sample. What is
    left afterwards has to be the same spectrum."""
    cube, _angles, _energy = make_scan([0.0, 0.0], flux=[1.0, 10.0], noise=False)
    out = kzmap.normalise_totals(cube)
    assert np.allclose(out[0], out[1], rtol=1e-9)


def test_an_empty_spectrum_is_left_alone_rather_than_made_nan():
    cube, _angles, _energy = make_scan([0.0, 0.0], noise=False)
    cube[1] = 0.0
    out = kzmap.normalise_totals(cube)
    assert np.isfinite(out).all()
    assert np.allclose(out[1], 0.0)


# --------------------------------------------------------------------------
# All of it
# --------------------------------------------------------------------------
def test_the_whole_thing_recovers_a_drifting_fermi_level():
    true_ef = 0.25 * np.sin(np.linspace(0, 2.0, 11)) - 0.05
    flux = np.linspace(1.0, 9.0, 11)
    cube, _angles, energy = make_scan(true_ef, flux=flux)
    result = kzmap.process_kz_map(cube, energy, whole(cube),
                                  temperature=TEMPERATURE, normalise=True)
    assert result.ok.all()
    assert np.abs(result.ef - true_ef).max() < 0.01
    assert result.spread == pytest.approx(true_ef.max() - true_ef.min(), abs=0.01)
    assert np.allclose(result.cube.sum(axis=(1, 2)), 1.0)
    assert np.isfinite(result.cube).all()
    for index in range(cube.shape[0]):
        fit = fit_fermi_edge(result.energy, result.cube[index].sum(axis=0),
                             temperature=TEMPERATURE)
        assert abs(fit.values["ef"]) < 0.01


def test_normalising_can_be_turned_off():
    flux = [1.0, 5.0, 10.0]
    cube, _angles, energy = make_scan([-0.1, 0.0, 0.1], flux=flux)
    result = kzmap.process_kz_map(cube, energy, whole(cube),
                                  temperature=TEMPERATURE, normalise=False)
    totals = result.cube.sum(axis=(1, 2))
    assert not np.allclose(totals, 1.0)
    assert totals[2] / totals[0] == pytest.approx(10.0, rel=0.1)
    assert result.normalised is False


def test_normalising_happens_after_the_crop():
    """Each total has to be a sum over the same energy range, or the
    normalisation carries the misalignment it is meant to be blind to.

    Checked by the order's own signature: the spectra are shifted by
    different amounts, so the crop removes a different fraction of each
    one's original intensity. Totals of exactly 1 afterwards can only
    happen if the sums were taken after the crop -- normalising first would
    leave each spectrum at 1 minus its own cropped fraction, and those
    differ.
    """
    cube, _angles, energy = make_scan([-0.2, 0.0, 0.2], noise=False)
    result = kzmap.process_kz_map(cube, energy, whole(cube),
                                  temperature=TEMPERATURE, normalise=True)
    assert np.allclose(result.cube.sum(axis=(1, 2)), 1.0)

    # Show the fractions really do differ, so the check above has teeth.
    aligned, _axis, _trim = kzmap.align(cube, energy, result.ef)
    fractions = aligned.sum(axis=(1, 2)) / cube.sum(axis=(1, 2))
    assert fractions.max() - fractions.min() > 0.01


def test_the_summary_says_what_happened():
    cube, _angles, energy = make_scan(np.linspace(-0.1, 0.1, 4))
    result = kzmap.process_kz_map(cube, energy, whole(cube),
                                  temperature=TEMPERATURE)
    text = result.summary()
    assert "4 spectra" in text
    assert "spread" in text
    assert "trimmed" in text
    assert "total intensity" in text


def test_progress_is_reported_once_per_spectrum():
    cube, _angles, energy = make_scan([0.0] * 6)
    seen = []
    kzmap.fit_levels(cube, energy, whole(cube), temperature=TEMPERATURE,
                     progress=lambda done, total: seen.append((done, total)))
    assert [d for d, _ in seen] == list(range(6))
    assert {t for _, t in seen} == {6}
