"""Tests for tools/cutops.py -- arithmetic between two cuts."""
import numpy as np
import pytest

from tools import cutops as C


ANG = ("Angle (deg)", "Energy (eV)")


def grid(nx=41, ny=31):
    return [np.linspace(-10, 10, nx), np.linspace(-1.0, 0.2, ny)]


def band(axes, amplitude=100.0, offset=5.0):
    x, y = np.meshgrid(axes[0], axes[1], indexing="ij")
    return offset + amplitude * np.exp(-((y + 0.1 * x ** 2 / 10 + 0.3) / 0.05) ** 2)


# -- the operations ----------------------------------------------------------
@pytest.mark.parametrize("operation,expected", [
    ("difference", lambda a, b: a - b),
    ("sum", lambda a, b: a + b),
    ("ratio", None),
    ("asymmetry", lambda a, b: (a - b) / (a + b)),
])
def test_operations(operation, expected):
    axes = grid()
    a = band(axes, 120.0)
    b = band(axes, 80.0)
    r = C.combine(a, axes, b, axes, operation=operation, a_labels=ANG, b_labels=ANG)
    assert r.values.shape == a.shape
    assert not r.resampled and r.overlap == 1.0 and r.scale == 1.0
    if expected is not None:
        np.testing.assert_allclose(r.values, expected(a, b))
    else:
        np.testing.assert_allclose(r.values, a / (b / b.mean()))


def test_bad_operation_and_shape():
    axes = grid()
    a = band(axes)
    with pytest.raises(ValueError):
        C.combine(a, axes, a, axes, operation="product")
    with pytest.raises(ValueError):
        C.combine(a, axes, a, axes, operation="difference", normalise="max")
    with pytest.raises(ValueError):
        C.combine(a[:-1], axes, a, axes, operation="difference")


# -- normalisation -----------------------------------------------------------
def test_total_normalisation_recovers_flux_ratio():
    axes = grid()
    a = band(axes)
    r = C.combine(a, axes, 0.37 * a, axes, operation="difference",
                  normalise="total")
    assert r.scale == pytest.approx(1 / 0.37)
    np.testing.assert_allclose(r.values, 0.0, atol=1e-9)
    np.testing.assert_allclose(r.b_used, a)


def test_region_normalisation_uses_only_the_box():
    axes = grid()
    a = band(axes)
    b = a.copy()
    # B is doubled only below -0.5 eV: a box above that sees scale 1.
    b[:, axes[1] < -0.5] *= 2.0
    r = C.combine(a, axes, b, axes, operation="difference", normalise="region",
                  region=(-10, -0.3, 10, 0.2))
    assert r.scale == pytest.approx(1.0)
    with pytest.raises(ValueError):
        C.combine(a, axes, b, axes, operation="difference", normalise="region")
    with pytest.raises(ValueError):
        C.combine(a, axes, b, axes, operation="difference", normalise="region",
                  region=(50, 5, 60, 6))


def test_ratio_is_never_rescaled():
    axes = grid()
    a = band(axes)
    r = C.combine(a, axes, 10 * a, axes, operation="ratio", normalise="total")
    assert r.scale == 1.0


# -- grids -------------------------------------------------------------------
def test_resampling_and_partial_overlap():
    def plane(ax):       # bilinear, so bilinear interpolation is exact
        x, y = np.meshgrid(*ax, indexing="ij")
        return 50 + 2 * x + 30 * y + 0.5 * x * y

    axes = grid()
    a = plane(axes)
    b_axes = [np.linspace(0, 20, 61), np.linspace(-1.0, 0.2, 45)]
    b = plane(b_axes)
    r = C.combine(a, axes, b, b_axes, operation="difference", a_labels=ANG,
                  b_labels=ANG)
    assert r.resampled
    assert 0.45 < r.overlap < 0.55
    assert np.all(np.isnan(r.values[axes[0] < 0]))
    inside = axes[0] >= 0
    np.testing.assert_allclose(r.values[inside], 0.0, atol=1e-9)
    assert r.sigma is None        # interpolated data are not Poisson counts


def test_descending_axis_is_handled():
    axes = grid()
    a = band(axes)
    b_axes = [axes[0][::-1], axes[1]]
    b = a[::-1].copy()
    r = C.combine(a, axes, b, b_axes, operation="difference")
    np.testing.assert_allclose(r.values, 0.0, atol=1e-8)


def test_unit_mismatch_is_refused():
    axes = grid()
    a = band(axes)
    with pytest.raises(ValueError, match="deg"):
        C.combine(a, axes, a, axes, operation="difference", a_labels=ANG,
                  b_labels=("k (Å⁻¹)", "Energy (eV)"))


def test_no_overlap_is_refused():
    axes = grid()
    a = band(axes)
    far = [axes[0] + 100, axes[1]]
    with pytest.raises(ValueError, match="overlap"):
        C.combine(a, axes, a, far, operation="difference")


# -- dividing by a reference -------------------------------------------------
def test_reference_profiles_and_floor():
    axes = grid()
    x, y = np.meshgrid(*axes, indexing="ij")
    sensitivity = 1.0 + 0.5 * np.cos(np.radians(9 * x))     # angle-dependent
    truth = band(axes)
    a = truth * sensitivity
    gold = 1000.0 * sensitivity * (1 + 0.0 * y)
    r = C.combine(a, axes, gold, axes, operation="ratio",
                  reference_shape="angle_profile")
    corrected = r.values / np.nanmean(r.values)
    np.testing.assert_allclose(corrected, truth / truth.mean(), rtol=1e-9)

    # energy profile: B varying only in energy removes only that
    fermi = 1 / (np.exp(y / 0.02) + 1)
    r = C.combine(truth * fermi, axes, 50 * fermi, axes, operation="ratio",
                  reference_shape="energy_profile", reference_floor=0.05)
    assert r.masked > 0                           # above E_F the gold is ~0
    ok = np.isfinite(r.values)
    ratio = r.values[ok] / truth[ok]
    np.testing.assert_allclose(ratio, ratio[0], rtol=1e-9)


def test_empty_reference_is_refused():
    axes = grid()
    with pytest.raises(ValueError, match="nothing to divide"):
        C.combine(band(axes), axes, np.zeros((41, 31)), axes, operation="ratio")


# -- masks and statistics ----------------------------------------------------
def test_intensity_floor_masks_the_background():
    axes = grid()
    a = band(axes, offset=0.0)
    r = C.combine(a, axes, 0.5 * a, axes, operation="asymmetry",
                  intensity_floor=0.02)
    assert 0 < r.masked < 1
    assert np.all(np.isnan(r.values[(a + 0.5 * a) < 0.02 * (1.5 * a).max()]))
    np.testing.assert_allclose(r.values[np.isfinite(r.values)], 1 / 3)


@pytest.mark.parametrize("operation", ["difference", "asymmetry"])
def test_poisson_sigma_matches_monte_carlo(operation):
    rng = np.random.default_rng(1)
    axes = [np.arange(20.0), np.arange(15.0)]
    mean_a = np.full((20, 15), 400.0)
    mean_b = np.full((20, 15), 250.0)
    draws = []
    for _ in range(400):
        a = rng.poisson(mean_a).astype(float)
        b = rng.poisson(mean_b).astype(float)
        draws.append(C.combine(a, axes, b, axes, operation=operation).values)
    spread = np.std(np.array(draws), axis=0).mean()
    a = rng.poisson(mean_a).astype(float)
    b = rng.poisson(mean_b).astype(float)
    r = C.combine(a, axes, b, axes, operation=operation)
    assert r.sigma is not None
    assert r.median_sigma == pytest.approx(spread, rel=0.06)
    assert r.significant is not None and 0.9 < r.significant <= 1.0


def test_poisson_sigma_with_scaling():
    rng = np.random.default_rng(2)
    axes = [np.arange(30.0), np.arange(20.0)]
    mean_a, mean_b = 300.0, 150.0
    a = rng.poisson(mean_a, (30, 20)).astype(float)
    b = rng.poisson(mean_b, (30, 20)).astype(float)
    r = C.combine(a, axes, b, axes, operation="difference", normalise="total")
    assert r.scale == pytest.approx(2.0, rel=0.02)
    # var = A + s^2 B = 300 + 4 * 150 = 900
    assert r.median_sigma == pytest.approx(30.0, rel=0.03)


def test_low_counts_note():
    rng = np.random.default_rng(3)
    axes = [np.arange(30.0), np.arange(20.0)]
    a = rng.poisson(5.0, (30, 20)).astype(float)
    b = rng.poisson(4.0, (30, 20)).astype(float)
    r = C.combine(a, axes, b, axes, operation="asymmetry")
    assert any("underestimate" in n for n in r.notes)
    assert "underestimate" in r.summary()


def test_no_sigma_for_processed_data():
    axes = grid()
    a = band(axes) + 0.123
    r = C.combine(a, axes, a * 0.9, axes, operation="difference")
    assert r.sigma is None and r.median_sigma is None


def test_looks_like_counts():
    assert C.looks_like_counts(np.array([0, 1, 5, 7.0]))
    assert not C.looks_like_counts(np.array([0.5, 1.2]))
    assert not C.looks_like_counts(np.array([-1.0, 3.0]))
    assert not C.looks_like_counts(np.zeros(4))
    assert not C.looks_like_counts(np.array([np.nan]))


def test_axis_unit():
    assert C.axis_unit("Angle (deg)") == "deg"
    assert C.axis_unit("Angle (°)") == "deg"
    assert C.axis_unit("k (Å⁻¹)") == C.axis_unit("k (A^-1)") == C.axis_unit("k (1/A)")
    assert C.axis_unit("Energy (eV)") == "ev"
    assert C.axis_unit("") == ""
    assert C.axis_unit("no unit") == ""


# -- metadata ----------------------------------------------------------------
def info(pol, hv=80.0, theta=0.0, x=1.0):
    return {"cassiopee.polarisation": pol, "photon_energy_eV": hv,
            "cassiopee.sample_theta_deg": theta, "cassiopee.sample_X_mm": x,
            "lens_mode": "A30"}


def test_metadata_differences():
    rows = {r[0]: r for r in C.metadata_differences(info("LH"), info("LV", theta=0.2,
                                                                      x=1.002))}
    assert rows["Polarisation"][3] is True
    assert rows["Photon energy (eV)"][3] is False
    assert rows["Theta (deg)"][3] is True
    assert rows["X (mm)"][3] is False          # 2 um is within 5 um
    assert rows["Lens mode"][3] is False
    assert "Temperature (K)" not in rows       # neither recorded it
    one_sided = {r[0]: r for r in C.metadata_differences(info("LH"), {})}
    assert one_sided["Polarisation"][3] is None


def test_metadata_pattern_keys():
    rows = {r[0]: r for r in C.metadata_differences(
        {"beamline.polarization": "LH", "manipulator.theta": 1.0},
        {"beamline.polarization": "LV", "manipulator.theta": 1.0})}
    assert rows["Polarisation"][3] is True
    assert rows["Theta (deg)"][3] is False


def test_polarisation_warnings():
    assert C.polarisation_warnings(info("LH"), info("LV"), "linear_dichroism") == []
    swapped = C.polarisation_warnings(info("LV"), info("LH"), "linear_dichroism")
    assert swapped and "Swap" in swapped[0]
    same = C.polarisation_warnings(info("LH"), info("LH"), "linear_dichroism")
    assert same and "both" in same[0]
    wrong = C.polarisation_warnings(info("CR"), info("CL"), "linear_dichroism")
    assert wrong and "linear" in wrong[0]
    wrong = C.polarisation_warnings(info("LH"), info("LV"), "circular_dichroism")
    assert wrong and "circular" in wrong[0]
    assert C.polarisation_warnings(info("LCP"), info("RCP"), "circular_dichroism") == []
    assert C.polarisation_warnings(info("LH"), info("LH"), "reference") == []
    assert C.polarisation_warnings({}, info("LH"), "linear_dichroism") == []
