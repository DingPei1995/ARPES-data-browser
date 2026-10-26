"""The map conversions after the speed / memory round.

* ``convert_map`` samples each energy plane with ``map_coordinates`` at
  fractional indices instead of building a ``RegularGridInterpolator`` per
  slice. Same arithmetic, so the result must match the interpolator to
  rounding, NaN for NaN -- on uneven and on descending angle axes too.
* Neither ``convert_map`` nor ``to_kz_cube`` copies the whole cube to
  float64 any more, and a float32 measurement comes back as float32.
* ``_resample_energy`` is vectorised; it must equal one ``np.interp`` per
  spectrum.
* The MDC/EDC fitter's linear solve skips the bounded solver when the
  plain least-squares heights are already non-negative.
* A cut shifted to E - E_F is refused with a reason that names the energy
  offset, not the angles.
"""
import numpy as np
import pytest
from scipy import optimize
from scipy.interpolate import RegularGridInterpolator

from tools.momentum import kspace, kzconv
from tools.spectra import peaks
from tools.momentum.cutk import convert_cut


def _reference_convert_map(theta, phi, energy, cube, **kw):
    """The previous implementation, one interpolator per energy slice."""
    R = kspace.origin_rotation(kw.get("theta_offset_deg", 0.0),
                               kw.get("phi_offset_deg", 0.0),
                               kw.get("azimuth_deg", 0.0))
    kx_axis, ky_axis, energy_out, _ = kspace.convert_map(theta, phi, energy,
                                                         cube, **kw)
    t_order, p_order = np.argsort(theta), np.argsort(phi)
    KX, KY = np.meshgrid(kx_axis, ky_axis, indexing="ij")
    k0 = kspace.k0_of_energy(energy_out)
    out = np.full(KX.shape + (energy_out.size,), np.nan)
    for m in range(energy_out.size):
        kxu, kyu = KX / k0[m], KY / k0[m]
        radial = kxu ** 2 + kyu ** 2
        inside = radial <= 1.0
        kzu = np.sqrt(np.clip(1.0 - radial, 0.0, None))
        v = R.T @ np.stack([kxu.ravel(), kyu.ravel(), kzu.ravel()])
        phi_q = np.degrees(np.arcsin(np.clip(v[1], -1, 1))).reshape(KX.shape)
        theta_q = np.degrees(np.arctan2(v[0], v[2])).reshape(KX.shape)
        interp = RegularGridInterpolator(
            (theta[t_order], phi[p_order]),
            np.asarray(cube, float)[np.ix_(t_order, p_order, [m])][:, :, 0],
            bounds_error=False, fill_value=np.nan)
        plane = np.full(KX.shape, np.nan)
        plane[inside] = interp(np.stack([theta_q[inside], phi_q[inside]], -1))
        out[:, :, m] = plane
    return out


@pytest.mark.parametrize("descending", [False, True])
def test_convert_map_matches_the_interpolator(descending):
    rng = np.random.default_rng(5)
    theta = np.sort(rng.uniform(-14, 14, 31))          # uneven steps
    phi = np.linspace(-15, 15, 47)
    energy = np.linspace(50.0, 51.0, 9)
    cube = rng.random((theta.size, phi.size, energy.size))
    cube[4, 6, 2] = np.nan
    if descending:
        theta, cube = theta[::-1], cube[::-1]
    # A grid at least as fine as the data in both directions: coarser
    # grids average first (tested below), which no plain interpolator does.
    kw = dict(theta_offset_deg=1.5, phi_offset_deg=-2.0, azimuth_deg=12.0,
              n_kx=60, n_ky=70)
    got = kspace.convert_map(theta, phi, energy, cube, **kw)[3]
    want = _reference_convert_map(theta, phi, energy, cube, **kw)
    assert np.array_equal(np.isnan(got), np.isnan(want))
    np.testing.assert_allclose(got[np.isfinite(got)], want[np.isfinite(want)],
                               rtol=0, atol=1e-10)


def test_convert_map_keeps_float32_and_does_not_upcast_the_input():
    cube = np.random.default_rng(1).random((11, 13, 5)).astype(np.float32)
    out = kspace.convert_map(np.linspace(-5, 5, 11), np.linspace(-6, 6, 13),
                             np.linspace(30, 31, 5), cube, n_kx=8, n_ky=8)[3]
    assert out.dtype == np.float32
    out64 = kspace.convert_map(np.linspace(-5, 5, 11), np.linspace(-6, 6, 13),
                               np.linspace(30, 31, 5), cube.astype(float),
                               n_kx=8, n_ky=8)[3]
    assert out64.dtype == np.float64
    np.testing.assert_allclose(out[np.isfinite(out)], out64[np.isfinite(out64)],
                               rtol=1e-6)


def test_resample_energy_equals_interp_per_spectrum():
    rng = np.random.default_rng(2)
    energy = np.sort(rng.uniform(10, 12, 23))[::-1]     # uneven, descending
    cube = rng.random((4, 5, energy.size))
    out, target = kspace._resample_energy(cube, energy, 37)
    order = np.argsort(energy)
    for i in range(4):
        for j in range(5):
            np.testing.assert_allclose(
                out[i, j], np.interp(target, energy[order], cube[i, j, order]),
                rtol=0, atol=1e-12)


def test_kz_cube_keeps_float32():
    hv = np.linspace(30, 60, 7)
    angle = np.linspace(-10, 10, 9)
    energy = np.linspace(-0.5, 0.0, 4)
    cube = np.ones((hv.size, angle.size, energy.size), dtype=np.float32)
    out = kzconv.to_kz_cube(hv, angle, energy, cube, inner_potential=12.0,
                            work_function=4.4, n_kz=6, n_kpar=6)[3]
    assert out.dtype == np.float32
    assert np.nanmax(np.abs(out - 1.0)) < 1e-6


def test_linear_solve_fast_path_agrees_with_the_bounded_solver():
    rng = np.random.default_rng(3)
    design = np.column_stack([rng.random(50), rng.random(50), np.ones(50)])
    sigma = np.full(50, 0.5)
    for heights in ([2.0, 3.0], [2.0, -1.0]):          # bound inactive / active
        y = design @ np.array(heights + [0.4]) + rng.normal(0, 0.01, 50)
        got = peaks._solve_linear(design, y, sigma, 2)
        lower = np.array([0.0, 0.0, -np.inf])
        want = optimize.lsq_linear(design / sigma[:, None], y / sigma,
                                   bounds=(lower, np.full(3, np.inf)),
                                   method="bvls").x
        np.testing.assert_allclose(got, want, atol=1e-8)
        assert np.all(got[:2] >= 0)


def test_a_cut_on_e_minus_ef_is_refused_with_the_offset_named():
    slit = np.linspace(-10, 10, 21)
    energy = np.linspace(-0.5, 0.1, 13)                # E - E_F, not kinetic
    with pytest.raises(ValueError, match="energy offset"):
        convert_cut(slit, energy, np.ones((21, 13)), deflector_deg=0.0)
    # ...and adding E_F back is all it takes
    k, e, values, _ = convert_cut(slit, energy, np.ones((21, 13)),
                                  deflector_deg=0.0, energy_offset_eV=55.6)
    assert np.isfinite(values).any()


def _small_map(scale=1.0):
    rng = np.random.default_rng(4)
    cube = (rng.random((11, 13, 5)) * scale).astype(np.float64)
    return (np.linspace(-5, 5, 11), np.linspace(-6, 6, 13),
            np.linspace(30, 31, 5), cube)


def test_float16_output_matches_float32_to_its_own_rounding():
    theta, phi, energy, cube = _small_map(4000.0)
    half = kspace.convert_map(theta, phi, energy, cube, n_kx=8, n_ky=8,
                              precision="float16")[3]
    single = kspace.convert_map(theta, phi, energy, cube, n_kx=8, n_ky=8,
                                precision="float32")[3]
    assert half.dtype == np.float16 and single.dtype == np.float32
    assert np.array_equal(np.isnan(half), np.isnan(single))
    ok = np.isfinite(single)
    # one float16 rounding: half a step of 2^-10 relative
    np.testing.assert_allclose(half[ok].astype(float), single[ok], rtol=2.0 ** -11)


@pytest.mark.parametrize("scale", [1e6, 1e-4])
def test_float16_falls_back_when_the_values_do_not_fit(scale):
    # 1e6 overflows float16; 1e-4 would sink its weak tail below float16's
    # normal range. Both come back as float32 rather than damaged.
    theta, phi, energy, cube = _small_map(scale)
    out = kspace.convert_map(theta, phi, energy, cube, n_kx=8, n_ky=8,
                             precision="float16")[3]
    assert out.dtype == np.float32
    assert np.isfinite(out).any() and np.nanmax(out) > 0


def test_float16_holds_the_range_it_claims():
    assert kspace.float16_holds(0.0, 4091.0)
    assert kspace.float16_holds(-1.0, kspace.FLOAT16_MAX)
    assert not kspace.float16_holds(0.0, 7e4)
    assert not kspace.float16_holds(0.0, 1e-3)
    assert kspace.float16_holds(np.nan, np.nan)


def test_no_precision_keeps_following_the_measurement():
    theta, phi, energy, cube = _small_map()
    assert kspace.convert_map(theta, phi, energy, cube, n_kx=4, n_ky=4)[3].dtype == np.float64
    with pytest.raises(ValueError, match="precision"):
        kspace.convert_map(theta, phi, energy, cube, n_kx=4, n_ky=4, precision="int8")


# -- reducing conversions average instead of picking samples ---------------
from tools.momentum import rebin


def test_box_kernel_covers_its_width():
    np.testing.assert_allclose(rebin.box_kernel(3.0), [1 / 3] * 3)
    # 2.5 wide: the centre sample whole, three quarters of each neighbour
    np.testing.assert_allclose(rebin.box_kernel(2.5), np.array([0.75, 1, 0.75]) / 2.5)
    np.testing.assert_allclose(rebin.box_kernel(4.0), np.array([0.5, 1, 1, 1, 0.5]) / 4)
    assert rebin.box_kernel(1.0).tolist() == [1.0]


def test_bin_average_conserves_the_mean_and_averages_whole_bins():
    energy = np.arange(12, dtype=float)                   # step 1
    values = np.arange(12, dtype=float)[None, :] * np.ones((2, 1))
    target = np.linspace(0.5, 10.5, 6)                    # step 2: pairs
    out = rebin.bin_average(values, energy, target, axis=1)
    np.testing.assert_allclose(out[0], [0.5, 2.5, 4.5, 6.5, 8.5, 10.5])
    flat = rebin.bin_average(np.full((3, 50), 7.0), np.linspace(1, 2, 50),
                             np.linspace(1, 2, 17), axis=1)
    np.testing.assert_allclose(flat, 7.0)


def test_bin_average_leaves_nan_out():
    values = np.array([[1.0, np.nan, 3.0, 5.0]])
    out = rebin.bin_average(values, np.arange(4.0), np.array([0.5, 2.5]), axis=1)
    np.testing.assert_allclose(out[0], [1.0, 4.0])


def test_fewer_energy_points_are_bin_averages_not_samples():
    rng = np.random.default_rng(7)
    energy = np.linspace(40.0, 41.0, 301)
    cube = rng.poisson(100.0, (3, 4, energy.size)).astype(float)
    out, target = kspace._resample_energy(cube, energy, 101)    # 3x fewer
    # three channels per point: the noise drops by about sqrt(3)
    assert out.std() < 0.7 * cube.std()
    assert abs(out.mean() - cube.mean()) < 0.01 * cube.mean()
    # within the threshold of the measured step it is still plain interpolation
    same, target = kspace._resample_energy(cube, energy, 290)
    np.testing.assert_allclose(same[1, 2], np.interp(target, energy, cube[1, 2]))


def test_a_coarse_k_grid_averages_and_the_measured_grid_does_not():
    rng = np.random.default_rng(8)
    theta, phi = np.linspace(-10, 10, 101), np.linspace(-12, 12, 121)
    energy = np.linspace(60.0, 60.1, 2)
    cube = rng.poisson(100.0, (theta.size, phi.size, energy.size)).astype(float)
    coarse = kspace.convert_map(theta, phi, energy, cube, n_kx=34, n_ky=40)[3]
    fine = kspace.convert_map(theta, phi, energy, cube, n_kx=101, n_ky=121)[3]
    ok = np.isfinite(coarse)
    # three steps per pixel each way: about nine samples per point
    assert np.nanstd(coarse) < 0.5 * np.nanstd(fine)
    assert abs(np.nanmean(coarse) - 100.0) < 1.0
    # the measured sampling is not touched: same as with averaging disabled
    original = rebin.box_filter
    try:
        rebin.box_filter = lambda plane, widths: plane
        plain = kspace.convert_map(theta, phi, energy, cube, n_kx=101, n_ky=121)[3]
    finally:
        rebin.box_filter = original
    np.testing.assert_array_equal(np.isnan(plain), np.isnan(fine))
    np.testing.assert_allclose(fine[np.isfinite(fine)], plain[np.isfinite(plain)])
    assert ok.any()


def test_a_cut_on_fewer_k_points_averages():
    rng = np.random.default_rng(9)
    slit = np.linspace(-15, 15, 601)
    energy = np.linspace(60.0, 60.2, 3)
    values = rng.poisson(100.0, (slit.size, energy.size)).astype(float)
    coarse = convert_cut(slit, energy, values, deflector_deg=0.0, n_k=150)[2]
    fine = convert_cut(slit, energy, values, deflector_deg=0.0, n_k=601)[2]
    assert np.nanstd(coarse) < 0.6 * np.nanstd(fine)
    assert abs(np.nanmean(coarse) - 100.0) < 1.0


def test_a_small_rotation_does_not_average_the_fine_slit_direction():
    # A deflector Map: 0.2 deg deflector steps against 0.045 deg slit steps.
    # Rotated by 2 deg, a fifth of the coarse step projects onto the slit
    # axis; that must not count as a reduction at the measured sampling.
    rng = np.random.default_rng(10)
    theta, phi = np.linspace(-4, 4, 41), np.linspace(-9, 9, 401)
    energy = np.linspace(50.0, 50.1, 2)
    cube = rng.poisson(100.0, (theta.size, phi.size, energy.size)).astype(float)
    kw = dict(azimuth_deg=-2.0, theta_offset_deg=0.6, phi_offset_deg=-0.95,
              n_kx=theta.size, n_ky=phi.size)
    got = kspace.convert_map(theta, phi, energy, cube, **kw)[3]
    original = rebin.box_filter
    try:
        rebin.box_filter = lambda plane, widths: plane
        plain = kspace.convert_map(theta, phi, energy, cube, **kw)[3]
    finally:
        rebin.box_filter = original
    np.testing.assert_array_equal(got, plain)
