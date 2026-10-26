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
    kw = dict(theta_offset_deg=1.5, phi_offset_deg=-2.0, azimuth_deg=12.0,
              n_kx=40, n_ky=44)
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
