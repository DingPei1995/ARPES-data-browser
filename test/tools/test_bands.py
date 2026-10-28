"""Calculated band structures (tools/crystal/bands.py, bandfiles.py) and
their comparison with a measurement (tools/crystal/compare.py): the files
read back what was written, the cut on the k_z sphere, the contours, the
Z and shift fit, the measurement assembled along a path."""
import os

import numpy as np
import pytest

from tools.crystal import bandfiles
from tools.crystal import bands as B
from tools.crystal import compare as C
from tools.crystal.lattice import LatticeParams

A, CC = 3.3, 6.0
LATTICE = np.array([[A, 0, 0], [-A / 2, A * np.sqrt(3) / 2, 0], [0, 0, CC]])


def _hoppings(scale=1.0, shift=0.0):
    """Two orbitals on the triangular lattice, one dispersing in k_z,
    coupled on site."""
    table = {}

    def add(R, i, j, v):
        table.setdefault(R, np.zeros((2, 2), complex))[i, j] += v / scale

    for o, (e0, t, kz) in enumerate(((-0.75, 0.12, 0.10), (-0.25, -0.10, 0.0))):
        add((0, 0, 0), o, o, e0 + 6 * t - shift)
        for R in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (1, 1, 0), (-1, -1, 0)):
            add(R, o, o, -t)
        if kz:
            add((0, 0, 1), o, o, kz / 2)
            add((0, 0, -1), o, o, kz / 2)
    add((0, 0, 0), 0, 1, 0.05)
    add((0, 0, 0), 1, 0, 0.05)
    return table


def _model(scale=1.0, shift=0.0):
    table = _hoppings(scale, shift)
    return B.WannierModel(LATTICE, np.array(list(table), float),
                          np.array(list(table.values())), orbital_names=["a", "b"])


def _closed_form(k):
    """Band 0 of orbital 0 alone (no coupling): e0 + t h + kz cos(kz c)."""
    kx, ky, kz = k.T
    f = 2 * np.cos(kx * A) + 4 * np.cos(kx * A / 2) * np.cos(np.sqrt(3) * ky * A / 2)
    return -0.75 + 0.12 * (6 - f) + 0.10 * np.cos(kz * CC)


def test_the_hamiltonian_is_the_tight_binding_band():
    table = _hoppings()
    table[(0, 0, 0)][0, 1] = table[(0, 0, 0)][1, 0] = 0      # uncoupled
    model = B.WannierModel(LATTICE, np.array(list(table), float),
                           np.array(list(table.values())))
    k = np.random.default_rng(0).uniform(-1.5, 1.5, (50, 3))
    energies = model.energies(k)
    expected = _closed_form(k)
    assert np.min(np.abs(energies - expected[None]), axis=0).max() < 1e-10


def test_hr_and_bxsf_files_read_back(tmp_path):
    model = _model()
    hr = str(tmp_path / "x_hr.dat")
    bandfiles.write_hr(hr, LATTICE, _hoppings(), 2, fermi=0.1, orbital_names=["a", "b"])
    back = bandfiles.read(hr)
    k = np.random.default_rng(1).uniform(-1, 1, (40, 3))
    assert np.allclose(back.energies(k), model.energies(k) - 0.1, atol=1e-5)
    assert back.orbital_names == ["a", "b"]
    assert os.path.basename(bandfiles.companion_files(hr)[0]) == "x.win"
    bxsf = str(tmp_path / "x.bxsf")
    bandfiles.write_bxsf(bxsf, model, n=(30, 30, 12))
    grid = bandfiles.read(bxsf)
    assert np.abs(grid.energies(k) - model.energies(k)).max() < 0.03


def test_named_path_through_the_hexagonal_zone():
    params = LatticeParams(A, A, CC, 90, 90, 120, space_group=191)
    k, distance, ticks = B.named_path(params, "G-M-K-G-A")
    assert [t[1] for t in ticks] == ["Γ", "M", "K", "Γ", "A"]
    assert np.isclose(ticks[1][0], 2 * np.pi / (np.sqrt(3) * A))
    assert np.isclose(ticks[-1][0] - ticks[-2][0], np.pi / CC)
    assert np.all(np.diff(distance) >= 0)


def test_a_jump_splits_the_drawn_bands():
    params = LatticeParams(A, A, CC, 90, 90, 120, space_group=191)
    k, distance, ticks = B.named_path(params, "G-M|K-G")
    segments = C.band_segments(distance, _model().energies(k), ticks)
    assert len(segments) == 2 * 2                     # two bands, two pieces


def test_the_fit_finds_z_and_the_shift():
    calc = _model(0.8, -0.05)                          # measured = 0.8 calc - 0.05
    truth = _model()
    k = np.linspace(-1.0, 1.0, 120)
    measured = B.bands_along_cut(truth, k, 90.0, kz=0.0)[0]
    result = B.fit_alignment(k, measured,
                             lambda kk: B.bands_along_cut(calc, kk, 90.0, kz=0.0))
    assert result.z == pytest.approx(0.8, abs=0.01)
    assert result.shift == pytest.approx(-0.05, abs=0.005)


def test_the_cut_on_the_sphere_takes_kz_from_the_photon_energy():
    model = _model()
    hv, v0, w = 60.0, 12.0, 4.4
    k = np.linspace(-0.5, 0.5, 21)
    on_sphere = B.bands_along_cut(model, k, 0.0, hv=hv, inner_potential=v0,
                                  work_function=w)
    kz = B.kz_of(0.0, on_sphere[0][10], hv, v0, w)
    flat = B.bands_along_cut(model, np.array([0.0]), 0.0, kz=float(kz))
    assert on_sphere[0][10] == pytest.approx(flat[0][0], abs=1e-6)


def test_contours_lie_on_the_energy():
    model = _model()
    x = np.linspace(-1.2, 1.2, 121)
    lines = B.contours(model, x, x, -0.3, kz=0.0)
    assert lines
    for band, xs, ys in lines:
        k = np.column_stack([xs, ys, np.zeros(len(xs))])
        assert np.abs(model.energies(k)[band] + 0.3).max() < 0.01


def test_the_path_is_assembled_from_the_maps_that_hold_it():
    params = LatticeParams(A, A, CC, 90, 90, 120, space_group=191)
    k, distance, ticks = B.named_path(params, "G-M-K-G-A")
    energy = np.linspace(-1, 0, 11)
    x = np.linspace(-0.4, 1.4, 61)                    # one side of the zone only
    plane = np.ones((x.size, x.size, energy.size))
    kz_axis = np.linspace(3.5, 5.0, 31)
    period = 2 * np.pi / CC
    scan = np.full((kz_axis.size, 21, energy.size), 2.0)
    sources = [C.PlaneMap(plane, x, x, energy, kz=4 * period),
               C.KzMap(scan, kz_axis, np.linspace(-0.3, 0.3, 21), energy)]
    image, owner = C.experiment_along_path(k, energy, sources, fold=6,
                                           kz_period=period)
    # within the plane tolerance (6 % of the period) of k_z = 0 the Γ map serves
    vertical = np.isclose(k[:, 0], 0) & np.isclose(k[:, 1], 0) & (k[:, 2] > 0.065)
    assert np.all(owner[~vertical] == 0)              # in plane: the map
    assert np.all(owner[vertical] == 1)               # Γ-A: the k_z scan
    assert np.allclose(image[vertical], 2.0)


def test_band_points_trace_a_measured_band():
    k = np.linspace(-1, 1, 161)
    e = np.linspace(-1, 0.1, 111)
    band = -0.8 + 0.7 * k ** 2
    values = 1 / ((e[None] - band[:, None]) ** 2 + 0.03 ** 2)
    ks, es = C.band_points(values, k, e)
    assert ks.size > 50
    assert np.abs(es - (-0.8 + 0.7 * ks ** 2)).max() < 0.02


def test_the_node_line_follows_one_crossing_through_kz():
    """Two uncoupled bands with opposite in-plane slopes and different k_z
    terms: the crossing's energy follows the difference."""
    table = {}

    def add(R, i, j, v):
        table.setdefault(R, np.zeros((2, 2), complex))[i, j] += v

    for o, (e0, t, kz) in enumerate(((-0.25, -0.10, -0.06), (-1.10, 0.15, -0.08))):
        add((0, 0, 0), o, o, e0 + 6 * t)
        for R in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (1, 1, 0), (-1, -1, 0)):
            add(R, o, o, -t)
        add((0, 0, 1), o, o, kz / 2)
        add((0, 0, -1), o, o, kz / 2)
    model = B.WannierModel(LATTICE, np.array(list(table), float),
                           np.array(list(table.values())))
    kz = np.linspace(0, np.pi / CC, 5)
    k, e, gap = B.node_line(model, np.linspace(0.2, 0.8, 601), kz, 90.0)
    assert np.all(gap < 0.01)
    # E = -0.25 - 0.1 h - 0.06 c with h = (0.85 + 0.02 c) / 0.25, c = cos(kz c)
    c = np.cos(kz * CC)
    assert np.allclose(e, -0.25 - 0.1 * (0.85 + 0.02 * c) / 0.25 - 0.06 * c, atol=0.005)


def test_two_fitted_branches_cross_at_the_dirac_point():
    from tools.spectra import dispersion as DISP

    rng = np.random.default_rng(3)
    k = np.linspace(0.3, 0.7, 30)
    up = DISP.fit_dispersion(k, 1.5 * (k - 0.5) - 0.6 + rng.normal(0, 0.003, k.size),
                             energy_error=np.full(k.size, 0.003), order=1)
    down = DISP.fit_dispersion(k, -1.2 * (k - 0.5) - 0.6 + rng.normal(0, 0.003, k.size),
                               energy_error=np.full(k.size, 0.003), order=1)
    ck, ce, sk, se = DISP.line_crossing(up, down)
    assert ck == pytest.approx(0.5, abs=3 * sk + 1e-3)
    assert ce == pytest.approx(-0.6, abs=3 * se + 1e-3)
    assert 0 < se < 0.005


def test_a_collapsed_peak_error_does_not_break_the_fit():
    from tools.spectra import dispersion as DISP

    k = np.linspace(0, 1, 10)
    errors = np.full(10, 0.01)
    errors[3] = 1e-13
    fit = DISP.fit_dispersion(k, 0.5 * k - 0.2, k_error=errors, order=1)
    assert fit.covariance is not None and np.all(np.isfinite(fit.errors))


def test_a_stack_of_slices_averages_its_window():
    from tools.processing import volume as V

    x = np.linspace(-1, 1, 21)
    y = np.linspace(-1, 1, 11)
    e = np.linspace(-1, 0, 51)
    cube = np.broadcast_to(e[None, None, :], (21, 11, 51)).copy()
    faces = V.stack_faces(cube, (x, y, e), 2, [-0.5, -0.2], width=0.1)
    assert [f.position for f in faces] == [-0.5, -0.2]
    assert np.allclose(np.nanmean(faces[0].image), -0.5, atol=1e-6)
    assert faces[0].outward == 0 and faces[0].axis == 2
