"""Tests for the ANTARES photon-energy (PhD) scan reader.

The synthetic files copy the layout of a real one (``NbSCl_..._.nxs``, from
its h5info dump): one top-level ``Phd_NNNN`` entry per photon energy, each
holding a single MBS frame in ``MBS_0_<n>`` with its own scales and ``hv``,
and a first entry with the beamline state and no frame at all.
"""
import h5py
import numpy as np
import pytest

import arpes_viewer.loader.registry as L
from arpes_viewer.loader.nxs_file import list_datasets, load_soleil_nxs

WORK_FUNCTION = 4.5
N_ANGLE, N_ENERGY = 40, 61
E_STEP = 0.01


def _frame(hv, e_min):
    """A Fermi edge at E_kin = hv - W across every angle, stored the way the
    file stores it: (energy, angle), the transpose of how it is read."""
    energy = e_min + E_STEP * np.arange(N_ENERGY)
    edge = 1.0 / (1.0 + np.exp((energy - (hv - WORK_FUNCTION)) / 0.02))
    return (1000 * edge[:, None] * np.ones(N_ANGLE)[None, :]).astype(np.int32)


def _write_phd(path, photon_energies, window=lambda hv: hv - WORK_FUNCTION - 0.4):
    with h5py.File(path, "w") as f:
        first = f.create_group("Phd_0001")
        first["title"] = "NbSCl"
        first.create_group("ANTARES/i12-m-c04-op-mono1")["energy"] = [photon_energies[0]]
        for index, hv in enumerate(photon_energies):
            entry = f.create_group(f"Phd_{index + 2:04d}")
            entry["title"] = "NbSCl"
            entry["start_time"] = "2026-10-21T23:18:59"
            entry.create_group("ANTARES/i12-m-c04-op-mono1")["energy"] = [hv]
            entry.create_group("User")["name"] = "someone"
            mbs = entry.create_group(f"MBS_0_{index}")
            e_min = window(hv)
            mbs["Image32_1"] = _frame(hv, e_min)
            mbs["Escale_Min_1"] = [e_min]
            mbs["Escale_Mult_1"] = [E_STEP]
            mbs["Escale_Max_1"] = [e_min + E_STEP * (N_ENERGY - 1)]
            mbs["Xscale_Min_1"] = [-10.0]
            mbs["Xscale_Mult_1"] = [0.5]
            mbs["Xscale_Max_1"] = [-10.0 + 0.5 * (N_ANGLE - 1)]
            mbs["hv"] = [hv]
            mbs["Theta"] = [2.0]
            mbs["Phi"] = [float(index)]
            mbs["Xproj_1"] = np.zeros(N_ANGLE)
            mbs["Yproj_1"] = np.zeros(N_ENERGY)
    return str(path)


def _edge_position(scan):
    """Where the edge sits on the cube's energy axis, per photon energy."""
    edc = np.asarray(scan.value, dtype=float).mean(axis=1)
    return np.array([np.interp(500, row[::-1], scan.z[::-1]) for row in edc])


@pytest.fixture()
def phd_file(tmp_path):
    return _write_phd(tmp_path / "NbSCl.nxs", [60.0, 62.0, 64.0, 66.0])


def test_the_frames_are_listed_as_one_kz_map(phd_file):
    rows = list_datasets(phd_file)
    assert len(rows) == 1
    assert rows[0]["entry"] == "Phd_*"
    assert rows[0]["kind"] == "kz map"
    assert rows[0]["title"] == "NbSCl"


def test_the_antares_loader_recognises_it(phd_file):
    assert L.detect(phd_file).name == "SOLEIL ANTARES"


def test_it_loads_as_an_hv_angle_energy_cube(phd_file):
    scan = L.load(phd_file, "Phd_*")
    assert scan.kind == "kz_map"
    assert scan.value.shape == (4, N_ANGLE, N_ENERGY)
    np.testing.assert_allclose(scan.x, [60, 62, 64, 66])
    np.testing.assert_allclose(scan.k, -10 + 0.5 * np.arange(N_ANGLE))
    assert scan.info["axis0.role"] == "photon_energy"
    assert scan.labels["x"] == "Photon energy (eV)"
    assert scan.info["antares.phd_members"] == 4
    assert scan.info["work_function_eV"] == WORK_FUNCTION


def test_kinetic_windows_that_follow_hv_are_stacked_untouched(phd_file):
    scan = load_soleil_nxs(phd_file, "Phd_*")
    assert "antares.energy_resampled" not in scan.info
    np.testing.assert_allclose(scan.z[0], -0.4, atol=1e-9)
    np.testing.assert_allclose(_edge_position(scan), 0.0, atol=2e-3)


def test_drifting_windows_are_resampled_onto_their_overlap(tmp_path):
    drift = {60.0: -0.3, 62.0: -0.137, 64.0: -0.25, 66.0: -0.401}
    path = _write_phd(tmp_path / "drift.nxs", list(drift),
                      window=lambda hv: hv - WORK_FUNCTION + drift[hv])
    scan = load_soleil_nxs(path, "Phd_*")
    assert "antares.energy_resampled" in scan.info
    assert scan.z[0] >= -0.137 - 1e-9
    assert scan.z[-1] <= -0.401 + E_STEP * (N_ENERGY - 1) + 1e-9
    np.testing.assert_allclose(_edge_position(scan), 0.0, atol=3e-3)


def test_no_entry_named_still_reads_the_whole_scan(phd_file):
    assert load_soleil_nxs(phd_file).kind == "kz_map"


def test_one_frame_opens_on_its_own_as_a_cut(phd_file):
    scan = load_soleil_nxs(phd_file, "Phd_0003")
    assert scan.kind == "cut"
    assert scan.value.shape == (N_ANGLE, N_ENERGY)
    np.testing.assert_allclose(scan.y[0], 62 - WORK_FUNCTION - 0.4)


def test_frames_out_of_order_are_sorted_by_photon_energy(tmp_path):
    path = _write_phd(tmp_path / "shuffled.nxs", [64.0, 60.0, 66.0, 62.0])
    scan = load_soleil_nxs(path, "Phd_*")
    np.testing.assert_allclose(scan.x, [60, 62, 64, 66])
    np.testing.assert_allclose(_edge_position(scan), 0.0, atol=2e-3)


def test_a_binding_energy_axis_is_left_as_it_is(tmp_path):
    path = _write_phd(tmp_path / "be.nxs", [60.0, 62.0],
                      window=lambda hv: -0.4)
    scan = load_soleil_nxs(path, "Phd_*")
    np.testing.assert_allclose(scan.z[0], -0.4)
    assert scan.labels["z"] == "Energy (eV)"
    assert "work_function_eV" not in scan.info
