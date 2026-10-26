"""The Elettra APE (Scienta SES) and ASTRID2 SGM4 (HDF5) readers.

No files from either beamline were available, so these are built as the
peaks readers they were ported from expect them: APE's SES files are
Bloch's with ``Location`` saying APE and ``polar tilt X Y Z`` motors; an
SGM4 file is ``Entry/Data/ScanDetails`` (the axes as start, step, length)
beside ``Entry/Data/TransformedData``, one row per point with the first
slow axis fastest.
"""
import h5py
import numpy as np
import pytest

from loader import registry, sgm4_loader
from test.test_cassiopee import write_spectrum
from test.test_maxiv_loader import _cube, _note, _ses_zip, _write_ibw

NA, NE = 5, 6


# --------------------------------------------------------------------------
# APE
# --------------------------------------------------------------------------
def _ape_note(**kwargs):
    return _note(**kwargs).replace("Location=Bloch", "Location=APE")


def test_ape_cut_is_told_from_bloch(tmp_path):
    ape, bloch = str(tmp_path / "ape.ibw"), str(tmp_path / "bloch.ibw")
    data = _cube((NA, NE))
    for path, note in ((ape, _ape_note()), (bloch, _note())):
        _write_ibw(path, data, (38.2, -16.9), (0.003, 0.0474),
                   ("Kinetic Energy [eV]", "Y-Scale [deg]"), note)
    assert registry.detect(ape).name == "Elettra APE"
    assert registry.detect(bloch).name == "MAX IV Bloch"
    scan = registry.load(ape)
    assert scan.kind == "cut" and scan.info["facility"] == "Elettra APE"
    np.testing.assert_array_equal(np.asarray(scan.value), data)
    # Chosen by hand, the Bloch file opens with the APE reader too.
    forced = registry.load(bloch, options=registry.LoadOptions(
        loader="Elettra APE"))
    assert forced.info["facility"] == "Elettra APE"


def test_ape_spatial_scan_with_its_motor_names(tmp_path):
    xs, zs = np.arange(4) * 0.1, 1 + np.arange(3) * 0.2
    rows = ["Name=Manipulator Scan",
            "\x0b".join(["Position [unit]", "Point [unit]", "polar", "tilt",
                         "X", "Y", "Z"])]
    rows += ["\x0b".join(f"{v:.3f}" for v in (z, x, 0, 0, x, 5, z))
             for z in zs for x in xs]
    path = str(tmp_path / "map.ibw")
    data = _cube((3, 4, NA, NE))
    _write_ibw(path, data, (38.2, -16.9, 0.0, 1.0), (0.0015, 0.0474, 0.1, 0.2),
               ("Kinetic Energy [eV]", "Y-Scale [deg]", "Point [unit]",
                "Position [unit]"), _ape_note(run_mode=rows))
    scan = registry.load(path)
    assert scan.kind == "spem_4d"
    assert scan.labels["x"] == "X (mm)" and scan.labels["y"] == "Z (mm)"
    np.testing.assert_array_equal(np.asarray(scan.value4d), data)


def test_ape_polar_scan_is_a_convertible_map(tmp_path):
    polar = np.arange(4) * 2.0
    rows = ["Name=Manipulator Scan",
            "\x0b".join(["Point [unit]", "polar", "tilt", "X", "Y", "Z"])]
    rows += ["\x0b".join(f"{v:.3f}" for v in (p, p, 0, 0, 0, 0))
             for p in polar]
    path = str(tmp_path / "polar.ibw")
    _write_ibw(path, _cube((4, NA, NE)), (38.2, -16.9, 0.0),
               (0.0015, 0.0474, 2.0), ("Kinetic Energy [eV]",
                                       "Y-Scale [deg]", "Point [unit]"),
               _ape_note(run_mode=rows))
    scan = registry.load(path)
    assert scan.kind == "map" and registry.role_is_angle(scan)
    assert scan.labels["x"].startswith("Polar angle")


def test_ape_zip(tmp_path):
    path = str(tmp_path / "fs.zip")
    data = _cube((4, NA, NE))
    _ses_zip(path, {"fs": data})
    import zipfile
    with zipfile.ZipFile(path) as z:
        members = {n: z.read(n) for n in z.namelist()}
    members["fs.ini"] = members["fs.ini"].replace(b"Location=Bloch",
                                                  b"Location=APE")
    with zipfile.ZipFile(path, "w") as z:
        for name, content in members.items():
            z.writestr(name, content)
    assert registry.detect(path).name == "Elettra APE"
    np.testing.assert_array_equal(np.asarray(registry.load(path).value), data)


def test_ape_text_spectrum(tmp_path):
    angle, energy = np.linspace(-15, 15, NA), np.linspace(20, 21, NE)
    values = np.arange(NA * NE, dtype=float).reshape(NA, NE)
    path = write_spectrum(str(tmp_path / "Ag_001.txt"), energy, angle, values)
    with open(path, encoding="latin-1") as f:
        content = f.read()
    soleil = str(tmp_path / "soleil.txt")
    with open(soleil, "w", encoding="latin-1", newline="") as f:
        f.write(content)
    with open(path, "w", encoding="latin-1", newline="") as f:
        f.write(content.replace("Location=Soleil", "Location=Elettra APE"))
    assert registry.detect(path).name == "Elettra APE"
    assert registry.detect(soleil).name == "SOLEIL CASSIOPEE"
    [entry] = registry.list_entries(path)
    assert entry["kind"] == "Cut"
    scan = registry.load(path, entry["entry"])
    assert scan.kind == "cut" and scan.info["facility"] == "Elettra APE"
    np.testing.assert_allclose(np.asarray(scan.value), values)
    np.testing.assert_allclose(scan.y, energy, atol=1e-4)


# --------------------------------------------------------------------------
# SGM4
# --------------------------------------------------------------------------
def _write_sgm4(path, slow, fast=None, points=None, energy_first=True):
    """``slow``: ``[(name, start, step, length)]`` as ScanDetails lists
    them; the spectra are ``np.arange`` on the full grid, ``(slow[-1], ...,
    slow[0], angle, energy)``, flattened to ``points`` rows."""
    fast = fast or [("Kinetic Energy", 20.0, 0.01, NE),
                    ("OrdinateRange", -15.0, 1.0, NA)]
    if not energy_first:
        fast = fast[::-1]
    grid = tuple(n for *_, n in reversed(slow))
    cube = _cube(grid + (NA, NE))
    with h5py.File(path, "w") as f:
        details = f.create_group("Entry/Data/ScanDetails")
        for which, axes in (("SlowAxis", slow), ("FastAxis", fast)):
            details[f"{which}_names"] = np.array([a[0].encode()
                                                  for a in axes])
            details[f"{which}_start"] = np.array([a[1] for a in axes], float)
            details[f"{which}_step"] = np.array([a[2] for a in axes], float)
            details[f"{which}_length"] = np.array([a[3] for a in axes])
        rows = cube.reshape((-1, NA, NE))
        if energy_first:
            rows = np.swapaxes(rows, 1, 2)
        if points is not None:
            rows = rows[:points]
        f["Entry/Data/TransformedData"] = rows
        f["Entry/Data/Command"] = b"scan FAT 20"
        f["Entry/Data/Timestamp"] = np.array([b"2024-05-01 10:00:00"])
        f["Entry/Instrument/Detector/Info"] = (
            b"Kinetic Start Energy: 20\r\nKinetic End Energy: 20.05\r\n"
            b"Samples: 3\r\nStep Width: 0.01\r\nPass Energy: 10\r\n"
            b"Dwell Time: 0.1\r\nLens Mode: WideAngleMode")
        f["Entry/Instrument/Monochromator/Photon Energy"] = np.array([47.0])
        f["Entry/Instrument/Temperature/Sample Stage"] = np.array([30.0, 32.0])
        f["Entry/Instrument/Positioner/SamPolar"] = np.array([46000.0])
    return cube


@pytest.mark.parametrize("energy_first", [True, False])
def test_sgm4_deflector_map(tmp_path, energy_first):
    path = str(tmp_path / "map.h5")
    cube = _write_sgm4(path, [("ShiftX", -5.0, 2.5, 5)],
                       energy_first=energy_first)
    assert registry.detect(path).name == "ASTRID2 SGM4"
    [entry] = registry.list_entries(path)
    assert entry["kind"] == "Map"
    scan = registry.load(path)
    assert scan.kind == "map" and registry.role_is_angle(scan)
    np.testing.assert_allclose(scan.x, np.linspace(-5, 5, 5))
    np.testing.assert_allclose(scan.k, -15 + np.arange(NA))
    np.testing.assert_allclose(scan.z, 20 + 0.01 * np.arange(NE))
    np.testing.assert_array_equal(np.asarray(scan.value), cube)
    np.testing.assert_array_equal(scan.value[3, 1:4, 2], cube[3, 1:4, 2])
    assert scan.info["photon_energy_eV"] == 47.0
    assert scan.info["analyser.pass_energy"] == 10.0
    assert scan.info["analyser.lens_mode"] == "WideAngleMode"
    assert scan.info["SampleTemperature_K"] == 31.0
    assert scan.info["manipulator.polar_deg"] == pytest.approx(1.0)
    scan.close()


def test_sgm4_spatial_scan_first_slow_axis_fastest(tmp_path):
    path = str(tmp_path / "xy.h5")
    cube = _write_sgm4(path, [("FSamX", 0.0, 1.0, 4), ("FSamY", 10.0, 2.0, 3)])
    scan = registry.load(path)
    assert scan.kind == "spem_4d"
    np.testing.assert_allclose(scan.x, np.arange(4.0))
    np.testing.assert_allclose(scan.y, 10 + 2 * np.arange(3))
    assert scan.labels["x"] == "X (µm)" and scan.labels["y"] == "Y (µm)"
    np.testing.assert_array_equal(np.asarray(scan.value4d), cube)
    np.testing.assert_array_equal(scan.value4d[2, 1], cube[2, 1])
    np.testing.assert_array_equal(scan.value4d[:, 3, 2, -1], cube[:, 3, 2, -1])
    scan.close()


def test_sgm4_stopped_scan_is_padded_with_nan(tmp_path):
    path = str(tmp_path / "part.h5")
    cube = _write_sgm4(path, [("FSamX", 0.0, 1.0, 4), ("FSamY", 0.0, 1.0, 3)],
                       points=6)
    scan = registry.load(path)
    full = np.asarray(scan.value4d)
    np.testing.assert_array_equal(full[0], cube[0])
    np.testing.assert_array_equal(full[1, :2], cube[1, :2])
    assert np.isnan(full[1, 2:]).all() and np.isnan(full[2]).all()
    assert "6 of 12" in scan.info["sgm4.incomplete"]
    scan.close()


def test_sgm4_photon_energy_scan(tmp_path):
    """One binding window: each hv's spectra stack as they are, on
    E_kin - hv0 + W with SGM4's own work function."""
    path = str(tmp_path / "hv.h5")
    cube = _write_sgm4(path, [("MEG", 40.0, 2.0, 3)])
    [entry] = registry.list_entries(path)
    assert entry["kind"] == "kz map"
    scan = registry.load(path)
    assert scan.kind == "kz_map"
    np.testing.assert_allclose(scan.x, [40, 42, 44])
    w = sgm4_loader.WORK_FUNCTION_EV
    assert scan.info["work_function_eV"] == w
    np.testing.assert_allclose(scan.z, 20 + 0.01 * np.arange(NE) - 40 + w)
    np.testing.assert_array_equal(scan.value, cube)


def test_sgm4_repeats_are_listed_one_per_loop(tmp_path):
    path = str(tmp_path / "loops.h5")
    cube = _write_sgm4(path, [("ShiftX", 0.0, 1.0, 4), ("Loop", 0, 1, 2)])
    entries = registry.list_entries(path)
    assert len(entries) == 2
    scan = registry.load(path, entries[1]["entry"])
    assert scan.kind == "map"
    np.testing.assert_array_equal(np.asarray(scan.value), cube[1])
    scan.close()


def test_sgm4_knife_edge_is_listed_but_not_opened(tmp_path):
    path = str(tmp_path / "knife.h5")
    with h5py.File(path, "w") as f:
        details = f.create_group("Entry/Data/ScanDetails")
        details["SlowAxis_names"] = np.array([b"SamX"])
        details["SlowAxis_start"] = np.array([0.0])
        details["SlowAxis_step"] = np.array([1.0])
        details["SlowAxis_length"] = np.array([5])
        f["Entry/Process/SumData"] = np.arange(5.0)
    assert registry.detect(path).name == "ASTRID2 SGM4"
    [entry] = registry.list_entries(path)
    assert entry["kind"] == "Knife edge"
    with pytest.raises(NotImplementedError, match="knife-edge"):
        registry.load(path)
