"""The ALS MAESTRO reader.

The files are built here, in the layout TensorSpec's tests describe: the
three top-level groups, a ``Low_Level_Scan`` table of motors, and the
spectra as a run of ``(angle, energy)`` planes. Every buffer holds
``np.arange``, so a value says exactly which point and pixel it came from,
and each test checks the opened array against the point order the scan kind
is documented to use -- the one thing a shape check cannot catch.
"""
import h5py
import numpy as np
import pytest

from loader.beamlines import maestro as maestro_loader
from loader import registry
from loader.nxs_file import HANDLES


def _row(tag, value):
    return (tag, tag, value, "")


def _plan_rows(name, loops):
    """``loops`` is a list of loops (inner first), each a list of
    ``(motor name, unit, start, end, n)``."""
    rows = [_row("lwlvnm", f"'{name}'"), _row("scanpar", "F"),
            _row("lwlvlpn", str(len(loops)))]
    for i, loop in enumerate(loops):
        rows.append(_row(f"nmsbdv{i}", str(len(loop))))
        for j, (motor, unit, start, end, n) in enumerate(loop):
            key = f"{i}_{j}"
            rows += [_row(f"nm_{key}", f"'{motor}'"),
                     _row(f"un_{key}", f"'{unit}'"),
                     _row(f"st_{key}", str(start)),
                     _row(f"en_{key}", str(end)),
                     _row(f"n_{key}", str(n))]
    return rows


def _write_fixed(path, loops, n_points, n_a=5, n_e=4, points_last=True,
                 recorded=None, name="Scan"):
    """A fixed-mode file; returns the buffer as ``(point, angle, energy)``."""
    recorded = n_points if recorded is None else recorded
    flat = np.arange(recorded * n_a * n_e, dtype=np.float32).reshape(
        recorded, n_a, n_e)
    stored = np.transpose(flat, (1, 2, 0)) if points_last else flat
    with h5py.File(path, "w") as f:
        f.create_group("0D_Data")
        headers = f.create_group("Headers")
        headers.create_dataset("DAQ_Fixed", data=np.array([]))
        headers.create_dataset("Low_Level_Scan",
                               data=np.array(_plan_rows(name, loops),
                                             dtype=object))
        ds = f.create_group("2D_Data").create_dataset(
            "Fixed_Spectra1", data=np.ascontiguousarray(stored))
        ds.attrs["unitNames"] = [b"eV", b"pixels"]
        ds.attrs["scaleOffset"] = [80.0, 0.0]
        ds.attrs["scaleDelta"] = [0.25, 1.0]
    return flat


def _load(path, entry=None):
    return registry.load(str(path), entry,
                         registry.LoadOptions(loader="ALS MAESTRO"))


@pytest.fixture(autouse=True)
def _no_handles_left():
    yield
    for path in list(HANDLES._open):
        HANDLES._open.pop(path)[0].close()


# --------------------------------------------------------------------------
def test_detected_by_structure_and_not_claimed_by_antares(tmp_path):
    path = tmp_path / "scan_00001.h5"
    _write_fixed(path, [[("Slit Defl.", "Deg", -1, 1, 3)]], 3)
    assert registry.detect(str(path)).name == "ALS MAESTRO"
    assert not registry.get_loader("SOLEIL ANTARES").can_open(str(path))

    other = tmp_path / "other.h5"
    with h5py.File(other, "w") as f:
        f.create_group("2D_Data")
    assert not maestro_loader.MaestroLoader().can_open(str(other))


def test_low_level_scan_parsing():
    plan = maestro_loader.parse_low_level_scan(_plan_rows(
        "Focus XY Fine",
        [[("Slit Defl.", "Deg", -1.0, 1.0, 3)],
         [("Scan X", "um", 0, 1, 2), ("Scan Y", "um", 10, 12, 3)]]))
    assert plan.mode_name == "Focus XY Fine"
    assert plan.expected_cycles == 18
    assert [m.name for m in plan.motors_slow_to_fast()] == [
        "Scan Y", "Scan X", "Slit Defl."]
    x, y = plan.xy_motors()
    assert (x.name, y.n) == ("Scan X", 3)
    np.testing.assert_allclose(plan.angle_motors()[0].axis, [-1, 0, 1])


@pytest.mark.parametrize("points_last", [True, False])
def test_deflector_scan_is_a_map(tmp_path, points_last):
    path = tmp_path / "fs.h5"
    flat = _write_fixed(path, [[("Slit Defl.", "Deg", -2, 2, 3)]], 3,
                        points_last=points_last)
    entries = maestro_loader.MaestroLoader().list_entries(str(path))
    assert [(e["entry"], e["kind"]) for e in entries] == [("scan", "Map")]

    scan = _load(path, "scan")
    assert scan.kind == "map"
    assert scan.info["axis0.role"] == "angle"
    np.testing.assert_allclose(scan.x, [-2, 0, 2])
    np.testing.assert_allclose(scan.z, 80 + 0.25 * np.arange(4))
    # Pixels become degrees about the detector's centre.
    np.testing.assert_allclose(scan.k, (np.arange(5) - 2) * 0.048)
    assert scan.labels["x"] == "Slit Defl. (deg)"
    np.testing.assert_array_equal(np.asarray(scan.value), flat)
    scan.close()


def test_xy_mesh_is_a_spatial_scan_with_x_fastest(tmp_path):
    path = tmp_path / "xy.h5"
    nx, ny = 3, 2
    flat = _write_fixed(path, [[("Scan X", "um", 0, 2, nx),
                                ("Scan Y", "um", 10, 11, ny)]], nx * ny)
    scan = _load(path)
    assert scan.kind == "spem_4d"
    np.testing.assert_allclose(scan.x, [0, 1, 2])
    np.testing.assert_allclose(scan.y, [10, 11])
    assert scan.labels["x"] == "Scan X (µm)"
    cube = scan.value4d
    assert cube.shape == (ny, nx, 5, 4)
    for yi in range(ny):
        for xi in range(nx):
            np.testing.assert_array_equal(cube[yi, xi], flat[yi * nx + xi])
    np.testing.assert_array_equal(np.asarray(cube),
                                  flat.reshape(ny, nx, 5, 4))
    scan.close()


def test_two_motor_line_opens_one_map_per_x(tmp_path):
    path = tmp_path / "line.h5"
    nd, nx = 3, 4
    flat = _write_fixed(path, [[("Slit Defl.", "Deg", -1, 1, nd),
                                ("Scan X", "um", 0, 3, nx)]], nd * nx)
    entries = maestro_loader.MaestroLoader().list_entries(str(path))
    assert len(entries) == nx
    assert entries[2]["entry"] == "point 2"
    assert "Scan X = 2 um" in entries[2]["name"]

    scan = _load(path, "point 2")
    assert scan.kind == "map"
    np.testing.assert_allclose(scan.x, [-1, 0, 1])
    # TensorSpec: in this kind the motor listed first steps slowest.
    expected = flat.reshape(nd, nx, 5, 4)[:, 2]
    np.testing.assert_array_equal(np.asarray(scan.value), expected)
    assert scan.fourd_info["Scan X"] == "2 um"
    scan.close()


def test_focus_scan_opens_one_spatial_scan_per_deflector_step(tmp_path):
    path = tmp_path / "focus.h5"
    nd, nx, ny = 3, 2, 2
    flat = _write_fixed(path, [[("Slit Defl.", "Deg", -1, 1, nd)],
                               [("Scan X", "um", 0, 1, nx),
                                ("Scan Y", "um", 10, 11, ny)]],
                        nd * nx * ny, points_last=False)
    entries = maestro_loader.MaestroLoader().list_entries(str(path))
    assert [e["kind"] for e in entries] == ["SPEM"] * nd

    scan = _load(path, "point 1")
    assert scan.kind == "spem_4d"
    # Y outer, X, deflector innermost.
    expected = flat.reshape(ny, nx, nd, 5, 4)[:, :, 1]
    np.testing.assert_array_equal(np.asarray(scan.value4d), expected)
    np.testing.assert_array_equal(scan.value4d[1, 0], expected[1, 0])
    scan.close()


def test_aborted_xy_mesh_keeps_its_complete_rows(tmp_path):
    path = tmp_path / "aborted.h5"
    nx, ny = 3, 4
    flat = _write_fixed(path, [[("Scan X", "um", 0, 2, nx),
                                ("Scan Y", "um", 0, 3, ny)]], nx * ny,
                        recorded=7)
    scan = _load(path)
    assert scan.kind == "spem_4d"
    np.testing.assert_allclose(scan.y, [0, 1])
    np.testing.assert_array_equal(np.asarray(scan.value4d),
                                  flat[:6].reshape(2, nx, 5, 4))
    assert "7 of 12" in scan.info["maestro.aborted"]
    scan.close()


def test_aborted_deflector_scan_is_a_list_of_points(tmp_path):
    path = tmp_path / "aborted_fs.h5"
    flat = _write_fixed(path, [[("Slit Defl.", "Deg", -1, 1, 9)]], 9,
                        recorded=6)
    scan = _load(path)
    assert scan.kind == "map"
    assert scan.info["axis0.role"] == "other"
    np.testing.assert_allclose(scan.x, np.arange(6))
    np.testing.assert_array_equal(np.asarray(scan.value), flat)
    scan.close()


def test_lazy_indexing_matches_numpy(tmp_path):
    path = tmp_path / "xy.h5"
    _write_fixed(path, [[("Scan X", "um", 0, 4, 5),
                         ("Scan Y", "um", 0, 3, 4)]], 20,
                 points_last=True, n_a=6, n_e=7)
    scan = _load(path)
    cube = scan.value4d
    full = np.asarray(cube)
    for key in [(2,), (slice(None), 3), (1, slice(1, 4)),
                (slice(None, None, 2), slice(None), slice(1, 5), 3),
                (Ellipsis, 2), (slice(None), slice(None), slice(None, None, -1)),
                (3, 4, [0, 2, 5]), (-1, -1, -1, -1), (slice(1, 3), [0, 4])]:
        np.testing.assert_array_equal(cube[key], full[key], err_msg=str(key))
    scan.close()


# --------------------------------------------------------------------------
# Process_000 (swept mode)
# --------------------------------------------------------------------------
def _write_process000(path, raw, motors=None, plan=None, fixed=False,
                      units=(b"eV", b"deg")):
    with h5py.File(path, "w") as f:
        zero_d = f.create_group("0D_Data")
        for name, values in (motors or {}).items():
            zero_d.create_dataset(name, data=values)
        headers = f.create_group("Headers")
        headers.create_dataset("DAQ_Fixed" if fixed else "DAQ_Swept",
                               data=np.array([(b"", b"Mode", b"Swept", b"")],
                                             dtype="S32"))
        if plan is not None:
            headers.create_dataset("Low_Level_Scan",
                                   data=np.array(plan, dtype=object))
        ds = f.create_group("2D_Data").create_dataset("Process_000", data=raw)
        ds.attrs["unitNames"] = list(units)
        ds.attrs["scaleOffset"] = [1.0, -2.0]
        ds.attrs["scaleDelta"] = [0.5, 1.0]


def test_single_swept_spectrum_is_a_cut(tmp_path):
    path = tmp_path / "cut.h5"
    raw = np.arange(12, dtype=np.float32).reshape(3, 4)    # (energy, angle)
    _write_process000(path, raw)
    scan = _load(path)
    assert scan.kind == "cut"
    np.testing.assert_allclose(scan.x, [-2, -1, 0, 1])
    np.testing.assert_allclose(scan.y, [1, 1.5, 2])
    np.testing.assert_array_equal(np.asarray(scan.value), raw.T)
    assert scan.info["maestro.mode"] == "Swept"
    scan.close()


def test_process000_motor_grid_with_y_fastest(tmp_path):
    path = tmp_path / "grid.h5"
    raw = np.arange(6 * 3 * 4, dtype=np.float32).reshape(6, 3, 4)
    motors = {"Scan X": [0, 0, 1, 1, 2, 2],
              "Scan Y": [10, 11, 10, 11, 10, 11],
              "Beam Current": [500.0, 499.9, 499.7, 499.8, 499.6, 499.5]}
    plan = _plan_rows("XY", [[("Scan Y", "um", 10, 11, 2),
                              ("Scan X", "um", 0, 2, 3)]])
    _write_process000(path, raw, motors, plan)
    scan = _load(path)
    assert scan.kind == "spem_4d"
    np.testing.assert_allclose(scan.x, [0, 1, 2])
    np.testing.assert_allclose(scan.y, [10, 11])
    # Stored X slow, Y fast; shown (y, x, angle, energy).
    expected = np.transpose(raw.reshape(3, 2, 3, 4), (1, 0, 3, 2))
    np.testing.assert_array_equal(np.asarray(scan.value4d), expected)
    np.testing.assert_array_equal(scan.value4d[1, 2], expected[1, 2])
    scan.close()


def test_process000_irregular_motors_become_points(tmp_path):
    path = tmp_path / "points.h5"
    raw = np.arange(5 * 3 * 4, dtype=np.float32).reshape(5, 3, 4)
    _write_process000(path, raw, {"Sample Theta": [0, 1, 3, 2, 5]})
    scan = _load(path)
    assert scan.kind == "map"
    assert scan.info["axis0.role"] == "other"
    np.testing.assert_array_equal(np.asarray(scan.value),
                                  np.transpose(raw, (0, 2, 1)))
    scan.close()


def test_process000_one_angle_motor_is_a_convertible_map(tmp_path):
    path = tmp_path / "theta.h5"
    raw = np.arange(4 * 3 * 5, dtype=np.float32).reshape(4, 3, 5)
    _write_process000(path, raw, {"Sample Theta": [-3, -1, 1, 3],
                                  "Sample Z": [7, 7, 7, 7]})
    scan = _load(path)
    assert scan.kind == "map"
    assert scan.info["axis0.role"] == "angle"
    np.testing.assert_allclose(scan.x, [-3, -1, 1, 3])
    assert scan.info["maestro.Sample Z"] == 7
    scan.close()


def test_file_closes_when_the_scan_does(tmp_path):
    path = tmp_path / "fs.h5"
    _write_fixed(path, [[("Slit Defl.", "Deg", -1, 1, 3)]], 3)
    scan = _load(path)
    assert HANDLES.is_open(str(path))
    scan.close()
    assert not HANDLES.is_open(str(path))
