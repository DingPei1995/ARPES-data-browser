"""SPEM display: ANTARES orientation, plain axis units, and the memory helpers.

* An ANTARES spatial scan opens the way the beamline draws it: coarse
  (ST/SZ) with both axes reversed, fine (PIX/PIY) with only X reversed.
  Other data keeps the plotting default.
* Axis ticks are in the unit the title names -- no "(x0.001)" rescaling.
* tools.memory reports something and never raises.

Needs Qt; run headless with QT_QPA_PLATFORM=offscreen.
"""
import os

import h5py
import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PyQt5")
from PyQt5.QtWidgets import QApplication                          # noqa: E402

from loader.nxs_file import (antares_spatial_orientation, load_soleil_nxs,  # noqa: E402
                             _antares_stage)
from tools import memory as memory_tools                           # noqa: E402
from ui.widgets import NxsData, MemoryData                         # noqa: E402
from ui import windows as W                                        # noqa: E402


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


def _write_spem(path, names, x_range, y_range):
    """A minimal ANTARES real-space scan: two trajectory actuators, the
    energy/slit scales and a small 4-D cube."""
    nx = int(round((x_range[1] - x_range[0]) / x_range[2])) + 1
    ny = int(round((y_range[1] - y_range[0]) / y_range[2])) + 1
    nk, ne = 5, 7
    with h5py.File(path, "w") as f:
        g = f.create_group("scan")
        sd = g.create_group("scan_data")
        sd["data_12"] = np.random.default_rng(0).random((ny, nx, nk, ne))
        for name, value in (("data_04", 10.0), ("data_05", 0.5), ("data_06", 13.0),
                            ("data_07", -2.0), ("data_08", 1.0), ("data_09", 2.0)):
            sd[name] = np.array([value])
        traj = g.create_group("scan_config/trajectory")
        for i, (name, (lo, hi, step)) in enumerate(zip(names, (x_range, y_range)), 1):
            act = traj.create_group(f"actuator_{i}_1")
            act["name"] = np.bytes_(name)
            act["from"] = np.array([lo])
            act["to"] = np.array([hi])
            act["delta"] = np.array([step])


def test_stage_names():
    assert _antares_stage("ST") == "st"
    assert _antares_stage("i12-m-cx1-ex-sample-mt_sz") == "sz"
    assert _antares_stage("PIX") == "pix"
    assert _antares_stage("piy") == "piy"
    assert _antares_stage("theta") is None
    assert _antares_stage(None) is None


def test_orientation_table():
    assert antares_spatial_orientation("ST", "SZ") == {
        "Spatial.invert_x": True, "Spatial.invert_y": True}
    assert antares_spatial_orientation("PIX", "PIY") == {
        "Spatial.invert_x": True, "Spatial.invert_y": False}
    # not ANTARES stages, or not where they usually are: left alone
    assert antares_spatial_orientation("motor1", "motor2") == {}
    assert antares_spatial_orientation("SZ", "ST") == {}


@pytest.mark.parametrize("names, expected", [
    (("ST", "SZ"), (True, True)),
    (("PIX", "PIY"), (True, False)),
])
def test_antares_window_orientation(app, tmp_path, names, expected):
    path = str(tmp_path / "spem.nxs")
    _write_spem(path, names, (-0.9, 0.1, 0.05), (39.6, 40.2, 0.05))
    data = NxsData.acquire(path)
    window = W.SpatialScanWindow(data, "spem", "gray", False)
    try:
        view = window.spatial_view
        assert (view.axis_inverted(0), view.axis_inverted(1)) == expected
        assert window.view_bar.invert_boxes[0].isChecked() == expected[0]
        # swapping the display carries each axis's direction with it
        window.set_swap_xy(True)
        assert (view.axis_inverted(0), view.axis_inverted(1)) == expected[::-1]
    finally:
        window.close()


def test_other_spem_keeps_default(app):
    x, y = np.linspace(0, 1, 6), np.linspace(0, 2, 5)
    k, e = np.linspace(-1, 1, 4), np.linspace(0, 1, 3)
    data = MemoryData("spem_4d", (x, y, k, e),
                      np.random.default_rng(1).random((y.size, x.size, k.size, e.size)),
                      {"x": "X (mm)", "y": "Y (mm)", "k": "k", "z": "E"},
                      source_label="s")
    window = W.SpatialScanWindow(data, "spem", "gray", False)
    try:
        view = window.spatial_view
        assert (view.axis_inverted(0), view.axis_inverted(1)) == (False, False)
    finally:
        window.close()


def test_axis_ticks_in_labelled_unit(app, tmp_path):
    path = str(tmp_path / "spem.nxs")
    _write_spem(path, ("ST", "SZ"), (-0.9, 0.1, 0.05), (39.6, 40.2, 0.05))
    scan = load_soleil_nxs(path)
    try:
        assert scan.labels["x"] == "X (mm)"
    finally:
        scan.close()
    data = NxsData.acquire(path)
    window = W.SpatialScanWindow(data, "spem", "gray", False)
    try:
        axis = window.spatial_view.view.getAxis("bottom")
        app.processEvents()
        assert axis.autoSIPrefixScale == 1.0
        assert "x0.001" not in axis.labelString()
    finally:
        window.close()


def test_memory_helpers():
    rss, _peak = memory_tools.process_memory()
    assert rss is None or rss > 0
    total, available = memory_tools.system_memory()
    assert total is None or total > 0
    assert memory_tools.format_bytes(None) == "n/a"
    assert memory_tools.format_bytes(1536) == "1.50 KB"
    assert memory_tools.release_memory() >= 0
