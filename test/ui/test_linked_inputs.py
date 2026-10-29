"""Greyed-out inputs show what is actually used, never a stale number.

Where a dialog offers two ways to give the same thing -- points or step
size, evenly spaced or typed positions, the bottom layer's own lattice or
the top one's -- the half that is not being typed is disabled. It used to
keep whatever it held when the dialog opened (0.01 Å⁻¹, 3.0 Å ...), which
read as the result. Now it follows: these tests change the live half and
read the greyed one back.
"""
import numpy as np
import pytest

pytest.importorskip("PyQt5")
from PyQt5.QtWidgets import QApplication                      # noqa: E402

from test.ui.test_defaults import _cut, _map, DEFL, SLIT, EKIN  # noqa: E402


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


def test_map_points_drive_the_greyed_steps_and_back(app):
    from arpes_viewer.ui.viewers import windows as W
    window = W.ContourWindow(_map(), "m", "gray", False)
    dialog = W.KConversionDialog(window)
    kx_lo, kx_hi, ky_lo, ky_hi = dialog.k_extent()
    assert not dialog.dkx_box.isEnabled()
    dialog.nkx_box.setValue(11)
    assert dialog.dkx_box.value() == pytest.approx((kx_hi - kx_lo) / 10, abs=1e-5)
    assert dialog.dky_box.value() == pytest.approx((ky_hi - ky_lo) / (SLIT.size - 1), abs=1e-5)
    span_meV = (EKIN[-1] - EKIN[0]) * 1000
    assert dialog.de_box.value() == pytest.approx(span_meV / (EKIN.size - 1), abs=1e-3)
    # the origin moves the k range, and the step a count gives with it
    before = dialog.dkx_box.value()
    dialog.theta_box.setValue(4.0)
    assert dialog.dkx_box.value() != pytest.approx(before, abs=1e-6)

    dialog.by_resolution.setChecked(True)
    assert not dialog.nkx_box.isEnabled() and dialog.dkx_box.isEnabled()
    # switching keeps the grid: the steps shown give back the same counts
    assert dialog.grid_counts()[0] == 11
    dialog.dkx_box.setValue(dialog.dkx_box.value() / 2)
    assert dialog.nkx_box.value() == dialog.grid_counts()[0] == 21
    assert dialog.settings()["n_kx"] == 21
    dialog.close()
    window.close()


def test_cut_points_drive_the_greyed_steps_and_back(app):
    from arpes_viewer.ui.viewers import windows as W
    window = W.CutWindow(_cut(), "c", "gray", False)
    dialog = W.CutKConversionDialog(window)
    lo, hi = dialog.k_extent()
    dialog.nk_box.setValue(31)
    assert dialog.dk_box.value() == pytest.approx((hi - lo) / 30, abs=1e-5)
    dialog.by_resolution.setChecked(True)
    dialog.dk_box.setValue((hi - lo) / 60)
    assert dialog.nk_box.value() == 61 == dialog.settings()["n_k"]
    dialog.close()
    window.close()


def test_the_step_boxes_are_called_steps(app):
    from PyQt5.QtWidgets import QLabel
    from arpes_viewer.ui.viewers import windows as W
    window = W.ContourWindow(_map(), "m", "gray", False)
    dialog = W.KConversionDialog(window)
    texts = [label.text() for label in dialog.findChildren(QLabel)]
    assert "kx step (Å⁻¹)" in texts and "E step (meV)" in texts
    assert not any("resolution (" in t for t in texts)
    assert dialog.by_resolution.text() == "by step size"
    dialog.close()
    window.close()


def test_slice_positions_show_each_other(app):
    from arpes_viewer.ui.panels.figure import SliceSeriesDialog
    axes = (DEFL, SLIT, EKIN)
    dialog = SliceSeriesDialog(None, axes, ("x", "y", "E"))
    dialog.count_box.setValue(3)
    dialog.from_box.setValue(55.0)
    dialog.to_box.setValue(55.5)
    assert dialog.values_edit.text() == "55, 55.25, 55.5"
    dialog.mode_typed.setChecked(True)
    dialog.values_edit.setText("55.6, 54.8")
    dialog._refresh()
    assert dialog.count_box.value() == 2
    assert dialog.from_box.value() == pytest.approx(54.8)
    assert dialog.to_box.value() == pytest.approx(55.6)
    dialog.close()


def test_switching_value_and_index_rewrites_the_bounds(app):
    from arpes_viewer.ui.viewers.windows import DataOperationsDialog
    dialog = DataOperationsDialog([("c", _cut())])
    lo, hi = dialog.trunc_boxes[0]
    axis = np.asarray(dialog.axes[0], dtype=float)
    lo.setText(f"{axis[4]:g}")
    hi.setText("")
    dialog.trunc_by_index.setChecked(True)
    assert lo.text() == "5" and hi.text() == ""
    dialog.trunc_by_value.setChecked(True)
    assert float(lo.text()) == pytest.approx(axis[4], abs=1e-4)
    dialog.close()


def test_the_bottom_layer_shows_the_top_lattice_while_the_same(app, tmp_path, monkeypatch):
    from arpes_viewer.ui.viewers import windows as W
    from arpes_viewer.ui.viewers.widgets import MemoryData
    monkeypatch.setattr(W.BrillouinZoneDialog, "last_used_path",
                        classmethod(lambda cls: str(tmp_path / "bz.json")))
    k = np.linspace(-1.2, 1.2, 21)
    data = MemoryData("k_map", (k, k, np.linspace(-0.5, 0.1, 5)),
                      np.ones((21, 21, 5)), {"x": "kx", "y": "ky", "z": "E"},
                      source_label="t")
    contour = W.ContourWindow(data, "t", "gray", False)
    dialog = contour.open_brillouin_zone()
    dialog.mode_2d.setChecked(True)
    dialog.moire_cb.setChecked(True)
    dialog.moire_same_cb.setChecked(True)
    dialog.a2_box.setValue(2.46)
    dialog.gamma2_box.setValue(60.0)
    assert not dialog.a2b_box.isEnabled()
    assert dialog.a2b_box.value() == pytest.approx(2.46)
    assert dialog.gamma2b_box.value() == pytest.approx(60.0)
    dialog.moire_same_cb.setChecked(False)
    assert dialog.a2b_box.isEnabled() and dialog.a2b_box.value() == pytest.approx(2.46)
    dialog.close()
    contour.close()
