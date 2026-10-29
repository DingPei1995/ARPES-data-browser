"""The settings each tool opens with are the conventional ones.

* Self-normalisation starts on the usual choice for the kind of data --
  each EDC of a cut, each cut (slit × energy) of a map or kz map -- and
  says in words what the ticked axes do. It used to start with nothing
  ticked, and "along" was read the wrong way round.
* The energy window of a normalisation can be on any of the summed axes,
  not only the first one ticked.
* The 3-D panel's smoothing widths start non-zero and its operations are
  named after the plane's real axes.
* The 2-D mirror line starts vertical (k -> -k), and the Fermi cut-off
  starts at an E_F on the energy axis.
* The map k conversion keeps the measured number of angles, and calls the
  origin "Γ deflector/slit angle" like the cut conversion does.
* Joining pieces starts on the axis they differ along.

Needs Qt for the dialog tests; run headless with QT_QPA_PLATFORM=offscreen.
"""
import os

import numpy as np
import pytest

from arpes_viewer.tools.processing import combine as C
from arpes_viewer.tools.processing.dataops import normalise_presets, self_normalize

DEFL = np.linspace(-8, 8, 33)
SLIT = np.linspace(-12, 12, 61)
EKIN = np.linspace(54.6, 55.7, 45)
EF = 55.6


# -- no Qt ------------------------------------------------------------------
def test_normalise_presets_start_on_the_conventional_choice():
    assert normalise_presets("cut", 2)[0][1] == (1,)            # each EDC
    for kind in ("map", "k_map", "kz_map", "kz_map_k"):
        assert normalise_presets(kind, 3)[0][1] == (1, 2)       # each cut
    assert normalise_presets("spem_4d", 4)[0][1] == (2, 3)
    assert normalise_presets("something new", 3)[0][1] == (2,)


def test_each_cut_of_a_map_is_divided_by_its_own_total():
    rng = np.random.default_rng(1)
    cube = rng.random((4, 5, 6)) + 1
    cube *= np.arange(1, 5)[:, None, None]                      # flux drift
    axes = [np.arange(n, dtype=float) for n in cube.shape]
    out = self_normalize(cube, axes, (1, 2))
    np.testing.assert_allclose(out.sum(axis=(1, 2)), 1.0)


def test_the_window_can_be_on_the_second_summed_axis():
    cube = np.ones((2, 3, 10))
    axes = [np.arange(2.0), np.arange(3.0), np.arange(10.0)]
    out = self_normalize(cube, axes, (1, 2), window=(0, 4), window_dim=2)
    np.testing.assert_allclose(out, 1.0 / 15)                   # 3 × 5 points
    with pytest.raises(ValueError):
        self_normalize(cube, axes, (2,), window=(0, 4), window_dim=0)


def test_the_join_axis_guess_follows_what_differs():
    angle, energy = np.linspace(-5, 5, 11), np.linspace(0, 1, 21)
    cuts = [(angle, energy[:12]), (angle, energy[9:])]
    assert C.guess_join_axis(cuts) == 1
    maps = [(np.arange(0, 10.0), angle, energy),
            (np.arange(8, 20.0), angle, energy)]
    assert C.guess_join_axis(maps) == 0
    assert C.guess_join_axis([(angle, energy), (angle, energy)]) == 0


# -- dialogs ----------------------------------------------------------------
@pytest.fixture(scope="module")
def app():
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    pytest.importorskip("PyQt5")
    from PyQt5.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def _cut():
    from arpes_viewer.ui.viewers.widgets import MemoryData
    rng = np.random.default_rng(7)
    edge = 1.0 / (np.exp((EKIN - EF) / 0.01) + 1.0)
    values = rng.random((SLIT.size, EKIN.size)) * 0.1 + 1 + 10 * edge
    return MemoryData("cut", (SLIT, EKIN), values,
                      {"x": "Angle, slit (deg)", "y": "Kinetic energy (eV)"},
                      source_label="c")


def _map(kind="map"):
    from arpes_viewer.ui.viewers.widgets import MemoryData
    rng = np.random.default_rng(3)
    cube = rng.random((DEFL.size, SLIT.size, EKIN.size)) + 1
    return MemoryData(kind, (DEFL, SLIT, EKIN), cube,
                      {"x": "Angle, deflector (deg)", "k": "Angle, slit (deg)",
                       "z": "Kinetic energy (eV)"}, source_label="m")


def test_the_normalise_tab_starts_on_each_edc_of_a_cut(app):
    from arpes_viewer.ui.viewers.windows import DataOperationsDialog
    dialog = DataOperationsDialog([("c", _cut())])
    assert [c.isChecked() for c in dialog.norm_checks] == [False, True]
    assert "EDC" in dialog.norm_preset.currentText()
    assert "Kinetic energy" in dialog.norm_explain.text()
    dialog.tabs.setCurrentIndex(1)
    created = dialog.apply()
    np.testing.assert_allclose(created[0].scan.value.sum(axis=1), 1.0)
    dialog.close()


def test_the_normalise_tab_starts_on_each_cut_of_a_map(app):
    from arpes_viewer.ui.viewers.windows import DataOperationsDialog
    for kind in ("map", "kz_map"):
        dialog = DataOperationsDialog([("m", _map(kind))])
        assert [c.isChecked() for c in dialog.norm_checks] == [False, True, True]
        # The window goes on the energy, the usual thing to restrict.
        assert dialog.norm_window_axis.currentData() == 2
        dialog.norm_checks[1].setChecked(False)               # by hand: EDCs
        assert "EDC" in dialog.norm_preset.currentText()
        dialog.norm_checks[0].setChecked(True)                # matches nothing
        assert dialog.norm_preset.currentText().startswith("custom")
        dialog.close()


def test_the_3d_panel_starts_with_a_smoothing_that_does_something(app):
    from arpes_viewer.ui.panels.volume import VolumeProcessDialog
    dialog = VolumeProcessDialog([("m", _map())])
    assert dialog.plane_width_u.value() > 0 and dialog.plane_width_v.value() > 0
    texts = [dialog.plane_operation.itemText(i)
             for i in range(dialog.plane_operation.count())]
    assert "Second derivative along Angle, deflector" in texts
    assert "Normalise each line along Angle, slit" in texts
    assert not any("first axis" in t for t in texts)
    dialog.plane_normal.setCurrentIndex(0)                   # slit–energy cuts
    assert "Curvature 1D along Kinetic energy" in [
        dialog.plane_operation.itemText(i)
        for i in range(dialog.plane_operation.count())]
    dialog.close()


def test_the_2d_panel_mirrors_k_and_finds_the_fermi_level(app):
    from arpes_viewer.ui.panels.process import ProcessDialog
    dialog = ProcessDialog([("c", _cut())])
    assert dialog.sym_angle.value() == 90.0
    assert dialog.sym_cx.value() == 0.0
    assert dialog.fd_ef.value() == pytest.approx(EF, abs=0.03)
    dialog.close()


def test_the_map_conversion_keeps_the_measured_angles(app):
    from arpes_viewer.ui.viewers import windows as W
    window = W.ContourWindow(_map(), "m", "gray", False)
    dialog = W.KConversionDialog(window)
    assert dialog.nkx_box.value() == DEFL.size
    assert dialog.nky_box.value() == SLIT.size
    from PyQt5.QtWidgets import QLabel
    texts = [label.text() for label in dialog.findChildren(QLabel)]
    assert "Γ deflector angle (deg)" in texts and "kx points" in texts
    assert not any("offset" in t.lower() and "Theta" in t for t in texts)
    dialog.close()
    window.close()


def test_the_map_conversion_stores_float16_unless_asked_otherwise(app):
    from arpes_viewer.ui.viewers import windows as W
    window = W.ContourWindow(_map(), "m", "gray", False)
    dialog = W.KConversionDialog(window)
    from PyQt5.QtWidgets import QLabel
    texts = " ".join(label.text() for label in dialog.findChildren(QLabel))
    assert "float16 keeps about 3 significant digits" in texts
    settings = dialog.settings()
    assert settings["precision"] == "float16"
    kdata = window.run_k_conversion(settings, name="k", background=False)
    assert kdata.scan.value.dtype == np.float16
    dialog.precision_combo.setCurrentIndex(1)
    settings = dialog.settings()
    assert settings["precision"] == "float32"
    assert window.run_k_conversion(settings, name="k32",
                                   background=False).scan.value.dtype == np.float32
    dialog.close()
    window.close()
