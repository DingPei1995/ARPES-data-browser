"""Fixes from running the program end to end for the user guides.

* A cut whose energy axis was shifted to E - E_F (Fermi level -> Offset
  energy axis) converts to k: the conversion adds E_F back by default and
  the result stays on E - E_F. It used to fail with "the cut's k axis is
  not monotonic", blaming the angles.
* An image panel with EDC/MDC curves gives the whole panel to the image
  while the readout cursor is off (it used to keep an empty strip right of
  and under it).
* The readout label flips to the inside of the plot near its top or right
  edge instead of being cut off.
* The processing previews pick "Equal axes" for a plane whose two axes
  share a unit (angle vs angle, k vs k) until the user picks a shape.
* The MDC/EDC fit panel builds without a signal reaching a timer that did
  not exist yet.
* A listed dataset whose entry has the row's own name is not called
  "name · name".

Needs Qt; run headless with QT_QPA_PLATFORM=offscreen.
"""
import os
import sys

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PyQt5")
from PyQt5.QtWidgets import QApplication                      # noqa: E402

from ui.widgets import MemoryData                             # noqa: E402
from ui import windows as W                                   # noqa: E402

SLIT = np.linspace(-12, 12, 61)
EKIN = np.linspace(54.6, 55.7, 45)
EF = 55.6


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


def _cut(energy=EKIN, label="Kinetic energy (eV)", info=None):
    rng = np.random.default_rng(7)
    return MemoryData("cut", (SLIT, energy), rng.random((SLIT.size, energy.size)) + 1,
                      {"x": "Angle, slit (deg)", "y": label}, source_label="c",
                      source_info=info)


def _shifted_cut():
    return _cut(EKIN - EF, "E - E_F (eV)", {"fitEF.ef_subtracted": EF})


def test_subtracted_fermi_level_is_read_back(app):
    assert W.subtracted_fermi_level(_shifted_cut()) == pytest.approx(EF)
    assert W.subtracted_fermi_level(_cut()) == 0.0


def test_a_shifted_cut_converts_and_stays_on_e_minus_ef(app):
    window = W.CutWindow(_shifted_cut(), "c [E-Ef]", "gray", False)
    made = []
    window.datasetCreated.connect(made.append)
    dialog = W.CutKConversionDialog(window)
    assert dialog.energy_offset_box.value() == pytest.approx(EF, abs=1e-4)
    assert "E − E_F" in dialog.file_note.text()
    result = window.run_cut_k_conversion(dialog.settings(), name="k")
    assert result is not None and made
    np.testing.assert_allclose(result.scan.y, EKIN - EF, atol=1e-9)
    assert result.scan.labels["y"] == "E - E_F (eV)"
    assert np.isfinite(result.scan.value).any()
    dialog.close()
    window.close()


def test_a_kinetic_cut_keeps_a_zero_offset(app):
    window = W.CutWindow(_cut(), "c", "gray", False)
    dialog = W.CutKConversionDialog(window)
    assert dialog.energy_offset_box.value() == 0.0
    dialog.close()
    window.close()


def test_the_image_takes_the_curves_space_while_the_cursor_is_off(app):
    window = W.CutWindow(_cut(), "c", "gray", False)
    window.resize(900, 700)
    window.show()
    app.processEvents()
    grid = window.panel._plot_grid
    assert grid.columnStretch(1) == 0 and grid.rowStretch(1) == 0
    window.view.set_readout_cursor_visible(True)
    app.processEvents()
    assert grid.columnStretch(1) == 1 and grid.rowStretch(1) == 1
    window.view.set_readout_cursor_visible(False)
    assert grid.columnStretch(1) == 0 and grid.rowStretch(1) == 0
    window.close()


def test_the_readout_label_stays_inside_the_plot(app):
    window = W.CutWindow(_cut(), "c", "gray", False)
    window.show()
    app.processEvents()
    view = window.view
    view.set_readout_cursor_visible(True)
    view.move_readout_to(SLIT[1], EKIN[2])              # lower left: up and right
    assert tuple(view.readout_label.anchor) == (0, 1)
    view.move_readout_to(SLIT[-2], EKIN[-2])            # upper right: down and left
    assert tuple(view.readout_label.anchor) == (1, 0)
    window.close()


def test_previews_choose_equal_axes_for_matching_units(app):
    from ui.process import PreviewPair
    pair = PreviewPair()
    pair.set_labels("Angle, deflector (deg)", "Angle, slit (deg)")
    assert pair.aspect.currentText() == "Equal axes"
    pair.set_labels("k (Å⁻¹)", "Energy (eV)")
    assert pair.aspect.currentText() == "Square pixels"
    pair.aspect.setCurrentText("Fill the box")
    pair.aspect.activated.emit(pair.aspect.currentIndex())   # as a click does
    pair.set_labels("kx (Å⁻¹)", "ky (Å⁻¹)")
    assert pair.aspect.currentText() == "Fill the box"      # the user's choice


def test_the_fit_panel_builds_cleanly(app, monkeypatch):
    from ui.fit import FitPanel
    errors = []
    monkeypatch.setattr(sys, "excepthook", lambda *exc: errors.append(exc))
    k = np.linspace(-1, 1, 61)
    data = MemoryData("cut", (k, EKIN - EF), np.ones((61, EKIN.size)),
                      {"x": "k (Å⁻¹)", "y": "E - E_F (eV)"},
                      source_label="k")
    panel = FitPanel(data, "k")
    app.processEvents()
    assert not errors
    panel.close()


def test_no_name_repeated_in_a_row_label():
    import ARPES_viewer as av
    saved = dict(av.loaded_items)
    try:
        av.loaded_items.clear()
        av.loaded_items[("/d/cuts.nxs", "cut_LH")] = {"name": "cut_LH", "entry": "cut_LH"}
        av.loaded_items[("/d/cuts.nxs", "cut_LV")] = {"name": "cuts.nxs", "entry": "cut_LV"}
        assert av.label_for(("/d/cuts.nxs", "cut_LH")) == "cut_LH"
        assert av.label_for(("/d/cuts.nxs", "cut_LV")) == "cuts.nxs · cut_LV"
    finally:
        av.loaded_items.clear()
        av.loaded_items.update(saved)


@pytest.mark.parametrize("lo, hi, target", [(-15, 15, 2), (0.31, 0.38, 2),
                                            (54.6, 55.72, 2), (-1, 1, 1)])
def test_an_axis_always_gets_at_least_two_tick_labels(lo, hi, target):
    from tools.export import tick_values
    values, step = tick_values(lo, hi, target)
    assert values.size >= 2
    assert np.all((values >= lo - 1e-9) & (values <= hi + 1e-9))
    np.testing.assert_allclose(np.diff(values), step)


def test_a_wide_x_label_is_moved_in_rather_than_cut_off():
    from ui.figure import FigurePainter
    centre = FigurePainter.label_centre
    assert centre(50.0, 20.0, 100.0) == 50.0              # fits: left alone
    assert centre(95.0, 20.0, 100.0) == 90.0              # right edge: moved in
    assert centre(4.0, 20.0, 100.0) == 10.0               # left edge: moved in
    assert centre(50.0, 120.0, 100.0) == 50.0             # wider than the page


def test_a_kz_map_cut_names_its_slider_in_ev(app):
    hv = np.linspace(30, 60, 5)
    data = MemoryData("kz_map", (hv, SLIT, EKIN - EF), np.ones((5, SLIT.size, EKIN.size)),
                      {"x": "Photon energy (eV)", "k": "Angle, slit (deg)",
                       "z": "E - E_F (eV)"}, source_label="kz")
    contour = W.ContourWindow(data, "kz", "gray", False)
    cut = contour.open_cut("slit")
    app.processEvents()
    assert cut.control.unit == "eV"
    assert "deg" not in cut.control.readout.text()
    contour.close()
