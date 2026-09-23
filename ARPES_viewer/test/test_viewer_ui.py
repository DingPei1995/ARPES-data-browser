"""The viewers' layout and the map's shared readout cursor.

* Tools are in a Functions menu; the window itself keeps a map's two cut
  buttons, the colormap and the view controls, on one row.
* One readout cursor for a contour and its cut windows (ui.cursorlink):
  switched on and off together, moved together, integration widths kept
  (and drawn) only in the window they are set in, typed positions, and a
  refusal to leave the data.

Needs Qt; run headless with QT_QPA_PLATFORM=offscreen.
"""
import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PyQt5")
from PyQt5.QtWidgets import QApplication, QPushButton, QToolButton   # noqa: E402

from ui.widgets import MemoryData, AXIS_COLORS                     # noqa: E402
from ui import windows as W                                        # noqa: E402

DEFL = np.linspace(-10, 10, 21)
SLIT = np.linspace(-15, 15, 31)
ENERGY = np.linspace(-1.0, 0.2, 61)
LABELS = {"x": "Angle, deflector (deg)", "k": "Angle, slit (deg)", "z": "Energy (eV)"}


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


def _map(kind="map"):
    rng = np.random.default_rng(3)
    return MemoryData(kind, (DEFL, SLIT, ENERGY),
                      rng.random((DEFL.size, SLIT.size, ENERGY.size)) + 1,
                      LABELS, source_label="m")


def _cut():
    rng = np.random.default_rng(4)
    return MemoryData("cut", (SLIT, ENERGY), rng.random((SLIT.size, ENERGY.size)),
                      {"x": "Angle, slit (deg)", "y": "Energy (eV)"}, source_label="c")


@pytest.fixture()
def linked(app):
    window = W.ContourWindow(_map(), "m", "gray", False)
    window.show()
    slit = window.open_cut("slit")
    defl = window.open_cut("deflector")
    app.processEvents()
    yield window, slit, defl
    window.close()


def _visible_texts(window):
    return [a.text() for a in window.function_actions() if a.isVisible()]


def _toolbar_buttons(window):
    """Push buttons on the window's top row (the Functions button is a
    QToolButton and is counted separately)."""
    bar = window.colormap_combo.parentWidget()
    return [b.text() for b in bar.findChildren(QPushButton)
            if b.parentWidget() is bar]


# -- layout -----------------------------------------------------------------
def test_map_keeps_only_its_cut_buttons(app):
    window = W.ContourWindow(_map(), "m", "gray", False)
    assert _toolbar_buttons(window) == ["Deflector cut", "Slit cut"]
    texts = _visible_texts(window)
    for expected in ("Arbitrary cut...", "Map k conversion...", "De-grid map...",
                     "Slice figure...", "Save slice to the main list",
                     "Open slice in a new panel"):
        assert expected in texts
    assert "Brillouin zone..." not in texts          # angle map: not offered
    assert "kz map processing..." not in texts
    window.close()


def test_colormap_and_view_controls_share_one_row(app):
    window = W.CutWindow(_cut(), "c", "gray", False)
    bar = window.colormap_combo.parentWidget()
    assert window.view_bar.parentWidget() is bar
    assert window.functions_button.parentWidget() is bar
    assert isinstance(window.functions_button, QToolButton)
    window.close()


def test_cut_viewer_keeps_no_tool_buttons(app):
    window = W.CutWindow(_cut(), "c", "gray", False)
    assert _toolbar_buttons(window) == []
    texts = _visible_texts(window)
    for expected in ("Fermi level...", "MDC / EDC fit...", "FS correction...",
                     "Cut k conversion...", "Cut arithmetic...", "De-grid..."):
        assert expected in texts
    # an angle cut cannot be band-fitted: the entry says why instead of vanishing
    assert not window.fit_action.isEnabled()
    assert window.fit_action.toolTip().startswith("Not available here")
    window.close()


def test_menu_sections_and_hidden_entries(app):
    window = W.ContourWindow(_map("kz_map"), "kz", "gray", False)
    window._rebuild_functions_menu()
    shown = [a.text() for a in window.functions_menu.actions()]
    assert shown[0] == "Data operations"
    assert "kz map processing..." in shown and "kz -> momentum..." in shown
    assert "Map k conversion..." not in shown
    assert shown.index("Slice") > shown.index("Visualization")
    k_map = W.ContourWindow(_map("k_map"), "k", "gray", False)
    assert "Brillouin zone..." in _visible_texts(k_map)
    assert "De-grid map..." not in _visible_texts(k_map)
    window.close()
    k_map.close()


def test_a_menu_entry_runs_its_tool(app):
    window = W.ContourWindow(_map(), "m", "gray", False)
    window.arbcut_action.trigger()
    assert getattr(window, "_arbcut_dialog", None) is not None or \
        any(type(c).__name__ == "ArbitraryCutDialog" for c in window.children())
    window.close()


# -- one cursor ----------------------------------------------------------------
def test_switching_on_in_the_contour_switches_on_everywhere(linked):
    contour, slit, defl = linked
    contour.view.set_readout_cursor_visible(True)
    assert all(w.view.readout_cursor is not None for w in (contour, slit, defl))
    assert slit.panel.cursor_widget.isVisibleTo(slit)
    slit.view.set_readout_cursor_visible(False)          # ...and off from a cut
    assert all(w.view.readout_cursor is None for w in (contour, slit, defl))


def test_switching_on_in_a_cut_switches_on_the_contour(linked):
    contour, slit, defl = linked
    defl.view.set_readout_cursor_visible(True)
    assert contour.view.readout_cursor is not None
    assert slit.view.readout_cursor is not None


def test_moving_on_the_contour_reslices_the_cuts(linked):
    contour, slit, defl = linked
    contour.view.set_readout_cursor_visible(True)
    assert contour.view.move_readout_to(x=4.0, y=-6.0) is None
    contour.cursor_link.flush()
    d, s = int(np.argmin(abs(DEFL - 4.0))), int(np.argmin(abs(SLIT + 6.0)))
    e = contour.e_control.index
    assert slit.control.index == d                     # slit cut at this deflector
    assert defl.control.index == s                     # deflector cut at this slit
    assert slit.view.readout_indices() == (s, e)
    assert defl.view.readout_indices() == (d, e)
    np.testing.assert_allclose(slit.view.last_frame, contour.cube[d])


def test_moving_in_energy_on_a_cut_moves_the_contour(linked):
    contour, slit, defl = linked
    slit.view.set_readout_cursor_visible(True)
    before = contour.view.readout_indices()
    assert slit.view.move_readout_to(y=-0.5) is None
    contour.cursor_link.flush()
    e = int(np.argmin(abs(ENERGY + 0.5)))
    assert contour.e_control.index == e
    assert defl.view.readout_indices()[1] == e
    assert contour.view.readout_indices() == before    # its in-plane point is unchanged


def test_the_energy_slider_moves_the_cuts_cursors(linked):
    contour, slit, defl = linked
    contour.view.set_readout_cursor_visible(True)
    contour.e_control.set_index(12)
    contour.cursor_link.flush()
    assert slit.view.readout_indices()[1] == 12
    assert defl.view.readout_indices()[1] == 12


def test_typing_a_position(linked):
    contour, slit, defl = linked
    contour.view.set_readout_cursor_visible(True)
    edit = slit.panel.cursor_edits[0]
    edit.setText("7")
    assert slit.panel._on_cursor_typed(0) is None
    contour.cursor_link.flush()
    s = int(np.argmin(abs(SLIT - 7.0)))
    assert contour.view.readout_indices()[1] == s
    assert defl.control.index == s


def test_a_position_outside_the_data_is_refused(linked):
    contour, slit, defl = linked
    contour.view.set_readout_cursor_visible(True)
    where = slit.view.readout_position()
    slit.panel.cursor_edits[1].setText("5")          # energy runs to 0.2 eV
    message = slit.panel._on_cursor_typed(1)
    assert message and "outside the data" in message
    assert slit.view.readout_position() == where
    assert "outside" in slit.panel.cursor_message.text()
    slit.panel.cursor_edits[0].setText("abc")
    assert "not a number" in slit.panel._on_cursor_typed(0)


def test_integration_widths_are_independent(linked):
    contour, slit, defl = linked
    contour.view.set_readout_cursor_visible(True)
    contour.e_control.width_spin.setValue(0.1)
    slit.control.width_spin.setValue(2.0)
    defl.panel.curves.edc_width.setValue(3.0)
    slit.panel.curves.mdc_width.setValue(0.05)
    assert slit.panel.curves.mdc_width.value() == pytest.approx(0.05)
    assert defl.panel.curves.mdc_width.value() == pytest.approx(0.0)
    assert contour.e_control.width_spin.value() == pytest.approx(0.1)
    assert slit.control.width_spin.value() == pytest.approx(2.0)
    assert defl.panel.curves.edc_width.value() == pytest.approx(3.0)
    assert slit.panel.curves.edc_width.value() == pytest.approx(0.0)


def test_widths_are_drawn_only_where_they_are_set(linked):
    contour, slit, defl = linked
    contour.view.set_readout_cursor_visible(True)
    contour.e_control.width_spin.setValue(0.1)       # the contour's energy window
    slit.control.width_spin.setValue(2.0)            # slit cut: over the deflector
    defl.control.width_spin.setValue(4.0)            # deflector cut: over the slit
    slit.panel.curves.edc_width.setValue(1.0)        # the slit cut's own EDC window
    slit.panel.curves.mdc_width.setValue(0.05)
    cv, sv, dv = contour.view, slit.view, defl.view
    # the contour's cursor shows no width at all
    assert not cv.readout_vband.isVisible() and not cv.readout_hband.isVisible()
    # a cut shows its own EDC and MDC windows, and nothing from the others
    assert sv.band_half_width(sv.readout_vband) == pytest.approx(1.0)
    assert sv.band_half_width(sv.readout_hband) == pytest.approx(0.05)
    assert not dv.readout_vband.isVisible() and not dv.readout_hband.isVisible()
    lo, hi = sv.readout_vband.getRegion()
    assert hi - lo == pytest.approx(2.0)


def test_axis_colours(linked):
    contour, slit, defl = linked
    contour.view.set_readout_cursor_visible(True)

    def colour(line):
        return line.pen.color().name()

    assert colour(contour.view.readout_vline) == AXIS_COLORS["defl"]
    assert colour(contour.view.readout_hline) == AXIS_COLORS["slit"]
    assert colour(slit.view.readout_vline) == AXIS_COLORS["slit"]
    assert colour(slit.view.readout_hline) == AXIS_COLORS["energy"]
    assert colour(defl.view.readout_vline) == AXIS_COLORS["defl"]
    assert slit.panel.curves.edc_curve.opts["pen"].color().name() == AXIS_COLORS["slit"]


def test_a_cut_opened_later_joins_in(app):
    contour = W.ContourWindow(_map(), "m", "gray", False)
    contour.show()
    contour.view.set_readout_cursor_visible(True)
    contour.view.move_readout_to(x=-3.0, y=2.0)
    slit = contour.open_cut("slit")
    assert slit.view.readout_cursor is not None
    d = int(np.argmin(abs(DEFL + 3.0)))
    assert slit.control.index == d
    slit.close()
    assert len(contour.cursor_link.members) == 1
    contour.view.move_readout_to(x=3.0)            # no error with the cut gone
    contour.cursor_link.flush()
    contour.close()


def test_a_lone_cut_has_the_cursor_row_too(app):
    window = W.CutWindow(_cut(), "c", "gray", False)
    window.view.set_readout_cursor_visible(True)
    window.panel.cursor_edits[0].setText("3")
    assert window.panel._on_cursor_typed(0) is None
    assert window.view.readout_position()[0] == pytest.approx(3.0)
    window.panel.curves.edc_width.setValue(1.0)
    assert window.view.readout_vband.isVisible()
    window.close()


# -- the curve viewer ---------------------------------------------------------------
def _curve(kind="edc", channels=1):
    rng = np.random.default_rng(5)
    x = np.linspace(-1, 0.2, 80)
    values = rng.poisson(100, size=(80, channels)).astype(float)
    info = {}
    if kind == "spin_edc":
        info["curve.channels"] = "|".join(f"C{i}" for i in range(channels))
    return MemoryData(kind, (x, np.arange(channels, dtype=float)), values,
                      {"x": "Energy (eV)"}, source_label="curve",
                      source_info=info)


def test_curve_viewer_has_a_functions_menu(app):
    from ui.curves import CurveWindow
    window = CurveWindow(_curve(), "edc")
    bar_buttons = [b.text() for b in window.findChildren(QPushButton)
                   if b.isVisibleTo(window)]
    assert "Curve fit..." not in bar_buttons and "As a figure..." not in bar_buttons
    texts = _visible_texts(window)
    for expected in ("Curve fit...", "As a figure...", "Bin...", "Normalise...",
                     "Crop to the range...", "Add counting errors (√N)"):
        assert expected in texts
    assert "Spin analysis..." not in texts               # not a spin EDC
    window._rebuild_functions_menu()
    shown = [a.text() for a in window.functions_menu.actions()]
    assert shown.index("Analysis") < shown.index("Data operations") < shown.index("Visualization")
    window.close()


def test_an_operation_from_the_menu_opens_the_strip_on_it(app):
    from ui.curves import CurveWindow
    window = CurveWindow(_curve(), "edc")
    assert not window.operations_box.isVisibleTo(window)
    window.operation_actions["bin"].trigger()
    assert window.operations_box.isVisibleTo(window)
    assert window.operation.currentData() == "bin"
    window.close()


def test_spin_analysis_is_offered_for_a_spin_edc(app):
    from ui.curves import CurveWindow
    window = CurveWindow(_curve("spin_edc", 4), "spin")
    assert window.spin_action.isVisible()
    window.close()


def test_an_unseen_contour_closes_with_its_last_cut(app):
    """A cut opened straight from the main list keeps its contour off
    screen; the contour goes once the last cut has, unless it was shown."""
    contour = W.ContourWindow(_map(), "m", "gray", False)     # never shown
    closed = []
    contour.closed.connect(lambda *_: closed.append(True))
    slit = contour.open_cut("slit")
    defl = contour.open_cut("deflector")
    slit.close()
    assert not closed                          # the deflector cut still needs it
    defl.close()
    assert closed
    shown = W.ContourWindow(_map(), "m", "gray", False)
    cut = shown.open_cut("slit")
    cut.show_map_action.trigger()             # "Show the map"
    assert shown.isVisible()
    cut.close()
    assert shown.isVisible() and shown.cursor_link.members
    shown.close()
