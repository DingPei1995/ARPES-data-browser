"""The launcher's layout: the tool panel beside the list, the buttons under
it, the right-click menu, and buttons that light up only for a selection
they can take."""
import numpy as np
import pytest
from PyQt5.QtCore import QPoint, QTimer
from PyQt5.QtWidgets import QApplication, QDialogButtonBox, QLabel, QPushButton

from ui.viewers.widgets import MemoryData

SLIT = np.linspace(-10, 10, 21)
EKIN = np.linspace(54, 56, 31)
DEFL = np.linspace(-5, 5, 11)


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


@pytest.fixture
def launcher(app, monkeypatch):
    import ARPES_viewer as av
    from ui.launcher import main_window as main_window_ui
    monkeypatch.setattr(av, "loaded_items", {})
    win = av._MainWindow()
    ui = main_window_ui.Ui_MainWindow()
    ui.setupUi(win)
    monkeypatch.setattr(av, "win", win)
    monkeypatch.setattr(av, "ui", ui)
    monkeypatch.setattr(av, "autosave", lambda data, name: None)
    av.SetConnect()
    cut = MemoryData("cut", (SLIT, EKIN), np.ones((SLIT.size, EKIN.size)),
                     {"x": "Angle (deg)", "y": "Kinetic energy (eV)"},
                     source_label="c")
    cube = MemoryData("map", (DEFL, SLIT, EKIN),
                      np.ones((DEFL.size, SLIT.size, EKIN.size)),
                      {"x": "Angle, deflector (deg)", "k": "Angle, slit (deg)",
                       "z": "Kinetic energy (eV)"}, source_label="m")
    av.add_memory_dataset(cut, "c")
    av.add_memory_dataset(cube, "m")
    yield av
    win.hide()
    win.deleteLater()


def _select(lv, *texts):
    lw = lv.ui.FilePathListWidget
    lw.clearSelection()
    for i in range(lw.count()):
        if any(lw.item(i).text().startswith(t + " ") for t in texts):
            lw.item(i).setSelected(True)


def test_tool_buttons_follow_the_selection(launcher):
    lv = launcher
    _select(lv)
    assert not any(b.isEnabled() for b in lv.ui.ToolButtons.values())
    assert "Select a dataset" in lv.ui.ToolButtons["smooth"].toolTip()
    _select(lv, "c")
    tools = lv.ui.ToolButtons
    assert tools["edcfit"].isEnabled() and tools["background"].isEnabled()
    assert not tools["mdcfit"].isEnabled()            # an angle cut
    assert "momentum" in tools["mdcfit"].toolTip()
    assert not tools["slices3d"].isEnabled()
    _select(lv, "m")
    assert tools["slices3d"].isEnabled() and tools["notch3d"].isEnabled()
    assert tools["projection3d"].isEnabled() and tools["smooth"].isEnabled()
    assert not tools["stack"].isEnabled()
    assert "a cut" in tools["stack"].toolTip()
    _select(lv, "c", "m")
    assert not tools["smooth"].isEnabled()
    assert lv.ui.ListButtons["truncate"].isEnabled()
    assert not lv.ui.ListButtons["combine"].isEnabled()   # two kinds


def test_memory_controls_are_in_the_tool_panel(launcher):
    panel = launcher.ui.ToolPanel
    for button in (launcher.ui.FreeMemoryButton, launcher.ui.CloseViewersButton,
                   launcher.ui.OpenAutosaveButton):
        assert panel.isAncestorOf(button)
    assert not panel.isAncestorOf(launcher.ui.AboutLabel)
    assert list(launcher.ui.ListButtons) == ["combine", "truncate", "compress",
                                             "normalise"]


def test_the_right_click_menu_holds_only_list_actions_in_order(launcher):
    lv = launcher
    _select(lv, "m")
    QTimer.singleShot(0, lambda: lv.last_list_menu.close())
    lv.show_list_menu(QPoint(5, 5))
    texts = [a.text() or "---" for a in lv.last_list_menu.actions()]
    assert texts == ["Open", "Open slit cut", "Open scan-axis cut", "---",
                     "Combine datasets...", "Save dataset...",
                     "Remove from list", "---", "Rename...",
                     "Show information", "Session log..."]


def test_each_direct_operation_opens_alone(launcher):
    lv = launcher
    _select(lv, "c")
    lv.ui.ListButtons["compress"].click()
    dialog = lv.data_ops_dialog
    assert dialog.windowTitle().startswith("Compress")
    assert dialog.tabs.count() == 1 and dialog.tabs.tabBar().isHidden()
    dialog.compress_boxes[1].setValue(2)
    created = dialog.apply()
    assert created[0].scan.value.shape == (SLIT.size, EKIN.size // 2)
    dialog.close()


def test_a_filter_button_opens_its_own_tab(launcher):
    lv = launcher
    _select(lv, "c")
    lv.ui.ToolButtons["curvature"].click()
    assert lv.process_dialog.tabs.tabText(
        lv.process_dialog.tabs.currentIndex()) == "Curvature"
    lv.process_dialog.close()
    _select(lv, "m")
    lv.ui.ToolButtons["derivative"].click()
    dialog = lv.process_dialog
    assert "Second derivative" in dialog.plane_operation.currentText()
    assert not any(b.text() == "Open the 3D view"
                   for b in dialog.findChildren(QPushButton))
    dialog.close()


def test_close_ends_every_button_row(app):
    from ui.viewers.widgets import ButtonRow, add_button
    row = ButtonRow()
    add_button(row, "Close", QDialogButtonBox.RejectRole)
    row = ButtonRow()
    apply_ = add_button(row, "Apply", QDialogButtonBox.ApplyRole)
    close = add_button(row, "Close", QDialogButtonBox.RejectRole)
    layout = row.layout()
    order = [layout.itemAt(i).widget() for i in range(layout.count())]
    assert order.index(apply_) < order.index(close) == len(order) - 1
    rejected = []
    row.rejected.connect(lambda: rejected.append(1))
    close.click()
    assert rejected == [1]


def test_a_k_cut_lights_the_mdc_fit(launcher):
    lv = launcher
    k = np.linspace(-0.5, 0.5, 21)
    kcut = MemoryData("cut", (k, EKIN), np.ones((k.size, EKIN.size)),
                      {"x": "k (Å⁻¹)", "y": "Kinetic energy (eV)"},
                      source_label="k")
    lv.add_memory_dataset(kcut, "k")
    _select(lv, "k")
    assert lv.ui.ToolButtons["mdcfit"].isEnabled()
    lv.ui.ToolButtons["mdcfit"].click()
    window = lv.viewer_windows[-1]
    assert window.windowTitle().startswith("MDC fit (k-cut only)")
    assert window.direction.currentIndex() == 0 and window.fermi_group.isHidden()
    window.close()
    _select(lv, "c")
    lv.ui.ToolButtons["edcfit"].click()
    window = lv.viewer_windows[-1]
    assert window.windowTitle().startswith("EDC fit")
    assert window.direction.currentIndex() == 1
    assert not window.dispersion_button.isEnabled()    # angle cut: no v_F
    window.close()


@pytest.mark.parametrize("operation, title", [
    ("smooth", "Smooth"), ("derivative", "Second derivative"),
    ("curvature", "Curvature"), ("background", "Background subtraction"),
    ("symmetry", "Symmetrise")])
def test_each_filter_opens_alone(launcher, operation, title):
    lv = launcher
    _select(lv, "c")
    lv.ui.ToolButtons[operation].click()
    dialog = lv.process_dialog
    assert dialog.windowTitle() == f"{title} (2D)"
    assert dialog.tabs.count() == 1 and dialog.tabs.tabBar().isHidden()
    tip = "smoothed first" in " ".join(
        label.text() for label in dialog.findChildren(QLabel))
    assert tip == (operation in ("derivative", "curvature"))
    dialog.close()
    _select(lv, "m")
    if operation == "background":
        assert not lv.ui.ToolButtons[operation].isEnabled()
        return
    lv.ui.ToolButtons[operation].click()
    dialog = lv.process_dialog
    assert dialog.windowTitle() == f"{title} (3D)" and dialog.tabs.count() == 1
    if operation != "symmetry":
        texts = [dialog.plane_operation.itemText(i)
                 for i in range(dialog.plane_operation.count())]
        word = {"smooth": ("smooth", "spikes"), "derivative": ("derivative",),
                "curvature": ("Curvature",)}[operation]
        assert all(any(w in t for w in word) for t in texts)
        dialog.refresh()
        assert "<span style='color:#b00'>" not in dialog.note.text()
    dialog.close()
