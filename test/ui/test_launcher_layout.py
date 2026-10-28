"""The launcher's layout: the tool panel beside the list, the buttons under
it, the right-click menu, and buttons that light up only for a selection
they can take."""
import numpy as np
import pytest
from PyQt5.QtCore import QPoint, QTimer
from PyQt5.QtWidgets import (QApplication, QDialogButtonBox, QLabel,
                             QMessageBox, QPushButton)

from ui.viewers import windows as W
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
        if any(lw.item(i).text().startswith(t + "   [") for t in texts):
            lw.item(i).setSelected(True)


def test_tool_buttons_follow_the_selection(launcher):
    lv = launcher
    _select(lv)
    # all but the figure composer, which opens a blank page
    assert not any(b.isEnabled() for name, b in lv.ui.ToolButtons.items()
                   if name != "figure")
    assert lv.ui.ToolButtons["figure"].isEnabled()
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


def test_the_right_click_menu_holds_only_the_ways_of_opening(launcher):
    lv = launcher
    _select(lv, "m")
    QTimer.singleShot(0, lambda: lv.last_list_menu.close())
    lv.show_list_menu(QPoint(5, 5))
    texts = [a.text() or "---" for a in lv.last_list_menu.actions()]
    assert texts == ["Open", "Open slit cut", "Open scan-axis cut"]


def test_the_rest_of_the_menu_is_the_datasets_box_above_memory(launcher):
    ui = launcher.ui
    panel = ui.ToolPanel
    assert [b.text() for b in ui.DatasetButtons.values()] == [
        "Experiment info", "Log format...", "Auto Exp. log", "Rename...",
        "Save datasets...", "Remove from list", "Session log..."]
    # the log's format and the log itself, side by side
    fmt, log = ui.DatasetButtons["explogformat"], ui.DatasetButtons["explog"]
    rows = [ui.DatasetBox.layout().itemAt(i).layout()
            for i in range(ui.DatasetBox.layout().count())]
    pair = next(r for r in rows if r.indexOf(fmt) >= 0)
    assert pair.indexOf(fmt) == 0 and pair.indexOf(log) == 1
    layout = panel.layout()
    order = [layout.indexOf(w) for w in (ui.ToolTabs, ui.DatasetBox,
                                         ui.MemoryBox)]
    assert order == sorted(order) and -1 not in order
    # combining is under the list only
    assert "combine" not in ui.DatasetButtons and "combine" in ui.ListButtons
    _select(launcher)
    assert not ui.DatasetButtons["info"].isEnabled()
    assert ui.DatasetButtons["explogformat"].isEnabled()
    assert not ui.DatasetButtons["explog"].isEnabled()
    _select(launcher, "c")
    assert ui.DatasetButtons["info"].isEnabled()
    _select(launcher, "c", "m")
    assert not ui.DatasetButtons["rename"].isEnabled()
    assert ui.DatasetButtons["save"].isEnabled()


def test_experiment_info_opens_the_information_window(launcher, monkeypatch):
    lv = launcher
    # the window belongs to this test's launcher, which goes with it
    monkeypatch.setattr(lv, "info_window", None)
    _select(lv, "c")
    lv.ui.DatasetButtons["info"].click()
    assert lv.info_window.isVisible()
    assert lv.info_window.windowTitle() == "Experiment info"
    lv.info_window.close()


def _set_up_log(lv, folder, name, template, keep=False, then=None):
    """Click Log format..., fill it in, and save it."""
    def fill():
        dialog = lv.last_exp_log_dialog
        dialog.folder_edit.setText(str(folder))
        dialog.name_edit.setText(name)
        dialog._sync_to_file()
        if template is not None:
            dialog.set_template(template)
        dialog.auto_box.setChecked(keep)
        if then is not None:
            then(dialog)
        dialog.accept()
    QTimer.singleShot(0, fill)
    lv.ui.DatasetButtons["explogformat"].click()


def test_the_log_is_set_up_once_then_one_click_logs(launcher, tmp_path,
                                                    monkeypatch):
    import csv
    from tools.output import explog as EL
    lv = launcher
    monkeypatch.setattr("ui.panels.explog.excel_conventions",
                        lambda: (",", "."))
    _select(lv, "m")

    def tick_scan(dialog):
        # the selected map is the reference; its fields are offered
        fields = [dialog.field_list.item(i).data(0x0100)
                  for i in range(dialog.field_list.count())]
        assert fields[:len(EL.DERIVED)] == list(EL.DERIVED)
        item = dialog.field_list.item(fields.index("@scan1"))
        item.setCheckState(2)
        assert dialog.template.columns[-1].field == "@scan1"
        dialog.table.selectRow(2)
        dialog.move_column(-5)
        assert [c.field for c in dialog.template.columns] == [
            "@scan1", "@file", "@kind"]
    _set_up_log(lv, tmp_path, "beamtime",
                EL.Template([EL.Column("File name", "@file"),
                             EL.Column("Type", "@kind")]), then=tick_scan)
    assert not (tmp_path / "beamtime.csv").exists()     # nothing logged yet
    # one click, no window
    lv.ui.DatasetButtons["explog"].click()
    _select(lv, "c", "m")
    lv.ui.DatasetButtons["explog"].click()              # "m" is in it
    with open(tmp_path / "beamtime.csv", encoding="utf-8-sig") as handle:
        rows = list(csv.reader(handle))
    assert rows == [["Scan axis", "File name", "Type"],
                    ["Angle, deflector (deg): -5/5/1", "m", "Map"],
                    ["", "c", "Cut"]]

    # the format changed later: the same file, still one header row
    def add_notes(dialog):
        assert dialog.template.headers() == rows[0]   # the log's own format
        dialog.add_blank()
        dialog.table.selectRow(3)
        dialog.move_column(-3)
    monkeypatch.setattr("ui.panels.explog.QMessageBox.question",
                        lambda *a, **k: QMessageBox.Yes)
    _set_up_log(lv, tmp_path, "beamtime", None, then=add_notes)
    with open(tmp_path / "beamtime.csv", encoding="utf-8-sig") as handle:
        rows = list(csv.reader(handle))
    assert rows == [["Comments", "Scan axis", "File name", "Type"],
                    ["", "Angle, deflector (deg): -5/5/1", "m", "Map"],
                    ["", "", "c", "Cut"]]


def test_the_first_click_asks_for_the_format(launcher, tmp_path, monkeypatch):
    from tools.output import explog as EL
    lv = launcher
    _select(lv, "c")
    asked = []
    monkeypatch.setattr(lv, "open_exp_log_format",
                        lambda: asked.append(1) or False)
    assert lv.record_exp_log() == 0 and asked == [1]
    EL.save_settings(EL.Settings(str(tmp_path), "x.csv", configured=True))
    assert lv.record_exp_log() == 1 and asked == [1]


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


def test_degrid_is_a_tool_for_cuts_and_angle_maps(launcher):
    lv = launcher
    _select(lv, "c")
    assert lv.ui.ToolButtons["degrid"].isEnabled()
    lv.ui.ToolButtons["degrid"].click()
    dialog = lv.degrid_dialog
    assert dialog.windowTitle().startswith("De-grid") and not dialog.is_map
    dialog.close()
    _select(lv, "m")
    lv.ui.ToolButtons["degrid"].click()
    assert lv.degrid_dialog.is_map
    lv.degrid_dialog.close()
    _select(lv, "c", "m")
    assert not lv.ui.ToolButtons["degrid"].isEnabled()


def test_degrid_says_why_not_for_a_map_in_momentum(app, monkeypatch):
    from PyQt5.QtWidgets import QMessageBox
    from ui.panels.degrid import open_degrid
    said = []
    monkeypatch.setattr(QMessageBox, "information",
                        staticmethod(lambda parent, title, text, *a: said.append(text)))
    kmap = MemoryData("k_map", (DEFL, SLIT, EKIN),
                      np.ones((DEFL.size, SLIT.size, EKIN.size)),
                      {"x": "kx (Å⁻¹)", "k": "ky (Å⁻¹)", "z": "E (eV)"},
                      source_label="k")
    assert open_degrid(None, kmap, "k") is None
    assert said and "momentum" in said[0]


def test_fermi_fit_and_kz_alignment_are_tools(launcher):
    lv = launcher
    _select(lv, "c")
    assert lv.ui.ToolButtons["fermi"].isEnabled()
    assert not lv.ui.ToolButtons["kzalign"].isEnabled()
    lv.ui.ToolButtons["fermi"].click()
    viewer = lv.viewer_windows[-1]
    dialog = viewer._fermi_dialog
    assert isinstance(viewer, W.CutWindow) and not dialog.slices
    dialog.close()
    QApplication.processEvents()
    assert not viewer.isVisible()                      # tied together

    _select(lv, "m")
    lv.ui.ToolButtons["fermi"].click()
    cut = next(w for w in QApplication.topLevelWidgets()
               if isinstance(w, W.MapCutWindow) and w.isVisible())
    assert cut._fermi_dialog.slices and not cut.contour.isVisible()
    cut.close()
    QApplication.processEvents()

    hv = np.linspace(30, 60, 7)
    rng = np.random.default_rng(2)
    kz = MemoryData("kz_map", (hv, SLIT, EKIN),
                    rng.random((hv.size, SLIT.size, EKIN.size)) + 1,
                    {"x": "Photon energy (eV)", "k": "Angle, slit (deg)",
                     "z": "Kinetic energy (eV)"}, source_label="kz")
    lv.add_memory_dataset(kz, "kz")
    _select(lv, "kz")
    assert lv.ui.ToolButtons["kzalign"].isEnabled()
    dialog = lv.open_kz_alignment()
    assert dialog is not None and dialog.cut_window.isVisible()
    dialog.close()
    QApplication.processEvents()
    assert not dialog.cut_window.isVisible()


def test_settings_are_remembered_per_size(tmp_path):
    from tools.runtime import remembered
    assert remembered.recall("t", (3, 4)) is None
    remembered.remember("t", (3, 4), {"a": 1.5})
    remembered.remember("t", (5, 4), {"a": 2.0})
    assert remembered.recall("t", (3, 4)) == {"a": 1.5}
    assert remembered.recall("t", (4, 3)) is None
    assert remembered.recall("other", (3, 4)) is None


def test_the_panel_has_a_basic_and_an_advanced_section(launcher):
    ui = launcher.ui
    basic, advanced = ui.ToolSections["Basic"], ui.ToolSections["Advanced"]
    for name in ("fermi", "degrid", "background", "kcut", "corelevel", "kmap",
                 "kzalign", "kzconv"):
        assert basic.isAncestorOf(ui.ToolButtons[name]), name
    for name in ("figure", "smooth", "derivative", "curvature", "symmetry",
                 "stack", "mdcfit", "edcfit", "arithmetic", "slices3d",
                 "notch3d", "projection3d"):
        assert advanced.isAncestorOf(ui.ToolButtons[name]), name


def test_map_k_conversion_opens_a_bare_contour_and_remembers(launcher, monkeypatch):
    lv = launcher
    run = W.ContourWindow.run_k_conversion
    monkeypatch.setattr(W.ContourWindow, "run_k_conversion",
                        lambda self, s, name=None, background=True:
                        run(self, s, name=name, background=False))
    _select(lv, "m")
    assert lv.ui.ToolButtons["kmap"].isEnabled()
    assert not lv.ui.ToolButtons["kzconv"].isEnabled()
    lv.ui.ToolButtons["kmap"].click()
    contour = lv.viewer_windows[-1]
    dialog = contour._kconv_dialog
    assert isinstance(contour, W.ContourWindow) and dialog is not None
    for hidden in (contour.functions_button, contour.plot_tools_button,
                   contour.slit_cut_button, contour.defl_cut_button):
        assert not hidden.isVisibleTo(contour)
    assert not dialog.restored
    dialog.theta_box.setValue(1.5)
    dialog.azimuth_box.setValue(10.0)
    rows = lv.ui.FilePathListWidget.count()
    dialog.accept()
    QApplication.processEvents()
    assert lv.ui.FilePathListWidget.count() == rows + 1      # the k-map listed
    assert not contour.isVisible()                            # closed with it

    _select(lv, "m")                                          # same size again
    lv.ui.ToolButtons["kmap"].click()
    dialog = lv.viewer_windows[-1]._kconv_dialog
    assert dialog.restored
    assert dialog.theta_box.value() == pytest.approx(1.5)
    assert dialog.azimuth_box.value() == pytest.approx(10.0)
    contour = lv.viewer_windows[-1]
    dialog.reject()
    QApplication.processEvents()
    assert not contour.isVisible()


def test_cut_k_conversion_picks_gamma_with_a_line_and_remembers(launcher):
    lv = launcher
    _select(lv, "c")
    lv.ui.ToolButtons["kcut"].click()
    window = lv.viewer_windows[-1]
    dialog = window._cut_k_dialog
    assert not window.functions_button.isVisibleTo(window)
    assert dialog.gamma_line.getViewBox() is window.view.view_box
    dialog.gamma_line.setValue(2.5)                  # dragged
    assert dialog.gamma_slit.value() == pytest.approx(2.5)
    geometry = dialog.geometry()
    assert geometry["gamma_deflector_deg"] == geometry["deflector_deg"]
    assert geometry["azimuth_deg"] == 0.0
    assert not dialog.advanced.isChecked()
    dialog.advanced.setChecked(True)
    dialog.gamma_defl.setValue(3.0)
    assert dialog.geometry()["gamma_deflector_deg"] == 3.0
    dialog.advanced.setChecked(False)
    dialog.accept()
    QApplication.processEvents()
    assert not window.isVisible()
    assert dialog.gamma_line.scene() is None             # taken off the cut

    _select(lv, "c")
    lv.ui.ToolButtons["kcut"].click()
    dialog = lv.viewer_windows[-1]._cut_k_dialog
    assert dialog.restored and dialog.gamma_slit.value() == pytest.approx(2.5)
    assert dialog.gamma_line.value() == pytest.approx(2.5)
    dialog.reject()


def test_a_k_cut_is_not_converted_again(launcher):
    lv = launcher
    k = np.linspace(-0.5, 0.5, 21)
    lv.add_memory_dataset(MemoryData("cut", (k, EKIN), np.ones((k.size, EKIN.size)),
                                     {"x": "k (Å⁻¹)", "y": "Kinetic energy (eV)"},
                                     source_label="kk"), "kk")
    _select(lv, "kk")
    assert not lv.ui.ToolButtons["kcut"].isEnabled()
    assert "already in momentum" in lv.ui.ToolButtons["kcut"].toolTip()


def test_kz_conversion_remembers(app, monkeypatch):
    from ui.panels import kzconv
    hv = np.linspace(30, 60, 7)
    rng = np.random.default_rng(2)
    data = MemoryData("kz_map", (hv, SLIT, EKIN),
                      rng.random((hv.size, SLIT.size, EKIN.size)) + 1,
                      {"x": "Photon energy (eV)", "k": "Angle, slit (deg)",
                       "z": "Kinetic energy (eV)"}, source_label="kz")
    contour = W.ContourWindow(data, "kz", "gray", False)
    monkeypatch.setattr(kzconv.nxs_jobs, "run_job", lambda *a, **k: True)
    dialog = kzconv.KzConversionDialog(contour)
    assert not dialog.restored
    dialog.inner_potential.setValue(15.5)
    dialog.lat_c.setValue(7.25)
    dialog.convert()
    dialog.close()
    again = kzconv.KzConversionDialog(contour)
    assert again.restored
    assert again.inner_potential.value() == pytest.approx(15.5)
    assert again.lat_c.value() == pytest.approx(7.25)
    again.close()
    contour.close()


def test_the_panel_opens_on_its_basic_tab(launcher):
    tabs = launcher.ui.ToolTabs
    assert [tabs.tabText(i) for i in range(tabs.count())] == ["Basic", "Advanced"]
    assert tabs.currentIndex() == 0
    basic = launcher.ui.ToolSections["Basic"]
    assert tabs.currentWidget().isAncestorOf(basic) or tabs.currentWidget().widget() is basic
    assert launcher.ui.ToolSections["Basic"].isAncestorOf(launcher.ui.ToolButtons["bz"])


def _k_map():
    kx = np.linspace(-1.2, 1.2, 25)
    ky = np.linspace(-1.0, 1.0, 21)
    data = MemoryData("k_map", (kx, ky, EKIN),
                      np.random.default_rng(4).random((kx.size, ky.size, EKIN.size)) + 1,
                      {"x": "kx (Å⁻¹)", "k": "ky (Å⁻¹)", "z": "Kinetic energy (eV)"},
                      source_label="km")
    # a hexagonal 2-D layer, so a plot needs no further typing
    data.scan.info.update({"bz.mode": "2d", "bz.a2": 3.0, "bz.b2": 3.0,
                           "bz.gamma2": 120.0, "bz.tile": 1})
    return data


def test_a_brillouin_zone_becomes_part_of_the_dataset(launcher, tmp_path):
    lv = launcher
    lv.add_memory_dataset(_k_map(), "km")
    _select(lv, "m")
    assert not lv.ui.ToolButtons["bz"].isEnabled()          # angle map
    _select(lv, "km")
    assert lv.ui.ToolButtons["bz"].isEnabled()
    lv.ui.ToolButtons["bz"].click()
    contour = lv.viewer_windows[-1]
    dialog = contour._bz_dialog
    assert not contour.functions_button.isVisibleTo(contour)
    assert dialog.plot() is not None
    key = lv.selected_key()
    assert lv.loaded_items[key]["properties"]["bz.overlay"]
    dialog.close()
    QApplication.processEvents()
    assert not contour.isVisible()                         # tied together

    # a viewer opened later offers it, off until ticked
    viewer = lv.open_selected()
    assert viewer.bz_check.isVisibleTo(viewer) and not viewer.bz_check.isChecked()
    assert viewer.view.bz_overlay_state() is None
    viewer.bz_check.setChecked(True)
    state = viewer.view.bz_overlay_state()
    assert state and state["layers"]

    # and its exports can draw it
    from ui.viewers.widgets import ExportDialog
    export = ExportDialog(viewer.view, ExportDialog.TAB_AXES)
    assert export.eps_zone.isVisibleTo(export)
    path = str(tmp_path / "axes.eps")
    export.save_axes(path)
    assert "setrgbcolor" in open(path, encoding="latin-1").read()
    export.eps_zone.setChecked(False)
    export.save_axes(path)
    assert "setrgbcolor" not in open(path, encoding="latin-1").read()
    export.with_axes.setChecked(False)
    with_zone = export.rendered()
    export.with_zone.setChecked(False)
    assert not np.array_equal(with_zone, export.rendered())
    export.close()
    viewer.close()

    # and it is written to the file the row is saved to
    item, problem = lv._dataset_for_saving(key)
    assert problem is None and item["info"]["bz.overlay"]


def test_a_tall_dialog_scrolls_and_keeps_its_buttons(app):
    from PyQt5.QtWidgets import QDialog, QVBoxLayout, QLabel, QScrollArea
    from ui.common.screen import fit_to_screen
    from ui.viewers.widgets import ButtonRow, add_button
    dialog = QDialog()
    layout = QVBoxLayout(dialog)
    for i in range(80):
        layout.addWidget(QLabel(f"row {i}"))
    buttons = ButtonRow()
    close = add_button(buttons, "Close", QDialogButtonBox.RejectRole)
    layout.addWidget(buttons)
    tall = dialog.minimumSizeHint().height()
    fit_to_screen(dialog)
    area = dialog.findChild(QScrollArea)
    assert area is not None and not area.widget().isAncestorOf(close)
    assert dialog.isAncestorOf(close)                     # kept below it
    assert dialog.minimumSizeHint().height() < tall / 3   # can be made short
    screen = QApplication.primaryScreen().availableGeometry()
    assert dialog.height() <= screen.height()


def test_the_tall_tool_windows_fit_the_screen(app):
    """The windows found taller than a laptop screen open no taller than
    the screen, and can be made shorter still."""
    from PyQt5.QtWidgets import QScrollArea
    import tools.spectra.fermi as F
    energy = np.linspace(55.0, 56.0, 120)
    cut = MemoryData("cut", (SLIT, energy),
                     100 * F.fermi_function(energy, 55.6, 20.0)[None, :]
                     * np.ones((SLIT.size, 1)) + 1,
                     {"x": "Angle (deg)", "y": "Kinetic energy (eV)"},
                     source_label="c")
    window = W.CutWindow(cut, "c", "gray", False)
    dialog = window.open_fermi_fit()
    screen = QApplication.primaryScreen().availableGeometry()
    assert dialog.findChild(QScrollArea) is not None
    assert dialog.height() <= screen.height()
    assert dialog.minimumSizeHint().height() < 300
    dialog.close()
    window.close()


def test_the_slice_series_figure_is_a_tool(launcher, monkeypatch):
    from PyQt5.QtWidgets import QDialog
    from ui.panels import figure
    lv = launcher
    _select(lv, "m")
    assert lv.ui.ToolButtons["slicefig"].text() == "Slice series figure..."
    assert lv.ui.ToolSections["Advanced"].isAncestorOf(lv.ui.ToolButtons["slicefig"])
    monkeypatch.setattr(figure.SliceSeriesDialog, "exec_", lambda self: QDialog.Accepted)
    window = lv.open_slice_series()
    assert window is not None and window.windowTitle().startswith("Slice series")
    assert len(window.figure.panels) >= 2
    assert not any(isinstance(w, W.ContourWindow) and w.isVisible()
                   for w in lv.viewer_windows)          # no contour left behind
    window.close()


def test_the_tool_panel_uses_the_main_panels_font(launcher):
    """No size of its own: the tool buttons and tabs read like the list's
    buttons beside them."""
    ui = launcher.ui
    main = ui.ListButtons["truncate"]
    main.ensurePolished()
    for button in ui.ToolButtons.values():
        button.ensurePolished()
        assert button.font().pointSizeF() == main.font().pointSizeF()
    bar = ui.ToolTabs.tabBar()
    bar.ensurePolished()
    assert bar.font().pointSizeF() == main.font().pointSizeF()


def test_tool_buttons_fit_and_do_not_overlap_in_fusion(app, monkeypatch):
    """macOS runs the launcher in Fusion (its native buttons overlapped in
    the tool panel); there every tool button keeps its full height, sits
    clear of its neighbours, and the panel is wide enough for its label."""
    from ui.common import platform_style
    from ui.launcher import main_window as main_window_ui
    from PyQt5.QtWidgets import QMainWindow
    old_style = app.style().objectName()
    monkeypatch.setattr(platform_style, "IS_MACOS", True)
    platform_style.apply_style(app)
    try:
        win = QMainWindow()
        ui = main_window_ui.Ui_MainWindow()
        ui.setupUi(win)
        win.show()
        for index in range(ui.ToolTabs.count()):
            ui.ToolTabs.setCurrentIndex(index)
            app.processEvents()
            page = ui.ToolTabs.widget(index).widget()
            buttons = [b for b in page.findChildren(QPushButton)]
            buttons.sort(key=lambda b: b.y())
            for above, below in zip(buttons, buttons[1:]):
                assert above.geometry().bottom() < below.geometry().top()
            for button in buttons:
                assert button.height() >= button.sizeHint().height()
                assert button.width() >= button.sizeHint().width()
        win.close()
    finally:
        app.setStyle(old_style)


def test_modifier_key_is_named_for_the_platform(monkeypatch):
    import importlib
    import sys
    from ui.common import platform_style
    monkeypatch.setattr(sys, "platform", "darwin")
    assert importlib.reload(platform_style).MOD_KEY == "⌘"
    monkeypatch.setattr(sys, "platform", "win32")
    assert importlib.reload(platform_style).MOD_KEY == "Ctrl"


def test_keep_recording_logs_each_file_added_to_the_list(launcher, tmp_path,
                                                         monkeypatch):
    import csv
    from loader.nxs_file import save_dataset
    from tools.output import explog as EL
    lv = launcher
    monkeypatch.setattr("ui.panels.explog.excel_conventions",
                        lambda: (";", ","))
    monkeypatch.setattr(lv, "exp_log_keep_recording", False)
    _select(lv, "c")
    _set_up_log(lv, tmp_path, "auto.csv",
                EL.Template([EL.Column("File name", "@file"),
                             EL.Column("hv", "info:photon_energy_eV")]),
                keep=True)
    assert lv.exp_log_keep_recording
    path = str(tmp_path / "new_cut.nxs")
    save_dataset(path, [{"name": "new_cut", "kind": "cut",
                         "axes": (SLIT, EKIN),
                         "labels": {"x": "Angle (deg)", "y": "Energy (eV)"},
                         "value": np.ones((SLIT.size, EKIN.size)),
                         "info": {"photon_energy_eV": 95.5}, "motors": {}}])
    lv.add_files([path], background=False)
    with open(tmp_path / "auto.csv", encoding="utf-8-sig") as handle:
        rows = list(csv.reader(handle, delimiter=";"))
    assert rows == [["File name", "hv"], ["new_cut", "95,5"]]
