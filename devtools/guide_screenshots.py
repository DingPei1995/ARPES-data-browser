"""Regenerate the screenshots in ``docs/images/`` for the user guides.

Runs the real program offscreen on the synthetic data of
``devtools/demo_data.py``, opens each window the guides describe, and saves
what it shows. Nothing here is part of the program; run it by hand after a
change to the windows it pictures.

Usage:  python devtools/guide_screenshots.py [output folder] [shot,shot,...]
        (default: docs/images/ and every shot; the demo data is written to
        a temporary folder, as are the program's config and session files)
"""
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
# Keep the program's config and session folders out of the real home.
_HOME = os.path.join(tempfile.gettempdir(), "arpes_guide")
os.makedirs(_HOME, exist_ok=True)
os.environ["HOME"] = _HOME
os.environ["USERPROFILE"] = _HOME

import numpy as np                                        # noqa: E402
from PyQt5 import QtCore                                  # noqa: E402
from PyQt5.QtGui import QFont                             # noqa: E402
from PyQt5.QtWidgets import QApplication, QPushButton, QWidget  # noqa: E402

OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "docs", "images")
ONLY = set(sys.argv[2].split(",")) if len(sys.argv) > 2 else None


def pump(n=6):
    for _ in range(n):
        QApplication.processEvents()
        QtCore.QThread.msleep(15)


def shot(widget: QWidget, name: str, size=None):
    """Save what ``widget`` shows as ``<OUT>/<name>.png``."""
    if size is not None:
        widget.resize(*size)
    widget.show()
    pump()
    widget.grab().save(os.path.join(OUT, f"{name}.png"))
    print("saved", name)


def grab_modal(name, then="reject", setup=None, size=None, tries=40):
    """Capture the modal dialog the next call is about to ``exec_()``, then
    close it with ``then`` ("accept" or "reject")."""
    def attempt(left=tries):
        dialog = QApplication.activeModalWidget()
        if dialog is None:
            if left:
                QtCore.QTimer.singleShot(50, lambda: attempt(left - 1))
            return
        try:
            if setup is not None:
                setup(dialog)
            if size is not None:
                dialog.resize(*size)
            pump()
            dialog.grab().save(os.path.join(OUT, f"{name}.png"))
            print("saved", name)
        finally:
            getattr(dialog, then)()
    QtCore.QTimer.singleShot(100, attempt)


def wait_for(condition, seconds=60.0):
    """Pump events until ``condition()`` is true (a worker-thread job has
    finished, say)."""
    import time
    end = time.time() + seconds
    while time.time() < end:
        QApplication.processEvents()
        if condition():
            return True
        QtCore.QThread.msleep(20)
    raise TimeoutError("gave up waiting")


def _quiet_message_boxes():
    """No one is there to click a message box: print it and answer Yes."""
    from PyQt5.QtWidgets import QMessageBox

    def say(kind):
        def box(parent, title, text, *args, **kwargs):
            print(f"[{kind}] {title}: {text}")
            return QMessageBox.Yes
        return staticmethod(box)
    for kind in ("information", "warning", "question", "critical"):
        setattr(QMessageBox, kind, say(kind))


def main():
    os.makedirs(OUT, exist_ok=True)
    app = QApplication(sys.argv[:1])
    # The virtual screen is 800 x 600; the windows are shown as they open
    # on a real one.
    from ui.common import screen
    screen.LIMIT_TO_SCREEN = False
    _quiet_message_boxes()
    app.setFont(QFont("DejaVu Sans", 9))
    from ui.common import gcguard
    gcguard.install(app)

    import ARPES_viewer as av
    from ui.launcher import main_window as main_window_ui
    from devtools import demo_data

    data_dir = os.path.join(_HOME, "demo_data")
    paths = demo_data.write_all(data_dir)

    av.win = av._MainWindow()
    av.ui = main_window_ui.Ui_MainWindow()
    av.ui.setupUi(av.win)
    av.InitializeUIWidgets()
    av.SetConnect()
    av.ui.FilePathLE.setText(data_dir)
    av.update_memory_status()
    av.win.show()
    pump()

    shots = Shots(av, app, paths, data_dir)
    import traceback
    for name in shots.ORDER:
        if ONLY and name not in ONLY and name != "main_window":
            continue
        try:
            getattr(shots, name)()
        except Exception:                                   # noqa: BLE001
            print("FAILED", name)
            traceback.print_exc()
        sys.stdout.flush()
    os._exit(0)


def demo_EF():
    from devtools import demo_data
    return demo_data.EF_KIN


class Shots:
    ORDER = ["main_window", "loader", "information", "exp_log", "map_viewer",
             "cut_viewer",
             "fermi", "fermi_map", "cut_k", "mdc_fit", "map_k", "brillouin",
             "arbitrary_cut",
             "fs_correction", "process", "stack", "volume", "cut_arithmetic",
             "kz_processing", "kz_conversion", "degrid", "curves", "spin",
             "figure", "export", "data_operations", "spatial", "slice_figure",
             "plot_tools", "combine"]

    def __init__(self, av, app, paths, data_dir):
        self.av, self.app, self.paths, self.data_dir = av, app, paths, data_dir

    def _select(self, text):
        lw = self.av.ui.FilePathListWidget
        lw.clearSelection()
        for i in range(lw.count()):
            if text in lw.item(i).text():
                lw.setCurrentItem(lw.item(i))
                lw.item(i).setSelected(True)
                return lw.item(i)
        raise KeyError(text)

    def main_window(self):
        self.av.add_files(self.paths, background=False)
        self._select("FS_map_hv60")
        self.av.update_memory_status()
        shot(self.av.win, "main_window", (860, 900))
        self.av.ui.ToolTabs.setCurrentIndex(1)
        pump()
        shot(self.av.win, "main_window_advanced", (860, 900))
        self.av.ui.ToolTabs.setCurrentIndex(0)

    def _open(self, text):
        self._select(text)
        window = self.av.open_selected()
        pump()
        return window

    def loader(self):
        from ui.launcher.loader_dialog import LoaderDialog
        dialog = LoaderDialog(self.av.win, self.data_dir)
        dialog.set_paths([p for p in self.paths if "map" in p or "cuts" in p])
        shot(dialog, "loader_dialog", (640, 720))
        dialog.reject()

    def information(self):
        self._select("Gamma_cut_LH")
        window = self.av.show_information()
        shot(self.av.win, "main_window_info", (860, 900))
        shot(window, "information_window", (520, 520))
        window.close()

    def exp_log(self):
        from tools.output import explog as EL
        from ui.panels.explog import ExpLogFormatDialog
        self._select("FS_map_hv60")
        dialog = ExpLogFormatDialog(
            self.av.win, entries=self.av._log_entries, load=self.av._log_load,
            selected=self.av.selected_keys, folder=self.data_dir,
            system=(";", ","))
        dialog.set_template(EL.Template([
            EL.Column("File name", "@file"), EL.Column("Data Type", "@kind"),
            EL.Column("Scan", "@scan1"), EL.Column("Energy", "@energy"),
            EL.Column("Comments", "@blank")]))
        dialog.folder_edit.setText(self.data_dir)
        dialog.name_edit.setText("exp_log.csv")
        dialog._sync_to_file()
        shot(dialog, "exp_log_dialog", (1000, 760))
        dialog.reject()

    def map_viewer(self):
        contour = self._open("FS_map_hv60")
        contour.e_control.set_index(int(np.argmin(abs(contour.E - (demo_EF() - 0.05)))))
        contour.view.set_readout_cursor_visible(True)
        contour.view.move_readout_to(2.0, 4.0)
        pump()
        shot(contour, "map_contour", (760, 820))
        slit = contour.open_cut("slit")
        shot(slit, "map_slit_cut", (820, 760))
        defl = contour.open_cut("deflector")
        shot(defl, "map_scan_axis_cut", (820, 760))
        self.contour = contour

    def _newest(self, kind=None):
        """The viewer opened most recently (of a class named ``kind``)."""
        windows = [w for w in self.av.viewer_windows
                   if kind is None or type(w).__name__ == kind]
        return windows[-1]

    def _open_listed(self, text):
        """Open the newest list row whose label contains ``text``."""
        lw = self.av.ui.FilePathListWidget
        for i in reversed(range(lw.count())):
            if text in lw.item(i).text():
                lw.clearSelection()
                lw.setCurrentItem(lw.item(i))
                lw.item(i).setSelected(True)
                window = self.av.open_selected()
                pump()
                return window
        raise KeyError(text)

    def cut_viewer(self):
        cut = self._open("Gamma_cut_LH")
        cut.view.set_readout_cursor_visible(True)
        cut.view.move_readout_to(-8.5, demo_EF() - 0.08)
        cut.panel.edc_width.setValue(0.3) if hasattr(cut.panel, "edc_width") else None
        pump()
        shot(cut, "cut_viewer", (1000, 760))
        self.cut = cut

    def fermi(self):
        self._select("Gamma_cut_LH")
        dialog = self.av.open_fermi_fit()          # the tool panel's route
        cut = dialog.window_
        pump()
        dialog.box_button.setChecked(True)
        dialog.energy_from.setValue(demo_EF() - 0.30)
        dialog.energy_to.setValue(demo_EF() + 0.10)
        dialog.angle_from.setValue(-6.0)
        dialog.angle_to.setValue(6.0)
        pump()
        dialog.estimate()
        dialog.run_fit()
        pump()
        shot(cut, "fermi_cut_box", (1000, 760))
        shot(dialog, "fermi_dialog")
        dialog.offset_energy_axis()
        pump()
        dialog.close()
        self.ef_cut = self._open_listed("Gamma_cut_LH")

    def fermi_map(self):
        """A map's Fermi edge fitted cut by cut along the deflector."""
        self._select("FS_map_hv60")
        dialog = self.av.open_fermi_fit()
        cut = dialog.window_
        pump()
        dialog.energy_from.setValue(demo_EF() - 0.25)
        dialog.energy_to.setValue(demo_EF() + 0.12)
        dialog.estimate()
        dialog.run_fit()
        dialog.offset_energy_axis(whole_map=False)
        size = cut.control.slider.maximum()
        for index in (size // 6, size - size // 6):
            cut.control.slider.setValue(index)
            pump()
            dialog.estimate()
            dialog.run_fit()
        pump(5)
        shot(dialog, "fermi_map_cuts", (1050, 900))
        shot(cut, "fermi_map_slit_cut", (820, 700))
        dialog.close()

    def cut_k(self):
        self._select("Gamma_cut_LH [E-Ef]")
        dialog = self.av.open_k_conversion("cut")   # the tool panel's route
        pump()
        dialog.gamma_line.setValue(0.0)
        pump()
        shot(dialog.window_, "cut_k_line", (1000, 700))
        shot(dialog, "cut_k_dialog", (520, 760))
        dialog.accept()
        pump()
        self.k_cut = self._open_listed("[k]") if any(
            "[k]" in self.av.ui.FilePathListWidget.item(i).text()
            for i in range(self.av.ui.FilePathListWidget.count())) else self._newest()
        shot(self.k_cut, "cut_in_k", (1000, 760))

    def mdc_fit(self):
        from tools.spectra import peaks as PK
        self._select("Gamma_cut_LH [E-Ef] [k]")
        panel = self.av.open_curve_fit("mdc")       # the tool panel's route
        pump()
        right = panel.bands[0] if panel.bands else panel.add_band()
        left = panel.add_band()
        from devtools import demo_data
        for e in (-0.02, -0.1, -0.2, -0.3):
            k = float(demo_data.pocket_k(e))
            right.seeds.append(PK.Seed(e, k, 0.04, 400.0))
            left.seeds.append(PK.Seed(e, -k, 0.04, 400.0))
        panel.sync_bands()
        panel.refresh_line()
        pump()
        panel.run_fit()
        panel.position_box.setValue(-0.1)
        pump()
        panel.try_line()
        pump()
        guide = panel.start_guide()             # the guided fit, at its seeds
        guide.go(4)
        pump()
        shot(panel, "mdc_fit_panel", (1250, 900))
        guide.stop()
        self.fit_panel = panel
        disp = panel.open_dispersion()
        pump()
        if disp is not None:
            shot(disp, "dispersion_window")


    def _rows(self):
        lw = self.av.ui.FilePathListWidget
        return [lw.item(i).text() for i in range(lw.count())]

    def map_k(self):
        self._select("FS_map_hv60")
        before = len(self._rows())
        dialog = self.av.open_k_conversion("map")   # a contour of its own
        contour = dialog.contour
        contour.e_control.set_index(int(np.argmin(abs(contour.E - (demo_EF() - 0.05)))))
        contour.view.move_readout_to(0.0, 0.0)
        pump()
        shot(contour, "map_k_contour", (700, 780))
        shot(dialog, "map_k_dialog", (480, 820))
        dialog.accept()
        wait_for(lambda: len(self._rows()) > before)
        self.k_map = self._open_listed("[k-map]")
        self.k_map.e_control.set_index(len(self.k_map.E) - 22)
        pump()
        shot(self.k_map, "k_map_contour", (760, 820))

    def brillouin(self):
        self._select("[k-map]")
        dialog = self.av.open_brillouin_zone()      # a contour of its own
        window = dialog.contour
        window.e_control.set_index(len(window.E) - 22)
        dialog.mode_3d.setChecked(True)
        dialog.spacegroup_box.setValue(191)          # P6/mmm, hexagonal
        dialog.a_box.setValue(3.3)
        dialog.c_box.setValue(6.0)
        dialog.tile_cb.setChecked(True)
        pump()
        dialog.plot()
        pump()
        shot(dialog, "brillouin_dialog", (430, 760))
        dialog.close()
        pump()
        # Kept with the dataset: its viewer shows it with the tick box.
        self.k_map.bz_check.setVisible(True)
        self.k_map.close()
        self.k_map = self._open_listed("[k-map]")
        self.k_map.e_control.set_index(len(self.k_map.E) - 22)
        self.k_map.bz_check.setChecked(True)
        pump()
        shot(self.k_map, "k_map_with_bz", (760, 820))
        dialog = self.k_map.view._open_export_dialog(2)
        pump()
        shot(dialog, "export_dialog_bz")
        dialog.close()
        self.k_map.bz_check.setChecked(False)

    def arbitrary_cut(self):
        contour = self.k_map
        contour.open_arbitrary_cut()
        pump()
        dialog = contour._arbcut_dialog
        k, m = 4 * np.pi / (3 * 3.3), 2 * np.pi / (np.sqrt(3) * 3.3)
        # Gamma -> M -> K -> Gamma
        for x, y in ((0.0, 0.0), (0.0, m), (k * np.sqrt(3) / 2, k / 2), (0.0, 0.0)):
            dialog.add_point(x, y)
        pump()
        shot(contour, "arbitrary_cut_path", (760, 820))
        shot(dialog, "arbitrary_cut_dialog")
        before = len(self._rows())
        dialog.accept()
        wait_for(lambda: len(self._rows()) > before)
        window = self._open_listed("[Cut]")
        shot(window, "arbitrary_cut_result", (1000, 700))

    def fs_correction(self):
        cut = self.cut
        dialog = cut.open_fs_correction()
        pump()
        if dialog is None:
            dialog = cut._fs_dialog
        shot(cut, "fs_correction_cut", (1000, 760))
        shot(dialog, "fs_correction_dialog")
        dialog.reject()

    def process(self):
        self._select("Gamma_cut_LH")
        dialog = self.av.open_processing("curvature")
        pump()
        shot(dialog, "process_dialog", (1100, 900))
        self.process_dialog = dialog

    def stack(self):
        self._select("Gamma_cut_LH")
        self.av.open_stack_plot()
        pump()
        window = [w for w in QApplication.topLevelWidgets()
                  if type(w).__name__ == "StackWindow" and w.isVisible()][-1]
        shot(window, "stack_plot", (900, 760))

    def volume(self):
        self._select("FS_map_hv60")
        window = self.av.open_volume_view("notch")
        pump(20)
        shot(window, "volume_view", (1150, 860))
        grab_modal("volume_range_dialog", then="reject")
        window.change_sampling()
        window.close()
        window = self.av.open_volume_view("slices")
        window.style_target.setCurrentIndex(0)
        window.opacity.setValue(75.0)
        window.fade.setValue(25.0)
        pump(20)
        shot(window, "volume_slices", (1150, 860))
        window.close()
        self._select("FS_map_hv60")
        dialog = self.av.open_processing("smooth")
        pump()
        shot(dialog, "volume_process_dialog", (1100, 880))

    def cut_arithmetic(self):
        lw = self.av.ui.FilePathListWidget
        lw.clearSelection()
        for i in range(lw.count()):
            if any(n in lw.item(i).text() for n in ("Gamma_cut_LH   ", "Gamma_cut_LV   ")):
                lw.item(i).setSelected(True)
        dialog = self.av.open_cut_arithmetic()
        pump(10)
        if dialog is None:
            dialog = [w for w in QApplication.topLevelWidgets()
                      if type(w).__name__ == "CutArithmeticDialog" and w.isVisible()][-1]
        shot(dialog, "cut_arithmetic", (1250, 900))
        before = len(self._rows())
        dialog.to_list.click()
        wait_for(lambda: len(self._rows()) > before)


    def kz_processing(self):
        contour = self._open("hv_scan_30-90eV")
        contour.e_control.set_index(len(contour.E) - 20)
        pump()
        shot(contour, "kz_map_contour", (760, 820))
        self._select("hv_scan_30-90eV")
        dialog = self.av.open_kz_alignment()        # the tool panel's route
        contour = dialog.contour
        pump()
        before = len(self._rows())
        dialog.run()
        wait_for(lambda: len(self._rows()) > before)
        pump()
        shot(dialog, "kz_processing_dialog", (1000, 820))
        cut = contour.cut_windows.get("slit")
        if cut is not None:
            shot(cut, "kz_processing_cut", (820, 700))
        dialog.close()
        self.kz_aligned = self._open_listed("[kz map]")

    def kz_conversion(self):
        self.kz_aligned.close()
        self._select("[kz map]")
        lw = self.av.ui.FilePathListWidget
        for i in reversed(range(lw.count())):        # the aligned one
            if "[kz map]" in lw.item(i).text():
                lw.clearSelection()
                lw.item(i).setSelected(True)
                break
        dialog = self.av.open_kz_conversion()       # the tool panel's route
        pump(10)
        shot(dialog, "kz_conversion_dialog", (1250, 860))
        before = len(self._rows())
        dialog.convert()
        wait_for(lambda: len(self._rows()) > before)
        kz_k = self._open_listed("kz map (k)")
        kz_k.e_control.set_index(len(kz_k.E) - 20)
        pump()
        shot(kz_k, "kz_map_in_k", (760, 820))

    def degrid(self):
        self._select("FS_map_mesh")
        dialog = self.av.open_degrid()              # the tool panel's route
        pump()
        dialog.run(background=False)
        dialog.slice_box.setValue(48)
        pump(10)
        shot(dialog, "degrid_dialog", (1250, 900))
        dialog.close()

    def curves(self):
        cut = self.ef_cut
        cut.view.set_readout_cursor_visible(True)
        cut.view.move_readout_to(0.0, -0.25)
        cut.panel.curves.edc_width.setValue(1.0)
        pump()
        before = len(self._rows())
        buttons = [b for b in cut.findChildren(QPushButton) if b.text() == "EDC \u2192 list"]
        buttons[0].click()
        wait_for(lambda: len(self._rows()) > before)
        window = self._open_listed("[EDC]")
        shot(window, "curve_viewer", (1000, 700))
        fit = window.open_fit()
        pump()
        fit.lo.setValue(-0.95)
        fit.hi.setValue(0.10)
        fit.fermi_cut.setChecked(True)
        fit.resolution.setValue(0.010)
        fit.find_peaks(max_peaks=3)
        pump()
        fit.run_fit()
        pump()
        shot(fit, "curve_fit_dialog", (1100, 860))

    def spin(self):
        window = self._open("spin_EDC")
        shot(window, "spin_edc_viewer", (1000, 700))
        dialog = window.open_spin_analysis()
        pump()
        dialog.refresh()
        pump()
        shot(dialog, "spin_analysis", (1150, 900))

    def figure(self):
        lw = self.av.ui.FilePathListWidget
        lw.clearSelection()
        for i in range(lw.count()):
            text = lw.item(i).text()
            if text.startswith(("Gamma_cut_LH   ", "Gamma_cut_LV   ")) or "LD" in text:
                lw.item(i).setSelected(True)
        window = self.av.plot_selected_as_figure()
        pump(10)
        if window is None:
            window = [w for w in QApplication.topLevelWidgets()
                      if type(w).__name__ == "FigureWindow" and w.isVisible()][-1]
        # A page laid out by hand: a zoom-in inset on (a), an arrow and a note.
        import tools.output.figure as nf
        window.figure.freeze()
        x, y, w, h = window.figure.panels[0].frame
        window.canvas._make_zoom(x + 0.35 * w, y + 0.1 * h, 0.3 * w, 0.3 * h)
        window.figure.shapes.append(nf.Shape(kind="arrow", x0=x + 0.1 * w,
                                             y0=y + 0.8 * h, x1=x + 0.3 * w,
                                             y1=y + 0.55 * h, color="#ffd400"))
        window.figure.shapes.append(nf.Shape(kind="text", x0=x + 0.05 * w,
                                             y0=y + 0.92 * h, text="band",
                                             color="#ffd400"))
        window.canvas.select("shape", 0)
        window.refresh()
        shot(window, "figure_composer", (1320, 860))

    def export(self):
        view = self.cut.view
        dialog = view._open_export_dialog(0)
        pump()
        shot(dialog, "export_dialog")
        dialog.close()

    def data_operations(self):
        self._select("Gamma_cut_LH")
        dialog = self.av.open_data_operations("normalise")  # the one with choices
        pump()
        shot(dialog, "data_operations")
        dialog.close()

    def spatial(self):
        window = self._open("SPEM_flake")
        pump(10)
        shot(window, "spatial_scan", (1250, 760))

    def slice_figure(self):
        self._select("FS_map_hv60")
        grab_modal("slice_figure_dialog", then="accept")
        window = self.av.open_slice_series()       # the tool panel's route
        pump(10)
        shot(window, "slice_figure", (1250, 760))

    def plot_tools(self):
        from ui.viewers.plot_tools import open_plot_tools
        cut = self.cut
        cut.view.set_readout_cursor_visible(False)
        dialog = open_plot_tools(cut)
        kinds = dialog.kind_combo
        for kind, setup, colour in (
                ("horizontal", lambda: (dialog.y_box.setValue(demo_EF()),
                                        dialog.text_edit.setText("E_F")),
                 "#ff3b30"),
                ("vertical", lambda: (dialog.x_box.setValue(0.0),
                                      dialog.text_edit.setText("k = 0")),
                 "#ffd60a"),
                ("angled", lambda: (dialog.x_box.setValue(-4.0),
                                    dialog.y_box.setValue(demo_EF() - 0.3),
                                    dialog.angle_box.setValue(60),
                                    dialog.text_edit.setText("")),
                 "#34c759")):
            dialog.list.setCurrentRow(-1)
            kinds.setCurrentIndex(kinds.findData(kind))
            setup()
            dialog.set_colour(colour)
            dialog.add()
        dialog.list.setCurrentRow(0)
        pump(10)
        shot(cut, "plot_tools_viewer", (1000, 760))
        shot(dialog, "plot_tools", (520, 680))
        dialog.clear()
        dialog.close()

    def combine(self):
        from ui.panels.stackcuts import CombineCutsDialog
        contour = self.contour
        x, k, z, _cube = contour.data.angle_cube
        x = np.asarray(x, dtype=float)
        half = x.size // 2
        labels = [contour.data.scan.labels.get(s, s) for s in ("x", "k", "z")]
        pieces = [((x[:half + 3], k, z), labels), ((x[half:], k, z), labels)]
        dialog = CombineCutsDialog(["FS_map_part1", "FS_map_part2"],
                                   reorder=True, parent=self.av.win,
                                   pieces=pieces, can_stack=False)
        shot(dialog, "combine_dialog", (520, 560))
        dialog.reject()

        cut = self.cut.data
        cut_axes = (np.asarray(cut.scan.x, dtype=float),
                    np.asarray(cut.scan.y, dtype=float))
        cut_labels = [cut.scan.labels.get("x", "x"), cut.scan.labels.get("y", "y")]
        names = [f"Gamma_cut_T{t}K" for t in (20, 60, 100, 140)]
        dialog = CombineCutsDialog(names, reorder=True, parent=self.av.win,
                                   pieces=[(cut_axes, cut_labels)] * 4)
        dialog.role.setCurrentIndex(dialog.role.findData("temperature"))
        dialog.start.setValue(20)
        dialog.step.setValue(40)
        shot(dialog, "combine_cuts_dialog", (520, 600))
        dialog.reject()


if __name__ == "__main__":
    main()
