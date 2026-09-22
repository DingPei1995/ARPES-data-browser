"""
ui/windows.py
=====================
The per-file viewer windows. The main ``ARPES_viewer`` window is only a
launcher -- a file list plus metadata -- and double-clicking a file opens
one of these on top of it:

- :class:`SpatialScanWindow` (``spem_4d`` / ``spem_1d``): the spatial map
  and, beside it, the E-vs-k spectrum at the cursor.
- :class:`CutWindow` (``cut``): the single E-vs-k spectrum.
- :class:`ContourWindow` (``map``): the constant-energy contour, which is
  what a deflector map is normally opened on. Its two buttons open the two
  orthogonal cuts through the cube, each in its own :class:`MapCutWindow`.

Why separate windows rather than the stacked pages this replaced: one
window per thing you are looking at means several files (or several cuts
of one file) can be on screen at once and arranged however you like, which
a single fixed layout cannot do.

Exporting belongs to the panel rather than to the window -- a window can
hold two pictures -- so it lives in each panel's own right-click menu (see
``ExportDialog`` in ``ui.widgets``): the data as ``.nxs``, the
picture as an image at a chosen size, and the axes alone as an ``.eps`` to
edit in Illustrator. The launcher keeps the metadata table and its CSV
export, since metadata belongs to the file rather than to any one view.
"""
from __future__ import annotations

import os

import numpy as np
import pyqtgraph as pg
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QPainterPath
from PyQt5.QtWidgets import (QWidget, QMainWindow, QVBoxLayout, QHBoxLayout,
                              QLabel, QPushButton, QFileDialog, QMessageBox,
                              QComboBox, QCheckBox, QDialog, QFormLayout,
                              QDoubleSpinBox, QSpinBox, QDialogButtonBox,
                              QApplication, QLineEdit, QInputDialog, QMenu,
                              QTableWidget, QTableWidgetItem, QHeaderView,
                              QScrollArea, QFrame, QGraphicsPathItem,
                              QTabWidget, QRadioButton, QGroupBox,
                              QProgressDialog)

from tools.analysis import arbitrary_cut, fit_feature, fs_correction
from tools.fermi import (PARAMETERS, PARAMETER_LABELS, K_B, initial_guess,
                       fermi_edge_model, fit_fermi_edge, fit_channels,
                       divide_fermi)
from tools.dataops import (truncate, self_normalize, compress, ARRAY_AXES,
                        CONSTRUCTOR_AXES)
from loader.nxs_file import CUBE_KINDS, energy_slot

from ui.widgets import (SpatialImageView, FrameImageView, ImagePanel,
                               FrameWindow, KMapData, MemoryData, _SliceControl,
                               ViewOptionsBar, apply_colormap, separator,
                               section_label, scroll_strip, strip_stock_menu,
                               COLORMAP_NAMES)


class ViewerWindow(QMainWindow):
    """Shared scaffolding for every viewer: title, colormap propagation,
    export of whatever this window currently shows, and deregistration on
    close."""

    closed = pyqtSignal(object)
    #: a dataset computed in this window (a saved slice, an arbitrary cut, a
    #: corrected map), for the launcher to list
    datasetCreated = pyqtSignal(object)

    def __init__(self, data, filename: str, colormap: str, flip: bool):
        super().__init__()
        self.data = data
        self.popout_windows = []
        self.filename = filename
        self.colormap = colormap
        self.flip = flip
        self.setWindowTitle(f"{filename}  [{data.kind}]")
        self.resize(1150, 820)

        self._central = QWidget()
        self.setCentralWidget(self._central)
        self.root = QVBoxLayout(self._central)
        self.root.setContentsMargins(4, 4, 4, 4)

        self.statusBar()

    # -- subclasses implement these ------------------------------------
    def image_panels(self):
        """Every ImagePanel this window owns, for colormap and display
        options."""
        return ()

    def exportable_panels(self):
        """``(suffix, array2d, x_axis, y_axis, x_label, y_label)`` for each
        panel worth writing out; arrays in (x, y) order matching the axes."""
        return []

    # -- shared behaviour ----------------------------------------------
    def build_toolbar(self, extra_widgets=()):
        """The two control rows along the top of every viewer.

        Row 1, left to right: any window-specific buttons (the cut and
        conversion buttons on a map), this window's colormap, and what to do
        with the slice on screen. Vertical rules separate the blocks so the
        row reads as groups rather than one strip of buttons.

        Exporting is deliberately *not* here: it belongs to the panel being
        looked at, not to the window, and a window can hold two of them. It
        lives in each panel's own right-click menu instead.

        Row 2 is the :class:`ViewOptionsBar`: axis ranges, invert, grid and
        alpha, which used to be buried in pyqtgraph's right-click menu.

        The colormap lives here rather than in the launcher so two windows can
        be tinted differently -- useful when comparing a faint map against a
        bright one, where one shared colormap suits neither.
        """
        bar_widget = QWidget()
        bar = QHBoxLayout(bar_widget)
        bar.setContentsMargins(0, 0, 0, 0)
        bar.setSpacing(6)
        for widget in extra_widgets:
            bar.addWidget(widget)
        if extra_widgets:
            bar.addWidget(separator())

        bar.addWidget(section_label("Colormap"))
        self.colormap_combo = QComboBox()
        self.colormap_combo.addItems(COLORMAP_NAMES)
        self.colormap_combo.setCurrentText(self.colormap)
        self.colormap_combo.setToolTip("Colormap for this window only.")
        self.colormap_combo.setMaximumWidth(132)
        bar.addWidget(self.colormap_combo)
        self.flip_cb = QCheckBox("Flip")
        self.flip_cb.setChecked(self.flip)
        self.flip_cb.setToolTip("Reverse the colormap (invert the intensity scale).")
        bar.addWidget(self.flip_cb)
        self.colormap_combo.currentTextChanged.connect(self._on_colormap_control_changed)
        self.flip_cb.toggled.connect(self._on_colormap_control_changed)

        bar.addWidget(separator())
        bar.addWidget(section_label("Slice"))
        self.save_slice_button = QPushButton("Save to the main list")
        self.save_slice_button.setToolTip(
            "Put the slice on screen in the launcher's list as a dataset of its "
            "own, to reopen, compare or save to a file later.")
        self.save_slice_button.clicked.connect(self.save_slice_to_list)
        self.popout_button = QPushButton("Open in a new panel")
        self.popout_button.setToolTip(
            "Open the slice on screen in its own frozen panel, for comparing "
            "two positions side by side.")
        self.popout_button.clicked.connect(self.pop_out_slice)
        bar.addWidget(self.save_slice_button)
        bar.addWidget(self.popout_button)

        bar.addStretch(1)

        # In a scroll strip like the row below it: a map window carries three
        # buttons plus the colormap and the slice actions, and without this
        # their combined width becomes the window's minimum -- which would
        # defeat opening a tall-narrow map at its own aspect ratio.
        self.root.addWidget(scroll_strip(bar_widget))

        # Directly under the colormap, as one visible row rather than three
        # levels of right-click menu.
        self.view_bar = ViewOptionsBar(self)
        self.root.addWidget(scroll_strip(self.view_bar))

    # -- the slice on screen, as a dataset or its own window --------------
    def slice_label(self) -> str:
        """How the slice currently on screen is identified -- the energy or
        angle it sits at. Windows that slice a cube override this; for a
        window showing one fixed image there is nothing to say."""
        return ""

    def _choose_panel(self, title: str):
        """Which panel to act on. One panel needs no question; two (a
        spatial scan) do, and the answer names them as the panels do."""
        panels = self.exportable_panels()
        if not panels:
            QMessageBox.warning(self, title, "Nothing displayed yet.")
            return None
        if len(panels) == 1:
            return panels[0]
        names = [entry[0] for entry in panels]
        choice, ok = QInputDialog.getItem(self, title, "Which panel?", names, 0, False)
        if not ok:
            return None
        return panels[names.index(choice)]

    def _slice_name(self, suffix: str) -> str:
        label = self.slice_label()
        stem = f"{self.filename} [{suffix}]"
        return f"{stem} {label}" if label else stem

    def save_slice_to_list(self):
        """Put the displayed slice in the launcher's list as a dataset of
        its own, the way the MATLAB tool's "Save Slice" puts it back in the
        workspace. From there it can be reopened, compared with others and
        saved to a file like any measurement."""
        entry = self._choose_panel("Save slice")
        if entry is None:
            return None
        suffix, array, xaxis, yaxis, xlabel, ylabel = entry
        default = self._slice_name(suffix)
        name, ok = QInputDialog.getText(self, "Save slice", "Name for the slice:",
                                        text=default)
        if not ok or not name.strip():
            return None
        data = MemoryData(
            "cut", (np.asarray(xaxis, dtype=float), np.asarray(yaxis, dtype=float)),
            np.array(array, dtype=float, copy=True),
            {"x": xlabel, "y": ylabel},
            source_label=name.strip(),
            parameters={"source": self.filename, "panel": suffix,
                        "position": self.slice_label() or "n/a"},
            prefix="slice",
            source_path=getattr(self.data, "path", ""),
            source_info=dict(getattr(self.data.scan, "info", {})))
        self.datasetCreated.emit(data)
        self.statusBar().showMessage(f"Saved slice as {name.strip()}")
        return data

    def pop_out_slice(self):
        """Open the displayed slice in a window of its own, frozen: the
        MATLAB tool's "Plot in new figure". Several can be open at once,
        which is the point -- two energies, or two cuts, side by side while
        this window moves on."""
        entry = self._choose_panel("Pop out")
        if entry is None:
            return None
        suffix, array, xaxis, yaxis, xlabel, ylabel = entry
        label = self.slice_label()
        title = f"{self.filename} · {suffix}" + (f" · {label}" if label else "")
        window = FrameWindow(title, np.asarray(array, dtype=float),
                             np.asarray(xaxis, dtype=float),
                             np.asarray(yaxis, dtype=float),
                             x_label=xlabel, y_label=ylabel,
                             colormap=self.colormap, flip=self.flip)
        self._arm_figure_source(window)
        window.closed.connect(lambda w: self.popout_windows.remove(w)
                              if w in self.popout_windows else None)
        self.popout_windows.append(window)
        window.show()
        self.statusBar().showMessage(f"Popped out {suffix} {label}".strip())
        return window

    # -- Fermi-surface correction ----------------------------------------
    def fs_angle_axis(self):
        """The axis the correction is fitted along -- the angle axis of the
        image on screen. Windows that offer the correction provide it."""
        return None

    def open_fs_correction(self):
        """Collect the points, fit, and straighten. Modeless, because the
        points are picked on the image behind the dialog."""
        existing = getattr(self, "_fs_dialog", None)
        if existing is not None and existing.isVisible():
            existing.raise_()
            existing.activateWindow()
            return None
        entry = self.exportable_panels()
        if not entry:
            QMessageBox.warning(self, "Fermi-surface correction",
                                "Nothing displayed yet.")
            return None
        _, _, _, _, x_label, y_label = entry[0]
        dialog = FSCorrectionDialog(self, self.image_panels()[0].view, x_label, y_label)
        self._fs_dialog = dialog
        self.statusBar().showMessage(
            "Click along the feature that should be flat, then press Correct.")

        def on_accepted():
            coeffs, order = dialog.coefficients(), dialog.order_box.value()
            points, name = dialog.points(), dialog.name_box.text().strip()
            self._fs_dialog = None
            self.run_fs_correction(coeffs, name=name, order=order, points=points)

        dialog.accepted.connect(on_accepted)
        dialog.rejected.connect(lambda: setattr(self, "_fs_dialog", None))
        dialog.show()
        dialog.raise_()
        return dialog

    def run_fs_correction(self, coeffs, name=None, order=2, points=(),
                          open_window: bool = True):
        """Apply the fit to this window's data and list the result.

        Subclasses say *what* is corrected (one cut, or the whole cube a cut
        was taken from) through :meth:`fs_correction_target`; the shifting
        itself is the same operation either way.
        """
        target = self.fs_correction_target()
        if target is None or coeffs is None:
            return None
        values, angle_axis, energy_axis, angle_dim, energy_dim, kind, axes_of = target
        try:
            corrected, new_energy = fs_correction(
                values, angle_axis, energy_axis, coeffs,
                angle_dim=angle_dim, energy_dim=energy_dim)
        except Exception as exc:
            QMessageBox.warning(self, "Fermi-surface correction",
                                f"The correction failed:\n{exc}")
            return None

        scan = self.data.scan
        source_info = {k: v for k, v in scan.info.items() if not k.startswith("fscorr.")}
        data = MemoryData(
            kind, axes_of(new_energy), corrected, dict(scan.labels),
            source_label=name or f"{self.filename} [FS corr]",
            parameters={"order": order, "coefficients": list(map(float, coeffs)),
                        "points": [list(map(float, p)) for p in points],
                        "angle_axis": self.fs_angle_label()},
            prefix="fscorr",
            source_path=getattr(self.data, "path", ""),
            source_info=source_info,
            source_motors=dict(scan.fourd_info))
        self.datasetCreated.emit(data)
        grew = new_energy.size - np.asarray(energy_axis).size
        self.statusBar().showMessage(
            f"Fermi-surface correction applied (order {order}); the energy axis "
            f"grew by {grew} bins -- added to the file list")
        if open_window:
            self.open_computed(data)
        return data

    def fs_correction_target(self):
        """What the correction acts on, or None when this window does not
        offer it."""
        return None

    def fs_angle_label(self) -> str:
        return ""

    # -- Fermi-edge fitting ------------------------------------------------
    def open_fermi_fit(self):
        """Fit the Fermi edge of what this window shows. Modeless, because
        the region it fits is the selection box on the image behind it."""
        existing = getattr(self, "_fermi_dialog", None)
        if existing is not None and existing.isVisible():
            existing.raise_()
            existing.activateWindow()
            return None
        try:
            dialog = FermiFitDialog(self)
        except ValueError as exc:
            QMessageBox.warning(self, "Fermi level fitting", str(exc))
            return None
        self._fermi_dialog = dialog
        dialog.finished.connect(lambda *_: setattr(self, "_fermi_dialog", None))
        dialog.show()
        dialog.raise_()
        return dialog

    def offset_energy_axis(self, ef: float, info: dict = None):
        """List a copy of this dataset with the energy axis shifted so the
        fitted E_F is zero. The axis is relabelled ``E - E_F`` -- after the
        shift the numbers are no longer the analyser's kinetic energies, and
        a label that still said "Energy" would be a lie."""
        data = self.data
        scan = data.scan
        axes_names = CONSTRUCTOR_AXES.get(data.kind)
        energy_name = energy_slot(data.kind)
        if axes_names is None or energy_name is None:
            QMessageBox.warning(self, "Offset energy axis",
                                 f"{data.kind} has no energy axis to shift.")
            return None
        values = scan.value4d if data.kind == "spem_4d" else scan.value
        if hasattr(values, "materialise"):
            values = values.materialise()
        axes = []
        for name in axes_names:
            axis = np.asarray(getattr(scan, name), dtype=float)
            axes.append(axis - float(ef) if name == energy_name else axis)
        labels = dict(scan.labels)
        unit = _unit_of(labels.get(energy_name, "")) or "eV"
        labels[energy_name] = f"E - E_F ({unit})"

        name = f"{self.filename} [E-Ef]"
        naming = getattr(self, "existing_names", None)
        taken = set(naming()) if callable(naming) else set()
        if name in taken:
            n = 2
            while f"{name} {n}" in taken:
                n += 1
            name = f"{name} {n}"
        parameters = dict(info or {})
        parameters["ef_subtracted"] = float(ef)
        result = MemoryData(data.kind, tuple(axes), np.asarray(values, dtype=float),
                            labels, source_label=name, parameters=parameters,
                            prefix="fitEF", source_path=getattr(data, "path", ""),
                            source_info=dict(scan.info),
                            source_motors=dict(scan.fourd_info))
        self.datasetCreated.emit(result)
        self.statusBar().showMessage(
            f"Energy axis shifted by {ef:.6g} eV -- added to the file list")
        return result

    def divide_fermi_cutoff(self, fit, cutoff_kt: float, info: dict = None):
        """Divide this dataset by the resolution-broadened Fermi function and
        list the result (the MATLAB tool's "Corr FS")."""
        data = self.data
        scan = data.scan
        axes_names = CONSTRUCTOR_AXES.get(data.kind)
        energy_name = energy_slot(data.kind)
        if axes_names is None or energy_name is None:
            return None
        values = scan.value4d if data.kind == "spem_4d" else scan.value
        if hasattr(values, "materialise"):
            values = values.materialise()
        values = np.asarray(values, dtype=float)
        array_names = ARRAY_AXES[data.kind]
        energy_dim = array_names.index(energy_name)
        energy = np.asarray(getattr(scan, energy_name), dtype=float)

        divided = divide_fermi(
            energy, values, ef=fit.values["ef"],
            temperature=fit.values["temperature"],
            resolution=fit.values["resolution"], energy_dim=energy_dim,
            cutoff_kt=float(cutoff_kt),
            background=(fit.values["bkg0"], fit.values["bkg1"]))

        name = f"{self.filename} [dFD]"
        parameters = dict(info or {})
        parameters["cutoff_kT"] = float(cutoff_kt)
        result = MemoryData(
            data.kind,
            tuple(np.asarray(getattr(scan, n), dtype=float) for n in axes_names),
            divided, dict(scan.labels), source_label=name,
            parameters=parameters, prefix="dFD",
            source_path=getattr(data, "path", ""),
            source_info=dict(scan.info), source_motors=dict(scan.fourd_info))
        self.datasetCreated.emit(result)
        kept = float(np.isfinite(divided).mean() * 100)
        self.statusBar().showMessage(
            f"Fermi cut-off divided out above {cutoff_kt:g} kT "
            f"({kept:.0f}% of the cube kept) -- added to the file list")
        return result

    def fit_channels_here(self, *, temperature, fixed, window, half_width,
                          step=1, start=None):
        """Fit the edge channel by channel on this window's own data and draw
        the result over the image."""
        try:
            angles, energy, frame, _, _ = reference_frame(self.data)
        except ValueError as exc:
            QMessageBox.warning(self, "Fermi level fitting", str(exc))
            return None
        return run_channel_fit(self, angles, energy, frame,
                               temperature=temperature, fixed=fixed,
                               window=window, half_width=half_width,
                               step=step, start=start, draw_on=self.view)

    def adopt_child(self, window):
        """Pass this window's launcher hooks to one it just opened.

        A cut window is opened by its contour, not by the launcher, so it
        would otherwise have no way to list a dataset it computes (the whole
        point of the Fermi correction living there). Its results are
        forwarded up through this window, which the launcher does listen to.
        """
        window._name_source = getattr(self, "_name_source", None)
        window._dataset_opener = getattr(self, "_dataset_opener", None)
        window.datasetCreated.connect(self.datasetCreated.emit)
        return window

    def open_computed(self, data):
        """Show a dataset this window has just computed.

        Goes through the launcher when there is one, so the new window is
        registered, named and coloured exactly like a double-clicked one;
        standing alone (in tests) it opens the window directly.
        """
        opener = getattr(self, "_dataset_opener", None)
        if callable(opener):
            return opener(data)
        window = open_viewer(data, getattr(data, "source_label", "computed"),
                             self.colormap, self.flip)
        if window is not None:
            self.popout_windows.append(window)
            window.closed.connect(lambda w: self.popout_windows.remove(w)
                                  if w in self.popout_windows else None)
            window.show()
        return window

    def _on_colormap_control_changed(self, *_):
        self.set_colormap(self.colormap_combo.currentText(), self.flip_cb.isChecked())

    def set_colormap(self, name: str, flip: bool):
        self.colormap, self.flip = name, flip
        combo = getattr(self, "colormap_combo", None)
        if combo is not None and combo.currentText() != name:
            combo.blockSignals(True)
            combo.setCurrentText(name)
            combo.blockSignals(False)
        flip_cb = getattr(self, "flip_cb", None)
        if flip_cb is not None and flip_cb.isChecked() != flip:
            flip_cb.blockSignals(True)
            flip_cb.setChecked(flip)
            flip_cb.blockSignals(False)
        for panel in self.image_panels():
            apply_colormap(panel.view, name, flip)

    def apply_display_options(self):
        for panel in self.image_panels():
            panel.apply_display_options()
        # The panels only exist after build_toolbar has run, so the view bar
        # is pointed at them here -- every window ends its __init__ with this
        # call, which is exactly the moment they are all in place.
        view_bar = getattr(self, "view_bar", None)
        if view_bar is not None and not view_bar.panels:
            view_bar.bind(self.image_panels())
        self.arm_panel_exports()

    def _arm_figure_source(self, window):
        """Let a popped-out panel's figure composer reach the launcher's
        list (for further panels) and this dataset's own metadata (for an
        E_F line taken from the fit)."""
        window.source_info = dict(getattr(self.data.scan, "info", {}) or {})
        entries = getattr(self, "_dataset_entries", None)
        loader = getattr(self, "_dataset_loader", None)
        if callable(entries) and callable(loader):
            window.figure_source = {"entries": entries, "loader": loader}

    def panel_export_names(self):
        """A short name for each :meth:`image_panels` entry, in that order --
        what an export of that panel is called.

        The names already exist as the suffixes of :meth:`exportable_panels`
        ("Cut", "SpatialMap", ...), which are far better filenames than the
        panel's on-screen title; this default reuses them whenever the two
        lists line up, and a window whose lists are in a different order
        (the spatial scan) says so itself.
        """
        panels = self.exportable_panels()
        if panels and len(panels) == len(self.image_panels()):
            return tuple(entry[0] for entry in panels)
        return ()

    def arm_panel_exports(self):
        """Tell each panel what an export of it should be called and what
        metadata goes with it.

        Only plain values are handed over -- no bound method of this window,
        which would make the panel keep the window alive. The slice position
        is read back the other way round, by the panel asking its window at
        export time.
        """
        info = dict(getattr(getattr(self.data, "scan", None), "info", {}) or {})
        names = self.panel_export_names()
        for index, panel in enumerate(self.image_panels()):
            view = getattr(panel, "view", None)
            if view is None:
                continue
            view.export_prefix = self.export_stem()
            view.export_info = info
            view.export_source_path = getattr(self.data, "path", "") or ""
            if index < len(names):
                view.export_title = names[index]

    def existing_names(self):
        """Names already in the launcher's list, so anything computed here
        can be given one that is free. Empty when no launcher is
        listening."""
        receiver = getattr(self, "_name_source", None)
        return list(receiver()) if callable(receiver) else []

    def export_stem(self) -> str:
        """A filesystem-safe basename for this window's exports.

        The window is titled after the dataset, which for a file holding
        several is "scan.nxs \u00b7 DeflX_0003" -- not something to feed a
        file dialog: the separator is not a legal Windows filename character,
        and splitext would cut at the wrong dot. So the entry, when there is
        one, is appended to the file's stem with an underscore.
        """
        name = self.filename
        entry = None
        if " \u00b7 " in name:
            name, entry = name.split(" \u00b7 ", 1)
        stem = os.path.splitext(name)[0]
        if entry:
            stem = f"{stem}_{entry}"
        return "".join(ch for ch in stem if ch.isalnum() or ch in "._- ") or "export"

    def size_to_data_aspect(self, view, extra_width=60, extra_height=260,
                            min_w=520, max_w=1500, min_h=520, max_h=1000):
        """Shape the window so an equal-aspect image fills it.

        With a 1:1 lock the view has to stretch whichever axis is short for
        the panel's shape, which is what leaves broad empty margins beside a
        tall-narrow contour in a wide window. Nothing inside the plot can fix
        that -- the cure is to open the window at the data's own aspect, so
        the lock has nothing to stretch. ``extra_height`` is everything above
        and below the image (two control rows, title, level bar, options,
        status bar), so it has to grow whenever a row is added. The user is free to resize
        afterwards (margins come back, as geometry demands) or to untick the
        ratio box for a tight fill on both axes.
        """
        extent = view.data_extent()
        if extent is None:
            return
        x0, x1, y0, y1 = extent
        dx, dy = abs(x1 - x0), abs(y1 - y0)
        if dx <= 0 or dy <= 0:
            return
        height = min(max(self.height(), min_h), max_h)
        plot_height = max(height - extra_height, 100)
        width = int(round(plot_height * dx / dy)) + extra_width
        self.resize(int(min(max(width, min_w), max_w)), int(height))

    def closeEvent(self, event):
        # Each window owns its own NxsData (and therefore its own open file
        # handle), so closing the window releases it. Snapshot windows hold
        # plain numpy copies and are unaffected.
        for window in list(getattr(self, "figure_windows", ())):
            window.close()
        try:
            self.data.close()
        except Exception:
            pass
        self.closed.emit(self)
        super().closeEvent(event)


# --------------------------------------------------------------------------
# Real-space scans
# --------------------------------------------------------------------------
class SpatialScanWindow(ViewerWindow):
    """Spatial map plus the E-vs-k spectrum at the cursor -- double-clicking
    a spatial scan opens both together, since neither is much use without
    the other."""

    def __init__(self, data, filename, colormap, flip):
        super().__init__(data, filename, colormap, flip)
        scan = data.scan
        labels = scan.labels

        hint = QLabel(
            "Drag the red target to move the readout pixel. Right-click a panel for the "
            "readout cursor (EDC/MDC), a selection box and its integration, or a new "
            "window at this position.")
        hint.setWordWrap(True)
        self.root.addWidget(hint)
        self.build_toolbar()

        panels = QHBoxLayout()
        self.spatial_view = SpatialImageView()
        self.spatial_panel = ImagePanel("Spatial overview", self.spatial_view,
                                        lock_ratio=True, parent=self)
        self.frame_view = FrameImageView()
        self.frame_panel = ImagePanel("E vs k at cursor / selection", self.frame_view,
                                      curves=True, parent=self)
        panels.addWidget(self.spatial_panel)
        panels.addWidget(self.frame_panel)
        self.root.addLayout(panels, stretch=1)

        self.frame_label = ""   # what the spectrum panel currently holds
        #: Whether the spatial panel shows the scan's Y axis horizontally.
        #: The file says which stage was scanned first, and that is what this
        #: program calls X -- but a lab whose own convention is the other way
        #: round would see its maps transposed, so the view can be flipped
        #: without touching the data. Only the *display* swaps: the cursor
        #: readout, every export and every saved slice keep the axes the
        #: measurement had.
        self.swap_xy = False
        self._spatial_map = None      # the current map in (y, x) data order

        if data.kind == "spem_4d":
            self.frame_view.set_axis_labels(labels.get("k"), labels.get("z"))
            self._set_spatial_map(data.spatial_overview)
            self.spatial_view.enable_swap_action(True)
            self.spatial_view.swapToggled.connect(self.set_swap_xy)
            self.spatial_view.enable_selection_action(
                True, "Integrate selection -> E vs k panel")
            self.frame_view.enable_selection_action(
                True, "Integrate selection -> spatial panel")
            self._show_frame_at(len(scan.y) // 2, len(scan.x) // 2)
        else:  # spem_1d
            self.spatial_view.set_axis_labels(labels.get("x"), labels.get("y"))
            self.frame_view.set_axis_labels(labels.get("y"), labels.get("z"))
            self.spatial_view.set_data(data.line_kmap.T, scan.x, scan.y)
            # The two-way region integration needs a 2D spatial map to write
            # back onto, which a line scan does not have.
            self.spatial_view.enable_selection_action(False)
            self.frame_view.enable_selection_action(False)
            self._show_frame_at_x(len(scan.x) // 2)

        self.spatial_view.enable_popout_action(True, "Open this position in a new window")
        self.spatial_view.pixelChanged.connect(self._on_cursor_moved)
        self.spatial_panel.selectionApplied.connect(self._on_spatial_selection)
        self.frame_panel.selectionApplied.connect(self._on_frame_selection)
        self.spatial_view.popoutRequested.connect(self.open_popout)

        self.set_colormap(colormap, flip)
        self.apply_display_options()

    def image_panels(self):
        return (self.spatial_panel, self.frame_panel)

    def panel_export_names(self):
        # exportable_panels() lists the spectrum first (it is the panel the
        # window exists for); image_panels() lists the map first, because
        # that is the one on the left. Spelt out rather than derived.
        return ("SpatialMap", "EvsK")

    def slice_label(self) -> str:
        return self.frame_label

    # -- which way round the spatial map is drawn -------------------------
    def _set_spatial_map(self, map_yx, reset_cursor: bool = True):
        """Draw a spatial map given in the data's own (y, x) order, the way
        round the swap setting asks for."""
        scan = self.data.scan
        labels = scan.labels
        self._spatial_map = np.asarray(map_yx)
        if self.swap_xy:
            self.spatial_view.set_axis_labels(labels.get("y"), labels.get("x"))
            self.spatial_view.set_data(self._spatial_map.T, scan.y, scan.x,
                                       reset_cursor=reset_cursor)
        else:
            self.spatial_view.set_axis_labels(labels.get("x"), labels.get("y"))
            self.spatial_view.set_data(self._spatial_map, scan.x, scan.y,
                                       reset_cursor=reset_cursor)

    def set_swap_xy(self, on: bool):
        """Exchange the spatial panel's two axes, keeping the cursor on the
        same measured pixel."""
        on = bool(on)
        if on == self.swap_xy or self._spatial_map is None:
            return
        cursor = self.spatial_view.cursor.pos()
        self.swap_xy = on
        self._set_spatial_map(self._spatial_map, reset_cursor=False)
        # The pixel the cursor was on has not moved; its coordinates have
        # swapped places, so the cursor follows them.
        self.spatial_view.cursor.setPos((float(cursor.y()), float(cursor.x())))
        self.spatial_panel.sync_level_range()
        self.set_colormap(self.colormap, self.flip)
        self.spatial_view.fit_to_data()
        action = self.spatial_view.action_swap
        if action.isChecked() != on:
            action.blockSignals(True)
            action.setChecked(on)
            action.blockSignals(False)
        self.statusBar().showMessage(
            "Spatial map drawn with Y horizontal" if on
            else "Spatial map drawn with X horizontal")

    def _data_pixel(self, row: int, col: int):
        """(y index, x index) for a pixel the *display* reported."""
        return (col, row) if self.swap_xy else (row, col)

    # -- cursor / frames -------------------------------------------------
    def _show_frame_at(self, row, col):
        scan = self.data.scan
        self.frame_view.set_frame(self.data.frame_at(row, col), scan.k, scan.z)
        self.frame_panel.sync_curve_source()
        self.frame_panel.sync_level_range()
        self.frame_label = f"x={scan.x[col]:.4g}, y={scan.y[row]:.4g}"
        self.statusBar().showMessage(f"pixel ({self.frame_label})")

    def _show_frame_at_x(self, xi):
        scan = self.data.scan
        self.frame_view.set_frame(self.data.frame_at_x(xi), scan.y, scan.z)
        self.frame_panel.sync_curve_source()
        self.frame_panel.sync_level_range()
        self.frame_label = f"x={scan.x[xi]:.4g}"
        self.statusBar().showMessage(self.frame_label)

    def _on_cursor_moved(self, row, col):
        if self.data.kind == "spem_4d":
            self._show_frame_at(*self._data_pixel(row, col))
        else:
            self._show_frame_at_x(col)

    # -- region integration ----------------------------------------------
    def _on_spatial_selection(self):
        if self.data.kind != "spem_4d":
            return
        ranges = self.spatial_view.selection_index_ranges()
        if ranges is None:
            self.statusBar().showMessage("Draw a selection box first "
                                          "(right-click -> Selection box)")
            return
        (ia0, ia1), (ib0, ib1) = ranges
        # The ranges come back in display order; swapped, the horizontal one
        # is the scan's Y.
        (ix0, ix1), (iy0, iy1) = (((ib0, ib1), (ia0, ia1)) if self.swap_xy
                                  else ((ia0, ia1), (ib0, ib1)))
        scan = self.data.scan
        frame = self.data.frame_over_region(slice(iy0, iy1 + 1), slice(ix0, ix1 + 1))
        self.frame_view.set_frame(frame, scan.k, scan.z)
        self.frame_panel.sync_curve_source()
        self.frame_panel.sync_level_range()
        self.frame_label = (f"x[{scan.x[ix0]:.4g}..{scan.x[ix1]:.4g}], "
                            f"y[{scan.y[iy0]:.4g}..{scan.y[iy1]:.4g}] integrated")
        self.set_colormap(self.colormap, self.flip)
        self.statusBar().showMessage(
            f"E vs k integrated over spatial {self.spatial_view.selection_summary()}")

    def _on_frame_selection(self):
        if self.data.kind != "spem_4d":
            return
        ranges = self.frame_view.selection_index_ranges()
        if ranges is None:
            self.statusBar().showMessage("Draw a selection box first "
                                          "(right-click -> Selection box)")
            return
        (ik0, ik1), (ie0, ie1) = ranges
        scan = self.data.scan
        smap = self.data.spatial_over_region(slice(ik0, ik1 + 1), slice(ie0, ie1 + 1))
        # The cursor stays put: this changes what the map shows, not which
        # pixel the user was looking at.
        self._set_spatial_map(smap, reset_cursor=False)
        self.spatial_panel.sync_level_range()
        self.set_colormap(self.colormap, self.flip)
        self.statusBar().showMessage(
            f"Spatial map integrated over {self.frame_view.selection_summary()}")

    # -- comparison snapshots --------------------------------------------
    def open_popout(self):
        frame = self.frame_view.last_frame
        if frame is None:
            return
        scan = self.data.scan
        x_axis = scan.k if self.data.kind == "spem_4d" else scan.y
        window = FrameWindow(f"{self.filename}  @ {self.frame_label or 'current view'}",
                             frame, x_axis, scan.z,
                             x_label=self.frame_view.x_label,
                             y_label=self.frame_view.y_label,
                             colormap=self.colormap, flip=self.flip)
        self._arm_figure_source(window)
        window.closed.connect(lambda w: self.popout_windows.remove(w)
                              if w in self.popout_windows else None)
        self.popout_windows.append(window)
        window.show()
        self.statusBar().showMessage(f"Opened snapshot ({len(self.popout_windows)} open)")

    def exportable_panels(self):
        scan = self.data.scan
        labels = scan.labels
        out = []
        if self.frame_view.last_frame is not None:
            x_axis = scan.k if self.data.kind == "spem_4d" else scan.y
            x_label = labels.get("k" if self.data.kind == "spem_4d" else "y", "k")
            out.append(("EvsK", self.frame_view.last_frame, x_axis, scan.z,
                        x_label, labels.get("z", "E")))
        if self._spatial_map is not None:
            # Always the measurement's own orientation, whichever way the
            # panel happens to be drawn: an exported map should not depend on
            # how it was being looked at.
            out.append(("SpatialMap", np.asarray(self._spatial_map).T, scan.x, scan.y,
                        labels.get("x", "x"), labels.get("y", "y")))
        elif self.spatial_view._image2d is not None:
            out.append(("SpatialMap", self.spatial_view._image2d.T, scan.x, scan.y,
                        labels.get("x", "x"), labels.get("y", "y")))
        return out


# --------------------------------------------------------------------------
# Single cut
# --------------------------------------------------------------------------
class CutWindow(ViewerWindow):
    """One E-vs-k spectrum (a scan taken at a single deflector position)."""

    def __init__(self, data, filename, colormap, flip):
        super().__init__(data, filename, colormap, flip)
        scan = data.scan
        hint = QLabel("Right-click for the readout cursor: it reports the coordinates and "
                      "value under it and draws the EDC and MDC through that point.")
        hint.setWordWrap(True)
        self.root.addWidget(hint)
        self.fermi_button = QPushButton("Fermi level")
        self.fermi_button.setToolTip(
            "Fit the Fermi edge of this cut: E_F, the resolution, and -- in its "
            "Advanced corner -- dividing the cut-off out and fitting E_F channel "
            "by channel.")
        self.fermi_button.clicked.connect(self.open_fermi_fit)
        self.fs_button = QPushButton("FS correction")
        self.fs_button.setToolTip(
            "Straighten a curved feature: click along it, fit a polynomial, and "
            "shift every angle column in energy so it comes out flat.")
        self.fs_button.clicked.connect(self.open_fs_correction)
        self.kconv_button = QPushButton("Cut k conversion")
        self.kconv_button.setToolTip(
            "Convert this cut from degrees to Å⁻¹. Needs to be "
            "told where Γ is, since a cut is a line that generally misses "
            "it -- inherit that from a converted map, or type the angles.")
        self.kconv_button.clicked.connect(self.open_cut_k_conversion)
        blocked = cut_k_conversion_blocked(data)
        if blocked:
            self.kconv_button.setEnabled(False)
            self.kconv_button.setToolTip(
                f"Not available here: {blocked}.")
        self.fit_button = QPushButton("MDC / EDC fit")
        self.fit_button.setToolTip(
            "Fit peaks line by line and read the band off the fitted centres: "
            "Fermi velocity, effective mass, self-energy.")
        self.fit_button.clicked.connect(self.open_mdc_edc_fit)
        # Deliberately only for converted cuts. A band fitted in degrees has
        # a slope in eV/degree, and every quantity built on it -- v_F, m*,
        # the self-energy -- would silently be in the wrong units. The button
        # says why rather than disappearing, so the route to it is obvious.
        from ui.fit import momentum_cut_reason
        reason = momentum_cut_reason(data)
        if reason:
            self.fit_button.setEnabled(False)
            # Only the first sentence: the rest of the explanation is the
            # message box's job, and a three-paragraph tooltip is unreadable.
            self.fit_button.setToolTip(
                f"Not available here: {reason.split(chr(10))[0]}")
        self.build_toolbar((self.fermi_button, self.fs_button, self.kconv_button,
                            self.fit_button))

        self.view = FrameImageView()
        # Titled from the axes themselves: a cut read from a file is an
        # E-vs-angle spectrum, but one computed along a path is energy
        # against distance, and the title should say so.
        title = (f"{scan.labels.get('y', 'E')} vs {scan.labels.get('x', 'k')}"
                 if scan.labels.get("x") else "E vs k (single cut)")
        self.panel = ImagePanel(title, self.view, curves=True, parent=self)
        self.root.addWidget(self.panel, stretch=1)

        self.view.set_axis_labels(scan.labels.get("x"), scan.labels.get("y"))
        self.view.set_frame(data.cut_frame, scan.x, scan.y)
        self.view.enable_selection_action(False)
        self._mark_path_corners()
        self.set_colormap(colormap, flip)
        self.apply_display_options()

    def image_panels(self):
        return (self.panel,)

    def exportable_panels(self):
        scan = self.data.scan
        if self.view.last_frame is None:
            return []
        return [("Cut", self.view.last_frame, scan.x, scan.y,
                 scan.labels.get("x", "k"), scan.labels.get("y", "E"))]

    # -- angle -> momentum -------------------------------------------------
    def open_cut_k_conversion(self):
        """The conversion dialog, modeless like the map's: the cut stays
        readable while the angles are being settled."""
        existing = getattr(self, "_cut_k_dialog", None)
        if existing is not None and existing.isVisible():
            existing.raise_()
            existing.activateWindow()
            return None
        dialog = CutKConversionDialog(self)
        self._cut_k_dialog = dialog

        def on_accepted():
            settings, name = dialog.settings(), dialog.output_name()
            self._cut_k_dialog = None
            self.run_cut_k_conversion(settings, name=name)

        dialog.accepted.connect(on_accepted)
        dialog.rejected.connect(lambda: setattr(self, "_cut_k_dialog", None))
        dialog.show()
        dialog.raise_()
        return dialog

    def open_mdc_edc_fit(self):
        """The MDC/EDC fit panel on this cut. Modeless and a window of its
        own, because a fit is worked at over several minutes and the cut
        behind it is what the seeds are read from."""
        from ui.fit import FitPanel, momentum_cut_reason

        reason = momentum_cut_reason(self.data)
        if reason:
            QMessageBox.information(self, "MDC / EDC fit", reason + ".")
            return None
        existing = getattr(self, "_fit_panel", None)
        if existing is not None and existing.isVisible():
            existing.raise_()
            existing.activateWindow()
            return existing
        panel = FitPanel(self.data, self.filename, self,
                         self.colormap, self.flip)
        self._fit_panel = panel
        panel.closed.connect(lambda *_: setattr(self, "_fit_panel", None))
        panel.datasetsCreated.connect(
            lambda made: [self.datasetCreated.emit(one) for one in made])
        panel.show()
        panel.raise_()
        return panel

    def run_cut_k_conversion(self, settings: dict, name: str = None):
        """Convert and list the result. Separate from the dialog so the
        conversion can be driven without one, which is how it is tested."""
        from tools.cutk import convert_cut, K_LABEL, K_RADIAL_LABEL

        scan = self.data.scan
        self.statusBar().showMessage("Converting the cut to k-space...")
        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            k_axis, energy, values, report = convert_cut(
                scan.x, scan.y, self.data.cut_frame, **settings)
        except Exception as exc:
            QMessageBox.warning(self, "Cut k conversion", f"Conversion failed:\n{exc}")
            self.statusBar().showMessage("Cut k conversion failed")
            return None
        finally:
            QApplication.restoreOverrideCursor()

        parameters = dict(settings)
        # The distance from Gamma is the one number the axes cannot carry,
        # and the whole reason this conversion needs a reference: it travels
        # with the result so a cut reopened next week still says where it ran.
        parameters.update({key: report[key] for key in
                           ("k_perp", "k_perp_min", "k_perp_max", "bow",
                            "energy_drift", "k_step", "measured_fraction")})
        source_info = {key: value for key, value in scan.info.items()
                       if not key.startswith("kcut.")}
        result = MemoryData(
            "cut", (np.asarray(k_axis, dtype=float), np.asarray(energy, dtype=float)),
            np.asarray(values, dtype=float),
            {"x": K_RADIAL_LABEL if settings.get("radial") else K_LABEL,
             "y": scan.labels.get("y", "Energy (eV)")},
            source_label=name or f"{self.filename} [k]",
            parameters=parameters, prefix="kcut",
            source_path=getattr(self.data, "path", ""),
            source_info=source_info,
            source_motors=dict(scan.fourd_info))
        self.datasetCreated.emit(result)
        self.statusBar().showMessage(
            f"Converted to k: {values.shape[0]}x{values.shape[1]}, "
            f"{report['measured_fraction'] * 100:.0f}% measured, "
            f"k⊥ = {abs(report['k_perp']):.4g} Å⁻¹ "
            f"-- added to the file list")
        return result

    def _mark_path_corners(self):
        """A cut taken along a path of several segments records where one
        segment ends and the next begins; mark them, as the MATLAB dialog's
        "boundary" option does, so a kink in a band is not mistaken for
        physics when it is only a corner in the path."""
        joints = self.data.scan.info.get("arbcut.joints")
        if joints is None:
            return
        for position in np.atleast_1d(joints):
            line = pg.InfiniteLine(pos=float(position), angle=90,
                                   pen=pg.mkPen("#2f7a3f", width=2, style=Qt.DashLine))
            line.setZValue(25)
            self.view.view_box.addItem(line, ignoreBounds=True)

    # -- Fermi-surface correction -----------------------------------------
    def fs_angle_axis(self):
        return self.data.scan.x

    def fs_angle_label(self) -> str:
        return self.data.scan.labels.get("x", "angle")

    def fs_correction_target(self):
        """The cut itself: (angle, energy), corrected into a new cut."""
        scan = self.data.scan
        return (np.asarray(self.data.cut_frame, dtype=float), scan.x, scan.y,
                0, 1, "cut", lambda new_energy: (np.asarray(scan.x, dtype=float),
                                                 new_energy))


# --------------------------------------------------------------------------
# Deflector map: contour first, cuts on demand
# --------------------------------------------------------------------------
class KConversionDialog(QDialog):
    """Settings for converting a Map from angle space to momentum space.

    Opens with the k-space origin already filled in from the readout cursor
    on the constant-energy contour -- picking the point that should become
    k = (0, 0) is the one thing the conversion genuinely needs from the user,
    and reading it off the contour is easier than typing two angles. The
    boxes stay editable, and keep following the cursor until the moment the
    user types into one of them (after that their value is theirs, and the
    "Read cursor" button is how they ask for the cursor's again).

    The remaining fields mirror the MATLAB dialog's, minus the parts that do
    not apply here: no lattice constant (the output is in A^-1, not pi/a),
    no k_z offset and no Ek (the 2-D-only box). A single Cut has its own
    dialog (:class:`CutKConversionDialog`), because it needs something a
    map does not: to be told where Gamma is.
    """

    def __init__(self, contour: "ContourWindow"):
        super().__init__(contour)
        self.contour = contour
        self.setWindowTitle("Map k conversion")
        self._user_edited = False

        layout = QVBoxLayout(self)
        intro = QLabel(
            "Pick the point on the constant-energy contour that should become "
            "the k-space origin (drag the readout cursor); its angles appear "
            "below and can be edited. Press OK to convert.")
        intro.setWordWrap(True)
        layout.addWidget(intro)

        form = QFormLayout()
        self.name_box = QLineEdit(self._default_name())
        self.name_box.setToolTip(
            "Name for the converted dataset in the file list. Converting the "
            "same map twice would otherwise give two identically-named rows, "
            "so this is pre-filled with a name that is not in use yet.")
        form.addRow("Name", self.name_box)

        self.theta_box = self._angle_box("Angle of the origin along the "
                                          "deflector axis (contour x).")
        self.phi_box = self._angle_box("Angle of the origin along the slit "
                                        "axis (contour y).")
        self.azimuth_box = self._angle_box(
            "Sample rotation about the surface normal. Leave at 0 unless the "
            "sample's axes are rotated against the analyser's.")
        form.addRow("Theta offset (deg)", self.theta_box)
        form.addRow("Phi offset (deg)", self.phi_box)
        form.addRow("Sample rotation (deg)", self.azimuth_box)

        self.rotation_button = QPushButton("Set rotation from contour...")
        self.rotation_button.setToolTip(
            "Pick a direction on the contour that should end up vertical: click "
            "two points along it (or one, and the line from the origin to it is "
            "used), and the rotation that stands it upright is filled in above.")
        self.rotation_button.setCheckable(True)
        self.rotation_button.toggled.connect(self._toggle_rotation_picking)
        form.addRow("", self.rotation_button)
        self.rotation_note = QLabel("")
        self.rotation_note.setWordWrap(True)
        self.rotation_note.setMinimumHeight(34)
        form.addRow("", self.rotation_note)
        self._rotation_points = []

        self.energy_offset_box = QDoubleSpinBox()
        self.energy_offset_box.setDecimals(4)
        self.energy_offset_box.setRange(-10000.0, 10000.0)
        self.energy_offset_box.setSingleStep(0.1)
        self.energy_offset_box.setToolTip(
            "Added to the kinetic energy before converting, for when the axis "
            "is not already a true kinetic energy.")
        form.addRow("Energy offset (eV)", self.energy_offset_box)

        # The output grid, asked for either way round: how many points, or
        # how far apart. Counts are the default because they say exactly how
        # big the result will be; a resolution says what it is worth
        # resolving, which is the question when comparing scans. One of the
        # two is always disabled, so there is never a pair of numbers
        # disagreeing about the same axis.
        n_e = len(contour.E)
        self.by_count = QRadioButton("by element number")
        self.by_resolution = QRadioButton("by resolution")
        self.by_count.setChecked(True)
        mode_row = QWidget()
        mode_box = QHBoxLayout(mode_row)
        mode_box.setContentsMargins(0, 0, 0, 0)
        mode_box.addWidget(self.by_count)
        mode_box.addWidget(self.by_resolution)
        mode_box.addStretch(1)
        form.addRow("Output grid", mode_row)

        self.nkx_box = self._count_box(100, "Number of kx points in the output.")
        self.nky_box = self._count_box(100, "Number of ky points in the output.")
        self.ne_box = self._count_box(
            n_e, f"Number of energy points. Defaults to the measured {n_e}; "
                 f"lower it to shrink the result.")
        form.addRow("kx element no.", self.nkx_box)
        form.addRow("ky element no.", self.nky_box)
        form.addRow("E element no.", self.ne_box)

        self.dkx_box = self._resolution_box(0.01, 4, "\u00c5\u207b\u00b9",
                                            "Spacing between kx points.")
        self.dky_box = self._resolution_box(0.01, 4, "\u00c5\u207b\u00b9",
                                            "Spacing between ky points.")
        energy_step = abs(float(contour.E[1] - contour.E[0])) * 1000 if n_e > 1 else 10.0
        self.de_box = self._resolution_box(
            max(round(energy_step, 3), 0.001), 3, "meV",
            "Spacing between energy points. The measurement's own step is "
            f"{energy_step:.4g} meV.")
        form.addRow("kx resolution (\u00c5\u207b\u00b9)", self.dkx_box)
        form.addRow("ky resolution (\u00c5\u207b\u00b9)", self.dky_box)
        form.addRow("E resolution (meV)", self.de_box)

        self.grid_note = QLabel("")
        self.grid_note.setWordWrap(True)
        # Two lines' worth: the note is the only place the other half of the
        # grid setting is visible, and a clipped one would be worse than none.
        self.grid_note.setMinimumHeight(34)
        form.addRow("", self.grid_note)
        layout.addLayout(form)

        for radio in (self.by_count, self.by_resolution):
            radio.toggled.connect(self._update_grid_mode)
        for box in (self.dkx_box, self.dky_box, self.de_box,
                    self.nkx_box, self.nky_box, self.ne_box):
            box.valueChanged.connect(self._update_grid_note)
        self._update_grid_mode()

        self.note = QLabel("")
        self.note.setWordWrap(True)
        layout.addWidget(self.note)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.button(QDialogButtonBox.Ok).setText("Convert")
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

        for box in (self.theta_box, self.phi_box):
            box.valueChanged.connect(self._mark_edited)

    def _default_name(self) -> str:
        """A name for the result that is not already taken.

        The launcher owns the list, so it is asked what is in use; when the
        dialog is exercised without one (in tests), the plain name is fine.
        """
        base = f"{self.contour.filename} [k]"
        taken = set()
        naming = getattr(self.contour, "existing_names", None)
        if callable(naming):
            taken = set(naming())
        if base not in taken:
            return base
        n = 2
        while f"{base} {n}" in taken:
            n += 1
        return f"{base} {n}"

    def output_name(self) -> str:
        return self.name_box.text().strip() or self._default_name()

    # -- the output grid ---------------------------------------------------
    @staticmethod
    def _resolution_box(value, decimals, unit, tip):
        box = QDoubleSpinBox()
        box.setDecimals(decimals)
        box.setRange(10.0 ** -decimals, 1e6)
        box.setSingleStep(10.0 ** -(decimals - 1))
        box.setValue(value)
        box.setToolTip(f"{tip} In {unit}.")
        return box

    def _update_grid_mode(self, *_):
        by_count = self.by_count.isChecked()
        for box in (self.nkx_box, self.nky_box, self.ne_box):
            box.setEnabled(by_count)
        for box in (self.dkx_box, self.dky_box, self.de_box):
            box.setEnabled(not by_count)
        self._update_grid_note()

    def k_extent(self):
        """The (kx, ky) box the current settings would produce, so a
        resolution can be turned into a point count. Returns None if the
        settings do not describe a conversion at all."""
        from tools.kspace import k_extent
        try:
            return k_extent(self.contour.angle_defl, self.contour.angle_slit,
                            self.contour.E,
                            theta_offset_deg=self.theta_box.value(),
                            phi_offset_deg=self.phi_box.value(),
                            azimuth_deg=self.azimuth_box.value(),
                            energy_offset_eV=self.energy_offset_box.value())
        except Exception:
            return None

    def grid_counts(self):
        """(n_kx, n_ky, n_E), whichever way the user asked for them."""
        if self.by_count.isChecked():
            return (int(self.nkx_box.value()), int(self.nky_box.value()),
                    int(self.ne_box.value()))
        extent = self.k_extent()
        if extent is None:
            raise ValueError("the conversion settings do not define a k range yet")
        kx_lo, kx_hi, ky_lo, ky_hi = extent
        n_kx = int(round((kx_hi - kx_lo) / self.dkx_box.value())) + 1
        n_ky = int(round((ky_hi - ky_lo) / self.dky_box.value())) + 1
        energy = np.asarray(self.contour.E, dtype=float)
        span_meV = abs(float(energy[-1] - energy[0])) * 1000.0
        n_e = int(round(span_meV / self.de_box.value())) + 1
        return max(n_kx, 2), max(n_ky, 2), max(n_e, 2)

    def _update_grid_note(self, *_):
        try:
            n_kx, n_ky, n_e = self.grid_counts()
        except Exception as exc:
            self.grid_note.setText(str(exc))
            return
        if self.by_count.isChecked():
            extent = self.k_extent()
            if extent is None:
                self.grid_note.setText("")
                return
            kx_lo, kx_hi, ky_lo, ky_hi = extent
            energy = np.asarray(self.contour.E, dtype=float)
            span_meV = abs(float(energy[-1] - energy[0])) * 1000.0
            self.grid_note.setText(
                f"\u2192 {(kx_hi - kx_lo) / max(n_kx - 1, 1):.4g} / "
                f"{(ky_hi - ky_lo) / max(n_ky - 1, 1):.4g} \u00c5\u207b\u00b9 and "
                f"{span_meV / max(n_e - 1, 1):.4g} meV per point")
        else:
            self.grid_note.setText(f"\u2192 {n_kx} x {n_ky} x {n_e} points")

    # -- sample rotation from the contour ---------------------------------
    def _toggle_rotation_picking(self, on: bool):
        if on:
            self._rotation_points = []
            self.contour.view.clear_overlays()
            self.contour.view.start_point_picking(self._add_rotation_point)
            self.rotation_note.setText(
                "Click one or two points on the contour along the direction "
                "that should end up vertical.")
        else:
            self.contour.view.stop_point_picking()
            self.contour.view.clear_overlays()
            if not self._rotation_points:
                self.rotation_note.setText("")
            # The readout cursor's own picking is not restored: it never used
            # clicks, only dragging, so nothing was taken away.

    def origin_point(self):
        """Where k = (0, 0) currently sits on the contour, in its own angles.

        That is the theta/phi offsets -- the point picked with the readout
        cursor -- not the contour's (0, 0). A single rotation point is
        measured from here, so the direction it defines is the one that will
        actually run through the origin of the converted map.
        """
        return float(self.theta_box.value()), float(self.phi_box.value())

    def _add_rotation_point(self, x: float, y: float):
        from tools.kspace import points_to_azimuth
        self._rotation_points.append((x, y))
        if len(self._rotation_points) > 2:
            self._rotation_points = self._rotation_points[-1:]
        xs = [p[0] for p in self._rotation_points]
        ys = [p[1] for p in self._rotation_points]
        self.contour.view.show_picked_points(xs, ys)
        if len(self._rotation_points) == 1:
            ox, oy = self.origin_point()
            pair = [(ox, oy), (xs[0], ys[0])]
        else:
            pair = list(self._rotation_points)
        try:
            rotation = points_to_azimuth(pair)
        except ValueError as exc:
            self.rotation_note.setText(str(exc))
            return
        self.azimuth_box.setValue(rotation)
        if len(self._rotation_points) == 1:
            ox, oy = pair[0]
            self.contour.view.show_overlay_path([ox, xs[0]], [oy, ys[0]])
            self.rotation_note.setText(
                f"One point: the line from the k-space origin "
                f"({ox:.3g}, {oy:.3g}) to it is taken as the direction "
                f"\u2192 {rotation:.3f}\u00b0. Click a second point to use the "
                f"line between the two instead.")
        else:
            self.contour.view.show_overlay_path(xs, ys)
            self.rotation_note.setText(
                f"Two points \u2192 {rotation:.3f}\u00b0 stands that direction "
                f"vertical (along ky).")
        self._update_grid_note()

    def _angle_box(self, tip):
        box = QDoubleSpinBox()
        box.setDecimals(4)
        box.setRange(-180.0, 180.0)
        box.setSingleStep(0.1)
        box.setToolTip(tip)
        return box

    @staticmethod
    def _count_box(value, tip):
        box = QSpinBox()
        box.setRange(2, 20000)
        box.setValue(int(value))
        box.setToolTip(tip)
        return box

    def _mark_edited(self, *_):
        self._user_edited = True

    def follow_cursor(self):
        """Track the contour's readout cursor, until the user types a value
        of their own."""
        if not self._user_edited:
            self.read_cursor(mark=False)

    def read_cursor(self, mark=True):
        """Copy the readout cursor's position into the offset boxes.

        There is no button for this any more: the boxes follow the cursor by
        themselves for as long as the user has not typed a value of their
        own, which is every moment the button would have been useful in.
        """
        position = self.contour.readout_position()
        if position is None:
            self.note.setText("No readout cursor on the contour yet -- "
                              "right-click it and choose \"Readout cursor\".")
            return
        theta, phi = position
        for box, value in ((self.theta_box, theta), (self.phi_box, phi)):
            box.blockSignals(True)
            box.setValue(float(value))
            box.blockSignals(False)
        self.note.setText(f"Origin taken from the cursor: "
                          f"theta {theta:.4g} deg, phi {phi:.4g} deg.")
        if mark:
            self._user_edited = True

    def settings(self) -> dict:
        n_kx, n_ky, n_energy = self.grid_counts()
        return {
            "theta_offset_deg": self.theta_box.value(),
            "phi_offset_deg": self.phi_box.value(),
            "azimuth_deg": self.azimuth_box.value(),
            "energy_offset_eV": self.energy_offset_box.value(),
            "n_kx": n_kx,
            "n_ky": n_ky,
            "n_energy": n_energy,
        }

    def closeEvent(self, event):
        # Picking must not outlive the dialog: a contour that still swallowed
        # clicks afterwards would be a puzzle with no visible cause.
        if self.rotation_button.isChecked():
            self.rotation_button.setChecked(False)
        else:
            self.contour.view.stop_point_picking()
            self.contour.view.clear_overlays()
        super().closeEvent(event)


class _OrbitViewBox(pg.ViewBox):
    """A ViewBox whose left-drag turns a 3-D scene around instead of panning
    the 2-D picture it is really showing.

    The Brillouin-zone preview is a hand-projected wireframe rather than a
    real 3-D scene graph (see :class:`CutPlanePreviewDialog` for why), so
    "rotate the view" has to mean "change the two angles the projection is
    computed from and redraw" -- this is where the drag gets turned into
    those two numbers. Right-drag is left alone, so the inherited zoom still
    works.
    """

    orbited = pyqtSignal(float, float)      # delta azimuth, delta elevation

    def mouseDragEvent(self, ev, axis=None):
        if ev.button() != Qt.LeftButton:
            super().mouseDragEvent(ev, axis)
            return
        ev.accept()
        delta = ev.screenPos() - ev.lastScreenPos()
        self.orbited.emit(-0.4 * delta.x(), 0.4 * delta.y())


class CutPlanePreviewDialog(QDialog):
    """A turnable 3-D view of the Brillouin zone with the cutting plane in
    it, so that "(h k l) = (1 0 1), offset 0.2" can be *seen* rather than
    imagined.

    Drawn by projecting the zone's own face loops onto the screen by hand --
    an orthographic camera in a dozen lines -- rather than through
    ``pyqtgraph.opengl``, which would add PyOpenGL to a program that
    currently installs with pip and no compiler, for one preview window. The
    faces pointing away from the camera are drawn faint and dashed and those
    facing it solid, which gives the shape enough depth to read without any
    hidden-surface machinery.

    What is drawn: the zone (or the irreducible wedge, whichever the dialog
    is set to), the plane as a dashed square, the polygon where the plane
    actually cuts the zone -- filled, and it is this polygon that becomes
    the overlay on the contour -- the plane's normal as an arrow from Gamma,
    and the three reciprocal-lattice vectors for orientation.
    """

    def __init__(self, parent: "BrillouinZoneDialog"):
        super().__init__(parent)
        self.bz_dialog = parent
        self.setModal(False)
        self.setWindowTitle("Cut plane in 3D")
        self.resize(560, 620)
        self._azimuth, self._elevation = 35.0, 22.0

        layout = QVBoxLayout(self)
        intro = QLabel(
            "Drag to turn the zone; right-drag to zoom. The filled polygon is "
            "the cut that goes onto the contour.")
        intro.setWordWrap(True)
        layout.addWidget(intro)

        self.view_box = _OrbitViewBox(lockAspect=True, enableMenu=False)
        self.view_box.orbited.connect(self._orbit)
        self.plot = pg.PlotWidget(viewBox=self.view_box)
        self.plot.hideAxis("left")
        self.plot.hideAxis("bottom")
        layout.addWidget(self.plot, stretch=1)

        controls = QHBoxLayout()
        self.azimuth_box = self._angle_box(self._azimuth, -360.0, 360.0)
        self.elevation_box = self._angle_box(self._elevation, -89.0, 89.0)
        for label, box in (("Azimuth", self.azimuth_box),
                           ("Elevation", self.elevation_box)):
            controls.addWidget(QLabel(label))
            controls.addWidget(box)
            box.valueChanged.connect(self._angles_typed)
        reset = QPushButton("Reset view")
        reset.clicked.connect(self._reset_view)
        controls.addWidget(reset)
        controls.addStretch(1)
        layout.addLayout(controls)

        self.note = QLabel("")
        self.note.setWordWrap(True)
        self.note.setMinimumHeight(34)
        layout.addWidget(self.note)

        buttons = QDialogButtonBox(QDialogButtonBox.Close)
        buttons.rejected.connect(self.close)
        layout.addWidget(buttons)

        self._frame_the_cut()
        for box, value in ((self.azimuth_box, self._azimuth),
                           (self.elevation_box, self._elevation)):
            box.blockSignals(True)
            box.setValue(value)
            box.blockSignals(False)
        self.refresh()

    @staticmethod
    def _angle_box(value, low, high):
        box = QDoubleSpinBox()
        box.setDecimals(1)
        box.setRange(low, high)
        box.setSingleStep(5.0)
        box.setValue(value)
        box.setKeyboardTracking(False)
        return box

    def _orbit(self, d_azimuth: float, d_elevation: float):
        self._azimuth = (self._azimuth + d_azimuth) % 360.0
        self._elevation = float(np.clip(self._elevation + d_elevation, -89.0, 89.0))
        for box, value in ((self.azimuth_box, self._azimuth),
                           (self.elevation_box, self._elevation)):
            box.blockSignals(True)
            box.setValue(value)
            box.blockSignals(False)
        self.refresh()

    def _angles_typed(self, *_):
        self._azimuth = self.azimuth_box.value()
        self._elevation = self.elevation_box.value()
        self.refresh()

    def _reset_view(self):
        self._frame_the_cut()
        for box, value in ((self.azimuth_box, self._azimuth),
                           (self.elevation_box, self._elevation)):
            box.blockSignals(True)
            box.setValue(value)
            box.blockSignals(False)
        self.plot.enableAutoRange()
        self.refresh()

    def _frame_the_cut(self):
        """Point the camera about 55 degrees away from the cut plane's own
        normal.

        A fixed starting view is fine until the plane's normal happens to
        point at the camera -- a (1 1 1) cut of a cubic crystal very nearly
        does -- and then the cut is seen face-on or edge-on, which is the
        one thing this window exists to avoid. Tilting off the normal by a
        fixed angle keeps every orientation legible without the user having
        to drag the view around first.
        """
        try:
            *_, normal, _point, _u, _v = self.bz_dialog.cut_geometry()
        except Exception:
            self._azimuth, self._elevation = 35.0, 22.0
            return
        normal_azimuth = float(np.degrees(np.arctan2(normal[1], normal[0])))
        normal_elevation = float(np.degrees(np.arcsin(np.clip(normal[2], -1.0, 1.0))))
        tilt = -55.0 if normal_elevation >= 0.0 else 55.0
        self._azimuth = (normal_azimuth + 25.0) % 360.0
        self._elevation = float(np.clip(normal_elevation + tilt, -89.0, 89.0))

    # -- projection --------------------------------------------------------
    def _camera(self):
        azimuth = np.radians(self._azimuth)
        elevation = np.radians(self._elevation)
        eye = np.array([np.cos(elevation) * np.cos(azimuth),
                        np.cos(elevation) * np.sin(azimuth),
                        np.sin(elevation)])
        up_reference = np.array([0.0, 0.0, 1.0])
        right = np.cross(up_reference, eye)
        if np.linalg.norm(right) < 1e-9:       # looking straight down z
            right = np.array([1.0, 0.0, 0.0])
        right = right / np.linalg.norm(right)
        return eye, right, np.cross(eye, right)

    def _project(self, points, camera=None):
        eye, right, up = camera or self._camera()
        points = np.atleast_2d(np.asarray(points, dtype=float))
        return points @ right, points @ up

    # -- drawing -----------------------------------------------------------
    def refresh(self):
        """Rebuild the picture from the dialog's current settings."""
        self.plot.clear()
        try:
            scene = self.bz_dialog.cut_geometry()
        except Exception as exc:
            self.note.setText(f"Nothing to show: {exc}")
            return

        zone, _b, b_conventional, normal, point, u, v = scene

        # A viewpoint chosen for one cut plane is not necessarily a useful
        # one for the next: swinging the normal round to (1 1 1) can leave
        # the old camera staring straight down it, which shows nothing. So
        # the view re-frames itself when the *normal* changes, and otherwise
        # leaves alone whatever the user has dragged it to.
        previous = getattr(self, "_framed_normal", None)
        if previous is None or not np.allclose(previous, normal, atol=1e-9):
            self._framed_normal = np.array(normal, dtype=float)
            self._frame_the_cut()
            for box, value in ((self.azimuth_box, self._azimuth),
                               (self.elevation_box, self._elevation)):
                box.blockSignals(True)
                box.setValue(value)
                box.blockSignals(False)
            self.plot.enableAutoRange()

        camera = self._camera()
        extent = float(np.max(np.linalg.norm(zone.vertices, axis=1))) or 1.0

        # the zone itself, front faces solid and back faces faint
        eye = camera[0]
        for loop in zone.faces:
            centre = loop[:-1].mean(axis=0)
            facing = float(np.dot(centre, eye)) >= 0.0
            xs, ys = self._project(loop, camera)
            pen = pg.mkPen("#7fd4ff" if facing else "#2f5a70",
                           width=1.6 if facing else 1.0,
                           style=Qt.SolidLine if facing else Qt.DashLine)
            self.plot.addItem(pg.PlotDataItem(xs, ys, pen=pen))

        # the plane, as a dashed square a little larger than the zone
        from tools.bz3d import plane_basis, cut_points_3d
        u_axis, v_axis = plane_basis(normal)
        half = 1.05 * extent
        corners = np.array([point + su * half * u_axis + sv * half * v_axis
                            for su, sv in ((1, 1), (-1, 1), (-1, -1), (1, -1), (1, 1))])
        xs, ys = self._project(corners, camera)
        self.plot.addItem(pg.PlotDataItem(
            xs, ys, pen=pg.mkPen("#ff3ce0", width=1.0, style=Qt.DashLine)))

        # where it actually cuts, filled: this polygon is the overlay
        if u is not None:
            cut = cut_points_3d(u, v, normal, point)
            xs, ys = self._project(cut, camera)
            path = QPainterPath()
            path.moveTo(float(xs[0]), float(ys[0]))
            for x, y in zip(xs[1:], ys[1:]):
                path.lineTo(float(x), float(y))
            path.closeSubpath()
            patch = QGraphicsPathItem(path)
            patch.setBrush(pg.mkBrush(255, 60, 224, 70))
            patch.setPen(pg.mkPen("#ff3ce0", width=2.0))
            self.plot.addItem(patch)

        # the normal from Gamma to the plane, and the reciprocal axes
        tip = point + normal * 0.55 * extent
        xs, ys = self._project(np.array([np.zeros(3), point, tip]), camera)
        self.plot.addItem(pg.PlotDataItem(
            xs, ys, pen=pg.mkPen("#ffd166", width=2.0)))
        self.plot.addItem(pg.PlotDataItem(
            [xs[-1]], [ys[-1]], symbol="o", symbolSize=7,
            symbolBrush="#ffd166", pen=None))
        # the *conventional* reciprocal axes, because those are the ones the
        # (h k l) boxes are read against (see BrillouinZoneDialog._zone_3d)
        for vector, name, colour in zip(b_conventional, ("a*", "b*", "c*"),
                                        ("#ff6b6b", "#7bd88f", "#8fa2ff")):
            scaled = vector / np.linalg.norm(vector) * 0.5 * extent
            xs, ys = self._project(np.array([np.zeros(3), scaled]), camera)
            self.plot.addItem(pg.PlotDataItem(xs, ys, pen=pg.mkPen(colour, width=1.4)))
            label = pg.TextItem(name, color=colour, anchor=(0.5, 0.5))
            label.setPos(float(xs[-1]), float(ys[-1]))
            self.plot.addItem(label)
        gamma = pg.TextItem("Γ", color="#dddddd", anchor=(1.2, 1.2))
        gamma.setPos(*[float(c[0]) for c in self._project(np.zeros(3), camera)])
        self.plot.addItem(gamma)

        self.note.setText(self.bz_dialog.cut_description())

    def closeEvent(self, event):
        self.bz_dialog.preview_closed()
        super().closeEvent(event)


class BrillouinZoneDialog(QDialog):
    """Overlay a Brillouin zone -- one crystal, or two twisted/mismatched
    layers' moire zone -- on a k-space constant-energy contour.

    Modeless, like the other contour tools that read points off the picture
    (:class:`KConversionDialog`, :class:`ArbitraryCutDialog`): the "Pick
    direction from contour..." button needs the contour to stay live behind
    it. Nothing is drawn until **Plot** is pressed, though: the settings
    here describe a crystal, and a half-typed lattice constant is not one,
    so redrawing on every keystroke would spend its time showing zones
    nobody asked for. What does update as you type is the reference text --
    which Bravais lattice the space group implies, how far along its normal
    the cut plane has been slid -- and the 3-D preview window, if it is
    open, since looking at the cut is the whole reason it exists.

    The drawn zone is deliberately left in place when this dialog is closed
    (see :meth:`_InteractiveImageBase.set_bz_layer`); "Clear overlay" is the
    explicit way to take it off again.

    Only ever opened on a k-space map (:attr:`ContourWindow.bz_button` is
    hidden otherwise, matching ``kconv_button``'s own gating the other way
    round): a Brillouin zone is a statement about the crystal's momentum-
    space periodicity, so it only means something once the contour's axes
    are actually momenta.

    Two modes and two kinds of zone
    -------------------------------
    **3D crystal** starts from the space-group number, and the space group
    then decides which of the six cell parameters are the user's to choose:
    a cubic group offers one length and no angles, a hexagonal group two
    lengths and no angles, P1 all six. The rest are disabled and follow
    along, so a cell that contradicts its own space group cannot be typed in
    the first place (which is a better answer than the warning this used to
    print underneath it). A rhombohedral group is offered in hexagonal axes
    only -- see :func:`tools.lattice.free_parameters`.

    **2D layer** is three numbers for a surface lattice with no third axis,
    and is also where the moire mode lives: a moire pattern is two stacked
    *layers*, so asking for it from a 3-D crystal would mean guessing which
    2-D sublattice the chosen cut plane exposes, which is a different (and
    unrequested) calculation.

    Either mode draws the **conventional** zone -- the first Brillouin zone,
    i.e. the Wigner-Seitz cell of the reciprocal lattice -- or the
    **irreducible** one, the wedge of it that the crystal's own symmetry
    repeats into the whole zone (:func:`tools.bz3d.irreducible_wedge`,
    :func:`tools.bz2d.irreducible_cell_2d`). In 3-D the symmetry comes from
    the space group's point group plus inversion (the Laue class, since time
    reversal gives ``E(k) = E(-k)`` even in a non-centrosymmetric crystal);
    in 2-D there is no space group to ask, so the plane lattice's own
    holohedry is used and the dialog says so, because a layer whose basis is
    less symmetric than its lattice has a larger wedge than the one drawn.

    Only the first zone is ever drawn, repeated by translation across the
    field of view if asked -- never a second- or higher-order zone.
    """

    #: flat "bz.<key>" keys in scan.info, matching the app's own
    #: "kconv.*"/"arbcut.*" convention -- save_dataset only ever writes
    #: scalars (see loader.nxs_file.save_dataset), so settings are flattened
    #: rather than kept as one nested dict.
    _PREFIX = "bz."

    def __init__(self, contour: "ContourWindow"):
        super().__init__(contour)
        self.contour = contour
        self.setModal(False)
        self.setWindowTitle("Brillouin zone")
        self._direction_points = []
        self._zone_cache = None
        self._preview = None

        layout = QVBoxLayout(self)

        # -- mode ---------------------------------------------------------
        mode_row = QWidget()
        mode_box = QHBoxLayout(mode_row)
        mode_box.setContentsMargins(0, 0, 0, 0)
        self.mode_3d = QRadioButton("3D crystal (space group)")
        self.mode_2d = QRadioButton("2D layer (surface lattice)")
        self.mode_3d.setChecked(True)
        mode_box.addWidget(self.mode_3d)
        mode_box.addWidget(self.mode_2d)
        mode_box.addStretch(1)
        layout.addWidget(mode_row)

        # -- 3D crystal: the space group first, then what it leaves free ---
        self.group_3d = QGroupBox("Crystal")
        form3d = QFormLayout(self.group_3d)
        self.spacegroup_box = QSpinBox()
        self.spacegroup_box.setRange(1, 230)
        self.spacegroup_box.setValue(1)
        self.spacegroup_box.setKeyboardTracking(False)
        self.spacegroup_box.setToolTip(
            "International Tables (1-230) space-group number. Choose this "
            "first: it fixes which cell parameters below are yours to set.")
        form3d.addRow("Space group no.", self.spacegroup_box)
        self.spacegroup_label = QLabel("")
        self.spacegroup_label.setWordWrap(True)
        self.spacegroup_label.setMinimumHeight(34)
        form3d.addRow("", self.spacegroup_label)

        self.a_box = self._length_box(3.0)
        self.b_box = self._length_box(3.0)
        self.c_box = self._length_box(3.0)
        self.alpha_box = self._angle_box_180(90.0)
        self.beta_box = self._angle_box_180(90.0)
        self.gamma_box = self._angle_box_180(90.0)
        for text, box in (("a (Å)", self.a_box), ("b (Å)", self.b_box),
                          ("c (Å)", self.c_box),
                          ("alpha (deg)", self.alpha_box),
                          ("beta (deg)", self.beta_box),
                          ("gamma (deg)", self.gamma_box)):
            box.setToolTip("Greyed out means the space group fixes this one.")
            form3d.addRow(text, box)
        layout.addWidget(self.group_3d)

        # -- 3D crystal: the cut plane -------------------------------------
        self.group_cut = QGroupBox("Cut plane")
        form_cut = QFormLayout(self.group_cut)
        normal_row = QWidget()
        normal_box = QHBoxLayout(normal_row)
        normal_box.setContentsMargins(0, 0, 0, 0)
        self.h_box = self._miller_box(0)
        self.k_box = self._miller_box(0)
        self.l_box = self._miller_box(1)
        for text, box in (("h", self.h_box), ("k", self.k_box), ("l", self.l_box)):
            normal_box.addWidget(QLabel(text))
            normal_box.addWidget(box)
        normal_box.addStretch(1)
        normal_row.setToolTip(
            "The plane's normal as Miller indices of the reciprocal lattice: "
            "h*b1 + k*b2 + l*b3. (0, 0, 1) is the usual surface-normal cut.")
        form_cut.addRow("Normal (hkl)", normal_row)

        self.offset_box = QDoubleSpinBox()
        self.offset_box.setDecimals(4)
        self.offset_box.setRange(-20.0, 20.0)
        self.offset_box.setSingleStep(0.01)
        self.offset_box.setKeyboardTracking(False)
        self.offset_box.setToolTip(
            "How far to slide the plane along its own normal, from Gamma. "
            "0 cuts through Gamma. This is a geometric offset in A^-1: "
            "turning a photon energy into a kz that belongs here needs an "
            "inner potential, which is the separate kz feature's job.")
        form_cut.addRow("Offset along normal (Å⁻¹)", self.offset_box)
        self.offset_label = QLabel("")
        self.offset_label.setWordWrap(True)
        self.offset_label.setMinimumHeight(34)
        form_cut.addRow("", self.offset_label)

        self.preview_button = QPushButton("Show cut plane in 3D...")
        self.preview_button.setToolTip(
            "Open a turnable 3-D view of the zone with this plane in it.")
        self.preview_button.clicked.connect(self.open_preview)
        form_cut.addRow("", self.preview_button)
        layout.addWidget(self.group_cut)

        # -- 2D layer -------------------------------------------------------
        self.group_2d = QGroupBox("Layer lattice")
        form2d = QFormLayout(self.group_2d)
        self.a2_box = self._length_box(3.0)
        self.b2_box = self._length_box(3.0)
        self.gamma2_box = self._angle_box_180(120.0)
        form2d.addRow("a (Å)", self.a2_box)
        form2d.addRow("b (Å)", self.b2_box)
        form2d.addRow("gamma (deg)", self.gamma2_box)
        layout.addWidget(self.group_2d)

        # -- which zone ------------------------------------------------------
        zone_group = QGroupBox("Zone")
        zone_form = QFormLayout(zone_group)
        zone_row = QWidget()
        zone_box = QHBoxLayout(zone_row)
        zone_box.setContentsMargins(0, 0, 0, 0)
        self.zone_conventional = QRadioButton("Conventional")
        self.zone_irreducible = QRadioButton("Irreducible")
        self.zone_conventional.setChecked(True)
        self.zone_conventional.setToolTip(
            "The first Brillouin zone: the Wigner-Seitz cell of the "
            "reciprocal lattice.")
        self.zone_irreducible.setToolTip(
            "The wedge of the zone that the crystal's symmetry repeats into "
            "the whole of it.")
        zone_box.addWidget(self.zone_conventional)
        zone_box.addWidget(self.zone_irreducible)
        zone_box.addStretch(1)
        zone_form.addRow("", zone_row)
        self.zone_label = QLabel("")
        self.zone_label.setWordWrap(True)
        self.zone_label.setMinimumHeight(34)
        zone_form.addRow("", self.zone_label)
        layout.addWidget(zone_group)

        # -- orientation ----------------------------------------------------
        orient_group = QGroupBox("Orientation")
        orient_form = QFormLayout(orient_group)
        self.azimuth_box = QDoubleSpinBox()
        self.azimuth_box.setDecimals(3)
        self.azimuth_box.setRange(-360.0, 360.0)
        self.azimuth_box.setSingleStep(1.0)
        self.azimuth_box.setKeyboardTracking(False)
        self.azimuth_box.setToolTip(
            "Rotates the drawn zone (counter-clockwise) about k = (0, 0) to "
            "line it up with the data -- the zone's own axes have no fixed "
            "relation to kx/ky until this is set.")
        orient_form.addRow("Azimuth (deg)", self.azimuth_box)
        self.direction_button = QPushButton("Pick direction from contour...")
        self.direction_button.setCheckable(True)
        self.direction_button.setToolTip(
            "Click one or two points on the contour along a direction whose "
            "crystallographic azimuth you already know. The angle shown is "
            "only reported, not applied -- set Azimuth above by hand to "
            "whatever value should put that direction there.")
        self.direction_button.toggled.connect(self._toggle_direction_picking)
        orient_form.addRow("", self.direction_button)
        self.direction_note = QLabel("")
        self.direction_note.setWordWrap(True)
        self.direction_note.setMinimumHeight(34)
        orient_form.addRow("", self.direction_note)
        self.tile_cb = QCheckBox("Repeat across the field of view")
        self.tile_cb.setChecked(True)
        self.tile_cb.setToolTip(
            "Draw every copy of the same zone that overlaps the visible k "
            "range, not only the one at the origin -- never a higher-order "
            "zone, only repeats of this one.")
        orient_form.addRow("", self.tile_cb)
        layout.addWidget(orient_group)

        # -- moire (2D layers only) -------------------------------------------
        self.moire_cb = QCheckBox("Moire (two twisted/mismatched layers)")
        self.moire_cb.setToolTip(
            "Uses the layer lattice above as the top layer and asks for a "
            "bottom layer below. Layers only: a moire pattern is two stacked "
            "2-D lattices, so this is not offered for a 3-D crystal.")
        layout.addWidget(self.moire_cb)
        self.moire_group = QGroupBox("Second (bottom) layer")
        moire_form = QFormLayout(self.moire_group)
        self.moire_same_cb = QCheckBox("Same lattice, twisted")
        self.moire_same_cb.setChecked(True)
        self.moire_same_cb.setToolTip(
            "On: the bottom layer is the same a/b/gamma as the top, just "
            "rotated by Twist -- the common case (e.g. twisted bilayer "
            "graphene). Off: give the bottom layer its own a/b/gamma below.")
        moire_form.addRow("", self.moire_same_cb)
        self.a2b_box = self._length_box(3.0)
        self.b2b_box = self._length_box(3.0)
        self.gamma2b_box = self._angle_box_180(120.0)
        moire_form.addRow("bottom a (Å)", self.a2b_box)
        moire_form.addRow("bottom b (Å)", self.b2b_box)
        moire_form.addRow("bottom gamma (deg)", self.gamma2b_box)
        self.twist_box = QDoubleSpinBox()
        self.twist_box.setDecimals(3)
        self.twist_box.setRange(-180.0, 180.0)
        self.twist_box.setSingleStep(0.1)
        self.twist_box.setKeyboardTracking(False)
        self.twist_box.setToolTip("Bottom layer's rotation relative to the top.")
        moire_form.addRow("Twist (deg)", self.twist_box)
        self.show_top_cb = QCheckBox("Show top-layer zone")
        self.show_top_cb.setChecked(True)
        self.show_bottom_cb = QCheckBox("Show bottom-layer zone")
        self.show_bottom_cb.setChecked(True)
        self.show_moire_cb = QCheckBox("Show moire zone")
        self.show_moire_cb.setChecked(True)
        for box in (self.show_top_cb, self.show_bottom_cb, self.show_moire_cb):
            moire_form.addRow("", box)
        layout.addWidget(self.moire_group)

        self.note = QLabel("")
        self.note.setWordWrap(True)
        self.note.setMinimumHeight(34)
        layout.addWidget(self.note)

        buttons_row = QHBoxLayout()
        self.plot_button = QPushButton("Plot")
        self.plot_button.setDefault(True)
        self.plot_button.setToolTip(
            "Work out the zone from these settings and draw it on the "
            "contour. Nothing is drawn until this is pressed.")
        self.plot_button.clicked.connect(self.plot)
        self.clear_button = QPushButton("Clear overlay")
        self.clear_button.setToolTip("Remove the drawn zone(s) from the contour.")
        self.clear_button.clicked.connect(self._clear_overlay)
        buttons_row.addWidget(self.plot_button)
        buttons_row.addWidget(self.clear_button)
        buttons_row.addStretch(1)
        layout.addLayout(buttons_row)
        buttons = QDialogButtonBox(QDialogButtonBox.Close)
        buttons.rejected.connect(self.close)
        layout.addWidget(buttons)

        # -- wiring ------------------------------------------------------------
        # Note what is *not* wired here: nothing recomputes the overlay. These
        # only keep the dialog's own reference text (and the 3-D preview, if
        # it is open) honest while the user types.
        self.spacegroup_box.valueChanged.connect(self._space_group_changed)
        self.a_box.valueChanged.connect(self._sync_mirrored_lengths)
        self.mode_3d.toggled.connect(self._toggle_mode)
        self.moire_cb.toggled.connect(self._toggle_moire)
        self.moire_same_cb.toggled.connect(self._toggle_moire_same)
        for box in (self.a_box, self.b_box, self.c_box, self.alpha_box,
                    self.beta_box, self.gamma_box, self.h_box, self.k_box,
                    self.l_box, self.offset_box):
            box.valueChanged.connect(self._refresh_info)
        for box in (self.zone_conventional, self.zone_irreducible):
            box.toggled.connect(self._refresh_info)

        self._load_saved_settings()
        self._apply_constraints()
        self._toggle_mode()
        self._toggle_moire()
        self._toggle_moire_same()
        self._refresh_info()

    # -- small widget factories --------------------------------------------
    @staticmethod
    def _length_box(value):
        box = QDoubleSpinBox()
        box.setDecimals(5)
        box.setRange(0.01, 1000.0)
        box.setSingleStep(0.01)
        box.setValue(value)
        box.setKeyboardTracking(False)
        return box

    @staticmethod
    def _angle_box_180(value):
        box = QDoubleSpinBox()
        box.setDecimals(3)
        box.setRange(1.0, 179.0)
        box.setSingleStep(0.5)
        box.setValue(value)
        box.setKeyboardTracking(False)
        return box

    @staticmethod
    def _miller_box(value):
        box = QSpinBox()
        box.setRange(-20, 20)
        box.setValue(value)
        box.setKeyboardTracking(False)
        return box

    # -- what the space group leaves free -----------------------------------
    def _length_boxes(self):
        return {"a": self.a_box, "b": self.b_box, "c": self.c_box}

    def _angle_boxes(self):
        return {"alpha": self.alpha_box, "beta": self.beta_box,
                "gamma": self.gamma_box}

    def _apply_constraints(self):
        """Disable and drive the cell parameters the space group fixes.

        The alternative -- letting anything be typed and printing a warning
        afterwards -- puts the user in the position of having to know the
        constraint anyway, so the boxes that are not theirs to set are
        switched off and follow the ones that are.
        """
        from tools.lattice import free_parameters
        constraints = free_parameters(int(self.spacegroup_box.value()))
        self._mirrors = dict(constraints.mirrors)

        lengths = self._length_boxes()
        for name, box in lengths.items():
            source = constraints.mirrors.get(name)
            box.setEnabled(source is None)
        for name, box in self._angle_boxes().items():
            fixed = constraints.fixed_angles.get(name)
            box.setEnabled(fixed is None)
            if fixed is not None:
                box.blockSignals(True)
                box.setValue(float(fixed))
                box.blockSignals(False)
        self._sync_mirrored_lengths()

    def _sync_mirrored_lengths(self, *_):
        lengths = self._length_boxes()
        for target, source in getattr(self, "_mirrors", {}).items():
            box = lengths[target]
            box.blockSignals(True)
            box.setValue(lengths[source].value())
            box.blockSignals(False)

    def _space_group_changed(self, *_):
        self._apply_constraints()
        self._refresh_info()

    # -- mode/section visibility ---------------------------------------------
    def _toggle_mode(self, *_):
        is_3d = self.mode_3d.isChecked()
        self.group_3d.setVisible(is_3d)
        self.group_cut.setVisible(is_3d)
        self.group_2d.setVisible(not is_3d)
        # A moire pattern is two stacked layers, so it belongs to the 2-D
        # mode only; hiding rather than disabling keeps the 3-D panel from
        # implying there is a 3-D moire calculation being withheld.
        self.moire_cb.setVisible(not is_3d)
        self.moire_group.setVisible(not is_3d and self.moire_cb.isChecked())
        if is_3d is False and self._preview is not None:
            self._preview.close()
        self._refresh_info()

    def _toggle_moire(self, *_):
        self.moire_group.setVisible(self.moire_cb.isChecked()
                                    and not self.mode_3d.isChecked())

    def _toggle_moire_same(self, *_):
        same = self.moire_same_cb.isChecked()
        for box in (self.a2b_box, self.b2b_box, self.gamma2b_box):
            box.setEnabled(not same)

    def zone_kind(self) -> str:
        return "irreducible" if self.zone_irreducible.isChecked() else "conventional"

    # -- reference text, and the 3-D preview ---------------------------------
    def _refresh_info(self, *_):
        """Keep the labels (and the preview window) current. Draws nothing on
        the contour -- that only happens on Plot."""
        from tools.spacegroups import spacegroup_info
        from tools.lattice import free_parameters

        if self.mode_3d.isChecked():
            try:
                number = int(self.spacegroup_box.value())
                info = spacegroup_info(number)
                self.spacegroup_label.setText(
                    f"{info.symbol} (#{info.number}) — {info.crystal_system}, "
                    f"Bravais {info.bravais_symbol}, point group "
                    f"{info.point_group}\n{free_parameters(number).note}")
            except Exception as exc:
                self.spacegroup_label.setText(str(exc))
            self.offset_label.setText(self.cut_description())
        else:
            self.spacegroup_label.setText("")
            self.offset_label.setText("")

        self.zone_label.setText(self._zone_description())
        if self._preview is not None and self._preview.isVisible():
            self._preview.refresh()

    def _zone_description(self) -> str:
        if self.zone_kind() == "conventional":
            return ("The first Brillouin zone: the Wigner-Seitz cell of the "
                    "reciprocal lattice.")
        if self.mode_3d.isChecked():
            try:
                from tools.lattice import point_group_operations
                order = len(point_group_operations(self.lattice_params()))
            except Exception as exc:
                return f"Irreducible zone unavailable: {exc}"
            return (f"1/{order} of the zone: the wedge the Laue class repeats "
                    "into the whole of it (the point group plus inversion, "
                    "since time reversal makes E(k) = E(-k)).")
        return ("The wedge of the zone repeated by the plane lattice's own "
                "symmetry. With no space group for a layer, the lattice's "
                "holohedry is assumed -- a layer whose basis is less "
                "symmetric than its lattice has a larger wedge than this.")

    def open_preview(self):
        """The 3-D view of the cut plane, in its own window."""
        if self._preview is not None and self._preview.isVisible():
            self._preview.raise_()
            self._preview.activateWindow()
            return self._preview
        self._preview = CutPlanePreviewDialog(self)
        self._preview.show()
        self._preview.raise_()
        return self._preview

    def preview_closed(self):
        self._preview = None

    # -- geometry ------------------------------------------------------------
    @staticmethod
    def _rotate_xy(x, y, deg):
        r = np.radians(deg)
        c, s = np.cos(r), np.sin(r)
        x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
        return x * c - y * s, x * s + y * c

    def _view_window(self):
        """A square half-extent ``(-half, half)`` for both axes, ``half``
        being the farthest any point of the contour's actual (kx, ky) range
        can be from the origin.

        Used as the tiling window in the zone's own (pre-azimuth-rotation)
        frame: a square of half-width R centred on the origin contains every
        point within distance R of the origin regardless of how that square
        is itself rotated (a disk of radius R fits inside any square of
        half-width R around the same centre) -- so tiling against this
        square, then rotating the result by the azimuth, is guaranteed to
        cover the true (kx, ky) viewport, whatever the azimuth turns out to
        be. It can only over-cover (a few extra tiles get computed and
        clipped away by the plot itself), never under-cover.
        """
        kx = np.asarray(self.contour.angle_defl, dtype=float)
        ky = np.asarray(self.contour.angle_slit, dtype=float)
        half = float(np.hypot(np.max(np.abs(kx)), np.max(np.abs(ky))))
        return (-max(half, 1e-6), max(half, 1e-6)), (-max(half, 1e-6), max(half, 1e-6))

    def lattice_params(self):
        from tools.lattice import LatticeParams
        return LatticeParams(
            self.a_box.value(), self.b_box.value(), self.c_box.value(),
            self.alpha_box.value(), self.beta_box.value(), self.gamma_box.value(),
            space_group=int(self.spacegroup_box.value()))

    def _zone_3d(self):
        """``(zone, b_primitive, b_conventional)`` for the current crystal
        and zone kind, cached.

        Two reciprocal bases, because they answer two different questions
        and confusing them is a real trap for a centred lattice:

        * ``b_primitive`` is the reciprocal of the *primitive* cell. It is
          the lattice the Brillouin zone is the Wigner-Seitz cell of, and
          the one whose translations tile copies of that zone across the
          field of view.
        * ``b_conventional`` is the reciprocal of the *conventional* cell,
          and it is what Miller indices mean. ``(0 0 1)`` on an FCC crystal
          is the cubic face everyone calls (001); read against the primitive
          basis the same three numbers would instead point along a <111>,
          which is not what anyone typing them intends.

        Cached because three different things ask for it while the user
        works -- the offset readout, the 3-D preview, and Plot itself -- and
        building a Wigner-Seitz cell (and then carving a wedge out of it) is
        the one genuinely non-trivial computation here.
        """
        from tools.lattice import (primitive_vectors, conventional_vectors,
                                 reciprocal_vectors, point_group_operations)
        from tools.bz3d import wigner_seitz_cell, irreducible_wedge

        params = self.lattice_params()
        kind = self.zone_kind()
        key = (params.a, params.b, params.c, params.alpha, params.beta,
               params.gamma, params.space_group, kind)
        if self._zone_cache is not None and self._zone_cache[0] == key:
            return self._zone_cache[1], self._zone_cache[2], self._zone_cache[3]

        b = reciprocal_vectors(primitive_vectors(params))
        b_conventional = reciprocal_vectors(conventional_vectors(params))
        zone = wigner_seitz_cell(b)
        if kind == "irreducible":
            zone = irreducible_wedge(zone, point_group_operations(params))
        self._zone_cache = (key, zone, b, b_conventional)
        return zone, b, b_conventional

    def cut_geometry(self):
        """``(zone, b, b_conventional, normal, point, u, v)`` for the
        current 3-D settings: the zone, the primitive reciprocal vectors
        that tile it, the conventional ones the Miller indices are read
        against, the cut plane's unit normal and a point on it, and the cut
        polygon in the plane's own coordinates (``None`` if the plane misses
        the zone).

        Shared by Plot and the 3-D preview so that what the preview shows is
        by construction the polygon that gets drawn on the contour, rather
        than a second calculation that could drift from it.
        """
        from tools.bz3d import plane_cut
        if not self.mode_3d.isChecked():
            raise ValueError("there is no cut plane in 2D-layer mode")
        zone, b, b_conventional = self._zone_3d()
        normal = (self.h_box.value() * b_conventional[0]
                  + self.k_box.value() * b_conventional[1]
                  + self.l_box.value() * b_conventional[2])
        length = float(np.linalg.norm(normal))
        if length < 1e-9:
            raise ValueError("the cut plane normal (h, k, l) cannot be (0, 0, 0)")
        normal = normal / length
        point = normal * float(self.offset_box.value())
        u, v = plane_cut(zone.faces, normal, point)
        return zone, b, b_conventional, normal, point, u, v

    @staticmethod
    def _distance_to_boundary(zone, normal) -> float:
        """How far the plane can slide along ``normal`` before it leaves the
        zone -- the number that makes an offset in A^-1 mean something."""
        from tools.bz3d import face_planes
        normals, offsets = face_planes(zone.faces)
        along = normals @ normal
        ahead = along > 1e-12
        if not np.any(ahead):
            return float("inf")
        return float(np.min(offsets[ahead] / along[ahead]))

    def cut_description(self) -> str:
        """One line about where the cut plane currently sits, for the dialog
        and for the preview window's own caption."""
        try:
            zone, _b, _b_conv, normal, _point, u, _v = self.cut_geometry()
        except Exception as exc:
            return str(exc)
        offset = float(self.offset_box.value())
        reach = self._distance_to_boundary(zone, normal)
        extent = float(np.max(np.linalg.norm(zone.vertices, axis=1)))
        text = (f"({self.h_box.value()} {self.k_box.value()} {self.l_box.value()}) "
                f"normal, offset {offset:.4g} Å⁻¹")
        if not np.isfinite(reach) or reach <= 1e-9 * max(extent, 1.0):
            # A wedge has its apex at Gamma, so along a normal that points
            # straight out of one of its own faces there is no room at all
            # to slide the plane -- worth saying, rather than printing a
            # fraction of something that is numerically zero.
            text += (" — this zone does not extend from Γ along that "
                     "normal at all (an irreducible wedge has its apex "
                     "there), so only offset 0 meets it")
        else:
            text += (f" = {offset / reach:.3f} of the {reach:.4g} "
                     f"Å⁻¹ from Γ to the zone edge along it")
        if u is None:
            text += " — the plane misses the zone entirely."
        else:
            text += f"; the cut is a {len(u) - 1}-sided polygon."
        return text

    def _zone_2d(self, a, b_length, gamma_deg):
        """``(polygons, g1, g2)`` for a plane lattice: the zone (or its
        irreducible wedge), tiled and rotated ready to draw."""
        from tools.bz2d import (reciprocal_vectors_2d, wigner_seitz_cell_2d,
                              tile_2d, lattice_point_group_2d,
                              irreducible_cell_2d)
        a1 = np.array([a, 0.0])
        a2 = np.array([b_length * np.cos(np.radians(gamma_deg)),
                       b_length * np.sin(np.radians(gamma_deg))])
        g1, g2 = reciprocal_vectors_2d(a1, a2)
        polygon = wigner_seitz_cell_2d(g1, g2)
        if self.zone_kind() == "irreducible":
            polygon = irreducible_cell_2d(polygon, lattice_point_group_2d(g1, g2))
        if self.tile_cb.isChecked():
            x_range, y_range = self._view_window()
            tiles = tile_2d(polygon, g1, g2, x_range, y_range)
        else:
            tiles = [polygon]
        azimuth = self.azimuth_box.value()
        return ([self._rotate_xy(t[:, 0], t[:, 1], azimuth) for t in tiles],
                g1, g2)

    # -- picking a direction off the contour ----------------------------------
    def _toggle_direction_picking(self, on: bool):
        if on:
            self._direction_points = []
            self.contour.view.clear_overlays()
            self.contour.view.start_point_picking(self._add_direction_point)
            self.direction_note.setText(
                "Click one or two points on the contour along the direction "
                "you want the angle of.")
        else:
            self.contour.view.stop_point_picking()
            self.contour.view.clear_overlays()

    @staticmethod
    def _angle_of(p0, p1) -> float:
        dx, dy = p1[0] - p0[0], p1[1] - p0[1]
        if abs(dx) < 1e-12 and abs(dy) < 1e-12:
            raise ValueError("the two points coincide")
        return float(np.degrees(np.arctan2(dy, dx)))

    def _add_direction_point(self, x: float, y: float):
        self._direction_points.append((x, y))
        if len(self._direction_points) > 2:
            self._direction_points = self._direction_points[-1:]
        xs = [p[0] for p in self._direction_points]
        ys = [p[1] for p in self._direction_points]
        self.contour.view.show_picked_points(xs, ys)
        if len(self._direction_points) == 1:
            pair = [(0.0, 0.0), (xs[0], ys[0])]
        else:
            pair = list(self._direction_points)
        try:
            angle = self._angle_of(pair[0], pair[1])
        except ValueError as exc:
            self.direction_note.setText(str(exc))
            return
        if len(self._direction_points) == 1:
            self.contour.view.show_overlay_path([0.0, xs[0]], [0.0, ys[0]])
            self.direction_note.setText(
                f"One point: the line from k = (0, 0) to it sits at "
                f"{angle:.3f}° (counter-clockwise from kx). Click a "
                "second point to use the line between the two instead.")
        else:
            self.contour.view.show_overlay_path(xs, ys)
            self.direction_note.setText(
                f"Two points: {angle:.3f}° (counter-clockwise from kx).")

    # -- drawing, on request only ----------------------------------------------
    def plot(self):
        """Work out the zone and draw it. The only thing here that touches
        the contour."""
        try:
            if self.moire_cb.isChecked() and not self.mode_3d.isChecked():
                message = self._plot_moire()
            else:
                message = self._plot_single()
        except Exception as exc:
            # The old drawing is taken away rather than left behind: an
            # overlay that no longer matches the settings beside it is worse
            # than no overlay.
            self.contour.view.clear_bz_overlay()
            self.note.setText(f"Could not plot: {exc}")
            return None
        self.note.setText(message)
        self.contour.statusBar().showMessage(message)
        self._save_settings()
        return message

    def _plot_single(self) -> str:
        if self.mode_3d.isChecked():
            from tools.bz3d import tile_and_cut
            zone, b, _b_conv, normal, point, u, v = self.cut_geometry()
            if self.tile_cb.isChecked():
                u_range, v_range = self._view_window()
                cuts = tile_and_cut(zone, b, normal, point, u_range, v_range)
            else:
                cuts = [] if u is None else [(u, v)]
            if not cuts:
                raise ValueError(
                    "the cut plane does not meet the zone -- reduce the "
                    "offset, or check the (h k l) normal")
            azimuth = self.azimuth_box.value()
            polygons = [self._rotate_xy(cu, cv, azimuth) for cu, cv in cuts]
        else:
            polygons, _g1, _g2 = self._zone_2d(
                self.a2_box.value(), self.b2_box.value(), self.gamma2_box.value())

        for name in ("bz_top", "bz_bottom", "bz_moire"):
            self.contour.view.set_bz_layer(name, [])
        self.contour.view.set_bz_layer("bz_main", polygons, color="#39c2ff")
        return (f"{self.zone_kind()} zone: {len(polygons)} polygon(s) drawn "
                f"on the contour.")

    def _plot_moire(self) -> str:
        from tools.bz2d import (reciprocal_vectors_2d, wigner_seitz_cell_2d,
                              tile_2d, lattice_point_group_2d,
                              irreducible_cell_2d)
        from tools.moire import hex_moire_lattice_fast, moire_bz

        a_top, b_top, gamma_top = (self.a2_box.value(), self.b2_box.value(),
                                   self.gamma2_box.value())
        if self.moire_same_cb.isChecked():
            a_bot, b_bot, gamma_bot = a_top, b_top, gamma_top
        else:
            a_bot, b_bot, gamma_bot = (self.a2b_box.value(), self.b2b_box.value(),
                                       self.gamma2b_box.value())
        twist = self.twist_box.value()
        azimuth = self.azimuth_box.value()
        irreducible = self.zone_kind() == "irreducible"

        top_polygons, g1_top, g2_top = self._zone_2d(a_top, b_top, gamma_top)

        # The bottom layer is built rotated by the twist: that is the
        # *relative* angle between the layers, while the dialog's azimuth
        # turns the finished picture as a whole, so the two must not be
        # rolled into one rotation.
        a1_bot = np.array([a_bot, 0.0])
        a2_bot = np.array([b_bot * np.cos(np.radians(gamma_bot)),
                           b_bot * np.sin(np.radians(gamma_bot))])
        radians = np.radians(twist)
        rotation = np.array([[np.cos(radians), -np.sin(radians)],
                             [np.sin(radians), np.cos(radians)]])
        a1_bot, a2_bot = rotation @ a1_bot, rotation @ a2_bot
        g1_bot, g2_bot = reciprocal_vectors_2d(a1_bot, a2_bot)

        def drawn(polygon, g1, g2):
            if irreducible:
                polygon = irreducible_cell_2d(polygon, lattice_point_group_2d(g1, g2))
            if self.tile_cb.isChecked():
                x_range, y_range = self._view_window()
                tiles = tile_2d(polygon, g1, g2, x_range, y_range)
            else:
                tiles = [polygon]
            return [self._rotate_xy(t[:, 0], t[:, 1], azimuth) for t in tiles]

        bottom_polygons = drawn(wigner_seitz_cell_2d(g1_bot, g2_bot), g1_bot, g2_bot)
        gm1, gm2, moire_polygon = moire_bz(g1_top, g2_top, g1_bot, g2_bot)
        moire_polygons = drawn(moire_polygon, gm1, gm2)

        am1, _am2 = reciprocal_vectors_2d(gm1, gm2)
        moire_period = float(np.linalg.norm(am1))

        self.contour.view.set_bz_layer(
            "bz_top", top_polygons if self.show_top_cb.isChecked() else [],
            color="#39c2ff")
        self.contour.view.set_bz_layer(
            "bz_bottom", bottom_polygons if self.show_bottom_cb.isChecked() else [],
            color="#ff9c39")
        self.contour.view.set_bz_layer(
            "bz_moire", moire_polygons if self.show_moire_cb.isChecked() else [],
            color="#ff3ce0")
        self.contour.view.set_bz_layer("bz_main", [])

        message = (f"moire real-space period ≈ {moire_period:.4g} Å "
                   f"({len(top_polygons)} top, {len(bottom_polygons)} bottom, "
                   f"{len(moire_polygons)} moire polygon(s))")
        hexagonal = (lambda length, other, angle:
                     abs(length - other) < 1e-6 * max(length, 1e-9)
                     and abs(angle - 120.0) < 0.5)
        if (hexagonal(a_top, b_top, gamma_top)
                and hexagonal(a_bot, b_bot, gamma_bot)):
            try:
                fast, _rotation = hex_moire_lattice_fast(a_top, a_bot, twist)
                message += f"; closed-form (hexagonal) fast path: {fast:.4g} Å"
                if abs(twist) > 30.0 and np.isclose(a_top, a_bot):
                    message += (" (only exact for |twist| <= 30 deg; the general "
                                "method used for the drawing folds correctly "
                                "beyond that, see tools/moire.py)")
            except Exception:
                pass
        return message

    def _clear_overlay(self):
        self.contour.view.clear_bz_overlay()
        self.note.setText("Overlay cleared.")

    # -- persistence: data.scan.info["bz.*"], flat scalars like kconv/arbcut -
    def _settings_widgets(self):
        return (self.mode_3d, self.mode_2d, self.a_box, self.b_box, self.c_box,
                self.alpha_box, self.beta_box, self.gamma_box,
                self.spacegroup_box, self.h_box, self.k_box, self.l_box,
                self.offset_box, self.a2_box, self.b2_box, self.gamma2_box,
                self.azimuth_box, self.tile_cb, self.moire_cb,
                self.moire_same_cb, self.a2b_box, self.b2b_box,
                self.gamma2b_box, self.twist_box, self.show_top_cb,
                self.show_bottom_cb, self.show_moire_cb,
                self.zone_conventional, self.zone_irreducible)

    def _load_saved_settings(self):
        widgets = self._settings_widgets()
        for widget in widgets:
            widget.blockSignals(True)
        try:
            info = self.contour.data.scan.info

            def saved(key, default):
                return info.get(self._PREFIX + key, default)

            # Not setChecked(False) on whichever radio is currently checked:
            # for an exclusive pair, Qt refuses to leave both unchecked, so
            # unchecking the sole checked one is a silent no-op -- checking
            # the *other* one is what actually switches the pair.
            if saved("mode", "3d") == "2d":
                self.mode_2d.setChecked(True)
            else:
                self.mode_3d.setChecked(True)
            if saved("zone", "conventional") == "irreducible":
                self.zone_irreducible.setChecked(True)
            else:
                self.zone_conventional.setChecked(True)
            self.spacegroup_box.setValue(int(saved("space_group", 1)))
            self.a_box.setValue(float(saved("a", 3.0)))
            self.b_box.setValue(float(saved("b", 3.0)))
            self.c_box.setValue(float(saved("c", 3.0)))
            self.alpha_box.setValue(float(saved("alpha", 90.0)))
            self.beta_box.setValue(float(saved("beta", 90.0)))
            self.gamma_box.setValue(float(saved("gamma", 90.0)))
            self.h_box.setValue(int(saved("h", 0)))
            self.k_box.setValue(int(saved("k", 0)))
            self.l_box.setValue(int(saved("l", 1)))
            self.offset_box.setValue(float(saved("offset", 0.0)))
            self.a2_box.setValue(float(saved("a2", 3.0)))
            self.b2_box.setValue(float(saved("b2", 3.0)))
            self.gamma2_box.setValue(float(saved("gamma2", 120.0)))
            self.azimuth_box.setValue(float(saved("azimuth_deg", 0.0)))
            self.tile_cb.setChecked(bool(int(saved("tile", 1))))
            self.moire_cb.setChecked(bool(int(saved("moire", 0))))
            self.moire_same_cb.setChecked(bool(int(saved("moire_same", 1))))
            self.a2b_box.setValue(float(saved("a2b", 3.0)))
            self.b2b_box.setValue(float(saved("b2b", 3.0)))
            self.gamma2b_box.setValue(float(saved("gamma2b", 120.0)))
            self.twist_box.setValue(float(saved("twist_deg", 0.0)))
            self.show_top_cb.setChecked(bool(int(saved("show_top", 1))))
            self.show_bottom_cb.setChecked(bool(int(saved("show_bottom", 1))))
            self.show_moire_cb.setChecked(bool(int(saved("show_moire", 1))))
        except Exception:
            pass    # a first-ever use of the tool has nothing to load yet
        finally:
            for widget in widgets:
                widget.blockSignals(False)

    def _save_settings(self):
        info = self.contour.data.scan.info
        values = {
            "mode": "3d" if self.mode_3d.isChecked() else "2d",
            "zone": self.zone_kind(),
            "space_group": int(self.spacegroup_box.value()),
            "a": self.a_box.value(), "b": self.b_box.value(), "c": self.c_box.value(),
            "alpha": self.alpha_box.value(), "beta": self.beta_box.value(),
            "gamma": self.gamma_box.value(),
            "h": int(self.h_box.value()), "k": int(self.k_box.value()),
            "l": int(self.l_box.value()),
            "offset": self.offset_box.value(),
            "a2": self.a2_box.value(), "b2": self.b2_box.value(),
            "gamma2": self.gamma2_box.value(),
            "azimuth_deg": self.azimuth_box.value(),
            "tile": int(self.tile_cb.isChecked()),
            "moire": int(self.moire_cb.isChecked()),
            "moire_same": int(self.moire_same_cb.isChecked()),
            "a2b": self.a2b_box.value(), "b2b": self.b2b_box.value(),
            "gamma2b": self.gamma2b_box.value(),
            "twist_deg": self.twist_box.value(),
            "show_top": int(self.show_top_cb.isChecked()),
            "show_bottom": int(self.show_bottom_cb.isChecked()),
            "show_moire": int(self.show_moire_cb.isChecked()),
        }
        for key, value in values.items():
            info[self._PREFIX + key] = value

    def closeEvent(self, event):
        # Same rule as the other picking dialogs: picking must not outlive
        # the dialog. The drawn zone itself is untouched -- it is meant to
        # stay on the contour after this closes.
        if self.direction_button.isChecked():
            self.direction_button.setChecked(False)
        else:
            self.contour.view.stop_point_picking()
            self.contour.view.clear_overlays()
        if self._preview is not None:
            self._preview.close()
        super().closeEvent(event)


class ContourWindow(ViewerWindow):
    """The constant-energy contour of a deflector map -- what a map file
    opens on by default, since that is the view you navigate the cube from.

    The two buttons open the orthogonal cuts, each in its own window, so a
    contour and one or both cuts can sit side by side rather than competing
    for space in one layout. Re-clicking a button raises the window that is
    already open instead of making a duplicate.

    "Map k conversion" turns the measurement into a momentum-space map (see
    ``tools.kspace``); the result is emitted on :attr:`kMapCreated` for the
    launcher to list, rather than replacing what this window shows.
    """

    #: a converted k-space map, for the launcher to add to its list
    kMapCreated = pyqtSignal(object)

    def __init__(self, data, filename, colormap, flip):
        super().__init__(data, filename, colormap, flip)
        self.angle_defl, self.angle_slit, self.E, self.cube = data.angle_cube
        labels = data.scan.labels
        self.defl_label = labels.get("x", "angle (deflector)")
        self.slit_label = labels.get("k", "angle (along slit)")
        self.energy_label = labels.get("z", "Energy (eV)")
        self.cut_windows = {}

        self.defl_cut_button = QPushButton("Deflector cut")
        self.defl_cut_button.setToolTip(
            "Open the deflector-vs-energy cut, integrated over the slit-angle window.")
        self.slit_cut_button = QPushButton("Slit cut")
        self.slit_cut_button.setToolTip(
            "Open the slit-vs-energy cut, integrated over the deflector-angle window.")
        self.kconv_button = QPushButton("Map k conversion")
        self.arbcut_button = QPushButton("Arbitrary cut")
        self.arbcut_button.setToolTip(
            "Cut along a path of up to six points picked on this contour, "
            "instead of along one of the two axes.")
        self.defl_cut_button.setMinimumWidth(100)
        self.slit_cut_button.setMinimumWidth(90)
        self.kconv_button.setToolTip(
            "Convert this map from angle space to momentum space (A^-1). "
            "Pick the point on the contour that should become k = (0, 0).")
        self.defl_cut_button.clicked.connect(lambda: self.open_cut("deflector"))
        self.slit_cut_button.clicked.connect(lambda: self.open_cut("slit"))
        self.kconv_button.clicked.connect(self.open_k_conversion)
        self.arbcut_button.clicked.connect(self.open_arbitrary_cut)
        # Already in momentum space: converting again would be meaningless.
        # A kz map is hidden too, for a different reason -- its first axis is
        # a photon energy, and the in-plane formula does not apply to it.
        # (The refusal in open_k_conversion() is still there as the backstop,
        # and explains why; this just stops the button being offered.)
        self.kconv_button.setVisible(data.kind not in ("k_map", "kz_map"))
        self.bz_button = QPushButton("Brillouin zone...")
        self.bz_button.setToolTip(
            "Overlay a Brillouin zone (or a moire zone) on this contour. "
            "Needs momentum, not angle -- only available once a map has "
            "been converted to k-space.")
        self.bz_button.clicked.connect(self.open_brillouin_zone)
        # A Brillouin zone is a statement about momentum-space periodicity,
        # so it only means something once the axes are actually k, not angle
        # -- the mirror image of kconv_button's own gating above.
        self.bz_button.setVisible(data.kind == "k_map")
        self.kz_button = QPushButton("kz map processing...")
        self.kz_button.setToolTip(
            "Fit the Fermi edge of every spectrum in this photon-energy "
            "scan, shift each one to put its own edge at zero, and crop to "
            "the energy range they all still cover.\n\n"
            "The scan arrives stacked as measured, because nothing in the "
            "files says where each spectrum's Fermi level is. This measures "
            "it from the spectra themselves.")
        self.kz_button.clicked.connect(self.open_kz_processing)
        # Only a photon-energy scan needs this: it is the one kind whose
        # members were each measured against a different reference.
        self.kz_button.setVisible(data.kind == "kz_map")
        self.slices_button = QPushButton("Slice figure...")
        self.slices_button.setToolTip(
            "A page of slices through this cube as one figure for a paper -- "
            "a row of constant-energy contours, or a row of cuts -- with a "
            "shared colour scale and labels only on the outer edges.")
        self.slices_button.clicked.connect(self.open_slice_figure)
        self.build_toolbar((self.defl_cut_button, self.slit_cut_button,
                            self.arbcut_button, self.kconv_button,
                            self.bz_button, self.kz_button,
                            self.slices_button))
        self.figure_windows = []

        self.e_control = _SliceControl("Energy (eV)", "eV", 0.05, self)
        self.root.addWidget(self.e_control)

        self.view = FrameImageView()
        self.panel = ImagePanel(
            f"Constant-E contour: {self.defl_label} vs {self.slit_label}",
            self.view, equal_ratio_option=True, lock_ratio=True, parent=self)
        self.root.addWidget(self.panel, stretch=1)

        self.view.set_axis_labels(self.defl_label, self.slit_label)
        self.view.enable_selection_action(False)
        self.e_control.configure(len(self.E), self.E)
        self.e_control.changed.connect(self.refresh_contour)
        self.refresh_contour()
        self.set_colormap(colormap, flip)
        self.apply_display_options()
        # Open at the contour's own aspect so the 1:1 lock has nothing to
        # stretch; without this a tall-narrow map leaves wide empty margins.
        self.size_to_data_aspect(self.view)
        self.view.fit_to_data()

    def image_panels(self):
        return (self.panel,)

    @staticmethod
    def index_window(axis, centre_idx, half_width):
        centre_idx = int(np.clip(centre_idx, 0, len(axis) - 1))
        if half_width <= 0:
            return np.array([centre_idx])
        centre = axis[centre_idx]
        idxs = np.nonzero((axis >= centre - half_width) & (axis <= centre + half_width))[0]
        return idxs if idxs.size else np.array([centre_idx])

    def full_cube(self):
        """The whole cube as a real numpy array.

        ``self.cube`` is usually a :class:`loader.nxs_file.LazyCube` reading
        from the file as it is sliced, which is what makes opening a map
        instant and lets a browsed-but-unused map cost nothing. An algorithm
        that needs every point -- a k conversion, an arbitrary cut, a page
        of slices -- has to have it read, and says so here rather than
        relying on a lazy array quacking like an ndarray in every last
        respect (it forwards indexing and ``np.asarray``, not ndarray's own
        methods).
        """
        return np.asarray(self.cube)

    def refresh_contour(self):
        ei = self.e_control.index
        idxs = self.index_window(self.E, ei, self.e_control.half_width)
        self.e_control.show_value(self.E[ei], idxs, self.E)
        self.view.set_frame(self.e_control.combine(self.cube[:, :, idxs], 2),
                            self.angle_defl, self.angle_slit)
        self.panel.sync_level_range()
        self.set_colormap(self.colormap, self.flip)

    def open_kz_processing(self):
        """Calibrate a photon-energy scan against its own Fermi edges.

        Modeless, and for a reason that is not the usual one: the box it
        works from is dragged on the *slit cut* window, which this opens
        alongside. A modal dialog would make that impossible.
        """
        from ui.kzmap import KzMapProcessDialog

        existing = getattr(self, "_kz_dialog", None)
        if existing is not None:
            existing.raise_()
            existing.activateWindow()
            return existing
        dialog = KzMapProcessDialog(self)
        dialog.datasetsCreated.connect(
            lambda made: [self.datasetCreated.emit(one) for one in made])
        self._kz_dialog = dialog
        dialog.finished.connect(lambda *_: setattr(self, "_kz_dialog", None))
        dialog.show()
        return dialog

    def open_cut(self, which: str):
        window = self.cut_windows.get(which)
        if window is not None:
            window.raise_()
            window.activateWindow()
            return
        window = MapCutWindow(self, which)
        self.adopt_child(window)
        window.closed.connect(lambda w, key=which: self.cut_windows.pop(key, None))
        self.cut_windows[which] = window
        window.set_colormap(self.colormap, self.flip)
        window.show()

    def exportable_panels(self):
        if self.view.last_frame is None:
            return []
        return [("ConstE_contour", self.view.last_frame, self.angle_defl, self.angle_slit,
                 self.defl_label, self.slit_label)]

    def slice_label(self) -> str:
        idx = self.e_control.index
        return f"{self.E[idx]:.4g} eV (ind {idx + 1})"

    # -- arbitrary-direction cut -------------------------------------------
    # -- a page of slices, as one figure ----------------------------------
    def slice_axes(self):
        """``(axes, labels)`` of the cube in the order it is stored, which
        is what the slicing dialog and :func:`slice_panels` work in."""
        labels = (self.defl_label, self.slit_label, self.energy_label)
        return (self.angle_defl, self.angle_slit, self.E), labels

    def open_slice_figure(self):
        """The MATLAB tool's ``massplott_3D``: cut the cube along one axis
        at a series of positions and lay the results out as a figure."""
        from ui.figure import SliceSeriesDialog, FigureWindow, slice_panels

        axes, labels = self.slice_axes()
        dialog = SliceSeriesDialog(self, axes, labels, kind=self.data.kind)
        if dialog.exec_() != QDialog.Accepted:
            return None
        settings = dialog.settings()
        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            panels = slice_panels(
                self.full_cube(), axes, labels, axis=settings["axis"],
                values=settings["values"], half_width=settings["half_width"],
                titles=settings["titles"], colormap=self.colormap,
                flip=self.flip)
        except Exception as exc:
            QMessageBox.warning(self, "Slice figure", f"Could not slice it:\n{exc}")
            return None
        finally:
            QApplication.restoreOverrideCursor()
        if not panels:
            return None

        window = FigureWindow(panels, f"Figure — {self.filename}", parent=None,
                              dataset_source=self._figure_source())
        window.figure.fit_grid(settings["cols"])
        window.figure.letter_panels()
        window.figure.style.shared_levels = settings["shared"]
        # Both axes of a constant-energy slice are angles or momenta, so an
        # equal scale is the honest default there and nowhere else.
        window.figure.style.equal_aspect = settings["axis"] == 2
        window.refresh()
        window.closed.connect(lambda w: self.figure_windows.remove(w)
                              if w in self.figure_windows else None)
        self.figure_windows.append(window)
        window.show()
        self.statusBar().showMessage(
            f"{len(panels)} slices laid out as a figure")
        return window

    def _figure_source(self):
        entries = getattr(self, "_dataset_entries", None)
        loader = getattr(self, "_dataset_loader", None)
        if callable(entries) and callable(loader):
            return {"entries": entries, "loader": loader}
        return None

    def open_arbitrary_cut(self):
        """Collect the path, then cut along it. Modeless, like the other
        picking dialogs: the points are picked on the contour behind it."""
        existing = getattr(self, "_arbcut_dialog", None)
        if existing is not None and existing.isVisible():
            existing.raise_()
            existing.activateWindow()
            return None

        dialog = ArbitraryCutDialog(self)
        self._arbcut_dialog = dialog
        self.statusBar().showMessage(
            "Click the contour to place the path's points, then press Plot cut.")

        def on_accepted():
            points, name = dialog.points(), dialog.name_box.text().strip()
            separate = dialog.separate_cb.isChecked()
            self._arbcut_dialog = None
            self.run_arbitrary_cut(points, name=name, separate=separate)

        dialog.accepted.connect(on_accepted)
        dialog.rejected.connect(lambda: setattr(self, "_arbcut_dialog", None))
        dialog.show()
        dialog.raise_()
        return dialog

    def run_arbitrary_cut(self, points, name: str = None, separate: bool = False,
                          open_window: bool = True):
        """Sample the cube along ``points`` and list the result.

        Split from the dialog so it can be driven directly (which is how the
        tests exercise it). Returns the datasets it created.
        """
        name = name or f"{self.filename} [arb cut]"
        unit = _unit_of(self.defl_label)
        distance_label = f"Distance along path ({unit})" if unit else "Distance along path"
        try:
            distance, energy, values, joints = arbitrary_cut(
                self.angle_defl, self.angle_slit, self.E, self.full_cube(), points)
        except Exception as exc:
            QMessageBox.warning(self, "Arbitrary cut", f"The cut failed:\n{exc}")
            return []

        source_info = {k: v for k, v in self.data.scan.info.items()
                       if not k.startswith("arbcut.")}
        created = []
        if separate:
            # One dataset per segment, each starting from distance 0 -- the
            # MATLAB dialog's "separate" option.
            edges = list(joints) + [float(distance[-1])]
            start = 0.0
            for i, end in enumerate(edges):
                mask = (distance >= start - 1e-9) & (distance <= end + 1e-9)
                seg_points = [points[i], points[i + 1]]
                created.append(self._make_cut_dataset(
                    f"{name} {i + 1}", distance[mask] - distance[mask][0], energy,
                    values[mask], distance_label,
                    {"points": seg_points, "segment": i + 1}, source_info))
                start = end
        else:
            created.append(self._make_cut_dataset(
                name, distance, energy, values, distance_label,
                {"points": list(points), "joints": list(map(float, joints))},
                source_info))

        for data in created:
            self.datasetCreated.emit(data)
        self.statusBar().showMessage(
            f"Arbitrary cut: {len(points)} points, {values.shape[0]} samples "
            f"-- added {len(created)} dataset(s) to the file list")
        if open_window and created:
            self.open_computed(created[0])
        return created

    def _make_cut_dataset(self, name, distance, energy, values, distance_label,
                          parameters, source_info):
        return MemoryData(
            "cut", (np.asarray(distance, dtype=float), np.asarray(energy, dtype=float)),
            np.asarray(values, dtype=float),
            {"x": distance_label, "y": self.energy_label},
            source_label=name, parameters=parameters, prefix="arbcut",
            source_path=getattr(self.data, "path", ""),
            source_info=source_info,
            source_motors=dict(self.data.scan.fourd_info))

    # -- angle -> momentum ------------------------------------------------
    def readout_position(self):
        """The readout cursor's (deflector angle, slit angle), or None if it
        is not switched on. This is what the conversion dialog reads to place
        the k-space origin."""
        cursor = getattr(self.view, "readout_cursor", None)
        if cursor is None:
            return None
        x, y = cursor.pos()
        return float(x), float(y)

    def open_k_conversion(self):
        """Ask for the conversion settings, then convert this map to k-space.

        The dialog is deliberately **modeless**: its whole premise is that
        the user drags the readout cursor on the contour behind it to choose
        the origin, and a modal dialog would swallow those clicks -- the
        contour could only be touched by closing the dialog first. So it is
        shown rather than exec'd, and the conversion runs from its
        ``accepted`` signal.

        The readout cursor is switched on first if it is not already, since
        an empty pair of offset boxes would be a poor way to ask for a point.
        """
        existing = getattr(self, "_kconv_dialog", None)
        if existing is not None and existing.isVisible():
            existing.raise_()
            existing.activateWindow()
            return None

        if getattr(self.view, "readout_cursor", None) is None:
            self.view.set_readout_cursor_visible(True)
        self.statusBar().showMessage(
            "Drag the readout cursor to the point that should become "
            "k = (0, 0), then press Convert.")

        dialog = KConversionDialog(self)
        dialog.setModal(False)
        self._kconv_dialog = dialog
        dialog.follow_cursor()

        # Keep the boxes on the cursor while the user drags it around, until
        # they type a value of their own.
        follow = lambda *_: dialog.follow_cursor()
        self.view.readoutMoved.connect(follow)

        def cleanup(clear_cursor: bool = False):
            try:
                self.view.readoutMoved.disconnect(follow)
            except TypeError:
                pass
            self._kconv_dialog = None
            # The picking overlays and the cursor exist to set up *this*
            # conversion; once it has run they are leftovers on the contour,
            # so the picture goes back to what it was.
            self.view.stop_point_picking()
            self.view.clear_overlays()
            if clear_cursor:
                self.view.set_readout_cursor_visible(False)
                action = getattr(self.view, "action_readout", None)
                if action is not None and action.isChecked():
                    action.blockSignals(True)
                    action.setChecked(False)
                    action.blockSignals(False)

        def on_accepted():
            settings = dialog.settings()
            name = dialog.output_name()
            cleanup(clear_cursor=True)
            self.run_k_conversion(settings, name=name)

        def on_rejected():
            cleanup()
            self.statusBar().showMessage("k conversion cancelled")

        dialog.accepted.connect(on_accepted)
        dialog.rejected.connect(on_rejected)
        dialog.show()
        dialog.raise_()
        return dialog

    def first_axis_is_an_angle(self) -> bool:
        """Whether this map's first axis is an emission angle.

        The same acquisition records a photon-energy, temperature or
        gate-voltage series, and the loader asks which it is. Converting any
        of the others to momentum would produce a picture that looks
        entirely plausible and means nothing, so the k conversion checks.
        """
        from loader import registry
        return registry.role_is_angle(self.data.scan)

    def run_k_conversion(self, settings: dict, name: str = None,
                         background: bool = True):
        """Do the conversion and hand the result to whoever is listening.

        Split from :meth:`open_k_conversion` so the conversion can be driven
        without a dialog (which is how the tests exercise it, with
        ``background=False`` so the answer is back when the call returns).

        The conversion itself is the longest thing this program does on a
        map -- it interpolates the whole cube onto a new grid -- so it goes
        to a worker thread. The cube is read *before* the thread starts: the
        GUI keeps slicing the same file to repaint, and HDF5 here is not
        built for two threads at once.
        """
        from ui import jobs
        from tools.kspace import convert_map

        if not self.first_axis_is_an_angle():
            role = self.data.scan.info.get("axis0.role", "not an angle")
            QMessageBox.information(
                self, "k conversion",
                f"This map's first axis was loaded as “{role}”, not an "
                f"emission angle, so converting it to momentum would not "
                f"mean anything.\n\nIf that is wrong, re-open the file from "
                f"the loader and set the first axis to an angle.")
            return None

        self.statusBar().showMessage("Reading the cube...")
        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            cube = self.full_cube()
        finally:
            QApplication.restoreOverrideCursor()

        defl, slit, energy_axis = self.angle_defl, self.angle_slit, self.E

        def work(report):
            report(None, "Converting to k-space...")
            return convert_map(defl, slit, energy_axis, cube, **settings)

        if background:
            started = jobs.run_job(
                self, "k conversion", work,
                on_done=lambda result: self._finish_k_conversion(result, settings, name),
                on_error=lambda exc: self._failed_k_conversion(exc),
                on_cancel=lambda: self.statusBar().showMessage("k conversion cancelled"))
            if started:
                self.statusBar().showMessage("Converting to k-space...")
                return None

        try:
            result = work(lambda *a, **k: None)
        except Exception as exc:
            self._failed_k_conversion(exc)
            return None
        return self._finish_k_conversion(result, settings, name)

    def _failed_k_conversion(self, exc):
        QMessageBox.warning(self, "k conversion", f"Conversion failed:\n{exc}")
        self.statusBar().showMessage("k conversion failed")

    def _finish_k_conversion(self, result, settings: dict, name: str = None):
        """Wrap the converted cube up as a dataset and announce it. Runs on
        the GUI thread, whether the conversion ran here or on a worker."""
        kx, ky, energy, cube = result

        source_info = {k: v for k, v in self.data.scan.info.items()
                       if not k.startswith("kconv.")}
        kdata = KMapData(kx, ky, energy, cube,
                         source_label=name or f"{self.filename} [k]",
                         parameters=settings,
                         source_path=getattr(self.data, "path", ""),
                         source_info=source_info,
                         source_motors=dict(self.data.scan.fourd_info))
        self.kMapCreated.emit(kdata)
        finite = float(np.isfinite(cube).mean() * 100.0)
        self.statusBar().showMessage(
            f"Converted to k-space: {cube.shape[0]}x{cube.shape[1]}x{cube.shape[2]}, "
            f"{finite:.0f}% inside the light cone -- added to the file list")
        return kdata

    # -- Brillouin-zone overlay ---------------------------------------------
    def open_brillouin_zone(self):
        """Open the Brillouin-zone dialog, or raise it if already open.

        Modeless and left running, like ``open_arbitrary_cut``'s dialog:
        there is no result to collect on close, since every change already
        redraws the contour directly.
        """
        existing = getattr(self, "_bz_dialog", None)
        if existing is not None and existing.isVisible():
            existing.raise_()
            existing.activateWindow()
            return existing

        dialog = BrillouinZoneDialog(self)
        self._bz_dialog = dialog
        dialog.finished.connect(lambda *_: setattr(self, "_bz_dialog", None))
        dialog.show()
        dialog.raise_()
        return dialog


class MapCutWindow(ViewerWindow):
    """One orthogonal cut through a deflector map, opened from the contour
    window's buttons. ``which`` is "deflector" (angle-deflector vs E,
    integrated over slit angle) or "slit" (angle-slit vs E, integrated over
    deflector angle)."""

    def __init__(self, contour: ContourWindow, which: str):
        # Take a reference of its own: the cut window and the contour share
        # one dataset and either may be closed first. Only file-backed data
        # is reference counted -- a converted k-map holds no file, so there
        # is nothing to keep alive.
        if hasattr(contour.data, "_refs"):
            contour.data._refs += 1
        super().__init__(contour.data, contour.filename, contour.colormap, contour.flip)
        self.contour = contour
        self.which = which

        if which == "deflector":
            self.x_axis, self.x_label = contour.angle_defl, contour.defl_label
            self.sum_axis, sum_label = contour.angle_slit, contour.slit_label
            title = f"{contour.defl_label} vs {contour.energy_label}"
        else:
            self.x_axis, self.x_label = contour.angle_slit, contour.slit_label
            self.sum_axis, sum_label = contour.angle_defl, contour.defl_label
            title = f"{contour.slit_label} vs {contour.energy_label}"
        self.setWindowTitle(f"{contour.filename}  [{title}]")
        self.resize(900, 760)

        # The correction is fitted along the analyser's slit, because that is
        # the direction its curvature lives in; a deflector-vs-E cut is a
        # different axis entirely, so the button is only offered here.
        self.fs_button = QPushButton("FS correction")
        self.fs_button.setToolTip(
            "Straighten a curved feature along the slit and apply the same "
            "shift to the whole cube, giving a corrected map in the file list.")
        self.fs_button.clicked.connect(self.open_fs_correction)
        self.build_toolbar((self.fs_button,) if which == "slit" else ())
        self.control = _SliceControl(f"Integrate over {sum_label}", "deg", 0.1, self)
        self.root.addWidget(self.control)

        self.view = FrameImageView()
        # No 1:1 default, even for a cut of a converted k-map: this panel is
        # momentum against *energy*, and an equal scale between A^-1 and eV
        # is not a shape, it is a coincidence of the numbers. (The contour
        # it was cut from is kx against ky and does open at 1:1, which is
        # where that default belongs.)
        self.panel = ImagePanel(title, self.view, curves=True, parent=self)
        self.root.addWidget(self.panel, stretch=1)

        self.view.set_axis_labels(self.x_label, contour.energy_label)
        self.view.enable_selection_action(False)
        self.control.configure(len(self.sum_axis), self.sum_axis)
        self.control.changed.connect(self.refresh_cut)
        self.refresh_cut()
        self.apply_display_options()

    def image_panels(self):
        return (self.panel,)

    def refresh_cut(self):
        idx = self.control.index
        idxs = ContourWindow.index_window(self.sum_axis, idx, self.control.half_width)
        self.control.show_value(self.sum_axis[idx], idxs, self.sum_axis)
        cube = self.contour.cube
        # cube is (deflector, slit, E): collapse whichever angle axis this
        # cut integrates over.
        frame = (self.control.combine(cube[:, idxs, :], 1) if self.which == "deflector"
                 else self.control.combine(cube[idxs, :, :], 0))
        self.view.set_frame(frame, self.x_axis, self.contour.E)
        self.panel.sync_curve_source()
        self.panel.sync_level_range()
        self.set_colormap(self.colormap, self.flip)

    def slice_label(self) -> str:
        idx = self.control.index
        return f"{self.sum_axis[idx]:.4g} (ind {idx + 1})"

    # -- Fermi-surface correction -----------------------------------------
    def fs_angle_axis(self):
        return self.x_axis

    def fs_angle_label(self) -> str:
        return self.x_label

    def fs_correction_target(self):
        """The whole cube, corrected along the slit axis.

        The fit is made on this window's slit-vs-E image (integrated over the
        deflector), and the resulting shift is applied to every deflector
        position -- which is the point: one fit straightens the entire map,
        exactly as ``correction.m`` corrects the whole 3D array from a single
        slice.
        """
        if self.which != "slit":
            return None
        contour = self.contour
        scan = contour.data.scan
        return (np.asarray(contour.cube, dtype=float), contour.angle_slit, contour.E,
                1, 2, contour.data.kind,
                lambda new_energy: (np.asarray(contour.angle_defl, dtype=float),
                                    np.asarray(contour.angle_slit, dtype=float),
                                    new_energy))

    def exportable_panels(self):
        if self.view.last_frame is None:
            return []
        suffix = "Deflector_vs_E" if self.which == "deflector" else "Slit_vs_E"
        return [(suffix, self.view.last_frame, self.x_axis, self.contour.E,
                 self.x_label, self.contour.energy_label)]


# --------------------------------------------------------------------------
# Picking points off a picture: the arbitrary cut and the FS correction
# --------------------------------------------------------------------------
def _unit_of(label: str) -> str:
    """The unit in an axis label, e.g. "Angle, deflector (deg)" -> "deg"."""
    if label and "(" in label and label.rstrip().endswith(")"):
        return label[label.rfind("(") + 1:-1]
    return ""


class _PickerDialog(QDialog):
    """Shared plumbing for the two dialogs that collect points off an image.

    Both are **modeless** for the same reason the k-conversion dialog is:
    their whole premise is that the user works on the picture behind them.
    While one is open, left-clicking the image adds a point, the picked
    points are drawn on it, and closing the dialog (however it closes)
    takes the overlays away again.
    """

    def __init__(self, window, view, title: str):
        super().__init__(window)
        self.window_ = window
        self.view = view
        self.setModal(False)
        self.setWindowTitle(title)
        self.finished.connect(lambda *_: self._release())

    def _begin_picking(self):
        self.view.start_point_picking(self.add_point)

    def _release(self):
        try:
            self.view.stop_point_picking()
            self.view.clear_overlays()
        except RuntimeError:        # the window went first
            pass

    def add_point(self, x: float, y: float):
        raise NotImplementedError


class ArbitraryCutDialog(_PickerDialog):
    """Up to six points defining a path across the constant-energy contour;
    the cut follows that path, segment by segment.

    Same shape as ``arbi_cut_plot_demo.m``: six numbered point rows, each
    with its own enable box, a "plot" that joins the segments end to end,
    and the choice of keeping the result as one dataset or one per segment.
    What is different is where the points come from -- MATLAB's ``getline``
    freezes the figure until you press Enter, while here the contour stays
    live: click it to fill the next free row, then drag the numbers if the
    click was off.
    """

    MAX_POINTS = 6

    def __init__(self, contour):
        super().__init__(contour, contour.view, "Arbitrary cut")
        self.contour = contour

        layout = QVBoxLayout(self)
        intro = QLabel(
            "Click the contour to place up to six points; the cut runs "
            "through them in order. Values can be typed or corrected here.")
        intro.setWordWrap(True)
        layout.addWidget(intro)

        unit = _unit_of(contour.defl_label) or "axis units"
        grid = QFormLayout()
        self.rows = []
        for i in range(self.MAX_POINTS):
            row = QWidget()
            box = QHBoxLayout(row)
            box.setContentsMargins(0, 0, 0, 0)
            enable = QCheckBox()
            enable.setToolTip("Include this point in the path.")
            x_box, y_box = QDoubleSpinBox(), QDoubleSpinBox()
            for spin, axis in ((x_box, contour.angle_defl), (y_box, contour.angle_slit)):
                spin.setDecimals(4)
                spin.setRange(float(np.min(axis)), float(np.max(axis)))
                spin.setKeyboardTracking(False)
                spin.setButtonSymbols(QDoubleSpinBox.NoButtons)
                spin.setMaximumWidth(92)
                spin.valueChanged.connect(self._refresh_overlay)
            enable.toggled.connect(self._refresh_overlay)
            box.addWidget(enable)
            box.addWidget(QLabel("x"))
            box.addWidget(x_box)
            box.addWidget(QLabel("y"))
            box.addWidget(y_box)
            box.addStretch(1)
            grid.addRow(f"P{i + 1}", row)
            self.rows.append((enable, x_box, y_box))
        layout.addLayout(grid)
        layout.addWidget(QLabel(f"Coordinates are in {unit}."))

        buttons_row = QHBoxLayout()
        self.clear_button = QPushButton("Clear points")
        self.clear_button.clicked.connect(self.clear_points)
        buttons_row.addWidget(self.clear_button)
        buttons_row.addStretch(1)
        layout.addLayout(buttons_row)

        form = QFormLayout()
        self.name_box = QLineEdit(self._default_name())
        form.addRow("Name", self.name_box)
        layout.addLayout(form)

        self.separate_cb = QCheckBox("One dataset per segment")
        self.separate_cb.setToolTip(
            "Off: the segments are laid end to end as a single cut, with the "
            "corners recorded so the viewer can mark them. On: each segment "
            "is listed on its own, starting from distance 0.")
        layout.addWidget(self.separate_cb)

        self.note = QLabel("")
        self.note.setWordWrap(True)
        layout.addWidget(self.note)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.button(QDialogButtonBox.Ok).setText("Plot cut")
        buttons.accepted.connect(self._on_accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

        self._begin_picking()
        self._refresh_overlay()

    # -- points ----------------------------------------------------------
    def _default_name(self) -> str:
        base = f"{self.contour.filename} [arb cut]"
        taken = set(self.contour.existing_names())
        if base not in taken:
            return base
        n = 2
        while f"{base} {n}" in taken:
            n += 1
        return f"{base} {n}"

    def add_point(self, x: float, y: float):
        """A click on the contour fills the first row that is still free; a
        seventh click is ignored rather than silently overwriting one."""
        for enable, x_box, y_box in self.rows:
            if not enable.isChecked():
                for box, value in ((x_box, x), (y_box, y)):
                    box.blockSignals(True)
                    box.setValue(value)
                    box.blockSignals(False)
                enable.setChecked(True)      # triggers the overlay refresh
                return
        self.note.setText("All six points are in use -- clear one to add another.")

    def points(self):
        return [(x_box.value(), y_box.value())
                for enable, x_box, y_box in self.rows if enable.isChecked()]

    def clear_points(self):
        for enable, _, _ in self.rows:
            enable.blockSignals(True)
            enable.setChecked(False)
            enable.blockSignals(False)
        self.note.setText("")
        self._refresh_overlay()

    def _refresh_overlay(self, *_):
        pts = self.points()
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        self.view.show_picked_points(xs, ys)
        self.view.show_overlay_path(xs if len(pts) > 1 else None,
                                    ys if len(pts) > 1 else None)
        if len(pts) < 2:
            self.note.setText("At least two points are needed.")
        else:
            length = float(np.sum(np.hypot(np.diff(xs), np.diff(ys))))
            self.note.setText(f"{len(pts)} points, {len(pts) - 1} segment(s), "
                              f"{length:.4g} {_unit_of(self.contour.defl_label)} long.")

    def _on_accept(self):
        if len(self.points()) < 2:
            QMessageBox.warning(self, "Arbitrary cut",
                                "Place at least two points on the contour first.")
            return
        self.accept()


class FSCorrectionDialog(_PickerDialog):
    """Points along a feature that should be flat, and the polynomial fitted
    through them.

    ``correction.m`` asks for the points with ``getpts`` and fits a
    parabola; here the points are clicked on the live image, listed in a
    table where they can be corrected, and the fit is drawn over the data as
    it is built, so a bad point is obvious before anything is computed.
    """

    def __init__(self, window, view, angle_label: str, energy_label: str):
        super().__init__(window, view, "Fermi-surface correction")
        self.angle_label = angle_label
        self.energy_label = energy_label

        layout = QVBoxLayout(self)
        intro = QLabel(
            "Click along the feature that should be flat (a Fermi edge, a band "
            "bottom). The fitted curve is drawn as you go; every energy column "
            "is then shifted to straighten it.")
        intro.setWordWrap(True)
        layout.addWidget(intro)

        self.table = QTableWidget(0, 2)
        self.table.setHorizontalHeaderLabels([angle_label, energy_label])
        header = self.table.horizontalHeader()
        header.setSectionResizeMode(QHeaderView.Stretch)
        self.table.setMaximumHeight(190)
        self.table.itemChanged.connect(lambda *_: self._refresh_overlay())
        layout.addWidget(self.table)

        row = QHBoxLayout()
        self.reference_button = QPushButton("From a reference (Au)...")
        self.reference_button.setToolTip(
            "Measure the curvature instead of clicking it: fit the Fermi edge "
            "of a gold reference channel by channel and use those points.")
        self.reference_button.clicked.connect(self.from_reference)
        self.remove_button = QPushButton("Remove selected")
        self.remove_button.clicked.connect(self.remove_selected)
        self.clear_button = QPushButton("Clear")
        self.clear_button.clicked.connect(self.clear_points)
        row.addWidget(self.reference_button)
        row.addWidget(self.remove_button)
        row.addWidget(self.clear_button)
        row.addStretch(1)
        row.addWidget(QLabel("Fit order"))
        self.order_box = QSpinBox()
        self.order_box.setRange(1, 4)
        self.order_box.setValue(2)
        self.order_box.setToolTip(
            "Polynomial order. 2 (a parabola) is what the MATLAB tool uses and "
            "what analyser curvature normally looks like; an order-n fit needs "
            "at least n+1 points.")
        self.order_box.valueChanged.connect(lambda *_: self._refresh_overlay())
        row.addWidget(self.order_box)
        layout.addLayout(row)

        form = QFormLayout()
        self.name_box = QLineEdit(self._default_name())
        form.addRow("Name", self.name_box)
        layout.addLayout(form)

        self.note = QLabel("")
        self.note.setWordWrap(True)
        layout.addWidget(self.note)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.button(QDialogButtonBox.Ok).setText("Correct")
        buttons.accepted.connect(self._on_accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

        self._begin_picking()
        self._refresh_overlay()

    def _default_name(self) -> str:
        base = f"{self.window_.filename} [FS corr]"
        taken = set(self.window_.existing_names()) if hasattr(self.window_, "existing_names") else set()
        if base not in taken:
            return base
        n = 2
        while f"{base} {n}" in taken:
            n += 1
        return f"{base} {n}"

    # -- points ----------------------------------------------------------
    def add_point(self, x: float, y: float):
        row = self.table.rowCount()
        self.table.blockSignals(True)
        self.table.insertRow(row)
        for column, value in enumerate((x, y)):
            self.table.setItem(row, column, QTableWidgetItem(f"{value:.5g}"))
        self.table.blockSignals(False)
        self._refresh_overlay()

    def points(self):
        out = []
        for row in range(self.table.rowCount()):
            try:
                x = float(self.table.item(row, 0).text())
                y = float(self.table.item(row, 1).text())
            except (AttributeError, ValueError):
                continue
            out.append((x, y))
        return out

    def remove_selected(self):
        rows = sorted({item.row() for item in self.table.selectedItems()}, reverse=True)
        self.table.blockSignals(True)
        for row in rows:
            self.table.removeRow(row)
        self.table.blockSignals(False)
        self._refresh_overlay()

    def clear_points(self):
        self.table.blockSignals(True)
        self.table.setRowCount(0)
        self.table.blockSignals(False)
        self._refresh_overlay()

    # -- points measured on a reference, rather than clicked ---------------
    def from_reference(self):
        """Fit a reference sample's edge channel by channel and put the
        result in the points table, where it is smoothed by the same
        polynomial as hand-picked points -- and can be inspected and pruned
        the same way."""
        entries = getattr(self.window_, "_dataset_entries", None)
        loader = getattr(self.window_, "_dataset_loader", None)
        if not callable(entries) or not callable(loader):
            QMessageBox.information(self, "Fermi surface from a reference",
                                     "No dataset list is available here.")
            return None
        listed = entries()
        if not listed:
            QMessageBox.information(self, "Fermi surface from a reference",
                                     "Load the reference measurement first.")
            return None

        chooser = AuReferenceDialog(self, listed, loader,
                                     metadata_temperature(self.window_.data))
        if chooser.exec_() != QDialog.Accepted:
            return None
        data = chooser.selected_data()
        settings = chooser.settings()
        if data is None:
            return None
        try:
            angles, energy, frame, _a_label, _e_label = reference_frame(data)
        except ValueError as exc:
            QMessageBox.warning(self, "Fermi surface from a reference", str(exc))
            return None

        result = run_channel_fit(
            self, angles, energy, frame,
            temperature=settings["temperature"], window=(None, None),
            half_width=settings["half_width"], step=settings["step"],
            draw_on=None,
            title=f"Fitting {settings['label']} channel by channel")
        if result is None:
            return None
        angles, ef, ok = result
        if not ok.any():
            QMessageBox.warning(self, "Fermi surface from a reference",
                                 "No channel of that reference could be fitted.")
            return None

        target = np.asarray(self.window_.fs_angle_axis(), dtype=float)
        if (angles.size != target.size
                or not np.allclose(angles, target, rtol=1e-3, atol=1e-6)):
            # A reference taken on a different angle grid still describes the
            # same analyser, so its points are kept as they are -- the
            # polynomial is evaluated on whatever axis the data being
            # corrected has.
            self.note.setText(
                f"note: the reference's {angles.size} channels do not match this "
                f"data's {target.size}; the fitted curve is used through the "
                f"polynomial, which is defined everywhere.")
        self.clear_points()
        self.table.blockSignals(True)
        for angle, value in zip(angles[ok], ef[ok]):
            row = self.table.rowCount()
            self.table.insertRow(row)
            self.table.setItem(row, 0, QTableWidgetItem(f"{angle:.5g}"))
            self.table.setItem(row, 1, QTableWidgetItem(f"{value:.6g}"))
        self.table.blockSignals(False)
        self._refresh_overlay()
        spread = float(np.nanmax(ef[ok]) - np.nanmin(ef[ok]))
        self.window_.statusBar().showMessage(
            f"{int(ok.sum())} channels fitted on {settings['label']}; "
            f"E_F varies by {spread * 1000:.3g} meV across the detector")
        return result

    def coefficients(self):
        order = self.order_box.value()
        pts = self.points()
        if len(pts) < order + 1:
            return None
        return fit_feature(pts, order=order)

    def _refresh_overlay(self, *_):
        pts = self.points()
        self.view.show_picked_points([p[0] for p in pts], [p[1] for p in pts])
        coeffs = self.coefficients()
        if coeffs is None:
            self.view.show_overlay_path(None, None)
            need = self.order_box.value() + 1
            self.note.setText(f"{len(pts)} point(s); an order-{self.order_box.value()} "
                              f"fit needs {need}.")
            return
        axis = np.asarray(self.window_.fs_angle_axis(), dtype=float)
        xs = np.linspace(float(axis.min()), float(axis.max()), 200)
        self.view.show_overlay_path(xs, np.polyval(coeffs, xs), dashed=True)
        fitted = np.polyval(coeffs, axis)
        span = float(np.max(fitted) - np.min(fitted))
        self.note.setText(
            f"{len(pts)} points fitted; the feature spans {span:.4g} "
            f"{_unit_of(self.energy_label) or 'eV'} and will be flattened.")

    def _on_accept(self):
        if self.coefficients() is None:
            QMessageBox.warning(self, "Fermi-surface correction",
                                f"An order-{self.order_box.value()} fit needs at least "
                                f"{self.order_box.value() + 1} points.")
            return
        self.accept()


# --------------------------------------------------------------------------
# Fermi-edge fitting
# --------------------------------------------------------------------------
def reference_frame(data):
    """An (angle, energy) frame and its axes from a dataset, for fitting the
    edge channel by channel.

    A cut is already that; a map is summed over the deflector, because the
    curvature being measured is the analyser's along its slit and every
    deflector position sees the same one.
    """
    scan = data.scan
    if data.kind == "cut":
        return (np.asarray(scan.x, dtype=float), np.asarray(scan.y, dtype=float),
                np.asarray(scan.value, dtype=float),
                scan.labels.get("x", "angle"), scan.labels.get("y", "E"))
    if data.kind in CUBE_KINDS:
        cube = np.asarray(scan.value, dtype=float)
        return (np.asarray(scan.k, dtype=float), np.asarray(scan.z, dtype=float),
                np.nansum(cube, axis=0),
                scan.labels.get("k", "angle"), scan.labels.get("z", "E"))
    raise ValueError(f"{data.kind} has no angle-vs-energy frame to fit")


def metadata_temperature(data, default: float = 30.0) -> float:
    """The sample temperature the file recorded, for holding T at.

    The fields are named differently from one beamline version to the next,
    so anything that looks like a sample temperature in kelvin is accepted;
    a reading far outside a cryostat's range is ignored rather than fitted
    from, since a wrong temperature quietly becomes a wrong resolution.
    """
    info = dict(getattr(data.scan, "info", {}) or {})
    info.update(getattr(data.scan, "fourd_info", {}) or {})
    for key, value in info.items():
        name = str(key).lower()
        if "temp" not in name or "setpoint" in name:
            continue
        try:
            number = float(np.asarray(value).reshape(-1)[0])
        except (TypeError, ValueError, IndexError):
            continue
        if 1.0 <= number <= 1000.0:
            return number
    return float(default)


# --------------------------------------------------------------------------
# k conversion for a single Cut
# --------------------------------------------------------------------------
def cut_deflector_angle(data, default: float = 0.0) -> float:
    """The deflector position a cut was taken at, as the file recorded it.

    Half the geometry a cut's k conversion needs comes for free: the file
    knows which deflector position this was, even though the axis was
    dropped when the single-position scan became a plain 2-D cut. The other
    half -- where Gamma is -- cannot come from the file, because where
    normal emission sits is a beamline convention, not a measurement.
    """
    info = dict(getattr(data.scan, "info", {}) or {})
    for key in ("DeflectorAngle_deg", "kcut.deflector_deg"):
        if key in info:
            try:
                return float(np.asarray(info[key]).reshape(-1)[0])
            except (TypeError, ValueError, IndexError):
                pass
    return float(default)


def cut_k_conversion_blocked(data):
    """Why this cut cannot be converted from angle to momentum, or None.

    The conversion assumes what a measured Cut is: intensity against the
    angle along the slit, at one deflector position. Several things in this
    program are also "cuts" without being that -- one already in momentum,
    one taken along an arbitrary path through a map, one sliced out of a map
    along the deflector axis -- and converting those would produce a picture
    whose axis means nothing. Checked from the metadata the producers leave
    behind rather than guessed from the numbers.
    """
    info = dict(getattr(data.scan, "info", {}) or {})
    if any(str(key).startswith("kcut.") for key in info):
        return "this cut is already in momentum"
    if any(str(key).startswith("arbcut.") for key in info):
        return ("this cut was taken along a path through a map, not at one "
                "deflector position — convert the map instead, then cut it")
    panel = str(info.get("slice.panel", ""))
    if panel and "slit" not in panel.lower() and "cut" not in panel.lower():
        return f"this is a slice through a map ({panel}), not a slit cut"
    label = str(getattr(data.scan, "labels", {}).get("x", ""))
    if "Å" in label:
        return "this cut's angle axis is already a momentum"
    return None


def sample_rotations(data) -> dict:
    """The manipulator angles the file recorded, for comparing a cut with
    the map whose alignment it is about to borrow."""
    motors = dict(getattr(data.scan, "fourd_info", {}) or {})
    return {name: float(motors[name]) for name in ("SRn", "SRz")
            if name in motors and isinstance(motors[name], (int, float, np.floating))}


class CutKConversionDialog(QDialog):
    """Convert one Cut from angle to momentum.

    The conversion needs two things the cut itself cannot supply: which
    deflector position it was taken at (the file knows) and at which angles
    Gamma sits (the file cannot know -- normal emission is a beamline
    convention). The second is what the three buttons at the top are for,
    in descending order of trustworthiness: inherit it from a Map already
    converted in this session, take the file's own geometry as a starting
    point, or type it.

    Deliberately *not* offered: picking a "centre line" on the cut. That is
    what the lab's MATLAB dialog does, and it conflates where normal
    emission is (which enters inside the sine, and so sets the k axis's
    scale) with where Gamma is (a translation applied afterwards). Getting
    the first from the second stretches one half of the cut and squashes the
    other. Here they are separate: geometry fixes the whole trajectory, and
    where the axis reads zero is chosen below.
    """

    def __init__(self, window: "CutWindow"):
        super().__init__(window)
        self.window_ = window
        self.setWindowTitle("Cut k conversion")
        self.setModal(False)
        scan = window.data.scan
        self.slit = np.asarray(scan.x, dtype=float)
        self.energy = np.asarray(scan.y, dtype=float)

        # Three groups of settings plus a readout is a tall dialog, and a
        # laptop screen is not. In a scroll area it opens at whatever height
        # there is instead of squashing its own rows until the spin boxes
        # are unreadable.
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        scroller = QScrollArea()
        scroller.setWidgetResizable(True)
        scroller.setFrameShape(QFrame.NoFrame)
        content = QWidget()
        scroller.setWidget(content)
        outer.addWidget(scroller, 1)

        layout = QVBoxLayout(content)
        intro = QLabel(
            "A cut is a line in the (kx, ky) plane and generally does not "
            "pass through Γ, so the conversion needs to be told where "
            "Γ is. Inherit it from a converted map, or type the angles.")
        intro.setWordWrap(True)
        layout.addWidget(intro)

        # -- where Gamma is --------------------------------------------------
        gamma_group = QGroupBox("Reference Γ")
        gamma = QFormLayout(gamma_group)
        self.map_button = QPushButton("From a converted map...")
        self.map_button.setToolTip(
            "Take the angles from a Map already converted in this session. "
            "Exact if the cut was measured in the same alignment -- the "
            "manipulator angles of both are compared and any difference is "
            "reported.")
        self.map_button.clicked.connect(self.from_converted_map)
        gamma.addRow("", self.map_button)

        self.gamma_defl = self._angle_box(
            "Deflector angle at which Γ was observed. This is the Map "
            "conversion's \"theta offset\".")
        self.gamma_slit = self._angle_box(
            "Angle along the slit at which Γ was observed. This is the "
            "Map conversion's \"phi offset\".")
        self.azimuth_box = self._angle_box(
            "Sample rotation about the surface normal, as in the Map "
            "conversion.")
        gamma.addRow("Γ deflector angle (deg)", self.gamma_defl)
        gamma.addRow("Γ slit angle (deg)", self.gamma_slit)
        gamma.addRow("Sample rotation (deg)", self.azimuth_box)
        self.source_note = QLabel("Typed by hand — nothing inherited yet.")
        self.source_note.setWordWrap(True)
        self.source_note.setMinimumHeight(34)
        gamma.addRow("", self.source_note)
        layout.addWidget(gamma_group)

        # -- what this cut is ------------------------------------------------
        cut_group = QGroupBox("This cut")
        cut = QFormLayout(cut_group)
        self.deflector_box = self._angle_box(
            "The deflector position this cut was taken at. Read from the "
            "file; edit it only if the file recorded it wrongly.")
        recorded = cut_deflector_angle(window.data, default=float("nan"))
        self.deflector_box.setValue(0.0 if not np.isfinite(recorded) else recorded)
        cut.addRow("Deflector angle (deg)", self.deflector_box)

        self.energy_offset_box = QDoubleSpinBox()
        self.energy_offset_box.setDecimals(4)
        self.energy_offset_box.setRange(-10000.0, 10000.0)
        self.energy_offset_box.setSingleStep(0.1)
        self.energy_offset_box.setToolTip(
            "Added to the energy axis before converting, for when it is not "
            "already a true kinetic energy.")
        cut.addRow("Energy offset (eV)", self.energy_offset_box)
        self.file_note = QLabel(
            "Deflector angle read from the file." if np.isfinite(recorded)
            else "This file records no deflector angle — type it, or "
                 "leave it at 0 if the cut was taken at the deflector's zero.")
        self.file_note.setWordWrap(True)
        cut.addRow("", self.file_note)
        layout.addWidget(cut_group)

        # -- the output -------------------------------------------------------
        out_group = QGroupBox("Output")
        out = QFormLayout(out_group)
        self.name_box = QLineEdit(self._default_name())
        out.addRow("Name", self.name_box)

        self.by_count = QRadioButton("by element number")
        self.by_resolution = QRadioButton("by resolution")
        self.by_count.setChecked(True)
        mode_row = QWidget()
        mode_box = QHBoxLayout(mode_row)
        mode_box.setContentsMargins(0, 0, 0, 0)
        mode_box.addWidget(self.by_count)
        mode_box.addWidget(self.by_resolution)
        mode_box.addStretch(1)
        out.addRow("Output grid", mode_row)

        n_e = self.energy.size
        self.nk_box = self._count_box(
            max(self.slit.size, 2),
            f"Number of k points. Defaults to the measured {self.slit.size} "
            f"angles.")
        self.ne_box = self._count_box(
            n_e, f"Number of energy points. Defaults to the measured {n_e}.")
        out.addRow("k element no.", self.nk_box)
        out.addRow("E element no.", self.ne_box)

        self.dk_box = self._resolution_box(0.01, 4, "Å⁻¹",
                                           "Spacing between k points.")
        step_meV = abs(float(self.energy[1] - self.energy[0])) * 1000 if n_e > 1 else 10.0
        self.de_box = self._resolution_box(
            max(round(step_meV, 3), 0.001), 3, "meV",
            f"Spacing between energy points. The measurement's own step is "
            f"{step_meV:.4g} meV.")
        out.addRow("k resolution (Å⁻¹)", self.dk_box)
        out.addRow("E resolution (meV)", self.de_box)

        self.radial_box = QCheckBox("Measure from Γ itself")
        self.radial_box.setToolTip(
            "Off: the axis is the momentum along the cut, zeroed where the "
            "cut passes closest to Γ -- the coordinate a band disperses "
            "in.\nOn: the axis is the distance from Γ, which leaves a "
            "real gap of ±k⊥ around zero because the cut never "
            "reaches Γ. Honest, but half the picture is then empty.")
        out.addRow("", self.radial_box)

        self.trim_box = QCheckBox("Trim to the range every energy covers")
        self.trim_box.setToolTip(
            "k scales with sqrt(E), so the lowest-energy rows do not reach as "
            "far as the highest. Off, the ends are missing for those rows; "
            "on, the output is cut back to the range all of them share.")
        out.addRow("", self.trim_box)

        self.grid_note = QLabel("")
        self.grid_note.setWordWrap(True)
        self.grid_note.setMinimumHeight(34)
        out.addRow("", self.grid_note)
        layout.addWidget(out_group)

        self.geometry_note = QLabel("")
        self.geometry_note.setWordWrap(True)
        self.geometry_note.setMinimumHeight(52)
        self.geometry_note.setTextInteractionFlags(Qt.TextSelectableByMouse)
        layout.addWidget(self.geometry_note)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.button(QDialogButtonBox.Ok).setText("Convert")
        buttons.button(QDialogButtonBox.Cancel).setText("Close")
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        # Outside the scroll area: Convert must not be something to scroll
        # down for.
        buttons.setContentsMargins(8, 4, 8, 8)
        outer.addWidget(buttons)
        self.resize(self.sizeHint().width() + 40, 700)

        for radio in (self.by_count, self.by_resolution):
            radio.toggled.connect(self._update_grid_mode)
        for box in (self.gamma_defl, self.gamma_slit, self.azimuth_box,
                    self.deflector_box, self.energy_offset_box,
                    self.nk_box, self.ne_box, self.dk_box, self.de_box):
            box.valueChanged.connect(self._refresh)
        for box in (self.radial_box, self.trim_box):
            box.toggled.connect(self._refresh)
        self._update_grid_mode()

    # -- small builders ----------------------------------------------------
    @staticmethod
    def _angle_box(tip: str):
        box = QDoubleSpinBox()
        box.setDecimals(3)
        box.setRange(-180.0, 180.0)
        box.setSingleStep(0.5)
        box.setToolTip(tip)
        return box

    @staticmethod
    def _count_box(value: int, tip: str):
        box = QSpinBox()
        box.setRange(2, 20000)
        box.setValue(int(value))
        box.setToolTip(tip)
        return box

    @staticmethod
    def _resolution_box(value, decimals, unit, tip):
        box = QDoubleSpinBox()
        box.setDecimals(decimals)
        box.setRange(10.0 ** -decimals, 1e6)
        box.setSingleStep(10.0 ** -(decimals - 1))
        box.setValue(value)
        box.setToolTip(f"{tip} In {unit}.")
        return box

    def _default_name(self) -> str:
        base = f"{self.window_.filename} [k]"
        taken = set(self.window_.existing_names())
        if base not in taken:
            return base
        n = 2
        while f"{base} {n}" in taken:
            n += 1
        return f"{base} {n}"

    def output_name(self) -> str:
        return self.name_box.text().strip() or self._default_name()

    # -- inheriting an alignment -------------------------------------------
    def from_converted_map(self):
        """Take Γ's angles from a Map converted earlier in this session.

        The Map conversion's theta/phi offsets *are* Γ's angles -- the
        user picked the point that became k = (0, 0) -- so they transfer to
        any cut taken in the same alignment without being measured again.
        """
        entries = getattr(self.window_, "_dataset_entries", None)
        loader = getattr(self.window_, "_dataset_loader", None)
        if not callable(entries) or not callable(loader):
            QMessageBox.information(self, "Cut k conversion",
                                     "No dataset list is available here.")
            return None
        candidates = [(label, key) for label, key in entries()
                      if label.rstrip().endswith("[k_map]")]
        if not candidates:
            QMessageBox.information(
                self, "Cut k conversion",
                "No converted map is in the list yet. Convert a map of this "
                "sample first ('Map k conversion'), then come back — or "
                "type Γ's angles below.")
            return None
        labels = [label for label, _ in candidates]
        choice, ok = QInputDialog.getItem(self, "Inherit Γ from a map",
                                           "Converted map:", labels, 0, False)
        if not ok:
            return None
        key = candidates[labels.index(choice)][1]
        try:
            data = loader(key)
        except Exception as exc:
            QMessageBox.warning(self, "Cut k conversion",
                                 f"That map could not be read:\n{exc}")
            return None
        return self.adopt_alignment(data, choice)

    def adopt_alignment(self, data, label: str = ""):
        """Copy a converted map's k-space origin into the boxes, and say how
        far the two measurements' manipulator angles agree."""
        info = dict(getattr(data.scan, "info", {}) or {})
        try:
            theta = float(info["kconv.theta_offset_deg"])
            phi = float(info["kconv.phi_offset_deg"])
        except (KeyError, TypeError, ValueError):
            QMessageBox.warning(
                self, "Cut k conversion",
                "That dataset does not record a k conversion, so there is no "
                "Γ to inherit from it.")
            return None
        azimuth = float(info.get("kconv.azimuth_deg", 0.0))
        for box, value in ((self.gamma_defl, theta), (self.gamma_slit, phi),
                           (self.azimuth_box, azimuth)):
            box.blockSignals(True)
            box.setValue(value)
            box.blockSignals(False)

        mine = sample_rotations(self.window_.data)
        theirs = sample_rotations(data)
        shared = sorted(set(mine) & set(theirs))
        drift = {name: mine[name] - theirs[name] for name in shared
                 if abs(mine[name] - theirs[name]) > 0.01}
        note = f"Γ taken from {label or 'a converted map'}."
        if drift:
            note += (" ⚠ The manipulator moved between the two: "
                     + ", ".join(f"{name} by {value:+.3f}°"
                                 for name, value in drift.items())
                     + ". The inherited angles are only right if that motion "
                       "does not change where Γ sits — check, or "
                       "correct the angles by hand.")
        elif shared:
            note += (" The manipulator was in the same position for both ("
                     + ", ".join(f"{name} {mine[name]:.3f}°"
                                 for name in shared) + ").")
        self.source_note.setText(note)
        self._refresh()
        return theta, phi, azimuth

    # -- the output grid ----------------------------------------------------
    def _update_grid_mode(self, *_):
        by_count = self.by_count.isChecked()
        for box in (self.nk_box, self.ne_box):
            box.setEnabled(by_count)
        for box in (self.dk_box, self.de_box):
            box.setEnabled(not by_count)
        self._refresh()

    def geometry(self) -> dict:
        return {
            "deflector_deg": self.deflector_box.value(),
            "gamma_deflector_deg": self.gamma_defl.value(),
            "gamma_slit_deg": self.gamma_slit.value(),
            "azimuth_deg": self.azimuth_box.value(),
            "energy_offset_eV": self.energy_offset_box.value(),
        }

    def k_extent(self):
        from tools.cutk import cut_extent
        try:
            return cut_extent(self.slit, self.energy,
                              radial=self.radial_box.isChecked(),
                              trim=self.trim_box.isChecked(), **self.geometry())
        except Exception:
            return None

    def grid_counts(self):
        """(n_k, n_E), whichever way the user asked for them."""
        if self.by_count.isChecked():
            return int(self.nk_box.value()), int(self.ne_box.value())
        extent = self.k_extent()
        if extent is None:
            raise ValueError("the settings do not define a k range yet")
        lo, hi = extent
        n_k = int(round((hi - lo) / self.dk_box.value())) + 1
        span_meV = abs(float(self.energy[-1] - self.energy[0])) * 1000.0
        n_e = int(round(span_meV / self.de_box.value())) + 1
        return max(n_k, 2), max(n_e, 2)

    def settings(self) -> dict:
        n_k, n_energy = self.grid_counts()
        settings = dict(self.geometry())
        settings.update({"n_k": n_k, "n_energy": n_energy,
                         "radial": self.radial_box.isChecked(),
                         "trim": self.trim_box.isChecked()})
        return settings

    # -- what the settings would give --------------------------------------
    def _refresh(self, *_):
        try:
            n_k, n_e = self.grid_counts()
            extent = self.k_extent()
        except Exception as exc:
            self.grid_note.setText(str(exc))
            self.geometry_note.setText("")
            return
        if extent is None:
            self.grid_note.setText("")
            self.geometry_note.setText("")
            return
        lo, hi = extent
        if self.by_count.isChecked():
            span_meV = abs(float(self.energy[-1] - self.energy[0])) * 1000.0
            self.grid_note.setText(
                f"→ {(hi - lo) / max(n_k - 1, 1):.4g} Å⁻¹ and "
                f"{span_meV / max(n_e - 1, 1):.4g} meV per point")
        else:
            self.grid_note.setText(f"→ {n_k} x {n_e} points")
        self.geometry_note.setText(self.describe_geometry(lo, hi))

    def describe_geometry(self, lo: float, hi: float) -> str:
        """Everything the axes cannot say: how far off Γ this cut runs,
        and how much of an idealisation "a line at constant k⊥" is for it.

        The conversion itself follows the exact arc, so the variation below
        is not an error in the result -- it is how much the *description* of
        the cut as one number is rounded. Said that way round, because the
        obvious reading of a number like that is "my conversion is wrong by
        this much", and it is not.
        """
        from tools.cutk import cut_momenta
        try:
            _k_par, k_perp = cut_momenta(self.slit, self.energy, **self.geometry())
        except Exception as exc:
            return str(exc)
        perp = np.abs(k_perp)
        try:
            n_k = self.grid_counts()[0]
        except Exception:
            n_k = max(self.slit.size, 2)
        pixel = abs(hi - lo) / max(n_k - 1, 1)

        if perp.max() < 1e-6:
            return (f"k from {lo:+.4g} to {hi:+.4g} Å⁻¹. This cut "
                    f"runs through Γ, so k∥ is the distance from it.")
        spread = float(perp.max() - perp.min())
        text = (f"k from {lo:+.4g} to {hi:+.4g} Å⁻¹. The cut runs "
                f"{perp.mean():.4g} Å⁻¹ from Γ, varying by "
                f"{spread:.4g} Å⁻¹ along it and across the energy "
                f"window — the conversion follows the exact arc, so that is "
                f"how much \"a line at constant k⊥\" is an idealisation, not "
                f"an error. One output point is {pixel:.4g} Å⁻¹.")
        if spread > 0.05:
            text += (" ⚠ At that size the cut is a visibly curved path "
                     "through the zone; compare it with a calculation along "
                     "the same arc, not a straight line.")
        return text

    def closeEvent(self, event):
        self.window_._cut_k_dialog = None
        super().closeEvent(event)


class FermiFitDialog(_PickerDialog):
    """Fit the Fermi edge of the cut on screen.

    The EDC is the cut summed over an angle range, and the fit runs over an
    energy range: both are the selection box on the image, so the region can
    be dragged there or typed here, whichever is quicker.

    The temperature is **held** by default at what the file recorded. On a
    single edge the temperature and the resolution are not separable -- the
    measured width is sqrt((3.53kT)^2 + FWHM^2) -- so fitting both gives two
    numbers that mean nothing apart. Unticking "hold" is allowed, and the
    summary then prints their correlation, which is how you can tell.

    The two rarely-needed operations live in "Advanced" at the bottom:
    dividing the Fermi cut-off out (superconducting gaps, mostly) and
    fitting the edge channel by channel.
    """

    def __init__(self, window):
        super().__init__(window, window.view, "Fermi level fitting")
        self.window_ = window
        self.fit_result = None
        self._angle_axis, self._energy_axis, self._frame, \
            self._angle_label, self._energy_label = reference_frame(window.data)

        layout = QVBoxLayout(self)
        intro = QLabel(
            "The EDC is this cut summed over an angle range; the fit runs over "
            "an energy range. Drag the selection box on the image or type the "
            "numbers here -- they follow each other.")
        intro.setWordWrap(True)
        layout.addWidget(intro)

        # -- the region -----------------------------------------------------
        region = QFormLayout()
        self.angle_from, self.angle_to = self._range_boxes(self._angle_axis)
        self.energy_from, self.energy_to = self._range_boxes(self._energy_axis)
        region.addRow(f"Integrate {self._angle_label}", self._pair(self.angle_from,
                                                                   self.angle_to))
        region.addRow(f"Fit window {self._energy_label}",
                      self._pair(self.energy_from, self.energy_to))
        self.box_button = QPushButton("Show the selection box")
        self.box_button.setCheckable(True)
        self.box_button.setToolTip(
            "Put the box on the image; its horizontal extent is the angle range "
            "and its vertical extent the energy window.")
        self.box_button.toggled.connect(self._toggle_box)
        region.addRow("", self.box_button)
        layout.addLayout(region)

        # -- parameters -----------------------------------------------------
        self.table = QTableWidget(len(PARAMETERS), 3)
        self.table.setHorizontalHeaderLabels(["Value", "Hold", "Uncertainty"])
        self.table.setVerticalHeaderLabels([PARAMETER_LABELS[p] for p in PARAMETERS])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.table.verticalHeader().setDefaultSectionSize(26)
        # Tall enough for all seven rows: a table showing two of them and a
        # scrollbar hides half the model from the person judging the fit.
        self.table.setMinimumHeight(7 * 26 + 30)
        self.table.setMaximumHeight(7 * 26 + 34)
        self.value_boxes, self.hold_boxes, self.error_labels = {}, {}, {}
        for row, name in enumerate(PARAMETERS):
            box = QDoubleSpinBox()
            box.setDecimals(6)
            box.setRange(-1e9, 1e9)
            box.setButtonSymbols(QDoubleSpinBox.NoButtons)
            self.table.setCellWidget(row, 0, box)
            self.value_boxes[name] = box

            hold = QCheckBox()
            hold.setChecked(name == "temperature")
            holder = QWidget()
            hbox = QHBoxLayout(holder)
            hbox.setContentsMargins(0, 0, 0, 0)
            hbox.addStretch(1)
            hbox.addWidget(hold)
            hbox.addStretch(1)
            self.table.setCellWidget(row, 1, holder)
            self.hold_boxes[name] = hold

            label = QLabel("")
            label.setAlignment(Qt.AlignCenter)
            self.table.setCellWidget(row, 2, label)
            self.error_labels[name] = label
        self.hold_boxes["temperature"].setToolTip(
            "Held by default: on one edge the temperature cannot be told apart "
            "from the resolution. Untick it only if the resolution is known "
            "independently -- the summary will print their correlation.")
        layout.addWidget(self.table)

        row = QHBoxLayout()
        self.estimate_button = QPushButton("Estimate")
        self.estimate_button.setToolTip("Read starting values off the data.")
        self.estimate_button.clicked.connect(self.estimate)
        self.fit_button = QPushButton("Fit")
        self.fit_button.clicked.connect(self.run_fit)
        # Next to Fit, because it is what you press straight afterwards: the
        # fit and the thing the fit is for belong on one row.
        self.offset_button = QPushButton("Offset energy axis")
        self.offset_button.setToolTip(
            "List a copy of this dataset with the energy axis shifted so the "
            "fitted E_F sits at zero.")
        self.offset_button.clicked.connect(self.offset_energy_axis)
        self.offset_button.setEnabled(False)
        row.addWidget(self.estimate_button)
        row.addWidget(self.fit_button)
        row.addWidget(self.offset_button)
        row.addStretch(1)
        layout.addLayout(row)

        # -- the picture ----------------------------------------------------
        self.plot = pg.PlotWidget()
        self.plot.setMinimumHeight(220)
        self.plot.setLabel("bottom", self._energy_label)
        self.plot.setLabel("left", "counts")
        strip_stock_menu(self.plot.getPlotItem())
        self.data_curve = self.plot.plot(pen=pg.mkPen("#1f6f8b", width=2), name="EDC")
        self.model_curve = self.plot.plot(pen=pg.mkPen("#bd4921", width=2))
        self.ef_line = pg.InfiniteLine(angle=90, pen=pg.mkPen("#2f7a3f", width=2,
                                                              style=Qt.DashLine))
        self.ef_line.setVisible(False)
        self.plot.addItem(self.ef_line, ignoreBounds=True)
        layout.addWidget(self.plot, stretch=1)

        self.summary = QLabel("")
        self.summary.setWordWrap(True)
        self.summary.setTextInteractionFlags(Qt.TextSelectableByMouse)
        layout.addWidget(self.summary)

        # -- what to do with the answer -------------------------------------
        actions = QHBoxLayout()
        actions.addStretch(1)
        close = QPushButton("Close")
        close.clicked.connect(self.reject)
        actions.addWidget(close)

        advanced = QGroupBox("Advanced -- superconducting gaps and reference samples")
        advanced.setCheckable(True)
        advanced.setChecked(False)
        adv = QFormLayout(advanced)
        self.cutoff_box = QDoubleSpinBox()
        self.cutoff_box.setRange(0.5, 20.0)
        self.cutoff_box.setValue(4.0)
        self.cutoff_box.setSingleStep(0.5)
        self.cutoff_box.setToolTip(
            "Everything above E_F + this many kT becomes missing: the divisor "
            "there is a small number known only to the noise.")
        self.divide_button = QPushButton("Divide the Fermi cut-off out")
        self.divide_button.setToolTip(
            "Divide by the resolution-broadened Fermi function only -- not by "
            "the fitted density of states, and with the background taken out "
            "first so it is not inflated near E_F.")
        self.divide_button.clicked.connect(self.divide_out)
        adv.addRow("Cut off above (kT)", self.cutoff_box)
        adv.addRow("", self.divide_button)

        self.channel_width = QSpinBox()
        self.channel_width.setRange(0, 200)
        self.channel_width.setValue(2)
        self.channel_width.setToolTip(
            "Add this many neighbouring channels on each side before fitting, "
            "for a reference too thin to fit channel by channel.")
        self.channel_step = QSpinBox()
        self.channel_step.setRange(1, 100)
        self.channel_step.setValue(1)
        self.channel_step.setToolTip(
            "Fit every Nth channel. The curvature is smooth, so there is "
            "rarely anything to gain from fitting all of them.")
        self.channels_button = QPushButton("Fit E_F channel by channel")
        self.channels_button.setToolTip(
            "Fit the edge in each detector channel of this cut and plot E_F "
            "against angle -- the curvature of the Fermi surface, measured "
            "rather than clicked.")
        self.channels_button.clicked.connect(self.fit_by_channel)
        adv.addRow("Neighbours (+/- channels)", self.channel_width)
        adv.addRow("Fit every Nth channel", self.channel_step)
        adv.addRow("", self.channels_button)
        layout.addWidget(advanced)
        layout.addLayout(actions)

        for box in (self.angle_from, self.angle_to, self.energy_from, self.energy_to):
            box.valueChanged.connect(self._on_range_typed)
        window.view.selectionChanged.connect(self._on_box_moved)
        self.value_boxes["temperature"].setValue(metadata_temperature(window.data))
        self.estimate()

    # -- small builders ---------------------------------------------------
    @staticmethod
    def _pair(first, second):
        holder = QWidget()
        box = QHBoxLayout(holder)
        box.setContentsMargins(0, 0, 0, 0)
        box.addWidget(first)
        box.addWidget(QLabel("to"))
        box.addWidget(second)
        box.addStretch(1)
        return holder

    @staticmethod
    def _range_boxes(axis):
        lo, hi = float(np.min(axis)), float(np.max(axis))
        boxes = []
        for value in (lo, hi):
            box = QDoubleSpinBox()
            box.setDecimals(4)
            box.setRange(lo, hi)
            box.setValue(value)
            box.setKeyboardTracking(False)
            box.setButtonSymbols(QDoubleSpinBox.NoButtons)
            box.setMaximumWidth(96)
            boxes.append(box)
        return boxes

    # -- the region and the EDC -------------------------------------------
    def _toggle_box(self, on: bool):
        view = self.window_.view
        view.set_selection_visible(bool(on))
        if on:
            view.set_selection_corners(self.angle_from.value(), self.energy_from.value(),
                                       self.angle_to.value(), self.energy_to.value())
        self.box_button.setText("Hide the selection box" if on
                                else "Show the selection box")

    def _on_box_moved(self):
        corners = self.window_.view.selection_corners()
        if corners is None:
            return
        x0, y0, x1, y1 = corners
        for box, value in ((self.angle_from, x0), (self.angle_to, x1),
                           (self.energy_from, y0), (self.energy_to, y1)):
            box.blockSignals(True)
            box.setValue(float(np.clip(value, box.minimum(), box.maximum())))
            box.blockSignals(False)
        self.refresh_edc()

    def _on_range_typed(self, *_):
        if self.box_button.isChecked():
            self.window_.view.set_selection_corners(
                self.angle_from.value(), self.energy_from.value(),
                self.angle_to.value(), self.energy_to.value())
        self.refresh_edc()

    def edc(self):
        """(energy, intensity) for the angle range on screen."""
        lo, hi = sorted((self.angle_from.value(), self.angle_to.value()))
        inside = (self._angle_axis >= lo) & (self._angle_axis <= hi)
        if not inside.any():
            inside = np.zeros_like(self._angle_axis, dtype=bool)
            inside[int(np.argmin(np.abs(self._angle_axis - lo)))] = True
        return self._energy_axis, np.nansum(self._frame[inside], axis=0)

    def window_range(self):
        lo, hi = sorted((self.energy_from.value(), self.energy_to.value()))
        return (lo, hi)

    def windowed_edc(self):
        """The part of the EDC the fit actually sees -- what the plot shows.

        Drawing the whole energy axis would put the edge in a corner and let
        the intensity scale be set by data nobody is fitting; restricted to
        the window, the edge fills the panel and the residual is visible.
        """
        energy, intensity = self.edc()
        lo, hi = self.window_range()
        inside = (energy >= lo) & (energy <= hi)
        if inside.sum() < 2:
            return energy, intensity
        return energy[inside], intensity[inside]

    def refresh_edc(self):
        energy, intensity = self.windowed_edc()
        self.data_curve.setData(energy, intensity)
        if energy.size:
            self.plot.setXRange(float(energy.min()), float(energy.max()), padding=0.02)
        if self.fit_result is not None:
            self.model_curve.setData(energy, self.fit_result.model(energy))
        else:
            self.model_curve.setData([], [])

    # -- fitting ----------------------------------------------------------
    def parameters(self) -> dict:
        return {name: box.value() for name, box in self.value_boxes.items()}

    def held(self) -> tuple:
        return tuple(name for name, box in self.hold_boxes.items() if box.isChecked())

    def _show_values(self, values: dict, errors: dict = None):
        for name, box in self.value_boxes.items():
            box.blockSignals(True)
            box.setValue(float(values[name]))
            box.blockSignals(False)
            error = (errors or {}).get(name)
            self.error_labels[name].setText(
                "" if error is None or not np.isfinite(error)
                else ("held" if name in self.held() and error == 0 else f"± {error:.4g}"))

    def estimate(self):
        """Starting values read off the data."""
        energy, intensity = self.edc()
        lo, hi = self.window_range()
        inside = (energy >= lo) & (energy <= hi)
        try:
            guess = initial_guess(energy[inside], intensity[inside],
                                  temperature=self.value_boxes["temperature"].value())
        except ValueError as exc:
            self.summary.setText(str(exc))
            return None
        guess["temperature"] = self.value_boxes["temperature"].value()
        self._show_values(guess)
        self.refresh_edc()
        shown = self.windowed_edc()[0]
        self.model_curve.setData(shown, fermi_edge_model(shown, **guess))
        self.summary.setText("Starting values estimated from the data. Press Fit.")
        return guess

    def run_fit(self):
        energy, intensity = self.edc()
        try:
            result = fit_fermi_edge(
                energy, intensity, start=self.parameters(), fixed=self.held(),
                temperature=self.value_boxes["temperature"].value(),
                window=self.window_range())
        except Exception as exc:
            QMessageBox.warning(self, "Fermi level fitting", f"The fit failed:\n{exc}")
            return None
        self.fit_result = result
        self._show_values(result.values, result.errors)
        self.refresh_edc()
        self.ef_line.setPos(result.values["ef"])
        self.ef_line.setVisible(True)
        self.summary.setText(result.summary())
        self.offset_button.setEnabled(True)
        self.window_.statusBar().showMessage(
            f"E_F = {result.values['ef']:.5g} eV, resolution "
            f"{result.values['resolution'] * 1000:.3g} meV FWHM")
        return result

    def fit_info(self) -> dict:
        """The fit as metadata, to travel with anything it produces."""
        if self.fit_result is None:
            return {}
        result = self.fit_result
        info = {f"{name}": float(result.values[name]) for name in PARAMETERS}
        info.update({f"{name}_err": float(result.errors.get(name, float("nan")))
                     for name in PARAMETERS})
        info["held"] = ", ".join(result.fixed) or "(nothing)"
        info["reduced_chi2"] = float(result.reduced_chi2)
        info["combined_width_eV"] = float(result.combined_width)
        info["window_eV"] = f"{result.window[0]:.6g}..{result.window[1]:.6g}"
        info["angle_range"] = (f"{min(self.angle_from.value(), self.angle_to.value()):.6g}"
                               f"..{max(self.angle_from.value(), self.angle_to.value()):.6g}")
        pair = result.correlation.get(("temperature", "resolution"))
        if pair is not None:
            info["corr_T_resolution"] = float(pair)
        return info

    # -- what comes out of it ---------------------------------------------
    def offset_energy_axis(self):
        if self.fit_result is None:
            return None
        return self.window_.offset_energy_axis(self.fit_result.values["ef"],
                                               self.fit_info())

    def divide_out(self):
        if self.fit_result is None:
            QMessageBox.information(self, "Fermi level fitting",
                                     "Fit the edge first.")
            return None
        return self.window_.divide_fermi_cutoff(self.fit_result,
                                                 self.cutoff_box.value(),
                                                 self.fit_info())

    def fit_by_channel(self):
        curve = self.window_.fit_channels_here(
            temperature=self.value_boxes["temperature"].value(),
            fixed=self.held(), window=self.window_range(),
            half_width=self.channel_width.value(), step=self.channel_step.value(),
            start=self.parameters())
        if curve is None:
            return None
        angles, ef, ok = curve
        if not ok.any():
            self.summary.setText("No channel could be fitted.")
            return None
        spread = float(np.nanmax(ef[ok]) - np.nanmin(ef[ok]))
        self.summary.setText(
            f"{int(ok.sum())} of {angles.size} channels fitted; E_F varies by "
            f"{spread * 1000:.3g} meV across the detector. The curve is drawn "
            f"over the cut -- use “FS correction” to flatten it.")
        return curve


def run_channel_fit(parent, angles, energy, frame, *, temperature, fixed=("temperature",),
                    window=(None, None), half_width=0, step=1, start=None,
                    draw_on=None, title="Fitting the edge channel by channel"):
    """Fit every (Nth) channel of an angle-vs-energy frame, with a progress
    dialog that can be cancelled, and optionally draw the result.

    Hundreds of fits take tens of seconds, so this is the one operation in
    the program that has to say how far along it is; cancelling keeps
    whatever has been fitted so far, which is usually enough to see whether
    the settings were right.

    This one still runs on the GUI thread, unlike the processing panels and
    the k conversion (see ``ui.jobs``), because it hands its answer back as
    a return value and both callers use it on the next line -- turning that
    into a callback would restructure the Fermi-level workflow for a fit
    that is already interruptible. What it does not do any more is leave the
    door open while it pumps events: the dialog is **application** modal, so
    ``processEvents`` below cannot dispatch a click on some other window
    into a second copy of this same code. That was the real hazard in the
    pattern, rather than the pumping itself.
    """
    progress = QProgressDialog(title, "Cancel", 0, 100, parent)
    progress.setWindowModality(Qt.ApplicationModal)
    progress.setMinimumDuration(400)
    progress.setValue(0)

    def report(done, total):
        progress.setMaximum(int(total))
        progress.setValue(int(done))
        QApplication.processEvents()
        return not progress.wasCanceled()

    try:
        ef, err, ok = fit_channels(angles, energy, frame, half_width=half_width,
                                   step=step, temperature=temperature, fixed=fixed,
                                   window=window, start=start, progress=report)
    except Exception as exc:
        progress.close()
        QMessageBox.warning(parent, "Fermi level fitting",
                             f"Channel fitting failed:\n{exc}")
        return None
    progress.close()

    if draw_on is not None and ok.any():
        draw_on.show_picked_points(angles[ok], ef[ok])
        order = np.argsort(angles[ok])
        draw_on.show_overlay_path(angles[ok][order], ef[ok][order])
    return angles, ef, ok


class AuReferenceDialog(QDialog):
    """Measure the Fermi surface's curvature on a reference sample.

    A gold (or freshly-scraped polycrystalline) reference has a Fermi edge
    at every angle with no band structure to confuse it, so fitting it
    channel by channel gives the analyser's curvature directly -- measured,
    with an error bar per channel, instead of clicked by eye.

    The reference may be a Cut or a Map (a map is summed over the deflector,
    since every deflector position sees the same slit curvature).
    """

    def __init__(self, parent, entries, loader, default_temperature: float = 30.0):
        super().__init__(parent)
        self.setWindowTitle("Fermi surface from a reference")
        self.entries = list(entries)
        self.loader = loader

        layout = QVBoxLayout(self)
        intro = QLabel(
            "Pick the reference measurement (gold, or anything with a clean "
            "edge at every angle). Its Fermi edge is fitted channel by channel "
            "and the result becomes the correction curve.")
        intro.setWordWrap(True)
        layout.addWidget(intro)

        form = QFormLayout()
        self.dataset_combo = QComboBox()
        for label, _key in self.entries:
            self.dataset_combo.addItem(label)
        form.addRow("Reference", self.dataset_combo)

        self.temperature_box = QDoubleSpinBox()
        self.temperature_box.setRange(0.1, 2000.0)
        self.temperature_box.setDecimals(2)
        self.temperature_box.setValue(float(default_temperature))
        self.temperature_box.setToolTip(
            "Held during the fit. It cannot be separated from the resolution "
            "on a single edge, and it does not affect where the edge is.")
        form.addRow("Temperature (K)", self.temperature_box)

        self.width_box = QSpinBox()
        self.width_box.setRange(0, 200)
        self.width_box.setValue(3)
        self.width_box.setToolTip(
            "Channels added on each side before fitting. A reference too thin "
            "to fit alone becomes fittable; the channel's position does not move.")
        form.addRow("Neighbours (+/- channels)", self.width_box)

        self.step_box = QSpinBox()
        self.step_box.setRange(1, 200)
        self.step_box.setValue(1)
        self.step_box.setToolTip(
            "Fit every Nth channel. The curvature is smooth, so sampling it and "
            "letting the polynomial interpolate is usually enough -- and a "
            "thousand fits take minutes.")
        form.addRow("Fit every Nth channel", self.step_box)
        layout.addLayout(form)

        self.note = QLabel("")
        self.note.setWordWrap(True)
        layout.addWidget(self.note)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.button(QDialogButtonBox.Ok).setText("Fit the reference")
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)
        self.dataset_combo.currentIndexChanged.connect(self._describe)
        self._describe()

    def _describe(self, *_):
        data = self.selected_data()
        if data is None:
            self.note.setText("")
            return
        try:
            angles, energy, _frame, angle_label, _e_label = reference_frame(data)
        except ValueError as exc:
            self.note.setText(str(exc))
            return
        self.step_box.setValue(max(1, int(np.ceil(angles.size / 120))))
        self.temperature_box.setValue(metadata_temperature(data,
                                                            self.temperature_box.value()))
        self.note.setText(f"{angles.size} channels along {angle_label}, "
                          f"{energy.size} energy points.")

    def selected_data(self):
        index = self.dataset_combo.currentIndex()
        if not (0 <= index < len(self.entries)):
            return None
        try:
            return self.loader(self.entries[index][1])
        except Exception as exc:
            self.note.setText(f"could not be read: {exc}")
            return None

    def settings(self) -> dict:
        return {"temperature": self.temperature_box.value(),
                "half_width": self.width_box.value(),
                "step": self.step_box.value(),
                "label": self.dataset_combo.currentText()}


# --------------------------------------------------------------------------
# Bulk operations on the datasets picked in the launcher
# --------------------------------------------------------------------------
class DataOperationsDialog(QDialog):
    """Truncate, self-normalise or compress the selected datasets.

    One dialog for the three MATLAB tools (``td_demo``,
    ``self_normalization_demo``, ``data_comb_resamp_demo``) because they are
    the same action from the user's side: pick datasets, set a few numbers
    per axis, get new datasets. Each result is listed under the original's
    name with the MATLAB suffix (``_tk``, ``_s_nor``, ``_comb``), so a chain
    of operations reads back off the names.

    Every tab drives **all** the selected datasets at once. That only means
    anything if they have the same shape -- one set of bounds cannot describe
    two different axes -- so the launcher checks that before opening this and
    says so if they disagree, rather than producing a batch of quietly wrong
    results.

    Resampling, the other half of ``data_comb_resamp_demo``, is deliberately
    not here: it keeps every Nth point and discards the rest, which throws
    away most of the counts in a photon-starved scan. Compressing sums each
    block instead, so the counts are kept and the noise improves.
    """

    #: the datasets this dialog produced, for the launcher to list
    datasetsCreated = pyqtSignal(list)

    def __init__(self, datasets, parent=None):
        super().__init__(parent)
        self.datasets = list(datasets)          # [(name, data), ...]
        self.setWindowTitle("Data operations")
        self.setModal(False)

        kind = self.datasets[0][1].kind
        self.kind = kind
        self.axis_names = ARRAY_AXES.get(kind, ())
        scan = self.datasets[0][1].scan
        self.axes = [np.asarray(getattr(scan, name), dtype=float)
                     for name in self.axis_names]
        self.axis_labels = [scan.labels.get(name, name) for name in self.axis_names]

        layout = QVBoxLayout(self)
        shape = " x ".join(str(a.size) for a in self.axes)
        names = ", ".join(name for name, _ in self.datasets[:4])
        if len(self.datasets) > 4:
            names += f", ... ({len(self.datasets)} in all)"
        header = QLabel(f"<b>{len(self.datasets)} dataset(s)</b> of kind "
                        f"<b>{kind}</b>, {shape} points<br>{names}")
        header.setWordWrap(True)
        layout.addWidget(header)

        self.tabs = QTabWidget()
        self.tabs.addTab(self._truncate_tab(), "Truncate")
        self.tabs.addTab(self._normalise_tab(), "Self-normalise")
        self.tabs.addTab(self._compress_tab(), "Compress")
        layout.addWidget(self.tabs)

        self.note = QLabel("")
        self.note.setWordWrap(True)
        layout.addWidget(self.note)

        buttons = QDialogButtonBox(QDialogButtonBox.Apply | QDialogButtonBox.Close)
        buttons.button(QDialogButtonBox.Apply).clicked.connect(self.apply)
        buttons.button(QDialogButtonBox.Close).clicked.connect(self.reject)
        layout.addWidget(buttons)

    # -- tabs -------------------------------------------------------------
    def _axis_caption(self, index: int) -> str:
        axis = self.axes[index]
        return (f"{self.axis_labels[index]}  "
                f"[{axis.min():.4g} .. {axis.max():.4g}, {axis.size} pts]")

    def _truncate_tab(self) -> QWidget:
        page = QWidget()
        form = QFormLayout(page)
        self.trunc_boxes = []
        for i in range(len(self.axes)):
            row = QWidget()
            box = QHBoxLayout(row)
            box.setContentsMargins(0, 0, 0, 0)
            lo, hi = QLineEdit(), QLineEdit()
            for edit, which in ((lo, "from"), (hi, "to")):
                edit.setPlaceholderText(which)
                edit.setMaximumWidth(90)
                edit.setToolTip("Leave blank for the end of the data.")
                box.addWidget(edit)
            box.addStretch(1)
            form.addRow(self._axis_caption(i), row)
            self.trunc_boxes.append((lo, hi))

        self.trunc_by_value = QRadioButton("by value")
        self.trunc_by_index = QRadioButton("by index (1-based)")
        self.trunc_by_value.setChecked(True)
        mode = QWidget()
        mode_box = QHBoxLayout(mode)
        mode_box.setContentsMargins(0, 0, 0, 0)
        mode_box.addWidget(self.trunc_by_value)
        mode_box.addWidget(self.trunc_by_index)
        mode_box.addStretch(1)
        form.addRow("Bounds mean", mode)
        return page

    def _normalise_tab(self) -> QWidget:
        page = QWidget()
        form = QFormLayout(page)
        self.norm_checks = []
        directions = QWidget()
        box = QHBoxLayout(directions)
        box.setContentsMargins(0, 0, 0, 0)
        for i, label in enumerate(self.axis_labels):
            check = QCheckBox(label.split(" (")[0])
            check.setToolTip(
                "Tick one axis to divide every line along it by that line's own "
                "total; tick two to divide every plane by the plane's total.")
            box.addWidget(check)
            self.norm_checks.append(check)
        box.addStretch(1)
        form.addRow("Along", directions)

        window = QWidget()
        wbox = QHBoxLayout(window)
        wbox.setContentsMargins(0, 0, 0, 0)
        self.norm_from, self.norm_to = QLineEdit(), QLineEdit()
        for edit, which in ((self.norm_from, "from"), (self.norm_to, "to")):
            edit.setPlaceholderText(which)
            edit.setMaximumWidth(90)
            edit.setToolTip(
                "Optional: only these points of the first ticked axis are summed "
                "to make the normalisation. Blank means all of them.")
            wbox.addWidget(edit)
        self.norm_by_value = QRadioButton("by value")
        self.norm_by_index = QRadioButton("by index")
        self.norm_by_value.setChecked(True)
        wbox.addWidget(self.norm_by_value)
        wbox.addWidget(self.norm_by_index)
        wbox.addStretch(1)
        form.addRow("Window", window)

        self.norm_to_peak = QCheckBox("normalise to the peak instead of the sum")
        form.addRow("", self.norm_to_peak)
        return page

    def _compress_tab(self) -> QWidget:
        page = QWidget()
        form = QFormLayout(page)
        self.compress_boxes = []
        for i, axis in enumerate(self.axes):
            spin = QSpinBox()
            spin.setRange(1, max(1, int(axis.size)))
            spin.setValue(1)
            spin.setMaximumWidth(80)
            spin.setToolTip(
                "Points summed into one. The axis takes each block's mean "
                "position; a remainder too short for a whole block is dropped.")
            spin.valueChanged.connect(self._update_compress_note)
            form.addRow(self._axis_caption(i), spin)
            self.compress_boxes.append(spin)
        self.compress_note = QLabel("")
        form.addRow("Result", self.compress_note)
        self._update_compress_note()
        return page

    def _update_compress_note(self, *_):
        sizes = [axis.size // spin.value()
                 for axis, spin in zip(self.axes, self.compress_boxes)]
        before = int(np.prod([axis.size for axis in self.axes]))
        after = int(np.prod(sizes)) if all(sizes) else 0
        self.compress_note.setText(
            " x ".join(str(s) for s in sizes) +
            (f"  ({100.0 * after / before:.1f}% of the points)" if before else ""))

    # -- applying ---------------------------------------------------------
    @staticmethod
    def _number(edit):
        text = edit.text().strip()
        if not text:
            return None
        try:
            return float(text)
        except ValueError:
            raise ValueError(f"{text!r} is not a number")

    def _values_and_axes(self, data):
        scan = data.scan
        values = scan.value4d if data.kind == "spem_4d" else scan.value
        if hasattr(values, "materialise"):
            values = values.materialise()
        axes = [np.asarray(getattr(scan, name), dtype=float)
                for name in self.axis_names]
        return np.asarray(values, dtype=float), axes

    def apply(self):
        tab = self.tabs.currentIndex()
        try:
            if tab == 0:
                operation, suffix, prefix = self._apply_truncate, "_tk", "trunc"
            elif tab == 1:
                operation, suffix, prefix = self._apply_normalise, "_s_nor", "snorm"
            else:
                operation, suffix, prefix = self._apply_compress, "_comb", "compress"
            settings = operation(None, probe=True)
        except ValueError as exc:
            QMessageBox.warning(self, "Data operations", str(exc))
            return

        if not self._confirm_size():
            return

        created = []
        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            for name, data in self.datasets:
                try:
                    values, axes = self._values_and_axes(data)
                    new_values, new_axes = operation((values, axes))
                except Exception as exc:
                    QApplication.restoreOverrideCursor()
                    QMessageBox.warning(self, "Data operations",
                                        f"{name}: {exc}")
                    return
                by_name = dict(zip(self.axis_names, new_axes))
                order = CONSTRUCTOR_AXES.get(data.kind, self.axis_names)
                created.append(MemoryData(
                    data.kind, tuple(by_name[n] for n in order), new_values,
                    dict(data.scan.labels), source_label=f"{name}{suffix}",
                    parameters=settings, prefix=prefix,
                    source_path=getattr(data, "path", ""),
                    source_info=dict(data.scan.info),
                    source_motors=dict(data.scan.fourd_info)))
        finally:
            QApplication.restoreOverrideCursor()

        self.datasetsCreated.emit(created)
        shapes = " x ".join(str(a.size) for a in created[0].scan_axes()) \
            if hasattr(created[0], "scan_axes") else ""
        self.note.setText(f"Created {len(created)} dataset(s) ending in "
                          f"“{suffix}”. {shapes}")
        return created

    def _confirm_size(self) -> bool:
        """A spatial scan's cube has to be in memory to be operated on, and
        a real one is GB-sized -- so say so before spending the minutes and
        the RAM, rather than appearing to hang."""
        if self.kind != "spem_4d":
            return True
        points = int(np.prod([a.size for a in self.axes])) * len(self.datasets)
        gigabytes = points * 8 / 1e9
        if gigabytes < 2.0:
            return True
        answer = QMessageBox.question(
            self, "Data operations",
            f"This reads about {gigabytes:.1f} GB of cube into memory "
            f"({len(self.datasets)} spatial scan(s)). Continue?",
            QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        return answer == QMessageBox.Yes

    def _apply_truncate(self, payload, probe: bool = False):
        bounds = [(self._number(lo), self._number(hi)) for lo, hi in self.trunc_boxes]
        by_index = self.trunc_by_index.isChecked()
        if probe:
            if all(lo is None and hi is None for lo, hi in bounds):
                raise ValueError("Set at least one bound to truncate.")
            return {"bounds": [list(b) for b in bounds],
                    "by": "index" if by_index else "value"}
        values, axes = payload
        return truncate(values, axes, bounds, by_index=by_index)

    def _apply_normalise(self, payload, probe: bool = False):
        dims = [i for i, check in enumerate(self.norm_checks) if check.isChecked()]
        window = (self._number(self.norm_from), self._number(self.norm_to))
        by_index = self.norm_by_index.isChecked()
        to_peak = self.norm_to_peak.isChecked()
        if probe:
            if not dims:
                raise ValueError("Tick the axis (or two axes) to normalise along.")
            if len(dims) > 2:
                raise ValueError("Normalise along one axis (lines) or two (planes), "
                                 "not three -- that would divide the data by itself.")
            return {"along": [self.axis_names[d] for d in dims],
                    "window": list(window), "by": "index" if by_index else "value",
                    "to_peak": bool(to_peak)}
        values, axes = payload
        return (self_normalize(values, axes, dims, window=window,
                               by_index=by_index, to_peak=to_peak), axes)

    def _apply_compress(self, payload, probe: bool = False):
        factors = [spin.value() for spin in self.compress_boxes]
        if probe:
            if all(f == 1 for f in factors):
                raise ValueError("Set a factor above 1 on at least one axis.")
            return {"factors": list(factors)}
        values, axes = payload
        return compress(values, axes, factors)


# --------------------------------------------------------------------------
def open_viewer(data, filename: str, colormap: str, flip: bool):
    """Build the right viewer for a loaded file, or return None for kinds
    that have nothing to show."""
    if data.kind in ("spem_4d", "spem_1d"):
        return SpatialScanWindow(data, filename, colormap, flip)
    if data.kind == "cut":
        return CutWindow(data, filename, colormap, flip)
    if data.kind in CUBE_KINDS:
        # A converted k-map is displayed exactly like the map it came from --
        # same contour, same two cuts -- with the axes labelled in A^-1.
        return ContourWindow(data, filename, colormap, flip)
    return None
