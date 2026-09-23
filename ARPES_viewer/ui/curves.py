"""
ui/curves.py
============
The curve viewer -- one-dimensional data -- and the two panels that open
from it: curve fitting, and spin analysis.

What it shows
-------------
EDCs, MDCs and spin EDCs, from a file or saved out of an image viewer's
readout cursor: every channel of the dataset as a curve, with its error bars
wherever the dataset carries an uncertainty for it (:mod:`tools.curves`
explains the ``σ <name>`` convention). Channels can be hidden, offset for a
waterfall, and read off with the cursor.

What it does
------------
* **Operations** -- crop, bin, normalise, subtract a background, shift the
  axis (E − E_F), add counting errors. Each makes a new dataset in the main
  list, carrying its history, exactly as the image processing does. They
  are the curve-aware versions: an uncertainty channel is propagated, not
  averaged.
* **Curve fit...** -- peaks (Lorentzian, Gaussian, Voigt) on a background,
  optionally cut off by a Fermi edge; or the Fermi edge itself. The solver
  is the MDC/EDC fitter's (:func:`tools.peaks.fit_line`) and the edge is the
  Fermi-level dialog's (:func:`tools.fermi.fit_fermi_edge`), so a peak
  fitted here and one fitted on a cut are the same number.
* **Spin analysis...** (spin EDCs) -- polarisation along an axis from the
  channels, with the instrumental asymmetry measured and, where the
  channels allow, cancelled; the spin-resolved spectra; all with counting
  errors. :mod:`tools.spin` has the model.

Nothing here computes anything the ``tools`` modules do not; this file is
drawing and wiring.
"""
from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtWidgets import (QAbstractItemView, QApplication, QCheckBox,
                             QComboBox, QDialog, QDialogButtonBox,
                             QDoubleSpinBox, QFormLayout, QGroupBox,
                             QHBoxLayout, QHeaderView, QLabel, QListWidget,
                             QListWidgetItem, QMainWindow, QMessageBox,
                             QPushButton, QSpinBox, QSplitter,
                             QStackedWidget, QTableWidget, QTableWidgetItem,
                             QVBoxLayout, QWidget)

from tools import curves as C
from tools import process as P
from tools import spin as S
from ui.widgets import MemoryData, _is_energy_label

__all__ = ["CurveWindow", "CurveFitDialog", "SpinAnalysisDialog",
           "curve_dataset", "curve_panel", "CHANNEL_COLOURS"]

#: Distinguishable in print and for the common colour-vision deficiencies;
#: the first two are the image viewers' EDC and MDC colours.
CHANNEL_COLOURS = ("#1f5aa6", "#c23b22", "#2a8a3e", "#8a4fb0", "#d08a00",
                   "#2aa1a8", "#6b6b6b", "#b0306e")


def _colour(index: int) -> str:
    return CHANNEL_COLOURS[index % len(CHANNEL_COLOURS)]


def _legend(text) -> str:
    """A legend entry. pyqtgraph renders legends as HTML, so a spin
    channel's "<0,0>" would be swallowed as a tag without this."""
    import html
    return html.escape(str(text))


def _spin(value=0.0, lo=-1e9, hi=1e9, decimals=4, step=0.01, width=90):
    box = QDoubleSpinBox()
    box.setRange(lo, hi)
    box.setDecimals(decimals)
    box.setSingleStep(step)
    box.setValue(value)
    box.setKeyboardTracking(False)
    box.setMaximumWidth(width)
    return box


def _unique(name: str, taken) -> str:
    taken = set(taken or ())
    if name not in taken:
        return name
    n = 2
    while f"{name} ({n})" in taken:
        n += 1
    return f"{name} ({n})"


def curve_dataset(kind, x, values, names, *, x_label, value_label, name,
                  step, parameters=None, source_info=None, source_path="",
                  prefix="curve", source=""):
    """A curve table as a :class:`MemoryData`, with its channels named and
    the operation that made it recorded."""
    values = np.asarray(values, dtype=float)
    if values.ndim == 1:
        values = values[:, None]
    info = P.record_step(dict(source_info or {}),
                         P.Step(step, dict(parameters or {}),
                                source=source or name))
    info.update(C.curve_info(names, value_label))
    return MemoryData(kind, (np.asarray(x, dtype=float),
                             np.arange(values.shape[1], dtype=float)),
                      values, {"x": x_label, "y": "Channel"},
                      source_label=name, parameters=parameters, prefix=prefix,
                      source_path=source_path, source_info=info)


def curve_panel(x, columns, *, x_label, y_label, name, errors=None):
    """A figure panel of curves: ``columns`` is ``[(label, y, colour)]``;
    ``errors`` an optional ``{label: σ}`` drawn as dotted ±σ lines."""
    import tools.figure as F

    x = np.asarray(x, dtype=float)
    overlays, ys = [], []
    for label, y, colour in columns:
        y = np.asarray(y, dtype=float)
        ys.append(y)
        overlays.append(F.Overlay(xs=x.tolist(), ys=y.tolist(), color=colour,
                                  width_pt=1.2, name=str(label)))
        sigma = (errors or {}).get(label)
        if sigma is not None:
            for sign in (-1, 1):
                edge = y + sign * np.asarray(sigma, dtype=float)
                ys.append(edge)
                overlays.append(F.Overlay(xs=x.tolist(), ys=edge.tolist(),
                                          color=colour, width_pt=0.5,
                                          style="dotted",
                                          name=f"{label} {'+' if sign > 0 else '−'}σ"))
    stacked = np.concatenate([y[np.isfinite(y)] for y in ys]) if ys else np.zeros(1)
    lo = float(np.min(stacked)) if stacked.size else 0.0
    hi = float(np.max(stacked)) if stacked.size else 1.0
    if hi <= lo:
        hi = lo + 1.0
    panel = F.Panel(data=F.PanelData(
        array=np.full((2, 2), np.nan),
        x=np.array([float(np.nanmin(x)), float(np.nanmax(x))]),
        y=np.array([lo, hi]), x_label=x_label, y_label=y_label, name=name,
        kind="stack"))
    panel.overlays = overlays
    return panel


# ==========================================================================
# The viewer
# ==========================================================================
OPERATIONS = (
    ("crop", "Crop to the range"),
    ("bin", "Bin"),
    ("normalise", "Normalise"),
    ("background", "Subtract a background"),
    ("shift", "Shift the axis"),
    ("poisson", "Add counting errors (√N)"),
)


class CurveWindow(QMainWindow):
    """One-dimensional data: every channel as a curve."""

    closed = pyqtSignal(object)
    datasetCreated = pyqtSignal(object)

    def __init__(self, data, filename: str, colormap: str = "gray",
                 flip: bool = False):
        super().__init__()
        self.data = data
        self.filename = filename
        # Read back by the launcher when it saves its settings; a curve
        # window has no colormap of its own but keeps what it was given.
        self.colormap, self.flip = colormap, flip
        self.kind = data.kind
        scan = data.scan
        self.x = np.asarray(scan.x, dtype=float).reshape(-1)
        values = np.asarray(scan.value, dtype=float)
        self.values = values[:, None] if values.ndim == 1 else values
        if self.values.shape[0] != self.x.size and self.values.shape[1] == self.x.size:
            self.values = self.values.T
        self.info = dict(scan.info or {})
        self.names = C.channel_names(self.info, self.values.shape[1], self.kind)
        self.x_label = (scan.labels or {}).get("x") or (
            "Energy (eV)" if self.kind != "mdc" else "Angle (deg)")
        self.value_label = C.value_label(self.info)
        self.setWindowTitle(f"{filename}  [{self.kind}]")
        self.resize(1000, 680)
        self._dialogs = []
        self._items = []
        self._build()
        self.redraw()

    # -- construction --------------------------------------------------------
    def _build(self):
        central = QWidget()
        self.setCentralWidget(central)
        root = QVBoxLayout(central)
        root.setContentsMargins(4, 4, 4, 4)

        bar = QHBoxLayout()
        self.fit_button = QPushButton("Curve fit...")
        self.fit_button.setToolTip(
            "Peaks on a background (optionally cut off by a Fermi edge), or "
            "the Fermi edge itself, with uncertainties on every parameter.")
        self.fit_button.clicked.connect(self.open_fit)
        bar.addWidget(self.fit_button)
        self.spin_button = QPushButton("Spin analysis...")
        self.spin_button.setToolTip(
            "Polarisation from the spin channels, the instrumental "
            "asymmetry, and the spin-resolved spectra.")
        self.spin_button.clicked.connect(self.open_spin_analysis)
        self.spin_button.setVisible(self.kind == "spin_edc")
        bar.addWidget(self.spin_button)
        self.figure_button = QPushButton("As a figure...")
        self.figure_button.clicked.connect(self.to_figure)
        bar.addWidget(self.figure_button)
        bar.addSpacing(12)
        self.errors_box = QCheckBox("Error bars")
        self.errors_box.setChecked(True)
        self.errors_box.toggled.connect(self.redraw)
        bar.addWidget(self.errors_box)
        bar.addWidget(QLabel("Offset"))
        self.offset_box = _spin(0.0, 0.0, 1e9, 4, 0.1, 80)
        self.offset_box.setToolTip(
            "Shift each shown channel up by this much more than the one "
            "before, for a waterfall.")
        self.offset_box.valueChanged.connect(self.redraw)
        bar.addWidget(self.offset_box)
        self.range_box = QCheckBox("Range")
        self.range_box.setToolTip(
            "A draggable x range: what Crop keeps, what a region "
            "normalisation or background is taken from, and the fit window.")
        self.range_box.toggled.connect(self._toggle_range)
        bar.addWidget(self.range_box)
        self.range_lo = _spin(0.0, width=80)
        self.range_hi = _spin(0.0, width=80)
        for box in (self.range_lo, self.range_hi):
            box.editingFinished.connect(self._range_typed)
            box.setEnabled(False)
            bar.addWidget(box)
        bar.addStretch(1)
        self.readout = QLabel("")
        bar.addWidget(self.readout)
        root.addLayout(bar)

        splitter = QSplitter(Qt.Horizontal)
        self.channel_list = QListWidget()
        self.channel_list.setMaximumWidth(260)
        for index in C.data_channels(self.names):
            item = QListWidgetItem(self.names[index])
            item.setData(Qt.UserRole, index)
            item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
            item.setCheckState(Qt.Checked)
            item.setForeground(pg.mkColor(_colour(index)))
            self.channel_list.addItem(item)
        self.channel_list.itemChanged.connect(lambda *_: self.redraw())
        splitter.addWidget(self.channel_list)

        self.plot = pg.PlotWidget()
        self.plot.setBackground("w")
        self.plot.showGrid(x=True, y=True, alpha=0.25)
        self.plot.setLabel("bottom", self.x_label)
        self.plot.setLabel("left", self.value_label)
        self.legend = self.plot.addLegend(offset=(-10, 10))
        self.plot.scene().sigMouseMoved.connect(self._mouse_moved)
        splitter.addWidget(self.plot)
        splitter.setStretchFactor(1, 1)
        root.addWidget(splitter, stretch=1)

        lo, hi = float(np.nanmin(self.x)), float(np.nanmax(self.x))
        span = hi - lo
        self.region = pg.LinearRegionItem(values=(lo + 0.25 * span,
                                                  lo + 0.75 * span))
        self.region.setZValue(-10)
        self.region.sigRegionChanged.connect(self._range_dragged)
        self._range_dragged()

        root.addWidget(self._operations())
        self.statusBar()

    def _operations(self) -> QWidget:
        box = QGroupBox("Operations  (each makes a new dataset in the list)")
        row = QHBoxLayout(box)
        self.operation = QComboBox()
        for key, text in OPERATIONS:
            self.operation.addItem(text, key)
        row.addWidget(self.operation)
        self.op_pages = QStackedWidget()
        # crop
        self.op_pages.addWidget(QLabel("keeps what is inside the Range"))
        # bin
        page = QWidget(); lay = QHBoxLayout(page); lay.setContentsMargins(0, 0, 0, 0)
        self.bin_factor = QSpinBox(); self.bin_factor.setRange(2, 1000)
        self.bin_factor.setValue(2)
        self.bin_how = QComboBox()
        self.bin_how.addItem("sum (counts stay counts)", "sum")
        self.bin_how.addItem("mean", "mean")
        lay.addWidget(QLabel("factor")); lay.addWidget(self.bin_factor)
        lay.addWidget(self.bin_how); lay.addStretch(1)
        self.op_pages.addWidget(page)
        # normalise
        page = QWidget(); lay = QHBoxLayout(page); lay.setContentsMargins(0, 0, 0, 0)
        self.norm_how = QComboBox()
        self.norm_how.addItem("to the maximum", "max")
        self.norm_how.addItem("to the area", "area")
        self.norm_how.addItem("to the mean over the Range", "region")
        self.norm_together = QCheckBox("one factor for all channels")
        self.norm_together.setChecked(True)
        self.norm_together.setToolTip(
            "Keeps the channels' relative sizes -- what you want for a "
            "spin-up and a spin-down spectrum.")
        lay.addWidget(self.norm_how); lay.addWidget(self.norm_together)
        lay.addStretch(1)
        self.op_pages.addWidget(page)
        # background
        page = QWidget(); lay = QHBoxLayout(page); lay.setContentsMargins(0, 0, 0, 0)
        self.bg_how = QComboBox()
        self.bg_how.addItem("constant: the mean over the Range", "constant")
        self.bg_how.addItem("linear: a line through the Range", "linear")
        self.bg_how.addItem("Shirley", "shirley")
        lay.addWidget(self.bg_how); lay.addStretch(1)
        self.op_pages.addWidget(page)
        # shift
        page = QWidget(); lay = QHBoxLayout(page); lay.setContentsMargins(0, 0, 0, 0)
        self.shift_value = _spin(0.0, decimals=5, step=0.001, width=100)
        self.shift_label = QCheckBox("relabel as E − E_F")
        self.shift_label.setChecked(_is_energy_label(self.x_label))
        self.shift_label.setEnabled(_is_energy_label(self.x_label))
        lay.addWidget(QLabel("subtract")); lay.addWidget(self.shift_value)
        lay.addWidget(self.shift_label); lay.addStretch(1)
        self.op_pages.addWidget(page)
        # poisson
        self.op_pages.addWidget(QLabel("σ = √N for every channel that is raw "
                                       "counts and has no σ yet"))
        self.operation.currentIndexChanged.connect(self.op_pages.setCurrentIndex)
        row.addWidget(self.op_pages, stretch=1)
        apply = QPushButton("Apply")
        apply.clicked.connect(lambda: self.apply_operation())
        row.addWidget(apply)
        return box

    # -- the range -----------------------------------------------------------
    def _toggle_range(self, on: bool):
        if on:
            self.plot.addItem(self.region)
        else:
            self.plot.removeItem(self.region)
        for box in (self.range_lo, self.range_hi):
            box.setEnabled(on)

    def _range_dragged(self, *_):
        lo, hi = self.region.getRegion()
        for box, value in ((self.range_lo, lo), (self.range_hi, hi)):
            box.blockSignals(True)
            box.setValue(float(value))
            box.blockSignals(False)

    def _range_typed(self):
        self.region.setRegion((self.range_lo.value(), self.range_hi.value()))

    def x_range(self):
        """The Range, or None when it is switched off."""
        if not self.range_box.isChecked():
            return None
        lo, hi = self.region.getRegion()
        return (float(min(lo, hi)), float(max(lo, hi)))

    def set_x_range(self, lo: float, hi: float):
        self.range_box.setChecked(True)
        self.region.setRegion((float(lo), float(hi)))

    # -- drawing -------------------------------------------------------------
    def shown_channels(self):
        out = []
        for row in range(self.channel_list.count()):
            item = self.channel_list.item(row)
            if item.checkState() == Qt.Checked:
                out.append(int(item.data(Qt.UserRole)))
        return out

    def redraw(self, *_):
        for item in self._items:
            self.plot.removeItem(item)
        self._items = []
        self.legend.clear()
        offset = self.offset_box.value()
        for order, index in enumerate(self.shown_channels()):
            y = self.values[:, index] + order * offset
            colour = _colour(index)
            curve = self.plot.plot(self.x, y, pen=pg.mkPen(colour, width=1.6),
                                   symbol="o", symbolSize=4,
                                   symbolBrush=colour, symbolPen=None,
                                   name=_legend(self.names[index]))
            self._items.append(curve)
            sigma_index = C.sigma_of(self.names, index)
            if sigma_index is not None and self.errors_box.isChecked():
                sigma = self.values[:, sigma_index]
                good = np.isfinite(sigma) & np.isfinite(y)
                bars = pg.ErrorBarItem(x=self.x[good], y=y[good],
                                       height=2 * sigma[good],
                                       pen=pg.mkPen(colour, width=1))
                self.plot.addItem(bars)
                self._items.append(bars)

    def _mouse_moved(self, position):
        box = self.plot.getPlotItem().vb
        if not box.sceneBoundingRect().contains(position):
            return
        point = box.mapSceneToView(position)
        self.readout.setText(f"x = {point.x():.5g}   y = {point.y():.5g}")

    # -- operations ----------------------------------------------------------
    def existing_names(self):
        receiver = getattr(self, "_name_source", None)
        return list(receiver()) if callable(receiver) else []

    def _emit(self, x, values, names, *, suffix, step, parameters,
              x_label=None, value_label=None, kind=None):
        name = _unique(f"{self.filename} [{suffix}]", self.existing_names())
        data = curve_dataset(kind or self.kind, x, values, names,
                             x_label=x_label or self.x_label,
                             value_label=value_label or self.value_label,
                             name=name, step=step, parameters=parameters,
                             source_info=self.info, source=self.filename,
                             source_path=getattr(self.data, "path", ""))
        self.datasetCreated.emit(data)
        self.statusBar().showMessage(f"Added “{name}” to the list")
        return data

    def apply_operation(self, key: str = None):
        """Run the selected (or named) operation and list the result.
        Returns the new dataset, or None if the operation was refused."""
        key = key or self.operation.currentData()
        x, values, names = self.x, self.values, list(self.names)
        region = self.x_range()
        try:
            if key == "crop":
                if region is None:
                    raise ValueError("switch the Range on and drag it first")
                x2, v2, n2 = C.crop(x, values, names, *region)
                return self._emit(x2, v2, n2, suffix="crop", step="crop",
                                  parameters={"range": list(region)})
            if key == "bin":
                factor, how = self.bin_factor.value(), self.bin_how.currentData()
                x2, v2, n2 = C.rebin(x, values, names, factor, how)
                return self._emit(x2, v2, n2, suffix=f"bin{factor}",
                                  step="bin",
                                  parameters={"factor": factor, "how": how})
            if key == "normalise":
                how = self.norm_how.currentData()
                if how == "region" and region is None:
                    raise ValueError("switch the Range on and drag it first")
                together = self.norm_together.isChecked()
                x2, v2, n2, factors = C.normalise(x, values, names, how,
                                                  region, together)
                return self._emit(
                    x2, v2, n2, suffix="norm", step="normalise",
                    parameters={"how": how, "together": together,
                                "region": list(region) if region else "",
                                "factors": ", ".join(f"{names[i]}={f:.6g}"
                                                     for i, f in factors.items())},
                    value_label=f"{self.value_label} (normalised)")
            if key == "background":
                how = self.bg_how.currentData()
                if how != "shirley" and region is None:
                    raise ValueError("switch the Range on and drag it over "
                                     "the background first")
                x2, v2, n2, _ = C.subtract_background(x, values, names, how,
                                                      region)
                return self._emit(x2, v2, n2, suffix="bg", step="background",
                                  parameters={"how": how,
                                              "region": list(region) if region else ""})
            if key == "shift":
                delta = self.shift_value.value()
                label = self.x_label
                if self.shift_label.isChecked() and self.shift_label.isEnabled():
                    label = "E − E_F (eV)"
                return self._emit(C.shift_x(x, delta), values, names,
                                  suffix="shifted", step="shift_x",
                                  parameters={"subtracted": delta},
                                  x_label=label)
            if key == "poisson":
                v2, n2, added = C.with_poisson_sigma(values, names)
                if not added:
                    raise ValueError("every channel either has a σ already or "
                                     "is not raw counts")
                return self._emit(x, v2, n2, suffix="σ", step="poisson_sigma",
                                  parameters={"added": ", ".join(added)})
            raise ValueError(f"unknown operation {key!r}")
        except ValueError as exc:
            QMessageBox.information(self, "Operation", str(exc))
            return None

    # -- the panels ------------------------------------------------------------
    def _keep(self, dialog):
        self._dialogs.append(dialog)
        dialog.finished.connect(lambda *_: self._dialogs.remove(dialog)
                                if dialog in self._dialogs else None)
        if hasattr(dialog, "datasetsCreated"):
            dialog.datasetsCreated.connect(
                lambda made: [self.datasetCreated.emit(one) for one in made])
        dialog.show()
        return dialog

    def open_fit(self):
        return self._keep(CurveFitDialog(self))

    def open_spin_analysis(self):
        if self.kind != "spin_edc":
            QMessageBox.information(self, "Spin analysis",
                                    "Only a spin EDC has spin channels.")
            return None
        try:
            dialog = SpinAnalysisDialog(self)
        except ValueError as exc:
            QMessageBox.information(self, "Spin analysis", str(exc))
            return None
        return self._keep(dialog)

    def to_figure(self):
        from ui.figure import FigureWindow

        offset = self.offset_box.value()
        columns, errors = [], {}
        for order, index in enumerate(self.shown_channels()):
            y = self.values[:, index] + order * offset
            columns.append((self.names[index], y, _colour(index)))
            s = C.sigma_of(self.names, index)
            if s is not None and self.errors_box.isChecked():
                errors[self.names[index]] = self.values[:, s]
        if not columns:
            return None
        panel = curve_panel(self.x, columns, x_label=self.x_label,
                            y_label=self.value_label, name=self.filename,
                            errors=errors)
        window = FigureWindow([panel], parent=self)
        window.show()
        self._dialogs.append(window)
        return window

    def closeEvent(self, event):
        # Everything is copied into memory at construction, so the file
        # reference the launcher handed over goes back now -- once, however
        # many close events arrive (see ViewerWindow.closeEvent).
        if getattr(self, "_released", False):
            super().closeEvent(event)
            return
        self._released = True
        for dialog in list(self._dialogs):
            dialog.close()
        try:
            self.data.close()
        except Exception:
            pass
        self.closed.emit(self)
        super().closeEvent(event)


# ==========================================================================
# Curve fitting
# ==========================================================================
class CurveFitDialog(QDialog):
    """Fit one channel of the viewer's curve."""

    datasetsCreated = pyqtSignal(list)

    def __init__(self, viewer: CurveWindow):
        super().__init__(viewer)
        from tools.peaks import BACKGROUNDS, PEAK_SHAPES

        self.viewer = viewer
        self.setWindowTitle(f"Curve fit — {viewer.filename}")
        self.setModal(False)
        self.resize(980, 820)
        self.last = None
        self.edge_found = False
        self._picking = False
        layout = QVBoxLayout(self)

        top = QFormLayout()
        self.channel = QComboBox()
        for index in C.data_channels(viewer.names):
            self.channel.addItem(viewer.names[index], index)
        self.channel.currentIndexChanged.connect(self._channel_changed)
        top.addRow("Channel", self.channel)
        self.model = QComboBox()
        self.model.addItem("Peaks on a background", "peaks")
        self.model.addItem("Fermi edge", "fermi")
        self.model.currentIndexChanged.connect(self._model_changed)
        top.addRow("Model", self.model)
        span = QHBoxLayout()
        lo, hi = float(np.nanmin(viewer.x)), float(np.nanmax(viewer.x))
        self.lo, self.hi = _spin(lo, decimals=5), _spin(hi, decimals=5)
        span.addWidget(self.lo); span.addWidget(QLabel("to")); span.addWidget(self.hi)
        from_viewer = QPushButton("From the viewer's Range")
        from_viewer.clicked.connect(self._range_from_viewer)
        span.addWidget(from_viewer); span.addStretch(1)
        top.addRow("Fit window", span)
        self.weighting = QComboBox()
        top.addRow("Weighting", self.weighting)
        layout.addLayout(top)

        # -- peaks -----------------------------------------------------------
        self.peaks_box = QGroupBox("Peaks")
        peaks = QVBoxLayout(self.peaks_box)
        self.table = QTableWidget(0, 4)
        self.table.setHorizontalHeaderLabels(["Shape", "Centre", "FWHM",
                                              "Hold FWHM"])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.table.setMaximumHeight(150)
        peaks.addWidget(self.table)
        buttons = QHBoxLayout()
        self.pick_button = QPushButton("Add by clicking the plot")
        self.pick_button.setCheckable(True)
        self.pick_button.toggled.connect(self._toggle_picking)
        find = QPushButton("Find peaks")
        find.clicked.connect(self.find_peaks)
        remove = QPushButton("Remove selected")
        remove.clicked.connect(self._remove_rows)
        for widget in (self.pick_button, find, remove):
            buttons.addWidget(widget)
        buttons.addStretch(1)
        peaks.addLayout(buttons)
        options = QHBoxLayout()
        self.background = QComboBox()
        self.background.addItems(BACKGROUNDS)
        self.background.setCurrentText("linear")
        options.addWidget(QLabel("Background")); options.addWidget(self.background)
        self.resolution = _spin(0.0, 0.0, 1e3, 5, 0.001, 90)
        self.resolution.setToolTip("Instrumental Gaussian FWHM, in x units; "
                                   "0 = none. Peaks are convolved with it.")
        options.addWidget(QLabel("Resolution")); options.addWidget(self.resolution)
        self.fermi_cut = QCheckBox("× Fermi edge at")
        self.fermi_cut.setToolTip(
            "Multiply the peaks by the Fermi–Dirac occupation before the "
            "resolution smears them -- for EDC peaks near E_F.")
        self.cut_ef = _spin(0.0, decimals=5, step=0.001)
        self.cut_t = _spin(30.0, 0.1, 2000, 1, 1, 70)
        options.addWidget(self.fermi_cut); options.addWidget(self.cut_ef)
        options.addWidget(QLabel("T (K)")); options.addWidget(self.cut_t)
        options.addStretch(1)
        peaks.addLayout(options)
        self.PEAK_SHAPES = PEAK_SHAPES
        layout.addWidget(self.peaks_box)

        # -- Fermi edge -------------------------------------------------------
        self.fermi_box = QGroupBox("Fermi edge")
        fermi = QHBoxLayout(self.fermi_box)
        self.edge_t = _spin(30.0, 0.1, 2000, 1, 1, 70)
        self.hold_t = QCheckBox("hold T")
        self.hold_t.setChecked(True)
        self.hold_t.setToolTip("On one edge T and the resolution trade off; "
                               "hold the one you know.")
        fermi.addWidget(QLabel("T (K)")); fermi.addWidget(self.edge_t)
        fermi.addWidget(self.hold_t)
        self.shift_button = QPushButton("List a copy with E_F = 0")
        self.shift_button.setEnabled(False)
        self.shift_button.clicked.connect(self.shift_by_ef)
        fermi.addWidget(self.shift_button)
        fermi.addStretch(1)
        layout.addWidget(self.fermi_box)

        # -- plots ------------------------------------------------------------
        self.plot = pg.PlotWidget()
        self.plot.setBackground("w")
        self.plot.showGrid(x=True, y=True, alpha=0.25)
        self.plot.addLegend(offset=(-10, 10))
        self.plot.setLabel("left", viewer.value_label)
        self.plot.scene().sigMouseClicked.connect(self._clicked)
        self.residual_plot = pg.PlotWidget()
        self.residual_plot.setBackground("w")
        self.residual_plot.setMaximumHeight(140)
        self.residual_plot.showGrid(x=True, y=True, alpha=0.25)
        self.residual_plot.setLabel("left", "residual / σ")
        self.residual_plot.setLabel("bottom", viewer.x_label)
        self.residual_plot.setXLink(self.plot)
        layout.addWidget(self.plot, stretch=3)
        layout.addWidget(self.residual_plot, stretch=1)

        self.results = QTableWidget(0, 3)
        self.results.setHorizontalHeaderLabels(["Parameter", "Value", "± σ"])
        self.results.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.results.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.results.setMaximumHeight(170)
        layout.addWidget(self.results)
        self.report = QLabel("")
        self.report.setWordWrap(True)
        self.report.setTextInteractionFlags(Qt.TextSelectableByMouse)
        layout.addWidget(self.report)

        box = QDialogButtonBox()
        fit = box.addButton("Fit", QDialogButtonBox.ActionRole)
        fit.clicked.connect(lambda: self.run_fit())
        self.to_list = box.addButton("Fit to list", QDialogButtonBox.ActionRole)
        self.to_list.clicked.connect(self.export)
        self.to_list.setEnabled(False)
        copy = box.addButton("Copy results", QDialogButtonBox.ActionRole)
        copy.clicked.connect(self.copy_results)
        box.addButton("Close", QDialogButtonBox.RejectRole).clicked.connect(
            self.reject)
        layout.addWidget(box)

        if viewer.x_range() is not None:
            self._range_from_viewer()
        self._channel_changed()
        self._model_changed()

    # -- inputs ---------------------------------------------------------------
    def _channel_changed(self, *_):
        index = self.channel.currentData()
        self.weighting.clear()
        if C.sigma_of(self.viewer.names, index) is not None:
            self.weighting.addItem("the channel's own σ", "sigma")
        y = self.viewer.values[:, index]
        if C.looks_like_counts(y):
            self.weighting.addItem("Poisson (√N)", "poisson")
        self.weighting.addItem("uniform", "none")
        self.draw_data()

    def _model_changed(self, *_):
        peaks = self.model.currentData() == "peaks"
        self.peaks_box.setVisible(peaks)
        self.fermi_box.setVisible(not peaks)

    def _range_from_viewer(self):
        region = self.viewer.x_range()
        if region is None:
            QMessageBox.information(self, "Fit window",
                                    "Switch the viewer's Range on first.")
            return
        self.lo.setValue(region[0])
        self.hi.setValue(region[1])
        self.draw_data()

    def fit_data(self):
        """``(x, y, σ or None)`` inside the fit window."""
        lo, hi = sorted((self.lo.value(), self.hi.value()))
        x = self.viewer.x
        keep = (x >= lo) & (x <= hi)
        index = self.channel.currentData()
        y = self.viewer.values[:, index]
        keep &= np.isfinite(y)
        sigma = None
        s = C.sigma_of(self.viewer.names, index)
        if s is not None:
            sigma = self.viewer.values[:, s][keep]
        return x[keep], y[keep], sigma

    # -- peaks table ------------------------------------------------------------
    def add_peak(self, centre: float, width: float = None, shape="lorentzian",
                 hold=False):
        from tools.peaks import suggest_seed

        x, y, _ = self.fit_data()
        if width is None and x.size > 3:
            _, width = suggest_seed(x, y, centre)
        row = self.table.rowCount()
        self.table.insertRow(row)
        combo = QComboBox()
        combo.addItems(self.PEAK_SHAPES)
        combo.setCurrentText(shape)
        self.table.setCellWidget(row, 0, combo)
        self.table.setItem(row, 1, QTableWidgetItem(f"{centre:.6g}"))
        self.table.setItem(row, 2, QTableWidgetItem(f"{(width or 0.1):.6g}"))
        check = QTableWidgetItem("")
        check.setFlags(Qt.ItemIsUserCheckable | Qt.ItemIsEnabled)
        check.setCheckState(Qt.Checked if hold else Qt.Unchecked)
        self.table.setItem(row, 3, check)

    def peaks(self):
        out = []
        for row in range(self.table.rowCount()):
            try:
                centre = float(self.table.item(row, 1).text())
                width = abs(float(self.table.item(row, 2).text()))
            except (AttributeError, ValueError):
                raise ValueError(f"peak {row + 1}: centre and FWHM must be numbers")
            out.append({"shape": self.table.cellWidget(row, 0).currentText(),
                        "centre": centre, "width": width,
                        "hold": self.table.item(row, 3).checkState() == Qt.Checked})
        return out

    def _remove_rows(self):
        for row in sorted({i.row() for i in self.table.selectedIndexes()},
                          reverse=True):
            self.table.removeRow(row)

    def _toggle_picking(self, on: bool):
        self._picking = on

    def _clicked(self, event):
        if not self._picking:
            return
        box = self.plot.getPlotItem().vb
        if not box.sceneBoundingRect().contains(event.scenePos()):
            return
        self.add_peak(float(box.mapSceneToView(event.scenePos()).x()))

    def find_peaks(self, max_peaks: int = 6):
        """Seed the table from the local maxima that stand out of the
        noise: prominence above 5% of the window's range."""
        from scipy.signal import find_peaks

        x, y, _ = self.fit_data()
        if x.size < 5:
            return 0
        kernel = np.ones(3) / 3.0
        smooth = np.convolve(y, kernel, mode="same")
        found, props = find_peaks(smooth, prominence=0.05 * float(np.ptp(smooth)))
        order = np.argsort(props["prominences"])[::-1][:max_peaks]
        for index in sorted(found[order]):
            self.add_peak(float(x[index]))
        return len(order)

    # -- fitting -----------------------------------------------------------------
    def run_fit(self):
        """Fit with the current settings. Returns the fit, or None (and says
        why) if it could not be done."""
        from tools.fermi import fit_fermi_edge
        from tools.peaks import Band, FitSettings, fit_line

        x, y, sigma = self.fit_data()
        weighting = self.weighting.currentData()
        try:
            if self.model.currentData() == "fermi":
                from tools.fermi import initial_guess, steepest_drop
                from tools.kzmap import is_an_edge

                fixed = ("temperature",) if self.hold_t.isChecked() else ()
                start = initial_guess(x, y, temperature=self.edge_t.value())
                start["ef"] = steepest_drop(x, y)
                fit = fit_fermi_edge(
                    x, y, start=start, temperature=self.edge_t.value(),
                    fixed=fixed,
                    weighting="uniform" if weighting == "none" else "poisson")
                step = float(np.median(np.diff(np.sort(x)))) if x.size > 1 else 0.0
                self.edge_found = is_an_edge(fit, float(x.min()),
                                             float(x.max()), 2.0 * step)
                self.last = ("fermi", fit, x, y, sigma)
            else:
                peaks = self.peaks()
                if not peaks:
                    raise ValueError("add at least one peak (click the plot, "
                                     "or Find peaks)")
                settings = FitSettings(
                    background=self.background.currentText(),
                    resolution=self.resolution.value(),
                    weighting="none" if weighting == "none" else "poisson",
                    fermi=self.fermi_cut.isChecked(), ef=self.cut_ef.value(),
                    temperature=self.cut_t.value())
                bands, guesses = [], {}
                for i, peak in enumerate(peaks):
                    name = f"P{i + 1}"
                    bands.append(Band(name, shape=peak["shape"],
                                      fix_width=peak["hold"]))
                    height = float(np.interp(peak["centre"], x, y))
                    guesses[name] = (peak["centre"], peak["width"], height)
                fit = fit_line(x, y, bands, settings, guesses=guesses,
                               sigma=sigma if weighting == "sigma" else None)
                self.last = ("peaks", fit, x, y, sigma)
                # Start the next fit from this one.
                for row in range(self.table.rowCount()):
                    self.table.item(row, 1).setText(f"{fit.centres[row]:.6g}")
                    self.table.item(row, 2).setText(f"{fit.widths[row]:.6g}")
        except (ValueError, RuntimeError, np.linalg.LinAlgError) as exc:
            self.last = None
            self.report.setText(f"<span style='color:#b00'>{exc}</span>")
            self.to_list.setEnabled(False)
            return None
        self.to_list.setEnabled(True)
        found = self.last[0] == "fermi" and getattr(self, "edge_found", False)
        self.shift_button.setEnabled(found)
        self._show_result()
        if self.last[0] == "fermi" and not found:
            self.report.setText(
                "<span style='color:#b00'><b>No edge found.</b> The fit "
                "converged, but E_F sits at the end of the window or the step "
                "is not significant -- it has fitted something else (a band, "
                "the secondary tail). Narrow the window around the edge and "
                "fit again.</span><br>" + self.report.text())
        return self.last[1]

    def rows(self):
        """``[(name, value, error)]`` for the last fit."""
        if self.last is None:
            return []
        kind, fit = self.last[0], self.last[1]
        if kind == "fermi":
            from tools.fermi import PARAMETERS

            return [(name, fit.values[name], fit.errors.get(name, float("nan")))
                    for name in PARAMETERS]
        rows = []
        for i, name in enumerate(fit.names):
            rows += [(f"{name} centre", fit.centres[i], fit.centre_errors[i]),
                     (f"{name} FWHM", fit.widths[i], fit.width_errors[i]),
                     (f"{name} height", fit.heights[i], fit.height_errors[i]),
                     (f"{name} area", fit.areas[i], float("nan"))]
        return rows

    def _curves(self):
        """``(model, components)`` of the last fit on its own x."""
        kind, fit, x = self.last[0], self.last[1], self.last[2]
        if kind == "fermi":
            return fit.model(x), []
        components = [(name, fit.component(i)) for i, name in enumerate(fit.names)]
        return fit.model, components

    def _show_result(self):
        kind, fit, x, y, sigma = self.last
        model, components = self._curves()
        self.draw_data()
        order = np.argsort(x) if kind == "fermi" else slice(None)
        xs = x[order] if kind == "fermi" else fit.x
        self.plot.plot(xs, np.asarray(model)[order] if kind == "fermi" else model,
                       pen=pg.mkPen("#c23b22", width=2), name=_legend("fit"))
        for i, (name, comp) in enumerate(components):
            self.plot.plot(fit.x, comp, pen=pg.mkPen(_colour(i + 2), width=1,
                                                      style=Qt.DashLine),
                           name=_legend(name))
        if kind == "peaks":
            self.plot.plot(fit.x, fit.background, pen=pg.mkPen("#888", width=1,
                                                               style=Qt.DotLine),
                           name=_legend("background"))
        data_y = y if kind == "fermi" else fit.y
        data_x = x if kind == "fermi" else fit.x
        model_y = np.asarray(model)
        if sigma is not None and self.weighting.currentData() == "sigma":
            scale = np.interp(data_x, x, sigma)
        elif self.weighting.currentData() == "none":
            scale = np.full_like(data_y, float(np.std(data_y - model_y)) or 1.0)
        else:
            scale = np.sqrt(np.maximum(np.abs(data_y), 1.0))
        self.residual_plot.clear()
        self.residual_plot.plot(data_x, (data_y - model_y) / scale,
                                pen=None, symbol="o", symbolSize=4,
                                symbolBrush="#444")
        self.residual_plot.addItem(pg.InfiniteLine(0, angle=0,
                                                   pen=pg.mkPen("#c23b22")))

        rows = self.rows()
        self.results.setRowCount(len(rows))
        for r, (name, value, error) in enumerate(rows):
            for c, text in enumerate((name, f"{value:.6g}",
                                      "—" if not np.isfinite(error) else f"{error:.2g}")):
                self.results.setItem(r, c, QTableWidgetItem(text))
        if kind == "fermi":
            self.report.setText(fit.summary().replace("\n", "<br>"))
        else:
            self.report.setText(
                f"reduced χ² = {fit.chi2:.4g}, R² = {fit.r_squared:.5f}"
                + ("" if fit.success else f"<br>solver: {fit.message}")
                + ("<br>Errors are scaled by √χ²: they assume the model is "
                   "right and the scatter is what it is."))

    def draw_data(self):
        self.plot.clear()
        x, y, sigma = self.fit_data()
        colour = _colour(self.channel.currentData() or 0)
        self.plot.plot(x, y, pen=None, symbol="o", symbolSize=5,
                       symbolBrush=colour, symbolPen=None, name=_legend("data"))
        if sigma is not None:
            self.plot.addItem(pg.ErrorBarItem(x=x, y=y, height=2 * sigma,
                                              pen=pg.mkPen(colour)))

    # -- leaving -------------------------------------------------------------
    def fit_parameters(self) -> dict:
        kind = self.last[0]
        out = {"model": kind, "channel": self.channel.currentText(),
               "window": [self.lo.value(), self.hi.value()],
               "weighting": self.weighting.currentData()}
        if kind == "fermi":
            out["reduced_chi2"] = float(self.last[1].reduced_chi2)
        else:
            out.update(background=self.background.currentText(),
                       resolution=self.resolution.value(),
                       reduced_chi2=float(self.last[1].chi2))
        for name, value, error in self.rows():
            key = name.replace(" ", "_")
            out[key] = float(value)
            if np.isfinite(error):
                out[f"{key}_err"] = float(error)
        return out

    def export(self):
        """The data, the fit, the residual and each component, as one curve
        dataset in the list."""
        if self.last is None:
            return None
        kind, fit, x, y, sigma = self.last
        model, components = self._curves()
        if kind == "fermi":
            order = np.argsort(x)
            xs, ys, model = x[order], y[order], np.asarray(model)[order]
            sigma = sigma[order] if sigma is not None else None
        else:
            xs, ys = fit.x, fit.y
            if sigma is not None:
                sigma = np.interp(xs, x, sigma)
        channel = self.channel.currentText()
        columns = [(channel, ys)]
        if sigma is not None:
            columns.append((C.sigma_name(channel), sigma))
        columns += [("fit", model), ("residual", ys - model)]
        columns += [(name, comp) for name, comp in components]
        if kind == "peaks":
            columns.append(("background", fit.background))
        xs, values, names = C.table(xs, columns)
        viewer = self.viewer
        target_kind = "mdc" if viewer.kind == "mdc" else "edc"
        name = _unique(f"{viewer.filename} [fit]", viewer.existing_names())
        data = curve_dataset(target_kind, xs, values, names,
                             x_label=viewer.x_label,
                             value_label=viewer.value_label, name=name,
                             step="curve_fit", parameters=self.fit_parameters(),
                             source_info=viewer.info, prefix="fit",
                             source=viewer.filename,
                             source_path=getattr(viewer.data, "path", ""))
        self.datasetsCreated.emit([data])
        self.report.setText(self.report.text()
                            + f"<br><b>Added “{name}” to the list.</b>")
        return data

    def shift_by_ef(self):
        """A copy of the viewer's whole dataset with the fitted E_F at 0."""
        if self.last is None or self.last[0] != "fermi":
            return None
        ef = float(self.last[1].values["ef"])
        viewer = self.viewer
        viewer.shift_value.setValue(ef)
        viewer.shift_label.setChecked(True)
        return viewer.apply_operation("shift")

    def copy_results(self):
        lines = ["parameter\tvalue\terror"]
        lines += [f"{n}\t{v:.8g}\t{e:.3g}" for n, v, e in self.rows()]
        QApplication.clipboard().setText("\n".join(lines))


# ==========================================================================
# Spin analysis
# ==========================================================================
class SpinAnalysisDialog(QDialog):
    """Polarisation, instrumental asymmetry and spin-resolved spectra from
    the viewer's spin channels."""

    datasetsCreated = pyqtSignal(list)

    def __init__(self, viewer: CurveWindow):
        super().__init__(viewer)
        self.viewer = viewer
        self.setWindowTitle(f"Spin analysis — {viewer.filename}")
        self.setModal(False)
        self.resize(1000, 900)
        self.last = None
        counts = viewer.values[:, C.data_channels(viewer.names)]
        if counts.shape[1] < 2:
            raise ValueError("a spin analysis needs at least two channels")
        self.counts = counts
        self.channels = S.channels_from_info(viewer.info, counts.shape[1])
        layout = QVBoxLayout(self)

        # -- the channels, as read -- and overridable -------------------------
        self.table = QTableWidget(len(self.channels), 5)
        self.table.setHorizontalHeaderLabels(
            ["Channel (as recorded)", "Setting", "Target", "Axis", "Counts as"])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.table.setMaximumHeight(40 + 30 * len(self.channels))
        self.table.setToolTip(
            "Read from the file's SpinComp labels. If a label was not "
            "understood, or is wrong, set the axis and sign here.")
        for row, ch in enumerate(self.channels):
            for col, text in enumerate((f"C{ch.index} {ch.label}", ch.setting,
                                        ch.magnetisation)):
                item = QTableWidgetItem(text)
                item.setFlags(item.flags() & ~Qt.ItemIsEditable)
                self.table.setItem(row, col, item)
            axis = QComboBox()
            axis.addItems(["—", "X", "Y", "Z"])
            axis.setCurrentText(ch.axis or "—")
            axis.currentIndexChanged.connect(self._channels_edited)
            sign = QComboBox()
            sign.addItem("+ (up)", 1)
            sign.addItem("− (down)", -1)
            sign.setCurrentIndex(0 if ch.sign >= 0 else 1)
            sign.currentIndexChanged.connect(self._channels_edited)
            self.table.setCellWidget(row, 3, axis)
            self.table.setCellWidget(row, 4, sign)
        layout.addWidget(self.table)

        form = QHBoxLayout()
        self.axis = QComboBox()
        self.axis.currentIndexChanged.connect(self._axis_changed)
        form.addWidget(QLabel("Axis")); form.addWidget(self.axis)
        self.method = QComboBox()
        self.method.currentIndexChanged.connect(lambda *_: self.refresh())
        form.addWidget(QLabel("Method")); form.addWidget(self.method)
        self.sherman = _spin(S.DEFAULT_SHERMAN, 0.01, 1.0, 3, 0.01, 70)
        self.sherman.setToolTip(
            "Effective Sherman function of the detector. 0.2 is the value "
            "published for this end station's FERRUM VLEED detector; use "
            "your own calibration if you have one. P scales as 1/S, and so "
            "does its error bar.")
        self.sherman.valueChanged.connect(lambda *_: self.refresh())
        form.addWidget(QLabel("S_eff")); form.addWidget(self.sherman)
        self.bin = QSpinBox()
        self.bin.setRange(1, 100)
        self.bin.setToolTip("Sum this many neighbouring points first.")
        self.bin.valueChanged.connect(lambda *_: self.refresh())
        form.addWidget(QLabel("Bin")); form.addWidget(self.bin)
        form.addStretch(1)
        layout.addLayout(form)

        zero = QHBoxLayout()
        self.zero = QCheckBox("Zero reference: take P = 0 between")
        self.zero.setToolTip(
            "Remove the asymmetry measured over a window you know to be "
            "unpolarised. It will also remove real polarisation if the window "
            "has any -- the cross ratio already cancels the instrumental "
            "asymmetry when the channels allow it.")
        self.zero.toggled.connect(lambda *_: self.refresh())
        lo, hi = float(viewer.x.min()), float(viewer.x.max())
        self.zero_lo = _spin(hi - 0.1 * (hi - lo), decimals=4)
        self.zero_hi = _spin(hi, decimals=4)
        for box in (self.zero_lo, self.zero_hi):
            box.valueChanged.connect(lambda *_: self.refresh())
        take = QPushButton("From the viewer's Range")
        take.clicked.connect(self._zero_from_viewer)
        for widget in (self.zero, self.zero_lo, QLabel("and"), self.zero_hi, take):
            zero.addWidget(widget)
        zero.addStretch(1)
        layout.addLayout(zero)

        # -- plots ------------------------------------------------------------
        self.raw_plot, self.p_plot, self.ud_plot = (pg.PlotWidget() for _ in range(3))
        for plot, label in ((self.raw_plot, "Counts"),
                            (self.p_plot, "Polarisation"),
                            (self.ud_plot, "Spin-resolved")):
            plot.setBackground("w")
            plot.showGrid(x=True, y=True, alpha=0.25)
            plot.setLabel("left", label)
            # A polarisation is a pure number: no "×0.001" on its axis.
            plot.getAxis("left").enableAutoSIPrefix(False)
            plot.addLegend(offset=(-10, 10))
            layout.addWidget(plot, stretch=1)
        self.p_plot.setXLink(self.raw_plot)
        self.ud_plot.setXLink(self.raw_plot)
        self.ud_plot.setLabel("bottom", viewer.x_label)

        self.report = QLabel("")
        self.report.setWordWrap(True)
        self.report.setTextInteractionFlags(Qt.TextSelectableByMouse)
        layout.addWidget(self.report)

        box = QDialogButtonBox()
        self.p_button = box.addButton("P to list", QDialogButtonBox.ActionRole)
        self.p_button.clicked.connect(self.export_polarisation)
        self.ud_button = box.addButton("I↑ / I↓ to list",
                                       QDialogButtonBox.ActionRole)
        self.ud_button.clicked.connect(self.export_spin_resolved)
        figure = box.addButton("As a figure", QDialogButtonBox.ActionRole)
        figure.clicked.connect(self.to_figure)
        box.addButton("Close", QDialogButtonBox.RejectRole).clicked.connect(
            self.reject)
        layout.addWidget(box)
        self._channels_edited()

    # -- inputs --------------------------------------------------------------
    def current_channels(self):
        out = []
        for row, ch in enumerate(self.channels):
            axis = self.table.cellWidget(row, 3).currentText()
            sign = self.table.cellWidget(row, 4).currentData()
            out.append(S.SpinChannel(ch.index, ch.label, ch.setting,
                                     ch.magnetisation,
                                     "" if axis == "—" else axis, int(sign)))
        return out

    def _channels_edited(self, *_):
        keep = self.axis.currentText()
        self.axis.blockSignals(True)
        self.axis.clear()
        self.axis.addItems(S.axes_available(self.current_channels()))
        if keep and self.axis.findText(keep) >= 0:
            self.axis.setCurrentText(keep)
        self.axis.blockSignals(False)
        self._axis_changed()

    def _axis_changed(self, *_):
        keep = self.method.currentData()
        self.method.blockSignals(True)
        self.method.clear()
        axis = self.axis.currentText()
        channels = self.current_channels()
        if axis:
            report = S.design_report(channels, axis)
            if report["n_plus"] == report["n_minus"] and report["n_plus"]:
                self.method.addItem("Cross ratio (all channels)", "cross")
            for setting, _p, _m in S.pairs(channels, axis):
                self.method.addItem(f"One pair: <{setting}>", f"pair:{setting}")
        if keep and self.method.findData(keep) >= 0:
            self.method.setCurrentIndex(self.method.findData(keep))
        self.method.blockSignals(False)
        self.refresh()

    def _zero_from_viewer(self):
        region = self.viewer.x_range()
        if region is None:
            QMessageBox.information(self, "Zero reference",
                                    "Switch the viewer's Range on first.")
            return
        self.zero_lo.setValue(region[0])
        self.zero_hi.setValue(region[1])
        self.zero.setChecked(True)

    def settings(self) -> dict:
        return {"axis": self.axis.currentText(),
                "method": self.method.currentData(),
                "sherman": self.sherman.value(),
                "bin_factor": self.bin.value(),
                "zero_region": ((self.zero_lo.value(), self.zero_hi.value())
                                if self.zero.isChecked() else None)}

    # -- computing -------------------------------------------------------------
    def compute(self) -> S.SpinResult:
        s = self.settings()
        if not s["axis"] or not s["method"]:
            raise ValueError("no axis has channels on both sides -- set the "
                             "axes and signs in the table")
        return S.analyse(self.viewer.x, self.counts, self.current_channels(),
                         s["axis"], sherman=s["sherman"], method=s["method"],
                         zero_region=s["zero_region"],
                         bin_factor=s["bin_factor"])

    def refresh(self, *_):
        if not hasattr(self, "report"):
            return
        try:
            result = self.compute()
        except ValueError as exc:
            self.last = None
            self.report.setText(f"<span style='color:#b00'>{exc}</span>")
            for plot in (self.p_plot, self.ud_plot):
                plot.clear()
            self.p_button.setEnabled(False)
            self.ud_button.setEnabled(False)
            self._draw_raw()
            return
        self.last = result
        self.p_button.setEnabled(True)
        self.ud_button.setEnabled(True)
        self._draw_raw()
        self._draw_result(result)
        self.report.setText(result.summary().replace("\n", "<br>"))

    def _draw_raw(self):
        self.raw_plot.clear()
        x = self.viewer.x
        for ch in self.channels:
            self.raw_plot.plot(x, self.counts[:, ch.index],
                               pen=pg.mkPen(_colour(ch.index), width=1.4),
                               name=_legend(f"C{ch.index} {ch.label}"))

    def _draw_result(self, r):
        axis = r.axis
        self.p_plot.clear()
        self.p_plot.addItem(pg.InfiniteLine(0, angle=0, pen=pg.mkPen("#999")))
        self.p_plot.plot(r.x, r.polarisation, pen=pg.mkPen("#1f5aa6", width=2),
                         symbol="o", symbolSize=5, symbolBrush="#1f5aa6",
                         name=_legend(f"P_{axis}"))
        self.p_plot.addItem(pg.ErrorBarItem(x=r.x, y=r.polarisation,
                                            height=2 * r.sigma_polarisation,
                                            pen=pg.mkPen("#1f5aa6")))
        for i, (setting, (a, _sa)) in enumerate(r.pair_asymmetries.items()):
            self.p_plot.plot(r.x, a / r.sherman,
                             pen=pg.mkPen(_colour(i + 3), width=1,
                                          style=Qt.DashLine),
                             name=_legend(f"<{setting}> alone"))
        self.p_plot.setLabel("left", f"P_{axis}")
        self.ud_plot.clear()
        for y, s, colour, name in ((r.up, r.sigma_up, "#c23b22", f"I↑ (+{axis})"),
                                   (r.down, r.sigma_down, "#1f5aa6",
                                    f"I↓ (−{axis})")):
            self.ud_plot.plot(r.x, y, pen=pg.mkPen(colour, width=2), name=_legend(name))
            self.ud_plot.addItem(pg.ErrorBarItem(x=r.x, y=y, height=2 * s,
                                                 pen=pg.mkPen(colour)))

    # -- leaving -----------------------------------------------------------------
    def _parameters(self, r) -> dict:
        s = self.settings()
        out = {"axis": r.axis, "method": r.method, "sherman": r.sherman,
               "bin_factor": s["bin_factor"],
               "channels_used": ", ".join(f"C{i}" for i in r.channels_used),
               "cancels_transmission": bool(r.design.get("transmission")),
               "cancels_reflectivity": bool(r.design.get("reflectivity"))}
        if s["zero_region"] is not None:
            out["zero_region"] = list(s["zero_region"])
            out["zero_asymmetry"] = float(np.tanh(r.zero_offset / 2))
        if r.instrumental is not None:
            out["instrumental_asymmetry"] = r.instrumental[0]
            out["instrumental_asymmetry_err"] = r.instrumental[1]
        return out

    def _emit(self, columns, suffix, value_label, step):
        viewer = self.viewer
        r = self.last
        x, values, names = C.table(r.x, columns)
        name = _unique(f"{viewer.filename} [{suffix}]", viewer.existing_names())
        data = curve_dataset("edc", x, values, names, x_label=viewer.x_label,
                             value_label=value_label, name=name, step=step,
                             parameters=self._parameters(r),
                             source_info=viewer.info, prefix="spin",
                             source=viewer.filename,
                             source_path=getattr(viewer.data, "path", ""))
        self.datasetsCreated.emit([data])
        self.report.setText(self.report.text()
                            + f"<br><b>Added “{name}” to the list.</b>")
        return data

    def export_polarisation(self):
        r = self.last
        if r is None:
            return None
        label = f"P_{r.axis}"
        return self._emit([(label, r.polarisation),
                           (C.sigma_name(label), r.sigma_polarisation)],
                          f"P{r.axis}", f"Spin polarisation {label}",
                          "spin_polarisation")

    def export_spin_resolved(self):
        r = self.last
        if r is None:
            return None
        up, down = f"I↑ (+{r.axis})", f"I↓ (−{r.axis})"
        return self._emit([(up, r.up), (down, r.down),
                           (C.sigma_name(up), r.sigma_up),
                           (C.sigma_name(down), r.sigma_down)],
                          f"spin {r.axis}",
                          f"Counts (sum of {len(r.channels_used)} channels)",
                          "spin_resolved")

    def to_figure(self):
        from ui.figure import FigureWindow

        r = self.last
        if r is None:
            return None
        panels = [
            curve_panel(r.x, [(f"P_{r.axis}", r.polarisation, "#1f5aa6")],
                        x_label=self.viewer.x_label, y_label=f"P_{r.axis}",
                        name="polarisation",
                        errors={f"P_{r.axis}": r.sigma_polarisation}),
            curve_panel(r.x, [(f"I↑ (+{r.axis})", r.up, "#c23b22"),
                              (f"I↓ (−{r.axis})", r.down, "#1f5aa6")],
                        x_label=self.viewer.x_label, y_label="Counts",
                        name="spin-resolved",
                        errors={f"I↑ (+{r.axis})": r.sigma_up,
                                f"I↓ (−{r.axis})": r.sigma_down}),
        ]
        window = FigureWindow(panels, parent=self)
        window.show()
        return window
