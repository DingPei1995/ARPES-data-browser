"""
ui/panels/fit.py
================
The MDC fit and EDC fit panels and what comes out of them: the dispersion
window (v_F, m*, the window-sensitivity scan) and the self-energy window.

The panel is built around **seeding**, the one idea worth keeping from
``fit_MDC_demo``: mark a band at four or five lines and the rest is
interpolated. Automatic peak finding fails where it matters -- crossings,
weak features, merging peaks -- and pointing at the band takes seconds.

Two panels, one class
---------------------
:class:`FitPanel` is opened as an **MDC fit** (``mode="mdc"``) or an **EDC
fit** (``mode="edc"``): the same seeding and fitting, along one direction
each. An EDC fit takes any cut -- an EDC's peak energy means the same in
degrees or in Å⁻¹ -- but its dispersion window, which turns the band into
v_F and m*, opens only for a momentum cut, for the reason below.

Momentum cuts only, for the MDC fit
-----------------------------------
The MDC fit opens for a cut whose first axis is a momentum: one converted by
"Cut k conversion", or a slice of a k-map. Not for an angle cut. A velocity is
eV*A, an effective mass is hbar^2/2m_e over eV*A^2, and a self-energy needs
both -- feed degrees into any of them and the numbers come out wrong by a
factor nobody can reconstruct afterwards, while looking perfectly
reasonable. Converting first is one button, so the panel asks for it rather
than carrying a unit through every formula and hoping.
"""
from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PyQt5.QtCore import Qt, pyqtSignal, QTimer
from PyQt5.QtWidgets import (QWidget, QMainWindow, QVBoxLayout, QHBoxLayout,
                             QLabel, QPushButton, QComboBox, QCheckBox,
                             QFormLayout, QDoubleSpinBox, QSpinBox, QLineEdit,
                             QApplication, QTabWidget, QGroupBox, QMessageBox,
                             QProgressDialog, QSplitter, QTableWidget,
                             QTableWidgetItem, QHeaderView, QFileDialog,
                             QAbstractItemView, QSizePolicy, QDialog,
                             QDialogButtonBox)

from ui.common.screen import fits_screen, wheel_scrolls_window

import tools.spectra.peaks as PK
import tools.spectra.dispersion as DISP
from ui.viewers.widgets import (MemoryData, apply_colormap, strip_stock_menu,
                                plain_image_view, show_frame, fit_frame_view,
                                COLORMAP_NAMES)
from ui.panels.process import spin, whole

#: Anything whose axis label carries one of these is a momentum.
MOMENTUM_MARKS = ("Å⁻¹", "A^-1", "1/Å", "Å-1", "AA^-1")
#: Colours for the bands, in order. Chosen to stay apart on a grey image.
BAND_COLOURS = ("#e8413c", "#2d9bf0", "#20a464", "#f5a623", "#a855c8",
                "#00b2b2", "#d64ea0", "#8a8f00")


def momentum_cut_reason(data) -> str:
    """Empty if this dataset can be fitted here, otherwise why not."""
    if getattr(data, "kind", "") != "cut":
        return ("The fit panel works on a 2-D cut. Take a slice of a map "
                "first.")
    labels = getattr(getattr(data, "scan", None), "labels", {}) or {}
    x_label = str(labels.get("x", ""))
    info = getattr(getattr(data, "scan", None), "info", {}) or {}
    converted = any(str(key).startswith("kcut.") for key in info)
    if converted or any(mark in x_label for mark in MOMENTUM_MARKS):
        return ""
    return (f"This cut's first axis is “{x_label or 'unnamed'}”, which is "
            f"not a momentum.\n\nFit only k-converted cuts: a velocity is "
            f"eV·Å and an effective mass is ħ²/2mₑ over eV·Å², so an "
            f"angle axis makes every number here wrong by a factor that "
            f"cannot be recovered later — while still looking reasonable.\n\n"
            f"Use “Convert to momentum (k)” in the cut viewer first, or take the "
            f"slice from a converted k-map.")


def edc_cut_reason(data) -> str:
    """Empty if the EDC fit can work on this dataset, otherwise why not:
    any cut will do."""
    if getattr(data, "kind", "") != "cut":
        return ("The EDC fit works on a 2-D cut. Take a slice of a map "
                "first.")
    return ""


#: the panel's two modes, by the title each window carries
FIT_MODES = {"mdc": "MDC fit (k-cut only)", "edc": "EDC fit"}


def fermi_defaults(data):
    """``(ef, temperature, resolution)`` read off the dataset if it says.

    A cut whose energy axis has already had E_F subtracted -- the Fermi
    round writes ``fitEF.ef_subtracted`` and relabels the axis -- has E_F at
    zero, and saying so saves the user typing a number the file already
    knows.
    """
    info = getattr(getattr(data, "scan", None), "info", {}) or {}
    labels = getattr(getattr(data, "scan", None), "labels", {}) or {}
    ef, temperature, resolution = None, None, None
    for key, value in info.items():
        tail = str(key).rsplit(".", 1)[-1].lower()
        try:
            number = float(value)
        except (TypeError, ValueError):
            continue
        if tail == "ef_subtracted":
            ef = 0.0
        elif tail == "ef" and ef is None:
            ef = number
        elif tail == "temperature" and temperature is None:
            temperature = number
        elif tail == "resolution" and resolution is None:
            resolution = number
    y_label = str(labels.get("y", ""))
    if ef is None and ("E_F" in y_label or "E-Ef" in y_label or "E − E_F" in y_label):
        ef = 0.0
    return (0.0 if ef is None else ef,
            30.0 if temperature is None else temperature,
            0.0 if resolution is None else resolution)


# ==========================================================================
# The panel
# ==========================================================================
@fits_screen
class FitPanel(QMainWindow):
    """Seed the bands on the picture, fit every line, read off the band."""

    datasetsCreated = pyqtSignal(list)
    closed = pyqtSignal(object)
    #: a seed was placed (the guide moves on to the next line)
    seedPlaced = pyqtSignal()
    #: the whole series was fitted
    seriesFitted = pyqtSignal()

    def __init__(self, data, label: str, parent=None, colormap="gray",
                 flip=False, mode: str = "mdc"):
        super().__init__(parent)
        if mode not in FIT_MODES:
            raise ValueError(f"unknown fit {mode!r}")
        self.mode = mode
        self.data = data
        #: whether the cut is in momentum, which v_F and m* need
        self.momentum = momentum_cut_reason(data) == ""
        self.label = label
        scan = data.scan
        self.values = np.asarray(scan.value, dtype=float)
        self.axes = [np.asarray(scan.x, dtype=float),
                     np.asarray(scan.y, dtype=float)]
        self.labels = [scan.labels.get("x", "k"), scan.labels.get("y", "E")]
        self.colormap, self.flip = colormap, flip
        self.bands = []
        self.series = None
        self._seeding = False

        # Made before the widgets: building the tabs fires the signals that
        # queue a redraw, and the timer has to exist by then.
        self._redraw = QTimer(self)
        self._redraw.setSingleShot(True)
        self._redraw.setInterval(90)
        self._redraw.timeout.connect(self.refresh_line)

        self.setWindowTitle(f"{FIT_MODES[mode]} — {label}")
        self.resize(1240, 860)

        central = QWidget()
        root = QVBoxLayout(central)
        root.setContentsMargins(6, 6, 6, 6)

        # The guided fit, first in sight: most users do not know where a fit
        # starts, and the button says it.
        start_row = QWidget()
        start_bar = QHBoxLayout(start_row)
        start_bar.setContentsMargins(0, 0, 0, 0)
        self.guide_button = QPushButton("\u25b6  Guided fit (step by step)")
        self.guide_button.setToolTip(
            "Walks through the fit one step at a time -- range, bands, "
            "seeds, a test line, the series -- and says what each step is "
            "for.")
        self.guide_button.setStyleSheet("font-weight: bold; padding: 4px 10px;")
        self.guide_button.clicked.connect(self.start_guide)
        start_bar.addWidget(self.guide_button)
        hint = QLabel("<i>New to this panel? The guide takes you from the "
                      "cut to the fitted band.  Ctrl + wheel zooms a plot.</i>")
        hint.setStyleSheet("color: #666;")
        start_bar.addWidget(hint, 1)
        root.addWidget(start_row)

        # Plots above, controls below, with a handle between them: the
        # controls are never pushed out of the window by the pictures, and
        # the user decides how much room each gets.
        self.splitter = QSplitter(Qt.Vertical)
        self.splitter.setChildrenCollapsible(False)
        top = QSplitter(Qt.Horizontal)
        top.addWidget(self._image_side())
        top.addWidget(self._line_side())
        top.setSizes([620, 620])
        self.splitter.addWidget(top)

        # The guide sits beside the tabs rather than above them, so opening
        # it takes no height from the pictures.
        controls = QWidget()
        side_by_side = QHBoxLayout(controls)
        side_by_side.setContentsMargins(0, 0, 0, 0)
        self.guide = FitGuide(self)
        self.guide.hide()
        side_by_side.addWidget(self.guide, 2)
        column = QWidget()
        below = QVBoxLayout(column)
        below.setContentsMargins(0, 0, 0, 0)
        side_by_side.addWidget(column, 3)
        self.tabs = QTabWidget()
        self.tabs.addTab(self._bands_tab(), "Bands")
        self.tabs.addTab(self._settings_tab(), "Fit settings")
        self.tabs.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Maximum)
        below.addWidget(self.tabs, 0)

        self.note = QLabel("Add a band, then click on the curve to place a "
                           "seed at this line. Four or five seeds down the "
                           "band are enough -- or press “Guided fit” above.")
        self.note.setWordWrap(True)
        below.addWidget(self.note)
        below.addStretch(1)
        self.splitter.addWidget(controls)
        self.splitter.setStretchFactor(0, 1)
        self.splitter.setStretchFactor(1, 0)
        self.splitter.setSizes([480, 320])
        root.addWidget(self.splitter, 1)
        self.setCentralWidget(central)

        self.show_image()
        self.add_band()
        self.refresh_line()
        wheel_scrolls_window(self.image.ui.graphicsView, self.curve_plot,
                             self.residual_plot)

    # -- layout ------------------------------------------------------------
    def _image_side(self) -> QWidget:
        box = QGroupBox("Cut")
        self.cut_box = box
        layout = QVBoxLayout(box)
        layout.setContentsMargins(4, 4, 4, 4)
        self.image = plain_image_view(self.labels[0], self.labels[1])
        self.image.setMinimumHeight(180)
        layout.addWidget(self.image)

        # The line being fitted, dragged with the mouse.
        self.cursor = pg.InfiniteLine(angle=0, movable=True,
                                      pen=pg.mkPen("#1f6f8b", width=3),
                                      hoverPen=pg.mkPen("#f5a623", width=4))
        self.cursor.setZValue(20)
        self.image.view.addItem(self.cursor)
        self.cursor.sigPositionChanged.connect(self._queue_line)

        self.seed_markers = {}      # band name -> ScatterPlotItem on the image
        self.fit_markers = {}
        row = QWidget()
        bar = QHBoxLayout(row)
        bar.setContentsMargins(0, 0, 0, 0)
        bar.addWidget(QLabel("Line at"))
        self.position_box = QDoubleSpinBox()
        self.position_box.setDecimals(5)
        self.position_box.setKeyboardTracking(False)
        self.position_box.setMaximumWidth(110)
        self.position_box.valueChanged.connect(self._position_typed)
        bar.addWidget(self.position_box)
        self.position_label = QLabel("")
        bar.addWidget(self.position_label)
        bar.addStretch(1)
        layout.addWidget(row)
        return box

    def _line_side(self) -> QWidget:
        box = QGroupBox("Line")
        self.line_box = box
        layout = QVBoxLayout(box)
        layout.setContentsMargins(4, 4, 4, 4)
        self.curve_plot = pg.PlotWidget()
        self.curve_plot.setMinimumHeight(140)
        self.curve_plot.addLegend(offset=(-10, 10))
        self.curve_plot.setLabel("left", "Intensity")
        strip_stock_menu(self.curve_plot.getPlotItem())
        self.curve_plot.scene().sigMouseClicked.connect(self._curve_clicked)
        layout.addWidget(self.curve_plot, 3)

        self.residual_plot = pg.PlotWidget()
        self.residual_plot.setLabel("left", "residual")
        self.residual_plot.setMaximumHeight(150)
        self.residual_plot.setMinimumHeight(60)
        self.residual_plot.setXLink(self.curve_plot)
        strip_stock_menu(self.residual_plot.getPlotItem())
        layout.addWidget(self.residual_plot, 1)

        row = QWidget()
        bar = QHBoxLayout(row)
        bar.setContentsMargins(0, 0, 0, 0)
        self.seed_button = QPushButton("Place a seed")
        self.seed_button.setCheckable(True)
        self.seed_button.setToolTip(
            "Click the top of the band on this curve. The height and width "
            "are read off the data; only the position is yours to give.")
        self.seed_button.toggled.connect(self._arm_seeding)
        bar.addWidget(self.seed_button)
        preview = QPushButton("Try this line")
        self.try_button = preview
        preview.setToolTip("Fit this one line, to check the settings before "
                           "committing to the whole series.")
        preview.clicked.connect(self.try_line)
        bar.addWidget(preview)
        bar.addStretch(1)
        self.line_caption = QLabel("")
        bar.addWidget(self.line_caption)
        layout.addWidget(row)
        return box

    def _bands_tab(self) -> QWidget:
        page = QWidget()
        layout = QHBoxLayout(page)
        self.band_table = QTableWidget(0, 5)
        self.band_table.setHorizontalHeaderLabels(
            ["Band", "Shape", "Same width as", "Seeds", "Range"])
        self.band_table.verticalHeader().setVisible(False)
        self.band_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.band_table.setSelectionMode(QAbstractItemView.SingleSelection)
        self.band_table.horizontalHeader().setSectionResizeMode(
            QHeaderView.ResizeToContents)
        self.band_table.setMaximumHeight(170)
        self.band_table.setMinimumHeight(90)
        self.band_table.itemSelectionChanged.connect(self._queue_line)
        layout.addWidget(self.band_table, 1)

        side = QWidget()
        column = QVBoxLayout(side)
        column.setContentsMargins(0, 0, 0, 0)
        for text, slot, tip in (
                ("Add band", self.add_band, "A new band to seed."),
                ("Remove band", self.remove_band, ""),
                ("Clear its seeds", self.clear_seeds, ""),
                ("Undo last seed", self.undo_seed, "")):
            button = QPushButton(text)
            if tip:
                button.setToolTip(tip)
            button.clicked.connect(slot)
            column.addWidget(button)
        column.addStretch(1)
        layout.addWidget(side)
        return page

    def _settings_tab(self) -> QWidget:
        page = QWidget()
        outer = QHBoxLayout(page)

        left = QWidget()
        form = QFormLayout(left)
        form.setFieldGrowthPolicy(QFormLayout.FieldsStayAtSizeHint)
        self.direction = QComboBox()
        self.direction.addItems([f"MDCs: lines along {self._short(0)}",
                                 f"EDCs: lines along {self._short(1)}"])
        self.direction.setMaximumWidth(230)
        # Fixed by the mode: the MDC fit and the EDC fit are two panels.
        self.direction.setCurrentIndex(0 if self.mode == "mdc" else 1)
        self.direction.currentIndexChanged.connect(self._direction_changed)
        form.addRow("Fit", self.direction)
        self.direction.setVisible(False)
        form.labelForField(self.direction).setVisible(False)
        self.shape = QComboBox()
        self.shape.addItems(PK.PEAK_SHAPES)
        self.shape.setMaximumWidth(140)
        self.shape.currentIndexChanged.connect(self._shape_changed)
        form.addRow("Line shape", self.shape)
        self.background = QComboBox()
        self.background.addItems(PK.BACKGROUNDS)
        self.background.setCurrentText("linear")
        self.background.setMaximumWidth(140)
        self.background.currentIndexChanged.connect(self._queue_line)
        form.addRow("Background", self.background)
        self.resolution = spin(0.0, 0.0, 10.0, 0.005, 4, " eV")
        self.resolution.setToolTip(
            "Gaussian FWHM of the instrument. With a Voigt shape this is "
            "held fixed and the Lorentzian width is fitted, which is the "
            "only way a width becomes an intrinsic one.")
        self.resolution.valueChanged.connect(self._queue_line)
        form.addRow("Resolution", self.resolution)
        self.combine = whole(1, 1, 99)
        self.combine.setToolTip("Sum this many neighbouring lines per fit.")
        form.addRow("Combine lines", self.combine)
        self.step = whole(1, 1, 99)
        form.addRow("Fit every", self.step)
        self.from_edit, self.to_edit = QLineEdit(), QLineEdit()
        for edit, which in ((self.from_edit, "from"), (self.to_edit, "to")):
            edit.setPlaceholderText(which)
            edit.setMaximumWidth(90)
        span = QWidget()
        bar = QHBoxLayout(span)
        bar.setContentsMargins(0, 0, 0, 0)
        bar.addWidget(self.from_edit)
        bar.addWidget(self.to_edit)
        bar.addStretch(1)
        form.addRow("Over", span)
        self.refine = QCheckBox("Second pass from the first result")
        self.refine.setChecked(True)
        self.refine.setToolTip(
            "The first pass follows the seeds, which are sparse; the second "
            "follows the fitted dispersion, which is dense.")
        form.addRow("", self.refine)
        outer.addWidget(left)

        self.fermi_group = QGroupBox("Fermi cut-off (EDC only)")
        self.fermi_group.setCheckable(True)
        self.fermi_group.setChecked(False)
        fermi = QFormLayout(self.fermi_group)
        fermi.setFieldGrowthPolicy(QFormLayout.FieldsStayAtSizeHint)
        ef, temperature, resolution = fermi_defaults(self.data)
        self.ef = spin(ef, -1e4, 1e4, 0.005, 5, " eV")
        self.temperature = spin(temperature, 0.1, 1000.0, 5.0, 1, " K")
        for widget, name in ((self.ef, "E_F"), (self.temperature, "Temperature")):
            widget.valueChanged.connect(self._queue_line)
            fermi.addRow(name, widget)
        self.fermi_group.toggled.connect(self._queue_line)
        note = QLabel(
            "<i>An EDC near E_F is (A·f)⊗R, not A: the occupation multiplies "
            "the spectrum before the resolution smears it. Left out, a peak "
            "10 meV above E_F comes back 20 meV too high. Dividing the "
            "cut-off out of the data first is worse — it divides the noise "
            "up too.</i>")
        note.setWordWrap(True)
        note.setMaximumWidth(330)
        fermi.addRow("", note)
        if resolution > 0:
            self.resolution.setValue(resolution)
        self.fermi_group.setTitle("Fermi cut-off")
        self.fermi_group.setVisible(self.mode == "edc")
        outer.addWidget(self.fermi_group)

        buttons = QWidget()
        column = QVBoxLayout(buttons)
        column.setContentsMargins(0, 0, 0, 0)
        run = QPushButton("Fit the series")
        self.run_button = run
        run.clicked.connect(self.run_fit)
        column.addWidget(run)
        self.dispersion_button = QPushButton("Dispersion...")
        self.dispersion_button.setEnabled(False)
        if not self.momentum:
            self.dispersion_button.setToolTip(
                "v_F and m* need a cut in momentum: this one's first axis is "
                "an angle. Convert the cut to k first; the fitted peaks can "
                "still be saved as a table here.")
        self.dispersion_button.clicked.connect(self.open_dispersion)
        column.addWidget(self.dispersion_button)
        self.table_button = QPushButton("Save table...")
        self.table_button.setEnabled(False)
        self.table_button.clicked.connect(self.save_table)
        column.addWidget(self.table_button)
        column.addStretch(1)
        outer.addWidget(buttons)
        return page

    def _short(self, dim: int) -> str:
        return self.labels[dim].split(" (")[0]

    # -- the image ---------------------------------------------------------
    def show_image(self):
        finite = self.values[np.isfinite(self.values)]
        levels = (np.percentile(finite, [1, 99]) if finite.size else (0, 1))
        show_frame(self.image, self.values, self.axes[0], self.axes[1],
                   levels=(float(levels[0]), float(levels[1])))
        apply_colormap(self.image, self.colormap, self.flip)
        self._fit_view()
        self._set_cursor_axis()

    def _fit_view(self):
        fit_frame_view(self.image, self.axes[0], self.axes[1])

    def showEvent(self, event):
        """Fit the view once the window actually has a size.

        The first ``_fit_view`` runs in the constructor, where the view is
        still 0 x 0 and the ranges it sets do not survive the first layout.
        """
        super().showEvent(event)
        if not getattr(self, "_shown", False):
            self._shown = True
            QTimer.singleShot(0, self._fit_view)

    def _set_cursor_axis(self):
        """Point the line, and the position box, at the axis the lines are
        spaced along -- which swaps when MDC becomes EDC."""
        mdc = self.direction.currentIndex() == 0 if hasattr(self, "direction") else True
        axis = self.axes[1] if mdc else self.axes[0]
        self.cursor.setAngle(0 if mdc else 90)
        middle = float(axis[axis.size // 2])
        self.cursor.blockSignals(True)
        self.cursor.setPos(middle)
        self.cursor.blockSignals(False)
        self.position_box.blockSignals(True)
        self.position_box.setRange(float(axis.min()), float(axis.max()))
        step = abs(float(axis[-1] - axis[0])) / max(axis.size - 1, 1)
        self.position_box.setSingleStep(step or 0.01)
        self.position_box.setValue(middle)
        self.position_box.blockSignals(False)
        self.position_label.setText(self.labels[1] if mdc else self.labels[0])

    def current_position(self) -> float:
        return float(self.cursor.value())

    def _position_typed(self, value: float):
        self.cursor.blockSignals(True)
        self.cursor.setPos(float(value))
        self.cursor.blockSignals(False)
        self._queue_line()

    def _queue_line(self, *_):
        self._redraw.start()

    def _direction_changed(self, *_):
        for band in self.bands:
            band.seeds.clear()
        self._set_cursor_axis()
        self.series = None
        self.dispersion_button.setEnabled(False)
        self.table_button.setEnabled(False)
        self.sync_bands()
        self.note.setText(
            "Direction changed, so the seeds were cleared: a seed is a point "
            "on one line, and the lines now run the other way.")
        self.refresh_line()

    def _shape_changed(self, *_):
        for band in self.bands:
            band.shape = self.shape.currentText()
        self._queue_line()

    # -- bands -------------------------------------------------------------
    def add_band(self):
        name = f"Band {len(self.bands) + 1}"
        taken = {b.name for b in self.bands}
        n = len(self.bands) + 1
        while name in taken:
            n += 1
            name = f"Band {n}"
        self.bands.append(PK.Band(name, [], shape=self.shape.currentText()
                                  if hasattr(self, "shape") else "lorentzian",
                                  colour=len(self.bands)))
        self.sync_bands()
        self.band_table.selectRow(len(self.bands) - 1)
        return self.bands[-1]

    def remove_band(self):
        index = self.selected_band_index()
        if index is None or len(self.bands) <= 1:
            return
        name = self.bands[index].name
        for band in self.bands:
            if band.share_width == name:
                band.share_width = ""
        self.bands.pop(index)
        for item in (self.seed_markers, self.fit_markers):
            marker = item.pop(name, None)
            if marker is not None:
                self.image.view.removeItem(marker)
        self.sync_bands()
        self.refresh_line()

    def clear_seeds(self):
        band = self.selected_band()
        if band is not None:
            band.seeds.clear()
            self.sync_bands()
            self.refresh_line()

    def undo_seed(self):
        band = self.selected_band()
        if band is not None and band.seeds:
            band.seeds.pop()
            self.sync_bands()
            self.refresh_line()

    def selected_band_index(self):
        rows = self.band_table.selectionModel().selectedRows()
        if not rows:
            return None
        row = rows[0].row()
        return row if 0 <= row < len(self.bands) else None

    def selected_band(self):
        index = self.selected_band_index()
        return None if index is None else self.bands[index]

    def sync_bands(self):
        self.band_table.blockSignals(True)
        keep = self.selected_band_index()
        self.band_table.setRowCount(len(self.bands))
        for row, band in enumerate(self.bands):
            span = band.span()
            cells = [band.name, band.shape, band.share_width or "—",
                     str(len(band.seeds)),
                     "—" if span is None else f"{span[0]:.4g} .. {span[1]:.4g}"]
            for column, text in enumerate(cells):
                item = QTableWidgetItem(text)
                if column == 0:
                    item.setForeground(pg.mkColor(
                        BAND_COLOURS[band.colour % len(BAND_COLOURS)]))
                item.setFlags(item.flags() & ~Qt.ItemIsEditable)
                self.band_table.setItem(row, column, item)
        self.band_table.blockSignals(False)
        if keep is not None and keep < len(self.bands):
            self.band_table.selectRow(keep)
        elif self.bands:
            self.band_table.selectRow(0)
        self._draw_seed_markers()

    def _draw_seed_markers(self):
        """Seeds on the image, in the band's colour."""
        mdc = self.direction.currentIndex() == 0 if hasattr(self, "direction") else True
        for band in self.bands:
            marker = self.seed_markers.get(band.name)
            if marker is None:
                marker = pg.ScatterPlotItem(
                    size=12, symbol="+",
                    pen=pg.mkPen(BAND_COLOURS[band.colour % len(BAND_COLOURS)],
                                 width=2), brush=None)
                marker.setZValue(30)
                self.image.view.addItem(marker)
                self.seed_markers[band.name] = marker
            if band.seeds:
                centres = [s.centre for s in band.seeds]
                positions = [s.position for s in band.seeds]
                xs, ys = (centres, positions) if mdc else (positions, centres)
                marker.setData(xs, ys)
            else:
                marker.setData([], [])

    # -- seeding -----------------------------------------------------------
    def _arm_seeding(self, on: bool):
        self._seeding = bool(on)
        self.curve_plot.setCursor(Qt.CrossCursor if on else Qt.ArrowCursor)
        if on:
            self.note.setText("Click the top of the band on the curve.")

    def _curve_clicked(self, event):
        if not self._seeding or event.button() != Qt.LeftButton:
            return
        view_box = self.curve_plot.getPlotItem().getViewBox()
        if not view_box.sceneBoundingRect().contains(event.scenePos()):
            return
        band = self.selected_band()
        if band is None:
            return
        point = view_box.mapSceneToView(event.scenePos())
        event.accept()
        along, line = self.current_line()
        height, width = PK.suggest_seed(along, line, float(point.x()))
        band.seeds = [s for s in band.seeds
                      if abs(s.position - self.current_position()) > 1e-12]
        band.seeds.append(PK.Seed(self.current_position(), float(point.x()),
                                  width, height))
        self.seed_button.setChecked(False)
        self.sync_bands()
        self.note.setText(
            f"{band.name}: seed at {self._short(0)} = {point.x():.5g}, "
            f"width {width:.4g}. {len(band.seeds)} seed(s) on this band.")
        self.refresh_line()
        self.seedPlaced.emit()

    # -- the current line --------------------------------------------------
    def current_line(self):
        """``(axis, values)`` of the line under the cursor, binned."""
        mdc = self.direction.currentIndex() == 0
        along = self.axes[0] if mdc else self.axes[1]
        across = self.axes[1] if mdc else self.axes[0]
        index = int(np.argmin(np.abs(across - self.current_position())))
        half = max(0, (int(self.combine.value()) - 1) // 2)
        first, last = max(0, index - half), min(across.size, index + half + 1)
        block = (self.values[:, first:last] if mdc
                 else self.values[first:last, :])
        return along, np.nansum(block, axis=1 if mdc else 0)

    def settings(self) -> PK.FitSettings:
        mdc = self.direction.currentIndex() == 0
        return PK.FitSettings(
            background=self.background.currentText(),
            resolution=self.resolution.value(),
            fermi=bool(self.fermi_group.isChecked() and not mdc),
            ef=self.ef.value(), temperature=self.temperature.value())

    def refresh_line(self):
        along, line = self.current_line()
        self.curve_plot.clear()
        self.residual_plot.clear()
        self.curve_plot.setLabel("bottom",
                                 self.labels[0] if self.direction.currentIndex() == 0
                                 else self.labels[1])
        self.curve_plot.plot(along, line, pen=None, symbol="o", symbolSize=3,
                             symbolBrush="#888", symbolPen=None, name="data")
        position = self.current_position()
        self.position_box.blockSignals(True)
        self.position_box.setValue(position)
        self.position_box.blockSignals(False)

        active = [b for b in self.bands if b.seeds]
        if not active:
            self.line_caption.setText("no seeds yet")
            return
        # The starting guess, drawn as it will be handed to the fitter.
        for band in active:
            centre, width, height = (float(v[0]) for v in band.guess([position]))
            colour = BAND_COLOURS[band.colour % len(BAND_COLOURS)]
            profile = height * PK.peak_profile(band.shape, along, centre, width,
                                               self.resolution.value())
            self.curve_plot.plot(along, profile + float(np.nanmin(line)),
                                 pen=pg.mkPen(colour, width=1,
                                              style=Qt.DashLine),
                                 name=f"{band.name} (guess)")
        self.line_caption.setText(
            f"{len(active)} band(s) here · guess shown dashed")

    def try_line(self):
        """Fit the line under the cursor and draw it."""
        position = self.current_position()
        along, line = self.current_line()
        tolerance = self._tolerance()
        active = [b for b in self.bands if b.covers(position, tolerance)]
        if not active:
            QMessageBox.information(self, "Fit",
                                    "No band is seeded at this line.")
            return None
        try:
            fit = PK.fit_line(along, line, active, self.settings(),
                              position=position)
        except Exception as exc:
            QMessageBox.warning(self, "Fit", str(exc))
            return None
        self.draw_fit(fit)
        return fit

    def _tolerance(self) -> float:
        across = self.axes[1] if self.direction.currentIndex() == 0 else self.axes[0]
        return abs(float(across[1] - across[0])) if across.size > 1 else 0.0

    def draw_fit(self, fit: PK.LineFit):
        self.refresh_line()
        self.curve_plot.plot(fit.x, fit.model, pen=pg.mkPen("#bd4921", width=2),
                             name="fit")
        self.curve_plot.plot(fit.x, fit.background,
                             pen=pg.mkPen("#999", width=1, style=Qt.DashLine),
                             name="background")
        for index, name in enumerate(fit.names):
            band = next((b for b in self.bands if b.name == name), None)
            colour = BAND_COLOURS[(band.colour if band else index)
                                  % len(BAND_COLOURS)]
            self.curve_plot.plot(fit.x, fit.component(index),
                                 pen=pg.mkPen(colour, width=1), name=name)
        self.residual_plot.plot(fit.x, fit.residual,
                                pen=pg.mkPen("#444", width=1))
        self.residual_plot.addLine(y=0, pen=pg.mkPen("#bbb", style=Qt.DashLine))
        pieces = [f"χ² = {fit.chi2:.3g}", f"R² = {fit.r_squared:.5f}"]
        for index, name in enumerate(fit.names):
            pieces.append(f"{name}: {fit.centres[index]:.5g} "
                          f"± {fit.centre_errors[index]:.2g}, "
                          f"FWHM {fit.widths[index]:.4g}")
        self.line_caption.setText("   ·   ".join(pieces))

    # -- the whole series --------------------------------------------------
    def run_fit(self):
        bands = [b for b in self.bands if b.seeds]
        if not bands:
            QMessageBox.information(
                self, "Fit", "Seed at least one band first: select it, press "
                             "“Place a seed” and click the band on the curve.")
            return None
        if any(len(b.seeds) < 2 for b in bands):
            answer = QMessageBox.question(
                self, "Fit",
                "A band with a single seed is fitted at that one line only, "
                "because its range comes from its seeds. Carry on?",
                QMessageBox.Yes | QMessageBox.No, QMessageBox.Yes)
            if answer != QMessageBox.Yes:
                return None

        bounds = (self._number(self.from_edit), self._number(self.to_edit))
        direction = "mdc" if self.direction.currentIndex() == 0 else "edc"
        # Application modal, not window modal: this loop pumps events to keep
        # the bar moving (the fit hands its answer straight back, so it stays
        # on the GUI thread -- see ui.common.jobs for the ones that do not), and a
        # window-modal dialog would still let a click on another viewer
        # re-enter this code on top of a half-finished fit.
        progress = QProgressDialog("Fitting lines...", "Stop", 0, 100, self)
        progress.setWindowModality(Qt.ApplicationModal)
        progress.setMinimumDuration(300)

        def tick(done, total):
            progress.setMaximum(total)
            progress.setValue(done)
            QApplication.processEvents()
            return not progress.wasCanceled()

        try:
            series = PK.fit_series(
                self.values, self.axes, bands, self.settings(),
                direction=direction, step=self.step.value(),
                combine=self.combine.value(), bounds=bounds,
                refine=self.refine.isChecked(), progress=tick)
        except Exception as exc:
            progress.close()
            QMessageBox.warning(self, "Fit", str(exc))
            return None
        finally:
            progress.close()

        series.x_label, series.y_label = self.labels
        self.series = series
        self.dispersion_button.setEnabled(self.momentum)
        self.table_button.setEnabled(True)
        self._draw_fitted_positions()
        done = int(np.isfinite(series.centres).any(axis=1).sum())
        self.note.setText(
            f"Fitted {done} of {len(series.positions)} lines, median "
            f"χ² = {np.nanmedian(series.chi2):.3g}. The fitted positions are "
            f"on the picture, sized by amplitude. “Dispersion...” for v_F "
            f"and m*.")
        self.seriesFitted.emit()
        return series

    def _draw_fitted_positions(self):
        """Fitted centres on the image, **sized by amplitude**.

        Borrowed from ``fit_MDC_demo``, and the cheapest honest quality
        indicator there is: where the band fades out the markers shrink, so
        a dispersion held up by three counts looks like what it is.
        """
        if self.series is None:
            return
        mdc = self.series.direction == "mdc"
        heights = np.asarray(self.series.heights, dtype=float)
        top = float(np.nanmax(heights)) if np.isfinite(heights).any() else 1.0
        for index, name in enumerate(self.series.names):
            marker = self.fit_markers.get(name)
            band = next((b for b in self.bands if b.name == name), None)
            colour = BAND_COLOURS[(band.colour if band else index)
                                  % len(BAND_COLOURS)]
            if marker is None:
                marker = pg.ScatterPlotItem(pen=None,
                                            brush=pg.mkBrush(colour))
                marker.setZValue(25)
                self.image.view.addItem(marker)
                self.fit_markers[name] = marker
            centres = self.series.centres[:, index]
            good = np.isfinite(centres)
            sizes = 3.0 + 9.0 * np.sqrt(np.clip(heights[good, index] / top, 0, 1))
            xs = centres[good] if mdc else self.series.positions[good]
            ys = self.series.positions[good] if mdc else centres[good]
            marker.setData(xs, ys, size=sizes)

    @staticmethod
    def _number(edit):
        text = edit.text().strip()
        if not text:
            return None
        try:
            return float(text)
        except ValueError:
            return None

    def open_dispersion(self):
        if self.series is None:
            return None
        window = DispersionWindow(self.series, self.label, self)
        window.datasetsCreated.connect(self.datasetsCreated)
        window.show()
        self._dispersion = window
        return window

    def save_table(self):
        if self.series is None:
            return None
        path, _ = QFileDialog.getSaveFileName(
            self, "Save the fitted parameters", f"{self.label}_fits.csv",
            "CSV (*.csv)")
        if not path:
            return None
        header = ["position", "chi2"]
        for name in self.series.names:
            header += [f"{name}_centre", f"{name}_centre_err",
                       f"{name}_fwhm", f"{name}_fwhm_err",
                       f"{name}_height", f"{name}_area"]
        rows = []
        for i, position in enumerate(self.series.positions):
            row = [position, self.series.chi2[i]]
            for j in range(len(self.series.names)):
                row += [self.series.centres[i, j], self.series.centre_errors[i, j],
                        self.series.widths[i, j], self.series.width_errors[i, j],
                        self.series.heights[i, j], self.series.areas[i, j]]
            rows.append(row)
        np.savetxt(path, np.array(rows, dtype=float), delimiter=",",
                   header=",".join(header), comments="")
        self.note.setText(f"Saved {len(rows)} rows to {path}")
        return path

    # -- the guided fit ----------------------------------------------------
    def start_guide(self):
        """Open the step-by-step guide at its first step."""
        self.guide.start()
        return self.guide

    def lines_axis(self):
        """The axis the lines are spaced along (energy for MDCs)."""
        return self.axes[1] if self.direction.currentIndex() == 0 else self.axes[0]

    def move_line_to(self, position: float):
        """Put the fitted line (the cursor) at ``position`` and redraw."""
        axis = self.lines_axis()
        position = float(np.clip(position, np.nanmin(axis), np.nanmax(axis)))
        self.cursor.blockSignals(True)
        self.cursor.setPos(position)
        self.cursor.blockSignals(False)
        self.refresh_line()

    def fit_range(self):
        """``(low, high)`` the series will cover: the Over boxes, or the
        whole axis where they are empty."""
        axis = self.lines_axis()
        low, high = self._number(self.from_edit), self._number(self.to_edit)
        low = float(np.nanmin(axis)) if low is None else low
        high = float(np.nanmax(axis)) if high is None else high
        return (min(low, high), max(low, high))

    def closeEvent(self, event):
        self.closed.emit(self)
        super().closeEvent(event)


# ==========================================================================
# The guided fit
# ==========================================================================
class FitGuide(QGroupBox):
    """The fit, one step at a time, inside the fit panel.

    Each step says in two or three sentences what it is for, lights up the
    part of the panel it uses, and offers the buttons that do it; *Next*
    waits until the step has what the fit needs (two seeds on a band, a
    fitted series). Nothing here fits anything by itself: every button
    calls the panel's own, so a user who has been through the guide once
    knows the panel.
    """

    #: how many lines the seeding step visits, spread over the fit range
    SEED_LINES = 5
    HIGHLIGHT = "2px solid #f5a623"

    def __init__(self, panel: "FitPanel"):
        super().__init__("Guided fit")
        self.panel = panel
        self.setMinimumWidth(380)
        self.step = 0
        self._lit = []
        self._seed_index = 0
        self.setStyleSheet("QGroupBox { font-weight: bold; }")
        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 4, 8, 6)
        self.title = QLabel("")
        self.title.setStyleSheet("font-size: 11pt; color: #1f6f8b;")
        layout.addWidget(self.title)
        self.text = QLabel("")
        self.text.setWordWrap(True)
        self.text.setTextFormat(Qt.RichText)
        self.text.setStyleSheet("font-weight: normal;")
        layout.addWidget(self.text)
        self.status = QLabel("")
        self.status.setWordWrap(True)
        self.status.setTextFormat(Qt.RichText)
        self.status.setStyleSheet("font-weight: normal; color: #2f6b2f;")
        layout.addWidget(self.status)

        self.actions = QHBoxLayout()
        self.actions.setContentsMargins(0, 0, 0, 0)
        layout.addLayout(self.actions)

        nav = QHBoxLayout()
        nav.setContentsMargins(0, 0, 0, 0)
        self.progress = QLabel("")
        self.progress.setStyleSheet("font-weight: normal; color: #666;")
        nav.addWidget(self.progress)
        nav.addStretch(1)
        self.back_button = QPushButton("\u2190 Back")
        self.back_button.clicked.connect(lambda: self.go(self.step - 1))
        nav.addWidget(self.back_button)
        self.next_button = QPushButton("Next \u2192")
        self.next_button.setDefault(True)
        self.next_button.clicked.connect(lambda: self.go(self.step + 1))
        nav.addWidget(self.next_button)
        close = QPushButton("Close guide")
        close.clicked.connect(self.stop)
        nav.addWidget(close)
        layout.addLayout(nav)

        panel.seedPlaced.connect(self._seed_placed)
        panel.seriesFitted.connect(self._update)
        self.steps = self._steps()

    # -- the steps ---------------------------------------------------------
    def _steps(self):
        """``(title, text, enter, ready)`` per step: ``enter`` builds the
        step's buttons and lights its widgets, ``ready`` says whether Next
        may be pressed."""
        p = self.panel
        mdc = p.mode == "mdc"
        curve, across = ("MDC", "energy") if mdc else ("EDC", "momentum")
        along = p._short(0) if mdc else p._short(1)
        return [
            ("What the fit does",
             f"An <b>{curve}</b> is the intensity along <b>{along}</b> at one "
             f"{across} -- one horizontal line of the cut"
             + ("" if mdc else " turned on its side") +
             f". Where a band crosses that line the curve has a peak. "
             f"Fitting a peak to every line, one after the other, follows the "
             f"band down the cut. You tell the fit roughly where the band is "
             f"(the <b>seeds</b>); it does the rest.",
             self._enter_intro, lambda: True),
            ("Line shape and background",
             "Choose what one peak looks like. <b>Lorentzian</b> is the "
             "natural shape of a band in an MDC; <b>Voigt</b> adds the "
             "instrument's Gaussian resolution, so the width that comes out "
             "is the intrinsic one. The <b>background</b> is the slowly "
             "varying intensity under the peaks -- linear is a safe start. "
             "The defaults are fine for a first try.",
             self._enter_settings, lambda: True),
            ("Where to fit",
             f"The blue line on the cut is the {curve} being looked at; drag "
             f"it with the mouse. Put it where the band <b>starts</b> and "
             f"press <i>Start here</i>, then where it <b>ends</b> and press "
             f"<i>End here</i>. Outside that range nothing is fitted -- leave "
             f"out the region where the band is too weak or crosses another.",
             self._enter_range, lambda: True),
            ("The bands",
             "Each <b>band</b> is one peak followed from line to line. "
             "<b>Band 1</b> is already there. If two bands cross the range "
             "(two peaks in one curve), add a second band now: they are "
             "fitted together, so each keeps to its own peak.",
             self._enter_bands, lambda: bool(p.bands)),
            ("Seeds: show the fit where the band is",
             "The guide moves the line through the range. At each stop, "
             "<b>click the top of the band's peak on the curve at the "
             "right</b> -- that click is a seed: the position is yours, the "
             "height and width are read off the data. Between seeds the "
             "band is interpolated, so 4-5 seeds spread along it are enough.",
             self._enter_seeds,
             lambda: any(len(b.seeds) >= 2 for b in p.bands)),
            ("Check one line",
             "Before fitting everything, try the line under the cursor. The "
             "dashed curves are the starting guess built from the seeds; "
             "<i>Try this line</i> fits it and draws the result with its "
             "residual below. A residual without structure and peaks on the "
             "data mean the settings are good; if not, go back and move a "
             "seed or change the line shape.",
             self._enter_try, lambda: True),
            ("Fit the series",
             f"Now every {curve} in the range is fitted, each starting from "
             f"the seeds' interpolation (and, on the second pass, from the "
             f"fitted neighbours). It takes a few seconds.",
             self._enter_fit, lambda: p.series is not None),
            ("The result",
             "The fitted peak positions are drawn on the cut, <b>sized by "
             "amplitude</b>: where the markers shrink, the band fades and "
             "the fit is held up by few counts. "
             + ("<i>Dispersion...</i> turns the band into v_F and m*; "
                if p.momentum else
                "v_F and m* need a cut in momentum, which this is not; ")
             + "<i>Save table...</i> writes every fitted parameter to a CSV "
             "file. Seeds and settings stay as they are: change one and fit "
             "again at any time.",
             self._enter_result, lambda: True),
        ]

    # -- navigation --------------------------------------------------------
    def start(self):
        self.show()
        self.panel.guide_button.setEnabled(False)
        self.go(0)

    def stop(self):
        self._unlight()
        self._clear_actions()
        self.panel.seed_button.setChecked(False)
        self.hide()
        self.panel.guide_button.setEnabled(True)

    def go(self, index: int):
        if index >= len(self.steps):
            self.stop()
            return
        self.step = max(0, index)
        self._unlight()
        self._clear_actions()
        self.status.setText("")
        title, text, enter, _ready = self.steps[self.step]
        self.title.setText(f"Step {self.step + 1} of {len(self.steps)}: {title}")
        self.text.setText(text)
        enter()
        self._update()

    def _update(self, *_):
        ready = self.steps[self.step][3]()
        last = self.step == len(self.steps) - 1
        self.next_button.setText("Finish" if last else "Next \u2192")
        self.next_button.setEnabled(bool(ready))
        self.back_button.setEnabled(self.step > 0)
        self.progress.setText("\u25cf " * (self.step + 1)
                              + "\u25cb " * (len(self.steps) - self.step - 1))
        if self.step == 6 and self.panel.series is not None:
            self.status.setText("Fitted -- press Next to read the result.")

    def _button(self, text, slot, tip=""):
        button = QPushButton(text)
        button.setStyleSheet("font-weight: normal;")
        if tip:
            button.setToolTip(tip)
        button.clicked.connect(slot)
        self.actions.addWidget(button)
        return button

    def _clear_actions(self):
        while self.actions.count():
            item = self.actions.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.hide()
                widget.setParent(None)
                widget.deleteLater()

    def _light(self, *widgets):
        """Outline the widgets a step uses, in the guide's colour."""
        for widget in widgets:
            if widget is None:
                continue
            name = widget.objectName() or f"guide_{id(widget)}"
            widget.setObjectName(name)
            self._lit.append((widget, widget.styleSheet()))
            kind = type(widget).__name__
            extra = ""
            if isinstance(widget, QGroupBox):
                # A group box's border carries its title: keep room for it.
                kind = "QGroupBox"
                extra = (f" margin-top: 1.1em; }} QGroupBox#{name}::title {{ "
                         f"subcontrol-origin: margin; left: 6px; "
                         f"color: #b86e00;")
            elif isinstance(widget, (QTabWidget, QTableWidget)):
                kind = "QWidget"
            widget.setStyleSheet(widget.styleSheet()
                                 + f" {kind}#{name} {{ border: {self.HIGHLIGHT}; "
                                   f"border-radius: 3px;{extra} }}")

    def _unlight(self):
        for widget, sheet in self._lit:
            try:
                widget.setStyleSheet(sheet)
            except RuntimeError:          # deleted with its window
                pass
        self._lit = []

    # -- what each step does when it opens ----------------------------------
    def _enter_intro(self):
        self._light(self.panel.cut_box, self.panel.line_box)

    def _enter_settings(self):
        p = self.panel
        p.tabs.setCurrentIndex(1)
        self._light(p.shape, p.background, p.resolution)
        if p.mode == "edc":
            self.text.setText(self.text.text() + (
                " For an EDC near E<sub>F</sub>, tick <b>Fermi cut-off</b> so "
                "the peak is fitted under the occupation."))
            self._light(p.fermi_group)

    def _enter_range(self):
        p = self.panel
        p.tabs.setCurrentIndex(1)
        self._light(p.cut_box, p.from_edit, p.to_edit)

        def take(edit, which):
            edit.setText(f"{p.current_position():.6g}")
            low, high = p.fit_range()
            self.status.setText(f"{which} set. The fit covers "
                                f"{low:.5g} .. {high:.5g}.")
            self._seed_index = 0

        def whole():
            p.from_edit.clear()
            p.to_edit.clear()
            low, high = p.fit_range()
            self.status.setText(f"The whole cut: {low:.5g} .. {high:.5g}.")
            self._seed_index = 0

        self._button("Start here", lambda: take(p.from_edit, "Start"),
                     "The line's position becomes the first line fitted.")
        self._button("End here", lambda: take(p.to_edit, "End"),
                     "The line's position becomes the last line fitted.")
        self._button("Whole cut", whole, "Fit every line of the cut.")
        self.actions.addStretch(1)
        low, high = p.fit_range()
        self.status.setText(f"Now: {low:.5g} .. {high:.5g}.")

    def _enter_bands(self):
        p = self.panel
        p.tabs.setCurrentIndex(0)
        self._light(p.band_table)

        def another():
            band = p.add_band()
            self.status.setText(f"Added {band.name}. {len(p.bands)} band(s).")
            self._update()

        self._button("Add another band", another,
                     "Only when a second peak runs through the range.")
        self.actions.addStretch(1)
        self.status.setText(f"{len(p.bands)} band(s): "
                            + ", ".join(b.name for b in p.bands))

    def seed_positions(self):
        """The lines the seeding step visits, spread over the fit range and
        kept a little inside its ends, where a band is usually weakest."""
        low, high = self.panel.fit_range()
        n = self.SEED_LINES
        margin = 0.06 * (high - low)
        return list(np.linspace(low + margin, high - margin, n))

    def _enter_seeds(self):
        p = self.panel
        p.tabs.setCurrentIndex(0)
        self._light(p.line_box, p.band_table)
        if p.selected_band() is None and p.bands:
            p.band_table.selectRow(0)
        self._button("\u2190 Previous line", lambda: self._visit(self._seed_index - 1))
        self._button("Skip this line \u2192", lambda: self._visit(self._seed_index + 1),
                     "The band is not visible here: go to the next line.")
        self._button("Undo last seed", self._undo)
        self.band_choice = QComboBox()
        self.band_choice.setStyleSheet("font-weight: normal;")
        for band in p.bands:
            self.band_choice.addItem(band.name)
        index = p.selected_band_index()
        self.band_choice.setCurrentIndex(0 if index is None else index)
        self.band_choice.currentIndexChanged.connect(self._band_chosen)
        if len(p.bands) > 1:
            self.actions.addWidget(QLabel("Seeding"))
            self.actions.addWidget(self.band_choice)
        else:
            self.band_choice.hide()
            self.actions.addWidget(self.band_choice)
        self.actions.addStretch(1)
        self._visit(0)

    def _band_chosen(self, index):
        if 0 <= index < len(self.panel.bands):
            self.panel.band_table.selectRow(index)
            self._visit(0)

    def _visit(self, index: int):
        p = self.panel
        positions = self.seed_positions()
        band = p.selected_band()
        if index >= len(positions):
            # Round done for this band: the next band, if one has fewer seeds.
            others = [i for i, b in enumerate(p.bands)
                      if len(b.seeds) < 2 and b is not band]
            if others:
                self.band_choice.setCurrentIndex(others[0])
                return
            p.seed_button.setChecked(False)
            self.status.setText(
                f"Every line visited: "
                + ", ".join(f"{b.name} {len(b.seeds)} seed(s)" for b in p.bands)
                + ". Press Next -- or Previous line to change one.")
            self._seed_index = len(positions)
            self._update()
            return
        self._seed_index = max(0, index)
        p.move_line_to(positions[self._seed_index])
        p.seed_button.setChecked(True)
        name = band.name if band is not None else "the band"
        self.status.setText(
            f"Line {self._seed_index + 1} of {len(positions)} "
            f"({p.position_label.text()} = {positions[self._seed_index]:.5g}): "
            f"click the peak of <b>{name}</b> on the curve at the right. "
            f"{name}: {len(band.seeds) if band else 0} seed(s).")
        self._update()

    def _seed_placed(self):
        if not self.isVisible() or self.step != 4:
            return
        self._visit(self._seed_index + 1)

    def _undo(self):
        self.panel.undo_seed()
        self._visit(max(0, self._seed_index - 1))

    def _enter_try(self):
        p = self.panel
        p.seed_button.setChecked(False)
        self._light(p.line_box, p.try_button)
        positions = self.seed_positions()
        p.move_line_to(positions[len(positions) // 2])

        def attempt():
            fit = p.try_line()
            if fit is not None:
                self.status.setText(f"Fitted this line: χ² = {fit.chi2:.3g}, "
                                    f"R² = {fit.r_squared:.4f}.")

        self._button("Try this line", attempt)
        self.actions.addStretch(1)
        self.status.setText("Drag the blue line anywhere in the range to try "
                            "another one.")

    def _enter_fit(self):
        p = self.panel
        p.tabs.setCurrentIndex(1)
        self._light(p.run_button)
        self._button("Fit the series now", p.run_fit)
        self.actions.addStretch(1)

    def _enter_result(self):
        p = self.panel
        p.tabs.setCurrentIndex(1)
        self._light(p.cut_box, p.dispersion_button, p.table_button)
        if p.momentum:
            self._button("Dispersion...", p.open_dispersion)
        self._button("Save table...", p.save_table)
        self.actions.addStretch(1)


# ==========================================================================
# The dispersion
# ==========================================================================
class DispersionWindow(QMainWindow):
    """The fitted band, what it says about v_F and m*, and how much that
    depends on the window it was measured over."""

    datasetsCreated = pyqtSignal(list)

    def __init__(self, series: PK.BandSeries, label: str, parent=None):
        super().__init__(parent)
        self.series = series
        self.label = label
        self.fit = None
        self.setWindowTitle(f"Dispersion — {label}")
        self.resize(1080, 700)

        central = QWidget()
        root = QVBoxLayout(central)
        split = QSplitter(Qt.Horizontal)

        left = QWidget()
        column = QVBoxLayout(left)
        column.setContentsMargins(0, 0, 0, 0)
        self.band_plot = pg.PlotWidget()
        self.band_plot.setLabel("bottom", series.x_label)
        self.band_plot.setLabel("left", series.y_label)
        strip_stock_menu(self.band_plot.getPlotItem())
        column.addWidget(self.band_plot)
        split.addWidget(left)

        right = QWidget()
        column2 = QVBoxLayout(right)
        column2.setContentsMargins(0, 0, 0, 0)
        self.scan_plot = pg.PlotWidget()
        self.scan_plot.setLabel("bottom", "half-width of the fitting window")
        self.scan_plot.setLabel("left", "value")
        strip_stock_menu(self.scan_plot.getPlotItem())
        column2.addWidget(self.scan_plot)
        split.addWidget(right)
        split.setSizes([540, 540])
        root.addWidget(split, 1)

        form = QHBoxLayout()
        self.band_choice = QComboBox()
        self.band_choice.addItems(series.names)
        self.band_choice.currentIndexChanged.connect(self.refresh)
        form.addWidget(QLabel("Band"))
        form.addWidget(self.band_choice)
        self.quantity = QComboBox()
        self.quantity.addItems(["v_F (straight line)", "m* (parabola)"])
        self.quantity.currentIndexChanged.connect(self.refresh)
        form.addWidget(QLabel("Measure"))
        form.addWidget(self.quantity)
        form.addWidget(QLabel("about"))
        self.centre = spin(0.0, -1e4, 1e4, 0.005, 5)
        self.centre.setToolTip(
            "E_F for a velocity, the band extremum for a mass. Both are "
            "local quantities, so this is where they are measured.")
        self.centre.valueChanged.connect(self.refresh)
        form.addWidget(self.centre)
        find = QPushButton("Guess")
        find.setToolTip("E_F = 0 for a velocity; the turning point of a "
                        "parabola fitted to everything, for a mass.")
        find.clicked.connect(self.guess_centre)
        form.addWidget(find)
        form.addStretch(1)
        root.addLayout(form)

        self.summary = QLabel("")
        self.summary.setTextInteractionFlags(Qt.TextSelectableByMouse)
        self.summary.setWordWrap(True)
        root.addWidget(self.summary)

        buttons = QHBoxLayout()
        for text, slot, tip in (
                ("Band to list", self.export_band,
                 "Add the fitted positions and widths to the main list."),
                ("Self-energy...", self.open_self_energy,
                 "Re Sigma and Im Sigma against a bare band, with the "
                 "Kramers-Kronig check.")):
            button = QPushButton(text)
            button.setToolTip(tip)
            button.clicked.connect(slot)
            buttons.addWidget(button)
        buttons.addStretch(1)
        root.addLayout(buttons)
        self.setCentralWidget(central)
        self.guess_centre()

    def _band_data(self):
        index = self.band_choice.currentIndex()
        k, energy, error = self.series.band(index)
        on_k = self.series.error_axis() == "x"
        return k, energy, (error if on_k else None), (None if on_k else error)

    def guess_centre(self):
        k, energy, k_error, energy_error = self._band_data()
        if self.quantity.currentIndex() == 0:
            self.centre.setValue(0.0)
        else:
            try:
                rough = DISP.fit_dispersion(k, energy, k_error=k_error,
                                            energy_error=energy_error, order=2)
                extremum = rough.band_extremum()[0]
                self.centre.setValue(float(extremum)
                                     if np.isfinite(extremum) else float(np.mean(k)))
            except Exception:
                self.centre.setValue(float(np.mean(k)))
        self.refresh()

    def refresh(self):
        k, energy, k_error, energy_error = self._band_data()
        self.band_plot.clear()
        self.scan_plot.clear()
        if k.size < 3:
            self.summary.setText("not enough converged fits for this band")
            return

        self.band_plot.plot(k, energy, pen=None, symbol="o", symbolSize=5,
                            symbolBrush="#bd4921", symbolPen=None)
        error = k_error if k_error is not None else energy_error
        bars = (pg.ErrorBarItem(x=k, y=energy,
                                width=2 * np.nan_to_num(error, nan=0.0))
                if k_error is not None else
                pg.ErrorBarItem(x=k, y=energy,
                                height=2 * np.nan_to_num(error, nan=0.0)))
        bars.setOpts(pen=pg.mkPen("#bd4921"))
        self.band_plot.addItem(bars)

        velocity = self.quantity.currentIndex() == 0
        order = 1 if velocity else 2
        try:
            self.fit = DISP.fit_dispersion(k, energy, k_error=k_error,
                                           energy_error=energy_error, order=order)
            grid = np.linspace(k.min(), k.max(), 200)
            self.band_plot.plot(grid, self.fit.value(grid),
                                pen=pg.mkPen("#0d7377", width=2))
        except Exception as exc:
            self.summary.setText(f"<span style='color:#b00'>{exc}</span>")
            return

        scan = DISP.window_scan(
            k, energy, k_error=k_error, energy_error=energy_error,
            quantity="velocity" if velocity else "mass",
            centre=self.centre.value())
        good = np.isfinite(scan.values)
        if good.any():
            self.scan_plot.plot(scan.widths[good], scan.values[good],
                                pen=pg.mkPen("#0d7377", width=2), symbol="o",
                                symbolSize=5, symbolBrush="#0d7377",
                                symbolPen=None)
            self.scan_plot.addItem(pg.ErrorBarItem(
                x=scan.widths[good], y=scan.values[good],
                height=2 * np.nan_to_num(scan.errors[good], nan=0.0),
                pen=pg.mkPen("#0d7377")))
        self.scan_plot.setLabel(
            "left", "dE/dk (eV·Å)" if velocity else "m*/mₑ")
        self.scan_plot.setLabel(
            "bottom", f"window half-width ({'eV' if velocity else 'Å⁻¹'})")

        plateau = scan.plateau_value()
        lines = [self.fit.summary().replace("\n", "<br>")]
        if plateau is None:
            lines.append(
                "<b>No plateau.</b> The answer moves with the window from the "
                "narrowest one onwards, so there is no window-independent "
                "value to quote here — which is itself a result about the "
                "band, not a failure of the fit.")
        else:
            value, error, widest = plateau
            first, last = scan.plateau(0.05)
            self.scan_plot.addItem(pg.LinearRegionItem(
                values=(scan.widths[first], scan.widths[last]), movable=False,
                brush=pg.mkBrush(13, 115, 119, 40)))
            unit = "eV·Å" if velocity else ""
            extra = (f" = {value * DISP.VELOCITY_FACTOR:.3g} m/s"
                     if velocity else "")
            lines.append(
                f"<b>Stable out to a window of {widest:.4g}"
                f"{' eV' if velocity else ' Å⁻¹'}</b>: "
                f"{value:.4g} ± {error:.2g} {unit}{extra}")
        lines.append(f"<i>Weighted by the error on "
                     f"{'k' if k_error is not None else 'E'}, which is the "
                     f"axis the fit actually measured.</i>")
        self.summary.setText("<br>".join(lines))

    def export_band(self):
        index = self.band_choice.currentIndex()
        k, energy, error = self.series.band(index)
        if k.size == 0:
            return None
        widths = self.series.widths[:, index]
        widths = widths[np.isfinite(self.series.centres[:, index])]
        values = np.column_stack([energy, error, widths[:k.size]])
        data = MemoryData(
            "cut", (k, np.array([0.0, 1.0, 2.0])), values,
            {"x": self.series.x_label, "y": "energy / error / FWHM"},
            source_label=f"{self.label}_{self.series.names[index]}",
            parameters={"band": self.series.names[index],
                        "direction": self.series.direction},
            prefix="proc.band")
        self.datasetsCreated.emit([data])
        return data

    def open_self_energy(self):
        if self.fit is None:
            return None
        if self.series.direction != "mdc":
            QMessageBox.information(
                self, "Self-energy",
                "A self-energy comes from MDC widths. Refit as MDCs.")
            return None
        try:
            result = DISP.self_energy(self.series, self.fit,
                                      band=self.band_choice.currentIndex())
        except Exception as exc:
            QMessageBox.warning(self, "Self-energy", str(exc))
            return None
        window = SelfEnergyWindow(result, self.label, self)
        window.show()
        self._sigma = window
        return window


class SelfEnergyWindow(QMainWindow):
    """``Re Sigma`` and ``Im Sigma``, with the Kramers-Kronig check."""

    def __init__(self, result: DISP.SelfEnergy, label: str, parent=None):
        super().__init__(parent)
        self.result = result
        self.setWindowTitle(f"Self-energy — {label}")
        self.resize(760, 560)
        central = QWidget()
        layout = QVBoxLayout(central)

        real = pg.PlotWidget()
        real.setLabel("left", "Re Σ (eV)")
        real.setLabel("bottom", "E − E_F (eV)")
        real.plot(result.energy, result.real, pen=None, symbol="o", symbolSize=4,
                  symbolBrush="#bd4921", symbolPen=None)
        real.addItem(pg.ErrorBarItem(
            x=result.energy, y=result.real,
            height=2 * np.nan_to_num(result.real_error), pen=pg.mkPen("#bd4921")))
        if result.kk_real is not None:
            real.plot(result.energy, result.kk_real,
                      pen=pg.mkPen("#0d7377", width=2, style=Qt.DashLine))
        strip_stock_menu(real.getPlotItem())
        layout.addWidget(real)

        imaginary = pg.PlotWidget()
        imaginary.setLabel("left", "Im Σ (eV)")
        imaginary.setLabel("bottom", "E − E_F (eV)")
        imaginary.setXLink(real)
        imaginary.plot(result.energy, result.imaginary, pen=None, symbol="o",
                       symbolSize=4, symbolBrush="#0d7377", symbolPen=None)
        imaginary.addItem(pg.ErrorBarItem(
            x=result.energy, y=result.imaginary,
            height=2 * np.nan_to_num(result.imaginary_error),
            pen=pg.mkPen("#0d7377")))
        strip_stock_menu(imaginary.getPlotItem())
        layout.addWidget(imaginary)

        consistency = result.consistency()
        note = QLabel(
            f"Dashed: the Kramers-Kronig transform of Im Σ, shifted to the "
            f"mean of Re Σ. Mismatch {consistency:.2f} "
            f"({'consistent' if consistency < 0.3 else 'inconsistent'}) — the "
            f"two are not independent, so a bare band chosen to enlarge a "
            f"kink fails this.")
        note.setWordWrap(True)
        layout.addWidget(note)
        self.setCentralWidget(central)
