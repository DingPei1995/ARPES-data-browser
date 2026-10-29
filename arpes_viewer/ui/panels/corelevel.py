"""
ui/panels/corelevel.py
======================
The core-level check of a cut, from the main list's right-click menu: what
element and orbital each peak of the angle-integrated EDC is.

Three steps, each the user can correct before the next one runs:

1. :class:`ElementsDialog` -- which elements the sample is made of (a
   formula is fine: ``NbSCl``), and how strict to be.
2. :class:`EnergyReferenceDialog` -- the photon energy, the work function
   and the Fermi level, found from the data where the window holds the
   edge and otherwise put at hv - W. These are what the binding energies
   hang on, so they are shown and editable rather than silently applied.
3. :class:`CoreLevelWindow` -- the EDC on the binding-energy scale with
   every peak marked, a table of what each one is (or might be), and what
   is missing. The assignment column is the user's to overwrite. The count
   axis can be logarithmic, so weak lines next to a strong one show; a peak
   the search missed can be added by hand (click the plot, or type its
   energy) and a false one removed; an added peak is named with the levels
   within ±2 eV. Two helpers open from it:

   * :class:`PeakTableDialog` -- every tabulated level, filtered by element
     and binding-energy range, to look a peak up by eye; a level picked
     there is shown on the plot and can be added as a peak.
   * :class:`ZoomWindow` -- a range of the EDC in a window of its own, with
     a cursor, a local fit and the second derivative, to place a peak.

The analysis itself is :mod:`tools.spectra.corelevel`, with no Qt in it.
"""
from __future__ import annotations

import json
import os

import numpy as np
import pyqtgraph as pg
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QColor, QGuiApplication
from PyQt5.QtWidgets import (QDialog, QMainWindow, QWidget, QVBoxLayout,
                             QHBoxLayout, QFormLayout, QLabel, QLineEdit,
                             QPushButton, QCheckBox, QComboBox, QAbstractItemView,
                             QDialogButtonBox, QTableWidget, QTableWidgetItem,
                             QHeaderView, QSplitter, QPlainTextEdit,
                             QMessageBox, QGroupBox)

import arpes_viewer.tools.spectra.corelevel as CL
from arpes_viewer.ui.panels.process import spin
from arpes_viewer.ui.viewers.widgets import strip_stock_menu

#: Plot symbol of a found peak and of one added by hand.
FOUND_SYMBOL, MANUAL_SYMBOL = "o", "d"

#: The level picked in the peak table, drawn on the plots.
PREVIEW_COLOUR = "#b0228c"

#: Marker colour of each status, on the plot and in the table.
STATUS_COLOURS = {
    "sample": "#2f7a3f",
    "ambiguous": "#1f6f8b",
    "uncertain": "#c47f00",
    "unexplained": "#bd4921",
}

LAST_USED_NAME = "corelevel.json"


def plotted_counts(values, intensity, log: bool):
    """Intensities as plotted: log10 on a log axis, where what is not
    positive is put just below the smallest positive count of
    ``intensity``."""
    values = np.asarray(values, dtype=float)
    if not log:
        return values
    intensity = np.asarray(intensity, dtype=float)
    positive = intensity[np.isfinite(intensity) & (intensity > 0)]
    floor = 0.5 * float(positive.min()) if positive.size else 1e-3
    with np.errstate(invalid="ignore"):
        return np.log10(np.where(values > floor, values, floor))


def _last_used_path() -> str:
    folder = os.environ.get("ARPES_VIEWER_CONFIG_DIR") or os.path.join(
        os.path.expanduser("~"), ".arpes_viewer")
    return os.path.join(folder, LAST_USED_NAME)


def load_last_used() -> dict:
    """The settings the check was last run with, or ``{}``."""
    try:
        with open(_last_used_path(), encoding="utf-8") as handle:
            settings = json.load(handle)
        return settings if isinstance(settings, dict) else {}
    except (OSError, ValueError):
        return {}


def save_last_used(settings: dict) -> None:
    path = _last_used_path()
    try:
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w", encoding="utf-8") as handle:
            json.dump(settings, handle, indent=1)
    except OSError:
        pass


# --------------------------------------------------------------------------
# Step 1: the elements
# --------------------------------------------------------------------------
class ElementsDialog(QDialog):
    """Ask what the sample is made of, and how strict the comparison is."""

    def __init__(self, parent=None, *, elements: str = None, settings=None):
        super().__init__(parent)
        self.setWindowTitle("Core-level check: elements")
        settings = dict(load_last_used() if settings is None else settings)
        layout = QVBoxLayout(self)
        intro = QLabel(
            "Which elements does the sample contain? Symbols in any order "
            "(Nb, S, Cl) or a formula (NbSCl, Bi2Se3). Each peak of the "
            "angle-integrated EDC is compared with their tabulated binding "
            "energies (X-ray Data Booklet, Table 1-1).")
        intro.setWordWrap(True)
        layout.addWidget(intro)

        form = QFormLayout()
        self.elements_box = QLineEdit(elements if elements is not None
                                      else settings.get("elements", ""))
        self.elements_box.setPlaceholderText("e.g. NbSCl  or  Nb S Cl")
        self.elements_box.textChanged.connect(self._update_parsed)
        form.addRow("Elements", self.elements_box)
        self.parsed = QLabel("")
        self.parsed.setWordWrap(True)
        form.addRow("", self.parsed)
        self.tolerance = spin(float(settings.get("tolerance", CL.TOLERANCE_EV)),
                              0.05, 5.0, 0.05, 2, " eV")
        self.tolerance.setToolTip(
            "A peak this close to a tabulated level of the sample is taken "
            "as that level.")
        form.addRow("Match within", self.tolerance)
        self.uncertain = spin(float(settings.get("uncertain", CL.UNCERTAIN_EV)),
                              0.1, 10.0, 0.1, 2, " eV")
        self.uncertain.setToolTip(
            "Further than the tolerance but within this, a peak may still be "
            "the sample's, chemically shifted or charged: it is marked "
            "'uncertain' and left to you.")
        form.addRow("Uncertain within", self.uncertain)
        self.chemical = QCheckBox("Compare typical chemical-state values too")
        self.chemical.setChecked(bool(settings.get("chemical_states", True)))
        self.chemical.setToolTip(
            "Oxides, chlorides, sulfides, adventitious carbon... The table "
            "gives O, N, F and Cl as free gases against the vacuum level, "
            "several eV from where they sit in a solid.")
        form.addRow("", self.chemical)
        layout.addLayout(form)

        self.buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        self.buttons.button(QDialogButtonBox.Ok).setText("Next")
        self.buttons.accepted.connect(self._accept)
        self.buttons.rejected.connect(self.reject)
        layout.addWidget(self.buttons)
        self._update_parsed()

    def elements(self):
        return CL.parse_elements(self.elements_box.text())

    def settings(self) -> dict:
        return {"elements": self.elements_box.text().strip(),
                "tolerance": self.tolerance.value(),
                "uncertain": self.uncertain.value(),
                "chemical_states": self.chemical.isChecked()}

    def _update_parsed(self):
        elements, unknown = self.elements()
        text = ("Read as: " + ", ".join(elements)) if elements else \
            "No element yet."
        if unknown:
            text += "  —  not an element: " + ", ".join(unknown)
        self.parsed.setText(text)
        self.buttons.button(QDialogButtonBox.Ok).setEnabled(bool(elements)
                                                            and not unknown)

    def _accept(self):
        elements, unknown = self.elements()
        if not elements or unknown:
            return
        save_last_used({**load_last_used(), **self.settings()})
        self.accept()


# --------------------------------------------------------------------------
# Step 2: the energy reference
# --------------------------------------------------------------------------
class EnergyReferenceDialog(QDialog):
    """Show, and let the user change, what the binding energies hang on."""

    def __init__(self, energy, intensity, energy_label: str, info: dict,
                 parent=None, *, temperature: float = 30.0, reference=None):
        super().__init__(parent)
        self.setWindowTitle("Core-level check: energy reference")
        self.energy = np.asarray(energy, dtype=float)
        self.intensity = np.asarray(intensity, dtype=float)
        self.temperature = float(temperature)
        self.estimate = None
        hv = CL.photon_energy_of(info)
        work = CL.work_function_of(info)
        scale = CL.detect_scale(self.energy, energy_label, hv)
        if reference is not None:
            hv, work, scale = reference.photon_energy, reference.work_function, \
                reference.scale

        layout = QVBoxLayout(self)
        intro = QLabel(
            "The binding energies are measured from the Fermi level. E_F was "
            "looked for in the integrated EDC; where the window does not hold "
            "the edge it is put at hν − W. Check the numbers, change any that "
            "are wrong, then press Analyse.")
        intro.setWordWrap(True)
        layout.addWidget(intro)

        form = QFormLayout()
        self.scale_box = QComboBox()
        for name in CL.SCALES:
            self.scale_box.addItem(CL.SCALE_LABELS[name], name)
        self.scale_box.setCurrentIndex(CL.SCALES.index(scale))
        self.scale_box.setToolTip(f"The axis is labelled “{energy_label}”.")
        form.addRow("Energy axis is", self.scale_box)
        self.hv_box = spin(hv or 0.0, 0.0, 20000.0, 0.1, 3, " eV")
        self.hv_box.setSpecialValueText("unknown")
        self.hv_box.setToolTip(
            "From the file. Needed to place E_F when the window does not "
            "hold the edge, and to recognise higher-order light. 0 = unknown.")
        form.addRow("Photon energy hν", self.hv_box)
        self.work_box = spin(work, 0.0, 10.0, 0.05, 3, " eV")
        self.work_box.setToolTip("The analyser work function W: E_F (kinetic) "
                                 "= hν − W.")
        form.addRow("Work function W", self.work_box)
        self.ef_box = spin(0.0, -1e6, 1e6, 0.01, 4, " eV")
        self.ef_box.setToolTip("Where the Fermi level is on the data's own "
                               "energy axis.")
        form.addRow("E_F on the axis", self.ef_box)
        layout.addLayout(form)

        row = QHBoxLayout()
        self.find_button = QPushButton("Find E_F in the data")
        self.find_button.clicked.connect(self.find_edge)
        self.expected_button = QPushButton("E_F = hν − W")
        self.expected_button.setToolTip("Put E_F where the photon energy and "
                                        "the work function say it is.")
        self.expected_button.clicked.connect(self.use_expected)
        row.addWidget(self.find_button)
        row.addWidget(self.expected_button)
        row.addStretch(1)
        layout.addLayout(row)

        self.message = QLabel("")
        self.message.setWordWrap(True)
        self.message.setTextInteractionFlags(Qt.TextSelectableByMouse)
        layout.addWidget(self.message)

        self.plot = pg.PlotWidget()
        self.plot.setMinimumHeight(200)
        self.plot.setLabel("bottom", energy_label)
        self.plot.setLabel("left", "counts (angle-integrated)")
        strip_stock_menu(self.plot.getPlotItem())
        self.curve = self.plot.plot(self.energy, self.intensity,
                                    pen=pg.mkPen("#1f6f8b", width=1.5))
        self.model = self.plot.plot(pen=pg.mkPen("#bd4921", width=2))
        self.ef_line = pg.InfiniteLine(angle=90, movable=True,
                                       pen=pg.mkPen("#2f7a3f", width=2,
                                                    style=Qt.DashLine))
        self.ef_line.setToolTip("E_F: drag it, or type the value above.")
        self.plot.addItem(self.ef_line, ignoreBounds=True)
        self.ef_line.sigPositionChanged.connect(
            lambda line: self.ef_box.setValue(float(line.value())))
        self.ef_box.valueChanged.connect(self._ef_changed)
        layout.addWidget(self.plot, stretch=1)

        self.check_label = QLabel("")
        self.check_label.setWordWrap(True)
        layout.addWidget(self.check_label)
        self.hv_box.valueChanged.connect(self._update_check)
        self.work_box.valueChanged.connect(self._update_check)
        self.scale_box.currentIndexChanged.connect(lambda *_: self.find_edge())

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.button(QDialogButtonBox.Ok).setText("Analyse")
        buttons.accepted.connect(self._accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)
        self.resize(560, 620)

        if reference is not None:
            self.ef_box.setValue(reference.fermi_level)
            self.message.setText("The reference the last analysis used.")
            self._zoom()
        else:
            self.find_edge()

    # -- values ------------------------------------------------------------
    def scale(self) -> str:
        return self.scale_box.currentData()

    def photon_energy(self):
        value = self.hv_box.value()
        return value if value > 0 else None

    def reference(self) -> CL.EnergyReference:
        return CL.EnergyReference(self.scale(), self.ef_box.value(),
                                  self.photon_energy(), self.work_box.value())

    # -- actions -----------------------------------------------------------
    def find_edge(self):
        self.estimate = CL.locate_fermi_edge(
            self.energy, self.intensity, self.scale(),
            photon_energy=self.photon_energy(),
            work_function=self.work_box.value(),
            temperature=self.temperature)
        level = self.estimate.level
        if level is None or not np.isfinite(level):
            level = float(self.energy[-1]) if self.scale() == "kinetic" else 0.0
        self.ef_box.setValue(level)
        fit = self.estimate.fit
        if fit is not None and self.scale() != "binding":
            lo, hi = fit.window
            e = self.energy[(self.energy >= lo) & (self.energy <= hi)]
            self.model.setData(e, fit.model(e))
        else:
            self.model.setData([], [])
        prefix = "" if self.estimate.found else "⚠ "
        self.message.setText(prefix + self.estimate.message)
        self._zoom()
        return self.estimate

    def use_expected(self):
        expected = CL.expected_fermi_level(self.scale(), self.photon_energy(),
                                           self.work_box.value())
        if expected is None:
            QMessageBox.information(self, "Energy reference",
                                    "Type the photon energy first.")
            return
        self.ef_box.setValue(expected)
        self.message.setText(f"E_F put at {expected:.3f} eV "
                             + ("(hν − W)." if self.scale() == "kinetic"
                                else "(the axis is already referred to E_F)."))
        self._zoom()

    def _ef_changed(self, value):
        if abs(self.ef_line.value() - value) > 1e-9:
            self.ef_line.blockSignals(True)
            self.ef_line.setValue(value)
            self.ef_line.blockSignals(False)
        self._update_check()

    def _zoom(self):
        """Show the edge when it is in the window, the whole EDC when not."""
        level = self.ef_box.value()
        lo, hi = float(self.energy.min()), float(self.energy.max())
        if lo - 1 <= level <= hi + 1 and (hi - lo) > 8:
            self.plot.setXRange(level - 3, level + 1.5, padding=0)
            near = (self.energy >= level - 3) & (self.energy <= level + 1.5)
            if near.any():
                values = self.intensity[near]
                self.plot.setYRange(float(np.nanmin(values)),
                                    float(np.nanmax(values)), padding=0.05)
        else:
            self.plot.autoRange()

    def _update_check(self, *_):
        reference = self.reference()
        texts = []
        implied = reference.implied_work_function()
        if implied is not None:
            texts.append(f"hν − E_F = {implied:.3f} eV (should be the work "
                         f"function, typically 4–5 eV).")
            if not 3.5 <= implied <= 5.5:
                texts.append("⚠ That is not a plausible work function: the "
                             "photon energy or E_F is off.")
        b = reference.binding(self.energy)
        texts.append(f"Binding-energy window: {np.nanmin(b):.2f} to "
                     f"{np.nanmax(b):.2f} eV.")
        if reference.max_binding is not None and np.nanmax(b) > reference.max_binding + 1:
            texts.append(f"The window goes deeper than hν − W = "
                         f"{reference.max_binding:.1f} eV, which only "
                         f"higher-order light reaches.")
        self.check_label.setText(" ".join(texts))

    def _accept(self):
        if not np.isfinite(self.ef_box.value()):
            return
        self.accept()


# --------------------------------------------------------------------------
# Step 3: the results
# --------------------------------------------------------------------------
class CoreLevelWindow(QMainWindow):
    """The EDC on the binding-energy scale, its peaks and what they are."""

    datasetsCreated = pyqtSignal(list)
    closed = pyqtSignal(object)
    #: after every analysis, so the zoom windows follow
    analysed = pyqtSignal()
    #: (binding energy or None, label): the level picked in the peak table
    previewChanged = pyqtSignal(object, str)

    COLUMNS = ("#", "E_B (eV)", "FWHM (eV)", "σ", "Status",
               "Assignment (editable)", "Candidates and notes")

    def __init__(self, energy, intensity, energy_label: str, elements,
                 reference: CL.EnergyReference, label: str, parent=None, *,
                 settings=None, info=None, source_path="", temperature=30.0):
        super().__init__(parent)
        self.setWindowTitle(f"Core-level check — {label}")
        self.energy = np.asarray(energy, dtype=float)
        self.intensity = np.asarray(intensity, dtype=float)
        self.energy_label = energy_label
        self.elements = list(elements)
        self.reference = reference
        self.label = label
        self.info = dict(info or {})
        self.source_path = source_path
        self.temperature = temperature
        self.report = None
        settings = dict(settings or {})
        self._items = []
        #: peaks added by hand and found peaks removed, as positions on the
        #: data's own axis so they follow a change of the energy reference
        self.manual = []
        self.ignored = []
        self.zooms = []
        self.peak_table = None
        self.preview_energy = None

        central = QWidget()
        self.setCentralWidget(central)
        layout = QVBoxLayout(central)

        # -- controls -------------------------------------------------------
        controls = QGroupBox("Comparison")
        row = QHBoxLayout(controls)
        self.elements_label = QLabel("")
        row.addWidget(self.elements_label)
        self.elements_button = QPushButton("Elements...")
        self.elements_button.clicked.connect(self.change_elements)
        row.addWidget(self.elements_button)
        self.reference_button = QPushButton("Energy reference...")
        self.reference_button.clicked.connect(self.change_reference)
        row.addWidget(self.reference_button)
        row.addSpacing(12)
        row.addWidget(QLabel("Match within"))
        self.tolerance = spin(float(settings.get("tolerance", CL.TOLERANCE_EV)),
                              0.05, 5.0, 0.05, 2, " eV")
        row.addWidget(self.tolerance)
        row.addWidget(QLabel("Uncertain within"))
        self.uncertain = spin(float(settings.get("uncertain", CL.UNCERTAIN_EV)),
                              0.1, 10.0, 0.1, 2, " eV")
        row.addWidget(self.uncertain)
        row.addWidget(QLabel("Peak threshold"))
        self.sensitivity = spin(float(settings.get("sensitivity", CL.SENSITIVITY)),
                                1.0, 100.0, 0.5, 1, " σ")
        self.sensitivity.setToolTip(
            "How far above the noise a maximum must stand to count as a "
            "peak. Lower it for weak lines, raise it if noise is marked.")
        row.addWidget(self.sensitivity)
        row.addWidget(QLabel("Min height"))
        self.min_relative = spin(
            100 * float(settings.get("min_relative", CL.MIN_RELATIVE)),
            0.0, 20.0, 0.1, 2, " %")
        self.min_relative.setToolTip(
            "Smallest peak kept, in % of the strongest one. The noise "
            "threshold already rejects noise; raise this only to hide weak "
            "structure next to a strong line.")
        row.addWidget(self.min_relative)
        self.chemical = QCheckBox("Chemical states")
        self.chemical.setChecked(bool(settings.get("chemical_states", True)))
        self.chemical.setToolTip("Compare typical solid-state values (oxides, "
                                 "chlorides, carbon species...) too.")
        row.addWidget(self.chemical)
        self.analyse_button = QPushButton("Re-analyse")
        self.analyse_button.clicked.connect(self.analyse)
        row.addWidget(self.analyse_button)
        row.addStretch(1)
        layout.addWidget(controls)

        # -- the peaks by hand ------------------------------------------------
        manual = QGroupBox("Peaks by hand")
        rows = QVBoxLayout(manual)
        row = QHBoxLayout()
        rows.addLayout(row)
        self.log_scale = QCheckBox("Log counts")
        self.log_scale.setChecked(bool(settings.get("log_scale", False)))
        self.log_scale.setToolTip(
            "Logarithmic count axis: weak lines next to a strong one show.")
        self.log_scale.toggled.connect(lambda *_: self._draw())
        row.addWidget(self.log_scale)
        row.addSpacing(12)
        self.click_add = QPushButton("Click to add peak")
        self.click_add.setCheckable(True)
        self.click_add.setToolTip(
            "While pressed, a click on the plot puts a peak there (a "
            "double-click always does). Its position is refined by a local "
            "fit, and it is named with the levels within "
            f"±{CL.HINT_EV:g} eV.")
        self.click_add.toggled.connect(
            lambda on: on and self.range_button.setChecked(False))
        row.addWidget(self.click_add)
        row.addWidget(QLabel("or at"))
        self.add_energy = spin(0.0, -100.0, 20000.0, 0.1, 2, " eV")
        self.add_energy.setToolTip("Binding energy of a peak to add.")
        row.addWidget(self.add_energy)
        self.add_button = QPushButton("Add")
        self.add_button.clicked.connect(
            lambda: self.add_peak(self.add_energy.value()))
        row.addWidget(self.add_button)
        row.addSpacing(12)
        self.remove_button = QPushButton("Remove selected")
        self.remove_button.setToolTip(
            "Drop the selected peaks: one added by hand goes, a found one "
            "is ignored until Reset.")
        self.remove_button.clicked.connect(self.remove_selected)
        row.addWidget(self.remove_button)
        self.reset_button = QPushButton("Reset peaks")
        self.reset_button.setToolTip("Forget the peaks added by hand and "
                                     "bring back the ones removed.")
        self.reset_button.clicked.connect(self.reset_peaks)
        row.addWidget(self.reset_button)
        row.addStretch(1)

        row = QHBoxLayout()
        rows.addLayout(row)
        self.table_button = QPushButton("Peak table...")
        self.table_button.setToolTip(
            "Every tabulated level, filtered by element and binding-energy "
            "range: look a peak up, show a level on the plot, add it.")
        self.table_button.clicked.connect(lambda: self.open_peak_table())
        row.addWidget(self.table_button)
        row.addSpacing(12)
        self.range_button = QPushButton("Select range")
        self.range_button.setCheckable(True)
        self.range_button.setToolTip(
            "Show a band on the plot; drag it and its edges over the range "
            "to look at, then open it in a zoom window.")
        self.range_button.toggled.connect(self._range_toggled)
        row.addWidget(self.range_button)
        self.zoom_button = QPushButton("Zoom window...")
        self.zoom_button.setToolTip(
            "Open the selected range (or, with no range shown, what the "
            "plot shows now) in a window of its own.")
        self.zoom_button.clicked.connect(lambda: self.open_zoom())
        row.addWidget(self.zoom_button)
        row.addStretch(1)
        layout.addWidget(manual)
        self.hint = QLabel("")
        self.hint.setWordWrap(True)
        self.hint.setTextInteractionFlags(Qt.TextSelectableByMouse)
        self.hint.setStyleSheet("color: #555555;")
        layout.addWidget(self.hint)

        # -- plot, table, notes -------------------------------------------
        splitter = QSplitter(Qt.Vertical)
        self.plot = pg.PlotWidget()
        self.plot.setLabel("bottom", "Binding energy (eV)")
        self.plot.setLabel("left", "counts (angle-integrated)")
        self.plot.invertX(True)          # the XPS convention: deeper to the left
        strip_stock_menu(self.plot.getPlotItem())
        self.curve = self.plot.plot(pen=pg.mkPen("#33415c", width=1.3))
        self.markers = pg.ScatterPlotItem(size=10)
        self.plot.addItem(self.markers)
        self.highlight = pg.InfiniteLine(angle=90, pen=pg.mkPen("#888888",
                                                                style=Qt.DotLine))
        self.highlight.setVisible(False)
        self.plot.addItem(self.highlight, ignoreBounds=True)
        self.preview = pg.InfiniteLine(angle=90, pen=pg.mkPen(
            PREVIEW_COLOUR, width=1.5, style=Qt.DashLine))
        self.preview.setVisible(False)
        self.plot.addItem(self.preview, ignoreBounds=True)
        self.region = pg.LinearRegionItem(brush=pg.mkBrush(31, 111, 139, 40))
        self.region.setZValue(-5)
        self.region.setVisible(False)
        self.plot.addItem(self.region, ignoreBounds=True)
        self.plot.scene().sigMouseClicked.connect(self._mouse_clicked)
        self.show_lines = QCheckBox("Show the sample's tabulated levels")
        self.show_lines.setChecked(True)
        self.show_lines.toggled.connect(lambda *_: self._draw())
        splitter.addWidget(self.plot)

        self.table = QTableWidget(0, len(self.COLUMNS))
        self.table.setHorizontalHeaderLabels(self.COLUMNS)
        header = self.table.horizontalHeader()
        for column in range(len(self.COLUMNS) - 1):
            header.setSectionResizeMode(column, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(len(self.COLUMNS) - 1, QHeaderView.Stretch)
        # one line per peak; the whole story of the selected one is in the
        # pane below, and in each cell's tooltip
        self.table.setWordWrap(False)
        self.table.verticalHeader().setVisible(False)
        self.table.setSelectionBehavior(QTableWidget.SelectRows)
        self.table.itemSelectionChanged.connect(self._row_selected)
        splitter.addWidget(self.table)

        self.notes = QPlainTextEdit()
        self.notes.setReadOnly(True)
        splitter.addWidget(self.notes)
        splitter.setSizes([440, 230, 150])
        layout.addWidget(splitter, stretch=1)

        buttons = QHBoxLayout()
        buttons.addWidget(self.show_lines)
        buttons.addStretch(1)
        self.shift_button = QPushButton("Apply suggested shift")
        self.shift_button.setToolTip(
            "Move E_F by the rigid shift the check found, as for a charged "
            "sample, and analyse again.")
        self.shift_button.clicked.connect(self.apply_shift)
        buttons.addWidget(self.shift_button)
        self.list_button = QPushButton("Add EDC to list")
        self.list_button.setToolTip("List the integrated EDC, on the binding-"
                                    "energy scale, as a dataset of its own.")
        self.list_button.clicked.connect(self.edc_to_list)
        buttons.addWidget(self.list_button)
        self.copy_button = QPushButton("Copy report")
        self.copy_button.clicked.connect(self.copy_report)
        buttons.addWidget(self.copy_button)
        layout.addLayout(buttons)

        self.statusBar()
        self.resize(1100, 860)
        self.analyse()

    # -- analysis ----------------------------------------------------------
    def binding(self):
        return self.reference.binding(self.energy)

    def analyse(self):
        own = self._own_assignments()
        self.report = CL.check(self.binding(), self.intensity, self.elements,
                               self.reference,
                               sensitivity=self.sensitivity.value(),
                               tolerance=self.tolerance.value(),
                               uncertain=self.uncertain.value(),
                               chemical_states=self.chemical.isChecked(),
                               min_relative=self.min_relative.value() / 100,
                               manual=self.reference.binding(self.manual),
                               ignored=self.reference.binding(self.ignored))
        self.elements_label.setText("Sample: " + ", ".join(self.elements))
        self._draw()
        self._fill_table(own)
        self.notes.setPlainText(self._notes_text())
        self.shift_button.setEnabled(self.report.shift is not None)
        counts = {s: sum(a.status == s for a in self.report.assignments)
                  for s in CL.STATUSES}
        by_hand = sum(p.manual for p in self.report.peaks)
        self.statusBar().showMessage(
            f"{len(self.report.peaks)} peak(s): {counts['sample']} sample, "
            f"{counts['ambiguous']} ambiguous, {counts['uncertain']} uncertain, "
            f"{counts['unexplained']} not the sample"
            + (f"; {by_hand} added by hand" if by_hand else "")
            + (f"; {len(self.ignored)} removed" if self.ignored else ""))
        self.analysed.emit()
        return self.report

    def _own_assignments(self) -> dict:
        """Binding energy -> what the user typed, for every row whose
        assignment is no longer the automatic one, so a re-analysis keeps
        it."""
        own = {}
        if self.report is None:
            return own
        for row, a in enumerate(self.report.assignments):
            item = self.table.item(row, 5)
            if item is not None and item.text() != self._assignment_text(a):
                own[a.peak.binding] = item.text()
        return own

    def _notes_text(self) -> str:
        report, rows = self.report, []
        ref = report.reference
        rows.append(f"E_F at {ref.fermi_level:.3f} eV on the axis "
                    f"({CL.SCALE_LABELS[ref.scale]})"
                    + (f", hν = {ref.photon_energy:g} eV" if ref.photon_energy
                       else ", hν unknown")
                    + f", W = {ref.work_function:g} eV.")
        if report.missing:
            rows.append("Sample levels in the window with no peak: "
                        + "; ".join(f"{line.label} {line.energy:g} eV{why}"
                                    for line, why in report.missing))
        for symbol, text in report.out_of_window.items():
            rows.append(f"{symbol}: {text}")
        rows += report.notes
        rows.append("Peaks marked ◆ were added by hand. Peak table... lists "
                    "every tabulated level; Zoom window... looks at a range "
                    "closely.")
        rows.append("Status: sample = within the tolerance of a level of the "
                    "sample; uncertain = within a chemical shift, check it; "
                    "not the sample = the candidates listed are common "
                    "elements (contamination?), higher-order light or Auger "
                    "lines. A core level stays put in binding energy when hν "
                    "changes; an Auger line and higher-order light do not.")
        return "\n".join(rows)

    def _y(self, values):
        return plotted_counts(values, self.intensity, self.log_scale.isChecked())

    def _draw(self):
        for item in self._items:
            self.plot.removeItem(item)
        self._items = []
        log = self.log_scale.isChecked()
        self.plot.getAxis("left").setLogMode(log)
        self.plot.setLabel("left", "counts (log scale)" if log
                           else "counts (angle-integrated)")
        b = self.binding()
        order = np.argsort(b)
        self.curve.setData(b[order], self._y(self.intensity[order]))
        report = self.report
        if report is None:
            return
        spots = []
        for n, a in enumerate(report.assignments, 1):
            colour = STATUS_COLOURS[a.status]
            y = float(self._y(a.peak.height))
            spots.append({"pos": (a.peak.binding, y),
                          "brush": pg.mkBrush(colour), "pen": pg.mkPen("w"),
                          "symbol": MANUAL_SYMBOL if a.peak.manual
                          else FOUND_SYMBOL,
                          "size": 12 if a.peak.manual else 10,
                          "data": n - 1})
            text = pg.TextItem(f"{n}: {a.label}", color=colour, anchor=(0.5, 1.2))
            text.setPos(a.peak.binding, y)
            self.plot.addItem(text, ignoreBounds=True)
            self._items.append(text)
        self.markers.setData(spots)
        finite = self._y(self.intensity[np.isfinite(self.intensity)])
        if finite.size:
            # room above the tallest peak for its label
            lo, hi = float(finite.min()), float(finite.max())
            self.plot.setYRange(lo, hi + 0.15 * (hi - lo), padding=0.02)
        if self.show_lines.isChecked():
            lo, hi = report.window
            for line in CL.reference_lines(self.elements, chemical_states=False,
                                           window=(lo, hi)):
                item = pg.InfiniteLine(
                    line.energy, angle=90,
                    pen=pg.mkPen("#9aa5b8", width=1, style=Qt.DashLine))
                # named on hover: printed labels pile up where a doublet's
                # two lines and a peak's own label all sit
                item.setToolTip(f"{line.name}: {line.energy:g} eV (Table 1-1)")
                item.setZValue(-10)
                self.plot.addItem(item, ignoreBounds=True)
                self._items.append(item)

    def _fill_table(self, own=None):
        report = self.report
        own = dict(own or {})
        self.table.blockSignals(True)
        self.table.setRowCount(len(report.assignments))
        for row, a in enumerate(report.assignments):
            p = a.peak
            details = self._candidates_text(a)
            assigned = self._assignment_text(a)
            if own:
                # what the user typed for this peak before the re-analysis
                previous = min(own, key=lambda e: abs(e - p.binding))
                if abs(previous - p.binding) <= 0.15:
                    assigned = own.pop(previous)
            values = (f"{row + 1} ◆" if p.manual else str(row + 1),
                      f"{p.binding:.2f}"
                      + (f" ± {p.error:.2f}" if np.isfinite(p.error) else ""),
                      f"{p.fwhm:.2f}" if np.isfinite(p.fwhm) else "—",
                      f"{p.significance:.0f}" if np.isfinite(p.significance)
                      else "—",
                      CL.STATUS_LABELS[a.status], assigned,
                      details.replace("\n", "  |  "))
            for column, value in enumerate(values):
                item = QTableWidgetItem(value)
                if column == 0 and p.manual:
                    item.setToolTip("Added by hand")
                if column == 6 and details:
                    item.setToolTip(details)
                if column != 5:
                    item.setFlags(item.flags() & ~Qt.ItemIsEditable)
                if column == 4:
                    item.setForeground(QColor(STATUS_COLOURS[a.status]))
                self.table.setItem(row, column, item)
        self.table.blockSignals(False)

    @staticmethod
    def _assignment_text(a) -> str:
        if a.status == "sample":
            m = a.matches[0]
            return f"{m.line.label} (Δ {m.delta:+.2f})"
        if a.status == "ambiguous":
            return " or ".join(f"{m.line.label} (Δ {m.delta:+.2f})"
                               for m in a.matches[:3])
        return ""

    @staticmethod
    def _candidates_text(a) -> str:
        parts = []
        if a.status in ("sample", "ambiguous"):
            doublet = [m for m in a.matches if m.partner is not None]
            if doublet:
                parts.append(f"doublet partner: peak {doublet[0].partner + 1}")
            extra = a.matches[1:] if a.status == "sample" else []
            if extra:
                parts.append("also within tolerance: "
                             + ", ".join(m.describe() for m in extra[:2]))
        if a.nearby:
            parts.append("nearest sample levels: "
                         + "; ".join(m.describe() for m in a.nearby[:3]))
        if a.contaminants:
            parts.append("possible contamination: "
                         + "; ".join(m.describe() for m in a.contaminants))
        if a.higher_order:
            parts.append("higher-order light: "
                         + "; ".join(m.describe() for m in a.higher_order))
        parts += a.notes
        sample, others = a.suggestions
        if sample:
            parts.append(f"sample within ±{CL.HINT_EV:g} eV: "
                         + "; ".join(m.describe() for m in sample))
        if others:
            parts.append(f"other elements within ±{CL.HINT_EV:g} eV: "
                         + "; ".join(m.describe() for m in others))
        return "\n".join(parts)

    def _row_selected(self):
        rows = {index.row() for index in self.table.selectedIndexes()}
        if len(rows) != 1 or self.report is None:
            self.highlight.setVisible(False)
            self.notes.setPlainText(self._notes_text())
            return
        row = rows.pop()
        a = self.report.assignments[row]
        self.highlight.setValue(a.peak.binding)
        self.highlight.setVisible(True)
        p = a.peak
        head = (f"Peak {row + 1} at {p.binding:.2f} eV (FWHM {p.fwhm:.2f} eV, "
                f"{p.significance:.0f}σ): {CL.STATUS_LABELS[a.status]}"
                + (" — added by hand" if p.manual else ""))
        assigned = self._assignment_text(a)
        body = ([f"Assignment: {assigned}"] if assigned else []) + \
            [self._candidates_text(a) or "Nothing near it in the tables."]
        if not p.manual:
            body.append(f"Within ±{CL.HINT_EV:g} eV: " + CL.describe_candidates(
                p.binding, self.elements,
                chemical_states=self.chemical.isChecked()))
        self.notes.setPlainText("\n".join([head] + body + ["", "—" * 20,
                                                          self._notes_text()]))

    # -- the peaks by hand --------------------------------------------------
    def _view_x(self, scene_pos):
        """The binding energy under a scene position, or None off the plot."""
        box = self.plot.getPlotItem().vb
        if not box.sceneBoundingRect().contains(scene_pos):
            return None
        return float(box.mapSceneToView(scene_pos).x())

    def _mouse_clicked(self, event):
        if event.button() != Qt.LeftButton:
            return
        if not (self.click_add.isChecked() or event.double()):
            return
        energy = self._view_x(event.scenePos())
        if energy is not None:
            event.accept()
            self.add_peak(energy)

    def add_peak(self, binding: float, label: str = None):
        """Put a peak at ``binding`` eV and analyse again; the new peak's
        row is selected, so its hints show. ``label`` (a level picked in
        the peak table) is written in as its assignment."""
        binding = float(binding)
        b = self.binding()
        if not np.nanmin(b) <= binding <= np.nanmax(b):
            self.statusBar().showMessage(
                f"{binding:.2f} eV is outside the window "
                f"({np.nanmin(b):.2f} to {np.nanmax(b):.2f} eV)")
            return None
        position = float(self.reference.axis(binding))
        if any(abs(position - e) < 0.05 for e in self.manual):
            return self.report      # the second click of a double-click
        self.manual.append(position)
        report = self.analyse()
        rows = [n for n, a in enumerate(report.assignments) if a.peak.manual]
        if rows:
            row = min(rows, key=lambda n: abs(
                report.assignments[n].peak.binding - binding))
            self.table.selectRow(row)
            a = report.assignments[row]
            if label:
                self.table.item(row, 5).setText(label)
            self.hint.setText(f"Added peak {row + 1} at {a.peak.binding:.2f} eV "
                              "— " + CL.describe_candidates(
                                  a.peak.binding, self.elements,
                                  chemical_states=self.chemical.isChecked()))
        return report

    def remove_selected(self):
        """Drop the selected peaks: one added by hand is forgotten, a found
        one is ignored until :meth:`reset_peaks`."""
        if self.report is None:
            return None
        rows = sorted({index.row() for index in self.table.selectedIndexes()})
        if not rows:
            self.statusBar().showMessage("Select the peaks to remove first")
            return None
        for row in rows:
            p = self.report.assignments[row].peak
            if p.manual:
                if self.manual:
                    positions = self.reference.binding(self.manual)
                    del self.manual[int(np.argmin(np.abs(positions - p.binding)))]
            else:
                self.ignored.append(float(self.reference.axis(p.binding)))
        return self.analyse()

    def reset_peaks(self):
        self.manual, self.ignored = [], []
        return self.analyse()

    # -- the peak table and the zoom windows ---------------------------------
    def view_range(self):
        """The binding-energy range the plot shows now, low to high."""
        lo, hi = self.plot.getPlotItem().vb.viewRange()[0]
        return min(lo, hi), max(lo, hi)

    def selected_range(self):
        """The band's range when it is shown, else what the plot shows."""
        if self.region.isVisible():
            lo, hi = self.region.getRegion()
            return min(lo, hi), max(lo, hi)
        return self.view_range()

    def _range_toggled(self, on):
        if on:
            self.click_add.setChecked(False)
            lo, hi = self.view_range()
            middle, half = 0.5 * (lo + hi), 0.05 * (hi - lo)
            self.region.setRegion((middle - half, middle + half))
        self.region.setVisible(on)

    def open_peak_table(self, window=None):
        """The peak table, filtered to ``window`` (binding, lo-hi) when
        given; one table per results window, raised if already open."""
        if self.peak_table is None:
            self.peak_table = PeakTableDialog(
                self.elements, self, window=window or self.window_range(),
                chemical_states=self.chemical.isChecked())
            self.peak_table.addRequested.connect(self.add_peak)
            self.peak_table.previewRequested.connect(self.show_preview)
        elif window is not None:
            self.peak_table.set_range(*window)
        self.peak_table.show()
        self.peak_table.raise_()
        return self.peak_table

    def window_range(self):
        b = self.binding()
        return float(np.nanmin(b)), float(np.nanmax(b))

    def show_preview(self, energy, label: str = ""):
        """Mark a tabulated level (from the peak table) on the plots."""
        self.preview_energy = None if energy is None else float(energy)
        if energy is None:
            self.preview.setVisible(False)
        else:
            self.preview.setValue(float(energy))
            self.preview.setToolTip(f"{label}: {energy:g} eV")
            self.preview.setVisible(True)
        self.previewChanged.emit(self.preview_energy, label)

    def open_zoom(self, lo: float = None, hi: float = None):
        """A zoom window on ``lo``-``hi`` eV binding, by default the
        selected range."""
        if lo is None or hi is None:
            lo, hi = self.selected_range()
        zoom = ZoomWindow(self, lo, hi)
        self.zooms.append(zoom)
        zoom.closed.connect(lambda z: z in self.zooms and self.zooms.remove(z))
        zoom.show()
        return zoom

    # -- the buttons --------------------------------------------------------
    def change_elements(self):
        dialog = ElementsDialog(self, elements=" ".join(self.elements),
                                settings=self._settings())
        if dialog.exec_() != QDialog.Accepted:
            return None
        self.elements = dialog.elements()[0]
        self.tolerance.setValue(dialog.tolerance.value())
        self.uncertain.setValue(dialog.uncertain.value())
        self.chemical.setChecked(dialog.chemical.isChecked())
        return self.analyse()

    def change_reference(self):
        dialog = EnergyReferenceDialog(self.energy, self.intensity,
                                       self.energy_label, self.info, self,
                                       temperature=self.temperature,
                                       reference=self.reference)
        if dialog.exec_() != QDialog.Accepted:
            return None
        self.reference = dialog.reference()
        return self.analyse()

    def apply_shift(self):
        if self.report is None or self.report.shift is None:
            return None
        shift = self.report.shift[0]
        ref = self.reference
        # E_B = E_F - E (or E - E_F): moving E_F by -shift (+shift) moves
        # every binding energy by -shift
        level = ref.fermi_level - shift if ref.scale != "binding" \
            else ref.fermi_level + shift
        self.reference = CL.EnergyReference(ref.scale, level, ref.photon_energy,
                                            ref.work_function)
        report = self.analyse()
        self.statusBar().showMessage(
            f"Every binding energy moved by {-shift:+.2f} eV; "
            + self.statusBar().currentMessage())
        return report

    def edc_to_list(self):
        from arpes_viewer.ui.viewers.curves import curve_dataset

        b = self.binding()
        order = np.argsort(b)
        ref = self.reference
        name = f"{self.label} [core-level EDC]"
        data = curve_dataset(
            "edc", b[order], self.intensity[order][:, None], ["EDC"],
            x_label="Binding energy (eV)",
            value_label="Intensity (angle-integrated)", name=name,
            step="core_level_edc",
            parameters={"elements": ", ".join(self.elements),
                        "scale": ref.scale, "fermi_level": ref.fermi_level,
                        "photon_energy": ref.photon_energy or 0.0,
                        "work_function": ref.work_function},
            prefix="corelevel", source_info=self.info,
            source_path=self.source_path, source=self.label)
        self.datasetsCreated.emit([data])
        self.statusBar().showMessage(f"Added “{name}” to the list")
        return data

    def copy_report(self):
        text = self.report_text()
        QGuiApplication.clipboard().setText(text)
        self.statusBar().showMessage("Report copied to the clipboard")
        return text

    def report_text(self) -> str:
        text = self.report.text()
        own = []
        for row, a in enumerate(self.report.assignments):
            item = self.table.item(row, 5)
            if item is not None and item.text() != self._assignment_text(a):
                own.append(f"Peak {row + 1}: {item.text()}")
        if self.ignored:
            text += ("\n\nRemoved by hand (found, but not peaks): "
                     + ", ".join(f"{e:.2f} eV" for e in
                                 self.reference.binding(self.ignored)))
        if own:
            text += "\n\nAssigned by hand:\n" + "\n".join(own)
        return text

    def _settings(self) -> dict:
        return {"tolerance": self.tolerance.value(),
                "uncertain": self.uncertain.value(),
                "chemical_states": self.chemical.isChecked(),
                "sensitivity": self.sensitivity.value(),
                "min_relative": self.min_relative.value() / 100,
                "log_scale": self.log_scale.isChecked()}

    def closeEvent(self, event):
        for zoom in list(self.zooms):
            zoom.close()
        if self.peak_table is not None:
            self.peak_table.close()
        save_last_used({**load_last_used(), **self._settings()})
        self.closed.emit(self)
        super().closeEvent(event)


# --------------------------------------------------------------------------
# The peak table
# --------------------------------------------------------------------------
class _NumberItem(QTableWidgetItem):
    """A cell that sorts by its number, not its text."""

    def __init__(self, value: float, text: str):
        super().__init__(text)
        self.value = float(value)

    def __lt__(self, other):
        if isinstance(other, _NumberItem):
            return self.value < other.value
        return super().__lt__(other)


class PeakTableDialog(QDialog):
    """Every tabulated level (Table 1-1 and the typical chemical states),
    filtered by element and binding-energy range, for the user to look a
    peak up by eye. Picking a row marks the level on the plots; *Add as
    peak* (or a double-click) adds a peak there, refined by a local fit and
    with the level written in as its assignment."""

    #: (binding energy, label)
    addRequested = pyqtSignal(float, str)
    #: (binding energy or None, label)
    previewRequested = pyqtSignal(object, str)

    COLUMNS = ("Element", "Level", "E_B (eV)", "Source", "In sample")

    def __init__(self, sample, parent=None, *, window=(0.0, 1500.0),
                 chemical_states: bool = True):
        super().__init__(parent)
        self.setWindowTitle("Core-level table")
        self.setModal(False)
        self.sample = list(sample)
        self.lines = []
        layout = QVBoxLayout(self)

        form = QFormLayout()
        row = QHBoxLayout()
        self.elements_box = QLineEdit("")
        self.elements_box.setPlaceholderText(
            "all elements  —  or e.g. Nb S Cl, NbSCl")
        self.elements_box.textChanged.connect(lambda *_: self.refresh())
        row.addWidget(self.elements_box, stretch=1)
        sample_button = QPushButton("Sample")
        sample_button.setToolTip("Only the sample's elements.")
        sample_button.clicked.connect(
            lambda: self.elements_box.setText(" ".join(self.sample)))
        row.addWidget(sample_button)
        all_button = QPushButton("All")
        all_button.clicked.connect(lambda: self.elements_box.setText(""))
        row.addWidget(all_button)
        form.addRow("Elements", row)

        row = QHBoxLayout()
        self.lo_box = spin(0.0, -100.0, 20000.0, 0.5, 2, " eV")
        self.hi_box = spin(0.0, -100.0, 20000.0, 0.5, 2, " eV")
        self.lo_box.valueChanged.connect(lambda *_: self.refresh())
        self.hi_box.valueChanged.connect(lambda *_: self.refresh())
        row.addWidget(self.lo_box)
        row.addWidget(QLabel("to"))
        row.addWidget(self.hi_box)
        self._window = tuple(sorted(window))
        whole = QPushButton("Whole window")
        whole.setToolTip("The binding-energy range of the data.")
        whole.clicked.connect(lambda: self.set_range(*self._window))
        row.addWidget(whole)
        if parent is not None and hasattr(parent, "view_range"):
            shown = QPushButton("Plot view")
            shown.setToolTip("What the results plot shows now (or its "
                             "selected range).")
            shown.clicked.connect(lambda: self.set_range(*parent.selected_range()))
            row.addWidget(shown)
        row.addStretch(1)
        form.addRow("Binding energy", row)
        self.chemical = QCheckBox("Typical chemical states too")
        self.chemical.setChecked(chemical_states)
        self.chemical.toggled.connect(lambda *_: self.refresh())
        form.addRow("", self.chemical)
        layout.addLayout(form)

        self.count = QLabel("")
        layout.addWidget(self.count)
        self.table = QTableWidget(0, len(self.COLUMNS))
        self.table.setHorizontalHeaderLabels(self.COLUMNS)
        self.table.verticalHeader().setVisible(False)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table.setSelectionBehavior(QTableWidget.SelectRows)
        self.table.setSelectionMode(QTableWidget.SingleSelection)
        header = self.table.horizontalHeader()
        for column in range(len(self.COLUMNS) - 2):
            header.setSectionResizeMode(column, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(3, QHeaderView.Stretch)
        header.setSectionResizeMode(4, QHeaderView.ResizeToContents)
        self.table.itemSelectionChanged.connect(self._selected)
        self.table.cellDoubleClicked.connect(lambda *_: self.add_selected())
        layout.addWidget(self.table, stretch=1)

        buttons = QHBoxLayout()
        self.add_button = QPushButton("Add as peak")
        self.add_button.setToolTip(
            "Add a peak at the selected level (refined by a local fit) with "
            "the level as its assignment. A double-click does the same.")
        self.add_button.clicked.connect(self.add_selected)
        buttons.addWidget(self.add_button)
        buttons.addStretch(1)
        close = QPushButton("Close")
        close.clicked.connect(self.close)
        buttons.addWidget(close)
        layout.addLayout(buttons)
        self.resize(560, 640)
        self.set_range(*self._window)

    def set_range(self, lo: float, hi: float):
        lo, hi = sorted((float(lo), float(hi)))
        for box, value in ((self.lo_box, lo), (self.hi_box, hi)):
            box.blockSignals(True)
            box.setValue(value)
            box.blockSignals(False)
        self.refresh()

    def refresh(self):
        text = self.elements_box.text().strip()
        if text:
            elements, unknown = CL.parse_elements(text)
        else:
            elements = [s for z, s in sorted(CL.ELEMENTS.items())]
            unknown = []
        lo, hi = sorted((self.lo_box.value(), self.hi_box.value()))
        self.lines = CL.reference_lines(elements,
                                        chemical_states=self.chemical.isChecked(),
                                        window=(lo, hi))
        self.table.setSortingEnabled(False)
        self.table.clearSelection()
        self.table.setRowCount(len(self.lines))
        bold = self.font()
        bold.setBold(True)
        for n, line in enumerate(self.lines):
            mine = line.element in self.sample
            cells = (QTableWidgetItem(line.element),
                     QTableWidgetItem(line.orbital),
                     _NumberItem(line.energy, f"{line.energy:g}"),
                     QTableWidgetItem("Table 1-1" if line.source == "table"
                                      else line.source),
                     QTableWidgetItem("yes" if mine else ""))
            cells[0].setData(Qt.UserRole, n)
            for column, cell in enumerate(cells):
                if mine:
                    cell.setForeground(QColor(STATUS_COLOURS["sample"]))
                    cell.setFont(bold)
                self.table.setItem(n, column, cell)
        self.table.setSortingEnabled(True)
        self.table.sortItems(2, Qt.AscendingOrder)
        message = f"{len(self.lines)} level(s) between {lo:g} and {hi:g} eV"
        if unknown:
            message += "  —  not an element: " + ", ".join(unknown)
        self.count.setText(message)
        self.add_button.setEnabled(False)
        self.previewRequested.emit(None, "")

    def selected_line(self):
        rows = {index.row() for index in self.table.selectedIndexes()}
        if len(rows) != 1:
            return None
        item = self.table.item(rows.pop(), 0)
        return self.lines[item.data(Qt.UserRole)]

    def select_line(self, line) -> bool:
        """Select the row of ``line``; True if it is in the table."""
        for row in range(self.table.rowCount()):
            if self.lines[self.table.item(row, 0).data(Qt.UserRole)] == line:
                self.table.selectRow(row)
                return True
        return False

    def _selected(self):
        line = self.selected_line()
        self.add_button.setEnabled(line is not None)
        if line is None:
            self.previewRequested.emit(None, "")
        else:
            self.previewRequested.emit(line.energy, line.label)

    def add_selected(self):
        line = self.selected_line()
        if line is None:
            return None
        self.addRequested.emit(line.energy, line.label)
        return line

    def closeEvent(self, event):
        self.previewRequested.emit(None, "")
        super().closeEvent(event)


# --------------------------------------------------------------------------
# A zoom window
# --------------------------------------------------------------------------
class ZoomWindow(QMainWindow):
    """A range of the results window's EDC on its own, to place a peak by
    eye: a movable cursor, the smoothed curve and its negative second
    derivative (a shoulder with no maximum of its own is a maximum there),
    a Gaussian fit at the cursor, and the peaks and sample levels of the
    results window, which it follows."""

    closed = pyqtSignal(object)

    def __init__(self, owner: CoreLevelWindow, lo: float, hi: float):
        super().__init__(owner)
        self.owner = owner
        self.setWindowTitle(f"Zoom — {owner.label}")
        self._items = []
        central = QWidget()
        self.setCentralWidget(central)
        layout = QVBoxLayout(central)

        row = QHBoxLayout()
        row.addWidget(QLabel("From"))
        self.lo_box = spin(0.0, -100.0, 20000.0, 0.1, 2, " eV")
        row.addWidget(self.lo_box)
        row.addWidget(QLabel("to"))
        self.hi_box = spin(0.0, -100.0, 20000.0, 0.1, 2, " eV")
        row.addWidget(self.hi_box)
        apply_button = QPushButton("Show")
        apply_button.clicked.connect(lambda: self.set_range(
            self.lo_box.value(), self.hi_box.value()))
        # Enter in From / to shows the range, as Show does
        self.enter_redraws = lambda: self.set_range(self.lo_box.value(),
                                                    self.hi_box.value())
        row.addWidget(apply_button)
        row.addSpacing(12)
        self.log_scale = QCheckBox("Log counts")
        self.log_scale.setChecked(owner.log_scale.isChecked())
        self.log_scale.toggled.connect(lambda *_: self.redraw())
        row.addWidget(self.log_scale)
        row.addWidget(QLabel("Smoothing"))
        self.smoothing = spin(0.1, 0.0, 5.0, 0.02, 2, " eV")
        self.smoothing.setToolTip("Gaussian smoothing of the dark curve and "
                                  "the second derivative; 0 = raw data.")
        self.smoothing.valueChanged.connect(lambda *_: self.redraw())
        row.addWidget(self.smoothing)
        self.derivative = QCheckBox("−d²I/dE²")
        self.derivative.setChecked(True)
        self.derivative.setToolTip(
            "The negative second derivative of the smoothed curve: every "
            "peak, shoulders included, is a maximum of it.")
        self.derivative.toggled.connect(lambda *_: self.redraw())
        row.addWidget(self.derivative)
        self.relative = QCheckBox("÷ I")
        self.relative.setChecked(True)
        self.relative.setToolTip(
            "Divide the second derivative by the smoothed intensity, so a "
            "weak peak shows as clearly as a strong one (for a Gaussian it "
            "is 1/σ² at the centre whatever the height).")
        self.relative.toggled.connect(lambda *_: self.redraw())
        row.addWidget(self.relative)
        table_button = QPushButton("Peak table...")
        table_button.setToolTip("The tabulated levels in this range.")
        table_button.clicked.connect(lambda: owner.open_peak_table(self.range()))
        row.addWidget(table_button)
        row.addStretch(1)
        layout.addLayout(row)

        self.plot = pg.PlotWidget()
        self.plot.invertX(True)
        self.plot.setLabel("left", "counts")
        strip_stock_menu(self.plot.getPlotItem())
        self.raw = self.plot.plot(pen=pg.mkPen("#9aa5b8", width=1))
        self.smooth = self.plot.plot(pen=pg.mkPen("#33415c", width=1.5))
        self.markers = pg.ScatterPlotItem(size=10)
        self.plot.addItem(self.markers)
        self.fit_curve = self.plot.plot(pen=pg.mkPen("#bd4921", width=1.5))
        self.preview = pg.InfiniteLine(angle=90, pen=pg.mkPen(
            PREVIEW_COLOUR, width=1.5, style=Qt.DashLine))
        self.preview.setVisible(False)
        self.plot.addItem(self.preview, ignoreBounds=True)
        self.cursor = pg.InfiniteLine(
            angle=90, movable=True, pen=pg.mkPen("#1f6f8b", width=2),
            label="{value:.2f} eV", labelOpts={"position": 0.95,
                                              "color": "#1f6f8b"})
        self.cursor.setToolTip("Drag it onto a peak; a click on the plot "
                               "puts it there, a double-click adds a peak.")
        self.plot.addItem(self.cursor, ignoreBounds=True)
        self.cursor.sigPositionChanged.connect(self._cursor_moved)
        self.plot.scene().sigMouseClicked.connect(self._mouse_clicked)

        self.d2_plot = pg.PlotWidget()
        self.d2_plot.invertX(True)
        self.d2_plot.setLabel("left", "−d²I/dE²")
        self.d2_plot.setLabel("bottom", "Binding energy (eV)")
        self.d2_plot.setXLink(self.plot)
        self.d2_plot.setMaximumHeight(170)
        strip_stock_menu(self.d2_plot.getPlotItem())
        self.d2_curve = self.d2_plot.plot(pen=pg.mkPen("#2f7a3f", width=1.3))
        self.d2_plot.addItem(pg.InfiniteLine(0, angle=0, pen=pg.mkPen(
            "#bbbbbb")), ignoreBounds=True)
        self.d2_cursor = pg.InfiniteLine(angle=90, pen=pg.mkPen(
            "#1f6f8b", width=1, style=Qt.DotLine))
        self.d2_plot.addItem(self.d2_cursor, ignoreBounds=True)
        splitter = QSplitter(Qt.Vertical)
        splitter.addWidget(self.plot)
        splitter.addWidget(self.d2_plot)
        splitter.setSizes([480, 160])
        layout.addWidget(splitter, stretch=1)

        row = QHBoxLayout()
        self.readout = QLabel("")
        self.readout.setTextInteractionFlags(Qt.TextSelectableByMouse)
        row.addWidget(self.readout, stretch=1)
        fit_button = QPushButton("Fit at cursor")
        fit_button.setToolTip("A Gaussian on a sloping background over 1 eV "
                              "either side of the cursor.")
        fit_button.clicked.connect(self.fit_at_cursor)
        row.addWidget(fit_button)
        add_button = QPushButton("Add peak at cursor")
        add_button.clicked.connect(self.add_at_cursor)
        row.addWidget(add_button)
        layout.addLayout(row)

        owner.analysed.connect(self.redraw)
        owner.previewChanged.connect(self._preview)
        self.resize(900, 700)
        self.set_range(lo, hi)
        self._preview(owner.preview_energy, "")

    # -- data ----------------------------------------------------------------
    def range(self):
        return self._range

    def _data(self):
        """Binding energy and counts inside the range, ascending."""
        b = self.owner.binding()
        y = self.owner.intensity
        lo, hi = self._range
        inside = np.isfinite(b) & np.isfinite(y) & (b >= lo) & (b <= hi)
        order = np.argsort(b[inside])
        return b[inside][order], y[inside][order]

    def smoothed(self):
        """``(grid, smoothed, -d2)`` of the range on a uniform grid."""
        from scipy.ndimage import gaussian_filter1d

        b, y = self._data()
        if b.size < 5:
            return b, y, np.zeros_like(y)
        step = float(np.median(np.diff(b)))
        grid = np.arange(b[0], b[-1] + 0.5 * step, step)
        raw = np.interp(grid, b, y)
        sigma = max(self.smoothing.value() / step, 1.0)
        smooth = gaussian_filter1d(raw, sigma, mode="nearest") \
            if self.smoothing.value() > 0 else raw
        d2 = -gaussian_filter1d(raw, sigma, order=2, mode="nearest") / step ** 2
        if self.relative.isChecked():
            level = gaussian_filter1d(raw, sigma, mode="nearest")
            positive = level[level > 0]
            floor = 0.01 * float(np.median(positive)) if positive.size else 1.0
            d2 = d2 / np.maximum(level, floor)
        return grid, smooth, d2

    # -- drawing -------------------------------------------------------------
    def set_range(self, lo: float, hi: float):
        lo, hi = sorted((float(lo), float(hi)))
        if hi - lo < 1e-3:
            hi = lo + 1.0
        self._range = (lo, hi)
        for box, value in ((self.lo_box, lo), (self.hi_box, hi)):
            box.blockSignals(True)
            box.setValue(value)
            box.blockSignals(False)
        self.setWindowTitle(f"Zoom {lo:.1f}–{hi:.1f} eV — {self.owner.label}")
        self.fit_curve.setData([], [])
        self.redraw()
        self.cursor.setValue(0.5 * (lo + hi))
        self._cursor_moved()

    def redraw(self):
        for item in self._items:
            self.plot.removeItem(item)
        self._items = []
        log = self.log_scale.isChecked()
        intensity = self.owner.intensity

        def y(values):
            return plotted_counts(values, intensity, log)
        self.plot.getAxis("left").setLogMode(log)
        self.plot.setLabel("left", "counts (log scale)" if log else "counts")
        b, counts = self._data()
        self.raw.setData(b, y(counts))
        grid, smooth, d2 = self.smoothed()
        self.smooth.setData(grid, y(smooth))
        self.d2_plot.setVisible(self.derivative.isChecked())
        self.relative.setEnabled(self.derivative.isChecked())
        self.d2_plot.setLabel("left", "−d²I/dE² ÷ I" if self.relative.isChecked()
                              else "−d²I/dE²")
        self.d2_curve.setData(grid, d2)
        lo, hi = self._range
        spots = []
        report = self.owner.report
        for n, a in enumerate(report.assignments if report else [], 1):
            if not lo <= a.peak.binding <= hi:
                continue
            colour = STATUS_COLOURS[a.status]
            h = float(y(a.peak.height))
            spots.append({"pos": (a.peak.binding, h), "brush": pg.mkBrush(colour),
                          "pen": pg.mkPen("w"), "size": 12 if a.peak.manual else 10,
                          "symbol": MANUAL_SYMBOL if a.peak.manual
                          else FOUND_SYMBOL})
            text = pg.TextItem(f"{n}: {a.label}", color=colour, anchor=(0.5, 1.2))
            text.setPos(a.peak.binding, h)
            self.plot.addItem(text, ignoreBounds=True)
            self._items.append(text)
        self.markers.setData(spots)
        for line in CL.reference_lines(self.owner.elements, chemical_states=False,
                                       window=(lo, hi)):
            item = pg.InfiniteLine(line.energy, angle=90, pen=pg.mkPen(
                "#9aa5b8", width=1, style=Qt.DashLine),
                label=line.name, labelOpts={"position": 0.85, "rotateAxis": (1, 0),
                                            "color": "#77809a"})
            item.setZValue(-10)
            self.plot.addItem(item, ignoreBounds=True)
            self._items.append(item)
        self.plot.setXRange(lo, hi, padding=0.01)
        finite = y(counts)
        if finite.size:
            bottom, top = float(np.nanmin(finite)), float(np.nanmax(finite))
            self.plot.setYRange(bottom, top + 0.15 * (top - bottom), padding=0.02)
        if d2.size:
            self.d2_plot.setYRange(float(np.nanmin(d2)), float(np.nanmax(d2)),
                                   padding=0.05)

    def _preview(self, energy, label):
        if energy is None:
            self.preview.setVisible(False)
        else:
            self.preview.setValue(float(energy))
            self.preview.setToolTip(f"{label}: {energy:g} eV")
            self.preview.setVisible(True)

    # -- the cursor ------------------------------------------------------------
    def cursor_energy(self) -> float:
        return float(self.cursor.value())

    def _cursor_moved(self, *_):
        energy = self.cursor_energy()
        self.d2_cursor.setValue(energy)
        b, counts = self._data()
        text = f"Cursor {energy:.2f} eV"
        if b.size:
            text += f", {np.interp(energy, b, counts):.0f} counts"
        self.readout.setText(text)

    def _mouse_clicked(self, event):
        if event.button() != Qt.LeftButton:
            return
        box = self.plot.getPlotItem().vb
        if not box.sceneBoundingRect().contains(event.scenePos()):
            return
        energy = float(box.mapSceneToView(event.scenePos()).x())
        event.accept()
        self.cursor.setValue(energy)
        if event.double():
            self.add_at_cursor()

    def fit_at_cursor(self):
        """Fit a Gaussian at the cursor, draw it and move the cursor to its
        centre; the peak, or None when the fit did not settle."""
        b, counts = self._data()
        peak = CL.peak_at(b, counts, self.cursor_energy())
        if not np.isfinite(peak.fwhm):
            self.fit_curve.setData([], [])
            self.readout.setText(f"No peak fits at {self.cursor_energy():.2f} eV.")
            return None
        self.cursor.setValue(peak.binding)
        sigma = peak.fwhm / 2.3548
        near = (b >= peak.binding - 3 * sigma) & (b <= peak.binding + 3 * sigma)
        x = b[near]
        base = peak.height - peak.prominence
        model = base + peak.prominence * np.exp(-0.5 * ((x - peak.binding)
                                                        / sigma) ** 2)
        self.fit_curve.setData(x, plotted_counts(model, self.owner.intensity,
                                                 self.log_scale.isChecked()))
        self.readout.setText(
            f"Fit: {peak.binding:.2f} ± {peak.error:.2f} eV, FWHM "
            f"{peak.fwhm:.2f} eV, {peak.significance:.0f}σ")
        return peak

    def add_at_cursor(self):
        return self.owner.add_peak(self.cursor_energy())

    def closeEvent(self, event):
        try:
            self.owner.analysed.disconnect(self.redraw)
            self.owner.previewChanged.disconnect(self._preview)
        except TypeError:
            pass
        self.closed.emit(self)
        super().closeEvent(event)


# --------------------------------------------------------------------------
# From a dataset
# --------------------------------------------------------------------------
def cut_edc(data):
    """``(energy, intensity, energy label)`` of a cut, summed over angle."""
    if data.kind != "cut":
        raise ValueError("the core-level check works on a cut")
    scan = data.scan
    energy, intensity = CL.integrated_edc(scan.value, scan.y, energy_dim=1)
    return energy, intensity, scan.labels.get("y", "Energy (eV)")


def start(data, label: str, parent=None, *, temperature: float = 30.0,
          elements_dialog=None, reference_dialog=None):
    """Run the three steps on a cut; the results window, or None if the
    user stopped on the way. The two ``*_dialog`` arguments let a test
    stand in for the modal dialogs."""
    energy, intensity, energy_label = cut_edc(data)
    info = dict(getattr(data.scan, "info", {}) or {})

    dialog = elements_dialog or ElementsDialog(parent)
    if elements_dialog is None and dialog.exec_() != QDialog.Accepted:
        return None
    elements = dialog.elements()[0]
    # the results window's own settings (threshold, log axis...) as last used
    settings = {**load_last_used(), **dialog.settings()}

    ref_dialog = reference_dialog or EnergyReferenceDialog(
        energy, intensity, energy_label, info, parent, temperature=temperature)
    if reference_dialog is None and ref_dialog.exec_() != QDialog.Accepted:
        return None
    reference = ref_dialog.reference()

    return CoreLevelWindow(energy, intensity, energy_label, elements,
                           reference, label, parent, settings=settings,
                           info=info, source_path=getattr(data, "path", ""),
                           temperature=temperature)
