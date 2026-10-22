"""
ui/corelevel.py
===============
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
   is missing. The assignment column is the user's to overwrite.

The analysis itself is :mod:`tools.corelevel`, with no Qt in it.
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
                             QPushButton, QCheckBox, QComboBox,
                             QDialogButtonBox, QTableWidget, QTableWidgetItem,
                             QHeaderView, QSplitter, QPlainTextEdit,
                             QMessageBox, QGroupBox)

import tools.corelevel as CL
from ui.process import spin
from ui.widgets import strip_stock_menu

#: Marker colour of each status, on the plot and in the table.
STATUS_COLOURS = {
    "sample": "#2f7a3f",
    "ambiguous": "#1f6f8b",
    "uncertain": "#c47f00",
    "unexplained": "#bd4921",
}

LAST_USED_NAME = "corelevel.json"


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
        self.report = CL.check(self.binding(), self.intensity, self.elements,
                               self.reference,
                               sensitivity=self.sensitivity.value(),
                               tolerance=self.tolerance.value(),
                               uncertain=self.uncertain.value(),
                               chemical_states=self.chemical.isChecked())
        self.elements_label.setText("Sample: " + ", ".join(self.elements))
        self._draw()
        self._fill_table()
        self.notes.setPlainText(self._notes_text())
        self.shift_button.setEnabled(self.report.shift is not None)
        counts = {s: sum(a.status == s for a in self.report.assignments)
                  for s in CL.STATUSES}
        self.statusBar().showMessage(
            f"{len(self.report.peaks)} peak(s): {counts['sample']} sample, "
            f"{counts['ambiguous']} ambiguous, {counts['uncertain']} uncertain, "
            f"{counts['unexplained']} not the sample")
        return self.report

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
        rows.append("Status: sample = within the tolerance of a level of the "
                    "sample; uncertain = within a chemical shift, check it; "
                    "not the sample = the candidates listed are common "
                    "elements (contamination?), higher-order light or Auger "
                    "lines. A core level stays put in binding energy when hν "
                    "changes; an Auger line and higher-order light do not.")
        return "\n".join(rows)

    def _draw(self):
        for item in self._items:
            self.plot.removeItem(item)
        self._items = []
        b = self.binding()
        order = np.argsort(b)
        self.curve.setData(b[order], self.intensity[order])
        report = self.report
        if report is None:
            return
        spots = []
        for n, a in enumerate(report.assignments, 1):
            colour = STATUS_COLOURS[a.status]
            spots.append({"pos": (a.peak.binding, a.peak.height),
                          "brush": pg.mkBrush(colour), "pen": pg.mkPen("w"),
                          "data": n - 1})
            text = pg.TextItem(f"{n}: {a.label}", color=colour, anchor=(0.5, 1.2))
            text.setPos(a.peak.binding, a.peak.height)
            self.plot.addItem(text, ignoreBounds=True)
            self._items.append(text)
        self.markers.setData(spots)
        finite = self.intensity[np.isfinite(self.intensity)]
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

    def _fill_table(self):
        report = self.report
        self.table.blockSignals(True)
        self.table.setRowCount(len(report.assignments))
        for row, a in enumerate(report.assignments):
            p = a.peak
            details = self._candidates_text(a)
            values = (str(row + 1), f"{p.binding:.2f}"
                      + (f" ± {p.error:.2f}" if np.isfinite(p.error) else ""),
                      f"{p.fwhm:.2f}", f"{p.significance:.0f}",
                      CL.STATUS_LABELS[a.status], self._assignment_text(a),
                      details.replace("\n", "  |  "))
            for column, value in enumerate(values):
                item = QTableWidgetItem(value)
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
                f"{p.significance:.0f}σ): {CL.STATUS_LABELS[a.status]}")
        assigned = self._assignment_text(a)
        body = ([f"Assignment: {assigned}"] if assigned else []) + \
            [self._candidates_text(a) or "Nothing near it in the tables."]
        self.notes.setPlainText("\n".join([head] + body + ["", "—" * 20,
                                                          self._notes_text()]))

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
        from ui.curves import curve_dataset

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
        for row in range(self.table.rowCount()):
            item = self.table.item(row, 5)
            automatic = self._assignment_text(self.report.assignments[row])
            if item is not None and item.text() != automatic:
                own.append(f"Peak {row + 1}: {item.text()}")
        if own:
            text += "\n\nAssigned by hand:\n" + "\n".join(own)
        return text

    def _settings(self) -> dict:
        return {"tolerance": self.tolerance.value(),
                "uncertain": self.uncertain.value(),
                "chemical_states": self.chemical.isChecked(),
                "sensitivity": self.sensitivity.value()}

    def closeEvent(self, event):
        save_last_used({**load_last_used(), **self._settings()})
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
    settings = dialog.settings()

    ref_dialog = reference_dialog or EnergyReferenceDialog(
        energy, intensity, energy_label, info, parent, temperature=temperature)
    if reference_dialog is None and ref_dialog.exec_() != QDialog.Accepted:
        return None
    reference = ref_dialog.reference()

    return CoreLevelWindow(energy, intensity, energy_label, elements,
                           reference, label, parent, settings=settings,
                           info=info, source_path=getattr(data, "path", ""),
                           temperature=temperature)
