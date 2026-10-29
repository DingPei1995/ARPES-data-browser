"""
ui/panels/repeats.py
====================
The "Sum repeated scans" window: a measurement taken in batches, checked
and added back together, step by step.

Reached from the tool panel's Datasets box with two or more cuts (or maps
of one kind) selected. Three steps:

1. **Are these repeats?** -- the grids and the recorded conditions of
   every scan against the first.
2. **Has the sample changed?** -- each scan against the sum of the others:
   chi2 per point, count rate and energy shift, with the changed ones
   flagged and their EDC drawn against the rest.
3. **Sum** -- which scans go in, whether the energy shifts are undone, sum
   or mean, and the name.

The arithmetic is :mod:`tools.processing.repeats`, which imports no Qt.
"""
from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QColor
from PyQt5.QtWidgets import (QAbstractItemView, QCheckBox, QComboBox, QDialog,
                             QFormLayout, QHBoxLayout, QHeaderView, QLabel,
                             QLineEdit, QPushButton, QTableWidget,
                             QTableWidgetItem, QVBoxLayout, QWidget)

from arpes_viewer.loader.nxs_file import axis_slots
from arpes_viewer.tools.processing import cutops
from arpes_viewer.tools.processing import process as P
from arpes_viewer.tools.processing import repeats as R
from arpes_viewer.ui.common.screen import fits_screen
from arpes_viewer.ui.common.stepper import Step, Stepper, guard_window
from arpes_viewer.ui.viewers.widgets import MemoryData

__all__ = ["RepeatsDialog"]

_CHANGED = QColor(255, 214, 196)
_DIFFERS = QColor(255, 226, 180)


@fits_screen
class RepeatsDialog(QDialog):
    """Sum repeats of one measurement. ``datasets`` is ``[(name, data),
    ...]`` in the list's order, all of one kind."""

    datasetsCreated = pyqtSignal(list)

    def __init__(self, datasets, parent=None, existing_names=None):
        super().__init__(parent)
        self.setWindowTitle("Sum repeated scans")
        self.resize(980, 640)
        self.datasets = list(datasets)
        self.names = [name for name, _ in self.datasets]
        self.kind = self.datasets[0][1].kind
        if any(data.kind != self.kind for _, data in self.datasets):
            raise ValueError("Repeated scans are one measurement taken again: "
                             "select datasets of one kind only.")
        self._existing = existing_names or (lambda: [])
        self.slots = axis_slots(self.kind, "constructor")
        self.axes = [tuple(np.asarray(getattr(data.scan, s), dtype=float)
                           for s in self.slots) for _, data in self.datasets]
        self.labels = dict(self.datasets[0][1].scan.labels or {})
        self.energy_dim = len(self.slots) - 1
        self.diagnosis = None
        self.result = None
        try:
            self.plan = R.grid_plan(self.axes)
            self.plan_problem = ""
        except ValueError as exc:
            self.plan, self.plan_problem = None, str(exc)

        layout = QVBoxLayout(self)
        self.stepper = Stepper([self._step_check(), self._step_compare(),
                                self._step_sum()], self,
                               finish_text="Sum to list")
        self.stepper.finished.connect(self.sum_to_list)
        layout.addWidget(self.stepper)
        guard_window(self, self.stepper, "sum_repeats", "Sum repeated scans")

    # -- data --------------------------------------------------------------
    def values(self, i):
        data = self.datasets[i][1]
        return np.asarray(data.cut_frame if self.kind == "cut"
                          else data.scan.value, dtype=float)

    # -- step 1 ------------------------------------------------------------
    def _step_check(self):
        page = QWidget()
        column = QVBoxLayout(page)
        self.check_table = QTableWidget()
        self.check_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        column.addWidget(self.check_table, 1)
        self.grid_label = QLabel("")
        self.grid_label.setWordWrap(True)
        column.addWidget(self.grid_label)
        self._fill_check_table()
        text = ("A fine cut or map taken as many short scans is one "
                "measurement, and only adds up if every scan really is the "
                "same one. The table lists each scan's grid and what it "
                "recorded; <b>orange</b> marks a value that differs from the "
                "first scan (photon energy, polarisation, angles, position, "
                "temperature, pass energy, lens mode). A small difference "
                "in the grid is fine -- the scans are interpolated onto the "
                "first one's points -- but a different position or photon "
                "energy means a different measurement: leave that scan out "
                "in step 3.")
        return Step("Are these repeats?", text, page,
                    ready=lambda: self.plan_problem)

    def _fill_check_table(self):
        infos = [dict(data.scan.info or {}) for _, data in self.datasets]
        quantities = []
        for row in cutops.metadata_differences(infos[0], infos[0]):
            quantities.append(row[0])
        for info in infos[1:]:
            for row in cutops.metadata_differences(infos[0], info):
                if row[0] not in quantities:
                    quantities.append(row[0])
        headers = ["Scan", "Points"] + quantities
        table = self.check_table
        table.setColumnCount(len(headers))
        table.setHorizontalHeaderLabels(headers)
        table.setRowCount(len(self.datasets))
        for i, (name, _data) in enumerate(self.datasets):
            table.setItem(i, 0, QTableWidgetItem(name))
            points = " × ".join(str(a.size) for a in self.axes[i])
            table.setItem(i, 1, QTableWidgetItem(points))
            rows = {r[0]: r for r in cutops.metadata_differences(infos[0], infos[i])}
            for j, quantity in enumerate(quantities):
                row = rows.get(quantity)
                value = "" if row is None or row[2] is None else str(row[2])
                item = QTableWidgetItem(value)
                if row is not None and row[3]:
                    item.setBackground(_DIFFERS)
                    item.setToolTip(f"First scan: {row[1]}")
                table.setItem(i, j + 2, item)
        table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        if self.plan_problem:
            self.grid_label.setText(f"<b>Cannot be summed:</b> {self.plan_problem}")
        elif self.plan["exact"]:
            self.grid_label.setText("All scans are on the same grid, point for "
                                    "point: they are added as measured.")
        else:
            self.grid_label.setText("<b>Grids differ:</b> "
                                    + "; ".join(self.plan["notes"]) + ".")

    # -- step 2 ------------------------------------------------------------
    def _step_compare(self):
        page = QWidget()
        column = QVBoxLayout(page)
        top = QHBoxLayout()
        self.compare_button = QPushButton("Compare the scans")
        self.compare_button.setStyleSheet("font-weight: bold;")
        self.compare_button.clicked.connect(lambda: self.compare())
        top.addWidget(self.compare_button)
        self.compare_status = QLabel("")
        self.compare_status.setWordWrap(True)
        top.addWidget(self.compare_status, 1)
        column.addLayout(top)
        row = QHBoxLayout()
        self.diag_table = QTableWidget(0, 5)
        self.diag_table.setHorizontalHeaderLabels(
            ["Scan", "Counts (rel.)", "χ² / point", "E shift (meV)",
             "Verdict"])
        self.diag_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.diag_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.diag_table.setSelectionMode(QAbstractItemView.SingleSelection)
        self.diag_table.itemSelectionChanged.connect(self._show_scan)
        self.diag_table.horizontalHeader().setSectionResizeMode(
            QHeaderView.ResizeToContents)
        row.addWidget(self.diag_table, 3)
        plots = QVBoxLayout()
        self.trend_plot = pg.PlotWidget(background="w")
        self.trend_plot.setLabel("bottom", "scan")
        self.trend_plot.setLabel("left", "χ² / point  (o),  counts rel. (□)")
        self.trend_plot.setMinimumHeight(160)
        plots.addWidget(self.trend_plot)
        self.edc_plot = pg.PlotWidget(background="w")
        self.edc_plot.setLabel("bottom", self.labels.get(self.slots[-1], "Energy"))
        self.edc_plot.setLabel("left", "EDC (normalised)")
        self.edc_plot.addLegend(offset=(5, 5))
        self.edc_plot.setMinimumHeight(160)
        plots.addWidget(self.edc_plot)
        row.addLayout(plots, 4)
        column.addLayout(row, 1)
        text = ("Press <b>Compare the scans</b>. Each scan is compared with "
                "the sum of all the others, allowing for counting noise: "
                "<b>χ² / point</b> is about the same for scans that "
                "differ only by noise and grows when the spectrum itself "
                "changed -- a band that moved or broadened as the sample "
                "aged. Scans well above the rest are marked "
                "<b>changed</b>. A falling <b>count rate</b> with an unchanged "
                "spectrum is the beam, not the sample. The <b>energy "
                "shift</b> is how far each scan's EDC sits from the others' "
                "(a drifting photon energy or surface potential). Click a "
                "row to see its EDC against the rest.")
        return Step("Has the sample changed?", text, page,
                    ready=lambda: "" if self.diagnosis is not None
                    else "Press “Compare the scans” first.")

    def compare(self, background: bool = True):
        if self.plan is None:
            return None

        def work(report=None):
            return R.diagnose(self.values, len(self.datasets), self.axes,
                              self.energy_dim, self.plan, progress=report)

        if background:
            from arpes_viewer.ui.common import jobs
            try:
                self.diagnosis = jobs.run_blocking(self, "Comparing the scans",
                                                   lambda report: work(report))
            except jobs.JobCancelled:
                return None
        else:
            self.diagnosis = work()
        self._fill_diagnosis()
        self.stepper.refresh()
        return self.diagnosis

    def _fill_diagnosis(self):
        d = self.diagnosis
        table = self.diag_table
        table.setRowCount(len(self.datasets))
        for i, name in enumerate(self.names):
            cells = [name, f"{d.relative_totals[i]:.3f}", f"{d.chi2[i]:.2f}",
                     f"{d.shifts[i] * 1000:+.1f}",
                     "changed" if d.flags[i] else "consistent"]
            for j, text in enumerate(cells):
                item = QTableWidgetItem(text)
                if d.flags[i]:
                    item.setBackground(_CHANGED)
                table.setItem(i, j, item)
        flagged = sum(d.flags)
        self.compare_status.setText(
            (f"{flagged} of {len(self.names)} scans look changed; they are "
             f"left out in step 3 unless you tick them back in."
             if flagged else "No scan stands out from the others.")
            + ("" if not d.notes else " Note: " + "; ".join(d.notes) + "."))
        self.trend_plot.clear()
        x = np.arange(1, len(self.names) + 1)
        self.trend_plot.plot(x, d.chi2, pen=pg.mkPen("#1f6f8b", width=2),
                             symbol="o", symbolBrush="#1f6f8b")
        self.trend_plot.plot(x, d.relative_totals, pen=pg.mkPen("#999", width=1),
                             symbol="s", symbolBrush="#bbb")
        for i in np.flatnonzero(d.flags):
            self.trend_plot.plot([x[i]], [d.chi2[i]], pen=None, symbol="o",
                                 symbolSize=14, symbolBrush=None,
                                 symbolPen=pg.mkPen("#d9534f", width=2))
        for i, check in enumerate(getattr(self, "include_boxes", [])):
            check.setChecked(not d.flags[i])
        if self.names:
            table.selectRow(int(np.argmax(np.where(np.isfinite(d.chi2), d.chi2, -1))))

    def _show_scan(self):
        rows = self.diag_table.selectionModel().selectedRows()
        if not rows or self.diagnosis is None:
            return
        i = rows[0].row()
        target = self.plan["target"]
        profiles = [R.edc_profile(R.onto_grid(self.values(j), self.axes[j], target),
                                  self.energy_dim) for j in range(len(self.names))]
        mine = profiles[i]
        rest = np.sum(profiles, axis=0) - mine
        energy = target[self.energy_dim]
        self.edc_plot.clear()
        scale = lambda c: c / (np.nanmax(c) or 1.0)          # noqa: E731
        self.edc_plot.plot(energy, scale(rest), pen=pg.mkPen("#888", width=2),
                           name="the others")
        self.edc_plot.plot(energy, scale(mine), pen=pg.mkPen("#d9534f", width=2),
                           name=self.names[i])

    # -- step 3 ------------------------------------------------------------
    def _step_sum(self):
        page = QWidget()
        column = QVBoxLayout(page)
        self.include_boxes = []
        includes = QVBoxLayout()
        for name in self.names:
            box = QCheckBox(name)
            box.setChecked(True)
            box.toggled.connect(lambda *_: self._summary())
            includes.addWidget(box)
            self.include_boxes.append(box)
        column.addLayout(includes)
        form = QFormLayout()
        self.align_box = QCheckBox("Undo the energy shifts before summing")
        self.align_box.setToolTip(
            "Move every scan back by its measured shift; the energy axis is "
            "then cropped to the range all of them still cover.")
        self.align_box.toggled.connect(lambda *_: self._summary())
        form.addRow(self.align_box)
        self.mode_combo = QComboBox()
        self.mode_combo.addItems(["Sum (counts stay counts)", "Mean"])
        form.addRow("Combine as", self.mode_combo)
        self.name_edit = QLineEdit(self._default_name())
        form.addRow("Name", self.name_edit)
        column.addLayout(form)
        self.sum_summary = QLabel("")
        self.sum_summary.setWordWrap(True)
        column.addWidget(self.sum_summary)
        column.addStretch(1)
        text = ("Tick the scans to add up -- the ones step 2 marked as "
                "changed start unticked. If step 2 found energy shifts "
                "larger than about a third of a point, <b>undo the energy "
                "shifts</b>; otherwise leave the scans as measured. "
                "<b>Sum</b> keeps counts as counts (what every fit and error "
                "bar expects); <b>Mean</b> keeps the scale of one scan. "
                "<b>Sum to list</b> adds the result as a new dataset.")
        return Step("Sum", text, page, ready=self._sum_ready,
                    enter=self._enter_sum)

    def _default_name(self):
        base = f"{self.names[0]} ×{len(self.names)} summed"
        taken = set(self._existing())
        name, n = base, 2
        while name in taken:
            name = f"{base} {n}"
            n += 1
        return name

    def _enter_sum(self):
        if self.diagnosis is not None:
            step = R._step(self.plan["target"][self.energy_dim])
            kept = [abs(s) for s, flagged in zip(self.diagnosis.shifts,
                                                 self.diagnosis.flags)
                    if not flagged] or [0.0]
            big = np.nanmax(kept) > step / 3
            self.align_box.setChecked(bool(big))
        self._summary()

    def included(self):
        return [i for i, box in enumerate(self.include_boxes) if box.isChecked()]

    def _sum_ready(self):
        return "" if self.included() else "Tick at least one scan."

    def _summary(self):
        chosen = self.included()
        left = [self.names[i] for i in range(len(self.names)) if i not in chosen]
        self.sum_summary.setText(
            f"{len(chosen)} of {len(self.names)} scans go in."
            + (f" Left out: {', '.join(left)}." if left else ""))
        self.stepper.refresh()

    # -- the result --------------------------------------------------------
    def build(self):
        """The summed dataset (not yet listed)."""
        include = self.included()
        shifts = None
        if self.align_box.isChecked() and self.diagnosis is not None:
            shifts = self.diagnosis.shifts
        mode = "mean" if self.mode_combo.currentIndex() == 1 else "sum"
        axes, values, notes = R.sum_repeats(
            self.values, self.axes, include, self.energy_dim, shifts=shifts,
            mode=mode)
        first = self.datasets[include[0]][1]
        name = self.name_edit.text().strip() or self._default_name()
        members = [self.names[i] for i in include]
        left = [self.names[i] for i in range(len(self.names)) if i not in include]
        parameters = {"members": members, "left_out": left, "mode": mode,
                      "aligned": shifts is not None,
                      "notes": " ".join(notes)}
        if self.diagnosis is not None:
            parameters["chi2"] = np.asarray(self.diagnosis.chi2, dtype=float)
            parameters["shift_eV"] = np.asarray(self.diagnosis.shifts, dtype=float)
        info = P.record_step(dict(first.scan.info or {}), P.Step(
            "sum_repeats", {"n": len(members), "mode": mode,
                            "aligned": shifts is not None},
            source=", ".join(members)))
        return MemoryData(self.kind, axes, values, dict(first.scan.labels or {}),
                          source_label=name, parameters=parameters,
                          prefix="repeats",
                          source_path=getattr(first, "path", ""),
                          source_info=info,
                          source_motors=dict(getattr(first.scan, "fourd_info", {})
                                             or {}))

    def sum_to_list(self):
        try:
            self.result = self.build()
        except ValueError as exc:
            self.stepper.missing.setText(str(exc))
            return None
        self.datasetsCreated.emit([self.result])
        self.stepper.missing.setText(f"Added {self.result.source_label} to the list.")
        return self.result
