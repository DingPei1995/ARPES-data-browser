"""
ui/panels/ridge.py
==================
The ridge-detection window: bands traced as ridges of the intensity.

Opened from the tool panel (Advanced -> Ridge detection...) with one or more
cuts, or maps, selected. The source on the left with the ridge points on
it, the ridge image on the right; the smoothing widths are in each axis'
own unit. A map is shown one plane at a time and processed plane by plane
on the planes you choose (constant-energy maps, or cuts along either
angle).

What comes out: the ridge image as a new dataset (each selected dataset,
same kind and grid), and the sub-pixel ridge points as a CSV file -- band
positions ready for a dispersion fit or the tight-binding fit.

The algorithm is :mod:`tools.processing.ridge`, which imports no Qt.
"""
from __future__ import annotations

import csv

import numpy as np
import pyqtgraph as pg
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtWidgets import (QComboBox, QDialog, QDoubleSpinBox, QFileDialog,
                             QFormLayout, QHBoxLayout, QLabel, QPushButton,
                             QSlider, QVBoxLayout, QWidget, QCheckBox,
                             QGroupBox)

from loader.nxs_file import axis_slots
from tools.processing import process as P
from tools.processing import ridge as RD
from ui.common.screen import fits_screen
from ui.panels.process import PreviewPair
from ui.viewers.widgets import MemoryData

__all__ = ["RidgeDialog", "PLANES"]

#: For a cube: (what the planes are, the axis held constant).
PLANES = (("Constant-energy planes", 2),
          ("Cuts along the analyser slit (one per scan step)", 0),
          ("Cuts along the scan axis (one per slit angle)", 1))


def _unit(label: str) -> str:
    label = str(label or "")
    return label[label.rfind("(") + 1:label.rfind(")")] if "(" in label else ""


@fits_screen
class RidgeDialog(QDialog):
    """Ridge detection on ``datasets`` (``[(name, data), ...]``, cuts or
    cubes; processed in one go if they are all of one kind)."""

    datasetsCreated = pyqtSignal(list)

    def __init__(self, datasets, parent=None, colormap="Greys",
                 existing_names=None):
        super().__init__(parent)
        self.setWindowTitle("Ridge detection")
        self.resize(1100, 700)
        self.datasets = list(datasets)
        self.colormap = colormap
        self._existing = existing_names or (lambda: [])
        self.result = None
        self.points = np.zeros((0, 3))

        layout = QVBoxLayout(self)
        hint = QLabel(
            "<b>How to use:</b> 1. pick the dataset (and, for a map, the "
            "planes and one to look at). 2. Set the <b>smoothing</b> along "
            "each axis, in its own unit, to about the band's width -- too "
            "little and noise makes short ridges, too much and close bands "
            "merge. 3. Raise the <b>threshold</b> until only the bands are "
            "left. The yellow points on the left are the ridge points, each "
            "placed to a fraction of a pixel. 4. <b>Result to list</b> for "
            "the ridge image, <b>Save ridge points</b> for the points (a "
            "CSV the tight-binding fit and any plotting program read).")
        hint.setWordWrap(True)
        hint.setStyleSheet("background: #eef4f8; border: 1px solid #c9dbe6; "
                           "border-radius: 4px; padding: 6px;")
        layout.addWidget(hint)

        self.preview = PreviewPair(self)
        layout.addWidget(self.preview, 1)
        self.scatter = pg.ScatterPlotItem(size=3, pen=None,
                                          brush=pg.mkBrush(255, 210, 0, 220))
        self.preview.before["view"].view.addItem(self.scatter)

        controls = QHBoxLayout()
        what = QGroupBox("Data")
        form = QFormLayout(what)
        self.dataset_combo = QComboBox()
        self.dataset_combo.addItems([name for name, _ in self.datasets])
        self.dataset_combo.currentIndexChanged.connect(self._dataset_changed)
        form.addRow("Dataset", self.dataset_combo)
        self.plane_combo = QComboBox()
        self.plane_combo.addItems([text for text, _ in PLANES])
        self.plane_combo.currentIndexChanged.connect(self._plane_changed)
        form.addRow("Planes", self.plane_combo)
        self.slice_slider = QSlider(Qt.Horizontal)
        self.slice_slider.valueChanged.connect(self._slice_changed)
        self.slice_label = QLabel("")
        slice_row = QHBoxLayout()
        slice_row.addWidget(self.slice_slider, 1)
        slice_row.addWidget(self.slice_label)
        form.addRow("Shown", slice_row)
        controls.addWidget(what, 1)

        settings = QGroupBox("Settings")
        sform = QFormLayout(settings)
        self.sigma_x = QDoubleSpinBox()
        self.sigma_y = QDoubleSpinBox()
        for box in (self.sigma_x, self.sigma_y):
            box.setDecimals(5)
            box.setRange(1e-6, 1e6)
            box.setKeyboardTracking(False)
            box.valueChanged.connect(self.update_preview)
        self.sigma_x_label, self.sigma_y_label = QLabel(""), QLabel("")
        sform.addRow(self.sigma_x_label, self.sigma_x)
        sform.addRow(self.sigma_y_label, self.sigma_y)
        self.threshold = QDoubleSpinBox()
        self.threshold.setRange(0.0, 100.0)
        self.threshold.setSuffix(" % of the sharpest")
        self.threshold.setValue(5.0)
        self.threshold.setKeyboardTracking(False)
        self.threshold.valueChanged.connect(self.update_preview)
        sform.addRow("Threshold", self.threshold)
        self.gamma = QDoubleSpinBox()
        self.gamma.setRange(0.1, 5.0)
        self.gamma.setSingleStep(0.1)
        self.gamma.setValue(0.5)
        self.gamma.setToolTip("Contrast of the ridge image only: below 1 "
                              "brings out weak ridges. The points are not "
                              "affected.")
        self.gamma.valueChanged.connect(self.update_preview)
        sform.addRow("Image gamma", self.gamma)
        self.kind_combo = QComboBox()
        self.kind_combo.addItems(["Maxima (bands)", "Minima (gaps, valleys)"])
        self.kind_combo.currentIndexChanged.connect(self.update_preview)
        sform.addRow("Trace", self.kind_combo)
        self.show_points = QCheckBox("Show ridge points on the source")
        self.show_points.setChecked(True)
        self.show_points.toggled.connect(self.update_preview)
        sform.addRow(self.show_points)
        controls.addWidget(settings, 1)

        buttons = QVBoxLayout()
        self.to_list_button = QPushButton("Result to list")
        self.to_list_button.setToolTip(
            "The ridge image of every selected dataset as new datasets (a "
            "map plane by plane, on the planes chosen).")
        self.to_list_button.clicked.connect(self.results_to_list)
        buttons.addWidget(self.to_list_button)
        self.save_points_button = QPushButton("Save ridge points...")
        self.save_points_button.setToolTip(
            "Every ridge point of this dataset (a map: of every plane) to a "
            "CSV file: the coordinates and the ridge's sharpness.")
        self.save_points_button.clicked.connect(lambda: self.save_points())
        buttons.addWidget(self.save_points_button)
        self.status = QLabel("")
        self.status.setWordWrap(True)
        buttons.addWidget(self.status)
        buttons.addStretch(1)
        controls.addLayout(buttons, 1)
        layout.addLayout(controls)

        self._dataset_changed(0)

    # -- the data ----------------------------------------------------------
    @property
    def data(self):
        return self.datasets[max(0, self.dataset_combo.currentIndex())][1]

    def _axes_labels(self, data):
        slots = axis_slots(data.kind, "constructor")
        labels = dict(data.scan.labels or {})
        axes = [np.asarray(getattr(data.scan, s), dtype=float) for s in slots]
        return axes, [labels.get(s, s) for s in slots]

    def _values(self, data):
        return np.asarray(data.cut_frame if data.kind == "cut"
                          else data.scan.value, dtype=float)

    def plane_axes(self):
        """The shown plane's two axes and their labels, and the held one."""
        axes, labels = self._axes_labels(self.data)
        if self.data.kind == "cut":
            return (axes[0], axes[1]), (labels[0], labels[1]), None
        hold = PLANES[self.plane_combo.currentIndex()][1]
        keep = [d for d in range(3) if d != hold]
        return ((axes[keep[0]], axes[keep[1]]),
                (labels[keep[0]], labels[keep[1]]), hold)

    def plane(self):
        values = self._values(self.data)
        _axes, _labels, hold = self.plane_axes()
        if hold is None:
            return values
        index = [slice(None)] * 3
        index[hold] = self.slice_slider.value()
        return values[tuple(index)]

    def _dataset_changed(self, _index):
        cube = self.data.kind != "cut"
        self.plane_combo.setEnabled(cube)
        self.slice_slider.setEnabled(cube)
        self._plane_changed()

    def _plane_changed(self, *_):
        (x, y), (lx, ly), hold = self.plane_axes()
        if hold is not None:
            axes, _ = self._axes_labels(self.data)
            self.slice_slider.blockSignals(True)
            self.slice_slider.setRange(0, axes[hold].size - 1)
            self.slice_slider.setValue(axes[hold].size // 2)
            self.slice_slider.blockSignals(False)
        else:
            self.slice_label.setText("")
        self.sigma_x_label.setText(f"Smoothing along {lx.split('(')[0].strip()}")
        self.sigma_y_label.setText(f"Smoothing along {ly.split('(')[0].strip()}")
        sx, sy = RD.default_sigmas(x, y)
        for box, value, label in ((self.sigma_x, sx, lx), (self.sigma_y, sy, ly)):
            box.blockSignals(True)
            box.setSingleStep(value / 2)
            box.setValue(value)
            box.setSuffix(f" {_unit(label)}" if _unit(label) else "")
            box.blockSignals(False)
        self.preview.set_labels(lx, ly)
        self._slice_changed()

    def _slice_changed(self, *_):
        _axes, _labels, hold = self.plane_axes()
        if hold is not None:
            axes, labels = self._axes_labels(self.data)
            value = axes[hold][self.slice_slider.value()]
            self.slice_label.setText(f"{value:.4g} {_unit(labels[hold])}")
        self.update_preview()

    # -- the result --------------------------------------------------------
    def options(self):
        return {"threshold": self.threshold.value() / 100.0,
                "gamma": self.gamma.value(),
                "kind": "min" if self.kind_combo.currentIndex() == 1 else "max"}

    def update_preview(self, *_):
        (x, y), _labels, _hold = self.plane_axes()
        values = self.plane()
        result = RD.ridges(values, x, y, self.sigma_x.value(),
                           self.sigma_y.value(), **self.options())
        self.points = result.points
        self.preview.show_source(values, (x, y), self.colormap, False)
        self.preview.show_result(result.strength, (x, y), "Greys", True)
        if self.show_points.isChecked() and result.points.size:
            self.scatter.setData(result.points[:, 0], result.points[:, 1])
        else:
            self.scatter.setData([], [])
        self.preview.report(f"{len(result.points)} ridge points")
        return result

    def settings(self):
        """The window's settings, read once (on the GUI thread)."""
        return {"sigmas": (self.sigma_x.value(), self.sigma_y.value()),
                "options": self.options(),
                "hold": PLANES[self.plane_combo.currentIndex()][1]}

    def compute(self, data, settings=None):
        """``(ridge values, points)`` for a whole dataset. Touches no widget
        when ``settings`` are given, so it can run off the GUI thread."""
        settings = settings or self.settings()
        values = self._values(data)
        axes, _labels = self._axes_labels(data)
        sigmas = settings["sigmas"]
        if data.kind == "cut":
            result = RD.ridges(values, axes[0], axes[1], *sigmas, **settings["options"])
            return result.strength, result.points
        return RD.ridge_cube(values, axes, settings["hold"], sigmas,
                             **settings["options"])

    def _name(self, base):
        taken = set(self._existing())
        name, n = f"{base}_rd", 2
        while name in taken:
            name = f"{base}_rd{n}"
            n += 1
        return name

    def results_to_list(self):
        from ui.common.jobs import JobCancelled, busy, run_blocking

        made = []
        datasets = list(self.datasets)
        settings = self.settings()

        def work(report):
            out = []
            for index, (name, data) in enumerate(datasets):
                report(index / len(datasets),
                       f"Ridges of {name} ({index + 1} of {len(datasets)})")
                out.append(self.compute(data, settings))
            return out

        try:
            # a ridge search through every plane of a map: the slow part
            computed = (work(lambda *_a: None) if busy()
                        else run_blocking(self, "Ridge detection", work))
        except JobCancelled:
            return []
        for (name, data), (values, points) in zip(datasets, computed):
            axes, _labels = self._axes_labels(data)
            options = self.options()
            parameters = {"sigma_x": self.sigma_x.value(),
                          "sigma_y": self.sigma_y.value(),
                          "threshold": options["threshold"],
                          "kind": options["kind"], "points": int(len(points))}
            if data.kind != "cut":
                parameters["plane"] = PLANES[self.plane_combo.currentIndex()][0]
            info = P.record_step(dict(data.scan.info or {}),
                                 P.Step("ridge", parameters, source=name))
            made.append(MemoryData(
                data.kind, tuple(axes), values, dict(data.scan.labels or {}),
                source_label=self._name(name), parameters=parameters,
                prefix="ridge", source_path=getattr(data, "path", ""),
                source_info=info,
                source_motors=dict(getattr(data.scan, "fourd_info", {}) or {})))
        self.result = made
        self.datasetsCreated.emit(made)
        self.status.setText(f"Added {len(made)} ridge image(s) to the list.")
        return made

    def save_points(self, path=None):
        data = self.data
        _values, points = self.compute(data)
        axes, labels = self._axes_labels(data)
        if path is None:
            path, _ = QFileDialog.getSaveFileName(
                self, "Save ridge points", f"{self.dataset_combo.currentText()}"
                f"_ridges.csv", "CSV (*.csv)")
            if not path:
                return None
        if data.kind == "cut":
            header = [labels[0], labels[1], "strength"]
        else:
            hold = PLANES[self.plane_combo.currentIndex()][1]
            keep = [d for d in range(3) if d != hold]
            header = [labels[hold], labels[keep[0]], labels[keep[1]], "strength"]
        with open(path, "w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(header)
            for row in points:
                writer.writerow([f"{v:.8g}" for v in row])
        self.status.setText(f"Wrote {len(points)} ridge points.")
        return path
