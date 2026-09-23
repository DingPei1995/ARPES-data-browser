"""
ui/cutops.py
============
The cut arithmetic window: A and B, and what comes of combining them.

One window for every way two cuts are combined -- linear and circular
dichroism, dividing by a reference, plain differences, ratios and sums. It
replaces the old right-click "Compare the two..." dialog, which already computed a
difference, a ratio and an asymmetry under a name that promised only a
look; two implementations of the same arithmetic would sooner or later have
disagreed about a scale factor.

It is reached two ways and is the same window either way: the
**Cut arithmetic...** button in a cut's own viewer (that cut is A,
and the next cut clicked in the main list is B), or the main list's
right-click menu with two cuts selected.

The arithmetic is :mod:`tools.cutops`, which imports no Qt; this file shows
A and B on a common scale, the result on a scale centred on zero where the
result is signed, what the two measurements recorded, and what the
combination did.
"""
from __future__ import annotations

import numpy as np
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QColor
from PyQt5.QtWidgets import (QComboBox, QDialog, QDialogButtonBox,
                             QDoubleSpinBox, QFormLayout, QGroupBox,
                             QHBoxLayout, QLabel, QPushButton, QTableWidget,
                             QTableWidgetItem, QVBoxLayout, QWidget,
                             QAbstractItemView, QHeaderView)

from tools import cutops
from tools import process as P
from ui.process import PreviewPair
from ui.widgets import (add_button, MemoryData, apply_colormap, fit_frame_view,
                        plain_image_view, show_frame)

__all__ = ["CutArithmeticDialog"]

#: Short names for the result, by operation and preset.
_SUFFIX = {
    "linear_dichroism": "LD",
    "circular_dichroism": "CD",
    "reference": "norm",
    "difference": "diff",
    "asymmetry": "asym",
    "ratio": "ratio",
    "sum": "sum",
}


class CutArithmeticDialog(QDialog):
    """Combine two cuts. ``a`` and ``b`` are ``(name, data)`` pairs.

    ``region_source`` is an optional callable returning the selection box of
    A's viewer as ``(x0, y0, x1, y1)``, or None; when given, the
    normalisation region can be taken straight from it.
    """

    datasetsCreated = pyqtSignal(list)

    def __init__(self, a, b, parent=None, colormap="gray", flip=False,
                 region_source=None, existing_names=None):
        super().__init__(parent)
        self.setWindowTitle("Cut arithmetic")
        self.setModal(False)
        self.resize(980, 900)
        self.colormap, self.flip = colormap, flip
        self.region_source = region_source
        self.existing_names = existing_names or (lambda: [])
        self.last = None
        self._set_pair(a, b)

        layout = QVBoxLayout(self)

        # -- who is A and who is B ---------------------------------------------
        header = QHBoxLayout()
        self.names_label = QLabel()
        self.names_label.setWordWrap(True)
        header.addWidget(self.names_label, stretch=1)
        self.swap_button = QPushButton("Swap A ↔ B")
        self.swap_button.setToolTip(
            "Which one is A decides the sign of a difference and the "
            "numerator of a ratio: LH − LV and LV − LH are the same "
            "map with opposite colours.")
        self.swap_button.clicked.connect(self.swap)
        header.addWidget(self.swap_button)
        layout.addLayout(header)

        # -- what the two measurements recorded --------------------------------
        self.metadata = QTableWidget(0, 3)
        self.metadata.setHorizontalHeaderLabels(["", "A", "B"])
        self.metadata.verticalHeader().setVisible(False)
        self.metadata.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.metadata.horizontalHeader().setSectionResizeMode(
            QHeaderView.Stretch)
        self.metadata.setMaximumHeight(150)
        self.metadata.setToolTip(
            "What each measurement recorded. Rows in orange differ by more "
            "than the tolerance -- a dichroism map between two sample "
            "positions, two temperatures or two photon energies is a map of "
            "that difference, not of the polarisation.")
        layout.addWidget(self.metadata)
        self.warnings = QLabel()
        self.warnings.setWordWrap(True)
        self.warnings.setStyleSheet("QLabel { color: #8a4b00; }")
        layout.addWidget(self.warnings)

        # -- the settings ------------------------------------------------------
        layout.addWidget(self._settings())

        # -- A, B and the result -----------------------------------------------
        self.preview = PreviewPair()
        self.preview.set_labels(*self.labels)
        self.preview.before["box"].setTitle("A")
        self.preview.after["box"].setTitle("B (after scaling, on A's grid)")
        layout.addWidget(self.preview, stretch=1)
        self.result_view = plain_image_view(*self.labels)
        self.result_view.setMinimumHeight(240)
        layout.addWidget(self.result_view, stretch=1)

        self.report = QLabel()
        self.report.setWordWrap(True)
        self.report.setTextInteractionFlags(Qt.TextSelectableByMouse)
        layout.addWidget(self.report)

        buttons = QDialogButtonBox()
        self.to_list = add_button(buttons,"Result to list",
                                         QDialogButtonBox.ActionRole)
        self.to_list.clicked.connect(self.export_result)
        self.sigma_to_list = add_button(buttons,"Uncertainty to list",
                                               QDialogButtonBox.ActionRole)
        self.sigma_to_list.setToolTip(
            "The propagated Poisson uncertainty of the result, as a cut of "
            "its own -- to mask with, or to show next to the result. Only "
            "available when both cuts are raw counts on the same grid.")
        self.sigma_to_list.clicked.connect(self.export_sigma)
        figure = add_button(buttons,"All three as a figure",
                                   QDialogButtonBox.ActionRole)
        figure.clicked.connect(self.to_figure)
        add_button(buttons,"Close", QDialogButtonBox.RejectRole).clicked.connect(
            self.reject)
        layout.addWidget(buttons)

        self._fill_names()
        self._fill_metadata()
        self.apply_preset()
        self.preview._apply_aspect()

    # -- the pair ----------------------------------------------------------
    def _set_pair(self, a, b):
        (self.name_a, self.a), (self.name_b, self.b) = a, b
        for name, data in (a, b):
            if getattr(data, "kind", None) != "cut":
                raise ValueError(
                    f"{name} is a {getattr(data, 'kind', 'dataset')}; cut "
                    f"arithmetic works on two cuts. Take a slice of a map "
                    f"first.")
        sa, sb = self.a.scan, self.b.scan
        self.a_axes = (np.asarray(sa.x, float), np.asarray(sa.y, float))
        self.b_axes = (np.asarray(sb.x, float), np.asarray(sb.y, float))
        self.a_values = np.asarray(self.a.cut_frame, dtype=float)
        self.b_values = np.asarray(self.b.cut_frame, dtype=float)
        self.labels = (sa.labels.get("x", "x"), sa.labels.get("y", "y"))
        self.b_labels = (sb.labels.get("x", "x"), sb.labels.get("y", "y"))

    def swap(self):
        self._set_pair((self.name_b, self.b), (self.name_a, self.a))
        self.preview.set_labels(*self.labels)
        self._fill_names()
        self._fill_metadata()
        self._region_defaults()
        self.refresh()

    def _fill_names(self):
        self.names_label.setText(
            f"<b>A</b> &nbsp;{self.name_a}<br><b>B</b> &nbsp;{self.name_b}")

    def _fill_metadata(self):
        rows = cutops.metadata_differences(self.a.scan.info or {},
                                           self.b.scan.info or {})
        self.metadata.setRowCount(len(rows))
        for index, (name, va, vb, differs) in enumerate(rows):
            cells = (QTableWidgetItem(name),
                     QTableWidgetItem("—" if va is None else _fmt(va)),
                     QTableWidgetItem("—" if vb is None else _fmt(vb)))
            # Polarisation is the thing that is *supposed* to differ in a
            # dichroism measurement, so it is shown but never flagged.
            if differs and name != "Polarisation":
                for cell in cells:
                    cell.setBackground(QColor("#ffe2b8"))
            for column, cell in enumerate(cells):
                self.metadata.setItem(index, column, cell)
        self.metadata.setVisible(bool(rows))
        self._differing = [name for name, _a, _b, differs in rows
                           if differs and name != "Polarisation"]

    # -- the settings ------------------------------------------------------
    def _settings(self) -> QWidget:
        box = QGroupBox("Operation")
        columns = QHBoxLayout(box)

        left = QFormLayout()
        self.preset = QComboBox()
        for key, spec in cutops.PRESETS.items():
            self.preset.addItem(spec["label"], key)
        self.preset.setToolTip(
            "A starting point for each measurement; every setting below "
            "stays editable.")
        self.preset.currentIndexChanged.connect(self.apply_preset)
        left.addRow("Purpose", self.preset)

        self.operation = QComboBox()
        for key, text in cutops.OPERATIONS.items():
            self.operation.addItem(text, key)
        self.operation.currentIndexChanged.connect(self._operation_changed)
        left.addRow("Result", self.operation)

        self.normalise = QComboBox()
        for key, text in cutops.NORMALISATIONS.items():
            self.normalise.addItem(text, key)
        self.normalise.setToolTip(
            "How B is scaled to A before they are combined. Two polarisations "
            "come off the undulator with different flux, so without this a "
            "linear-dichroism map is mostly the flux ratio.")
        self.normalise.currentIndexChanged.connect(self._normalise_changed)
        left.addRow("Scale B to A", self.normalise)
        columns.addLayout(left, stretch=1)

        middle = QFormLayout()
        self.region_boxes = []
        row = QHBoxLayout()
        for _ in range(4):
            spin = QDoubleSpinBox()
            spin.setDecimals(4)
            spin.setRange(-1e6, 1e6)
            spin.setKeyboardTracking(False)
            spin.setButtonSymbols(QDoubleSpinBox.NoButtons)
            spin.setMaximumWidth(78)
            spin.valueChanged.connect(self.refresh)
            self.region_boxes.append(spin)
            row.addWidget(spin)
        holder = QWidget()
        holder.setLayout(row)
        row.setContentsMargins(0, 0, 0, 0)
        middle.addRow("Region x0 y0 x1 y1", holder)
        self.from_box = QPushButton("Take the box from A's viewer")
        self.from_box.setToolTip(
            "Use the selection box on the cut viewer this was opened from. "
            "A region above E_F, or a band known not to be dichroic, makes a "
            "better reference than the whole cut.")
        self.from_box.clicked.connect(self._region_from_viewer)
        self.from_box.setEnabled(self.region_source is not None)
        middle.addRow("", self.from_box)
        columns.addLayout(middle, stretch=1)

        right = QFormLayout()
        self.reference_shape = QComboBox()
        for key, text in cutops.REFERENCE_SHAPES.items():
            self.reference_shape.addItem(text, key)
        self.reference_shape.setToolTip(
            "What of B is divided by. For a gold reference the angular "
            "profile is usual: it carries the detector's channel-to-channel "
            "sensitivity without printing gold's own Fermi edge and noise "
            "into the result.")
        self.reference_shape.currentIndexChanged.connect(self.refresh)
        right.addRow("Divide by", self.reference_shape)

        self.reference_floor = self._percent(5.0,
            "Hide the ratio wherever the reference is below this fraction "
            "of its maximum -- dividing by almost nothing is not a "
            "correction.")
        right.addRow("Reference floor", self.reference_floor)

        self.intensity_floor = self._percent(2.0,
            "Hide the result wherever A + B is below this fraction of its "
            "maximum. An asymmetry is a ratio, and in the background it is "
            "a ratio of two noises.")
        right.addRow("Hide where A+B below", self.intensity_floor)
        columns.addLayout(right, stretch=1)

        self._region_defaults()
        return box

    def _percent(self, value, tip):
        spin = QDoubleSpinBox()
        spin.setRange(0.0, 100.0)
        spin.setDecimals(1)
        spin.setSuffix(" %")
        spin.setValue(value)
        spin.setKeyboardTracking(False)
        spin.setToolTip(tip)
        spin.valueChanged.connect(self.refresh)
        return spin

    def _region_defaults(self):
        """The whole of A, until something better is chosen."""
        x, y = self.a_axes
        values = (float(x.min()), float(y.min()), float(x.max()), float(y.max()))
        for spin, value in zip(getattr(self, "region_boxes", ()), values):
            spin.blockSignals(True)
            spin.setValue(value)
            spin.blockSignals(False)

    def _region_from_viewer(self):
        corners = self.region_source() if self.region_source else None
        if corners is None:
            self.report.setText("A's viewer has no selection box on it. Turn "
                                "one on there, drag it, then press this again.")
            return
        for spin, value in zip(self.region_boxes, corners):
            spin.blockSignals(True)
            spin.setValue(float(value))
            spin.blockSignals(False)
        self.normalise.setCurrentIndex(self.normalise.findData("region"))
        self.refresh()

    def region(self):
        return tuple(float(spin.value()) for spin in self.region_boxes)

    # -- presets and dependencies ------------------------------------------
    def apply_preset(self, *_):
        spec = cutops.PRESETS[self.preset.currentData()]
        for combo, key in ((self.operation, "operation"),
                           (self.normalise, "normalise"),
                           (self.reference_shape, "reference_shape")):
            if key in spec:
                combo.blockSignals(True)
                combo.setCurrentIndex(combo.findData(spec[key]))
                combo.blockSignals(False)
        self._operation_changed()

    def _operation_changed(self, *_):
        ratio = self.operation.currentData() == "ratio"
        # A ratio divides by a normalised reference, so scaling B to A first
        # would change nothing; the reference options belong to it alone.
        self.normalise.setEnabled(not ratio)
        self.reference_shape.setEnabled(ratio)
        self.reference_floor.setEnabled(ratio)
        self.intensity_floor.setEnabled(not ratio)
        self._normalise_changed()

    def _normalise_changed(self, *_):
        regional = (self.normalise.currentData() == "region"
                    and self.normalise.isEnabled())
        for spin in self.region_boxes:
            spin.setEnabled(regional)
        self.refresh()

    # -- computing and drawing ---------------------------------------------
    def settings(self) -> dict:
        return dict(
            operation=self.operation.currentData(),
            normalise=self.normalise.currentData(),
            region=self.region(),
            reference_shape=self.reference_shape.currentData(),
            reference_floor=self.reference_floor.value() / 100.0,
            intensity_floor=self.intensity_floor.value() / 100.0)

    def compute(self) -> cutops.CombineResult:
        return cutops.combine(self.a_values, self.a_axes, self.b_values,
                              self.b_axes, a_labels=self.labels,
                              b_labels=self.b_labels, **self.settings())

    def refresh(self, *_):
        if not hasattr(self, "result_view"):
            return                  # still being built
        notes = cutops.polarisation_warnings(self.a.scan.info or {},
                                             self.b.scan.info or {},
                                             self.preset.currentData())
        if self._differing:
            notes.append("A and B differ in " + ", ".join(self._differing)
                         + " -- see the table above.")
        self.warnings.setText("<br>".join(notes))
        self.warnings.setVisible(bool(notes))

        try:
            result = self.compute()
        except (ValueError, ImportError) as exc:
            self.last = None
            self.report.setText(f"<span style='color:#b00'>{exc}</span>")
            self.result_view.clear()
            self.to_list.setEnabled(False)
            self.sigma_to_list.setEnabled(False)
            return
        self.last = result
        self.to_list.setEnabled(True)
        self.sigma_to_list.setEnabled(result.sigma is not None)
        self._draw_inputs(result)
        self._draw_result(result)

    def _draw_inputs(self, result):
        a, b = self.a_values, result.b_used
        pool = np.concatenate([a[np.isfinite(a)].ravel(),
                               b[np.isfinite(b)].ravel()])
        levels = np.percentile(pool, [1, 99]) if pool.size else (0.0, 1.0)
        for target, values in ((self.preview.before, a),
                               (self.preview.after, b)):
            self.preview._show(target, values, list(self.a_axes),
                               self.colormap, self.flip)
            target["view"].setLevels(float(levels[0]), float(levels[1]))
            target["caption"].setText(
                f"{values.shape[0]} × {values.shape[1]}   "
                f"common scale [{levels[0]:.4g} .. {levels[1]:.4g}]")

    def _draw_result(self, result):
        values = result.values
        finite = values[np.isfinite(values)]
        operation = self.operation.currentData()
        signed = operation in ("difference", "asymmetry")
        # A signed result wants a window symmetric about zero and a diverging
        # map: an asymmetric window makes a small positive change look like a
        # large negative one.
        if not finite.size:
            lo, hi = -1.0, 1.0
        elif signed:
            reach = float(np.percentile(np.abs(finite), 98)) or 1.0
            lo, hi = -reach, reach
        else:
            lo, hi = (float(v) for v in np.percentile(finite, [2, 98]))
        x, y = self.a_axes
        show_frame(self.result_view, values, x, y, levels=(lo, hi))
        apply_colormap(self.result_view,
                       "redwhiteblue" if signed else self.colormap,
                       False if signed else self.flip)
        fit_frame_view(self.result_view, x, y)

        text = [f"<b>{self.operation.currentText()}</b>: shown on "
                f"[{lo:.4g} .. {hi:.4g}]"]
        summary = result.summary()
        if summary:
            text += summary.split("\n")
        self.report.setText("<br>".join(text))

    # -- leaving -----------------------------------------------------------
    def _result_name(self, suffix: str) -> str:
        base = f"{self.name_a} {suffix}"
        taken = set(self.existing_names() or ())
        if base not in taken:
            return base
        index = 2
        while f"{base} ({index})" in taken:
            index += 1
        return f"{base} ({index})"

    def _parameters(self, result) -> dict:
        settings = self.settings()
        operation = settings["operation"]
        params = {"a": self.name_a, "b": self.name_b,
                  "preset": self.preset.currentData(),
                  "operation": operation,
                  "scale_applied_to_b": float(result.scale),
                  "b_resampled": bool(result.resampled),
                  "overlap": float(result.overlap)}
        if operation == "ratio":
            params.update(reference_shape=settings["reference_shape"],
                          reference_floor=settings["reference_floor"])
        else:
            params.update(normalise=settings["normalise"],
                          intensity_floor=settings["intensity_floor"])
            if settings["normalise"] == "region":
                params["region"] = list(settings["region"])
        if result.median_sigma is not None:
            params["median_sigma"] = float(result.median_sigma)
            params["fraction_beyond_2sigma"] = float(result.significant)
        return params

    def _dataset(self, values, suffix, params):
        info = P.record_step(dict(self.a.scan.info or {}),
                             P.Step("cut_arithmetic", params,
                                    source=f"{self.name_a} ∘ {self.name_b}"))
        return MemoryData(
            "cut", tuple(self.a_axes), values, dict(self.a.scan.labels),
            source_label=self._result_name(suffix), parameters=params,
            prefix="cutops", source_path=getattr(self.a, "path", ""),
            source_info=info)

    def export_result(self):
        result = self.last
        if result is None:
            return None
        preset = self.preset.currentData()
        suffix = _SUFFIX.get(preset if preset != "custom" else
                             self.operation.currentData(), "calc")
        data = self._dataset(result.values, suffix, self._parameters(result))
        self.datasetsCreated.emit([data])
        self.report.setText(self.report.text()
                            + f"<br><b>Added “{data.source_label}” "
                              f"to the list.</b>")
        return data

    def export_sigma(self):
        result = self.last
        if result is None or result.sigma is None:
            return None
        params = self._parameters(result)
        params["quantity"] = "Poisson standard deviation of the result"
        data = self._dataset(result.sigma, "σ", params)
        self.datasetsCreated.emit([data])
        self.report.setText(self.report.text()
                            + f"<br><b>Added “{data.source_label}” "
                              f"to the list.</b>")
        return data

    def to_figure(self):
        from ui.figure import FigureWindow, panel_from_arrays

        result = self.last
        if result is None:
            return None
        x, y = self.a_axes
        panels = [
            panel_from_arrays(self.a_values, x, y, *self.labels, self.name_a),
            panel_from_arrays(result.b_used, x, y, *self.labels,
                              self.name_b),
            panel_from_arrays(result.values, x, y, *self.labels,
                              self.operation.currentText()),
        ]
        if self.operation.currentData() in ("difference", "asymmetry"):
            panels[2].colormap = "redwhiteblue"
        window = FigureWindow(panels, parent=self)
        window.figure.cols = 3
        window.show()
        return window


def _fmt(value) -> str:
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)
