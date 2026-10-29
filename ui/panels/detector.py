"""
ui/panels/detector.py
=====================
The detector corrections, step by step: a trapezoidal (keystone) image made
rectangular, and dead or bad regions patched -- on a cut, or on every
detector image of a map at once.

Reached from the tool panel (Basic -> Detector corrections..., under Remove
detector grid) with one cut, map or kz map selected that is still on the
detector's pixels. Three steps:

1. **Keystone** -- drag the four corners of the trapezoid the lit area
   forms; the preview shows it stretched to a rectangle.
2. **Bad regions** -- draw rectangles over dead or hot spots (or let them
   be found), and choose how they are filled.
3. **Apply** -- the corrected dataset into the list.

The corrections are :mod:`tools.processing.detector`, which imports no Qt.
"""
from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PyQt5.QtCore import pyqtSignal
from PyQt5.QtWidgets import (QCheckBox, QComboBox, QDialog, QDoubleSpinBox,
                             QFormLayout, QHBoxLayout, QLabel, QLineEdit,
                             QListWidget, QPushButton, QVBoxLayout, QWidget)

from loader.nxs_file import axis_slots
from tools.processing import detector as DT
from tools.processing import process as P
from ui.common.screen import fits_screen
from ui.common.stepper import Step, Stepper, guard_window
from ui.panels.process import PreviewPair
from ui.viewers.widgets import MemoryData

__all__ = ["DetectorDialog"]


@fits_screen
class DetectorDialog(QDialog):
    """Keystone and bad-region corrections. ``dataset`` is ``(name,
    data)``: a cut, map or kz map as measured."""

    datasetsCreated = pyqtSignal(list)

    def __init__(self, dataset, parent=None, colormap="Greys",
                 existing_names=None):
        super().__init__(parent)
        self.name, self.data = dataset
        reason = DT.not_on_detector(dict(self.data.scan.info or {}))
        if self.data.kind not in ("cut", "map", "kz_map"):
            reason = reason or "it is in momentum, resampled off the detector's pixels"
        if reason:
            raise ValueError(
                f"{self.name} cannot be corrected: {reason}. These corrections "
                f"belong on the detector image as measured -- correct the "
                f"original, then do the rest.")
        self.setWindowTitle(f"Detector corrections -- {self.name}")
        self.resize(1080, 720)
        self.colormap = colormap
        self._existing = existing_names or (lambda: [])
        self.result = None
        self._cache = {}

        slots = axis_slots(self.data.kind, "constructor")
        labels = dict(self.data.scan.labels or {})
        self.axes = [np.asarray(getattr(self.data.scan, s), dtype=float)
                     for s in slots]
        fx, fy = DT.detector_frames(self.data.kind)
        self.fx, self.fy = self.axes[fx], self.axes[fy]
        self.labels = (labels.get(slots[fx], slots[fx]),
                       labels.get(slots[fy], slots[fy]))
        self.frame = self._mean_frame()

        layout = QVBoxLayout(self)
        self.stepper = Stepper([self._step_keystone(), self._step_regions(),
                                self._step_apply()], self,
                               finish_text="Apply to list")
        self.stepper.finished.connect(lambda: self.apply_to_list())
        layout.addWidget(self.stepper)
        guard_window(self, self.stepper, "detector", "Detector corrections")
        self.update_keystone()
        self.update_regions()

    def _values(self):
        return (self.data.cut_frame if self.data.kind == "cut"
                else self.data.scan.value)

    def _mean_frame(self):
        """The detector image: the cut itself, or a map's slices averaged
        (one slice at a time, so a map read lazily stays on disk)."""
        values = self._values()
        if self.data.kind == "cut":
            return np.asarray(values, dtype=float)
        total = None
        for i in range(values.shape[0]):
            plane = np.asarray(values[i], dtype=float)
            total = plane if total is None else total + plane
        return total / max(1, values.shape[0])

    # ==================================================================
    # Step 1: keystone
    # ==================================================================
    def _step_keystone(self):
        page = QWidget()
        column = QVBoxLayout(page)
        row = QHBoxLayout()
        self.keystone_box = QCheckBox("Correct a trapezoidal image")
        self.keystone_box.toggled.connect(self.update_keystone)
        row.addWidget(self.keystone_box)
        row.addWidget(QLabel("Stretch along"))
        self.along_combo = QComboBox()
        self.along_combo.addItems([f"{self.labels[0]} (edges change with "
                                   f"{self.labels[1]})",
                                   f"{self.labels[1]} (edges change with "
                                   f"{self.labels[0]})"])
        self.along_combo.currentIndexChanged.connect(self._reset_trapezoid)
        row.addWidget(self.along_combo)
        self.conserve_box = QCheckBox("Keep each line's total counts")
        self.conserve_box.setChecked(True)
        self.conserve_box.toggled.connect(self.update_keystone)
        row.addWidget(self.conserve_box)
        reset = QPushButton("Reset corners")
        reset.clicked.connect(self._reset_trapezoid)
        row.addWidget(reset)
        row.addStretch(1)
        column.addLayout(row)
        self.keystone_preview = PreviewPair(page)
        self.keystone_preview.set_labels(*self.labels)
        column.addWidget(self.keystone_preview, 1)
        self.keystone_label = QLabel("")
        column.addWidget(self.keystone_label)
        self.trapezoid = pg.PolyLineROI([[0, 0], [1, 0], [1, 1], [0, 1]],
                                        closed=True, pen=pg.mkPen("#ffcc00", width=2))
        self.keystone_preview.before["view"].view.addItem(self.trapezoid)
        self.trapezoid.sigRegionChangeFinished.connect(self.update_keystone)
        self._reset_trapezoid()
        text = ("Some lens modes light a <b>trapezoid</b> of the detector "
                "rather than a rectangle: the angular range shrinks or grows "
                "with the kinetic energy. Tick <b>Correct</b>, then drag the "
                "four yellow corners onto the lit area's corners in the "
                "picture on the left. The right shows it stretched so every "
                "line spans the extent of the trapezoid's first edge. Leave "
                "it unticked if the image is already rectangular.")
        return Step("Keystone", text, page, ready=self._keystone_problem)

    def _reset_trapezoid(self, *_):
        x0, x1 = float(self.fx.min()), float(self.fx.max())
        y0, y1 = float(self.fy.min()), float(self.fy.max())
        mx, my = 0.15 * (x1 - x0), 0.15 * (y1 - y0)
        if self.along_combo.currentIndex() == 0:
            corners = [[x0 + mx, y1 - my], [x1 - mx, y1 - my],
                       [x1 - 2 * mx, y0 + my], [x0 + 2 * mx, y0 + my]]
        else:
            corners = [[x1 - mx, y0 + my], [x1 - mx, y1 - my],
                       [x0 + mx, y1 - 2 * my], [x0 + mx, y0 + 2 * my]]
        self.set_corners(corners)

    def set_corners(self, corners):
        self.trapezoid.blockSignals(True)
        self.trapezoid.setPos(0, 0)
        self.trapezoid.setPoints([list(map(float, c)) for c in corners],
                                 closed=True)
        self.trapezoid.blockSignals(False)
        self.update_keystone()

    def corners(self):
        return [(p.x(), p.y()) for p in
                (self.trapezoid.mapToParent(pos) for _name, pos in
                 self.trapezoid.getLocalHandlePositions())]

    def keystone(self):
        if not self.keystone_box.isChecked():
            return None
        return DT.Keystone.from_corners(self.corners(),
                                        along=self.along_combo.currentIndex())

    def _keystone_problem(self):
        try:
            self.keystone()
        except ValueError as exc:
            return str(exc)
        return ""

    def update_keystone(self, *_):
        preview = self.keystone_preview
        preview.show_source(self.frame, (self.fx, self.fy), self.colormap, False)
        problem = self._keystone_problem()
        keystone = None if problem else self.keystone()
        if keystone is None:
            preview.show_result(self.frame, (self.fx, self.fy), self.colormap, False)
            self.keystone_label.setText(problem or "No keystone correction.")
        else:
            out = DT.keystone_correct(self.frame, self.fx, self.fy, keystone,
                                      conserve=self.conserve_box.isChecked())
            preview.show_result(out, (self.fx, self.fy), self.colormap, False)
            s0, e0 = keystone.start_a, keystone.end_a
            s1, e1 = keystone.start_b, keystone.end_b
            self.keystone_label.setText(
                f"Edge at {keystone.level_a:.4g}: {s0:.4g} .. {e0:.4g}; at "
                f"{keystone.level_b:.4g}: {s1:.4g} .. {e1:.4g} -> stretched to "
                f"{s0:.4g} .. {e0:.4g} everywhere.")
        self.trapezoid.setVisible(self.keystone_box.isChecked())
        if hasattr(self, "stepper"):
            self.stepper.refresh()

    # ==================================================================
    # Step 2: bad regions
    # ==================================================================
    def _step_regions(self):
        page = QWidget()
        column = QVBoxLayout(page)
        row = QHBoxLayout()
        add = QPushButton("Add a region")
        add.setToolTip("A rectangle in the middle of the view: drag it over "
                       "the bad spot, and its corner to size it.")
        add.clicked.connect(lambda: self.add_region())
        row.addWidget(add)
        find = QPushButton("Find bad pixels")
        find.setToolTip("Pixels that differ from their 5 x 5 neighbourhood's "
                        "median by more than the threshold (on a map: in the "
                        "image averaged over the scan, where the "
                        "photoemission is smeared out and a faulty pixel is "
                        "not).")
        find.clicked.connect(self.find_regions)
        row.addWidget(find)
        self.find_threshold = QDoubleSpinBox()
        self.find_threshold.setRange(5.0, 500.0)
        self.find_threshold.setValue(30.0)
        self.find_threshold.setSuffix(" % off")
        row.addWidget(self.find_threshold)
        remove = QPushButton("Remove selected")
        remove.clicked.connect(self.remove_selected_region)
        row.addWidget(remove)
        clear = QPushButton("Clear")
        clear.clicked.connect(self.clear_regions)
        row.addWidget(clear)
        row.addWidget(QLabel("Fill"))
        self.method_combo = QComboBox()
        self.method_combo.addItems(DT.METHODS)
        self.method_combo.setCurrentIndex(2)
        self.method_combo.setToolTip(
            "linear along x / y: blend the pixels just outside the region on "
            "either side (the MATLAB tool's way).\nsmooth: the smoothest "
            "surface meeting the data all round the region.")
        self.method_combo.currentIndexChanged.connect(self.update_regions)
        row.addWidget(self.method_combo)
        row.addStretch(1)
        column.addLayout(row)
        body = QHBoxLayout()
        self.region_preview = PreviewPair(page)
        self.region_preview.set_labels(*self.labels)
        body.addWidget(self.region_preview, 4)
        self.region_list = QListWidget()
        self.region_list.setMaximumWidth(230)
        body.addWidget(self.region_list, 1)
        column.addLayout(body, 1)
        self.region_rois = []
        text = ("Dead or hot spots of the detector -- a dead channel, a "
                "burnt-in spot, a damaged MCP region -- are in every image "
                "it takes. <b>Add a region</b> and drag it over one, or press "
                "<b>Find bad pixels</b> to have them found. They are filled "
                "from their surroundings: <b>smooth</b> is the best guess "
                "for a spot, <b>linear along</b> an axis for a stripe that "
                "runs across the other. Skip this step if there is nothing "
                "to patch.")
        return Step("Bad regions", text, page, enter=self.update_regions)

    def add_region(self, region: DT.Region = None):
        if region is None:
            x0, x1 = float(self.fx.min()), float(self.fx.max())
            y0, y1 = float(self.fy.min()), float(self.fy.max())
            w, h = 0.08 * (x1 - x0), 0.08 * (y1 - y0)
            region = DT.Region((x0 + x1) / 2 - w / 2, (x0 + x1) / 2 + w / 2,
                               (y0 + y1) / 2 - h / 2, (y0 + y1) / 2 + h / 2)
        roi = pg.RectROI([region.x0, region.y0],
                         [region.x1 - region.x0, region.y1 - region.y0],
                         pen=pg.mkPen("#ff5050", width=2))
        roi.sigRegionChangeFinished.connect(self.update_regions)
        self.region_preview.before["view"].view.addItem(roi)
        self.region_rois.append(roi)
        self.update_regions()
        return roi

    def regions(self):
        out = []
        for roi in self.region_rois:
            pos, size = roi.pos(), roi.size()
            out.append(DT.Region(pos.x(), pos.x() + size.x(),
                                 pos.y(), pos.y() + size.y()))
        return out

    def find_regions(self):
        found = DT.find_bad_regions(self.frame, self.fx, self.fy,
                                    threshold=self.find_threshold.value() / 100.0)
        for region in found:
            self.add_region(region)
        self.stepper.missing.setText(f"Found {len(found)} bad region(s)." if found
                                     else "No pixel stands out at this threshold.")
        return found

    def remove_selected_region(self):
        row = self.region_list.currentRow()
        if 0 <= row < len(self.region_rois):
            roi = self.region_rois.pop(row)
            self.region_preview.before["view"].view.removeItem(roi)
            self.update_regions()

    def clear_regions(self):
        for roi in self.region_rois:
            self.region_preview.before["view"].view.removeItem(roi)
        self.region_rois = []
        self.update_regions()

    def update_regions(self, *_):
        preview = self.region_preview
        preview.show_source(self.frame, (self.fx, self.fy), self.colormap, False)
        regions = self.regions()
        self.region_list.clear()
        for r in regions:
            self.region_list.addItem(f"{r.x0:.4g}..{r.x1:.4g}, {r.y0:.4g}..{r.y1:.4g}")
        try:
            out = DT.patch_regions(self.frame, self.fx, self.fy, regions,
                                   self.method_combo.currentText(), self._cache)
        except ValueError as exc:
            out = self.frame
            self.stepper.missing.setText(str(exc)) if hasattr(self, "stepper") else None
        preview.show_result(out, (self.fx, self.fy), self.colormap, False)

    # ==================================================================
    # Step 3: apply
    # ==================================================================
    def _step_apply(self):
        page = QWidget()
        form = QFormLayout(page)
        self.name_edit = QLineEdit(self._default_name())
        form.addRow("Name", self.name_edit)
        self.summary = QLabel("")
        self.summary.setWordWrap(True)
        form.addRow(self.summary)
        what = "every slice of the map" if self.data.kind != "cut" else "the cut"
        text = (f"The corrections are applied to {what}: the bad regions "
                f"patched first, on the pixels as they were, then the "
                f"keystone stretched out. <b>Apply to list</b> adds the "
                f"corrected dataset; do de-gridding, the Fermi level and the "
                f"k conversion on it afterwards.")
        return Step("Apply", text, page, enter=self._summary,
                    ready=lambda: "" if (self.keystone_box.isChecked()
                                         or self.region_rois)
                    else "Nothing to apply: set a keystone or add a region.")

    def _default_name(self):
        taken = set(self._existing())
        name, n = f"{self.name}_det", 2
        while name in taken:
            name, n = f"{self.name}_det{n}", n + 1
        return name

    def _summary(self):
        parts = []
        if self.keystone_box.isChecked():
            parts.append("keystone corrected")
        if self.region_rois:
            parts.append(f"{len(self.region_rois)} region(s) patched "
                         f"({self.method_combo.currentText()})")
        self.summary.setText(("; ".join(parts) or "nothing set") + ".")

    def build(self, progress=None):
        keystone = self.keystone()
        regions = self.regions()
        method = self.method_combo.currentText()
        values = DT.apply_to(self._values(), self.axes, self.data.kind,
                             keystone=keystone, regions=regions, method=method,
                             conserve=self.conserve_box.isChecked(),
                             progress=progress)
        parameters = {"regions": np.asarray([r.as_tuple() for r in regions],
                                            dtype=float).reshape(-1, 4),
                      "fill": method}
        if keystone is not None:
            parameters.update({f"keystone_{k}": v
                               for k, v in keystone.as_dict().items()})
            parameters["keystone_conserve"] = bool(self.conserve_box.isChecked())
        info = P.record_step(dict(self.data.scan.info or {}), P.Step(
            "detector", {"keystone": keystone is not None,
                         "regions": len(regions), "fill": method},
            source=self.name))
        return MemoryData(self.data.kind, tuple(self.axes), np.nan_to_num(values),
                          dict(self.data.scan.labels or {}),
                          source_label=self.name_edit.text().strip()
                          or self._default_name(),
                          parameters=parameters, prefix="detcorr",
                          source_path=getattr(self.data, "path", ""),
                          source_info=info,
                          source_motors=dict(getattr(self.data.scan,
                                                     "fourd_info", {}) or {}))

    def apply_to_list(self, background: bool = True):
        if self.data.kind != "cut" and background:
            from ui.common import jobs
            try:
                self.result = jobs.run_blocking(self, "Correcting the map",
                                                lambda report: self.build(report))
            except jobs.JobCancelled:
                return None
        else:
            self.result = self.build()
        self.datasetsCreated.emit([self.result])
        self.stepper.missing.setText(f"Added {self.result.source_label} to the list.")
        return self.result
