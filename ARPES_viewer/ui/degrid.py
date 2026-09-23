"""
ui/degrid.py
============
The De-grid window: remove the detector's grid from a cut or a whole map.

Opened from a cut viewer (**De-grid...**, that cut only) or from a map's
contour or slit-cut window (**De-grid map...**, the whole map -- angle maps
and kz maps alike). The defaults are the configuration that did best on the
real data it was developed on; the Advanced box is there for when it does
not, and can stay closed. :mod:`tools.degrid` has the method and the
numbers.

What the window shows, so the result can be judged rather than trusted:

* the slice before and after, on the same colour scale, optionally as
  −∂²I/∂E² -- a second derivative multiplies a grid by its frequency
  squared, which makes it the most sensitive view there is;
* the grid pattern that was found, and where it sits in k-space;
* for a map, the contrast and the displacement fitted to every slice;
* the grid-to-background power in the grid regions before and after.

The result goes to the main list as a new dataset. For a map, the grid
pattern can go there too (``[grid]``): a cut taken later with the same
detector settings can then be de-gridded with it, with no reference
measurement.
"""
from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtWidgets import (QCheckBox, QComboBox, QDialog, QDialogButtonBox,
                             QDoubleSpinBox, QFormLayout, QGridLayout,
                             QGroupBox, QHBoxLayout, QLabel,
                             QPushButton, QSpinBox, QVBoxLayout, QWidget)
from scipy.ndimage import gaussian_filter

from tools import degrid as DG
from tools import process as P
from ui.widgets import MemoryData, apply_colormap, plain_image_view, show_frame

__all__ = ["DegridDialog", "grid_candidates", "default_source"]

INTRO_MAP = (
    "Finds the detector's grid from the whole map -- it stays on the same "
    "pixels while the photoemission moves from slice to slice -- and divides "
    "it out of every slice, with each slice's own contrast and sub-pixel "
    "shift fitted. The defaults are the best found on real maps; they "
    "should not need changing.")
INTRO_CUT = (
    "Removes the detector's grid from this cut. Best: a grid measured on a "
    "map taken with the same detector settings (de-grid that map first and "
    "keep its [grid]). Otherwise the cut's own grid peaks are notched out, "
    "which also removes a little of the photoemission at the same "
    "frequencies.")


def grid_candidates(entries, loader, frame_shape, info):
    """``[(label, pattern, box, matches)]`` for every ``[grid]`` dataset in
    the list with this detector frame. ``matches`` says whether its lens
    mode and pass energy are this cut's as well."""
    out = []
    for label, key in (entries() if callable(entries) else []):
        if "[grid]" not in str(label):
            continue
        try:
            data = loader(key)
        except Exception:                                   # noqa: BLE001
            continue
        scan = data.scan
        if not scan.info.get("degrid.is_grid"):
            continue
        pattern = np.asarray(scan.value, dtype=float)
        if pattern.shape != tuple(frame_shape):
            continue
        box = _parse_box(scan.info.get("degrid.box"))
        # the list's labels end in "   [Kind]"; the name is what goes on record
        name = str(label).split("   [")[0].strip()
        out.append((name, pattern, box, settings_match(scan.info, info)))
    return out


#: Where each reader records the lens mode and the pass energy.
_SETTINGS_KEYS = ("lens_mode", "pass_energy_eV", "MBS.lens_mode",
                  "MBS.passenergy", "scienta.Lens Mode", "scienta.Pass Energy",
                  "mbs.Lens Mode", "mbs.Pass Energy")


def settings_match(a: dict, b: dict) -> bool:
    """Whether two measurements used the same lens mode and pass energy, as
    far as their metadata says. With nothing to compare, the answer is no:
    a grid from unknown settings should be offered, not chosen."""
    compared = False
    for key in _SETTINGS_KEYS:
        va, vb = (a or {}).get(key), (b or {}).get(key)
        if va is None or vb is None:
            continue
        compared = True
        if str(va).strip() != str(vb).strip():
            return False
    return compared


def default_source(candidates) -> int:
    """Index into ``["notch"] + candidates`` to start on: a grid from the
    same detector settings if there is one, the notch otherwise. A grid
    from other settings is offered but not chosen: on a real pair (pass
    energy 20 vs 50 eV) it removed ~88 % of the grid's power, where the
    notch removes all of it."""
    for i, candidate in enumerate(candidates):
        if candidate[3]:
            return i + 1
    return 0


def _parse_box(text):
    try:
        a, e = str(text).split(",")
        a0, a1 = (int(v) for v in a.split(":"))
        e0, e1 = (int(v) for v in e.split(":"))
        return slice(a0, a1), slice(e0, e1)
    except (ValueError, AttributeError):
        return None


def _box_text(box):
    return f"{box[0].start}:{box[0].stop},{box[1].start}:{box[1].stop}"


def second_derivative(image, sigma=1.5):
    """−∂²I/∂E², energy along the second axis."""
    smooth = gaussian_filter(np.asarray(image, dtype=float), sigma)
    return -np.gradient(np.gradient(smooth, axis=1), axis=1)


class DegridDialog(QDialog):
    """``window`` is the viewer it came from; ``data`` a cut or a map."""

    datasetsCreated = pyqtSignal(list)

    def __init__(self, window, data, name, *, candidates=(), colormap="gray",
                 flip=False, existing_names=None):
        super().__init__(window)
        self.window_ = window
        self.data = data
        self.name = name
        self.is_map = data.kind != "cut"
        self.colormap, self.flip = colormap, flip
        self.existing_names = existing_names or (lambda: [])
        self.candidates = list(candidates)
        self.result = None
        self.setWindowTitle(f"De-grid — {name}")
        self.setModal(False)
        self.resize(1150, 950)

        scan = data.scan
        if self.is_map:
            self.axes = (np.asarray(scan.x, float), np.asarray(scan.k, float),
                         np.asarray(scan.z, float))
            self.labels = (scan.labels.get("x", "x"), scan.labels.get("k", "angle"),
                           scan.labels.get("z", "energy"))
        else:
            self.axes = (np.asarray(scan.x, float), np.asarray(scan.y, float))
            self.labels = (scan.labels.get("x", "angle"), scan.labels.get("y", "energy"))

        layout = QVBoxLayout(self)
        intro = QLabel(INTRO_MAP if self.is_map else INTRO_CUT)
        intro.setWordWrap(True)
        layout.addWidget(intro)

        # -- what to use and how ----------------------------------------------
        top = QHBoxLayout()
        self.source = QComboBox()
        if not self.is_map:
            self.source.addItem("This cut alone (FFT notch)")
            for label, _p, _b, same in self.candidates:
                self.source.addItem(f"{label}" + ("" if same else
                                    "   (settings differ or unknown)"))
            self.source.setCurrentIndex(default_source(self.candidates))
            top.addWidget(QLabel("Grid from"))
            top.addWidget(self.source, stretch=1)
        self.run_button = QPushButton("Find and remove the grid")
        self.run_button.clicked.connect(lambda: self.run())
        top.addWidget(self.run_button)
        layout.addLayout(top)

        self.advanced = QGroupBox("Advanced (the defaults are the tested best)")
        self.advanced.setCheckable(True)
        self.advanced.setChecked(False)
        form = QFormLayout()
        adv_widget = QWidget()
        adv_widget.setLayout(form)
        box_layout = QVBoxLayout(self.advanced)
        box_layout.addWidget(adv_widget)
        self.advanced.toggled.connect(adv_widget.setVisible)
        adv_widget.setVisible(False)
        d = DG.DEFAULTS
        self.tiles = QSpinBox(); self.tiles.setRange(0, 12); self.tiles.setValue(d["tiles"])
        self.tiles.setToolTip("Local refinement of the contrast and shift over "
                              "n × n tiles; 0 fits one of each per slice.")
        self.iterations = QSpinBox(); self.iterations.setRange(0, 4)
        self.iterations.setValue(d["iterations"])
        self.iterations.setToolTip("Passes that re-estimate the map's grid with "
                                   "every slice shifted back into register.")
        self.smooth = QDoubleSpinBox(); self.smooth.setRange(2, 50); self.smooth.setValue(d["smooth"])
        self.smooth.setToolTip("Gaussian σ (pixels) of the smooth image the grid "
                               "is measured against. Must be well above the "
                               "grid's period.")
        self.threshold = QDoubleSpinBox(); self.threshold.setRange(3, 1000)
        self.threshold.setValue(d["seed_threshold"])
        self.threshold.setToolTip("How far above its surroundings a k-space "
                                  "peak must stand to be called grid.")
        form.addRow("Local tiles per side", self.tiles)
        if self.is_map:
            form.addRow("Registration passes", self.iterations)
        form.addRow("Smoothing σ (px)", self.smooth)
        form.addRow("Peak threshold", self.threshold)
        layout.addWidget(self.advanced)

        # -- views ---------------------------------------------------------------
        view_bar = QHBoxLayout()
        self.derivative = QCheckBox("Show −∂²I/∂E² (makes any grid obvious)")
        self.derivative.toggled.connect(lambda *_: self.redraw())
        view_bar.addWidget(self.derivative)
        self.slice_box = QSpinBox()
        self.slice_box.valueChanged.connect(lambda *_: self.redraw())
        if self.is_map:
            self.slice_box.setRange(0, int(np.asarray(data.angle_cube[3].shape)[0]) - 1)
            view_bar.addWidget(QLabel("Slice"))
            view_bar.addWidget(self.slice_box)
            self.slice_label = QLabel("")
            view_bar.addWidget(self.slice_label)
        view_bar.addStretch(1)
        layout.addLayout(view_bar)

        grid = QGridLayout()
        xl, yl = (self.labels[1], self.labels[2]) if self.is_map else self.labels
        self.before_view = plain_image_view(xl, yl)
        self.after_view = plain_image_view(xl, yl)
        self.pattern_view = plain_image_view("pixel", "pixel")
        self.kspace_view = plain_image_view("cycles/px (angle)", "cycles/px (energy)")
        for view, (r, c), title in ((self.before_view, (0, 0), "Before"),
                                    (self.after_view, (0, 1), "After"),
                                    (self.pattern_view, (2, 0), "The grid found (central 160 × 160 px)"),
                                    (self.kspace_view, (2, 1), "k-space: grid regions in red")):
            caption = QLabel(f"<b>{title}</b>")
            grid.addWidget(caption, r, c)
            grid.addWidget(view, r + 1, c)
            view.setMinimumHeight(220)
        layout.addLayout(grid, stretch=3)

        self.fit_plot = pg.PlotWidget()
        self.fit_plot.setBackground("w")
        self.fit_plot.setMaximumHeight(170)
        self.fit_plot.showGrid(x=True, y=True, alpha=0.25)
        self.fit_plot.addLegend(offset=(-10, 5))
        self.fit_plot.setLabel("bottom", "slice")
        self.fit_plot.setVisible(self.is_map)
        layout.addWidget(self.fit_plot, stretch=1)

        self.report = QLabel("Press “Find and remove the grid”.")
        self.report.setWordWrap(True)
        self.report.setTextInteractionFlags(Qt.TextSelectableByMouse)
        layout.addWidget(self.report)

        buttons = QDialogButtonBox()
        self.to_list = buttons.addButton("Result to list", QDialogButtonBox.ActionRole)
        self.to_list.clicked.connect(self.export)
        self.to_list.setEnabled(False)
        self.with_grid = QCheckBox("Also list the grid pattern (for cuts taken "
                                   "later with the same settings)")
        self.with_grid.setChecked(self.is_map)
        self.with_grid.setVisible(self.is_map)
        buttons.addButton("Close", QDialogButtonBox.RejectRole).clicked.connect(self.reject)
        bottom = QHBoxLayout()
        bottom.addWidget(self.with_grid)
        bottom.addStretch(1)
        bottom.addWidget(buttons)
        layout.addLayout(bottom)

    # -- running -------------------------------------------------------------------
    def settings(self) -> dict:
        return {"tiles": self.tiles.value(), "iterations": self.iterations.value(),
                "smooth": self.smooth.value(),
                "seed_threshold": self.threshold.value()}

    def cube(self):
        return self.data.angle_cube[3]

    def run(self, background: bool = True):
        """Compute. ``background=False`` runs on this thread (tests)."""
        s = self.settings()
        try:
            if self.is_map:
                cube = self.cube()

                def work(report=None):
                    def progress(done, total, label):
                        if report is not None:
                            report(done / max(total, 1), f"{label} ({done} of {total})")
                    return DG.degrid_map(cube, settings=s, progress=progress)
                if background:
                    from ui.jobs import run_blocking, JobCancelled
                    try:
                        result = run_blocking(self, "De-gridding the map", work)
                    except JobCancelled:
                        self.report.setText("Cancelled.")
                        return None
                else:
                    result = work()
            else:
                frame = np.asarray(self.data.cut_frame, dtype=float)
                index = self.source.currentIndex()
                if index <= 0:
                    result = DG.degrid_cut_notch(frame, settings=s)
                else:
                    _label, pattern, box, _same = self.candidates[index - 1]
                    result = DG.degrid_cut_with_grid(frame, pattern, box=box, settings=s)
        except DG.GridNotFound as exc:
            self.result = None
            self.to_list.setEnabled(False)
            self.report.setText(
                f"<b>No grid found</b> ({exc}). Nothing to remove -- or, in a "
                f"single cut with few counts, a grid hidden in the noise: "
                f"de-grid a map taken with the same settings and use its [grid].")
            return None
        except ValueError as exc:
            self.result = None
            self.to_list.setEnabled(False)
            self.report.setText(f"<span style='color:#b00'>{exc}</span>")
            return None
        self.result = result
        self.to_list.setEnabled(True)
        if self.is_map:
            self.slice_box.blockSignals(True)
            self.slice_box.setValue(result.preview_index)
            self.slice_box.blockSignals(False)
        self.report.setText(result.summary().replace("\n", "<br>"))
        self.redraw()
        self._draw_fits()
        return result

    # -- drawing -------------------------------------------------------------------
    def _images(self):
        if self.is_map:
            i = self.slice_box.value()
            before = np.asarray(self.cube()[i], dtype=float)
            after = (np.asarray(self.result.values[i], dtype=float)
                     if self.result is not None else None)
            self.slice_label.setText(f"{self.labels[0]} = {self.axes[0][i]:.4g}")
            axes = (self.axes[1], self.axes[2])
        else:
            before = np.asarray(self.data.cut_frame, dtype=float)
            after = (np.asarray(self.result.values, dtype=float)
                     if self.result is not None else None)
            axes = self.axes
        return before, after, axes

    def redraw(self):
        before, after, (x, y) = self._images()
        if self.derivative.isChecked():
            before = second_derivative(before)
            after = second_derivative(after) if after is not None else None
        finite = before[np.isfinite(before)]
        lo, hi = (np.percentile(finite, [1, 99.5]) if finite.size else (0, 1))
        if self.derivative.isChecked():
            lo = 0.0
        for view, image in ((self.before_view, before), (self.after_view, after)):
            if image is None:
                view.clear()
                continue
            show_frame(view, image, x, y, levels=(float(lo), float(hi)))
            apply_colormap(view, self.colormap, self.flip)
        if self.result is not None:
            self._draw_grid()

    def _draw_grid(self):
        model = self.result.model
        g = model.pattern()
        c0, c1 = g.shape[0] // 2, g.shape[1] // 2
        patch = g[max(c0 - 80, 0):c0 + 80, max(c1 - 80, 0):c1 + 80]
        reach = float(np.percentile(np.abs(patch), 99)) or 1.0
        show_frame(self.pattern_view, patch, np.arange(patch.shape[0]),
                   np.arange(patch.shape[1]), levels=(-reach, reach))
        apply_colormap(self.pattern_view, "gray", False)
        # k-space: the preview image's own spectrum, the grid regions in red
        before, _after, _ = self._images()
        image = before[model.box]
        ratio = image / np.maximum(gaussian_filter(image, 8.0), 1e-9) - 1
        n0, n1 = ratio.shape
        power = np.abs(np.fft.fft2(ratio * np.outer(np.hanning(n0), np.hanning(n1)))) ** 2
        logp = np.fft.fftshift(np.log10(power + 1e-12))
        lo, hi = np.percentile(logp, [50, 99.9])
        grey = np.clip((logp - lo) / max(hi - lo, 1e-9), 0, 1)
        rgb = np.stack([grey] * 3, axis=-1)
        mask = np.fft.fftshift(model.region)
        rgb[mask] = 0.35 * rgb[mask] + 0.65 * np.array([1.0, 0.1, 0.1])
        f0 = np.fft.fftshift(np.fft.fftfreq(n0))
        f1 = np.fft.fftshift(np.fft.fftfreq(n1))
        view = self.kspace_view
        item = view.getImageItem()
        # An RGB picture: no lookup table, no level window, and the same
        # (row = vertical axis) orientation show_frame uses for everything.
        item.setLookupTable(None)
        item.setImage((np.transpose(rgb, (1, 0, 2)) * 255).astype(np.uint8),
                      autoLevels=False, levels=(0, 255))
        item.setRect(pg.QtCore.QRectF(float(f0[0]), float(f1[0]),
                                      float(f0[-1] - f0[0]), float(f1[-1] - f1[0])))
        for axis in ("left", "bottom"):
            view.getView().getAxis(axis).enableAutoSIPrefix(False)
        view.getView().autoRange()

    def _draw_fits(self):
        self.fit_plot.clear()
        if not self.is_map or self.result is None:
            return
        fits = self.result.fits
        idx = np.arange(len(fits))
        beta = np.array([f.beta if not f.empty else np.nan for f in fits])
        shift = np.array([f.shift if not f.empty else (np.nan, np.nan) for f in fits])
        self.fit_plot.plot(idx, beta, pen=pg.mkPen("#1f5aa6", width=2),
                           name="contrast β")
        self.fit_plot.plot(idx, shift[:, 0], pen=pg.mkPen("#c23b22"),
                           name="shift, angle (px)")
        self.fit_plot.plot(idx, shift[:, 1], pen=pg.mkPen("#2a8a3e"),
                           name="shift, energy (px)")

    # -- leaving ---------------------------------------------------------------------
    def _unique(self, base):
        taken = set(self.existing_names() or ())
        if base not in taken:
            return base
        n = 2
        while f"{base} ({n})" in taken:
            n += 1
        return f"{base} ({n})"

    def parameters(self) -> dict:
        r = self.result
        m = r.model
        params = {"method": r.method, **self.settings(),
                  "grid_peaks": len(m.peaks),
                  "grid_kspace_fraction": float(m.region.mean()),
                  "grid_rms": m.rms, "box": _box_text(m.box),
                  "contrast_before": r.contrast_before,
                  "contrast_after": r.contrast_after}
        fits = [f for f in r.fits if not f.empty]
        if len(fits) > 1:
            betas = np.array([f.beta for f in fits])
            shifts = np.array([f.shift for f in fits])
            params.update(beta_median=float(np.median(betas)),
                          beta_min=float(betas.min()), beta_max=float(betas.max()),
                          shift_std_angle_px=float(shifts[:, 0].std()),
                          shift_std_energy_px=float(shifts[:, 1].std()))
        if not self.is_map and self.source.currentIndex() > 0:
            params["grid_source"] = self.candidates[self.source.currentIndex() - 1][0]
        return params

    def export(self):
        if self.result is None:
            return []
        scan = self.data.scan
        params = self.parameters()
        info = P.record_step(dict(scan.info), P.Step("degrid", params, source=self.name))
        made = []
        if self.is_map:
            made.append(MemoryData(
                self.data.kind, self.axes, self.result.values, dict(scan.labels),
                source_label=self._unique(f"{self.name} [degrid]"),
                parameters=params, prefix="degrid",
                source_path=getattr(self.data, "path", ""), source_info=info,
                source_motors=dict(getattr(scan, "fourd_info", {}) or {})))
            if self.with_grid.isChecked():
                grid_info = dict(scan.info)
                grid_info["degrid.is_grid"] = 1
                grid_info["degrid.box"] = _box_text(self.result.model.box)
                made.append(MemoryData(
                    "cut", (self.axes[1], self.axes[2]), self.result.model.full(),
                    {"x": self.labels[1], "y": self.labels[2]},
                    source_label=self._unique(f"{self.name} [grid]"),
                    parameters={"is_grid": 1, "box": _box_text(self.result.model.box),
                                "grid_rms": self.result.model.rms},
                    prefix="degrid", source_path=getattr(self.data, "path", ""),
                    source_info=grid_info))
        else:
            made.append(MemoryData(
                "cut", self.axes, self.result.values, dict(scan.labels),
                source_label=self._unique(f"{self.name} [degrid]"),
                parameters=params, prefix="degrid",
                source_path=getattr(self.data, "path", ""), source_info=info,
                source_motors=dict(getattr(scan, "fourd_info", {}) or {})))
        self.datasetsCreated.emit(made)
        self.report.setText(self.report.text() + "<br><b>Added "
                            + ", ".join(f"“{d.source_label}”" for d in made)
                            + " to the list.</b>")
        return made
