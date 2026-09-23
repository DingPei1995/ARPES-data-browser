"""
ui/process.py
====================
The data-processing panels: one for 2-D data (:class:`ProcessDialog`), one
for 3-D (:class:`VolumeProcessDialog`), plus the displays that belong with
them -- the EDC/MDC stack plot (:class:`StackWindow`). Arithmetic between
two cuts lives in :mod:`ui.cutops`; MDC/EDC curve fitting has its own
panel, :mod:`ui.fit`.

How this differs from the MATLAB tools it replaces
--------------------------------------------------
Those tools process **in place**: ``Curvature.m`` overwrites ``<var>_cur``
in the workspace, and the parameters that produced it are gone the moment
the window closes. Here every operation produces a **new dataset in the
launcher's list**, carrying the whole chain of steps that made it
(``proc.step.1``, ``proc.step.2``, ...) in its metadata, which is saved
into the ``.nxs`` and shown by "Show information". Six months later the
figure can still say how it was made.

Nothing is computed until "Apply": the preview runs on a decimated copy so
that dragging a slider stays interactive on a 2000 x 2000 cut, and the full
resolution is only paid for once, when the result is wanted.
"""
from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PyQt5.QtCore import Qt, pyqtSignal, QTimer
from PyQt5.QtWidgets import (QSizePolicy, QWidget, QMainWindow, QVBoxLayout, QHBoxLayout,
                             QLabel, QPushButton, QComboBox, QCheckBox, QDialog,
                             QFormLayout, QDoubleSpinBox, QSpinBox, QLineEdit,
                             QDialogButtonBox, QApplication, QTabWidget,
                             QGroupBox, QMessageBox, QProgressDialog)

import tools.process as P
import tools.volume as V
from ui.widgets import (add_button, MemoryData, apply_colormap, separator,
                               section_label, COLORMAP_NAMES, strip_stock_menu,
                               plain_image_view, show_frame, fit_frame_view)

#: Above this many points the preview works on a decimated copy. Chosen so
#: that the slowest operation here (symmetrisation, which resamples the
#: whole image once per symmetry element) still redraws inside ~0.2 s.
PREVIEW_POINTS = 160_000


def spin(value, lo, hi, step, decimals=3, suffix=""):
    box = QDoubleSpinBox()
    box.setRange(lo, hi)
    box.setSingleStep(step)
    box.setDecimals(decimals)
    box.setValue(value)
    if suffix:
        box.setSuffix(suffix)
    box.setMinimumWidth(96)
    return box


def whole(value, lo, hi, step=1, suffix=""):
    box = QSpinBox()
    box.setRange(lo, hi)
    box.setSingleStep(step)
    box.setValue(value)
    if suffix:
        box.setSuffix(suffix)
    box.setMinimumWidth(80)
    return box


#: Above this luminance (0-255) a line is too pale to read against the
#: plot's light background. Several colormaps end at pure white -- ``gray``
#: and ``hot`` do, and ``jet`` passes through a luminance of 202 in the
#: middle -- so a stack coloured straight from the table loses whichever
#: curves land there. On ``gray``, which is the usual default, that is the
#: topmost curve: white on white, and invisible.
#:
#: 140 rather than something higher because that is where the contrast
#: against a near-white background reaches about 3:1, the usual floor for a
#: line you are meant to read a shape off. A cap of 170 still leaves the
#: palest curves at roughly 2:1, which looks fixed in a screenshot and is
#: still hard to follow on a projector.
MAX_LINE_LUMINANCE = 140.0


def curve_colour(lut, fraction: float):
    """A readable line colour from a colormap, as ``(r, g, b)``.

    Darkened towards black when it would otherwise be too pale, which keeps
    the hue -- and so the ordering the colour is carrying -- while making
    every curve visible. Scaling all three channels rather than clipping any
    one of them is what preserves the hue.
    """
    index = int(np.clip(fraction, 0.0, 1.0) * (len(lut) - 1))
    colour = np.asarray(lut[index], dtype=float)
    luminance = float(0.299 * colour[0] + 0.587 * colour[1] + 0.114 * colour[2])
    if luminance > MAX_LINE_LUMINANCE:
        colour = colour * (MAX_LINE_LUMINANCE / luminance)
    return tuple(int(round(c)) for c in colour)


def _decimate(values, axes, budget: int = PREVIEW_POINTS):
    """Thin an image down to roughly ``budget`` points for the preview.

    Slicing rather than binning: the preview only has to *look* like the
    result, and binning would smooth the data before the smoothing being
    previewed, which would make every slider look less effective than it is.
    """
    values = np.asarray(values, dtype=float)
    total = values.size
    if total <= budget:
        return values, [np.asarray(a, dtype=float) for a in axes]
    factor = int(np.ceil(np.sqrt(total / budget)))
    sl = tuple(slice(None, None, factor) for _ in range(values.ndim))
    return values[sl], [np.asarray(a, dtype=float)[::factor] for a in axes]


# ==========================================================================
# A small before/after preview
# ==========================================================================
class PreviewPair(QWidget):
    """The source on the left, the result on the right, same colour scale
    unless the operation changes what the numbers mean.

    Showing both is not decoration: curvature, a second derivative and a
    background subtraction all produce pictures that look plausible whatever
    the parameters, and the only way to see that a setting has eaten the
    band is to have the original beside it.

    The two views are locked together -- pan or zoom one and the other
    follows -- because a comparison between two pictures showing different
    parts of the data is worse than no comparison at all.

    **Shape.** The default is one array element drawn square, which is the
    only shape that shows the data undistorted; stretching the image to fill
    whatever box the dialog happens to give it (the obvious thing, and what
    this did at first) makes a band look steeper or shallower than it is.
    "Equal axes" is there for a k-vs-k map, where one inverse angstrom
    really should be the same length on both axes, and "Fill the box" is
    kept because it is the largest picture and sometimes that is what you
    want.
    """

    ASPECTS = ("Square pixels", "Equal axes", "Fill the box")

    def __init__(self, parent=None):
        super().__init__(parent)
        root = QVBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(3)

        row = QWidget()
        layout = QHBoxLayout(row)
        layout.setContentsMargins(0, 0, 0, 0)
        self.before = self._view("Source")
        self.after = self._view("Result")
        layout.addWidget(self.before["box"], 1)
        layout.addWidget(self.after["box"], 1)
        root.addWidget(row, 1)

        # The two viewboxes share a range, so one reset serves both and
        # dragging either moves both. Linked ViewBox to ViewBox rather than
        # PlotItem to PlotItem: the link lives on the ViewBox either way,
        # and going through it keeps this readable next to _box().
        self._box(self.after).setXLink(self._box(self.before))
        self._box(self.after).setYLink(self._box(self.before))

        controls = QWidget()
        bar = QHBoxLayout(controls)
        bar.setContentsMargins(0, 0, 0, 0)
        bar.addWidget(QLabel("Shape"))
        self.aspect = QComboBox()
        self.aspect.addItems(self.ASPECTS)
        self.aspect.setToolTip(
            "Square pixels draws the data undistorted.\n"
            "Equal axes gives one data unit the same length on both axes -- "
            "right for a k-vs-k map, meaningless for k against energy.\n"
            "Fill the box is the biggest picture, at the cost of the shape.")
        self.aspect.currentIndexChanged.connect(self._apply_aspect)
        bar.addWidget(self.aspect)
        reset = QPushButton("Reset view")
        reset.setToolTip("Undo any panning and zooming.")
        reset.clicked.connect(self.reset_view)
        bar.addWidget(reset)
        bar.addWidget(QLabel("<i>drag to pan, wheel to zoom</i>"))
        bar.addStretch(1)
        # Where the cursor is, in the data's own units. Without this the
        # previews are pictures with no numbers on them, and a centre or a
        # mirror line can only be guessed at and typed.
        self.readout = QLabel("")
        self.readout.setMinimumWidth(230)
        self.readout.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.readout.setStyleSheet("color: #555;")
        bar.addWidget(self.readout)
        root.addWidget(controls)

        self._labels = ("x", "y")
        self._pick_callback = None
        self._marker = None
        self._guide = None
        for target in (self.before, self.after):
            target["view"].scene.sigMouseMoved.connect(
                lambda position, t=target: self._on_move(t, position))
        self.before["view"].scene.sigMouseClicked.connect(self._on_click)

        # A floor, not a target: the previews take whatever the layout has
        # left over (they are the only stretching widget), and a minimum
        # large enough to be comfortable is large enough to push the tab
        # pane off the bottom of a short dialog.
        self.setMinimumHeight(240)
        self._steps = (1.0, 1.0)
        self._axes = None

    # -- coordinates, picking and guides -----------------------------------
    def set_labels(self, x_label: str, y_label: str):
        self._labels = (x_label, y_label)
        for target in (self.before, self.after):
            target["view"].view.setLabel("bottom", x_label)
            target["view"].view.setLabel("left", y_label)

    @staticmethod
    def _box(target):
        """The ViewBox of a preview.

        ``view.view`` is a PlotItem here, so that the previews carry their
        axes. A PlotItem's own ``sceneBoundingRect`` covers the axes as
        well as the data, so asking it whether it contains the pointer
        would report coordinates for a click on the axis labels.
        """
        return getattr(target["view"].view, "vb", target["view"].view)

    def _on_move(self, target, position):
        view_box = self._box(target)
        if not view_box.sceneBoundingRect().contains(position):
            return
        point = view_box.mapSceneToView(position)
        hint = "  — click to place" if self._pick_callback else ""
        self.readout.setText(
            f"{self._labels[0]} = {point.x():.5g}   "
            f"{self._labels[1]} = {point.y():.5g}{hint}")

    def start_picking(self, callback):
        """Take the next left click on the **source** preview.

        The source, not the result: the point being chosen is a property of
        the data going in, and on a symmetrised or differentiated picture
        the feature it should sit on may no longer be there to aim at.
        """
        self._pick_callback = callback
        for target in (self.before, self.after):
            target["view"].setCursor(Qt.CrossCursor)

    def stop_picking(self):
        self._pick_callback = None
        for target in (self.before, self.after):
            target["view"].unsetCursor()

    def _on_click(self, event):
        callback = self._pick_callback
        if callback is None or event.button() != Qt.LeftButton:
            return
        view_box = self._box(self.before)
        if not view_box.sceneBoundingRect().contains(event.scenePos()):
            return
        point = view_box.mapSceneToView(event.scenePos())
        event.accept()
        self.stop_picking()
        callback(float(point.x()), float(point.y()))

    def set_marker(self, x, y):
        """Show the chosen point on both previews."""
        if self._marker is None:
            self._marker = []
            for target in (self.before, self.after):
                item = pg.ScatterPlotItem(
                    size=13, symbol="+", pen=pg.mkPen("#f5a623", width=2),
                    brush=None)
                item.setZValue(50)
                target["view"].view.addItem(item)
                self._marker.append(item)
        for item in self._marker:
            item.setData([float(x)], [float(y)])

    def set_guide(self, x, y, angle):
        """Draw the mirror line itself through the chosen point.

        Typing an angle and a point and then working out from the result
        where the line went is the slow way round; drawing it is one line of
        code and removes the guess.
        """
        if self._guide is None:
            self._guide = []
            for target in (self.before, self.after):
                line = pg.InfiniteLine(angle=0.0, movable=False,
                                       pen=pg.mkPen("#f5a623", width=2,
                                                    style=Qt.DashLine))
                line.setZValue(49)
                target["view"].view.addItem(line)
                self._guide.append(line)
        for line in self._guide:
            line.setAngle(float(angle))
            line.setPos(pg.Point(float(x), float(y)))

    def clear_overlays(self):
        for group, adder in ((self._marker, None), (self._guide, None)):
            if not group:
                continue
            for target, item in zip((self.before, self.after), group):
                target["view"].view.removeItem(item)
        self._marker = None
        self._guide = None

    def _view(self, title):
        box = QGroupBox(title)
        inner = QVBoxLayout(box)
        inner.setContentsMargins(4, 4, 4, 4)
        view = plain_image_view()
        inner.addWidget(view)
        caption = QLabel("")
        caption.setAlignment(Qt.AlignCenter)
        caption.setStyleSheet("color: #666;")
        inner.addWidget(caption)
        return {"box": box, "view": view, "caption": caption}

    def _apply_aspect(self, *_):
        """Lock both viewboxes to the chosen shape.

        pyqtgraph's ratio is xScale/yScale -- screen pixels per x unit over
        screen pixels per y unit -- so one array element comes out square
        when that ratio is the y step over the x step.
        """
        dx, dy = self._steps
        mode = self.aspect.currentText()
        for target in (self.before, self.after):
            view_box = self._box(target)
            if mode == "Fill the box":
                view_box.setAspectLocked(False)
            elif mode == "Equal axes":
                view_box.setAspectLocked(True, 1.0)
            else:
                view_box.setAspectLocked(True, abs(dy / dx) if dx else 1.0)
        self.reset_view()

    def reset_view(self):
        for target in (self.before, self.after):
            if self._axes is None:
                self._box(target).autoRange(padding=0)
            else:
                fit_frame_view(target["view"], *self._axes)

    def _show(self, target, values, axes, colormap, flip, percentile=1.0):
        values = np.asarray(values, dtype=float)
        finite = values[np.isfinite(values)]
        if finite.size:
            lo, hi = np.percentile(finite, [percentile, 100.0 - percentile])
            if hi <= lo:
                lo, hi = float(finite.min()), float(finite.max()) or 1.0
        else:
            lo, hi = 0.0, 1.0
        x, y = np.asarray(axes[0], float), np.asarray(axes[1], float)
        self._axes = (x, y)
        # The view is never re-ranged here: the user's pan and zoom must
        # survive a parameter change, or every slider nudge throws away
        # where they were looking. reset_view() is how they get back.
        self._steps = show_frame(target["view"], values, x, y,
                                 levels=(float(lo), float(hi)))
        apply_colormap(target["view"], colormap, flip)
        target["caption"].setText(
            f"{values.shape[0]} × {values.shape[1]}   "
            f"[{lo:.4g} .. {hi:.4g}]")

    def show_source(self, values, axes, colormap, flip):
        previous = self._steps
        self._show(self.before, values, axes, colormap, flip)
        # Re-fit when the source is genuinely different data rather than the
        # same data with a parameter changed. The 3-D panel swaps between
        # xy, xz and yz planes, whose steps differ, and a view fitted to one
        # of them shows a sliver of another.
        if (self.before["view"].image is None or previous == (1.0, 1.0)
                or not np.allclose(previous, self._steps, rtol=1e-9)):
            self._apply_aspect()

    def show_result(self, values, axes, colormap, flip):
        self._show(self.after, values, axes, colormap, flip)

    def report(self, text: str):
        self.after["caption"].setText(text)


# ==========================================================================
# The 2-D processing panel
# ==========================================================================
class ProcessDialog(QDialog):
    """Smooth, differentiate, symmetrise, subtract, normalise, repair and
    fit a 2-D dataset.

    One panel for all of them because from the user's side they are the same
    action -- choose a dataset, set a few numbers, watch the preview, get a
    new dataset -- and because a chain like "smooth, then curvature, then
    normalise" should not mean opening three windows. Each result is named
    with the operation's suffix, so the chain reads back off the list.
    """

    #: datasets this panel produced, for the launcher to list
    datasetsCreated = pyqtSignal(list)

    SUFFIXES = {"smooth": "_sm", "derivative": "_d2", "curvature": "_cur",
                "gradient": "_grad", "background": "_bg", "fermi": "_fd",
                "despike": "_fix", "symmetry": "_symm"}

    def __init__(self, datasets, parent=None, colormap="gray", flip=False):
        super().__init__(parent)
        self.datasets = list(datasets)          # [(label, data), ...]
        self.colormap, self.flip = colormap, flip
        self.setWindowTitle("Data processing (2D)")
        self.setModal(False)
        self.resize(1000, 900)

        label, data = self.datasets[0]
        scan = data.scan
        self.values = np.asarray(scan.value, dtype=float)
        self.axes = [np.asarray(scan.x, dtype=float), np.asarray(scan.y, dtype=float)]
        self.labels = [scan.labels.get("x", "x"), scan.labels.get("y", "y")]
        self.preview_values, self.preview_axes = _decimate(self.values, self.axes)

        layout = QVBoxLayout(self)
        names = ", ".join(name for name, _ in self.datasets[:3])
        if len(self.datasets) > 3:
            names += f", ... ({len(self.datasets)} in all)"
        header = QLabel(
            f"<b>{names}</b><br>{self.values.shape[0]} × {self.values.shape[1]} "
            f"points &nbsp;·&nbsp; {self.labels[0]} × {self.labels[1]}")
        header.setWordWrap(True)
        layout.addWidget(header)

        self.preview = PreviewPair()
        self.preview.set_labels(*self.labels)
        layout.addWidget(self.preview, 1)

        self.tabs = QTabWidget()
        self.tabs.addTab(self._smooth_tab(), "Smooth")
        self.tabs.addTab(self._derivative_tab(), "Derivative")
        self.tabs.addTab(self._curvature_tab(), "Curvature")
        self.tabs.addTab(self._background_tab(), "Background")
        self.tabs.addTab(self._symmetry_tab(), "Symmetrise")
        self.tabs.currentChanged.connect(self.refresh)
        self.tabs.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Maximum)
        layout.addWidget(self.tabs, 0)

        self.note = QLabel("")
        self.note.setWordWrap(True)
        layout.addWidget(self.note)

        buttons = QDialogButtonBox()
        self.live = QCheckBox("Live preview")
        self.live.setChecked(True)
        self.live.setToolTip(
            "Recompute the right-hand picture whenever a setting changes.\n"
            "Turn off for a very large dataset, then use Preview.")
        buttons.addButton(self.live, QDialogButtonBox.ResetRole)
        preview_button = add_button(buttons,"Preview", QDialogButtonBox.ActionRole)
        preview_button.clicked.connect(self.refresh)
        apply_button = add_button(buttons,"Apply", QDialogButtonBox.AcceptRole)
        apply_button.clicked.connect(self.apply)
        apply_button.setToolTip("Add the result to the main list as a new dataset.")
        close = add_button(buttons,"Close", QDialogButtonBox.RejectRole)
        close.clicked.connect(self.reject)
        layout.addWidget(buttons)

        # Coalesce a burst of spinbox changes into one recompute, so holding
        # an arrow key does not queue fifty full previews.
        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(120)
        self._timer.timeout.connect(self.refresh)

        self.preview.show_source(self.preview_values, self.preview_axes,
                                 self.colormap, self.flip)
        self.refresh()

    # -- wiring ------------------------------------------------------------
    def _watch(self, *widgets):
        for widget in widgets:
            for signal in ("valueChanged", "currentIndexChanged", "toggled",
                           "editingFinished"):
                if hasattr(widget, signal):
                    getattr(widget, signal).connect(self._queue)
                    break
        return widgets[0] if len(widgets) == 1 else widgets

    def _queue(self, *_):
        if self.live.isChecked():
            self._timer.start()

    def _span(self, dim: int) -> float:
        axis = self.axes[dim]
        return float(abs(axis[-1] - axis[0])) or 1.0

    # -- tabs --------------------------------------------------------------
    def _smooth_tab(self):
        page = QWidget()
        form = QFormLayout(page)
        self.smooth_kind = QComboBox()
        self.smooth_kind.addItems(["Gaussian", "Savitzky-Golay", "Moving average"])
        self.smooth_kind.setToolTip(
            "Gaussian: the honest blur, width in the axis's own units.\n"
            "Savitzky-Golay: keeps peak height and width -- use before a fit.\n"
            "Moving average: a rebin that keeps the sampling.")
        form.addRow("Method", self._watch(self.smooth_kind))

        self.smooth_x = spin(self._span(0) / 100, 0.0, self._span(0), self._span(0) / 200, 5)
        self.smooth_y = spin(self._span(1) / 100, 0.0, self._span(1), self._span(1) / 200, 5)
        form.addRow(f"Width along {self.labels[0]}", self._watch(self.smooth_x))
        form.addRow(f"Width along {self.labels[1]}", self._watch(self.smooth_y))
        form.addRow("", QLabel("<i>Gaussian: sigma. Others: full window width. "
                               "In the axis's own units; 0 leaves that axis "
                               "alone.</i>"))

        self.despike_on = QCheckBox("Remove spikes first")
        self.despike_on.setToolTip(
            "Cosmic rays and dead pixels, found by a median test against the "
            "local counting noise.")
        self.despike_sigma = spin(6.0, 2.0, 20.0, 0.5, 1)
        self.despike_scale = QComboBox()
        self.despike_scale.addItems(["poisson", "global", "local"])
        row = QWidget()
        box = QHBoxLayout(row)
        box.setContentsMargins(0, 0, 0, 0)
        box.addWidget(self._watch(self.despike_on))
        box.addWidget(QLabel("above"))
        box.addWidget(self._watch(self.despike_sigma))
        box.addWidget(QLabel("σ, noise from"))
        box.addWidget(self._watch(self.despike_scale))
        box.addStretch(1)
        form.addRow("Repair", row)
        return page

    def _derivative_tab(self):
        page = QWidget()
        form = QFormLayout(page)
        self.deriv_kind = QComboBox()
        self.deriv_kind.addItems([
            f"-d²/d{self._short(1)}²  (sharpens EDC peaks)",
            f"-d²/d{self._short(0)}²  (sharpens MDC peaks)",
            "-(both), each normalised",
            f"-d/d{self._short(1)}  (first derivative)"])
        form.addRow("Quantity", self._watch(self.deriv_kind))
        self.deriv_smooth_x = spin(self._span(0) / 60, 0.0, self._span(0),
                                   self._span(0) / 200, 5)
        self.deriv_smooth_y = spin(self._span(1) / 60, 0.0, self._span(1),
                                   self._span(1) / 200, 5)
        form.addRow(f"Pre-smooth along {self.labels[0]}", self._watch(self.deriv_smooth_x))
        form.addRow(f"Pre-smooth along {self.labels[1]}", self._watch(self.deriv_smooth_y))
        self.deriv_window = spin(0.0, 0.0, max(self._span(0), self._span(1)),
                                 0.001, 5)
        form.addRow("Fit window (0 = 7 points)", self._watch(self.deriv_window))
        note = QLabel(
            "<i>A second derivative multiplies the noise by roughly the "
            "square of the sampling rate, so it is only readable after "
            "smoothing. The border is fitted, not forced to zero.</i>")
        note.setWordWrap(True)
        form.addRow("", note)
        return page

    def _curvature_tab(self):
        page = QWidget()
        form = QFormLayout(page)
        self.curv_mode = QComboBox()
        self.curv_mode.addItems(["2D (isotropic)",
                                 f"1D along {self._short(1)}",
                                 f"1D along {self._short(0)}",
                                 "Gradient enhancement I/|∇I|"])
        self.curv_mode.setToolTip(
            "2D for a Fermi surface or constant-energy map; 1D along the "
            "energy axis for a dispersion.")
        form.addRow("Mode", self._watch(self.curv_mode))

        self.curv_a0 = QDoubleSpinBox()
        self.curv_a0.setDecimals(6)
        self.curv_a0.setRange(1e-6, 1e6)
        self.curv_a0.setValue(1.0)
        self.curv_a0.setSingleStep(0.1)
        suggest = QPushButton("Suggest")
        suggest.setToolTip("Set a0 to the mean square slope of the normalised "
                           "image -- the scale at which the denominator starts "
                           "to matter.")
        suggest.clicked.connect(self._suggest_a0)
        row = QWidget()
        box = QHBoxLayout(row)
        box.setContentsMargins(0, 0, 0, 0)
        box.addWidget(self._watch(self.curv_a0))
        box.addWidget(suggest)
        box.addStretch(1)
        form.addRow("a₀ (dimensionless)", row)

        self.curv_ratio = spin(1.0, 0.01, 100.0, 0.1, 3)
        self.curv_ratio.setToolTip(
            "1 keeps the two axes balanced, which is what the [0,1] scaling "
            "is for. Raise it to weight the second axis less.")
        form.addRow("Anisotropy ratio", self._watch(self.curv_ratio))
        self.curv_smooth_x = spin(self._span(0) / 60, 0.0, self._span(0),
                                  self._span(0) / 200, 5)
        self.curv_smooth_y = spin(self._span(1) / 60, 0.0, self._span(1),
                                  self._span(1) / 200, 5)
        form.addRow(f"Pre-smooth along {self.labels[0]}", self._watch(self.curv_smooth_x))
        form.addRow(f"Pre-smooth along {self.labels[1]}", self._watch(self.curv_smooth_y))
        note = QLabel(
            "<i>Both axes and the intensity are mapped to [0,1] first, so a₀ "
            "means the same thing on the next dataset. Negative curvature is "
            "kept -- clip it in the display levels if you do not want it.</i>")
        note.setWordWrap(True)
        form.addRow("", note)
        return page

    def _background_tab(self):
        page = QWidget()
        form = QFormLayout(page)
        self.bkg_kind = QComboBox()
        self.bkg_kind.addItems(["Shirley (per line)", "Tougaard (per line)",
                                "Angle-independent (percentile)",
                                "2-D polynomial",
                                "Divide out the Fermi cut-off"])
        form.addRow("Method", self._watch(self.bkg_kind))
        self.bkg_dim = QComboBox()
        self.bkg_dim.addItems([f"along {self.labels[1]}", f"along {self.labels[0]}"])
        form.addRow("Direction", self._watch(self.bkg_dim))
        self.bkg_percentile = spin(5.0, 0.0, 50.0, 1.0, 1, " %")
        form.addRow("Percentile", self._watch(self.bkg_percentile))
        self.bkg_order = whole(2, 0, 6)
        form.addRow("Polynomial order", self._watch(self.bkg_order))
        self.bkg_clip = QCheckBox("Hold the result at zero")
        self.bkg_clip.setToolTip(
            "Physically a count cannot be negative -- but negative pixels are "
            "also the only sign that the background was over-subtracted.")
        form.addRow("", self._watch(self.bkg_clip))

        group = QGroupBox("Fermi cut-off")
        fd = QFormLayout(group)
        self.fd_ef = spin(0.0, -1e4, 1e4, 0.01, 4, " eV")
        self.fd_temp = spin(30.0, 0.1, 1000.0, 5.0, 1, " K")
        self.fd_res = spin(0.015, 0.0, 1.0, 0.005, 4, " eV")
        self.fd_cut = spin(4.0, 1.0, 12.0, 0.5, 1, " kT")
        fd.addRow("E_F", self._watch(self.fd_ef))
        fd.addRow("Temperature", self._watch(self.fd_temp))
        fd.addRow("Resolution (FWHM)", self._watch(self.fd_res))
        fd.addRow("Stop above", self._watch(self.fd_cut))
        fd.addRow("", QLabel("<i>Divides by the occupation only, and returns "
                             "NaN where the divisor is smaller than its own "
                             "noise.</i>"))
        form.addRow(group)
        return page

    def _symmetry_tab(self):
        """Mirror the data about one line.

        A line, and nothing else. Rotational symmetrisation belongs to a
        constant-energy map, which is a plane of a cube and so is handled by
        the 3-D panel; on a 2-D cut -- momentum against energy -- a rotation
        is not a symmetry of anything, so offering the order, the inversion
        and a sector wedge here only invited a meaningless choice.

        A line needs a point and a direction, so the two boxes that stay are
        the centre and the angle.
        """
        page = QWidget()
        form = QFormLayout(page)

        self.sym_angle = spin(0.0, -180.0, 180.0, 15.0, 2, "°")
        self.sym_angle.setToolTip(
            "Direction of the mirror line, measured from the x axis. "
            "0° mirrors top to bottom about a horizontal line; 90° mirrors "
            "left to right.")
        form.addRow("Mirror line", self._watch(self.sym_angle))

        self.sym_cx = spin(0.0, -1e4, 1e4, 0.01, 5)
        self.sym_cy = spin(0.0, -1e4, 1e4, 0.01, 5)
        row = QWidget()
        box = QHBoxLayout(row)
        box.setContentsMargins(0, 0, 0, 0)
        box.addWidget(QLabel(self._short(0)))
        box.addWidget(self._watch(self.sym_cx))
        box.addWidget(QLabel(self._short(1)))
        box.addWidget(self._watch(self.sym_cy))
        self.sym_pick = QPushButton("Pick on the image")
        self.sym_pick.setCheckable(True)
        self.sym_pick.setToolTip(
            "Click a point on the Source preview to put the line through it. "
            "The dashed line shows where it lies.")
        self.sym_pick.toggled.connect(self._arm_pick)
        box.addWidget(self.sym_pick)
        box.addStretch(1)
        form.addRow("The line passes through", row)

        note = QLabel("<i>Mirroring hides a real broken symmetry as readily "
                      "as it fills a gap. The note below the tabs reports how "
                      "far the two halves actually differed — read it "
                      "before believing the result.</i>")
        note.setWordWrap(True)
        form.addRow("", note)
        return page

    def _short(self, dim: int) -> str:
        return self.labels[dim].split(" (")[0]

    # -- the operation for the current tab ---------------------------------
    def operation(self):
        """``(key, parameters, function)`` for the tab that is showing.

        One definition serves both the preview and Apply, so what is shown
        is by construction what is produced.

        Dispatched on the tab's **name**, not its index. An index-based
        table is silently wrong the moment a tab is added or removed -- the
        preview then shows one operation while Apply performs another, and
        nothing about that failure looks like a failure.
        """
        by_name = {
            "Smooth": self._smooth_operation,
            "Derivative": self._derivative_operation,
            "Curvature": self._curvature_operation,
            "Background": self._background_operation,
            "Symmetrise": self._symmetry_operation,
        }
        name = self.tabs.tabText(self.tabs.currentIndex())
        try:
            return by_name[name]()
        except KeyError:
            raise ValueError(f"no operation is defined for the {name!r} tab")

    def _smooth_operation(self):
        kind = ["gaussian", "savgol", "box"][self.smooth_kind.currentIndex()]
        widths = (self.smooth_x.value(), self.smooth_y.value())
        spikes = self.despike_on.isChecked()
        threshold = self.despike_sigma.value()
        scale = self.despike_scale.currentText()
        params = {"method": kind, "width_x": widths[0], "width_y": widths[1]}
        if spikes:
            params.update({"despike": threshold, "despike_scale": scale})

        def run(values, axes):
            out = values
            if spikes:
                out, _ = P.despike(out, threshold=threshold, scale=scale)
            if any(w > 0 for w in widths):
                out = P.SMOOTHERS[kind](out, axes, widths)
            return out, axes
        return "smooth", params, run

    def _derivative_operation(self):
        choice = self.deriv_kind.currentIndex()
        widths = (self.deriv_smooth_x.value(), self.deriv_smooth_y.value())
        window = self.deriv_window.value() or None
        params = {"quantity": ["d2_y", "d2_x", "both", "d1_y"][choice],
                  "smooth_x": widths[0], "smooth_y": widths[1]}
        if window:
            params["window"] = window

        def run(values, axes):
            out = values
            if any(w > 0 for w in widths):
                out = P.gaussian_smooth(out, axes, widths)
            if choice == 0:
                out = P.derivative(out, axes, 1, 2, window=window)
            elif choice == 1:
                out = P.derivative(out, axes, 0, 2, window=window)
            elif choice == 2:
                out = P.laplacian(out, axes, window=window)
            else:
                out = P.derivative(out, axes, 1, 1, window=window)
            return out, axes
        return "derivative", params, run

    def _curvature_operation(self):
        choice = self.curv_mode.currentIndex()
        a0 = self.curv_a0.value()
        ratio = self.curv_ratio.value()
        widths = (self.curv_smooth_x.value(), self.curv_smooth_y.value())
        mode = ["2d", "y", "x"][choice] if choice < 3 else None
        params = ({"mode": mode, "a0": a0, "ratio": ratio,
                   "smooth_x": widths[0], "smooth_y": widths[1]} if mode
                  else {"smooth_x": widths[0], "smooth_y": widths[1]})

        def run(values, axes):
            out = values
            if any(w > 0 for w in widths):
                out = P.gaussian_smooth(out, axes, widths)
            if mode is None:
                return P.gradient_enhance(out, axes), axes
            return P.curvature(out, axes, mode=mode, a0=a0, ratio=ratio), axes
        return ("curvature" if mode else "gradient"), params, run

    def _background_operation(self):
        choice = self.bkg_kind.currentIndex()
        dim = 1 if self.bkg_dim.currentIndex() == 0 else 0
        clip = self.bkg_clip.isChecked()
        percentile = self.bkg_percentile.value()
        order = self.bkg_order.value()
        if choice == 4:
            ef, temp = self.fd_ef.value(), self.fd_temp.value()
            res, cut = self.fd_res.value(), self.fd_cut.value()
            params = {"ef": ef, "temperature": temp, "resolution": res,
                      "cutoff_kt": cut, "axis": dim}

            def run(values, axes):
                return P.divide_fermi_edge(values, axes, ef=ef, temperature=temp,
                                           resolution=res, dim=dim,
                                           cutoff_kt=cut), axes
            return "fermi", params, run

        method = ["shirley", "tougaard", "percentile", "polynomial"][choice]
        params = {"method": method, "axis": dim, "clip": clip}
        if method == "percentile":
            params["percentile"] = percentile
        if method == "polynomial":
            params["order"] = order

        def run(values, axes):
            out, _ = P.subtract_background(values, axes, method, dim=dim,
                                           clip=clip, percentile=percentile,
                                           order=order)
            return out, axes
        return "background", params, run

    def _symmetry_operation(self):
        angle = self.sym_angle.value()
        centre = (self.sym_cx.value(), self.sym_cy.value())
        params = {"mirror": angle, "centre": list(centre)}

        def run(values, axes):
            result = P.symmetrise(values, axes[0], axes[1], fold=1,
                                  centre=centre, mirror_angles=(angle,))
            # How far the two halves disagreed is the one number that says
            # whether the result should be trusted, so it is reported rather
            # than left for the user to work out from the picture.
            self._symmetry_residual = result.residual
            return result.values, axes
        return "symmetry", params, run

    @staticmethod
    def _number(edit):
        text = edit.text().strip()
        if not text:
            return None
        try:
            return float(text)
        except ValueError:
            raise ValueError(f"{text!r} is not a number")

    # -- preview -----------------------------------------------------------
    def refresh(self):
        try:
            key, params, run = self.operation()
            values, axes = run(self.preview_values, self.preview_axes)
        except Exception as exc:
            self.note.setText(f"<span style='color:#b00'>{exc}</span>")
            return
        self.preview.show_result(values, axes, self.colormap, self.flip)
        if key == "symmetry":
            self._show_mirror_guide()
        else:
            self.preview.clear_overlays()

        pieces = [f"{key}: " + ", ".join(f"{k}={v}" for k, v
                                         in sorted(params.items()))]
        finite = np.isfinite(values)
        missing = 100.0 * (1.0 - finite.mean()) if finite.size else 0.0
        if missing > 0.05:
            pieces.append(f"{missing:.1f}% missing")
        residual = getattr(self, "_symmetry_residual", None)
        if key == "symmetry" and residual is not None and np.isfinite(residual):
            pieces.append(f"the two halves differed by {100 * residual:.1f}%")
        if self.preview_values.size < self.values.size:
            pieces.append("<i>preview on a decimated copy</i>")
        self.note.setText("  ·  ".join(pieces))

    def _suggest_a0(self):
        choice = self.curv_mode.currentIndex()
        mode = ["2d", "y", "x"][choice] if choice < 3 else "2d"
        widths = (self.curv_smooth_x.value(), self.curv_smooth_y.value())
        values = self.preview_values
        if any(w > 0 for w in widths):
            values = P.gaussian_smooth(values, self.preview_axes, widths)
        self.curv_a0.setValue(P.suggest_a0(values, self.preview_axes, mode=mode))
        self.refresh()

    def _arm_pick(self, on: bool):
        """Take the mirror line's point from a click on the source preview."""
        if not on:
            self.preview.stop_picking()
            return
        self.note.setText("Click a point on the Source preview to put the "
                          "mirror line through it.")
        self.preview.start_picking(self._picked)

    def _picked(self, x: float, y: float):
        self.sym_pick.setChecked(False)
        for box, value in ((self.sym_cx, x), (self.sym_cy, y)):
            box.blockSignals(True)
            box.setValue(value)
            box.blockSignals(False)
        self.refresh()

    def _show_mirror_guide(self):
        """Keep the marker and the dashed line on the chosen point."""
        x, y = self.sym_cx.value(), self.sym_cy.value()
        self.preview.set_marker(x, y)
        self.preview.set_guide(x, y, self.sym_angle.value())

    # -- apply -------------------------------------------------------------
    def apply(self, background: bool = True):
        """Run the chosen operation over every selected dataset.

        The arrays are read here, on the GUI thread, and only the arithmetic
        goes to the worker: a dataset read from a file is backed by h5py,
        which is not safe to touch from two threads, and the GUI keeps
        slicing the same file to repaint.

        ``background=False`` runs it inline, which is how a script (or a
        test) gets the result as the return value rather than through
        ``datasetsCreated``.
        """
        from ui import jobs

        try:
            key, params, run = self.operation()
        except Exception as exc:
            QMessageBox.warning(self, "Data processing", str(exc))
            return None

        try:
            prepared = [
                (name, data,
                 np.asarray(data.scan.value, dtype=float),
                 [np.asarray(data.scan.x, dtype=float),
                  np.asarray(data.scan.y, dtype=float)])
                for name, data in self.datasets]
        except Exception as exc:
            QMessageBox.warning(self, "Data processing", f"Could not read: {exc}")
            return None

        def work(report):
            made = []
            for index, (name, data, values, axes) in enumerate(prepared):
                report(index / max(len(prepared), 1), name)
                out, new_axes = run(values, axes)
                made.append(self._wrap(data, name, key, params, out, new_axes))
            return made

        if background:
            started = jobs.run_job(
                self, "Processing", work,
                on_done=lambda made: self._finish_apply(made, key),
                on_error=lambda exc: QMessageBox.warning(
                    self, "Data processing", str(exc)),
                on_cancel=lambda: self.note.setText("Processing cancelled."))
            if started:
                self.note.setText("Processing...")
                return None

        try:
            created = work(lambda *a, **k: None)
        except Exception as exc:
            QMessageBox.warning(self, "Data processing", str(exc))
            return None
        return self._finish_apply(created, key)

    def _finish_apply(self, created, key):
        """Hand the results to the list. On the GUI thread either way."""
        if created:
            self.datasetsCreated.emit(created)
            self.note.setText(
                f"Added {len(created)} dataset(s) ending in "
                f"“{self.SUFFIXES.get(key, '_proc')}” to the list.")
        return created

    def _wrap(self, data, name, key, params, values, axes):
        suffix = self.SUFFIXES.get(key, "_proc")
        info = P.record_step(dict(data.scan.info),
                             P.Step(key, params, source=name))
        return MemoryData("cut", tuple(axes), values, dict(data.scan.labels),
                          source_label=f"{name}{suffix}", parameters=params,
                          prefix=f"proc.{key}",
                          source_path=getattr(data, "path", ""),
                          source_info=info,
                          source_motors=dict(data.scan.fourd_info))


# ==========================================================================
# The EDC/MDC stack plot  (EDC_MDC_plot.m's display half, stack_plot.m)
# ==========================================================================
class StackWindow(QMainWindow):
    """A waterfall of curves taken across a 2-D dataset.

    Three things the MATLAB pair did not do. The offset can follow the
    curve's **real position** rather than its index, so twenty EDCs taken at
    unevenly spaced angles stack in proportion to the angle rather than
    evenly; the colour can run along a colormap keyed to that position, which
    is the only readable option past about eight curves where
    ``stack_plot.m``'s eight-colour cycle starts repeating; and each curve
    can be normalised, so a stack does not simply show the intensity falling
    off with binding energy.
    """

    datasetsCreated = pyqtSignal(list)

    def __init__(self, values, axes, labels, source_name: str, parent=None,
                 colormap="jet"):
        super().__init__(parent)
        self.values = np.asarray(values, dtype=float)
        self.axes = [np.asarray(a, dtype=float) for a in axes]
        self.labels = list(labels)
        self.source_name = source_name
        self.setWindowTitle(f"Stack plot — {source_name}")
        self.resize(900, 700)

        central = QWidget()
        layout = QHBoxLayout(central)
        self.plot = pg.PlotWidget()
        strip_stock_menu(self.plot.getPlotItem())
        layout.addWidget(self.plot, 1)

        side = QWidget()
        side.setMinimumWidth(380)
        side.setMaximumWidth(400)
        form = QFormLayout(side)

        self.direction = QComboBox()
        self.direction.addItems([f"EDCs: along {self.labels[1]}",
                                 f"MDCs: along {self.labels[0]}"])
        form.addRow("Curves", self.direction)

        self.count = whole(12, 2, 400)
        form.addRow("How many", self.count)
        self.combine = whole(1, 1, 99)
        self.combine.setToolTip("Sum this many neighbouring lines per curve.")
        form.addRow("Combine", self.combine)

        self.from_edit, self.to_edit = QLineEdit(), QLineEdit()
        for edit, which in ((self.from_edit, "from"), (self.to_edit, "to")):
            edit.setPlaceholderText(which)
            edit.setMaximumWidth(90)
        row = QWidget()
        box = QHBoxLayout(row)
        box.setContentsMargins(0, 0, 0, 0)
        box.addWidget(self.from_edit)
        box.addWidget(self.to_edit)
        form.addRow("Curves from", row)

        # Where each curve starts and stops along its own axis. The pair
        # above chooses *which* curves; this pair trims *each* of them, and
        # a selection box on the viewer sets both at once.
        self.limit_from, self.limit_to = QLineEdit(), QLineEdit()
        for edit, which in ((self.limit_from, "from"), (self.limit_to, "to")):
            edit.setPlaceholderText(which)
            edit.setMaximumWidth(90)
        limit_row = QWidget()
        limit_box = QHBoxLayout(limit_row)
        limit_box.setContentsMargins(0, 0, 0, 0)
        limit_box.addWidget(self.limit_from)
        limit_box.addWidget(self.limit_to)
        form.addRow("Each curve over", limit_row)

        self.normalise = QComboBox()
        self.normalise.addItems(["none", "area", "max"])
        form.addRow("Normalise each", self.normalise)

        self.offset_mode = QComboBox()
        self.offset_mode.addItems(["by position (true spacing)", "evenly by index"])
        self.offset_mode.setToolTip(
            "By position keeps the vertical spacing proportional to where "
            "each curve was taken, so an unevenly sampled series does not "
            "pretend to be even.")
        form.addRow("Offset", self.offset_mode)
        self.offset = spin(1.0, 0.0, 1e6, 0.1, 4)
        form.addRow("Offset scale", self.offset)
        self.x_offset = spin(0.0, -1e6, 1e6, 0.001, 5)
        form.addRow("Shear along x", self.x_offset)

        self.colour_mode = QComboBox()
        self.colour_mode.addItems(["colormap by position", "single colour"])
        form.addRow("Colour", self.colour_mode)
        self.colormap = QComboBox()
        self.colormap.addItems(COLORMAP_NAMES)
        self.colormap.setCurrentText(colormap if colormap in COLORMAP_NAMES else "jet")
        form.addRow("Colormap", self.colormap)
        self.width = spin(1.2, 0.2, 8.0, 0.2, 1)
        form.addRow("Line width", self.width)
        self.annotate = QCheckBox("Label each curve")
        self.annotate.setChecked(True)
        form.addRow("", self.annotate)
        self.fill = QCheckBox("Fill under the curves")
        self.fill.setToolTip("Hides the curves behind, which is what makes a "
                             "dense stack readable.")
        form.addRow("", self.fill)

        # The direction gets its own handler, which has to run *before* the
        # redraw: the two range pairs mean different axes in the two
        # directions, so they have to be swapped first.
        self.direction.currentIndexChanged.connect(self._swap_ranges)

        for widget in (self.direction, self.count, self.combine, self.normalise,
                       self.offset_mode, self.offset, self.x_offset,
                       self.colour_mode, self.colormap, self.width,
                       self.annotate, self.fill):
            for name in ("valueChanged", "currentIndexChanged", "toggled"):
                if hasattr(widget, name):
                    getattr(widget, name).connect(self.redraw)
                    break
        for edit in (self.from_edit, self.to_edit,
                     self.limit_from, self.limit_to):
            edit.editingFinished.connect(self.redraw)

        buttons = QWidget()
        bbox = QVBoxLayout(buttons)
        bbox.setContentsMargins(0, 0, 0, 0)
        to_figure = QPushButton("Send to figure composer")
        to_figure.setToolTip("One panel with these curves, ready for export.")
        to_figure.clicked.connect(self.to_figure)
        bbox.addWidget(to_figure)
        form.addRow("", buttons)

        self.note = QLabel("")
        self.note.setWordWrap(True)
        form.addRow("", self.note)
        layout.addWidget(side)
        self.setCentralWidget(central)
        self.redraw()

    def _swap_ranges(self):
        """Exchange the two range pairs when the direction changes.

        "Curves from" is read against the axis the curves are spaced along
        and "Each curve over" against the axis they run along -- and those
        two axes trade places when EDC becomes MDC. Left alone, the boxes
        keep the *other* axis's numbers: switching a stack of EDCs taken
        over k = 0.1 .. 0.45 to MDCs asked for energies of 0.1 and 0.45 eV,
        both past the top of the energy axis, so both clamped to the same
        index and the stack collapsed to a single curve.
        """
        pairs = ((self.from_edit, self.limit_from), (self.to_edit, self.limit_to))
        for a, b in pairs:
            a_text, b_text = a.text(), b.text()
            a.setText(b_text)
            b.setText(a_text)

    def curves(self):
        """``(x, list of (position, y))`` for the current settings."""
        edc = self.direction.currentIndex() == 0
        along = self.axes[1] if edc else self.axes[0]
        across = self.axes[0] if edc else self.axes[1]
        lo, hi = 0, across.size - 1
        for edit, default, pick in ((self.from_edit, 0, min),
                                    (self.to_edit, across.size - 1, max)):
            text = edit.text().strip()
            if not text:
                continue
            try:
                value = float(text)
            except ValueError:
                continue
            index = int(np.argmin(np.abs(across - value)))
            if edit is self.from_edit:
                lo = index
            else:
                hi = index
        if lo > hi:
            lo, hi = hi, lo
        wanted = np.unique(np.linspace(lo, hi, int(self.count.value())).round()
                           .astype(int))
        half = max(0, (int(self.combine.value()) - 1) // 2)

        # Trim every curve to the same window along its own axis. Done
        # before the normalisation on purpose: normalising to the area of
        # the part being shown is what makes a set of curves comparable, and
        # normalising to the area of a tail that is then cut off is not.
        keep = np.ones(along.size, dtype=bool)
        for edit in (self.limit_from, self.limit_to):
            text = edit.text().strip()
            if not text:
                continue
            try:
                value = float(text)
            except ValueError:
                continue
            keep &= (along >= value) if edit is self.limit_from else (along <= value)
        if not keep.any():
            keep = np.ones(along.size, dtype=bool)
        along = along[keep]

        out = []
        for i in wanted:
            first, last = max(0, i - half), min(across.size, i + half + 1)
            block = (self.values[first:last, :] if edc
                     else self.values[:, first:last])
            line = np.nansum(block, axis=0 if edc else 1)[keep]
            mode = self.normalise.currentText()
            if mode == "area":
                total = np.nansum(line) * P.axis_step(along)
                line = line / total if abs(total) > 1e-30 else line
            elif mode == "max":
                top = np.nanmax(line)
                line = line / top if abs(top) > 1e-30 else line
            out.append((float(across[i]), line))
        return along, out

    def redraw(self):
        from tools import colormaps
        self.plot.clear()
        try:
            along, curves = self.curves()
        except Exception as exc:
            self.note.setText(f"<span style='color:#b00'>{exc}</span>")
            return
        if not curves:
            return
        edc = self.direction.currentIndex() == 0
        self.plot.setLabel("bottom", self.labels[1] if edc else self.labels[0])
        self.plot.setLabel("left", "Intensity (offset)")

        positions = np.array([p for p, _ in curves], dtype=float)
        span = float(positions[-1] - positions[0]) or 1.0
        scale = float(self.offset.value())
        # A "unit" offset means one typical curve height, so the scale box
        # behaves the same whether the data is in counts or normalised.
        heights = [np.nanmax(y) - np.nanmin(y) for _, y in curves]
        unit = float(np.nanmedian(heights)) or 1.0

        lut = colormaps.get_lut(self.colormap.currentText())
        by_position = self.offset_mode.currentIndex() == 0
        pen_width = float(self.width.value())

        for i, (position, y) in enumerate(curves):
            fraction = ((position - positions[0]) / span if by_position
                        else i / max(len(curves) - 1, 1))
            shift = fraction * scale * unit * max(len(curves) - 1, 1)
            x = along + fraction * float(self.x_offset.value()) * \
                max(len(curves) - 1, 1)
            if self.colour_mode.currentIndex() == 0:
                colour = curve_colour(lut, fraction)
                pen = pg.mkPen(colour, width=pen_width)
                brush = pg.mkBrush(*colour, 90)
            else:
                pen = pg.mkPen("#bd4921", width=pen_width)
                brush = pg.mkBrush(189, 73, 33, 90)
            item = self.plot.plot(x, y + shift, pen=pen)
            if self.fill.isChecked():
                item.setFillLevel(shift)
                item.setBrush(brush)
            if self.annotate.isChecked():
                text = pg.TextItem(f"{position:.4g}", anchor=(0, 0.5),
                                   color=pen.color())
                text.setPos(float(x[-1]), float(y[-1] + shift))
                self.plot.addItem(text)
        self.plot.getPlotItem().enableAutoRange()
        asked = int(self.count.value())
        if len(curves) < asked:
            # A range covering fewer points than the number of curves asked
            # for silently collapses to however many distinct lines exist.
            # Saying so is the difference between "the range is narrow" and
            # "the tool is broken".
            self._shortfall = (
                f"  ·  only {len(curves)} of the {asked} asked for: the "
                f"range covers {len(curves)} point(s)")
        else:
            self._shortfall = ""
        if self.annotate.isChecked():
            # The labels sit past the end of each curve, so the view has to
            # leave room for them or they are cut off by the frame -- the
            # autoRange only knows about the curves themselves.
            self.plot.getPlotItem().getViewBox().updateAutoRange()
            (x0, x1), _ = self.plot.getPlotItem().getViewBox().viewRange()
            self.plot.setXRange(x0, x1 + 0.16 * (x1 - x0), padding=0)
        self.note.setText(f"{len(curves)} curves, "
                          f"{positions[0]:.4g} to {positions[-1]:.4g}"
                          + getattr(self, "_shortfall", ""))

    def to_figure(self):
        from ui.figure import FigureWindow
        import tools.figure as F
        from tools import colormaps

        along, curves = self.curves()
        if not curves:
            return None
        edc = self.direction.currentIndex() == 0
        positions = np.array([p for p, _ in curves], dtype=float)
        span = float(positions[-1] - positions[0]) or 1.0
        heights = [np.nanmax(y) - np.nanmin(y) for _, y in curves]
        unit = float(np.nanmedian(heights)) or 1.0
        scale = float(self.offset.value())
        lut = colormaps.get_lut(self.colormap.currentText())

        lines = []
        for i, (position, y) in enumerate(curves):
            fraction = ((position - positions[0]) / span
                        if self.offset_mode.currentIndex() == 0
                        else i / max(len(curves) - 1, 1))
            shift = fraction * scale * unit * max(len(curves) - 1, 1)
            colour = curve_colour(lut, fraction)
            lines.append(F.Overlay(
                xs=np.asarray(along, dtype=float).tolist(),
                ys=np.asarray(y + shift, dtype=float).tolist(),
                color="#%02x%02x%02x" % colour,
                width_pt=float(self.width.value()), name=f"{position:.4g}"))

        # A stack has no image, so the panel carries an empty one of the right
        # extent and the curves go in as overlays -- the polyline primitive
        # the figure model already has. The painter then handles ticks,
        # labels, panel letters and export exactly as it does for a picture,
        # which is the whole reason for routing a curve plot through it.
        ys = np.concatenate([np.asarray(line.ys, dtype=float) for line in lines])
        blank = np.full((2, 2), np.nan)
        panel = F.Panel(data=F.PanelData(
            array=blank,
            x=np.array([float(along[0]), float(along[-1])]),
            y=np.array([float(np.nanmin(ys)), float(np.nanmax(ys))]),
            x_label=self.labels[1] if edc else self.labels[0],
            y_label="Intensity (offset)", name=self.source_name, kind="stack"))
        panel.overlays = lines
        window = FigureWindow([panel], parent=self)
        window.show()
        self._figure = window
        return window

