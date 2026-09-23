"""
ui/volume.py
===================
The 3-D half of the data tools: :class:`VolumeProcessDialog`, which runs
the plane-wise operations and the 3-D symmetrisation over a whole cube, and
:class:`VolumeWindow`, which draws it -- orthogonal slices, the notched
cube, an isosurface and the volume projections.

The renderer is a plain ``QPainter``. Every surface in these views is flat,
and the parallel projection of a rectangle is a parallelogram, which is
exactly an affine image transform; so a painter's algorithm over
depth-sorted faces gives a pixel-exact picture with no OpenGL driver and no
PyOpenGL dependency. It also means the view can be painted at any size,
which is what lets a notched cube be exported as a 300-dpi figure instead
of a screen-sized bitmap -- ``volume_3d_plot_w_notch.m``'s output could
only ever be as large as the window it was drawn in.
"""
from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PyQt5.QtCore import Qt, QPointF, QRectF, pyqtSignal, QTimer
from PyQt5.QtGui import (QImage, QPainter, QPixmap, QTransform, QColor, QPen,
                         QPolygonF, QFont)
from PyQt5.QtWidgets import (QSizePolicy, QWidget, QMainWindow, QVBoxLayout, QHBoxLayout,
                             QLabel, QPushButton, QComboBox, QCheckBox, QDialog,
                             QFormLayout, QDoubleSpinBox, QSpinBox, QLineEdit,
                             QDialogButtonBox, QApplication, QTabWidget,
                             QGroupBox, QMessageBox, QProgressDialog, QSlider,
                             QFileDialog, QScrollArea)

from tools import colormaps
import tools.process as P
import tools.volume as V
from ui.widgets import add_button, MemoryData, COLORMAP_NAMES
from ui.process import spin, whole, _decimate


# ==========================================================================
# Painting a set of faces
# ==========================================================================
def _colour_face(image, lut, levels, gamma: float = 1.0, shade: float = 1.0):
    """One face's samples as an ARGB image, missing points transparent.

    Transparency rather than a background colour: a face that reaches past
    the measured volume should show the face behind it, not a block of grey
    that reads as data.
    """
    values = np.asarray(image, dtype=float)
    lo, hi = float(levels[0]), float(levels[1])
    span = (hi - lo) or 1.0
    scaled = np.clip((values - lo) / span, 0.0, 1.0)
    if gamma != 1.0:
        scaled = np.power(scaled, 1.0 / max(float(gamma), 1e-6))
    # NaN casts to an arbitrary integer, so it is mapped to 0 first; those
    # pixels are then made transparent by the alpha channel below anyway.
    index = np.clip(np.nan_to_num(scaled * (len(lut) - 1), nan=0.0)
                    .astype(np.int32), 0, len(lut) - 1)
    rgb = (lut[index] * float(np.clip(shade, 0.0, 1.0))).astype(np.uint8)
    alpha = np.where(np.isfinite(values), 255, 0).astype(np.uint8)

    # QImage wants (row, column) = (v, u), so the sampled (u, v) array is
    # transposed once here rather than at every call site.
    height, width = rgb.shape[1], rgb.shape[0]
    buffer = np.empty((height, width, 4), dtype=np.uint8)
    buffer[..., 0] = rgb[..., 2].T          # BGRA, which is what Format_ARGB32 is
    buffer[..., 1] = rgb[..., 1].T
    buffer[..., 2] = rgb[..., 0].T
    buffer[..., 3] = alpha.T
    qimage = QImage(buffer.data, width, height, 4 * width, QImage.Format_ARGB32)
    clone = qimage.copy()                    # own the memory, the buffer dies here
    return clone


def _face_transform(screen, nu: int, nv: int) -> QTransform:
    """The affine map from a face's image pixels to its projected corners.

    ``screen`` holds the four projected corners in the order
    ``(origin, +u, +u+v, +v)``, so the u and v edges give the two basis
    vectors directly -- no solving, and exact for a parallel projection.
    """
    origin = np.asarray(screen[0], dtype=float)
    du = (np.asarray(screen[1], dtype=float) - origin) / max(nu, 1)
    dv = (np.asarray(screen[3], dtype=float) - origin) / max(nv, 1)
    return QTransform(du[0], du[1], dv[0], dv[1], origin[0], origin[1])


class VolumeRangeDialog(QDialog):
    """What to load into the 3-D views, before anything is computed.

    Opening a whole measured cube in the 3-D views is usually the wrong
    thing: most of it is the part of the detector nobody is looking at, and
    every one of the four views resamples the full array before it can draw
    a single frame. Choosing the momentum and energy window -- and how
    finely to sample it -- is the difference between a view that turns as
    the mouse moves and one that stalls on every drag.

    The sampling is deliberately in **points per axis**, not a "quality"
    setting: it is the number the cost actually depends on, and the dialog
    shows the resulting array size as it is changed.
    """

    def __init__(self, values, axes, labels, parent=None):
        super().__init__(parent)
        self.values = np.asarray(values, dtype=float)
        self.axes = [np.asarray(a, dtype=float) for a in axes]
        self.labels = list(labels)
        self.setWindowTitle("What to show in 3D")
        self.setModal(True)

        layout = QVBoxLayout(self)
        layout.addWidget(QLabel(
            "<b>Choose the region and the sampling.</b><br>"
            "The views resample the whole cube for every frame, so this is "
            "what decides whether they keep up with the mouse."))

        form = QFormLayout()
        self.lo_boxes, self.hi_boxes, self.point_boxes = [], [], []
        for dim, axis in enumerate(self.axes):
            row = QWidget()
            box = QHBoxLayout(row)
            box.setContentsMargins(0, 0, 0, 0)
            lo = QDoubleSpinBox()
            hi = QDoubleSpinBox()
            step = abs(float(axis[-1] - axis[0])) / max(axis.size - 1, 1)
            for widget, value in ((lo, float(axis.min())), (hi, float(axis.max()))):
                widget.setDecimals(5)
                widget.setRange(float(axis.min()), float(axis.max()))
                widget.setSingleStep(step or 0.01)
                widget.setValue(value)
                widget.setKeyboardTracking(False)
                widget.setMaximumWidth(110)
                widget.valueChanged.connect(self._update_estimate)
                box.addWidget(widget)
            points = QSpinBox()
            points.setRange(8, 2000)
            points.setValue(int(min(axis.size, 160)))
            points.setKeyboardTracking(False)
            points.setMaximumWidth(90)
            points.setToolTip("How many points to keep on this axis. "
                              "Neighbours are averaged, not thrown away.")
            points.valueChanged.connect(self._update_estimate)
            box.addWidget(QLabel("points"))
            box.addWidget(points)
            box.addStretch(1)
            form.addRow(self.labels[dim], row)
            self.lo_boxes.append(lo)
            self.hi_boxes.append(hi)
            self.point_boxes.append(points)
        layout.addLayout(form)

        presets = QWidget()
        bar = QHBoxLayout(presets)
        bar.setContentsMargins(0, 0, 0, 0)
        bar.addWidget(QLabel("Sampling"))
        for name, target in (("Fast (80)", 80), ("Balanced (160)", 160),
                             ("Fine (300)", 300), ("Everything", None)):
            button = QPushButton(name)
            button.clicked.connect(lambda _, t=target: self._preset(t))
            bar.addWidget(button)
        reset = QPushButton("Whole range")
        reset.clicked.connect(self._reset_ranges)
        bar.addWidget(reset)
        bar.addStretch(1)
        layout.addWidget(presets)

        self.estimate = QLabel("")
        self.estimate.setWordWrap(True)
        layout.addWidget(self.estimate)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)
        self._update_estimate()

    def _preset(self, target):
        for dim, box in enumerate(self.point_boxes):
            box.blockSignals(True)
            box.setValue(self.axes[dim].size if target is None
                         else int(min(target, self.axes[dim].size)))
            box.blockSignals(False)
        self._update_estimate()

    def _reset_ranges(self):
        for dim, axis in enumerate(self.axes):
            for widget, value in ((self.lo_boxes[dim], float(axis.min())),
                                  (self.hi_boxes[dim], float(axis.max()))):
                widget.blockSignals(True)
                widget.setValue(value)
                widget.blockSignals(False)
        self._update_estimate()

    def selection(self):
        ranges = [(self.lo_boxes[d].value(), self.hi_boxes[d].value())
                  for d in range(3)]
        targets = [self.point_boxes[d].value() for d in range(3)]
        return ranges, targets

    def _kept(self):
        """Points per axis after cropping and binning, without doing it."""
        ranges, targets = self.selection()
        out = []
        for axis, (lo, hi), target in zip(self.axes, ranges, targets):
            inside = int(((axis >= min(lo, hi)) & (axis <= max(lo, hi))).sum())
            factor = max(1, round(inside / max(target, 1)))
            out.append(max(inside // factor, 1))
        return out

    def _update_estimate(self):
        kept = self._kept()
        total = int(np.prod(kept))
        before = int(self.values.size)
        self.estimate.setText(
            f"{' × '.join(str(n) for n in kept)} = {total:,} points "
            f"({100.0 * total / max(before, 1):.1f}% of the {before:,} "
            f"measured), about {total * 8 / 1e6:.1f} MB.")

    def result(self):
        """The cropped, binned cube and its axes."""
        ranges, targets = self.selection()
        return V.reduce_volume(self.values, self.axes, ranges, targets)


class AxisPosition(QWidget):
    """A position along one axis: a short slider plus the value, typed.

    Two problems with a bare slider. It answers "roughly where" but never
    "exactly where", and a plane at *exactly* E_F or k = 0 is the one anyone
    actually wants; and a slider stretched across the panel emits a value
    for every pixel it is dragged through, so a cube that takes a fraction
    of a second to resample turns a single drag into a queue of rebuilds
    that the interface cannot keep up with.

    So the slider is kept short, it only reports when the drag **ends**
    (``setTracking(False)``), and the spin box beside it takes the axis's
    own units to full precision.
    """

    changed = pyqtSignal()

    def __init__(self, axis, parent=None):
        super().__init__(parent)
        self.axis = np.asarray(axis, dtype=float)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(4)

        self.slider = QSlider(Qt.Horizontal)
        self.slider.setRange(0, self.axis.size - 1)
        self.slider.setValue(self.axis.size // 2)
        # Only on release: dragging through 200 positions must not ask for
        # 200 rebuilds of the scene.
        self.slider.setTracking(False)
        self.slider.setMaximumWidth(120)
        self.slider.valueChanged.connect(self._from_slider)
        layout.addWidget(self.slider)

        step = abs(float(self.axis[-1] - self.axis[0])) / max(self.axis.size - 1, 1)
        self.box = QDoubleSpinBox()
        self.box.setDecimals(5)
        self.box.setRange(float(min(self.axis[0], self.axis[-1])),
                          float(max(self.axis[0], self.axis[-1])))
        self.box.setSingleStep(step or 0.01)
        self.box.setValue(float(self.axis[self.axis.size // 2]))
        self.box.setKeyboardTracking(False)
        self.box.setMaximumWidth(100)
        self.box.valueChanged.connect(self._from_box)
        layout.addWidget(self.box)
        layout.addStretch(1)

    def index(self) -> int:
        return int(self.slider.value())

    def position(self) -> float:
        return float(self.axis[self.index()])

    def set_position(self, value: float):
        index = int(np.argmin(np.abs(self.axis - float(value))))
        self._apply(index)

    def _apply(self, index: int):
        index = int(np.clip(index, 0, self.axis.size - 1))
        for widget, value in ((self.slider, index),
                              (self.box, float(self.axis[index]))):
            widget.blockSignals(True)
            widget.setValue(value)
            widget.blockSignals(False)
        self.changed.emit()

    def _from_slider(self, index: int):
        self._apply(index)

    def _from_box(self, value: float):
        self._apply(int(np.argmin(np.abs(self.axis - float(value)))))


class SceneView(QWidget):
    """Draws a list of :class:`V.Face` objects, and lets the mouse
    turn them."""

    rotated = pyqtSignal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.faces = []
        self.bounds = [(0, 1)] * 3
        self.camera = V.Camera()
        self.levels = (0.0, 1.0)
        self.gamma = 1.0
        self.colormap = "gray"
        self.flip = False
        self.show_edges = True
        self.show_axes = True
        self.axis_labels = ("x", "y", "z")
        self.background = "#ffffff"
        self.cull = True
        self.image = None                    # a flat projected image, if any
        self.setMinimumSize(360, 320)
        self.setMouseTracking(True)
        self._drag = None

    def set_scene(self, faces, bounds):
        self.faces = list(faces)
        self.bounds = list(bounds)
        self.update()

    # -- painting ----------------------------------------------------------
    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing, True)
        painter.setRenderHint(QPainter.SmoothPixmapTransform, True)
        painter.fillRect(self.rect(), QColor(self.background))
        self.render_into(painter, self.width(), self.height())
        painter.end()

    def render_into(self, painter, width: int, height: int):
        """Paint the scene into any painter at any size.

        Kept separate from ``paintEvent`` so the same code produces the
        widget, a PNG at 600 dpi and a PDF page -- the view on screen and the
        file are then the same picture by construction rather than by care.
        """
        if not self.faces and self.image is None:
            return
        lut = colormaps.get_lut(self.colormap, flip=self.flip)
        scale, offset = self._fit_for(width, height)
        if self.image is not None:
            self._paint_image(painter, width, height)
        visible = V.visible(self.faces, self.camera, self.bounds, cull=self.cull)
        # Which of each face's edges belong to the shape rather than to the
        # way it was cut up for drawing. Computed over the visible set, so a
        # seam whose other half is hidden is still recognised as a seam.
        real_edges = dict(zip((id(f) for f in visible),
                              V.outline_edges(visible)))
        for face, screen in V.sort_faces(visible, self.camera, self.bounds):
            pixels = self._pixels_for(screen, scale, offset, height)
            nu, nv = face.image.shape
            qimage = _colour_face(face.image, lut, self.levels, self.gamma,
                                  self._shade(face))
            painter.save()
            painter.setTransform(_face_transform(pixels, nu, nv), True)
            painter.drawImage(QRectF(0, 0, nu, nv), qimage)
            painter.restore()
            if self.show_edges:
                painter.setPen(QPen(QColor("#33000000"), 1.0))
                for i, draw in enumerate(real_edges[id(face)]):
                    if draw:
                        painter.drawLine(QPointF(*pixels[i]),
                                         QPointF(*pixels[(i + 1) % 4]))
        if self.show_axes:
            self._paint_axes(painter, width, height)

    def _fit_for(self, width, height):
        corners = np.array([[self.bounds[0][i], self.bounds[1][j], self.bounds[2][k]]
                            for i in (0, 1) for j in (0, 1) for k in (0, 1)],
                           dtype=float)
        screen, _ = self.camera.project(corners, self.bounds)
        lo, hi = screen.min(axis=0), screen.max(axis=0)
        span = np.maximum(hi - lo, 1e-9)
        margin = 0.12
        usable = np.array([width * (1 - 2 * margin), height * (1 - 2 * margin)])
        scale = float(np.min(usable / span))
        centre = (lo + hi) / 2.0
        offset = np.array([width / 2.0, height / 2.0]) - centre * scale
        return scale, offset

    def _pixels_for(self, screen, scale, offset, height):
        out = np.asarray(screen, dtype=float) * scale + offset
        out[:, 1] = height - out[:, 1]
        return out

    def _shade(self, face):
        """A little directional shading, so the three walls of a notch are
        told apart even in a single-hue colormap."""
        if not face.outward:
            return 1.0
        normal = face.normal()
        light = np.array([0.35, 0.25, 0.90])
        light /= np.linalg.norm(light)
        return 0.78 + 0.22 * float(abs(np.dot(normal, light)))

    def _paint_image(self, painter, width, height):
        values = np.asarray(self.image, dtype=float)
        lut = colormaps.get_lut(self.colormap, flip=self.flip)
        qimage = _colour_face(values, lut, self.levels, self.gamma)
        side = min(width, height) * 0.86
        left = (width - side) / 2.0
        top = (height - side) / 2.0
        painter.drawImage(QRectF(left, top, side, side), qimage)

    def _paint_axes(self, painter, width, height):
        """A labelled tripod in the corner of the view.

        Anchored to the widget rather than to a corner of the box on
        purpose: a tripod drawn at a box corner is hidden inside the solid
        for half of all viewing angles, and a notched cube is normally
        looked at from exactly those angles.
        """
        basis = self.camera.basis()
        arm = min(width, height) * 0.09
        origin = np.array([width - arm * 1.9, height - arm * 1.5])
        font = QFont()
        font.setPointSizeF(max(7.0, min(11.0, arm * 0.28)))
        painter.setFont(font)
        # Draw the axis pointing away from the viewer first, so the two in
        # front overlap it rather than the other way round.
        order = np.argsort([basis[2][d] for d in range(3)])
        for d in order:
            direction = np.zeros(3)
            direction[d] = 1.0
            tip = origin + np.array([float(direction @ basis[0]),
                                     -float(direction @ basis[1])]) * arm
            colour = ("#c0392b", "#27ae60", "#2471a3")[d]
            painter.setPen(QPen(QColor(colour), 1.8))
            painter.drawLine(QPointF(*origin), QPointF(*tip))
            painter.drawText(QPointF(*(tip + np.sign(tip - origin) *
                                       np.array([3.0, 4.0]))),
                             self.axis_labels[d])

    # -- mouse -------------------------------------------------------------
    def mousePressEvent(self, event):
        self._drag = (event.pos(), self.camera.azimuth, self.camera.elevation)

    def mouseMoveEvent(self, event):
        if self._drag is None:
            return
        start, azimuth, elevation = self._drag
        self.camera.azimuth = azimuth - 0.5 * (event.pos().x() - start.x())
        self.camera.elevation = float(np.clip(
            elevation + 0.5 * (event.pos().y() - start.y()), -89.0, 89.0))
        self.update()
        self.rotated.emit()

    def mouseReleaseEvent(self, event):
        self._drag = None


# ==========================================================================
# The 3-D viewer
# ==========================================================================
class VolumeWindow(QMainWindow):
    """Orthogonal slices, the notched cube, an isosurface and the volume
    projections, all of one cube.

    The four are tabs of one window rather than four windows because they
    answer the same question at different depths: where is the intensity in
    this volume. Switching between them keeps the camera, the levels and the
    colormap, so the comparison is honest.
    """

    closed = pyqtSignal(object)

    def __init__(self, values, axes, labels, source_name: str, parent=None,
                 colormap="gray", flip=False):
        super().__init__(parent)
        self.values = np.asarray(values, dtype=float)
        self.axes = [np.asarray(a, dtype=float) for a in axes]
        self.labels = list(labels)
        self.source_name = source_name
        self.setWindowTitle(f"3D view — {source_name}")
        self.resize(1020, 720)

        finite = self.values[np.isfinite(self.values)]
        self._auto_levels = ((float(np.percentile(finite, 1)),
                              float(np.percentile(finite, 99.5)))
                             if finite.size else (0.0, 1.0))

        # Built before any control that can fire _queue: a slider's initial
        # setValue emits valueChanged, and construction order is not a thing
        # a reader should have to hold in their head.
        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(140)
        self._timer.timeout.connect(self.rebuild)

        central = QWidget()
        layout = QHBoxLayout(central)
        self.scene = SceneView()
        self.scene.colormap, self.scene.flip = colormap, flip
        self.scene.levels = self._auto_levels
        self.scene.axis_labels = tuple(l.split(" (")[0] for l in self.labels)
        self.scene.camera.aspect = (1.0, 1.0, 0.75)
        self.scene.rotated.connect(self._camera_moved)
        layout.addWidget(self.scene, 1)

        side = QScrollArea()
        side.setWidgetResizable(True)
        side.setMinimumWidth(400)
        side.setMaximumWidth(430)
        panel = QWidget()
        side.setWidget(panel)
        form = QVBoxLayout(panel)

        self.tabs = QTabWidget()
        self.tabs.addTab(self._slice_tab(), "Slices")
        self.tabs.addTab(self._notch_tab(), "Notched cube")
        self.tabs.addTab(self._projection_tab(), "Projection")
        self.tabs.currentChanged.connect(self.rebuild)
        form.addWidget(self.tabs)
        form.addWidget(self._view_group())
        form.addStretch(1)

        self.note = QLabel("")
        self.note.setWordWrap(True)
        form.addWidget(self.note)
        layout.addWidget(side)
        self.setCentralWidget(central)
        self.rebuild()

    # -- tabs --------------------------------------------------------------
    def _axis_slider(self, dim: int):
        control = AxisPosition(self.axes[dim])
        control.changed.connect(self._queue)
        return control

    def _slice_tab(self):
        page = QWidget()
        form = QFormLayout(page)
        form.setFieldGrowthPolicy(QFormLayout.FieldsStayAtSizeHint)
        self.slice_on, self.slice_at = [], []
        for dim in range(3):
            check = QCheckBox(self.labels[dim].split(" (")[0])
            check.setChecked(True)
            check.toggled.connect(self._queue)
            control = self._axis_slider(dim)
            check.setMinimumWidth(56)
            form.addRow(check, control)
            self.slice_on.append(check)
            self.slice_at.append(control)
        hint = QLabel("<i>Type a position or drag the slider; the view "
                      "rebuilds when the drag ends. Drag the picture to turn "
                      "it.</i>")
        hint.setWordWrap(True)
        form.addRow("", hint)
        return page

    def _notch_tab(self):
        page = QWidget()
        form = QFormLayout(page)
        form.setFieldGrowthPolicy(QFormLayout.FieldsStayAtSizeHint)
        self.notch_at = []
        for dim in range(3):
            control = self._axis_slider(dim)
            control.set_position(float(self.axes[dim][int(self.axes[dim].size * 0.6)]))
            form.addRow(self.labels[dim].split(" (")[0], control)
            self.notch_at.append(control)
        self.notch_corner = QComboBox()
        for signs, text in (((1, 1, 1), "+ + +"), ((-1, 1, 1), "− + +"),
                            ((1, -1, 1), "+ − +"), ((-1, -1, 1), "− − +"),
                            ((1, 1, -1), "+ + −"), ((-1, 1, -1), "− + −"),
                            ((1, -1, -1), "+ − −"),
                            ((-1, -1, -1), "− − −")):
            self.notch_corner.addItem(text, signs)
        self.notch_corner.currentIndexChanged.connect(self._queue)
        self.notch_corner.setMaximumWidth(120)
        form.addRow("Corner removed", self.notch_corner)
        face = QPushButton("Look into it")
        face.setToolTip("Point the camera at the notched corner.")
        face.setMaximumWidth(140)
        face.clicked.connect(self._face_the_notch)
        form.addRow("", face)
        hint = QLabel("<i>Any of the eight corners, at any position inside "
                      "the box — a notch pushed all the way to an edge "
                      "simply drops the faces it squeezes to nothing.</i>")
        hint.setWordWrap(True)
        form.addRow("", hint)
        return page

    def _projection_tab(self):
        page = QWidget()
        form = QFormLayout(page)
        form.setFieldGrowthPolicy(QFormLayout.FieldsStayAtSizeHint)
        self.projection_mode = QComboBox()
        self.projection_mode.setMaximumWidth(190)
        self.projection_mode.addItems([
            "Maximum intensity", "Sum along the view", "Alpha composite",
            "Isosurface (shaded)"])
        self.projection_mode.currentIndexChanged.connect(self._queue)
        form.addRow("Mode", self.projection_mode)
        self.projection_samples = whole(140, 40, 400, 10)
        self.projection_samples.valueChanged.connect(self._queue)
        form.addRow("Ray samples", self.projection_samples)
        self.iso_level = spin(0.5, 0.0, 1.0, 0.02, 3)
        self.iso_level.setToolTip("As a fraction of the level window.")
        self.iso_level.valueChanged.connect(self._queue)
        form.addRow("Isosurface level", self.iso_level)
        self.alpha_opacity = spin(0.06, 0.005, 1.0, 0.01, 3)
        self.alpha_opacity.valueChanged.connect(self._queue)
        form.addRow("Opacity per sample", self.alpha_opacity)
        note = QLabel("<i>Maximum intensity is the only one of these in which "
                      "every pixel is a value that was actually measured; the "
                      "others are renderings, and should be labelled as "
                      "such.</i>")
        note.setWordWrap(True)
        form.addRow("", note)
        return page

    def _view_group(self):
        group = QGroupBox("View")
        form = QFormLayout(group)
        form.setFieldGrowthPolicy(QFormLayout.FieldsStayAtSizeHint)
        self.colormap = QComboBox()
        self.colormap.setMaximumWidth(160)
        self.colormap.addItems(COLORMAP_NAMES)
        self.colormap.setCurrentText(self.scene.colormap)
        self.colormap.currentIndexChanged.connect(self._restyle)
        form.addRow("Colormap", self.colormap)
        self.flip = QCheckBox("reversed")
        self.flip.setChecked(self.scene.flip)
        self.flip.toggled.connect(self._restyle)
        form.addRow("", self.flip)

        self.level_lo = spin(self._auto_levels[0], -1e12, 1e12,
                             max(abs(self._auto_levels[1]) / 50, 1e-6), 4)
        self.level_hi = spin(self._auto_levels[1], -1e12, 1e12,
                             max(abs(self._auto_levels[1]) / 50, 1e-6), 4)
        for box in (self.level_lo, self.level_hi):
            box.valueChanged.connect(self._restyle)
        self.level_lo.setMaximumWidth(140)
        self.level_lo.setMaximumWidth(140)
        self.level_lo.setMaximumWidth(140)
        self.level_lo.setMaximumWidth(140)
        self.level_lo.setMaximumWidth(140)
        self.level_lo.setMaximumWidth(140)
        form.addRow("Level low", self.level_lo)
        self.level_hi.setMaximumWidth(140)
        form.addRow("Level high", self.level_hi)
        self.gamma = spin(1.0, 0.1, 5.0, 0.1, 2)
        self.gamma.valueChanged.connect(self._restyle)
        form.addRow("Gamma", self.gamma)

        self.azimuth = spin(45.0, -360.0, 360.0, 5.0, 1, "°")
        self.elevation = spin(25.0, -89.0, 89.0, 5.0, 1, "°")
        for box in (self.azimuth, self.elevation):
            box.valueChanged.connect(self._camera_typed)
        form.addRow("Azimuth", self.azimuth)
        form.addRow("Elevation", self.elevation)
        self.z_aspect = spin(0.75, 0.05, 5.0, 0.05, 2)
        self.z_aspect.setToolTip(
            "How tall the third axis is drawn. An inverse angstrom and an "
            "electronvolt have no common length, so there is no correct "
            "value — only a legible one.")
        self.z_aspect.valueChanged.connect(self._camera_typed)
        form.addRow("Depth aspect", self.z_aspect)
        self.edges = QCheckBox("Outline the faces")
        self.edges.setChecked(True)
        self.edges.toggled.connect(self._restyle)
        form.addRow("", self.edges)

        save = QPushButton("Save the view...")
        save.clicked.connect(self.save_view)
        form.addRow("", save)
        return group

    # -- behaviour ---------------------------------------------------------
    def _queue(self, *_):
        self._timer.start()

    def _restyle(self, *_):
        self.scene.colormap = self.colormap.currentText()
        self.scene.flip = self.flip.isChecked()
        lo, hi = self.level_lo.value(), self.level_hi.value()
        self.scene.levels = (lo, hi) if hi > lo else self._auto_levels
        self.scene.gamma = self.gamma.value()
        self.scene.show_edges = self.edges.isChecked()
        self.scene.update()

    def _camera_typed(self, *_):
        self.scene.camera.azimuth = self.azimuth.value()
        self.scene.camera.elevation = self.elevation.value()
        self.scene.camera.aspect = (1.0, 1.0, self.z_aspect.value())
        if self.tabs.currentIndex() == 2:
            self._queue()               # a projection depends on the direction
        else:
            self.scene.update()

    def _camera_moved(self):
        for box, value in ((self.azimuth, self.scene.camera.azimuth),
                           (self.elevation, self.scene.camera.elevation)):
            box.blockSignals(True)
            box.setValue(value)
            box.blockSignals(False)
        if self.tabs.currentIndex() == 2:
            self._queue()

    def _face_the_notch(self):
        signs = self.notch_corner.currentData()
        self.azimuth.setValue(float(np.degrees(np.arctan2(signs[1], signs[0]))))
        self.elevation.setValue(28.0 * signs[2])

    def rebuild(self):
        tab = self.tabs.currentIndex()
        bounds = V._bounds(self.axes)
        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            if tab == 0:
                positions = [self.slice_at[d].position()
                             if self.slice_on[d].isChecked() else None
                             for d in range(3)]
                self.scene.image = None
                self.scene.cull = False
                faces = V.slice_faces(self.values, self.axes, positions,
                                      density=360)
                self.scene.set_scene(faces, bounds)
                self.note.setText(f"{len(faces)} plane(s) at " + ", ".join(
                    f"{self.labels[d].split(' (')[0]} = {p:.4g}"
                    for d, p in enumerate(positions) if p is not None))
            elif tab == 1:
                notch = [self.notch_at[d].position() for d in range(3)]
                self.scene.image = None
                self.scene.cull = True
                faces = V.notched_box(self.values, self.axes, notch,
                                      corner=self.notch_corner.currentData(),
                                      density=360)
                self.scene.set_scene(faces, bounds)
                self.note.setText(f"{len(faces)} faces built, "
                                  f"{len(V.visible(faces, self.scene.camera, bounds))}"
                                  f" facing the camera")
            else:
                self._build_projection(bounds)
        except Exception as exc:
            self.note.setText(f"<span style='color:#b00'>{exc}</span>")
        finally:
            QApplication.restoreOverrideCursor()

    def _build_projection(self, bounds):
        samples = int(self.projection_samples.value())
        cube = V.view_volume(self.values, self.axes, self.scene.camera,
                             samples=samples, bounds=bounds)
        mode = self.projection_mode.currentIndex()
        lo, hi = self.scene.levels
        if mode == 3:
            level = lo + float(self.iso_level.value()) * (hi - lo)
            depth, shade = V.isosurface_depth(cube, level)
            self.scene.image = shade
            hits = float(np.isfinite(depth).mean())
            self.note.setText(
                f"Isosurface at {level:.4g}: {100 * hits:.1f}% of the view "
                f"hits it. Shading is Lambert on the depth map's gradient.")
            self.scene.set_scene([], bounds)
            self._iso_levels = (0.0, 1.0)
            self.scene.levels = (0.0, 1.0)
            return
        kind = ["mip", "sum", "alpha"][mode]
        image = V.project_volume(cube, mode=kind, levels=(lo, hi),
                                 alpha=float(self.alpha_opacity.value()))
        self.scene.image = image
        self.scene.set_scene([], bounds)
        if kind == "alpha":
            self.scene.levels = (0.0, float(np.nanmax(image)) or 1.0)
        elif kind == "sum":
            self.scene.levels = (float(np.nanmin(image)), float(np.nanmax(image)))
        else:
            self.scene.levels = (lo, hi)
        self.note.setText(f"{self.projection_mode.currentText()} over "
                          f"{samples} samples along the view")

    def save_view(self):
        path, _ = QFileDialog.getSaveFileName(
            self, "Save the 3D view", f"{self.source_name}_3d.png",
            "PNG (*.png);;PDF (*.pdf);;SVG (*.svg)")
        if not path:
            return None
        lower = path.lower()
        width, height = 2000, int(2000 * self.scene.height() / max(self.scene.width(), 1))
        if lower.endswith(".pdf"):
            from PyQt5.QtGui import QPdfWriter
            from PyQt5.QtCore import QSizeF
            writer = QPdfWriter(path)
            writer.setResolution(600)
            writer.setPageSizeMM(QSizeF(160.0, 160.0 * height / width))
            painter = QPainter(writer)
            self.scene.render_into(painter, writer.width(), writer.height())
            painter.end()
        elif lower.endswith(".svg"):
            from PyQt5.QtSvg import QSvgGenerator
            generator = QSvgGenerator()
            generator.setFileName(path)
            generator.setSize(pg.QtCore.QSize(width, height))
            generator.setViewBox(QRectF(0, 0, width, height))
            painter = QPainter(generator)
            painter.fillRect(QRectF(0, 0, width, height),
                             QColor(self.scene.background))
            self.scene.render_into(painter, width, height)
            painter.end()
        else:
            pixmap = QPixmap(width, height)
            pixmap.fill(QColor(self.scene.background))
            painter = QPainter(pixmap)
            painter.setRenderHint(QPainter.Antialiasing, True)
            painter.setRenderHint(QPainter.SmoothPixmapTransform, True)
            self.scene.render_into(painter, width, height)
            painter.end()
            pixmap.save(path)
        self.note.setText(f"Saved {width} × {height} to {path}")
        return path

    def closeEvent(self, event):
        self.closed.emit(self)
        super().closeEvent(event)


# ==========================================================================
# 3-D processing
# ==========================================================================
class VolumeProcessDialog(QDialog):
    """The same operations as the 2-D panel, over a whole cube.

    Split from :class:`ui.process.ProcessDialog` rather than folded
    into it because the questions are different. On a cube the first
    question is always *which plane is this defined on* -- a curvature of a
    constant-energy map and a curvature of a dispersion are different
    quantities that happen to share a formula -- and answering it on every
    tab of a shared dialog would be worse than having two.
    """

    datasetsCreated = pyqtSignal(list)

    SUFFIXES = {"smooth": "_sm", "derivative": "_d2", "curvature": "_cur",
                "normalise": "_nor", "symmetry": "_symm", "despike": "_fix"}

    def __init__(self, datasets, parent=None, colormap="gray", flip=False):
        super().__init__(parent)
        self.datasets = list(datasets)
        self.colormap, self.flip = colormap, flip
        self.setWindowTitle("Data processing (3D)")
        self.setModal(False)
        self.resize(1000, 880)

        label, data = self.datasets[0]
        x, k, z, cube = data.angle_cube
        self.values = np.asarray(cube, dtype=float)
        self.axes = [np.asarray(x, float), np.asarray(k, float),
                     np.asarray(z, float)]
        scan = data.scan
        self.labels = [scan.labels.get("x", "x"), scan.labels.get("k", "y"),
                       scan.labels.get("z", "z")]

        layout = QVBoxLayout(self)
        shape = " × ".join(str(a.size) for a in self.axes)
        layout.addWidget(QLabel(
            f"<b>{label}</b><br>{shape} points &nbsp;·&nbsp; "
            + " × ".join(self.labels)))

        # The preview is one plane of the cube, at a position the user picks:
        # a 3-D operation cannot be previewed as a volume, and showing the
        # plane it is about to be applied to is more use than a thumbnail.
        self.preview_plane = QSlider(Qt.Horizontal)
        self.preview_plane.setRange(0, self.axes[2].size - 1)
        self.preview_plane.setValue(self.axes[2].size // 2)
        self.preview_plane.valueChanged.connect(self._queue)
        row = QWidget()
        box = QHBoxLayout(row)
        box.setContentsMargins(0, 0, 0, 0)
        box.addWidget(QLabel("Preview at"))
        box.addWidget(self.preview_plane, 1)
        self.plane_caption = QLabel("")
        box.addWidget(self.plane_caption)
        layout.addWidget(row)

        from ui.process import PreviewPair
        self.preview = PreviewPair()
        layout.addWidget(self.preview, 1)

        self.tabs = QTabWidget()
        self.tabs.addTab(self._plane_tab(), "Plane-wise")
        self.tabs.addTab(self._symmetry_tab(), "Symmetrise")
        self.tabs.currentChanged.connect(self.refresh)
        self.tabs.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Maximum)
        layout.addWidget(self.tabs, 0)

        self.note = QLabel("")
        self.note.setWordWrap(True)
        layout.addWidget(self.note)

        buttons = QDialogButtonBox()
        preview = add_button(buttons,"Preview", QDialogButtonBox.ActionRole)
        preview.clicked.connect(self.refresh)
        apply_button = add_button(buttons,"Apply to the whole cube",
                                         QDialogButtonBox.AcceptRole)
        apply_button.clicked.connect(self.apply)
        view = add_button(buttons,"Open the 3D view", QDialogButtonBox.ActionRole)
        view.clicked.connect(self.open_view)
        add_button(buttons,"Close", QDialogButtonBox.RejectRole).clicked.connect(
            self.reject)
        layout.addWidget(buttons)

        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(140)
        self._timer.timeout.connect(self.refresh)
        self.refresh()

    # -- tabs --------------------------------------------------------------
    def _watch(self, *widgets):
        for widget in widgets:
            for name in ("valueChanged", "currentIndexChanged", "toggled"):
                if hasattr(widget, name):
                    getattr(widget, name).connect(self._queue)
                    break
        return widgets[0] if len(widgets) == 1 else widgets

    def _queue(self, *_):
        self._timer.start()

    def _plane_tab(self):
        page = QWidget()
        form = QFormLayout(page)
        self.plane_normal = QComboBox()
        for dim in range(3):
            others = [self.labels[d].split(" (")[0] for d in range(3) if d != dim]
            self.plane_normal.addItem(
                f"{others[0]} – {others[1]} planes "
                f"(at each {self.labels[dim].split(' (')[0]})", dim)
        self.plane_normal.setCurrentIndex(2)
        form.addRow("Work on", self._watch(self.plane_normal))

        self.plane_operation = QComboBox()
        self.plane_operation.addItems([
            "Gaussian smooth", "Savitzky-Golay smooth", "Remove spikes",
            "Second derivative (first axis)", "Second derivative (second axis)",
            "Curvature 2D", "Curvature 1D (first axis)",
            "Curvature 1D (second axis)", "Normalise each line"])
        form.addRow("Operation", self._watch(self.plane_operation))

        self.plane_width_u = spin(0.0, 0.0, 1e4, 0.001, 5)
        self.plane_width_v = spin(0.0, 0.0, 1e4, 0.001, 5)
        form.addRow("Width, first axis", self._watch(self.plane_width_u))
        form.addRow("Width, second axis", self._watch(self.plane_width_v))
        self.plane_a0 = QDoubleSpinBox()
        self.plane_a0.setDecimals(6)
        self.plane_a0.setRange(1e-6, 1e6)
        self.plane_a0.setValue(1.0)
        form.addRow("a₀", self._watch(self.plane_a0))
        self.plane_norm_mode = QComboBox()
        self.plane_norm_mode.addItems(["area", "max", "mean"])
        form.addRow("Normalise to", self._watch(self.plane_norm_mode))
        note = QLabel("<i>Applied plane by plane, so a curvature of "
                      "constant-energy maps and a curvature of dispersions "
                      "are separate choices rather than the same button.</i>")
        note.setWordWrap(True)
        form.addRow("", note)
        return page

    def _symmetry_tab(self):
        page = QWidget()
        form = QFormLayout(page)
        self.sym_fold = whole(6, 1, 12)
        form.addRow("Rotation order", self._watch(self.sym_fold))
        self.sym_mirrors = QLineEdit()
        self.sym_mirrors.setPlaceholderText("e.g. 0, 30")
        self.sym_mirrors.editingFinished.connect(self._queue)
        form.addRow("Mirror lines (°)", self.sym_mirrors)
        self.sym_inversion = QCheckBox("also k → −k")
        form.addRow("", self._watch(self.sym_inversion))
        self.sym_cx = spin(0.0, -1e4, 1e4, 0.01, 5)
        self.sym_cy = spin(0.0, -1e4, 1e4, 0.01, 5)
        row = QWidget()
        box = QHBoxLayout(row)
        box.setContentsMargins(0, 0, 0, 0)
        box.addWidget(self._watch(self.sym_cx))
        box.addWidget(self._watch(self.sym_cy))
        self.sym_pick = QPushButton("Pick on the image")
        self.sym_pick.setCheckable(True)
        self.sym_pick.setToolTip(
            "Click the rotation centre on the Source preview. One centre is "
            "used for the whole cube -- a per-plane centre would shear the "
            "volume -- so pick it on the plane that shows the symmetry most "
            "clearly.")
        self.sym_pick.toggled.connect(self._arm_pick)
        box.addWidget(self.sym_pick)
        form.addRow("Centre", row)
        self.sym_output = QComboBox()
        self.sym_output.addItems(["symmetrised", "coverage map"])
        form.addRow("Result", self._watch(self.sym_output))
        note = QLabel("<i>The missing-data mask is taken per plane, not once "
                      "from the first one: for a converted k-map or a photon "
                      "energy scan the valid region differs in every "
                      "plane.</i>")
        note.setWordWrap(True)
        form.addRow("", note)
        return page

    # -- operations --------------------------------------------------------
    def _plane_axes(self, dim):
        return [self.axes[d] for d in range(3) if d != dim]

    def plane_operation_function(self):
        dim = self.plane_normal.currentData()
        axes = self._plane_axes(dim)
        choice = self.plane_operation.currentIndex()
        widths = (self.plane_width_u.value(), self.plane_width_v.value())
        a0 = self.plane_a0.value()
        mode = self.plane_norm_mode.currentText()
        params = {"planes": self.labels[dim],
                  "operation": self.plane_operation.currentText()}

        def run(plane):
            if choice == 0:
                return P.gaussian_smooth(plane, axes, widths)
            if choice == 1:
                return P.savgol_smooth(plane, axes, widths)
            if choice == 2:
                return P.despike(plane)[0]
            if choice == 3:
                return P.derivative(plane, axes, 0, 2)
            if choice == 4:
                return P.derivative(plane, axes, 1, 2)
            if choice == 5:
                return P.curvature(plane, axes, mode="2d", a0=a0)
            if choice == 6:
                return P.curvature(plane, axes, mode="x", a0=a0)
            if choice == 7:
                return P.curvature(plane, axes, mode="y", a0=a0)
            return P.normalise(plane, axes, 0, mode=mode)

        if choice in (0, 1):
            params.update({"width_u": widths[0], "width_v": widths[1]})
        if choice in (5, 6, 7):
            params["a0"] = a0
        if choice == 8:
            params["mode"] = mode
        key = {0: "smooth", 1: "smooth", 2: "despike", 3: "derivative",
               4: "derivative", 5: "curvature", 6: "curvature",
               7: "curvature", 8: "normalise"}[choice]
        return key, params, dim, run

    def _mirror_angles(self):
        text = self.sym_mirrors.text().strip()
        if not text:
            return ()
        return tuple(float(part) for part in text.replace(";", ",").split(",")
                     if part.strip())

    # -- preview -----------------------------------------------------------
    def _current_plane(self):
        """The plane the preview shows, and the two axes it lives on."""
        if self.tabs.currentIndex() == 0:
            dim = self.plane_normal.currentData()
        else:
            dim = 2
        index = int(np.clip(self.preview_plane.value(), 0, self.axes[dim].size - 1))
        plane = np.take(self.values, index, axis=dim)
        return dim, index, plane, self._plane_axes(dim)

    def refresh(self):
        dim, index, plane, axes = self._current_plane()
        # The preview slider always addresses the axis being stepped over,
        # which changes when the plane orientation does.
        if self.preview_plane.maximum() != self.axes[dim].size - 1:
            self.preview_plane.blockSignals(True)
            self.preview_plane.setRange(0, self.axes[dim].size - 1)
            self.preview_plane.setValue(self.axes[dim].size // 2)
            self.preview_plane.blockSignals(False)
            dim, index, plane, axes = self._current_plane()
        self.plane_caption.setText(
            f"{self.labels[dim].split(' (')[0]} = {self.axes[dim][index]:.4g}")

        # The two axes of the preview change with the plane orientation, so
        # the readout's names have to follow or it will label a kz axis "ky".
        others = [self.labels[d] for d in range(3) if d != dim]
        self.preview.set_labels(*others)

        small, small_axes = _decimate(plane, axes)
        self.preview.show_source(small, small_axes, self.colormap, self.flip)
        try:
            if self.tabs.currentIndex() == 0:
                key, params, _, run = self.plane_operation_function()
                result = run(small)
            else:
                key, params = "symmetry", {}
                result = P.symmetrise(
                    small, small_axes[0], small_axes[1],
                    fold=self.sym_fold.value(),
                    centre=(self.sym_cx.value(), self.sym_cy.value()),
                    mirror_angles=self._mirror_angles(),
                    inversion=self.sym_inversion.isChecked())
                result = (result.coverage if self.sym_output.currentIndex() == 1
                          else result.values)
        except Exception as exc:
            self.note.setText(f"<span style='color:#b00'>{exc}</span>")
            return
        self.preview.show_result(result, small_axes, self.colormap, self.flip)
        if key == "symmetry":
            self.preview.set_marker(self.sym_cx.value(), self.sym_cy.value())
        else:
            self.preview.clear_overlays()
        self.note.setText(f"{key}: " + ", ".join(f"{k}={v}" for k, v
                                                 in sorted(params.items())))

    def _arm_pick(self, on: bool):
        """Take the rotation centre from a click on the source preview."""
        if not on:
            self.preview.stop_picking()
            return
        self.note.setText("Click the rotation centre on the Source preview.")
        self.preview.start_picking(self._picked)

    def _picked(self, x: float, y: float):
        self.sym_pick.setChecked(False)
        for box, value in ((self.sym_cx, x), (self.sym_cy, y)):
            box.blockSignals(True)
            box.setValue(value)
            box.blockSignals(False)
        self.refresh()

    # -- apply -------------------------------------------------------------
    def apply(self, background: bool = True):
        """Run the chosen 3-D operation over every selected cube.

        Plane-wise work over a cube is the longest thing this panel does --
        hundreds of 2-D operations in a row -- so it goes to the worker
        thread, and its own per-plane progress callback becomes the job's
        progress report, which is also what makes Cancel take effect between
        planes rather than only between datasets.

        The cubes are read here first: a cube read from a file is backed by
        h5py, which two threads must not touch at once.
        """
        from ui import jobs

        try:
            prepared = []
            for name, data in self.datasets:
                x, k, z, cube = data.angle_cube
                prepared.append((name, data, np.asarray(cube, dtype=float),
                                 [np.asarray(x, float), np.asarray(k, float),
                                  np.asarray(z, float)]))
        except Exception as exc:
            QMessageBox.warning(self, "Data processing (3D)",
                                f"Could not read the cube: {exc}")
            return None

        # Every widget is read here, before the worker starts. Reading a
        # spin box from another thread mostly appears to work, which is the
        # worst kind of unsafe; and a value that changed halfway through a
        # run would silently process half the datasets differently.
        try:
            plane_wise = self.tabs.currentIndex() == 0
            if plane_wise:
                key, params, dim, run = self.plane_operation_function()
                symmetry = None
            else:
                key, dim, run = "symmetry", None, None
                symmetry = {
                    "fold": self.sym_fold.value(),
                    "centre": (self.sym_cx.value(), self.sym_cy.value()),
                    "mirrors": self._mirror_angles(),
                    "inversion": self.sym_inversion.isChecked(),
                    "coverage": self.sym_output.currentIndex() == 1,
                }
                params = {"fold": symmetry["fold"],
                          "centre": list(symmetry["centre"]),
                          "inversion": symmetry["inversion"],
                          "output": self.sym_output.currentText()}
                if symmetry["mirrors"]:
                    params["mirrors"] = list(symmetry["mirrors"])
        except Exception as exc:
            QMessageBox.warning(self, "Data processing (3D)", str(exc))
            return None

        def work(report):
            created = []
            for index, (name, data, values, axes) in enumerate(prepared):
                base = index / max(len(prepared), 1)
                span = 1.0 / max(len(prepared), 1)

                def tick(done, total, _base=base, _span=span, _name=name):
                    # The job runner turns a cancelled job into an exception
                    # raised out of report(), so returning True here is
                    # enough: nothing else has to check a flag.
                    report(_base + _span * (done / max(total, 1)), _name)
                    return True
                if plane_wise:
                    out = V.apply_plane_wise(values, run, axis=dim, progress=tick)
                else:
                    out, coverage, _centre = V.symmetrise_volume(
                        values, axes[0], axes[1], axes[2],
                        fold=symmetry["fold"], centre=symmetry["centre"],
                        mirror_angles=symmetry["mirrors"],
                        inversion=symmetry["inversion"], progress=tick)
                    if symmetry["coverage"]:
                        out = coverage
                info = P.record_step(dict(data.scan.info),
                                     P.Step(key, params, source=name))
                created.append(MemoryData(
                    data.kind, tuple(axes), out, dict(data.scan.labels),
                    source_label=f"{name}{self.SUFFIXES.get(key, '_proc')}",
                    parameters=params, prefix=f"proc.{key}",
                    source_path=getattr(data, "path", ""), source_info=info,
                    source_motors=dict(data.scan.fourd_info)))
            return created

        if background:
            started = jobs.run_job(
                self, "Processing the cube", work,
                on_done=self._finish_apply,
                on_error=lambda exc: QMessageBox.warning(
                    self, "Data processing (3D)", str(exc)),
                on_cancel=lambda: self.note.setText("Processing cancelled."))
            if started:
                self.note.setText("Processing the cube...")
                return None

        try:
            created = work(lambda *a, **k: None)
        except Exception as exc:
            QMessageBox.warning(self, "Data processing (3D)", str(exc))
            return None
        return self._finish_apply(created)

    def _finish_apply(self, created):
        """Hand the results to the list. On the GUI thread either way."""
        if created:
            self.datasetsCreated.emit(created)
            self.note.setText(f"Added {len(created)} dataset(s) to the list.")
        return created

    def open_view(self):
        label, data = self.datasets[0]
        chooser = VolumeRangeDialog(self.values, self.axes, self.labels, self)
        if chooser.exec_() != QDialog.Accepted:
            return None
        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            values, axes = chooser.result()
        except Exception as exc:
            QMessageBox.warning(self, "3D view", str(exc))
            return None
        finally:
            QApplication.restoreOverrideCursor()
        window = VolumeWindow(values, axes, self.labels, label, self,
                              self.colormap, self.flip)
        window.show()
        self._view = window
        return window
