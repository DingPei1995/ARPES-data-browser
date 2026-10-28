"""
ui/panels/volume.py
===================
The 3-D half of the data tools: :class:`VolumeProcessDialog`, which runs
the plane-wise operations and the 3-D symmetrisation over a whole cube, and
:class:`VolumeWindow`, which draws it -- one window each for the
orthogonal slices, the notched cube, and the volume projections (with the
isosurface).

The renderer is a plain ``QPainter``. Every surface in these views is flat,
and the parallel projection of a rectangle is a parallelogram, which is
exactly an affine image transform; so a painter's algorithm over
depth-sorted faces gives a pixel-exact picture with no OpenGL driver and no
PyOpenGL dependency. It also means the view can be painted at any size,
which is what lets a notched cube be exported as a 300-dpi figure instead
of a screen-sized bitmap -- ``volume_3d_plot_w_notch.m``'s output could
only ever be as large as the window it was drawn in. The figure composer
draws a 3-D panel with the same function (:func:`paint_scene`), handing it
the figure's font and line widths, so the axes of a 3-D panel are labelled
like those of every other panel.
"""
from __future__ import annotations

import copy
from dataclasses import dataclass, field, replace

import numpy as np
import pyqtgraph as pg
from PyQt5.QtCore import Qt, QPointF, QRectF, pyqtSignal, QTimer
from PyQt5.QtGui import (QImage, QPainter, QPixmap, QTransform, QColor, QPen,
                         QPolygonF, QFont, QFontMetricsF)
from PyQt5.QtWidgets import (QSizePolicy, QWidget, QMainWindow, QVBoxLayout, QHBoxLayout,
                             QLabel, QPushButton, QComboBox, QCheckBox, QDialog,
                             QFormLayout, QDoubleSpinBox, QSpinBox, QLineEdit,
                             QDialogButtonBox, QApplication, QTabWidget,
                             QGroupBox, QMessageBox, QProgressDialog, QSlider,
                             QFileDialog, QScrollArea, QInputDialog)

from ui.common.screen import fits_screen
from tools.output import colormaps
import tools.processing.process as P
import tools.processing.volume as V
from ui.viewers.widgets import add_button, ButtonRow, MemoryData, COLORMAP_NAMES
from ui.panels.process import spin, whole, _decimate


# ==========================================================================
# Painting a set of faces
# ==========================================================================
def _colour_face(image, lut, levels, gamma: float = 1.0, shade: float = 1.0,
                 opacity: float = 1.0, fade: float = 0.0):
    """One face's samples as an ARGB image, missing points transparent.

    Transparency rather than a background colour: a face that reaches past
    the measured volume should show the face behind it, not a block of grey
    that reads as data.

    ``opacity`` (0..1) is the whole face's; ``fade`` (0..1 of the level
    window) makes the weakest intensities transparent too, rising from
    nothing at the low level to ``opacity`` at that fraction of the window
    -- so a slice can be looked *through* where it holds only background.
    """
    values = np.asarray(image, dtype=float)
    lo, hi = float(levels[0]), float(levels[1])
    span = (hi - lo) or 1.0
    scaled = np.clip((values - lo) / span, 0.0, 1.0)
    # NaN casts to an arbitrary integer, so it is mapped to 0 first; those
    # pixels are then made transparent by the alpha channel below anyway.
    scaled = np.nan_to_num(scaled, nan=0.0)
    coloured = scaled
    if gamma != 1.0:
        coloured = np.power(scaled, 1.0 / max(float(gamma), 1e-6))
    index = np.clip((coloured * (len(lut) - 1)).astype(np.int32), 0, len(lut) - 1)
    rgb = (lut[index] * float(np.clip(shade, 0.0, 1.0))).astype(np.uint8)
    alpha = np.where(np.isfinite(values), float(np.clip(opacity, 0.0, 1.0)), 0.0)
    if fade > 0:
        alpha = alpha * np.clip(scaled / float(fade), 0.0, 1.0)
    alpha = np.rint(alpha * 255.0).astype(np.uint8)

    # QImage wants (row, column) = (v, u), so the sampled (u, v) array is
    # transposed once here rather than at every call site.
    height, width = rgb.shape[1], rgb.shape[0]
    buffer = np.empty((height, width, 4), dtype=np.uint8)
    buffer[..., 0] = rgb[..., 2].T          # BGRA, which is what Format_ARGB32 is
    buffer[..., 1] = rgb[..., 1].T
    buffer[..., 2] = rgb[..., 0].T
    buffer[..., 3] = alpha.T
    # Premultiplied would be wrong here: the colour is straight, the alpha
    # separate, which is what Format_ARGB32 means.
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


@dataclass
class FaceStyle:
    """How the faces normal to one axis are coloured.

    One per axis rather than one for the scene, because the three planes of
    a cube are different pictures -- a constant-energy map and two
    dispersions -- and the window that brings out one is rarely right for
    the others.
    """
    colormap: str = "gray"
    flip: bool = False
    levels: tuple = (0.0, 1.0)
    gamma: float = 1.0
    opacity: float = 1.0            # 0..1, the whole face
    fade: float = 0.0               # 0..1 of the level window; 0 = off

    def opaque(self) -> bool:
        return self.opacity >= 0.999 and self.fade <= 0.0


@dataclass
class AxesSpec:
    """The coordinate axes drawn into the scene: which, on which edge of the
    box, which way. ``edges`` None places them automatically, on the
    outline edges facing the viewer (:func:`V.auto_edges`), and they follow
    the camera. See :func:`tools.processing.volume.axis_segments`."""
    shown: list = field(default_factory=lambda: [True, True, True])
    edges: list = None
    directions: list = field(default_factory=lambda: [1, 1, 1])
    titles: tuple = ("x", "y", "z")


#: what a snapshot of a scene keeps: everything :func:`paint_scene` reads
SCENE_FIELDS = ("faces", "bounds", "camera", "levels", "gamma", "colormap",
                "flip", "show_edges", "show_axes", "axis_labels", "background",
                "cull", "image", "styles", "axes", "marks", "frame_color")


class SceneSnapshot:
    """A 3-D view frozen as it is now, for a figure panel: what the view
    shows, independent of the window, which can go on turning or close."""

    def __init__(self, source):
        for name in SCENE_FIELDS:
            setattr(self, name, getattr(source, name, None))
        self.faces = list(self.faces or [])
        self.bounds = list(self.bounds)
        self.camera = copy.deepcopy(self.camera)
        self.styles = ({axis: replace(style) for axis, style in self.styles.items()}
                       if self.styles else None)
        self.axes = copy.deepcopy(self.axes)
        self.marks = copy.deepcopy(self.marks or [])
        if self.image is not None:
            self.image = np.array(self.image, copy=True)


def _fit_for(state, width, height, furniture: float = 0.0):
    """Scale and offset that fit the box into a ``width`` x ``height`` view.

    ``furniture`` is the room, in pixels, the axes' ticks, values and
    titles take beyond the box's edges; it is kept free on every side, so
    the labels stay inside the view (and the figure panel) at any size."""
    bounds = state.bounds
    corners = np.array([[bounds[0][i], bounds[1][j], bounds[2][k]]
                        for i in (0, 1) for j in (0, 1) for k in (0, 1)],
                       dtype=float)
    screen, _ = state.camera.project(corners, bounds)
    lo, hi = screen.min(axis=0), screen.max(axis=0)
    span = np.maximum(hi - lo, 1e-9)
    margin = 0.04 if furniture > 0 else 0.12
    usable = np.array([max(width * (1 - 2 * margin) - 2 * furniture, 1.0),
                       max(height * (1 - 2 * margin) - 2 * furniture, 1.0)])
    scale = float(np.min(usable / span))
    centre = (lo + hi) / 2.0
    offset = np.array([width / 2.0, height / 2.0]) - centre * scale
    return scale, offset


def _shade(face):
    """A little directional shading, so the three walls of a notch are
    told apart even in a single-hue colormap."""
    if not face.outward:
        return 1.0
    normal = face.normal()
    light = np.array([0.35, 0.25, 0.90])
    light /= np.linalg.norm(light)
    return 0.78 + 0.22 * float(abs(np.dot(normal, light)))


def _text_beside(painter, metrics, text, point, normal, angle: float = 0.0,
                 measure=None):
    """Draw ``text`` just beyond ``point`` in the direction ``normal`` (a
    unit screen vector), so it never overlaps the line it labels, turned by
    ``angle`` degrees. Returns how far it reaches along ``normal``."""
    width = metrics.horizontalAdvance(text)
    height = metrics.height()
    theta = np.radians(angle)
    # the text box's extent along the normal, once turned
    along_x = abs(np.cos(theta)) * width + abs(np.sin(theta)) * height
    along_y = abs(np.sin(theta)) * width + abs(np.cos(theta)) * height
    reach = abs(normal[0]) * along_x + abs(normal[1]) * along_y
    centre = np.asarray(point, dtype=float) + normal * (reach / 2.0)
    if measure is not None:
        # where it would be drawn, not drawn: for fitting the picture
        measure.append(QRectF(centre[0] - along_x / 2.0, centre[1] - along_y / 2.0,
                              along_x, along_y))
        return reach
    painter.save()
    painter.translate(QPointF(*centre))
    painter.rotate(angle)
    if "_" in text or "^" in text:
        _draw_markup_centred(painter, text, width, metrics)
    else:
        painter.drawText(QRectF(-width / 2.0, -height / 2.0, width, height),
                         Qt.AlignCenter, text)
    painter.restore()
    return reach


def _draw_markup_centred(painter, text, width, metrics):
    """``text`` with "E_F"-style sub- and superscripts, centred on the
    painter's origin (the figure's own markup rules)."""
    from tools.output.figure import parse_markup

    base = QFont(painter.font())
    small = QFont(base)
    small.setPointSizeF(base.pointSizeF() * 0.7) if base.pointSizeF() > 0 else \
        small.setPixelSize(max(1, int(base.pixelSize() * 0.7)))
    runs = parse_markup(text)
    total = sum(QFontMetricsF(base if kind == "base" else small).horizontalAdvance(run)
                for run, kind in runs)
    x = -total / 2.0
    baseline = 0.5 * (metrics.ascent() - metrics.descent())
    for run, kind in runs:
        font = base if kind == "base" else small
        painter.setFont(font)
        shift = {"sub": 0.3, "super": -0.45}.get(kind, 0.0) * metrics.height()
        painter.drawText(QPointF(x, baseline + shift), run)
        x += QFontMetricsF(font).horizontalAdvance(run)
    painter.setFont(base)


def paint_scene(state, painter, rect, *, font=None, title_font=None,
                line_px: float = 1.0, tick_px: float = None, color=None):
    """Paint a 3-D scene into ``rect`` of any painter.

    ``state`` is a :class:`SceneView` or a :class:`SceneSnapshot` -- the
    window and a figure panel are drawn by this one function, so a figure
    shows exactly what the window did. ``font``, ``title_font``,
    ``line_px``, ``tick_px`` and ``color`` size and colour the axes and
    their labels: screen defaults in the window, the figure's own style
    (font, point sizes, axis line width, tick length) in a figure.
    """
    if not state.faces and state.image is None:
        return
    rect = QRectF(rect)
    width, height = rect.width(), rect.height()
    font = QFont(font) if font is not None else QFont()
    title_font = QFont(title_font) if title_font is not None else QFont(font)
    color = QColor(color) if color is not None else QColor("#202020")
    has_axes = getattr(state, "axes", None) is not None
    scale, offset = _fit_for(state, width, height, 1.0 if has_axes else 0.0)
    fit = {"scale": scale, "offset": offset}

    def to_pixels(screen):
        out = np.asarray(screen, dtype=float) * fit["scale"] + fit["offset"]
        out[:, 1] = height - out[:, 1]
        out[:, 0] += rect.x()
        out[:, 1] += rect.y()
        return out

    def world(points):
        screen, _ = state.camera.project(np.atleast_2d(points), state.bounds)
        return to_pixels(screen)

    if has_axes:
        _fit_with_labels(state, rect, fit, world, font, title_font, tick_px)

    if state.image is not None:
        values = np.asarray(state.image, dtype=float)
        lut = colormaps.get_lut(state.colormap, flip=state.flip)
        qimage = _colour_face(values, lut, state.levels, state.gamma)
        side = min(width, height) * 0.86
        painter.drawImage(QRectF(rect.x() + (width - side) / 2.0,
                                 rect.y() + (height - side) / 2.0, side, side),
                          qimage)

    styles = state.styles or {}
    # Faces turned away cannot be seen through an opaque solid, but they
    # can through a translucent one.
    opaque = all(style.opaque() for style in styles.values())
    visible = V.visible(state.faces, state.camera, state.bounds,
                        cull=bool(state.cull) and opaque)
    # Which of each face's edges belong to the shape rather than to the
    # way it was cut up for drawing. Computed over the visible set, so a
    # seam whose other half is hidden is still recognised as a seam.
    real_edges = dict(zip((id(f) for f in visible), V.outline_edges(visible)))
    luts = {}
    for face, screen in V.sort_faces(visible, state.camera, state.bounds):
        pixels = to_pixels(screen)
        nu, nv = face.image.shape
        style = styles.get(face.axis)
        if style is None:
            style = FaceStyle(state.colormap, state.flip, tuple(state.levels),
                              state.gamma)
        key = (style.colormap, style.flip)
        if key not in luts:
            luts[key] = colormaps.get_lut(style.colormap, flip=style.flip)
        qimage = _colour_face(face.image, luts[key], style.levels, style.gamma,
                              _shade(face), style.opacity, style.fade)
        painter.save()
        painter.setTransform(_face_transform(pixels, nu, nv), True)
        painter.drawImage(QRectF(0, 0, nu, nv), qimage)
        painter.restore()
        frame = getattr(state, "frame_color", "")
        if frame:
            # a stacked slice: its whole outline, so each layer reads as a
            # sheet even where it holds only background
            painter.setPen(QPen(QColor(frame), line_px))
            painter.setBrush(Qt.NoBrush)
            painter.drawPolygon(QPolygonF([QPointF(*p) for p in pixels]))
        elif state.show_edges:
            painter.setPen(QPen(QColor(0, 0, 0, 0x33), line_px))
            for i, draw in enumerate(real_edges[id(face)]):
                if draw:
                    painter.drawLine(QPointF(*pixels[i]),
                                     QPointF(*pixels[(i + 1) % 4]))

    for mark in getattr(state, "marks", None) or []:
        _paint_mark(mark, painter, world, font, line_px)

    if getattr(state, "axes", None) is not None:
        _paint_data_axes(state, painter, world, font, title_font, line_px,
                         tick_px, color)
    if state.show_axes:
        _paint_tripod(state, painter, rect, font)


def _fit_with_labels(state, rect, fit, world, font, title_font, tick_px):
    """Scale the box so that it *and its axis labels* fill ``rect``.

    The labels are a fixed size in pixels while the box scales, so a
    fixed margin for them is either wasted on a large panel or eats a
    small one whole. Here they are measured where they fall: the picture
    is ``scale x box + labels`` across, solved for the scale, then
    centred."""
    bounds = state.bounds
    corners = np.array([[bounds[0][i], bounds[1][j], bounds[2][k]]
                        for i in (0, 1) for j in (0, 1) for k in (0, 1)], float)
    margin = 0.03
    for _round in range(3):
        rects = []
        _paint_data_axes(state, None, world, font, title_font, 1.0, tick_px,
                         QColor("black"), measure=rects)
        box = world(corners)
        lo, hi = box.min(axis=0), box.max(axis=0)
        all_lo, all_hi = lo.copy(), hi.copy()
        for r in rects:
            all_lo = np.minimum(all_lo, [r.left(), r.top()])
            all_hi = np.maximum(all_hi, [r.right(), r.bottom()])
        span = np.maximum(hi - lo, 1e-9)
        extra = (all_hi - all_lo) - span
        room = np.array([rect.width(), rect.height()]) * (1 - 2 * margin) - extra
        factor = float(np.min(np.maximum(room, 0.1 * span) / span))
        fit["scale"] *= factor
        # centre what is drawn (box and labels) in the rect
        box = world(corners)
        rects = []
        _paint_data_axes(state, None, world, font, title_font, 1.0, tick_px,
                         QColor("black"), measure=rects)
        all_lo, all_hi = box.min(axis=0), box.max(axis=0)
        for r in rects:
            all_lo = np.minimum(all_lo, [r.left(), r.top()])
            all_hi = np.maximum(all_hi, [r.right(), r.bottom()])
        centre = 0.5 * (all_lo + all_hi)
        want = np.array([rect.center().x(), rect.center().y()])
        shift = want - centre
        fit["offset"] = fit["offset"] + np.array([shift[0], -shift[1]])
        if abs(factor - 1.0) < 0.01:
            break


def _paint_mark(mark, painter, world, font, line_px):
    """A curve through points of the volume (the Dirac points of a stack of
    cuts, a band followed through the layers), drawn on top, with its dots
    and a label at its first point. ``mark`` is a dict: ``points`` (n, 3)
    in data units, ``color``, ``width`` (x the axis line), ``dots``,
    ``text``."""
    points = np.asarray(mark.get("points") or [], dtype=float).reshape(-1, 3)
    if not len(points):
        return
    pixels = world(points)
    colour = QColor(mark.get("color", "#d62728"))
    pen = QPen(colour, line_px * float(mark.get("width", 2.0)))
    pen.setCapStyle(Qt.RoundCap)
    pen.setJoinStyle(Qt.RoundJoin)
    if mark.get("style") == "dashed":
        pen.setStyle(Qt.DashLine)
    painter.setPen(pen)
    painter.setBrush(Qt.NoBrush)
    if len(pixels) > 1:
        painter.drawPolyline(QPolygonF([QPointF(*p) for p in pixels]))
    if mark.get("dots", True):
        painter.setBrush(colour)
        radius = 1.6 * pen.widthF()
        for p in pixels:
            painter.drawEllipse(QPointF(*p), radius, radius)
    text = str(mark.get("text", "") or "")
    if text:
        painter.setFont(font)
        painter.setPen(colour)
        metrics = QFontMetricsF(font)
        at = pixels[int(mark.get("text_at", 0)) % len(pixels)]
        painter.drawText(QPointF(at[0] + 0.5 * metrics.height(),
                                 at[1] - 0.4 * metrics.height()), text)


def _paint_data_axes(state, painter, world, font, title_font, line_px,
                     tick_px, color, measure=None):
    """The coordinate axes, with round ticks and their values, along the
    lines :attr:`state.axes` chooses. Ticks and labels go on the side of
    each axis away from the middle of the box, where they cannot run into
    the picture they label."""
    from tools.output.export import tick_format, tick_values

    spec = state.axes
    metrics = QFontMetricsF(font)
    title_metrics = QFontMetricsF(title_font)
    tick = float(tick_px) if tick_px else 0.45 * metrics.height()
    gap = 0.3 * metrics.height()
    centre = world([np.mean(b) for b in state.bounds])[0]
    drawing = measure is None
    if drawing:
        painter.setPen(QPen(color, line_px))
    edges = spec.edges or V.auto_edges(state.camera, state.bounds)
    for d, start, end, ticks, step in V.axis_segments(
            state.bounds, edges, spec.directions, spec.shown):
        p0, p1 = world(np.array([start, end]))
        along = p1 - p0
        length = float(np.hypot(*along))
        if length < 2.0:
            continue                         # pointing straight at the viewer
        along /= length
        normal = np.array([-along[1], along[0]])
        middle = (p0 + p1) / 2.0
        if np.dot(normal, middle - centre) < 0:
            normal = -normal
        # as many round ticks as there is room for their numbers
        fmt = tick_format(step)
        label = max(metrics.horizontalAdvance(fmt % v) for v in ticks) if len(ticks) else 0
        room = abs(along[0]) * label + abs(along[1]) * metrics.height()
        target = int(np.clip(length / max(1.6 * room, 1.0), 2, 5))
        lo_v, hi_v = sorted((start[d], end[d]))
        ticks, step = tick_values(lo_v, hi_v, target)
        fmt = tick_format(step)
        if drawing:
            painter.setFont(font)
            painter.drawLine(QPointF(*p0), QPointF(*p1))
        reach = 0.0
        for value in ticks:
            point = start.copy()
            point[d] = value
            p = world(point)[0]
            if drawing:
                painter.drawLine(QPointF(*p), QPointF(*(p + normal * tick)))
            reach = max(reach, _text_beside(painter, metrics, fmt % value,
                                            p + normal * (tick + gap), normal,
                                            measure=measure))
        # The title runs along its axis, never upside down, so a long
        # label beside an upright axis does not push out of the picture.
        angle = float(np.degrees(np.arctan2(along[1], along[0])))
        if angle > 90.0:
            angle -= 180.0
        elif angle < -90.0:
            angle += 180.0
        if drawing:
            painter.setFont(title_font)
        _text_beside(painter, title_metrics, spec.titles[d],
                     middle + normal * (tick + 2 * gap + reach), normal, angle,
                     measure=measure)


def _paint_tripod(state, painter, rect, font):
    """A small labelled tripod in the corner of the view, for orientation.

    Anchored to the view rather than to a corner of the box on purpose: a
    tripod drawn at a box corner is hidden inside the solid for half of all
    viewing angles, and a notched cube is normally looked at from exactly
    those angles.
    """
    width, height = rect.width(), rect.height()
    basis = state.camera.basis()
    arm = min(width, height) * 0.09
    origin = np.array([rect.x() + width - arm * 1.9, rect.y() + height - arm * 1.5])
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
                         state.axis_labels[d])


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
        self.show_axes = True                # the orientation tripod
        self.axis_labels = ("x", "y", "z")
        self.background = "#ffffff"
        self.cull = True
        self.image = None                    # a flat projected image, if any
        self.styles = None                   # {axis: FaceStyle}, or one look
        self.axes = None                     # AxesSpec, or no axes
        self.marks = []                      # curves and labels on top
        self.frame_color = ""                # outline every face (a stack)
        self.setMinimumSize(360, 320)
        self.setMouseTracking(True)
        self._drag = None

    def set_scene(self, faces, bounds):
        self.faces = list(faces)
        self.bounds = list(bounds)
        self.update()

    def snapshot(self) -> SceneSnapshot:
        return SceneSnapshot(self)

    # -- painting ----------------------------------------------------------
    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing, True)
        painter.setRenderHint(QPainter.SmoothPixmapTransform, True)
        painter.fillRect(self.rect(), QColor(self.background))
        self.render_into(painter, self.width(), self.height())
        painter.end()

    def render_into(self, painter, width: int, height: int):
        """Paint the scene into any painter at ``width`` x ``height``."""
        paint_scene(self, painter, QRectF(0, 0, width, height), font=self.font())

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

    def __init__(self, values, axes, labels, parent=None, ranges=None,
                 targets=None):
        super().__init__(parent)
        # Kept as given: a lazy cube is read only by result(), and only
        # where the ranges ask.
        self.values = values if hasattr(values, "shape") else np.asarray(values, float)
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
            points.setValue(int(targets[dim]) if targets is not None
                            else int(min(axis.size, BALANCED_POINTS)))
            points.setKeyboardTracking(False)
            points.setMaximumWidth(90)
            points.setToolTip("How many points to keep on this axis. "
                              "Neighbours are averaged, not thrown away.")
            points.valueChanged.connect(self._update_estimate)
            box.addWidget(QLabel("points"))
            box.addWidget(points)
            box.addStretch(1)
            form.addRow(self.labels[dim], row)
            if ranges is not None and ranges[dim] is not None:
                for widget, value in zip((lo, hi), ranges[dim]):
                    widget.blockSignals(True)
                    widget.setValue(float(value))
                    widget.blockSignals(False)
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
        before = int(np.prod(self.values.shape))
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

    def set_axis(self, axis):
        """A new axis -- the cube was resampled -- keeping the position as
        close as the new points allow. Quietly: the caller redraws."""
        position = self.position()
        self.axis = np.asarray(axis, dtype=float)
        step = abs(float(self.axis[-1] - self.axis[0])) / max(self.axis.size - 1, 1)
        index = int(np.argmin(np.abs(self.axis - position)))
        for widget in (self.slider, self.box):
            widget.blockSignals(True)
        self.slider.setRange(0, self.axis.size - 1)
        self.slider.setValue(index)
        self.box.setRange(float(min(self.axis[0], self.axis[-1])),
                          float(max(self.axis[0], self.axis[-1])))
        self.box.setSingleStep(step or 0.01)
        self.box.setValue(float(self.axis[index]))
        for widget in (self.slider, self.box):
            widget.blockSignals(False)

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


# ==========================================================================
# The 3-D viewers
# ==========================================================================
#: points per axis a 3-D view opens with: the sampling dialog's "Balanced"
BALANCED_POINTS = 160

#: the three views, by the name the launcher opens them with
MODES = {"slices": "3D slices", "notch": "Notched cube",
         "projection": "Volume projection", "stack": "Stacked slices"}

#: the eight corners of the box, as the signs of their coordinates
CORNERS = (((-1, -1, -1), "− − −"), ((1, -1, -1), "+ − −"),
           ((-1, 1, -1), "− + −"), ((1, 1, -1), "+ + −"),
           ((-1, -1, 1), "− − +"), ((1, -1, 1), "+ − +"),
           ((-1, 1, 1), "− + +"), ((1, 1, 1), "+ + +"))


class VolumeWindow(QMainWindow):
    """One of the three 3-D views of a cube: orthogonal slices, the notched
    cube, or a volume projection -- ``mode`` says which.

    Three windows, not three tabs of one: they are different pictures with
    different controls, and a tab hid two of them behind a button labelled
    with the first. Each opens straight away, sampled at the "Balanced"
    setting over the whole cube, so the data is on screen before any
    question is asked; the sampling is shown under the view's buttons, and
    *Change sampling...* resamples from the full cube and redraws.

    ``values``/``axes`` are the cube at full resolution -- a lazy,
    file-backed one is fine, and is only read where the sampling asks.
    """

    closed = pyqtSignal(object)

    def __init__(self, values, axes, labels, source_name: str, parent=None,
                 colormap="gray", flip=False, mode: str = "slices",
                 sampling=None):
        super().__init__(parent)
        if mode not in MODES:
            raise ValueError(f"unknown 3-D view {mode!r}")
        self.mode = mode
        self.source_values = values
        self.source_axes = [np.asarray(a, dtype=float) for a in axes]
        self.labels = list(labels)
        self.short_labels = [l.split(" (")[0] for l in self.labels]
        self.source_name = source_name
        self.setWindowTitle(f"{MODES[mode]} — {source_name}")
        self.resize(1080, 760)

        if sampling is None:
            sampling = ([(float(a.min()), float(a.max())) for a in self.source_axes],
                        [int(min(BALANCED_POINTS, a.size)) for a in self.source_axes])
        self._resample(*sampling)

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
        self.scene.axis_labels = tuple(self.short_labels)
        self.scene.camera.aspect = (1.0, 1.0, 0.75)
        self.scene.rotated.connect(self._camera_moved)
        layout.addWidget(self.scene, 1)

        side = QScrollArea()
        side.setWidgetResizable(True)
        side.setMinimumWidth(460)
        side.setMaximumWidth(490)
        panel = QWidget()
        side.setWidget(panel)
        form = QVBoxLayout(panel)

        if mode == "projection":
            form.addWidget(self._projection_group())
            form.addWidget(self._flat_colour_group())
        else:
            form.addWidget({"slices": self._slice_group,
                            "notch": self._notch_group,
                            "stack": self._stack_group}[mode]())
            self.scene.styles = {d: FaceStyle(colormap, flip, self._auto_levels)
                                 for d in range(3)}
            if mode == "stack":
                # a layer is looked through where it holds only background
                for style in self.scene.styles.values():
                    style.fade = 0.35
                self.scene.camera.elevation = 18.0
                self.scene.frame_color = "#8a8a8a"
            self.scene.axes = AxesSpec(titles=tuple(self.labels))
            form.addWidget(self._colour_group())
            if mode == "stack":
                form.addWidget(self._marks_group())
            form.addWidget(self._axes_group())
        form.addWidget(self._camera_group())
        form.addLayout(self._output_row())
        form.addWidget(self._sampling_group())
        form.addStretch(1)

        self.note = QLabel("")
        self.note.setWordWrap(True)
        form.addWidget(self.note)
        layout.addWidget(side)
        self.setCentralWidget(central)
        self._update_sampling_label()
        self.rebuild()

    # -- sampling ------------------------------------------------------------
    def _resample(self, ranges, targets):
        """Crop and bin the full cube to ``ranges`` and ``targets`` (points
        per axis) -- the cube every view of this window is drawn from."""
        values, axes = V.reduce_volume(self.source_values, self.source_axes,
                                       ranges, targets)
        self.values = np.asarray(values, dtype=float)
        self.axes = [np.asarray(a, dtype=float) for a in axes]
        self.sampling = ([(float(a.min()), float(a.max())) if r is None else tuple(r)
                          for a, r in zip(self.source_axes, ranges)],
                         [a.size if t is None else int(t)
                          for a, t in zip(self.source_axes, targets)])
        finite = self.values[np.isfinite(self.values)]
        self._auto_levels = ((float(np.percentile(finite, 1)),
                              float(np.percentile(finite, 99.5)))
                             if finite.size else (0.0, 1.0))

    def _sampling_group(self):
        group = QGroupBox("Sampling")
        box = QVBoxLayout(group)
        self.sampling_label = QLabel("")
        self.sampling_label.setWordWrap(True)
        self.sampling_label.setTextInteractionFlags(Qt.TextSelectableByMouse)
        box.addWidget(self.sampling_label)
        self.sampling_button = QPushButton("Change sampling...")
        self.sampling_button.setToolTip(
            "Choose the region of the cube and how many points per axis to "
            "keep, then redraw. More points is finer and slower.")
        self.sampling_button.setMaximumWidth(170)
        self.sampling_button.clicked.connect(self.change_sampling)
        box.addWidget(self.sampling_button)
        return group

    def _update_sampling_label(self):
        measured = int(np.prod([a.size for a in self.source_axes]))
        kept = int(self.values.size)
        rows = []
        for d, axis in enumerate(self.axes):
            full = self.source_axes[d]
            whole = (abs(axis.min() - full.min()) <= abs(full[1] - full[0]) if
                     full.size > 1 else True) and \
                    (abs(axis.max() - full.max()) <= abs(full[1] - full[0]) if
                     full.size > 1 else True)
            rows.append(f"{self.short_labels[d]}: {axis.size} of {full.size} pts, "
                        f"{axis.min():.4g} … {axis.max():.4g}"
                        + ("" if whole else " (cropped)"))
        self.sampling_label.setText(
            f"<b>{' × '.join(str(a.size) for a in self.axes)}</b> = "
            f"{kept:,} points, {100.0 * kept / max(measured, 1):.1f}% of the "
            f"{measured:,} measured ({kept * 8 / 1e6:.1f} MB)<br>"
            + "<br>".join(rows))

    def change_sampling(self):
        """Ask for a new region and sampling, resample the full cube to it,
        and redraw -- keeping the positions, colours and camera."""
        ranges, targets = self.sampling
        dialog = VolumeRangeDialog(self.source_values, self.source_axes,
                                   self.labels, self, ranges=ranges,
                                   targets=targets)
        if dialog.exec_() != QDialog.Accepted:
            return False
        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            self._resample(*dialog.selection())
        except Exception as exc:
            QApplication.restoreOverrideCursor()
            QMessageBox.warning(self, "Change sampling", str(exc))
            return False
        QApplication.restoreOverrideCursor()
        for control in getattr(self, "slice_at", []) + getattr(self, "notch_at", []):
            control.set_axis(self.axes[self._control_axis(control)])
        if hasattr(self, "axes_edge"):
            self._axes_edited()
        self._update_sampling_label()
        self.rebuild()
        return True

    def _control_axis(self, control):
        for controls in (getattr(self, "slice_at", []),
                         getattr(self, "notch_at", [])):
            if control in controls:
                return controls.index(control)
        return 0

    # -- the view's own controls ---------------------------------------------
    def _axis_slider(self, dim: int):
        control = AxisPosition(self.axes[dim])
        control.changed.connect(self._queue)
        return control

    def _slice_group(self):
        group = QGroupBox("Slices")
        form = QFormLayout(group)
        form.setFieldGrowthPolicy(QFormLayout.FieldsStayAtSizeHint)
        self.slice_on, self.slice_at = [], []
        for dim in range(3):
            check = QCheckBox(self.short_labels[dim])
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
        return group

    def _notch_group(self):
        group = QGroupBox("Notch")
        form = QFormLayout(group)
        form.setFieldGrowthPolicy(QFormLayout.FieldsStayAtSizeHint)
        self.notch_at = []
        for dim in range(3):
            control = self._axis_slider(dim)
            control.set_position(float(self.axes[dim][int(self.axes[dim].size * 0.6)]))
            form.addRow(self.short_labels[dim], control)
            self.notch_at.append(control)
        self.notch_corner = QComboBox()
        for signs, text in CORNERS[::-1]:
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
        return group

    def _stack_group(self):
        group = QGroupBox("Stack")
        form = QFormLayout(group)
        form.setFieldGrowthPolicy(QFormLayout.FieldsStayAtSizeHint)
        self.stack_axis = QComboBox()
        for d in range(3):
            self.stack_axis.addItem(f"along {self.short_labels[d]}", d)
        self.stack_axis.setCurrentIndex(2)
        self.stack_axis.setToolTip(
            "The axis the layers are stacked along: the energy for "
            "constant-energy maps one above the other, k_z (or the photon "
            "energy) for cuts one behind the other.")
        self.stack_axis.currentIndexChanged.connect(lambda *_: (
            self.even_positions(), self._camera_typed() if hasattr(self, "z_aspect") else None))
        form.addRow("Layers", self.stack_axis)
        self.stack_count = QSpinBox()
        self.stack_count.setRange(1, 60)
        self.stack_count.setValue(8)
        self.stack_count.valueChanged.connect(lambda *_: self.even_positions())
        form.addRow("How many", self.stack_count)
        self.stack_positions = QLineEdit("")
        self.stack_positions.setToolTip(
            "Where each layer is, in the axis's units, separated by commas. "
            "Filled in evenly from the number above; type your own.")
        self.stack_positions.editingFinished.connect(self._queue)
        form.addRow("At", self.stack_positions)
        self.stack_width = spin(0.0, 0.0, 1e6, 0.01, 4)
        self.stack_width.setToolTip("Each layer is the average over this width "
                                    "round its position (0: the nearest "
                                    "sample only).")
        self.stack_width.valueChanged.connect(self._queue)
        form.addRow("Integrate", self.stack_width)
        self.stack_frame = QComboBox()
        self.stack_frame.addItems(["grey", "black", "none"])
        self.stack_frame.setToolTip("An outline round every layer.")
        self.stack_frame.currentIndexChanged.connect(lambda *_: self._restyle())
        form.addRow("Frames", self.stack_frame)
        self.stack_spread = spin(1.0, 0.2, 10.0, 0.5, 1)
        self.stack_spread.setToolTip(
            "Stretch the box along the stacking axis, so the layers stand "
            "apart and each shows whole from an oblique view.")
        self.stack_spread.valueChanged.connect(lambda *_: self._camera_typed())
        form.addRow("Spread", self.stack_spread)
        hint = QLabel("<i>Fade low values (Colour) makes each layer's "
                      "background transparent, so the layers behind show "
                      "through.</i>")
        hint.setWordWrap(True)
        form.addRow("", hint)
        self.even_positions()
        return group

    def even_positions(self):
        """Spread the layers evenly over the stacking axis."""
        axis = self.axes[self.stack_axis.currentData()]
        n = self.stack_count.value()
        lo, hi = float(axis.min()), float(axis.max())
        if n == 1:
            values = [0.5 * (lo + hi)]
        else:
            pad = 0.5 * (hi - lo) / n
            values = np.linspace(lo + pad, hi - pad, n)
        self.stack_positions.setText(", ".join(f"{v:.4g}" for v in values))
        self._queue()

    def set_stack(self, axis: int, positions, width: float = None):
        """Stack along ``axis`` at ``positions`` (for scripts and tests)."""
        self.stack_axis.blockSignals(True)
        self.stack_axis.setCurrentIndex(int(axis))
        self.stack_axis.blockSignals(False)
        self.stack_positions.setText(", ".join(f"{float(v):.6g}" for v in positions))
        if width is not None:
            self.stack_width.setValue(float(width))
        if hasattr(self, "z_aspect"):
            self._camera_typed()
        self.rebuild()

    def _positions(self):
        values = []
        for part in self.stack_positions.text().replace(";", ",").split(","):
            try:
                values.append(float(part))
            except ValueError:
                continue
        return values

    def _marks_group(self):
        group = QGroupBox("Mark a curve through the layers")
        column = QVBoxLayout(group)
        from PyQt5.QtWidgets import QTableWidget, QHeaderView
        self.mark_table = QTableWidget(0, 3)
        self.mark_table.setHorizontalHeaderLabels(self.short_labels)
        self.mark_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.mark_table.setMaximumHeight(150)
        self.mark_table.setToolTip("Points in the data's own units -- a Dirac "
                                   "point on each layer -- joined by a curve "
                                   "drawn over the stack.")
        self.mark_table.itemChanged.connect(lambda *_: self._marks_changed())
        column.addWidget(self.mark_table)
        row = QHBoxLayout()
        add = QPushButton("Add point")
        add.clicked.connect(lambda: self._add_mark_row())
        remove = QPushButton("Remove")
        remove.clicked.connect(lambda: self._remove_mark_row())
        row.addWidget(add)
        row.addWidget(remove)
        self.mark_colour = QComboBox()
        self.mark_colour.addItems(["red", "black", "white", "blue", "green"])
        self.mark_colour.currentIndexChanged.connect(lambda *_: self._marks_changed())
        row.addWidget(self.mark_colour)
        column.addLayout(row)
        self.mark_text = QLineEdit("")
        self.mark_text.setPlaceholderText("Label at the first point (e.g. DP)")
        self.mark_text.editingFinished.connect(self._marks_changed)
        column.addWidget(self.mark_text)
        return group

    def _add_mark_row(self, values=None):
        from PyQt5.QtWidgets import QTableWidgetItem
        row = self.mark_table.rowCount()
        self.mark_table.blockSignals(True)
        self.mark_table.insertRow(row)
        if values is None:
            values = [float(np.mean(a)) for a in self.axes]
        for d, v in enumerate(values):
            self.mark_table.setItem(row, d, QTableWidgetItem(f"{float(v):.6g}"))
        self.mark_table.blockSignals(False)
        self._marks_changed()

    def _remove_mark_row(self):
        row = self.mark_table.currentRow()
        if row < 0:
            row = self.mark_table.rowCount() - 1
        if row >= 0:
            self.mark_table.removeRow(row)
            self._marks_changed()

    def set_marks(self, points, text: str = "", colour: str = "red"):
        """The curve through ``points`` ((n, 3), data units), labelled."""
        self.mark_table.setRowCount(0)
        for p in points:
            self._add_mark_row(p)
        self.mark_text.setText(text)
        self.mark_colour.setCurrentText(colour)
        self._marks_changed()

    def _marks_changed(self):
        points = []
        for row in range(self.mark_table.rowCount()):
            try:
                points.append([float(self.mark_table.item(row, d).text())
                               for d in range(3)])
            except (AttributeError, ValueError):
                continue
        colour = {"red": "#d62728", "black": "#000000", "white": "#ffffff",
                  "blue": "#1f4e9c", "green": "#2ca02c"}[self.mark_colour.currentText()]
        self.scene.marks = ([{"points": points, "color": colour, "width": 2.0,
                              "text": self.mark_text.text()}] if points else [])
        self.scene.update()

    def _projection_group(self):
        group = QGroupBox("Projection")
        form = QFormLayout(group)
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
        return group

    def _level_box(self, value):
        box = spin(value, -1e12, 1e12, max(abs(self._auto_levels[1]) / 50, 1e-6), 4)
        box.setMaximumWidth(140)
        return box

    def _flat_colour_group(self):
        """The projection's colours: one picture, one look."""
        group = QGroupBox("Colour")
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
        self.level_lo = self._level_box(self._auto_levels[0])
        self.level_hi = self._level_box(self._auto_levels[1])
        for box in (self.level_lo, self.level_hi):
            box.valueChanged.connect(self._restyle)
        form.addRow("Level low", self.level_lo)
        form.addRow("Level high", self.level_hi)
        self.gamma = spin(1.0, 0.1, 5.0, 0.1, 2)
        self.gamma.valueChanged.connect(self._restyle)
        form.addRow("Gamma", self.gamma)
        return group

    # -- colours, one look per orientation -----------------------------------
    def _colour_group(self):
        """Colours for the faces normal to each axis, or all of them.

        "All" edits every orientation at once -- only the setting that was
        changed, so a colormap chosen for all keeps each slice's own levels.
        """
        group = QGroupBox("Colour")
        form = QFormLayout(group)
        form.setFieldGrowthPolicy(QFormLayout.FieldsStayAtSizeHint)
        what = "slices" if self.mode == "slices" else "faces"
        self.style_target = QComboBox()
        self.style_target.addItem(f"All {what}", None)
        for d in range(3):
            self.style_target.addItem(
                f"{self.short_labels[d]} = const ({'slice' if what == 'slices' else 'faces'})", d)
        self.style_target.setToolTip(
            "Which planes the settings below apply to. A slice at constant "
            "energy and one at constant angle are different pictures, and "
            "rarely want the same window.")
        self.style_target.currentIndexChanged.connect(self._load_style_controls)
        form.addRow("Apply to", self.style_target)

        self.colormap = QComboBox()
        self.colormap.setMaximumWidth(160)
        self.colormap.addItems(COLORMAP_NAMES)
        form.addRow("Colormap", self.colormap)
        self.flip = QCheckBox("reversed")
        form.addRow("", self.flip)
        self.level_lo = self._level_box(self._auto_levels[0])
        self.level_hi = self._level_box(self._auto_levels[1])
        form.addRow("Level low", self.level_lo)
        form.addRow("Level high", self.level_hi)
        levels_row = QWidget()
        row = QHBoxLayout(levels_row)
        row.setContentsMargins(0, 0, 0, 0)
        auto = QPushButton("Auto")
        auto.setToolTip("Levels from the data of the planes chosen above "
                        "(1st to 99.5th percentile); for all of them, one "
                        "window over all.")
        auto.clicked.connect(self.auto_levels)
        each = QPushButton("Each its own")
        each.setToolTip("Every orientation its own levels, from its own data.")
        each.clicked.connect(self.auto_levels_each)
        row.addWidget(auto)
        row.addWidget(each)
        row.addStretch(1)
        form.addRow("", levels_row)
        self.gamma = spin(1.0, 0.1, 5.0, 0.1, 2)
        form.addRow("Gamma", self.gamma)
        self.opacity = spin(100.0, 0.0, 100.0, 5.0, 0, " %")
        self.opacity.setToolTip("How opaque the planes are. Below 100 % the "
                                "planes behind show through.")
        form.addRow("Opacity", self.opacity)
        self.fade = spin(0.0, 0.0, 100.0, 5.0, 0, " %")
        self.fade.setToolTip(
            "Make the weakest intensities transparent: from fully transparent "
            "at the low level to the opacity above at this share of the level "
            "window. 0 % = off.")
        form.addRow("Fade low values", self.fade)

        for widget, name in ((self.colormap, "colormap"), (self.flip, "flip"),
                             (self.level_lo, "lo"), (self.level_hi, "hi"),
                             (self.gamma, "gamma"), (self.opacity, "opacity"),
                             (self.fade, "fade")):
            signal = (widget.currentIndexChanged if isinstance(widget, QComboBox)
                      else widget.toggled if isinstance(widget, QCheckBox)
                      else widget.valueChanged)
            signal.connect(lambda *_a, n=name: self._style_edited(n))
        self._load_style_controls()
        return group

    def _style_axes(self):
        target = self.style_target.currentData()
        return [0, 1, 2] if target is None else [int(target)]

    def _load_style_controls(self, *_):
        """Show the chosen planes' settings (the first of them, for all)."""
        style = self.scene.styles[self._style_axes()[0]]
        values = ((self.colormap, style.colormap), (self.flip, style.flip),
                  (self.level_lo, style.levels[0]), (self.level_hi, style.levels[1]),
                  (self.gamma, style.gamma), (self.opacity, style.opacity * 100.0),
                  (self.fade, style.fade * 100.0))
        for widget, value in values:
            widget.blockSignals(True)
            if isinstance(widget, QComboBox):
                widget.setCurrentText(value)
            elif isinstance(widget, QCheckBox):
                widget.setChecked(bool(value))
            else:
                widget.setValue(float(value))
            widget.blockSignals(False)

    def _style_edited(self, name):
        for d in self._style_axes():
            style = self.scene.styles[d]
            if name == "colormap":
                style.colormap = self.colormap.currentText()
            elif name == "flip":
                style.flip = self.flip.isChecked()
            elif name == "lo":
                style.levels = (self.level_lo.value(), style.levels[1])
            elif name == "hi":
                style.levels = (style.levels[0], self.level_hi.value())
            elif name == "gamma":
                style.gamma = self.gamma.value()
            elif name == "opacity":
                style.opacity = self.opacity.value() / 100.0
            elif name == "fade":
                style.fade = self.fade.value() / 100.0
        self.scene.update()

    def _plane_levels(self, axis=None):
        """1st to 99.5th percentile of what the planes normal to ``axis``
        show (every plane, for None); the cube's, if none is drawn."""
        images = [np.asarray(f.image, dtype=float).ravel() for f in self.scene.faces
                  if axis is None or f.axis == axis]
        values = np.concatenate(images) if images else self.values.ravel()
        values = values[np.isfinite(values)]
        if values.size == 0:
            return self._auto_levels
        lo, hi = float(np.percentile(values, 1)), float(np.percentile(values, 99.5))
        return (lo, hi) if hi > lo else self._auto_levels

    def auto_levels(self):
        axes = self._style_axes()
        levels = self._plane_levels(None if len(axes) > 1 else axes[0])
        for d in axes:
            self.scene.styles[d].levels = levels
        self._load_style_controls()
        self.scene.update()

    def auto_levels_each(self):
        for d in range(3):
            self.scene.styles[d].levels = self._plane_levels(d)
        self._load_style_controls()
        self.scene.update()

    # -- coordinate axes -------------------------------------------------------
    def _axes_group(self):
        """Where the coordinate axes go: for each axis, the edge of the box
        it is drawn on and the way it runs -- which together say where it
        starts. Automatic puts them on the outline edges facing the viewer,
        and moves them as the view turns."""
        group = QGroupBox("Axes")
        form = QFormLayout(group)
        form.setFieldGrowthPolicy(QFormLayout.FieldsStayAtSizeHint)
        self.axes_on = QCheckBox("Draw the axes, with ticks and values")
        self.axes_on.setChecked(True)
        self.axes_on.toggled.connect(self._axes_edited)
        form.addRow(self.axes_on)

        self.axes_auto = QComboBox()
        self.axes_auto.addItem("Automatic: outer edges facing you", True)
        self.axes_auto.addItem("As chosen below", False)
        self.axes_auto.setToolTip(
            "Automatic keeps each axis on an outline edge of the box, at the "
            "bottom or the left as the view turns. Choose an edge or a "
            "direction below to place them yourself.")
        self.axes_auto.currentIndexChanged.connect(self._axes_edited)
        form.addRow("Placement", self.axes_auto)

        self.axes_show, self.axes_edge, self.axes_direction = [], [], []
        for d in range(3):
            others = [o for o in range(3) if o != d]
            row = QWidget()
            box = QHBoxLayout(row)
            box.setContentsMargins(0, 0, 0, 0)
            show = QCheckBox()
            show.setChecked(True)
            show.setToolTip(f"Draw the {self.short_labels[d]} axis.")
            edge = QComboBox()
            for signs in V.EDGES:
                edge.addItem(", ".join(
                    f"{self._tiny(o)} {'max' if sign > 0 else 'min'}"
                    for o, sign in zip(others, signs)), signs)
            edge.setToolTip(
                f"The edge of the box the {self.short_labels[d]} axis is drawn "
                f"on: where the other two coordinates are at their low (min) "
                f"or high (max) end.")
            direction = QComboBox()
            direction.addItem("from min", 1)
            direction.addItem("from max", -1)
            direction.setToolTip(
                f"Which end of the edge the {self.short_labels[d]} axis starts "
                f"from: its start point is that corner of the box.")
            box.addWidget(show)
            box.addWidget(edge)
            box.addWidget(direction)
            box.addStretch(1)
            form.addRow(self.short_labels[d], row)
            show.toggled.connect(self._axes_edited)
            edge.activated.connect(self._axes_chosen)
            direction.activated.connect(self._axes_chosen)
            self.axes_show.append(show)
            self.axes_edge.append(edge)
            self.axes_direction.append(direction)
        self._axes_edited()
        return group

    def _tiny(self, d):
        """An axis's name cut to a word, for the edge lists."""
        words = self.short_labels[d].replace(",", "").split()
        return words[-1] if words else "xyz"[d]

    def _axes_chosen(self, *_):
        """An edge or a direction picked by hand: placement is manual now."""
        self.axes_auto.blockSignals(True)
        self.axes_auto.setCurrentIndex(1)
        self.axes_auto.blockSignals(False)
        self._axes_edited()

    def _show_auto_edges(self):
        """With automatic placement, the edge lists are greyed out but show
        the edge each axis is on now."""
        if self.scene.axes is None or self.scene.axes.edges is not None:
            return
        for combo, edge in zip(self.axes_edge,
                               V.auto_edges(self.scene.camera, self.scene.bounds)):
            combo.blockSignals(True)
            combo.setCurrentIndex(combo.findData(tuple(edge)))
            combo.blockSignals(False)

    def _axes_edited(self, *_):
        automatic = bool(self.axes_auto.currentData())
        for combo in self.axes_edge:
            combo.setEnabled(not automatic and self.axes_on.isChecked())
        for widget in self.axes_show + self.axes_direction + [self.axes_auto]:
            widget.setEnabled(self.axes_on.isChecked())
        if not self.axes_on.isChecked():
            self.scene.axes = None
        else:
            self.scene.axes = AxesSpec(
                shown=[box.isChecked() for box in self.axes_show],
                edges=None if automatic else
                [tuple(combo.currentData()) for combo in self.axes_edge],
                directions=[combo.currentData() for combo in self.axes_direction],
                titles=tuple(self.labels))
        self._show_auto_edges()
        self.scene.update()

    # -- camera ----------------------------------------------------------------
    def _camera_group(self):
        group = QGroupBox("Camera")
        form = QFormLayout(group)
        form.setFieldGrowthPolicy(QFormLayout.FieldsStayAtSizeHint)
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
        self.edges.setVisible(self.mode != "projection")
        form.addRow("", self.edges)
        self.tripod = QCheckBox("Orientation tripod in the corner")
        self.tripod.setChecked(self.mode == "projection")
        self.tripod.toggled.connect(self._restyle)
        form.addRow("", self.tripod)
        self.scene.show_axes = self.tripod.isChecked()
        return group

    def _output_row(self):
        row = QHBoxLayout()
        save = QPushButton("Save the view...")
        save.setToolTip(
            "Write this view as a figure -- PDF, SVG, PNG or TIFF -- with the "
            "figure composer's default style: 89 mm wide, 600 dpi, its font "
            "and line widths for the axes. For another style, use To figure "
            "composer.")
        # Through a lambda: clicked() passes its `checked` flag, which would
        # land in save_view's `path` (and to_figure_composer's `target`).
        save.clicked.connect(lambda: self.save_view())
        compose = QPushButton("To figure composer...")
        compose.setToolTip(
            "Add this view, as it is now, as a panel of a figure: a new one, "
            "or one already open. There it takes the figure's style and can "
            "be lettered, titled and annotated beside other panels.")
        compose.clicked.connect(lambda: self.to_figure_composer())
        row.addWidget(save)
        row.addWidget(compose)
        row.addStretch(1)
        return row

    # -- behaviour ---------------------------------------------------------
    def _queue(self, *_):
        self._timer.start()

    def _restyle(self, *_):
        if self.mode == "projection":
            self.scene.colormap = self.colormap.currentText()
            self.scene.flip = self.flip.isChecked()
            lo, hi = self.level_lo.value(), self.level_hi.value()
            self.scene.levels = (lo, hi) if hi > lo else self._auto_levels
            self.scene.gamma = self.gamma.value()
            self._queue()                    # the projection uses the levels
        self.scene.show_edges = self.edges.isChecked()
        self.scene.show_axes = self.tripod.isChecked()
        if self.mode == "stack":
            self.scene.frame_color = {"grey": "#8a8a8a", "black": "#000000",
                                      "none": ""}[self.stack_frame.currentText()]
        self.scene.update()

    def _camera_typed(self, *_):
        self.scene.camera.azimuth = self.azimuth.value()
        self.scene.camera.elevation = self.elevation.value()
        self.scene.camera.aspect = self._camera_aspect()
        if self.mode == "projection":
            self._queue()               # a projection depends on the direction
        else:
            if hasattr(self, "axes_edge"):
                self._show_auto_edges()
            self.scene.update()

    def _camera_aspect(self):
        """The box's proportions: the energy axis's from *z aspect*, and
        in a stack the stacking axis stretched by *Spread*, so the layers
        stand apart instead of covering one another."""
        aspect = [1.0, 1.0, self.z_aspect.value()]
        if self.mode == "stack" and hasattr(self, "stack_spread"):
            aspect[self.stack_axis.currentData()] *= self.stack_spread.value()
        return tuple(aspect)

    def _camera_moved(self):
        for box, value in ((self.azimuth, self.scene.camera.azimuth),
                           (self.elevation, self.scene.camera.elevation)):
            box.blockSignals(True)
            box.setValue(value)
            box.blockSignals(False)
        if self.mode == "projection":
            self._queue()
        elif hasattr(self, "axes_edge"):
            self._show_auto_edges()

    def _face_the_notch(self):
        signs = self.notch_corner.currentData()
        self.azimuth.setValue(float(np.degrees(np.arctan2(signs[1], signs[0]))))
        self.elevation.setValue(28.0 * signs[2])

    def rebuild(self):
        bounds = V._bounds(self.axes)
        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            if self.mode == "slices":
                positions = [self.slice_at[d].position()
                             if self.slice_on[d].isChecked() else None
                             for d in range(3)]
                self.scene.image = None
                self.scene.cull = False
                faces = V.slice_faces(self.values, self.axes, positions,
                                      density=360)
                self.scene.set_scene(faces, bounds)
                self.note.setText(f"{len(faces)} plane(s) at " + ", ".join(
                    f"{self.short_labels[d]} = {p:.4g}"
                    for d, p in enumerate(positions) if p is not None))
            elif self.mode == "notch":
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
            elif self.mode == "stack":
                axis = self.stack_axis.currentData()
                self.scene.image = None
                self.scene.cull = False
                faces = V.stack_faces(self.values, self.axes, axis,
                                      self._positions(),
                                      width=self.stack_width.value(), density=300)
                self.scene.set_scene(faces, bounds)
                self.note.setText(f"{len(faces)} layer(s) along "
                                  f"{self.short_labels[axis]}")
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
        lo, hi = self.level_lo.value(), self.level_hi.value()
        if hi <= lo:
            lo, hi = self._auto_levels
        if mode == 3:
            level = lo + float(self.iso_level.value()) * (hi - lo)
            depth, shade = V.isosurface_depth(cube, level)
            self.scene.image = shade
            hits = float(np.isfinite(depth).mean())
            self.note.setText(
                f"Isosurface at {level:.4g}: {100 * hits:.1f}% of the view "
                f"hits it. Shading is Lambert on the depth map's gradient.")
            self.scene.set_scene([], bounds)
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

    # -- output ------------------------------------------------------------
    def figure_panel(self):
        """This view as it is now, as a figure-composer panel."""
        from tools.output.figure import Panel, PanelData

        snapshot = self.scene.snapshot()
        style = (snapshot.styles or {}).get(2)
        colormap = style.colormap if style else snapshot.colormap
        flip = style.flip if style else snapshot.flip
        levels = tuple(style.levels if style else snapshot.levels)
        data = PanelData(array=np.array([levels, levels], dtype=float),
                         x=np.array([0.0, 1.0]), y=np.array([0.0, 1.0]),
                         name=f"{self.source_name} ({MODES[self.mode]})",
                         kind="scene3d", scene=snapshot)
        return Panel(data=data, colormap=colormap, flip=flip, levels=levels,
                     colorbar_label="Intensity")

    def _aspect(self) -> float:
        """Height / width for a figure panel of this view: the shape of the
        box as the camera sees it, not of the window, so the page has no
        bands of white above and below it."""
        bounds = self.scene.bounds
        corners = np.array([[bounds[0][i], bounds[1][j], bounds[2][k]]
                            for i in (0, 1) for j in (0, 1) for k in (0, 1)],
                           dtype=float)
        screen, _ = self.scene.camera.project(corners, bounds)
        span = np.maximum(screen.max(axis=0) - screen.min(axis=0), 1e-9)
        return float(np.clip(span[1] / span[0], 0.45, 1.6))

    def save_view(self, path: str = None):
        """Write the view as a one-panel figure, by the composer's own
        export: the same style, fonts and resolution as any figure."""
        from tools.output.figure import Figure, FigureStyle
        from ui.panels.figure import FigurePainter, write_figure

        if path is None:
            path, _ = QFileDialog.getSaveFileName(
                self, "Save the view", f"{self.source_name}_3d.pdf",
                "PDF (*.pdf);;SVG (*.svg);;PNG (*.png);;TIFF (*.tif)")
            if not path:
                return None
        fmt = path.lower().rsplit(".", 1)[-1] if "." in path else "png"
        fmt = {"tiff": "tif"}.get(fmt, fmt)
        if fmt not in ("pdf", "svg", "png", "tif"):
            fmt, path = "png", path + ".png"
        figure = Figure([self.figure_panel()],
                        style=FigureStyle(panel_aspect=self._aspect()))
        figure.panels[0].label = ""
        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            written = write_figure(FigurePainter(figure), path, fmt)
        except Exception as exc:
            QApplication.restoreOverrideCursor()
            QMessageBox.warning(self, "Save the view", f"Could not write it:\n{exc}")
            return None
        QApplication.restoreOverrideCursor()
        self.note.setText(f"Saved {', '.join(written)}")
        return path

    def to_figure_composer(self, target=None):
        """Put this view into a figure composer: ``target`` (a window), or
        one the user picks among those open, or a new one."""
        from ui.panels.figure import FigureWindow, open_figure_windows
        from ui.viewers.windows import keep_alive

        panel = self.figure_panel()
        if target is None:
            windows = open_figure_windows()
            if windows:
                names = ["A new figure"] + [w.windowTitle() for w in windows]
                choice, ok = QInputDialog.getItem(
                    self, "To figure composer", "Add this view to:", names,
                    0, False)
                if not ok:
                    return None
                index = names.index(choice)
                target = windows[index - 1] if index > 0 else None
        if target is None:
            target = FigureWindow([panel], f"Figure — {self.source_name}",
                                  parent=None)
            target.figure.style.panel_aspect = self._aspect()
            target.refresh()
            keep_alive(target)
            target.show()
        else:
            target.add_panels([panel])
            target.raise_()
            target.activateWindow()
        self.note.setText(f"Added to “{target.windowTitle()}”")
        return target

    def closeEvent(self, event):
        self.closed.emit(self)
        super().closeEvent(event)


# ==========================================================================
# 3-D processing
# ==========================================================================
@fits_screen
class VolumeProcessDialog(QDialog):
    """The same operations as the 2-D panel, over a whole cube.

    Split from :class:`ui.panels.process.ProcessDialog` rather than folded
    into it because the questions are different. On a cube the first
    question is always *which plane is this defined on* -- a curvature of a
    constant-energy map and a curvature of a dispersion are different
    quantities that happen to share a formula -- and answering it on every
    tab of a shared dialog would be worse than having two.
    """

    datasetsCreated = pyqtSignal(list)

    SUFFIXES = {"smooth": "_sm", "derivative": "_d2", "curvature": "_cur",
                "normalise": "_nor", "symmetry": "_symm", "despike": "_fix"}

    #: what each of the launcher's tool buttons opens this panel as: the
    #: plane-wise operations it offers (indices of PLANE_OPERATIONS), or
    #: the symmetrisation, and the window's title
    STARTS = {"smooth": ((0, 1, 2), "Smooth"),
              "derivative": ((3, 4), "Second derivative"),
              "curvature": ((5, 6, 7), "Curvature"),
              "symmetry": (None, "Symmetrise")}

    def __init__(self, datasets, parent=None, colormap="gray", flip=False,
                 start=None):
        super().__init__(parent)
        self.datasets = list(datasets)
        self.colormap, self.flip = colormap, flip
        self.start = start if start in self.STARTS else None
        self.setWindowTitle("Smooth / derivative (3D)")
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

        from ui.panels.process import PreviewPair
        self.preview = PreviewPair()
        layout.addWidget(self.preview, 1)

        self.tabs = QTabWidget()
        self._plane_page = self._plane_tab()
        self._symmetry_page = self._symmetry_tab()
        # One operation per window when a tool button opened it, as in the
        # 2-D panel: the plane-wise page then lists only that operation's
        # variants, and only the settings they use.
        if self.start is None:
            self.tabs.addTab(self._plane_page, "Plane-wise")
            self.tabs.addTab(self._symmetry_page, "Symmetrise")
        else:
            choices, title = self.STARTS[self.start]
            if choices is None:
                self.tabs.addTab(self._symmetry_page, "Symmetrise")
            else:
                self.tabs.addTab(self._plane_page, "Plane-wise")
            self.tabs.tabBar().hide()
            self.setWindowTitle(f"{title} (3D)")
        self.tabs.currentChanged.connect(self.refresh)
        self.tabs.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Maximum)
        layout.addWidget(self.tabs, 0)
        from ui.panels.process import fit_tabs_to_current_page
        fit_tabs_to_current_page(self.tabs)

        self.note = QLabel("")
        self.note.setWordWrap(True)
        layout.addWidget(self.note)

        buttons = ButtonRow()
        preview = add_button(buttons,"Preview", QDialogButtonBox.ActionRole)
        preview.clicked.connect(self.refresh)
        apply_button = add_button(buttons,"Apply to the whole cube",
                                         QDialogButtonBox.AcceptRole)
        # Not connected directly: clicked() passes `checked` (False), which
        # would become apply(background=False) and freeze the window.
        apply_button.clicked.connect(lambda: self.apply())
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
        # Named after the plane's real axes (filled in by _plane_changed),
        # not "first axis"/"second axis", which left the user to work out
        # which of the plane's two axes the program counts first. Each
        # entry carries its index in PLANE_OPERATIONS.
        choices = (self.STARTS[self.start][0] if self.start is not None
                   else None) or range(len(self.PLANE_OPERATIONS))
        for index in choices:
            self.plane_operation.addItem("", index)
        form.addRow("Operation", self._watch(self.plane_operation))
        if self.start in ("derivative", "curvature"):
            from ui.panels.process import smooth_first_hint
            form.insertRow(0, smooth_first_hint(
                "A second derivative" if self.start == "derivative"
                else "The curvature", presmooth=False,
                smoothed=all(P.is_smoothed(d) for _, d in self.datasets)))

        self.plane_width_u = spin(0.0, 0.0, 1e4, 0.001, 5)
        self.plane_width_v = spin(0.0, 0.0, 1e4, 0.001, 5)
        form.addRow("Width along u", self._watch(self.plane_width_u))
        form.addRow("Width along v", self._watch(self.plane_width_v))
        self._plane_form = form
        self.plane_a0 = QDoubleSpinBox()
        self.plane_a0.setDecimals(6)
        self.plane_a0.setRange(1e-6, 1e6)
        self.plane_a0.setValue(1.0)
        form.addRow("a₀", self._watch(self.plane_a0))
        self.plane_norm_mode = QComboBox()
        self.plane_norm_mode.addItems(["area", "max", "mean"])
        form.addRow("Normalise to", self._watch(self.plane_norm_mode))
        self.plane_normal.currentIndexChanged.connect(self._plane_changed)
        self._plane_changed()
        # Only the settings this window's operations use.
        used = {"smooth": (self.plane_width_u, self.plane_width_v),
                "derivative": (), "curvature": (self.plane_a0,)}.get(self.start)
        if used is not None:
            for widget in (self.plane_width_u, self.plane_width_v,
                           self.plane_a0, self.plane_norm_mode):
                if widget not in used:
                    widget.setVisible(False)
                    form.labelForField(widget).setVisible(False)
        note = QLabel("<i>Applied plane by plane, so a curvature of "
                      "constant-energy maps and a curvature of dispersions "
                      "are separate choices rather than the same button.</i>")
        note.setWordWrap(True)
        form.addRow("", note)
        return page

    #: The plane-wise operations, in menu order; ``{u}``/``{v}`` are the
    #: plane's two axes.
    PLANE_OPERATIONS = (
        "Gaussian smooth", "Savitzky-Golay smooth", "Remove spikes",
        "Second derivative along {u}", "Second derivative along {v}",
        "Curvature 2D", "Curvature 1D along {u}", "Curvature 1D along {v}",
        "Normalise each line along {u}", "Normalise each line along {v}")

    def _plane_changed(self, *_):
        """Relabel the operations and widths with the new plane's axes, and
        start the widths at 1/100 of each axis -- the 2-D panel's default --
        so that a smooth applied straight away does something."""
        dim = self.plane_normal.currentData()
        u, v = [self.labels[d].split(" (")[0] for d in range(3) if d != dim]
        for row in range(self.plane_operation.count()):
            text = self.PLANE_OPERATIONS[self.plane_operation.itemData(row)]
            self.plane_operation.setItemText(row, text.format(u=u, v=v))
        for box, name, axis in ((self.plane_width_u, u, self._plane_axes(dim)[0]),
                                (self.plane_width_v, v, self._plane_axes(dim)[1])):
            self._plane_form.labelForField(box).setText(f"Width along {name}")
            span = float(abs(axis[-1] - axis[0])) if axis.size > 1 else 1.0
            box.blockSignals(True)
            box.setMaximum(max(span, 1e-6))
            box.setSingleStep(max(span / 200, 1e-5))
            box.setValue(span / 100)
            box.blockSignals(False)

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
        choice = int(self.plane_operation.currentData())
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
            return P.normalise(plane, axes, 0 if choice == 8 else 1, mode=mode)

        if choice in (0, 1):
            params.update({"width_u": widths[0], "width_v": widths[1]})
        if choice in (5, 6, 7):
            params["a0"] = a0
        if choice in (8, 9):
            params["mode"] = mode
        key = {0: "smooth", 1: "smooth", 2: "despike", 3: "derivative",
               4: "derivative", 5: "curvature", 6: "curvature",
               7: "curvature", 8: "normalise", 9: "normalise"}[choice]
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
        if self.tabs.currentWidget() is self._plane_page:
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
            if self.tabs.currentWidget() is self._plane_page:
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
        from ui.common import jobs

        try:
            prepared = []
            for name, data in self.datasets:
                x, k, z, cube = data.angle_cube
                prepared.append((name, data, np.asarray(cube, dtype=float),
                                 [np.asarray(x, float), np.asarray(k, float),
                                  np.asarray(z, float)]))
        except Exception as exc:
            QMessageBox.warning(self, "Smooth / derivative (3D)",
                                f"Could not read the cube: {exc}")
            return None

        # Every widget is read here, before the worker starts. Reading a
        # spin box from another thread mostly appears to work, which is the
        # worst kind of unsafe; and a value that changed halfway through a
        # run would silently process half the datasets differently.
        try:
            plane_wise = self.tabs.currentWidget() is self._plane_page
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
            QMessageBox.warning(self, "Smooth / derivative (3D)", str(exc))
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
                    self, "Smooth / derivative (3D)", str(exc)),
                on_cancel=lambda: self.note.setText("Processing cancelled."))
            if started:
                self.note.setText("Processing the cube...")
                return None

        try:
            created = work(lambda *a, **k: None)
        except Exception as exc:
            QMessageBox.warning(self, "Smooth / derivative (3D)", str(exc))
            return None
        return self._finish_apply(created)

    def _finish_apply(self, created):
        """Hand the results to the list. On the GUI thread either way."""
        if created:
            self.datasetsCreated.emit(created)
            self.note.setText(f"Added {len(created)} dataset(s) to the list.")
        return created
