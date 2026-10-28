"""
ui/panels/figure.py
===================
The figure composer: the window a frozen slice's "Figure composer..."
button leads to, the figure tools that dress a figure up, and the painter
that draws it. (Not to be confused with a viewer's "Plot tools", which draw
guide lines on the live plot -- see ``ui.viewers.plot_tools``.)

One painter, every output
-------------------------
:class:`FigurePainter` draws a :class:`nf.Figure` onto any Qt paint
device. The preview on screen, the PNG, the TIFF, the vector PDF and the
SVG are the same call with a different device, so a figure cannot look
right in the window and wrong in the file. Text stays text in the PDF and
the SVG; the data itself is a raster embedded at the figure's own dpi,
which is what a journal wants -- an ARPES panel is an image, and pretending
otherwise makes a 40 MB vector file of a million little rectangles.

Everything is laid out in millimetres by ``tools.output.figure`` and multiplied by a
single scale factor here, so "89 mm wide" means 89 mm wide in the PDF and
the same figure at 600 dpi is 2102 px.

What the tools act on
---------------------
:class:`FigureToolsDialog` edits the *selected panel* of the figure, or every
panel at once when "Apply to all panels" is ticked -- which is the honest
form of the MATLAB original's Apply2All, where each button had to
re-implement the loop. Its tabs follow what a figure is made of rather than
what Qt calls things: Panels, Axes, Text, Colour, Overlays, Export.
"""
from __future__ import annotations

import copy
import os
import time
import uuid

import numpy as np
from PyQt5.QtCore import Qt, QRectF, QPointF, QSize, pyqtSignal
from PyQt5.QtGui import (QColor, QCursor, QFont, QFontMetricsF, QImage,
                         QKeySequence, QPainter, QPainterPath, QPen, QBrush,
                         QPolygonF, QPdfWriter, QPageSize)
from PyQt5.QtWidgets import (QWidget, QMainWindow, QVBoxLayout, QHBoxLayout,
                             QLabel, QPushButton, QComboBox, QCheckBox,
                             QDoubleSpinBox, QSpinBox, QLineEdit, QFormLayout,
                             QTabWidget, QDialog, QDialogButtonBox, QGroupBox,
                             QFileDialog, QMessageBox, QListWidget, QScrollArea,
                             QFrame, QApplication, QRadioButton, QInputDialog,
                             QSizePolicy, QColorDialog, QListWidgetItem,
                             QMenu, QShortcut, QSplitter, QGridLayout)

from tools.output import colormaps
from tools.output import export
from tools.output import figfile
import tools.output.figure as nf
from tools.output.figure import (Figure, FigureStyle, Panel, PanelData, Line,
                                 TextNote, Arrow, ScaleBar, Inset, Overlay,
                                 JOURNAL_PRESETS, parse_markup, pt_to_mm)

#: line styles, as the tools name them
PEN_STYLES = {"solid": Qt.SolidLine, "dashed": Qt.DashLine,
              "dotted": Qt.DotLine, "dashdot": Qt.DashDotLine}

#: where a thing can sit inside a panel
CORNERS = ("top left", "top right", "bottom left", "bottom right",
           "top center", "bottom center")


# --------------------------------------------------------------------------
# Painting
# --------------------------------------------------------------------------
class FigurePainter:
    """Draws a :class:`nf.Figure`.

    ``scale`` is device units per millimetre: ``dpi / 25.4`` for a raster,
    ``72 / 25.4`` for a PDF in points. Fonts are sized in device units from
    their point size for the same reason -- one conversion, applied once.
    """

    def __init__(self, figure: Figure):
        self.figure = figure
        #: device units per mm to generate the data raster at; None means
        #: the figure's own dpi (see :meth:`_raster`)
        self._raster_scale = None
        #: draw the empty slots (on screen only: never in a file)
        self.placeholders = False

    # -- sizes -------------------------------------------------------------
    def size_mm(self):
        width, height, _rects = self.figure.layout()
        return width, height

    def pixel_size(self, dpi: float = None):
        dpi = float(dpi or self.figure.style.dpi)
        width, height = self.size_mm()
        scale = dpi / nf.MM_PER_INCH
        return max(int(round(width * scale)), 1), max(int(round(height * scale)), 1)

    # -- small helpers ------------------------------------------------------
    def _font(self, size_pt: float, scale: float, bold: bool = False) -> QFont:
        font = QFont(self.figure.style.font)
        font.setPixelSize(max(1, int(round(pt_to_mm(size_pt) * scale))))
        font.setBold(bool(bold))
        return font

    def _pen(self, color, width_pt: float, scale: float, style="solid") -> QPen:
        pen = QPen(QColor(color))
        pen.setWidthF(max(pt_to_mm(width_pt) * scale, 0.1))
        pen.setStyle(PEN_STYLES.get(style, Qt.SolidLine))
        pen.setCapStyle(Qt.FlatCap)
        return pen

    @staticmethod
    def _runs_width(painter, runs, base: QFont, small: QFont) -> float:
        total = 0.0
        for text, kind in runs:
            painter.setFont(base if kind == "base" else small)
            total += QFontMetricsF(painter.font()).horizontalAdvance(text)
        return total

    def _draw_markup(self, painter, x, y, text, size_pt, scale, color,
                     anchor="left", bold=False, rotate=0.0):
        """Draw a label with ``_``/``^`` honoured, anchored at ``(x, y)``
        on its baseline. Returns the width drawn."""
        runs = parse_markup(text)
        if not runs:
            return 0.0
        base = self._font(size_pt, scale, bold)
        small = self._font(size_pt * 0.72, scale, bold)
        painter.save()
        painter.setPen(QPen(QColor(color)))
        if rotate:
            painter.translate(x, y)
            painter.rotate(rotate)
            x = y = 0.0
        width = self._runs_width(painter, runs, base, small)
        cursor = x - (width / 2 if anchor == "center" else width if anchor == "right" else 0)
        shift = 0.34 * base.pixelSize()
        for run, kind in runs:
            painter.setFont(base if kind == "base" else small)
            offset = 0.0 if kind == "base" else (shift if kind == "sub" else -shift * 1.15)
            painter.drawText(QPointF(cursor, y + offset), run)
            cursor += QFontMetricsF(painter.font()).horizontalAdvance(run)
        painter.restore()
        return width

    def _markup_width(self, painter, text, size_pt, scale, bold=False) -> float:
        runs = parse_markup(text)
        return self._runs_width(painter, runs, self._font(size_pt, scale, bold),
                                self._font(size_pt * 0.72, scale, bold))

    @staticmethod
    def label_centre(centre: float, width: float, page_width: float = None) -> float:
        """Where to centre an axis label of ``width``: under its panel, but
        moved in rather than cut off when it is wider than the panel and the
        panel is at the edge of the page -- "Angle, deflector (de" on the
        right-hand column of a narrow figure."""
        if not page_width or width >= page_width:
            return centre
        return float(min(max(centre, width / 2.0), page_width - width / 2.0))

    @staticmethod
    def _ticks(lo, hi, step, target=4):
        if step and step > 0:
            first = np.ceil(lo / step - 1e-9) * step
            n = int(np.floor((hi - first) / step + 1e-9)) + 1
            values = first + step * np.arange(max(n, 0))
            return values[(values >= lo - 1e-9) & (values <= hi + 1e-9)], step
        return export.tick_values(lo, hi, target)

    # -- the whole figure ---------------------------------------------------
    def paint(self, painter: QPainter, scale: float, selection: int = None):
        """Draw at ``scale`` device units per millimetre. ``selection``
        outlines one panel; exports pass None so the marker never reaches a
        file."""
        figure = self.figure
        style = figure.style
        width, height, rects = figure.layout()
        self._page_width = width * scale
        painter.setRenderHint(QPainter.Antialiasing, True)
        painter.setRenderHint(QPainter.TextAntialiasing, True)
        painter.setRenderHint(QPainter.SmoothPixmapTransform, True)
        painter.fillRect(QRectF(0, 0, width * scale, height * scale),
                         QColor(style.background))

        self._align_labels(rects, scale)
        for rect in rects:
            self._paint_panel(painter, figure.panels[rect["index"]], rect,
                              scale, rect["index"])
        self._draw_zoom_marks(painter, rects, scale)
        for shape in getattr(figure, "shapes", ()):
            self.paint_shape(painter, shape, scale)
        if figure.title:
            self._draw_markup(
                painter, width * scale / 2,
                (style.margin_top_mm * 0.9) * scale,
                figure.title, (style.label_pt or style.font_pt) + 1, scale,
                style.axis_color, anchor="center", bold=True)
        if selection is not None and 0 <= selection < len(rects):
            rect = rects[selection]
            pen = QPen(QColor("#2ca6cd"))
            pen.setWidthF(max(1.2, 0.25 * scale))
            pen.setStyle(Qt.DashLine)
            painter.setPen(pen)
            painter.setBrush(Qt.NoBrush)
            painter.drawRect(QRectF((rect["x"] - 0.6) * scale, (rect["y"] - 0.6) * scale,
                                    (rect["w"] + 1.2) * scale, (rect["h"] + 1.2) * scale))

    # -- one panel ----------------------------------------------------------
    def _paint_panel(self, painter, panel: Panel, rect, scale, index):
        """One panel, in its own font size when it has one (the figure's
        sizes are swapped for the panel's while it is drawn)."""
        style = self.figure.style
        if not panel.font_pt:
            return self._paint_panel_body(painter, panel, rect, scale, index)
        kept = (style.font_pt, style.tick_pt, style.label_pt,
                style.annotation_pt, style.panel_label_pt)
        style.font_pt = float(panel.font_pt)
        style.tick_pt = style.label_pt = style.annotation_pt = 0.0
        style.panel_label_pt = kept[4] or (kept[0] + 1)
        try:
            return self._paint_panel_body(painter, panel, rect, scale, index)
        finally:
            (style.font_pt, style.tick_pt, style.label_pt,
             style.annotation_pt, style.panel_label_pt) = kept

    def _paint_panel_body(self, painter, panel: Panel, rect, scale, index):
        figure, style = self.figure, self.figure.style
        (x0, x1), (y0, y1) = figure.panel_ranges(panel)
        levels = figure.panel_levels(panel)
        box = QRectF(rect["x"] * scale, rect["y"] * scale,
                     rect["w"] * scale, rect["h"] * scale)

        if panel.data.kind == "empty":
            self._paint_empty(painter, panel, box, scale)
            return
        if panel.data.kind == "picture":
            self._paint_picture(painter, panel, box, scale)
            return
        if panel.data.kind == "scene3d":
            self._paint_scene_panel(painter, panel, box, scale, levels)
            return
        self._draw_image(painter, panel, box, (x0, x1), (y0, y1), levels, scale)
        painter.setClipRect(box)
        self._draw_overlays(painter, panel, box, (x0, x1), (y0, y1), scale)
        painter.setClipping(False)

        self._draw_frame_and_ticks(painter, panel, rect, box, scale, index,
                                   (x0, x1), (y0, y1))
        if panel.inset is not None:
            self._draw_inset(painter, panel, box, scale, levels)
        if panel.scale_bar is not None:
            self._draw_scale_bar(painter, panel, box, scale, (x0, x1))
        if panel.colorbar:
            self._draw_colorbar(painter, panel, rect, box, scale, levels)
        self._draw_panel_label(painter, panel, box, scale)

    def _paint_empty(self, painter, panel, box, scale):
        """An empty slot: a dashed outline with its letter on screen, and
        nothing at all in a file."""
        if not self.placeholders:
            return
        pen = QPen(QColor("#9aa3ad"))
        pen.setWidthF(max(1.0, 0.25 * scale))
        pen.setStyle(Qt.DashLine)
        painter.setPen(pen)
        painter.setBrush(QBrush(QColor(240, 243, 247, 160)))
        painter.drawRect(box)
        size = 7.0 if box.width() / scale > 18 else 5.0
        text = ("empty — double-click to load" if box.width() / scale > 45
                else "empty")
        self._draw_markup(painter, box.center().x(),
                          box.center().y() + pt_to_mm(size) * scale / 3,
                          text, size, scale, "#7a828b", anchor="center")
        self._draw_panel_label(painter, panel, box, scale)

    def _paint_picture(self, painter, panel, box, scale):
        """A picture from a file, as large as fits the frame without
        changing its shape, centred in it."""
        rgba = np.ascontiguousarray(panel.data.array, dtype=np.uint8)
        if rgba.ndim == 3 and rgba.shape[0] > 0 and rgba.shape[1] > 0:
            height, width = rgba.shape[:2]
            image = QImage(rgba.data, width, height, width * 4,
                           QImage.Format_RGBA8888)
            fit = min(box.width() / width, box.height() / height)
            target = QRectF(0, 0, width * fit, height * fit)
            target.moveCenter(box.center())
            painter.drawImage(target, image)
        self._draw_panel_label(painter, panel, box, scale)

    def _draw_zoom_marks(self, painter, rects, scale):
        """On each panel an inset magnifies, a dashed box round the region
        the inset shows."""
        figure, style = self.figure, self.figure.style
        boxes = {figure.panels[r["index"]].uid: r for r in rects}
        for panel in figure.panels:
            if panel.role != "inset" or not panel.zoom_of:
                continue
            host = figure.panel_by_uid(panel.zoom_of)
            rect = boxes.get(panel.zoom_of)
            if host is None or rect is None or not host.data.is_plot():
                continue
            box = QRectF(rect["x"] * scale, rect["y"] * scale,
                         rect["w"] * scale, rect["h"] * scale)
            xr, yr = figure.panel_ranges(host)
            (a0, a1), (b0, b1) = panel.ranges()
            corner0 = self._map_point(box, host, xr, yr, a0, b0)
            corner1 = self._map_point(box, host, xr, yr, a1, b1)
            painter.save()
            painter.setClipRect(box)
            painter.setPen(self._pen(style.axis_color if style.axis_color != "#000000"
                                     else "#ffffff", style.axis_width_pt, scale,
                                     "dashed"))
            painter.setBrush(Qt.NoBrush)
            painter.drawRect(QRectF(corner0, corner1).normalized())
            painter.restore()

    def paint_shape(self, painter, shape, scale):
        """One arrow, line, box, ellipse or note of the page (mm -> device
        units by ``scale``)."""
        style = self.figure.style
        start = QPointF(shape.x0 * scale, shape.y0 * scale)
        end = QPointF(shape.x1 * scale, shape.y1 * scale)
        painter.save()
        if shape.kind == "text":
            size = shape.size_pt or style.annotation_pt or style.font_pt
            self._draw_markup(painter, start.x(), start.y(), shape.text or "text",
                              size, scale, shape.color, bold=shape.bold)
            painter.restore()
            return
        pen = self._pen(shape.color, shape.width_pt, scale, shape.style)
        pen.setCapStyle(Qt.RoundCap)
        painter.setPen(pen)
        if shape.kind in ("rect", "ellipse"):
            painter.setBrush(QBrush(QColor(shape.fill)) if shape.fill else Qt.NoBrush)
            area = QRectF(start, end).normalized()
            if shape.kind == "rect":
                painter.drawRect(area)
            else:
                painter.drawEllipse(area)
        elif shape.kind == "line":
            painter.drawLine(start, end)
        elif shape.kind == "arrow":
            self._page_arrow(painter, start, end, shape, scale)
        painter.restore()

    def _page_arrow(self, painter, start, end, shape, scale):
        vector = np.array([end.x() - start.x(), end.y() - start.y()])
        length = float(np.hypot(*vector))
        if length < 1e-6:
            return
        unit = vector / length
        head = min(shape.head_mm * scale, 0.45 * length)
        normal = np.array([-unit[1], unit[0]]) * head * 0.38
        tips = [(end, unit)] + ([(start, -unit)] if shape.double else [])
        a = QPointF(start.x() + unit[0] * head, start.y() + unit[1] * head) \
            if shape.double else start
        b = QPointF(end.x() - unit[0] * head, end.y() - unit[1] * head)
        painter.drawLine(a, b)
        painter.setPen(Qt.NoPen)
        painter.setBrush(QBrush(QColor(shape.color)))
        for tip, direction in tips:
            base = QPointF(tip.x() - direction[0] * head, tip.y() - direction[1] * head)
            painter.drawPolygon(QPolygonF([
                tip, QPointF(base.x() + normal[0], base.y() + normal[1]),
                QPointF(base.x() - normal[0], base.y() - normal[1])]))

    def text_box(self, shape, scale):
        """Where a note is drawn, in device units, for picking it."""
        style = self.figure.style
        size = shape.size_pt or style.annotation_pt or style.font_pt
        font = self._font(size, scale, shape.bold)
        metrics = QFontMetricsF(font)
        width = max(metrics.horizontalAdvance(nf.plain_text(shape.text or "text")), 4.0)
        return QRectF(shape.x0 * scale, shape.y0 * scale - metrics.ascent(),
                      width, metrics.height())

    def _paint_scene_panel(self, painter, panel, box, scale, levels):
        """A 3-D view (ui.panels.volume) as a panel: drawn by the same code
        as on screen, with this figure's font, sizes and line width, so its
        axes are labelled like every other panel's. No frame -- the axes are
        in the picture -- but notes, arrows, a colour bar, the title and the
        panel letter as usual."""
        from ui.panels.volume import paint_scene

        style = self.figure.style
        tick_pt = style.tick_pt or style.font_pt
        label_pt = style.label_pt or style.font_pt
        painter.save()
        painter.setClipRect(box)
        drawn = paint_scene(panel.data.scene, painter, box,
                    font=self._font(tick_pt, scale),
                    title_font=self._font(label_pt, scale),
                    line_px=max(pt_to_mm(style.axis_width_pt) * scale, 0.5),
                    tick_px=style.tick_len_mm * scale,
                    color=QColor(style.axis_color))
        painter.restore()
        ranges = ((0.0, 1.0), (0.0, 1.0))
        painter.setClipRect(box)
        self._draw_overlays(painter, panel, box, *ranges, scale)
        painter.setClipping(False)
        if panel.title:
            metrics = QFontMetricsF(self._font(label_pt, scale))
            self._draw_markup(painter, box.center().x(),
                              box.top() - metrics.descent(), panel.title,
                              label_pt, scale, style.axis_color, anchor="center")
        if panel.colorbar:
            # against the cube as drawn, not the whole frame round it: as
            # tall as the cube beside it, as wide as the cube under the frame
            # (below the cube's own axis labels)
            against = box
            if drawn is not None:
                drawn = drawn.intersected(box)
                if getattr(panel, "colorbar_position", "right") == "below":
                    against = QRectF(drawn.left(), box.top(), drawn.width(), box.height())
                else:
                    against = QRectF(box.left(), drawn.top(), box.width(), drawn.height())
            self._draw_colorbar(painter, panel, None, against, scale, levels)
        self._draw_panel_label(painter, panel, box, scale)

    #: no raster wider or taller than this, whatever the dpi asks for
    MAX_RASTER = 8000

    def _raster(self, panel, box, x_range, y_range, levels, scale):
        """The panel's data as pixels.

        Sized from the panel's size in *millimetres* times the figure's dpi,
        not from the device units it is being drawn at: that is what puts a
        600-dpi image inside a vector PDF whose own units are points, and
        what keeps a preview cheap without changing the file.
        """
        target = self._raster_scale or (self.figure.style.dpi / nf.MM_PER_INCH)
        width = max(int(round(box.width() / scale * target)), 1)
        height = max(int(round(box.height() / scale * target)), 1)
        width, height = min(width, self.MAX_RASTER), min(height, self.MAX_RASTER)
        try:
            lut = colormaps.get_lut(panel.colormap, flip=panel.flip)
        except Exception:
            lut = colormaps.get_lut("gray")
        return export.render_rgba(
            panel.data.array, panel.data.x, panel.data.y,
            x_range=x_range, y_range=y_range, lut=lut, levels=levels,
            gamma=panel.gamma, width=width, height=height,
            interpolate=panel.smooth, invert_x=panel.invert_x,
            invert_y=panel.invert_y)

    def _draw_image(self, painter, panel, box, x_range, y_range, levels, scale):
        rgba = self._raster(panel, box, x_range, y_range, levels, scale)
        buffer = np.ascontiguousarray(rgba)
        image = QImage(buffer.data, buffer.shape[1], buffer.shape[0],
                       buffer.shape[1] * 4, QImage.Format_RGBA8888)
        painter.drawImage(box, image)

    # -- the axis furniture -------------------------------------------------
    def _draw_frame_and_ticks(self, painter, panel, rect, box, scale, index,
                              x_range, y_range):
        figure, style = self.figure, self.figure.style
        tick_pt = style.tick_pt or style.font_pt
        label_pt = style.label_pt or style.font_pt
        inward = -1.0 if style.tick_direction == "in" else 1.0
        tick = style.tick_len_mm * scale
        gap = 0.35 * pt_to_mm(tick_pt) * scale
        painter.setBrush(Qt.NoBrush)
        painter.setPen(self._pen(style.axis_color, style.axis_width_pt, scale))
        painter.drawRect(box)

        all_sides = style.tick_sides == "all"
        (x0, x1), (y0, y1) = x_range, y_range

        def px(value):
            frac = (value - x0) / (x1 - x0) if x1 > x0 else 0.0
            if panel.invert_x:
                frac = 1.0 - frac
            return box.left() + frac * box.width()

        def py(value):
            frac = (value - y0) / (y1 - y0) if y1 > y0 else 0.0
            if panel.invert_y:
                frac = 1.0 - frac
            return box.bottom() - frac * box.height()

        target_x, target_y = figure.tick_targets(rect["w"], rect["h"])
        x_ticks, x_step = self._ticks(x0, x1, style.x_tick_step, target_x)
        named = [(float(v), str(t)) for v, t in (panel.x_tick_labels or ())]
        if named:
            # A band-structure path: one tick per labelled point, no
            # numbers and no minor ticks between them.
            x_ticks = np.array([v for v, _t in named])
            x_step = 0.0
        y_ticks, y_step = self._ticks(y0, y1, style.y_tick_step, target_y)
        x_fmt, y_fmt = export.tick_format(x_step), export.tick_format(y_step)

        def draw_ticks(values, length):
            for value in values:
                if value < x0 - 1e-9 or value > x1 + 1e-9:
                    continue
                x = px(value)
                painter.drawLine(QPointF(x, box.bottom()),
                                 QPointF(x, box.bottom() + inward * length))
                if all_sides:
                    painter.drawLine(QPointF(x, box.top()),
                                     QPointF(x, box.top() - inward * length))
            return values

        def draw_yticks(values, length):
            for value in values:
                if value < y0 - 1e-9 or value > y1 + 1e-9:
                    continue
                y = py(value)
                painter.drawLine(QPointF(box.left(), y),
                                 QPointF(box.left() - inward * length, y))
                if all_sides:
                    painter.drawLine(QPointF(box.right(), y),
                                     QPointF(box.right() + inward * length, y))

        draw_ticks(x_ticks, tick)
        draw_yticks(y_ticks, tick)
        if style.minor_ticks > 0:
            minor = style.minor_ticks
            sub_x = np.concatenate([x_ticks + x_step * (i + 1) / (minor + 1)
                                    for i in range(minor)]
                                   + [x_ticks - x_step * (i + 1) / (minor + 1)
                                      for i in range(minor)]) if (
                x_ticks.size and not named) else np.zeros(0)
            sub_y = np.concatenate([y_ticks + y_step * (i + 1) / (minor + 1)
                                    for i in range(minor)]
                                   + [y_ticks - y_step * (i + 1) / (minor + 1)
                                      for i in range(minor)]) if y_ticks.size else y_ticks
            draw_ticks(sub_x, tick * 0.55)
            draw_yticks(sub_y, tick * 0.55)

        metrics = QFontMetricsF(self._font(tick_pt, scale))
        outward = tick if style.tick_direction == "out" else 0.0
        if figure.shows_x_labels(index):
            baseline = box.bottom() + outward + gap + metrics.ascent()
            names = dict(named)
            for value in x_ticks:
                text = names.get(float(value), None) if named else None
                self._draw_markup(painter, px(value), baseline,
                                  x_fmt % value if text is None else text,
                                  tick_pt, scale, style.axis_color, anchor="center")
            label_width = self._markup_width(painter, panel.label_x(), label_pt, scale)
            self._draw_markup(painter,
                              self.label_centre(box.center().x(), label_width,
                                                getattr(self, "_page_width", None)),
                              baseline + 1.25 * metrics.height(),
                              panel.label_x(), label_pt, scale, style.axis_color,
                              anchor="center")
        if figure.shows_y_labels(index):
            widest = 0.0
            for value in y_ticks:
                text = y_fmt % value
                widest = max(widest, self._markup_width(painter, text, tick_pt, scale))
                self._draw_markup(painter, box.left() - outward - gap,
                                  py(value) + metrics.ascent() / 2 - metrics.descent() / 2,
                                  text, tick_pt, scale, style.axis_color, anchor="right")
            self._draw_markup(
                painter, box.left() - outward - gap - widest - 0.45 * metrics.height(),
                box.center().y(), panel.label_y(), label_pt, scale,
                style.axis_color, anchor="center", rotate=-90.0)
        if panel.title:
            self._draw_markup(painter, box.center().x(),
                              box.top() - outward - gap - metrics.descent(),
                              panel.title, label_pt, scale, style.axis_color,
                              anchor="center")
        self._draw_secondary(painter, panel, box, scale, x_range, y_range,
                             tick_pt, label_pt, inward, tick, gap)

    def _draw_secondary(self, painter, panel, box, scale, x_range, y_range,
                        tick_pt, label_pt, inward, tick, gap):
        """A second scale on the far edge: binding against kinetic energy,
        Å⁻¹ against degrees, k in units of π/a."""
        style = self.figure.style
        metrics = QFontMetricsF(self._font(tick_pt, scale))
        painter.setPen(self._pen(style.axis_color, style.axis_width_pt, scale))

        axis = panel.top_axis
        if axis is not None and axis.enabled and axis.scale:
            x0, x1 = x_range
            lo, hi = sorted((axis.scale * x0 + axis.offset,
                             axis.scale * x1 + axis.offset))
            values, step = export.tick_values(lo, hi, 4)
            fmt = export.tick_format(step)
            for value in values:
                primary = (value - axis.offset) / axis.scale
                frac = (primary - x0) / (x1 - x0) if x1 > x0 else 0.0
                if panel.invert_x:
                    frac = 1.0 - frac
                if not (-1e-9 <= frac <= 1 + 1e-9):
                    continue
                x = box.left() + frac * box.width()
                painter.drawLine(QPointF(x, box.top()),
                                 QPointF(x, box.top() - inward * tick))
                self._draw_markup(painter, x, box.top() - gap - metrics.descent(),
                                  fmt % value, tick_pt, scale, style.axis_color,
                                  anchor="center")
            if axis.label:
                self._draw_markup(painter, box.center().x(),
                                  box.top() - gap - 1.35 * metrics.height(),
                                  axis.label, label_pt, scale, style.axis_color,
                                  anchor="center")

        axis = panel.right_axis
        if axis is not None and axis.enabled and axis.scale:
            y0, y1 = y_range
            lo, hi = sorted((axis.scale * y0 + axis.offset,
                             axis.scale * y1 + axis.offset))
            values, step = export.tick_values(lo, hi, 4)
            fmt = export.tick_format(step)
            widest = 0.0
            for value in values:
                primary = (value - axis.offset) / axis.scale
                frac = (primary - y0) / (y1 - y0) if y1 > y0 else 0.0
                if panel.invert_y:
                    frac = 1.0 - frac
                if not (-1e-9 <= frac <= 1 + 1e-9):
                    continue
                y = box.bottom() - frac * box.height()
                painter.drawLine(QPointF(box.right(), y),
                                 QPointF(box.right() + inward * tick, y))
                text = fmt % value
                widest = max(widest, self._markup_width(painter, text, tick_pt, scale))
                self._draw_markup(painter, box.right() + gap,
                                  y + metrics.ascent() / 2 - metrics.descent() / 2,
                                  text, tick_pt, scale, style.axis_color)
            if axis.label:
                self._draw_markup(painter, box.right() + gap + widest
                                  + 0.9 * metrics.height(), box.center().y(),
                                  axis.label, label_pt, scale, style.axis_color,
                                  anchor="center", rotate=-90.0)

    # -- what is drawn on top -----------------------------------------------
    def _map_point(self, box, panel, x_range, y_range, x, y):
        (x0, x1), (y0, y1) = x_range, y_range
        fx = (x - x0) / (x1 - x0) if x1 > x0 else 0.0
        fy = (y - y0) / (y1 - y0) if y1 > y0 else 0.0
        if panel.invert_x:
            fx = 1.0 - fx
        if panel.invert_y:
            fy = 1.0 - fy
        return QPointF(box.left() + fx * box.width(), box.bottom() - fy * box.height())

    def _draw_overlays(self, painter, panel, box, x_range, y_range, scale):
        style = self.figure.style
        (x0, x1), (y0, y1) = x_range, y_range
        for line in panel.lines:
            painter.setPen(self._pen(line.color, line.width_pt, scale, line.style))
            if line.axis == "x":
                point = self._map_point(box, panel, x_range, y_range, line.value, y0)
                painter.drawLine(QPointF(point.x(), box.top()),
                                 QPointF(point.x(), box.bottom()))
            else:
                point = self._map_point(box, panel, x_range, y_range, x0, line.value)
                painter.drawLine(QPointF(box.left(), point.y()),
                                 QPointF(box.right(), point.y()))
            if line.label:
                size = style.annotation_pt or style.font_pt
                if line.axis == "x":
                    self._draw_markup(painter, point.x() + 0.5 * scale,
                                      box.top() + 1.6 * pt_to_mm(size) * scale,
                                      line.label, size, scale, line.color)
                else:
                    self._draw_markup(painter, box.left() + 0.6 * scale,
                                      point.y() - 0.5 * scale, line.label, size,
                                      scale, line.color)
        for overlay in panel.overlays:
            self._draw_overlay(painter, panel, box, x_range, y_range, overlay,
                               scale)
        for arrow in panel.arrows:
            self._draw_arrow(painter, panel, box, x_range, y_range, arrow, scale)
        for note in panel.texts:
            size = note.size_pt or style.annotation_pt or style.font_pt
            point = self._map_point(box, panel, x_range, y_range, note.x, note.y)
            self._draw_markup(painter, point.x(), point.y(), note.text, size,
                              scale, note.color, anchor=note.anchor,
                              bold=note.bold)

    def _draw_overlay(self, painter, panel, box, x_range, y_range, overlay,
                      scale):
        """One polyline: broken wherever a point is NaN, with its symbols
        and error bars when it has them."""
        xs = np.asarray(overlay.xs, dtype=float)
        ys = np.asarray(overlay.ys, dtype=float)
        if xs.size == 0 or xs.size != ys.size:
            return
        good = np.isfinite(xs) & np.isfinite(ys)
        points = [self._map_point(box, panel, x_range, y_range, x, y) if ok else None
                  for x, y, ok in zip(xs, ys, good)]
        pen = self._pen(overlay.color, overlay.width_pt, scale, overlay.style)
        painter.setPen(pen)
        painter.setBrush(Qt.NoBrush)
        if getattr(overlay, "line", True) and xs.size >= 2:
            path = QPainterPath()
            drawing = False
            for point in points:
                if point is None:
                    drawing = False
                elif drawing:
                    path.lineTo(point)
                else:
                    path.moveTo(point)
                    drawing = True
            if overlay.closed and all(p is not None for p in points):
                path.closeSubpath()
            painter.drawPath(path)
        cap = 0.5 * max(getattr(overlay, "marker_mm", 1.0), 0.6) * scale
        solid = self._pen(overlay.color, overlay.width_pt, scale, "solid")
        for errors, vertical in ((getattr(overlay, "y_errors", None), True),
                                 (getattr(overlay, "x_errors", None), False)):
            if not errors or len(errors) != xs.size:
                continue
            painter.setPen(solid)
            for x, y, e, ok in zip(xs, ys, errors, good):
                if not ok or not np.isfinite(e) or e <= 0:
                    continue
                if vertical:
                    a = self._map_point(box, panel, x_range, y_range, x, y - e)
                    b = self._map_point(box, panel, x_range, y_range, x, y + e)
                    ends = ((QPointF(a.x() - cap, a.y()), QPointF(a.x() + cap, a.y())),
                            (QPointF(b.x() - cap, b.y()), QPointF(b.x() + cap, b.y())))
                else:
                    a = self._map_point(box, panel, x_range, y_range, x - e, y)
                    b = self._map_point(box, panel, x_range, y_range, x + e, y)
                    ends = ((QPointF(a.x(), a.y() - cap), QPointF(a.x(), a.y() + cap)),
                            (QPointF(b.x(), b.y() - cap), QPointF(b.x(), b.y() + cap)))
                painter.drawLine(a, b)
                for p, q in ends:
                    painter.drawLine(p, q)
        marker = getattr(overlay, "marker", "")
        if marker:
            size = max(float(getattr(overlay, "marker_mm", 1.0)), 0.2) * scale
            painter.setPen(Qt.NoPen)
            painter.setBrush(QBrush(QColor(overlay.color)))
            for point in points:
                if point is None:
                    continue
                self._draw_marker(painter, point, marker, size)
            painter.setBrush(Qt.NoBrush)

    @staticmethod
    def _draw_marker(painter, point, marker, size):
        half = size / 2.0
        x, y = point.x(), point.y()
        if marker == "s":
            painter.drawRect(QRectF(x - half, y - half, size, size))
        elif marker in ("^", "v", "d"):
            if marker == "^":
                corners = [(x, y - half), (x + half, y + half), (x - half, y + half)]
            elif marker == "v":
                corners = [(x, y + half), (x + half, y - half), (x - half, y - half)]
            else:
                corners = [(x, y - half), (x + half, y), (x, y + half), (x - half, y)]
            painter.drawPolygon(QPolygonF([QPointF(*c) for c in corners]))
        else:
            painter.drawEllipse(point, half, half)

    def _draw_arrow(self, painter, panel, box, x_range, y_range, arrow, scale):
        start = self._map_point(box, panel, x_range, y_range, arrow.x0, arrow.y0)
        end = self._map_point(box, panel, x_range, y_range, arrow.x1, arrow.y1)
        painter.setPen(self._pen(arrow.color, arrow.width_pt, scale))
        vector = np.array([end.x() - start.x(), end.y() - start.y()])
        length = float(np.hypot(*vector))
        if length < 1e-6:
            return
        unit = vector / length
        head = arrow.head_mm * scale
        base = QPointF(end.x() - unit[0] * head, end.y() - unit[1] * head)
        painter.drawLine(start, base)
        normal = np.array([-unit[1], unit[0]]) * head * 0.38
        painter.setBrush(QBrush(QColor(arrow.color)))
        painter.setPen(Qt.NoPen)
        painter.drawPolygon(QPolygonF([
            end,
            QPointF(base.x() + normal[0], base.y() + normal[1]),
            QPointF(base.x() - normal[0], base.y() - normal[1])]))
        painter.setBrush(Qt.NoBrush)

    def _draw_scale_bar(self, painter, panel, box, scale, x_range):
        bar = panel.scale_bar
        style = self.figure.style
        x0, x1 = x_range
        span = abs(x1 - x0)
        length = bar.length if bar.length > 0 else nf.nice_bar_length(span)
        fraction = length / span if span > 0 else 0.2
        width = fraction * box.width()
        margin = 0.06 * min(box.width(), box.height())
        thickness = pt_to_mm(bar.thickness_pt) * scale
        left = (box.left() + margin if "left" in bar.position
                else box.center().x() - width / 2 if "center" in bar.position
                else box.right() - margin - width)
        size = style.annotation_pt or style.font_pt
        text_h = pt_to_mm(size) * scale
        top = (box.top() + margin + (text_h * 1.3 if bar.show_label else 0)
               if "top" in bar.position else box.bottom() - margin - thickness)
        painter.fillRect(QRectF(left, top, width, thickness), QColor(bar.color))
        if bar.show_label:
            unit = bar.unit or nf.unit_of(panel.label_x())
            text = f"{length:g} {unit}".strip()
            baseline = top - 0.35 * text_h if "top" not in bar.position else top - 0.35 * text_h
            self._draw_markup(painter, left + width / 2, baseline, text, size,
                              scale, bar.color, anchor="center")

    def _draw_inset(self, painter, panel, box, scale, levels):
        inset = panel.inset
        style = self.figure.style
        x_range = (min(inset.x0, inset.x1), max(inset.x0, inset.x1))
        y_range = (min(inset.y0, inset.y1), max(inset.y0, inset.y1))
        if x_range[1] <= x_range[0] or y_range[1] <= y_range[0]:
            return
        short = min(box.width(), box.height())
        width = max(inset.size * short, 8.0)
        height = width * (box.height() / box.width()) \
            if box.width() > 0 else width
        margin = 0.04 * short
        left = (box.left() + margin if "left" in inset.position
                else box.right() - margin - width)
        top = (box.top() + margin if "top" in inset.position
               else box.bottom() - margin - height)
        target = QRectF(left, top, width, height)
        if inset.indicate:
            painter.setPen(self._pen(style.axis_color, style.axis_width_pt * 0.9,
                                     scale, "dashed"))
            painter.setBrush(Qt.NoBrush)
            corner0 = self._map_point(box, panel, panel.ranges()[0], panel.ranges()[1],
                                      x_range[0], y_range[0])
            corner1 = self._map_point(box, panel, panel.ranges()[0], panel.ranges()[1],
                                      x_range[1], y_range[1])
            painter.drawRect(QRectF(corner0, corner1).normalized())
        self._draw_image(painter, panel, target, x_range, y_range, levels, scale)
        if inset.frame:
            painter.setPen(self._pen(style.axis_color, style.axis_width_pt, scale))
            painter.setBrush(Qt.NoBrush)
            painter.drawRect(target)

    def _draw_colorbar(self, painter, panel, rect, box, scale, levels):
        style = self.figure.style
        if getattr(panel, "colorbar_position", "right") == "below":
            return self._draw_colorbar_below(painter, panel, box, scale, levels)
        tick_pt = style.tick_pt or style.font_pt
        left = box.right() + style.colorbar_gap_mm * scale
        width = style.colorbar_mm * scale
        fraction = float(np.clip(getattr(panel, "colorbar_length", 1.0) or 1.0, 0.1, 1.0))
        strip = QRectF(left, box.bottom() - fraction * box.height(), width,
                       fraction * box.height())
        try:
            lut = np.asarray(colormaps.get_lut(panel.colormap, flip=panel.flip))
        except Exception:
            lut = np.asarray(colormaps.get_lut("gray"))
        ramp = np.linspace(0.0, 1.0, 256)
        if panel.gamma and panel.gamma != 1.0:
            shown = ramp ** panel.gamma
        else:
            shown = ramp
        index = np.rint(shown * (lut.shape[0] - 1)).astype(np.intp)
        column = np.empty((256, 1, 4), dtype=np.uint8)
        column[..., :3] = lut[index][:, None, :3][::-1]
        column[..., 3] = 255
        buffer = np.ascontiguousarray(column)
        image = QImage(buffer.data, 1, 256, 4, QImage.Format_RGBA8888)
        painter.drawImage(strip, image)
        painter.setPen(self._pen(style.axis_color, style.axis_width_pt, scale))
        painter.setBrush(Qt.NoBrush)
        painter.drawRect(strip)

        lo, hi = levels
        metrics = QFontMetricsF(self._font(tick_pt, scale))
        gap = 0.3 * pt_to_mm(tick_pt) * scale
        widest = 0.0
        ends = getattr(panel, "colorbar_ends", "values")
        for value, y in ((hi, strip.top()), (lo, strip.bottom())):
            if ends == "none":
                continue
            text = (("High" if value == hi else "Low") if ends == "high-low"
                    else f"{value:.4g}")
            widest = max(widest, self._markup_width(painter, text, tick_pt, scale))
            baseline = (y + metrics.ascent() if value == hi and y == strip.top()
                        else y)
            self._draw_markup(painter, strip.right() + gap, baseline, text,
                              tick_pt, scale, style.axis_color)
        if panel.colorbar_label:
            self._draw_markup(
                painter, strip.right() + gap + widest + 0.9 * metrics.height(),
                strip.center().y(), panel.colorbar_label,
                style.label_pt or style.font_pt, scale, style.axis_color,
                anchor="center", rotate=-90.0)

    def _colorbar_column(self, panel, length: int = 256) -> QImage:
        """The colour map as an image ``length`` px long, low to high."""
        try:
            lut = np.asarray(colormaps.get_lut(panel.colormap, flip=panel.flip))
        except Exception:
            lut = np.asarray(colormaps.get_lut("gray"))
        ramp = np.linspace(0.0, 1.0, length)
        if panel.gamma and panel.gamma != 1.0:
            ramp = ramp ** panel.gamma
        index = np.rint(ramp * (lut.shape[0] - 1)).astype(np.intp)
        row = np.empty((1, length, 4), dtype=np.uint8)
        row[..., :3] = lut[index][None, :, :3]
        row[..., 3] = 255
        buffer = np.ascontiguousarray(row)
        return QImage(buffer.data, length, 1, 4 * length,
                      QImage.Format_RGBA8888).copy()

    def _draw_colorbar_below(self, painter, panel, box, scale, levels):
        """A lying colour bar under the panel -- under its tick labels and
        axis title when it has them -- centred, ``colorbar_length`` of the
        panel's width, its ends (or its levels) at either side."""
        style = self.figure.style
        tick_pt = style.tick_pt or style.font_pt
        label_pt = style.label_pt or style.font_pt
        below = 0.8 * scale
        if panel.data.is_plot() and panel.x_labels:
            below += (1.3 * pt_to_mm(tick_pt) + 1.6 * pt_to_mm(label_pt)) * scale
        fraction = float(np.clip(getattr(panel, "colorbar_length", 1.0) or 1.0, 0.1, 1.0))
        width = fraction * box.width()
        strip = QRectF(box.center().x() - width / 2, box.bottom() + below, width,
                       style.colorbar_mm * scale)
        painter.drawImage(strip, self._colorbar_column(panel))
        painter.setPen(self._pen(style.axis_color, style.axis_width_pt, scale))
        painter.setBrush(Qt.NoBrush)
        painter.drawRect(strip)
        ends = getattr(panel, "colorbar_ends", "values")
        if ends != "none":
            lo, hi = levels
            low, high = (("Low", "High") if ends == "high-low"
                         else (f"{lo:.4g}", f"{hi:.4g}"))
            gap = 0.6 * pt_to_mm(tick_pt) * scale
            baseline = strip.center().y() + 0.35 * pt_to_mm(tick_pt) * scale
            self._draw_markup(painter, strip.left() - gap, baseline, low, tick_pt,
                              scale, style.axis_color, anchor="right")
            self._draw_markup(painter, strip.right() + gap, baseline, high, tick_pt,
                              scale, style.axis_color, anchor="left")
        if panel.colorbar_label:
            metrics = QFontMetricsF(self._font(label_pt, scale))
            self._draw_markup(painter, strip.center().x(),
                              strip.bottom() + 1.1 * metrics.ascent(),
                              panel.colorbar_label, label_pt, scale,
                              style.axis_color, anchor="center")

    def _draw_panel_label(self, painter, panel, box, scale):
        style = self.figure.style
        if not panel.label:
            return
        size = style.panel_label_pt or (style.font_pt + 1)
        margin = 0.05 * min(box.width(), box.height())
        text_h = pt_to_mm(size) * scale
        inside = style.panel_label_inside
        position = style.panel_label_position
        if "left" in position:
            x, anchor = (box.left() + margin if inside else box.left()), "left"
        elif "right" in position:
            x, anchor = (box.right() - margin if inside else box.right()), "right"
        else:
            x, anchor = box.center().x(), "center"
        if "top" in position:
            # Outside, the letter has to clear whatever else is above the
            # box -- a title, or a second axis with its own numbers and
            # label -- or it lands on top of them.
            above = 0.0
            if not inside:
                label_mm = pt_to_mm(style.label_pt or style.font_pt)
                tick_mm = pt_to_mm(style.tick_pt or style.font_pt)
                if panel.title:
                    above += 1.4 * label_mm * scale
                if panel.top_axis is not None and panel.top_axis.enabled:
                    above += (1.1 * tick_mm + 1.35 * label_mm) * scale
            y = box.top() + margin + text_h if inside else box.top() - margin - above
        else:
            y = box.bottom() - margin if inside else box.bottom() + margin + text_h
        aligned = getattr(self, "_label_at", {}).get(id(panel))
        if aligned is not None and not inside:
            x = aligned[0] if aligned[0] is not None else x
            y = aligned[1] if aligned[1] is not None else y
        dx, dy = getattr(panel, "label_offset", (0.0, 0.0)) or (0.0, 0.0)
        x, y = x + float(dx) * scale, y + float(dy) * scale
        color = style.panel_label_color if inside else style.axis_color
        self._draw_markup(painter, x, y, panel.label, size, scale, color,
                          anchor=anchor, bold=style.panel_label_bold)

    def _natural_label_point(self, panel, box, scale):
        """Where :meth:`_draw_panel_label` puts a letter on its own (outside,
        top left): ``(x, y)`` in device units."""
        style = self.figure.style
        margin = 0.05 * min(box.width(), box.height())
        label_mm = pt_to_mm(style.label_pt or style.font_pt)
        tick_mm = pt_to_mm(style.tick_pt or style.font_pt)
        above = 0.0
        if panel.title:
            above += 1.4 * label_mm * scale
        if panel.top_axis is not None and panel.top_axis.enabled:
            above += (1.1 * tick_mm + 1.35 * label_mm) * scale
        return box.left(), box.top() - margin - above

    def _align_labels(self, rects, scale):
        """The letters of a row of panels at one height (the highest any of
        them needs), and of a column at one left edge -- so (a), (b), (c)
        line up however the panels' titles and axes differ."""
        self._label_at = {}
        style = self.figure.style
        if (not getattr(style, "panel_label_aligned", True) or style.panel_label_inside
                or "top" not in style.panel_label_position
                or "left" not in style.panel_label_position):
            return
        boxes, natural = {}, {}
        for rect in rects:
            panel = self.figure.panels[rect["index"]]
            if not panel.label or panel.role == "inset":
                continue
            boxes[id(panel)] = (rect["x"], rect["y"], rect["w"], rect["h"])
            box = QRectF(rect["x"] * scale, rect["y"] * scale,
                         rect["w"] * scale, rect["h"] * scale)
            natural[id(panel)] = self._natural_label_point(panel, box, scale)
        rows, columns = Figure.label_groups(boxes)
        point = {key: [None, None] for key in boxes}
        for row in rows:
            y = min(natural[key][1] for key in row)
            for key in row:
                point[key][1] = y
        for column in columns:
            x = min(natural[key][0] for key in column)
            for key in column:
                point[key][0] = x
        self._label_at = {key: tuple(value) for key, value in point.items()}

    # -- outputs -------------------------------------------------------------
    def to_image(self, dpi: float = None, selection: int = None) -> QImage:
        dpi = float(dpi or self.figure.style.dpi)
        width, height = self.pixel_size(dpi)
        image = QImage(width, height, QImage.Format_RGBA8888)
        image.fill(QColor(self.figure.style.background))
        painter = QPainter(image)
        try:
            self.paint(painter, dpi / nf.MM_PER_INCH, selection)
        finally:
            painter.end()
        return image

    def preview(self, max_width_px: int, selection: int = None) -> QImage:
        """A screen-sized render. The data raster is generated at the preview
        scale rather than the figure's dpi, which is the difference between
        a preview that redraws while you drag a spin box and one that does
        not."""
        width_mm, _height = self.size_mm()
        scale = max(min(max_width_px / max(width_mm, 1.0), 20.0), 1.0)
        self._raster_scale = scale
        try:
            width = max(int(round(width_mm * scale)), 1)
            height = max(int(round(self.size_mm()[1] * scale)), 1)
            image = QImage(width, height, QImage.Format_RGBA8888)
            image.fill(QColor(self.figure.style.background))
            painter = QPainter(image)
            try:
                self.paint(painter, scale, selection)
            finally:
                painter.end()
        finally:
            self._raster_scale = None
        return image

    def to_pdf(self, path: str):
        """A vector PDF: text is text, the axes are paths, the data is an
        image embedded at the figure's dpi."""
        width_mm, height_mm = self.size_mm()
        writer = QPdfWriter(path)
        writer.setPageSize(QPageSize(
            QSizeF_mm(width_mm, height_mm), QPageSize.Millimeter,
            "figure", QPageSize.ExactMatch))
        writer.setPageMargins(QMarginsF_zero())
        writer.setResolution(1200)
        painter = QPainter(writer)
        try:
            self.paint(painter, writer.resolution() / nf.MM_PER_INCH, None)
        finally:
            painter.end()
        return path

    def to_svg(self, path: str):
        from PyQt5.QtSvg import QSvgGenerator
        width_mm, height_mm = self.size_mm()
        scale = 72.0 / nf.MM_PER_INCH          # points
        generator = QSvgGenerator()
        generator.setFileName(path)
        generator.setSize(QSize(int(round(width_mm * scale)),
                                int(round(height_mm * scale))))
        generator.setViewBox(QRectF(0, 0, width_mm * scale, height_mm * scale))
        generator.setResolution(72)
        painter = QPainter(generator)
        try:
            self.paint(painter, scale, None)
        finally:
            painter.end()
        return path


def QSizeF_mm(width, height):
    from PyQt5.QtCore import QSizeF
    return QSizeF(float(width), float(height))


def QMarginsF_zero():
    from PyQt5.QtCore import QMarginsF
    return QMarginsF(0.0, 0.0, 0.0, 0.0)


# --------------------------------------------------------------------------
# The canvas
# --------------------------------------------------------------------------
class FigureCanvas(QWidget):
    """The page itself, edited with the mouse.

    The figure is drawn by :class:`FigurePainter`, the same call that
    writes the file, so the page on screen is the page in the PDF. On top
    of that picture this widget draws only what is not part of the figure:
    the selection's handles, the snapping guides, and the shape being
    drawn.

    The tools (``tool``):

    * ``select`` -- click a panel, an inset or a drawing to select it; drag
      it to move it, drag a handle to resize it. Panels snap to each
      other's edges, and to a 0.5 mm grid.
    * ``slot`` -- drag out a new empty panel (inside another panel: an
      inset).
    * ``zoom`` -- drag out a region of a panel: an inset magnifying it.
    * ``arrow``, ``line``, ``rect``, ``ellipse`` -- drag to draw one.
    * ``text`` -- click where a note goes.

    The page is laid out in millimetres; this widget only converts.
    """

    panelPicked = pyqtSignal(int)
    pointPicked = pyqtSignal(int, float, float)
    #: the selection changed (a panel, a drawing, or nothing)
    selectionChanged = pyqtSignal()
    #: about to change the figure: the window keeps an undo step
    aboutToChange = pyqtSignal()
    #: the figure was changed by the mouse or the keyboard
    changed = pyqtSignal()
    #: a drawing tool finished, and the canvas is back to ``select``
    toolFinished = pyqtSignal()
    #: a panel was double-clicked
    panelActivated = pyqtSignal(int)
    #: right click: (what, index, global position), what in panel/shape/page
    menuRequested = pyqtSignal(str, int, object)

    TOOLS = ("select", "slot", "zoom", "arrow", "line", "rect", "ellipse", "text")
    #: handles, by the edges of the frame each one moves
    HANDLES = {"nw": "lt", "n": "t", "ne": "rt", "e": "r",
               "se": "rb", "s": "b", "sw": "lb", "w": "l"}
    MARGIN_PX = 14

    def __init__(self, figure: Figure, parent=None):
        super().__init__(parent)
        self.figure = figure
        self.painter = FigurePainter(figure)
        self.painter.placeholders = True
        self.selection = None          # ("panel", index) / ("shape", index)
        #: more panels picked with Ctrl / Shift + click, for arranging
        self.extra = []
        self.tool = "select"
        self.snap = True
        #: None fits the page to the window; otherwise screen px per mm
        self.zoom = None
        self._image = None
        self._scale = 1.0
        self._origin = (0.0, 0.0)
        self._rect = None              # where the page sits: x, y, w, h (px)
        self._picking = None           # a callback while a point is picked
        self._drag = None              # what the mouse is doing
        self._guides = []              # snapping guides: ("v"/"h", mm)
        self._draft = None             # the shape or box being drawn
        self.setMinimumSize(280, 220)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.setFocusPolicy(Qt.StrongFocus)
        self.setMouseTracking(True)

    # -- the selection, as the tools see it ---------------------------------
    @property
    def selected(self) -> int:
        """The selected panel's index (the first panel when a drawing or
        nothing is selected), for the tools that work on one panel."""
        if self.selection and self.selection[0] == "panel":
            return self.selection[1]
        return 0

    @selected.setter
    def selected(self, index):
        if 0 <= int(index) < len(self.figure.panels):
            self.selection = ("panel", int(index))
        self.update()

    def selected_shape(self):
        if self.selection and self.selection[0] == "shape" \
                and 0 <= self.selection[1] < len(self.figure.shapes):
            return self.figure.shapes[self.selection[1]]
        return None

    def selected_panel(self):
        if self.selection and self.selection[0] == "panel" \
                and 0 <= self.selection[1] < len(self.figure.panels):
            return self.figure.panels[self.selection[1]]
        return None

    def selected_indices(self):
        """Every selected panel: the one clicked first, then those added
        with Ctrl or Shift."""
        first = [self.selection[1]] if self.selection and self.selection[0] == "panel" else []
        return first + [i for i in self.extra if i not in first
                        and 0 <= i < len(self.figure.panels)]

    def select(self, what, index=None):
        self.extra = []
        self.selection = None if what is None else (what, int(index))
        self.update()
        self.selectionChanged.emit()
        if what == "panel":
            self.panelPicked.emit(int(index))

    def set_tool(self, tool: str):
        if tool not in self.TOOLS:
            raise ValueError(f"unknown tool {tool!r}")
        self.tool = tool
        self._drag = self._draft = None
        self.setCursor(Qt.ArrowCursor if tool == "select" else Qt.CrossCursor)
        self.update()

    # -- rendering ---------------------------------------------------------
    def refresh(self):
        self._image = None
        self.update()

    def _page_scale(self):
        width_mm, height_mm = self.painter.size_mm()
        if self.zoom:
            return float(self.zoom)
        available_w = max(self.width() - 2 * self.MARGIN_PX, 60)
        available_h = max(self.height() - 2 * self.MARGIN_PX, 60)
        return min(available_w / max(width_mm, 1.0),
                   available_h / max(height_mm, 1.0))

    def _ensure_image(self):
        if self._image is not None:
            return
        width_mm, height_mm = self.painter.size_mm()
        scale = self._page_scale()
        if self.zoom:
            self.setMinimumSize(int(width_mm * scale) + 2 * self.MARGIN_PX,
                                int(height_mm * scale) + 2 * self.MARGIN_PX)
        else:
            self.setMinimumSize(280, 220)
        self._image = self.painter.preview(max(int(width_mm * scale), 60), None)
        self._scale = self._image.width() / max(width_mm, 1.0)
        x = max((self.width() - self._image.width()) // 2, self.MARGIN_PX)
        y = max((self.height() - self._image.height()) // 2, self.MARGIN_PX)
        self._origin = (float(x), float(y))
        self._rect = (x, y, self._image.width(), self._image.height())

    def to_mm(self, pos):
        return ((pos.x() - self._origin[0]) / self._scale,
                (pos.y() - self._origin[1]) / self._scale)

    def to_px(self, x_mm, y_mm):
        return QPointF(self._origin[0] + x_mm * self._scale,
                       self._origin[1] + y_mm * self._scale)

    def paintEvent(self, event):
        self._ensure_image()
        painter = QPainter(self)
        painter.fillRect(self.rect(), QColor("#e9e6ef"))
        if self._image is not None:
            x, y, w, h = self._rect
            painter.fillRect(QRectF(x + 3, y + 3, w, h), QColor(0, 0, 0, 40))
            painter.drawImage(x, y, self._image)
        painter.setRenderHint(QPainter.Antialiasing, True)
        self._paint_draft(painter)
        self._paint_selection(painter)
        self._paint_guides(painter)
        painter.end()

    def _paint_draft(self, painter):
        draft = self._draft
        if draft is None:
            return
        if isinstance(draft, nf.Shape):
            painter.save()
            painter.translate(*self._origin)
            self.painter.paint_shape(painter, draft, self._scale)
            painter.restore()
            return
        x0, y0, x1, y1 = draft
        pen = QPen(QColor("#2ca6cd"))
        pen.setWidthF(1.5)
        pen.setStyle(Qt.DashLine)
        painter.setPen(pen)
        painter.setBrush(QBrush(QColor(44, 166, 205, 40)))
        painter.drawRect(QRectF(self.to_px(x0, y0), self.to_px(x1, y1)).normalized())

    def _paint_selection(self, painter):
        if self.selection is None:
            return
        pen = QPen(QColor("#2ca6cd"))
        pen.setWidthF(1.2)
        pen.setStyle(Qt.DashLine)
        painter.setPen(pen)
        painter.setBrush(Qt.NoBrush)
        shape = self.selected_shape()
        if shape is not None and shape.kind == "text":
            painter.drawRect(self.painter.text_box(shape, self._scale)
                             .translated(*self._origin).adjusted(-2, -2, 2, 2))
        elif shape is None:
            for index in self.selected_indices():
                frame = self._frame_of(self.figure.panels[index])
                if frame is not None:
                    x, y, w, h = frame
                    painter.drawRect(QRectF(self.to_px(x, y), self.to_px(x + w, y + h)))
        painter.setPen(QPen(QColor("#1d7fa0"), 1.0))
        painter.setBrush(QBrush(QColor("#ffffff")))
        for _name, point in self._handles():
            painter.drawRect(QRectF(point.x() - 4, point.y() - 4, 8, 8))

    def _paint_guides(self, painter):
        if not self._guides or self._rect is None:
            return
        pen = QPen(QColor("#e0389b"))
        pen.setWidthF(1.0)
        painter.setPen(pen)
        x, y, w, h = self._rect
        for direction, value in self._guides:
            point = self.to_px(value, value)
            if direction == "v":
                painter.drawLine(QPointF(point.x(), y), QPointF(point.x(), y + h))
            else:
                painter.drawLine(QPointF(x, point.y()), QPointF(x + w, point.y()))

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self._image = None

    # -- geometry ------------------------------------------------------------
    def _frame_of(self, panel):
        """A panel's data box in mm, whether the page is free or a grid."""
        if panel is None:
            return None
        if panel.frame is not None:
            return panel.frame
        _w, _h, rects = self.figure.layout()
        index = self.figure.panels.index(panel)
        for rect in rects:
            if rect["index"] == index:
                return (rect["x"], rect["y"], rect["w"], rect["h"])
        return None

    def _handles(self):
        """``[(name, point in px)]`` for the selection."""
        shape = self.selected_shape()
        if shape is not None:
            if shape.kind in ("arrow", "line"):
                return [("p0", self.to_px(shape.x0, shape.y0)),
                        ("p1", self.to_px(shape.x1, shape.y1))]
            if shape.kind in ("rect", "ellipse"):
                left, top, right, bottom = shape.bounds()
                return [("nw", self.to_px(left, top)), ("ne", self.to_px(right, top)),
                        ("se", self.to_px(right, bottom)), ("sw", self.to_px(left, bottom))]
            return []
        frame = self._frame_of(self.selected_panel())
        if frame is None:
            return []
        x, y, w, h = frame
        points = {"nw": (x, y), "n": (x + w / 2, y), "ne": (x + w, y),
                  "e": (x + w, y + h / 2), "se": (x + w, y + h),
                  "s": (x + w / 2, y + h), "sw": (x, y + h), "w": (x, y + h / 2)}
        return [(name, self.to_px(*xy)) for name, xy in points.items()]

    def _tolerance_mm(self):
        return 5.0 / max(self._scale, 1e-6)

    def _hit_shape(self, shape, x, y, tol) -> bool:
        if shape.kind == "text":
            box = self.painter.text_box(shape, self._scale)
            return QRectF(box.x() / self._scale - tol, box.y() / self._scale - tol,
                          box.width() / self._scale + 2 * tol,
                          box.height() / self._scale + 2 * tol).contains(QPointF(x, y))
        if shape.kind in ("arrow", "line"):
            return _segment_distance(x, y, shape.x0, shape.y0, shape.x1, shape.y1) <= tol
        left, top, right, bottom = shape.bounds()
        inside = left - tol <= x <= right + tol and top - tol <= y <= bottom + tol
        if not inside:
            return False
        if shape.fill:
            return True
        if shape.kind == "rect":
            return (min(abs(x - left), abs(x - right)) <= tol
                    or min(abs(y - top), abs(y - bottom)) <= tol)
        cx, cy = (left + right) / 2, (top + bottom) / 2
        rx, ry = max((right - left) / 2, 1e-6), max((bottom - top) / 2, 1e-6)
        radius = np.hypot((x - cx) / rx, (y - cy) / ry)
        return abs(radius - 1.0) * min(rx, ry) <= tol

    def item_at(self, pos):
        """``("handle", name)``, ``("shape", i)``, ``("panel", i)`` or None
        under a widget position: the selection's handles first, then the
        drawings, the insets and the panels, topmost first."""
        for name, point in self._handles():
            if abs(point.x() - pos.x()) <= 6 and abs(point.y() - pos.y()) <= 6:
                return ("handle", name)
        x, y = self.to_mm(pos)
        tol = self._tolerance_mm()
        for index in range(len(self.figure.shapes) - 1, -1, -1):
            if self._hit_shape(self.figure.shapes[index], x, y, tol):
                return ("shape", index)
        panels = self.figure.panels
        order = ([i for i, p in enumerate(panels) if p.role == "inset"][::-1]
                 + [i for i, p in enumerate(panels) if p.role != "inset"][::-1])
        for index in order:
            frame = self._frame_of(panels[index])
            if frame is None:
                continue
            fx, fy, fw, fh = frame
            if fx - tol <= x <= fx + fw + tol and fy - tol <= y <= fy + fh + tol:
                return ("panel", index)
        return None

    def _locate(self, pos):
        """Which panel was clicked, and where in its data."""
        if self._rect is None:
            return None
        mm_x, mm_y = self.to_mm(pos)
        _w, _h, rects = self.figure.layout()
        for rect in reversed(rects):
            if (rect["x"] <= mm_x <= rect["x"] + rect["w"]
                    and rect["y"] <= mm_y <= rect["y"] + rect["h"]):
                panel = self.figure.panels[rect["index"]]
                (x0, x1), (y0, y1) = self.figure.panel_ranges(panel)
                fx = (mm_x - rect["x"]) / rect["w"]
                fy = 1.0 - (mm_y - rect["y"]) / rect["h"]
                if panel.invert_x:
                    fx = 1.0 - fx
                if panel.invert_y:
                    fy = 1.0 - fy
                return rect["index"], x0 + fx * (x1 - x0), y0 + fy * (y1 - y0)
        return None

    # -- snapping ------------------------------------------------------------
    def _snap_value(self, value):
        return round(value * 2.0) / 2.0 if self.snap else value

    def _snap_frame(self, moving, frame, edges):
        """Snap ``frame``'s ``edges`` (some of "l", "r", "t", "b") to the
        other panels' edges within a few pixels, else to the 0.5 mm grid.
        Returns the frame, with the guides that apply."""
        x, y, w, h = frame
        self._guides = []
        if not self.snap:
            return frame
        threshold = max(6.0 / self._scale, 0.6)
        xs, ys = [0.0, self.figure.style.width_mm], [0.0]
        for panel in self.figure.panels:
            if panel is moving or panel.frame is None:
                continue
            if panel.role == "inset" and moving.role != "inset":
                continue          # they move with their panel
            px, py, pw, ph = panel.frame
            xs += [px, px + pw]
            ys += [py, py + ph]

        def best(value, candidates):
            near = [(abs(value - c), c) for c in candidates if abs(value - c) <= threshold]
            return min(near)[1] if near else None

        moving_whole = set(edges) == {"l", "r", "t", "b"}
        # x
        if moving_whole:
            left, right = best(x, xs), best(x + w, xs)
            if left is not None and (right is None or abs(left - x) <= abs(right - x - w)):
                x = left
                self._guides.append(("v", left))
            elif right is not None:
                x = right - w
                self._guides.append(("v", right))
            else:
                x = self._snap_value(x)
            top, bottom = best(y, ys), best(y + h, ys)
            if top is not None and (bottom is None or abs(top - y) <= abs(bottom - y - h)):
                y = top
                self._guides.append(("h", top))
            elif bottom is not None:
                y = bottom - h
                self._guides.append(("h", bottom))
            else:
                y = self._snap_value(y)
            return (x, y, w, h)
        left, top, right, bottom = x, y, x + w, y + h
        for edge in edges:
            if edge == "l":
                found = best(left, xs)
                left = found if found is not None else self._snap_value(left)
                if found is not None:
                    self._guides.append(("v", found))
            elif edge == "r":
                found = best(right, xs)
                right = found if found is not None else self._snap_value(right)
                if found is not None:
                    self._guides.append(("v", found))
            elif edge == "t":
                found = best(top, ys)
                top = found if found is not None else self._snap_value(top)
                if found is not None:
                    self._guides.append(("h", found))
            elif edge == "b":
                found = best(bottom, ys)
                bottom = found if found is not None else self._snap_value(bottom)
                if found is not None:
                    self._guides.append(("h", found))
        return (left, top, max(right - left, 2.0), max(bottom - top, 2.0))

    # -- the mouse -------------------------------------------------------------
    def start_point_picking(self, callback):
        """The next click reports ``(panel index, x, y)`` in data
        coordinates -- how an annotation gets placed where it is wanted
        rather than typed in blind."""
        self._picking = callback
        self.setCursor(Qt.CrossCursor)

    def stop_point_picking(self):
        self._picking = None
        self.setCursor(Qt.ArrowCursor if self.tool == "select" else Qt.CrossCursor)

    def _begin_change(self):
        """Before the first change of a gesture: an undo step, and a page
        laid out by hand from here on."""
        self.aboutToChange.emit()
        self.figure.freeze()

    def mousePressEvent(self, event):
        self.setFocus()
        self._ensure_image()
        if event.button() == Qt.RightButton:
            found = self.item_at(event.pos())
            if found is not None and found[0] in ("panel", "shape"):
                self.select(*found)
                self.menuRequested.emit(found[0], found[1], event.globalPos())
            else:
                self.menuRequested.emit("page", -1, event.globalPos())
            return
        if event.button() != Qt.LeftButton:
            return
        if self._picking is not None:
            found = self._locate(event.pos())
            if found is None:
                return
            callback, self._picking = self._picking, None
            self.stop_point_picking()
            callback(*found)
            return
        x, y = self.to_mm(event.pos())
        if self.tool == "select":
            found = self.item_at(event.pos())
            if found is None:
                self.select(None)
                return
            adding = event.modifiers() & (Qt.ControlModifier | Qt.ShiftModifier)
            if (adding and found[0] == "panel" and self.selection
                    and self.selection[0] == "panel"):
                # Ctrl / Shift: add the panel to the selection (or take it out)
                index = int(found[1])
                if index in self.extra:
                    self.extra.remove(index)
                elif index != self.selection[1]:
                    self.extra.append(index)
                self.update()
                self.selectionChanged.emit()
                return
            if found[0] != "handle":
                self.select(*found)
            self._drag = {"kind": found[0], "name": found[1] if found[0] == "handle" else None,
                          "start": (x, y), "started": False,
                          "origin": self._origin_geometry()}
            return
        if self.tool == "text":
            self._place_text(self._snap_value(x), self._snap_value(y))
            return
        sx, sy = self._snap_value(x), self._snap_value(y)
        self._drag = {"kind": "draw", "start": (sx, sy), "started": True}
        if self.tool in ("arrow", "line", "rect", "ellipse"):
            self._draft = nf.Shape(kind=self.tool, x0=sx, y0=sy, x1=sx, y1=sy,
                                   color=self.default_colour())
        else:
            self._draft = (sx, sy, sx, sy)
        self.update()

    def default_colour(self) -> str:
        """Drawings start in the figure's axis colour."""
        return self.figure.style.axis_color

    def _origin_geometry(self):
        shape = self.selected_shape()
        if shape is not None:
            return (shape.x0, shape.y0, shape.x1, shape.y1)
        frame = self._frame_of(self.selected_panel())
        return frame

    def mouseMoveEvent(self, event):
        drag = self._drag
        x, y = self.to_mm(event.pos())
        if drag is None:
            if self.tool == "select":
                found = self.item_at(event.pos())
                if found and found[0] == "handle":
                    self.setCursor(_HANDLE_CURSORS.get(found[1], Qt.SizeAllCursor))
                elif found:
                    self.setCursor(Qt.SizeAllCursor)
                else:
                    self.setCursor(Qt.ArrowCursor)
            return
        dx, dy = x - drag["start"][0], y - drag["start"][1]
        if drag["kind"] == "draw":
            sx, sy = drag["start"]
            ex, ey = self._snap_value(x), self._snap_value(y)
            if isinstance(self._draft, nf.Shape):
                if event.modifiers() & Qt.ShiftModifier and self._draft.kind in ("arrow", "line"):
                    # Shift: horizontal, vertical or 45°
                    angle = np.arctan2(ey - sy, ex - sx)
                    step = np.pi / 4
                    angle = round(angle / step) * step
                    length = np.hypot(ex - sx, ey - sy)
                    ex, ey = sx + length * np.cos(angle), sy + length * np.sin(angle)
                self._draft.x1, self._draft.y1 = ex, ey
            else:
                self._draft = (sx, sy, ex, ey)
            self.update()
            return
        if not drag["started"]:
            if abs(dx) * self._scale < 3 and abs(dy) * self._scale < 3:
                return
            drag["started"] = True
            self._begin_change()
            drag["origin"] = self._origin_geometry()
        origin = drag["origin"]
        if origin is None:
            return
        shape = self.selected_shape()
        if shape is not None:
            x0, y0, x1, y1 = origin
            name = drag["name"]
            if drag["kind"] == "shape":
                sdx, sdy = self._snap_value(x0 + dx) - x0, self._snap_value(y0 + dy) - y0
                shape.x0, shape.y0, shape.x1, shape.y1 = x0 + sdx, y0 + sdy, x1 + sdx, y1 + sdy
            elif name == "p0":
                shape.x0, shape.y0 = self._snap_value(x), self._snap_value(y)
            elif name == "p1":
                shape.x1, shape.y1 = self._snap_value(x), self._snap_value(y)
            else:
                left, top, right, bottom = min(x0, x1), min(y0, y1), max(x0, x1), max(y0, y1)
                sx, sy = self._snap_value(x), self._snap_value(y)
                if "w" in name:
                    left = sx
                if "e" in name:
                    right = sx
                if "n" in name:
                    top = sy
                if "s" in name:
                    bottom = sy
                shape.x0, shape.y0, shape.x1, shape.y1 = left, top, right, bottom
        else:
            panel = self.selected_panel()
            if panel is None:
                return
            fx, fy, fw, fh = origin
            if drag["kind"] == "panel":
                frame = self._snap_frame(panel, (fx + dx, fy + dy, fw, fh),
                                         ("l", "r", "t", "b"))
            else:
                edges = self.HANDLES[drag["name"]]
                left, top, right, bottom = fx, fy, fx + fw, fy + fh
                if "l" in edges:
                    left = min(fx + dx, right - 2.0)
                if "r" in edges:
                    right = max(fx + fw + dx, left + 2.0)
                if "t" in edges:
                    top = min(fy + dy, bottom - 2.0)
                if "b" in edges:
                    bottom = max(fy + fh + dy, top + 2.0)
                if (event.modifiers() & Qt.ShiftModifier) and len(edges) == 2:
                    # Shift on a corner: keep the panel's shape
                    ratio = fh / fw if fw else 1.0
                    width = right - left
                    if "t" in edges:
                        top = bottom - width * ratio
                    else:
                        bottom = top + width * ratio
                frame = self._snap_frame(panel, (left, top, right - left, bottom - top),
                                         tuple(edges))
            self.figure.move_panel(panel, frame)
        self.refresh()
        self.changed.emit()

    def mouseReleaseEvent(self, event):
        drag, self._drag = self._drag, None
        self._guides = []
        if drag is None:
            return
        if drag["kind"] == "draw":
            draft, self._draft = self._draft, None
            self._finish_drawing(draft)
            return
        if drag.get("started"):
            if drag["kind"] in ("panel", "handle") and self.selected_shape() is None:
                self.figure.auto_edge_labels()
            self.refresh()
            self.changed.emit()
        self.update()

    def _finish_drawing(self, draft):
        tool = self.tool
        if isinstance(draft, nf.Shape):
            if np.hypot(draft.x1 - draft.x0, draft.y1 - draft.y0) < 1.0:
                self.update()
                return
            self._begin_change()
            self.figure.shapes.append(draft)
            self._done(("shape", len(self.figure.shapes) - 1))
            return
        if draft is None:
            return
        x0, y0, x1, y1 = draft
        left, top = min(x0, x1), min(y0, y1)
        width, height = abs(x1 - x0), abs(y1 - y0)
        if width < 2.0 or height < 2.0:
            self.update()
            return
        if tool == "slot":
            self._begin_change()
            panel = Panel(data=nf.empty_data(), frame=(left, top, width, height))
            host = self._panel_containing(left + width / 2, top + height / 2)
            if host is not None and host.frame is not None \
                    and width < host.frame[2] and height < host.frame[3]:
                panel.role = "inset"
                panel.x_labels = panel.y_labels = False
            self.figure.panels.append(panel)
            if panel.role != "inset":
                self.figure.letter_panels()
            self._done(("panel", len(self.figure.panels) - 1))
        elif tool == "zoom":
            self._make_zoom(left, top, width, height)

    def _panel_containing(self, x, y, plots_only=False):
        for panel in self.figure.main_panels():
            frame = self._frame_of(panel)
            if frame is None or (plots_only and not panel.data.is_plot()):
                continue
            fx, fy, fw, fh = frame
            if fx <= x <= fx + fw and fy <= y <= fy + fh:
                return panel
        return None

    def _make_zoom(self, left, top, width, height):
        """An inset magnifying the region dragged out on a panel, placed in
        the panel's top right corner."""
        host = self._panel_containing(left + width / 2, top + height / 2,
                                      plots_only=True)
        if host is None:
            self.update()
            return None
        self._begin_change()
        fx, fy, fw, fh = host.frame
        (x0, x1), (y0, y1) = self.figure.panel_ranges(host)

        def data_x(mm):
            frac = (mm - fx) / fw
            frac = 1.0 - frac if host.invert_x else frac
            return x0 + frac * (x1 - x0)

        def data_y(mm):
            frac = 1.0 - (mm - fy) / fh
            frac = 1.0 - frac if host.invert_y else frac
            return y0 + frac * (y1 - y0)

        xs = sorted((data_x(left), data_x(left + width)))
        ys = sorted((data_y(top), data_y(top + height)))
        inset = copy.copy(host)
        for name in ("lines", "texts", "arrows", "overlays", "x_tick_labels"):
            setattr(inset, name, copy.deepcopy(getattr(host, name)))
        inset.uid = uuid.uuid4().hex[:10]
        inset.role, inset.zoom_of, inset.label = "inset", host.uid, ""
        inset.x_range, inset.y_range = tuple(xs), tuple(ys)
        inset.colorbar, inset.scale_bar, inset.inset = False, None, None
        inset.title = ""
        inset.top_axis, inset.right_axis = nf.SecondaryAxis(), nf.SecondaryAxis()
        inset.x_labels = inset.y_labels = False
        # The region's own shape, about 40 % of the panel, in its corner.
        size = 0.42 * fw
        aspect = (height / width) if width else 1.0
        w_in, h_in = size, min(size * aspect, 0.6 * fh)
        if h_in < size * aspect:
            w_in = h_in / max(aspect, 1e-6)
        inset.frame = (fx + fw - w_in - 0.04 * fw, fy + 0.04 * fh, w_in, h_in)
        self.figure.panels.append(inset)
        self._done(("panel", len(self.figure.panels) - 1))
        return inset

    def _place_text(self, x, y):
        text, ok = QInputDialog.getText(self, "Text", "Text (k_x, E_F, Å^-1 ...):")
        if not ok or not text.strip():
            return None
        self._begin_change()
        self.figure.shapes.append(nf.Shape(kind="text", x0=x, y0=y, x1=x, y1=y,
                                           text=text, color=self.default_colour()))
        self._done(("shape", len(self.figure.shapes) - 1))
        return self.figure.shapes[-1]

    def _done(self, selection):
        self.selection = selection
        self.set_tool("select")
        self.refresh()
        self.changed.emit()
        self.selectionChanged.emit()
        self.toolFinished.emit()

    def mouseDoubleClickEvent(self, event):
        found = self.item_at(event.pos())
        if found is None or found[0] == "handle":
            return
        if found[0] == "shape":
            shape = self.figure.shapes[found[1]]
            if shape.kind == "text":
                text, ok = QInputDialog.getText(self, "Text", "Text:", text=shape.text)
                if ok:
                    self.aboutToChange.emit()
                    shape.text = text
                    self.refresh()
                    self.changed.emit()
            return
        self.panelActivated.emit(found[1])

    # -- the keyboard -----------------------------------------------------------
    def delete_selection(self):
        if self.selection is None:
            return False
        what, index = self.selection
        self.aboutToChange.emit()
        if what == "shape" and 0 <= index < len(self.figure.shapes):
            self.figure.shapes.pop(index)
        elif what == "panel" and 0 <= index < len(self.figure.panels):
            panel = self.figure.panels[index]
            if len(self.figure.main_panels()) <= 1 and panel.role != "inset":
                return False              # a figure keeps one panel
            self.figure.freeze()
            self.figure.panels.pop(index)
            # Insets go with the panel they magnify.
            self.figure.panels = [p for p in self.figure.panels
                                  if p.zoom_of != panel.uid]
            if panel.role != "inset":
                self.figure.letter_panels()
        self.selection = None
        self.refresh()
        self.changed.emit()
        self.selectionChanged.emit()
        return True

    def nudge(self, dx, dy):
        if self.selection is None:
            return
        self._begin_change()
        shape = self.selected_shape()
        if shape is not None:
            shape.move(dx, dy)
        else:
            panel = self.selected_panel()
            x, y, w, h = panel.frame
            self.figure.move_panel(panel, (x + dx, y + dy, w, h))
        self.refresh()
        self.changed.emit()

    def keyPressEvent(self, event):
        key = event.key()
        step = 5.0 if event.modifiers() & Qt.ShiftModifier else 0.5
        if key in (Qt.Key_Delete, Qt.Key_Backspace):
            self.delete_selection()
        elif key == Qt.Key_Escape:
            self._drag = self._draft = None
            self.set_tool("select")
            self.toolFinished.emit()
        elif key == Qt.Key_Left:
            self.nudge(-step, 0)
        elif key == Qt.Key_Right:
            self.nudge(step, 0)
        elif key == Qt.Key_Up:
            self.nudge(0, -step)
        elif key == Qt.Key_Down:
            self.nudge(0, step)
        else:
            super().keyPressEvent(event)


_HANDLE_CURSORS = {"nw": Qt.SizeFDiagCursor, "se": Qt.SizeFDiagCursor,
                   "ne": Qt.SizeBDiagCursor, "sw": Qt.SizeBDiagCursor,
                   "n": Qt.SizeVerCursor, "s": Qt.SizeVerCursor,
                   "e": Qt.SizeHorCursor, "w": Qt.SizeHorCursor,
                   "p0": Qt.CrossCursor, "p1": Qt.CrossCursor}


def _segment_distance(x, y, x0, y0, x1, y1) -> float:
    """Distance from a point to a line segment."""
    dx, dy = x1 - x0, y1 - y0
    length = dx * dx + dy * dy
    t = 0.0 if length <= 0 else max(0.0, min(1.0, ((x - x0) * dx + (y - y0) * dy) / length))
    return float(np.hypot(x - (x0 + t * dx), y - (y0 + t * dy)))


# --------------------------------------------------------------------------
# Small shared widgets
# --------------------------------------------------------------------------
class ColorButton(QPushButton):
    """A button that shows, and picks, a colour."""

    colorChanged = pyqtSignal(str)

    def __init__(self, color="#ffffff", parent=None):
        super().__init__(parent)
        self._color = color
        self.setFixedWidth(52)
        self.clicked.connect(self._choose)
        self._refresh()

    def _refresh(self):
        self.setStyleSheet(f"background-color: {self._color};")
        self.setText("")

    def color(self):
        return self._color

    def set_color(self, color):
        self._color = str(color)
        self._refresh()

    def _choose(self):
        chosen = QColorDialog.getColor(QColor(self._color), self, "Colour")
        if chosen.isValid():
            self.set_color(chosen.name())
            self.colorChanged.emit(self._color)


def spin(value, lo, hi, decimals=2, step=0.1, tip="", width=90):
    box = QDoubleSpinBox()
    box.setDecimals(decimals)
    box.setRange(lo, hi)
    box.setSingleStep(step)
    box.setValue(value)
    box.setKeyboardTracking(False)
    box.setMaximumWidth(width)
    if tip:
        box.setToolTip(tip)
    return box


def whole(value, lo, hi, tip="", width=90):
    box = QSpinBox()
    box.setRange(lo, hi)
    box.setValue(int(value))
    box.setMaximumWidth(width)
    if tip:
        box.setToolTip(tip)
    return box


# --------------------------------------------------------------------------
# The tools
# --------------------------------------------------------------------------
class FigureToolsDialog(QDialog):
    """Everything that dresses a figure up, in six tabs.

    Each control writes into the model and asks the canvas to redraw, so
    there is no "apply" step and nothing can be set in the dialog but
    missing from the file. "Apply to all panels" makes a control act on
    every panel at once -- the MATLAB original had to re-implement that loop
    inside each of its forty callbacks.
    """

    def __init__(self, window: "FigureWindow", embedded: bool = False):
        super().__init__(window)
        self.window_ = window
        self.embedded = embedded
        self.setWindowTitle("Figure tools")
        self.setModal(False)
        if embedded:
            # Part of the composer's window, beside the page.
            self.setWindowFlags(Qt.Widget)
        self._syncing = False

        outer = QVBoxLayout(self)
        top = QHBoxLayout()
        top.addWidget(QLabel("Panel"))
        self.panel_combo = QComboBox()
        self.panel_combo.currentIndexChanged.connect(self._panel_chosen)
        top.addWidget(self.panel_combo, 1)
        self.all_box = QCheckBox("Apply to all panels")
        self.all_box.setToolTip(
            "Every control below then acts on every panel, which is how a "
            "figure's panels stay consistent with one another.")
        top.addWidget(self.all_box)
        outer.addLayout(top)

        self.tabs = QTabWidget()
        for title, builder in (("Panels", self._build_panels),
                               ("Axes", self._build_axes),
                               ("Text", self._build_text),
                               ("Colour", self._build_colour),
                               ("Overlays", self._build_overlays),
                               ("Export", self._build_export)):
            page = QWidget()
            scroller = QScrollArea()
            scroller.setWidgetResizable(True)
            scroller.setFrameShape(QFrame.NoFrame)
            scroller.setWidget(page)
            builder(page)
            self.tabs.addTab(scroller, title)
        outer.addWidget(self.tabs, 1)

        self.status = QLabel("")
        self.status.setWordWrap(True)
        outer.addWidget(self.status)
        if not embedded:
            close = QDialogButtonBox(QDialogButtonBox.Close)
            close.rejected.connect(self.reject)
            outer.addWidget(close)
            self.resize(560, 700)
        else:
            outer.setContentsMargins(0, 0, 0, 0)
        self.setMinimumWidth(400)
        self.sync()

    @property
    def figure(self):
        # The window's, always: undo puts a new style object on it.
        return self.window_.figure

    # -- plumbing ----------------------------------------------------------
    def panel(self):
        index = self.panel_combo.currentIndex()
        if 0 <= index < len(self.figure.panels):
            return self.figure.panels[index]
        return None

    def targets(self):
        if self.all_box.isChecked():
            return list(self.figure.panels)
        panel = self.panel()
        return [panel] if panel is not None else []

    def apply(self, setter):
        """Run ``setter(panel)`` on whatever the dialog is pointed at."""
        if self._syncing:
            return
        self._remember()
        for panel in self.targets():
            setter(panel)
        self.redraw()

    def _remember(self):
        remember = getattr(self.window_, "remember", None)
        if remember is not None:
            remember("tools")

    def redraw(self):
        self.window_.refresh()
        width, height = self.window_.painter.size_mm()
        pixels = self.window_.painter.pixel_size()
        self.status.setText(
            f"{width:.1f} × {height:.1f} mm — "
            f"{pixels[0]} × {pixels[1]} px at {self.figure.style.dpi} dpi")

    def sync(self):
        """Put the model's values into the widgets."""
        self._syncing = True
        try:
            current = self.panel_combo.currentIndex()
            self.panel_combo.clear()
            for i, panel in enumerate(self.figure.panels):
                name = panel.data.name or ("empty" if panel.data.kind == "empty"
                                           else f"panel {i + 1}")
                if panel.role == "inset":
                    letter = "inset"
                else:
                    letter = panel.label or f"{i + 1}."
                self.panel_combo.addItem(f"{letter}  {name}")
            if 0 <= current < self.panel_combo.count():
                self.panel_combo.setCurrentIndex(current)
            elif self.panel_combo.count():
                self.panel_combo.setCurrentIndex(self.window_.canvas.selected)
            self._sync_figure()
            self._sync_panel()
        finally:
            self._syncing = False
        self.redraw()

    def _panel_chosen(self, index):
        if self._syncing:
            return
        self.window_.canvas.selected = max(0, index)
        self.window_.refresh()
        self._syncing = True
        try:
            self._sync_panel()
        finally:
            self._syncing = False

    # ------------------------------------------------------------------
    # Panels tab: the figure itself
    # ------------------------------------------------------------------
    def _build_panels(self, page):
        style = self.figure.style
        outer = QVBoxLayout(page)

        group = QGroupBox("Page")
        form = QFormLayout(group)
        self.preset_combo = QComboBox()
        self.preset_combo.addItem("(custom)")
        for label, width, font in JOURNAL_PRESETS:
            self.preset_combo.addItem(f"{label} — {width:.0f} mm")
        self.preset_combo.setToolTip(
            "A journal's column width and a font size that suits it. The "
            "figure is built at that size, so nothing is scaled afterwards "
            "and the type comes out the size it says.")
        self.preset_combo.currentIndexChanged.connect(self._apply_preset)
        form.addRow("Journal", self.preset_combo)

        self.width_box = spin(style.width_mm, 20, 500, 1, 1.0, "Page width.")
        self.height_box = spin(style.height_mm, 0, 800, 1, 1.0,
                               "Page height; 0 lets the panels choose it.")
        self.dpi_box = whole(style.dpi, 72, 2400,
                             "Resolution of the data raster inside the figure.")
        for box, name in ((self.width_box, "width_mm"), (self.height_box, "height_mm")):
            box.valueChanged.connect(lambda v, n=name: self._set_style(n, v))
        self.dpi_box.valueChanged.connect(lambda v: self._set_style("dpi", int(v)))
        form.addRow("Width (mm)", self.width_box)
        form.addRow("Height (mm)", self.height_box)
        form.addRow("Resolution (dpi)", self.dpi_box)
        outer.addWidget(group)

        group = QGroupBox("Arrange as a grid")
        form = QFormLayout(group)
        grid_row = QWidget()
        grid_box = QHBoxLayout(grid_row)
        grid_box.setContentsMargins(0, 0, 0, 0)
        self.rows_box = whole(self.figure.rows, 1, 12, "Rows of panels.", 60)
        self.cols_box = whole(self.figure.cols, 1, 12, "Columns of panels.", 60)
        self.rows_box.valueChanged.connect(self._set_grid)
        self.cols_box.valueChanged.connect(self._set_grid)
        grid_box.addWidget(QLabel("rows"))
        grid_box.addWidget(self.rows_box)
        grid_box.addWidget(QLabel("cols"))
        grid_box.addWidget(self.cols_box)
        grid_box.addStretch(1)
        form.addRow("Grid", grid_row)
        self.gap_x = spin(style.gap_x_mm, 0, 40, 1, 0.5, "Horizontal gap between panels.")
        self.gap_y = spin(style.gap_y_mm, 0, 40, 1, 0.5, "Vertical gap between panels.")
        self.gap_x.valueChanged.connect(lambda v: self._set_style("gap_x_mm", v))
        self.gap_y.valueChanged.connect(lambda v: self._set_style("gap_y_mm", v))
        form.addRow("Column gap (mm)", self.gap_x)
        form.addRow("Row gap (mm)", self.gap_y)
        self.aspect_box = spin(style.panel_aspect, 0, 8, 2, 0.05,
                               "Panel height / width; 0 chooses it.")
        self.aspect_box.valueChanged.connect(lambda v: self._set_style("panel_aspect", v))
        form.addRow("Panel aspect", self.aspect_box)
        self.fill_slots = QCheckBox("Fill the grid with empty panels")
        self.fill_slots.setChecked(True)
        self.fill_slots.setToolTip(
            "Rows × columns slots in all: the ones without a dataset stay "
            "empty, to be filled from the list or with a picture.")
        form.addRow("", self.fill_slots)
        arrange = QPushButton("Arrange as grid")
        arrange.setToolTip(
            "Put the panels on a regular grid of this many rows and "
            "columns, in reading order. Insets move with their panel.")
        arrange.clicked.connect(self._arrange)
        auto = QPushButton("Auto rows / cols")
        auto.clicked.connect(self._auto_grid)
        form.addRow("", self._pair(arrange, auto))
        outer.addWidget(group)

        group = QGroupBox("All panels")
        form = QFormLayout(group)
        self.shared_levels = QCheckBox("Same colour scale everywhere")
        self.shared_levels.setToolTip(
            "One colour window over all the panels, so their intensities can "
            "be compared by eye instead of each being stretched to its own "
            "maximum.")
        self.shared_ranges = QCheckBox("Same axis ranges everywhere")
        self.edge_labels = QCheckBox("Tick labels on outer edges only")
        self.edge_labels.setToolTip(
            "The usual look for a grid of panels: the y axis is labelled on "
            "the left column and the x axis on the bottom row.")
        self.equal_aspect = QCheckBox("Keep the data's aspect (1:1)")
        for box, name in ((self.shared_levels, "shared_levels"),
                          (self.shared_ranges, "shared_ranges"),
                          (self.edge_labels, "edge_labels_only"),
                          (self.equal_aspect, "equal_aspect")):
            box.toggled.connect(lambda v, n=name: self._set_style(n, bool(v)))
            form.addRow("", box)
        outer.addWidget(group)

        group = QGroupBox("Selected panel")
        form = QFormLayout(group)
        self.content_label = QLabel("")
        self.content_label.setWordWrap(True)
        form.addRow("Content", self.content_label)
        loads = QWidget()
        row = QHBoxLayout(loads)
        row.setContentsMargins(0, 0, 0, 0)
        for text, slot in (("From list...", lambda: self.window_.load_from_list()),
                           ("Picture...", lambda: self.window_.load_picture()),
                           ("Empty", lambda: self.window_.empty_panel())):
            button = QPushButton(text)
            button.clicked.connect(slot)
            row.addWidget(button)
        row.addStretch(1)
        form.addRow("Load", loads)
        self.pos_x = spin(0, -500, 1000, 2, 0.5, "Left edge of the data box (mm).", 80)
        self.pos_y = spin(0, -500, 1000, 2, 0.5, "Top edge of the data box (mm).", 80)
        self.pos_w = spin(10, 2, 1000, 2, 0.5, "Width of the data box (mm).", 80)
        self.pos_h = spin(10, 2, 1000, 2, 0.5, "Height of the data box (mm).", 80)
        for box in (self.pos_x, self.pos_y, self.pos_w, self.pos_h):
            box.valueChanged.connect(self._set_frame)
        form.addRow("x / y (mm)", self._pair(self.pos_x, self.pos_y))
        form.addRow("w / h (mm)", self._pair(self.pos_w, self.pos_h))
        self.role_combo = QComboBox()
        self.role_combo.addItems(["panel (lettered)", "inset (no letter)"])
        self.role_combo.currentIndexChanged.connect(self._set_role)
        form.addRow("Kind", self.role_combo)
        self.xlab_box = QCheckBox("x numbers and label")
        self.ylab_box = QCheckBox("y numbers and label")
        self.xlab_box.toggled.connect(lambda v: self._set_shown("x_labels", v))
        self.ylab_box.toggled.connect(lambda v: self._set_shown("y_labels", v))
        form.addRow("Show", self._pair(self.xlab_box, self.ylab_box))
        self.panel_font = spin(0.0, 0, 40, 1, 0.5,
                               "This panel's font size; 0 follows the figure "
                               "(an inset usually wants a smaller one).", 80)
        self.panel_font.valueChanged.connect(
            lambda v: self.apply(lambda p: setattr(p, "font_pt", float(v))))
        form.addRow("Font size (pt)", self.panel_font)

        buttons = QWidget()
        row = QHBoxLayout(buttons)
        row.setContentsMargins(0, 0, 0, 0)
        add = QPushButton("Add panel...")
        add.setToolTip("A new panel holding a dataset from the main list.")
        add.clicked.connect(self.window_.add_panel_from_list)
        remove = QPushButton("Remove")
        remove.clicked.connect(self._remove_panel)
        letter = QPushButton("Re-letter")
        letter.setToolTip("(a), (b), (c)... in reading order.")
        letter.clicked.connect(self._letter)
        row.addWidget(add)
        row.addWidget(remove)
        row.addWidget(letter)
        row.addStretch(1)
        form.addRow("", buttons)
        outer.addWidget(group)
        outer.addStretch(1)

    def _arrange(self):
        """Arrange the panels as the grid the boxes say, adding empty
        panels to fill it when asked."""
        figure = self.figure
        self._remember()
        rows, cols = int(self.rows_box.value()), int(self.cols_box.value())
        figure.freeze()
        if self.fill_slots.isChecked():
            while len(figure.main_panels()) < rows * cols:
                slot = Panel(data=nf.empty_data())
                slot.frame = (0.0, 1e6 + len(figure.panels), 1.0, 1.0)  # last in reading order
                figure.panels.append(slot)
        figure.arrange_grid(rows, cols)
        figure.letter_panels()
        self.sync()

    def _set_shown(self, name, value):
        """Tick numbers on or off for a panel: a page laid out by hand from
        then on, since a grid decides them by itself."""
        if self._syncing:
            return
        self.figure.freeze()
        self.apply(lambda p: setattr(p, name, bool(value)))

    def _set_frame(self, *_):
        if self._syncing:
            return
        panel = self.panel()
        if panel is None:
            return
        self._remember()
        self.figure.freeze()
        self.figure.move_panel(panel, (self.pos_x.value(), self.pos_y.value(),
                                       self.pos_w.value(), self.pos_h.value()))
        self.redraw()

    def _set_role(self, index):
        if self._syncing:
            return
        panel = self.panel()
        if panel is None:
            return
        self._remember()
        self.figure.freeze()
        panel.role = "inset" if index == 1 else "panel"
        if panel.role == "panel":
            panel.zoom_of = ""
        self.figure.letter_panels()
        self.sync()

    def _set_style(self, name, value):
        if self._syncing:
            return
        self._remember()
        setattr(self.figure.style, name, value)
        if name == "edge_labels_only" and self.figure.free:
            self.figure.auto_edge_labels()
            self._resync()
        self.redraw()

    def _set_grid(self, *_):
        if self._syncing or self.figure.free:
            return          # a page laid out by hand changes on "Arrange"
        self._remember()
        self.figure.rows = int(self.rows_box.value())
        self.figure.cols = int(self.cols_box.value())
        self.redraw()

    def _auto_grid(self):
        count = len(self.figure.main_panels())
        cols = min(count, max(1, int(round(np.sqrt(count * 1.4))))) if count else 1
        if self.figure.free:
            self._syncing = True
            self.cols_box.setValue(cols)
            self.rows_box.setValue(int(np.ceil(max(count, 1) / cols)))
            self._syncing = False
            return
        self.figure.fit_grid()
        self._syncing = True
        self.rows_box.setValue(self.figure.rows)
        self.cols_box.setValue(self.figure.cols)
        self._syncing = False
        self.redraw()

    def _apply_preset(self, index):
        if self._syncing or index <= 0:
            return
        _label, width, font = JOURNAL_PRESETS[index - 1]
        self.figure.style.width_mm = width
        self.figure.style.font_pt = font
        self._syncing = True
        self.width_box.setValue(width)
        self.font_size.setValue(font)
        self._syncing = False
        self.redraw()

    def _letter(self):
        self._remember()
        self.figure.letter_panels()
        self.sync()

    def _remove_panel(self):
        index = self.panel_combo.currentIndex()
        if len(self.figure.panels) <= 1 or index < 0:
            return
        if self.figure.free:
            self.window_.canvas.selection = ("panel", index)
            self.window_.canvas.delete_selection()
            return
        self._remember()
        self.figure.panels.pop(index)
        self.figure.fit_grid(self.figure.cols)
        self.window_.canvas.selected = min(index, len(self.figure.panels) - 1)
        self.sync()

    def _move(self, delta):
        index = self.panel_combo.currentIndex()
        target = index + delta
        if not (0 <= index < len(self.figure.panels)) or not (0 <= target < len(self.figure.panels)):
            return
        panels = self.figure.panels
        panels[index], panels[target] = panels[target], panels[index]
        self.window_.canvas.selected = target
        self.sync()
        self.panel_combo.setCurrentIndex(target)

    def _sync_figure(self):
        style = self.figure.style
        self.width_box.setValue(style.width_mm)
        self.height_box.setValue(style.height_mm)
        self.dpi_box.setValue(int(style.dpi))
        self.rows_box.setValue(self.figure.rows)
        self.cols_box.setValue(self.figure.cols)
        self.gap_x.setValue(style.gap_x_mm)
        self.gap_y.setValue(style.gap_y_mm)
        self.shared_levels.setChecked(style.shared_levels)
        self.shared_ranges.setChecked(style.shared_ranges)
        self.edge_labels.setChecked(style.edge_labels_only)
        self.equal_aspect.setChecked(style.equal_aspect)
        self.aspect_box.setValue(style.panel_aspect)
        self.font_size.setValue(style.font_pt)
        self.font_combo.setCurrentText(style.font)
        self.tick_dir.setCurrentText(style.tick_direction)
        self.tick_len.setValue(style.tick_len_mm)
        self.minor_box.setValue(int(style.minor_ticks))
        self.tick_sides.setCurrentText(style.tick_sides)
        self.axis_width.setValue(style.axis_width_pt)
        self.axis_color.set_color(style.axis_color)
        self.background.set_color(style.background)
        self.figure_title.setText(self.figure.title)
        self.label_position.setCurrentText(style.panel_label_position)
        self.label_inside.setChecked(style.panel_label_inside)
        self.label_aligned.setChecked(getattr(style, "panel_label_aligned", True))
        self.label_bold.setChecked(style.panel_label_bold)
        self.label_color.set_color(style.panel_label_color)

    # ------------------------------------------------------------------
    # Axes tab
    # ------------------------------------------------------------------
    def _build_axes(self, page):
        style = self.figure.style
        form = QFormLayout(page)

        self.x_from = spin(0, -1e9, 1e9, 4, 0.1, "Left edge of this panel.")
        self.x_to = spin(1, -1e9, 1e9, 4, 0.1, "Right edge.")
        self.y_from = spin(0, -1e9, 1e9, 4, 0.1, "Bottom edge.")
        self.y_to = spin(1, -1e9, 1e9, 4, 0.1, "Top edge.")
        for box in (self.x_from, self.x_to, self.y_from, self.y_to):
            box.valueChanged.connect(self._set_ranges)
        form.addRow("x from / to", self._pair(self.x_from, self.x_to))
        form.addRow("y from / to", self._pair(self.y_from, self.y_to))
        reset = QPushButton("Whole data")
        reset.clicked.connect(self._reset_ranges)
        form.addRow("", reset)

        self.invert_x = QCheckBox("Reverse x")
        self.invert_y = QCheckBox("Reverse y")
        self.invert_x.toggled.connect(
            lambda v: self.apply(lambda p: setattr(p, "invert_x", bool(v))))
        self.invert_y.toggled.connect(
            lambda v: self.apply(lambda p: setattr(p, "invert_y", bool(v))))
        form.addRow("", self.invert_x)
        form.addRow("", self.invert_y)

        self.tick_dir = QComboBox()
        self.tick_dir.addItems(["in", "out"])
        self.tick_dir.currentTextChanged.connect(
            lambda v: self._set_style("tick_direction", v))
        form.addRow("Tick direction", self.tick_dir)
        self.tick_len = spin(style.tick_len_mm, 0.1, 6, 2, 0.1, "Tick length.")
        self.tick_len.valueChanged.connect(lambda v: self._set_style("tick_len_mm", v))
        form.addRow("Tick length (mm)", self.tick_len)
        self.minor_box = whole(style.minor_ticks, 0, 9,
                               "Minor ticks between each labelled one.")
        self.minor_box.valueChanged.connect(
            lambda v: self._set_style("minor_ticks", int(v)))
        form.addRow("Minor ticks", self.minor_box)
        self.tick_sides = QComboBox()
        self.tick_sides.addItems(["all", "left-bottom"])
        self.tick_sides.setToolTip("Ticks on all four sides, or only where "
                                   "the numbers are.")
        self.tick_sides.currentTextChanged.connect(
            lambda v: self._set_style("tick_sides", v))
        form.addRow("Tick sides", self.tick_sides)

        self.x_step = spin(style.x_tick_step, 0, 1e6, 4, 0.1,
                           "Spacing of labelled x ticks; 0 chooses round ones.")
        self.y_step = spin(style.y_tick_step, 0, 1e6, 4, 0.1,
                           "Spacing of labelled y ticks; 0 chooses round ones.")
        self.x_step.valueChanged.connect(lambda v: self._set_style("x_tick_step", v))
        self.y_step.valueChanged.connect(lambda v: self._set_style("y_tick_step", v))
        form.addRow("x tick step", self.x_step)
        form.addRow("y tick step", self.y_step)

        self.axis_width = spin(style.axis_width_pt, 0.1, 4, 2, 0.1, "Frame line width.")
        self.axis_width.valueChanged.connect(
            lambda v: self._set_style("axis_width_pt", v))
        form.addRow("Axis width (pt)", self.axis_width)
        self.axis_color = ColorButton(style.axis_color)
        self.axis_color.colorChanged.connect(lambda c: self._set_style("axis_color", c))
        form.addRow("Axis colour", self.axis_color)
        self.background = ColorButton(style.background)
        self.background.colorChanged.connect(lambda c: self._set_style("background", c))
        form.addRow("Background", self.background)

        for side, attribute in (("Top", "top_axis"), ("Right", "right_axis")):
            group = QGroupBox(f"{side} axis (a second scale)")
            group.setToolTip(
                "A second scale on the far edge: binding against kinetic "
                "energy, Å⁻¹ against degrees, k in units of "
                "π/a. The second value is scale × first + offset.")
            inner = QFormLayout(group)
            enabled = QCheckBox("Show")
            scale_box = spin(1.0, -1e6, 1e6, 6, 0.1, "Second = scale × first + offset.")
            offset_box = spin(0.0, -1e9, 1e9, 6, 0.1, "")
            label_box = QLineEdit()
            inner.addRow("", enabled)
            inner.addRow("Scale", scale_box)
            inner.addRow("Offset", offset_box)
            inner.addRow("Label", label_box)
            form.addRow(group)

            def make(attr, en, sc, off, lab):
                def update(*_):
                    self.apply(lambda p: setattr(
                        p, attr, nf.SecondaryAxis(en.isChecked(), sc.value(),
                                                  off.value(), lab.text())))
                return update
            update = make(attribute, enabled, scale_box, offset_box, label_box)
            enabled.toggled.connect(update)
            scale_box.valueChanged.connect(update)
            offset_box.valueChanged.connect(update)
            label_box.editingFinished.connect(update)
            setattr(self, f"{attribute}_widgets",
                    (enabled, scale_box, offset_box, label_box))

    @staticmethod
    def _pair(first, second):
        holder = QWidget()
        box = QHBoxLayout(holder)
        box.setContentsMargins(0, 0, 0, 0)
        box.addWidget(first)
        box.addWidget(second)
        box.addStretch(1)
        return holder

    def _set_ranges(self, *_):
        if self._syncing:
            return
        x = (self.x_from.value(), self.x_to.value())
        y = (self.y_from.value(), self.y_to.value())
        if x[1] <= x[0] or y[1] <= y[0]:
            return
        self.apply(lambda p: (setattr(p, "x_range", x), setattr(p, "y_range", y)))

    def _reset_ranges(self):
        self.apply(lambda p: (setattr(p, "x_range", None), setattr(p, "y_range", None)))
        self._syncing = True
        self._sync_panel()
        self._syncing = False

    # ------------------------------------------------------------------
    # Text tab
    # ------------------------------------------------------------------
    def _build_text(self, page):
        style = self.figure.style
        form = QFormLayout(page)
        hint = QLabel("Sub- and superscripts with _ and ^ : "
                      "<tt>k_{||} (Å^-1)</tt>, <tt>E - E_F</tt>.")
        hint.setWordWrap(True)
        form.addRow(hint)

        self.title_edit = QLineEdit()
        self.xlabel_edit = QLineEdit()
        self.ylabel_edit = QLineEdit()
        self.title_edit.editingFinished.connect(
            lambda: self.apply(lambda p: setattr(p, "title", self.title_edit.text())))
        self.xlabel_edit.editingFinished.connect(
            lambda: self.apply(lambda p: setattr(p, "x_label", self.xlabel_edit.text())))
        self.ylabel_edit.editingFinished.connect(
            lambda: self.apply(lambda p: setattr(p, "y_label", self.ylabel_edit.text())))
        form.addRow("Panel title", self.title_edit)
        form.addRow("x label", self.xlabel_edit)
        form.addRow("y label", self.ylabel_edit)

        self.figure_title = QLineEdit()
        self.figure_title.editingFinished.connect(self._set_figure_title)
        form.addRow("Figure title", self.figure_title)

        self.font_combo = QComboBox()
        self.font_combo.setEditable(True)
        self.font_combo.addItems(["Helvetica", "Arial", "DejaVu Sans",
                                  "Times New Roman", "Nimbus Sans"])
        self.font_combo.currentTextChanged.connect(lambda v: self._set_style("font", v))
        form.addRow("Font", self.font_combo)
        self.font_size = spin(style.font_pt, 3, 40, 1, 0.5, "Base font size.")
        self.font_size.valueChanged.connect(lambda v: self._set_style("font_pt", v))
        form.addRow("Font size (pt)", self.font_size)
        self.label_size = spin(style.label_pt, 0, 40, 1, 0.5,
                               "Axis labels; 0 follows the base size.")
        self.tick_size = spin(style.tick_pt, 0, 40, 1, 0.5,
                              "Tick numbers; 0 follows the base size.")
        self.annot_size = spin(style.annotation_pt, 0, 40, 1, 0.5,
                               "Annotations; 0 follows the base size.")
        self.label_size.valueChanged.connect(lambda v: self._set_style("label_pt", v))
        self.tick_size.valueChanged.connect(lambda v: self._set_style("tick_pt", v))
        self.annot_size.valueChanged.connect(lambda v: self._set_style("annotation_pt", v))
        form.addRow("Label size (pt)", self.label_size)
        form.addRow("Tick size (pt)", self.tick_size)
        form.addRow("Annotation size", self.annot_size)

        group = QGroupBox("Panel letter")
        inner = QFormLayout(group)
        self.label_edit = QLineEdit()
        self.label_edit.editingFinished.connect(
            lambda: self.apply(lambda p: setattr(p, "label", self.label_edit.text())))
        inner.addRow("Text", self.label_edit)
        nudge = QHBoxLayout()
        self.label_dx = spin(0.0, -50, 50, 1, 0.5, "Move this panel's letter right (mm).")
        self.label_dy = spin(0.0, -50, 50, 1, 0.5, "Move it down (mm).")
        for box in (self.label_dx, self.label_dy):
            box.valueChanged.connect(lambda _v: self.apply(lambda p: setattr(
                p, "label_offset", (self.label_dx.value(), self.label_dy.value()))))
        nudge.addWidget(QLabel("x"))
        nudge.addWidget(self.label_dx)
        nudge.addWidget(QLabel("y"))
        nudge.addWidget(self.label_dy)
        inner.addRow("Nudge (mm)", nudge)
        self.label_position = QComboBox()
        self.label_position.addItems(list(CORNERS))
        self.label_position.currentTextChanged.connect(
            lambda v: self._set_style("panel_label_position", v))
        inner.addRow("Position", self.label_position)
        self.label_inside = QCheckBox("Inside the panel")
        self.label_inside.toggled.connect(
            lambda v: self._set_style("panel_label_inside", bool(v)))
        inner.addRow("", self.label_inside)
        self.label_aligned = QCheckBox("In line: one height per row, one edge per column")
        self.label_aligned.setToolTip(
            "Letters of panels side by side sit at the same height, of panels "
            "one above another at the same left edge, whatever titles or axes "
            "each panel has.")
        self.label_aligned.toggled.connect(
            lambda v: self._set_style("panel_label_aligned", bool(v)))
        inner.addRow("", self.label_aligned)
        self.label_bold = QCheckBox("Bold")
        self.label_bold.toggled.connect(
            lambda v: self._set_style("panel_label_bold", bool(v)))
        inner.addRow("", self.label_bold)
        self.label_pt = spin(style.panel_label_pt, 0, 40, 1, 0.5,
                             "0 follows the base size + 1.")
        self.label_pt.valueChanged.connect(
            lambda v: self._set_style("panel_label_pt", v))
        inner.addRow("Size (pt)", self.label_pt)
        self.label_color = ColorButton(style.panel_label_color)
        self.label_color.colorChanged.connect(
            lambda c: self._set_style("panel_label_color", c))
        inner.addRow("Colour (inside)", self.label_color)
        form.addRow(group)

    def _set_figure_title(self):
        self.figure.title = self.figure_title.text()
        self.redraw()

    # ------------------------------------------------------------------
    # Colour tab
    # ------------------------------------------------------------------
    def _build_colour(self, page):
        form = QFormLayout(page)
        self.cmap_combo = QComboBox()
        self.cmap_combo.addItems(colormaps.COLORMAP_NAMES)
        self.cmap_combo.currentTextChanged.connect(
            lambda v: self.apply(lambda p: setattr(p, "colormap", v)))
        form.addRow("Colormap", self.cmap_combo)
        self.flip_box = QCheckBox("Reverse")
        self.flip_box.toggled.connect(
            lambda v: self.apply(lambda p: setattr(p, "flip", bool(v))))
        form.addRow("", self.flip_box)

        self.level_lo = spin(0, -1e12, 1e12, 4, 1.0, "Bottom of the colour window.")
        self.level_hi = spin(1, -1e12, 1e12, 4, 1.0, "Top of the colour window.")
        for box in (self.level_lo, self.level_hi):
            box.valueChanged.connect(self._set_levels)
        form.addRow("Levels", self._pair(self.level_lo, self.level_hi))
        auto = QPushButton("From the data")
        auto.clicked.connect(lambda: (self.apply(lambda p: setattr(p, "levels", None)),
                                      self._resync()))
        form.addRow("", auto)

        self.gamma_box = spin(1.0, 0.05, 20, 2, 0.05,
                              "Below 1 lifts weak features out of the background.")
        self.gamma_box.valueChanged.connect(
            lambda v: self.apply(lambda p: setattr(p, "gamma", v)))
        form.addRow("Gamma", self.gamma_box)
        self.smooth_box = QCheckBox("Interpolate between data points")
        self.smooth_box.toggled.connect(
            lambda v: self.apply(lambda p: setattr(p, "smooth", bool(v))))
        form.addRow("", self.smooth_box)

        group = QGroupBox("Colour bar")
        inner = QFormLayout(group)
        self.cbar_box = QCheckBox("Show")
        self.cbar_box.toggled.connect(
            lambda v: self.apply(lambda p: setattr(p, "colorbar", bool(v))))
        inner.addRow("", self.cbar_box)
        self.cbar_label = QLineEdit()
        self.cbar_label.editingFinished.connect(
            lambda: self.apply(lambda p: setattr(p, "colorbar_label",
                                                 self.cbar_label.text())))
        inner.addRow("Label", self.cbar_label)
        self.cbar_ends = QComboBox()
        self.cbar_ends.addItems(["values", "high-low", "none"])
        self.cbar_ends.setToolTip("What the bar's two ends say: the levels, "
                                  "the words High and Low, or nothing.")
        self.cbar_ends.currentTextChanged.connect(
            lambda v: self.apply(lambda p: setattr(p, "colorbar_ends", v)))
        inner.addRow("Ends", self.cbar_ends)
        self.cbar_length = spin(100.0, 10.0, 100.0, 0, 5.0,
                                "The bar's length, % of the panel's height.")
        self.cbar_length.valueChanged.connect(
            lambda v: self.apply(lambda p: setattr(p, "colorbar_length", v / 100.0)))
        inner.addRow("Length (%)", self.cbar_length)
        self.cbar_position = QComboBox()
        self.cbar_position.addItems(["right", "below"])
        self.cbar_position.setToolTip("Upright beside the panel, or lying "
                                      "under it (centred).")
        self.cbar_position.currentTextChanged.connect(
            lambda v: self.apply(lambda p: setattr(p, "colorbar_position", v)))
        inner.addRow("Position", self.cbar_position)
        self.cbar_width = spin(self.figure.style.colorbar_mm, 0.5, 20, 1, 0.2, "")
        self.cbar_width.valueChanged.connect(
            lambda v: self._set_style("colorbar_mm", v))
        inner.addRow("Width (mm)", self.cbar_width)
        form.addRow(group)

        group = QGroupBox("Scale bar")
        group.setToolTip("For a real-space map, where the axes are usually "
                         "taken off and a bar of known length replaces them.")
        inner = QFormLayout(group)
        self.bar_box = QCheckBox("Show")
        self.bar_box.toggled.connect(self._set_scale_bar)
        inner.addRow("", self.bar_box)
        self.bar_length = spin(0, 0, 1e9, 4, 1.0,
                               "Length in data units; 0 picks a round one.")
        self.bar_position = QComboBox()
        self.bar_position.addItems(list(CORNERS))
        self.bar_color = ColorButton("#ffffff")
        self.bar_thickness = spin(2.5, 0.2, 20, 1, 0.2, "")
        self.bar_unit = QLineEdit()
        self.bar_unit.setPlaceholderText("from the x label")
        self.bar_label_box = QCheckBox("Write the length beside it")
        self.bar_label_box.setChecked(True)
        for widget, signal in ((self.bar_length, "valueChanged"),
                               (self.bar_thickness, "valueChanged"),
                               (self.bar_position, "currentTextChanged"),
                               (self.bar_color, "colorChanged"),
                               (self.bar_label_box, "toggled")):
            getattr(widget, signal).connect(self._set_scale_bar)
        self.bar_unit.editingFinished.connect(self._set_scale_bar)
        inner.addRow("Length", self.bar_length)
        inner.addRow("Position", self.bar_position)
        inner.addRow("Colour", self.bar_color)
        inner.addRow("Thickness (pt)", self.bar_thickness)
        inner.addRow("Unit", self.bar_unit)
        inner.addRow("", self.bar_label_box)
        form.addRow(group)

    def _set_levels(self, *_):
        if self._syncing:
            return
        lo, hi = self.level_lo.value(), self.level_hi.value()
        if hi <= lo:
            return
        self.apply(lambda p: setattr(p, "levels", (lo, hi)))

    def _resync(self):
        self._syncing = True
        try:
            self._sync_panel()
        finally:
            self._syncing = False
        self.redraw()

    def _set_scale_bar(self, *_):
        if self._syncing:
            return
        if not self.bar_box.isChecked():
            self.apply(lambda p: setattr(p, "scale_bar", None))
            return
        bar = ScaleBar(length=self.bar_length.value(),
                       position=self.bar_position.currentText(),
                       color=self.bar_color.color(),
                       unit=self.bar_unit.text(),
                       thickness_pt=self.bar_thickness.value(),
                       show_label=self.bar_label_box.isChecked())
        self.apply(lambda p: setattr(p, "scale_bar", ScaleBar(**vars(bar))))

    # ------------------------------------------------------------------
    # Overlays tab
    # ------------------------------------------------------------------
    def _build_overlays(self, page):
        layout = QVBoxLayout(page)

        group = QGroupBox("Reference lines")
        inner = QFormLayout(group)
        self.line_axis = QComboBox()
        self.line_axis.addItems(["horizontal (a y value)", "vertical (an x value)"])
        self.line_value = spin(0, -1e9, 1e9, 5, 0.05, "Where the line goes.")
        self.line_color = ColorButton("#ffffff")
        self.line_style = QComboBox()
        self.line_style.addItems(list(PEN_STYLES))
        self.line_style.setCurrentText("dashed")
        self.line_width = spin(0.8, 0.1, 6, 2, 0.1, "")
        self.line_label = QLineEdit()
        self.line_label.setPlaceholderText("optional, e.g. E_F")
        inner.addRow("Direction", self.line_axis)
        inner.addRow("Value", self.line_value)
        inner.addRow("Colour", self.line_color)
        inner.addRow("Style", self.line_style)
        inner.addRow("Width (pt)", self.line_width)
        inner.addRow("Label", self.line_label)
        row = QWidget()
        box = QHBoxLayout(row)
        box.setContentsMargins(0, 0, 0, 0)
        add = QPushButton("Add")
        add.clicked.connect(self._add_line)
        ef = QPushButton("From the fit")
        ef.setToolTip("Uses the fitted E_F this dataset carries, if it was "
                      "produced by the Fermi-level fit.")
        ef.clicked.connect(self._add_ef_line)
        clear = QPushButton("Clear")
        clear.clicked.connect(lambda: self.apply(lambda p: p.lines.clear()))
        box.addWidget(add)
        box.addWidget(ef)
        box.addWidget(clear)
        inner.addRow("", row)
        layout.addWidget(group)

        note = QLabel(
            "Arrows, lines, boxes and text are drawn straight on the page "
            "with the tools on the left of the window; a magnified region "
            "is the <b>Zoom-in inset</b> tool. The lines here sit at a "
            "value of the data (E_F, k = 0) and follow the panel's axes.")
        note.setWordWrap(True)
        layout.addWidget(note)

        group = QGroupBox("Curves drawn on the panel")
        inner = QVBoxLayout(group)
        note = QLabel(
            "Polylines in data coordinates — a Brillouin-zone boundary, "
            "a calculated band, a fitted dispersion. A Brillouin-zone panel "
            "will hand them over through this list; anything can, by "
            "appending nf.Overlay objects to the panel.")
        note.setWordWrap(True)
        inner.addWidget(note)
        self.overlay_list = QListWidget()
        self.overlay_list.setMaximumHeight(110)
        self.overlay_list.setSelectionMode(QListWidget.ExtendedSelection)
        self.overlay_list.setToolTip("Select curves to restyle them below "
                                     "(none selected: all of them).")
        inner.addWidget(self.overlay_list)
        style_row = QWidget()
        form = QFormLayout(style_row)
        form.setContentsMargins(0, 0, 0, 0)
        self.curve_colour = ColorButton("#000000")
        self.curve_width = spin(0.8, 0.1, 8, 2, 0.1, "Line width (pt).")
        self.curve_style = QComboBox()
        self.curve_style.addItems(list(PEN_STYLES))
        self.curve_marker = QComboBox()
        for text, code in (("none", ""), ("● circle", "o"), ("■ square", "s"),
                           ("▲ triangle", "^"), ("▼ triangle", "v"),
                           ("◆ diamond", "d")):
            self.curve_marker.addItem(text, code)
        self.curve_marker_size = spin(1.0, 0.2, 6, 2, 0.1, "Symbol size (mm).")
        self.curve_joined = QCheckBox("joined by a line")
        self.curve_joined.setChecked(True)
        self.curve_bars = QCheckBox("± σ as error bars")
        self.curve_bars.setToolTip(
            "A curve whose uncertainty came in as two dotted ± σ curves gets "
            "them as error bars at its points instead (and the dotted "
            "curves go).")
        marker_row = QHBoxLayout()
        marker_row.addWidget(self.curve_marker)
        marker_row.addWidget(self.curve_marker_size)
        form.addRow("Colour", self.curve_colour)
        form.addRow("Width (pt)", self.curve_width)
        form.addRow("Style", self.curve_style)
        form.addRow("Symbols", marker_row)
        form.addRow("", self.curve_joined)
        form.addRow("", self.curve_bars)
        restyle = QPushButton("Apply to the curves")
        restyle.clicked.connect(lambda: self.restyle_curves())
        form.addRow("", restyle)
        inner.addWidget(style_row)
        row = QWidget()
        box = QHBoxLayout(row)
        box.setContentsMargins(0, 0, 0, 0)
        copy_all = QPushButton("Copy to all panels")
        copy_all.setToolTip("The MATLAB tool's \"copy lines to all axes\".")
        copy_all.clicked.connect(self._copy_overlays)
        clear = QPushButton("Clear")
        clear.clicked.connect(lambda: (self.apply(lambda p: p.overlays.clear()),
                                       self._resync()))
        box.addWidget(copy_all)
        box.addWidget(clear)
        box.addStretch(1)
        inner.addWidget(row)
        layout.addWidget(group)
        layout.addStretch(1)

    def _add_line(self):
        axis = "y" if self.line_axis.currentIndex() == 0 else "x"
        line = Line(value=self.line_value.value(), axis=axis,
                    color=self.line_color.color(),
                    width_pt=self.line_width.value(),
                    style=self.line_style.currentText(),
                    label=self.line_label.text())
        self.apply(lambda p: p.lines.append(Line(**vars(line))))

    def _add_ef_line(self):
        panel = self.panel()
        if panel is None:
            return
        info = dict(panel.data.info or {})
        for key in ("fitEF.ef", "kcut.ef", "ef"):
            if key in info:
                try:
                    value = float(info[key])
                except (TypeError, ValueError):
                    continue
                self.line_value.setValue(value)
                self.line_label.setText("E_F")
                self._add_line()
                return
        QMessageBox.information(
            self, "Reference line",
            "This dataset does not carry a fitted E_F. Fit the Fermi level "
            "first, or type the value above.")

    def restyle_curves(self, rows=None):
        """Colour, width, dash, symbols and error bars for the selected
        curves of the panel (all of them when none is selected)."""
        panel = self.panel()
        if panel is None or not panel.overlays:
            return
        if rows is None:
            rows = sorted(i.row() for i in self.overlay_list.selectedIndexes())
        rows = rows or list(range(len(panel.overlays)))
        chosen = [panel.overlays[i] for i in rows if i < len(panel.overlays)]
        self._remember()
        removed = []
        for overlay in chosen:
            overlay.color = self.curve_colour.color()
            overlay.width_pt = self.curve_width.value()
            overlay.style = self.curve_style.currentText()
            overlay.marker = self.curve_marker.currentData() or ""
            overlay.marker_mm = self.curve_marker_size.value()
            overlay.line = self.curve_joined.isChecked()
            if self.curve_bars.isChecked() and not overlay.y_errors:
                plus = next((o for o in panel.overlays
                             if o.name == f"{overlay.name} +σ"), None)
                minus = next((o for o in panel.overlays
                              if o.name == f"{overlay.name} −σ"), None)
                if plus is not None and len(plus.ys) == len(overlay.ys):
                    overlay.y_errors = [float(p) - float(y)
                                        for p, y in zip(plus.ys, overlay.ys)]
                    removed += [plus] + ([minus] if minus is not None else [])
            elif not self.curve_bars.isChecked():
                overlay.y_errors = []
        panel.overlays = [o for o in panel.overlays if o not in removed]
        self._resync()
        self.redraw()

    def _copy_overlays(self):
        panel = self.panel()
        if panel is None:
            return
        source = [Overlay(**vars(o)) for o in panel.overlays]
        for other in self.figure.panels:
            if other is panel:
                continue
            other.overlays = [Overlay(**vars(o)) for o in source]
        self.redraw()

    # ------------------------------------------------------------------
    # Export tab
    # ------------------------------------------------------------------
    def _build_export(self, page):
        form = QFormLayout(page)
        note = QLabel(
            "The preview and the file come from the same painter, so what "
            "is on screen is what is written.")
        note.setWordWrap(True)
        form.addRow(note)

        self.format_combo = QComboBox()
        self.format_combo.addItems([
            "PDF (vector, data embedded at the dpi above)",
            "SVG (vector)",
            "PNG (raster, transparent background allowed)",
            "TIFF (raster)",
            "Two files: EPS axes + PNG data (for Illustrator)"])
        form.addRow("Format", self.format_combo)
        save = QPushButton("Save figure...")
        save.clicked.connect(self._save)
        form.addRow("", save)

        clipboard = QPushButton("Copy figure to the clipboard")
        clipboard.setToolTip("At the dpi above, for pasting into a draft.")
        clipboard.clicked.connect(self.window_.copy_to_clipboard)
        form.addRow("", clipboard)

        group = QGroupBox("Style")
        inner = QFormLayout(group)
        hint = QLabel("Save how this figure looks and apply it to the next "
                      "one, so every figure in a paper matches. The data is "
                      "not part of it.")
        hint.setWordWrap(True)
        inner.addRow(hint)
        row = QWidget()
        box = QHBoxLayout(row)
        box.setContentsMargins(0, 0, 0, 0)
        save_style = QPushButton("Save style...")
        save_style.clicked.connect(self._save_style)
        load_style = QPushButton("Load style...")
        load_style.clicked.connect(self._load_style)
        box.addWidget(save_style)
        box.addWidget(load_style)
        inner.addRow("", row)
        form.addRow(group)

    def _save(self):
        index = self.format_combo.currentIndex()
        self.window_.export_figure(("pdf", "svg", "png", "tif", "layers")[index])

    def _save_style(self):
        path, _ = QFileDialog.getSaveFileName(
            self, "Save the figure style", "figure_style.json", "JSON (*.json)")
        if not path:
            return
        if not path.lower().endswith(".json"):
            path += ".json"
        self.figure.save_style(path)
        self.status.setText(f"Style written to {os.path.basename(path)}")

    def _load_style(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Load a figure style", "", "JSON (*.json)")
        if not path:
            return
        try:
            self.figure.load_style(path)
        except Exception as exc:
            QMessageBox.warning(self, "Style", f"Could not read that style:\n{exc}")
            return
        self.sync()

    # ------------------------------------------------------------------
    def _sync_panel(self):
        panel = self.panel()
        if panel is None:
            return
        (x0, x1), (y0, y1) = self.figure.panel_ranges(panel)
        self.x_from.setValue(x0)
        self.x_to.setValue(x1)
        self.y_from.setValue(y0)
        self.y_to.setValue(y1)
        self.invert_x.setChecked(panel.invert_x)
        self.invert_y.setChecked(panel.invert_y)
        self.title_edit.setText(panel.title)
        self.xlabel_edit.setText(panel.label_x())
        self.ylabel_edit.setText(panel.label_y())
        self.label_edit.setText(panel.label)
        dx, dy = getattr(panel, "label_offset", (0.0, 0.0)) or (0.0, 0.0)
        self.label_dx.setValue(float(dx))
        self.label_dy.setValue(float(dy))
        self.cmap_combo.setCurrentText(panel.colormap)
        self.flip_box.setChecked(panel.flip)
        lo, hi = self.figure.panel_levels(panel)
        self.level_lo.setValue(lo)
        self.level_hi.setValue(hi)
        self.gamma_box.setValue(panel.gamma)
        self.smooth_box.setChecked(panel.smooth)
        self.cbar_box.setChecked(panel.colorbar)
        self.cbar_label.setText(panel.colorbar_label)
        self.cbar_ends.setCurrentText(getattr(panel, "colorbar_ends", "values"))
        self.cbar_length.setValue(100.0 * float(getattr(panel, "colorbar_length", 1.0)))
        self.cbar_position.setCurrentText(getattr(panel, "colorbar_position", "right"))
        bar = panel.scale_bar
        self.bar_box.setChecked(bar is not None)
        if bar is not None:
            self.bar_length.setValue(bar.length)
            self.bar_position.setCurrentText(bar.position)
            self.bar_color.set_color(bar.color)
            self.bar_thickness.setValue(bar.thickness_pt)
            self.bar_unit.setText(bar.unit)
            self.bar_label_box.setChecked(bar.show_label)
        kind = panel.data.kind
        self.content_label.setText(
            "empty -- load a dataset or a picture" if kind == "empty"
            else f"picture: {panel.data.name}" if kind == "picture"
            else f"{panel.data.name or kind}"
            + (" (magnified region)" if panel.zoom_of else ""))
        frame = self.window_.canvas._frame_of(panel) \
            if hasattr(self.window_, "canvas") else panel.frame
        if frame is not None:
            for box, value in zip((self.pos_x, self.pos_y, self.pos_w, self.pos_h), frame):
                box.setValue(float(value))
        self.role_combo.setCurrentIndex(1 if panel.role == "inset" else 0)
        index = self.figure.panels.index(panel)
        self.xlab_box.setChecked(self.figure.shows_x_labels(index))
        self.ylab_box.setChecked(self.figure.shows_y_labels(index))
        self.panel_font.setValue(panel.font_pt)
        for attribute in ("top_axis", "right_axis"):
            axis = getattr(panel, attribute) or nf.SecondaryAxis()
            enabled, scale_box, offset_box, label_box = getattr(
                self, f"{attribute}_widgets")
            enabled.setChecked(axis.enabled)
            scale_box.setValue(axis.scale)
            offset_box.setValue(axis.offset)
            label_box.setText(axis.label)
        self.overlay_list.clear()
        for overlay in panel.overlays:
            self.overlay_list.addItem(
                overlay.name or f"{len(overlay.xs)} points")

    def closeEvent(self, event):
        if not self.embedded:
            self.window_._tools = None
        super().closeEvent(event)


# --------------------------------------------------------------------------
# Turning this program's datasets into panels
# --------------------------------------------------------------------------
def panel_from_arrays(array, x, y, x_label="", y_label="", name="", info=None,
                      colormap="gray", flip=False, smooth=True,
                      sources=None) -> Panel:
    """A panel of an image. ``sources`` names what it was made from
    (:func:`tools.output.figfile.source_record`); by default the origin
    ``info`` records."""
    from tools.output.figfile import sources_from_info

    data = PanelData(array=np.asarray(array, dtype=float),
                     x=np.asarray(x, dtype=float), y=np.asarray(y, dtype=float),
                     x_label=x_label, y_label=y_label, name=name,
                     info=dict(info or {}))
    data.sources = list(sources) if sources is not None else sources_from_info(info)
    return Panel(data=data, colormap=colormap, flip=flip, smooth=smooth)


def panel_from_dataset(data, name="") -> Panel:
    """A Cut (or anything with a 2-D ``scan.value``) as a panel; a curve
    dataset as its curves, with ±σ where it carries uncertainties."""
    scan = data.scan
    from loader.nxs_file import CURVE_KINDS
    if getattr(data, "kind", None) in CURVE_KINDS:
        from tools.spectra import curves as C
        from ui.viewers.curves import curve_panel, _colour

        values = np.asarray(scan.value, dtype=float)
        values = values[:, None] if values.ndim == 1 else values
        names = C.channel_names(scan.info, values.shape[1], data.kind)
        columns, errors = [], {}
        for index in C.data_channels(names):
            columns.append((names[index], values[:, index], _colour(index)))
            sigma = C.sigma_of(names, index)
            if sigma is not None:
                errors[names[index]] = values[:, sigma]
        labels = getattr(scan, "labels", {}) or {}
        panel = curve_panel(scan.x, columns, x_label=labels.get("x", ""),
                            y_label=C.value_label(scan.info),
                            name=name or getattr(data, "source_label", "") or "curve",
                            errors=errors)
        panel.data.sources = _dataset_sources(data)
        return panel
    if getattr(scan, "value", None) is None or np.asarray(scan.value).ndim != 2:
        raise ValueError("only a two-dimensional dataset can be a panel on "
                         "its own; take a slice of a map first")
    labels = getattr(scan, "labels", {}) or {}
    return panel_from_arrays(scan.value, scan.x, scan.y,
                             labels.get("x", ""), labels.get("y", ""),
                             name or getattr(data, "source_label", "") or "panel",
                             info=getattr(scan, "info", {}),
                             sources=_dataset_sources(data))


def _dataset_sources(data) -> list:
    """The original file a dataset began as (its info's origin, else the
    file it was read from), for a panel made of it."""
    from tools.output.figfile import sources_from_info

    info = getattr(getattr(data, "scan", None), "info", None) or {}
    return sources_from_info(info, getattr(data, "path", "") or "")


class SliceSeriesDialog(QDialog):
    """A page of slices through a cube, which is what a Fermi-surface figure
    or an energy-dependence figure is made of.

    The MATLAB tool's ``massplott_3D``, with its three choices kept -- which
    axis to slice, where, and how far to integrate -- and its guesswork
    dropped: the positions can be typed, or taken as *n* evenly spaced
    values, and the titles are written from the axis's own label and unit
    rather than a hard-coded list of four units.
    """

    def __init__(self, parent, axes, labels, kind="map"):
        super().__init__(parent)
        self.setWindowTitle("Slice series figure")
        self.axes = axes            # (x, y, z) arrays
        self.labels = labels        # (x, y, z) labels
        layout = QVBoxLayout(self)
        intro = QLabel(
            "One panel per slice through the cube. The other two axes become "
            "each panel's picture.")
        intro.setWordWrap(True)
        layout.addWidget(intro)

        form = QFormLayout()
        self.direction = QComboBox()
        for index, label in enumerate(labels):
            self.direction.addItem(f"along {label}")
        self.direction.setCurrentIndex(2)
        self.direction.currentIndexChanged.connect(self._direction_changed)
        form.addRow("Slice", self.direction)

        self.mode_even = QRadioButton("evenly spaced")
        self.mode_typed = QRadioButton("at these values")
        self.mode_even.setChecked(True)
        mode_row = QWidget()
        mode_box = QHBoxLayout(mode_row)
        mode_box.setContentsMargins(0, 0, 0, 0)
        mode_box.addWidget(self.mode_even)
        mode_box.addWidget(self.mode_typed)
        mode_box.addStretch(1)
        form.addRow("Positions", mode_row)

        self.count_box = whole(6, 1, 64, "How many slices.")
        self.count_box.valueChanged.connect(self._refresh)
        form.addRow("Number", self.count_box)
        self.from_box = spin(0, -1e9, 1e9, 4, 0.05, "First slice.")
        self.to_box = spin(1, -1e9, 1e9, 4, 0.05, "Last slice.")
        self.from_box.valueChanged.connect(self._refresh)
        self.to_box.valueChanged.connect(self._refresh)
        form.addRow("From / to", FigureToolsDialog._pair(self.from_box, self.to_box))
        self.values_edit = QLineEdit()
        self.values_edit.setPlaceholderText("e.g. 0, -0.1, -0.25")
        self.values_edit.editingFinished.connect(self._refresh)
        form.addRow("Values", self.values_edit)

        self.width_box = spin(0, 0, 1e9, 4, 0.01,
                              "Integrate this far either side of each slice; "
                              "0 takes the single nearest plane.")
        self.width_box.valueChanged.connect(self._refresh)
        form.addRow("Integrate ±", self.width_box)

        self.cols_box = whole(3, 1, 12, "Columns in the figure.")
        form.addRow("Columns", self.cols_box)
        self.title_box = QCheckBox("Title each panel with its position")
        self.title_box.setChecked(True)
        form.addRow("", self.title_box)
        self.shared_box = QCheckBox("Same colour scale everywhere")
        self.shared_box.setChecked(True)
        form.addRow("", self.shared_box)
        layout.addLayout(form)

        self.note = QLabel("")
        self.note.setWordWrap(True)
        layout.addWidget(self.note)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.button(QDialogButtonBox.Ok).setText("Make the figure")
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)
        for radio in (self.mode_even, self.mode_typed):
            radio.toggled.connect(self._refresh)
        self._direction_changed()

    def _direction_changed(self, *_):
        axis = np.asarray(self.axes[self.direction.currentIndex()], dtype=float)
        self.from_box.setValue(float(axis.min()))
        self.to_box.setValue(float(axis.max()))
        self._refresh()

    def _refresh(self, *_):
        even = self.mode_even.isChecked()
        self.count_box.setEnabled(even)
        self.from_box.setEnabled(even)
        self.to_box.setEnabled(even)
        self.values_edit.setEnabled(not even)
        values = self.values()
        # The greyed half shows what the other one gives -- the positions
        # the even spacing works out to, or the number and range of the
        # typed ones -- never what was left there from before.
        if even:
            self.values_edit.setText(", ".join(f"{v:.6g}" for v in values))
        else:
            for box, value in ((self.count_box, len(values)),
                               (self.from_box, min(values)), (self.to_box, max(values))):
                box.blockSignals(True)
                box.setValue(value)
                box.blockSignals(False)
        axis = np.asarray(self.axes[self.direction.currentIndex()], dtype=float)
        step = abs(float(axis[1] - axis[0])) if axis.size > 1 else 0.0
        planes = max(int(round(2 * self.width_box.value() / step)) + 1, 1) if step else 1
        self.note.setText(
            f"{len(values)} panel(s); each sums {planes} plane(s) of the "
            f"{axis.size} along {self.labels[self.direction.currentIndex()]}.")

    def values(self):
        if self.mode_even.isChecked():
            return list(np.linspace(self.from_box.value(), self.to_box.value(),
                                    max(int(self.count_box.value()), 1)))
        out = []
        for chunk in self.values_edit.text().replace(";", ",").split(","):
            chunk = chunk.strip()
            if not chunk:
                continue
            try:
                out.append(float(chunk))
            except ValueError:
                continue
        return out or [float(np.mean(self.axes[self.direction.currentIndex()]))]

    def settings(self):
        return {"axis": self.direction.currentIndex(),
                "values": self.values(),
                "half_width": self.width_box.value(),
                "cols": int(self.cols_box.value()),
                "titles": self.title_box.isChecked(),
                "shared": self.shared_box.isChecked()}


def slice_panels(cube, axes, labels, *, axis: int, values, half_width: float = 0.0,
                 titles: bool = True, colormap="gray", flip=False):
    """Cut a cube into panels along one axis.

    ``cube`` is ``(x, y, z)``; slicing along one axis leaves the other two
    as the panel's picture, in the order they appear -- so a constant-energy
    slice of a map comes out as the two angles, and a slice along an angle
    comes out as the other angle against energy.
    """
    cube = np.asarray(cube, dtype=float)
    axes = [np.asarray(a, dtype=float) for a in axes]
    slicing = axes[axis]
    others = [i for i in range(3) if i != axis]
    unit = nf.unit_of(labels[axis])
    panels = []
    for value in values:
        centre = int(np.argmin(np.abs(slicing - value)))
        if half_width > 0:
            inside = np.abs(slicing - slicing[centre]) <= half_width
            indices = np.flatnonzero(inside)
            if indices.size == 0:
                indices = np.array([centre])
        else:
            indices = np.array([centre])
        block = np.take(cube, indices, axis=axis)
        frame = np.nanmean(block, axis=axis)
        panel = panel_from_arrays(
            frame, axes[others[0]], axes[others[1]],
            labels[others[0]], labels[others[1]],
            name=f"{labels[axis].split('(')[0].strip()} = {slicing[centre]:.4g}",
            colormap=colormap, flip=flip)
        if titles:
            panel.title = f"{slicing[centre]:.4g} {unit}".strip()
        panels.append(panel)
    return panels


# --------------------------------------------------------------------------
# The window
# --------------------------------------------------------------------------
class ShapeEditor(QGroupBox):
    """The selected drawing's colour, line, fill and text -- above the
    figure tools while a drawing is selected."""

    def __init__(self, window: "FigureWindow"):
        super().__init__("Selected drawing")
        self.window_ = window
        self.shape = None
        self._syncing = False
        form = QFormLayout(self)
        self.kind = QLabel("")
        form.addRow("Kind", self.kind)
        self.color = ColorButton("#000000")
        self.color.colorChanged.connect(lambda c: self._set("color", c))
        form.addRow("Colour", self.color)
        self.width = spin(0.8, 0.1, 12, 2, 0.1, "Line width (pt).")
        self.width.valueChanged.connect(lambda v: self._set("width_pt", v))
        form.addRow("Width (pt)", self.width)
        self.style = QComboBox()
        self.style.addItems(list(PEN_STYLES))
        self.style.currentTextChanged.connect(lambda v: self._set("style", v))
        form.addRow("Line", self.style)
        self.head = spin(1.6, 0.3, 10, 1, 0.1, "Size of the arrow head (mm).")
        self.head.valueChanged.connect(lambda v: self._set("head_mm", v))
        form.addRow("Head (mm)", self.head)
        self.double = QCheckBox("Heads at both ends")
        self.double.toggled.connect(lambda v: self._set("double", bool(v)))
        form.addRow("", self.double)
        fill_row = QWidget()
        box = QHBoxLayout(fill_row)
        box.setContentsMargins(0, 0, 0, 0)
        self.filled = QCheckBox("Fill")
        self.fill = ColorButton("#ffffff")
        self.filled.toggled.connect(self._set_fill)
        self.fill.colorChanged.connect(self._set_fill)
        box.addWidget(self.filled)
        box.addWidget(self.fill)
        box.addStretch(1)
        form.addRow("", fill_row)
        self.text = QLineEdit()
        self.text.setPlaceholderText("k_x, E_F, Å^-1 ...")
        self.text.editingFinished.connect(lambda: self._set("text", self.text.text()))
        form.addRow("Text", self.text)
        self.size = spin(0.0, 0, 60, 1, 0.5, "Text size (pt); 0 follows the figure.")
        self.size.valueChanged.connect(lambda v: self._set("size_pt", v))
        form.addRow("Size (pt)", self.size)
        self.bold = QCheckBox("Bold")
        self.bold.toggled.connect(lambda v: self._set("bold", bool(v)))
        form.addRow("", self.bold)
        delete = QPushButton("Delete")
        delete.clicked.connect(self.window_.canvas.delete_selection)
        form.addRow("", delete)
        self._rows = {"line": (self.width, self.style), "arrow": (self.head, self.double),
                      "fill": (fill_row,), "text": (self.text, self.size, self.bold)}
        self._form = form

    def show_shape(self, shape):
        self.shape = shape
        self.setVisible(shape is not None)
        if shape is None:
            return
        self._syncing = True
        try:
            self.kind.setText({"arrow": "Arrow", "line": "Line", "rect": "Box",
                               "ellipse": "Ellipse", "text": "Text"}[shape.kind])
            self.color.set_color(shape.color)
            self.width.setValue(shape.width_pt)
            self.style.setCurrentText(shape.style)
            self.head.setValue(shape.head_mm)
            self.double.setChecked(shape.double)
            self.filled.setChecked(bool(shape.fill))
            if shape.fill:
                self.fill.set_color(shape.fill)
            self.text.setText(shape.text)
            self.size.setValue(shape.size_pt)
            self.bold.setChecked(shape.bold)
        finally:
            self._syncing = False
        text = shape.kind == "text"
        visible = {"line": not text, "arrow": shape.kind == "arrow",
                   "fill": shape.kind in ("rect", "ellipse"), "text": text}
        for group, widgets in self._rows.items():
            for widget in widgets:
                widget.setVisible(visible[group])
                label = self._form.labelForField(widget)
                if label is not None:
                    label.setVisible(visible[group])

    def _set(self, name, value):
        if self._syncing or self.shape is None:
            return
        self.window_.remember("shape")
        setattr(self.shape, name, value)
        self.window_.refresh()

    def _set_fill(self, *_):
        self._set("fill", self.fill.color() if self.filled.isChecked() else "")


class FigureWindow(QMainWindow):
    """Where a figure for a paper is put together, like the page it will
    be printed on.

    It starts as a page of empty panels -- (a), (b), (c), (d) -- or with
    the panels it was opened with. Every panel can be dragged and resized
    on the page, a small inset put on top of one (a magnified region of it,
    or another dataset), and each slot filled from the main list or with a
    picture from a file. Arrows, lines, boxes and notes are drawn straight
    on the page with the tools on the left. The figure tools on the right
    set the axes, fonts, colours and the export of the selected panel.
    """

    closed = pyqtSignal(object)

    #: the drawing tools on the left: (tool, label, what it does)
    TOOL_BUTTONS = (
        ("select", "↖  Select / move",
         "Click a panel or a drawing to select it; drag it to move it, drag "
         "a handle to resize it. Del removes it, the arrow keys nudge it "
         "(Shift: 5 mm)."),
        ("slot", "▭  Panel",
         "Drag out a new empty panel. Drawn inside another panel, it becomes "
         "an inset (no letter)."),
        ("zoom", "⊕  Zoom-in inset",
         "Drag out a region of a panel: an inset magnifying it appears in "
         "the panel's corner, and the region is outlined."),
        ("arrow", "→  Arrow", "Drag from the tail to the head. Shift: "
         "horizontal, vertical or 45°."),
        ("line", "╱  Line", "Drag from one end to the other. Shift: "
         "horizontal, vertical or 45°."),
        ("rect", "□  Box", "Drag out a box round a feature."),
        ("ellipse", "○  Ellipse", "Drag out an ellipse."),
        ("text", "T  Text", "Click where the text goes. _ and ^ make sub- "
         "and superscripts: k_x, Å^-1."),
    )

    def __init__(self, panels, title="Figure", parent=None, dataset_source=None):
        super().__init__(parent)
        self.setWindowTitle(title)
        panels = list(panels)
        blank = not panels
        self.figure = Figure(panels, style=FigureStyle())
        if blank:
            self.figure.panels = [Panel(data=nf.empty_data()) for _ in range(4)]
            self.figure.fit_grid(2)
            self.figure.arrange_grid(2, 2)
        else:
            self.figure.fit_grid()
        self.figure.letter_panels()
        #: ``() -> [(label, key)]`` and ``key -> data``, so a panel can be
        #: taken from the launcher's list
        self.dataset_source = dataset_source
        self._undo, self._redo = [], []
        self._last_memory = (None, 0.0)

        central = QWidget()
        self.setCentralWidget(central)
        root = QVBoxLayout(central)
        root.setContentsMargins(4, 4, 4, 4)

        bar = QHBoxLayout()
        self.add_button = QPushButton("Add panel...")
        self.add_button.setToolTip("A new panel holding a dataset from the main list.")
        self.add_button.clicked.connect(self.add_panel_from_list)
        bar.addWidget(self.add_button)
        self.undo_button = QPushButton("↶ Undo")
        self.undo_button.clicked.connect(self.undo)
        self.redo_button = QPushButton("↷ Redo")
        self.redo_button.clicked.connect(self.redo)
        bar.addWidget(self._separator())
        bar.addWidget(self.undo_button)
        bar.addWidget(self.redo_button)
        bar.addWidget(self._separator())
        #: the .afig file this figure was opened from or last saved to
        self.file_path = ""
        for text, slot, tip in (
                ("Open...", lambda: self.open_file(),
                 "Open a figure saved by this program (.afig) in a new window."),
                ("Save", lambda: self.save_file(self.file_path or None),
                 "Save the figure to edit later (.afig), with the data it was "
                 "made from in a “data for …” folder beside it."),
                ("Save as...", lambda: self.save_file(),
                 "Save the figure under a new name (.afig), with its data "
                 "folder.")):
            button = QPushButton(text)
            button.setToolTip(tip)
            button.clicked.connect(slot)
            bar.addWidget(button)
        bar.addWidget(self._separator())
        for text, fmt in (("PDF", "pdf"), ("PNG", "png"), ("SVG", "svg")):
            button = QPushButton(text)
            button.setMaximumWidth(56)
            button.clicked.connect(lambda _c, f=fmt: self.export_figure(f))
            bar.addWidget(button)
        copy_button = QPushButton("Copy")
        copy_button.setMaximumWidth(56)
        copy_button.setToolTip("Put the figure on the clipboard.")
        copy_button.clicked.connect(self.copy_to_clipboard)
        bar.addWidget(copy_button)
        bar.addWidget(self._separator())
        bar.addWidget(QLabel("Zoom"))
        self.zoom_combo = QComboBox()
        self.zoom_combo.addItems(["Fit", "100 %", "150 %", "200 %", "300 %"])
        self.zoom_combo.setToolTip("100 % is the printed size on a 96-dpi screen.")
        self.zoom_combo.currentIndexChanged.connect(self._zoom_changed)
        bar.addWidget(self.zoom_combo)
        self.snap_box = QCheckBox("Snap")
        self.snap_box.setChecked(True)
        self.snap_box.setToolTip("Panels snap to each other's edges and to a "
                                 "0.5 mm grid while they are dragged.")
        bar.addWidget(self.snap_box)
        bar.addStretch(1)
        self.tools_button = QPushButton("Figure tools")
        self.tools_button.setCheckable(True)
        self.tools_button.setChecked(True)
        self.tools_button.setToolTip("Show or hide the figure tools on the right.")
        bar.addWidget(self.tools_button)
        root.addLayout(bar)

        body = QSplitter(Qt.Horizontal)
        body.addWidget(self._tool_palette())
        self.canvas = FigureCanvas(self.figure, self)
        self.painter = self.canvas.painter
        self.snap_box.toggled.connect(lambda on: setattr(self.canvas, "snap", bool(on)))
        scroller = QScrollArea()
        scroller.setWidgetResizable(True)
        scroller.setWidget(self.canvas)
        self.scroller = scroller
        body.addWidget(scroller)

        side = QWidget()
        side_layout = QVBoxLayout(side)
        side_layout.setContentsMargins(0, 0, 0, 0)
        self.shape_editor = ShapeEditor(self)
        self.shape_editor.hide()
        side_layout.addWidget(self.shape_editor)
        self._tools = FigureToolsDialog(self, embedded=True)
        side_layout.addWidget(self._tools, 1)
        self.side = side
        body.addWidget(side)
        body.setStretchFactor(0, 0)
        body.setStretchFactor(1, 1)
        body.setStretchFactor(2, 0)
        body.setSizes([170, 700, 440])
        self.tools_button.toggled.connect(side.setVisible)
        root.addWidget(body, 1)

        self.hint = QLabel("")
        self.hint.setWordWrap(True)
        root.addWidget(self.hint)

        self.canvas.panelPicked.connect(self._panel_picked)
        self.canvas.selectionChanged.connect(self._selection_changed)
        self.canvas.aboutToChange.connect(lambda: self.remember("canvas", force=True))
        self.canvas.changed.connect(self._canvas_changed)
        self.canvas.toolFinished.connect(self._tool_finished)
        self.canvas.panelActivated.connect(self._panel_activated)
        self.canvas.menuRequested.connect(self._context_menu)
        for keys, slot in ((QKeySequence.Undo, self.undo), (QKeySequence.Redo, self.redo),
                           ("Ctrl+Shift+Z", self.redo)):
            shortcut = QShortcut(QKeySequence(keys), self)
            shortcut.activated.connect(slot)

        self.statusBar()
        self.resize(1320, 820)
        self._update_undo()
        self.set_tool("select")
        self.refresh()

    @staticmethod
    def _separator():
        line = QFrame()
        line.setFrameShape(QFrame.VLine)
        line.setFrameShadow(QFrame.Sunken)
        return line

    def _tool_palette(self) -> QWidget:
        palette = QWidget()
        column = QVBoxLayout(palette)
        column.setContentsMargins(2, 2, 2, 2)
        heading = QLabel("<b>Draw</b>")
        column.addWidget(heading)
        self.tool_buttons = {}
        for tool, text, tip in self.TOOL_BUTTONS:
            button = QPushButton(text)
            button.setCheckable(True)
            button.setToolTip(tip)
            button.setStyleSheet("text-align: left; padding: 4px 8px;")
            button.clicked.connect(lambda _c, t=tool: self.set_tool(t))
            column.addWidget(button)
            self.tool_buttons[tool] = button
        column.addSpacing(8)
        column.addWidget(QLabel("<b>Selected panel</b>"))
        for text, slot, tip in (
                ("Load from list...", lambda: self.load_from_list(),
                 "Put a dataset from the main list into the selected panel."),
                ("Load picture...", lambda: self.load_picture(),
                 "Put a picture from a file (PNG, JPEG, TIFF...) into the "
                 "selected panel."),
                ("Empty it", lambda: self.empty_panel(),
                 "Take the content out; the panel stays as a slot."),
                ("Delete", lambda: self.canvas.delete_selection(),
                 "Remove the selected panel or drawing (Del).")):
            button = QPushButton(text)
            button.setToolTip(tip)
            button.setStyleSheet("text-align: left; padding: 4px 8px;")
            button.clicked.connect(slot)
            column.addWidget(button)
        column.addSpacing(8)
        column.addWidget(QLabel("<b>Page</b>"))
        letter = QPushButton("Re-letter (a), (b)...")
        letter.setToolTip("Letter the panels again in reading order: top to "
                          "bottom, left to right. Insets take no letter.")
        letter.setStyleSheet("text-align: left; padding: 4px 8px;")
        letter.clicked.connect(self.reletter)
        column.addWidget(letter)
        column.addSpacing(8)
        column.addWidget(QLabel("<b>Arrange</b> <i>(Ctrl+click: several)</i>"))
        grid = QGridLayout()
        grid.setSpacing(2)
        for n, (text, operation, tip) in enumerate((
                ("⇤ Left", "left", "Line up the selected panels' left edges."),
                ("Right ⇥", "right", "Line up their right edges."),
                ("⤒ Top", "top", "Line up their top edges."),
                ("Bottom ⤓", "bottom", "Line up their bottom edges."),
                ("↔ Centres", "hcenter", "Their centres on one vertical line "
                                          "(the first-selected panel's)."),
                ("↕ Centres", "vcenter", "Their centres on one horizontal line."),
                ("Space ↔", "space-x", "Equal gaps between them, left to right "
                                       "(three or more)."),
                ("Space ↕", "space-y", "Equal gaps between them, top to bottom."),
                ("Same width", "width", "The first-selected panel's width for all."),
                ("Same height", "height", "Its height for all."))):
            button = QPushButton(text)
            button.setToolTip(tip)
            button.setStyleSheet("padding: 3px 4px;")
            button.clicked.connect(lambda _c, o=operation: self.arrange(o))
            grid.addWidget(button, n // 2, n % 2)
        column.addLayout(grid)
        column.addStretch(1)
        palette.setMinimumWidth(150)
        palette.setMaximumWidth(220)
        return palette

    # -- tools ---------------------------------------------------------------
    def set_tool(self, tool: str):
        self.canvas.set_tool(tool)
        for name, button in self.tool_buttons.items():
            button.setChecked(name == tool)
        tip = next(t for name, _x, t in self.TOOL_BUTTONS if name == tool)
        self.hint.setText(tip if tool != "select" else
                          "Click a panel or a drawing to select it, drag to move, "
                          "drag a handle to resize. Double-click an empty panel to "
                          "fill it. Right-click for more.")

    def _tool_finished(self):
        self.set_tool(self.canvas.tool)

    def _zoom_changed(self, index):
        self.canvas.zoom = None if index == 0 else (
            (100, 150, 200, 300)[index - 1] / 100.0 * 96.0 / nf.MM_PER_INCH)
        self.canvas.refresh()

    # -- state -------------------------------------------------------------
    def refresh(self):
        self.canvas.refresh()
        width, height = self.painter.size_mm()
        mains = len(self.figure.main_panels())
        insets = len(self.figure.panels) - mains
        self.statusBar().showMessage(
            f"{mains} panel(s)" + (f", {insets} inset(s)" if insets else "")
            + (f", {len(self.figure.shapes)} drawing(s)" if self.figure.shapes else "")
            + f" — {width:.1f} x {height:.1f} mm at {self.figure.style.dpi} dpi")

    def _canvas_changed(self):
        self.refresh()
        if self._tools is not None:
            self._tools.sync()

    def _panel_picked(self, index):
        if self._tools is not None:
            self._tools.panel_combo.setCurrentIndex(index)

    def _selection_changed(self):
        self.shape_editor.show_shape(self.canvas.selected_shape())

    def open_tools(self):
        """Show the figure tools (they live on the right of the window)."""
        self.tools_button.setChecked(True)
        self.side.setVisible(True)
        return self._tools

    def pick_point(self, message, callback):
        self.hint.setText(message)
        self.canvas.start_point_picking(
            lambda index, x, y: (self.hint.setText(""), callback(index, x, y)))

    # -- undo ----------------------------------------------------------------
    def _snapshot(self):
        figure = self.figure
        return {"panels": [(p, {k: copy.deepcopy(v) for k, v in vars(p).items()
                                if k != "data"}, p.data) for p in figure.panels],
                "shapes": copy.deepcopy(figure.shapes),
                "style": copy.deepcopy(figure.style),
                "grid": (figure.rows, figure.cols), "title": figure.title}

    def _restore(self, state):
        figure = self.figure
        figure.panels = []
        for panel, fields, data in state["panels"]:
            for name, value in fields.items():
                setattr(panel, name, copy.deepcopy(value))
            panel.data = data
            figure.panels.append(panel)
        figure.shapes = copy.deepcopy(state["shapes"])
        figure.style = copy.deepcopy(state["style"])
        figure.rows, figure.cols = state["grid"]
        figure.title = state["title"]
        self.canvas.selection = None
        self._selection_changed()
        if self._tools is not None:
            self._tools.sync()
        self.refresh()

    def remember(self, source="property", force=False):
        """Keep an undo step before a change. A burst of changes from one
        control (a spin box held down) is one step."""
        last_source, when = self._last_memory
        now = time.monotonic()
        if not force and source == last_source and now - when < 1.0:
            self._last_memory = (source, now)
            return
        self._last_memory = (source, now)
        self._undo.append(self._snapshot())
        del self._undo[:-60]
        self._redo.clear()
        self._update_undo()

    def undo(self):
        if not self._undo:
            return False
        self._redo.append(self._snapshot())
        self._restore(self._undo.pop())
        self._last_memory = (None, 0.0)
        self._update_undo()
        return True

    def redo(self):
        if not self._redo:
            return False
        self._undo.append(self._snapshot())
        self._restore(self._redo.pop())
        self._last_memory = (None, 0.0)
        self._update_undo()
        return True

    def _update_undo(self):
        self.undo_button.setEnabled(bool(self._undo))
        self.redo_button.setEnabled(bool(self._redo))

    # -- panels ------------------------------------------------------------
    def add_panels(self, panels, cols=None):
        if self.figure.free:
            self.remember("add", force=True)
            for panel in panels:
                self.figure.place_new(panel)
                self.figure.panels.append(panel)
            self.figure.auto_edge_labels()
        else:
            self.figure.panels.extend(panels)
            self.figure.fit_grid(cols or self.figure.cols)
        self.figure.letter_panels()
        if self._tools is not None:
            self._tools.sync()
        self.refresh()

    def _choose_dataset(self, title):
        """``(data, name)`` picked from the launcher's list, or None."""
        source = self.dataset_source
        if not source or not callable(source.get("entries")):
            QMessageBox.information(self, title, "No dataset list is available here.")
            return None
        listed = source["entries"]()
        if not listed:
            QMessageBox.information(self, title, "The list is empty.")
            return None
        labels = [label for label, _key in listed]
        choice, ok = QInputDialog.getItem(self, title, "Dataset:", labels, 0, False)
        if not ok:
            return None
        key = listed[labels.index(choice)][1]
        try:
            return source["loader"](key), choice.split("   [")[0]
        except Exception as exc:
            QMessageBox.warning(self, title, f"Could not read it:\n{exc}")
            return None

    def add_panel_from_list(self):
        """Take another dataset from the launcher's list, as a new panel."""
        chosen = self._choose_dataset("Add panel")
        if chosen is None:
            return None
        try:
            panel = panel_from_dataset(*chosen)
        except Exception as exc:
            QMessageBox.warning(self, "Add panel",
                                f"That dataset cannot be a panel:\n{exc}")
            return None
        self.add_panels([panel])
        return panel

    def _target_panel(self, what):
        panel = self.canvas.selected_panel()
        if panel is None:
            QMessageBox.information(self, what, "Click a panel on the page first.")
        return panel

    def fill_panel(self, panel, source: Panel):
        """Put ``source``'s content (its data and how it is drawn) into
        ``panel``, which keeps its place, size, letter and role."""
        self.remember("fill", force=True)
        keep = {name: getattr(panel, name) for name in
                ("frame", "role", "label", "uid", "x_labels", "y_labels", "font_pt")}
        for name, value in vars(source).items():
            setattr(panel, name, value)
        for name, value in keep.items():
            setattr(panel, name, value)
        panel.zoom_of = ""
        # Insets that magnified the old content no longer show it.
        self.figure.panels = [p for p in self.figure.panels if p.zoom_of != panel.uid]
        if self._tools is not None:
            self._tools.sync()
        self.refresh()
        return panel

    def load_from_list(self, panel=None):
        """Fill the selected panel with a dataset from the main list."""
        panel = panel or self._target_panel("Load from list")
        if panel is None:
            return None
        chosen = self._choose_dataset("Load from list")
        if chosen is None:
            return None
        try:
            source = panel_from_dataset(*chosen)
        except Exception as exc:
            QMessageBox.warning(self, "Load from list",
                                f"That dataset cannot be a panel:\n{exc}")
            return None
        return self.fill_panel(panel, source)

    PICTURE_FILTER = "Pictures (*.png *.jpg *.jpeg *.tif *.tiff *.bmp *.gif *.svg)"

    def load_picture(self, panel=None, path=None):
        """Fill the selected panel with a picture from a file; the panel
        takes the picture's shape."""
        panel = panel or self._target_panel("Load picture")
        if panel is None:
            return None
        if path is None:
            path, _ = QFileDialog.getOpenFileName(self, "Load a picture", "",
                                                  self.PICTURE_FILTER)
            if not path:
                return None
        image = QImage(path)
        if image.isNull():
            QMessageBox.warning(self, "Load picture", f"Could not read {path}.")
            return None
        rgba = _qimage_array(image)
        source = Panel(data=nf.picture_data(rgba, os.path.basename(path)))
        source.data.sources = [figfile.source_record(path, kind="picture")]
        self.figure.freeze()
        self.fill_panel(panel, source)
        x, y, w, h = panel.frame
        height, width = rgba.shape[:2]
        panel.frame = (x, y, w, w * height / max(width, 1))
        panel.x_labels = panel.y_labels = False
        self.refresh()
        return panel

    def empty_panel(self, panel=None):
        panel = panel or self._target_panel("Empty it")
        if panel is None:
            return None
        return self.fill_panel(panel, Panel(data=nf.empty_data()))

    def reletter(self):
        self.remember("letter", force=True)
        self.figure.freeze()
        self.figure.letter_panels()
        if self._tools is not None:
            self._tools.sync()
        self.refresh()

    def _panel_activated(self, index):
        panel = self.figure.panels[index]
        if panel.data.kind != "empty":
            self.open_tools()
            return
        menu = QMenu(self)
        menu.addAction("Load from the list...", lambda: self.load_from_list(panel))
        menu.addAction("Load a picture...", lambda: self.load_picture(panel))
        menu.exec_(QCursor.pos())

    def _context_menu(self, what, index, position):
        menu = QMenu(self)
        if what == "panel":
            panel = self.figure.panels[index]
            menu.addAction("Load from the list...", lambda: self.load_from_list(panel))
            menu.addAction("Load a picture...", lambda: self.load_picture(panel))
            if panel.data.kind != "empty":
                menu.addAction("Empty it", lambda: self.empty_panel(panel))
            if panel.data.is_plot() and panel.role != "inset":
                menu.addAction("Zoom-in inset of a region...",
                               lambda: self.set_tool("zoom"))
            menu.addSeparator()
            role = "Make it a panel (lettered)" if panel.role == "inset" \
                else "Make it an inset (no letter)"
            menu.addAction(role, lambda: self._toggle_role(panel))
            menu.addAction("Bring to front", lambda: self._restack(panel, +1))
            menu.addAction("Send to back", lambda: self._restack(panel, -1))
            menu.addSeparator()
            menu.addAction("Delete", self.canvas.delete_selection)
        elif what == "shape":
            menu.addAction("Bring to front", lambda: self._restack_shape(index, +1))
            menu.addAction("Send to back", lambda: self._restack_shape(index, -1))
            menu.addAction("Delete", self.canvas.delete_selection)
        else:
            menu.addAction("Add an empty panel here",
                           lambda: self.set_tool("slot"))
            menu.addAction("Arrange as grid...", self.open_tools)
        menu.exec_(position)

    def _toggle_role(self, panel):
        self.remember("role", force=True)
        self.figure.freeze()
        panel.role = "panel" if panel.role == "inset" else "inset"
        if panel.role == "panel":
            panel.zoom_of = ""
        self.figure.letter_panels()
        self._canvas_changed()

    def _restack(self, panel, direction):
        self.remember("order", force=True)
        panels = self.figure.panels
        panels.remove(panel)
        if direction > 0:
            panels.append(panel)
        else:
            panels.insert(0, panel)
        self.canvas.selection = ("panel", panels.index(panel))
        self._canvas_changed()

    def _restack_shape(self, index, direction):
        self.remember("order", force=True)
        shapes = self.figure.shapes
        shape = shapes.pop(index)
        if direction > 0:
            shapes.append(shape)
        else:
            shapes.insert(0, shape)
        self.canvas.selection = ("shape", shapes.index(shape))
        self._canvas_changed()

    # -- output ------------------------------------------------------------
    def _default_name(self, suffix):
        stem = "".join(ch for ch in (self.windowTitle() or "figure")
                       if ch.isalnum() or ch in "._- ") or "figure"
        return stem.strip() + suffix

    def export_figure(self, fmt: str):
        filters = {"pdf": "PDF (*.pdf)", "svg": "SVG (*.svg)",
                   "png": "PNG (*.png)", "tif": "TIFF (*.tif)",
                   "layers": "PNG (*.png)"}
        suffix = {"pdf": ".pdf", "svg": ".svg", "png": ".png",
                  "tif": ".tif", "layers": ".png"}[fmt]
        path, _ = QFileDialog.getSaveFileName(
            self, "Save the figure", self._default_name(suffix), filters[fmt])
        if not path:
            return None
        if not path.lower().endswith(suffix):
            path += suffix
        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            written = self.write(path, fmt)
        except Exception as exc:
            QMessageBox.warning(self, "Save the figure",
                                 f"Could not write it:\n{exc}")
            return None
        finally:
            QApplication.restoreOverrideCursor()
        self.statusBar().showMessage("Wrote " + ", ".join(
            os.path.basename(p) for p in written))
        return written

    def arrange(self, operation: str, indices=None):
        """Align, space or size the selected panels (or ``indices``)."""
        indices = list(indices if indices is not None else self.canvas.selected_indices())
        if len(indices) < 2:
            self.statusBar().showMessage("Select two or more panels: click one, "
                                         "then Ctrl+click the others")
            return False
        self.remember("arrange", force=True)
        if operation in ("left", "right", "top", "bottom", "hcenter", "vcenter"):
            self.figure.align_panels(indices, operation)
        elif operation in ("space-x", "space-y"):
            self.figure.distribute_panels(indices, operation[-1])
        elif operation in ("width", "height"):
            self.figure.match_size(indices, operation)
        else:
            raise ValueError(f"unknown arrangement {operation!r}")
        self.refresh()
        return True

    # -- the program's own figure file ---------------------------------------
    def save_file(self, path=None, *, with_data: bool = True):
        """Save the figure as an .afig file to open and edit again, and
        (``with_data``) the data each panel was made from in "data for
        <name>" beside it. Returns the path, or None."""
        if path is None:
            start = self.file_path or self._default_name(figfile.EXTENSION)
            path, _ = QFileDialog.getSaveFileName(
                self, "Save the figure", start,
                f"ARPES viewer figure (*{figfile.EXTENSION})")
            if not path:
                return None
        if not path.lower().endswith(figfile.EXTENSION):
            path += figfile.EXTENSION
        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            from version import VERSION
            path = figfile.save(path, self.figure, flatten=self._flatten_panel,
                                app_version=VERSION)
            copied = (figfile.write_data_folder(path, self.figure,
                                                picture_writer=_write_png)
                      if with_data else {})
        except Exception as exc:
            QMessageBox.warning(self, "Save the figure", f"Could not save it:\n{exc}")
            return None
        finally:
            QApplication.restoreOverrideCursor()
        self.file_path = path
        # the next save finds the originals in this one's data folder
        self.figure.loaded_from = path
        self.setWindowTitle(os.path.splitext(os.path.basename(path))[0])
        files = sum(len(v) for v in copied.values())
        self.statusBar().showMessage(
            f"Saved {os.path.basename(path)}"
            + (f"; {files} file(s) in \u201c{os.path.basename(figfile.data_folder_of(path))}\u201d"
               if with_data else ""))
        return path

    def open_file(self, path=None):
        """Open an .afig figure in a new window (this one stays)."""
        if path is None:
            path, _ = QFileDialog.getOpenFileName(
                self, "Open a figure", os.path.dirname(self.file_path or ""),
                f"ARPES viewer figure (*{figfile.EXTENSION})")
            if not path:
                return None
        try:
            window = open_figure_file(path, dataset_source=self.dataset_source)
        except Exception as exc:
            QMessageBox.warning(self, "Open a figure", f"Could not open it:\n{exc}")
            return None
        keep_figure_window(window)
        return window

    def set_figure(self, figure):
        """Show ``figure`` (a loaded one) in this window: its panels, page,
        drawings and style replace what is here."""
        self.figure.panels = list(figure.panels)
        self.figure.shapes = list(figure.shapes)
        self.figure.style = figure.style
        self.figure.rows, self.figure.cols = figure.rows, figure.cols
        self.figure.title = figure.title
        self.figure.loaded_from = getattr(figure, "loaded_from", "")
        self._undo.clear()
        self._redo.clear()
        self._update_undo()
        self.canvas.selection = None
        self._selection_changed()
        if self._tools is not None:
            self._tools.sync()
        self.refresh()

    def _flatten_panel(self, panel):
        """A 3-D view panel as the picture it shows (an .afig keeps
        pictures and arrays, not live 3-D scenes)."""
        x, y, w, h = panel.frame or (0.0, 0.0, 60.0, 60.0)
        single = copy.copy(panel)
        single.label, single.colorbar, single.frame = "", False, (0.0, 0.0, w, h)
        page = Figure([single], style=copy.deepcopy(self.figure.style))
        page.style.margin_left_mm = page.style.margin_right_mm = 0.0
        page.style.margin_top_mm = page.style.margin_bottom_mm = 0.0
        page.style.width_mm, page.style.height_mm = w, h
        return _qimage_array(FigurePainter(page).to_image())

    def write(self, path: str, fmt: str):
        """Write the figure. Returns every file written."""
        if fmt == "layers":
            return self._write_layers(path)
        return write_figure(self.file_painter(), path, fmt)

    def file_painter(self) -> FigurePainter:
        """A painter for files and the clipboard: the canvas's own draws
        the empty slots, which never belong in the figure."""
        return FigurePainter(self.figure)

    def _write_layers(self, path: str):
        """The two-file route: the data as a raster and the axes as vector
        PostScript, for a figure that will be rearranged in Illustrator.

        One panel only -- a page of panels is what the PDF is for.
        """
        panel = self.figure.panels[self.canvas.selected]
        if not panel.data.is_plot():
            raise ValueError("only a panel of data has a separate data layer "
                             "and axes; save this one as PDF, SVG or PNG")
        (x0, x1), (y0, y1) = self.figure.panel_ranges(panel)
        _w, _h, rects = self.figure.layout()
        rect = next(r for r in rects if r["index"] == self.canvas.selected)
        dpi_scale = self.figure.style.dpi / nf.MM_PER_INCH
        rgba = self.painter._raster(
            panel, QRectF(0, 0, rect["w"], rect["h"]), (x0, x1), (y0, y1),
            self.figure.panel_levels(panel), 1.0)
        export.write_image(path, rgba, pil_format="PNG", keep_alpha=True,
                               dpi=self.figure.style.dpi)
        eps = os.path.splitext(path)[0] + "_axes.eps"
        export.axes_eps(
            eps, x_range=(x0, x1), y_range=(y0, y1),
            x_label=nf.plain_text(panel.label_x()),
            y_label=nf.plain_text(panel.label_y()),
            width_mm=rect["w"], height_mm=rect["h"],
            font_size=self.figure.style.label_pt or self.figure.style.font_pt,
            line_width=self.figure.style.axis_width_pt,
            minor_ticks=self.figure.style.minor_ticks,
            ticks_inward=self.figure.style.tick_direction == "in",
            invert_x=panel.invert_x, invert_y=panel.invert_y)
        return [path, eps]

    def copy_to_clipboard(self):
        QApplication.clipboard().setImage(self.file_painter().to_image())
        self.statusBar().showMessage("Figure copied to the clipboard")

    def closeEvent(self, event):
        if self._tools is not None:
            self._tools.close()
        self.closed.emit(self)
        super().closeEvent(event)


def write_figure(painter: FigurePainter, path: str, fmt: str):
    """Write the figure ``painter`` draws as PDF, SVG, PNG or TIFF -- the
    composer's own export, for a window that saves a figure without
    opening one (a 3-D view's "Save the view..."). Returns the files."""
    if fmt == "pdf":
        return [painter.to_pdf(path)]
    if fmt == "svg":
        return [painter.to_svg(path)]
    if fmt in ("png", "tif"):
        image = painter.to_image()
        export.write_image(path, _qimage_array(image),
                           pil_format="PNG" if fmt == "png" else "TIFF",
                           keep_alpha=True, dpi=painter.figure.style.dpi)
        return [path]
    raise ValueError(f"unknown format {fmt!r}")


def _write_png(rgba, path):
    export.write_image(path, np.asarray(rgba), pil_format="PNG", keep_alpha=True)


def open_figure_file(path: str, dataset_source=None, show: bool = True):
    """A composer window with the figure saved in ``path`` (.afig), ready
    to edit."""
    figure = figfile.load(path)
    window = FigureWindow([Panel(data=nf.empty_data())],
                          os.path.splitext(os.path.basename(path))[0],
                          dataset_source=dataset_source)
    window.set_figure(figure)
    window.file_path = os.path.abspath(path)
    if show:
        window.show()
    return window


def keep_figure_window(window):
    from ui.viewers.windows import keep_alive

    keep_alive(window)
    window.show()
    window.raise_()
    return window


def open_figure_windows():
    """The figure composers open now, newest last -- for a window that
    offers to add its picture to one of them."""
    return [w for w in QApplication.topLevelWidgets()
            if isinstance(w, FigureWindow) and w.isVisible()]


def _qimage_array(image: QImage) -> np.ndarray:
    image = image.convertToFormat(QImage.Format_RGBA8888)
    size = (image.sizeInBytes() if hasattr(image, "sizeInBytes")
            else image.byteCount())
    buffer = image.constBits()
    buffer.setsize(int(size))
    array = np.frombuffer(buffer, dtype=np.uint8).reshape(
        image.height(), image.bytesPerLine() // 4, 4)
    return np.array(array[:, :image.width(), :], copy=True)
