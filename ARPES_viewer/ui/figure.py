"""
ui/figure.py
===================
The figure composer: the window "Open in a new panel" leads to, the plot
tools that dress a figure up, and the painter that draws it.

One painter, every output
-------------------------
:class:`FigurePainter` draws a :class:`nf.Figure` onto any Qt paint
device. The preview on screen, the PNG, the TIFF, the vector PDF and the
SVG are the same call with a different device, so a figure cannot look
right in the window and wrong in the file. Text stays text in the PDF and
the SVG; the data itself is a raster embedded at the figure's own dpi,
which is what a journal wants -- an ARPES panel is an image, and pretending
otherwise makes a 40 MB vector file of a million little rectangles.

Everything is laid out in millimetres by ``tools.figure`` and multiplied by a
single scale factor here, so "89 mm wide" means 89 mm wide in the PDF and
the same figure at 600 dpi is 2102 px.

What the tools act on
---------------------
:class:`PlotToolsDialog` edits the *selected panel* of the figure, or every
panel at once when "Apply to all panels" is ticked -- which is the honest
form of the MATLAB original's Apply2All, where each button had to
re-implement the loop. Its tabs follow what a figure is made of rather than
what Qt calls things: Panels, Axes, Text, Colour, Overlays, Export.
"""
from __future__ import annotations

import os

import numpy as np
from PyQt5.QtCore import Qt, QRectF, QPointF, QSize, pyqtSignal
from PyQt5.QtGui import (QColor, QFont, QFontMetricsF, QImage, QPainter,
                         QPainterPath, QPen, QBrush, QPolygonF, QPdfWriter,
                         QPageSize)
from PyQt5.QtWidgets import (QWidget, QMainWindow, QVBoxLayout, QHBoxLayout,
                             QLabel, QPushButton, QComboBox, QCheckBox,
                             QDoubleSpinBox, QSpinBox, QLineEdit, QFormLayout,
                             QTabWidget, QDialog, QDialogButtonBox, QGroupBox,
                             QFileDialog, QMessageBox, QListWidget, QScrollArea,
                             QFrame, QApplication, QRadioButton, QInputDialog,
                             QSizePolicy, QColorDialog, QListWidgetItem)

from tools import colormaps
from tools import export
import tools.figure as nf
from tools.figure import (Figure, FigureStyle, Panel, PanelData, Line, TextNote,
                        Arrow, ScaleBar, Inset, Overlay, JOURNAL_PRESETS,
                        parse_markup, pt_to_mm)

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
        painter.setRenderHint(QPainter.Antialiasing, True)
        painter.setRenderHint(QPainter.TextAntialiasing, True)
        painter.setRenderHint(QPainter.SmoothPixmapTransform, True)
        painter.fillRect(QRectF(0, 0, width * scale, height * scale),
                         QColor(style.background))

        for rect in rects:
            self._paint_panel(painter, figure.panels[rect["index"]], rect,
                              scale, rect["index"])
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
        figure, style = self.figure, self.figure.style
        (x0, x1), (y0, y1) = figure.panel_ranges(panel)
        levels = figure.panel_levels(panel)
        box = QRectF(rect["x"] * scale, rect["y"] * scale,
                     rect["w"] * scale, rect["h"] * scale)

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
                                      for i in range(minor)]) if x_ticks.size else x_ticks
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
            for value in x_ticks:
                self._draw_markup(painter, px(value), baseline, x_fmt % value,
                                  tick_pt, scale, style.axis_color, anchor="center")
            self._draw_markup(painter, box.center().x(),
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
            if len(overlay.xs) < 2:
                continue
            painter.setPen(self._pen(overlay.color, overlay.width_pt, scale,
                                     overlay.style))
            path = QPainterPath()
            first = self._map_point(box, panel, x_range, y_range,
                                    overlay.xs[0], overlay.ys[0])
            path.moveTo(first)
            for x, y in zip(overlay.xs[1:], overlay.ys[1:]):
                path.lineTo(self._map_point(box, panel, x_range, y_range, x, y))
            if overlay.closed:
                path.closeSubpath()
            painter.setBrush(Qt.NoBrush)
            painter.drawPath(path)
        for arrow in panel.arrows:
            self._draw_arrow(painter, panel, box, x_range, y_range, arrow, scale)
        for note in panel.texts:
            size = note.size_pt or style.annotation_pt or style.font_pt
            point = self._map_point(box, panel, x_range, y_range, note.x, note.y)
            self._draw_markup(painter, point.x(), point.y(), note.text, size,
                              scale, note.color, anchor=note.anchor,
                              bold=note.bold)

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
        tick_pt = style.tick_pt or style.font_pt
        left = box.right() + style.colorbar_gap_mm * scale
        width = style.colorbar_mm * scale
        strip = QRectF(left, box.top(), width, box.height())
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
        for value, y in ((hi, strip.top()), (lo, strip.bottom())):
            text = f"{value:.4g}"
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
        color = style.panel_label_color if inside else style.axis_color
        self._draw_markup(painter, x, y, panel.label, size, scale, color,
                          anchor=anchor, bold=style.panel_label_bold)

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
    """Shows the rendered figure and lets a panel be picked by clicking it.

    Nothing here is a pyqtgraph plot: the figure is a picture drawn by
    :class:`FigurePainter`, which is the only way the preview and the file
    can be guaranteed to agree.
    """

    panelPicked = pyqtSignal(int)
    pointPicked = pyqtSignal(int, float, float)

    def __init__(self, figure: Figure, parent=None):
        super().__init__(parent)
        self.figure = figure
        self.painter = FigurePainter(figure)
        self.selected = 0
        self._image = None
        self._rect = None          # where the image sits in the widget
        self._picking = None       # a callback while a point is being picked
        self.setMinimumSize(280, 220)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.setCursor(Qt.ArrowCursor)

    # -- rendering ---------------------------------------------------------
    def refresh(self):
        self._image = None
        self.update()

    def _ensure_image(self):
        if self._image is not None:
            return
        width_mm, height_mm = self.painter.size_mm()
        margin = 12
        available_w = max(self.width() - 2 * margin, 60)
        available_h = max(self.height() - 2 * margin, 60)
        scale = min(available_w / max(width_mm, 1.0),
                    available_h / max(height_mm, 1.0))
        self._image = self.painter.preview(max(int(width_mm * scale), 60),
                                           self.selected)
        x = (self.width() - self._image.width()) // 2
        y = (self.height() - self._image.height()) // 2
        self._rect = (x, y, self._image.width(), self._image.height())

    def paintEvent(self, event):
        self._ensure_image()
        painter = QPainter(self)
        painter.fillRect(self.rect(), QColor("#e9e6ef"))
        if self._image is not None:
            painter.drawImage(self._rect[0], self._rect[1], self._image)
        painter.end()

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self._image = None

    # -- picking -----------------------------------------------------------
    def start_point_picking(self, callback):
        """The next click reports ``(panel index, x, y)`` in data
        coordinates -- how an annotation gets placed where it is wanted
        rather than typed in blind."""
        self._picking = callback
        self.setCursor(Qt.CrossCursor)

    def stop_point_picking(self):
        self._picking = None
        self.setCursor(Qt.ArrowCursor)

    def _locate(self, pos):
        """Which panel was clicked, and where in its data."""
        if self._rect is None:
            return None
        width_mm, height_mm = self.painter.size_mm()
        scale = self._rect[2] / max(width_mm, 1.0)
        mm_x = (pos.x() - self._rect[0]) / scale
        mm_y = (pos.y() - self._rect[1]) / scale
        _w, _h, rects = self.figure.layout()
        for rect in rects:
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

    def mousePressEvent(self, event):
        found = self._locate(event.pos())
        if found is None:
            return
        index, x, y = found
        if self._picking is not None:
            callback, self._picking = self._picking, None
            self.setCursor(Qt.ArrowCursor)
            callback(index, x, y)
            return
        self.selected = index
        self.refresh()
        self.panelPicked.emit(index)


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
class PlotToolsDialog(QDialog):
    """Everything that dresses a figure up, in six tabs.

    Each control writes into the model and asks the canvas to redraw, so
    there is no "apply" step and nothing can be set in the dialog but
    missing from the file. "Apply to all panels" makes a control act on
    every panel at once -- the MATLAB original had to re-implement that loop
    inside each of its forty callbacks.
    """

    def __init__(self, window: "FigureWindow"):
        super().__init__(window)
        self.window_ = window
        self.figure = window.figure
        self.setWindowTitle("Plot tools")
        self.setModal(False)
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
        close = QDialogButtonBox(QDialogButtonBox.Close)
        close.rejected.connect(self.reject)
        outer.addWidget(close)
        self.setMinimumWidth(470)
        self.resize(560, 700)
        self.sync()

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
        for panel in self.targets():
            setter(panel)
        self.redraw()

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
                name = panel.data.name or f"panel {i + 1}"
                self.panel_combo.addItem(f"{i + 1}. {name}")
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
        form = QFormLayout(page)

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

        self.width_box = spin(style.width_mm, 20, 500, 1, 1.0, "Figure width.")
        self.height_box = spin(style.height_mm, 0, 800, 1, 1.0,
                               "Figure height; 0 lets the panels choose it.")
        self.dpi_box = whole(style.dpi, 72, 2400,
                             "Resolution of the data raster inside the figure.")
        for box, name in ((self.width_box, "width_mm"), (self.height_box, "height_mm")):
            box.valueChanged.connect(lambda v, n=name: self._set_style(n, v))
        self.dpi_box.valueChanged.connect(lambda v: self._set_style("dpi", int(v)))
        form.addRow("Width (mm)", self.width_box)
        form.addRow("Height (mm)", self.height_box)
        form.addRow("Resolution (dpi)", self.dpi_box)

        grid_row = QWidget()
        grid_box = QHBoxLayout(grid_row)
        grid_box.setContentsMargins(0, 0, 0, 0)
        self.rows_box = whole(self.figure.rows, 1, 12, "Rows of panels.", 60)
        self.cols_box = whole(self.figure.cols, 1, 12, "Columns of panels.", 60)
        self.rows_box.valueChanged.connect(self._set_grid)
        self.cols_box.valueChanged.connect(self._set_grid)
        auto = QPushButton("Auto")
        auto.setMaximumWidth(56)
        auto.clicked.connect(self._auto_grid)
        grid_box.addWidget(QLabel("rows"))
        grid_box.addWidget(self.rows_box)
        grid_box.addWidget(QLabel("cols"))
        grid_box.addWidget(self.cols_box)
        grid_box.addWidget(auto)
        grid_box.addStretch(1)
        form.addRow("Grid", grid_row)

        self.gap_x = spin(style.gap_x_mm, 0, 40, 1, 0.5, "Horizontal gap between panels.")
        self.gap_y = spin(style.gap_y_mm, 0, 40, 1, 0.5, "Vertical gap between panels.")
        self.gap_x.valueChanged.connect(lambda v: self._set_style("gap_x_mm", v))
        self.gap_y.valueChanged.connect(lambda v: self._set_style("gap_y_mm", v))
        form.addRow("Column gap (mm)", self.gap_x)
        form.addRow("Row gap (mm)", self.gap_y)

        self.shared_levels = QCheckBox("Same colour scale everywhere")
        self.shared_levels.setToolTip(
            "One colour window over all the panels, so their intensities can "
            "be compared by eye instead of each being stretched to its own "
            "maximum.")
        self.shared_ranges = QCheckBox("Same axis ranges everywhere")
        self.edge_labels = QCheckBox("Tick labels on outer edges only")
        self.edge_labels.setToolTip(
            "The usual look for a grid of panels: the y axis is labelled on "
            "the left column and the x axis on the bottom row. The panels "
            "keep the same size either way, so the grid stays regular.")
        self.equal_aspect = QCheckBox("Keep the data's aspect (1:1)")
        for box, name in ((self.shared_levels, "shared_levels"),
                          (self.shared_ranges, "shared_ranges"),
                          (self.edge_labels, "edge_labels_only"),
                          (self.equal_aspect, "equal_aspect")):
            box.toggled.connect(lambda v, n=name: self._set_style(n, bool(v)))
            form.addRow("", box)

        self.aspect_box = spin(style.panel_aspect, 0, 8, 2, 0.05,
                               "Panel height / width; 0 chooses it.")
        self.aspect_box.valueChanged.connect(lambda v: self._set_style("panel_aspect", v))
        form.addRow("Panel aspect", self.aspect_box)

        buttons = QWidget()
        row = QHBoxLayout(buttons)
        row.setContentsMargins(0, 0, 0, 0)
        add = QPushButton("Add panel...")
        add.setToolTip("Put another dataset from the file list into the figure.")
        add.clicked.connect(self.window_.add_panel_from_list)
        remove = QPushButton("Remove")
        remove.clicked.connect(self._remove_panel)
        letter = QPushButton("Letter")
        letter.clicked.connect(self._letter)
        row.addWidget(add)
        row.addWidget(remove)
        row.addWidget(letter)
        form.addRow("", buttons)

        order = QWidget()
        order_row = QHBoxLayout(order)
        order_row.setContentsMargins(0, 0, 0, 0)
        up = QPushButton("Earlier")
        down = QPushButton("Later")
        up.clicked.connect(lambda: self._move(-1))
        down.clicked.connect(lambda: self._move(+1))
        order_row.addWidget(up)
        order_row.addWidget(down)
        form.addRow("", order)

    def _set_style(self, name, value):
        if self._syncing:
            return
        setattr(self.figure.style, name, value)
        self.redraw()

    def _set_grid(self, *_):
        if self._syncing:
            return
        self.figure.rows = int(self.rows_box.value())
        self.figure.cols = int(self.cols_box.value())
        self.redraw()

    def _auto_grid(self):
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
        self.figure.letter_panels()
        self.sync()

    def _remove_panel(self):
        index = self.panel_combo.currentIndex()
        if len(self.figure.panels) <= 1 or index < 0:
            return
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
        self.label_position = QComboBox()
        self.label_position.addItems(list(CORNERS))
        self.label_position.currentTextChanged.connect(
            lambda v: self._set_style("panel_label_position", v))
        inner.addRow("Position", self.label_position)
        self.label_inside = QCheckBox("Inside the panel")
        self.label_inside.toggled.connect(
            lambda v: self._set_style("panel_label_inside", bool(v)))
        inner.addRow("", self.label_inside)
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

        group = QGroupBox("Text and arrows")
        inner = QFormLayout(group)
        self.note_text = QLineEdit()
        self.note_text.setPlaceholderText("Γ, M, \"gap\" ...")
        self.note_color = ColorButton("#ffffff")
        self.note_bold = QCheckBox("Bold")
        inner.addRow("Text", self.note_text)
        inner.addRow("Colour", self.note_color)
        inner.addRow("", self.note_bold)
        row = QWidget()
        box = QHBoxLayout(row)
        box.setContentsMargins(0, 0, 0, 0)
        place = QPushButton("Place text...")
        place.setToolTip("Then click where it should go.")
        place.clicked.connect(self._place_text)
        arrow = QPushButton("Arrow...")
        arrow.setToolTip("Click the tail, then the head.")
        arrow.clicked.connect(self._place_arrow)
        clear = QPushButton("Clear")
        clear.clicked.connect(lambda: self.apply(
            lambda p: (p.texts.clear(), p.arrows.clear())))
        box.addWidget(place)
        box.addWidget(arrow)
        box.addWidget(clear)
        inner.addRow("", row)
        layout.addWidget(group)

        group = QGroupBox("Inset")
        inner = QFormLayout(group)
        self.inset_box = QCheckBox("Show a magnified region")
        self.inset_box.toggled.connect(self._set_inset)
        inner.addRow("", self.inset_box)
        self.inset_x0 = spin(0, -1e9, 1e9, 4, 0.05, "")
        self.inset_x1 = spin(1, -1e9, 1e9, 4, 0.05, "")
        self.inset_y0 = spin(0, -1e9, 1e9, 4, 0.05, "")
        self.inset_y1 = spin(1, -1e9, 1e9, 4, 0.05, "")
        self.inset_position = QComboBox()
        self.inset_position.addItems(list(CORNERS))
        self.inset_size = spin(0.38, 0.1, 0.9, 2, 0.02,
                               "Fraction of the panel's short side.")
        self.inset_frame = QCheckBox("Frame")
        self.inset_frame.setChecked(True)
        self.inset_indicate = QCheckBox("Outline the region on the panel")
        self.inset_indicate.setChecked(True)
        inner.addRow("x from / to", self._pair(self.inset_x0, self.inset_x1))
        inner.addRow("y from / to", self._pair(self.inset_y0, self.inset_y1))
        inner.addRow("Position", self.inset_position)
        inner.addRow("Size", self.inset_size)
        inner.addRow("", self.inset_frame)
        inner.addRow("", self.inset_indicate)
        for widget, signal in ((self.inset_x0, "valueChanged"), (self.inset_x1, "valueChanged"),
                               (self.inset_y0, "valueChanged"), (self.inset_y1, "valueChanged"),
                               (self.inset_position, "currentTextChanged"),
                               (self.inset_size, "valueChanged"),
                               (self.inset_frame, "toggled"),
                               (self.inset_indicate, "toggled")):
            getattr(widget, signal).connect(self._set_inset)
        pick = QPushButton("Pick the region...")
        pick.setToolTip("Click two opposite corners on the figure.")
        pick.clicked.connect(self._pick_inset)
        inner.addRow("", pick)
        layout.addWidget(group)

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
        self.overlay_list.setMaximumHeight(80)
        inner.addWidget(self.overlay_list)
        row = QWidget()
        box = QHBoxLayout(row)
        box.setContentsMargins(0, 0, 0, 0)
        copy = QPushButton("Copy to all panels")
        copy.setToolTip("The MATLAB tool's \"copy lines to all axes\".")
        copy.clicked.connect(self._copy_overlays)
        clear = QPushButton("Clear")
        clear.clicked.connect(lambda: (self.apply(lambda p: p.overlays.clear()),
                                       self._resync()))
        box.addWidget(copy)
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

    def _place_text(self):
        if not self.note_text.text().strip():
            QMessageBox.information(self, "Text", "Type the text first.")
            return
        text, color = self.note_text.text(), self.note_color.color()
        bold = self.note_bold.isChecked()

        def placed(index, x, y):
            self.figure.panels[index].texts.append(
                TextNote(x=x, y=y, text=text, color=color, bold=bold))
            self.redraw()
        self.window_.pick_point("Click where the text should go", placed)

    def _place_arrow(self):
        color = self.note_color.color()
        start = {}

        def first(index, x, y):
            start.update(index=index, x=x, y=y)
            self.window_.pick_point("Now click the head of the arrow", second)

        def second(index, x, y):
            if index != start.get("index"):
                return
            self.figure.panels[index].arrows.append(
                Arrow(x0=start["x"], y0=start["y"], x1=x, y1=y, color=color))
            self.redraw()
        self.window_.pick_point("Click the tail of the arrow", first)

    def _pick_inset(self):
        corner = {}

        def first(index, x, y):
            corner.update(index=index, x=x, y=y)
            self.window_.pick_point("Now the opposite corner", second)

        def second(index, x, y):
            if index != corner.get("index"):
                return
            self._syncing = True
            self.inset_x0.setValue(min(corner["x"], x))
            self.inset_x1.setValue(max(corner["x"], x))
            self.inset_y0.setValue(min(corner["y"], y))
            self.inset_y1.setValue(max(corner["y"], y))
            self.inset_box.setChecked(True)
            self._syncing = False
            self._set_inset()
        self.window_.pick_point("Click one corner of the region to magnify", first)

    def _set_inset(self, *_):
        if self._syncing:
            return
        if not self.inset_box.isChecked():
            self.apply(lambda p: setattr(p, "inset", None))
            return
        values = dict(x0=self.inset_x0.value(), x1=self.inset_x1.value(),
                      y0=self.inset_y0.value(), y1=self.inset_y1.value(),
                      position=self.inset_position.currentText(),
                      size=self.inset_size.value(),
                      frame=self.inset_frame.isChecked(),
                      indicate=self.inset_indicate.isChecked())
        self.apply(lambda p: setattr(p, "inset", Inset(**values)))

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
        self.cmap_combo.setCurrentText(panel.colormap)
        self.flip_box.setChecked(panel.flip)
        lo, hi = self.figure.panel_levels(panel)
        self.level_lo.setValue(lo)
        self.level_hi.setValue(hi)
        self.gamma_box.setValue(panel.gamma)
        self.smooth_box.setChecked(panel.smooth)
        self.cbar_box.setChecked(panel.colorbar)
        self.cbar_label.setText(panel.colorbar_label)
        bar = panel.scale_bar
        self.bar_box.setChecked(bar is not None)
        if bar is not None:
            self.bar_length.setValue(bar.length)
            self.bar_position.setCurrentText(bar.position)
            self.bar_color.set_color(bar.color)
            self.bar_thickness.setValue(bar.thickness_pt)
            self.bar_unit.setText(bar.unit)
            self.bar_label_box.setChecked(bar.show_label)
        inset = panel.inset
        self.inset_box.setChecked(inset is not None)
        if inset is not None:
            self.inset_x0.setValue(inset.x0)
            self.inset_x1.setValue(inset.x1)
            self.inset_y0.setValue(inset.y0)
            self.inset_y1.setValue(inset.y1)
            self.inset_position.setCurrentText(inset.position)
            self.inset_size.setValue(inset.size)
            self.inset_frame.setChecked(inset.frame)
            self.inset_indicate.setChecked(inset.indicate)
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
        self.window_._tools = None
        super().closeEvent(event)


# --------------------------------------------------------------------------
# Turning this program's datasets into panels
# --------------------------------------------------------------------------
def panel_from_arrays(array, x, y, x_label="", y_label="", name="", info=None,
                      colormap="gray", flip=False, smooth=True) -> Panel:
    return Panel(data=PanelData(array=np.asarray(array, dtype=float),
                                x=np.asarray(x, dtype=float),
                                y=np.asarray(y, dtype=float),
                                x_label=x_label, y_label=y_label, name=name,
                                info=dict(info or {})),
                 colormap=colormap, flip=flip, smooth=smooth)


def panel_from_dataset(data, name="") -> Panel:
    """A Cut (or anything with a 2-D ``scan.value``) as a panel; a curve
    dataset as its curves, with ±σ where it carries uncertainties."""
    scan = data.scan
    from loader.nxs_file import CURVE_KINDS
    if getattr(data, "kind", None) in CURVE_KINDS:
        from tools import curves as C
        from ui.curves import curve_panel, _colour

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
        return curve_panel(scan.x, columns, x_label=labels.get("x", ""),
                           y_label=C.value_label(scan.info),
                           name=name or getattr(data, "source_label", "") or "curve",
                           errors=errors)
    if getattr(scan, "value", None) is None or np.asarray(scan.value).ndim != 2:
        raise ValueError("only a two-dimensional dataset can be a panel on "
                         "its own; take a slice of a map first")
    labels = getattr(scan, "labels", {}) or {}
    return panel_from_arrays(scan.value, scan.x, scan.y,
                             labels.get("x", ""), labels.get("y", ""),
                             name or getattr(data, "source_label", "") or "panel",
                             info=getattr(scan, "info", {}))


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
        self.setWindowTitle("Multi-slice figure")
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
        form.addRow("From / to", PlotToolsDialog._pair(self.from_box, self.to_box))
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
class FigureWindow(QMainWindow):
    """Where a figure for a paper is put together.

    Seeded with one panel by "Plot tools" in a popped-out panel, or with a
    page of them by a map's multi-slice button. The canvas is the figure
    itself at the size it will be printed; everything else is in the tools.
    """

    closed = pyqtSignal(object)

    def __init__(self, panels, title="Figure", parent=None, dataset_source=None):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.figure = Figure(list(panels), style=FigureStyle())
        self.figure.fit_grid()
        self.figure.letter_panels()
        self.painter = FigurePainter(self.figure)
        #: ``() -> [(label, key)]`` and ``key -> data``, so a panel can be
        #: taken from the launcher's list
        self.dataset_source = dataset_source
        self._tools = None

        central = QWidget()
        self.setCentralWidget(central)
        root = QVBoxLayout(central)
        root.setContentsMargins(4, 4, 4, 4)

        bar = QHBoxLayout()
        self.tools_button = QPushButton("Plot tools")
        self.tools_button.setToolTip(
            "Everything that dresses the figure up: panels, axes, text, "
            "colour, annotations and export.")
        self.tools_button.clicked.connect(self.open_tools)
        bar.addWidget(self.tools_button)
        self.add_button = QPushButton("Add panel...")
        self.add_button.clicked.connect(self.add_panel_from_list)
        bar.addWidget(self.add_button)
        bar.addWidget(self._separator())
        for text, fmt in (("PDF", "pdf"), ("PNG", "png"), ("SVG", "svg")):
            button = QPushButton(text)
            button.setMaximumWidth(56)
            button.clicked.connect(lambda _c, f=fmt: self.export_figure(f))
            bar.addWidget(button)
        copy = QPushButton("Copy")
        copy.setMaximumWidth(56)
        copy.setToolTip("Put the figure on the clipboard.")
        copy.clicked.connect(self.copy_to_clipboard)
        bar.addWidget(copy)
        bar.addStretch(1)
        root.addLayout(bar)

        self.canvas = FigureCanvas(self.figure, self)
        self.canvas.panelPicked.connect(self._panel_picked)
        root.addWidget(self.canvas, 1)

        self.hint = QLabel("")
        self.hint.setWordWrap(True)
        root.addWidget(self.hint)
        self.statusBar()
        self.resize(900, 760)
        self.refresh()

    @staticmethod
    def _separator():
        line = QFrame()
        line.setFrameShape(QFrame.VLine)
        line.setFrameShadow(QFrame.Sunken)
        return line

    # -- state -------------------------------------------------------------
    def refresh(self):
        self.canvas.refresh()
        width, height = self.painter.size_mm()
        self.statusBar().showMessage(
            f"{len(self.figure.panels)} panel(s), "
            f"{self.figure.rows} x {self.figure.cols} — "
            f"{width:.1f} x {height:.1f} mm at {self.figure.style.dpi} dpi")

    def _panel_picked(self, index):
        if self._tools is not None:
            self._tools.panel_combo.setCurrentIndex(index)

    def open_tools(self):
        if self._tools is not None and self._tools.isVisible():
            self._tools.raise_()
            self._tools.activateWindow()
            return self._tools
        self._tools = PlotToolsDialog(self)
        self._tools.show()
        return self._tools

    def pick_point(self, message, callback):
        self.hint.setText(message)
        self.canvas.start_point_picking(
            lambda index, x, y: (self.hint.setText(""), callback(index, x, y)))

    # -- panels ------------------------------------------------------------
    def add_panels(self, panels, cols=None):
        self.figure.panels.extend(panels)
        self.figure.fit_grid(cols or self.figure.cols)
        self.figure.letter_panels()
        if self._tools is not None:
            self._tools.sync()
        self.refresh()

    def add_panel_from_list(self):
        """Take another dataset from the launcher's list."""
        source = self.dataset_source
        if not source or not callable(source.get("entries")):
            QMessageBox.information(self, "Add panel",
                                     "No dataset list is available here.")
            return None
        listed = source["entries"]()
        if not listed:
            QMessageBox.information(self, "Add panel", "The list is empty.")
            return None
        labels = [label for label, _key in listed]
        choice, ok = QInputDialog.getItem(self, "Add panel", "Dataset:",
                                           labels, 0, False)
        if not ok:
            return None
        key = listed[labels.index(choice)][1]
        try:
            data = source["loader"](key)
            panel = panel_from_dataset(data, choice.split("   [")[0])
        except Exception as exc:
            QMessageBox.warning(self, "Add panel",
                                 f"That dataset cannot be a panel:\n{exc}")
            return None
        self.add_panels([panel])
        return panel

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

    def write(self, path: str, fmt: str):
        """Write the figure. Returns every file written."""
        if fmt == "pdf":
            return [self.painter.to_pdf(path)]
        if fmt == "svg":
            return [self.painter.to_svg(path)]
        if fmt in ("png", "tif"):
            image = self.painter.to_image()
            export.write_image(path, _qimage_array(image),
                                   pil_format="PNG" if fmt == "png" else "TIFF",
                                   keep_alpha=True, dpi=self.figure.style.dpi)
            return [path]
        if fmt == "layers":
            return self._write_layers(path)
        raise ValueError(f"unknown format {fmt!r}")

    def _write_layers(self, path: str):
        """The two-file route: the data as a raster and the axes as vector
        PostScript, for a figure that will be rearranged in Illustrator.

        One panel only -- a page of panels is what the PDF is for.
        """
        panel = self.figure.panels[self.canvas.selected]
        (x0, x1), (y0, y1) = self.figure.panel_ranges(panel)
        _w, _h, rects = self.figure.layout()
        rect = rects[self.canvas.selected]
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
        QApplication.clipboard().setImage(self.painter.to_image())
        self.statusBar().showMessage("Figure copied to the clipboard")

    def closeEvent(self, event):
        if self._tools is not None:
            self._tools.close()
        self.closed.emit(self)
        super().closeEvent(event)


def _qimage_array(image: QImage) -> np.ndarray:
    image = image.convertToFormat(QImage.Format_RGBA8888)
    size = (image.sizeInBytes() if hasattr(image, "sizeInBytes")
            else image.byteCount())
    buffer = image.constBits()
    buffer.setsize(int(size))
    array = np.frombuffer(buffer, dtype=np.uint8).reshape(
        image.height(), image.bytesPerLine() // 4, 4)
    return np.array(array[:, :image.width(), :], copy=True)
