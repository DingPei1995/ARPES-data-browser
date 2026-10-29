"""
ui/viewers/plot_tools.py
========================
**Plot tools**: guide lines and labels drawn on a live viewer.

Every viewer (and every frozen slice) has a *Plot tools* button on its top
row. It opens a small, non-modal window that draws on the plot in front of
you -- not on a copy of it:

- vertical and horizontal lines (E_F, a high-symmetry point, k = 0),
- a line at any angle through a point (a band's slope, a mirror plane),
- text labels,

each in a colour, width and line style of its own. Lines can be dragged on
the plot; the window follows. Guides belong to the plot they were drawn on
and stay while its slice moves (an energy cut's E_F line stays at E_F), and
the panel's image export draws them too (see :func:`draw_guides`).

This is *not* the figure composer's **Figure tools**, which dress up a
figure for a paper in a window of its own (``ui.panels.figure``). Plot tools are
for looking at the data; Figure tools are for publishing it.

What a window offers is decided by :func:`plot_targets`: the image panels
of a viewer, or the plot of a curve viewer.
"""
from __future__ import annotations

import weakref
from dataclasses import dataclass, field

import numpy as np
import pyqtgraph as pg
from PyQt5.QtCore import QPointF, Qt
from PyQt5.QtGui import QColor, QFont, QImage, QPainter, QPen
from PyQt5.QtWidgets import (QCheckBox, QColorDialog, QComboBox, QDialog,
                             QDialogButtonBox, QDoubleSpinBox, QFormLayout,
                             QGridLayout, QGroupBox, QHBoxLayout, QLabel,
                             QLineEdit, QListWidget, QPushButton, QSpinBox,
                             QVBoxLayout)

from arpes_viewer.ui.common.screen import fits_screen


__all__ = ["Guide", "PlotToolsDialog", "plot_tools_button", "plot_targets",
           "guides_of", "draw_guides", "KINDS", "PALETTE"]

#: what can be drawn, as offered
KINDS = (("vertical", "Vertical line"),
         ("horizontal", "Horizontal line"),
         ("angled", "Line at an angle"),
         ("text", "Text label"))

#: one-click colours: distinct from one another and from the usual
#: colormaps' extremes, so a guide is visible on dark and light data alike
PALETTE = ("#ff3b30", "#ff9500", "#ffd60a", "#34c759", "#00c7be",
           "#0a84ff", "#bf5af2", "#ffffff", "#000000")

STYLES = {"solid": Qt.SolidLine, "dashed": Qt.DashLine,
          "dotted": Qt.DotLine, "dash-dot": Qt.DashDotLine}


@dataclass
class Guide:
    """One guide, in data coordinates. ``angle`` is in degrees from the
    horizontal axis, in data units: the line's slope is tan(angle) in
    y-units per x-unit (90 is vertical whatever the aspect)."""
    kind: str = "vertical"
    x: float = 0.0
    y: float = 0.0
    angle: float = 45.0
    color: str = PALETTE[0]
    width: float = 1.5
    style: str = "dashed"
    text: str = ""
    size: int = 12
    movable: bool = True
    item: object = field(default=None, repr=False, compare=False)

    def line_angle(self) -> float:
        return {"vertical": 90.0, "horizontal": 0.0}.get(self.kind, self.angle)

    def describe(self) -> str:
        what = dict(KINDS)[self.kind]
        if self.kind == "vertical":
            where = f"x = {self.x:.4g}"
        elif self.kind == "horizontal":
            where = f"y = {self.y:.4g}"
        elif self.kind == "angled":
            where = f"({self.x:.4g}, {self.y:.4g}) at {self.angle:.3g}°"
        else:
            where = f"“{self.text}” at ({self.x:.4g}, {self.y:.4g})"
        label = (f"  “{self.text}”"
                 if self.text and self.kind != "text" else "")
        return f"{what}: {where}{label}"


# --------------------------------------------------------------------------
# Where guides go
# --------------------------------------------------------------------------
def _plot_item_of(obj):
    """The pyqtgraph PlotItem behind a panel, an image view or a plot
    widget."""
    view = getattr(obj, "view", obj)
    plot = getattr(view, "plot", None)
    if isinstance(plot, pg.PlotItem):
        return plot
    if isinstance(obj, pg.PlotWidget):
        return obj.getPlotItem()
    if isinstance(obj, pg.PlotItem):
        return obj
    return None


def plot_targets(window):
    """``[(name, plot_item, image_view_or_None), ...]`` a window's Plot
    tools can draw on.

    A window may say for itself (``plot_tool_targets()``); otherwise its
    image panels are used, and failing those its ``plot`` attribute (the
    curve viewer).
    """
    own = getattr(window, "plot_tool_targets", None)
    if callable(own):
        return list(own())
    targets = []
    panels = []
    image_panels = getattr(window, "image_panels", None)
    if callable(image_panels):
        panels = list(image_panels())
    elif getattr(window, "panel", None) is not None:
        panels = [window.panel]
    for i, panel in enumerate(panels):
        plot = _plot_item_of(panel)
        if plot is None:
            continue
        title = getattr(getattr(panel, "title_label", None), "text", None)
        name = title() if callable(title) else f"Panel {i + 1}"
        targets.append((name or f"Panel {i + 1}", plot,
                        getattr(panel, "view", None)))
    if not targets and getattr(window, "plot", None) is not None:
        plot = _plot_item_of(window.plot)
        if plot is not None:
            targets.append(("Plot", plot, None))
    return targets


def guides_of(view) -> list:
    """The guides drawn on ``view`` (an image view), for its export."""
    return list(getattr(view, "_plot_guides", None) or [])


# --------------------------------------------------------------------------
# Drawing on the plot
# --------------------------------------------------------------------------
def _pen(guide: Guide):
    pen = pg.mkPen(QColor(guide.color), width=float(guide.width))
    pen.setStyle(STYLES.get(guide.style, Qt.SolidLine))
    pen.setCosmetic(True)
    return pen


class _DraggableText(pg.TextItem):
    """A text label that can be dragged with the mouse."""

    def __init__(self, *args, on_moved=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.on_moved = on_moved
        self.movable = True
        self._grab = None

    def mouseDragEvent(self, ev):
        if not self.movable or ev.button() != Qt.LeftButton:
            ev.ignore()
            return
        ev.accept()
        parent = self.parentItem()
        here = parent.mapFromScene(ev.scenePos()) if parent else ev.scenePos()
        if ev.isStart():
            self._grab = self.pos() - here
        if self._grab is None:
            return
        self.setPos(here + self._grab)
        if ev.isFinish():
            self._grab = None
        if self.on_moved is not None:
            self.on_moved(self)


def make_item(guide: Guide, on_moved=None):
    """A fresh pyqtgraph item drawing ``guide``."""
    if guide.kind == "text":
        item = _DraggableText(guide.text or "label", color=QColor(guide.color),
                              anchor=(0, 1), on_moved=on_moved)
        font = QFont()
        font.setPointSize(int(guide.size))
        item.setFont(font)
        item.setPos(guide.x, guide.y)
        item.movable = guide.movable
    else:
        options = {}
        if guide.text:
            options = {"label": guide.text,
                       "labelOpts": {"position": 0.92, "movable": True,
                                     "color": QColor(guide.color),
                                     "fill": (0, 0, 0, 90)}}
        item = pg.InfiniteLine(pos=QPointF(guide.x, guide.y),
                               angle=guide.line_angle(), pen=_pen(guide),
                               movable=guide.movable, **options)
        hover = _pen(guide)
        hover.setWidthF(float(guide.width) + 2)
        item.setHoverPen(hover)
        if on_moved is not None:
            item.sigPositionChanged.connect(on_moved)
    item.setZValue(40)
    return item


# --------------------------------------------------------------------------
# Drawing on an exported image
# --------------------------------------------------------------------------
def draw_guides(rgba, guides, *, x_range, y_range, invert_x=False,
                invert_y=False, scale: float = 1.0) -> np.ndarray:
    """``rgba`` (the inside of a panel's frame, as exported) with
    ``guides`` drawn on it. ``x_range`` / ``y_range`` are what the image
    covers; ``scale`` multiplies line widths and text sizes."""
    guides = [g for g in guides if g is not None]
    rgba = np.ascontiguousarray(np.asarray(rgba, dtype=np.uint8))
    if not guides:
        return rgba
    height, width = rgba.shape[:2]
    image = QImage(rgba.data, width, height, width * 4,
                   QImage.Format_RGBA8888).copy()
    x0, x1 = float(min(x_range)), float(max(x_range))
    y0, y1 = float(min(y_range)), float(max(y_range))

    def px(x, y):
        fx = (x - x0) / (x1 - x0) if x1 != x0 else 0.0
        fy = (y - y0) / (y1 - y0) if y1 != y0 else 0.0
        return (QPointF((1 - fx if invert_x else fx) * width,
                        (fy if invert_y else 1 - fy) * height))

    painter = QPainter(image)
    painter.setRenderHint(QPainter.Antialiasing, True)
    painter.setRenderHint(QPainter.TextAntialiasing, True)
    for g in guides:
        color = QColor(g.color)
        if g.kind == "text":
            font = QFont("Helvetica")
            font.setPixelSize(max(6, int(round(g.size * 1.33 * scale))))
            painter.setFont(font)
            painter.setPen(color)
            painter.drawText(px(g.x, g.y), g.text)
            continue
        pen = QPen(color)
        pen.setWidthF(max(float(g.width) * scale, 0.5))
        pen.setStyle(STYLES.get(g.style, Qt.SolidLine))
        painter.setPen(pen)
        if g.kind == "vertical":
            a, b = px(g.x, y0), px(g.x, y1)
        elif g.kind == "horizontal":
            a, b = px(x0, g.y), px(x1, g.y)
        else:
            # The angle is in data units, as pyqtgraph draws it: the line's
            # slope is tan(angle) in y-units per x-unit.
            t = np.radians(g.angle)
            reach = 4.0 * max(x1 - x0, y1 - y0, 1e-12) / max(
                min(abs(np.cos(t)), abs(np.sin(t))), 1e-3)
            a = px(g.x - np.cos(t) * reach, g.y - np.sin(t) * reach)
            b = px(g.x + np.cos(t) * reach, g.y + np.sin(t) * reach)
        painter.drawLine(a, b)
        if g.text:
            font = QFont("Helvetica")
            font.setPixelSize(max(6, int(round(12 * scale))))
            painter.setFont(font)
            painter.setPen(color)
            end = b if g.kind != "vertical" else a
            painter.drawText(QPointF(min(max(end.x(), 4), width - 60),
                                     min(max(end.y(), 14), height - 4)),
                             g.text)
    painter.end()
    out = image.convertToFormat(QImage.Format_RGBA8888)
    ptr = out.constBits()
    ptr.setsize(out.sizeInBytes() if hasattr(out, "sizeInBytes")
                else out.byteCount())
    return np.frombuffer(ptr, dtype=np.uint8).reshape(
        height, out.bytesPerLine() // 4, 4)[:, :width].copy()


# --------------------------------------------------------------------------
# The window
# --------------------------------------------------------------------------
@fits_screen
class PlotToolsDialog(QDialog):
    """Add, edit and remove guides on one of a window's plots.

    The controls edit the guide selected in the list, live; with nothing
    selected they set up the next one, which **Add** draws (at the middle
    of the view, or wherever **Place by clicking** is clicked).
    """

    def __init__(self, window):
        super().__init__(window)
        self._window = weakref.ref(window)
        self.setWindowTitle("Plot tools — guide lines and labels")
        self.setModal(False)
        self._syncing = False
        self._picking = None
        self.guides = {}           # id(plot) -> [Guide]
        self.targets = plot_targets(window)

        layout = QVBoxLayout(self)
        intro = QLabel(
            "Draw on the plot itself: lines in any colour to mark E<sub>F</sub>, "
            "k = 0, a band or a mirror plane, and text labels. Drag a line "
            "on the plot to move it. (For a figure to publish, use the "
            "frozen slice's <i>Figure composer</i>.)")
        intro.setWordWrap(True)
        layout.addWidget(intro)

        form = QFormLayout()
        self.target_combo = QComboBox()
        for name, _plot, _view in self.targets:
            self.target_combo.addItem(name)
        self.target_combo.currentIndexChanged.connect(self._target_changed)
        form.addRow("Draw on", self.target_combo)
        if len(self.targets) < 2:
            self.target_combo.setEnabled(False)
        layout.addLayout(form)

        # -- the guide's settings -------------------------------------------
        box = QGroupBox("Guide")
        grid = QFormLayout(box)
        self.kind_combo = QComboBox()
        for key, text in KINDS:
            self.kind_combo.addItem(text, key)
        grid.addRow("Type", self.kind_combo)

        self.x_box, self.y_box = self._number(), self._number()
        position = QHBoxLayout()
        self.x_label, self.y_label = QLabel("x"), QLabel("y")
        position.addWidget(self.x_label)
        position.addWidget(self.x_box, 1)
        position.addWidget(self.y_label)
        position.addWidget(self.y_box, 1)
        grid.addRow("Position", position)
        self.angle_box = QDoubleSpinBox()
        self.angle_box.setRange(-180.0, 180.0)
        self.angle_box.setDecimals(2)
        self.angle_box.setSuffix("°")
        self.angle_box.setValue(45.0)
        self.angle_box.setToolTip(
            "From the horizontal axis, counter-clockwise positive, in the "
            "plot's own units: the slope is tan(angle) y-units per x-unit, "
            "so on a plot with unequal axes it looks steeper or flatter.")
        grid.addRow("Angle", self.angle_box)
        self.text_edit = QLineEdit()
        self.text_edit.setPlaceholderText("optional for a line, e.g. E_F, Γ")
        grid.addRow("Label", self.text_edit)

        colours = QGridLayout()
        colours.setSpacing(2)
        self.colour_buttons = []
        for i, colour in enumerate(PALETTE):
            button = QPushButton()
            button.setFixedSize(22, 22)
            button.setToolTip(colour)
            button.setStyleSheet(f"background-color: {colour};"
                                 " border: 1px solid #888;")
            button.clicked.connect(lambda _=False, c=colour: self.set_colour(c))
            colours.addWidget(button, 0, i)
            self.colour_buttons.append(button)
        self.more_colour = QPushButton("More...")
        self.more_colour.clicked.connect(self._choose_colour)
        colours.addWidget(self.more_colour, 0, len(PALETTE))
        self.colour_swatch = QLabel()
        self.colour_swatch.setFixedSize(40, 22)
        colours.addWidget(self.colour_swatch, 0, len(PALETTE) + 1)
        grid.addRow("Colour", colours)
        self.colour = PALETTE[0]

        look = QHBoxLayout()
        self.width_box = QDoubleSpinBox()
        self.width_box.setRange(0.5, 12.0)
        self.width_box.setSingleStep(0.5)
        self.width_box.setValue(1.5)
        self.style_combo = QComboBox()
        self.style_combo.addItems(list(STYLES))
        self.style_combo.setCurrentText("dashed")
        self.size_box = QSpinBox()
        self.size_box.setRange(6, 48)
        self.size_box.setValue(12)
        self.size_box.setSuffix(" pt")
        look.addWidget(QLabel("Width"))
        look.addWidget(self.width_box)
        look.addWidget(QLabel("Style"))
        look.addWidget(self.style_combo)
        look.addWidget(QLabel("Text size"))
        look.addWidget(self.size_box)
        grid.addRow("Line", look)
        self.movable_box = QCheckBox("Can be dragged on the plot")
        self.movable_box.setChecked(True)
        grid.addRow("", self.movable_box)
        layout.addWidget(box)

        add_row = QHBoxLayout()
        self.add_button = QPushButton("Add")
        self.add_button.setToolTip("Draw a new guide with these settings, at "
                                   "the position above.")
        self.add_button.clicked.connect(lambda: self.add())
        self.pick_button = QPushButton("Place by clicking")
        self.pick_button.setCheckable(True)
        self.pick_button.setToolTip("Click on the plot to draw a new guide "
                                    "there with these settings.")
        self.pick_button.toggled.connect(self._toggle_picking)
        self.centre_button = QPushButton("Centre of view")
        self.centre_button.setToolTip("Put the position boxes at the middle "
                                      "of what the plot shows.")
        self.centre_button.clicked.connect(self._centre_position)
        add_row.addWidget(self.add_button)
        add_row.addWidget(self.pick_button)
        add_row.addWidget(self.centre_button)
        layout.addLayout(add_row)

        layout.addWidget(QLabel("On this plot (select one to edit it):"))
        self.list = QListWidget()
        self.list.currentRowChanged.connect(self._selected)
        layout.addWidget(self.list, 1)
        list_row = QHBoxLayout()
        self.deselect_button = QPushButton("New guide")
        self.deselect_button.setToolTip("Stop editing the selected guide; "
                                        "the controls then set up a new one.")
        self.deselect_button.clicked.connect(
            lambda: self.list.setCurrentRow(-1))
        self.remove_button = QPushButton("Remove")
        self.remove_button.clicked.connect(self.remove_selected)
        self.clear_button = QPushButton("Remove all")
        self.clear_button.clicked.connect(self.clear)
        list_row.addWidget(self.deselect_button)
        list_row.addWidget(self.remove_button)
        list_row.addWidget(self.clear_button)
        list_row.addStretch(1)
        layout.addLayout(list_row)

        self.status = QLabel("")
        self.status.setWordWrap(True)
        layout.addWidget(self.status)
        close = QDialogButtonBox(QDialogButtonBox.Close)
        close.rejected.connect(self.close)
        layout.addWidget(close)

        for signal in (self.x_box.valueChanged, self.y_box.valueChanged,
                       self.angle_box.valueChanged, self.width_box.valueChanged,
                       self.size_box.valueChanged,
                       self.style_combo.currentIndexChanged,
                       self.text_edit.textChanged, self.movable_box.toggled):
            signal.connect(self._edited)
        self.kind_combo.currentIndexChanged.connect(self._kind_changed)
        self.set_colour(PALETTE[0])
        self._kind_changed()
        self._centre_position()
        self.resize(520, 640)

    @staticmethod
    def _number():
        box = QDoubleSpinBox()
        box.setRange(-1e9, 1e9)
        box.setDecimals(4)
        box.setSingleStep(0.01)
        box.setKeyboardTracking(False)
        return box

    # -- which plot --------------------------------------------------------
    def plot(self):
        index = self.target_combo.currentIndex()
        if 0 <= index < len(self.targets):
            return self.targets[index][1]
        return None

    def image_view(self):
        index = self.target_combo.currentIndex()
        if 0 <= index < len(self.targets):
            return self.targets[index][2]
        return None

    def current_guides(self) -> list:
        plot = self.plot()
        return self.guides.setdefault(id(plot), []) if plot is not None else []

    def _publish(self):
        """Let the image view know its guides, for the export."""
        view = self.image_view()
        if view is not None:
            view._plot_guides = list(self.current_guides())

    def _target_changed(self, *_):
        self._stop_picking()
        self._refresh_list()
        self._centre_position()

    # -- the controls ------------------------------------------------------
    def set_colour(self, colour: str):
        self.colour = QColor(colour).name()
        self.colour_swatch.setStyleSheet(
            f"background-color: {self.colour}; border: 1px solid #888;")
        self._edited()

    def _choose_colour(self):
        colour = QColorDialog.getColor(QColor(self.colour), self,
                                       "Guide colour")
        if colour.isValid():
            self.set_colour(colour.name())

    def kind(self) -> str:
        return self.kind_combo.currentData()

    def _kind_changed(self, *_):
        kind = self.kind()
        self.x_box.setEnabled(kind != "horizontal")
        self.x_label.setEnabled(kind != "horizontal")
        self.y_box.setEnabled(kind != "vertical")
        self.y_label.setEnabled(kind != "vertical")
        self.angle_box.setEnabled(kind == "angled")
        self.width_box.setEnabled(kind != "text")
        self.style_combo.setEnabled(kind != "text")
        self.size_box.setEnabled(kind == "text")
        self._edited()

    def _centre_position(self):
        plot = self.plot()
        if plot is None:
            return
        (xa, xb), (ya, yb) = plot.vb.viewRange()
        self._syncing = True
        try:
            self.x_box.setValue((xa + xb) / 2)
            self.y_box.setValue((ya + yb) / 2)
        finally:
            self._syncing = False
        self._edited()

    def settings(self, **overrides) -> Guide:
        """A Guide from the controls."""
        guide = Guide(kind=self.kind(), x=self.x_box.value(),
                      y=self.y_box.value(), angle=self.angle_box.value(),
                      color=self.colour, width=self.width_box.value(),
                      style=self.style_combo.currentText(),
                      text=self.text_edit.text().strip(),
                      size=self.size_box.value(),
                      movable=self.movable_box.isChecked())
        for key, value in overrides.items():
            setattr(guide, key, value)
        return guide

    def _show(self, guide: Guide):
        self._syncing = True
        try:
            self.kind_combo.setCurrentIndex(self.kind_combo.findData(guide.kind))
            self.x_box.setValue(guide.x)
            self.y_box.setValue(guide.y)
            self.angle_box.setValue(guide.angle)
            self.width_box.setValue(guide.width)
            self.style_combo.setCurrentText(guide.style)
            self.text_edit.setText(guide.text)
            self.size_box.setValue(guide.size)
            self.movable_box.setChecked(guide.movable)
            self.colour = guide.color
            self.colour_swatch.setStyleSheet(
                f"background-color: {guide.color}; border: 1px solid #888;")
        finally:
            self._syncing = False
        self._kind_changed_quiet()

    def _kind_changed_quiet(self):
        self._syncing = True
        try:
            self._kind_changed()
        finally:
            self._syncing = False

    # -- the guides --------------------------------------------------------
    def add(self, guide: Guide = None) -> Guide:
        """Draw ``guide`` (or one from the controls) and select it."""
        plot = self.plot()
        if plot is None:
            return None
        guide = guide or self.settings()
        if guide.kind == "text" and not guide.text:
            guide.text = "label"
        self._draw(guide)
        guides = self.current_guides()
        guides.append(guide)
        self._publish()
        self._refresh_list()
        self.list.setCurrentRow(len(guides) - 1)
        self.status.setText(f"Added: {guide.describe()}")
        return guide

    def _draw(self, guide: Guide):
        plot = self.plot()
        if guide.item is not None:
            self._undraw(guide)
        guide.item = make_item(guide, on_moved=lambda item, g=guide:
                               self._dragged(g, item))
        plot.addItem(guide.item, ignoreBounds=True)

    def _undraw(self, guide: Guide):
        item = guide.item
        guide.item = None
        if item is None:
            return
        scene = item.scene()
        try:
            if item.getViewBox() is not None:
                item.getViewBox().removeItem(item)
            elif scene is not None:
                scene.removeItem(item)
        except RuntimeError:
            pass

    def _dragged(self, guide: Guide, item):
        if isinstance(item, pg.InfiniteLine):
            pos = item.pos()
            guide.x, guide.y = float(pos.x()), float(pos.y())
        else:
            pos = item.pos()
            guide.x, guide.y = float(pos.x()), float(pos.y())
        guides = self.current_guides()
        if guide in guides:
            row = guides.index(guide)
            if self.list.item(row) is not None:
                self.list.item(row).setText(guide.describe())
            if self.list.currentRow() == row:
                self._syncing = True
                try:
                    self.x_box.setValue(guide.x)
                    self.y_box.setValue(guide.y)
                finally:
                    self._syncing = False
        self._publish()

    def selected(self):
        row = self.list.currentRow()
        guides = self.current_guides()
        return guides[row] if 0 <= row < len(guides) else None

    def _selected(self, row):
        guide = self.selected()
        self.remove_button.setEnabled(guide is not None)
        self.deselect_button.setEnabled(guide is not None)
        if guide is not None:
            self._show(guide)

    def _edited(self, *_):
        """Push the controls into the selected guide, live."""
        if self._syncing:
            return
        guide = self.selected()
        if guide is None:
            return
        new = self.settings()
        for name in ("kind", "x", "y", "angle", "color", "width", "style",
                     "text", "size", "movable"):
            setattr(guide, name, getattr(new, name))
        if guide.kind == "text" and not guide.text:
            guide.text = "label"
        self._draw(guide)
        self._publish()
        item = self.list.currentItem()
        if item is not None:
            item.setText(guide.describe())

    def _refresh_list(self):
        self.list.blockSignals(True)
        self.list.clear()
        for guide in self.current_guides():
            self.list.addItem(guide.describe())
        self.list.blockSignals(False)
        self.list.setCurrentRow(-1)
        self._selected(-1)

    def remove_selected(self):
        guide = self.selected()
        if guide is None:
            return
        self._undraw(guide)
        self.current_guides().remove(guide)
        self._publish()
        self._refresh_list()

    def clear(self):
        for guide in self.current_guides():
            self._undraw(guide)
        self.current_guides().clear()
        self._publish()
        self._refresh_list()

    # -- placing by clicking -------------------------------------------------
    def _toggle_picking(self, on: bool):
        plot = self.plot()
        if not on or plot is None:
            self._stop_picking()
            return
        scene = plot.scene()
        if scene is None:
            self.pick_button.setChecked(False)
            return
        scene.sigMouseClicked.connect(self._clicked)
        self._picking = scene
        self.status.setText("Click on the plot to place the guide; "
                            "untick “Place by clicking” to stop.")

    def _stop_picking(self):
        if self._picking is not None:
            try:
                self._picking.sigMouseClicked.disconnect(self._clicked)
            except (TypeError, RuntimeError):
                pass
            self._picking = None
        if self.pick_button.isChecked():
            self.pick_button.blockSignals(True)
            self.pick_button.setChecked(False)
            self.pick_button.blockSignals(False)

    def _clicked(self, event):
        plot = self.plot()
        if plot is None or event.button() != Qt.LeftButton:
            return
        if not plot.vb.sceneBoundingRect().contains(event.scenePos()):
            return
        point = plot.vb.mapSceneToView(event.scenePos())
        event.accept()
        self.list.setCurrentRow(-1)
        self.add(self.settings(x=float(point.x()), y=float(point.y())))

    def closeEvent(self, event):
        # The guides stay on the plot; only the click-to-place stops.
        self._stop_picking()
        super().closeEvent(event)


def open_plot_tools(window):
    """Show ``window``'s Plot tools, creating them on first use (one per
    window, so its guides are kept while it is closed and reopened)."""
    dialog = getattr(window, "_plot_tools_dialog", None)
    if dialog is None:
        dialog = PlotToolsDialog(window)
        window._plot_tools_dialog = dialog
    dialog.show()
    dialog.raise_()
    dialog.activateWindow()
    return dialog


def plot_tools_button(window) -> QPushButton:
    """The **Plot tools** button for a window's top row."""
    button = QPushButton("Plot tools")
    button.setToolTip(
        "Draw on this plot: guide lines (vertical, horizontal or at any "
        "angle) and text labels, each in its own colour, width and style. "
        "Lines can be dragged; image exports include them.")
    button.clicked.connect(lambda *_: open_plot_tools(window))
    return button
