"""
tools/output/figure.py
======================
The model behind the figure composer: what a publication figure *is*, laid
out in millimetres, with no Qt in sight. The painting of it lives in
``ui.panels.figure``; everything here can be exercised headlessly.

Why a model at all, rather than a pile of buttons that poke at the plot
being displayed -- which is what the lab's ``plot_tools2_demo.m`` does with
``set(gca, ...)``. Three reasons:

* A figure in a paper is several panels that must agree with one another --
  same colour scale, same limits, labels only on the outside edges. That is
  a property of the *figure*, not of whichever axes MATLAB last touched.
* What you see has to be what you get. One model, one layout calculation,
  one painter: the preview on screen, the PNG, the PDF and the vector axes
  all come from the same arithmetic, so a figure cannot look right in the
  window and wrong in the file.
* A paper's figures should look alike. A model can be written to JSON and
  applied to the next figure (:meth:`Figure.style_dict`); a pile of
  button presses cannot.

Everything is in millimetres, because that is the unit a journal's column
width comes in. The painter multiplies by one scale factor at the end.
"""
from __future__ import annotations

import json
import math
import uuid
from dataclasses import dataclass, field, asdict

import numpy as np

from tools.output.export import (MM_PER_INCH, PT_PER_MM, tick_values,
                                 tick_format)

__all__ = [
    "Figure", "Panel", "PanelData", "FigureStyle", "Line", "TextNote",
    "Arrow", "ScaleBar", "Inset", "SecondaryAxis", "Overlay", "Shape",
    "empty_data", "picture_data", "JOURNAL_PRESETS", "parse_markup", "text_width_mm", "pt_to_mm",
]

#: Widths journals actually ask for, with a sensible default font size.
#: (label, width_mm, font_pt)
JOURNAL_PRESETS = (
    ("Nature, single column", 89.0, 7.0),
    ("Nature, double column", 183.0, 7.0),
    ("APS/PRB, single column", 3.375 * MM_PER_INCH, 8.0),
    ("APS/PRB, double column", 6.75 * MM_PER_INCH, 8.0),
    ("Science, single column", 55.0, 7.0),
    ("Elsevier, single column", 90.0, 8.0),
    ("Elsevier, double column", 190.0, 8.0),
    ("Slide / talk (16:9 half)", 120.0, 12.0),
)


def pt_to_mm(points: float) -> float:
    return float(points) / PT_PER_MM


# --------------------------------------------------------------------------
# Text: a little markup, because "k_x" and "A^-1" are most of what a label is
# --------------------------------------------------------------------------
def parse_markup(text: str):
    """Split a label into ``(run, kind)`` with kind in ``base/sub/super``.

    ``_`` and ``^`` take the next character, or a ``{...}`` group:
    ``k_x``, ``E_F``, ``Å^-1``, ``k_{∥}``. A backslash escapes
    either marker. This is not TeX and deliberately is not: the two things
    an ARPES axis ever needs are a subscript and a superscript, and a real
    TeX renderer is a dependency and a font problem.
    """
    runs, buffer, index = [], "", 0
    text = str(text or "")
    while index < len(text):
        char = text[index]
        if char == "\\" and index + 1 < len(text) and text[index + 1] in "_^\\":
            buffer += text[index + 1]
            index += 2
            continue
        if char in "_^" and index + 1 < len(text):
            if buffer:
                runs.append((buffer, "base"))
                buffer = ""
            kind = "sub" if char == "_" else "super"
            index += 1
            if text[index] == "{":
                end = text.find("}", index)
                end = len(text) if end < 0 else end
                runs.append((text[index + 1:end], kind))
                index = end + 1
            else:
                # TeX takes a single character; this takes an optional sign
                # and the word after it, because what people type is
                # "A^-1" and "k_max", and making them write "A^{-1}" to get
                # the obvious thing is a tax with no benefit here.
                end = index
                if text[end] in "+-−":
                    end += 1
                while end < len(text) and (text[end].isalnum() or text[end] == "."):
                    end += 1
                end = max(end, index + 1)
                runs.append((text[index:end], kind))
                index = end
            continue
        buffer += char
        index += 1
    if buffer:
        runs.append((buffer, "base"))
    return runs


def plain_text(text: str) -> str:
    """The markup stripped, for width estimates and file names."""
    return "".join(run for run, _ in parse_markup(text))


def text_width_mm(text: str, font_pt: float) -> float:
    """Roughly how wide a label is. Helvetica averages about 0.55 em over
    digits and lower case; sub/superscripts are drawn smaller. Only used to
    reserve margins, and the painter re-measures with real font metrics
    before drawing, so a few per cent out is harmless."""
    width = 0.0
    for run, kind in parse_markup(text):
        factor = 0.58 if kind == "base" else 0.58 * 0.72
        width += factor * len(run)
    return width * pt_to_mm(font_pt)


# --------------------------------------------------------------------------
# The things that can sit on a panel
# --------------------------------------------------------------------------
@dataclass
class Line:
    """A reference line across the panel: E_F, k = 0, a gap edge."""
    value: float = 0.0
    axis: str = "y"            # "y" = horizontal line at this y, "x" = vertical
    color: str = "#ffffff"
    width_pt: float = 0.8
    style: str = "dashed"      # solid / dashed / dotted / dashdot
    label: str = ""


@dataclass
class TextNote:
    """A label anchored in data coordinates (Γ, M, a band name)."""
    x: float = 0.0
    y: float = 0.0
    text: str = ""
    color: str = "#ffffff"
    size_pt: float = 0.0       # 0 = the figure's annotation size
    anchor: str = "center"     # left / center / right
    bold: bool = False


@dataclass
class Arrow:
    """A straight arrow in data coordinates, pointing at something."""
    x0: float = 0.0
    y0: float = 0.0
    x1: float = 0.0
    y1: float = 0.0
    color: str = "#ffffff"
    width_pt: float = 0.8
    head_mm: float = 1.6


@dataclass
class ScaleBar:
    """A bar of known length, for a real-space map where the axes are
    usually removed."""
    length: float = 0.0        # in data units; 0 = choose a round one
    position: str = "bottom right"
    color: str = "#ffffff"
    unit: str = ""             # taken from the x label when empty
    thickness_pt: float = 2.5
    show_label: bool = True


@dataclass
class Inset:
    """A magnified corner of the same panel."""
    x0: float = 0.0
    x1: float = 0.0
    y0: float = 0.0
    y1: float = 0.0
    position: str = "top right"
    size: float = 0.38         # fraction of the panel's short side
    frame: bool = True
    indicate: bool = True      # outline the region on the main panel


@dataclass
class SecondaryAxis:
    """A second scale on the opposite edge -- binding energy against kinetic,
    Å⁻¹ against degrees, k in units of π/a.

    ``value2 = scale * value + offset``, which covers every one of those.
    """
    enabled: bool = False
    scale: float = 1.0
    offset: float = 0.0
    label: str = ""


@dataclass
class Overlay:
    """A polyline in data coordinates.

    This is the hook a Brillouin-zone tool will fill: it hands over a list
    of these and the panel draws them. Nothing in this module knows or cares
    what the lines mean, so a BZ panel can be written separately without
    touching the figure code.
    """
    xs: list = field(default_factory=list)
    ys: list = field(default_factory=list)
    color: str = "#ff00ff"
    width_pt: float = 0.8
    style: str = "solid"
    closed: bool = False
    name: str = ""
    #: a symbol at every point -- "o", "s" (square), "^", "v", "d" -- for
    #: measured points rather than a curve; "" draws none
    marker: str = ""
    marker_mm: float = 1.0
    #: whether the points are joined (False: symbols only)
    line: bool = True
    #: error bars, one per point (empty: none), in data units
    y_errors: list = field(default_factory=list)
    x_errors: list = field(default_factory=list)


@dataclass
class Shape:
    """Something drawn on the page itself, in millimetres from its top
    left corner: an arrow from one panel to another, a guide line, a box
    round a feature, a note.

    Not tied to a panel's data, because on a manuscript figure these cross
    panels and margins freely -- the panel-bound annotations (a reference
    line at E_F, a label at Γ) are the panel's own lists.
    """
    kind: str = "arrow"        # arrow / line / rect / ellipse / text
    x0: float = 0.0
    y0: float = 0.0
    x1: float = 0.0
    y1: float = 0.0
    color: str = "#000000"
    width_pt: float = 0.8
    style: str = "solid"       # solid / dashed / dotted / dashdot
    fill: str = ""             # a box's or ellipse's fill; "" = none
    text: str = ""
    size_pt: float = 0.0       # 0 = the figure's annotation size
    bold: bool = False
    head_mm: float = 1.6       # an arrow's head
    #: an arrow's head at both ends
    double: bool = False

    #: the kinds the composer's drawing tools make
    KINDS = ("arrow", "line", "rect", "ellipse", "text")

    def bounds(self):
        """``(left, top, right, bottom)`` in mm (a text: its anchor only)."""
        return (min(self.x0, self.x1), min(self.y0, self.y1),
                max(self.x0, self.x1), max(self.y0, self.y1))

    def move(self, dx: float, dy: float):
        self.x0 += dx
        self.x1 += dx
        self.y0 += dy
        self.y1 += dy


@dataclass(eq=False)
class PanelData:
    """The measurement a panel draws. Not part of the saved style."""
    array: np.ndarray = None       # (nx, ny)
    x: np.ndarray = None
    y: np.ndarray = None
    x_label: str = ""
    y_label: str = ""
    name: str = ""
    kind: str = "cut"
    info: dict = field(default_factory=dict)
    #: for ``kind == "scene3d"``: the 3-D view to draw, a snapshot made by
    #: ui.panels.volume (faces, camera, colours, axes). The panel's own x
    #: and y then run 0..1 across the picture, for placing notes on it.
    scene: object = None
    #: what the panel was made from, traced back as far as it goes: the
    #: original measurement file (and entry), a calculation file, a picture
    #: file -- ``tools.output.figfile.source_record`` dicts. Saved with the
    #: figure, and copied into its "data for ..." folder.
    sources: list = field(default_factory=list)

    def is_plot(self) -> bool:
        """Whether this is data drawn on axes (not a slot or a picture)."""
        return self.kind not in ("empty", "picture", "scene3d")


def empty_data(name: str = "") -> PanelData:
    """An empty slot of a figure, waiting for a dataset or a picture."""
    return PanelData(array=np.zeros((2, 2)), x=np.array([0.0, 1.0]),
                     y=np.array([0.0, 1.0]), name=name, kind="empty")


def picture_data(rgba, name: str = "") -> PanelData:
    """A picture read from a file (a sketch, a photograph, a figure made
    elsewhere): ``rgba`` is ``(height, width, 4)`` uint8. Drawn as it is,
    without axes."""
    rgba = np.asarray(rgba, dtype=np.uint8)
    return PanelData(array=rgba, x=np.array([0.0, 1.0]),
                     y=np.array([0.0, 1.0]), name=name, kind="picture")


@dataclass(eq=False)
class Panel:
    """One picture in the figure, and everything drawn on it.

    Compared by identity (``eq=False``): a panel is itself, not its
    contents -- comparing two by value compares their arrays, and
    ``list.index`` could not find one."""
    data: PanelData = None
    colormap: str = "gray"
    flip: bool = False
    levels: tuple = None           # None -> from the data
    gamma: float = 1.0
    smooth: bool = True
    x_range: tuple = None          # None -> the data's own extent
    y_range: tuple = None
    invert_x: bool = False
    invert_y: bool = False
    title: str = ""
    x_label: str = None            # None -> the data's own label
    y_label: str = None
    label: str = ""                # "(a)"
    colorbar: bool = False
    colorbar_label: str = ""
    #: what the colour bar's ends say: "values" (the levels), "high-low"
    #: (the words, as a figure of raw intensity usually has), or "none"
    colorbar_ends: str = "values"
    #: the colour bar's length, as a fraction of the panel's height (it
    #: sits at the bottom of the panel's edge)
    colorbar_length: float = 1.0
    lines: list = field(default_factory=list)
    texts: list = field(default_factory=list)
    arrows: list = field(default_factory=list)
    overlays: list = field(default_factory=list)
    scale_bar: ScaleBar = None
    inset: Inset = None
    top_axis: SecondaryAxis = field(default_factory=SecondaryAxis)
    right_axis: SecondaryAxis = field(default_factory=SecondaryAxis)
    #: where the panel's data box sits on the page, ``(x, y, w, h)`` in mm
    #: from the top left; None while the figure is still a plain grid
    frame: tuple = None
    #: "panel" -- one of (a), (b), (c) -- or "inset", a small picture on
    #: top of another panel, which takes no letter
    role: str = "panel"
    #: an inset that magnifies part of another panel: that panel's uid
    zoom_of: str = ""
    #: tick numbers and axis label along x / y (on a page laid out by hand;
    #: a grid decides this from the panel's place in it)
    x_labels: bool = True
    y_labels: bool = True
    #: this panel's font size; 0 = the figure's
    font_pt: float = 0.0
    #: named ticks along x instead of numbers: ``[(value, text), ...]`` --
    #: the high-symmetry points of a band-structure path (Γ, M, K, ...)
    x_tick_labels: list = field(default_factory=list)
    uid: str = field(default_factory=lambda: uuid.uuid4().hex[:10])

    # -- what the painter asks for ---------------------------------------
    def label_x(self) -> str:
        return self.data.x_label if self.x_label is None else self.x_label

    def label_y(self) -> str:
        return self.data.y_label if self.y_label is None else self.y_label

    def extent(self):
        """The data's outer bounds: the axis range grown by half a pixel at
        each end, since a pixel is centred on its value."""
        x, y = np.asarray(self.data.x, float), np.asarray(self.data.y, float)
        if self.data.kind in ("scene3d", "picture", "empty"):
            return (0.0, 1.0, 0.0, 1.0)       # a picture, not pixels on axes
        if self.data.kind == "stack":
            # Curves only: the "image" is a blank placeholder of two pixels,
            # and growing it by half a pixel doubled both ranges -- the curves
            # sat in the middle half of the panel. A curve is not a pixel:
            # its x range is its own, and y gets a small margin so the lines
            # do not run along the frame.
            margin = 0.04 * (float(np.nanmax(y) - np.nanmin(y)) or 1.0)
            return (float(np.nanmin(x)), float(np.nanmax(x)),
                    float(np.nanmin(y) - margin), float(np.nanmax(y) + margin))
        dx = (x[1] - x[0]) if x.size > 1 else 1.0
        dy = (y[1] - y[0]) if y.size > 1 else 1.0
        return (float(x[0] - dx / 2), float(x[-1] + dx / 2),
                float(y[0] - dy / 2), float(y[-1] + dy / 2))

    def ranges(self):
        x0, x1, y0, y1 = self.extent()
        xr = self.x_range or (min(x0, x1), max(x0, x1))
        yr = self.y_range or (min(y0, y1), max(y0, y1))
        return (float(min(xr)), float(max(xr))), (float(min(yr)), float(max(yr)))

    def auto_levels(self):
        if not self.data.is_plot() and self.data.kind != "scene3d":
            return (0.0, 1.0)
        finite = np.asarray(self.data.array, float)
        finite = finite[np.isfinite(finite)]
        if finite.size == 0:
            return (0.0, 1.0)
        return (float(finite.min()), float(finite.max()))

    def effective_levels(self):
        return tuple(self.levels) if self.levels else self.auto_levels()


@dataclass
class FigureStyle:
    """Everything that should look the same across a paper's figures."""
    width_mm: float = 89.0
    height_mm: float = 0.0         # 0 -> from the panels' own shape
    dpi: int = 600
    font: str = "Helvetica"
    font_pt: float = 7.0
    label_pt: float = 0.0          # 0 -> font_pt
    tick_pt: float = 0.0           # 0 -> font_pt
    annotation_pt: float = 0.0     # 0 -> font_pt
    panel_label_pt: float = 0.0    # 0 -> font_pt + 1
    panel_label_bold: bool = True
    panel_label_position: str = "top left"
    #: outside by default: a letter inside the panel lands on the data, and
    #: whether it needs to be black or white then depends on what the data
    #: happens to look like in that corner
    panel_label_inside: bool = False
    panel_label_color: str = "#ffffff"
    axis_width_pt: float = 0.6
    axis_color: str = "#000000"
    background: str = "#ffffff"
    tick_direction: str = "in"
    tick_len_mm: float = 1.2
    minor_ticks: int = 1
    tick_sides: str = "all"        # all / left-bottom
    x_tick_step: float = 0.0       # 0 -> chosen automatically
    y_tick_step: float = 0.0
    margin_left_mm: float = 1.0
    margin_right_mm: float = 1.0
    margin_top_mm: float = 1.0
    margin_bottom_mm: float = 1.0
    gap_x_mm: float = 2.0
    gap_y_mm: float = 2.0
    colorbar_mm: float = 2.2
    colorbar_gap_mm: float = 1.0
    shared_levels: bool = False
    shared_ranges: bool = False
    edge_labels_only: bool = True
    equal_aspect: bool = False
    panel_aspect: float = 0.0      # 0 -> from the data/ranges


class Figure:
    """A grid of panels plus the style they share."""

    def __init__(self, panels=None, rows: int = 1, cols: int = 1,
                 style: FigureStyle = None, title: str = ""):
        self.panels = list(panels or [])
        self.rows = max(int(rows), 1)
        self.cols = max(int(cols), 1)
        self.style = style or FigureStyle()
        self.title = title
        #: arrows, lines, boxes and notes drawn on the page (see Shape)
        self.shapes = []

    # -- a page laid out by hand --------------------------------------------
    @property
    def free(self) -> bool:
        """Whether the panels sit where the user put them (each has a
        ``frame``) rather than in the automatic grid."""
        return bool(self.panels) and all(p.frame is not None for p in self.panels)

    def main_panels(self):
        return [p for p in self.panels if p.role != "inset"]

    def panel_by_uid(self, uid: str):
        return next((p for p in self.panels if p.uid == uid), None)

    def freeze(self):
        """Turn the automatic grid into frames the user can move.

        The page keeps the size the grid gave it, and each panel keeps its
        place, so nothing on the figure moves at the moment it becomes
        editable by hand. Tick labels stay where the grid had them.
        """
        if self.free:
            return self
        width, height, rects = self.layout()
        for rect in rects:
            panel = self.panels[rect["index"]]
            if panel.frame is None:
                panel.frame = (rect["x"], rect["y"], rect["w"], rect["h"])
                panel.x_labels = self.shows_x_labels(rect["index"])
                panel.y_labels = self.shows_y_labels(rect["index"])
        self.style.width_mm = float(width)
        self.style.height_mm = float(height)
        return self

    def grid_frames(self, rows: int, cols: int, count: int = None):
        """``count`` data boxes in a ``rows`` x ``cols`` grid filling the
        page, with room for the axis labels round each -- what *Arrange as
        grid* and a new blank figure lay out."""
        style = self.style
        saved = (self.rows, self.cols, list(self.panels))
        count = rows * cols if count is None else count
        try:
            self.rows, self.cols = max(int(rows), 1), max(int(cols), 1)
            if not self.panels:
                self.panels = [Panel(data=empty_data())]
            left, bottom, top, right = self.margins_mm()
        finally:
            self.rows, self.cols, self.panels = saved
        height = style.height_mm if style.height_mm > 0 else None
        inner_w = (style.width_mm - style.margin_left_mm - style.margin_right_mm
                   - (cols - 1) * style.gap_x_mm)
        cell_w = max(inner_w / cols, 4.0)
        data_w = max(cell_w - left - right, 2.0)
        if height:
            inner_h = (height - style.margin_top_mm - style.margin_bottom_mm
                       - (rows - 1) * style.gap_y_mm)
            cell_h = max(inner_h / rows, 4.0)
            data_h = max(cell_h - top - bottom, 2.0)
        else:
            data_h = max(data_w * (style.panel_aspect or self.DEFAULT_ASPECT), 2.0)
            cell_h = data_h + top + bottom
        frames = []
        for index in range(count):
            row, col = divmod(index, cols)
            frames.append((style.margin_left_mm + col * (cell_w + style.gap_x_mm) + left,
                           style.margin_top_mm + row * (cell_h + style.gap_y_mm) + top,
                           data_w, data_h))
        page_h = (cell_h * rows + (rows - 1) * style.gap_y_mm
                  + style.margin_top_mm + style.margin_bottom_mm)
        return frames, (height or page_h)

    def arrange_grid(self, rows: int, cols: int, fill_height: bool = False):
        """Put the panels (not the insets) back on a regular grid, in their
        reading order; insets travel with the panel they sit on.

        The panels take the figure's panel shape and the page's height
        follows -- or, with ``fill_height``, they share out the page's
        present height.
        """
        mains = self.reading_order()
        rows = max(int(rows), int(math.ceil(len(mains) / max(int(cols), 1))), 1)
        kept = self.style.height_mm
        if not fill_height:
            self.style.height_mm = 0.0
        try:
            frames, height = self.grid_frames(rows, cols, len(mains))
        finally:
            self.style.height_mm = kept
        self.style.height_mm = float(height)
        for panel, frame in zip(mains, frames):
            old = panel.frame
            panel.frame = frame
            if old is not None:
                self._carry_insets(panel, old, frame)
        self.rows, self.cols = rows, int(cols)
        self.auto_edge_labels()
        return self

    def _carry_insets(self, panel, old, new):
        """Move and scale the insets that sat on ``panel`` (whose frame
        was ``old``) along with it."""
        ox, oy, ow, oh = old
        nx, ny, nw, nh = new
        if ow <= 0 or oh <= 0:
            return
        for inset in self.panels:
            if inset.role != "inset" or inset.frame is None or inset is panel:
                continue
            if not self._sits_on(inset, panel, old):
                continue
            x, y, w, h = inset.frame
            inset.frame = (nx + (x - ox) * nw / ow, ny + (y - oy) * nh / oh,
                           w * nw / ow, h * nh / oh)

    @staticmethod
    def _sits_on(inset, panel, box) -> bool:
        if inset.zoom_of:
            return inset.zoom_of == panel.uid
        cx = inset.frame[0] + inset.frame[2] / 2
        cy = inset.frame[1] + inset.frame[3] / 2
        return box[0] <= cx <= box[0] + box[2] and box[1] <= cy <= box[1] + box[3]

    def host_of(self, inset):
        """The panel an inset sits on: the one its zoom names, or else the
        panel whose box holds the inset's centre."""
        if inset.zoom_of:
            host = self.panel_by_uid(inset.zoom_of)
            if host is not None:
                return host
        if inset.frame is None:
            return None
        for panel in self.main_panels():
            if panel.frame is not None and self._sits_on(inset, panel, panel.frame):
                return panel
        return None

    def move_panel(self, panel, frame):
        """Give ``panel`` a new frame; its insets follow."""
        old = panel.frame
        panel.frame = tuple(float(v) for v in frame)
        if old is not None and panel.role != "inset":
            self._carry_insets(panel, old, panel.frame)

    def reading_order(self):
        """The panels (not insets) top to bottom, then left to right: rows
        are panels whose tops are within half a panel of one another."""
        mains = self.main_panels()
        if not all(p.frame is not None for p in mains):
            return mains
        ordered = sorted(mains, key=lambda p: (p.frame[1], p.frame[0]))
        rows, current = [], []
        for panel in ordered:
            if current and panel.frame[1] > current[0].frame[1] + 0.5 * current[0].frame[3]:
                rows.append(current)
                current = []
            current.append(panel)
        if current:
            rows.append(current)
        return [p for row in rows for p in sorted(row, key=lambda p: p.frame[0])]

    def auto_edge_labels(self):
        """On a hand-made page, label only the outer edges (when the style
        says so): a panel keeps its y numbers unless another panel sits
        right beside it on the left at the same height, and its x numbers
        unless one sits right below it."""
        mains = [p for p in self.main_panels() if p.frame is not None]
        for panel in mains:
            if not self.style.edge_labels_only:
                panel.x_labels = panel.y_labels = True
                continue
            x, y, w, h = panel.frame
            left = below = False
            for other in mains:
                if other is panel:
                    continue
                ox, oy, ow, oh = other.frame
                overlap_v = min(y + h, oy + oh) - max(y, oy)
                overlap_h = min(x + w, ox + ow) - max(x, ox)
                if (ox + ow <= x + 0.5 and x - (ox + ow) < w
                        and overlap_v > 0.6 * min(h, oh)
                        and abs(oy - y) < 1.0 and abs(oh - h) < 1.0):
                    left = True
                if (oy >= y + h - 0.5 and oy - (y + h) < h
                        and overlap_h > 0.6 * min(w, ow)
                        and abs(ox - x) < 1.0 and abs(ow - w) < 1.0):
                    below = True
            panel.y_labels = not left
            panel.x_labels = not below

    def place_new(self, panel):
        """A frame for a panel added to a hand-made page: the size of the
        last panel, below everything there is; the page grows to hold it."""
        mains = [p for p in self.panels if p.frame is not None and p.role != "inset"]
        style = self.style
        if mains:
            w, h = mains[-1].frame[2], mains[-1].frame[3]
            left = min(p.frame[0] for p in mains)
        else:
            frames, _ = self.grid_frames(1, 1)
            left, _top, w, h = frames[0]
        bottom = max([p.frame[1] + p.frame[3] for p in self.panels
                      if p.frame is not None] + [0.0])
        top_margin = self.margins_mm((w, h))[2]
        y = bottom + self.margins_mm((w, h))[1] + style.gap_y_mm + top_margin
        panel.frame = (left, y, w, h)
        needed = y + h + self.margins_mm((w, h))[1] + style.margin_bottom_mm
        if needed > style.height_mm:
            style.height_mm = float(needed)
        return panel

    # -- the grid ---------------------------------------------------------
    def fit_grid(self, cols: int = None):
        """Choose a grid that holds every panel."""
        n = max(len(self.panels), 1)
        if cols:
            self.cols = max(int(cols), 1)
        else:
            self.cols = min(n, max(1, int(round(math.sqrt(n * 1.4)))))
        self.rows = int(math.ceil(n / self.cols))
        return self.rows, self.cols

    def cell(self, index: int):
        return divmod(int(index), self.cols)

    def shows_y_labels(self, index: int) -> bool:
        if self.free:
            return bool(self.panels[index].y_labels)
        if not self.style.edge_labels_only:
            return True
        return self.cell(index)[1] == 0

    def shows_x_labels(self, index: int) -> bool:
        if self.free:
            return bool(self.panels[index].x_labels)
        if not self.style.edge_labels_only:
            return True
        row, col = self.cell(index)
        # the bottom-most panel of this column, which is not always the last
        # row when the grid is not full
        below = (row + 1) * self.cols + col
        return row == self.rows - 1 or below >= len(self.panels)

    # -- shared settings --------------------------------------------------
    def common_levels(self):
        """One colour window over every panel, for the figure whose panels
        must be comparable by eye."""
        lows, highs = [], []
        for panel in self.panels:
            if not panel.data.is_plot():
                continue
            lo, hi = panel.auto_levels()
            lows.append(lo)
            highs.append(hi)
        if not lows:
            return (0.0, 1.0)
        return (float(min(lows)), float(max(highs)))

    def common_ranges(self):
        xs, ys = [], []
        for panel in self.panels:
            if not panel.data.is_plot() or panel.role == "inset":
                continue
            (x0, x1), (y0, y1) = panel.ranges()
            xs += [x0, x1]
            ys += [y0, y1]
        if not xs:
            return (0.0, 1.0), (0.0, 1.0)
        return (min(xs), max(xs)), (min(ys), max(ys))

    def panel_levels(self, panel: Panel):
        if self.style.shared_levels and panel.levels is None:
            return self.common_levels()
        return panel.effective_levels()

    def panel_ranges(self, panel: Panel):
        if (self.style.shared_ranges and panel.role != "inset"
                and panel.data.is_plot()):
            return self.common_ranges()
        return panel.ranges()

    def letter_panels(self, template: str = "({letter})", start: int = 0):
        """(a), (b), (c)... in reading order; insets take no letter."""
        for panel in self.panels:
            if panel.role == "inset":
                panel.label = ""
        for i, panel in enumerate(self.reading_order()):
            letter = chr(ord("a") + (i + start) % 26)
            panel.label = template.format(letter=letter, LETTER=letter.upper(),
                                          number=i + start + 1)

    # -- layout -----------------------------------------------------------
    def tick_targets(self, width_mm: float, height_mm: float):
        """Roughly how many labelled ticks fit on a panel this size.

        The painter and the margin calculation must agree on this, or the
        margin is reserved for four short numbers and five long ones get
        drawn -- which is how a y label ends up half off the page.
        """
        em = pt_to_mm(self.style.tick_pt or self.style.font_pt)
        x = int(np.clip(width_mm / (6.5 * em), 2, 8))
        y = int(np.clip(height_mm / (3.2 * em), 2, 8))
        return x, y

    def margins_mm(self, panel_size=None):
        """``(left, bottom, top, right)`` reserved inside every cell for the
        axis furniture.

        The same amount is reserved for every panel whether or not it draws
        labels, so that all the data rectangles come out identical and the
        grid lines up. Hiding a label frees no space -- that is the point:
        panels in a column must be the same size.
        """
        style = self.style
        tick_pt = style.tick_pt or style.font_pt
        label_pt = style.label_pt or style.font_pt
        gap = 0.35 * pt_to_mm(tick_pt)
        outward = style.tick_len_mm if style.tick_direction == "out" else 0.0
        # A 3-D view carries its axes inside the picture, so a figure made
        # only of those reserves nothing for a frame's ticks and labels.
        scenes_only = bool(self.panels) and all(
            getattr(panel.data, "kind", "") == "scene3d" for panel in self.panels)

        target = self.tick_targets(*(panel_size or (40.0, 30.0)))[1]
        widest = 0.0
        for index, panel in enumerate(self.panels):
            if not self.shows_y_labels(index) or not panel.data.is_plot():
                continue
            (_x, _x1), (y0, y1) = self.panel_ranges(panel)
            values, step = tick_values(y0, y1, target)
            fmt = tick_format(step)
            for value in values:
                widest = max(widest, text_width_mm(fmt % value, tick_pt))
        if widest == 0.0 and any(p.data.kind == "empty" for p in self.panels):
            # Slots still to be filled: room for numbers like "-0.25".
            widest = text_width_mm("-0.25", tick_pt)
        # A font's line box is about 1.17 em tall and the axis label sits a
        # line below the tick labels, so the reserve is counted in line
        # boxes rather than em -- get this wrong by 10% and the label is
        # clipped off the bottom of the figure, which is how it was found.
        left = outward + gap + widest + 1.7 * pt_to_mm(label_pt)
        bottom = (outward + gap + 1.25 * pt_to_mm(tick_pt)
                  + 1.9 * pt_to_mm(label_pt))
        top = outward
        right = outward
        if scenes_only:
            left = bottom = top = right = 0.0
        if any(panel.title for panel in self.panels):
            top += 1.4 * pt_to_mm(label_pt)
        if any(panel.top_axis and panel.top_axis.enabled for panel in self.panels):
            top += 1.0 * pt_to_mm(tick_pt) + 1.35 * pt_to_mm(label_pt) + gap
        if any(panel.right_axis and panel.right_axis.enabled for panel in self.panels):
            right += 2.0 * pt_to_mm(tick_pt) + 1.9 * pt_to_mm(label_pt) + gap
        if any(panel.colorbar for panel in self.panels):
            # A rotated label is drawn on its baseline, so it needs a line
            # box either side of it, not an em -- the descender of a "y"
            # falls off the edge of the figure otherwise.
            right += (style.colorbar_gap_mm + style.colorbar_mm
                      + 1.5 * pt_to_mm(tick_pt)
                      + (1.9 * pt_to_mm(label_pt)
                         if any(panel.colorbar_label for panel in self.panels) else 0.0))
        if self.style.panel_label_position and not self.style.panel_label_inside:
            top += 1.4 * pt_to_mm(style.panel_label_pt or style.font_pt + 1)
        return left, bottom, top, right

    #: what a panel looks like when nothing says otherwise
    DEFAULT_ASPECT = 0.75

    def data_aspect(self):
        """Height/width for a panel.

        The ratio of the data's own spans is used **only** when an equal
        aspect was asked for. Otherwise it is meaningless -- degrees against
        eV would make a panel 35 times wider than tall -- so a plain, chosen
        shape is the honest default, the same one for every panel so the
        grid is regular.
        """
        if self.style.panel_aspect > 0:
            return float(self.style.panel_aspect)
        if not self.style.equal_aspect or not self.panels:
            return self.DEFAULT_ASPECT
        (x0, x1), (y0, y1) = self.panel_ranges(self.panels[0])
        span_x, span_y = abs(x1 - x0), abs(y1 - y0)
        if span_x <= 0 or span_y <= 0:
            return self.DEFAULT_ASPECT
        return span_y / span_x

    def _panel_size(self, left, bottom, top, right):
        """The data rectangle these margins would give, for feeding back
        into the tick-count estimate."""
        style = self.style
        inner_w = (style.width_mm - style.margin_left_mm - style.margin_right_mm
                   - (self.cols - 1) * style.gap_x_mm)
        data_w = max(inner_w / self.cols - left - right, 2.0)
        if style.height_mm > 0:
            inner_h = (style.height_mm - style.margin_top_mm
                       - style.margin_bottom_mm - (self.rows - 1) * style.gap_y_mm)
            data_h = max(inner_h / self.rows - top - bottom, 2.0)
        else:
            data_h = max(data_w * self.data_aspect(), 2.0)
        return data_w, data_h

    def layout(self):
        """Every panel's data rectangle in millimetres, plus the figure size.

        Returns ``(width_mm, height_mm, rects)`` with one
        ``{"index", "row", "col", "x", "y", "w", "h"}`` per panel, ``y``
        measured from the top so the painter can use it directly.
        """
        style = self.style
        if self.free:
            # Drawn in order: the panels, then the insets on top of them.
            order = ([i for i, p in enumerate(self.panels) if p.role != "inset"]
                     + [i for i, p in enumerate(self.panels) if p.role == "inset"])
            rects = []
            for index in order:
                x, y, w, h = self.panels[index].frame
                rects.append({"index": index, "row": 0, "col": 0,
                              "x": float(x), "y": float(y),
                              "w": max(float(w), 1.0), "h": max(float(h), 1.0)})
            height = style.height_mm
            if height <= 0:
                height = max(r["y"] + r["h"] for r in rects) + 12.0
            return max(float(style.width_mm), 20.0), float(height), rects
        # Two passes: the margins depend on how wide the tick labels are,
        # which depends on how many there are, which depends on the panel
        # size the margins decide. One round of feedback settles it.
        left, bottom, top, right = self.margins_mm()
        for _ in range(2):
            size = self._panel_size(left, bottom, top, right)
            left, bottom, top, right = self.margins_mm(size)
        width = max(float(style.width_mm), 20.0)

        inner_w = (width - style.margin_left_mm - style.margin_right_mm
                   - (self.cols - 1) * style.gap_x_mm)
        cell_w = max(inner_w / self.cols, 4.0)
        data_w = max(cell_w - left - right, 2.0)

        if style.height_mm > 0:
            height = float(style.height_mm)
            inner_h = (height - style.margin_top_mm - style.margin_bottom_mm
                       - (self.rows - 1) * style.gap_y_mm)
            cell_h = max(inner_h / self.rows, 4.0)
            data_h = max(cell_h - top - bottom, 2.0)
        else:
            data_h = max(data_w * self.data_aspect(), 2.0)
            cell_h = data_h + top + bottom
            height = (cell_h * self.rows + (self.rows - 1) * style.gap_y_mm
                      + style.margin_top_mm + style.margin_bottom_mm)

        rects = []
        for index in range(len(self.panels)):
            row, col = self.cell(index)
            x = (style.margin_left_mm + col * (cell_w + style.gap_x_mm) + left)
            y = (style.margin_top_mm + row * (cell_h + style.gap_y_mm) + top)
            rects.append({"index": index, "row": row, "col": col,
                          "x": x, "y": y, "w": data_w, "h": data_h})
        return width, height, rects

    # -- saving the look, not the data ------------------------------------
    def style_dict(self) -> dict:
        """The figure's appearance, without any measurement in it.

        This is what makes every figure in a paper match: save it once,
        apply it to the next figure. Panel *contents* are deliberately left
        out -- only the settings that should travel.
        """
        return {
            "version": 1,
            "rows": self.rows,
            "cols": self.cols,
            "title": self.title,
            "style": asdict(self.style),
        }

    def apply_style_dict(self, payload: dict, geometry: bool = False):
        """Take the appearance back out. ``geometry`` also restores the grid
        shape, which is usually *not* wanted -- the next figure has its own
        number of panels."""
        style = dict(payload.get("style") or {})
        known = {f for f in FigureStyle().__dict__}
        for key, value in style.items():
            if key in known:
                setattr(self.style, key, value)
        if geometry:
            self.rows = max(int(payload.get("rows", self.rows)), 1)
            self.cols = max(int(payload.get("cols", self.cols)), 1)
        return self

    def save_style(self, path: str):
        with open(path, "w", encoding="utf-8") as handle:
            json.dump(self.style_dict(), handle, indent=2)
        return path

    def load_style(self, path: str, geometry: bool = False):
        with open(path, "r", encoding="utf-8") as handle:
            return self.apply_style_dict(json.load(handle), geometry=geometry)


# --------------------------------------------------------------------------
# Scale bars
# --------------------------------------------------------------------------
def nice_bar_length(span: float) -> float:
    """A round bar about a fifth of the panel wide: 1, 2, 5 x 10^n."""
    if span <= 0 or not np.isfinite(span):
        return 1.0
    raw = span / 5.0
    magnitude = 10.0 ** math.floor(math.log10(raw))
    for multiple in (1.0, 2.0, 5.0, 10.0):
        if raw <= multiple * magnitude * (1 + 1e-9):
            return multiple * magnitude
    return 10.0 * magnitude


def unit_of(label: str) -> str:
    """The unit out of an axis label: "X (µm)" -> "µm"."""
    label = str(label or "")
    if "(" in label and label.rstrip().endswith(")"):
        return label[label.rfind("(") + 1:label.rstrip().rfind(")")]
    return ""
