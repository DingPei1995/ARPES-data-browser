"""
ui/widgets.py
=====================
Qt widgets for browsing SOLEIL nano-ARPES ``*.nxs`` files, playing the same
role ``h5Loader_Widgets.py`` played for the old ``*.h5`` format.

What each data ``kind`` (see loader.nxs_file.NxsScan.kind) gets:

- ``spem_4d`` / ``spem_1d`` (real-space scan): :class:`SpatialImageView`
  (spatial overview, drag the red target to move the readout pixel) next to
  a :class:`FrameImageView` showing the **E-vs-k spectrum** there. Both
  panels support a drag-out rectangle selection: select a box on either
  one, right-click -> "Integrate selection into the other panel", and the
  opposite panel shows the sum over that box.
- ``cut``: a single :class:`FrameImageView` of the one E-vs-k spectrum.
- ``map``: :class:`KCubeExplorer`, three linked orthogonal slices through
  the (deflector angle, slit angle, E) cube, each with its own integration
  controls directly above it.

Shared interaction, on every image panel:
- right-click -> "Readout cursor" puts a crosshair you can drag, which
  reports the data coordinates and the value under it;
- right-click -> "Selection box" adds/removes the rectangle ROI;
- the stock pyqtgraph entries ("Mouse Mode", "X axis", "Y axis" and the
  PlotItem's "Plot Options") are removed from the right-click menu, which
  therefore holds only this app's own actions; what was worth keeping from
  them -- axis range, invert, grid, image alpha -- lives in the visible
  :class:`ViewOptionsBar` under each window's colormap (see
  :func:`strip_stock_menu`).

Colormaps (including the lab's Colormap.mat tables and MATLAB's
gray/jet/hsv/hot/copper) come from the generated ``tools.colormaps`` module;
:func:`apply_colormap` also handles the flip/reverse option.

All widgets are self-contained PyQt5/pyqtgraph -- no dependency on the lab's
``pyNanoScanSystem``/``tools_packages`` packages.
"""
from __future__ import annotations

import os
import weakref
import warnings

import numpy as np
import pyqtgraph as pg
from PyQt5.QtCore import Qt, QObject, QTimer, pyqtSignal
from PyQt5.QtGui import QColor, QFont, QFontMetricsF, QImage, QPainter, QPen
from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
                              QLabel, QSlider, QDoubleSpinBox, QCheckBox,
                              QPushButton, QComboBox, QFrame, QSizePolicy,
                              QScrollArea, QSpinBox, QDialog, QTabWidget,
                              QFormLayout, QLineEdit, QFileDialog, QMessageBox,
                              QMenu, QAction)

from tools import colormaps
from tools import export
from loader.nxs_file import (load_soleil_nxs, to_kspace_cube, NxsScan,
                             save_dataset, CUBE_KINDS)

pg.setConfigOption("background", (245, 241, 249))
pg.setConfigOption("foreground", (102, 126, 161))
pg.setConfigOption("antialias", True)
# Every image in the program is handed to pyqtgraph as (row, column) =
# (y, x). This is global state, so it is set once here, at import, rather
# than from a widget's constructor: set there, an image built before the
# first viewer happened to exist would be drawn with its axes exchanged,
# and which of the two it got would depend on the order the user opened
# windows in.
pg.setConfigOptions(imageAxisOrder="row-major")


def _register_view_without_name(self, name=None):
    """Stand-in for ``pg.ViewBox.register``: every view is still listed, but
    none is registered under a name.

    pyqtgraph names a view so that other views can link axes to it from the
    context menu, a feature this program never offers. The name costs a
    ``destroyed`` callback that, when the view goes, walks *every* other
    view to update its link menu. Every image viewer here is an
    ``ImageView`` and gets the same name, so two viewers closed and then
    reclaimed together -- a map and its slit cut, or any two windows freed
    in one garbage collection -- run that walk over views already half torn
    down, and the interpreter crashes (a segmentation fault, no traceback).
    """
    pg.ViewBox.AllViews[self] = None
    self.name = None


pg.ViewBox.register = _register_view_without_name

RED_PEN = pg.mkPen("#bd4921", width=3)
EDC_COLOR = "#bd4921"   # the EDC curve and the vertical line that feeds it
MDC_COLOR = "#0d7377"   # the MDC curve and the horizontal line that feeds it

#: One colour per axis of a map, used by every readout cursor: a line at
#: constant deflector angle is red wherever it is drawn, one at constant slit
#: angle green, one at constant energy blue -- so the same coordinate is the
#: same colour in the contour and in both of its cuts, and the EDC or MDC a
#: line feeds is drawn in that line's colour. A lone cut (slit angle against
#: energy) uses the slit and energy colours.
AXIS_COLORS = {"defl": "#d62728", "slit": "#2ca02c", "energy": "#1f77b4"}
#: Opacity (0-255) of the shaded EDC/MDC integration windows around a cursor.
BAND_ALPHA = 60


def _rgba(color: str, alpha: int):
    qcolor = QColor(color)
    return (qcolor.red(), qcolor.green(), qcolor.blue(), int(alpha))
# Selection rectangle: deliberately heavy so the box and its corner grips
# stand out against dense ARPES data, and easy to grab.
ROI_PEN = pg.mkPen("#1f6f8b", width=4)
ROI_HOVER_PEN = pg.mkPen("#2ca6cd", width=5)
ROI_HANDLE_PEN = pg.mkPen("#1f6f8b", width=3)
ROI_HANDLE_HOVER_PEN = pg.mkPen("#f5a623", width=4)
ROI_HANDLE_SIZE = 14   # pyqtgraph's default is 5 px, too small to hit reliably
CURSOR_LABEL_OPTS = {
    "position": 0.1,
    "color": (255, 200, 200),
    "fill": (63, 24, 11, 150),
}

COLORMAP_NAMES = colormaps.COLORMAP_NAMES
DEFAULT_COLORMAP = "gray"


def apply_colormap(view, name: str, flip: bool = False) -> bool:
    """Apply a named colormap (see ``colormaps.COLORMAP_NAMES``) to any
    pg.ImageView-based widget, optionally reversed. Returns False silently
    if the name is unknown or this pyqtgraph version rejects it, so callers
    can fire-and-forget across several views."""
    try:
        lut = colormaps.get_lut(name, flip=flip)
        positions = np.linspace(0.0, 1.0, len(lut))
        colors = np.column_stack([lut, np.full(len(lut), 255, dtype=np.uint8)])
        view.setColorMap(pg.ColorMap(positions, colors))
        # Remembered so an export can re-create exactly this colouring at a
        # different pixel size; pyqtgraph keeps the ColorMap but not the name
        # it came from, and the export needs the table, not the widget.
        view._colormap_name, view._colormap_flip = name, bool(flip)
        return True
    except Exception:
        return False


def plain_image_view(x_label: str = "", y_label: str = "") -> pg.ImageView:
    """A bare ``pg.ImageView`` set up the way every picture here is drawn.

    The panels that are not built on :class:`_InteractiveImageBase` -- the
    processing previews, the comparison's third panel, the fit panel's cut
    -- each used to construct their own ``pg.ImageView()`` and each got a
    different subset of this right. Two of pyqtgraph's defaults are wrong
    for a physics plot and neither announces itself:

    * ``ImageView`` inverts the y axis, because its default subject is a
      photograph, whose first row is the top one. On a cut of energy
      against momentum that puts the Fermi level at the bottom and the
      picture upside down.
    * ``ImageView`` also locks the aspect ratio, which squashes any plot
      whose two axes carry unrelated units (eV against A^-1).

    Building one here, once, is what keeps those decisions from being made
    three times. Use :func:`show_frame` to put data into it.
    """
    view = pg.ImageView(view=pg.PlotItem())
    view.ui.roiBtn.hide()
    view.ui.menuBtn.hide()
    view.ui.histogram.hide()
    strip_stock_menu(view.view)
    view.view.invertY(False)
    view.view.setAspectLocked(False)
    view.view.setDefaultPadding(0.0)
    if x_label:
        view.view.setLabel("bottom", x_label)
    if y_label:
        view.view.setLabel("left", y_label)
    return view


def show_frame(view, values, x_axis, y_axis, *, levels=None,
               auto_range: bool = False, auto_levels: bool = False):
    """Draw an ``(x, y)``-indexed array on an ImageView, in data coordinates.

    Every array in this program is indexed ``[ix, iy]`` -- momentum down the
    rows, energy across the columns -- because that is the order the axes
    come out of the file in. pyqtgraph is configured to read an image as
    ``(row, column) = (y, x)``, so the array is transposed **here**, on the
    way to the screen, and nowhere else: transposing the stored array
    instead would silently move every cursor, marker and fitted centre onto
    the wrong axis while the picture looked right.

    ``pos`` and ``scale`` place the image so each pixel is *centred* on its
    axis value rather than starting at it, which is what makes a readout at
    a pixel's centre report that pixel's coordinate.
    """
    values = np.asarray(values, dtype=float)
    x = np.asarray(x_axis, dtype=float).ravel()
    y = np.asarray(y_axis, dtype=float).ravel()
    dx = (x[-1] - x[0]) / max(x.size - 1, 1) if x.size else 1.0
    dy = (y[-1] - y[0]) / max(y.size - 1, 1) if y.size else 1.0
    view.setImage(values.T, autoRange=False, autoLevels=auto_levels,
                  levels=levels,
                  pos=(float(x[0] - dx / 2), float(y[0] - dy / 2)),
                  scale=(float(dx) or 1.0, float(dy) or 1.0))
    if auto_range:
        fit_frame_view(view, x, y)
    return float(dx) or 1.0, float(dy) or 1.0


def fit_frame_view(view, x_axis, y_axis) -> None:
    """Put the whole frame in the box, by explicit range.

    Not ``autoRange``: these views carry infinite cursor lines, which an
    automatic range has to treat as items without bounds, and a view that
    has not been laid out yet has no size to compute a padding against.
    Between them the picture arrives flattened.
    """
    x = np.asarray(x_axis, dtype=float).ravel()
    y = np.asarray(y_axis, dtype=float).ravel()
    if x.size == 0 or y.size == 0:
        return
    dx = abs((x[-1] - x[0]) / max(x.size - 1, 1))
    dy = abs((y[-1] - y[0]) / max(y.size - 1, 1))
    view.view.setXRange(float(x.min() - dx / 2), float(x.max() + dx / 2),
                        padding=0)
    view.view.setYRange(float(y.min() - dy / 2), float(y.max() + dy / 2),
                        padding=0)


#: stock right-click submenus removed from every plot here. "Mouse Mode"
#: only toggles pan-vs-rubber-band dragging and fights the rectangle
#: selection; "X axis"/"Y axis" are replaced by the visible
#: :class:`ViewOptionsBar`, which offers the same range/invert controls
#: without the options this app has no use for (visible-data-only, auto-pan,
#: axis linking, and a mouse-enabled switch that is now always on).
STOCK_SUBMENUS = ("mouse mode", "x axis", "y axis")


def strip_stock_menu(plot_item) -> None:
    """Cut pyqtgraph's stock entries out of a plot's right-click menu,
    leaving room for this app's own actions.

    Removed: the "X axis"/"Y axis" submenus and "Mouse Mode" from the
    ViewBox menu, and the PlotItem's whole "Plot Options" menu. Everything
    worth keeping from those (range, invert, grid, alpha) now lives in the
    :class:`ViewOptionsBar` under the window's colormap, where it is visible
    instead of buried three levels into a context menu.

    Takes the PlotItem when there is one (so "Plot Options" can go too) or a
    bare ViewBox. Handles both menu layouts: pyqtgraph <= 0.13 exposed the
    mouse-mode submenu as ``ViewBoxMenu.leftMenu``, while 0.14 dropped that
    attribute and only leaves the submenu action in place, so it has to be
    matched by title.
    """
    view_box = getattr(plot_item, "vb", plot_item)
    try:
        # "Plot Options" (transforms, downsample, average, alpha, grid,
        # points) belongs to the PlotItem, not the ViewBox. Disabling only
        # the PlotItem's menu leaves the ViewBox's -- and therefore this
        # app's own actions on it -- untouched.
        if hasattr(plot_item, "setMenuEnabled") and hasattr(plot_item, "vb"):
            plot_item.setMenuEnabled(False, None)
    except Exception:
        pass
    try:
        # The mouse is always draggable now, so the switch that could turn
        # it off is gone rather than left in a state nothing exposes.
        view_box.setMouseEnabled(x=True, y=True)
    except Exception:
        pass
    try:
        # pyqtgraph's scene appends its own "Export..." (its exporter dialog)
        # to every context menu. This program has its own Export submenu,
        # which writes the three things a figure is actually built from, so
        # the stock one is only a second door to a worse answer.
        scene = plot_item.scene()
        if scene is not None:
            scene.contextMenu = []
    except Exception:
        pass
    try:
        menu = view_box.menu
        left_menu = getattr(menu, "leftMenu", None)
        if left_menu is not None:
            menu.removeAction(left_menu.menuAction())
        for action in list(menu.actions()):
            title = action.text().replace("&", "").strip().lower()
            if action.menu() is not None and title in STOCK_SUBMENUS:
                menu.removeAction(action)
    except Exception:
        pass



def _menu_action(menu, text: str) -> "QAction":
    """``menu.addAction(text)``, but with the action made from Python.

    ``QMenu.addAction(text)`` has Qt create the action, and PyQt is not told
    when an object Qt created is deleted; a reference kept to one (every
    image view keeps its right-click entries) can outlive it. An action made
    here is PyQt's own and its deletion is tracked. Defensive only: the
    random crashes on opening and closing windows had another cause (see
    ui/gcguard.py).
    """
    action = QAction(text, menu)
    menu.addAction(action)
    return action



def add_button(box, text: str, role):
    """``box.addButton(text, role)`` with the button made from Python, for
    the same reason as :func:`_menu_action`."""
    button = QPushButton(text)
    box.addButton(button, role)
    return button


class AspectBox(QWidget):
    """Holds a plot and gives it the *shape* the data needs.

    This is the answer to a problem that has no good solution inside a
    ViewBox. An equal-aspect image in a panel that is not the image's shape
    can be handled two ways: stretch the picture, or widen the view range
    past the data. pyqtgraph does the second, which is why a square map in a
    wide window ends up with an axis running to +-40 degrees when the data
    stops at 12 -- the axes no longer sit against the data, and every export
    taken from that view carries the empty margins with it.

    The third way, and the right one for a measurement, is to leave the view
    range exactly on the data and make the *plot box* the shape the data is.
    So this widget sizes its child -- allowing for whatever the axes and
    their labels consume inside it -- until the plotting area has the pixel
    aspect the data asks for, and centres it. The empty space then lands
    outside the frame, where it belongs, and the axes hug the data at every
    window size. With the range no longer being padded, pyqtgraph's own
    aspect lock is switched off: the geometry already guarantees the scale.

    ``child`` is the widget that gets resized and ``view`` the image whose
    plotting area must end up the right shape. They are usually the same
    widget, but not on a panel that also carries EDC and MDC curves: there
    the image's edges are pinned to the curves beside it, so what gets
    resized is the whole block of three and the image's share of it comes
    out right because the layout divides the block in fixed proportions.
    Either way the arithmetic is the same -- everything between the block's
    edge and the image's plotting area is "margin".

    ``set_target(None)`` goes back to filling the whole box, which is what
    an unlocked ratio means.
    """

    #: corrections allowed per resize (see :meth:`_verify`)
    BUDGET = 4
    #: the plotting area never goes below this, in pixels
    MIN_PLOT = 48.0

    def __init__(self, child, view, parent=None):
        super().__init__(parent)
        child.setParent(self)
        self.child = child
        self.view = view
        self._target = None          # wanted inner height / inner width
        self._retries = 0
        self._busy = False
        self._clamped = False
        self._effective = None       # the wanted shape, capped to what fits
        self._margins = (60.0, 44.0)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.setMinimumSize(80, 80)

    def plot_margins(self):
        """Pixels spent on everything that is not the image's plotting area
        -- its axes and labels, and on a curves panel the EDC and MDC too.
        What is left is what the data is drawn in.

        Remembered between calls, and a reading is only believed when it is
        plausible. pyqtgraph re-lays its axes out on the next turn of the
        event loop rather than inside ``setGeometry``, so asking twice within
        one pass gives the old ViewBox size against the new widget size -- a
        difference that can come out negative. The last sane pair is a far
        better answer than a fresh nonsensical one.
        """
        try:
            box = self.view.view_box
            margin_x = float(self.child.width() - box.width())
            margin_y = float(self.child.height() - box.height())
            if (0 < margin_x < 0.95 * self.child.width()
                    and 0 < margin_y < 0.95 * self.child.height()):
                self._margins = (margin_x, margin_y)
        except Exception:
            pass
        return self._margins

    def set_target(self, ratio):
        ratio = None if not ratio or not np.isfinite(ratio) or ratio <= 0 else float(ratio)
        if ratio is not None and self._target is not None \
                and abs(ratio - self._target) < 1e-4:
            return
        self._target = ratio
        self._retries = self.BUDGET
        self.relayout()

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self._retries = self.BUDGET
        self.relayout()

    def recheck(self):
        """The child's plotting area changed size without us asking -- which
        happens whenever pyqtgraph decides its tick labels need a wider or
        narrower axis. If that has knocked the shape out, spend one of the
        corrections this resize is allowed."""
        self._verify()

    def relayout(self):
        """Re-solve the child's geometry, once.

        Deliberately one pass and no loop. The margins are read from the
        child's *current* layout, and pyqtgraph only re-lays its axes out on
        the next turn of the event loop -- so a second pass inside the same
        call would solve against a widget that has already been resized and a
        ViewBox that has not, and land further from the answer than it
        started. The iterating is done by :meth:`_verify` instead, one pass
        per event-loop turn, where every reading is real.
        """
        width, height = self.width(), self.height()
        if width <= 0 or height <= 0 or self._busy:
            return
        self._busy = True
        try:
            if self._target is None:
                self.child.setGeometry(0, 0, width, height)
                return
            margin_x, margin_y = self.plot_margins()
            inner_w = max(width - margin_x, self.MIN_PLOT)
            inner_h = max(height - margin_y, self.MIN_PLOT)

            # A ratio can ask for a shape this box cannot hold: 1:1 between
            # 28 degrees and 0.8 eV wants a strip 35 times wider than tall,
            # which on a panel carrying EDC and MDC curves leaves the image
            # nothing at all. The wanted shape is therefore capped to what
            # still leaves a usable plot, and the cap -- not the original --
            # is what the geometry then solves for, so it converges instead
            # of chasing an impossible target. The axes still sit exactly on
            # the data, which is the guarantee that matters; what is given
            # up is the last of the scale, on a setting that was already
            # asking for something unreadable.
            thinnest, tallest = self.MIN_PLOT / inner_w, inner_h / self.MIN_PLOT
            target = min(max(self._target, thinnest), tallest)
            self._clamped = abs(target - self._target) > 1e-9 * max(target, 1.0)
            self._effective = target

            if inner_w * target <= inner_h:
                inner_h = inner_w * target
            else:
                inner_w = inner_h / target
            child_w = int(round(inner_w + margin_x))
            child_h = int(round(inner_h + margin_y))
            self.child.setGeometry((width - child_w) // 2, (height - child_h) // 2,
                                   child_w, child_h)
        finally:
            self._busy = False
        self._verify()

    def _verify(self):
        """Check the plotting area really came out the shape asked for, and
        if not, try again after the event loop has let pyqtgraph finish
        laying its axes out.

        A correction budget, rather than a loop, is what makes this safe: a
        pass is spent only when the shape is actually wrong, and the budget
        is refilled only by a real resize. A good result leaves the budget
        alone, because the axis widths can still change afterwards -- which
        is exactly the case this exists for, a panel that comes out right and
        is then knocked out of shape when the tick labels get shorter.
        """
        if self._target is None or self._busy:
            return
        box = getattr(self.view, "view_box", None)
        if box is None or box.width() < 1 or box.height() < 1:
            return
        # Against the capped shape, not the asked-for one: an impossible
        # ratio would otherwise spend every correction and never arrive.
        wanted = self._effective or self._target
        achieved = float(box.height()) / float(box.width())
        if abs(achieved - wanted) <= 0.01 * wanted or self._retries <= 0:
            return
        self._retries -= 1
        QTimer.singleShot(0, self.relayout)


def strip_mouse_mode(view_box) -> None:
    """Backwards-compatible alias for :func:`strip_stock_menu`."""
    strip_stock_menu(view_box)


def _convolve1d(array: np.ndarray, kernel: np.ndarray, axis: int) -> np.ndarray:
    """Separable 1D convolution with edge padding, using only numpy (scipy
    is not a dependency of this app)."""
    pad = len(kernel) // 2
    widths = [(pad, pad) if i == axis else (0, 0) for i in range(array.ndim)]
    padded = np.pad(array, widths, mode="edge")
    out = np.zeros_like(array, dtype=float)
    for offset, weight in enumerate(kernel):
        sl = [slice(None)] * array.ndim
        sl[axis] = slice(offset, offset + array.shape[axis])
        out += weight * padded[tuple(sl)]
    return out


#: A displayed image is resampled up to about this many samples along each
#: axis before it is drawn, which is what turns a coarse scan from a wall of
#: rectangles into a continuous picture.
INTERP_TARGET = 500
#: ... but never by more than this, so a tiny array cannot blow up into a
#: huge texture.
MAX_INTERP_FACTOR = 16


def interp_factor(n: int, target: int = INTERP_TARGET) -> int:
    """How many display samples to put between each pair of data points.

    Data that is already dense gets 1 (no resampling): there is nothing to
    gain from interpolating a 500-point axis onto 500 points.
    """
    if n < 2:
        return 1
    return int(np.clip(int(np.ceil(target / n)), 1, MAX_INTERP_FACTOR))


def _upsample_axis(array: np.ndarray, factor: int, axis: int) -> np.ndarray:
    """Linear interpolation along one axis onto ``factor`` times as many
    samples, covering exactly the same extent.

    Sample j of the output sits at index ``-0.5 + (j + 0.5)/factor`` of the
    input, i.e. the output pixels tile the input's outer bounds (each input
    pixel is centred on its axis value and spans half a step either side).
    That keeps the image in exactly the same place on the axes, so only the
    ``scale`` passed to setImage changes.
    """
    if factor <= 1:
        return array
    n = array.shape[axis]
    pos = np.clip(-0.5 + (np.arange(n * factor) + 0.5) / factor, 0, n - 1)
    i0 = np.floor(pos).astype(int)
    i1 = np.minimum(i0 + 1, n - 1)
    w = (pos - i0).reshape([-1 if d == axis else 1 for d in range(array.ndim)])
    return np.take(array, i0, axis=axis) * (1 - w) + np.take(array, i1, axis=axis) * w


def upsample2d(array: np.ndarray, factors) -> np.ndarray:
    """Bilinear resampling of a 2D array, NaN-aware.

    This is what MATLAB's ``pcolor(...); shading interp`` does for the same
    data, and it is the reason the MATLAB plots look continuous where a
    plain image of the same scan looks like a mosaic: a deflector map is
    typically only a few tens of points across, so each data point covers a
    big block of screen. Blurring alone cannot fix that -- a Gaussian on the
    coarse grid just makes softer *blocks* -- so the array is resampled onto
    a finer grid first.

    Missing points are carried by the same weights as the values, so a NaN
    pulls in its neighbours instead of spreading a hole the size of the
    interpolation kernel.
    """
    fy, fx = (int(max(1, f)) for f in factors)
    if fy <= 1 and fx <= 1:
        return array
    data = np.asarray(array, dtype=float)
    finite = np.isfinite(data)
    if finite.all():
        out = _upsample_axis(_upsample_axis(data, fy, 0), fx, 1)
        return out
    filled = np.where(finite, data, 0.0)
    num = _upsample_axis(_upsample_axis(filled, fy, 0), fx, 1)
    den = _upsample_axis(_upsample_axis(finite.astype(float), fy, 0), fx, 1)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = num / den
    return np.where(den > 1e-9, out, np.nan)


def smooth2d(array: np.ndarray, sigma: float = 1.0) -> np.ndarray:
    """Gaussian-smooth a 2D array for *display*, to take the staircase edges
    off sparsely-sampled scans.

    NaN-aware: missing points are excluded from both the weighted sum and
    its normalisation, so a gap pulls in its neighbours' values instead of
    poisoning the whole neighbourhood (plain convolution would spread the
    NaN). Pixels with no finite neighbour at all stay NaN.

    This never touches stored data -- callers keep the raw array and smooth
    only what they hand to setImage, so exports stay unsmoothed.
    """
    data = np.asarray(array, dtype=float)
    if sigma <= 0 or data.ndim != 2:
        return data
    radius = max(1, int(np.ceil(2 * sigma)))
    offsets = np.arange(-radius, radius + 1)
    kernel = np.exp(-0.5 * (offsets / sigma) ** 2)
    kernel /= kernel.sum()

    finite = np.isfinite(data)
    filled = np.where(finite, data, 0.0)
    numerator = _convolve1d(_convolve1d(filled, kernel, 0), kernel, 1)
    denominator = _convolve1d(_convolve1d(finite.astype(float), kernel, 0), kernel, 1)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = numerator / denominator
    return np.where(denominator > 0, out, np.nan)


# --------------------------------------------------------------------------
# Data wrapper
# --------------------------------------------------------------------------
def _options_key(options):
    """A hashable identity for a set of load options, for the shared-handle
    registry. None and the default options are the same thing."""
    if options is None:
        return None
    key = (getattr(options, "loader", None),
           tuple(options.permutation) if getattr(options, "permutation", None) else None,
           getattr(options, "axis0_role", "angle") or "angle",
           getattr(options, "axis0_label", "") or "")
    return None if key == (None, None, "angle", "") else key


class NxsData:
    """Loads one data file through the loader registry (:mod:`loader.registry`)
    and exposes the GUI-friendly derived views on top of it. ``self.scan`` is
    the raw :class:`~loader.nxs_file.NxsScan`."""

    #: path -> the single open NxsData for it (see acquire()).
    _open_files: dict = {}

    def __init__(self, nxs_path: str, entry: str = None, options=None,
                 progress=None):
        self.path = nxs_path
        self.options = options
        # Through the loader registry rather than straight to the SOLEIL
        # reader: which beamline wrote this, and the axis choices made at
        # load time, are settings now (see loader.registry / ui.loader_dialog).
        # Falling back to the old direct call keeps a file opening even if
        # no loader recognises it but the SOLEIL reader can cope.
        from loader import registry
        try:
            self.scan: NxsScan = registry.load(nxs_path, entry=entry,
                                                  options=options,
                                                  progress=progress)
        except ValueError:
            if options is not None and getattr(options, "loader", None):
                raise
            self.scan = load_soleil_nxs(nxs_path, entry=entry)
        self.kind = self.scan.kind
        self._overview = None
        self._key = None
        self._refs = 1

    @classmethod
    def acquire(cls, nxs_path: str, entry: str = None, options=None,
                progress=None) -> "NxsData":
        """Get the shared NxsData for a file, opening it if needed.

        HDF5 keeps one underlying handle per file per process, so two
        independently opened h5py.File objects for the same path are *not*
        independent -- closing either one invalidates both, and a second
        window on the same scan would kill the first window's lazy reads.
        Callers therefore share one instance and each release it with
        :meth:`close`; the file is only really closed when the last user has
        let go.

        The load options are part of the identity: the same file read with
        its axes in a different order is a different dataset, and handing
        back the cached one would quietly ignore what was asked for.
        """
        key = (os.path.abspath(nxs_path), entry, _options_key(options))
        existing = cls._open_files.get(key)
        if existing is not None and existing.alive():
            existing._refs += 1
            return existing
        if existing is not None:
            # Its file handle died under it. Handing it out again is the
            # "identifier is not of specified type" crash on the first read;
            # drop it and open the file afresh instead.
            cls._open_files.pop(key, None)
            existing._key = None
        data = cls(nxs_path, entry, options=options, progress=progress)
        data._key = key
        cls._open_files[key] = data
        return data

    def retain(self):
        """Take one more reference (released by :meth:`close`) and return
        self. For a holder that did not get this object from
        :meth:`acquire` -- the memory budget handing out what it caches --
        so that every holder owns exactly one reference and closes exactly
        once."""
        self._refs += 1
        return self

    def alive(self) -> bool:
        """Whether this dataset can still be read: an in-memory one always
        can; a lazy one only while its HDF5 dataset is open."""
        for value in (getattr(self.scan, "value", None),
                      getattr(self.scan, "value4d", None)):
            dset = getattr(value, "_dset", None)
            if dset is not None:
                try:
                    if not dset.id.valid:
                        return False
                except Exception:                           # noqa: BLE001
                    return False
        return self._refs > 0

    def close(self):
        """Release one reference; the file closes when the last one goes."""
        self._refs -= 1
        if self._refs > 0:
            return
        if self._key is not None:
            NxsData._open_files.pop(self._key, None)
            self._key = None
        self.scan.close()

    # -- spem_4d --------------------------------------------------------
    def _overview_from_preview(self):
        """Spatial overview from one of the file's reduced preview cubes.

        ``data_11`` is the cube already summed over energy and ``data_01``
        over the slit angle, so summing either one's last axis gives the same
        map as summing the whole cube -- at a thousandth of the data volume,
        which is the difference between opening a scan instantly and reading
        a gigabyte off a network share.

        The result is checked against the real cube on a few pixels before
        being trusted; returns None if no preview is present or the check
        fails, and the caller streams the cube instead.
        """
        for name in ("data_11", "data_01"):
            entry = self.scan.previews.get(name)
            if entry is None:
                continue
            dset, perm = entry
            try:
                preview = np.transpose(dset[()], perm)       # (y, x, other)
                overview = preview.sum(axis=2, dtype=np.float64)
            except Exception as exc:
                warnings.warn(f"{name}: could not be read as a preview ({exc})")
                continue

            ny, nx = overview.shape
            samples = [(0, 0), (ny // 2, nx // 2), (ny - 1, nx - 1)]
            ok = True
            for yi, xi in samples:
                true_sum = float(np.asarray(self.scan.value4d[yi, xi], dtype=np.float64).sum())
                claimed = float(overview[yi, xi])
                scale = max(abs(true_sum), abs(claimed), 1.0)
                if abs(true_sum - claimed) / scale > 1e-3:
                    ok = False
                    break
            if ok:
                return overview
            warnings.warn(
                f"{name}: does not reproduce the cube's own sums, so it is not "
                f"a plain projection in this file; falling back to reading the "
                f"cube for the spatial overview.")
        return None

    def _overview_by_streaming(self, progress=None):
        """Fallback: sum the cube one spatial row at a time, so peak memory
        stays at one row rather than the whole gigabyte."""
        ny, nx = self.scan.value4d.shape[:2]
        overview = np.empty((ny, nx), dtype=np.float64)
        for yi in range(ny):
            row = np.asarray(self.scan.value4d[yi], dtype=np.float64)   # (x, k, E)
            overview[yi] = row.sum(axis=(1, 2))
            if progress is not None:
                progress(yi + 1, ny)
        return overview

    # -- spem_4d --------------------------------------------------------
    def spatial_overview_computed(self, progress=None) -> np.ndarray:
        """(y, x) map, summed over (k, E). Uses the file's reduced preview
        cube when it checks out, and otherwise streams the real cube."""
        assert self.kind == "spem_4d"
        if self._overview is None:
            overview = self._overview_from_preview()
            if overview is None:
                overview = self._overview_by_streaming(progress)
            self._overview = overview
        return self._overview

    @property
    def spatial_overview(self) -> np.ndarray:
        return self.spatial_overview_computed()

    def frame_at(self, yi: int, xi: int) -> np.ndarray:
        """(k, E) spectrum at spatial pixel (yi, xi)."""
        assert self.kind == "spem_4d"
        return self.scan.value4d[yi, xi]

    def frame_over_region(self, yslice, xslice) -> np.ndarray:
        """(k, E) spectrum summed over a rectangular spatial region, read one
        spatial row at a time so a large selection never has to fit in
        memory all at once."""
        assert self.kind == "spem_4d"
        rows = range(*yslice.indices(self.scan.value4d.shape[0]))
        total = None
        for yi in rows:
            block = np.asarray(self.scan.value4d[yi, xslice], dtype=np.float64)
            block = block.sum(axis=0)
            total = block if total is None else total + block
        return total

    def spatial_over_region(self, kslice, eslice, progress=None) -> np.ndarray:
        """(y, x) map summed over a rectangular (k, E) region -- the inverse
        selection: pick a feature in the spectrum, see where on the sample it
        comes from.

        This one genuinely has to touch every spatial pixel, but only inside
        the selected (k, E) window, and it reads a row at a time.
        """
        assert self.kind == "spem_4d"
        ny, nx = self.scan.value4d.shape[:2]
        out = np.empty((ny, nx), dtype=np.float64)
        for yi in range(ny):
            block = np.asarray(self.scan.value4d[yi, :, kslice, eslice], dtype=np.float64)
            out[yi] = block.sum(axis=(1, 2))
            if progress is not None:
                progress(yi + 1, ny)
        return out

    # -- spem_1d ----------------------------------------------------------
    @property
    def line_kmap(self) -> np.ndarray:
        """(x, k) map, summed over E, for a real-space *line* scan."""
        assert self.kind == "spem_1d"
        # np.asarray rather than .sum() straight off scan.value: a scan's
        # array may be a LazyCube/LazyArray, which forwards indexing and
        # np.asarray but not ndarray's own methods.
        return np.asarray(self.scan.value).sum(axis=2)

    def frame_at_x(self, xi: int) -> np.ndarray:
        assert self.kind == "spem_1d"
        return self.scan.value[xi]

    # -- cut ----------------------------------------------------------
    @property
    def cut_frame(self) -> np.ndarray:
        assert self.kind == "cut"
        return self.scan.value

    # -- map (angle/angle/E cube) ---------------------------------------
    @property
    def angle_cube(self):
        """The three axes and the cube, for the Map viewer.

        Named for the angle map it was written for, but a k-map saved to a
        file and read back arrives here too -- same shape, same viewer, axes
        that are momenta and labelled as such. Returns (x, y, E, cube).
        """
        assert self.kind in CUBE_KINDS, self.kind
        return self.scan.x, self.scan.k, self.scan.z, self.scan.value

    def kcube(self, kinetic_energy_eV=None):
        """Optional, NOT used by the GUI by default: converts the
        deflector-angle axis to momentum via loader.nxs_file.to_kspace_cube
        (small-angle free-electron approximation). Returns
        (kx, ky, E, cube). Available for scripted/offline analysis."""
        assert self.kind == "map"
        return to_kspace_cube(self.scan, kinetic_energy_eV=kinetic_energy_eV)


class _MemScan:
    """The subset of :class:`loader.nxs_file.NxsScan` that a computed dataset
    needs to expose.

    A converted k-map, an arbitrary-direction cut and a Fermi-corrected map
    are all computed rather than read: no file, no HDF5 entry, no lazy cube
    behind them. The viewers only ever touch these few attributes, so
    presenting them is enough for a computed dataset to be displayed,
    exported and saved exactly like a measured one.

    Which attributes carry the axes depends on the kind, matching the file
    parser: a 2D ``cut`` uses ``x``/``y``, a 3D ``map``/``k_map`` uses
    ``x``/``k``/``z``.
    """

    def __init__(self, kind, axes, value, labels, info, fourd_info=None):
        from loader.nxs_file import axis_slots

        self.kind = kind
        self.labels = labels
        self.info = info
        self.fourd_info = dict(fourd_info or {})
        self.previews = {}
        self.value = value
        self.value4d = None
        self.x = self.y = self.k = self.z = None

        # Which slots this kind fills comes from loader.nxs_file.AXIS_SLOTS,
        # the one table that says so -- this used to be one of five
        # hand-written copies of it.
        slots = axis_slots(kind, "constructor")
        if not slots:
            slots = ("x", "y")             # unknown 2-D-ish kind: best effort
        if len(axes) != len(slots):
            raise ValueError(
                f"a {kind} takes {len(slots)} axes {slots}, got {len(axes)}")
        for slot, values in zip(slots, axes):
            setattr(self, slot, values)

        if kind == "spem_4d":
            # The cube lives in value4d, where the viewer looks for it, and
            # value mirrors it so the generic paths (saving, the data
            # operations) need no special case.
            self.value4d = value
        elif kind in CUBE_KINDS:
            # A map's second axis is the analyser's, which the viewers reach
            # as either .k or .y depending on which of them is asking.
            self.y = self.k

    def close(self):
        pass


class MemoryData:
    """A dataset computed in this session, quacking like :class:`NxsData`.

    Everything the viewers use -- ``kind``, ``scan``, ``cut_frame`` /
    ``angle_cube``, ``path``, ``close`` -- is here, so a computed dataset
    opens through exactly the same window as the measurement it came from,
    is exported by the same code, and is written by ``save_dataset`` in the
    same format. ``info`` records where it came from and the settings that
    produced it, so a dataset saved now still says how it was made.
    """

    def __init__(self, kind, axes, value, labels, *, source_label,
                 parameters=None, prefix="calc", source_path="",
                 source_info=None, source_motors=None):
        info = dict(source_info or {})
        info["_kind"] = kind
        info["_source"] = source_label
        for key, val in (parameters or {}).items():
            info[f"{prefix}.{key}"] = val

        self.path = source_path
        self.kind = kind
        self.source_label = source_label
        self.parameters = dict(parameters or {})
        self.scan = _MemScan(kind, axes, value, labels, info, source_motors)

    # -- what the viewers ask for -----------------------------------------
    @property
    def cut_frame(self):
        assert self.kind == "cut", self.kind
        return self.scan.value

    # -- spatial scan, for a SPEM computed in this session ------------------
    def spatial_overview_computed(self, progress=None) -> np.ndarray:
        """(y, x) map, summed over (k, E). In memory the whole cube is here
        already, so there is nothing to stream or approximate."""
        assert self.kind == "spem_4d", self.kind
        if getattr(self, "_overview", None) is None:
            self._overview = np.asarray(self.scan.value4d, dtype=np.float64).sum(axis=(2, 3))
        return self._overview

    @property
    def spatial_overview(self) -> np.ndarray:
        return self.spatial_overview_computed()

    def frame_at(self, yi: int, xi: int) -> np.ndarray:
        assert self.kind == "spem_4d", self.kind
        return np.asarray(self.scan.value4d[yi, xi], dtype=float)

    def frame_over_region(self, yslice, xslice) -> np.ndarray:
        assert self.kind == "spem_4d", self.kind
        return np.asarray(self.scan.value4d[yslice, xslice],
                          dtype=np.float64).sum(axis=(0, 1))

    def spatial_over_region(self, kslice, eslice, progress=None) -> np.ndarray:
        assert self.kind == "spem_4d", self.kind
        return np.asarray(self.scan.value4d[:, :, kslice, eslice],
                          dtype=np.float64).sum(axis=(2, 3))

    @property
    def line_kmap(self) -> np.ndarray:
        assert self.kind == "spem_1d", self.kind
        return np.asarray(self.scan.value, dtype=float).sum(axis=2)

    def frame_at_x(self, xi: int) -> np.ndarray:
        assert self.kind == "spem_1d", self.kind
        return np.asarray(self.scan.value[xi], dtype=float)

    @property
    def angle_cube(self):
        """(x, y, E, cube) -- the Map viewer's accessor. Keeps the name for
        a k-map, whose axes are momenta: renaming it would mean forking the
        viewer for no gain, and the axis *labels* say A^-1, which is what
        the user reads."""
        assert self.kind in CUBE_KINDS, self.kind
        return self.scan.x, self.scan.k, self.scan.z, self.scan.value

    def close(self):
        pass

    def retain(self):
        """Nothing to count for an in-memory dataset; see NxsData.retain."""
        return self

    def alive(self) -> bool:
        return True


def KMapData(kx, ky, energy, cube, *, source_label: str, parameters: dict,
             source_path: str = "", source_info: dict = None,
             source_motors: dict = None) -> MemoryData:
    """A Map converted to momentum space, held in memory.

    A thin front for :class:`MemoryData` with the k-map's labels and the
    ``kconv.*`` parameter prefix the conversion round established (the `.m`
    file kept the same record in ``data_k.para``).
    """
    return MemoryData("k_map", (kx, ky, energy), cube,
                      {"x": "kx (\u00c5\u207b\u00b9)", "k": "ky (\u00c5\u207b\u00b9)",
                       "z": "Energy (eV)"},
                      source_label=source_label, parameters=parameters,
                      prefix="kconv", source_path=source_path,
                      source_info=source_info, source_motors=source_motors)


# --------------------------------------------------------------------------
# Shared interaction: readout cursor + rectangle selection
# --------------------------------------------------------------------------
class _InteractiveImageBase(pg.ImageView):
    """Common machinery for every image panel here: axis bookkeeping, an
    optional draggable readout crosshair, an optional rectangle selection,
    and a cleaned-up right-click menu.

    Subclasses call :meth:`_store_axes` whenever they set new image data so
    coordinate<->index conversion stays correct.
    """

    #: emitted when the user asks (via right-click) to apply the current
    #: rectangle selection; the owner decides what that means.
    selectionApplied = pyqtSignal()
    #: emitted whenever the selection rectangle moves or is created/removed,
    #: so the coordinate boxes below the panel can follow it.
    selectionChanged = pyqtSignal()
    #: emitted with the readout cursor's (ix, iy) array indices whenever it
    #: moves, so the EDC/MDC curves under the panel can follow it.
    readoutMoved = pyqtSignal(int, int)
    #: emitted when the readout cursor is switched on or off.
    readoutToggled = pyqtSignal(bool)
    #: emitted when the user asks to snapshot the current position into its
    #: own window (spatial panel only).
    popoutRequested = pyqtSignal()
    #: emitted when the user asks for the two displayed axes to be exchanged
    #: (spatial panel only).
    swapToggled = pyqtSignal(bool)

    def __init__(self, parent=None, name="view"):
        self.plot = pg.PlotItem()
        super().__init__(parent=parent, name=name, view=self.plot)
        self.ui.roiBtn.hide()
        self.ui.menuBtn.hide()
        # Box frame: all four axes drawn, values only on the bottom and left.
        # Without the top/right lines the data just fades out at the edge of
        # the panel, which makes it hard to see where the data actually ends.
        self.plot.showAxes(True, showValues=(True, False, False, True))
        # No padding: the axes sit tight against the data instead of leaving
        # pyqtgraph's default margin of blank space around every image.
        self.view_box.setDefaultPadding(0.0)
        # The panel's LevelBar replaces this: same job, far less screen space,
        # and no intensity distribution to draw.
        self.ui.histogram.hide()
        self.view.invertY(False)

        self.x_label = None   # remembered axis titles (units included)
        self.y_label = None
        self._xarray = None   # horizontal axis values
        self._yarray = None   # vertical axis values

        self.readout_cursor = None
        self.readout_label = None
        #: colours of the vertical (constant x) and horizontal (constant y)
        #: cursor lines and their bands; see AXIS_COLORS
        self._readout_colors = (AXIS_COLORS["slit"], AXIS_COLORS["energy"])
        self.readout_vline = self.readout_hline = None
        self.readout_vband = self.readout_hband = None
        self._band_half = {}        # band -> half-width it was last given
        self.selection_roi = None
        self._value_lookup = None   # callable(ix, iy) -> value, set by subclass
        self._smooth_on = False     # display interpolation + optional blur
        self._smooth_sigma = 0.0    # extra Gaussian blur, in data pixels
        #: how many display samples per data point the last _for_display
        #: produced, per array axis -- setImage's scale divides by these
        self._display_factors = (1, 1)
        self._gamma = 1.0           # display-only intensity exponent
        self._levels = None         # the colour window gamma is applied in
        self._aspect_ratio = None   # locked x:y scale, or None to stretch
        self._aspect_box = None     # the AspectBox shaping this view, if any
        self._fitting = False       # re-entrancy guard for fit_to_data
        self._pick_callback = None  # set while a dialog is collecting points
        self._pick_markers = None
        self._overlay_curve = None
        self._bz_layers = {}
        self._colormap_name, self._colormap_flip = DEFAULT_COLORMAP, False

        # What an export calls this panel and what metadata travels with it.
        # Set by whoever owns the view: ImagePanel supplies the title, the
        # viewer window the file name and the file's metadata. Plain data
        # only -- deliberately no reference back to the window, because a
        # view holding its own window alive turns the whole viewer into a
        # reference cycle, and a cycle is freed at whatever moment the
        # garbage collector chooses rather than when the window closes.
        # Qt objects destroyed at an arbitrary moment crash pyqtgraph.
        self.export_title = ""
        self.export_prefix = "export"
        self.export_info = {}
        self.export_source_path = ""
        self._export_dialog = None

        # Re-fitting on resize must not fight a user who has zoomed or
        # panned deliberately, so track whether they have.
        self._user_ranged = False
        self.view_box.sigRangeChangedManually.connect(self._on_manual_range)
        # The required range depends on the viewbox's pixel size whenever an
        # aspect is locked, and the final size is only known after the layout
        # has settled -- later than any showEvent we could hook. sigResized
        # fires on every geometry change, including that last pass.
        self.view_box.sigResized.connect(lambda *_: self.refit_if_untouched())
        self.view_box.sigResized.connect(self._on_plot_area_resized)
        # The plot box's shape follows the range it is showing, so a zoom
        # reshapes the box instead of re-introducing padding.
        self.view_box.sigRangeChanged.connect(self.update_aspect_box)

        strip_stock_menu(self.plot)
        self._add_context_actions()

    def _on_manual_range(self, *_):
        self._user_ranged = True

    def refit_if_untouched(self):
        """Re-apply the tight fit after a resize. An aspect lock stretches
        whichever axis is short *for the widget's current shape*, so the
        right range changes every time the panel is resized; this keeps the
        image tight and centred as the window changes, but leaves the view
        alone once the user has zoomed or panned."""
        if not self._user_ranged:
            self.fit_to_data()

    @property
    def view_box(self):
        """The ViewBox owning the right-click menu. ``self.view`` is the
        PlotItem we passed to ImageView, not the ViewBox itself, so the menu
        lives one level down in ``.vb``."""
        return getattr(self.view, "vb", self.view)

    # -- display options -------------------------------------------------
    def set_axis_labels(self, x_label: str = None, y_label: str = None):
        """Label the axes and remember the labels, so anything built on top
        of this view (the EDC/MDC curves, pop-out windows) can reuse the
        same names and units instead of re-deriving them."""
        if x_label:
            self.x_label = x_label
            self.view.setLabel("bottom", x_label)
        if y_label:
            self.y_label = y_label
            self.view.setLabel("left", y_label)

    def set_smoothing(self, enabled: bool, sigma: float = 0.0):
        """Turn display smoothing on or off, with ``sigma`` pixels of extra
        Gaussian blur on top of it (0 = interpolation only). Display only:
        the stored arrays and every export stay raw."""
        self._smooth_on = bool(enabled)
        self._smooth_sigma = max(0.0, float(sigma))
        self.redraw()

    def _for_display(self, array: np.ndarray) -> np.ndarray:
        """Turn a data array into what is actually drawn.

        Two steps, in this order because it is much the cheaper one:

        1. an optional Gaussian blur on the *data* grid, where the kernel is
           a few pixels rather than a few hundred;
        2. bilinear resampling onto a finer grid (:func:`upsample2d`), which
           is what removes the mosaic of rectangles a coarse scan otherwise
           draws -- the same thing MATLAB's ``shading interp`` does.

        Records the per-axis resampling factors so the caller can divide the
        image scale by them and leave the picture exactly where it was.
        """
        data = np.asarray(array, dtype=float)
        if not self._smooth_on:
            self._display_factors = (1, 1)
            return self._apply_gamma(data)
        if self._smooth_sigma > 0:
            data = smooth2d(data, self._smooth_sigma)
        factors = tuple(interp_factor(n) for n in data.shape[:2])
        self._display_factors = factors
        return self._apply_gamma(upsample2d(data, factors))

    def _apply_gamma(self, data: np.ndarray) -> np.ndarray:
        """Re-shape the intensities inside the colour window, the way
        ``PlotSlices.m``'s plot() does: clip to the window, normalise, raise
        to gamma, map back. Levels are left alone, so the numbers on the
        level bar still mean what they say."""
        if self._gamma == 1.0 or self._levels is None:
            return data
        lo, hi = self._levels
        if hi <= lo:
            return data
        with np.errstate(invalid="ignore"):
            norm = np.clip((data - lo) / (hi - lo), 0.0, 1.0) ** self._gamma
        return lo + norm * (hi - lo)

    def _reapply_levels(self):
        """setImage resets the levels whenever autoLevels runs; put the
        user's window back so a redraw never changes the contrast."""
        if self._levels is not None:
            self.setLevels(*self._levels)

    def redraw(self):
        """Re-render the current data with the present display options.
        Subclasses override to re-issue their own setImage call."""

    def set_aspect(self, locked: bool, ratio: float = 1.0):
        """Lock the on-screen x:y scale to ``ratio`` data units, or let the
        image stretch freely to fill its box.

        Where the panel gave us an :class:`AspectBox` -- which is every panel
        now, curves or not -- the lock is enforced by the plot box's *shape*,
        so the view range can stay exactly on the data at any window size.
        pyqtgraph's own lock is the fallback for a view nobody has given a
        box to, and it keeps the scale honest at the price of padding the
        range, which is what the box exists to avoid.
        """
        self._aspect_ratio = float(ratio) if (locked and ratio > 0) else None
        if self._aspect_box is not None:
            self.view.setAspectLocked(False)
            self.update_aspect_box()
        elif self._aspect_ratio is not None:
            self.view.setAspectLocked(True, ratio=self._aspect_ratio)
        else:
            self.view.setAspectLocked(False)

    def attach_aspect_box(self, box):
        self._aspect_box = box
        self.set_aspect(self._aspect_ratio is not None, self._aspect_ratio or 1.0)

    def aspect_target(self):
        """Wanted plot-area height/width for the range on screen, or None
        when the image is free to stretch."""
        if self._aspect_ratio is None:
            return None
        (x0, x1), (y0, y1) = self.view_box.viewRange()
        span_x, span_y = abs(x1 - x0), abs(y1 - y0)
        if span_x <= 0 or span_y <= 0:
            return None
        return span_y / (self._aspect_ratio * span_x)

    def update_aspect_box(self, *_):
        if self._aspect_box is not None:
            self._aspect_box.set_target(self.aspect_target())

    def _on_plot_area_resized(self, *_):
        """The plotting area changed size -- which, when we did not ask for
        it, means pyqtgraph has re-measured its axes and the shape may have
        drifted."""
        if self._aspect_box is not None:
            self._aspect_box.recheck()

    def export_aspect(self):
        """Height/width an export of this panel should have.

        With a locked ratio this comes from the data and the ratio, never
        from the window: an export must not change shape because the window
        was dragged wider. Unlocked, the image genuinely is whatever shape
        the panel is, so the panel's shape is the honest answer.
        """
        target = self.aspect_target()
        if target is not None:
            return target
        try:
            box = self.view_box
            if box.width() > 1 and box.height() > 1:
                return float(box.height()) / float(box.width())
        except Exception:
            pass
        return 0.75

    # -- what the ViewOptionsBar drives ---------------------------------
    def set_grid(self, on: bool, alpha: float = 0.25):
        """Draw (or drop) the x/y grid over the image."""
        try:
            self.plot.showGrid(x=bool(on), y=bool(on), alpha=float(alpha))
        except Exception:
            pass

    def set_gamma(self, gamma: float):
        """Intensity exponent applied inside the colour window, as in the
        MATLAB tool: the window is normalised to 0..1, raised to ``gamma``
        and mapped back. Below 1 it lifts weak features out of the
        background (the usual way to see a faint Fermi surface under a
        bright band), above 1 it suppresses them. Display only."""
        gamma = float(gamma)
        if gamma <= 0:
            return
        self._gamma = gamma
        self.redraw()

    def set_level_window(self, lo: float, hi: float):
        """The min/max the colours span. Kept on the view rather than only
        pushed into setLevels because gamma has to be applied *inside* this
        window -- it is what "0" and "1" mean for the exponent."""
        if not (np.isfinite(lo) and np.isfinite(hi)) or hi <= lo:
            return
        self._levels = (float(lo), float(hi))
        self.setLevels(float(lo), float(hi))
        if self._gamma != 1.0:
            self.redraw()

    def axis_range(self, axis: int):
        """Current (min, max) of axis 0 = x, 1 = y."""
        return tuple(float(v) for v in self.view_box.viewRange()[axis])

    def set_axis_range(self, axis: int, lo: float, hi: float):
        """Fix one axis to an exact range. Counts as a deliberate range, so
        the automatic re-fit on resize leaves it alone."""
        if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
            return False
        self.view_box.disableAutoRange()
        # Marked deliberate *before* the range moves, not after: setXRange
        # emits while it runs, the plot box reshapes itself to the new range,
        # and the resize that follows would otherwise see a view nobody had
        # touched yet and fit it straight back to the data.
        self._user_ranged = True
        if axis == 0:
            self.view_box.setXRange(lo, hi, padding=0)
        else:
            self.view_box.setYRange(lo, hi, padding=0)
        return True

    def set_axis_inverted(self, axis: int, inverted: bool):
        if axis == 0:
            self.view_box.invertX(bool(inverted))
        else:
            self.view_box.invertY(bool(inverted))

    def axis_inverted(self, axis: int) -> bool:
        return bool(self.view_box.state["yInverted" if axis else "xInverted"])

    def set_auto_range(self, on: bool):
        """"Auto" here means this app's tight fit to the data (see
        :meth:`fit_to_data`), not pyqtgraph's padded autoRange: the whole
        point of these panels is that the axes sit against the data."""
        if on:
            self.fit_to_data()
        else:
            self._user_ranged = True

    # -- context menu ---------------------------------------------------
    def _add_context_actions(self):
        menu = self.view_box.menu
        menu.addSeparator()
        self.action_readout = _menu_action(menu, "Readout cursor")
        self.action_readout.setCheckable(True)
        self.action_readout.toggled.connect(self.set_readout_cursor_visible)

        self.action_selection = _menu_action(menu, "Selection box")
        self.action_selection.setCheckable(True)
        self.action_selection.toggled.connect(self.set_selection_visible)

        self.action_apply = _menu_action(menu, "Integrate selection into the other panel")
        self.action_apply.setEnabled(False)
        self.action_apply.triggered.connect(self.selectionApplied.emit)

        # Only advertised on panels where it means something (the spatial
        # overview); enable_popout_action() turns it on.
        self.action_popout = _menu_action(menu, "Open this position in a new window")
        self.action_popout.setVisible(False)
        self.action_popout.triggered.connect(self.popoutRequested.emit)

        # Likewise: exchanging the two axes only means something on a map
        # whose axes are the same kind of quantity.
        self.action_swap = _menu_action(menu, "Swap X and Y")
        self.action_swap.setCheckable(True)
        self.action_swap.setVisible(False)
        self.action_swap.toggled.connect(self.swapToggled.emit)

        # A cut along any direction, with a settable integration width. The
        # EDC and MDC curves are the two axis-aligned special cases of this;
        # a band that runs at an angle needs the general one, and reading it
        # off by eye from a picture is how a dispersion gets misquoted.
        self.action_profile = _menu_action(menu, "Line profile...")
        self.action_profile.setToolTip(
            "Drag a line across the image and read the intensity along it, "
            "integrated over a width you choose.")
        self.action_profile.triggered.connect(self.open_line_profile)

        # A stack of curves across whatever is inside the selection box. The
        # box is the natural way to say "these energies, this momentum
        # range" -- it is already on screen and already framing the feature,
        # so re-typing four numbers into the stack window is work the user
        # has effectively already done.
        self.action_stack = _menu_action(menu, "Stack plot of the selection...")
        self.action_stack.setToolTip(
            "EDCs or MDCs across the selection box, as a waterfall.")
        self.action_stack.triggered.connect(self.open_stack_plot)
        menu.aboutToShow.connect(self._sync_stack_action)

        # Export lives here and nowhere else. It belongs to the panel being
        # looked at rather than to the window -- a window may hold two -- and
        # right-clicking the picture you want is the least ambiguous way to
        # say which one.
        menu.addSeparator()
        # Made here, not by addMenu(): see _menu_action.
        export_menu = QMenu("Export", menu)
        menu.addMenu(export_menu)
        self._export_menu = export_menu
        self.action_export_nxs = _menu_action(export_menu, "Data (.nxs)...")
        self.action_export_nxs.setToolTip(
            "Write this panel's data, axes and metadata as a dataset this "
            "program opens again.")
        self.action_export_nxs.triggered.connect(self.export_dataset)
        self.action_export_image = _menu_action(export_menu, "Image...")
        self.action_export_image.setToolTip(
            "Render what is on screen at a chosen size and format, with or "
            "without the axes around it.")
        self.action_export_image.triggered.connect(self.export_image)
        self.action_export_eps = _menu_action(export_menu, "Axes only (.eps)...")
        self.action_export_eps.setToolTip(
            "The frame, ticks and labels as vector PostScript, to edit in "
            "Illustrator with the image placed inside.")
        self.action_export_eps.triggered.connect(self.export_axes)

    # -- stack plot -------------------------------------------------------
    def _sync_stack_action(self):
        """The stack entry says what it will act on before it is clicked."""
        corners = self.selection_corners()
        if corners is None:
            self.action_stack.setText("Stack plot (whole panel)...")
        else:
            x0, y0, x1, y1 = corners
            self.action_stack.setText(
                f"Stack plot of the selection "
                f"({x0:.3g}–{x1:.3g}, {y0:.3g}–{y1:.3g})...")

    def open_stack_plot(self):
        """A waterfall of curves across the selection box, or the whole panel.

        The curves are taken from the **full** data and the box only sets
        the range, rather than the data being cropped first: a curve that
        runs a little past the box is still one curve, and cropping it would
        change its area and so its normalisation.
        """
        from ui.process import StackWindow

        arrays = self.export_arrays()
        if arrays is None:
            return None
        values, x_axis, y_axis = arrays
        labels = (getattr(self, "x_label", "") or "x",
                  getattr(self, "y_label", "") or "y")
        window = StackWindow(values, (x_axis, y_axis), labels,
                             self.export_name(), self.window(),
                             getattr(self, "_colormap_name", "jet"))

        corners = self.selection_corners()
        if corners is not None:
            x0, y0, x1, y1 = corners
            # The longer side of the box says which way the curves should
            # run: a box that is wide and short is asking for a few curves
            # across many energies, not the other way round.
            span_x = (x1 - x0) / (abs(x_axis[-1] - x_axis[0]) or 1.0)
            span_y = (y1 - y0) / (abs(y_axis[-1] - y_axis[0]) or 1.0)
            edc = span_y >= span_x
            window.direction.setCurrentIndex(0 if edc else 1)
            lo, hi = ((x0, x1) if edc else (y0, y1))
            window.from_edit.setText(f"{lo:.6g}")
            window.to_edit.setText(f"{hi:.6g}")
            window.limit_from.setText(f"{(y0 if edc else x0):.6g}")
            window.limit_to.setText(f"{(y1 if edc else x1):.6g}")
            window.redraw()

        window.show()
        self._stack_window = window
        return window

    # -- line profile -----------------------------------------------------
    def open_line_profile(self):
        """Open (or raise) this panel's line-profile window."""
        existing = getattr(self, "_profile_window", None)
        if existing is not None and existing.isVisible():
            existing.raise_()
            return existing
        arrays = self.export_arrays()
        if arrays is None:
            return None
        window = LineProfileWindow(self, self.window())
        window.show()
        self._profile_window = window
        return window

    # -- export ----------------------------------------------------------
    def export_arrays(self):
        """``(array_xy, x_axis, y_axis)`` for this panel, raw -- never the
        smoothed or gamma-shaped version that is drawn. Subclasses that hold
        an image provide it."""
        return None

    def export_name(self) -> str:
        """What the export is called by default: the file, this panel, and
        -- for a window whose slider moves through a cube -- where the
        slider is *now*.

        The position is asked of the window through ``self.window()`` rather
        than kept here, so nothing about exporting makes a panel hold its
        window alive (see :meth:`__init__`).
        """
        parts = [str(self.export_prefix or "export")]
        if self.export_title:
            parts.append(f"[{self.export_title}]")
        owner = self.window()
        position = owner.slice_label() if hasattr(owner, "slice_label") else ""
        if position:
            parts.append(str(position))
        return " ".join(parts)

    def export_lut(self):
        try:
            return colormaps.get_lut(self._colormap_name, flip=self._colormap_flip)
        except Exception:
            return colormaps.get_lut(DEFAULT_COLORMAP)

    def export_ranges(self):
        """The (x, y) ranges an export covers: what the view is showing.

        There used to be a choice here -- the data's own extent, or the view
        -- because an aspect lock made the view wider than the data and
        neither answer was right. The plot box is reshaped instead of the
        range now (see :class:`AspectBox`), so the view *is* the data until
        somebody zooms, and then the zoom is what they want exported.
        """
        return self.axis_range(0), self.axis_range(1)

    def _open_export_dialog(self, page: int):
        arrays = self.export_arrays()
        if arrays is None or arrays[0] is None:
            from PyQt5.QtWidgets import QMessageBox
            QMessageBox.warning(self, "Export", "Nothing displayed to export yet.")
            return None
        existing = self._export_dialog
        if existing is not None and existing.isVisible():
            existing.tabs.setCurrentIndex(page)
            existing.raise_()
            existing.activateWindow()
            return existing
        dialog = ExportDialog(self, page)
        self._export_dialog = dialog
        weak = weakref.ref(self)
        def forget(*_):
            view = weak()
            if view is not None:
                view._export_dialog = None
        dialog.finished.connect(forget)
        dialog.show()
        return dialog

    def export_dataset(self):
        return self._open_export_dialog(0)

    def export_image(self):
        return self._open_export_dialog(1)

    def export_axes(self):
        return self._open_export_dialog(2)

    def enable_popout_action(self, enabled: bool, text: str = None):
        self.action_popout.setVisible(enabled)
        if text:
            self.action_popout.setText(text)

    # -- picking points off the image ------------------------------------
    def start_point_picking(self, callback):
        """Left-clicking the image now reports (x, y) to ``callback``.

        This is how both tools that need points off a picture get them --
        the corners of an arbitrary cut and the points along a feature to
        straighten. The MATLAB originals call ``getpts``/``getline``, which
        take over the figure until the user presses Enter; here the window
        stays live throughout, so the picked points can be looked at,
        dragged and retyped while the dialog that collects them is open.
        """
        self.stop_point_picking()
        self._pick_callback = callback
        self._pick_connection = self.scene.sigMouseClicked.connect(self._on_pick_click)

    def stop_point_picking(self):
        if getattr(self, "_pick_callback", None) is None:
            return
        try:
            self.scene.sigMouseClicked.disconnect(self._on_pick_click)
        except (TypeError, RuntimeError):
            pass
        self._pick_callback = None

    def _on_pick_click(self, event):
        callback = getattr(self, "_pick_callback", None)
        if callback is None or event.button() != Qt.LeftButton:
            return
        point = self.view_box.mapSceneToView(event.scenePos())
        callback(float(point.x()), float(point.y()))

    def show_picked_points(self, xs, ys):
        """Draw (or clear) the markers for the points picked so far."""
        if not len(xs):
            if self._pick_markers is not None:
                self.view_box.removeItem(self._pick_markers)
                self._pick_markers = None
            return
        if self._pick_markers is None:
            # Bright, not black: these markers sit on top of the data, which
            # is dark wherever it matters -- a black cross on a band is
            # invisible exactly where it is being placed.
            self._pick_markers = pg.ScatterPlotItem(
                size=15, symbol="+", pen=pg.mkPen("#f5a623", width=3),
                brush=None)
            self._pick_markers.setZValue(30)
            self.view_box.addItem(self._pick_markers, ignoreBounds=True)
        self._pick_markers.setData(np.asarray(xs, dtype=float),
                                   np.asarray(ys, dtype=float))

    def show_overlay_path(self, xs, ys, dashed: bool = False):
        """Draw (or clear) a line over the image: the path an arbitrary cut
        will follow, or the curve fitted through picked points."""
        if xs is None or not len(xs):
            if self._overlay_curve is not None:
                self.view_box.removeItem(self._overlay_curve)
                self._overlay_curve = None
            return
        if self._overlay_curve is None:
            pen = pg.mkPen("#bd4921", width=2,
                           style=Qt.DashLine if dashed else Qt.SolidLine)
            self._overlay_curve = pg.PlotDataItem(pen=pen)
            self._overlay_curve.setZValue(29)
            self.view_box.addItem(self._overlay_curve, ignoreBounds=True)
        self._overlay_curve.setData(np.asarray(xs, dtype=float),
                                    np.asarray(ys, dtype=float))

    def clear_overlays(self):
        self.show_picked_points([], [])
        self.show_overlay_path(None, None)

    # -- Brillouin-zone overlay --------------------------------------------
    # Deliberately separate from show_overlay_path/clear_overlays above:
    # those are transient picking aids that every tool wipes on its own
    # close (see _PickerDialog._release), while a Brillouin-zone drawing is
    # something the user turns on and expects to stay on the contour while
    # they go use some other tool (an arbitrary cut, a k-conversion redo,
    # ...) alongside it.
    #
    # Kept as named layers (rather than one curve) so the moire mode can
    # show the top layer's, the bottom layer's and the moire cell's own
    # zone at once, each in its own colour, and update or drop any one of
    # them without touching the others.
    def set_bz_layer(self, name: str, polygons, color: str = "#39c2ff",
                     width: float = 1.4):
        """Draw (``polygons`` non-empty) or remove (empty/None) one named
        Brillouin-zone layer. Each polygon is ``(xs, ys)``; a layer with
        several (a tiled zone) is drawn as one
        :class:`pyqtgraph.PlotDataItem` with ``NaN`` separating each
        polygon -- pyqtgraph treats that as a break in the line, so this is
        one GPU-friendly item per layer rather than one per polygon, which
        matters once tiling hands this a few dozen of them.
        """
        if not polygons:
            item = self._bz_layers.pop(name, None)
            if item is not None:
                self.view_box.removeItem(item)
            return
        xs, ys = [], []
        for px, py in polygons:
            if xs:
                xs.append(np.nan)
                ys.append(np.nan)
            xs.extend(np.asarray(px, dtype=float))
            ys.extend(np.asarray(py, dtype=float))
        item = self._bz_layers.get(name)
        if item is None:
            item = pg.PlotDataItem(pen=pg.mkPen(color, width=width), connect="finite")
            item.setZValue(28)
            self.view_box.addItem(item, ignoreBounds=True)
            self._bz_layers[name] = item
        else:
            item.setPen(pg.mkPen(color, width=width))
        item.setData(np.asarray(xs), np.asarray(ys))

    def clear_bz_overlay(self):
        """Remove every Brillouin-zone layer set by :meth:`set_bz_layer`."""
        for name in list(self._bz_layers):
            self.set_bz_layer(name, [])

    def enable_swap_action(self, enabled: bool, text: str = None):
        self.action_swap.setVisible(enabled)
        if text:
            self.action_swap.setText(text)

    def enable_selection_action(self, enabled: bool, text: str = None):
        """Owners call this to advertise (or hide) the 'apply selection'
        menu entry, which only makes sense on some pages."""
        self.action_apply.setVisible(enabled)
        if text:
            self.action_apply.setText(text)

    # -- axes ------------------------------------------------------------
    def _store_axes(self, xarray, yarray):
        self._xarray = np.asarray(xarray)
        self._yarray = np.asarray(yarray)

    @staticmethod
    def _step(axis):
        return (axis[1] - axis[0]) if len(axis) > 1 else 1.0

    def data_extent(self):
        """The image's outer bounds, i.e. the axis range extended by half a
        pixel at each end (pixels are centred on their axis values)."""
        if self._xarray is None or self._yarray is None:
            return None
        dx, dy = self._step(self._xarray), self._step(self._yarray)
        return (float(self._xarray[0]) - dx / 2, float(self._xarray[-1]) + dx / 2,
                float(self._yarray[0]) - dy / 2, float(self._yarray[-1]) + dy / 2)

    def fit_to_data(self):
        """Set the view to exactly the image's bounds, with no padding.

        Done explicitly rather than through ``autoRange()``, which pads. On
        a panel shaped by an :class:`AspectBox` there is nothing more to it:
        no aspect is locked in the ViewBox, so both axes land exactly on the
        data at every window size.

        The other case is a panel with EDC/MDC curves pinned to the image's
        edges, which cannot be reshaped and so keeps pyqtgraph's own aspect
        lock. There one axis has to be widened past the data, and pyqtgraph
        does it around the *previous* view centre -- leaving the data off to
        one side and inflating both axes past what the lock needs. So the
        expansion is solved for here and applied symmetrically instead.
        """
        if self._fitting:
            return
        extent = self.data_extent()
        if extent is None:
            return
        x0, x1, y0, y1 = extent
        cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
        dx, dy = x1 - x0, y1 - y0
        vb = self.view_box
        vb.disableAutoRange()

        aspect = vb.state["aspectLocked"]
        if aspect:
            # pyqtgraph's lock enforces exactly
            #     (x_span / width_px) == aspect * (y_span / height_px)
            # so the expansion can be solved for directly. Letting setRange
            # discover it instead takes several passes and settles slightly
            # wide on *both* axes; solving it keeps one axis exactly on the
            # data and grows only the other, symmetrically about the centre.
            w_px = max(float(vb.width()), 1.0)
            h_px = max(float(vb.height()), 1.0)
            x_needed = aspect * dy * w_px / h_px
            if x_needed >= dx:
                dx = x_needed          # y stays tight
            else:
                dy = dx * h_px / (aspect * w_px)   # x stays tight

        # Guarded, because setRange emits while it runs: the plot box
        # reshapes, the child resizes, and the resize asks for a fit again.
        self._fitting = True
        try:
            vb.setRange(xRange=(cx - dx / 2, cx + dx / 2),
                        yRange=(cy - dy / 2, cy + dy / 2), padding=0)
        finally:
            self._fitting = False
        self._user_ranged = False

    def _image_origin(self):
        """Bottom-left corner for setImage so that each pixel is *centred* on
        its axis value.

        pyqtgraph places the first pixel's corner at ``pos``, so passing
        ``axis[0]`` shifts the whole image half a pixel against the tick
        labels and makes the image extend a full step past the last data
        point -- invisible while the view had padding, obvious once the axes
        sit tight against the data.
        """
        return (self._xarray[0] - self._step(self._xarray) / 2.0,
                self._yarray[0] - self._step(self._yarray) / 2.0)

    # -- readout cursor --------------------------------------------------
    def set_readout_cursor_visible(self, visible: bool):
        if visible and self.readout_cursor is None:
            self._build_readout_cursor()
        elif not visible and self.readout_cursor is not None:
            for item in (self.readout_cursor, self.readout_vline, self.readout_hline,
                         self.readout_vband, self.readout_hband):
                self.removeItem(item)
            self.plot.removeItem(self.readout_label)
            self.readout_cursor = None
            self.readout_label = None
            self.readout_vline = self.readout_hline = None
            self.readout_vband = self.readout_hband = None
            self._band_half = {}
        if self.action_readout.isChecked() != visible:
            self.action_readout.setChecked(visible)
        self.readoutToggled.emit(visible)

    def _build_readout_cursor(self):
        if self._xarray is None:
            return
        cx = float(self._xarray[len(self._xarray) // 2])
        cy = float(self._yarray[len(self._yarray) // 2])
        self.readout_cursor = pg.TargetItem(size=14, pos=(cx, cy),
                                             pen=pg.mkPen("#222222", width=2))

        # Crosshair lines, like the spatial map's, so it is obvious *which*
        # column and row the EDC and MDC are being taken from. Each is drawn
        # in the colour of the axis it holds constant (AXIS_COLORS), and the
        # curve it feeds -- the vertical line the EDC, the horizontal the
        # MDC -- in the same colour.
        self.readout_vline = pg.InfiniteLine(angle=90, movable=False)
        self.readout_hline = pg.InfiniteLine(angle=0, movable=False)
        # Shaded bands showing the "+/-" integration windows, so a widened
        # EDC/MDC shows on the image exactly how much it is summing.
        self.readout_vband = pg.LinearRegionItem(orientation="vertical", movable=False)
        self.readout_hband = pg.LinearRegionItem(orientation="horizontal", movable=False)
        for band in (self.readout_vband, self.readout_hband):
            # Over the image, translucent: at the old z of -10 the bands sat
            # *behind* the image and were never seen.
            band.setZValue(5)
            band.hide()
        self.readout_vline.setZValue(10)
        self.readout_hline.setZValue(10)
        self._apply_readout_colors()

        # ignoreBounds: these are overlays, not data -- an InfiniteLine or a
        # floating label must never take part in working out the view range.
        for item in (self.readout_vband, self.readout_hband,
                     self.readout_vline, self.readout_hline, self.readout_cursor):
            self.plot.addItem(item, ignoreBounds=True)
        self.readout_label = pg.TextItem(color="#111111", anchor=(0, 1),
                                         fill=pg.mkBrush(255, 255, 255, 190))
        self.readout_label.setZValue(30)
        self.readout_cursor.setZValue(20)
        self.plot.addItem(self.readout_label, ignoreBounds=True)
        self.readout_cursor.sigPositionChanged.connect(self._update_readout)
        self._update_readout()

    def set_readout_colors(self, x_color: str, y_color: str):
        """Colours for the vertical line (constant x) and the horizontal one
        (constant y), with their integration bands."""
        self._readout_colors = (x_color, y_color)
        self._apply_readout_colors()

    def _apply_readout_colors(self):
        if self.readout_cursor is None:
            return
        x_color, y_color = self._readout_colors
        self.readout_vline.setPen(pg.mkPen(x_color, width=2))
        self.readout_hline.setPen(pg.mkPen(y_color, width=2))
        # Translucent fill, and dashed edges in the same colour: a blue band
        # on the blue end of a colormap is invisible as a fill alone.
        for band, color in ((self.readout_vband, x_color), (self.readout_hband, y_color)):
            band.setBrush(pg.mkBrush(*_rgba(color, BAND_ALPHA)))
            for line in band.lines:
                line.setPen(pg.mkPen(color, width=1, style=Qt.DashLine))
            band.update()

    def readout_position(self):
        """The cursor's (x, y) in data units, or None when it is off."""
        if self.readout_cursor is None:
            return None
        x, y = self.readout_cursor.pos()
        return float(x), float(y)

    def readout_indices(self):
        """The data point under the cursor as (ix, iy), or None."""
        position = self.readout_position()
        if position is None or self._xarray is None:
            return None
        return (int(np.argmin(np.abs(self._xarray - position[0]))),
                int(np.argmin(np.abs(self._yarray - position[1]))))

    def axis_bounds(self, axis: int):
        """(lowest, highest) data value on axis 0 = x or 1 = y."""
        values = self._xarray if axis == 0 else self._yarray
        return float(np.nanmin(values)), float(np.nanmax(values))

    def move_readout_to(self, x=None, y=None):
        """Put the cursor at ``(x, y)`` in data units; None keeps that
        coordinate. Returns None on success, or -- if either value is outside
        the data -- a sentence saying so, and the cursor does not move at
        all (not even along the axis whose value was fine)."""
        if self.readout_cursor is None or self._xarray is None:
            return "The readout cursor is not on."
        current = self.readout_position()
        target = [current[0] if x is None else float(x),
                  current[1] if y is None else float(y)]
        for axis, value, given in ((0, target[0], x), (1, target[1], y)):
            if given is None:
                continue
            lo, hi = self.axis_bounds(axis)
            tol = 1e-9 * max(1.0, abs(hi - lo))
            if not np.isfinite(value) or value < lo - tol or value > hi + tol:
                name = (self.x_label if axis == 0 else self.y_label) or ("x", "y")[axis]
                return (f"{name} = {value:.6g} is outside the data "
                        f"({lo:.6g} … {hi:.6g}); the cursor stays where it is.")
        self.readout_cursor.setPos(target[0], target[1])
        return None

    def set_readout_index(self, ix=None, iy=None):
        """Put the cursor on a data point; None keeps that coordinate. Does
        nothing if it is already there, so a cursor being dragged between
        two points is not snapped back onto one."""
        indices = self.readout_indices()
        if indices is None:
            return
        ix = indices[0] if ix is None else int(np.clip(ix, 0, len(self._xarray) - 1))
        iy = indices[1] if iy is None else int(np.clip(iy, 0, len(self._yarray) - 1))
        if (ix, iy) == indices:
            return
        x, y = self.readout_position()
        if ix != indices[0]:
            x = float(self._xarray[ix])
        if iy != indices[1]:
            y = float(self._yarray[iy])
        self.readout_cursor.setPos(x, y)

    def refresh_readout_text(self):
        """Re-read the value under the cursor after new data has landed,
        without reporting a move."""
        if self.readout_cursor is None or self._xarray is None:
            return
        indices = self.readout_indices()
        self.readout_label.setText(self._readout_text(*indices))

    def _readout_text(self, ix, iy):
        text = f"x={self._xarray[ix]:.4g}\ny={self._yarray[iy]:.4g}"
        if self._value_lookup is not None:
            try:
                value = self._value_lookup(ix, iy)
                if value is not None:
                    text += f"\nvalue={value:.6g}"
            except Exception:
                pass
        return text

    def readout_value(self):
        """The data value under the cursor, or None."""
        indices = self.readout_indices()
        if indices is None or self._value_lookup is None:
            return None
        try:
            value = self._value_lookup(*indices)
            return None if value is None else float(value)
        except Exception:
            return None

    def set_readout_bands(self, x_half_width: float, y_half_width: float):
        """Show/hide the shaded integration windows around the crosshair.
        Called by the panel whenever an EDC/MDC "+/-" box changes."""
        self._set_bands(self.readout_vband, self.readout_hband,
                        x_half_width, y_half_width)

    def band_half_width(self, band) -> float:
        """The half-width a band is drawn with (0 when hidden)."""
        return float(self._band_half.get(band, 0.0)) if band is not None else 0.0

    def _set_bands(self, vband, hband, x_half_width, y_half_width):
        if self.readout_cursor is None or self._xarray is None:
            return
        x, y = self.readout_cursor.pos()
        for band, centre, half in ((vband, x, x_half_width), (hband, y, y_half_width)):
            half = float(half or 0.0)
            self._band_half[band] = half
            if half > 0:
                band.setRegion((centre - half, centre + half))
                band.show()
            else:
                band.hide()

    def _update_readout(self):
        if self.readout_cursor is None or self._xarray is None:
            return
        x, y = self.readout_cursor.pos()
        ix = int(np.argmin(np.abs(self._xarray - x)))
        iy = int(np.argmin(np.abs(self._yarray - y)))
        self.readout_label.setText(self._readout_text(ix, iy))
        self.readout_label.setPos(x, y)
        self.readout_vline.setPos(x)
        self.readout_hline.setPos(y)
        self._reposition_bands(x, y)
        self.readoutMoved.emit(ix, iy)

    def _reposition_bands(self, x, y):
        """Keep the integration bands centred on the crosshair as it moves,
        preserving each band's current width."""
        for band, centre in ((self.readout_vband, x), (self.readout_hband, y)):
            half = self._band_half.get(band, 0.0)
            if half > 0:
                band.setRegion((centre - half, centre + half))

    # -- rectangle selection ---------------------------------------------
    def set_selection_visible(self, visible: bool):
        if visible and self.selection_roi is None:
            self._build_selection_roi()
        elif not visible and self.selection_roi is not None:
            self.removeItem(self.selection_roi)
            self.selection_roi = None
            self.action_apply.setEnabled(False)
        if self.action_selection.isChecked() != visible:
            self.action_selection.setChecked(visible)
        self.selectionChanged.emit()

    def _build_selection_roi(self):
        if self._xarray is None:
            return
        # default box: the middle half of each axis, so it is immediately
        # visible and grabbable rather than a zero-size box in a corner
        x0, x1 = float(self._xarray[0]), float(self._xarray[-1])
        y0, y1 = float(self._yarray[0]), float(self._yarray[-1])
        w, h = (x1 - x0), (y1 - y0)
        self.selection_roi = pg.RectROI([x0 + 0.25 * w, y0 + 0.25 * h],
                                         [0.5 * w, 0.5 * h],
                                         pen=ROI_PEN, hoverPen=ROI_HOVER_PEN,
                                         handlePen=ROI_HANDLE_PEN,
                                         handleHoverPen=ROI_HANDLE_HOVER_PEN)
        # Bigger grab targets than pyqtgraph's 5px default -- the stock
        # handles are fiddly to hit, especially on a high-DPI screen.
        self.selection_roi.handleSize = ROI_HANDLE_SIZE
        # Corner handles on the two draggable corners whose coordinates the
        # panel below reports.
        self.selection_roi.addScaleHandle([1, 1], [0, 0])
        self.selection_roi.addScaleHandle([0, 0], [1, 1])
        self.plot.addItem(self.selection_roi, ignoreBounds=True)
        self.action_apply.setEnabled(True)
        self.selection_roi.sigRegionChanged.connect(lambda *_: self.selectionChanged.emit())

    def selection_corners(self):
        """The two draggable corners as ``(x0, y0, x1, y1)`` in data
        coordinates (lower-left then upper-right), or None."""
        if self.selection_roi is None:
            return None
        pos = self.selection_roi.pos()
        size = self.selection_roi.size()
        x0, x1 = sorted([float(pos[0]), float(pos[0] + size[0])])
        y0, y1 = sorted([float(pos[1]), float(pos[1] + size[1])])
        return x0, y0, x1, y1

    def set_selection_corners(self, x0, y0, x1, y1):
        """Move/resize the selection to the given data coordinates, creating
        it if necessary. Used by the coordinate boxes below the panel so a
        region can be typed in exactly instead of dragged by eye."""
        if self.selection_roi is None:
            self.set_selection_visible(True)
        if self.selection_roi is None:   # still nothing to place it on
            return
        x0, x1 = sorted([float(x0), float(x1)])
        y0, y1 = sorted([float(y0), float(y1)])
        # A zero-width box would select nothing; widen it to one pixel.
        if self._xarray is not None and x1 <= x0:
            x1 = x0 + abs(self._step(self._xarray))
        if self._yarray is not None and y1 <= y0:
            y1 = y0 + abs(self._step(self._yarray))
        self.selection_roi.setPos((x0, y0), finish=False)
        self.selection_roi.setSize((x1 - x0, y1 - y0), finish=True)
        self.selectionChanged.emit()

    def selection_index_ranges(self):
        """Current selection as ``((ix0, ix1), (iy0, iy1))`` inclusive index
        bounds into the stored axes, or None if there's no selection.
        Bounds are clipped to the data and always span at least one pixel,
        so a box dragged partly outside the image still yields a usable
        region rather than an empty sum."""
        if self.selection_roi is None or self._xarray is None:
            return None
        pos = self.selection_roi.pos()
        size = self.selection_roi.size()
        xlo, xhi = sorted([pos[0], pos[0] + size[0]])
        ylo, yhi = sorted([pos[1], pos[1] + size[1]])

        def bounds(axis, lo, hi):
            i0 = int(np.clip(np.searchsorted(axis, lo, side="left"), 0, len(axis) - 1))
            i1 = int(np.clip(np.searchsorted(axis, hi, side="right") - 1, 0, len(axis) - 1))
            if i1 < i0:
                i0, i1 = min(i0, i1), max(i0, i1)
            return i0, i1

        ascending_x = len(self._xarray) < 2 or self._xarray[1] >= self._xarray[0]
        ascending_y = len(self._yarray) < 2 or self._yarray[1] >= self._yarray[0]
        xr = bounds(self._xarray if ascending_x else self._xarray[::-1], xlo, xhi)
        yr = bounds(self._yarray if ascending_y else self._yarray[::-1], ylo, yhi)
        if not ascending_x:
            xr = (len(self._xarray) - 1 - xr[1], len(self._xarray) - 1 - xr[0])
        if not ascending_y:
            yr = (len(self._yarray) - 1 - yr[1], len(self._yarray) - 1 - yr[0])
        return xr, yr

    def selection_summary(self):
        """Human-readable description of the current selection, for the
        status bar."""
        ranges = self.selection_index_ranges()
        if ranges is None:
            return "no selection"
        (ix0, ix1), (iy0, iy1) = ranges
        return (f"x[{self._xarray[ix0]:.4g}..{self._xarray[ix1]:.4g}] "
                f"y[{self._yarray[iy0]:.4g}..{self._yarray[iy1]:.4g}] "
                f"({ix1 - ix0 + 1}x{iy1 - iy0 + 1} px)")


class SpatialImageView(_InteractiveImageBase):
    """Spatial overview (x, y) map for spem_4d/spem_1d. Drag the red target
    to move the readout pixel; the separate teal readout cursor and the
    selection box are opt-in from the right-click menu."""

    pixelChanged = pyqtSignal(int, int)

    def __init__(self, parent=None):
        super().__init__(parent, name="SpatialImageView")
        # Both axes are the same physical length unit here, so 1:1 is the
        # honest default -- a stretched spatial map misrepresents the shape
        # of whatever is on the sample.
        self.view.setAspectLocked(True, ratio=1.0)

        # Just the marker: full-width crosshair lines were visual noise here,
        # and the position is already reported in the status bar. (The EDC/MDC
        # crosshair on the spectrum panels is a different thing -- there the
        # lines show which row and column the curves come from.)
        self.cursor = pg.TargetItem(size=16, pos=(0, 0), pen=RED_PEN)
        self.plot.addItem(self.cursor, ignoreBounds=True)
        self.cursor.sigPositionChanged.connect(self._on_cursor_moved)
        self._image2d = None

    def set_data(self, array2d: np.ndarray, xarray: np.ndarray, yarray: np.ndarray,
                 reset_cursor: bool = True):
        """``array2d`` is (row, col) i.e. (len(yarray), len(xarray))."""
        self._store_axes(xarray, yarray)
        self._image2d = array2d
        self._value_lookup = lambda ix, iy: self._image2d[iy, ix]
        self._render(auto_range=True)
        if reset_cursor:
            cx = float(self._xarray[len(self._xarray) // 2])
            cy = float(self._yarray[len(self._yarray) // 2])
            self.cursor.setPos((cx, cy))

    def _render(self, auto_range=False):
        # _image2d is (row, col) = (y, x), so the display factors come back
        # in that order too.
        image = self._for_display(self._image2d)
        fy, fx = self._display_factors
        self.setImage(image,
                      pos=self._image_origin(),
                      scale=(self._step(self._xarray) / fx,
                             self._step(self._yarray) / fy),
                      autoRange=False, autoLevels=self._levels is None)
        self._reapply_levels()
        if auto_range:
            self.fit_to_data()

    def redraw(self):
        if self._image2d is not None and self._xarray is not None:
            self._render()

    def export_arrays(self):
        if self._image2d is None or self._xarray is None:
            return None
        # stored (row, col) = (y, x); exports are (x, y) throughout
        return np.asarray(self._image2d).T, self._xarray, self._yarray

    def _on_cursor_moved(self):
        if self._xarray is None:
            return
        x, y = self.cursor.pos()
        x = min(max(x, self._xarray[0]), self._xarray[-1])
        y = min(max(y, self._yarray[0]), self._yarray[-1])
        col = int(np.argmin(np.abs(self._xarray - x)))
        row = int(np.argmin(np.abs(self._yarray - y)))
        self.pixelChanged.emit(row, col)


class FrameImageView(_InteractiveImageBase):
    """2D frame display (E vs k, constant-E map, ...). Aspect starts
    unlocked, because the two axes usually carry unrelated physical units
    (eV vs degrees) where an equal on-screen scale would just squash the
    image; the panel's ratio controls can lock it when that does make
    sense (e.g. a constant-E map, whose axes are both degrees)."""

    def __init__(self, parent=None):
        super().__init__(parent, name="FrameImageView")
        self.view.setAspectLocked(False)
        self._last_frame = None

    def set_frame(self, frame: np.ndarray, x_axis: np.ndarray, y_axis: np.ndarray):
        """``frame`` has shape (len(x_axis), len(y_axis)) and is displayed
        with x_axis horizontal, y_axis vertical (transposed for row-major
        rendering). ``last_frame`` keeps the *raw* array, so smoothing stays
        a display effect and exports write unsmoothed data."""
        self._last_frame = frame
        self._store_axes(x_axis, y_axis)
        self._value_lookup = lambda ix, iy: self._last_frame[ix, iy]
        self._render(auto_range=True)
        # The cursor stays put; what is under it has changed.
        self.refresh_readout_text()

    def _render(self, auto_range=False):
        # last_frame is (x, y); it is transposed for row-major drawing, so
        # the first display factor is the x one.
        image = self._for_display(self._last_frame)
        fx, fy = self._display_factors
        self.setImage(image.T,
                      pos=self._image_origin(),
                      scale=(self._step(self._xarray) / fx,
                             self._step(self._yarray) / fy),
                      autoRange=False, autoLevels=self._levels is None)
        self._reapply_levels()
        if auto_range:
            self.fit_to_data()

    def redraw(self):
        if self._last_frame is not None and self._xarray is not None:
            self._render()

    @property
    def last_frame(self):
        return self._last_frame

    def export_arrays(self):
        if self._last_frame is None or self._xarray is None:
            return None
        return np.asarray(self._last_frame), self._xarray, self._yarray


# --------------------------------------------------------------------------
# Export
# --------------------------------------------------------------------------
def _qimage_to_rgba(image: QImage) -> np.ndarray:
    """A QImage as an (h, w, 4) uint8 array, so the writing is done by PIL
    for every format rather than depending on which Qt image plugins the
    machine happens to have."""
    image = image.convertToFormat(QImage.Format_RGBA8888)
    width, height = image.width(), image.height()
    buffer = image.constBits()
    # byteCount() was renamed sizeInBytes() in Qt 5.10; both are around in
    # the wild, so ask for whichever this build has.
    size = (image.sizeInBytes() if hasattr(image, "sizeInBytes")
            else image.byteCount())
    buffer.setsize(int(size))
    array = np.frombuffer(buffer, dtype=np.uint8).reshape(height,
                                                          image.bytesPerLine() // 4, 4)
    return np.array(array[:, :width, :], copy=True)


def render_framed(rgba: np.ndarray, *, x_range, y_range, x_label: str = "",
                  y_label: str = "", font_px: float = 0.0, line_width: float = 0.0,
                  invert_x: bool = False, invert_y: bool = False,
                  background=(255, 255, 255)) -> np.ndarray:
    """Put the rendered data inside a drawn axis frame, at the same pixel
    scale.

    The axes are drawn here rather than screen-grabbed from the panel: a
    grab is only ever as sharp as the window, whereas this draws the frame,
    the ticks and the text at the export's own size. It uses the same tick
    machinery as the EPS export (:func:`export.tick_values`), so a
    framed PNG and the vector axes beside it carry identical ticks.

    ``rgba`` becomes the inside of the box, unchanged; the returned array is
    larger by the margins the labels need.
    """
    rgba = np.asarray(rgba, dtype=np.uint8)
    box_h, box_w = rgba.shape[0], rgba.shape[1]
    if font_px <= 0:
        # Legible at any export size: bigger pictures get proportionally
        # bigger type rather than 9-pt text lost in a 3000-px image.
        font_px = max(11.0, round(min(box_w, box_h) / 34.0))
    if line_width <= 0:
        line_width = max(1.0, font_px / 8.0)

    font = QFont("Helvetica")
    font.setPixelSize(int(round(font_px)))
    metrics = QFontMetricsF(font)

    def text_width(text: str) -> float:
        # QFontMetrics.width() is deprecated from Qt 5.11 in favour of
        # horizontalAdvance(); both exist somewhere in the versions this has
        # to run on, so prefer the newer name where it is there.
        if hasattr(metrics, "horizontalAdvance"):
            return metrics.horizontalAdvance(text)
        return metrics.width(text)

    x_ticks, x_step = export.tick_values(*x_range, target=5)
    y_ticks, y_step = export.tick_values(*y_range, target=5)
    x_texts = [export.tick_format(x_step) % v for v in x_ticks]
    y_texts = [export.tick_format(y_step) % v for v in y_ticks]

    gap = 0.4 * font_px
    text_h = metrics.height()
    left = int(round(gap * 2 + max([text_width(t) for t in y_texts] or [0])
                     + (text_h * 1.2 if y_label else 0)))
    bottom = int(round(gap * 2 + text_h + (text_h * 1.1 if x_label else 0)))
    top = int(round(text_h * 0.6))
    right = int(round(max(text_h * 0.6,
                          (text_width(x_texts[-1]) / 2 if x_texts else 0))))

    image = QImage(left + box_w + right, top + box_h + bottom,
                   QImage.Format_RGBA8888)
    image.fill(QColor(*background))
    painter = QPainter(image)
    painter.setRenderHint(QPainter.Antialiasing, True)
    painter.setRenderHint(QPainter.TextAntialiasing, True)
    painter.setFont(font)

    # The buffer must outlive the QImage wrapping it -- a temporary here
    # would be collected while Qt is still reading from it.
    buffer = np.ascontiguousarray(rgba)
    data = QImage(buffer.data, box_w, box_h, box_w * 4, QImage.Format_RGBA8888)
    painter.drawImage(left, top, data)

    pen = QPen(QColor(0, 0, 0))
    pen.setWidthF(line_width)
    painter.setPen(pen)
    painter.drawRect(left, top, box_w, box_h)

    x0, x1 = float(min(x_range)), float(max(x_range))
    y0, y1 = float(min(y_range)), float(max(y_range))
    tick = max(3.0, font_px * 0.45)

    def to_px(value):
        frac = (value - x0) / (x1 - x0) if x1 > x0 else 0.0
        return left + (1.0 - frac if invert_x else frac) * box_w

    def to_py(value):
        frac = (value - y0) / (y1 - y0) if y1 > y0 else 0.0
        # screen y grows downward, so an un-inverted axis counts from the bottom
        return top + (frac if invert_y else 1.0 - frac) * box_h

    for value, text in zip(x_ticks, x_texts):
        x = to_px(value)
        painter.drawLine(int(x), top + box_h, int(x), int(top + box_h - tick))
        painter.drawLine(int(x), top, int(x), int(top + tick))
        painter.drawText(int(round(x - text_width(text) / 2)),
                         int(round(top + box_h + gap + metrics.ascent())), text)
    for value, text in zip(y_ticks, y_texts):
        y = to_py(value)
        painter.drawLine(left, int(y), int(left + tick), int(y))
        painter.drawLine(left + box_w, int(y), int(left + box_w - tick), int(y))
        painter.drawText(int(round(left - gap - text_width(text))),
                         int(round(y + metrics.ascent() / 2 - 1)), text)

    if x_label:
        painter.drawText(
            int(round(left + box_w / 2 - text_width(x_label) / 2)),
            int(round(top + box_h + gap + metrics.ascent() + text_h)), x_label)
    if y_label:
        painter.save()
        painter.translate(gap + metrics.ascent(), top + box_h / 2)
        painter.rotate(-90)
        painter.drawText(int(round(-text_width(y_label) / 2)), 0, y_label)
        painter.restore()
    painter.end()
    return _qimage_to_rgba(image)


class ExportDialog(QDialog):
    """The one place a panel is written out from.

    Three tabs, because a figure is assembled from three different files and
    they must agree with one another:

    * **Data (.nxs)** -- the numbers, reopenable here.
    * **Image** -- what is drawn, at a chosen pixel size, with or without
      the axes around it.
    * **Axes (.eps)** -- the frame and its annotation as vector PostScript,
      nothing inside, to edit in Illustrator with the image placed in the
      box.

    All three cover the same region -- whatever the panel is showing, which
    is the data itself unless it has been zoomed -- so an axes-only EPS
    always matches the image meant to go inside it. That is the whole reason
    these live in one dialog instead of three menu entries.
    """

    TAB_DATA, TAB_IMAGE, TAB_AXES = 0, 1, 2

    def __init__(self, view, page: int = 0):
        super().__init__(view.window())
        # Weakly, and through a property: the panel keeps a reference to this
        # dialog while it is open, so a strong one back would be a reference
        # cycle -- and a cycle of Qt objects is destroyed whenever the
        # garbage collector happens to run, which is how pyqtgraph gets torn
        # down in the middle of one of its own calls.
        self._view_ref = weakref.ref(view)
        self.setWindowTitle("Export")
        self.setModal(False)
        self.resize(430, 470)

        layout = QVBoxLayout(self)

        header = QLabel("All three exports cover what the panel is showing:")
        header.setWordWrap(True)
        layout.addWidget(header)

        self.range_label = QLabel("")
        self.range_label.setWordWrap(True)
        self.range_label.setStyleSheet("color: #667ea1;")
        layout.addWidget(self.range_label)

        self.tabs = QTabWidget()
        self.tabs.addTab(self._build_data_tab(), "Data (.nxs)")
        self.tabs.addTab(self._build_image_tab(), "Image")
        self.tabs.addTab(self._build_axes_tab(), "Axes (.eps)")
        self.tabs.setCurrentIndex(int(page))
        self.tabs.currentChanged.connect(lambda *_: self._refresh())
        layout.addWidget(self.tabs, 1)

        self.status = QLabel("")
        self.status.setWordWrap(True)
        self.status.setTextInteractionFlags(Qt.TextSelectableByMouse)
        layout.addWidget(self.status)

        buttons = QHBoxLayout()
        self.save_button = QPushButton("Save...")
        self.save_button.setDefault(True)
        self.save_button.clicked.connect(self.save)
        close = QPushButton("Close")
        close.clicked.connect(self.reject)
        buttons.addStretch(1)
        buttons.addWidget(self.save_button)
        buttons.addWidget(close)
        layout.addLayout(buttons)

        # Both size pairs start at the shape on screen rather than at the
        # spin boxes' defaults.
        self._width_changed()
        self._eps_width_changed()

    @property
    def view(self):
        return self._view_ref()

    # -- the tabs ---------------------------------------------------------
    def _build_data_tab(self):
        page = QWidget()
        form = QFormLayout(page)
        self.name_edit = QLineEdit(self.view.export_name())
        self.name_edit.setToolTip("The dataset's name inside the file, and "
                                  "what the list shows when it is reopened.")
        form.addRow("Dataset name", self.name_edit)
        note = QLabel(
            "Writes the panel's array, both axes with their labels and this "
            "file's metadata into a .nxs this program reads back as a dataset "
            "of its own.")
        note.setWordWrap(True)
        form.addRow(note)
        self.data_shape_label = QLabel("")
        form.addRow("Size", self.data_shape_label)
        return page

    def _build_image_tab(self):
        page = QWidget()
        form = QFormLayout(page)
        self.format_combo = QComboBox()
        for label, _, _, _ in export.IMAGE_FORMATS:
            self.format_combo.addItem(label)
        form.addRow("Format", self.format_combo)

        self.width_box = QSpinBox()
        self.width_box.setRange(16, 20000)
        self.width_box.setValue(1600)
        self.width_box.setSingleStep(100)
        self.height_box = QSpinBox()
        self.height_box.setRange(16, 20000)
        self.height_box.setValue(1200)
        self.keep_aspect = QCheckBox("Keep the picture's proportions")
        self.keep_aspect.setToolTip(
            "The height follows the width at the shape the panel is drawn at "
            "-- the data's own shape when the ratio is locked. Untick it only "
            "to stretch the export deliberately.")
        self.keep_aspect.setChecked(True)
        self.width_box.valueChanged.connect(self._width_changed)
        self.keep_aspect.toggled.connect(self._width_changed)
        form.addRow("Width (px)", self.width_box)
        form.addRow("Height (px)", self.height_box)
        form.addRow("", self.keep_aspect)

        self.dpi_box = QSpinBox()
        self.dpi_box.setRange(36, 2400)
        self.dpi_box.setValue(300)
        self.dpi_box.setToolTip(
            "Written into the file as metadata: the pixel count is what is "
            "set above, this is how big a layout program places it.")
        form.addRow("Resolution (dpi)", self.dpi_box)

        self.with_axes = QCheckBox("Include the axes, ticks and labels")
        self.with_axes.setChecked(True)
        self.with_axes.setToolTip(
            "Off gives the coloured data alone, the exact size set above and "
            "with everything outside the data transparent -- which is what "
            "goes inside the box of the .eps axes.")
        self.with_axes.toggled.connect(lambda *_: self._refresh())
        form.addRow("", self.with_axes)

        self.smooth_export = QCheckBox("Interpolate between data points")
        self.smooth_export.setChecked(bool(getattr(self.view, "_smooth_on", False)))
        self.smooth_export.setToolTip(
            "Follows the panel's Smooth box by default. Off draws each "
            "measured point as a hard rectangle, which is the honest choice "
            "for a coarse scan.")
        form.addRow("", self.smooth_export)

        self.image_size_label = QLabel("")
        self.image_size_label.setWordWrap(True)
        form.addRow("File", self.image_size_label)
        return page

    def _build_axes_tab(self):
        page = QWidget()
        form = QFormLayout(page)
        note = QLabel(
            "Frame, ticks and labels only. The image export (axes off) drops "
            "into the empty box this reports.")
        note.setWordWrap(True)
        form.addRow(note)

        self.eps_width = QDoubleSpinBox()
        self.eps_width.setRange(5.0, 500.0)
        self.eps_width.setValue(80.0)
        self.eps_width.setSuffix(" mm")
        self.eps_height = QDoubleSpinBox()
        self.eps_height.setRange(5.0, 500.0)
        self.eps_height.setValue(60.0)
        self.eps_height.setSuffix(" mm")
        self.eps_keep_aspect = QCheckBox("Keep the picture's proportions")
        self.eps_keep_aspect.setChecked(True)
        self.eps_width.valueChanged.connect(self._eps_width_changed)
        self.eps_keep_aspect.toggled.connect(self._eps_width_changed)
        form.addRow("Box width", self.eps_width)
        form.addRow("Box height", self.eps_height)
        form.addRow("", self.eps_keep_aspect)

        self.eps_font = QDoubleSpinBox()
        self.eps_font.setRange(3.0, 40.0)
        self.eps_font.setValue(9.0)
        self.eps_font.setSuffix(" pt")
        form.addRow("Font size", self.eps_font)
        self.eps_line = QDoubleSpinBox()
        self.eps_line.setRange(0.1, 5.0)
        self.eps_line.setSingleStep(0.1)
        self.eps_line.setValue(0.8)
        self.eps_line.setSuffix(" pt")
        form.addRow("Line width", self.eps_line)
        self.eps_minor = QSpinBox()
        self.eps_minor.setRange(0, 9)
        self.eps_minor.setValue(1)
        self.eps_minor.setToolTip("Minor ticks between each pair of labelled ones.")
        form.addRow("Minor ticks", self.eps_minor)
        self.eps_inward = QCheckBox("Ticks inside the box")
        self.eps_inward.setChecked(True)
        form.addRow("", self.eps_inward)
        return page

    # -- state ------------------------------------------------------------
    def ranges(self):
        return self.view.export_ranges()

    def _aspect(self) -> float:
        """Height / width of the export.

        Asked of the panel (:meth:`_InteractiveImageBase.export_aspect`): a
        locked ratio gives the data's own shape, an unlocked one the shape
        the panel is drawn at. Never the window's shape -- an export must not
        change proportions because a window was dragged wider.
        """
        return self.view.export_aspect()

    def _width_changed(self, *_):
        if self.keep_aspect.isChecked():
            self.height_box.blockSignals(True)
            self.height_box.setValue(max(16, int(round(self.width_box.value()
                                                        * self._aspect()))))
            self.height_box.blockSignals(False)
        self.height_box.setEnabled(not self.keep_aspect.isChecked())
        self._refresh()

    def _eps_width_changed(self, *_):
        if self.eps_keep_aspect.isChecked():
            self.eps_height.blockSignals(True)
            self.eps_height.setValue(max(5.0, self.eps_width.value() * self._aspect()))
            self.eps_height.blockSignals(False)
        self.eps_height.setEnabled(not self.eps_keep_aspect.isChecked())
        self._refresh()

    def _refresh(self):
        (x0, x1), (y0, y1) = self.ranges()
        xl = self.view.x_label or "x"
        yl = self.view.y_label or "y"
        self.range_label.setText(f"{xl}: {x0:.6g} to {x1:.6g}\n{yl}: {y0:.6g} to {y1:.6g}")
        cropped = self._cropped()
        if cropped is not None:
            array = cropped[0]
            self.data_shape_label.setText(f"{array.shape[0]} x {array.shape[1]} points")
        self.image_size_label.setText(
            f"{self.width_box.value()} x {self.height_box.value()} px of data"
            + (", plus the margins the labels need" if self.with_axes.isChecked()
               else ", transparent outside the data"))

    def _cropped(self):
        """The panel's arrays restricted to what is on screen. Unzoomed that
        is all of it; zoomed in, the .nxs holds the part being looked at,
        like the other two exports."""
        arrays = self.view.export_arrays()
        if arrays is None or arrays[0] is None:
            return None
        array, xaxis, yaxis = arrays
        array = np.asarray(array, dtype=float)
        xaxis, yaxis = np.asarray(xaxis, dtype=float), np.asarray(yaxis, dtype=float)
        (x0, x1), (y0, y1) = self.ranges()
        keep_x = (xaxis >= min(x0, x1)) & (xaxis <= max(x0, x1))
        keep_y = (yaxis >= min(y0, y1)) & (yaxis <= max(y0, y1))
        if not keep_x.any() or not keep_y.any():
            return array, xaxis, yaxis
        return array[np.ix_(keep_x, keep_y)], xaxis[keep_x], yaxis[keep_y]

    # -- writing ----------------------------------------------------------
    def _ask_path(self, title, suffix, file_filter):
        stem = "".join(ch for ch in self.view.export_name()
                       if ch.isalnum() or ch in "._- ") or "export"
        folder = os.path.dirname(self.view.export_source_path or "") or os.getcwd()
        path, _ = QFileDialog.getSaveFileName(self, title,
                                              os.path.join(folder, stem + suffix),
                                              file_filter)
        if not path:
            return None
        if not path.lower().endswith(suffix.lower()):
            path += suffix
        return path

    def save(self):
        index = self.tabs.currentIndex()
        try:
            if index == self.TAB_DATA:
                return self.save_dataset()
            if index == self.TAB_IMAGE:
                return self.save_image()
            return self.save_axes()
        except (OSError, PermissionError, ValueError) as exc:
            QMessageBox.warning(self, "Export", f"Could not write the file:\n{exc}")
            return None

    def save_dataset(self, path: str = None):
        cropped = self._cropped()
        if cropped is None:
            return None
        array, xaxis, yaxis = cropped
        path = path or self._ask_path("Export data", ".nxs", "NeXus files (*.nxs)")
        if path is None:
            return None
        name = self.name_edit.text().strip() or self.view.export_name()
        info = dict(self.view.export_info or {})
        info.update({"export.source": os.path.basename(self.view.export_source_path or ""),
                     "export.panel": self.view.export_title or "",
                     "export.region": "%.6g..%.6g, %.6g..%.6g" % (
                         *self.ranges()[0], *self.ranges()[1])})
        save_dataset(path, [{
            "name": name, "kind": "cut", "axes": (xaxis, yaxis),
            "labels": {"x": self.view.x_label or "x", "y": self.view.y_label or "y"},
            "value": array, "info": info, "motors": {}}])
        self.status.setText(f"Wrote {os.path.basename(path)} "
                            f"({array.shape[0]} x {array.shape[1]} points).")
        return path

    def rendered(self) -> np.ndarray:
        """The image as it will be written, framed or not."""
        arrays = self.view.export_arrays()
        if arrays is None or arrays[0] is None:
            raise ValueError("nothing displayed to export")
        array, xaxis, yaxis = arrays
        x_range, y_range = self.ranges()
        rgba = export.render_rgba(
            array, xaxis, yaxis, x_range=x_range, y_range=y_range,
            lut=self.view.export_lut(), levels=getattr(self.view, "_levels", None),
            gamma=getattr(self.view, "_gamma", 1.0),
            width=self.width_box.value(), height=self.height_box.value(),
            interpolate=self.smooth_export.isChecked(),
            invert_x=self.view.axis_inverted(0), invert_y=self.view.axis_inverted(1))
        if not self.with_axes.isChecked():
            return rgba
        return render_framed(rgba, x_range=x_range, y_range=y_range,
                             x_label=self.view.x_label or "",
                             y_label=self.view.y_label or "",
                             invert_x=self.view.axis_inverted(0),
                             invert_y=self.view.axis_inverted(1))

    def save_image(self, path: str = None):
        label, suffix, pil_format, alpha = export.IMAGE_FORMATS[
            self.format_combo.currentIndex()]
        path = path or self._ask_path("Export image", suffix,
                                      f"{pil_format} (*{suffix})")
        if path is None:
            return None
        rgba = self.rendered()
        export.write_image(path, rgba, pil_format=pil_format,
                               keep_alpha=alpha, dpi=self.dpi_box.value())
        self.status.setText(f"Wrote {os.path.basename(path)} "
                            f"({rgba.shape[1]} x {rgba.shape[0]} px).")
        return path

    def save_axes(self, path: str = None):
        path = path or self._ask_path("Export axes", ".eps", "EPS (*.eps)")
        if path is None:
            return None
        x_range, y_range = self.ranges()
        box = export.axes_eps(
            path, x_range=x_range, y_range=y_range,
            x_label=self.view.x_label or "", y_label=self.view.y_label or "",
            width_mm=self.eps_width.value(), height_mm=self.eps_height.value(),
            font_size=self.eps_font.value(), line_width=self.eps_line.value(),
            minor_ticks=self.eps_minor.value(),
            ticks_inward=self.eps_inward.isChecked(),
            invert_x=self.view.axis_inverted(0), invert_y=self.view.axis_inverted(1))
        self.status.setText(
            f"Wrote {os.path.basename(path)}. Place the image (axes off) at "
            f"x {box['x']:.1f} pt, y {box['y']:.1f} pt, "
            f"{box['width']:.1f} x {box['height']:.1f} pt -- the empty box.")
        return box


def _pin_axis_geometry(image_view, curve_readout, left_width: int = 62,
                       bottom_height: int = 30) -> None:
    """Give the image and its two curve panels identical axis extents.

    Linking the views (``setXLink``/``setYLink``) makes the *data* ranges
    match, but the panels only line up on screen if their plot areas start
    at the same pixel -- which means the left axes (image and MDC) and the
    bottom axes (image and EDC) must be exactly as wide/tall as each other.
    pyqtgraph sizes axes to their tick labels, so without pinning them the
    columns drift apart whenever the numbers change length.
    """
    image_view.view.getAxis("left").setWidth(left_width)
    image_view.view.getAxis("bottom").setHeight(bottom_height)
    curve_readout.mdc_plot.getPlotItem().getAxis("left").setWidth(left_width)
    curve_readout.edc_plot.getPlotItem().getAxis("bottom").setHeight(bottom_height)


class LineProfileWindow(QWidget):
    """Intensity along an arbitrary line across a panel, integrated over a
    width.

    The EDC and MDC readouts are the two axis-aligned special cases of this.
    A band that runs at an angle -- which is most of them away from a high
    symmetry direction -- has no axis-aligned cut that follows it, and the
    usual workaround is to read positions off the picture by eye. This gives
    the actual numbers, and can hand them to the figure composer or add them
    to the list as a dataset.

    The line lives on the panel it was opened from and is removed when this
    window closes, so a panel is never left with a stray ROI on it.
    """

    def __init__(self, view, parent=None):
        super().__init__(parent, Qt.Window)
        self.view = view
        self.setWindowTitle("Line profile")
        self.resize(640, 420)

        arrays = view.export_arrays()
        values, x_axis, y_axis = arrays
        self.x_axis = np.asarray(x_axis, dtype=float)
        self.y_axis = np.asarray(y_axis, dtype=float)

        layout = QVBoxLayout(self)
        self.plot = pg.PlotWidget()
        self.plot.setLabel("bottom", "distance along the line")
        self.plot.setLabel("left", "Intensity")
        strip_stock_menu(self.plot.getPlotItem())
        layout.addWidget(self.plot, 1)

        row = QHBoxLayout()
        row.addWidget(QLabel("Width"))
        self.width_box = QDoubleSpinBox()
        self.width_box.setDecimals(5)
        self.width_box.setRange(0.0, 1e6)
        span_y = float(abs(self.y_axis[-1] - self.y_axis[0])) or 1.0
        self.width_box.setSingleStep(span_y / 100)
        self.width_box.setToolTip(
            "Integrate this far either side of the line, in the y axis's "
            "units. Zero samples the line itself.")
        self.width_box.valueChanged.connect(self.refresh)
        row.addWidget(self.width_box)
        self.samples_box = QSpinBox()
        self.samples_box.setRange(16, 4000)
        self.samples_box.setValue(400)
        self.samples_box.valueChanged.connect(self.refresh)
        row.addWidget(QLabel("Samples"))
        row.addWidget(self.samples_box)
        self.readout = QLabel("")
        row.addWidget(self.readout, 1)
        to_list = QPushButton("To list")
        to_list.setToolTip("Add the profile to the main list as a dataset.")
        to_list.clicked.connect(self.to_list)
        row.addWidget(to_list)
        layout.addLayout(row)

        # Start across the middle of the image, which is visible whatever
        # the data is and needs no guessing about where the interesting
        # feature is.
        x0, x1 = float(self.x_axis[0]), float(self.x_axis[-1])
        y0, y1 = float(self.y_axis[0]), float(self.y_axis[-1])
        self.roi = pg.LineSegmentROI(
            [[x0 + 0.2 * (x1 - x0), y0 + 0.5 * (y1 - y0)],
             [x0 + 0.8 * (x1 - x0), y0 + 0.5 * (y1 - y0)]],
            pen=pg.mkPen("#f5a623", width=3),
            hoverPen=pg.mkPen("#ffcc66", width=4))
        for handle in self.roi.getHandles():
            handle.pen = pg.mkPen("#f5a623", width=3)
        self.view.view_box.addItem(self.roi)
        self.roi.sigRegionChanged.connect(self.refresh)
        self.refresh()

    def points(self):
        """The line's two ends in data coordinates."""
        handles = self.roi.getSceneHandlePositions()
        ends = [self.view.view_box.mapSceneToView(position)
                for _, position in handles]
        return ((float(ends[0].x()), float(ends[0].y())),
                (float(ends[1].x()), float(ends[1].y())))

    def profile(self):
        """``(distance, intensity)`` along the line."""
        arrays = self.view.export_arrays()
        if arrays is None:
            return None, None
        values, x_axis, y_axis = arrays
        values = np.asarray(values, dtype=float)
        x_axis = np.asarray(x_axis, dtype=float)
        y_axis = np.asarray(y_axis, dtype=float)
        (px0, py0), (px1, py1) = self.points()
        n = int(self.samples_box.value())
        xs = np.linspace(px0, px1, n)
        ys = np.linspace(py0, py1, n)

        length = float(np.hypot(px1 - px0, py1 - py0))
        distance = np.linspace(0.0, length, n)

        width = float(self.width_box.value())
        if width > 0 and length > 0:
            # Offsets are taken perpendicular to the line in *data* units,
            # which is the only definition that stays put when the window is
            # resized -- a perpendicular in pixels is not a perpendicular in
            # the data.
            nx, ny = -(py1 - py0) / length, (px1 - px0) / length
            steps = max(3, int(round(width / _axis_step_of(y_axis))) | 1)
            offsets = np.linspace(-width, width, steps)
        else:
            nx, ny, offsets = 0.0, 0.0, np.array([0.0])

        total = np.zeros(n)
        count = np.zeros(n)
        for offset in offsets:
            sample = _sample_grid(values, x_axis, y_axis,
                                  xs + nx * offset, ys + ny * offset)
            good = np.isfinite(sample)
            total = np.where(good, total + np.nan_to_num(sample), total)
            count = count + good
        with np.errstate(invalid="ignore", divide="ignore"):
            intensity = np.where(count > 0, total / np.maximum(count, 1e-12),
                                 np.nan)
        return distance, intensity

    def refresh(self):
        distance, intensity = self.profile()
        if distance is None:
            return
        self.plot.clear()
        self.plot.plot(distance, intensity, pen=pg.mkPen("#bd4921", width=2))
        (px0, py0), (px1, py1) = self.points()
        self.readout.setText(
            f"({px0:.4g}, {py0:.4g}) \u2192 ({px1:.4g}, {py1:.4g})   "
            f"length {distance[-1]:.4g}")

    def to_list(self):
        distance, intensity = self.profile()
        if distance is None:
            return None
        owner = self.view.window()
        emitter = getattr(owner, "datasetCreated", None)
        if emitter is None:
            QMessageBox.information(
                self, "Line profile",
                "This window cannot add datasets to the list. Export the "
                "panel instead, or open the profile from a viewer.")
            return None
        (px0, py0), (px1, py1) = self.points()
        data = MemoryData(
            "cut", (distance, np.array([0.0])),
            np.asarray(intensity, dtype=float).reshape(-1, 1),
            {"x": "distance along the line", "y": "Intensity"},
            source_label=f"{self.view.export_prefix or 'panel'}_profile",
            parameters={"from": [px0, py0], "to": [px1, py1],
                        "width": self.width_box.value()},
            prefix="proc.profile")
        emitter.emit(data)
        return data

    def closeEvent(self, event):
        try:
            self.view.view_box.removeItem(self.roi)
        except Exception:
            pass
        self.view._profile_window = None
        super().closeEvent(event)


def _axis_step_of(axis) -> float:
    axis = np.asarray(axis, dtype=float)
    if axis.size < 2:
        return 1.0
    step = abs(float(axis[-1] - axis[0])) / (axis.size - 1)
    return step if step > 0 else 1.0


def _sample_grid(values, x_axis, y_axis, xq, yq):
    """Bilinear samples of a regular grid, NaN-aware, outside -> NaN."""
    from scipy import ndimage

    x0, y0 = float(x_axis[0]), float(y_axis[0])
    dx = _axis_step_of(x_axis) * (1 if x_axis[-1] >= x_axis[0] else -1)
    dy = _axis_step_of(y_axis) * (1 if y_axis[-1] >= y_axis[0] else -1)
    ix = (np.asarray(xq, float) - x0) / dx
    iy = (np.asarray(yq, float) - y0) / dy
    # A handle dragged onto the edge of the image lands a rounding error
    # past it, and a strict bounds test then drops that one sample -- a NaN
    # at the end of every profile that happens to span the full width. The
    # tolerance is a millionth of a pixel, far too small to invent data.
    TOL = 1e-6
    inside = ((ix >= -TOL) & (ix <= values.shape[0] - 1 + TOL) &
              (iy >= -TOL) & (iy <= values.shape[1] - 1 + TOL))
    ix = np.clip(ix, 0.0, values.shape[0] - 1)
    iy = np.clip(iy, 0.0, values.shape[1] - 1)

    mask = np.isfinite(values)
    filled = np.where(mask, values, 0.0)
    coords = np.array([ix, iy])
    num = ndimage.map_coordinates(filled, coords, order=1, mode="constant", cval=0.0)
    den = ndimage.map_coordinates(mask.astype(float), coords, order=1,
                                  mode="constant", cval=0.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where((den > 0.5) & inside, num / np.maximum(den, 1e-12), np.nan)


class LevelBar(QWidget):
    """A slim two-handle min/max bar for the image's colour scaling, sitting
    above the image.

    This replaces pyqtgraph's HistogramLUTWidget, which drew the whole
    intensity distribution and a vertical gradient strip down the side of
    every panel -- a lot of screen space for what is normally a two-number
    decision. Here the distribution is drawn *inside* the bar itself, behind
    the two handles, so it costs no extra height.

    The distribution is worth having because choosing levels blind is
    guesswork: ARPES intensities are heavily skewed, a handful of hot pixels
    can hold the top of the scale ten times above everything else, and the
    only way to see that is to look at where the counts actually are. It is
    drawn on a square-root vertical scale, since on a linear one the
    background peak is the only thing visible and the interesting tail is a
    flat line. The **Clip** button snaps the handles to a percentile pair,
    which is the same judgement made numerically.
    """

    levelsChanged = pyqtSignal(float, float)
    gammaChanged = pyqtSignal(float)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._lo = 0.0
        self._hi = 1.0
        self._syncing = False
        self._sample = None         # the values behind the drawn histogram

        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(4)

        layout.addWidget(QLabel("Levels"))
        self.min_label = QLabel("-")
        self.min_label.setMinimumWidth(70)
        self.min_label.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        layout.addWidget(self.min_label)

        self.bar = pg.PlotWidget()
        self.bar.setFixedHeight(28)
        self.bar.setMenuEnabled(False)
        self.bar.hideAxis("left")
        self.bar.hideAxis("bottom")
        self.bar.setMouseEnabled(x=False, y=False)
        self.bar.getPlotItem().hideButtons()
        self.bar.setYRange(0, 1, padding=0)
        layout.addWidget(self.bar, stretch=1)

        #: the intensity distribution, drawn behind the handles
        self.histogram = pg.PlotCurveItem(
            pen=pg.mkPen("#9aa5ab", width=1), brush=pg.mkBrush(154, 165, 171, 110),
            fillLevel=0.0)
        self.histogram.setZValue(-10)
        self.bar.addItem(self.histogram)

        self.max_label = QLabel("-")
        self.max_label.setMinimumWidth(70)
        layout.addWidget(self.max_label)

        self.clip_button = QPushButton("Clip")
        self.clip_button.setMaximumWidth(46)
        self.clip_button.setToolTip(
            "Put the handles at the 1st and 99th percentile of the counts, "
            "which throws away the hot pixels that otherwise hold the whole "
            "top of the scale. Right-click for other pairs.")
        self.clip_button.clicked.connect(lambda: self.clip_percentiles(1.0, 99.0))
        self.clip_button.setContextMenuPolicy(Qt.CustomContextMenu)
        self.clip_button.customContextMenuRequested.connect(self._clip_menu)
        layout.addWidget(self.clip_button)

        # The same two numbers as a percentage of the data's own range, which
        # is how the MATLAB tool asks for them: "show me the bottom 40%" is a
        # setting that carries from one scan to the next, where an absolute
        # count does not.
        layout.addWidget(separator())
        self.min_pct = QDoubleSpinBox()
        self.max_pct = QDoubleSpinBox()
        for box, label, value, tip in (
                (self.min_pct, "Min%", 0.0, "Lower end of the colour scale, as a "
                 "percentage of this image's range."),
                (self.max_pct, "Max%", 100.0, "Upper end of the colour scale, as a "
                 "percentage of this image's range.")):
            layout.addWidget(QLabel(label))
            box.setRange(0.0, 100.0)
            box.setDecimals(2)
            box.setValue(value)
            box.setKeyboardTracking(False)
            box.setButtonSymbols(QDoubleSpinBox.NoButtons)
            box.setMaximumWidth(66)
            box.setToolTip(tip)
            box.editingFinished.connect(self._on_percent_typed)
            layout.addWidget(box)

        layout.addWidget(separator())
        layout.addWidget(QLabel("Gamma"))
        self.gamma_box = QDoubleSpinBox()
        self.gamma_box.setRange(0.05, 5.0)
        self.gamma_box.setSingleStep(0.1)
        self.gamma_box.setDecimals(2)
        self.gamma_box.setValue(1.0)
        self.gamma_box.setKeyboardTracking(False)
        self.gamma_box.setMaximumWidth(64)
        self.gamma_box.setToolTip(
            "Intensity exponent inside the colour window: below 1 lifts weak "
            "features out of the background, above 1 suppresses them. "
            "Display only -- exported data stays raw.")
        self.gamma_box.valueChanged.connect(self.gammaChanged.emit)
        layout.addWidget(self.gamma_box)

        self.region = pg.LinearRegionItem(values=(0.0, 1.0), orientation="vertical",
                                           brush=pg.mkBrush(13, 115, 119, 60),
                                           hoverBrush=pg.mkBrush(13, 115, 119, 100))
        for line in self.region.lines:
            line.setPen(pg.mkPen("#0d7377", width=3))
            line.setHoverPen(pg.mkPen("#bd4921", width=4))
        self.bar.addItem(self.region)
        self.region.sigRegionChanged.connect(self._on_region_changed)

    def set_data_range(self, lo: float, hi: float, keep_levels: bool = False):
        """Point the bar at a new image's value range. Unless ``keep_levels``,
        the handles reset to the full range (which is what autoLevels just
        did to the image)."""
        if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
            lo, hi = 0.0, 1.0
        self._lo, self._hi = float(lo), float(hi)
        self.bar.setXRange(self._lo, self._hi, padding=0)
        self.region.setBounds((self._lo, self._hi))
        if not keep_levels:
            self.region.blockSignals(True)
            self.region.setRegion((self._lo, self._hi))
            self.region.blockSignals(False)
        self._update_labels()

    def set_histogram(self, values, bins: int = 256):
        """Draw the intensity distribution of the current image in the bar.

        The vertical scale is the square root of the counts per bin,
        normalised: on a linear scale an ARPES image is one spike at the
        background level and nothing else visible, and the part that decides
        the levels is the tail.
        """
        values = np.asarray(values, dtype=float).ravel()
        values = values[np.isfinite(values)]
        if values.size < 2 or self._hi <= self._lo:
            self.histogram.setData(x=np.array([]), y=np.array([]))
            self._sample = None
            return
        # Keep the values themselves (subsampled if there are a lot of them)
        # rather than reading percentiles back off the drawn bins. An ARPES
        # image spans four decades, so 256 bins put the entire background
        # peak -- most of the pixels -- inside bin zero, and every percentile
        # below about 98% then comes back as the same number.
        self._sample = (values if values.size <= 400_000
                        else values[::int(np.ceil(values.size / 400_000))])
        counts, edges = np.histogram(values, bins=int(bins),
                                     range=(self._lo, self._hi))
        height = np.sqrt(counts.astype(float))
        top = float(height.max()) or 1.0
        centres = (edges[:-1] + edges[1:]) / 2.0
        # 0.9 rather than 1.0 so the tallest bin does not touch the handles
        # and make them hard to grab.
        self.histogram.setData(x=centres, y=0.9 * height / top)

    def clip_percentiles(self, low: float, high: float):
        """Move the handles to a percentile pair of the image's values."""
        sample = getattr(self, "_sample", None)
        if sample is None or sample.size < 2:
            return
        lo, hi = (float(v) for v in np.percentile(sample, [low, high]))
        if hi <= lo:
            lo, hi = self._lo, self._hi
        self.region.setRegion((lo, hi))

    def _clip_menu(self, position):
        menu = QMenu(self.clip_button)
        pairs = [(0.5, 99.5), (1.0, 99.0), (2.0, 98.0), (5.0, 95.0), (0.0, 100.0)]
        actions = {menu.addAction(f"{lo:g}% – {hi:g}%"): (lo, hi)
                   for lo, hi in pairs}
        chosen = menu.exec_(self.clip_button.mapToGlobal(position))
        if chosen in actions:
            self.clip_percentiles(*actions[chosen])

    def levels(self):
        lo, hi = self.region.getRegion()
        return float(lo), float(hi)

    def _on_region_changed(self):
        self._update_labels()
        lo, hi = self.levels()
        if hi > lo:
            self.levelsChanged.emit(lo, hi)

    # -- percentages ------------------------------------------------------
    def _to_percent(self, value: float) -> float:
        span = self._hi - self._lo
        if span <= 0:
            return 0.0
        return float(np.clip(100.0 * (value - self._lo) / span, 0.0, 100.0))

    def _from_percent(self, pct: float) -> float:
        return self._lo + (self._hi - self._lo) * float(pct) / 100.0

    def _on_percent_typed(self):
        if self._syncing:
            return
        lo = self._from_percent(self.min_pct.value())
        hi = self._from_percent(self.max_pct.value())
        if hi <= lo:                    # an inverted window has no meaning
            self._update_labels()
            return
        self.region.setRegion((lo, hi))

    def gamma(self) -> float:
        return float(self.gamma_box.value())

    def _update_labels(self):
        lo, hi = self.levels()
        self.min_label.setText(f"{lo:.4g}")
        self.max_label.setText(f"{hi:.4g}")
        self._syncing = True
        try:
            self.min_pct.setValue(self._to_percent(lo))
            self.max_pct.setValue(self._to_percent(hi))
        finally:
            self._syncing = False


class CurveReadout(QObject):
    """The EDC/MDC pair that goes with the readout cursor on an E-vs-k panel.

    - **EDC** (energy distribution curve): counts vs the vertical axis,
      summed over a window of the horizontal axis centred on the cursor.
      Drawn to the **right** of the image with counts horizontal and energy
      vertical, its energy axis linked to the image's, so features line up
      across the two panels.
    - **MDC** (momentum distribution curve): counts vs the horizontal axis,
      summed over a window of the vertical axis centred on the cursor. Drawn
      **below** the image, angle axis linked to the image's.

    This is a plain QObject, not a widget: the two plots and the width boxes
    are separate widgets that the owning :class:`ImagePanel` places into its
    grid, which is what makes the shared-axis alignment possible.

    Each curve has its own "+/-" half-width in physical units (degrees, eV
    -- whatever the axis carries); 0 means a single row/column, which is the
    plain single-coordinate readout. Integrating over a window is the usual
    way to get a usable curve out of a noisy single-pixel line.

    The names are literal only when the vertical axis really is energy; on
    other axis pairs the titles fall back to generic profiles.

    **EDC → list / MDC → list** put the curve on screen in the launcher's
    list as a one-dimensional dataset, to open in the curve viewer, fit, or
    save. :attr:`curveToList` carries it out as a plain dict; the window
    that owns the panel turns it into a dataset, since it knows the source.
    """

    #: ``{"which": "edc"|"mdc", "x", "y", "x_label", "position",
    #: "position_label", "half_width", "n_summed"}``
    curveToList = pyqtSignal(dict)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._frame = None          # (x, y) ordered, raw
        self._xarray = None
        self._yarray = None
        self._ix = self._iy = 0
        self._x_label = "x"
        self._y_label = "y"

        # -- EDC: counts (horizontal) vs energy (vertical), right of image --
        self.edc_plot = pg.PlotWidget()
        self.edc_plot.setMinimumWidth(90)
        self.edc_plot.showGrid(x=True, y=True, alpha=0.25)
        self.edc_plot.hideAxis("left")          # energy is already on the image
        strip_stock_menu(self.edc_plot.getPlotItem())
        self.edc_curve = self.edc_plot.plot(pen=pg.mkPen(EDC_COLOR, width=2))

        # -- MDC: angle (horizontal) vs counts (vertical), below image --
        self.mdc_plot = pg.PlotWidget()
        self.mdc_plot.setMinimumHeight(90)
        self.mdc_plot.showGrid(x=True, y=True, alpha=0.25)
        self.mdc_plot.hideAxis("bottom")        # angle is already on the image
        strip_stock_menu(self.mdc_plot.getPlotItem())
        self.mdc_curve = self.mdc_plot.plot(pen=pg.mkPen(MDC_COLOR, width=2))

        # -- the two +/- width boxes, placed in the panel's options row --
        self.controls = QWidget()
        row = QHBoxLayout(self.controls)
        row.setContentsMargins(0, 0, 0, 0)
        row.addWidget(QLabel("EDC +/-"))
        self.edc_width = QDoubleSpinBox()
        self.edc_width.setDecimals(3)
        self.edc_width.setRange(0.0, 1e6)
        self.edc_width.setSingleStep(0.1)
        self.edc_width.setMaximumWidth(80)
        self.edc_width.setToolTip(
            "Half-width on the horizontal axis to sum the EDC over "
            "(0 = a single column under the cursor).")
        row.addWidget(self.edc_width)
        row.addWidget(QLabel("MDC +/-"))
        self.mdc_width = QDoubleSpinBox()
        self.mdc_width.setDecimals(3)
        self.mdc_width.setRange(0.0, 1e6)
        self.mdc_width.setSingleStep(0.1)
        self.mdc_width.setMaximumWidth(80)
        self.mdc_width.setToolTip(
            "Half-width on the vertical axis to sum the MDC over "
            "(0 = a single row under the cursor).")
        row.addWidget(self.mdc_width)
        # The positions are reported here rather than as plot titles: a title
        # adds height to the top of the EDC only, which would push its plot
        # area 30 px below the image's and break the shared-axis alignment.
        self.info_label = QLabel("-")
        row.addWidget(self.info_label)
        self.edc_button = QPushButton("EDC → list")
        self.edc_button.setToolTip(
            "Put this EDC (summed over its ± window) in the main list as a "
            "curve of its own: open, fit, or save it from there.")
        self.edc_button.clicked.connect(lambda: self.export_curve("edc"))
        self.mdc_button = QPushButton("MDC → list")
        self.mdc_button.setToolTip(
            "Put this MDC (summed over its ± window) in the main list as a "
            "curve of its own.")
        self.mdc_button.clicked.connect(lambda: self.export_curve("mdc"))
        row.addWidget(self.edc_button)
        row.addWidget(self.mdc_button)
        row.addStretch(1)

        self.edc_width.valueChanged.connect(self.refresh)
        self.mdc_width.valueChanged.connect(self.refresh)

    def set_colors(self, edc_color: str, mdc_color: str):
        """The EDC in the colour of the vertical cursor line that feeds it,
        the MDC in the horizontal line's."""
        self.edc_curve.setPen(pg.mkPen(edc_color, width=2))
        self.mdc_curve.setPen(pg.mkPen(mdc_color, width=2))

    def widgets(self):
        """Everything that shows and hides with the readout cursor."""
        return (self.edc_plot, self.mdc_plot, self.controls)

    def link_to(self, image_view):
        """Share axes with the image: the EDC's vertical axis follows the
        image's vertical axis and the MDC's horizontal axis follows the
        image's horizontal one, so panning or zooming the image keeps all
        three aligned."""
        plot_item = image_view.view
        self.edc_plot.getPlotItem().setYLink(plot_item)
        self.mdc_plot.getPlotItem().setXLink(plot_item)

    def set_source(self, frame, xarray, yarray, x_label="x", y_label="y"):
        """Point the curves at a frame. ``frame`` is (x, y) ordered, the same
        convention ``FrameImageView.set_frame`` uses."""
        self._frame = None if frame is None else np.asarray(frame)
        self._xarray = None if xarray is None else np.asarray(xarray)
        self._yarray = None if yarray is None else np.asarray(yarray)
        self._x_label, self._y_label = x_label or "x", y_label or "y"
        self.edc_plot.setLabel("bottom", "counts")
        self.mdc_plot.setLabel("left", "counts")
        self.refresh()

    def set_cursor(self, ix: int, iy: int):
        self._ix, self._iy = int(ix), int(iy)
        self.refresh()

    def export_curve(self, which: str):
        """Send the EDC (``"edc"``) or MDC (``"mdc"``) now on screen out
        through :attr:`curveToList`. Returns the payload, or None if there
        is no curve yet."""
        payload = _curve_payload(self, which)
        if payload is not None:
            self.curveToList.emit(payload)
        return payload

    @staticmethod
    def _window(axis, centre_idx, half_width):
        """Indices within +/- half_width (axis units) of axis[centre_idx]."""
        centre_idx = int(np.clip(centre_idx, 0, len(axis) - 1))
        if half_width <= 0:
            return np.array([centre_idx])
        centre = axis[centre_idx]
        idxs = np.nonzero((axis >= centre - half_width) & (axis <= centre + half_width))[0]
        return idxs if idxs.size else np.array([centre_idx])

    def refresh(self, *_):
        if self._frame is None or self._xarray is None or self._yarray is None:
            self.edc_curve.setData([], [])
            self.mdc_curve.setData([], [])
            return
        xi = int(np.clip(self._ix, 0, self._frame.shape[0] - 1))
        yi = int(np.clip(self._iy, 0, self._frame.shape[1] - 1))
        energy_like = _is_energy_label(self._y_label)

        xw = self._window(self._xarray, xi, self.edc_width.value())
        edc = np.nansum(self._frame[xw, :], axis=0)
        # counts horizontal, energy vertical -- the axes are swapped relative
        # to the MDC so this panel can sit beside the image and share its
        # vertical axis.
        self.edc_curve.setData(np.asarray(edc, dtype=float),
                               np.asarray(self._yarray, dtype=float))

        yw = self._window(self._yarray, yi, self.mdc_width.value())
        mdc = np.nansum(self._frame[:, yw], axis=1)
        self.mdc_curve.setData(np.asarray(self._xarray, dtype=float),
                               np.asarray(mdc, dtype=float))

        edc_name = "EDC" if energy_like else "profile"
        mdc_name = "MDC" if energy_like else "profile"
        self.edc_button.setText(f"{edc_name} → list")
        self.mdc_button.setText(f"{mdc_name} → list")
        self.info_label.setText(
            f"{edc_name} @ {self._xarray[xi]:.4g} ({len(xw)} pt{'s' if len(xw) > 1 else ''})"
            f"   |   {mdc_name} @ {self._yarray[yi]:.4g} ({len(yw)} pt{'s' if len(yw) > 1 else ''})")


def _curve_payload(readout: "CurveReadout", which: str):
    frame, xs, ys = readout._frame, readout._xarray, readout._yarray
    if frame is None or xs is None or ys is None:
        return None
    xi = int(np.clip(readout._ix, 0, frame.shape[0] - 1))
    yi = int(np.clip(readout._iy, 0, frame.shape[1] - 1))
    if which == "edc":
        window = readout._window(xs, xi, readout.edc_width.value())
        return {"which": "edc", "x": np.asarray(ys, dtype=float).copy(),
                "y": np.nansum(np.asarray(frame, dtype=float)[window, :], axis=0),
                "x_label": readout._y_label, "position": float(xs[xi]),
                "position_label": readout._x_label,
                "half_width": float(readout.edc_width.value()),
                "n_summed": int(len(window))}
    window = readout._window(ys, yi, readout.mdc_width.value())
    return {"which": "mdc", "x": np.asarray(xs, dtype=float).copy(),
            "y": np.nansum(np.asarray(frame, dtype=float)[:, window], axis=1),
            "x_label": readout._x_label, "position": float(ys[yi]),
            "position_label": readout._y_label,
            "half_width": float(readout.mdc_width.value()),
            "n_summed": int(len(window))}



def separator() -> QFrame:
    """A thin vertical rule, for grouping a row of controls into blocks that
    read as blocks instead of one long undifferentiated strip."""
    line = QFrame()
    line.setFrameShape(QFrame.VLine)
    line.setFrameShadow(QFrame.Sunken)
    line.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Preferred)
    return line


def section_label(text: str) -> QLabel:
    """A small bold caption naming a block of controls."""
    label = QLabel(text)
    font = label.font()
    font.setBold(True)
    label.setFont(font)
    return label


def scroll_strip(widget: QWidget) -> QWidget:
    """Wrap a row of controls so a narrow window can still show it.

    A control row has a hard minimum width (the sum of its widgets), and a
    window can never be narrower than that. Several windows here size
    themselves to their data's aspect ratio and can legitimately want to be
    narrow, so the row goes in a fixed-height scroll area instead: it looks
    like an ordinary row whenever there is room, and grows a small
    horizontal scrollbar when there is not, rather than forcing the whole
    window wide.
    """
    area = QScrollArea()
    area.setWidget(widget)
    area.setWidgetResizable(True)
    area.setFrameShape(QFrame.NoFrame)
    area.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
    area.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
    area.setFixedHeight(widget.sizeHint().height() + 14)
    area.setMinimumWidth(0)
    area.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
    return area


class ViewOptionsBar(QWidget):
    """The axis / grid / alpha controls, as a visible row under the window's
    colormap rather than entries in the right-click menu.

    It replaces pyqtgraph's stock "X axis"/"Y axis" submenus and its "Plot
    Options" menu, keeping only what means something here:

    - per-axis **min/max** boxes (type a range, press Enter) and **Auto**,
      which is this app's tight fit to the data, not pyqtgraph's padded
      auto-range;
    - per-axis **Invert**;
    - **Grid**, which applies to every panel in the window (like the colormap
      above it), since it is how the window is drawn rather than where one
      panel is looking. Intensity scaling (Min%/Max%/Gamma) lives with each
      panel's own level bar instead, where the numbers it works in are.

    Dropped with the menus: "visible data only" and "auto pan only" (both
    only matter for live-updating curves), axis linking (these panels have
    different axes, so linking them is never right), and the mouse-enabled
    switch -- dragging is now always on.

    A window with more than one panel gets a panel chooser in front of the
    axis boxes; with a single panel there is nothing to choose and the
    chooser stays hidden.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.panels = []
        self._syncing = False

        row = QHBoxLayout(self)
        row.setContentsMargins(0, 0, 0, 0)
        row.setSpacing(6)

        row.addWidget(section_label("View"))
        self.panel_combo = QComboBox()
        self.panel_combo.setToolTip("Which panel the axis boxes act on.")
        self.panel_combo.currentIndexChanged.connect(lambda *_: self.sync_from_view())
        self.panel_combo.setVisible(False)
        row.addWidget(self.panel_combo)

        self.axis_boxes = []      # [(min, max), ...] per axis
        self.auto_boxes = []
        self.invert_boxes = []
        for axis, name in enumerate(("X", "Y")):
            row.addWidget(separator())
            row.addWidget(QLabel(f"{name}:"))
            spins = []
            for which in ("min", "max"):
                spin = QDoubleSpinBox()
                spin.setDecimals(4)
                spin.setRange(-1e9, 1e9)
                spin.setMaximumWidth(72)
                spin.setKeyboardTracking(False)
                # No up/down arrows: nudging a view range one step at a time
                # is not how anyone uses these, and the arrows cost more width
                # than the digits do.
                spin.setButtonSymbols(QDoubleSpinBox.NoButtons)
                spin.setToolTip(f"{name} axis {which}. Press Enter to apply.")
                spin.editingFinished.connect(
                    lambda axis=axis: self._apply_range(axis))
                row.addWidget(spin)
                spins.append(spin)
            self.axis_boxes.append(spins)

            auto = QCheckBox("Auto")
            auto.setChecked(True)
            auto.setToolTip("Keep this panel fitted tightly to its data. "
                            "Zooming or typing a range turns it off.")
            auto.toggled.connect(lambda on, axis=axis: self._apply_auto(on))
            row.addWidget(auto)
            self.auto_boxes.append(auto)

            invert = QCheckBox("Inv")
            invert.setToolTip(f"Reverse the {name} axis direction.")
            invert.toggled.connect(lambda on, axis=axis: self._apply_invert(axis, on))
            row.addWidget(invert)
            self.invert_boxes.append(invert)

        row.addWidget(separator())
        self.reset_button = QPushButton("Reset")
        self.reset_button.setToolTip("Fit the panel back to its data.")
        self.reset_button.clicked.connect(self._reset)
        row.addWidget(self.reset_button)

        row.addWidget(separator())
        self.grid_cb = QCheckBox("Grid")
        self.grid_cb.setToolTip("Draw a grid over every panel in this window.")
        self.grid_cb.toggled.connect(self._apply_grid)
        row.addWidget(self.grid_cb)

        row.addStretch(1)

    # -- binding ---------------------------------------------------------
    def bind(self, panels):
        """Point the bar at this window's panels. Called once the panels
        exist, which is after the toolbar is built."""
        self.panels = [p for p in panels if p is not None]
        self.panel_combo.blockSignals(True)
        self.panel_combo.clear()
        for panel in self.panels:
            self.panel_combo.addItem(panel.title_label.text())
        self.panel_combo.setCurrentIndex(0)
        self.panel_combo.blockSignals(False)
        self.panel_combo.setVisible(len(self.panels) > 1)
        for panel in self.panels:
            view = panel.view
            view.view_box.sigRangeChanged.connect(
                lambda *_, v=view: self._on_view_ranged(v))
            view.view_box.sigRangeChangedManually.connect(
                lambda *_, v=view: self._on_manual_range(v))
        self.sync_from_view()

    def current_view(self):
        if not self.panels:
            return None
        index = max(0, min(self.panel_combo.currentIndex(), len(self.panels) - 1))
        return self.panels[index].view

    # -- view -> boxes ---------------------------------------------------
    def sync_from_view(self):
        view = self.current_view()
        if view is None:
            return
        self._syncing = True
        try:
            for axis in range(2):
                lo, hi = view.axis_range(axis)
                for spin, value in zip(self.axis_boxes[axis], (lo, hi)):
                    spin.setValue(value)
                self.invert_boxes[axis].setChecked(view.axis_inverted(axis))
                self.auto_boxes[axis].setChecked(not view._user_ranged)
        finally:
            self._syncing = False

    def _on_view_ranged(self, view):
        if view is self.current_view() and not self._syncing:
            self._syncing = True
            try:
                for axis in range(2):
                    lo, hi = view.axis_range(axis)
                    for spin, value in zip(self.axis_boxes[axis], (lo, hi)):
                        spin.setValue(value)
            finally:
                self._syncing = False

    def _on_manual_range(self, view):
        """A drag or wheel zoom is a deliberate range, so Auto comes off --
        otherwise the next resize would silently undo what the user just
        did."""
        if view is self.current_view():
            self._syncing = True
            try:
                for box in self.auto_boxes:
                    box.setChecked(False)
            finally:
                self._syncing = False

    # -- boxes -> view ---------------------------------------------------
    def _apply_range(self, axis: int):
        view = self.current_view()
        if view is None or self._syncing:
            return
        lo = self.axis_boxes[axis][0].value()
        hi = self.axis_boxes[axis][1].value()
        current = view.axis_range(axis)
        # Leaving a box without changing it is not a request to fix the range,
        # so an unchanged value must not switch Auto off.
        if abs(lo - current[0]) < 1e-12 and abs(hi - current[1]) < 1e-12:
            return
        if view.set_axis_range(axis, lo, hi):
            self._syncing = True
            try:
                for box in self.auto_boxes:
                    box.setChecked(False)
            finally:
                self._syncing = False
        else:
            self.sync_from_view()

    def _apply_auto(self, on: bool):
        view = self.current_view()
        if view is None or self._syncing:
            return
        # The two axes cannot be fitted independently (an aspect lock ties
        # them together), so the pair moves as one.
        self._syncing = True
        try:
            for box in self.auto_boxes:
                box.setChecked(on)
        finally:
            self._syncing = False
        view.set_auto_range(on)
        self._on_view_ranged(view)

    def _apply_invert(self, axis: int, on: bool):
        view = self.current_view()
        if view is None or self._syncing:
            return
        view.set_axis_inverted(axis, on)
        if self.auto_boxes[axis].isChecked():
            view.fit_to_data()

    def _reset(self):
        view = self.current_view()
        if view is None:
            return
        view.fit_to_data()
        self._syncing = True
        try:
            for box in self.auto_boxes:
                box.setChecked(True)
        finally:
            self._syncing = False
        self._on_view_ranged(view)

    # -- whole-window options --------------------------------------------
    def _apply_grid(self, on: bool):
        for panel in self.panels:
            panel.view.set_grid(on)


def _is_energy_label(label: str) -> bool:
    """Whether an axis label denotes energy, so the curves can be called
    EDC/MDC rather than generic profiles."""
    return "eV" in (label or "") or "energy" in (label or "").lower()


class ImagePanel(QWidget):
    """An image view plus the controls that belong under it: a Smooth
    checkbox, the x/y ratio controls, and -- whenever a selection box is
    active -- the coordinates of its two draggable corners with a Refresh
    button.

    The coordinate boxes are two-way: they follow the rectangle while it is
    dragged, and typing values then pressing Refresh moves the rectangle
    there and re-runs the integration. That makes an exact region (e.g.
    x 110..120, y 110..120) reproducible, which dragging by eye is not.

    ``selectionApplied`` fires for both the Refresh button and the
    right-click menu entry, so the owner only wires one signal.
    """

    selectionApplied = pyqtSignal()
    #: an EDC or MDC to list (see :attr:`CurveReadout.curveToList`)
    curveToList = pyqtSignal(dict)

    def __init__(self, title: str, view: "_InteractiveImageBase",
                 equal_ratio_option: bool = False, curves: bool = False,
                 lock_ratio: bool = False, parent=None):
        super().__init__(parent)
        self.view = view
        view.setParent(self)
        # The panel's own title is what an export of it is called; the
        # window fills in the file name and the slice position around it.
        view.export_title = title

        layout = QVBoxLayout(self)
        layout.setContentsMargins(2, 2, 2, 2)
        layout.setSpacing(3)

        self.title_label = QLabel(title)
        self.title_label.setAlignment(Qt.AlignCenter)
        self.title_label.setWordWrap(True)
        layout.addWidget(self.title_label)

        # -- min/max level bar, directly above the image --
        self.level_bar = LevelBar(self)
        # In a scroll strip: its numeric boxes must not become the panel's
        # (and so the window's) minimum width -- a contour window has to stay
        # free to open narrow, at its data's aspect.
        layout.addWidget(scroll_strip(self.level_bar))
        self.level_bar.levelsChanged.connect(self._on_levels_changed)
        self.level_bar.gammaChanged.connect(self.view.set_gamma)

        # -- image, with the EDC to its right and the MDC underneath --
        # A grid, not nested boxes: the curves have to share a cell edge with
        # the image for their linked axes to line up on screen.
        self.curves = CurveReadout(self) if curves else None
        if self.curves is not None:
            self.curves.curveToList.connect(self.curveToList.emit)
        # The image and its curves live in a widget of their own so that the
        # whole block can be reshaped as one: a locked ratio is enforced by
        # giving the plotting area the data's shape (see AspectBox), and on a
        # curves panel the image cannot be resized alone without unpinning it
        # from the EDC and MDC beside it. The layout divides the block in
        # fixed proportions, so resizing the block gets the image right.
        plot_cell = QWidget()
        plot_grid = QGridLayout(plot_cell)
        plot_grid.setContentsMargins(0, 0, 0, 0)
        plot_grid.setSpacing(2)
        plot_grid.addWidget(view, 0, 0)
        plot_grid.setColumnStretch(0, 4)
        plot_grid.setRowStretch(0, 4)
        if self.curves is not None:
            plot_grid.addWidget(self.curves.edc_plot, 0, 1)
            plot_grid.addWidget(self.curves.mdc_plot, 1, 0)
            # The side panels get a quarter of the space, so the image itself
            # stays the thing you are actually looking at.
            plot_grid.setColumnStretch(1, 1)
            plot_grid.setRowStretch(1, 1)
            self.curves.link_to(view)
            _pin_axis_geometry(view, self.curves)
        self.aspect_box = AspectBox(plot_cell, view, self)
        view.attach_aspect_box(self.aspect_box)
        layout.addWidget(self.aspect_box, stretch=1)

        # -- the readout cursor's position, typed (shown while it is on) --
        self._build_cursor_row(layout)

        # -- display options row --
        options = QHBoxLayout()
        options.setContentsMargins(0, 0, 0, 0)
        options.setSpacing(6)
        self.smooth_cb = QCheckBox("Smooth")
        self.smooth_cb.setToolTip(
            "Draw the image interpolated between data points instead of as one "
            "rectangle per point -- the same thing MATLAB's \"shading interp\" does, "
            "and what a coarse deflector map needs to stop looking like a mosaic. "
            "Display only: exported data stays raw.")
        options.addWidget(self.smooth_cb)

        options.addWidget(QLabel("blur"))
        self.smooth_sigma = QDoubleSpinBox()
        self.smooth_sigma.setRange(0.0, 10.0)
        self.smooth_sigma.setSingleStep(0.2)
        self.smooth_sigma.setValue(0.0)
        self.smooth_sigma.setDecimals(1)
        self.smooth_sigma.setMaximumWidth(56)
        self.smooth_sigma.setToolTip(
            "Extra Gaussian blur on top of the interpolation, in data pixels. "
            "0 keeps every measured feature and only removes the blockiness; "
            "raise it to average down noise.")
        self.smooth_sigma.setEnabled(False)
        options.addWidget(self.smooth_sigma)

        if self.curves is not None:
            self.curves_separator = separator()
            options.addWidget(self.curves_separator)
            options.addWidget(self.curves.controls)
        options.addStretch(1)
        options.addWidget(separator())
        self.ratio_cb = QCheckBox("Equal ratio" if equal_ratio_option else "Lock x/y")
        self.ratio_cb.setToolTip(
            "Lock the on-screen x:y scale to the ratio on the right "
            "(1 = one x unit occupies the same space as one y unit). "
            "Unchecked, the image stretches to fill the panel.")
        options.addWidget(self.ratio_cb)

        self.ratio_spin = QDoubleSpinBox()
        self.ratio_spin.setRange(0.001, 1000.0)
        self.ratio_spin.setDecimals(3)
        self.ratio_spin.setSingleStep(0.1)
        self.ratio_spin.setValue(1.0)
        self.ratio_spin.setMaximumWidth(80)
        self.ratio_spin.setToolTip("x:y ratio applied when the lock is on.")
        self.ratio_spin.setEnabled(False)
        options.addWidget(self.ratio_spin)
        # Panels whose two axes carry the same physical unit (a spatial map,
        # a constant-E contour) start locked at 1:1 rather than stretched.
        if lock_ratio:
            self.ratio_cb.setChecked(True)
            self.ratio_spin.setEnabled(True)
        layout.addLayout(options)

        # -- selection coordinates row (hidden until a box exists) --
        self.selection_widget = QWidget()
        sel = QHBoxLayout(self.selection_widget)
        sel.setContentsMargins(0, 0, 0, 0)
        sel.setSpacing(6)
        sel.addWidget(section_label("Selection"))
        self.corner_spins = []
        for label in ("x0", "y0", "x1", "y1"):
            sel.addWidget(QLabel(label))
            spin = QDoubleSpinBox()
            spin.setDecimals(4)
            spin.setRange(-1e9, 1e9)
            spin.setMaximumWidth(90)
            sel.addWidget(spin)
            self.corner_spins.append(spin)
        self.refresh_button = QPushButton("Refresh")
        self.refresh_button.setToolTip(
            "Apply the coordinates above to the selection box and update the data.")
        sel.addWidget(self.refresh_button)
        sel.addStretch(1)
        layout.addWidget(self.selection_widget)
        self.selection_widget.setVisible(False)

        self.smooth_cb.toggled.connect(self._on_smooth_changed)
        self.smooth_sigma.valueChanged.connect(self._on_smooth_changed)
        self.ratio_cb.toggled.connect(self._on_ratio_changed)
        self.ratio_spin.valueChanged.connect(self._on_ratio_changed)
        self.refresh_button.clicked.connect(self._on_refresh_clicked)
        if self.curves is not None:
            for widget in self._curve_widgets():
                widget.setVisible(False)

        view.selectionChanged.connect(self._sync_from_selection)
        view.selectionApplied.connect(self.selectionApplied.emit)
        view.readoutToggled.connect(self._on_cursor_toggled)
        view.readoutMoved.connect(lambda *_: self._sync_cursor_row())
        self.set_axis_colors(*view._readout_colors)
        if self.curves is not None:
            view.readoutToggled.connect(self._on_readout_toggled)
            view.readoutMoved.connect(self.curves.set_cursor)
            view.readoutMoved.connect(lambda *_: self._sync_readout_bands())
            self.curves.edc_width.valueChanged.connect(self._sync_readout_bands)
            self.curves.mdc_width.valueChanged.connect(self._sync_readout_bands)

    # -- the cursor row -----------------------------------------------------
    def _build_cursor_row(self, layout):
        """``Cursor  <x name> [____]  <y name> [____]  value  <message>``.

        The two boxes follow the cursor as it is dragged; typing a value and
        pressing Enter moves the cursor there (and, on a map, the windows
        linked to it). A value outside the data leaves the cursor where it
        is and says so in the message, in red.
        """
        self.cursor_widget = QWidget()
        row = QHBoxLayout(self.cursor_widget)
        row.setContentsMargins(0, 0, 0, 0)
        row.setSpacing(6)
        row.addWidget(section_label("Cursor"))
        self.cursor_labels, self.cursor_edits = [], []
        for axis in (0, 1):
            label = QLabel(("x", "y")[axis])
            edit = QLineEdit()
            edit.setMaximumWidth(96)
            edit.setToolTip("Type a position and press Enter to move the cursor there.")
            edit.returnPressed.connect(lambda axis=axis: self._on_cursor_typed(axis))
            row.addWidget(label)
            row.addWidget(edit)
            self.cursor_labels.append(label)
            self.cursor_edits.append(edit)
        self.cursor_value = QLabel("")
        row.addWidget(self.cursor_value)
        self.cursor_message = QLabel("")
        self.cursor_message.setStyleSheet("color: #c0392b;")
        self.cursor_message.setWordWrap(True)
        row.addWidget(self.cursor_message, 1)
        layout.addWidget(self.cursor_widget)
        self.cursor_widget.setVisible(False)

    def set_axis_colors(self, x_color: str, y_color: str):
        """Colour the cursor lines, their bands, the EDC/MDC they feed and
        the names in the cursor row by axis (see ``AXIS_COLORS``)."""
        self.view.set_readout_colors(x_color, y_color)
        if self.curves is not None:
            self.curves.set_colors(x_color, y_color)
        for label, color in zip(self.cursor_labels, (x_color, y_color)):
            label.setStyleSheet(f"color: {color}; font-weight: bold;")

    def _on_cursor_toggled(self, visible: bool):
        self.cursor_widget.setVisible(visible)
        self.cursor_message.setText("")
        if visible:
            self._sync_cursor_row()

    def _sync_cursor_row(self):
        """Show where the cursor is -- the data point it reads, not the
        pixel the mouse let go of."""
        indices = self.view.readout_indices()
        if indices is None:
            return
        view = self.view
        names = (view.x_label or "x", view.y_label or "y")
        for axis, (label, edit) in enumerate(zip(self.cursor_labels, self.cursor_edits)):
            label.setText(names[axis])
            values = view._xarray if axis == 0 else view._yarray
            if not edit.hasFocus() or not edit.isModified():
                edit.setText(f"{float(values[indices[axis]]):.6g}")
                edit.setModified(False)
        value = view.readout_value()
        self.cursor_value.setText("" if value is None else f"value {value:.6g}")

    def _on_cursor_typed(self, axis: int):
        edit = self.cursor_edits[axis]
        text = edit.text().strip().replace(",", ".")
        try:
            value = float(text)
        except ValueError:
            message = f"“{edit.text()}” is not a number."
        else:
            message = self.view.move_readout_to(**{("x", "y")[axis]: value})
        edit.setModified(False)
        self.cursor_message.setText(message or "")
        if message:
            window = self.window()
            status = getattr(window, "statusBar", None)
            if callable(status):
                status().showMessage(message, 8000)
        self._sync_cursor_row()
        return message

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self.view.refit_if_untouched()

    def showEvent(self, event):
        super().showEvent(event)
        self.view.refit_if_untouched()

    # -- levels ----------------------------------------------------------
    def _on_levels_changed(self, lo: float, hi: float):
        self.view.set_level_window(lo, hi)

    def sync_level_range(self):
        """Point the level bar at the current image's value range. Called
        after new data lands, when setImage's autoLevels has just reset the
        image's own levels to the full range."""
        frame = getattr(self.view, "last_frame", None)
        if frame is None:
            frame = getattr(self.view, "_image2d", None)
        if frame is None:
            return
        finite = np.asarray(frame, dtype=float)
        finite = finite[np.isfinite(finite)]
        if finite.size == 0:
            return
        self.level_bar.set_data_range(float(finite.min()), float(finite.max()))
        self.level_bar.set_histogram(finite)
        # New data means a new window; tell the view, so gamma keeps being
        # applied inside the levels actually on screen.
        self.view.set_level_window(*self.level_bar.levels())

    # -- EDC / MDC -------------------------------------------------------
    def _on_readout_toggled(self, visible: bool):
        """The curves only mean something while the cursor that defines them
        is on screen, so they appear and disappear with it."""
        if self.curves is None:
            return
        if visible:
            self.sync_curve_source()
        for widget in self._curve_widgets():
            widget.setVisible(visible)
        if visible:
            self._sync_readout_bands()

    def _curve_widgets(self):
        """The EDC/MDC plots, their controls, and the rule that separates
        those controls from the rest of the row -- all of which appear and
        disappear with the readout cursor."""
        if self.curves is None:
            return ()
        return tuple(self.curves.widgets()) + (self.curves_separator,)

    def _sync_readout_bands(self, *_):
        """Mirror the EDC/MDC integration widths onto the image as shaded
        bands around the crosshair."""
        if self.curves is None:
            return
        self.view.set_readout_bands(self.curves.edc_width.value(),
                                    self.curves.mdc_width.value())

    def sync_curve_source(self):
        """Re-point the curves at whatever the view is currently showing.
        Called when the readout is switched on and after new data lands."""
        if self.curves is None:
            return
        frame = getattr(self.view, "last_frame", None)
        self.curves.set_source(frame, self.view._xarray, self.view._yarray,
                               self.view.x_label or "x", self.view.y_label or "y")

    # -- display options -------------------------------------------------
    def _on_smooth_changed(self, *_):
        on = self.smooth_cb.isChecked()
        self.smooth_sigma.setEnabled(on)
        self.view.set_smoothing(on, self.smooth_sigma.value())

    def _on_ratio_changed(self, *_):
        on = self.ratio_cb.isChecked()
        self.ratio_spin.setEnabled(on)
        self.view.set_aspect(on, self.ratio_spin.value())
        # Re-fit: locking stretches an axis and unlocking should give the
        # space straight back, rather than leaving the view as the previous
        # setting left it.
        self.view.fit_to_data()

    def apply_display_options(self):
        """Re-assert smoothing and aspect after new data is loaded (a fresh
        setImage with autoRange resets the view's aspect state), and re-point
        the EDC/MDC curves at the new frame."""
        self._on_smooth_changed()
        self._on_ratio_changed()
        self.sync_level_range()
        if self.curves is not None and self.curves.edc_plot.isVisibleTo(self):
            self.sync_curve_source()

    # -- selection -------------------------------------------------------
    def _sync_from_selection(self):
        corners = self.view.selection_corners()
        self.selection_widget.setVisible(corners is not None)
        if corners is None:
            return
        for spin, value in zip(self.corner_spins, corners):
            spin.blockSignals(True)
            spin.setValue(value)
            spin.blockSignals(False)

    def _on_refresh_clicked(self):
        values = [spin.value() for spin in self.corner_spins]
        self.view.set_selection_corners(*values)
        self.selectionApplied.emit()

    def set_selection_enabled(self, enabled: bool, text: str = None):
        """Forwarded to the view; also hides the coordinate row for panels
        where region integration isn't defined."""
        self.view.enable_selection_action(enabled, text)
        self.refresh_button.setVisible(enabled)


# --------------------------------------------------------------------------
# Angle-cube explorer (the "3D" visualization for `map` data)
# --------------------------------------------------------------------------
class _SliceControl(QWidget):
    """A labeled slider with typed **Pos** and **Ind** boxes and a "+/-"
    integration half-width, sized to sit directly above the plot panel it
    drives (so it inherits that panel's column width exactly).

    Pos and Ind are the two ways of naming the same slice, as in the MATLAB
    tool: type a value in data units (eV, degrees) and the nearest
    measured point is shown; type an index and that point is shown. Either
    box, and the ``<`` / ``>`` step buttons, move the slider, so there is
    one source of truth. Indices are **1-based**, matching the MATLAB
    tool's "Ind" box, so the same number means the same slice in both.

    The spin box is ``width_spin``, not ``width``: a ``width`` attribute
    would shadow ``QWidget.width()`` and silently break any layout code
    that asks this widget how wide it is.
    """

    changed = pyqtSignal()

    def __init__(self, title: str, unit: str, width_step: float, parent=None):
        super().__init__(parent)
        self.unit = unit
        self.axis = None
        layout = QVBoxLayout(self)
        layout.setContentsMargins(2, 2, 2, 2)
        layout.setSpacing(2)

        self.title_label = QLabel(title)
        self.title_label.setAlignment(Qt.AlignCenter)
        layout.addWidget(self.title_label)

        self.readout = QLabel("-")
        self.readout.setAlignment(Qt.AlignCenter)
        layout.addWidget(self.readout)

        row = QHBoxLayout()
        row.setContentsMargins(0, 0, 0, 0)
        row.setSpacing(4)
        self.slider = QSlider(Qt.Horizontal)
        row.addWidget(self.slider, stretch=1)

        row.addWidget(QLabel("Pos"))
        self.pos_box = QDoubleSpinBox()
        self.pos_box.setDecimals(4)
        self.pos_box.setRange(-1e9, 1e9)
        self.pos_box.setKeyboardTracking(False)
        self.pos_box.setButtonSymbols(QDoubleSpinBox.NoButtons)
        self.pos_box.setMaximumWidth(86)
        self.pos_box.setToolTip(
            f"Slice position in {unit}. Type a value and press Enter -- the "
            "nearest measured point is shown.")
        row.addWidget(self.pos_box)

        row.addWidget(QLabel("Ind"))
        self.prev_button = QPushButton("<")
        self.next_button = QPushButton(">")
        for button in (self.prev_button, self.next_button):
            button.setFixedWidth(24)
            button.setToolTip("Step one point along the axis.")
        self.ind_box = QSpinBox()
        self.ind_box.setRange(1, 1)
        self.ind_box.setKeyboardTracking(False)
        self.ind_box.setButtonSymbols(QSpinBox.NoButtons)
        self.ind_box.setMaximumWidth(64)
        self.ind_box.setToolTip(
            "Slice index, 1-based as in the MATLAB tool. Type and press Enter.")
        row.addWidget(self.prev_button)
        row.addWidget(self.ind_box)
        row.addWidget(self.next_button)

        row.addWidget(QLabel("+/-"))
        self.width_spin = QDoubleSpinBox()
        self.width_spin.setDecimals(3)
        self.width_spin.setSingleStep(width_step)
        self.width_spin.setRange(0.0, 1e6)
        self.width_spin.setValue(0.0)
        self.width_spin.setMaximumWidth(80)
        self.width_spin.setToolTip(
            f"Integration half-width in {unit}: every point within +/- this "
            "much of the position goes into the slice. The slider's range "
            "shrinks to match, so the window always stays inside the data.")
        row.addWidget(self.width_spin)

        self.reduce_combo = QComboBox()
        self.reduce_combo.addItems(["mean", "sum"])
        self.reduce_combo.setMaximumWidth(72)
        self.reduce_combo.setToolTip(
            "How the points in the window are combined. mean keeps the "
            "brightness the same however wide the window is, so the levels "
            "stay put (this is what the MATLAB tool does); sum gives total "
            "counts, which doubles when you double the window.")
        row.addWidget(self.reduce_combo)
        layout.addLayout(row)

        self.slider.valueChanged.connect(self._on_slider_moved)
        self.width_spin.valueChanged.connect(self._on_width_changed)
        self.reduce_combo.currentTextChanged.connect(lambda *_: self.changed.emit())
        self.pos_box.editingFinished.connect(self._on_pos_typed)
        self.ind_box.editingFinished.connect(self._on_index_typed)
        self.prev_button.clicked.connect(lambda: self.step(-1))
        self.next_button.clicked.connect(lambda: self.step(+1))

    def configure(self, n: int, axis=None):
        self.axis = None if axis is None else np.asarray(axis, dtype=float)
        self.slider.setRange(0, max(n - 1, 0))
        self.ind_box.setRange(1, max(n, 1))
        self.slider.setValue(n // 2)
        self._sync_boxes()

    @property
    def index(self) -> int:
        return self.slider.value()

    # -- keeping slider, Pos and Ind in step ------------------------------
    def _on_slider_moved(self, *_):
        self._sync_boxes()
        self.changed.emit()

    def _sync_boxes(self):
        """Show the slider's slice in both boxes without re-triggering a
        redraw -- they are a readout as much as an input."""
        for box in (self.pos_box, self.ind_box):
            box.blockSignals(True)
        try:
            self.ind_box.setValue(self.index + 1)
            if self.axis is not None and self.axis.size:
                idx = int(np.clip(self.index, 0, self.axis.size - 1))
                self.pos_box.setValue(float(self.axis[idx]))
        finally:
            for box in (self.pos_box, self.ind_box):
                box.blockSignals(False)

    def set_index(self, index: int):
        index = int(np.clip(index, self.slider.minimum(), self.slider.maximum()))
        if index == self.index:
            self._sync_boxes()      # snap a typed value back onto its point
            return
        self.slider.setValue(index)

    def step(self, delta: int):
        self.set_index(self.index + int(delta))

    def _on_pos_typed(self):
        """A typed position is matched to the nearest measured point: the
        data only exists at those points, so showing anything else would be
        a lie about which slice is on screen."""
        if self.axis is None or not self.axis.size:
            return
        target = self.pos_box.value()
        self.set_index(int(np.argmin(np.abs(self.axis - target))))

    def _on_index_typed(self):
        self.set_index(self.ind_box.value() - 1)

    @property
    def half_width(self) -> float:
        return self.width_spin.value()

    @property
    def reduce(self) -> str:
        return self.reduce_combo.currentText()

    def combine(self, values, axis: int):
        """Collapse the integration window the way the combo box says."""
        values = np.asarray(values)
        return values.mean(axis=axis) if self.reduce == "mean" else values.sum(axis=axis)

    # -- integration window ------------------------------------------------
    def _on_width_changed(self, *_):
        self._apply_width_limits()
        self.changed.emit()

    def _apply_width_limits(self):
        """Keep the integration window inside the measured range.

        With a half-width of w, a slice centred nearer than w to either end
        would silently be integrated over fewer points than one in the
        middle -- the same slice, quietly weaker. The MATLAB tool pulls the
        slider's own limits in by w instead, which is what this does; a
        half-width wider than half the axis is clamped rather than refused,
        so the box always shows a width that is actually in use.
        """
        if self.axis is None or self.axis.size < 2:
            return
        lo_v, hi_v = float(self.axis.min()), float(self.axis.max())
        width = self.width_spin.value()
        limit = (hi_v - lo_v) / 2.0
        if width > limit:
            width = limit
            self.width_spin.blockSignals(True)
            self.width_spin.setValue(width)
            self.width_spin.blockSignals(False)
        tol = 1e-9 * max(1.0, abs(hi_v - lo_v))
        valid = np.nonzero((self.axis >= lo_v + width - tol)
                           & (self.axis <= hi_v - width + tol))[0]
        if valid.size == 0:
            valid = np.array([int(np.argmin(np.abs(self.axis - (lo_v + hi_v) / 2)))])
        first, last = int(valid.min()), int(valid.max())
        self.slider.setRange(first, last)
        self.ind_box.setRange(first + 1, last + 1)
        self._sync_boxes()

    def show_value(self, center, idxs, axis):
        if len(idxs) == 1:
            self.readout.setText(f"{center:.4g} {self.unit}")
        else:
            self.readout.setText(
                f"{center:.4g} {self.unit} ({self.reduce} "
                f"{axis[idxs[0]]:.4g}~{axis[idxs[-1]]:.4g}, {len(idxs)} pts)")


class FrameWindow(QWidget):
    """A standalone window holding one frozen E-vs-k frame.

    Opened from the spatial panel's right-click menu to snapshot the
    spectrum at the cursor's current position. Several can be open at once
    and each keeps its own copy of the data, so they stay put while the main
    window's cursor moves on -- which is the whole point: comparing several
    sample positions side by side.

    The window carries the same ImagePanel controls as the main panels
    (Smooth, x/y ratio, readout cursor with EDC/MDC), so a snapshot is as
    inspectable as the live view.
    """

    #: emitted on close so the owner can drop its reference
    closed = pyqtSignal(object)

    def __init__(self, title: str, frame, x_axis, y_axis,
                 x_label=None, y_label=None, colormap=None, flip=False):
        super().__init__()
        self.setWindowTitle(title)
        self.resize(560, 620)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)

        # Its own colormap control: giving one snapshot a different scale from
        # another is exactly what makes two positions easy to tell apart.
        self.colormap = colormap or DEFAULT_COLORMAP
        self.flip = flip
        bar = QHBoxLayout()
        bar.setSpacing(6)
        bar.addWidget(section_label("Colormap"))
        self.colormap_combo = QComboBox()
        self.colormap_combo.addItems(COLORMAP_NAMES)
        self.colormap_combo.setCurrentText(self.colormap)
        self.flip_cb = QCheckBox("Flip")
        self.flip_cb.setChecked(self.flip)
        bar.addWidget(self.colormap_combo)
        bar.addWidget(self.flip_cb)
        bar.addWidget(separator())
        self.figure_button = QPushButton("Plot tools")
        self.figure_button.setToolTip(
            "Turn this panel into a figure for a paper: panel letters, a "
            "colour bar, editable labels, annotations, a second axis, and "
            "export at a journal's column width.")
        self.figure_button.clicked.connect(self.open_figure)
        bar.addWidget(self.figure_button)
        bar.addStretch(1)
        layout.addLayout(bar)
        #: {"entries": callable, "loader": callable} so the figure composer
        #: can take further panels from the launcher's list
        self.figure_source = None
        self.figure_windows = []

        self.view = FrameImageView()
        self.panel = ImagePanel(title, self.view, curves=True, parent=self)
        self.view.enable_selection_action(False)
        # Same axis/grid/alpha row as the main viewers, directly under the
        # colormap, so a snapshot is as adjustable as the live view.
        self.view_bar = ViewOptionsBar(self)
        layout.addWidget(scroll_strip(self.view_bar))
        layout.addWidget(self.panel)

        self.colormap_combo.currentTextChanged.connect(
            lambda *_: self.set_colormap(self.colormap_combo.currentText(),
                                          self.flip_cb.isChecked()))
        self.flip_cb.toggled.connect(
            lambda *_: self.set_colormap(self.colormap_combo.currentText(),
                                          self.flip_cb.isChecked()))

        self.view.set_axis_labels(x_label, y_label)
        # A frozen panel is named after itself: its title already carries the
        # file and the position it was taken at.
        self.view.export_prefix = title
        self.view.export_title = ""
        # a copy, not a reference: the main view's array is replaced whenever
        # the cursor moves, and this window must not follow it
        self.view.set_frame(np.array(frame, copy=True),
                            np.array(x_axis, copy=True), np.array(y_axis, copy=True))
        if colormap:
            apply_colormap(self.view, colormap, flip)
        self.panel.apply_display_options()
        self.view_bar.bind([self.panel])

    def set_colormap(self, name: str, flip: bool = False):
        self.colormap, self.flip = name, flip
        apply_colormap(self.view, name, flip)

    def open_figure(self):
        """Hand this panel to the figure composer.

        Imported here rather than at the top of the module: the composer
        pulls in Qt's PDF and SVG writers, and a viewer that never makes a
        figure should not pay for them.
        """
        from ui.figure import FigureWindow, panel_from_arrays

        panel = panel_from_arrays(
            self.view.last_frame, self.view._xarray, self.view._yarray,
            self.view.x_label or "", self.view.y_label or "",
            name=self.windowTitle(), info=dict(getattr(self, "source_info", {})),
            colormap=self.colormap, flip=self.flip,
            smooth=self.panel.smooth_cb.isChecked())
        panel.levels = self.view._levels
        panel.gamma = self.view._gamma
        window = FigureWindow([panel], f"Figure — {self.windowTitle()}",
                              dataset_source=self.figure_source)
        window.closed.connect(lambda w: self.figure_windows.remove(w)
                              if w in self.figure_windows else None)
        self.figure_windows.append(window)
        window.show()
        return window

    def closeEvent(self, event):
        for window in list(self.figure_windows):
            window.close()
        self.closed.emit(self)
        super().closeEvent(event)


class KCubeExplorer(QWidget):
    """Three linked orthogonal slice views of a (deflector angle, slit
    angle, E) cube, laid out as three equal columns:

    - column 1: constant-E map      (deflector angle vs slit angle)
    - column 2: deflector angle vs E, integrated over a slit-angle window
    - column 3: slit angle vs E, integrated over a deflector-angle window

    Each column's integration control sits directly above that column's
    plot and therefore shares its exact width.

    Deliberately labeled in raw **angles**, not kx/ky: converting the
    deflector axis to momentum needs a reference kinetic energy and a
    small-angle approximation (see loader.nxs_file.deflector_angle_to_k),
    a modeling choice better made explicitly during analysis than baked
    into the default viewer. Use ``NxsData.kcube()`` for the converted
    version in a script.

    The "+/-" box on each control is a half-width in physical units: 0
    gives a single-pixel slice, while e.g. 0.25 at a slider position of
    88 eV sums every channel in [87.75, 88.25] eV -- the usual way to cut
    noise in a constant-energy map.

    Needs no OpenGL. If ``pyqtgraph.opengl`` + PyOpenGL are importable, a
    volumetric rendering is offered additionally via
    :meth:`try_build_gl_volume`; the slice viewer alone already gives full
    access to the cube.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.angle_defl = self.angle_slit = self.E = None
        self.cube = None  # (len(angle_defl), len(angle_slit), len(E))

        grid = QGridLayout(self)
        for col in range(3):
            grid.setColumnStretch(col, 1)
        grid.setRowStretch(0, 0)   # slice controls: natural height
        grid.setRowStretch(1, 1)   # ImagePanels take the rest

        self.energy_map = FrameImageView()
        self.defl_e_cut = FrameImageView()
        self.slit_e_cut = FrameImageView()

        # Wrapped in ImagePanels so each column carries its own Smooth and
        # x/y-ratio controls under the plot. Only the constant-E map has two
        # axes in the same units (both angles), so it is the one that gets
        # the explicit "Equal ratio" tickbox; the E cuts mix degrees with eV,
        # where a 1:1 lock would be meaningless, so they get the general
        # ratio lock instead.
        self.energy_panel = ImagePanel(
            "Constant-E map: angle (deflector) vs angle (along slit)",
            self.energy_map, equal_ratio_option=True, parent=self)
        self.defl_panel = ImagePanel("Angle (deflector) vs E", self.defl_e_cut,
                                     curves=True, parent=self)
        self.slit_panel = ImagePanel("Angle (along slit) vs E", self.slit_e_cut,
                                     curves=True, parent=self)

        # The control above each column collapses the axis that column's
        # panel does NOT show, which is what makes the label meaningful.
        self.e_control = _SliceControl("Energy (eV)", "eV", 0.05, self)
        self.slit_control = _SliceControl("Integrate over angle (along slit)", "deg", 0.1, self)
        self.defl_control = _SliceControl("Integrate over angle (deflector)", "deg", 0.1, self)

        # Row 0 is the slice control, row 1 the ImagePanel (which carries its
        # own title and display controls), so control and plot share a column
        # and therefore an exact width.
        for col, (control, panel) in enumerate([
            (self.e_control, self.energy_panel),
            (self.slit_control, self.defl_panel),
            (self.defl_control, self.slit_panel),
        ]):
            grid.addWidget(control, 0, col)
            grid.addWidget(panel, 1, col)

        self.e_control.changed.connect(self._refresh_e)
        self.slit_control.changed.connect(self._refresh_defl_cut)
        self.defl_control.changed.connect(self._refresh_slit_cut)

        self.gl_view = None

    def set_cube(self, angle_defl, angle_slit, E, cube):
        self.angle_defl = np.asarray(angle_defl)
        self.angle_slit = np.asarray(angle_slit)
        self.E = np.asarray(E)
        self.cube = cube
        self.e_control.configure(len(self.E), self.E)
        self.slit_control.configure(len(self.angle_slit), self.angle_slit)
        self.defl_control.configure(len(self.angle_defl), self.angle_defl)
        self._refresh_e()
        self._refresh_defl_cut()
        self._refresh_slit_cut()
        # A fresh setImage with autoRange drops the aspect lock, so put the
        # user's Smooth/ratio choices back after reloading the cube.
        for panel in self.image_panels():
            panel.apply_display_options()
        self.try_build_gl_volume()

    @staticmethod
    def _index_range(axis: np.ndarray, center_idx: int, half_width: float) -> np.ndarray:
        """Indices to sum over for a center index +/- a physical half-width
        (same units as ``axis``); half_width<=0 gives just the center."""
        if half_width <= 0:
            return np.array([center_idx])
        center_val = axis[center_idx]
        idxs = np.nonzero((axis >= center_val - half_width) & (axis <= center_val + half_width))[0]
        return idxs if idxs.size else np.array([center_idx])

    def _refresh_e(self):
        if self.cube is None:
            return
        ei = self.e_control.index
        idxs = self._index_range(self.E, ei, self.e_control.half_width)
        self.e_control.show_value(self.E[ei], idxs, self.E)
        self.energy_map.set_frame(self.e_control.combine(self.cube[:, :, idxs], 2),
                                   self.angle_defl, self.angle_slit)
        self.energy_panel.sync_level_range()

    def _refresh_defl_cut(self):
        """Deflector-vs-E panel: collapses the slit-angle axis."""
        if self.cube is None:
            return
        si = self.slit_control.index
        idxs = self._index_range(self.angle_slit, si, self.slit_control.half_width)
        self.slit_control.show_value(self.angle_slit[si], idxs, self.angle_slit)
        self.defl_e_cut.set_frame(self.slit_control.combine(self.cube[:, idxs, :], 1),
                                  self.angle_defl, self.E)
        self.defl_panel.sync_curve_source()
        self.defl_panel.sync_level_range()

    def _refresh_slit_cut(self):
        """Slit-vs-E panel: collapses the deflector-angle axis."""
        if self.cube is None:
            return
        di = self.defl_control.index
        idxs = self._index_range(self.angle_defl, di, self.defl_control.half_width)
        self.defl_control.show_value(self.angle_defl[di], idxs, self.angle_defl)
        self.slit_e_cut.set_frame(self.defl_control.combine(self.cube[idxs, :, :], 0),
                                  self.angle_slit, self.E)
        self.slit_panel.sync_curve_source()
        self.slit_panel.sync_level_range()

    def current_energy_map_frame(self):
        """The (deflector, slit) frame currently shown in the constant-E
        panel -- used for 'save current view' in ARPES_viewer.py."""
        return self.energy_map.last_frame

    def set_axis_labels(self, defl_label: str, slit_label: str, energy_label: str):
        """Label all three panels from the scan's axis names/units. Each
        panel shows a different pair, so the labels are distributed rather
        than applied uniformly."""
        self.energy_map.set_axis_labels(defl_label, slit_label)
        self.defl_e_cut.set_axis_labels(defl_label, energy_label)
        self.slit_e_cut.set_axis_labels(slit_label, energy_label)

    def image_panels(self):
        return (self.energy_panel, self.defl_panel, self.slit_panel)

    def panels(self):
        return (self.energy_map, self.defl_e_cut, self.slit_e_cut)

    def set_colormap(self, name: str, flip: bool = False):
        for view in self.panels():
            apply_colormap(view, name, flip)

    def try_build_gl_volume(self):
        """Best-effort bonus: a real 3D volume rendering, if
        pyqtgraph.opengl + PyOpenGL happen to be installed. Skipped
        silently otherwise -- the slice viewer already gives full access."""
        try:
            import pyqtgraph.opengl as gl
        except ImportError:
            return
        if self.gl_view is None:
            self.gl_view = gl.GLViewWidget()
            self.gl_view.setWindowTitle("Deflector angle / slit angle / E volume (bonus 3D view)")
            self.gl_view.opts["distance"] = 200
        self.gl_view.clear()
        vol = np.nan_to_num(self.cube, nan=0.0)
        vol = vol - vol.min()
        vmax = vol.max()
        if vmax > 0:
            vol = vol / vmax
        rgba = np.zeros(vol.shape + (4,), dtype=np.ubyte)
        rgba[..., 0] = (vol * 255).astype(np.ubyte)
        rgba[..., 3] = (np.clip(vol, 0, 1) ** 1.5 * 180).astype(np.ubyte)
        item = gl.GLVolumeItem(rgba, sliceDensity=1)
        item.translate(-vol.shape[0] / 2, -vol.shape[1] / 2, -vol.shape[2] / 2)
        self.gl_view.addItem(item)
        self.gl_view.show()
