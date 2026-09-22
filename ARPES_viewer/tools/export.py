"""
tools/export.py
=============
What leaves the program, and in what shape.

Three things can be exported from a viewer panel, and they exist because a
figure is built from three different kinds of file:

``.nxs``
    The panel's own data, written by ``loader.nxs_file.save_dataset`` so this
    program opens it again as a dataset like any measured one. This is the
    archival export -- axes, labels and the file's metadata all travel with
    the numbers.

an image
    Exactly what the panel draws -- the same colormap, the same level
    window, the same gamma -- rendered at whatever pixel size is asked for
    rather than at the size the window happens to be.  :func:`render_rgba`
    resamples the data onto the requested pixel grid itself instead of
    scaling a screenshot, so a 3000-px export is genuinely 3000 px of data,
    not a blown-up 600-px one.  Points outside the data become transparent.

an ``.eps`` of the axes alone
    The frame, the ticks, the tick labels and the axis titles as vector
    PostScript, with nothing inside the box. This is the half of a figure
    that has to stay editable: in Illustrator the text is text and the
    ticks are paths, and the raster above drops into the empty box.  The
    two exports are generated from the same ranges, so they line up when
    the image is placed at the box :func:`axes_eps` reports.

Kept free of Qt (like ``tools.analysis``/``tools.fermi``) so the rendering and
the PostScript can be tested headlessly.
"""
from __future__ import annotations

import math
import time

import numpy as np

__all__ = ["nice_step", "tick_values", "tick_format", "render_rgba",
           "write_image", "axes_eps", "IMAGE_FORMATS", "MM_PER_INCH", "PT_PER_MM"]

MM_PER_INCH = 25.4
PT_PER_MM = 72.0 / 25.4

#: what the image export offers: (label, extension, PIL format, keeps alpha)
IMAGE_FORMATS = (
    ("PNG (lossless, transparent background)", ".png", "PNG", True),
    ("TIFF (lossless, for publishers)", ".tif", "TIFF", True),
    ("JPEG (small, no transparency)", ".jpg", "JPEG", False),
)


# --------------------------------------------------------------------------
# Ticks
# --------------------------------------------------------------------------
def nice_step(span: float, target: int = 5) -> float:
    """A round tick spacing giving roughly ``target`` intervals over ``span``.

    The 1/2/2.5/5/10 ladder, i.e. what any plotting library picks: ticks
    land on numbers a reader can hold in their head. Computed here rather
    than taken from pyqtgraph because the EPS is a different physical size
    from the panel on screen, so it wants its own tick density.

    The rung is the one *nearest* the ideal spacing on a log scale, not the
    first one above it: always rounding up turns a five-tick request into
    three ticks whenever the ideal spacing lands just past a rung.
    """
    span = abs(float(span))
    if span <= 0 or not np.isfinite(span):
        return 1.0
    raw = span / max(int(target), 1)
    magnitude = 10.0 ** math.floor(math.log10(raw))
    candidates = [multiple * magnitude for multiple in (1.0, 2.0, 2.5, 5.0, 10.0)]
    return min(candidates, key=lambda step: abs(math.log(step / raw)))


def tick_values(lo: float, hi: float, target: int = 5):
    """``(values, step)`` -- the round tick positions inside ``[lo, hi]``."""
    lo, hi = float(min(lo, hi)), float(max(lo, hi))
    step = nice_step(hi - lo, target)
    first = math.ceil(lo / step - 1e-9) * step
    n = int(math.floor((hi - first) / step + 1e-9)) + 1
    values = first + step * np.arange(max(n, 0))
    return values[(values >= lo - 1e-9) & (values <= hi + 1e-9)], step


def tick_format(step: float) -> str:
    """The shortest ``%``-format that still writes ``step`` exactly.

    Not simply ``-log10(step)``: the ladder in :func:`nice_step` includes
    2.5, and rounding 2.5 to no decimals turns a row of ticks into
    "-5 -2 0 2 5". So the number of decimals is the smallest that
    reproduces the step itself.
    """
    if step <= 0 or not np.isfinite(step):
        return "%g"
    if step >= 1e5 or step < 1e-4:
        return "%.1e"
    for decimals in range(7):
        if abs(round(step, decimals) - step) <= 1e-9 * step:
            return f"%.{decimals}f"
    return "%.6f"


# --------------------------------------------------------------------------
# The picture
# --------------------------------------------------------------------------
def _axis_to_index(axis: np.ndarray, coords: np.ndarray):
    """Fractional array index for each coordinate, and which ones are inside.

    "Inside" is the image's own extent -- the axis range grown by half a
    pixel at each end, since a pixel is centred on its axis value. This is
    the same convention the viewer draws with, so an export covers exactly
    the coloured area and not a half-pixel more.
    """
    axis = np.asarray(axis, dtype=float)
    n = axis.size
    if n == 0:
        return np.zeros_like(coords), np.zeros(coords.shape, dtype=bool)
    if n == 1:
        half = 0.5
        inside = np.abs(coords - axis[0]) <= half
        return np.zeros_like(coords), inside
    step = float(axis[1] - axis[0])
    index = (coords - float(axis[0])) / step if step else np.zeros_like(coords)
    inside = (index >= -0.5) & (index <= n - 0.5)
    return np.clip(index, 0.0, n - 1.0), inside


def _sample(array: np.ndarray, ix: np.ndarray, iy: np.ndarray, interpolate: bool):
    """Sample ``array`` (nx, ny) at fractional indices, NaN-aware.

    Bilinear where asked for, but with the weights renormalised over the
    finite corners only: a single missing pixel must not smear a hole four
    pixels wide across the export, and a point whose four neighbours are all
    missing stays missing.
    """
    nx, ny = array.shape
    if not interpolate:
        return array[np.rint(ix).astype(np.intp).clip(0, nx - 1)[:, None],
                     np.rint(iy).astype(np.intp).clip(0, ny - 1)[None, :]]

    x0 = np.floor(ix).astype(np.intp).clip(0, nx - 1)
    y0 = np.floor(iy).astype(np.intp).clip(0, ny - 1)
    x1 = np.minimum(x0 + 1, nx - 1)
    y1 = np.minimum(y0 + 1, ny - 1)
    tx = (ix - x0)[:, None]
    ty = (iy - y0)[None, :]

    out = np.zeros((ix.size, iy.size), dtype=float)
    weight = np.zeros_like(out)
    for xi, wx in ((x0, 1.0 - tx), (x1, tx)):
        for yi, wy in ((y0, 1.0 - ty), (y1, ty)):
            corner = array[xi[:, None], yi[None, :]]
            w = wx * wy
            good = np.isfinite(corner)
            out += np.where(good, corner * w, 0.0)
            weight += np.where(good, w, 0.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = np.where(weight > 1e-12, out / np.maximum(weight, 1e-12), np.nan)
    return out


def render_rgba(array_xy, x_axis, y_axis, *, x_range, y_range, lut,
                levels=None, gamma: float = 1.0, width: int = 1200,
                height: int = 900, interpolate: bool = True,
                invert_x: bool = False, invert_y: bool = False) -> np.ndarray:
    """Draw ``array_xy`` the way the panel draws it, at an exact pixel size.

    ``array_xy`` is (nx, ny) with ``x_axis`` horizontal; the result is an
    ``(height, width, 4)`` uint8 RGBA image covering exactly ``x_range`` x
    ``y_range``, row 0 at the top of the picture. Anything outside the data
    -- the margins an aspect lock leaves, and any NaN inside -- comes back
    fully transparent, so the export can be laid over a figure background
    without a white rectangle around it.

    ``levels`` is the colour window (the panel's Min/Max); ``gamma`` is
    applied inside that window exactly as on screen, so what is saved is
    what was being looked at.
    """
    array = np.asarray(array_xy, dtype=float)
    if array.ndim != 2:
        raise ValueError("render_rgba expects a 2D array")
    width, height = max(int(width), 1), max(int(height), 1)

    x0, x1 = float(min(x_range)), float(max(x_range))
    y0, y1 = float(min(y_range)), float(max(y_range))
    # pixel centres, not edges: a pixel covers a finite span and samples the
    # data at its middle
    xs = x0 + (np.arange(width) + 0.5) * (x1 - x0) / width
    ys = y0 + (np.arange(height) + 0.5) * (y1 - y0) / height
    if invert_x:
        xs = xs[::-1]
    # the top row of an image is the *last* row of a y-up plot
    if not invert_y:
        ys = ys[::-1]

    ix, in_x = _axis_to_index(x_axis, xs)
    iy, in_y = _axis_to_index(y_axis, ys)
    values = _sample(array, ix, iy, interpolate)          # (width, height)
    values = np.where(in_x[:, None] & in_y[None, :], values, np.nan)
    values = values.T                                      # (height, width)

    finite = np.isfinite(values)
    if levels is None:
        if finite.any():
            lo, hi = float(np.nanmin(values)), float(np.nanmax(values))
        else:
            lo, hi = 0.0, 1.0
    else:
        lo, hi = float(min(levels)), float(max(levels))
    if not np.isfinite(hi - lo) or hi <= lo:
        hi = lo + 1.0

    with np.errstate(invalid="ignore"):
        norm = np.clip((values - lo) / (hi - lo), 0.0, 1.0)
        if gamma and gamma != 1.0:
            norm = norm ** float(gamma)
    norm = np.nan_to_num(norm, nan=0.0)

    table = np.asarray(lut, dtype=np.uint8)
    if table.ndim != 2 or table.shape[1] < 3:
        raise ValueError("lut must be (N, 3) or (N, 4)")
    index = np.rint(norm * (table.shape[0] - 1)).astype(np.intp)
    rgba = np.zeros(values.shape + (4,), dtype=np.uint8)
    rgba[..., :3] = table[index][..., :3]
    rgba[..., 3] = np.where(finite, 255, 0)
    return rgba


def write_image(path: str, rgba: np.ndarray, *, pil_format: str = "PNG",
                keep_alpha: bool = True, dpi: float = None,
                background=(255, 255, 255)) -> str:
    """Write an RGBA array with PIL, flattening onto ``background`` for the
    formats that have no alpha (JPEG). ``dpi`` is only metadata -- the pixel
    count is what was asked for -- but it is what tells a layout program how
    big to place the picture."""
    from PIL import Image

    rgba = np.asarray(rgba, dtype=np.uint8)
    if keep_alpha and pil_format in ("PNG", "TIFF"):
        image = Image.fromarray(rgba, mode="RGBA")
    else:
        alpha = rgba[..., 3:4].astype(float) / 255.0
        flat = rgba[..., :3].astype(float) * alpha + np.asarray(
            background, dtype=float) * (1.0 - alpha)
        image = Image.fromarray(np.rint(flat).astype(np.uint8), mode="RGB")
    options = {}
    if dpi:
        options["dpi"] = (float(dpi), float(dpi))
    if pil_format == "JPEG":
        options["quality"] = 95
    image.save(path, format=pil_format, **options)
    return path


# --------------------------------------------------------------------------
# The axes, as vector PostScript
# --------------------------------------------------------------------------
_EPS_PROLOG = """/ctr { dup stringwidth pop 2 div neg 0 rmoveto show } bind def
/rgt { dup stringwidth pop neg 0 rmoveto show } bind def
/L { moveto lineto stroke } bind def
"""


def _label_width(text: str, font_size: float) -> float:
    """Helvetica is about 0.55 em per character averaged over digits and
    lower case -- close enough to reserve the right margin without parsing
    font metrics."""
    return 0.58 * font_size * max(len(text), 1)


def axes_eps(path: str, *, x_range, y_range, x_label: str = "", y_label: str = "",
             width_mm: float = 80.0, height_mm: float = 60.0,
             font_size: float = 9.0, line_width: float = 0.8,
             tick_length: float = 4.0, minor_ticks: int = 1,
             ticks_inward: bool = True, invert_x: bool = False,
             invert_y: bool = False, x_targets: int = 5, y_targets: int = 5,
             title: str = "") -> dict:
    """Write the axis frame alone as an EPS, and report where the box is.

    Nothing is drawn inside the frame: the point is to place the raster
    export (rendered over the same ``x_range``/``y_range``) in the box and
    keep the annotation editable. The returned dict gives the box in
    PostScript points relative to the file's lower-left corner --
    ``{"x", "y", "width", "height"}`` -- which is what the caller tells the
    user to type into Illustrator's transform palette.

    Ticks point inward by default, which is the convention in ARPES figures
    and keeps the box free of anything that would collide with a
    neighbouring panel.
    """
    x0, x1 = float(min(x_range)), float(max(x_range))
    y0, y1 = float(min(y_range)), float(max(y_range))
    box_w = max(float(width_mm), 5.0) * PT_PER_MM
    box_h = max(float(height_mm), 5.0) * PT_PER_MM

    x_ticks, x_step = tick_values(x0, x1, x_targets)
    y_ticks, y_step = tick_values(y0, y1, y_targets)
    x_fmt, y_fmt = tick_format(x_step), tick_format(y_step)
    x_texts = [x_fmt % v for v in x_ticks]
    y_texts = [y_fmt % v for v in y_ticks]

    gap = 0.4 * font_size
    outward = 0.0 if ticks_inward else float(tick_length)
    left = (outward + gap + max([_label_width(t, font_size) for t in y_texts] or [0])
            + (1.6 * font_size if y_label else 0.4 * font_size))
    bottom = (outward + gap + 1.0 * font_size
              + (1.5 * font_size if x_label else 0.3 * font_size))
    top = (1.8 * font_size if title else 0.5 * font_size) + outward
    right = 0.5 * font_size + max(_label_width(x_texts[-1], font_size) / 2
                                  if x_texts else 0.0, outward)

    total_w, total_h = left + box_w + right, bottom + box_h + top

    def px(value):
        frac = (value - x0) / (x1 - x0) if x1 > x0 else 0.0
        if invert_x:
            frac = 1.0 - frac
        return left + frac * box_w

    def py(value):
        frac = (value - y0) / (y1 - y0) if y1 > y0 else 0.0
        if invert_y:
            frac = 1.0 - frac
        return bottom + frac * box_h

    out = ["%!PS-Adobe-3.0 EPSF-3.0",
           f"%%BoundingBox: 0 0 {math.ceil(total_w)} {math.ceil(total_h)}",
           f"%%HiResBoundingBox: 0 0 {total_w:.3f} {total_h:.3f}",
           "%%Creator: ARPES viewer (axes only -- place the image export in the box)",
           f"%%CreationDate: {time.strftime('%Y-%m-%d %H:%M:%S')}",
           "%%DocumentData: Clean7Bit",
           "%%EndComments",
           f"% data box: x {left:.3f} y {bottom:.3f} w {box_w:.3f} h {box_h:.3f} pt",
           f"% x range: {x0:.9g} .. {x1:.9g}",
           f"% y range: {y0:.9g} .. {y1:.9g}",
           _EPS_PROLOG,
           f"{line_width:.3f} setlinewidth 0 setgray 1 setlinecap 1 setlinejoin",
           f"/Helvetica findfont {font_size:.3f} scalefont setfont"]

    # the frame
    out.append(f"newpath {left:.3f} {bottom:.3f} moveto "
               f"{box_w:.3f} 0 rlineto 0 {box_h:.3f} rlineto "
               f"{-box_w:.3f} 0 rlineto closepath stroke")

    sign = -1.0 if ticks_inward else 1.0

    def ticks(values, major: bool):
        length = tick_length * (1.0 if major else 0.55)
        for value in values:
            if value < x0 - 1e-9 or value > x1 + 1e-9:
                continue
            x = px(value)
            out.append(f"{x:.3f} {bottom:.3f} {x:.3f} {bottom - sign * length:.3f} L")
            out.append(f"{x:.3f} {bottom + box_h:.3f} "
                       f"{x:.3f} {bottom + box_h + sign * length:.3f} L")

    def yticks(values, major: bool):
        length = tick_length * (1.0 if major else 0.55)
        for value in values:
            if value < y0 - 1e-9 or value > y1 + 1e-9:
                continue
            y = py(value)
            out.append(f"{left:.3f} {y:.3f} {left - sign * length:.3f} {y:.3f} L")
            out.append(f"{left + box_w:.3f} {y:.3f} "
                       f"{left + box_w + sign * length:.3f} {y:.3f} L")

    ticks(x_ticks, True)
    yticks(y_ticks, True)
    if minor_ticks > 0:
        sub = np.concatenate([x_ticks + x_step * (i + 1) / (minor_ticks + 1)
                              for i in range(minor_ticks)] +
                             [x_ticks - x_step * (i + 1) / (minor_ticks + 1)
                              for i in range(minor_ticks)]) if x_ticks.size else x_ticks
        ticks(sub, False)
        sub = np.concatenate([y_ticks + y_step * (i + 1) / (minor_ticks + 1)
                              for i in range(minor_ticks)] +
                             [y_ticks - y_step * (i + 1) / (minor_ticks + 1)
                              for i in range(minor_ticks)]) if y_ticks.size else y_ticks
        yticks(sub, False)

    # tick labels
    baseline = bottom - outward - gap - font_size * 0.78
    for value, text in zip(x_ticks, x_texts):
        out.append(f"{px(value):.3f} {baseline:.3f} moveto ({_ps(text)}) ctr")
    for value, text in zip(y_ticks, y_texts):
        out.append(f"{left - outward - gap:.3f} {py(value) - font_size * 0.34:.3f} "
                   f"moveto ({_ps(text)}) rgt")

    # axis titles
    if x_label:
        out.append(f"{left + box_w / 2:.3f} {baseline - font_size * 1.25:.3f} "
                   f"moveto ({_ps(x_label)}) ctr")
    if y_label:
        out.append("gsave")
        out.append(f"{left - outward - gap - max([_label_width(t, font_size) for t in y_texts] or [0]) - font_size * 0.6:.3f} "
                   f"{bottom + box_h / 2:.3f} translate 90 rotate")
        out.append(f"0 0 moveto ({_ps(y_label)}) ctr")
        out.append("grestore")
    if title:
        out.append(f"{left + box_w / 2:.3f} {bottom + box_h + outward + font_size * 0.6:.3f} "
                   f"moveto ({_ps(title)}) ctr")

    out.append("showpage")
    out.append("%%EOF")
    with open(path, "w", encoding="latin-1", errors="replace") as fh:
        fh.write("\n".join(out) + "\n")
    return {"x": left, "y": bottom, "width": box_w, "height": box_h,
            "page_width": total_w, "page_height": total_h}


def _ps(text: str) -> str:
    """Escape a string for a PostScript literal. Non-Latin-1 characters --
    the degree sign survives, but an angstrom or a minus sign may not -- are
    replaced rather than written raw, because a broken EPS opens as nothing
    at all."""
    text = (str(text).replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)"))
    try:
        text.encode("latin-1")
    except UnicodeEncodeError:
        replacements = {"Å": "A", "−": "-", "°": "\\260",
                        "⁻": "-", "¹": "1", "²": "2",
                        "—": "-", "–": "-", "·": "."}
        text = "".join(replacements.get(ch, ch if ord(ch) < 256 else "?")
                       for ch in text)
    return text
