"""
tools/processing/detector.py
============================
Corrections made on the detector's own pixels, before anything else: a
trapezoidal (keystone) image made rectangular again, and dead or bad
regions patched from their neighbours.

Ported from the lab's Quick_Panel functions ``StretchDitherImage.m`` (the
keystone) and ``InterpRegionAlongX.m`` / ``InterpRegionAlongY.m`` (the
patch), with ``inpaint_nans.m``'s idea for a smoother patch, and these
changes:

* **The keystone as an exact inverse map.** The .m file pushed every pixel
  to its stretched position and re-gridded the scattered cloud with
  ``TriScatteredInterp`` (a Delaunay triangulation of every pixel of the
  image, for every image). Here each output pixel's source position is
  written down directly -- along each row the trapezoid's edges at that
  row are mapped onto the rectangle's -- and read off the regular grid by
  linear interpolation: exact, a single vectorised pass, and no
  triangulation artefacts at the edges.
* **Counts are kept.** Stretching a row spreads its counts over more
  pixels; the .m file left the values as they were, so a stretched row came
  out brighter in total than it went in. Here each row is multiplied by its
  stretch factor's inverse (optional), so the total is conserved.
* **The patch: from both sides, or smooth.** ``InterpRegionAlongX.m``
  blends the column just left of the region with the one just right of it
  (and ``...AlongY.m`` the rows above and below). That is kept ("along
  x" / "along y"), and a third way added: the harmonic fill, which solves
  Laplace's equation inside the region with the surrounding pixels as its
  boundary -- the smoothest patch that meets the data on every side, the
  idea of ``inpaint_nans.m``. Its sparse factorisation is made once and
  reused for every slice of a map.
* **Bad pixels can be found.** On a map, a dead or hot detector pixel is
  dead or hot in every slice while the photoemission moves; the mean
  detector image over the scan, divided by its own median filter, shows
  them, and :func:`find_bad_regions` returns them as regions to patch.

Applied to a map, everything is done on each slice's detector image --
(slit angle, energy) for every scan step -- with the same geometry.

Qt-free, so it can be tested on its own.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage, sparse
from scipy.sparse.linalg import splu

__all__ = ["Keystone", "keystone_correct", "Region", "patch_regions",
           "find_bad_regions", "detector_frames", "apply_to", "not_on_detector"]


def not_on_detector(info: dict):
    """Why this dataset is no longer on the detector's pixels, or None.

    The de-grid tool's test (:func:`tools.processing.degrid.not_pixel_locked`)
    except that de-gridding itself is allowed before these corrections: it
    does not move a pixel.
    """
    from tools.processing.degrid import not_pixel_locked

    kept = {}
    for key, value in (info or {}).items():
        if str(key).startswith("degrid."):
            continue
        if (str(key).startswith("proc.step.") and not str(key).endswith(".from")
                and str(value).split("(")[0].strip() == "degrid"):
            continue
        kept[key] = value
    return not_pixel_locked(kept)


def _index(axis, values):
    """Fractional indices of ``values`` on the (monotonic) ``axis``."""
    axis = np.asarray(axis, dtype=float)
    order = np.argsort(axis)
    return np.interp(values, axis[order], order.astype(float),
                     left=np.nan, right=np.nan)


# --------------------------------------------------------------------------
# Keystone
# --------------------------------------------------------------------------
class Keystone:
    """A trapezoid on the detector image, to be made the rectangle of its
    first edge.

    The image is ``[x, y]``; ``along`` (0 or 1) is the axis the image is
    stretched along -- the axis whose acceptance changes, normally the
    angle (0). The trapezoid's two edges along that axis run from
    ``(start_a, end_a)`` at ``y = level_a`` to ``(start_b, end_b)`` at
    ``y = level_b`` (data units; ``y`` meaning the other axis). After the
    correction every line between the two levels spans ``start_a ..
    end_a``, the extent at ``level_a``.
    """

    def __init__(self, start_a, end_a, level_a, start_b, end_b, level_b,
                 along: int = 0):
        self.start_a, self.end_a, self.level_a = float(start_a), float(end_a), float(level_a)
        self.start_b, self.end_b, self.level_b = float(start_b), float(end_b), float(level_b)
        self.along = int(along)
        if self.level_a == self.level_b:
            raise ValueError("the trapezoid's two edges are at the same level")
        if self.end_a == self.start_a:
            raise ValueError("the reference edge has no width")

    @classmethod
    def from_corners(cls, corners, along: int = 0):
        """From four (x, y) corners in any order: the two with the larger
        "level" coordinate form one edge, the other two the other."""
        corners = np.asarray(corners, dtype=float).reshape(4, 2)
        level = 1 - along
        order = np.argsort(corners[:, level])
        low, high = corners[order[:2]], corners[order[2:]]
        low = low[np.argsort(low[:, along])]
        high = high[np.argsort(high[:, along])]
        return cls(high[0, along], high[1, along], high[:, level].mean(),
                   low[0, along], low[1, along], low[:, level].mean(), along)

    def edges_at(self, level):
        """The trapezoid's two edges (start, end) at ``level``."""
        f = (np.asarray(level, dtype=float) - self.level_a) / (self.level_b - self.level_a)
        return (self.start_a + f * (self.start_b - self.start_a),
                self.end_a + f * (self.end_b - self.end_a))

    def as_dict(self):
        return {k: getattr(self, k) for k in ("start_a", "end_a", "level_a",
                                              "start_b", "end_b", "level_b",
                                              "along")}


def keystone_correct(image, x_axis, y_axis, keystone: Keystone, *,
                     conserve: bool = True):
    """The image with ``keystone``'s trapezoid stretched to a rectangle.

    For every output point ``(u, level)`` the source point along the
    stretched axis is ``s(level) + (u - start_a) * w(level) / w_a`` with
    ``s``, ``w`` the trapezoid's start and width at that level; read off
    the input by linear interpolation. Points whose source falls outside
    the image are NaN. ``conserve`` scales each line by ``w(level) / w_a``
    so that its total count is unchanged.
    """
    image = np.asarray(image, dtype=float)
    axes = [np.asarray(x_axis, dtype=float), np.asarray(y_axis, dtype=float)]
    along, level_dim = keystone.along, 1 - keystone.along
    u = axes[along]
    levels = axes[level_dim]
    start, end = keystone.edges_at(levels)
    width = end - start
    ratio = width / (keystone.end_a - keystone.start_a)
    uu = u[:, None] if along == 0 else u[None, :]
    ss = start[None, :] if along == 0 else start[:, None]
    rr = ratio[None, :] if along == 0 else ratio[:, None]
    source = ss + (uu - keystone.start_a) * rr
    source_index = _index(u, source)
    grid = np.indices(image.shape).astype(float)
    grid[along] = source_index
    out = ndimage.map_coordinates(np.nan_to_num(image), [np.nan_to_num(g, nan=-1e9)
                                                         for g in grid],
                                  order=1, mode="constant", cval=np.nan)
    out[np.isnan(source_index)] = np.nan
    if conserve:
        out = out * rr
    return out


# --------------------------------------------------------------------------
# Bad regions
# --------------------------------------------------------------------------
class Region:
    """A rectangle of the detector image, in data units:
    ``x0..x1`` by ``y0..y1``."""

    def __init__(self, x0, x1, y0, y1):
        self.x0, self.x1 = sorted((float(x0), float(x1)))
        self.y0, self.y1 = sorted((float(y0), float(y1)))

    def mask(self, x_axis, y_axis):
        x = np.asarray(x_axis, dtype=float)
        y = np.asarray(y_axis, dtype=float)
        return (((x >= self.x0) & (x <= self.x1))[:, None]
                & ((y >= self.y0) & (y <= self.y1))[None, :])

    def as_tuple(self):
        return (self.x0, self.x1, self.y0, self.y1)


def _linear_fill(image, mask, dim: int):
    """Blend the pixels just outside the masked runs, along ``dim``."""
    out = image.copy()
    moved = np.moveaxis(out, dim, 0)
    holes = np.moveaxis(mask, dim, 0)
    n = moved.shape[0]
    index = np.arange(n, dtype=float)
    for col in range(moved.shape[1]):
        hole = holes[:, col]
        if not hole.any():
            continue
        good = ~hole & np.isfinite(moved[:, col])
        if good.sum() == 0:
            continue
        moved[hole, col] = np.interp(index[hole], index[good], moved[good, col])
    return out


class _Harmonic:
    """The Laplace fill for one mask, factorised once."""

    def __init__(self, mask):
        self.mask = np.asarray(mask, bool)
        self.cells = np.flatnonzero(self.mask)
        shape = self.mask.shape
        number = -np.ones(self.mask.size, dtype=int)
        number[self.cells] = np.arange(self.cells.size)
        rows, cols, vals = [], [], []
        self.boundary = []          # (equation, flat index of a known pixel)
        ii, jj = np.unravel_index(self.cells, shape)
        for eq, (i, j) in enumerate(zip(ii, jj)):
            degree = 0
            for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                a, b = i + di, j + dj
                if not (0 <= a < shape[0] and 0 <= b < shape[1]):
                    continue
                degree += 1
                flat = a * shape[1] + b
                if self.mask.flat[flat]:
                    rows.append(eq)
                    cols.append(number[flat])
                    vals.append(-1.0)
                else:
                    self.boundary.append((eq, flat))
            rows.append(eq)
            cols.append(eq)
            vals.append(float(degree))
        matrix = sparse.csc_matrix((vals, (rows, cols)),
                                   shape=(self.cells.size, self.cells.size))
        self.solver = splu(matrix)
        self.boundary = np.array(self.boundary, dtype=int).reshape(-1, 2)

    def fill(self, image):
        out = np.array(image, dtype=float)
        flat = np.nan_to_num(out).ravel()
        rhs = np.zeros(self.cells.size)
        if self.boundary.size:
            np.add.at(rhs, self.boundary[:, 0], flat[self.boundary[:, 1]])
        out.ravel()[self.cells] = self.solver.solve(rhs)
        return out


METHODS = ("linear along x", "linear along y", "smooth (harmonic)")


def patch_regions(image, x_axis, y_axis, regions, method: str = METHODS[2],
                  _cache=None):
    """``image`` with every region replaced from its surroundings.

    ``method``: "linear along x" (the pixels just outside, before and after
    the region along x, blended -- ``InterpRegionAlongX.m``), "linear along
    y", or "smooth (harmonic)" (Laplace's equation inside, the surrounding
    pixels as its boundary). ``_cache`` (a dict) keeps the harmonic
    factorisation between calls with the same mask.
    """
    image = np.asarray(image, dtype=float)
    mask = np.zeros(image.shape, dtype=bool)
    for region in regions:
        mask |= region.mask(x_axis, y_axis)
    if not mask.any():
        return image.copy()
    if mask.all():
        raise ValueError("the regions cover the whole image: nothing is left "
                         "to patch them from")
    if method == METHODS[0]:
        return _linear_fill(image, mask, 0)
    if method == METHODS[1]:
        return _linear_fill(image, mask, 1)
    key = mask.tobytes()
    if _cache is not None and _cache.get("key") == key:
        solver = _cache["solver"]
    else:
        solver = _Harmonic(mask)
        if _cache is not None:
            _cache.update(key=key, solver=solver)
    return solver.fill(image)


def find_bad_regions(frame, x_axis, y_axis, *, threshold: float = 0.3,
                     size: int = 5, grow: int = 1, max_regions: int = 50,
                     sigmas: float = 6.0):
    """Pixels that stand out of their neighbourhood, as regions.

    ``frame`` is a mean detector image -- of a map, averaged over the scan,
    so that the photoemission (which moves) is smoothed out and a pixel's
    own fault (which does not) is left. A pixel is bad when it differs from
    the median of its ``size`` x ``size`` neighbourhood by more than
    ``threshold`` (a fraction of that median) and by more than ``sigmas``
    times the image's robust pixel-to-pixel scatter: dead pixels below, hot
    ones above. Neighbouring bad pixels are joined, grown by ``grow`` pixels,
    and returned as bounding rectangles, worst first.
    """
    frame = np.asarray(frame, dtype=float)
    local = ndimage.median_filter(np.nan_to_num(frame), size=size)
    floor = max(1e-12, 0.02 * float(np.nanmax(local)))
    relative = (frame - local) / np.maximum(local, floor)
    bad = np.abs(np.nan_to_num(relative)) > threshold
    bad &= local > floor                      # not the unlit edge
    # ...and well outside the image's own pixel-to-pixel scatter, so that
    # noise (which a short acquisition has plenty of) is not mistaken for a
    # faulty pixel.
    residual = np.nan_to_num(frame - local)
    spread = 1.4826 * float(np.median(np.abs(residual - np.median(residual))))
    if spread > 0:
        bad &= np.abs(residual) > sigmas * spread
    if grow:
        bad = ndimage.binary_dilation(bad, iterations=int(grow))
    labels, count = ndimage.label(bad)
    x = np.asarray(x_axis, dtype=float)
    y = np.asarray(y_axis, dtype=float)
    dx = float(np.median(np.abs(np.diff(x)))) / 2 if x.size > 1 else 0.0
    dy = float(np.median(np.abs(np.diff(y)))) / 2 if y.size > 1 else 0.0
    found = []
    for index, box in enumerate(ndimage.find_objects(labels), start=1):
        if box is None:
            continue
        strength = float(np.max(np.abs(relative[box][labels[box] == index])))
        sx, sy = box
        found.append((strength, Region(x[sx.start] - dx, x[sx.stop - 1] + dx,
                                       y[sy.start] - dy, y[sy.stop - 1] + dy)))
    found.sort(key=lambda item: -item[0])
    return [region for _s, region in found[:max_regions]]


# --------------------------------------------------------------------------
# Whole datasets
# --------------------------------------------------------------------------
def detector_frames(kind: str):
    """Which axes of a ``kind`` form the detector image: (0, 1) for a cut,
    (1, 2) -- slit angle and energy -- for a map or kz map."""
    if kind == "cut":
        return (0, 1)
    if kind in ("map", "kz_map"):
        return (1, 2)
    raise ValueError(f"a {kind} is not on the detector's pixels")


def apply_to(values, axes, kind: str, *, keystone: Keystone = None,
             regions=(), method: str = METHODS[2], conserve: bool = True,
             progress=None):
    """Both corrections on a cut, or on every slice of a map: the patch
    first (on the pixels as they are), then the keystone."""
    values = np.asarray(values, dtype=float)
    fx, fy = detector_frames(kind)
    x_axis, y_axis = axes[fx], axes[fy]
    cache = {}

    def one(frame):
        if regions:
            frame = patch_regions(frame, x_axis, y_axis, regions, method, cache)
        if keystone is not None:
            frame = keystone_correct(frame, x_axis, y_axis, keystone,
                                     conserve=conserve)
        return frame

    if kind == "cut":
        return one(values)
    out = np.empty_like(values)
    n = values.shape[0]
    for i in range(n):
        out[i] = one(values[i])
        if progress and (i % 8 == 0 or i == n - 1):
            progress((i + 1) / n, f"Slice {i + 1} of {n}")
    return out
