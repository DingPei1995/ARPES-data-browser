"""
tools/volume.py
=============
Three-dimensional display geometry, with no Qt and no OpenGL: orthogonal
slices, the notched cube of ``volume_3d_plot_w_notch.m``, a shaded
isosurface and volume projections.

Everything a 3-D view needs is produced here as plain arrays -- a list of
:class:`Face` objects, each an image plus the four corners of the
parallelogram it is drawn into, or a single projected image. The painter
that draws them is the ordinary 2-D one, which is the point: the same
figure pipeline that exports a cut as PDF exports a notched cube as PDF,
at whatever resolution the journal asks for, instead of leaving it trapped
in a screen-sized OpenGL bitmap.

Why not OpenGL
--------------
Every surface here is planar, and the parallel projection of a rectangle is
a parallelogram, which is exactly an affine image transform. So a painter's
algorithm over depth-sorted faces gives a pixel-accurate result with no
graphics driver, no PyOpenGL dependency, and -- unlike a screen grab of a
3-D canvas -- a vector-exportable one.

Two repairs to the MATLAB original
----------------------------------
``volume_3d_plot_w_notch.m`` writes its eight notch orientations out
longhand, about five hundred lines of near-identical code, and computes the
face sizes as ``round((xe - xb) / xstep)``, which is zero or negative when
the notch vertex reaches an edge -- an empty ``linspace`` and an error from
inside ``slice``. Here the orientation is three sign flips over one
construction (:func:`notched_box`), and a face that has been squeezed to
nothing is simply not produced.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import ndimage


# ==========================================================================
# Faces
# ==========================================================================
@dataclass
class Face:
    """One flat, textured quadrilateral of a 3-D scene.

    ``corners`` are the four world-space corners in the order
    ``(origin, +u, +u+v, +v)``, and ``image`` is the data sampled on a
    regular ``(nu, nv)`` grid over that rectangle. The painter needs only
    an affine map from image space to the projected corners.
    """
    image: np.ndarray
    corners: np.ndarray                  # (4, 3)
    name: str = ""
    axis: int = 0                        # which axis the plane is normal to
    position: float = 0.0                # where along it
    outward: int = 0                     # +1/-1 along `axis`; 0 = two-sided

    def centroid(self) -> np.ndarray:
        return np.asarray(self.corners, dtype=float).mean(axis=0)

    def normal(self) -> np.ndarray:
        """The outward unit normal, or the geometric one for a two-sided face.

        Which side of a face is outside is a fact about how the solid was
        built, not something to be recovered afterwards from the winding of
        the corners or from where the face sits relative to the box centre.
        Both of those guesses fail on the inner walls of a notch, so the
        builder records it and the renderer reads it.
        """
        n = np.zeros(3)
        if self.outward:
            n[self.axis] = float(np.sign(self.outward))
            return n
        c = np.asarray(self.corners, dtype=float)
        v = np.cross(c[1] - c[0], c[3] - c[0])
        length = float(np.linalg.norm(v))
        return v / length if length > 0 else v


def _index_coords(axis, values):
    """Fractional index of each value on a regular axis."""
    axis = np.asarray(axis, dtype=float)
    start = float(axis[0])
    step = (float(axis[-1]) - start) / max(axis.size - 1, 1)
    if step == 0:
        step = 1.0
    return (np.asarray(values, dtype=float) - start) / step


def sample_plane(values, axes, axis: int, position: float, u_range, v_range,
                 nu: int, nv: int) -> np.ndarray:
    """Sample the volume on a rectangle of the plane ``axis = position``.

    The two remaining axes, in their natural order, become u and v.
    Anything outside the measured volume comes back NaN rather than as the
    edge value: a face that reaches past the data should show that it does.
    """
    values = np.asarray(values, dtype=float)
    axis = int(axis) % 3
    others = [d for d in range(3) if d != axis]
    u = np.linspace(float(u_range[0]), float(u_range[1]), int(nu))
    v = np.linspace(float(v_range[0]), float(v_range[1]), int(nv))
    gu, gv = np.meshgrid(u, v, indexing="ij")

    coords = [None, None, None]
    coords[axis] = np.full(gu.shape, _index_coords(axes[axis], position))
    coords[others[0]] = _index_coords(axes[others[0]], gu)
    coords[others[1]] = _index_coords(axes[others[1]], gv)
    stacked = np.array([c.ravel() for c in coords])

    mask = np.isfinite(values)
    filled = np.where(mask, values, 0.0)
    num = ndimage.map_coordinates(filled, stacked, order=1, mode="constant", cval=0.0)
    den = ndimage.map_coordinates(mask.astype(float), stacked, order=1,
                                  mode="constant", cval=0.0)
    inside = np.ones(stacked.shape[1], dtype=bool)
    for d in range(3):
        inside &= (stacked[d] >= 0) & (stacked[d] <= values.shape[d] - 1)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = np.where((den > 0.5) & inside, num / np.maximum(den, 1e-12), np.nan)
    return out.reshape(gu.shape)


def _face_corners(axis: int, position: float, u_range, v_range) -> np.ndarray:
    """The four world corners of a plane rectangle, in (0, +u, +u+v, +v)."""
    axis = int(axis) % 3
    others = [d for d in range(3) if d != axis]
    corners = np.zeros((4, 3), dtype=float)
    corners[:, axis] = float(position)
    us = [u_range[0], u_range[1], u_range[1], u_range[0]]
    vs = [v_range[0], v_range[0], v_range[1], v_range[1]]
    corners[:, others[0]] = us
    corners[:, others[1]] = vs
    return corners


def _span(a, b):
    return (min(float(a), float(b)), max(float(a), float(b)))


def _resolution(span, limits, density: int) -> int:
    """How many samples a face of this length deserves, capped so that a
    thin sliver of a notch does not ask for a megapixel."""
    total = abs(float(limits[1]) - float(limits[0])) or 1.0
    fraction = abs(span[1] - span[0]) / total
    return int(np.clip(round(density * fraction), 2, density))


# ==========================================================================
# Orthogonal slices
# ==========================================================================
def _split(interval, cut):
    """One interval as one or two pieces, cut at ``cut`` if it falls inside."""
    lo, hi = float(interval[0]), float(interval[1])
    if cut is None or not (lo < float(cut) < hi):
        return [(lo, hi)]
    return [(lo, float(cut)), (float(cut), hi)]


def slice_faces(values, axes, positions, *, density: int = 400,
                limits=None) -> list:
    """The three orthogonal cut planes of ``slice_3d_plot.m``.

    ``positions`` is ``(x, y, z)``; any entry that is None is skipped, so
    one, two or three planes can be shown. ``limits`` optionally crops the
    box each plane is drawn over.

    Each plane is **cut where the others cross it**, into up to four
    rectangles. That is not cosmetic. Three orthogonal planes all pass
    through one point, so no ordering of them as whole planes is correct:
    a painter's algorithm has to draw one entirely in front of another, and
    the picture then shows a vertical plane floating on top of a horizontal
    one it should be half behind. Split at the intersections, the pieces are
    disjoint, and sorting by depth is exact.
    """
    values = np.asarray(values, dtype=float)
    bounds = _bounds(axes, limits)
    faces = []
    for axis, position in enumerate(positions):
        if position is None:
            continue
        others = [d for d in range(3) if d != axis]
        for u_range in _split(bounds[others[0]], positions[others[0]]):
            for v_range in _split(bounds[others[1]], positions[others[1]]):
                nu = _resolution(u_range, bounds[others[0]], density)
                nv = _resolution(v_range, bounds[others[1]], density)
                image = sample_plane(values, axes, axis, position,
                                     u_range, v_range, nu, nv)
                faces.append(Face(
                    image, _face_corners(axis, position, u_range, v_range),
                    name="xyz"[axis] + f" = {float(position):.4g}",
                    axis=axis, position=float(position)))
    return faces


def _bounds(axes, limits=None):
    out = []
    for d in range(3):
        axis = np.asarray(axes[d], dtype=float)
        lo, hi = float(axis.min()), float(axis.max())
        if limits is not None and limits[d] is not None:
            want = limits[d]
            if want[0] is not None:
                lo = max(lo, float(want[0]))
            if want[1] is not None:
                hi = min(hi, float(want[1]))
        if hi <= lo:
            lo, hi = float(axis.min()), float(axis.max())
        out.append((lo, hi))
    return out


# ==========================================================================
# The notched cube
# ==========================================================================
def notched_box(values, axes, notch, *, corner=(1, 1, 1), density: int = 400,
                limits=None) -> list:
    """The cube with one corner cut away -- the best single picture of a
    3-D ARPES volume, because it shows a constant-energy map and two
    dispersions at once.

    ``notch`` is the ``(x, y, z)`` position of the inner vertex; ``corner``
    is three signs saying which corner is removed (``+1`` means the high
    end of that axis). That replaces the MATLAB version's eight hand-written
    cases and its letter-labelled diagram: the corner is now just the corner
    the user clicked, and any of the eight follows from the signs.

    Twelve faces come out -- three whole outer faces, three L-shaped outer
    faces cut into two rectangles each, and the three inner faces of the
    notch -- minus any that the notch has squeezed to nothing.
    """
    values = np.asarray(values, dtype=float)
    bounds = _bounds(axes, limits)
    signs = [1 if s >= 0 else -1 for s in corner]

    near = [bounds[d][1] if signs[d] > 0 else bounds[d][0] for d in range(3)]
    far = [bounds[d][0] if signs[d] > 0 else bounds[d][1] for d in range(3)]
    vertex = [float(np.clip(notch[d], min(bounds[d]), max(bounds[d])))
              for d in range(3)]

    faces = []

    def add(axis, position, u_range, v_range, name, outward):
        u_range = _span(*u_range)
        v_range = _span(*v_range)
        if (u_range[1] - u_range[0]) <= 0 or (v_range[1] - v_range[0]) <= 0:
            return                       # the notch reaches the edge here
        nu = _resolution(u_range, bounds[[d for d in range(3) if d != axis][0]],
                         density)
        nv = _resolution(v_range, bounds[[d for d in range(3) if d != axis][1]],
                         density)
        image = sample_plane(values, axes, axis, position, u_range, v_range, nu, nv)
        faces.append(Face(image, _face_corners(axis, position, u_range, v_range),
                          name=name, axis=axis, position=float(position),
                          outward=outward))

    for axis in range(3):
        b, c = [d for d in range(3) if d != axis]
        letter = "xyz"[axis]
        s = signs[axis]

        # 1. The whole face at the far end of this axis: outside is further
        #    from the notched corner.
        add(axis, far[axis], bounds[b], bounds[c], f"{letter} outer", -s)

        # 2. The face at the notched end: an L, cut on the full quadrant
        #    grid rather than into the two rectangles it minimally needs.
        #    Two rectangles would overlap only *partially* along their
        #    common line -- one runs the whole width, the other stops at the
        #    notch -- and a seam like that cannot be told from a real edge
        #    by comparing edges, so it gets drawn as a line across the face
        #    that exists nowhere in the geometry. On the quadrant grid every
        #    shared edge coincides exactly, which :func:`outline_edges` can
        #    recognise and drop. The cost is one extra rectangle per face.
        for u_range in _split(bounds[b], vertex[b]):
            for v_range in _split(bounds[c], vertex[c]):
                inside_notch = (
                    min(abs(u - near[b]) for u in u_range) < 1e-12
                    and min(abs(v - near[c]) for v in v_range) < 1e-12)
                if inside_notch:
                    continue          # this quadrant is the removed corner
                add(axis, near[axis], u_range, v_range, f"{letter} near", s)

        # 3. The inner wall of the notch. The material sits on the far side
        #    of it -- the corner block in front has been removed -- so it
        #    faces the same way as the near face, towards the viewer who is
        #    looking into the notch.
        add(axis, vertex[axis], (vertex[b], near[b]), (vertex[c], near[c]),
            f"{letter} notch", s)
    return faces


# ==========================================================================
# Projection
# ==========================================================================
@dataclass
class Camera:
    """A parallel (orthographic) camera.

    Parallel rather than perspective on purpose: a figure is measured off
    the page, and a perspective view makes two equal momentum intervals
    print at different lengths.

    ``aspect`` scales the three axes before projection, which is what makes
    a k-k-E volume presentable -- ``axis equal`` is exactly wrong there,
    since an inverse angstrom and an electronvolt have no common length.
    """
    azimuth: float = 45.0                # degrees, around the z axis
    elevation: float = 25.0              # degrees, above the xy plane
    roll: float = 0.0
    aspect: tuple = (1.0, 1.0, 1.0)

    def basis(self) -> np.ndarray:
        """Rows: right, up, towards the viewer."""
        az = np.radians(self.azimuth)
        el = np.radians(self.elevation)
        view = np.array([np.cos(el) * np.cos(az), np.cos(el) * np.sin(az),
                         np.sin(el)])
        right = np.array([-np.sin(az), np.cos(az), 0.0])
        up = np.cross(view, right)
        if self.roll:
            r = np.radians(self.roll)
            right, up = (np.cos(r) * right + np.sin(r) * up,
                         -np.sin(r) * right + np.cos(r) * up)
        return np.array([right, up, view])

    def normalise(self, points, bounds):
        """World points into the unit cube, with ``aspect`` applied."""
        points = np.atleast_2d(np.asarray(points, dtype=float))
        out = np.empty_like(points)
        for d in range(3):
            lo, hi = bounds[d]
            span = (hi - lo) or 1.0
            out[:, d] = ((points[:, d] - lo) / span - 0.5) * float(self.aspect[d])
        return out

    def project(self, points, bounds):
        """``(xy, depth)``: screen coordinates, and distance towards the
        viewer for the painter's ordering."""
        local = self.normalise(points, bounds)
        basis = self.basis()
        screen = local @ basis[:2].T
        depth = local @ basis[2]
        return screen, depth


def _edge_key(a, b, places: int = 9):
    """A corner pair as a hashable key, independent of which end is first."""
    ends = tuple(sorted((tuple(round(float(v), places) for v in a),
                         tuple(round(float(v), places) for v in b))))
    return ends


def outline_edges(faces):
    """Which of each face's four edges are edges of the **solid**.

    Returns one list of four booleans per face, in the corner order
    ``(0-1, 1-2, 2-3, 3-0)``.

    Faces here are split for two reasons that have nothing to do with the
    shape being drawn: an L-shaped outer face is cut into two rectangles
    because a rectangle is what can be texture-mapped, and a slice plane is
    cut into quadrants because that is what makes depth sorting exact.
    Outlining each rectangle then draws the cuts as if they were edges --
    a line across the top of the cube and another down its side, neither of
    which is anywhere in the geometry.

    An edge shared by two faces in the same plane is such a cut; an edge
    that only one face has is a real boundary. That single test removes
    every seam and keeps every genuine edge, including the concave ones
    around a notch, where the two faces meet at a right angle rather than
    in one plane.
    """
    counts = {}
    for face in faces:
        corners = np.asarray(face.corners, dtype=float)
        for i in range(4):
            key = (face.axis, round(float(face.position), 9),
                   _edge_key(corners[i], corners[(i + 1) % 4]))
            counts[key] = counts.get(key, 0) + 1

    out = []
    for face in faces:
        corners = np.asarray(face.corners, dtype=float)
        flags = []
        for i in range(4):
            key = (face.axis, round(float(face.position), 9),
                   _edge_key(corners[i], corners[(i + 1) % 4]))
            flags.append(counts[key] == 1)
        out.append(flags)
    return out


def sort_faces(faces, camera: Camera, bounds):
    """Faces back to front, and their projected corners.

    Ordering by the centroid's depth is exact for the disjoint axis-aligned
    rectangles produced here: none of them intersect, so no face is
    partly in front of and partly behind another.
    """
    entries = []
    for face in faces:
        screen, depth = camera.project(face.corners, bounds)
        entries.append((float(np.mean(depth)), face, screen))
    entries.sort(key=lambda item: item[0])
    return [(face, screen) for _, face, screen in entries]


def visible(faces, camera: Camera, bounds, *, cull: bool = True):
    """Drop the faces that point away from the viewer.

    A face is kept when its recorded outward normal has a non-negative
    component towards the camera; a face seen exactly edge-on is kept rather
    than dropped, so a view straight down an axis does not lose a wall. A
    two-sided face (a free-standing slice plane) is always kept.
    """
    if not cull:
        return list(faces)
    basis = camera.basis()
    kept = []
    for face in faces:
        if not face.outward:
            kept.append(face)
            continue
        # The axis scaling is positive on every axis, so it cannot change
        # the sign of an axis-aligned normal's depth component.
        if float(np.dot(face.normal(), basis[2])) >= -1e-9:
            kept.append(face)
    return kept


# ==========================================================================
# View-aligned resampling, and what it makes possible
# ==========================================================================
def view_volume(values, axes, camera: Camera, *, samples: int = 180,
                bounds=None) -> np.ndarray:
    """Resample the volume onto a grid aligned with the camera.

    Returns ``(right, up, depth)`` with depth increasing **towards** the
    viewer, so a projection along the last axis is a projection along the
    line of sight. Doing this once means the maximum-intensity projection,
    the alpha composite and the isosurface all come from the same array and
    always agree with each other.
    """
    values = np.asarray(values, dtype=float)
    bounds = _bounds(axes) if bounds is None else bounds
    basis = camera.basis()
    # The unit cube's half-diagonal, so the view grid always contains it
    # whatever the angles.
    reach = 0.5 * float(np.linalg.norm(camera.aspect))
    grid = np.linspace(-reach, reach, int(samples))
    gr, gu, gd = np.meshgrid(grid, grid, grid, indexing="ij")
    local = (gr[..., None] * basis[0] + gu[..., None] * basis[1]
             + gd[..., None] * basis[2])

    coords = []
    for d in range(3):
        lo, hi = bounds[d]
        span = (hi - lo) or 1.0
        world = (local[..., d] / float(camera.aspect[d]) + 0.5) * span + lo
        coords.append(_index_coords(axes[d], world))
    stacked = np.array([c.ravel() for c in coords])

    mask = np.isfinite(values)
    filled = np.where(mask, values, 0.0)
    num = ndimage.map_coordinates(filled, stacked, order=1, mode="constant", cval=0.0)
    den = ndimage.map_coordinates(mask.astype(float), stacked, order=1,
                                  mode="constant", cval=0.0)
    inside = np.ones(stacked.shape[1], dtype=bool)
    for d in range(3):
        inside &= (stacked[d] >= 0) & (stacked[d] <= values.shape[d] - 1)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = np.where((den > 0.5) & inside, num / np.maximum(den, 1e-12), np.nan)
    return out.reshape(gr.shape)


def project_volume(cube, *, mode: str = "mip", levels=None, alpha: float = 0.06):
    """Collapse a view-aligned cube into one image.

    ``"mip"`` takes the maximum along the line of sight -- the simplest and
    most honest volume view, since every pixel is a real measured value.
    ``"sum"`` integrates. ``"alpha"`` composites front to back with an
    opacity proportional to intensity, which shows a Fermi surface inside
    the volume rather than only its brightest shell; ``levels`` sets the
    intensity window it maps opacity over.
    """
    cube = np.asarray(cube, dtype=float)
    if mode == "mip":
        with np.errstate(invalid="ignore"):
            return np.nanmax(np.where(np.isfinite(cube), cube, -np.inf), axis=2)
    if mode == "sum":
        return np.nansum(cube, axis=2)
    if mode != "alpha":
        raise ValueError(f"projection mode {mode!r} is not mip, sum or alpha")

    finite = cube[np.isfinite(cube)]
    if levels is None:
        levels = ((float(finite.min()), float(finite.max())) if finite.size
                  else (0.0, 1.0))
    lo, hi = float(levels[0]), float(levels[1])
    span = (hi - lo) or 1.0
    scaled = np.clip((np.nan_to_num(cube, nan=lo) - lo) / span, 0.0, 1.0)
    opacity = np.clip(scaled * float(alpha), 0.0, 1.0)

    # Front to back: the depth axis increases towards the viewer, so walk
    # it backwards and stop accumulating once the ray is opaque.
    out = np.zeros(cube.shape[:2])
    transmitted = np.ones(cube.shape[:2])
    for k in range(cube.shape[2] - 1, -1, -1):
        contribution = transmitted * opacity[:, :, k]
        out += contribution * scaled[:, :, k]
        transmitted *= (1.0 - opacity[:, :, k])
    return out


def isosurface_depth(cube, level: float):
    """Where each line of sight first crosses ``level``, and a shading.

    A full marching-cubes mesh is not needed to *look at* an isosurface, and
    it would add a dependency the lab may not have. The first crossing along
    each ray gives a depth map; its gradient gives a surface normal; Lambert
    shading of that normal is the picture. Returns ``(depth, shade)``, both
    NaN where no crossing was found.
    """
    cube = np.asarray(cube, dtype=float)
    filled = np.nan_to_num(cube, nan=-np.inf)
    above = filled >= float(level)
    # Depth increases towards the viewer, so the *last* index that is above
    # the level is the nearest surface the ray meets.
    any_hit = above.any(axis=2)
    nearest = cube.shape[2] - 1 - np.argmax(above[:, :, ::-1], axis=2)
    depth = np.where(any_hit, nearest.astype(float), np.nan)

    smooth = ndimage.gaussian_filter(np.nan_to_num(depth, nan=0.0), 1.5,
                                     mode="nearest")
    gx, gy = np.gradient(smooth)
    normal = np.stack([-gx, -gy, np.ones_like(smooth)], axis=-1)
    normal /= np.linalg.norm(normal, axis=-1, keepdims=True)
    light = np.array([-0.4, -0.5, 0.75])
    light = light / np.linalg.norm(light)
    shade = np.clip(normal @ light, 0.0, 1.0) * 0.8 + 0.2
    return depth, np.where(any_hit, shade, np.nan)


# ==========================================================================
# 3-D symmetrisation
# ==========================================================================
def symmetrise_volume(values, x_axis, y_axis, z_axis, *, fold: int = 1,
                      centre=(0.0, 0.0), mirror_angles=(), inversion: bool = False,
                      sector=None, progress=None):
    """Symmetrise every constant-``z`` plane of a volume about the same axis.

    The one real fix to ``data_symmetrization_pro.m``: it takes its missing
    -data mask from the first z plane alone and applies it to all of them,
    with a comment admitting the problem. For a photon-energy scan or a
    converted k-map the valid region is different in every plane -- that is
    what the conversion does -- so the mask is taken per plane here.

    The origin is the caller's and is held for every plane: a per-plane
    origin would make the volume shear.
    """
    from tools.process import symmetrise

    values = np.asarray(values, dtype=float)
    if values.ndim != 3:
        raise ValueError(f"this works on a 3-D volume, not {values.ndim}-D")
    z_axis = np.asarray(z_axis, dtype=float)

    out = np.empty_like(values)
    coverage = np.empty(values.shape, dtype=float)
    for k in range(values.shape[2]):
        result = symmetrise(values[:, :, k], x_axis, y_axis, fold=fold,
                            centre=centre, mirror_angles=mirror_angles,
                            inversion=inversion, sector=sector)
        out[:, :, k] = result.values
        coverage[:, :, k] = result.coverage
        if progress is not None and not progress(k + 1, values.shape[2]):
            out[:, :, k + 1:] = np.nan
            coverage[:, :, k + 1:] = 0.0
            break
    return out, coverage, (float(centre[0]), float(centre[1]))


def reduce_volume(values, axes, ranges=None, targets=None):
    """Crop a cube to ``ranges`` and bin it down to about ``targets`` points.

    ``ranges`` is one ``(lo, hi)`` per axis (either end may be None) in the
    axes' own units; ``targets`` is the number of points wanted on each axis
    afterwards, or None to leave that axis alone.

    Binning is a NaN-aware **mean** rather than a sum: this is for looking
    at, and a sum changes the levels every time the factor changes, so the
    colour scale would move under the user as they traded resolution for
    speed. Blocks are whole, and a remainder that does not fill one is
    dropped, which keeps every output point made of the same number of
    inputs.

    Returns ``(values, axes)``.
    """
    values = np.asarray(values, dtype=float)
    axes = [np.asarray(a, dtype=float) for a in axes]
    if values.ndim != len(axes):
        raise ValueError("one axis per dimension, please")

    if ranges is not None:
        index = []
        for axis, span in zip(axes, ranges):
            if span is None:
                index.append(slice(None))
                continue
            lo = -np.inf if span[0] is None else float(span[0])
            hi = np.inf if span[1] is None else float(span[1])
            keep = np.flatnonzero((axis >= min(lo, hi)) & (axis <= max(lo, hi)))
            if keep.size == 0:
                raise ValueError("that range leaves nothing on one of the axes")
            index.append(slice(int(keep[0]), int(keep[-1]) + 1))
        values = values[tuple(index)]
        axes = [axis[sl] for axis, sl in zip(axes, index)]

    if targets is not None:
        for dim, target in enumerate(targets):
            if target is None:
                continue
            factor = int(max(1, round(values.shape[dim] / max(int(target), 1))))
            if factor <= 1:
                continue
            blocks = values.shape[dim] // factor
            if blocks < 1:
                continue
            kept = blocks * factor
            values = np.take(values, np.arange(kept), axis=dim)
            shape = list(values.shape)
            shape[dim:dim + 1] = [blocks, factor]
            reshaped = values.reshape(shape)
            with np.errstate(invalid="ignore"):
                counts = np.isfinite(reshaped).sum(axis=dim + 1)
                totals = np.nansum(reshaped, axis=dim + 1)
                values = np.where(counts > 0, totals / np.maximum(counts, 1),
                                  np.nan)
            axes[dim] = axes[dim][:kept].reshape(blocks, factor).mean(axis=1)
    return values, axes


def apply_plane_wise(values, function, *, axis: int = 2, progress=None):
    """Run a 2-D operation on every plane of a volume.

    Smoothing, derivatives and curvature are all defined on a plane; a
    volume is processed by doing each plane in turn, which is also what
    makes the progress bar possible and the memory bounded.
    """
    values = np.asarray(values, dtype=float)
    axis = int(axis) % 3
    moved = np.moveaxis(values, axis, 0)
    out = np.empty_like(moved)
    for i in range(moved.shape[0]):
        out[i] = function(moved[i])
        if progress is not None and not progress(i + 1, moved.shape[0]):
            out[i + 1:] = np.nan
            break
    return np.moveaxis(out, 0, axis)
