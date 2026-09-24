"""
Where Gamma is on a constant-energy map -- a suggestion, never a fact.

In the viewer the k-space origin is picked by eye on the contour. A batch
run needs a number, so this offers one: the point about which the
constant-energy map is most nearly inversion-symmetric,
``I(c + d) ~ I(c - d)``. Photoemission intensity is only approximately
symmetric (matrix elements, polarisation), so the answer is reported with
a score between 0 and 1 and is meant to be looked at on the preview before
it is trusted. A map that does not contain Gamma has no such centre, and
its score says so.

The centre is found from the self-convolution of the image: ``(I * I)(m)``
peaks at ``m = 2c`` when ``I`` is symmetric about ``c``.
"""
from __future__ import annotations

import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.signal import fftconvolve


def energy_slice(cube, energy_axis, energy: float, width: float = 0.05):
    """Sum of the slices within ``energy +- width``."""
    energy_axis = np.asarray(energy_axis, dtype=float)
    inside = np.abs(energy_axis - float(energy)) <= float(width)
    if not inside.any():
        inside[int(np.argmin(np.abs(energy_axis - float(energy))))] = True
    return np.nansum(np.asarray(cube)[..., inside], axis=-1)


def _peak_subpixel(values, index):
    """Parabolic refinement of a peak at ``index`` along each axis."""
    out = []
    for dim, i in enumerate(index):
        if 0 < i < values.shape[dim] - 1:
            sl = list(index)
            sl[dim] = i - 1
            a = values[tuple(sl)]
            sl[dim] = i + 1
            c = values[tuple(sl)]
            b = values[tuple(index)]
            denom = a - 2 * b + c
            out.append(i + (0.5 * (a - c) / denom if denom else 0.0))
        else:
            out.append(float(i))
    return out


def symmetry_centre(image, axis0, axis1, *, smooth_px: float = 1.5):
    """``(c0, c1, score)``: the inversion centre in axis units, and how
    symmetric the image is about it (1 = perfectly)."""
    image = np.nan_to_num(np.asarray(image, dtype=float))
    if smooth_px:
        image = gaussian_filter(image, smooth_px)
    image = image - np.median(image)
    image = np.clip(image, 0, None)
    norm = float(np.sum(image ** 2))
    if norm <= 0:
        return float("nan"), float("nan"), 0.0
    conv = fftconvolve(image, image, mode="full")
    index = np.unravel_index(int(np.argmax(conv)), conv.shape)
    m0, m1 = _peak_subpixel(conv, index)
    score = float(conv[index] / norm)
    axis0, axis1 = np.asarray(axis0, dtype=float), np.asarray(axis1, dtype=float)
    c0 = float(np.interp(m0 / 2.0, np.arange(axis0.size), axis0))
    c1 = float(np.interp(m1 / 2.0, np.arange(axis1.size), axis1))
    return c0, c1, score


def suggest_centre(ds, *, energy: float = None, width: float = 0.05) -> dict:
    """For a map in angles: the deflector and slit angles of the symmetry
    centre at ``energy`` (default: 0.1 eV below the top of the energy axis
    after an E_F calibration, else the middle of the window)."""
    energy_axis = ds.axis("z")
    if energy is None:
        top = float(np.nanmax(energy_axis))
        energy = -0.1 if -0.5 < top < 1.0 else float(np.mean(energy_axis))
    cube = ds.array(np.float32)
    image = energy_slice(cube, energy_axis, energy, width)
    c0, c1, score = symmetry_centre(image, ds.axis("x"), ds.axis("k"))
    return {"theta_offset_deg": c0, "phi_offset_deg": c1, "score": score,
            "energy": float(energy), "width": float(width)}


def mirror_centre(image, axis_values, *, along: int = 0, smooth_px: float = 1.5):
    """``(centre, score)``: where ``image`` is most nearly mirror-symmetric
    along one axis -- normal emission on a cut's slit, or on the slit axis
    of a photon-energy scan. Each line across the other axis is
    self-convolved along ``along`` and the results summed, so every energy
    row votes for the same centre."""
    image = np.nan_to_num(np.asarray(image, dtype=float))
    image = np.moveaxis(image, along, 0)
    if smooth_px:
        image = gaussian_filter(image, smooth_px)
    image = np.clip(image - np.median(image), 0, None)
    norm = float(np.sum(image ** 2))
    if norm <= 0:
        return float("nan"), 0.0
    conv = np.zeros(2 * image.shape[0] - 1)
    for column in image.T:
        conv += fftconvolve(column, column, mode="full")
    index = int(np.argmax(conv))
    (m,) = _peak_subpixel(conv, (index,))
    axis_values = np.asarray(axis_values, dtype=float)
    centre = float(np.interp(m / 2.0, np.arange(axis_values.size), axis_values))
    return centre, float(conv[index] / norm)
