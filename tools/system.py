"""
tools/system.py
============
Small, dependency-light stand-ins for the lab-internal ``tools_packages``
helpers used by the original ``*.h5`` loader (``myTools.save_img``,
``myTools.load_yaml``/``dump_yaml``, ``Win10Notif``). They are used by
``ARPES_viewer.py`` / ``ui/widgets.py`` so the nano-ARPES ``*.nxs``
viewer runs standalone on any machine with PyQt5 + pyqtgraph + numpy, without
needing SOLEIL's internal packages installed.

If you run this on the lab PC where ``tools_packages``/``pyNanoScanSystem``
are already installed, you can freely swap these calls back out for the
originals (e.g. to get the toast-style Win10Notif popups and the
shared-memory "Sync with running scan" feature) -- the function signatures
here were kept close to the originals for that reason.
"""
from __future__ import annotations

import json
import os
import numpy as np


# --------------------------------------------------------------------------
# Thumbnail / image export
# --------------------------------------------------------------------------
def normalize_to_uint8(array: np.ndarray) -> np.ndarray:
    arr = np.asarray(array, dtype=float)
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        return np.zeros(arr.shape, dtype=np.uint8)
    lo, hi = np.nanpercentile(finite, [1, 99])
    if hi <= lo:
        lo, hi = finite.min(), finite.max()
    if hi <= lo:
        return np.zeros(arr.shape, dtype=np.uint8)
    scaled = np.clip((arr - lo) / (hi - lo), 0, 1)
    scaled = np.nan_to_num(scaled, nan=0.0)
    return (scaled * 255).astype(np.uint8)


def save_img(path, array: np.ndarray, lut=None):
    """Save (or just render, if ``path`` is None) ``array`` as a PNG.

    ``lut``: optional (N,3) uint8 RGB lookup table -- pass the one the GUI is
    currently displaying (see ``tools.colormaps.get_lut``) so an exported image
    matches what's on screen. Without it this falls back to matplotlib's
    'inferno', then to plain greyscale if matplotlib isn't installed.

    Returns a PIL.Image so the caller can also put it on the clipboard.
    Mirrors the (path, array) -> PIL.Image signature of the original
    ``myTools.save_img``.
    """
    from PIL import Image
    u8 = normalize_to_uint8(array)

    if lut is not None:
        lut = np.asarray(lut, dtype=np.uint8)
        # Map the 0..255 normalized values onto however many entries the
        # colormap actually has (the lab's tables are 64/100/101 long).
        idx = (u8.astype(np.uint16) * (len(lut) - 1) // 255).astype(np.intp)
        img = Image.fromarray(lut[idx], mode="RGB")
    else:
        try:
            import matplotlib
            try:
                colormap = matplotlib.colormaps["inferno"]  # matplotlib >= 3.5
            except AttributeError:
                import matplotlib.cm as cm
                colormap = cm.get_cmap("inferno")  # matplotlib < 3.9 fallback
            colored = (colormap(u8) * 255).astype(np.uint8)  # RGBA
            img = Image.fromarray(colored, mode="RGBA").convert("RGB")
        except ImportError:
            img = Image.fromarray(u8, mode="L")

    if path is not None:
        img.save(path)
    return img


# --------------------------------------------------------------------------
# Tiny config persistence (last-used folder, window prefs, ...)
# --------------------------------------------------------------------------
def load_config(path: str, defaults: dict) -> dict:
    if os.path.isfile(path):
        try:
            with open(path, "r", encoding="utf-8") as fh:
                data = json.load(fh)
            merged = dict(defaults)
            merged.update(data)
            return merged
        except (json.JSONDecodeError, OSError):
            pass
    return dict(defaults)


def dump_config(path: str, values: dict) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(values, fh, indent=2)
