"""Draw the icon of ``ARPES viewer.app`` (the macOS launcher).

A cut through a hole band and an electron pocket, with the Fermi cut-off,
in the "hot" colormap, on a rounded square as macOS icons are. Written
as ``ARPES viewer.app/Contents/Resources/AppIcon.icns`` (every size macOS
asks for, by Pillow). Nothing here is part of the program; run it by hand
to change the icon.

Usage:  python devtools/make_mac_icon.py
"""
import os
import sys

import numpy as np
from PIL import Image, ImageDraw

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from tools.output.colormaps import get_lut                 # noqa: E402

OUT = os.path.join(HERE, "ARPES viewer.app", "Contents", "Resources",
                   "AppIcon.icns")
SIZE = 1024


def spectrum(n: int = SIZE) -> np.ndarray:
    """I(k, E) on an n x n grid, scaled 0..1: E from -1 (bottom) to 0.25."""
    k = np.linspace(-1.0, 1.0, n)[None, :]
    energy = np.linspace(0.25, -1.0, n)[:, None]
    image = np.zeros((n, n))
    for centre, width in ((0.05 - 1.1 * k ** 2, 0.035),        # hole band
                          (-0.75 + 1.6 * k ** 2, 0.035)):      # electron pocket
        lorentz = width / ((energy - centre) ** 2 + width ** 2)
        image += lorentz * (1.0 + 2.0 * np.abs(energy))
    image *= 1.0 / (1.0 + np.exp(energy / 0.03))               # Fermi cut-off
    image = np.sqrt(image / image.max())
    return image


def icon() -> Image.Image:
    lut = np.asarray(get_lut("hot"))[:, :3]
    values = spectrum()
    rgb = lut[np.clip((values * (len(lut) - 1)).astype(int), 0, len(lut) - 1)]
    picture = Image.fromarray(rgb.astype(np.uint8), "RGB").convert("RGBA")
    # macOS icons: a rounded square, inset from the canvas edge.
    inset, radius = SIZE // 10, SIZE // 5
    mask = Image.new("L", (SIZE, SIZE), 0)
    ImageDraw.Draw(mask).rounded_rectangle(
        (inset, inset, SIZE - inset, SIZE - inset), radius, fill=255)
    picture = picture.resize((SIZE, SIZE))
    out = Image.new("RGBA", (SIZE, SIZE), (0, 0, 0, 0))
    out.paste(picture, (0, 0), mask)
    return out


def main():
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    icon().save(OUT, format="ICNS")
    print("saved", OUT)


if __name__ == "__main__":
    main()
