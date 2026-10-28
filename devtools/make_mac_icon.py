"""Make the icon of ``ARPES viewer.app`` (the macOS launcher).

The picture is a line drawing of a hemispherical analyser by the author
(``devtools/icon_drawing.png``: the electrons bent through the hemisphere and
landing on the detector, spread in energy and angle), set on a white rounded
square as macOS icons are. Written as
``ARPES viewer.app/Contents/Resources/AppIcon.icns``, in every size macOS
asks for (by Pillow). Nothing here is part of the program; run it by hand
after changing the drawing.

Usage:  python devtools/make_mac_icon.py
"""
import os

from PIL import Image, ImageDraw

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DRAWING = os.path.join(HERE, "devtools", "icon_drawing.png")
OUT = os.path.join(HERE, "ARPES viewer.app", "Contents", "Resources",
                   "AppIcon.icns")
SIZE = 1024
#: Apple's icon grid: an 824 px rounded square on a 1024 px canvas, with a
#: corner radius of about 185 px.
TILE, RADIUS = 824, 185
#: the drawing's margin inside the tile
MARGIN = 45


def icon() -> Image.Image:
    canvas = Image.new("RGBA", (SIZE, SIZE), (0, 0, 0, 0))
    inset = (SIZE - TILE) // 2
    box = (inset, inset, inset + TILE - 1, inset + TILE - 1)
    draw = ImageDraw.Draw(canvas)
    draw.rounded_rectangle(box, RADIUS, fill=(255, 255, 255, 255),
                           outline=(200, 204, 214, 255), width=6)
    drawing = Image.open(DRAWING).convert("RGB")
    side = TILE - 2 * MARGIN
    drawing.thumbnail((side, side), Image.LANCZOS)
    # The drawing's white ground is the tile's own: paste only the lines.
    lines = drawing.convert("L").point(lambda v: 255 - v)
    lines = lines.point(lambda v: min(255, v * 3))
    at = ((SIZE - drawing.width) // 2, (SIZE - drawing.height) // 2)
    canvas.paste(drawing, at, lines)
    return canvas


def main():
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    icon().save(OUT, format="ICNS")
    print("saved", OUT)


if __name__ == "__main__":
    main()
