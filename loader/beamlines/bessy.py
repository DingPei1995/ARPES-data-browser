"""
loader/beamlines/bessy.py
=========================
BESSY II's two Scienta end stations: **One Cube** (1^3, DA30-L) and **One
Square** (1^2, R8000). Both write SES files, read by
:mod:`loader.formats.ses`:

==============  ============================  ==============================
end station     file                          opens as
==============  ============================  ==============================
One Cube        ``.ibw`` cut                  ``cut``
                ``.ibw`` photon-energy scan   ``kz_map`` (binding energy)
                ``.zip`` deflector map        ``map`` (Thetay)
One Square      ``.ibw`` cut                  ``cut``
                ``.ibw`` polar map            ``map`` (``P-Axis``, a polar
                                              angle, convertible to k)
                ``.ibw`` photon-energy scan   ``kz_map`` (binding energy)
                ``.ibw`` manipulator scan     ``spem_4d`` (``Position`` = Z,
                                              ``Point`` = polar)
==============  ============================  ==============================

Written against the headers in ``for loader/BESSY``. One Square names
itself in the header (``Location=ARPES 1^2``); One Cube's acquisition PC
leaves SES's default ``Location=Scienta``, so it is also recognised by its
analyser (``Instrument=DA30L_8000-11623``). Either can be chosen by hand in
the Loader window for a file that says neither.

Nothing in here imports Qt.
"""
from __future__ import annotations

import re

from loader.formats.ses import SESLoader
from loader.registry import register

__all__ = ["OneCubeLoader", "OneSquareLoader", "is_one_cube",
           "is_one_square"]

_CUBE = re.compile(r"1\s*\^\s*3|1³|one\s*cube|1-cubed", re.IGNORECASE)
_SQUARE = re.compile(r"1\s*\^\s*2|1²|one\s*square|1-squared", re.IGNORECASE)
#: The analysers' serials as their headers state them.
CUBE_INSTRUMENTS = ("DA30L_8000-11623",)
SQUARE_INSTRUMENTS = ("R8000-8ES202",)


def is_one_cube(header: dict) -> bool:
    return bool(_CUBE.search(header.get("Location", ""))) \
        or header.get("Instrument", "").strip() in CUBE_INSTRUMENTS


def is_one_square(header: dict) -> bool:
    return bool(_SQUARE.search(header.get("Location", ""))) \
        or header.get("Instrument", "").strip() in SQUARE_INSTRUMENTS


class OneCubeLoader(SESLoader):
    name = "BESSY One Cube"
    description = ("Scienta SES files from the 1^3 end station at BESSY II: "
                   ".ibw cuts and photon-energy scans, .zip deflector maps.")
    patterns = ("*.ibw", "*.zip", "*.txt")
    reads_text = True
    priority = 22
    facility = "BESSY One Cube"
    prefix = "BESSY1C_"

    def claims(self, header: dict) -> bool:
        return is_one_cube(header)


class OneSquareLoader(SESLoader):
    name = "BESSY One Square"
    description = ("Scienta SES files from the 1^2 end station at BESSY II: "
                   ".ibw cuts, polar maps, photon-energy scans and "
                   "manipulator scans.")
    patterns = ("*.ibw", "*.zip", "*.txt")
    reads_text = True
    priority = 22
    facility = "BESSY One Square"
    prefix = "BESSY1S_"

    def claims(self, header: dict) -> bool:
        return is_one_square(header)


register(OneCubeLoader())
register(OneSquareLoader())
