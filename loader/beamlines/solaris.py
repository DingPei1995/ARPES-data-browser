"""
loader/beamlines/solaris.py
===========================
SOLARIS URANOS (Scienta DA30-L): SES files, read by
:mod:`loader.formats.ses`:

============================  =============================================
file                          opens as
============================  =============================================
``.ibw`` cut                  ``cut``
``.ibw`` manipulator scan     ``spem_4d``: ``Point`` steps ``X:``,
                              ``Position`` steps ``Y:``; SES numbers both
                              ``1, 2, 3 ...`` and the motor readbacks in the
                              run-mode table become the axes
``.zip`` deflector map        ``map`` (Thetay)
``.txt`` photon-energy scan   ``kz_map``: one ``[Region n]`` per photon
                              energy, each referred to E_F by its own hv
                              (``E_kin - hv + W``) and resampled onto the
                              range they share when the windows differ
============================  =============================================

Written against the headers in ``for loader/SOLARIS`` and a text scan cut
down to its first photon energy. Recognised by the header's ``Location``
(``URANOS`` in the ``.ibw``/``.zip``, ``Solaris`` in the text export).

Nothing in here imports Qt.
"""
from __future__ import annotations

import re

from loader.formats.ses import SESLoader
from loader.registry import register

__all__ = ["UranosLoader", "is_uranos"]

_URANOS = re.compile(r"uranos|solaris", re.IGNORECASE)


def is_uranos(header: dict) -> bool:
    return bool(_URANOS.search(header.get("Location", "")))


class UranosLoader(SESLoader):
    name = "SOLARIS URANOS"
    description = ("Scienta SES files from the URANOS beamline at SOLARIS: "
                   ".ibw cuts and spatial scans, .zip deflector maps, .txt "
                   "photon-energy scans.")
    patterns = ("*.ibw", "*.zip", "*.txt")
    reads_text = True
    # Above the CASSIOPEE readers, which would otherwise take the .txt.
    priority = 22
    facility = "SOLARIS URANOS"
    prefix = "SOLARIS_"

    def claims(self, header: dict) -> bool:
        return is_uranos(header)


register(UranosLoader())
