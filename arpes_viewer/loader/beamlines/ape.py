"""
loader/beamlines/ape.py
=======================
Elettra APE (the LE / spin-ARPES branch): Scienta SES files, as at MAX IV
Bloch.

Ported from the APE reader in peaks
(<https://github.com/phrgab/peaks>, ``peaks/core/fileIO/loaders/ape.py``),
which is the generic SES reader with APE's manipulator names -- ``polar``,
``tilt``, ``X``, ``Y``, ``Z``. Here too the reading is shared: ``.ibw``,
``.zip`` and ``.txt`` all go through :mod:`loader.formats.ses` (whose
run-mode table already knows APE's motor names, and whose text path uses
the Scienta parser of :mod:`loader.beamlines.cassiopee`).

The three formats are written the same way at every SES beamline, so the
facility is told by the ``Location`` the acquisition PC writes into the
header: this loader claims a file whose ``Location`` says APE or Elettra,
and the loader dialog can choose it for one that does not. Spin-resolved
data (the VLEED detectors) is not read, as in peaks.

Untested on real APE files: none were available when it was written.

Nothing in here imports Qt.
"""
from __future__ import annotations

import re

from arpes_viewer.loader.formats.ses import SESLoader
from arpes_viewer.loader.registry import register

__all__ = ["ApeLoader", "is_ape_location"]

_APE = re.compile(r"\bAPE\b|elettra", re.IGNORECASE)


def is_ape_location(location: str) -> bool:
    return bool(_APE.search(location or ""))


class ApeLoader(SESLoader):
    name = "Elettra APE"
    description = ("Scienta SES files from the APE beamline at Elettra: "
                   ".ibw and .txt spectra and scans, .zip deflector maps.")
    patterns = ("*.ibw", "*.zip", "*.txt")
    reads_text = True
    # Above the CASSIOPEE readers, which would otherwise take an APE .txt;
    # this one only claims a file whose header says APE.
    priority = 22
    facility = "Elettra APE"
    prefix = "APE_"

    def claims(self, header: dict) -> bool:
        return is_ape_location(header.get("Location", ""))


register(ApeLoader())
