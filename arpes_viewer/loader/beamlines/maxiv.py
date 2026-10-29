"""
loader/beamlines/maxiv.py
=========================
MAX IV Bloch (A branch, Scienta DA30-L): Scienta SES ``.ibw`` cuts,
photon-energy scans and manipulator scans, and ``.zip`` deflector maps.

The reading is :mod:`loader.formats.ses`; this is the facility's name on
it. As the lowest-priority SES reader it also opens any SES file no other
facility claims.

Nothing in here imports Qt.
"""
from __future__ import annotations

from arpes_viewer.loader.formats.ses import SESLoader
from arpes_viewer.loader.registry import register

__all__ = ["MaxIVLoader"]


class MaxIVLoader(SESLoader):
    """Bloch -- and, as the lowest-priority SES reader, any SES file no
    other facility claims."""
    name = "MAX IV Bloch"
    description = ("Scienta SES files from the Bloch beamline at MAX IV: "
                   ".ibw cuts, photon-energy scans and manipulator scans, "
                   "and .zip deflector maps.")
    priority = 12
    facility = "MAX IV Bloch"
    prefix = "MAXIV_"


register(MaxIVLoader())
