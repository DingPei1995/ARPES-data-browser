"""
tools/crystal/last_crystal.py
=============================
The crystal last typed into the Brillouin-zone tool, for the tools that
need the same crystal again: the tight-binding fit and the crystal-
structure viewer.

The Brillouin-zone dialog (``ui.viewers.windows.BrillouinZoneDialog``)
keeps its last settings in ``brillouin_zone.json`` beside the program's
config (``~/.arpes_viewer``, or ``ARPES_VIEWER_CONFIG_DIR``). This reads
that file and hands back the lattice in the forms the other tools use; it
never writes it.
"""
from __future__ import annotations

import json
import os

__all__ = ["path", "load", "lattice_params", "vectors_2d", "describe"]

FILE_NAME = "brillouin_zone.json"


def path() -> str:
    folder = os.environ.get("ARPES_VIEWER_CONFIG_DIR") or os.path.join(
        os.path.expanduser("~"), ".arpes_viewer")
    return os.path.join(folder, FILE_NAME)


def load() -> dict:
    """The saved settings, or ``{}`` when the tool has not been used."""
    try:
        with open(path(), encoding="utf-8") as handle:
            values = json.load(handle)
        return values if isinstance(values, dict) else {}
    except (OSError, ValueError):
        return {}


def lattice_params(values: dict = None):
    """The 3-D crystal as :class:`tools.crystal.lattice.LatticeParams`, or
    None when the last zone was a 2-D layer (or nothing was saved)."""
    from tools.crystal.lattice import LatticeParams

    values = load() if values is None else values
    if not values or values.get("mode", "3d") != "3d":
        return None
    return LatticeParams(float(values.get("a", 3.0)), float(values.get("b", 3.0)),
                         float(values.get("c", 3.0)),
                         float(values.get("alpha", 90.0)),
                         float(values.get("beta", 90.0)),
                         float(values.get("gamma", 90.0)),
                         space_group=int(values.get("space_group", 1)))


def vectors_2d(values: dict = None):
    """``(a, b, gamma_deg, azimuth_deg)`` of the in-plane lattice: the 2-D
    layer's own, or a 3-D crystal's a, b, gamma; None when nothing was
    saved."""
    values = load() if values is None else values
    if not values:
        return None
    azimuth = float(values.get("azimuth_deg", 0.0))
    if values.get("mode", "3d") == "2d":
        return (float(values.get("a2", 3.0)), float(values.get("b2", 3.0)),
                float(values.get("gamma2", 120.0)), azimuth)
    return (float(values.get("a", 3.0)), float(values.get("b", 3.0)),
            float(values.get("gamma", 90.0)), azimuth)


def describe(values: dict = None) -> str:
    """One line saying what was saved, for a button's tooltip."""
    values = load() if values is None else values
    if not values:
        return "The Brillouin-zone tool has not been used yet."
    if values.get("mode", "3d") == "2d":
        return (f"2-D layer: a = {values.get('a2')} Å, b = {values.get('b2')} Å, "
                f"gamma = {values.get('gamma2')}°")
    return (f"Space group {values.get('space_group')}: a = {values.get('a')}, "
            f"b = {values.get('b')}, c = {values.get('c')} Å, "
            f"alpha/beta/gamma = {values.get('alpha')}/{values.get('beta')}/"
            f"{values.get('gamma')}°")

