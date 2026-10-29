"""A min/max the user sets on the level bar holds while the slice moves.

Dragging the energy bar of a constant-E contour used to reset the level
bar to each new frame's full range. Now the window the user chose stays
for every position, and a pop-out snapshot opens on the same colour scale;
dragging the handles back to both ends (or 0%–100%)
goes back to following each frame's range.

Needs Qt; run headless with QT_QPA_PLATFORM=offscreen.
"""
import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PyQt5")
from PyQt5.QtWidgets import QApplication                           # noqa: E402

from arpes_viewer.ui.viewers.widgets import MemoryData  # noqa: E402
from arpes_viewer.ui.viewers import windows as W  # noqa: E402

DEFL = np.linspace(-10, 10, 21)
SLIT = np.linspace(-15, 15, 31)
ENERGY = np.linspace(-1.0, 0.2, 61)
LABELS = {"x": "Angle, deflector (deg)", "k": "Angle, slit (deg)", "z": "Energy (eV)"}


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


def _contour():
    rng = np.random.default_rng(5)
    # Intensity grows with energy, so every slice has a different range.
    cube = rng.random((DEFL.size, SLIT.size, ENERGY.size)) * np.arange(1, 62)
    return W.ContourWindow(MemoryData("map", (DEFL, SLIT, ENERGY), cube, LABELS,
                                      source_label="m"), "m", "gray", False)


def test_levels_survive_energy_moves(app):
    window = _contour()
    bar = window.panel.level_bar
    bar.region.setRegion((3.0, 7.0))
    for index in (0, 10, 60):
        window.e_control.set_index(index)
        assert bar.levels() == pytest.approx((3.0, 7.0))
        assert window.view._levels == pytest.approx((3.0, 7.0))
    window.close()


def test_full_range_follows_each_slice(app):
    window = _contour()
    bar = window.panel.level_bar
    bar.region.setRegion((3.0, 7.0))
    bar.clip_percentiles(0.0, 100.0)
    window.e_control.set_index(60)
    frame = window.view.last_frame
    assert bar.levels() == pytest.approx((float(frame.min()), float(frame.max())))
    window.close()


def test_pop_out_keeps_the_colour_scale(app):
    window = _contour()
    window.panel.level_bar.region.setRegion((3.0, 7.0))
    window.panel.level_bar.gamma_box.setValue(0.5)
    snapshot = window.pop_out_slice()
    assert snapshot.panel.level_bar.levels() == pytest.approx((3.0, 7.0))
    assert snapshot.view._levels == pytest.approx((3.0, 7.0))
    assert snapshot.view._gamma == pytest.approx(0.5)
    snapshot.close()
    window.close()
