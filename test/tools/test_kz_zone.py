"""The Brillouin zone as a k_z map sees it, for any surface normal.

A k_z map is the plane holding the surface normal and the slit direction.
The zone drawn on it used to be vertical lines at multiples of the period
along the normal -- Gamma points, not boundaries, and the same picture
whatever the surface. These pin the real cut, the dialogs that draw it, and
the lattice they build it from.
"""
import os

import numpy as np
import pytest

from tools.crystal.bz3d import kz_frame, kz_zone_cuts, wigner_seitz_cell
from tools.crystal.lattice import (LatticeParams, primitive_vectors,
                                   reciprocal_vectors)


def _zone(params):
    b = reciprocal_vectors(primitive_vectors(params))
    return wigner_seitz_cell(b), b


def _centres_on_axis(polygons):
    return sorted(round(float(xs[:-1].mean()), 4) for xs, ys in polygons
                  if abs(ys[:-1].mean()) < 1e-6)


def test_hexagonal_001_is_rectangles_with_period_2pi_over_c():
    params = LatticeParams(3.2, 3.2, 6.0, 90, 90, 120, space_group=194)
    zone, b = _zone(params)
    polygons = kz_zone_cuts(zone, b, kz_frame(params, (0, 0, 1)),
                            (0.0, 3.0), (-1.0, 1.0))
    centres = _centres_on_axis(polygons)
    assert np.allclose(np.diff(centres), 2 * np.pi / 6.0, atol=1e-4)
    assert all(len(xs) - 1 == 4 for xs, _ys in polygons)
    # half-width along k_par is Gamma-M, |a*| / 2
    xs, ys = polygons[0]
    assert np.isclose(np.max(np.abs(ys)), 2 * np.pi / 3.2 / np.sqrt(3), atol=1e-6)


def test_bcc_001_repeats_every_4pi_over_a_not_2pi():
    a = np.pi
    params = LatticeParams(a, a, a, space_group=229)
    zone, b = _zone(params)
    polygons = kz_zone_cuts(zone, b, kz_frame(params, (0, 0, 1)),
                            (-0.5, 9.0), (-0.5, 0.5))
    assert np.allclose(np.diff(_centres_on_axis(polygons)), 4 * np.pi / a)


def test_fcc_111_cuts_are_hexagons_along_the_111_period():
    a = np.pi
    params = LatticeParams(a, a, a, space_group=225)
    zone, b = _zone(params)
    polygons = kz_zone_cuts(zone, b, kz_frame(params, (1, 1, 1)),
                            (0.0, 8.0), (-1.0, 1.0))
    assert np.allclose(np.diff(_centres_on_axis(polygons)),
                       2 * np.pi * np.sqrt(3) / a, atol=1e-4)
    assert {len(xs) - 1 for xs, _ys in polygons} == {6}


def test_the_frame_follows_the_normal_and_the_kpar_choice():
    params = LatticeParams(4.0, 4.0, 4.0, space_group=221)
    normal, along, across = kz_frame(params, (1, 1, 0), (0, 0, 1))
    assert np.allclose(normal, [1 / np.sqrt(2), 1 / np.sqrt(2), 0])
    assert np.allclose(along, [0, 0, 1])
    assert np.allclose(np.cross(normal, along), across)
    _n, turned, _x = kz_frame(params, (1, 1, 0), (0, 0, 1), 90.0)
    assert np.allclose(turned, across)
    with pytest.raises(ValueError):
        kz_frame(params, (1, 1, 0), (2, 2, 0))
    with pytest.raises(ValueError):
        kz_frame(params, (0, 0, 0))


def test_untiled_keeps_the_copy_in_the_window():
    params = LatticeParams(3.0, 3.0, 3.0, space_group=221)
    zone, b = _zone(params)
    one = kz_zone_cuts(zone, b, kz_frame(params, (0, 0, 1)), (5.5, 7.0),
                       (-0.5, 0.5), tile=False)
    assert len(one) == 1
    assert np.isclose(one[0][0][:-1].mean(), 3 * 2 * np.pi / 3.0)


# -- the dialogs ------------------------------------------------------------
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


@pytest.fixture(scope="module")
def app():
    pytest.importorskip("PyQt5")
    from PyQt5.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def _kz_map(kind):
    from ui.viewers.widgets import MemoryData
    first = (np.linspace(40, 120, 21) if kind == "kz_map"
             else np.linspace(2.0, 6.0, 41))
    second = (np.linspace(-12, 12, 25) if kind == "kz_map"
              else np.linspace(-1.0, 1.0, 21))
    energy = np.linspace(-0.5, 0.1, 13)
    rng = np.random.default_rng(1)
    labels = {"x": "x", "k": "k", "z": "E - E_F (eV)"}
    return MemoryData(kind, (first, second, energy),
                      rng.random((first.size, second.size, energy.size)) + 1,
                      labels, source_label="kz")


def test_kz_conversion_window_draws_the_cut_for_any_normal(app):
    from ui.viewers import windows as W
    from ui.panels.kzconv import KzConversionDialog
    contour = W.ContourWindow(_kz_map("kz_map"), "kz", "gray", False)
    dialog = KzConversionDialog(contour)
    # the default hexagonal group gets its gamma = 120, not 90
    assert dialog.lattice().gamma == 120.0
    assert not dialog.lat_gamma.isEnabled() and not dialog.lat_b.isEnabled()
    assert dialog.surface_hkl() == (0, 0, 1)
    # a lattice edit keeps the chosen normal
    target = next(i for i in range(dialog.normal_index.count())
                  if tuple(dialog.normal_index.itemData(i)[0]) == (1, 0, 0))
    dialog.normal_index.setCurrentIndex(target)
    dialog.lat_c.setValue(6.5)
    assert dialog.surface_hkl() == (1, 0, 0)
    xs, ys = dialog.zone_item.getData()
    assert xs is not None and np.isfinite(xs).any()
    # cubic: b and c follow a
    dialog.space_group.setValue(225)
    dialog.lat_a.setValue(4.1)
    assert dialog.lat_b.value() == dialog.lat_c.value() == 4.1
    dialog.close()
    contour.close()


def test_brillouin_zone_dialog_is_offered_and_draws_on_a_kz_map(app):
    from ui.viewers import windows as W
    contour = W.ContourWindow(_kz_map("kz_map_k"), "kz", "gray", False)
    info = contour.data.scan.info
    info.update({"bz.space_group": 225, "bz.a": np.pi, "bz.h": 1,
                 "bz.k": 1, "bz.l": 1})
    dialog = contour.open_brillouin_zone()
    assert dialog._kz_map and not dialog.mode_2d.isVisible()
    message = dialog.plot()
    assert message and "polygon" in message, dialog.note.text()
    item = contour.view._bz_layers["bz_main"]
    xs, _ys = item.getData()
    finite = xs[np.isfinite(xs)]
    # drawn in k_z, inside (or at the edge of) the map's k_z range
    assert finite.max() > 2.0 and finite.min() < 6.0
    assert info["bz.kpar_h"] == 0
    dialog.close()
    contour.close()
