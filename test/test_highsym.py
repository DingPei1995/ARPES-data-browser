"""The named high-symmetry points, and the cut offsets that pass through them.

The table (tools.highsym.POINT_TABLE) is the part that can be wrong without
looking wrong, so the first check is geometric: every named point must lie
on the surface of the zone built independently by tools.bz3d. The rest pin
the rhombohedral case the tool is most often used for, with the cell of
space group 166, a = 5.416 Å, c = 13.413 Å, where the distances along c*
are known in closed form (L at pi/c, F at 2 pi/c, Z at 3 pi/c).
"""
import collections

import numpy as np
import pytest

from tools.bz3d import face_planes, tile_and_cut, plane_basis, plane_cut, \
    wigner_seitz_cell
from tools.highsym import (equivalent_points, high_symmetry_points,
                           lattice_variant, plane_offsets, points_in_frame,
                           points_on_plane)
from tools.lattice import (LatticeParams, conventional_vectors,
                           primitive_vectors, reciprocal_vectors)

R166 = LatticeParams(5.416, 5.416, 13.413, 90, 90, 120, space_group=166)

CELLS = {
    "hR, alpha < 90": (R166, "RHL1"),
    "hR, alpha > 90": (LatticeParams(5, 5, 3, 90, 90, 120, space_group=166), "RHL2"),
    "cF": (LatticeParams(4, 4, 4, space_group=225), "FCC"),
    "cI": (LatticeParams(3, 3, 3, space_group=229), "BCC"),
    "cP": (LatticeParams(3, 3, 3, space_group=221), "CUB"),
    "hP": (LatticeParams(3, 3, 5, 90, 90, 120, space_group=194), "HEX"),
    "hP trigonal": (LatticeParams(3, 3, 5, 90, 90, 120, space_group=164), "HEX"),
    "tP": (LatticeParams(3, 3, 5, space_group=123), "TET"),
    "tI, c < a": (LatticeParams(4, 4, 3, space_group=139), "BCT1"),
    "tI, c > a": (LatticeParams(3, 3, 8, space_group=139), "BCT2"),
    "oP": (LatticeParams(3, 4, 5, space_group=47), "ORC"),
    "mC": (LatticeParams(3, 4, 5, 90, 100, 90, space_group=12), None),
    "aP": (LatticeParams(3, 4, 5, 70, 80, 100, space_group=2), None),
}


def _zone(params):
    b = reciprocal_vectors(primitive_vectors(params))
    return wigner_seitz_cell(b), b


def _normal(params, hkl):
    normal = np.asarray(hkl, dtype=float) @ reciprocal_vectors(
        conventional_vectors(params))
    return normal / np.linalg.norm(normal)


@pytest.mark.parametrize("name", list(CELLS))
def test_every_named_point_lies_on_the_zone_surface(name):
    params, variant = CELLS[name]
    assert lattice_variant(params) == variant
    zone, _b = _zone(params)
    normals, offsets = face_planes(zone.faces)
    for point in high_symmetry_points(params):
        excess = normals @ point.k - offsets
        if point.label == "Γ":
            assert np.allclose(point.k, 0.0)
        else:
            assert abs(excess.max()) < 1e-7, (point.label, excess.max())


@pytest.mark.parametrize("name", list(CELLS))
def test_copies_are_distinct_and_one_name_per_place(name):
    params, _variant = CELLS[name]
    copies = equivalent_points(params)
    places = [tuple(np.round(p.k, 6)) for p in copies]
    assert len(places) == len(set(places))


def test_rhombohedral_166_points_and_their_counts():
    counts = collections.Counter(p.label for p in equivalent_points(R166))
    # the four kinds of time-reversal-invariant point: 1 Γ, 2 Z, 6 F, 6 L
    assert counts["Γ"] == 1 and counts["Z"] == 2
    assert counts["F"] == 6 and counts["L"] == 6
    # second names for the same places are not listed again
    assert "L1" not in counts and "B1" not in counts and "P2" not in counts


def test_rhombohedral_166_offsets_along_c_star_are_multiples_of_pi_over_c():
    c = R166.c
    presets = dict((round(o, 6), labels) for o, labels in
                   plane_offsets(equivalent_points(R166), _normal(R166, (0, 0, 1))))
    offsets = sorted(presets)
    assert offsets == pytest.approx([0.0, np.pi / c, 2 * np.pi / c, 3 * np.pi / c],
                                    abs=1e-6)
    assert "Γ" in presets[offsets[0]]
    assert "L" in presets[offsets[1]]
    assert "F" in presets[offsets[2]]
    assert "Z" in presets[offsets[3]]


def test_rhombohedral_166_z_is_the_top_face_and_a_cut_there_is_that_face():
    zone, _b = _zone(R166)
    normal = _normal(R166, (0, 0, 1))
    z = next(p for p in high_symmetry_points(R166) if p.label == "Z")
    assert np.allclose(z.k[:2], 0.0, atol=1e-12)
    assert z.k[2] == pytest.approx(3 * np.pi / R166.c)
    u, v = plane_cut(zone.faces, normal, z.k)
    assert u is not None and len(u) - 1 == 6          # the hexagonal top face
    # and just above it there is nothing
    u, _v = plane_cut(zone.faces, normal, z.k * (1 + 1e-6))
    assert u is None


@pytest.mark.parametrize("hkl", [(0, 0, 1), (1, 0, 0), (1, 1, 0), (1, 0, 4)])
def test_rhombohedral_166_every_offered_offset_holds_its_points(hkl):
    copies = equivalent_points(R166)
    zone, _b = _zone(R166)
    normal = _normal(R166, hkl)
    for offset, labels in plane_offsets(copies, normal):
        held = {p.label for p in points_on_plane(copies, normal, offset)}
        assert set(labels) == held
        u, _v = plane_cut(zone.faces, normal, normal * offset)
        # Inside the zone's extent the plane cuts a polygon; at the far end
        # it can only touch a corner, which is still a point worth offering
        # (the copies beside it are cut there) but has no polygon of its own.
        if offset < float(np.max(zone.vertices @ normal)) - 1e-6:
            assert u is not None


def test_rhombohedral_166_labels_in_the_plane_of_a_tiled_cut():
    copies = equivalent_points(R166)
    _zone_, b = _zone(R166)
    normal = _normal(R166, (0, 0, 1))
    u_axis, v_axis = plane_basis(normal)
    found = points_in_frame(copies, b, np.zeros(3), u_axis, v_axis,
                            (-1.5, 1.5), (-1.5, 1.5))
    labels = collections.Counter(label for label, _x, _y in found)
    assert labels["Γ"] >= 1
    # every labelled point really is in the plane and in the window
    for _label, x, y in found:
        assert -1.5 - 1e-6 <= x <= 1.5 + 1e-6 and -1.5 - 1e-6 <= y <= 1.5 + 1e-6
    # nearest the window's centre first, and the centre is Γ itself
    assert found[0][0] == "Γ" and np.hypot(found[0][1], found[0][2]) < 1e-9


def test_tiling_finds_every_copy_when_the_plane_is_slid_far_out():
    """The cut used to search lattice indices around Gamma only, so at a
    large offset it found too few copies -- or none."""
    for params in (R166, LatticeParams(4.14, 4.14, 28.64, 90, 90, 120,
                                       space_group=166)):
        zone, b = _zone(params)
        normal = _normal(params, (0, 0, 1))
        for offset in (0.0, 2.0, 6.0):
            point = normal * offset
            found = tile_and_cut(zone, b, normal, point, (-2, 2), (-2, 2))
            brute = 0
            for i in range(-12, 13):
                for j in range(-12, 13):
                    for k in range(-12, 13):
                        shift = i * b[0] + j * b[1] + k * b[2]
                        if abs(np.dot(shift, normal) - offset) > 1.0:
                            continue
                        u, v = plane_cut([f + shift for f in zone.faces],
                                         normal, point)
                        if u is None:
                            continue
                        if (u.max() < -2.2 or u.min() > 2.2
                                or v.max() < -2.2 or v.min() > 2.2):
                            continue
                        brute += 1
            assert len(found) == brute > 0, (params.c, offset)


@pytest.fixture(scope="module")
def qapp():
    pytest.importorskip("PyQt5")
    from PyQt5.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def test_the_dialog_offers_and_applies_the_offsets(qapp):
    from ui import windows as W
    from ui.widgets import MemoryData
    k = np.linspace(-1.2, 1.2, 41)
    data = MemoryData("k_map", (k, k, np.linspace(-0.5, 0.1, 7)),
                      np.random.default_rng(0).random((41, 41, 7)) + 1,
                      {"x": "kx", "y": "ky", "z": "E"}, source_label="t")
    contour = W.ContourWindow(data, "t", "gray", False)
    dialog = contour.open_brillouin_zone()
    dialog.spacegroup_box.setValue(166)
    dialog.a_box.setValue(5.416)
    dialog.c_box.setValue(13.413)
    combo = dialog.offset_points_combo
    texts = [combo.itemText(i) for i in range(combo.count())]
    assert any(t.startswith("Γ") for t in texts)
    z_index = next(i for i, t in enumerate(texts) if "Z" in t.split("(")[0])
    combo.setCurrentIndex(z_index)
    dialog._offset_point_chosen(z_index)
    # the exact offset, not the four-decimal box value, so the plane is
    # in the top face rather than a hair above it
    assert dialog.offset_value() == pytest.approx(3 * np.pi / 13.413, abs=1e-12)
    assert "Z" in {label for label, _k in dialog.points_in_plane()}
    preview = dialog.open_preview()
    assert "Z" in preview.note.text()
    message = dialog.plot()
    assert message and "Z" in message, dialog.note.text()
    assert contour.view._bz_labels
    dialog.labels_cb.setChecked(False)
    dialog.plot()
    assert not contour.view._bz_labels
    # typing a different offset by hand drops the preset
    dialog.offset_box.setValue(0.1)
    assert dialog.offset_value() == pytest.approx(0.1)
    preview.close()
    dialog.close()
    contour.close()


def _k_map_contour():
    from ui import windows as W
    from ui.widgets import MemoryData
    k = np.linspace(-1.2, 1.2, 21)
    data = MemoryData("k_map", (k, k, np.linspace(-0.5, 0.1, 5)),
                      np.ones((21, 21, 5)), {"x": "kx", "y": "ky", "z": "E"},
                      source_label="t")
    return W.ContourWindow(data, "t", "gray", False)


def test_the_dialog_starts_from_the_last_crystal_used(qapp):
    """Same experiment, same crystal: a second map opens the dialog on the
    space group and cell typed for the first, though not on its offset."""
    first = _k_map_contour()
    dialog = first.open_brillouin_zone()
    dialog.spacegroup_box.setValue(166)
    dialog.a_box.setValue(5.416)
    dialog.c_box.setValue(13.413)
    dialog.offset_box.setValue(0.2)
    dialog.close()
    first.close()

    second = _k_map_contour()
    dialog = second.open_brillouin_zone()
    assert dialog.spacegroup_box.value() == 166
    assert dialog.a_box.value() == pytest.approx(5.416)
    assert dialog.c_box.value() == pytest.approx(13.413)
    assert dialog.offset_box.value() == 0.0
    # a dataset's own saved settings still win
    dialog.close()
    second.close()
    third = _k_map_contour()
    third.data.scan.info.update({"bz.space_group": 225, "bz.a": 4.0})
    dialog = third.open_brillouin_zone()
    assert dialog.spacegroup_box.value() == 225
    assert dialog.a_box.value() == pytest.approx(4.0)
    dialog.close()
    third.close()
