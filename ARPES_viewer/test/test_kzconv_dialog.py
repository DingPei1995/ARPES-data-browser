"""The kz -> momentum window builds the cell the space group describes.

It used to take a and c only and set b = a and every angle to 90 degrees,
whatever the space group: for a hexagonal group -- 194, the window's
default -- that is a tetragonal cell, so every surface normal other than
(001) was listed with the period of the wrong reciprocal lattice.
"""
import numpy as np
import pytest

pytest.importorskip("PyQt5")
from PyQt5.QtWidgets import QApplication                           # noqa: E402

from tools import cleavage                                          # noqa: E402
from tools.lattice import LatticeParams, validate_lattice_parameters  # noqa: E402
from ui.widgets import MemoryData                                  # noqa: E402
from ui import windows as W                                        # noqa: E402

HV = np.linspace(40.0, 80.0, 21)
SLIT = np.linspace(-15, 15, 31)
ENERGY = np.linspace(-1.0, 0.2, 61)


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


@pytest.fixture()
def dialog(app):
    # Opened the way the program opens it, and closed with its window (which
    # closes the dialogs it opened) -- closing the dialog on its own first
    # left objects for the collector to reclaim later, in whatever window
    # test came next.
    rng = np.random.default_rng(5)
    data = MemoryData("kz_map", (HV, SLIT, ENERGY),
                      rng.random((HV.size, SLIT.size, ENERGY.size)) + 1,
                      {"x": "Photon energy (eV)", "k": "Angle, slit (deg)",
                       "z": "E - E_F (eV)"}, source_label="kz")
    window = W.ContourWindow(data, "kz", "gray", False)
    box = window.open_kz_conversion()
    app.processEvents()
    yield box
    window.close()
    app.processEvents()


def _listed(box):
    """{(h, k, l): period} as the surface-normal box offers them."""
    out = {}
    for index in range(box.normal_index.count()):
        entry = box.normal_index.itemData(index)
        if entry is not None:
            out[tuple(int(v) for v in entry[0])] = float(entry[1])
    return out


def test_a_hexagonal_group_gets_a_hexagonal_cell(dialog):
    assert dialog.space_group.value() == 194
    assert dialog.lat_gamma.value() == pytest.approx(120.0)
    assert not dialog.lat_gamma.isEnabled() and not dialog.lat_b.isEnabled()
    dialog.lat_a.setValue(3.16)
    assert dialog.lat_b.value() == pytest.approx(3.16)            # b follows a
    lattice = dialog.lattice()
    assert validate_lattice_parameters(lattice) == []
    # (1 0 0) of a hexagonal cell repeats every 4 pi / (sqrt(3) a), not 2 pi / a
    periods = _listed(dialog)
    assert periods[(1, 0, 0)] == pytest.approx(4 * np.pi / (np.sqrt(3) * 3.16), rel=1e-6)
    assert periods[(0, 0, 1)] == pytest.approx(2 * np.pi / dialog.lat_c.value(), rel=1e-6)


def test_the_list_matches_the_lattice_it_claims(dialog):
    for group, a, b, c in ((194, 3.16, 3.16, 12.3), (139, 3.96, 3.96, 13.02),
                           (62, 5.4, 7.6, 5.5), (166, 4.38, 4.38, 30.5)):
        dialog.space_group.setValue(group)
        dialog.lat_a.setValue(a)
        if dialog.lat_b.isEnabled():
            dialog.lat_b.setValue(b)
        dialog.lat_c.setValue(c)
        lattice = dialog.lattice()
        assert validate_lattice_parameters(lattice) == [], group
        expected = {tuple(h): L for h, L, _u in
                    cleavage.reciprocal_lengths(lattice, max_index=2)[:12]}
        assert _listed(dialog) == pytest.approx(expected), group


def test_the_constraints_follow_the_space_group(dialog):
    dialog.space_group.setValue(139)                   # I4/mmm: tetragonal
    assert dialog.lat_gamma.value() == pytest.approx(90.0)
    assert not dialog.lat_gamma.isEnabled() and not dialog.lat_b.isEnabled()
    # body-centred: along (001) the first reflection is (002), 4 pi / c
    assert _listed(dialog)[(0, 0, 2)] == pytest.approx(4 * np.pi / dialog.lat_c.value())
    dialog.space_group.setValue(62)                    # Pnma: orthorhombic
    assert dialog.lat_b.isEnabled()                    # b is its own
    dialog.lat_b.setValue(7.6)
    assert dialog.lattice().b == pytest.approx(7.6) != dialog.lattice().a
    dialog.space_group.setValue(12)                    # C2/m: monoclinic, beta free
    assert dialog.lat_beta.isEnabled() and not dialog.lat_alpha.isEnabled()
    assert isinstance(dialog.lattice(), LatticeParams)
