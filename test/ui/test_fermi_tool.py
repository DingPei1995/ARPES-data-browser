"""The Fermi edge fit as a tool (ui/viewers/windows.py FermiFitDialog,
ARPES_viewer.open_fermi_fit): a cut's fit, a map's fit on its slit cut --
one E_F for the whole map, or one per cut with the slices in between
interpolated -- and the fit and its window closing together."""
import numpy as np
import pytest
from PyQt5.QtWidgets import QApplication, QMessageBox

from tools.spectra.fermi import levels_along, fermi_function
from ui.viewers import windows as W
from ui.viewers.widgets import MemoryData

DEFL = np.linspace(-8, 8, 17)
SLIT = np.linspace(-10, 10, 41)
EKIN = np.linspace(55.0, 56.0, 201)


def _drift(defl):
    """The Fermi level of each slice: drifting by 40 meV across the map."""
    return 55.60 + 0.0025 * defl


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


@pytest.fixture(autouse=True)
def _no_modal_boxes(monkeypatch):
    def fail(parent, title, text, *args, **kwargs):
        raise AssertionError(f"{title}: {text}")
    for name in ("warning", "information", "critical"):
        monkeypatch.setattr(QMessageBox, name, staticmethod(fail))


def _map():
    rng = np.random.default_rng(5)
    cube = np.empty((DEFL.size, SLIT.size, EKIN.size))
    for i, d in enumerate(DEFL):
        cube[i] = 100 * fermi_function(EKIN, _drift(d), 20.0)[None, :] + 5
    cube += rng.normal(0, 0.3, cube.shape)
    return MemoryData("map", (DEFL, SLIT, EKIN), cube,
                      {"x": "Angle, deflector (deg)", "k": "Angle, slit (deg)",
                       "z": "Kinetic energy (eV)"}, source_label="m")


def _slit_fit(app):
    contour = W.ContourWindow(_map(), "m", "gray", False)
    cut = contour.open_cut("slit")
    made = []
    contour.datasetCreated.connect(made.append)
    dialog = cut.open_fermi_fit(bound=True)
    app.processEvents()
    return contour, cut, dialog, made


def _fit_at(cut, dialog, index):
    cut.control.set_index(index) if hasattr(cut.control, "set_index") else None
    if cut.control.index != index:
        cut.control.slider.setValue(index)
    assert cut.control.index == index
    result = dialog.run_fit()
    assert result is not None
    return result


def test_levels_between_fits_are_a_straight_line_held_at_the_ends():
    axis = np.array([-2.0, 0.0, 1.0, 2.0, 4.0, 6.0])
    levels = levels_along(axis, [4.0, 0.0], [1.4, 1.0])
    np.testing.assert_allclose(levels, [1.0, 1.0, 1.1, 1.2, 1.4, 1.4])
    with pytest.raises(ValueError):
        levels_along(axis, [], [])


def test_the_fit_follows_the_slit_cut_along_the_map(app):
    contour, cut, dialog, _ = _slit_fit(app)
    assert dialog.slices
    low = _fit_at(cut, dialog, 0).values["ef"]
    high = _fit_at(cut, dialog, DEFL.size - 1).values["ef"]
    assert low == pytest.approx(_drift(DEFL[0]), abs=3e-3)
    assert high == pytest.approx(_drift(DEFL[-1]), abs=3e-3)
    dialog.close()


def test_one_offset_for_the_whole_map(app):
    contour, cut, dialog, made = _slit_fit(app)
    _fit_at(cut, dialog, DEFL.size // 2)
    result = dialog.offset_energy_axis(whole_map=True)
    assert result is not None and made == [result]
    assert result.scan.labels["z"].startswith("E - E_F")
    np.testing.assert_allclose(result.scan.z, EKIN - dialog.fit_result.values["ef"])
    dialog.close()


def test_cut_by_cut_with_the_slices_between_interpolated(app):
    contour, cut, dialog, made = _slit_fit(app)
    _fit_at(cut, dialog, 2)
    assert dialog.offset_energy_axis(whole_map=False) is None
    assert dialog.per_cut and dialog.cuts_panel.isVisibleTo(dialog)
    _fit_at(cut, dialog, DEFL.size - 3)                 # added by itself
    assert [round(e["position"], 6) for e in dialog.cut_fits] == [DEFL[2], DEFL[-3]]
    # an edited E_F is used as typed
    item = dialog.cuts_table.item(0, 1)
    item.setText(f"{dialog.cut_fits[0]['ef']:.6f}")
    result = dialog.apply_cut_fits()
    assert result is not None and made[-1] is result
    energy = result.scan.z

    def edge(i):
        edc = result.scan.value[i].mean(axis=0)
        half = edc.min() + 0.5 * (edc.max() - edc.min())
        return energy[np.argmin(np.abs(edc - half))]

    # between the fits the line follows the drift: every edge at zero
    for i in (2, 5, 8, 12, DEFL.size - 3):
        assert abs(edge(i)) < 2e-3, (i, edge(i))
    # beyond them the nearest fit is held, so the drift past it remains
    for i, nearest in ((0, 2), (DEFL.size - 1, DEFL.size - 3)):
        expected = _drift(DEFL[i]) - _drift(DEFL[nearest])
        assert edge(i) == pytest.approx(expected, abs=2e-3)
    dialog.close()


def test_closing_either_closes_the_other(app):
    contour, cut, dialog, _ = _slit_fit(app)
    dialog.close()
    app.processEvents()
    assert not cut.isVisible()
    contour, cut, dialog, _ = _slit_fit(app)
    cut.close()
    app.processEvents()
    assert not dialog.isVisible()


def test_the_viewers_no_longer_offer_it(app):
    from test.ui.test_viewer_ui import _visible_texts
    contour = W.ContourWindow(MemoryData(
        "kz_map", (np.linspace(30, 60, 5), SLIT, EKIN),
        np.ones((5, SLIT.size, EKIN.size)),
        {"x": "Photon energy (eV)", "k": "Angle, slit (deg)",
         "z": "Kinetic energy (eV)"}, source_label="kz"), "kz", "gray", False)
    assert "Align Fermi levels (kz map)..." not in _visible_texts(contour)
    cut = W.CutWindow(MemoryData("cut", (SLIT, EKIN), np.ones((SLIT.size, EKIN.size)),
                                 {"x": "Angle (deg)", "y": "Kinetic energy (eV)"},
                                 source_label="c"), "c", "gray", False)
    assert "Fermi edge fit..." not in _visible_texts(cut)
    contour.close()
    cut.close()
