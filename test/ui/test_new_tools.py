"""The tool-panel round from the MATLAB comparison: the panel's new rows and
Calculation tab, and the windows behind the new buttons -- summing repeated
scans, ridge detection, detector corrections, the tight-binding fit, the
crystal structure -- plus Save for Igor."""
import os

import numpy as np
import pytest
from PyQt5.QtWidgets import QApplication

from tools.crystal import tightbinding as TB
from ui.launcher.list_actions import availability
from ui.viewers.widgets import MemoryData



@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


@pytest.fixture
def launcher(app, monkeypatch):
    import ARPES_viewer as av
    from ui.launcher import main_window as main_window_ui
    monkeypatch.setattr(av, "loaded_items", {})
    win = av._MainWindow()
    ui = main_window_ui.Ui_MainWindow()
    ui.setupUi(win)
    monkeypatch.setattr(av, "win", win)
    monkeypatch.setattr(av, "ui", ui)
    monkeypatch.setattr(av, "autosave", lambda data, name: None)
    av.SetConnect()
    yield av
    win.hide()
    win.deleteLater()


def _cut(name="c", value=None, info=None):
    k = np.linspace(-1, 1, 40)
    e = np.linspace(-0.5, 0.1, 61)
    value = np.ones((k.size, e.size)) if value is None else value
    return MemoryData("cut", (k, e), value, {"x": "Angle (deg)",
                                             "y": "Kinetic energy (eV)"},
                      source_label=name, source_info=dict(info or {}))


def _k_map():
    vectors = TB.rotate(TB.lattice_2d(3.5, 3.5, 90), 3)
    model = TB.TBFit(vectors, TB.shells(vectors, 2), 0.1, [-0.3, 0.05], 0,
                     [0, 0], 0, 0, False)
    kx = np.linspace(-1, 1, 121)
    energy = np.linspace(-0.6, 0.1, 71)
    band = TB.evaluate(model, kx, kx)
    cube = 0.03 / ((energy[None, None, :] - band[:, :, None]) ** 2 + 0.03 ** 2)
    return MemoryData("k_map", (kx, kx, energy), cube,
                      {"x": "kx (Å⁻¹)", "k": "ky (Å⁻¹)", "z": "E - E_F (eV)"},
                      source_label="km")


# -- the panel ---------------------------------------------------------------
def test_memory_buttons_share_a_row_with_short_names(launcher):
    ui = launcher.ui
    texts = [b.text() for b in (ui.FreeMemoryButton, ui.CloseViewersButton,
                                ui.OpenAutosaveButton)]
    assert texts == ["Free RAM", "Close all", "Autosaves"]
    rows = [ui.MemoryBox.layout().itemAt(i).layout()
            for i in range(ui.MemoryBox.layout().count())]
    row = next(r for r in rows if r is not None
               and r.indexOf(ui.FreeMemoryButton) >= 0)
    assert [row.indexOf(b) for b in (ui.FreeMemoryButton, ui.CloseViewersButton,
                                     ui.OpenAutosaveButton)] == [0, 1, 2]


def test_the_new_tools_sit_where_they_were_asked_for(launcher):
    ui = launcher.ui
    basic, advanced = ui.ToolSections["Basic"], ui.ToolSections["Advanced"]
    calculation = ui.ToolSections["Calculation"]
    assert basic.isAncestorOf(ui.ToolButtons["detcorr"])
    column = basic.layout()
    assert (column.indexOf(ui.ToolButtons["detcorr"])
            == column.indexOf(ui.ToolButtons["degrid"]) + 1)
    assert advanced.isAncestorOf(ui.ToolButtons["ridge"])
    assert calculation.isAncestorOf(ui.ToolButtons["lattice"])
    assert calculation.isAncestorOf(ui.ToolButtons["tbfit"])
    assert "sumrepeat" in ui.DatasetButtons


def test_availability_of_the_new_actions():
    cut, cube, kmap = ("a", "Cut"), ("b", "Map"), ("k", "k-map")
    assert availability("sumrepeat", [cut, cut])[0]
    assert not availability("sumrepeat", [cut])[0]
    assert not availability("sumrepeat", [cut, cube])[0]
    assert availability("ridge", [cube])[0] and availability("ridge", [cut])[0]
    assert availability("detcorr", [cube])[0]
    assert not availability("detcorr", [kmap])[0]
    assert availability("tbfit", [])[0] and availability("tbfit", [kmap])[0]
    assert not availability("tbfit", [cut])[0]
    assert availability("lattice", [cut, cube])[0]


# -- Save for Igor -------------------------------------------------------------
def test_save_writes_igor_text_and_keeps_the_rows(launcher, tmp_path):
    lv = launcher
    lv.add_memory_dataset(_cut("one"), "one")
    lv.add_memory_dataset(_cut("two"), "two")
    lv.ui.FilePathListWidget.selectAll()
    path = lv.save_selected(str(tmp_path / "both"), "itx")
    text = open(path, encoding="utf-8-sig").read()
    assert path.endswith(".itx")
    assert "'one'" in text and "'two'" in text
    assert all(r["memory"] is not None for r in lv.loaded_items.values())
    written = lv.save_selected(str(tmp_path / "w.ibw"), "ibw")
    assert os.path.exists(str(tmp_path / "w_one.ibw"))
    assert written.endswith(".ibw")


# -- the windows -----------------------------------------------------------------
def test_sum_repeats_window_walks_through_its_steps(app):
    from ui.panels.repeats import RepeatsDialog
    rng = np.random.default_rng(1)
    base = 100 + 50 * np.exp(-np.linspace(-1, 1, 40)[:, None] ** 2 / 0.1)
    scans = [(f"s{i}", _cut(f"s{i}", rng.poisson(base * np.ones((40, 61)))
                            .astype(float), {"photon_energy_eV": 60.0}))
             for i in range(3)]
    dialog = RepeatsDialog(scans)
    assert dialog.stepper.refresh() == ""
    dialog.stepper.go(1)
    assert dialog.stepper.refresh()                 # compare first
    dialog.compare(background=False)
    dialog.stepper.go(2)
    made = []
    dialog.datasetsCreated.connect(made.extend)
    dialog.stepper.next_button.click()
    assert made and made[0].scan.value.sum() == pytest.approx(
        sum(d.scan.value.sum() for _n, d in scans))
    assert made[0].scan.info["repeats.members"] == ["s0", "s1", "s2"]


def test_sum_repeats_refuses_two_kinds(app):
    from ui.panels.repeats import RepeatsDialog
    cube = MemoryData("map", (np.arange(2.0), np.arange(3.0), np.arange(4.0)),
                      np.ones((2, 3, 4)), {}, source_label="m")
    with pytest.raises(ValueError):
        RepeatsDialog([("c", _cut()), ("m", cube)])


def test_ridge_window_lists_ridge_images_and_saves_points(app, tmp_path):
    from ui.panels.ridge import RidgeDialog
    k = np.linspace(-0.6, 0.6, 121)
    e = np.linspace(-0.4, 0.05, 91)
    kk, ee = np.meshgrid(k, e, indexing="ij")
    image = 1 / (1 + ((ee - 1.2 * kk ** 2 + 0.3) / 0.02) ** 2)
    cut = MemoryData("cut", (k, e), image, {"x": "k (Å⁻¹)", "y": "E (eV)"},
                     source_label="c")
    dialog = RidgeDialog([("c", cut)])
    assert len(dialog.points) > 50
    made = dialog.results_to_list()
    assert made[0].kind == "cut" and made[0].scan.value.shape == image.shape
    path = dialog.save_points(str(tmp_path / "p.csv"))
    assert open(path, encoding="utf-8").readline().startswith("k (Å⁻¹),E (eV)")


def test_detector_window_corrects_a_map(app):
    from ui.panels.detector import DetectorDialog
    x = np.linspace(-15, 15, 61)
    y = np.linspace(50, 52, 41)
    frame = np.ones((61, 41))
    frame[30, 20] = 80
    cube = MemoryData("map", (np.arange(3.0), x, y), np.stack([frame] * 3),
                      {"x": "Deflector (deg)", "k": "Slit (deg)",
                       "z": "Kinetic energy (eV)"}, source_label="m")
    dialog = DetectorDialog(("m", cube))
    dialog.stepper.go(1)
    assert len(dialog.find_regions()) == 1
    dialog.stepper.go(2)
    result = dialog.apply_to_list(background=False)
    assert result.scan.value.shape == (3, 61, 41)
    assert np.allclose(result.scan.value, 1.0)
    assert result.scan.info["detcorr.fill"] == "smooth (harmonic)"


def test_detector_window_refuses_data_off_the_detector(app):
    from ui.panels.detector import DetectorDialog
    with pytest.raises(ValueError, match="cannot be corrected"):
        DetectorDialog(("c", _cut(info={"kcut.gamma": 0.0})))


def test_tight_binding_window_traces_fits_and_lists_a_model(app, tmp_path):
    from ui.panels.tightbinding import TightBindingDialog
    dialog = TightBindingDialog(("km", _k_map()))
    dialog.a2.setValue(3.5)
    dialog.b2.setValue(3.5)
    dialog.stepper.go(1)
    for energy in (0.0, -0.15, -0.3):
        dialog.energy_slider.setValue(dialog._nearest_energy(energy))
        assert dialog.trace_plane() > 50
    assert dialog.stepper.refresh() == ""
    dialog.stepper.go(2)
    dialog.fit_azimuth.setChecked(True)
    model = dialog.run_fit()
    assert model.azimuth_deg == pytest.approx(3.0, abs=0.2)
    assert model.t == pytest.approx([-0.3, 0.05], abs=0.005)
    dialog.stepper.go(3)
    result = dialog.model_to_list()
    assert result.kind == "k_map" and result.scan.value.shape == (121, 121, 71)
    assert os.path.exists(dialog.save_parameters(str(tmp_path / "tb.json")))


def test_tight_binding_from_a_file_alone(app, tmp_path):
    from ui.panels.tightbinding import TightBindingDialog
    vectors = TB.lattice_2d(3.0, 3.0, 120)
    model = TB.TBFit(vectors, TB.shells(vectors, 1), 0.0, [-0.5], 0, [0], 0, 0,
                     False)
    rng = np.random.default_rng(5)
    k = rng.uniform(-1, 1, (60, 2))
    path = tmp_path / "pts.csv"
    path.write_text("kx,ky,E\n" + "\n".join(
        f"{a},{b},{e}" for (a, b), e in zip(k, model.energy(k))))
    dialog = TightBindingDialog(None)
    dialog.a2.setValue(3.0)
    dialog.b2.setValue(3.0)
    dialog.gamma2.setValue(120.0)
    dialog.n_shells.setValue(1)
    assert dialog.read_file(str(path)) == 60
    fit = dialog.run_fit()
    assert fit.t[0] == pytest.approx(-0.5, abs=1e-6)
    assert dialog.model_to_list().kind == "k_map"


def test_tight_binding_takes_the_brillouin_zone_crystal(app, tmp_path, monkeypatch):
    import json
    from ui.panels.tightbinding import TightBindingDialog
    monkeypatch.setenv("ARPES_VIEWER_CONFIG_DIR", str(tmp_path))
    (tmp_path / "brillouin_zone.json").write_text(json.dumps(
        {"mode": "2d", "a2": 2.46, "b2": 2.46, "gamma2": 120.0,
         "azimuth_deg": 15.0}))
    dialog = TightBindingDialog(None)
    assert dialog.take_from_bz()
    assert dialog.a2.value() == pytest.approx(2.46)
    assert dialog.azimuth.value() == pytest.approx(15.0)


def test_crystal_window_reads_a_cif_draws_a_plane_and_exports(app, tmp_path):
    from test.tools.test_new_analysis import NACL_CIF
    from ui.panels.crystal import CrystalDialog
    cif = tmp_path / "nacl.cif"
    cif.write_text(NACL_CIF, encoding="utf-8")
    dialog = CrystalDialog()
    assert dialog.read_cif(str(cif)) is not None
    assert len(dialog.structure.atoms) == 8
    dialog.stepper.go(1)
    assert dialog.style_table.rowCount() == 2
    dialog.set_style("Na", gloss=0.1, alpha=0.5, size=1.5)
    dialog.stepper.go(2)
    dialog.plane_box.setChecked(True)
    for box, value in zip(dialog.hkl_boxes, (1, 1, 1)):
        box.setValue(value)
    assert dialog.layer_list.count() > 3
    dialog.set_plane_height(dialog.layer_list.item(2).data(0x0100))
    dialog.above_combo.setCurrentIndex(2)
    image = dialog.view.render_image(300, 200)
    assert (image.width(), image.height()) == (300, 200)
    dialog.stepper.go(3)
    assert os.path.getsize(dialog.save_picture(str(tmp_path / "x.png"))) > 1000
    window = dialog.to_figure_composer(ask=False)
    assert window.figure.panels[0].data.kind == "picture"
    window.close()


def test_crystal_window_typed_in_and_centred(app):
    from ui.panels.crystal import CrystalDialog
    dialog = CrystalDialog()
    dialog.space_group.setValue(229)          # Im-3m: body centred
    assert len(dialog.structure.atoms) == 2
    dialog.add_atom("O", 0.5, 0, 0)
    assert len(dialog.structure.atoms) == 4
