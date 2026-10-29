"""The figure composer's own file (.afig): saved, opened again and edited;
the "data for <name>" folder beside it, with each panel's processed data
and the original files it was traced back to; the launcher opening one;
and the Calculated-bands window's panels naming what they came from."""
import os

import numpy as np
import pytest
from PyQt5.QtWidgets import QApplication

import arpes_viewer.tools.output.figure as nf
from arpes_viewer.tools.output import figfile


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


def _measurement(tmp_path):
    """A cut written to a file, read back as the launcher reads it."""
    from arpes_viewer.loader.nxs_file import save_dataset
    from arpes_viewer.ui.viewers.widgets import NxsData

    path = str(tmp_path / "raw.nxs")
    k = np.linspace(-1, 1, 50)
    e = np.linspace(-1, 0.1, 40)
    values = np.exp(-((e[None] + 0.5 - 0.4 * k[:, None] ** 2) / 0.05) ** 2)
    save_dataset(path, [{"name": "cut1", "kind": "cut", "axes": (k, e),
                         "labels": {"x": "k (Å⁻¹)", "y": "E - E_F (eV)"},
                         "value": values, "info": {}}])
    return path, NxsData(path, "cut1")


def _processed(data):
    """A dataset made from ``data``, the way a processing step makes one."""
    from arpes_viewer.tools.processing.process import Step, record_step
    from arpes_viewer.ui.viewers.widgets import MemoryData

    scan = data.scan
    info = record_step(scan.info, Step("smooth", {"sigma": 0.02}, source="cut1"))
    return MemoryData("cut", (scan.x, scan.y), np.asarray(scan.value) * 2,
                      dict(scan.labels), source_label="cut1 smoothed",
                      source_info=info, source_path=data.path)


def _window(panels):
    from arpes_viewer.ui.panels.figure import FigureWindow
    return FigureWindow(panels, "fig")


def test_a_dataset_remembers_the_file_it_began_as(app, tmp_path):
    path, data = _measurement(tmp_path)
    derived = _processed(data)
    assert derived.scan.info["origin.path"] == os.path.abspath(path)
    assert derived.scan.info["origin.entry"] == "cut1"
    from arpes_viewer.ui.panels.figure import panel_from_dataset
    panel = panel_from_dataset(derived, "smoothed")
    source = panel.data.sources[0]
    assert source["path"] == os.path.abspath(path)
    assert source["history"] and source["history"][0].startswith("smooth")
    data.close()


def test_save_and_open_again_unchanged(app, tmp_path):
    from arpes_viewer.ui.panels.figure import open_figure_file, panel_from_dataset

    path, data = _measurement(tmp_path)
    panel = panel_from_dataset(_processed(data), "smoothed")
    panel.title = "LH"
    panel.lines.append(nf.Line(value=0.0, axis="y", color="#ff0000"))
    panel.texts.append(nf.TextNote(x=0.1, y=-0.2, text="α", color="#000000"))
    panel.overlays.append(nf.Overlay(xs=[0, 0.5, np.nan, 0.6], ys=[0, -0.2, 0, -0.3],
                                     name="band"))
    panel.x_tick_labels = [(-1.0, "M"), (0.0, "Γ")]
    panel.colorbar, panel.colorbar_ends = True, "high-low"
    window = _window([panel])
    window.figure.shapes.append(nf.Shape(kind="arrow", x0=1, y0=1, x1=5, y1=5))
    saved = window.save_file(str(tmp_path / "fig1"))
    assert saved.endswith(".afig")
    first = str(tmp_path / "a.png")
    window.write(first, "png")

    again = open_figure_file(saved, show=False)
    back = again.figure.panels[0]
    assert back.title == "LH" and back.colorbar_ends == "high-low"
    assert back.x_tick_labels == [(-1.0, "M"), (0.0, "Γ")]
    assert isinstance(back.lines[0], nf.Line) and back.texts[0].text == "α"
    assert np.isnan(back.overlays[0].xs[2])
    assert np.array_equal(back.data.array, panel.data.array)
    assert again.figure.shapes[0].kind == "arrow"
    second = str(tmp_path / "b.png")
    again.write(second, "png")
    from PIL import Image
    assert np.array_equal(np.asarray(Image.open(first)), np.asarray(Image.open(second)))
    data.close()


def test_the_data_folder_keeps_processed_and_original(app, tmp_path):
    from arpes_viewer.ui.panels.figure import panel_from_dataset

    path, data = _measurement(tmp_path)
    panel = panel_from_dataset(_processed(data), "smoothed")
    panel.overlays.append(nf.Overlay(xs=[0, 1], ys=[0, -1], name="calc"))
    picture = str(tmp_path / "sketch.png")
    from PyQt5.QtGui import QImage
    image = QImage(20, 10, QImage.Format_RGB32)
    image.fill(0xff3366)
    image.save(picture)
    window = _window([panel, nf.Panel(data=nf.empty_data())])
    window.figure.freeze()
    window.load_picture(window.figure.panels[1], picture)
    out = str(tmp_path / "out")
    saved = window.save_file(os.path.join(out, "fig1.afig"))
    folder = figfile.data_folder_of(saved)
    assert os.path.basename(folder) == "data for fig1"
    a = os.path.join(folder, "a")
    assert os.path.isfile(os.path.join(a, "original", "raw.nxs"))
    assert os.path.isfile(os.path.join(a, "processed", "smoothed.nxs"))
    assert os.path.isfile(os.path.join(a, "processed", "lines.csv"))
    text = open(os.path.join(a, "sources.txt"), encoding="utf-8").read()
    assert "entry: cut1" in text and "then: smooth" in text
    from arpes_viewer.loader import registry
    processed = registry.load(os.path.join(a, "processed", "smoothed.nxs"))
    assert np.allclose(processed.value, panel.data.array)
    b = os.path.join(folder, "b")
    assert os.path.isfile(os.path.join(b, "original", "sketch.png"))
    assert os.path.isfile(os.path.join(b, "processed", "sketch.png"))
    data.close()


def test_a_reopened_figure_keeps_its_originals_when_they_moved(app, tmp_path):
    from arpes_viewer.ui.panels.figure import open_figure_file, panel_from_dataset

    path, data = _measurement(tmp_path)
    window = _window([panel_from_dataset(data, "raw")])
    saved = window.save_file(str(tmp_path / "one" / "fig.afig"))
    data.close()
    os.remove(path)                                    # the original is gone
    again = open_figure_file(saved, show=False)
    resaved = again.save_file(str(tmp_path / "two" / "fig.afig"))
    kept = os.path.join(figfile.data_folder_of(resaved), "a", "original", "raw.nxs")
    assert os.path.isfile(kept)


def test_a_3d_view_is_saved_as_its_picture(app, tmp_path):
    from arpes_viewer.ui.panels.figure import open_figure_file
    from arpes_viewer.ui.panels.volume import VolumeWindow

    rng = np.random.default_rng(0)
    view = VolumeWindow(rng.random((20, 15, 12)), (np.arange(20.0), np.arange(15.0),
                                                   np.arange(12.0)),
                        ("x", "y", "E"), "cube", mode="slices")
    window = _window([view.figure_panel()])
    window.figure.freeze()
    saved = window.save_file(str(tmp_path / "v.afig"), with_data=False)
    back = open_figure_file(saved, show=False).figure.panels[0]
    assert back.data.kind == "picture" and back.data.array.ndim == 3
    view.close()


def test_the_launcher_opens_a_saved_figure(app, tmp_path, monkeypatch):
    import arpes_viewer.app as lv
    from arpes_viewer.ui.panels.figure import panel_from_arrays

    window = _window([panel_from_arrays(np.ones((3, 3)), [0, 1, 2], [0, 1, 2])])
    saved = window.save_file(str(tmp_path / "f.afig"), with_data=False)
    monkeypatch.setattr(lv, "figure_windows", [])
    monkeypatch.setattr(lv, "win", None)
    lv.add_files([saved], background=False)
    assert len(lv.figure_windows) == 1
    assert lv.figure_windows[0].file_path == saved
    lv.figure_windows[0].close()


def test_not_a_figure_is_refused(tmp_path):
    import zipfile
    import json
    bad = str(tmp_path / "x.afig")
    with zipfile.ZipFile(bad, "w") as archive:
        archive.writestr("figure.json", json.dumps({"format": "other"}))
        archive.writestr("arrays.npz", b"")
    with pytest.raises(Exception):
        figfile.load(bad)


def test_the_bands_window_fits_and_names_its_sources(app, tmp_path):
    """A calculation file and a measured cut in k: Z and the shift found,
    and the panels name the measurement's file and the calculation's."""
    from arpes_viewer.tools.crystal import bandfiles
    from arpes_viewer.tools.crystal import bands as B
    from arpes_viewer.ui.panels.bands import BandsDialog
    from arpes_viewer.ui.viewers.widgets import MemoryData

    a = 3.3
    lattice = np.array([[a, 0, 0], [-a / 2, a * np.sqrt(3) / 2, 0], [0, 0, 6.0]])
    hop = {(0, 0, 0): np.array([[-0.75 + 0.72]]) / 0.8 + 0.05 / 0.8}
    for R in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (1, 1, 0), (-1, -1, 0)):
        hop[R] = np.array([[-0.12 / 0.8]])
    calc = str(tmp_path / "c_hr.dat")
    bandfiles.write_hr(calc, lattice, hop, 1, fermi=0.0)
    truth = bandfiles.read(calc)
    k = np.linspace(-1.0, 1.0, 161)
    e = np.linspace(-1.0, 0.1, 111)
    band = 0.8 * B.bands_along_cut(truth, k, 90.0, kz=0.0)[0] - 0.05
    values = 1 / ((e[None] - band[:, None]) ** 2 + 0.02 ** 2)
    cut = MemoryData("cut", (k, e), values, {"x": "k (Å⁻¹)", "y": "E - E_F (eV)"},
                     source_label="cut", source_info={"origin.path": calc + ".raw"})
    dialog = BandsDialog([("cut", cut)])
    dialog.open_file(calc)
    dialog.stepper.go(2)
    dialog.compare_combo.setCurrentText("cut")
    dialog.direction_box.setValue(90.0)
    dialog.kz_plane_mode.setChecked(True)
    result = dialog.fit_alignment()
    assert result.z == pytest.approx(0.8, abs=0.02)
    assert result.shift == pytest.approx(-0.05, abs=0.01)
    panel = dialog.panel_compare()
    kinds = [s["kind"] for s in panel.data.sources]
    assert kinds == ["dataset", "calculation"]
    assert panel.data.sources[1]["companions"][0].endswith("c.win")
    spectrum = dialog.build_spectrum()
    assert spectrum.kind == "cut" and spectrum.scan.value.shape == (k.size, e.size)


def test_a_calculation_alone_goes_to_the_output(app, tmp_path):
    """No measurement selected: the FPLO tools' band.dat along a path and a
    Hamiltonian open, and their panels (the bands along the path, the
    contours at E_F) come out without a compare."""
    from arpes_viewer.tools.crystal import bandfiles
    from arpes_viewer.ui.panels.bands import BandsDialog

    a = 3.3
    lattice = np.array([[a, 0, 0], [-a / 2, a * np.sqrt(3) / 2, 0], [0, 0, 6.0]])
    hop = {(0, 0, 0): np.array([[0.5]])}
    for R in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (1, 1, 0), (-1, -1, 0)):
        hop[R] = np.array([[-0.15]])
    hr = str(tmp_path / "wannier90_hr.dat")
    bandfiles.write_hr(hr, lattice, hop, 1)
    path = str(tmp_path / "band.dat")
    with open(path, "w") as handle:
        handle.write("".join(f"  {d:.4f}  {-0.5 + d:.4f}\n" for d in np.linspace(0, 1, 11)))
        handle.write("\n")
    dialog = BandsDialog([])
    assert dialog.open_file(path) is not None
    assert not dialog.panel_boxes["path"].isChecked()
    dialog.stepper.go(3)
    panels, problems = dialog.panels()
    assert len(panels) == 1 and not problems
    dialog.open_file(hr)
    dialog.set_goal("contours")
    assert dialog.panel_boxes["contours"].isChecked()
    dialog.series_energies.setText("0")
    panels, problems = dialog.panels()
    assert len(panels) == 1 and not problems                 # the contours alone
    assert any(o.name.startswith("band") for o in panels[0].overlays)


def test_the_stacked_slices_view_draws_layers_and_a_mark(app):
    from arpes_viewer.ui.panels.volume import VolumeWindow

    rng = np.random.default_rng(1)
    axes = (np.linspace(-1, 1, 30), np.linspace(-1, 1, 30), np.linspace(-1, 0, 40))
    window = VolumeWindow(rng.random((30, 30, 40)), axes, ("k_x", "k_y", "E - E_F"),
                          "m", mode="stack")
    window.set_stack(2, [-0.8, -0.5, -0.2], width=0.05)
    assert len(window.scene.faces) == 3
    window.set_marks([[0, 0.5, -0.5], [0, 0.4, -0.2]], "DP")
    assert window.scene.marks[0]["text"] == "DP"
    image = window.scene.grab().toImage()
    assert not image.isNull()
    panel = window.figure_panel()
    assert panel.data.kind == "scene3d" and panel.data.scene.marks
    window.close()


def test_points_with_error_bars_are_drawn(app):
    from arpes_viewer.ui.panels.figure import FigureWindow
    from arpes_viewer.ui.viewers.curves import curve_panel

    panel = curve_panel([0, 1, 2], [("c", np.array([0.0, 1.0, 0.5]), "#000000")],
                        x_label="x", y_label="y", name="p")
    panel.overlays.append(nf.Overlay(xs=[0.5, 1.5], ys=[0.2, 0.8], color="#ff0000",
                                     marker="o", line=False, y_errors=[0.1, 0.2]))
    window = FigureWindow([panel], "f")
    image = window.painter.to_image()
    assert not image.isNull()
    tools = window._tools
    tools.overlay_list.setCurrentRow(0)
    tools.curve_marker.setCurrentIndex(2)        # squares
    tools.restyle_curves([0])
    assert panel.overlays[0].marker == "s"


def test_arrange_aligns_spaces_and_sizes_panels(app):
    from arpes_viewer.ui.panels.figure import FigureWindow, panel_from_arrays

    panels = [panel_from_arrays(np.ones((3, 3)), [0, 1, 2], [0, 1, 2]) for _ in range(3)]
    window = FigureWindow(panels, "f")
    figure = window.figure
    figure.freeze()
    for panel, frame in zip(figure.panels, ((5, 10, 30, 20), (45, 13, 25, 24),
                                            (100, 8, 20, 20))):
        panel.frame = frame
    assert window.arrange("top", [0, 1, 2])
    assert {p.frame[1] for p in figure.panels} == {8}
    window.arrange("space-x", [0, 1, 2])
    gaps = [figure.panels[i + 1].frame[0] - (figure.panels[i].frame[0] + figure.panels[i].frame[2])
            for i in range(2)]
    assert gaps[0] == pytest.approx(gaps[1])
    window.arrange("height", [0, 1, 2])
    assert {p.frame[3] for p in figure.panels} == {20}
    assert not window.arrange("left", [0])          # one panel: nothing to do


def test_letters_of_a_row_share_a_height(app):
    from arpes_viewer.ui.panels.figure import FigurePainter, FigureWindow, panel_from_arrays

    a = panel_from_arrays(np.ones((3, 3)), [0, 1, 2], [0, 1, 2])
    b = panel_from_arrays(np.ones((3, 3)), [0, 1, 2], [0, 1, 2])
    b.title = "with a title"
    window = FigureWindow([a, b], "f")
    figure = window.figure
    figure.freeze()
    a.frame, b.frame = (10, 12, 30, 30), (60, 13, 30, 30)
    painter = FigurePainter(figure)
    painter.to_image(dpi=100)
    ya, yb = (painter._label_at[id(p)][1] for p in (a, b))
    assert ya == pytest.approx(yb)
    figure.style.panel_label_aligned = False
    painter.to_image(dpi=100)
    assert painter._label_at == {}


def test_a_colour_bar_can_lie_under_the_panel(app):
    from arpes_viewer.ui.panels.figure import FigureWindow, panel_from_arrays

    panel = panel_from_arrays(np.arange(9.0).reshape(3, 3), [0, 1, 2], [0, 1, 2])
    panel.colorbar, panel.colorbar_position, panel.colorbar_length = True, "below", 0.5
    panel.colorbar_ends = "high-low"
    window = FigureWindow([panel], "f")
    image = window.painter.to_image(dpi=150)
    assert not image.isNull()
    _left, bottom, _top, _right = window.figure.margins_mm()
    panel.colorbar_position = "right"
    assert window.figure.margins_mm()[1] < bottom   # room kept below for it


def test_zone_labels_move_to_their_leftmost_equivalent(app, tmp_path):
    from arpes_viewer.ui.panels.crystal import CrystalDialog

    cif = tmp_path / "hex.cif"
    cif.write_text("""data_x
_cell_length_a 3.3
_cell_length_b 3.3
_cell_length_c 6.0
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 120
_space_group_IT_number 191
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Ti1 Ti 0 0 0
""")
    dialog = CrystalDialog()
    dialog.read_cif(str(cif))
    dialog.stepper.go(3)
    view = dialog.view
    right = view.camera.basis()[0]

    def label_x(name):
        items, _frame = view._zone_items()
        return next(float(i[1] @ right) for i in items if i[0] == "text" and i[2] == name)

    before = label_x("M")
    dialog.zone_left_edit.setText("M")
    dialog._zone_changed()
    assert label_x("M") < before
    items, _frame = view._zone_items()
    assert sum(1 for i in items if i[0] == "arrow") == 3      # k_x, k_y, k_z
    assert not view.render_image(300, 300).isNull()
