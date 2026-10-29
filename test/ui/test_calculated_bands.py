"""The Calculated-bands window as an experimentalist meets it: say what to
see first (dispersion, contours, 3-D Fermi surface), find the file with
"Which file?" (renamed and cut-up files included), a progress bar for the
long reads, the 3-D Fermi surface drawn, the calculation saved to the main
panel with its high-symmetry points -- and the step-by-step windows asking
before they close and keeping their settings."""
import os

import numpy as np
import pytest
from PyQt5.QtWidgets import QApplication

from arpes_viewer.tools.crystal import bandfiles
from arpes_viewer.tools.crystal import bands as B

A = 3.3
LATTICE = np.array([[A, 0, 0], [-A / 2, A * np.sqrt(3) / 2, 0], [0, 0, 6.0]])


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


def _hoppings():
    hop = {(0, 0, 0): np.array([[0.5]])}
    for R in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (1, 1, 0), (-1, -1, 0)):
        hop[R] = np.array([[-0.15]])
    hop[(0, 0, 1)] = hop[(0, 0, -1)] = np.array([[0.05]])
    return hop


def _path_file(folder):
    path = os.path.join(folder, "band.dat")
    with open(path, "w") as handle:
        for shift in (-0.5, 0.2):
            handle.write("".join(f"  {d:.4f}  {shift + d:.4f}\n"
                                 for d in np.linspace(0, 1, 11)))
            handle.write("\n")
    return path


# --------------------------------------------------------------------------
# the files, by content
# --------------------------------------------------------------------------
def test_a_renamed_hamiltonian_and_a_cut_up_file_are_found(tmp_path):
    hr = str(tmp_path / "x_hr.dat")
    bandfiles.write_hr(hr, LATTICE, _hoppings(), 1)
    os.rename(hr, tmp_path / "hamiltonian.txt")
    os.rename(tmp_path / "x.win", tmp_path / "hamiltonian.win")
    assert bandfiles.detect(str(tmp_path / "hamiltonian.txt")) == "hr"
    model = bandfiles.read(str(tmp_path / "hamiltonian.txt"))
    assert isinstance(model, B.WannierModel)
    rows = {r["name"]: r for r in bandfiles.survey(str(tmp_path))}
    assert "unusual name" in " ".join(rows["hamiltonian.txt"]["notes"])
    assert rows["hamiltonian.txt"]["goals"] == {"dispersion", "contours", "fs"}
    assert rows["hamiltonian.win"]["kind"] == "win"

    # a band.dat cut in two by rows, as `split` would to pass an upload limit
    lines = open(_path_file(str(tmp_path))).read().splitlines(keepends=True)
    os.remove(tmp_path / "band.dat")
    (tmp_path / "band.dat.1").write_text("".join(lines[:9]))
    (tmp_path / "band.dat.2").write_text("".join(lines[9:]))
    rows = {r["name"]: r for r in bandfiles.survey(str(tmp_path))}
    assert rows["band.dat.1"]["pieces"] and "piece 1 of 2" in " ".join(rows["band.dat.1"]["notes"])
    joined = bandfiles.join_pieces(rows["band.dat.1"]["pieces"], str(tmp_path / "band.dat"))
    assert bandfiles.read(joined).energies.shape == (2, 11)


def test_matlab_grids_read_in_both_versions(tmp_path):
    from scipy.io import savemat

    from arpes_viewer.tools.crystal.fermisurface import sampled

    model = B.WannierModel(LATTICE, np.array(list(_hoppings()), float),
                           np.array(list(_hoppings().values())))
    n = 13
    axis = np.linspace(-0.5, 0.5, n)
    frac = np.stack(np.meshgrid(axis, axis, axis, indexing="ij"), -1).reshape(-1, 3)
    values = model.energies(frac @ model.reciprocal).reshape(1, n, n, n)
    ek = np.moveaxis(values, 0, -1)                              # MATLAB: k1 k2 k3 band
    k = np.random.default_rng(0).uniform(-1, 1, (20, 3))
    savemat(tmp_path / "Ek_old.mat", {"Ek_3D": ek})
    import h5py

    with h5py.File(tmp_path / "Ek_new.mat", "w", userblock_size=512) as handle:
        handle["Ek_3D"] = ek.transpose(3, 2, 1, 0)               # as MATLAB writes it
    with open(tmp_path / "Ek_new.mat", "r+b") as handle:
        handle.write(b"MATLAB 7.3 MAT-file, written by a test")
    for name in ("Ek_old.mat", "Ek_new.mat"):
        assert bandfiles.detect(str(tmp_path / name)) == "matlab"
        grid = bandfiles.read(str(tmp_path / name), lattice=LATTICE)
        assert np.abs(grid.energies(k) - model.energies(k)).max() < 0.02
    with pytest.raises(ValueError, match="cell is needed"):
        bandfiles.read(str(tmp_path / "Ek_old.mat"))
    assert sampled(grid) is grid


def test_big_reads_report_their_progress(tmp_path):
    hr = str(tmp_path / "x_hr.dat")
    bandfiles.write_hr(hr, LATTICE, _hoppings(), 1)
    seen = []
    bandfiles.read(hr, progress=lambda f, m: seen.append(f))
    assert seen and seen[-1] == pytest.approx(1.0)


# --------------------------------------------------------------------------
# the window
# --------------------------------------------------------------------------
def test_the_goal_decides_what_can_be_opened_and_shown(app, tmp_path):
    from arpes_viewer.ui.panels.bands import BandsDialog

    path = _path_file(str(tmp_path))
    dialog = BandsDialog([])
    dialog.set_goal("contours")
    filters = bandfiles.formats_for("contours")
    assert "*_band.dat" not in filters[0][1] and "*.bxsf" in filters[0][1]
    dialog.open_file(path)
    assert "cannot give constant-energy contours" in dialog._calculation_ready()
    dialog.set_goal("dispersion")
    assert dialog._calculation_ready() == ""
    assert dialog.stepper.steps[1].title == "Along a path"
    dialog.set_goal("fs")
    assert dialog.stepper.steps[1].title == "3D Fermi surface"


def test_the_fermi_surface_is_drawn_and_made_a_figure(app, tmp_path):
    from arpes_viewer.ui.panels.bands import BandsDialog

    hr = str(tmp_path / "wannier90_hr.dat")
    bandfiles.write_hr(hr, LATTICE, _hoppings(), 1)
    dialog = BandsDialog([])
    dialog.set_goal("fs")
    dialog.open_file(hr)
    dialog.fs_detail.setCurrentIndex(0)
    dialog.stepper.go(1)                         # entering it draws
    surface = dialog._fs_drawn
    assert surface is not None and len(surface.surfaces) == 1
    assert dialog.fs_table.rowCount() == 1
    assert dialog.panel_boxes["fs"].isChecked()
    panels, problems = dialog.panels()
    assert not problems and panels[0].data.kind == "picture"
    image = dialog.fs_view.render_image(300, 260)
    assert not image.isNull()
    window = dialog.open_fs_window()
    assert window.view.surface is surface
    assert os.path.exists(window.save_picture(str(tmp_path / "fs.png")))
    window.close()


def test_the_saved_path_keeps_its_high_symmetry_points(app, tmp_path):
    from arpes_viewer.ui.panels.bands import BandsDialog
    from arpes_viewer.ui.panels.figure import panel_from_dataset
    from arpes_viewer.ui.viewers.windows import CutWindow, path_points

    with open(tmp_path / "KPOINTS", "w") as handle:
        handle.write("2 6   Direct Lattice Vectors(Ang)   Reciprocal\n")
        for row in LATTICE:
            handle.write(" ".join(f"{v:.6f}" for v in (*row, 0, 0, 0)) + "\n")
        handle.write("0 0 0 ! Γ\n0.5 0 0 ! M\n\n0.5 0 0 ! M\n0.333 0.333 0 ! K\n")
    with open(tmp_path / "band.dat", "w") as handle:
        d = np.r_[np.linspace(0, 1, 6), np.linspace(1, 1.6, 6)]
        handle.write("".join(f"  {x:.4f}  {x - 0.8:.4f}\n" for x in d))
    dialog = BandsDialog([])
    dialog.open_file(str(tmp_path / "band.dat"))
    created = []
    dialog.datasetsCreated.connect(created.extend)
    dialog.spectrum_to_list()
    data = created[0]
    points = path_points(data.scan.info)
    assert [label for _d, label in points] == ["Γ", "M", "K"]
    window = CutWindow(data, "calc", "gray", False)
    top = window.view.plot.getAxis("top")
    assert window.points_box.isChecked() and len(window._point_lines) == 3
    window.points_box.setChecked(False)
    assert not window._point_lines
    window.close()
    assert [label for _d, label in panel_from_dataset(data).x_tick_labels] == ["Γ", "M", "K"]
    del top


def test_which_file_lists_the_folder_and_opens_the_choice(app, tmp_path):
    from arpes_viewer.ui.panels.bandfinder import WhichFileDialog

    _path_file(str(tmp_path))
    (tmp_path / "notes.txt").write_text("from the theory side\n")
    dialog = WhichFileDialog(str(tmp_path), "dispersion")
    names = [dialog.table.item(i, 0).text() for i in range(dialog.table.rowCount())]
    assert names == ["band.dat", "notes.txt"]
    assert "not a file this program reads" in dialog.table.item(1, 1).text()
    chosen = []
    dialog.chosen.connect(chosen.append)
    dialog.open_selected()                       # band.dat was selected for us
    assert chosen == [str(tmp_path / "band.dat")]
    dialog = WhichFileDialog(str(tmp_path), "fs")
    assert "No file here can show" in dialog.missing_label.text()


def test_a_step_window_asks_before_closing_and_keeps_its_settings(app, tmp_path,
                                                                  monkeypatch):
    from arpes_viewer.ui.common import stepper
    from arpes_viewer.ui.panels.bands import BandsDialog

    path = _path_file(str(tmp_path))
    dialog = BandsDialog([])
    dialog.show()
    dialog.open_file(path)
    dialog.set_goal("fs")
    dialog.keep_lo.setValue(-3.5)
    monkeypatch.delenv("ARPES_VIEWER_NO_CONFIRM")
    monkeypatch.setattr(stepper, "confirm_close", lambda *a: False)
    dialog.stepper.go(3)
    dialog.stepper.next_button.click()           # "Close": asks, and is told no
    assert dialog.isVisible()
    assert not dialog.stepper.next_button.isDefault()   # Enter cannot press it
    monkeypatch.setattr(stepper, "confirm_close", lambda *a: True)
    dialog.close()
    assert not dialog.isVisible()
    again = BandsDialog([])
    assert again.goal() == "fs" and again.keep_lo.value() == -3.5
    assert not again.reopen_button.isHidden() and "band.dat" in again.reopen_button.text()
    again.stepper.defaults_button.click()
    assert again.goal() == "dispersion" and again.keep_lo.value() == -2.0


# --------------------------------------------------------------------------
# round 3: the plane, the cell, Enter, the bands' looks, the goal's steps
# --------------------------------------------------------------------------
def _two_band_hr(folder):
    R = [(0, 0, 0), (1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (1, 1, 0),
         (-1, -1, 0), (0, 0, 1), (0, 0, -1)]
    H = np.zeros((9, 2, 2), complex)
    H[0] = np.diag([0.3, -0.4])
    H[1:7, 0, 0], H[1:7, 1, 1] = -0.15, 0.12
    H[7:, 0, 0], H[7:, 1, 1] = 0.08, -0.05
    hr = os.path.join(folder, "x_hr.dat")
    bandfiles.write_hr(hr, LATTICE, {r: H[i] for i, r in enumerate(R)}, 2)
    return hr


def test_the_contour_plane_is_set_like_the_zone_tool(app, tmp_path):
    from arpes_viewer.ui.panels.bands import BandsDialog

    dialog = BandsDialog([])
    dialog.open_file(_two_band_hr(str(tmp_path)))
    dialog.set_goal("contours")
    dialog.stepper.go(1)
    items = [dialog.cont_through.itemText(i) for i in range(dialog.cont_through.count())]
    assert any(text.startswith("A, H, L") for text in items)
    top = next(i for i, t in enumerate(items) if t.startswith("A, H, L"))
    dialog.cont_through.setCurrentIndex(top)
    dialog._through_chosen(top)
    assert dialog.cont_offset.value() == pytest.approx(np.pi / 6.0, abs=1e-4)
    assert "holds A, H, L" in dialog.cont_label.text()
    dialog.cont_window.setValue(0.2)                     # along the (hkl) normal
    overlays, axis, image = dialog.plane_overlays(0.0)
    assert image is not None and image.shape == (axis.size, axis.size)
    first = dialog.current_plane().e1
    dialog.cont_rotation.setValue(30.0)                  # turns the drawn axes
    turned = dialog.current_plane()
    assert np.dot(first, turned.e1) == pytest.approx(np.cos(np.pi / 6))
    assert np.allclose(np.cross(turned.e1, turned.e2), turned.normal)
    dialog.cont_h.setValue(1)
    dialog.cont_l.setValue(0)
    assert abs(dialog.current_plane().normal[2]) < 1e-9
    zone = dialog.open_zone_window()
    assert zone.view.planes and "(1 0 0)" in zone.label.text()
    zone.close()


def test_the_cell_is_locked_until_asked(app, tmp_path):
    from arpes_viewer.ui.panels.bands import BandsDialog

    dialog = BandsDialog([])
    dialog.open_file(_two_band_hr(str(tmp_path)))
    assert not dialog.space_group.isEnabled() and not dialog.take_button.isEnabled()
    assert not any(box.isEnabled() for box in dialog.cell_boxes)
    k = np.array([[0.4, 0.1, 0.2]])
    before = dialog.model.energies(k)
    dialog.change_lattice.setChecked(True)
    assert dialog.space_group.isEnabled() and dialog.cell_boxes[0].isEnabled()
    for box in dialog.cell_boxes[:3]:
        box.setValue(2 * box.value())
    dialog._cell_typed()
    assert np.allclose(dialog.model.energies(k / 2), before)


def test_enter_in_a_field_takes_the_value_and_presses_nothing(app, tmp_path):
    from PyQt5.QtCore import Qt
    from PyQt5.QtTest import QTest

    from arpes_viewer.ui.panels.bands import BandsDialog

    dialog = BandsDialog([])
    dialog.show()
    dialog.open_file(_path_file(str(tmp_path)))
    dialog.options_toggle.setChecked(True)
    dialog.keep_lo.setFocus()
    app.processEvents()
    QTest.keyClick(QApplication.focusWidget(), Qt.Key_Return)
    app.processEvents()
    assert dialog.stepper.index == 0                 # Next was not pressed
    dialog.stepper.next_button.setFocus()
    QTest.keyClick(dialog.stepper.next_button, Qt.Key_Return)
    assert dialog.stepper.index == 1                 # a focused button still is
    dialog.close()


def test_bands_take_their_own_colours_and_opacity(app, tmp_path):
    from arpes_viewer.ui.panels.bands import BandsDialog

    dialog = BandsDialog([])
    dialog.set_goal("fs")
    dialog.open_file(_two_band_hr(str(tmp_path)))
    dialog.fs_detail.setCurrentIndex(0)
    dialog.fs_energy.setValue(0.05)
    dialog._fs_energy_entered()                      # what Enter does
    surface = dialog._fs_drawn
    assert surface is not None and surface.energy == pytest.approx(0.05)
    bands = dialog.fs_table.bands
    for band in bands:                               # one group, one colour
        dialog.fs_table.set_colour(band, "#123456")
    dialog.fs_table.cellWidget(0, 2).setValue(0.4)
    assert all(dialog.fs_view.style_of(i, b)["colour"] == "#123456"
               for i, b in enumerate(bands))
    assert dialog.fs_view.style_of(0, bands[0])["opacity"] == pytest.approx(0.4)
    assert not dialog.fs_view.render_image(200, 180).isNull()
    # the contours take the same colours
    dialog.set_goal("contours")
    dialog.stepper.go(1)
    overlays, _axis, _image = dialog.plane_overlays(0.05)
    assert {o.color for o in overlays if o.name.startswith("band")} == {"#123456"}


def test_the_goal_shapes_the_steps_and_what_is_saved(app, tmp_path):
    from arpes_viewer.ui.panels.bands import BandsDialog

    dialog = BandsDialog([])
    dialog.open_file(_two_band_hr(str(tmp_path)))
    dialog.set_goal("fs")
    assert 2 in dialog.stepper.skipped
    assert not dialog.save_button.isEnabled()
    dialog.fs_detail.setCurrentIndex(0)
    dialog.stepper.go(1)
    dialog.stepper.next_button.click()
    assert dialog.stepper.index == 3                 # past Compare
    assert dialog.stepper.title.text().startswith("Step 3 of 3")
    dialog.set_goal("contours")
    assert 2 not in dialog.stepper.skipped
    assert dialog.save_button.isEnabled() and "map" in dialog.save_button.text()
    created = []
    dialog.datasetsCreated.connect(created.extend)
    dialog.map_lo.setValue(-0.2)
    dialog.map_hi.setValue(0.1)
    dialog.spectrum_to_list()
    data = created[0]
    assert data.kind == "k_map" and data.scan.value.shape[2] == 31
    dialog.set_goal("dispersion")
    dialog.spectrum_to_list()
    assert created[1].kind == "cut"
