"""Paper-style figures from the demo data, made with the program's tools.

A worked example of the Calculation tab, end to end, and the source of the
figures in the user guide's "Calculated bands" section. From the demo
crystal (``devtools/demo_data.py``: a layered hexagonal crystal, its
"calculation" as a Wannier90 Hamiltonian, its CIF and an EDS spectrum):

* **Figure 1** -- (a) the crystal structure with its coordination
  octahedra, (b) the EDS spectrum with the elements' lines and the
  composition, (c) the 3-D Brillouin zone with its labelled points, the
  k_z = 0 plane and the zone-top plane and the path, (d) the measured
  intensity along Γ-M-K-Γ-A-L-H-A, assembled from the Γ-plane and
  zone-top maps and the photon-energy scan, with the calculation over it;
* **Figure 2** -- (a) a mirror-symmetrised constant-energy map with the
  calculated contours and the cut's position, (b, c) the cut in LH and LV
  polarisation with the fitted calculation, (d) the calculation alone
  along the same direction, (e) the cut's second derivative, with band
  labels and arrows.

Every step is a call to the program's own code -- the k conversions, the
symmetrisation, the second derivative, the Calculated-bands window, the
crystal window and the figure composer -- so the figures are what a user
gets by clicking through the same tools.

Usage:  python devtools/paper_figures.py [PNG folder] [.afig folder]
        (PNGs default to docs/images; with a second folder the figures are
        also saved there as .afig files, each with its "data for ..." folder)
"""
import os
import sys
import tempfile

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for folder in (HERE, os.path.join(HERE, "devtools")):
    if folder not in sys.path:
        sys.path.insert(0, folder)

from arpes_viewer.loader import registry                          # noqa: E402

import demo_data as D                                # noqa: E402


# --------------------------------------------------------------------------
# The measurement, converted to momentum as the program converts it
# --------------------------------------------------------------------------
def _scan(path, entry):
    """A measurement, read as the launcher reads it: its info names the
    file and entry it came from, which every dataset made from it keeps."""
    scan = registry.load(path, entry)
    scan.info.setdefault("origin.path", os.path.abspath(path))
    scan.info.setdefault("origin.entry", entry)
    return scan


def converted(folder):
    """The demo measurements in momentum: ``[(name, MemoryData)]`` --
    the Γ-plane and zone-top k-maps, the photon-energy scan as (k_z, k_par,
    E), and the cuts in k."""
    from arpes_viewer.tools.momentum.cutk import convert_cut
    from arpes_viewer.tools.momentum.kspace import convert_map
    from arpes_viewer.tools.momentum.kzconv import to_kz_cube
    from arpes_viewer.ui.viewers.widgets import MemoryData

    path = os.path.join(folder, "demo_bulk_bands.nxs")
    out = []
    for entry, label in ((f"FS_map_hv{D.HV_GAMMA:g}_Gamma_plane", "Γ-plane map"),
                         (f"FS_map_hv{D.HV_TOP:g}_A_plane", "A-plane map")):
        scan = _scan(path, entry)
        hv = float(scan.info["photon_energy_eV"])
        kx, ky, energy, cube = convert_map(scan.x, scan.k, scan.z, scan.value,
                                           n_kx=181, n_ky=181)
        energy = np.asarray(energy) - (hv - D.WORK_FUNCTION)
        out.append((label, MemoryData(
            "k_map", (kx, ky, energy), np.asarray(cube, dtype=np.float32),
            {"x": "k_x (Å⁻¹)", "k": "k_y (Å⁻¹)", "z": "E - E_F (eV)"},
            source_label=label, source_info=dict(scan.info))))
    scan = _scan(path, "hv_scan_GM_50-92eV")
    kz, kpar, energy, cube = to_kz_cube(scan.x, scan.k, scan.z, scan.value,
                                        inner_potential=D.INNER_POTENTIAL,
                                        work_function=D.WORK_FUNCTION,
                                        n_kz=181, n_kpar=121)[:4]
    out.append(("k_z scan", MemoryData(
        "kz_map_k", (kz, kpar, energy), np.asarray(cube, dtype=np.float32),
        {"x": "k_z (Å⁻¹)", "k": "k_∥ (Å⁻¹)", "z": "E - E_F (eV)"},
        source_label="k_z scan", source_info=dict(scan.info))))
    for entry in ("cut_GM_LH", "cut_GM_LV", "cut_GK_LH"):
        scan = _scan(path, entry)
        hv = float(scan.info["photon_energy_eV"])
        k, energy, values, _report = convert_cut(scan.x, scan.y, scan.value,
                                                 deflector_deg=0.0, n_k=301)
        energy = np.asarray(energy) - (hv - D.WORK_FUNCTION)
        out.append((entry, MemoryData(
            "cut", (k, energy), np.asarray(values, dtype=np.float32),
            {"x": "k (Å⁻¹)", "y": "E - E_F (eV)"}, source_label=entry,
            source_info=dict(scan.info))))
    return out



# --------------------------------------------------------------------------
# The windows, driven as a user drives them
# --------------------------------------------------------------------------
def bands_window(folder, datasets):
    """The Calculated-bands window with the demo calculation open, the
    k_z scan's direction set, and Z and the shift fitted on the Γ-M cut."""
    from arpes_viewer.ui.panels.bands import BandsDialog

    dialog = BandsDialog(datasets, colormap="blue")
    dialog.open_file(os.path.join(folder, "demo_calc_hr.dat"))
    dialog.set_source("k_z scan", value=90.0)    # the hv scan ran along Γ-M
    dialog.e_lo.setValue(-1.4)
    dialog.e_hi.setValue(0.1)
    dialog.stepper.go(2)
    dialog.compare_combo.setCurrentText("cut_GM_LH")
    dialog.direction_box.setValue(90.0)
    alignment = dialog.fit_alignment()
    print(f"fitted: Z = {alignment.z:.3f} ± {alignment.z_err:.3f}, shift = "
          f"{alignment.shift * 1000:+.1f} ± {alignment.shift_err * 1000:.1f} meV "
          f"(the demo crystal: Z = {D.CALC_Z}, shift = {D.CALC_SHIFT * 1000:+.0f} meV)")
    return dialog


def crystal_pictures(folder, size=900):
    """The crystal with its octahedra, and its Brillouin zone, from the
    crystal window: two RGBA pictures."""
    from arpes_viewer.tools.crystal.bands import parse_path               # noqa: F401
    from arpes_viewer.ui.panels.crystal import CrystalDialog
    from arpes_viewer.ui.panels.figure import _qimage_array

    dialog = CrystalDialog()
    dialog.read_cif(os.path.join(folder, "demo_crystal.cif"))
    dialog.stepper.go(1)
    for box, n in zip(dialog.repeat_boxes, (2, 2, 2)):
        box.setValue(n)
    dialog.poly_combo.setCurrentText("Ti")
    dialog.poly_colour.set_colour("#8fc3e8")
    dialog.scale_box.setValue(0.5)
    dialog.view.camera.azimuth, dialog.view.camera.elevation = -62.0, 12.0
    crystal = _trim(_qimage_array(dialog.view.render_image(size, int(size * 1.25))))
    dialog.stepper.go(3)
    dialog.path_edit.setText("Γ-M-K-Γ-A-L-H-A")
    dialog._zone_changed()
    dialog.view.camera.azimuth, dialog.view.camera.elevation = -32.0, 16.0
    dialog.view.font_px = 30
    zone = _trim(_qimage_array(dialog.view.render_image(size, size)))
    return crystal, zone


def _trim(rgba, pad=8):
    """A picture cut down to what is drawn on it (not the white round it)."""
    ink = np.any(rgba[..., :3] < 245, axis=2) & (rgba[..., 3] > 0)
    rows, cols = np.flatnonzero(ink.any(axis=1)), np.flatnonzero(ink.any(axis=0))
    if not rows.size:
        return rgba
    r0, r1 = max(rows[0] - pad, 0), min(rows[-1] + pad + 1, rgba.shape[0])
    c0, c1 = max(cols[0] - pad, 0), min(cols[-1] + pad + 1, rgba.shape[1])
    return np.ascontiguousarray(rgba[r0:r1, c0:c1])


def eds_panel(folder):
    from arpes_viewer.tools.output import figure as F
    from arpes_viewer.ui.panels.figure import panel_from_dataset
    from arpes_viewer.ui.viewers.widgets import MemoryData

    scan = _scan(os.path.join(folder, "demo_eds.nxs"), "EDS_spectrum")
    data = MemoryData("edc", (scan.x, np.zeros(1)), scan.value, dict(scan.labels),
                      source_label="EDS", source_info=dict(scan.info))
    panel = panel_from_dataset(data, "EDS")
    for overlay in panel.overlays:
        overlay.color, overlay.width_pt = "#1f4e9c", 0.6
    counts = np.asarray(scan.value, float).ravel()
    energy = np.asarray(scan.x, float)
    for name, centre, _height in D.EDS_LINES:
        if name.endswith("β") or name.startswith(("C", "O")):
            continue
        top = float(counts[np.abs(energy - centre) < 0.08].max())
        panel.texts.append(F.TextNote(x=centre + 0.25, y=top * 1.04,
                                      text=name, color="#000000", anchor="left"))
    panel.texts.append(F.TextNote(x=13.6, y=float(counts.max()) * 0.62,
                                  text="Ti : Se = 1 : 2", color="#000000",
                                  anchor="right"))
    panel.x_range = (0.0, 14.0)
    panel.y_range = (0.0, float(counts.max()) * 1.15)
    panel.y_label = "Intensity (counts)"
    return panel


def _place(figure, frames, width, height):
    for panel, frame in zip(figure.panels, frames):
        panel.frame = frame
    figure.style.width_mm, figure.style.height_mm = width, height
    figure.letter_panels()


def _figure_style(figure):
    style = figure.style
    style.font_pt = 7.0
    style.panel_label_inside = False
    style.panel_label_color = "#000000"
    style.tick_direction = "in"
    style.dpi = 300


def figure_1(folder, dialog, out, save_to=None):
    """Crystal, EDS, Brillouin zone, and ARPES along the path."""
    from arpes_viewer.tools.output import figure as F
    from arpes_viewer.ui.panels.figure import FigureWindow, Panel

    crystal, zone = crystal_pictures(folder)
    from arpes_viewer.tools.output.figfile import source_record

    pa = Panel(data=F.picture_data(crystal, "crystal"))
    pc = Panel(data=F.picture_data(zone, "Brillouin zone"))
    cif = os.path.join(folder, "demo_crystal.cif")
    pa.data.sources = [source_record(cif, kind="structure",
                                     note="the picture is drawn from this CIF")]
    pc.data.sources = [source_record(cif, kind="structure",
                                     note="the zone is that of this crystal")]
    for panel in (pa, pc):
        panel.x_labels = panel.y_labels = False
    pb = eds_panel(folder)
    dialog.weight_combo.setCurrentIndex(0)
    pd = dialog.panel_path(colormap="blue")
    pd.flip = True
    pd.colorbar = True
    pd.colorbar_ends = "high-low"
    window = FigureWindow([pa, pb, pc, pd], "Figure 1 — demo crystal")
    figure = window.figure
    _figure_style(figure)
    _place(figure, [(4, 6, 40, 50), (56, 6, 58, 42), (124, 6, 52, 48),
                    (12, 68, 150, 62)], 180, 142)
    window.refresh()
    written = window.write(os.path.join(out, "paper_figure_1.png"), "png")
    if save_to:
        window.save_file(os.path.join(save_to, "paper_figure_1.afig"))
    return written, window


def _second_derivative(values, k, energy):
    from arpes_viewer.tools.processing.process import derivative, gaussian_smooth

    smooth = gaussian_smooth(values, (k, energy), (0.012, 0.012))
    d2 = derivative(smooth, (k, energy), 1, 2, window=0.06)
    return np.clip(d2, 0, None)


def figure_2(folder, dialog, datasets, out, save_to=None):
    """Constant-energy map, LH and LV cuts with the calculation, the
    calculation alone, and the second derivative."""
    from arpes_viewer.tools.output import figure as F
    from arpes_viewer.tools.processing.process import symmetrise
    from arpes_viewer.ui.panels.figure import FigureWindow, panel_from_arrays
    from arpes_viewer.ui.viewers.curves import curve_panel

    data = dict(datasets)
    # (a) the Fermi surface on the Γ plane, symmetrised, 50 meV
    dialog.compare_combo.setCurrentText("Γ-plane map")
    dialog.integrate_box.setValue(0.05)
    dialog.direction_box.setValue(0.0)
    kx, ky, plane, energy = dialog.map_plane(data["Γ-plane map"], energy=0.0)
    sym = symmetrise(plane, kx, ky, fold=6, mirror_angles=(0.0, 30.0))
    pa = panel_from_arrays(sym.values, kx, ky, "k_x (Å⁻¹)", "k_y (Å⁻¹)",
                           "E_F", colormap="blue",
                           sources=dialog.provenance([data["Γ-plane map"]]))
    for band, xs, ys in dialog.map_contours(kx, ky, energy):
        pa.overlays.append(F.Overlay(xs=list(xs), ys=list(ys), color="#d62728",
                                     width_pt=0.6, style="dashed",
                                     name=f"calc {band + 1}"))
    gm = 2 * np.pi / (np.sqrt(3) * D.LATTICE_A)
    pa.overlays.append(F.Overlay(xs=[0.0, 0.0], ys=[-1.2 * gm, 1.2 * gm],
                                 color="#000000", width_pt=0.9, name="cut"))
    for text, x, y in (("Γ", 0.08, -0.08), ("M", 0.08, gm), ("K", 4 * np.pi / (3 * D.LATTICE_A), 0.08)):
        pa.texts.append(F.TextNote(x=x, y=y, text=text, color="#000000", anchor="left"))
    lim = 1.25
    pa.x_range = pa.y_range = (-lim, lim)
    # (b, c) the Γ-M cut in LH and LV, the fitted calculation over each
    dialog.direction_box.setValue(90.0)
    panels = []
    for name, title in (("cut_GM_LH", "LH"), ("cut_GM_LV", "LV")):
        panel = dialog.panel_compare(data[name], colormap="blue",
                                     colour="#d62728")
        panel.title = title
        panel.x_range = (-1.2, 1.2)
        panel.y_range = (-1.3, 0.1)
        panels.append(panel)
    pb, pc = panels
    # (d) the calculation alone along M-Γ-M
    k = np.linspace(-1.2, 1.2, 481)
    bands = dialog.cut_bands(k)
    pd = curve_panel(k, [], x_label="", y_label="E - E_F (eV)", name="calculation")
    pd.data.sources = dialog.provenance([])
    pd.data.x = np.array([-1.2, 1.2])
    pd.data.y = np.array([-1.3, 0.1])
    pd.y_range = (-1.3, 0.1)
    for band, row in enumerate(bands):
        pd.overlays.append(F.Overlay(xs=list(k), ys=list(row), color="#000000",
                                     width_pt=0.9, name=f"band {band + 1}"))
    pd.x_tick_labels = [(-gm, "M"), (0.0, "Γ"), (gm, "M")]
    pd.lines.append(F.Line(value=0.0, axis="y", color="#888888", width_pt=0.5,
                           style="dashed"))
    pd.title = "Calc."
    # (e) the LH cut's second derivative, with the bands named
    kk, ee, values = dialog.cut_arrays(data["cut_GM_LH"])
    pe = panel_from_arrays(_second_derivative(values, kk, ee), kk, ee, "k (Å⁻¹)",
                           "E - E_F (eV)", "LH 2nd derivative", colormap="blue",
                           sources=dialog.provenance([data["cut_GM_LH"]]))
    pe.x_range, pe.y_range = (-1.2, 1.2), (-1.3, 0.1)
    pe.title = "LH, d²I/dE²"
    calc_at_0 = dialog.cut_bands(np.array([0.0, 0.45]))
    for label, band, x in (("α", 1, 0.45), ("β", 0, 0.0), ("γ", 3, 0.0)):
        y = float(calc_at_0[band][0 if x == 0.0 else 1])
        pe.arrows.append(F.Arrow(x0=x + 0.35, y0=y + 0.12, x1=x + 0.06, y1=y + 0.02,
                                 color="#d62728", width_pt=0.8))
        pe.texts.append(F.TextNote(x=x + 0.38, y=y + 0.14, text=label,
                                   color="#d62728", anchor="left"))
    for panel in (pb, pc, pe):
        panel.colorbar = False
    pe.colorbar = True
    pe.colorbar_ends = "high-low"
    pe.gamma = 0.6
    for panel in (pa, pb, pc, pe):
        panel.flip = True
    window = FigureWindow([pa, pb, pc, pd, pe], "Figure 2 — demo crystal")
    figure = window.figure
    _figure_style(figure)
    h = 46
    _place(figure, [(12, 8, 40, 40), (66, 8, 26, h), (96, 8, 26, h),
                    (126, 8, 18, h), (148, 8, 26, h)], 186, 64)
    for panel in (pc, pd, pe):
        panel.y_labels = False
    window.refresh()
    written = window.write(os.path.join(out, "paper_figure_2.png"), "png")
    if save_to:
        window.save_file(os.path.join(save_to, "paper_figure_2.afig"))
    return written, window


def main(out=None, save_to=None):
    """The two figures as PNGs in ``out``; with ``save_to``, also as .afig
    files there, each with its "data for ..." folder."""
    from PyQt5.QtWidgets import QApplication

    app = QApplication.instance() or QApplication([])
    out = out or os.path.join(HERE, "docs", "images")
    os.makedirs(out, exist_ok=True)
    with tempfile.TemporaryDirectory() as folder:
        D.write_all(folder)
        datasets = converted(folder)
        dialog = bands_window(folder, datasets)
        written, _w1 = figure_1(folder, dialog, out, save_to)
        written2, _w2 = figure_2(folder, dialog, datasets, out, save_to)
    for path in written + written2:
        print(path)
    return app


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else None,
         sys.argv[2] if len(sys.argv) > 2 else None)
