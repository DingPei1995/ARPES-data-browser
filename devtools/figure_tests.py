"""The two example figures of the user guide's worked example (§30), rebuilt
from the demo data as a test of the figure tools.

The formats are those of a study of a nodal-line semimetal (HfSiS). Every
panel is made from the demo crystal of ``devtools/demo_data.py`` -- a
layered hexagonal crystal whose two bands cross along a line through the
zone -- and the demo tetragonal (ZrSiS-type) structure, with the program's
own tools:

* **Example figure 1** -- (a) the notched cube of a k-map (3-D view), (b)
  the constant-energy maps at four energies above the calculated contours
  at the same energies (Calculated bands), coloured by orbital weight, the
  zones dashed, (c) the constant-energy maps stacked up the energy axis
  (Stacked slices);
* **Example figure 2** -- (a) the layered crystal structure: glossy atoms,
  the cell's lattice constants, the cleavage plane; (b) the Brillouin zone
  with the measured constant-energy maps on its Γ and zone-top planes and
  the crossing line; (c) the cuts at four k_z stacked along k_z with the
  crossings joined (Stacked slices); (d) (i)-(iv) the four cuts with the
  calculated bands as dashed guides, (v) the MDC stack round the crossing
  (Stack plot), (vi) the fitted crossing energies against k_z with the
  calculated line.

Usage:  python devtools/figure_tests.py [PNG folder] [.afig folder]
        (PNGs default to docs/images)
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

from loader import registry                          # noqa: E402,F401

import demo_data as D                                # noqa: E402
import paper_figures as P                            # noqa: E402

PI_C = np.pi / D.LATTICE_C                            # k_z unit, π/c


def _style(figure, width, height):
    P._figure_style(figure)
    figure.style.width_mm, figure.style.height_mm = width, height


def _write(window, out, name, save_to):
    written = window.write(os.path.join(out, f"{name}.png"), "png")
    if save_to:
        window.save_file(os.path.join(save_to, f"{name}.afig"))
    return written


def _note(panel, text, x, y, colour="#000000", anchor="left", bold=False):
    from tools.output import figure as F
    panel.texts.append(F.TextNote(x=x, y=y, text=text, color=colour,
                                  anchor=anchor, bold=bold))


# --------------------------------------------------------------------------
# Test figure 1
# --------------------------------------------------------------------------
def notch_panel(data, name="Γ-plane map"):
    """(a): the notched cube, the corner in front cut away down to E_F."""
    from ui.panels.volume import VolumeWindow

    scan = data.scan
    labels = ("k_x (Å⁻¹)", "k_y (Å⁻¹)", "E - E_F (eV)")
    window = VolumeWindow(scan.value, (scan.x, scan.k, scan.z), labels, name,
                          colormap="blue", flip=True, mode="notch",
                          sampling=([(-1.25, 1.25), (-1.25, 1.25), (-1.3, 0.0)],
                                    [140, 140, 120]))
    for control, value in zip(window.notch_at, (0.0, 0.0, -1.3)):
        control.set_position(value)
    window.notch_corner.setCurrentIndex(
        [window.notch_corner.itemData(i) for i in range(8)].index((1, -1, 1)))
    window.azimuth.setValue(-52.0)
    window.elevation.setValue(24.0)
    window.z_aspect.setValue(1.1)
    for d in range(3):
        window.scene.styles[d].levels = window._auto_levels
    window.rebuild()
    panel = window.figure_panel()
    window.close()
    return panel


def stack_panel(data, energies, mark=None, name="Γ-plane map"):
    """(c): constant-energy maps stacked up the energy axis."""
    from ui.panels.volume import VolumeWindow

    scan = data.scan
    labels = ("k_x (Å⁻¹)", "k_y (Å⁻¹)", "E - E_F (eV)")
    window = VolumeWindow(scan.value, (scan.x, scan.k, scan.z), labels, name,
                          colormap="blue", flip=True, mode="stack",
                          sampling=([(-1.2, 1.2), (-1.2, 1.2), (-1.0, 0.05)],
                                    [150, 150, 120]))
    window.set_stack(2, energies, width=0.04)
    for d in range(3):
        window.scene.styles[d].fade = 0.45
    window.azimuth.setValue(-45.0)
    window.elevation.setValue(14.0)
    window.z_aspect.setValue(2.2)
    if mark is not None:
        window.set_marks(*mark)
    panel = window.figure_panel()
    window.close()
    return panel


def test_figure_1(folder, dialog, datasets, out, save_to=None):
    from tools.output import figure as F
    from ui.panels.figure import FigureWindow

    data = dict(datasets)
    plane = data["Γ-plane map"]
    pa = notch_panel(plane)
    pa.colorbar, pa.colorbar_ends, pa.colorbar_label = True, "high-low", ""
    pa.colorbar_position, pa.colorbar_length = "below", 0.6

    energies = [0.0, -0.3, -0.6, -0.9]
    dialog.compare_combo.setCurrentText("Γ-plane map")
    dialog.direction_box.setValue(0.0)
    dialog.integrate_box.setValue(0.05)
    dialog.weight_combo.setCurrentText("Se-pxy")       # coloured by its weight
    series = dialog.series_panels(plane, energies)
    dialog.weight_combo.setCurrentIndex(0)
    gm = 2 * np.pi / (np.sqrt(3) * D.LATTICE_A)
    gk = 4 * np.pi / (3 * D.LATTICE_A)
    for i, panel in enumerate(series):
        panel.x_range = panel.y_range = (-1.2, 1.2)
        panel.colormap, panel.flip = "blue", True
        _note(panel, f"({['i', 'ii', 'iii', 'iv', 'v', 'vi', 'vii', 'viii'][i]})",
              -1.12, 1.0, "#000000")
        panel.title = panel.title if i < 4 else ""
    top = series[0]
    _note(top, "Γ", 0.06, -0.12)
    _note(top, "M", 0.06, gm - 0.16)
    _note(top, "K", gk - 0.15, 0.08)
    for panel in series[4:]:
        panel.title = ""
    pc = stack_panel(plane, [0.0, -0.15, -0.3, -0.45, -0.6, -0.75, -0.9])

    window = FigureWindow([pa] + series + [pc], "Example figure 1")
    figure = window.figure
    _style(figure, 200, 80)
    size, gap = 20.0, 2.0
    frames = [(2, 13, 44, 50)]
    for i in range(8):
        row, col = divmod(i, 4)
        frames.append((66 + col * (size + gap), 12 + row * (size + 8.5), size, size))
    frames.append((154, 12, 44, 66))
    for panel, frame in zip(figure.panels, frames):
        panel.frame = frame
    # the composer's Arrange tools: (a), (b) and (c) start on one line
    window.arrange("top", [0, 1, 2, 3, 4, 9])
    labels = ["(a)", "(b)"] + [""] * 7 + ["(c)"]
    for panel, label in zip(figure.panels, labels):
        panel.label = label
    series[0].label_offset = (-6.0, 0.0)             # clear of its k_y axis
    for i, panel in enumerate(series):
        panel.x_labels = i >= 4
        panel.y_labels = i % 4 == 0
    pa.x_labels = pa.y_labels = True
    window.refresh()
    return _write(window, out, "example_figure_1", save_to), window


# --------------------------------------------------------------------------
# Test figure 2
# --------------------------------------------------------------------------
def sorted_curve(branch):
    """``(E, k)`` of a branch, sorted by energy (for interpolating k(E))."""
    k, e = (np.asarray(v, float) for v in branch)
    order = np.argsort(e)
    return e[order], k[order]


def node_branches(model, alignment, kz, k_node, span=0.25, step=0.05):
    """The calculation's two bands crossing at ``k_node`` on the plane
    ``kz`` (relative to Γ), along the cut, as aligned straight ``(k, E)``
    lines through the node with the bands' own slopes there. (Following
    the bands by index away from the node would hop onto other bands
    wherever they reorder.)"""
    from tools.crystal import bands as B

    at = alignment.apply(B.bands_along_cut(model, np.array([k_node - step, k_node]),
                                           90.0, kz=kz))
    node = at[:, 1]
    gaps = np.abs(node[:, None] - node[None, :])
    np.fill_diagonal(gaps, np.inf)
    i, j = np.unravel_index(np.argmin(gaps), gaps.shape)
    e_node = 0.5 * (node[i] + node[j])
    k = np.linspace(k_node - span, k_node + span, 201)
    lines = []
    for band in (i, j):
        slope = (node[band] - at[band, 0]) / step
        lines.append((k, e_node + slope * (k - k_node)))
    return tuple(lines)


def traced_bands(model, alignment, kz, k_node, k_range, n=241):
    """The two calculated bands crossing at ``k_node``, each followed by
    its character straight through the crossing (and past any band it
    only anticrosses), aligned: ``[(k, E), (k, E)]`` along the cut."""
    from tools.crystal import bands as B

    k = np.linspace(k_range[0], k_range[1], n)
    points = np.column_stack([np.zeros(n), k, np.full(n, kz)])   # along k_y
    start = min(int(np.argmin(np.abs(k - k_node))) + 3, n - 1)  # just past the node
    column = B.bands_along_cut(model, k[start:start + 1], 90.0, kz=kz)[:, 0]
    gaps = np.abs(column[:, None] - column[None, :])
    np.fill_diagonal(gaps, np.inf)
    pair = np.unravel_index(np.argmin(gaps), gaps.shape)
    out = []
    for b in pair:
        energies = alignment.apply(B.follow_band(model, points, start, int(b)))
        # where the trace jumps a small anticrossing, leave a gap rather
        # than draw the jump
        jump = np.abs(np.diff(energies)) > 4 * np.nanmedian(np.abs(np.diff(energies)))
        energies[1:][jump] = np.nan
        out.append((k, energies))
    return out


def kz_cut(kz_map, kz):
    """The cut through a k_z map at one k_z (its nearest k_z plane, as the
    3-D viewer's slice takes it)."""
    scan = kz_map.scan
    i = int(np.argmin(np.abs(np.asarray(scan.x) - kz)))
    return (np.asarray(scan.k, float), np.asarray(scan.z, float),
            np.asarray(scan.value, float)[i], float(scan.x[i]))


def fit_dirac_point(k, energy, values, branches, guess_e, half=0.14):
    """The two straight branches of the crossing, from MDC fits seeded on
    the calculation's two branches (``branches``: two ``(k, E)`` curves,
    aligned): ``(lines, (k, E, σk, σE))``."""
    from tools.spectra import dispersion as DISP
    from tools.spectra import peaks as PK

    window = (guess_e - half, guess_e + half)
    guess_k = float(np.mean([np.interp(guess_e, *sorted_curve(b)) for b in branches]))
    keep = (k > guess_k - 0.3) & (k < guess_k + 0.3)
    kk, vv = k[keep], values[keep]
    bands = []
    for name, branch in zip(("first", "second"), branches):
        band = PK.Band(name)
        for e in np.linspace(window[0], window[1], 4):
            band.seeds.append(PK.Seed(e, float(np.interp(e, *sorted_curve(branch))),
                                      0.05, 1.0))
        bands.append(band)
    series = PK.fit_series(vv, (kk, energy), bands, PK.FitSettings(),
                           direction="mdc", bounds=window)
    lines = []
    for index, branch in enumerate(branches):
        kb, eb, sk = series.band(index)
        # the crossing itself blurs, and a weak branch's peak can wander onto
        # a neighbour: keep the centres near the branch they were seeded on
        near = np.abs(kb - np.interp(eb, *sorted_curve(branch))) < 0.04
        keep = (np.abs(eb - guess_e) > 0.03) & near & np.isfinite(sk)
        lines.append(DISP.fit_dispersion(kb[keep], eb[keep], k_error=sk[keep], order=1))
    return lines, DISP.line_crossing(*lines)


def bz_picture(folder, datasets, size=900):
    """(a): the zone with the measured maps on its Γ and zone-top planes and
    the line node (K-H) in red."""
    from ui.panels.crystal import CrystalDialog
    from ui.panels.figure import _qimage_array

    maps = [(n, d) for n, d in datasets if d.kind == "k_map"]
    dialog = CrystalDialog(datasets=maps, colormap="blue", flip=True)
    dialog.read_cif(os.path.join(folder, "demo_crystal.cif"))
    dialog.stepper.go(3)
    dialog.path_edit.setText("K-H")
    dialog.zone_left_edit.setText("L, M")             # clear of the path and axes
    for (box, height, colour, source, energy), name, e in zip(
            dialog.plane_rows, ("Γ-plane map", "A-plane map"), (-0.65, -0.55)):
        source.setCurrentText(name)
        energy.setValue(e)
    dialog._zone_changed()
    dialog.view.camera.azimuth, dialog.view.camera.elevation = -28.0, 22.0
    dialog.view.font_px = 30
    image = P._trim(_qimage_array(dialog.view.render_image(size, size)))
    return image, os.path.join(folder, "demo_crystal.cif")


def kz_stack_panel(kz_map, kz_values, dps):
    """(c), right: the cuts at the four k_z, one behind the other along
    k_z, with the Dirac points joined."""
    from ui.panels.volume import VolumeWindow

    scan = kz_map.scan
    labels = ("k_z (Å⁻¹)", "k_∥ (Å⁻¹)", "E - E_F (eV)")
    window = VolumeWindow(scan.value, (scan.x, scan.k, scan.z), labels, "k_z",
                          colormap="blue", flip=True, mode="stack",
                          sampling=([(min(kz_values) - 0.03, max(kz_values) + 0.03),
                                     (0.1, 0.9), (-1.0, 0.05)], [60, 120, 120]))
    window.set_stack(0, kz_values, width=0.0)
    for d in range(3):
        window.scene.styles[d].fade = 0.3
    # each cut seen almost face on, the stack spread out along k_z so none
    # hides the next
    window.stack_spread.setValue(3.0)
    window.azimuth.setValue(40.0)
    window.elevation.setValue(12.0)
    window.z_aspect.setValue(1.0)
    window.set_marks([[kz, k, e] for kz, (k, e) in zip(kz_values, dps)],
                     "", "red")
    panel = window.figure_panel()
    window.close()
    return panel


def test_figure_2(folder, dialog, datasets, out, save_to=None):
    from tools.crystal import bands as B
    from tools.output import figure as F
    from tools.output.figfile import source_record
    from ui.panels.figure import FigureWindow, Panel, panel_from_arrays
    from ui.viewers.curves import curve_panel

    data = dict(datasets)
    kz_map = data["k_z scan"]
    gamma_kz = 8 * PI_C
    kz_values = gamma_kz + np.array([0.0, 0.3, 0.6, 1.0]) * PI_C

    # (a) the zone, the measured maps on its planes
    picture, cif = bz_picture(folder, datasets)
    pa = Panel(data=F.picture_data(picture, "Brillouin zone"))
    pa.x_labels = pa.y_labels = False
    pa.data.sources = (dialog.provenance([data["Γ-plane map"], data["A-plane map"]])
                       + [source_record(cif, kind="structure")])

    # (d) (i)-(iv): the four cuts, the branches fitted, the Dirac points
    model = dialog.model
    cuts, dps, fits = [], [], []
    measured_kz = [kz_cut(kz_map, kz)[3] for kz in kz_values]
    # one call through every k_z, starting on the Γ plane, so the same pair
    # of bands is followed (two other bands touch elsewhere on the cut)
    node_k, node_e, _gap = B.node_line(model, np.linspace(0.2, 0.8, 301),
                                       np.r_[0.0, np.array(measured_kz) - gamma_kz],
                                       90.0)
    node_k, node_e = node_k[1:], dialog.alignment.apply(node_e[1:])
    nk_cut, ne_cut = node_k, node_e
    for kz, nk, ne in zip(kz_values, node_k, node_e):
        k, energy, values, at = kz_cut(kz_map, kz)
        guess_e = float(ne)
        branches = node_branches(model, dialog.alignment, at - gamma_kz, float(nk))
        lines, crossing = fit_dirac_point(k, energy, values, branches, guess_e)
        cuts.append((k, energy, values, at))
        fits.append(lines)
        dps.append(crossing)
    panels_d = []
    for i, ((k, energy, values, at), lines, dp) in enumerate(zip(cuts, fits, dps)):
        panel = panel_from_arrays(values, k, energy, "k_∥ (Å⁻¹)", "E - E_F (eV)",
                                  f"k_z = {at / PI_C:.1f}π/c", colormap="blue",
                                  sources=dialog.provenance([kz_map]))
        panel.flip = True
        panel.x_range = (max(0.1, float(k.min())), min(0.9, float(k.max())))
        panel.y_range = (-1.0, 0.05)
        panel.title = f"k_z = {at / PI_C:.1f}π/c"
        # dashed guides on the two bands that cross, from the calculation
        for xs, ys in traced_bands(model, dialog.alignment, at - gamma_kz,
                                   float(nk_cut[i]), (0.15, 0.85)):
            panel.overlays.append(F.Overlay(xs=list(xs), ys=list(ys),
                                            color="#d62728", width_pt=0.8,
                                            style="dashed", name="band"))
        panel.lines.append(F.Line(value=0.0, axis="y", color="#d62728",
                                  width_pt=0.5, style="dashed"))
        _note(panel, f"({['i', 'ii', 'iii', 'iv'][i]})", 0.13, -0.1)
        if i == 0:
            box = (0.33, 0.67, ne_cut[0] - 0.15, ne_cut[0] + 0.15)
            panel.overlays.append(F.Overlay(
                xs=[box[0], box[1], box[1], box[0], box[0]],
                ys=[box[2], box[2], box[3], box[3], box[2]], color="#000000",
                width_pt=0.6, style="dashed", name="MDC region"))
        panels_d.append(panel)

    # (v) the MDC stack round the first crossing (the Stack plot window),
    # with the same dashed guides as the cut, through each MDC
    from ui.panels.process import StackWindow
    k, energy, values, at = cuts[0]
    half = 0.15
    keep_k = (k > 0.33) & (k < 0.67)
    keep_e = np.abs(energy - ne_cut[0]) <= half
    stack = StackWindow(values[keep_k][:, keep_e], (k[keep_k], energy[keep_e]),
                        ("k_∥ (Å⁻¹)", "E - E_F (eV)"), "MDCs",
                        sources=dialog.provenance([kz_map]))
    stack.direction.setCurrentIndex(1)                # MDCs
    stack.count.setValue(40)
    stack.combine.setValue(3)                         # each the sum of 3 rows
    stack.offset.setValue(0.32)
    stack.colour_mode.setCurrentIndex(0)
    stack.redraw()
    pv = stack.figure_panel()
    mdcs = list(pv.overlays)
    for overlay in mdcs:
        overlay.color, overlay.width_pt = "#000000", 0.45
    for k_band, e_band in traced_bands(model, dialog.alignment, at - gamma_kz,
                                       float(nk_cut[0]), (0.33, 0.67)):
        xs, ys = [], []
        for overlay in mdcs:
            e_curve = float(overlay.name)
            good = np.isfinite(e_band)
            order = np.argsort(e_band[good])
            kb = float(np.interp(e_curve, e_band[good][order], k_band[good][order],
                                 left=np.nan, right=np.nan))
            if np.isfinite(kb):
                xs.append(kb)
                ys.append(float(np.interp(kb, overlay.xs, overlay.ys)))
        pv.overlays.append(F.Overlay(xs=xs, ys=ys, color="#d62728", width_pt=0.8,
                                     style="dashed", name="band"))
    pv.y_label = "Intensity (arb. units)"
    pv.title = f"k_z = {cuts[0][3] / PI_C:.1f}π/c"
    stack.close()

    # (vi) the Dirac points against k_z, with the calculated line node
    kz_dense = np.linspace(kz_values[0] - 0.4 * PI_C, kz_values[-1] + 1.4 * PI_C, 121)
    _k, calc_e, _g = B.node_line(model, np.linspace(0.2, 0.8, 301),
                                 np.r_[0.0, kz_dense - gamma_kz], 90.0)
    calc_e = dialog.alignment.apply(calc_e[1:])
    measured = np.array([d[1] for d in dps])
    errors = [max(float(d[3]), 0.004) for d in dps]
    pw = curve_panel(kz_dense / PI_C, [("calculation", calc_e, "#d62728")],
                     x_label="k_z (π/c)", y_label="E_c - E_F (eV)", name="E_c")
    pw.overlays.append(F.Overlay(xs=list(np.array([c[3] for c in cuts]) / PI_C),
                                 ys=list(measured), color="#1f4e9c", width_pt=0.6,
                                 marker="o", marker_mm=1.1, line=False,
                                 y_errors=errors, name="measured"))
    pw.data.sources = dialog.provenance([kz_map])
    pw.y_range = (float(min(measured.min(), calc_e.min())) - 0.03,
                  float(max(measured.max(), calc_e.max())) + 0.03)
    pw.lines.append(F.Line(value=9.0, axis="x", color="#888888", width_pt=0.5,
                           style="dashed"))
    _note(pw, "(vi)", 7.75, pw.y_range[1] - 0.03)

    # (c), right: the cuts stacked along k_z, the line node through them
    pc = kz_stack_panel(kz_map, [c[3] for c in cuts], [(d[0], d[1]) for d in dps])

    pw.y_label = "E_c - E_F (eV)"                    # the crossing's energy
    crystal = crystal_panel(folder)
    window = FigureWindow([crystal, pa, pc] + panels_d + [pv, pw],
                          "Example figure 2")
    figure = window.figure
    _style(figure, 186, 112)
    w, h = 21.0, 34.0
    frames = [(4, 8, 34, 50), (46, 8, 48, 50), (96, 8, 88, 50)]
    for i in range(4):
        frames.append((12 + i * (w + 2.0), 70, w, h))
    frames.append((12 + 4 * (w + 2.0) + 12, 70, 22, h))
    frames.append((12 + 4 * (w + 2.0) + 50, 70, 26, h))
    for panel, frame in zip(figure.panels, frames):
        panel.frame = frame
    for panel, label in zip(figure.panels, ["(a)", "(b)", "(c)", "(d)", "", "", "",
                                            "", ""]):
        panel.label = label
    window.arrange("top", [0, 1, 2])                   # the first row
    window.arrange("height", [0, 1, 2])
    window.arrange("top", list(range(3, 9)))           # the second
    window.arrange("height", list(range(3, 9)))
    for i, panel in enumerate(panels_d):
        panel.y_labels = i == 0
    window.refresh()
    print("crossings: " + ", ".join(
        f"k_z = {c[3] / PI_C:.2f}π/c: E = {d[1]:.3f} ± {d[3]:.3f} eV"
        for c, d in zip(cuts, dps)))
    return _write(window, out, "example_figure_2", save_to), window


# --------------------------------------------------------------------------
# Figure 2 (a): the crystal structure
# --------------------------------------------------------------------------
def crystal_panel(folder, size=1000):
    """The layered structure as a picture panel: glossy atoms, the cell's
    lattice constants and the cleavage plane between the S layers."""
    from tools.output import figure as F
    from tools.output.figfile import source_record
    from ui.panels.crystal import CrystalDialog
    from ui.panels.figure import Panel, _qimage_array

    cif = os.path.join(folder, "demo_tetragonal_crystal.cif")
    dialog = CrystalDialog()
    dialog.read_cif(cif)
    dialog.stepper.go(1)
    for box, n in zip(dialog.repeat_boxes, (1, 1, 1)):
        box.setValue(n)
    dialog.scale_box.setValue(0.52)
    for element, colour in (("Hf", "#3fb4b0"), ("Si", "#1f3f9e"), ("S", "#d6df2c")):
        dialog.set_style(element, colour=colour, gloss=0.85, alpha=1.0)
    dialog.bond_grey.setChecked(True)
    dialog.lengths_box.setChecked(True)
    dialog.stepper.go(2)
    dialog.plane_box.setChecked(True)
    dialog.set_plane_height(0.5 * 8.0)                # between the S layers
    dialog.plane_colour.set_colour("#9a9a9a")
    dialog._plane_changed()
    dialog.view.camera.azimuth, dialog.view.camera.elevation = -68.0, 14.0
    dialog.view.font_px = 26
    dialog._look_changed()
    picture = P._trim(_qimage_array(dialog.view.render_image(int(size * 0.75), size)))
    panel = Panel(data=F.picture_data(picture, "crystal structure"))
    panel.x_labels = panel.y_labels = False
    panel.data.sources = [source_record(cif, kind="structure")]
    return panel


def main(out=None, save_to=None):
    from PyQt5.QtWidgets import QApplication

    app = QApplication.instance() or QApplication([])
    out = out or os.path.join(HERE, "docs", "images")
    os.makedirs(out, exist_ok=True)
    with tempfile.TemporaryDirectory() as folder:
        D.write_all(folder)
        datasets = P.converted(folder)
        dialog = P.bands_window(folder, datasets)
        written = []
        for make in (lambda: test_figure_1(folder, dialog, datasets, out, save_to),
                     lambda: test_figure_2(folder, dialog, datasets, out, save_to)):
            files, _window = make()
            written += files
    for path in written:
        print(path)
    return app


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else None,
         sys.argv[2] if len(sys.argv) > 2 else None)
