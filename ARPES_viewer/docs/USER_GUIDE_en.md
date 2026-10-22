# ARPES viewer — User Guide

*[中文版 / Chinese version](USER_GUIDE_zh.md) · [Version française](USER_GUIDE_fr.md)*

This guide walks through the program's main features with screenshots, in
the order you would normally use them: loading data, looking at it,
converting it to momentum, analysing it, and making figures. The technical
background for each step (formulas, test results, design choices) is in the
[README](../README.md) and [CHANGELOG](../CHANGELOG.md).

Every screenshot was taken from the running program, using the synthetic
data written by `devtools/demo_data.py`, so you can follow along without a
beamtime (see [Trying it with demo data](#2-trying-it-with-demo-data)).

---

## Contents

1. [Installing and starting](#1-installing-and-starting)
2. [Trying it with demo data](#2-trying-it-with-demo-data)
3. [The main window](#3-the-main-window)
4. [Loading data](#4-loading-data)
5. [What the viewers have in common](#5-what-the-viewers-have-in-common)
6. [Maps: the contour and its two cuts](#6-maps-the-contour-and-its-two-cuts)
7. [Cuts: from raw spectrum to band parameters](#7-cuts-from-raw-spectrum-to-band-parameters)
8. [From angles to momentum for a map, and the Brillouin zone](#8-from-angles-to-momentum-for-a-map-and-the-brillouin-zone)
9. [Arbitrary cuts and Fermi-surface correction](#9-arbitrary-cuts-and-fermi-surface-correction)
10. [Processing: smoothing, derivatives, curvature](#10-processing-smoothing-derivatives-curvature)
11. [Stack plots and the 3-D view](#11-stack-plots-and-the-3-d-view)
12. [Cut arithmetic: dichroism and reference division](#12-cut-arithmetic-dichroism-and-reference-division)
13. [Photon-energy scans (kz)](#13-photon-energy-scans-kz)
14. [Removing the detector grid](#14-removing-the-detector-grid)
15. [Curves: EDCs, MDCs and spin-resolved data](#15-curves-edcs-mdcs-and-spin-resolved-data)
16. [Spatially resolved scans (nano-ARPES)](#16-spatially-resolved-scans-nano-arpes)
17. [Figures and export](#17-figures-and-export)
18. [Autosave, memory and recovery](#18-autosave-memory-and-recovery)
19. [Tips and troubleshooting](#19-tips-and-troubleshooting)

---

## 1. Installing and starting

You need Python 3.9 or newer and five packages:

```bash
pip install PyQt5 pyqtgraph numpy h5py scipy
```

Then, from the `ARPES_viewer` folder:

```bash
python ARPES_viewer.py
```

The program reads:

| Beamline / source | Files |
|---|---|
| SOLEIL **ANTARES** | NeXus/HDF5 `*.nxs` (cuts, deflector maps, spatial scans) |
| SOLEIL **CASSIOPEE** | Scienta SES text files — single spectra, and numbered folders assembled into a map or a photon-energy scan |
| CASSIOPEE **spin end station** (MBS A-1) | `*.krx` and `*.txt` — cuts, deflector maps and spin-resolved EDCs |
| This program | its own saved `*.nxs` format |

## 2. Trying it with demo data

```bash
python devtools/demo_data.py            # writes ./demo_data/*.nxs
```

This writes six files of synthetic data: a Fermi-surface map, a pair of cuts
taken with LH and LV polarisation, a photon-energy scan, a spin-resolved EDC,
a small real-space scan and a map with a detector mesh printed on it. They
come from a toy band structure (a hexagonal electron pocket and a hole band
at Γ, hole pockets at K) with a Fermi cut-off, the analyser resolution and
Poisson noise, and open exactly like measured data.

To re-take every screenshot in this guide after changing a window:

```bash
python devtools/guide_screenshots.py    # writes docs/images/*.png
```

---

## 3. The main window

![The main window](images/main_window_info.png)

The main window is a **browser**: one row per *dataset*, not per file (one
`.nxs` file often holds several measurements, and each gets its own row). The
kind is written after the name: `[Cut]`, `[Map]`, `[kz map]`, `[k-map]`,
`[SPEM]`, `[EDC]`, `[Spin EDC]` …

| Action | How |
|---|---|
| Open a dataset in its own window | **double-click** the row, or **Open selected** |
| Everything else that works on datasets | **right-click** the row (or a selection of rows) |
| Read the axes: first/last value, points, step | right-click → **Show information** (fills *Data Information*) |
| Truncate / normalise / compress several datasets at once | **Data operations...** |
| Smooth, differentiate, curvature, backgrounds, symmetrise | **Process...** |
| See and free memory | the **Memory** box at the bottom |

The right-click menu only offers what fits the selection. A greyed entry says
why when you hover over it — for example, *Cut arithmetic* needs exactly two
cuts, the *3D view* needs a map.

**Show information** also opens the full metadata of the dataset in a window
of its own, with a filter box and **Export metadata (.csv)**:

![Dataset information](images/information_window.png)

**Clear list** asks first, and offers to save anything computed in this
session that has not been saved yet. Files on disk are never deleted.

---

## 4. Loading data

**Load data...** opens the loader. **Quick add...** skips it and adds the
chosen files with everything detected automatically.

![The loader](images/loader_dialog.png)

1. **Files** — pick one or many (*Select file(s)...*). The file filter is built
   from what the installed readers can open.
2. **Reader** — detected automatically, and the result is shown
   (*Detected: 3 × This program's own format*). Override it only if a file is
   being misread.
3. **Datasets found** — what is inside each file, read from its *structure*
   only. This is instant even for a multi-GB scan; nothing is read until you
   open a row.
4. **Axes** — *Axis order* reorders the array's dimensions once, at the door.
   *First axis of a map is* says what was scanned: a deflector angle (the
   default), a **photon energy**, a temperature, a gate voltage, a delay… This
   keeps the right unit everywhere afterwards, and stops a k conversion being
   offered for an axis that is not an angle.

**CASSIOPEE folders.** Select any file of a numbered series
(`name_1_ROI1_.txt`, `name_2_ROI1_.txt`, …) and the whole folder is listed
once, as one assembled cube. Whether it becomes an angle map or a kz map is
read from the `_i` parameter files next to the spectra.

### Data kinds

| Kind | Axes | What it is |
|---|---|---|
| `cut` | angle (or k) × E | one spectrum |
| `map` | angle × angle × E | a deflector or polar-angle scan |
| `k_map` | kx × ky × E | a map converted to momentum |
| `kz_map` | hν × angle × E | a photon-energy scan |
| `kz_map_k` | k_z × k∥ × E | that scan converted to momentum |
| `spem_4d`, `spem_1d` | x (× y) × angle × E | real-space scans |
| `edc`, `mdc`, `spin_edc` | E (or k) × channels | curves, with optional σ channels |

---

## 5. What the viewers have in common

Each dataset opens in its own window, so you can put a map, its cuts and
another file side by side.

![A cut viewer with the readout cursor](images/cut_viewer.png)

- **The top row**: the **Functions** menu (everything the window can do), the
  colormap (per window) with **Flip**, and the view controls — axis ranges,
  **Auto**, **Inv** (invert), **Reset**, **Grid**.
- **Above each image**: the **Levels** bar (drag the two handles, or type
  **Min%/Max%** percentiles; **Clip** applies them) and **Gamma**.
- **Under each image**: **Smooth** (display only, does not change the data)
  and **Lock x/y** (aspect ratio).
- **Readout cursor** — right-click the image → *Readout cursor*. It reports
  the value under it and, on a cut, draws the **EDC** (right) and **MDC**
  (below) through that point. **EDC ±** / **MDC ±** integrate over a window,
  shown as a shaded band on the image. Type a position into the **Cursor**
  row and press Enter to put it exactly there.
- **EDC → list** / **MDC → list** put the curve on screen in the main list as
  a dataset of its own (see [§15](#15-curves-edcs-mdcs-and-spin-resolved-data)).
- **Right-click the image → Export** saves what the panel shows — data
  (`.nxs`), an image, or vector axes (`.eps`); see [§17](#17-figures-and-export).
- **Functions → Slice → Save slice to the main list / Open slice in a new
  panel** keeps a copy of the slice on screen.

The **Functions** menu is organised in sections:

| Viewer | Analysis | Data operations | Visualization |
|---|---|---|---|
| map / kz map / k-map | | Arbitrary cut, Map k conversion, kz map processing, kz → momentum, De-grid map | Brillouin zone (k-maps), Slice figure |
| cut | Fermi level, MDC / EDC fit | FS correction, Cut k conversion, Cut arithmetic, De-grid | |
| curve | Curve fit, Spin analysis | Crop, Bin, Normalise, Background, Shift axis, Add √N errors | As a figure |

An entry that does not apply is not shown; one that applies but cannot run
yet is greyed, with the reason in its tooltip.

---

## 6. Maps: the contour and its two cuts

A map opens on its **constant-energy contour**. The slider at the top picks
the energy (**Pos** in eV, or **Ind** as an index; **±** integrates over a
window, **mean** or **sum**). **Deflector cut** and **Slit cut** open the two
orthogonal cuts, each in its own window.

![Constant-energy contour of a map](images/map_contour.png)

| Slit cut (slit angle vs energy) | Deflector cut (deflector angle vs energy) |
|---|---|
| ![Slit cut](images/map_slit_cut.png) | ![Deflector cut](images/map_deflector_cut.png) |

**One cursor for all three windows.** Switching the readout cursor on in any
of them switches it on in the others, and moving it anywhere re-slices the
others through the new point. The axes keep one colour each everywhere:
**deflector red, slit green, energy blue**. A cursor line has the colour of
the axis it holds constant, and the EDC or MDC it feeds is drawn in the same
colour. Integration widths (the ± boxes) belong to each window and are never
shared.

From the list you can also open a map's cut straight away (right-click →
**Open slit cut** / **Open deflector cut**); the contour then stays hidden
until you ask for it (*Functions → Show the map*).

---

## 7. Cuts: from raw spectrum to band parameters

The usual order for a single cut is: **Fermi level → offset the energy axis →
k conversion → MDC/EDC fit**.

### 7.1 Fermi level

*Functions → Fermi level...* fits a Fermi edge (Fermi–Dirac × linear DOS,
convolved with a Gaussian resolution) to the EDC summed over a box.

![Fermi level dialog](images/fermi_dialog.png)

1. Press **Show the selection box** and drag the box on the cut over the edge
   (or type the angle and energy ranges). The angle range says which channels
   are summed; the energy range is the fit window.

   ![The fit box on the cut](images/fermi_cut_box.png)

2. **Estimate**, then **Fit**. Tick **Hold** to keep a parameter fixed
   (the temperature is held at the recorded value by default). Every
   parameter is reported with its uncertainty.
3. **Offset energy axis** lists a copy of the cut on **E − E_F**.

The *Advanced* corner divides the Fermi cut-off out (for gaps and reference
samples) or fits E_F channel by channel.

### 7.2 Cut k conversion

*Functions → Cut k conversion...* converts the angle axis to k∥ (Å⁻¹). A cut
is a line in the (kx, ky) plane that generally misses Γ, so the dialog asks
where Γ is: type the angles, or **From a converted map...** to inherit them
from a map already converted.

![Cut k conversion](images/cut_k_dialog.png)

If the cut is already on E − E_F, the dialog fills in **Energy offset** with
the E_F that was subtracted (the note under it says so): the conversion uses
true kinetic energies, and the result stays on E − E_F.

![The cut in momentum](images/cut_in_k.png)

### 7.3 MDC / EDC fit

*Functions → MDC / EDC fit...* fits peaks line by line and reads the band
off the fitted centres. It only accepts a cut whose axis is already in Å⁻¹ —
a velocity in eV·deg would be meaningless.

![MDC fit panel](images/mdc_fit_panel.png)

1. **Bands** tab: **Add band** for each branch; select a band, press **Place
   a seed** and click the peak on the curve (right). Four or five seeds down a
   band are enough; **Line at** moves the line being looked at.
2. **Try this line** shows the fit of the current line (dashed: the guess).
3. **Fit settings** tab: MDC or EDC, peak shape (Lorentzian, Gaussian,
   Voigt), background, step, then **Fit the series**. Fitted positions are
   drawn on the image, sized by amplitude.
4. **Dispersion...** fits the band: **v_F** (straight line), **m\*** (parabola)
   — and checks how far out the fit is stable. **Self-energy...** gives Re Σ
   and Im Σ.

![Dispersion window](images/dispersion_window.png)

**Save table...** writes every fitted line; **Band to list** puts a band in
the main list.

### 7.4 Core-level check

A wide scan of the core levels is listed as a cut too. Right-click it in the
main list → **Core-level check...** to find out what its peaks are.

1. Type the sample's elements — `NbSCl` or `Nb, S, Cl` — and press **Next**.
   *Match within* (0.5 eV) is how close a peak must be to a tabulated level
   to count as the sample's.

   ![Elements](images/core_level_elements.png)

2. Check the energy reference. E_F is fitted when the window holds the
   Fermi edge, otherwise put at hν − W. Change the photon energy, the work
   function or E_F if they are wrong (the green line can be dragged); a
   hν − E_F far from 4–5 eV means the photon energy is off. Press
   **Analyse**.

   ![Energy reference](images/core_level_reference.png)

3. The results: the EDC on a binding-energy axis, each peak marked —
   green *sample* (orbital named), orange *uncertain* (a level of the
   sample within 2 eV: chemical shift or charging — your call), red *not
   the sample* (common elements that could explain it: contamination,
   higher-order light, Auger lines). Click a row for the full list of
   candidates. Type your own assignment in the *Assignment* column.

   ![Core-level check](images/core_level_check.png)

If every peak is off by the same amount, the notes say so and **Apply
suggested shift** corrects it (a charged sample, or a wrong E_F). **Add EDC
to list** lists the EDC on the binding-energy scale; **Copy report** copies
the whole result as text. Measuring at a second photon energy settles what
stays doubtful: core levels stay put in binding energy, Auger lines and
higher-order light move.

---

## 8. From angles to momentum for a map, and the Brillouin zone

### 8.1 Map k conversion

*Functions → Map k conversion...* on the contour window. Drag the readout
cursor to the point that should become **k = (0, 0)** — the dialog follows
it — then **Convert**. The result is a new `k-map` row in the list.

![Map k conversion](images/map_k_dialog.png)

- **Sample rotation** turns the map about the surface normal; **Set rotation
  from contour...** reads it off a direction you click on the contour (two
  points along a mirror line, say).
- **Output grid**: number of points, or resolution (Å⁻¹, meV).
- **Energy offset** is added to the energy axis before converting, for an
  axis that is not a true kinetic energy (filled in automatically for a map on
  E − E_F).

![The converted k-map](images/k_map_contour.png)

### 8.2 Brillouin zone

On a k-map: *Functions → Visualization → Brillouin zone...*.

![Brillouin-zone dialog](images/brillouin_dialog.png)

1. **Space group** first (1–230); the lattice parameters it fixes are greyed
   and follow the free ones. Or switch to **2D layer** for a surface lattice.
2. **Cut plane**: the surface normal (hkl) and an offset along it — the cut of
   the 3-D zone by that plane is drawn, for any surface, not only (001).
   **Through point** lists the offsets at which the plane passes through a
   high-symmetry point (Γ, Z, L, F ... for a rhombohedral crystal); choose
   one rather than working the number out. **Show cut plane in 3D...**
   previews it, with the named points marked and the ones in the plane
   ringed; **Label high-symmetry points** writes them on the contour, which
   is the quickest way to set the azimuth. The dialog opens on the last
   crystal you used, so a second map of the same sample needs no retyping.
3. **Conventional** or **Irreducible** zone; **Azimuth** (or **Pick
   direction from contour...**); **Repeat across the field of view**.
4. **Moiré** adds a second, twisted or mismatched layer and its moiré zone.
5. **Plot**.

![k-map with the zone overlaid](images/k_map_with_bz.png)

---

## 9. Arbitrary cuts and Fermi-surface correction

### 9.1 Arbitrary cut

*Functions → Arbitrary cut...* on any contour: click up to six points (or type
them) and the cut runs through them in order — here Γ → M → K → Γ on the
k-map. **One dataset per segment** lists each leg separately.

| Path on the contour | Points |
|---|---|
| ![Arbitrary cut path](images/arbitrary_cut_path.png) | ![Arbitrary cut dialog](images/arbitrary_cut_dialog.png) |

The result is a cut against the distance along the path; the corners of the
path are marked with dashed green lines, so a kink there is not mistaken for
physics.

![The cut along Γ–M–K–Γ](images/arbitrary_cut_result.png)

### 9.2 Fermi-surface (FS) correction

*Functions → FS correction...* on a cut (or a map's slit cut) straightens a
feature that should be flat — usually the Fermi edge, which curves across the
analyser angle. Click along the feature; a polynomial (order 2 by default) is
fitted and every angle column is shifted in energy. **From a reference
(Au)...** takes the curve from a gold spectrum in the list instead. On a
map's slit cut the same shift is applied to the whole cube.

![FS correction](images/fs_correction_dialog.png)

---

## 10. Processing: smoothing, derivatives, curvature

Select one or more cuts (or maps) in the list and press **Process...**.
Every result is a new dataset carrying the steps that made it.

![2-D processing](images/process_dialog.png)

- Tabs: **Smooth** (Gaussian, box, …; optional despiking), **Derivative**,
  **Curvature** (1-D and 2-D), **Background**, **Symmetrise**.
- **Source** and **Result** are shown side by side on the same scale, so you
  can see when a setting has eaten the band. They pan and zoom together.
- **Shape**: *Square pixels* (one data point is square — undistorted for an
  angle-vs-energy cut), *Equal axes* (one unit is the same length on both
  axes — right for angle-vs-angle or k-vs-k; chosen automatically when both
  axes have the same unit), or *Fill the box*.
- **Live preview** re-computes as you change a value; **Apply** lists the
  result.

For a map the panel works plane by plane, and asks **which planes** first —
constant-energy maps or dispersions — since the curvature of a contour and
the curvature of a cut are different quantities. **Preview at** picks the
plane shown; **Apply to the whole cube** runs it on every plane.

![3-D processing](images/volume_process_dialog.png)

**Data operations...** (main window) truncates, self-normalises or compresses
(bins) one or several datasets at once — as long as they have the same kind
and the same number of points on every axis.

![Data operations](images/data_operations.png)

---

## 11. Stack plots and the 3-D view

**Stack plot** (right-click a cut): a waterfall of EDCs or MDCs, with the
number of curves, binning, normalisation, offset, shear and colouring
controls; **Send to figure composer** makes it a figure panel.

![Stack plot](images/stack_plot.png)

**3D view** (right-click a map) first asks for the region and the sampling,
since every frame resamples the cube — *Fast*, *Balanced* or *Fine* decide
whether it keeps up with the mouse.

![3-D view: region and sampling](images/volume_range_dialog.png)

Then drag the picture to turn it. **Slices** shows three orthogonal planes at
the positions you set, **Notched cube** cuts a corner out of the cube, and
**Projection** sums along an axis. Colormap, levels, gamma, azimuth,
elevation and depth aspect are on the right; **Save the view...** exports
the picture.

![3-D view](images/volume_view.png)

---

## 12. Cut arithmetic: dichroism and reference division

Two cuts combined into a third: linear dichroism (LH − LV), circular
dichroism, dividing by a gold reference, or plain differences, ratios and
sums. Select two cuts in the list → right-click → **Cut arithmetic on the
two...** (the first selected is A), or in a cut's viewer *Functions → Cut
arithmetic...* and then click the other cut in the list.

![Cut arithmetic](images/cut_arithmetic.png)

- **The table at the top** compares what the two measurements recorded
  (polarisation, photon energy, angles, sample position, temperature, pass
  energy, lens mode). Rows in orange differ beyond a tolerance — a dichroism
  map between two spots on the sample is a map of *that* difference.
- **Purpose** presets: *Linear dichroism* (A − B, B scaled to the same total
  intensity), *Circular dichroism* ((A − B)/(A + B)), *Divide by a reference*,
  *Custom*. Warnings appear if the polarisations look swapped (**Swap A ↔ B**).
- **Scale B to A** by total intensity or by a region; **Hide where A+B below**
  blanks the background.
- A and B are drawn on one colour scale; the result on a red-white-blue scale
  symmetric about zero. With raw counts on both, the Poisson uncertainty is
  propagated and summarised under the plot.
- **Result to list**, **Uncertainty to list**, **All three as a figure**.

---

## 13. Photon-energy scans (kz)

A photon-energy scan is listed as a `kz map`; its contour is photon energy
against angle at one energy.

![A kz map](images/kz_map_contour.png)

### 13.1 kz map processing

Nothing in the files says where each spectrum's Fermi level really is.
*Functions → kz map processing...* measures it from the spectra themselves.

1. The slit cut opens with a box on it. Drag the box over a Fermi edge — the
   angle range sets the channels summed, the energy range the fit window.

   ![The box on the slit cut](images/kz_processing_cut.png)

2. **Fit and assemble**: the same box (by detector index) is used for every
   spectrum; each is shifted to put its own edge at zero; the cube is cropped
   to the range they all cover; and, by default, each spectrum is
   **normalised by its total intensity**, which divides out the photon flux.

![kz map processing](images/kz_processing_dialog.png)

E_F against photon energy is plotted. A smooth curve means the fits are
sound; points whose fit is not believed are marked and filled in from their
neighbours.

### 13.2 kz → momentum

*Functions → kz -> momentum...* converts the scan to (k_z, k∥).

![kz -> momentum](images/kz_conversion_dialog.png)

- **Inner potential V0** is the unknown: the preview reconverts one energy
  slice as you change it, in milliseconds.
- With the **lattice** (space group, a, b, c), the **surface normal** and the
  in-plane direction, the zone boundaries are drawn over the preview.
- **Scan V0** converts at a range of V0 values and reports where the k_z
  period matches the lattice — and how well that pins V0.
- **Measure a k_z period...**: click two points that are the same feature one
  zone apart; the matching reciprocal-lattice planes are listed.
- The yellow banner warns when the Fermi edge still bends across the analyser
  angle — correct it first (FS correction on the slit cut).
- **Convert the whole cube** lists a `kz map (k)`.

![The converted kz map](images/kz_map_in_k.png)

---

## 14. Removing the detector grid

MCP (hexagonal) and mesh (square) patterns print on every image. **Do it
first**: only data still on the detector's pixels can be de-gridded, so the
button refuses anything already converted, corrected or smoothed, and says
why.

*Functions → De-grid map...* on a map (or its slit cut), or *De-grid...* on a
cut. For a map, no reference is needed — the grid stays on the same pixels in
every slice while the photoemission moves. **Find and remove the grid**:

![De-grid](images/degrid_dialog.png)

- **Before / After** a chosen slice on one colour scale (tick **Show
  −∂²I/∂E²** to make any grid obvious).
- **The grid found**, its **k-space regions** (in red), and the **contrast
  and shift** fitted for every slice.
- The power left in the grid regions (1 = no grid left).
- **Also list the grid pattern** keeps a `[grid]` dataset, which a cut taken
  later with the same lens mode and pass energy will find and use.

---

## 15. Curves: EDCs, MDCs and spin-resolved data

**EDC → list** / **MDC → list** under any image, or a spin-resolved file,
give a curve dataset, which opens in the **curve viewer**:

![Curve viewer](images/curve_viewer.png)

Channels can be hidden, offset into a waterfall (**Offset**), shown with
**Error bars** (from a σ channel), and **Range** marks an x window used by
the operations. *Functions → Data operations* gives **Crop**, **Bin**,
**Normalise**, **Subtract a background** (constant, linear, Shirley),
**Shift the axis** and **Add counting errors (√N)**; each makes a new dataset
and propagates σ correctly.

### 15.1 Curve fit

*Functions → Curve fit...*: **Peaks on a background** (Lorentzian, Gaussian
or Voigt peaks placed by clicking or **Find peaks**; none / constant /
linear / Shirley background; instrumental resolution; optional **× Fermi
edge**), or the **Fermi edge** itself. Weighting is the channel's own σ,
√N, or uniform. Every parameter comes with its error; the residual is shown
in units of σ.

![Curve fit](images/curve_fit_dialog.png)

**Fit to list** keeps data, fit, residual and components; **Copy results**
copies a tab-separated table. A Fermi-edge fit can **List a copy with
E_F = 0**.

### 15.2 Spin analysis

A spin EDC opens with one curve per spin channel:

![Spin EDC](images/spin_edc_viewer.png)

*Functions → Spin analysis...*:

![Spin analysis](images/spin_analysis.png)

- The table reads each channel's label (manipulator setting, target
  magnetisation, axis and sign counted as "up"); **Axis** and **Counts as**
  can be overridden.
- **Cross ratio** (default when the channels allow it) cancels both the
  transmission and the target reflectivity exactly; **one pair** is available
  but carries the instrumental asymmetry, which is measured and reported (ε).
- **S_eff** (Sherman function, default 0.2), **Bin**, optional **Zero
  reference** window.
- Plots: raw channels, P with its error bars, and I↑ / I↓ with errors.
- **P to list**, **I↑ / I↓ to list**, **As a figure**.

---

## 16. Spatially resolved scans (nano-ARPES)

A real-space scan opens the **spatial map** (left, the total intensity at
each position) together with the **E vs k spectrum at the cursor** (right).
Drag the red target to read the spectrum at another point. Right-click a
panel for the two integrations:

- on the map, **Integrate selection → E vs k panel** sums the spectra over a
  region of the sample;
- on the spectrum, **Integrate selection → spatial panel** maps the intensity
  of a chosen (k, E) window across the sample — how a domain or a flake edge
  is found.

![Spatial scan](images/spatial_scan.png)

---

## 17. Figures and export

### 17.1 Figure composer

Select cuts (or curves) in the list → right-click → **Plot as a figure**. The
composer lays panels out for a paper: one size per column, labels only on
the outer edges, (a)/(b)/(c) letters, shared colour scales if you want them.

![Figure composer](images/figure_composer.png)

**Plot tools** sets everything (size in mm, fonts, ticks, colour bars, lines,
arrows, text, insets, scale bars); **Add panel...** adds another dataset;
**PDF**, **PNG**, **SVG** export, **Copy** puts it on the clipboard.

### 17.2 Slice figure

From a map: *Functions → Slice figure...* lays out a series of slices through
the cube — constant-energy contours, or cuts — evenly spaced or at chosen
values, with integration, one shared colour scale and titles.

| Settings | Result |
|---|---|
| ![Slice figure settings](images/slice_figure_dialog.png) | ![Slice figure](images/slice_figure.png) |

### 17.3 Exporting a panel

Right-click any image → **Export** → **Data (.nxs)...**, **Image...** or
**Axes only (.eps)...**. It always exports exactly what the panel shows
(range and all): the data as a dataset this program reads back; an image
with the colormap applied (PNG, TIFF or JPEG, at a chosen size and dpi, with
or without the axes); or the frame, ticks and labels as vector PostScript to
place the image in, in Illustrator.

![Export](images/export_dialog.png)

**Save** (right-click rows in the list) writes datasets to a `.nxs` file of
your own; several selected rows go into one file.

---

## 18. Autosave, memory and recovery

- **Autosave.** Every computed dataset (a k-map, a processed cut, a fit
  result…) is written, uncompressed, to this session's folder in
  `~/.arpes_viewer/sessions/` the moment it exists. Real-space scans are not
  (they are the measurement itself, and large).
- **A normal close deletes the autosaved copies.** When you close the main
  window and confirm (everything you want is saved, or *Close anyway*), the
  open viewers are closed and the session folder is deleted straight away.
  A dataset's copy also goes as soon as you save it to a file of your own.
- **Recovery after a crash.** Only if the program did not close normally (a
  crash, a kill) is the session folder kept: the next start offers its
  datasets back, and it is deleted after three days.
- **Operations log.** `~/.arpes_viewer/sessions/operations.log`, one line per
  store / save / discard — how to find what was where.
- **Memory box.** This program's memory, the system's, how many viewers and
  HDF5 files are open, and how much computed data is held. **Free memory**
  drops cached copies of computed data (they are re-read from their autosaved
  file when needed); **Close all viewers** closes every window and then frees
  memory.
- **Open auto-save folder** (below the Memory box) shows the session folders
  in the system's file browser, to delete old working copies by hand. The
  folder of the session now running (named in the status bar) is in use.
- **Long operations** (reading a big map, converting, de-gridding, assembling
  a folder) run off the GUI thread with a progress bar and a working Cancel.
- **Closing** with unsaved computed data asks first, and names the datasets:
  **Save to a file...**, **Close anyway** (the autosaved copies are deleted)
  or **Cancel** (nothing happens).

## 19. Tips and troubleshooting

- **"Not available here" in the Functions menu** — hover over the entry: the
  tooltip says why (e.g. the MDC/EDC fit needs a cut already in Å⁻¹).
- **A k conversion is not offered** — the first axis of the map is not an
  angle (a photon energy, temperature…). Check *First axis of a map is* in
  the loader.
- **Cut k conversion of an E − E_F cut** — the energy offset is filled in for
  you from the Fermi fit; if you shifted the axis some other way, type the
  kinetic energy of E_F (photon energy − work function) into *Energy offset*.
- **A file is misread** — choose the reader explicitly in the loader, and see
  *Assumptions worth knowing about* in the README (axis order in SOLEIL
  files).
- **Dragging the cursor on a big map is sluggish** — a map read lazily from
  its file re-slices in about 0.1 s per window, so the other windows wait for
  the mouse to pause. Once the map has been read into memory it follows live.
- **The machine is paging** — watch the Memory bar (amber past a third of the
  RAM, red past a half); use **Free memory** or **Close all viewers**.
