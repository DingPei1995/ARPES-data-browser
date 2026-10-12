# ARPES viewer — User Guide

*[中文版 / Chinese version](USER_GUIDE_zh.md) · [Version française](USER_GUIDE_fr.md)*

**Program version 2026-10-12** · source code:
<https://github.com/DingPei1995/ARPES-data-browser>

This guide has four parts:

- **Part I — Installation and setup**: the version, where the code is, what
  data the program reads, and how to set up Python and start it.
- **Part II — Looking at data during the experiment**: the main window, the
  right-click menus, and the viewers for SPEM, map and cut data. These are
  the tools for deciding quickly, at the beamline, where to measure next.
- **Part III — Advanced processing for publication**: de-gridding, the Fermi
  level, k conversion, the Brillouin zone, kz, curvature, dichroism, peak
  fits and figures, in the order a paper's data usually goes through them.
- **Part IV — Troubleshooting**: autosave and crash recovery, memory, and
  answers to common problems.

Every screenshot was taken from the running program on the synthetic data
written by `devtools/demo_data.py`, so you can follow along without a
beamtime ([§3](#3-trying-it-with-demo-data)). The formulas, tests and design
choices behind each step are in the [README](../README.md) and the
[CHANGELOG](../CHANGELOG.md).

---

## Contents

**Part I — Installation and setup**

1. [About the program](#1-about-the-program)
2. [Installing Python and starting the program](#2-installing-python-and-starting-the-program)
3. [Trying it with demo data](#3-trying-it-with-demo-data)

**Part II — Looking at data during the experiment**

4. [The main window and its right-click menu](#4-the-main-window-and-its-right-click-menu)
5. [Loading data](#5-loading-data)
6. [Image windows: controls and the right-click menu](#6-image-windows-controls-and-the-right-click-menu)
7. [SPEM (real-space) scans: finding the spot to measure](#7-spem-real-space-scans-finding-the-spot-to-measure)
8. [Maps: the contour and its two cuts](#8-maps-the-contour-and-its-two-cuts)
9. [Cuts: quick looks](#9-cuts-quick-looks)
10. [Photon-energy scans, curves and spin EDCs: quick looks](#10-photon-energy-scans-curves-and-spin-edcs-quick-looks)
11. [The 3-D view of a map](#11-the-3-d-view-of-a-map)

**Part III — Advanced processing for publication**

12. [The usual processing order](#12-the-usual-processing-order)
13. [Removing the detector grid](#13-removing-the-detector-grid)
14. [Fermi level and Fermi-surface correction](#14-fermi-level-and-fermi-surface-correction)
15. [From angles to momentum](#15-from-angles-to-momentum)
16. [The Brillouin zone and arbitrary cuts](#16-the-brillouin-zone-and-arbitrary-cuts)
17. [Photon-energy scans: kz processing and kz → momentum](#17-photon-energy-scans-kz-processing-and-kz--momentum)
18. [Processing: smoothing, derivatives, curvature](#18-processing-smoothing-derivatives-curvature)
19. [Cut arithmetic: dichroism and reference division](#19-cut-arithmetic-dichroism-and-reference-division)
20. [MDC / EDC fit: dispersion, v_F, m\*, self-energy](#20-mdc--edc-fit-dispersion-v_f-m-self-energy)
21. [Curve fit and spin analysis](#21-curve-fit-and-spin-analysis)
22. [Figures, export and saving](#22-figures-export-and-saving)

**Part IV — Troubleshooting**

23. [Autosave, crash recovery and memory](#23-autosave-crash-recovery-and-memory)
24. [Common problems](#24-common-problems)

---

# Part I — Installation and setup

## 1. About the program

| | |
|---|---|
| **Version** | 2026-10-12 (last modified 12 October 2026). The version is also at the top of the [README](../README.md); [CHANGELOG.md](../CHANGELOG.md) lists every change. |
| **Source code** | <https://github.com/DingPei1995/ARPES-data-browser> (the program is in the `ARPES_viewer/` folder) |
| **Written in** | Python 3, PyQt5 and pyqtgraph; numpy, scipy and h5py for the numerics and the files |
| **Platforms** | Windows, macOS and Linux — anywhere PyQt5 runs |

### Data the program reads

| Beamline / source | Files | Opens as |
|---|---|---|
| SOLEIL **ANTARES** | NeXus/HDF5 `*.nxs` | cuts, deflector maps, real-space (SPEM) scans |
| SOLEIL **CASSIOPEE** | Scienta SES text files: single spectra, and numbered folders | cuts; a folder becomes one map or one photon-energy (kz) scan |
| CASSIOPEE **spin end station** (MBS A-1) | `*.krx` and `*.txt` | cuts, deflector maps, spin-resolved EDCs |
| This program | its own saved `*.nxs` | everything it saves, with the steps that made it |

Adding a beamline means adding one reader file to `loader/` (see *Adding a
beamline* in the README).

### Data kinds

Every dataset in the program has one of these kinds. The kind decides which
viewer opens it and which functions it is offered.

| Kind | Listed as | Axes | What it is |
|---|---|---|---|
| `cut` | `[Cut]` | angle (or k) × E | one spectrum |
| `map` | `[Map]` | angle × angle × E | a deflector or polar-angle scan |
| `k_map` | `[k-map]` | kx × ky × E | a map converted to momentum |
| `kz_map` | `[kz map]` | hν × angle × E | a photon-energy scan |
| `kz_map_k` | `[kz map (k)]` | k_z × k∥ × E | that scan converted to momentum |
| `spem_4d`, `spem_1d` | `[SPEM]` | x (× y) × angle × E | real-space (nano-ARPES) scans |
| `edc`, `mdc`, `spin_edc` | `[EDC]`, `[MDC]`, `[Spin EDC]` | E (or k) × channels | curves, with optional σ channels |

## 2. Installing Python and starting the program

**Python.** Use Python **3.9 or newer** (it has been tested up to 3.13).
Any distribution works: python.org, Anaconda/Miniconda, or the system
Python on Linux.

**Packages.** Five are needed:

```bash
pip install PyQt5 pyqtgraph numpy h5py scipy
```

It is cleanest to install them into an environment of their own, so the
program does not depend on, or upset, other Python work on the machine:

```bash
# with venv
python -m venv arpes-env
source arpes-env/bin/activate          # Windows: arpes-env\Scripts\activate
pip install PyQt5 pyqtgraph numpy h5py scipy

# or with conda
conda create -n arpes python=3.12 pyqt pyqtgraph numpy h5py scipy -c conda-forge
conda activate arpes
```

Versions it was last tested with: PyQt5 5.15, pyqtgraph 0.14 (0.13 works
too), numpy 2.3, h5py 3.16, scipy 1.18.

**Getting the code.** Download it from the GitHub page (*Code → Download
ZIP*), or:

```bash
git clone https://github.com/DingPei1995/ARPES-data-browser.git
cd ARPES-data-browser/ARPES_viewer
```

**Starting.** From the `ARPES_viewer` folder:

```bash
python ARPES_viewer.py
```

**Settings and working files** are kept in the folder `.arpes_viewer` in
your home directory:

| File | What it holds |
|---|---|
| `~/.arpes_viewer/config.json` | the last data folder, and the colormap a new window starts with |
| `~/.arpes_viewer/sessions/` | the autosaved working copies of computed data ([§23](#23-autosave-crash-recovery-and-memory)) |
| `~/.arpes_viewer/sessions/operations.log` | one line for every autosave, save and deletion |

Deleting `config.json` resets the settings; nothing else is stored.

**Checking the installation.** The test suite runs without a display:

```bash
python -m pytest            # from the ARPES_viewer folder; needs pip install pytest
```

## 3. Trying it with demo data

```bash
python devtools/demo_data.py            # writes ./demo_data/*.nxs
```

This writes six files of synthetic data: a Fermi-surface map, a pair of cuts
taken with LH and LV polarisation, a photon-energy scan, a spin-resolved EDC,
a small real-space scan of a flake with two domains, and a map with a
detector mesh printed on it. They come from a toy band structure (a
hexagonal electron pocket and a hole band at Γ, hole pockets at K) with a
Fermi cut-off, the analyser resolution and Poisson noise, and open exactly
like measured data. Point the main window's folder box at `demo_data/` and
add them.

To re-take every screenshot in this guide after changing a window:

```bash
python devtools/guide_screenshots.py    # writes docs/images/*.png
python devtools/build_guide_pdfs.py     # rebuilds the PDFs (needs pandoc and playwright)
```

---

# Part II — Looking at data during the experiment

## 4. The main window and its right-click menu

![The main window](images/main_window_info.png)

The main window is a **browser**: one row per *dataset*, not per file. One
`.nxs` file often holds several measurements, and each gets its own row. The
kind is written after the name: `[Cut]`, `[Map]`, `[kz map]`, `[k-map]`,
`[SPEM]`, `[EDC]`, `[Spin EDC]` …

| To … | do this |
|---|---|
| open a dataset in its own window | **double-click** the row, or **Open selected** |
| work on one or several datasets | select them and **right-click** (menu below) |
| read the axes: first and last value, points, step | right-click → **Show information** (fills *Data Information*) |
| truncate, normalise or bin several datasets at once | **Data operations...** ([§18](#18-processing-smoothing-derivatives-curvature)) |
| smooth, differentiate, curvature, background, symmetrise | **Process...** ([§18](#18-processing-smoothing-derivatives-curvature)) |
| see and free memory | the **Memory** box at the bottom ([§23](#23-autosave-crash-recovery-and-memory)) |

### The list's right-click menu

Select one or more rows and right-click. The menu acts on the **whole
selection** (shift-click and ctrl-click select several rows):

![The list's right-click menu, with two cuts selected](images/list_menu.png)

| Entry | Works on | What it does |
|---|---|---|
| **Open** | one dataset | opens it in its viewer — the same as a double-click |
| **Open slit cut** / **Open deflector cut** | maps (one or several) | opens just that cut of each map, without the contour |
| **Show information** | one dataset | kind, axes and full metadata in the information window |
| **Rename...** | one dataset | changes the name it is listed and saved under (the file on disk keeps its name) |
| **Plot as a figure...** | cuts or curves | one panel per dataset in the figure composer ([§22](#22-figures-export-and-saving)) |
| **Stack plot (EDC / MDC)...** | one cut | a waterfall of EDCs or MDCs ([§9](#9-cuts-quick-looks)) |
| **MDC / EDC fit...** | one cut in k | the peak-fit panel ([§20](#20-mdc--edc-fit-dispersion-v_f-m-self-energy)) |
| **3D view...** | one map | slices, notched cube and projections ([§11](#11-the-3-d-view-of-a-map)) |
| **Cut arithmetic on the two...** | exactly two cuts | dichroism, ratios, differences ([§19](#19-cut-arithmetic-dichroism-and-reference-division)); the first selected is A |
| **Save dataset(s)...** | any | writes the selection to one `.nxs` file of your own |
| **Remove from list** | any | takes the rows out of the list, after asking; files on disk are never deleted |
| **Session log...** | — | what has been autosaved, saved and deleted, and where ([§23](#23-autosave-crash-recovery-and-memory)) |

The menu is decided from the selection alone, without reading any data, so it
opens at once however large the files are. An entry that does not fit the
selection is **greyed out**; hover over it and its tooltip says why (for
example, *Cut arithmetic* needs exactly two cuts, the *3D view* needs a map).

**Show information** also opens the full metadata of the dataset in a window
of its own, with a filter box and **Export metadata (.csv)**:

![Dataset information](images/information_window.png)

**Clear list** asks first, and offers to save anything computed in this
session that has not been saved yet. Files on disk are never deleted.

## 5. Loading data

**Load data...** opens the loader. **Quick add...** skips it and adds the
chosen files with everything detected automatically — the fast way during a
beamtime.

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

## 6. Image windows: controls and the right-click menu

Each dataset opens in its own window, so a map, its cuts and another file
can sit side by side.

![A cut viewer with the readout cursor](images/cut_viewer.png)

- **Top row**: the **Functions** menu (everything the window can do), the
  colormap (per window) with **Flip**, and the view controls — axis ranges,
  **Auto**, **Inv** (invert), **Reset**, **Grid**.
- **Above each image**: the **Levels** bar (drag the two handles, or type
  **Min%/Max%** percentiles and press **Clip**; *right-click* **Clip** for
  ready-made pairs: 0.5–99.5 %, 1–99 %, 2–98 %, 5–95 %, 0–100 %) and **Gamma**.
- **Under each image**: **Smooth** (display only; the data are not changed)
  and **Lock x/y** (aspect ratio).
- **Mouse**: drag to pan, wheel to zoom, **right-drag** to zoom each axis
  separately.

### The image right-click menu

Right-click on any image:

![The image right-click menu](images/image_menu.png)

| Entry | What it does |
|---|---|
| **View All** | fits the view back to the whole image |
| **Readout cursor** | a crosshair you drag. The **Cursor** row under the image shows its position and the value there; type a position and press Enter to put it exactly there. On a cut it also draws the **EDC** (right) and **MDC** (below) through that point; **EDC ±** / **MDC ±** integrate over a window, shown as a shaded band. |
| **Selection box** | a rectangle with two corner handles. Its corners appear in the **Selection** row under the image, where they can be typed (then **Refresh**). In a SPEM window the box is what you integrate over ([§7](#7-spem-real-space-scans-finding-the-spot-to-measure)); everywhere it sets the range of the stack plot. |
| **Line profile...** | the intensity along a line at any angle, with an integration width ([below](#line-profile)) |
| **Stack plot of the selection...** | a waterfall of EDCs or MDCs across the selection box — or across the whole panel when there is no box. The entry shows the range it will use. |
| **Export → Data (.nxs) / Image / Axes only (.eps)** | saves exactly what this panel shows ([§22](#22-figures-export-and-saving)) |

SPEM windows add three entries of their own, described in
[§7](#7-spem-real-space-scans-finding-the-spot-to-measure). pyqtgraph's own
menu items (mouse mode, axis submenus, plot options, its export dialog) are
removed on purpose: what was useful in them is on the visible view bar.

**EDC → list** / **MDC → list** (next to the curves) put the curve on screen
into the main list as a dataset of its own
([§10](#10-photon-energy-scans-curves-and-spin-edcs-quick-looks)).
**Functions → Slice → Save slice to the main list / Open slice in a new
panel** keeps a copy of the slice on screen.

### Line profile

**Right-click → Line profile...** puts an orange line across the image and
opens a window with the intensity along it. Drag the line, or its ends;
**Width** integrates that far either side (in the units of the vertical
axis), **Samples** sets the number of points, and **To list** keeps the
profile as a dataset. Closing the window removes the line. Use it for a band
that runs at an angle, where neither the EDC nor the MDC follows it.

| The line on the cut | The profile along it |
|---|---|
| ![Line on a cut](images/line_profile_cut.png) | ![Line profile](images/line_profile.png) |

### The Functions menu

| Viewer | Analysis | Data operations | Visualization |
|---|---|---|---|
| map / kz map / k-map | | Arbitrary cut, Map k conversion, kz map processing, kz → momentum, De-grid map | Brillouin zone (k-maps), Slice figure |
| cut | Fermi level, MDC / EDC fit | FS correction, Cut k conversion, Cut arithmetic, De-grid | |
| curve | Curve fit, Spin analysis | Crop, Bin, Normalise, Background, Shift axis, Add √N errors | As a figure |

An entry that does not apply is not shown; one that applies but cannot run
yet is greyed, with the reason in its tooltip.

## 7. SPEM (real-space) scans: finding the spot to measure

A real-space scan (nano-ARPES, `[SPEM]`) opens two panels side by side: the
**spatial map** on the left — the total intensity at each sample position —
and the **E vs k spectrum** on the right. The red target on the map chooses
the pixel whose spectrum is shown, and the right panel's title says which
pixel it is (*E vs k at x=20, y=24*).

![A SPEM scan: map and the spectrum at the red target](images/spatial_scan.png)

**Drag the red target** to read the spectrum anywhere on the sample. This is
the quickest way to see whether a flake shows bands at all, and how sharp
they are.

### Right-click on the spatial map

![Right-click on the spatial map](images/spem_spatial_menu.png)

| Entry | What it does |
|---|---|
| **Selection box** | draws a box on the map |
| **Integrate selection -> E vs k panel** | the spectrum **summed over every pixel in the box**: more counts, and the average over one region of the sample |
| **Spatial map over the whole spectrum** | appears after a spectrum box has been integrated onto the map (below): goes back to the map summed over all k and E |
| **Open this position in a new window** | a snapshot of the spectrum shown now, in a window of its own — keep several open to compare spots |
| **Swap X and Y** | draws the map transposed, for a lab whose convention has the other stage horizontal. The cursor and the selection box stay on the same sample pixels; only the drawing changes. |
| Readout cursor, Line profile, Stack plot, Export | as on every image ([§6](#6-image-windows-controls-and-the-right-click-menu)) |

**Integrating a region of the sample.** Right-click the map → **Selection
box**, drag the box over the region (or type its corners in the
**Selection** row and press **Refresh**), then right-click → **Integrate
selection -> E vs k panel**. The spectrum is now the sum over the box, and
its title says so (*E vs k, summed over x 14 to 26, y 14 to 30 (63 pixels)*):

![The spectrum summed over a box on the flake](images/spem_box_to_spectrum.png)

### Right-click on the spectrum: where on the sample does a feature come from?

![Right-click on the spectrum](images/spem_spectrum_menu.png)

The opposite direction. Right-click the spectrum → **Selection box**, frame a
feature in (k, E) — a band, a shoulder, the Fermi edge of one phase — and
right-click → **Integrate selection -> spatial panel**. The map is redrawn as
the intensity of **that window only** at every position, and its title says
which window. Here the box frames a band that only one domain of the flake
has, and the map lights up exactly that domain:

![The map of the band in the box: one domain lights up](images/spem_box_to_map.png)

This is how domains, flake edges, terraces or a second phase are found. To
get the full map back, right-click the map → **Spatial map over the whole
spectrum**.

**A practical routine at the beamline:**

1. Open the overview scan; drag the red target over the flake and look at
   the spectra.
2. Box a feature in the spectrum → integrate onto the map to see where that
   feature is strongest.
3. Box the best region of the map → integrate onto the spectrum to check the
   bands with good statistics.
4. **Open this position in a new window** for the candidates, compare them
   side by side, and read the coordinates off the cursor row for the next
   measurement.

A line scan (`spem_1d`) opens the same way, with position × angle on the left;
the box integration needs a 2-D map, so it is offered only for 2-D scans.

## 8. Maps: the contour and its two cuts

A map opens on its **constant-energy contour**. The slider at the top picks
the energy (**Pos** in eV, or **Ind** as an index; **±** integrates over a
window, **mean** or **sum**). **Deflector cut** and **Slit cut** open the two
orthogonal cuts, each in its own window.

![Constant-energy contour of a map](images/map_contour.png)

| Slit cut (slit angle vs energy) | Deflector cut (deflector angle vs energy) |
|---|---|
| ![Slit cut](images/map_slit_cut.png) | ![Deflector cut](images/map_deflector_cut.png) |

**One cursor for all three windows.** Switching the readout cursor on in any
of them (right-click → *Readout cursor*) switches it on in the others, and
moving it anywhere re-slices the others through the new point. The axes keep
one colour each everywhere: **deflector red, slit green, energy blue**. A
cursor line has the colour of the axis it holds constant, and the EDC or MDC
it feeds is drawn in the same colour. Integration widths (the ± boxes) belong
to each window and are never shared.

From the list you can also open a map's cut straight away (right-click →
**Open slit cut** / **Open deflector cut**); the contour then stays hidden
until you ask for it (*Functions → Show the map*).

During a beamtime this answers *is the sample aligned, where is Γ, which
direction is the cut going to follow* — before any conversion.

## 9. Cuts: quick looks

A cut opens in a single image panel ([§6](#6-image-windows-controls-and-the-right-click-menu)).
The right-click tools give fast answers:

- **Readout cursor** with **EDC ±** / **MDC ±**: the curves through any point,
  integrated over a window.
- **Line profile...**: the intensity along any line.
- **Stack plot**: right-click the cut (or its row in the list) → **Stack plot
  of the selection...**. The curves run across the selection box, if there is
  one: a box taller than it is wide gives EDCs across its angle range, a
  wide one MDCs across its energy range. The window sets the number of
  curves, binning, normalisation, offset, shear and colours; **Send to figure
  composer** makes it a figure panel.

![Stack plot](images/stack_plot.png)

- **Fermi level**: *Functions → Fermi level...* gives E_F and the resolution
  in a few seconds ([§14](#14-fermi-level-and-fermi-surface-correction)) —
  useful to check the energy calibration right after a photon-energy change.

## 10. Photon-energy scans, curves and spin EDCs: quick looks

**Photon-energy scans.** A photon-energy scan is listed as a `kz map`; its
contour is photon energy against angle at one energy, and its cuts open like
a map's. Before processing (§17) the Fermi level of each spectrum is not
aligned, so treat it as a quick look.

![A kz map](images/kz_map_contour.png)

**Curves.** **EDC → list** / **MDC → list** under any image, or a spin file,
give a curve dataset, which opens in the **curve viewer**:

![Curve viewer](images/curve_viewer.png)

Channels can be hidden, offset into a waterfall (**Offset**), shown with
**Error bars** (from a σ channel), and **Range** marks an x window used by
the operations. *Functions → Data operations* gives **Crop**, **Bin**,
**Normalise**, **Subtract a background** (constant, linear, Shirley),
**Shift the axis** and **Add counting errors (√N)**; each makes a new dataset
and propagates σ.

**Spin EDCs.** A spin-resolved EDC opens with one curve per spin channel, so
the count rates of all channels can be checked as they come in. The
polarisation itself is in *Functions → Spin analysis...*
([§21](#21-curve-fit-and-spin-analysis)).

![Spin EDC](images/spin_edc_viewer.png)

## 11. The 3-D view of a map

**3D view** (right-click a map in the list) first asks for the region and the
sampling, since every frame resamples the cube — *Fast*, *Balanced* or
*Fine* decide whether it keeps up with the mouse.

![3-D view: region and sampling](images/volume_range_dialog.png)

Then drag the picture to turn it. **Slices** shows three orthogonal planes at
the positions you set, **Notched cube** cuts a corner out of the cube, and
**Projection** sums along an axis. Colormap, levels, gamma, azimuth,
elevation and depth aspect are on the right; **Save the view...** exports
the picture.

![3-D view](images/volume_view.png)

---

# Part III — Advanced processing for publication

## 12. The usual processing order

Each step makes a **new dataset** in the list, which records the steps that
made it, so the raw data are never changed and every result can be traced.
The order matters, because some steps only work on raw detector data:

1. **De-grid** the raw map or cut, if a mesh or MCP pattern shows ([§13](#13-removing-the-detector-grid)).
2. **Fermi level**, and **FS correction** if the edge curves across the
   analyser angle ([§14](#14-fermi-level-and-fermi-surface-correction)).
3. **k conversion** of the cut or map ([§15](#15-from-angles-to-momentum)); for a
   photon-energy scan, **kz processing** then **kz → momentum** ([§17](#17-photon-energy-scans-kz-processing-and-kz--momentum)).
4. **Brillouin zone** overlay and **arbitrary cuts** along high-symmetry paths ([§16](#16-the-brillouin-zone-and-arbitrary-cuts)).
5. **Curvature / second derivative** to sharpen bands, **dichroism** ([§18](#18-processing-smoothing-derivatives-curvature), [§19](#19-cut-arithmetic-dichroism-and-reference-division)).
6. **Peak fits**: dispersion, v_F, m\*, self-energy ([§20](#20-mdc--edc-fit-dispersion-v_f-m-self-energy)); curve and spin fits ([§21](#21-curve-fit-and-spin-analysis)).
7. **Figures and export** ([§22](#22-figures-export-and-saving)).

## 13. Removing the detector grid

MCP (hexagonal) and mesh (square) patterns print on every image. **Do it
first**: only data still on the detector's pixels can be de-gridded, so the
function refuses anything already converted, corrected or smoothed, and says
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

## 14. Fermi level and Fermi-surface correction

### 14.1 Fermi level

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

### 14.2 Fermi-surface (FS) correction

*Functions → FS correction...* on a cut (or a map's slit cut) straightens a
feature that should be flat — usually the Fermi edge, which curves across the
analyser angle. Click along the feature; a polynomial (order 2 by default) is
fitted and every angle column is shifted in energy. **From a reference
(Au)...** takes the curve from a gold spectrum in the list instead. On a
map's slit cut the same shift is applied to the whole cube.

![FS correction](images/fs_correction_dialog.png)

## 15. From angles to momentum

### 15.1 Cut k conversion

*Functions → Cut k conversion...* converts the angle axis to k∥ (Å⁻¹). A cut
is a line in the (kx, ky) plane that generally misses Γ, so the dialog asks
where Γ is: type the angles, or **From a converted map...** to inherit them
from a map already converted.

![Cut k conversion](images/cut_k_dialog.png)

If the cut is already on E − E_F, the dialog fills in **Energy offset** with
the E_F that was subtracted (the note under it says so): the conversion uses
true kinetic energies, and the result stays on E − E_F.

![The cut in momentum](images/cut_in_k.png)

### 15.2 Map k conversion

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

## 16. The Brillouin zone and arbitrary cuts

### 16.1 Brillouin zone

On a k-map: *Functions → Visualization → Brillouin zone...*.

![Brillouin-zone dialog](images/brillouin_dialog.png)

1. **Space group** first (1–230); the lattice parameters it fixes are greyed
   and follow the free ones. Or switch to **2D layer** for a surface lattice.
2. **Cut plane**: the surface normal (hkl) and an offset along it — the cut of
   the 3-D zone by that plane is drawn, for any surface, not only (001).
   **Show cut plane in 3D...** previews it (drag to turn it, right-drag to
   zoom).
3. **Conventional** or **Irreducible** zone; **Azimuth** (or **Pick
   direction from contour...**); **Repeat across the field of view**.
4. **Moiré** adds a second, twisted or mismatched layer and its moiré zone.
5. **Plot**.

![k-map with the zone overlaid](images/k_map_with_bz.png)

### 16.2 Arbitrary cut

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

## 17. Photon-energy scans: kz processing and kz → momentum

### 17.1 kz map processing

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

### 17.2 kz → momentum

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

## 18. Processing: smoothing, derivatives, curvature

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

## 19. Cut arithmetic: dichroism and reference division

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

## 20. MDC / EDC fit: dispersion, v_F, m\*, self-energy

*Functions → MDC / EDC fit...* (or right-click the row → **MDC / EDC
fit...**) fits peaks line by line and reads the band off the fitted
centres. It only accepts a cut whose axis is already in Å⁻¹ — a velocity in
eV·deg would be meaningless — so convert the cut first ([§15](#15-from-angles-to-momentum)).

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

## 21. Curve fit and spin analysis

### 21.1 Curve fit

On a curve: *Functions → Curve fit...*. **Peaks on a background**
(Lorentzian, Gaussian or Voigt peaks placed by clicking or **Find peaks**;
none / constant / linear / Shirley background; instrumental resolution;
optional **× Fermi edge**), or the **Fermi edge** itself. Weighting is the
channel's own σ, √N, or uniform. Every parameter comes with its error; the
residual is shown in units of σ.

![Curve fit](images/curve_fit_dialog.png)

**Fit to list** keeps data, fit, residual and components; **Copy results**
copies a tab-separated table. A Fermi-edge fit can **List a copy with
E_F = 0**.

### 21.2 Spin analysis

On a spin EDC: *Functions → Spin analysis...*:

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

## 22. Figures, export and saving

### 22.1 Figure composer

Select cuts (or curves) in the list → right-click → **Plot as a figure**. The
composer lays panels out for a paper: one size per column, labels only on
the outer edges, (a)/(b)/(c) letters, shared colour scales if you want them.

![Figure composer](images/figure_composer.png)

**Plot tools** sets everything (size in mm, fonts, ticks, colour bars, lines,
arrows, text, insets, scale bars); **Add panel...** adds another dataset;
**PDF**, **PNG**, **SVG** export, **Copy** puts it on the clipboard. A stack
plot (**Send to figure composer**) and a dichroism result (**All three as a
figure**) come here too.

### 22.2 Slice figure

From a map: *Functions → Slice figure...* lays out a series of slices through
the cube — constant-energy contours, or cuts — evenly spaced or at chosen
values, with integration, one shared colour scale and titles.

| Settings | Result |
|---|---|
| ![Slice figure settings](images/slice_figure_dialog.png) | ![Slice figure](images/slice_figure.png) |

### 22.3 Exporting a panel

Right-click any image → **Export** → **Data (.nxs)...**, **Image...** or
**Axes only (.eps)...**. It always exports exactly what the panel shows
(range and all): the data as a dataset this program reads back; an image
with the colormap applied (PNG, TIFF or JPEG, at a chosen size and dpi, with
or without the axes); or the frame, ticks and labels as vector PostScript to
place the image in, in Illustrator. In a SPEM window, right-click the panel
you want — the map and the spectrum export separately.

![Export](images/export_dialog.png)

### 22.4 Saving datasets

Right-click rows in the list → **Save dataset(s)...** writes them to a `.nxs`
file of your own; several selected rows go into one file, each under its
listed name (use **Rename...** first). The saved file opens again with every
step that made each dataset, and the save is written to the session log.

---

# Part IV — Troubleshooting

## 23. Autosave, crash recovery and memory

- **Autosave.** Every computed dataset (a k-map, a processed cut, a fit
  result…) is written, uncompressed, to this session's folder in
  `~/.arpes_viewer/sessions/` the moment it exists. Real-space scans are not
  (they are the measurement itself, and large).
- **A normal close deletes the autosaved copies.** When you close the main
  window and confirm (everything you want is saved, or *Close anyway*), the
  open viewers are closed and the session folder is deleted straight away.
  A dataset's copy also goes as soon as you save it to a file of your own.
- **Recovery after a crash.** Only if the program did not close normally (a
  crash, a kill, a power cut) is the session folder kept: the next start
  offers its datasets back, and it is deleted after three days.
- **Operations log.** `~/.arpes_viewer/sessions/operations.log` (or
  right-click in the list → **Session log...**), one line per autosave,
  save and deletion, with paths — how to find what was where.
- **Closing** with unsaved computed data asks first, and names the datasets:
  **Save to a file...**, **Close anyway** (the autosaved copies are deleted)
  or **Cancel** (nothing happens).
- **Memory box.** This program's memory, the system's, how many viewers and
  HDF5 files are open, and how much computed data is held. **Free memory**
  drops cached copies of computed data (they are re-read from their autosaved
  file when needed); **Close all viewers** closes every window and then frees
  memory.
- **Long operations** (reading a big map, converting, de-gridding, assembling
  a folder) run off the GUI thread with a progress bar and a working Cancel.

## 24. Common problems

### Installation and start-up

- **`ModuleNotFoundError: No module named 'PyQt5'`** (or pyqtgraph, h5py…) —
  the packages are installed in a different Python from the one running the
  program. Activate the environment ([§2](#2-installing-python-and-starting-the-program))
  and check with `python -c "import PyQt5, pyqtgraph, h5py, scipy"`.
- **`qt.qpa.plugin: Could not load the Qt platform plugin "xcb"`** (Linux) —
  a system Qt library is missing; install it with your package manager
  (`sudo apt install libxcb-xinerama0 libxkbcommon-x11-0` on Ubuntu), or use
  the conda installation.
- **Fonts or windows too large or too small** (high-DPI screens) — the
  program follows the system scaling; set `QT_SCALE_FACTOR=0.9` (or another
  value) in the environment before starting it.
- **The settings look wrong** — delete `~/.arpes_viewer/config.json`.

### Reading data

- **A file is misread** — choose the reader explicitly in the loader, and see
  *Assumptions worth knowing about* in the README (axis order in SOLEIL
  files).
- **A k conversion is not offered** — the first axis of the map is not an
  angle (a photon energy, temperature…). Check *First axis of a map is* in
  the loader.
- **A CASSIOPEE folder is listed as separate spectra** — the `_i` parameter
  files next to the spectra are missing or renamed; they say what was
  stepped.

### Right-click menus

- **A list entry is greyed out** — hover over it: the tooltip says what the
  selection lacks (one row, two cuts, a map, a cut in k…).
- **"Integrate selection" is greyed** — draw a box first (right-click →
  **Selection box**).
- **"Integrate selection" is missing** — it exists only in SPEM windows with
  a 2-D map; on a line scan, a cut or a map the box is used for the stack
  plot and the coordinates.
- **The SPEM map stays dark or shows only part of the sample** — it is
  showing the intensity of a spectrum box (its title says so). Right-click
  the map → **Spatial map over the whole spectrum**.
- **The SPEM map looks transposed or mirrored** — **Swap X and Y** (right-
  click the map) and the **Inv** boxes in the view bar change only the
  drawing; ANTARES scans open the way the beamline draws them.
- **The line profile's line has gone** — it is removed when its window
  closes; open it again from the right-click menu.
- **The MDC / EDC fit refuses a cut** — its axis is still an angle; convert
  it first ([§15](#15-from-angles-to-momentum)).

### Processing

- **"Not available here" in the Functions menu** — hover over the entry: the
  tooltip says why.
- **De-grid refuses a dataset** — it has already been converted, corrected or
  smoothed; de-grid the raw data first ([§13](#13-removing-the-detector-grid)).
- **Cut k conversion of an E − E_F cut** — the energy offset is filled in for
  you from the Fermi fit; if you shifted the axis some other way, type the
  kinetic energy of E_F (photon energy − work function) into *Energy offset*.

### Speed and memory

- **Dragging the cursor on a big map is sluggish** — a map read lazily from
  its file re-slices in about 0.1 s per window, so the other windows wait for
  the mouse to pause. Once the map has been read into memory it follows live.
- **The machine is paging** — watch the Memory bar (amber past a third of the
  RAM, red past a half); use **Free memory** or **Close all viewers**.
- **After a crash** — restart the program: it offers the autosaved datasets
  back. The log (**Session log...**) says where each one is.
