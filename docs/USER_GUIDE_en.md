# ARPES viewer — User Guide

*[中文版 / Chinese version](USER_GUIDE_zh.md) · [Version française](USER_GUIDE_fr.md)*

*Last updated: 29 October 2026 — for ARPES viewer v1.2.2*

The version and its date are in the main window's title bar and at the bottom of
the main window. Bug reports, questions and the latest version:
<https://github.com/DingPei1995/ARPES-data-browser> (please give the version
number when reporting a bug).

This guide walks through the program's features with screenshots, in the
order of the tool panel: getting started (installing, loading, the main
window and the viewers), then the **Basic**, **Advanced** and
**Calculation** tabs, then figures — ending with two figures built step
by step in the figure composer. The technical
background for each step (formulas, test results, design choices) is in the
technical references beside the code — [reading data](../arpes_viewer/loader/README.md),
[the windows](../arpes_viewer/ui/README.md), [analysis](../arpes_viewer/tools/README.md) — and in the
[CHANGELOG](../CHANGELOG.md).

Every screenshot was taken from the running program, using the synthetic
data written by `devtools/demo_data.py`, so you can follow along without a
beamtime (see [Trying it with demo data](#2-trying-it-with-demo-data)).

---

## Contents

**[Part I — Getting started](#part-i--getting-started)**

1. [Installing and starting](#1-installing-and-starting)  
2. [Trying it with demo data](#2-trying-it-with-demo-data)  
3. [The main window](#3-the-main-window)  
4. [Loading data](#4-loading-data)  
5. [What the viewers have in common](#5-what-the-viewers-have-in-common)  
6. [Maps: the contour and its two cuts](#6-maps-the-contour-and-its-two-cuts)  
7. [Arbitrary cuts and Fermi-surface correction](#7-arbitrary-cuts-and-fermi-surface-correction)

**[Part II — The Basic tab](#part-ii--the-basic-tab)**

8. [Fermi edge fit](#8-fermi-edge-fit)  
9. [Removing the detector grid and detector corrections](#9-removing-the-detector-grid-and-detector-corrections)  
10. [Background subtraction](#10-background-subtraction)  
11. [Converting a cut to momentum](#11-converting-a-cut-to-momentum)  
12. [Core-level check](#12-core-level-check)  
13. [Converting a map to momentum](#13-converting-a-map-to-momentum)  
14. [Photon-energy scans (kz)](#14-photon-energy-scans-kz)  
15. [Brillouin zone](#15-brillouin-zone)

**[Part III — The Advanced tab](#part-iii--the-advanced-tab)**

16. [Smoothing, derivatives, curvature and symmetrising](#16-smoothing-derivatives-curvature-and-symmetrising)  
17. [Ridge detection](#17-ridge-detection)  
18. [Stack plots](#18-stack-plots)  
19. [MDC fit and EDC fit](#19-mdc-fit-and-edc-fit)  
20. [Cut arithmetic: dichroism and reference division](#20-cut-arithmetic-dichroism-and-reference-division)  
21. [Slice series figure](#21-slice-series-figure)  
22. [The 3-D views](#22-the-3-d-views)  
23. [Curves: EDCs, MDCs and spin-resolved data](#23-curves-edcs-mdcs-and-spin-resolved-data)  
24. [Spatially resolved scans (nano-ARPES)](#24-spatially-resolved-scans-nano-arpes)

**[Part IV — The Calculation tab](#part-iv--the-calculation-tab)**

25. [Crystal structure](#25-crystal-structure)  
26. [Tight-binding fit](#26-tight-binding-fit)  
27. [Calculated bands](#27-calculated-bands)

**[Part V — Figures and output](#part-v--figures-and-output)**

28. [The figure composer](#28-the-figure-composer)  
29. [Exporting a panel and saving data](#29-exporting-a-panel-and-saving-data)  
30. [Worked example: two figures in the figure composer](#30-worked-example-two-figures-in-the-figure-composer)

**[Part VI — Housekeeping](#part-vi--housekeeping)**

31. [Autosave, memory and recovery](#31-autosave-memory-and-recovery)  
32. [Tips and troubleshooting](#32-tips-and-troubleshooting)  

---

# Part I — Getting started

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

**On Windows, double-click `ARPES viewer.pyw`** in the program's folder,
not `ARPES_viewer.py`: Windows runs a `.pyw` with `pythonw`, which opens no
black console window -- closing that window by mistake used to close the
program with it. (A `.py` double-clicked all the same now hides the console
it was given; one started from a command prompt keeps its console.) The
program's messages go to `~/.arpes_viewer/launcher.log`; if it stops with an
error, a message box shows the end of that log.

**Errors are kept as a text file.** Whatever the console would show when
something goes wrong -- the error and its traceback -- is saved by itself in
**`error_logs/errors_<date>_<time>.txt` in the program's folder** (in
`~/.arpes_viewer/error_logs` if that folder cannot be written), only when an
error happens. A message says so, and the program **keeps running** instead
of closing. With a bug report, **upload that file** (on GitHub Issues: drag
it into the text box) **or paste its text** — the main panel's foot has
**Copy the newest error log** (its text to the clipboard) and **Open the
error logs** (the folder).

**Long operations show a progress bar.** Anything that can take more than a
moment on a big dataset or a slow computer — reading a whole map, k
conversion, FS correction, arbitrary cuts, data operations, joining
datasets, saving, ridge detection, the fits, and in Calculated bands reading
a calculation, the plane's bands and the 3-D Fermi surface — shows a
progress window after 0.4 s: a bar that moves, the seconds taken, and Cancel
where the step can stop. The program is working, not frozen.

**On a Mac you do not need Terminal to start it:** double-click
**`ARPES viewer.app`** in the program's folder (drag it to the Dock to keep
it there). It finds the Python that has the packages above -- the same one
as `python3` in Terminal, conda environments included -- and remembers it.
If no Python has them, a dialog gives the `pip install` command (already
copied to the clipboard: paste it into Terminal once) or lets you choose
a Python yourself.

- If macOS says the app *cannot be opened because it is from an unidentified
  developer* (the program was downloaded as a ZIP): open **System Settings →
  Privacy & Security** and click **Open Anyway** (older macOS: right-click
  the app → **Open**). Only the first time.
- A copy of the app moved out of the program's folder (to Applications, say)
  asks once where the folder is.
- If the program stops with an error, the dialog shows the end of
  `~/.arpes_viewer/launcher.log`; please send that file with a bug report.

The program reads:

| Beamline / source | Files |
|---|---|
| SOLEIL **ANTARES** | NeXus/HDF5 `*.nxs` (cuts, deflector maps, spatial scans) |
| SOLEIL **CASSIOPEE** | Scienta SES text files — single spectra, and numbered folders assembled into a map or a photon-energy scan |
| CASSIOPEE **spin end station** (MBS A-1) | `*.krx` and `*.txt` — cuts, deflector maps and spin-resolved EDCs |
| This program | its own saved `*.nxs` format |

## 2. Trying it with demo data

```bash
python devtools/demo_data.py            # writes ./demo_data/
```

This writes synthetic data that opens exactly like measured data:

- the original set -- a Fermi-surface map, a pair of cuts in LH and LV
  polarisation, a photon-energy scan, a spin-resolved EDC, a small
  real-space scan and a map with a detector mesh printed on it -- from a toy
  band structure (a hexagonal electron pocket and a hole band at Γ, hole
  pockets at K);
- **a layered crystal in three dimensions** (`demo_bulk_bands.nxs`): maps at
  the photon energies of the Γ plane (59.3 eV) and of the zone top (77 eV),
  cuts along Γ-M (LH and LV) and Γ-K, and a photon-energy scan along Γ-M
  through a whole period of k_z. Two of its bands cross along a line through
  the zone, whose energy changes with k_z;
- **a calculation of that crystal**, as a theory group would send it: a
  Wannier90 Hamiltonian (`demo_calc_hr.dat` with `demo_calc.win`) and the
  same bands on a grid (`demo_calc.bxsf`). It is the crystal before
  renormalisation -- its bands 1/0.8 times wider and 50 meV higher -- which is
  what fitting it to the data recovers ([§27](#27-calculated-bands));
- two crystal structures as CIF files (`demo_crystal.cif`, a 1T layered
  hexagonal crystal; `demo_tetragonal_crystal.cif`, a ZrSiS-type one) and an
  EDS spectrum (`demo_eds.nxs`).

All come with a Fermi cut-off, the analyser resolution and Poisson noise.

To re-take every screenshot in this guide after changing a window:

```bash
python devtools/guide_screenshots.py    # writes docs/images/*.png
```

The two example figures of [§30](#30-worked-example-two-figures-in-the-figure-composer)
are rebuilt, from the demo data and with the program's own windows, by

```bash
python devtools/figure_tests.py docs/images [folder for the .afig files]
```

## 3. The main window

![The main window](images/main_window_info.png)

The main window is a **browser**: one row per *dataset*, not per file (one
`.nxs` file often holds several measurements, and each gets its own row). The
kind is written after the name: `[Cut]`, `[Map]`, `[kz map]`, `[k-map]`,
`[SPEM]`, `[EDC]`, `[Spin EDC]` …

It has three parts:

- **The tool panel** (left) — every analysis and processing tool, as a
  button, in three tabs. **Basic** (shown at start) holds what nearly every
  dataset goes through: Fermi edge fit, detector-grid removal, background
  subtraction, the three k conversions, the core-level check, aligning a kz
  map's Fermi levels, the Brillouin zone and the detector corrections.
  **Advanced** holds the rest:
  figure composer, smoothing and the band-sharpening filters, symmetrising,
  ridge detection, stack plots, MDC and EDC fits, cut arithmetic, the slice
  series figure and the four 3-D views. **Calculation** holds what is
  worked out rather than measured: the crystal structure, the
  tight-binding fit and calculated band structures (see
  [Part IV](#part-iv--the-calculation-tab)).
  Within a tab the tools are grouped by the data they
  take (*2D & 3D — cuts, maps*, *2D only — cuts*, *3D only — maps* …).
- **The list** (right), with its buttons.
- **The Memory box** (bottom left) and the version line (bottom right).

![The Advanced tab of the tool panel](images/main_window_advanced.png)

**The tool buttons follow the selection.** A button lights up only when the
rows selected in the list are what it takes — a k conversion (cut) needs one
angle cut, cut arithmetic exactly two cuts, a 3-D view a map. A greyed button
says why when you hover over it.

| Action | How |
|---|---|
| Open a dataset in its own window | **double-click** the row (or right-click → **Open**) |
| Run a tool on the selected rows | its button in the **tool panel** |
| Stack several cuts into one map, or join a map / kz map measured in segments | **Combine...** under the list |
| Crop to a range, bin points, normalise each spectrum to itself | **Truncate...**, **Compress...**, **Self-normalise...** under the list |
| Open a map's slit or scan-axis cut | **right-click** the row (or a selection of rows) |
| Read the metadata, write the experiment log, rename, save, remove | the **Datasets** box of the tool panel, above Memory |
| See and free memory, close every viewer, open the autosave folder | the **Memory** box, bottom left |
| See which version you are running, report a bug, get the latest version | the title bar and the last line of the window (link to GitHub) |

The four buttons under the list act on the rows selected above them, and each
makes new rows — the originals are never changed.

**The right-click menu** holds only the ways of opening the rows: **Open**,
**Open slit cut** and **Open scan-axis cut** (for a map). An entry that does
not fit the selection is greyed, with the reason in its tooltip.

**The Datasets box** (tool panel, between the tools and Memory) holds the
rest of what is done to the rows: **Experiment info**, **Log format...** and
**Auto Exp. log** side by side (see [the experiment log](#the-experiment-log)),
**Rename**, **Save** and **Remove** on one row, **Sum repeated scans...**
and **Session log...**. Its buttons grey out like the tools.

### Where things are now (v1.2.2)

If you used an earlier version: the tools moved out of the right-click menu
and the viewers' *Functions* menus into the tool panel, so each has one way
in, always in sight.

| Now | Was |
|---|---|
| Tool panel → Basic → **Fermi edge fit...** (cuts and maps) | viewer *Functions → Fermi edge fit* |
| Tool panel → Basic → **Remove detector grid...** | viewer *Functions → Remove detector grid* |
| Tool panel → Basic → **k conversion (cut)** / **(map)**, **kz conversion** | viewer *Functions → Convert to momentum* |
| Tool panel → Basic → **Align Fermi levels (kz)...** | viewer *Functions → Align Fermi levels (kz map)* |
| Tool panel → Basic → **Brillouin zone...** | k-map viewer *Functions → Brillouin zone* |
| Tool panel → Basic → **Core-level check...** | right-click → Core-level check |
| Tool panel → Basic/Advanced → **Background subtraction**, **Smooth**, **Second derivative**, **Curvature**, **Symmetrise** (one window each) | **Smooth / derivative...** (one window with tabs) |
| **Truncate...**, **Compress...**, **Self-normalise...** under the list | **Crop / bin / normalise...** |
| Tool panel → Advanced → **MDC fit (k-cut only)**, **EDC fit** | viewer *Functions → MDC / EDC fit* |
| Tool panel → Advanced → **Cut arithmetic (two cuts)** | right-click / viewer *Functions → Cut arithmetic* |
| Tool panel → Advanced → **Slice series figure...** | map viewer *Functions → Multi-slice figure* |
| Tool panel → Advanced → **3D slices**, **Notched cube**, **Volume projection** | right-click → 3D view (one window, three modes) |
| Tool panel → Advanced → **Figure composer...** | right-click → Open in figure composer |
| Right-click → **Open scan-axis cut**, viewer button **Scan-axis cut** | Open deflector cut, Deflector cut |
| Memory box → **Free RAM**, **Close all**, **Autosaves** (one row) | under the list; then Free memory, Close all viewers, Open auto-save folder |
| Datasets box → **Experiment info** | right-click → Show information |
| Datasets box → **Rename**, **Save**, **Remove** (one row), **Session log...** | right-click; then Rename..., Save datasets..., Remove from list |
| **Combine...** under the list (only) | also right-click → Combine N datasets |

**Experiment info** also opens the full metadata of the dataset in a window
of its own, with a filter box and **Export metadata (.csv)**:

![Dataset information](images/information_window.png)

### The experiment log

Two buttons side by side in the **Datasets** box write the beamtime's log --
a CSV file, one row per measurement -- from the files' own metadata:
**Log format...** (left) sets it up once, **Auto Exp. log** (right) adds the
selected datasets to it with one click.

![The Log format window](images/exp_log_dialog.png)

**Log format...**

1. **Log file.** Choose the folder and type the file name. If a log of that
   name already exists, the window shows **that log's own format**, ready
   to adjust.
2. **Format.** *Separator* and *Decimal mark* are by default what Excel
   uses on this computer. Excel splits a CSV at the system's list
   separator -- a **semicolon** wherever the decimal mark is a comma (French,
   German, Italian... settings) -- so a comma-separated file opens there
   with everything in the first column. Leave them on *auto* unless the log
   is for a computer set up differently. *Ranges as:* writes a range
   `start/stop/step` or `start:step:stop`.
3. **Columns.** *Fields of:* picks a dataset of the main list -- choose one
   of the kind you log most -- and **Read fields** lists everything it has,
   with its value: the file name, data type, start time and scanned ranges
   worked out from the data, then every metadata field and motor. **Tick**
   a field to add it as a column. In the table on the right, double-click a
   heading to rename it (`hν/Δ`, `Polar/Δ`, `Epass`...), and use **Top /
   Up / Down / Bottom** to put the columns in any order. **Add empty
   column** gives one to fill in by hand (comments). **Save format as...** /
   **Load format...** keep one format per beamline.
4. **If scanned.** What a position or angle column says depends on the
   measurement:
   - a **cut**: the sample stood still -- the motor's value;
   - a **SPEM** scan: the scanned stage axis as `start/stop/step`
     (choose *SPEM: scanned X* / *Y* for the stage columns);
   - a **map** taken with the manipulator rather than the deflector: the
     angle as `start/stop/step` (*map: scanned angle*), and a sample
     position the beamline corrected along with the angle as
     `first..last`;
   - a **kz map**: the photon energy range (*kz map: photon energy*).

   *auto* (the default) recognises the scanned motor by name, and a
   position that moved by its values. Any other measurement gets the
   field's own value.
5. **Save format.** If the log already holds records and you changed its
   format, it is **rewritten under the new header**: each record moves with
   its columns (renamed or moved columns are followed), new columns are
   empty for the old records, and columns you took out are named first,
   since what they hold goes. The header stays the only first row.

*Skip datasets already in the log* leaves out files whose name is already in
the log's file-name column. **Keep recording** logs every file you add to the
list from then on, until the program is closed: during a beamtime, add each
file as it is measured and the log keeps itself.

**Auto Exp. log** adds the selected rows to the log, after its last record,
without opening anything; the status bar says how many were added. The very
first time, it opens *Log format...* first.

If the log is open in Excel, Windows will not let it be written: close it
there and click again.

### Combining datasets: stacking and joining

**Combine...** (under the list) does two different things, and asks which:

- **Stack** — *cuts only*. Each cut becomes one slice of a new map, along an
  axis you name: a temperature series, a gate-voltage series, one cut per
  photon energy. The result is looked at like any map — sliced at constant
  energy, or along the new axis.
- **Join** — *cuts, maps, kz maps or k-maps, all of one kind*. The pieces of
  one measurement that was taken in segments are put back end to end along
  one of their own axes: a map whose deflector range was scanned in two
  parts, a kz map whose photon-energy range was split over several runs, a
  cut whose energy window was taken in two halves.

**Stacking cuts**

1. Select the cuts in the list (Ctrl/Shift-click) and press
   **Combine...**. Stacking needs the same axes point for point (same number
   of points, same values, same units); if they differ, *Stack* is greyed
   with the first difference named, and *Join* is offered instead.
2. It asks whether to combine them **in the list's order**. *Yes* keeps that
   order; *No* lets you choose it.
3. Choose **Stack**, then set the new axis: **New axis is** (temperature,
   gate voltage, photon energy, delay, position, other — this fills in
   *Name* and *Unit*, which you can change), **Start** and **Step**. Each
   cut is listed with the value it gets. If you chose to reorder, **drag**
   the rows (or use **Move up** / **Move down** / **Reverse**). Press
   **Combine**.

![Combining cuts](images/combine_cuts_dialog.png)

A photon-energy series of angle cuts becomes a `kz map`, so *Align Fermi
levels (kz map)* and *Convert to momentum (kz, k∥)* take it (see
[§14](#14-photon-energy-scans-kz)); any other axis gives a `map` whose first
axis is marked as not an angle, so no k conversion is offered for it.

**Joining pieces of a measurement**

1. Select the pieces (all maps, or all kz maps, or all cuts …) and press
   **Combine...**. Maps and kz maps go straight to *Join*.
2. **Join along**: the axis the measurement was split along. It starts on
   the axis whose range differs between the pieces — a map's deflector
   axis, a kz map's photon energy, or the energy of cuts taken in two
   energy windows — so it is usually right already.
3. **Where pieces overlap**: *Average them*, or keep the earlier / later
   piece in the order shown (drag the rows to change it).
4. The dialog lists each piece's range and shows the **Result** — the joined
   range and the number of points — with notes on anything it had to do.
   Press **Combine**.

![Joining a map measured in two parts](images/combine_dialog.png)

How the axes are handled:

- Pieces on one common grid (same step, ranges that meet or overlap by whole
  points — the usual case) are **copied as measured**, nothing
  interpolated.
- A piece with a different step or a half-step offset is **interpolated**
  onto the common axis (the first piece's step), and the dialog says so.
- If the **other axes** differ — an energy window a few meV different
  between runs, a different number of slit channels — every piece is
  resampled onto the range they **all** cover, and the dialog says so.
- A stretch that no piece covers is a **gap**, left at zero and reported.
- A joined kz map is still a `kz map`: run *Align Fermi levels (kz map)*
  on the joined result rather than on the pieces (per-photon-energy fit
  results of the pieces are not carried over).

### Summing repeated scans

For a fine cut or map taken as many short scans so that an ageing sample can
be caught. Select the scans (two or more of one kind) → **Sum repeated
scans...**:

1. **Are these repeats?** Each scan's grid and recorded conditions (photon
   energy, polarisation, angles, position, temperature, pass energy, lens
   mode); orange marks a value that differs from the first scan.
2. **Has the sample changed?** Press **Compare the scans**. Each scan is
   compared with the sum of the others: *χ² / point* is about the same for
   scans that differ only by counting noise and grows when the spectrum
   changed; scans well above the rest are marked **changed**. The energy
   shift of each scan's EDC and its count rate are listed; click a row to
   see its EDC against the rest.
3. **Sum.** Tick the scans to keep (changed ones start unticked), choose
   whether to undo the energy shifts, **Sum** or **Mean**, and **Sum to
   list**.

### Step-by-step windows

The tools that take several decisions open as **numbered steps** (sum
repeated scans, detector corrections, tight-binding fit, crystal structure,
calculated bands): the title says *Step 2 of 4: …*, a blue box says what
the step is for and what to do, **Next →** stays greyed (with the reason
beside it) until the step is done, and **← Back** or the list of steps on
the left goes back.

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
read from the `_i` parameter files next to the spectra. The row is named
with the photon energy and polarisation read from those files —
`FineMap 117eV LH`, or `KzScan 40-120eV LH` for a photon-energy scan — so
repeats of the same scan at another energy or polarisation can be told
apart.

**ANTARES photon-energy (PhD) files.** A PhD `.nxs` file (one `Phd_NNNN`
entry per photon energy) is listed as one `kz map`: the frames are sorted by
photon energy and each kinetic-energy axis is referred to E_F with its own
hν (work function 4.5 eV, which the k_z conversion then uses). Each frame
still opens on its own as a cut.

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

## 5. What the viewers have in common

Each dataset opens in its own window, so you can put a map, its cuts and
another file side by side.

![A cut viewer with the readout cursor](images/cut_viewer.png)

- **The top row**: the **Functions** menu (everything the window can do),
  **Plot tools** (guide lines and labels on this plot — see below), the
  colormap (per window) with **Flip**, and the view controls — axis ranges,
  **Auto**, **Inv** (invert), **Reset**, **Grid**.
- **Above each image**: the **Levels** bar (drag the two handles, or type
  **Min%/Max%** percentiles; **Clip** applies them) and **Gamma**. A min/max
  you set is kept while you move the slice (the energy of a contour, the
  angle of a cut, the cursor of a spatial map), so every position is on the
  same colour scale; drag the handles back to both ends, or set 0%–100%, to
  follow each frame's own range again.
- **Under each image**: **Smooth** (display only, does not change the data)
  and **Lock x/y** (aspect ratio).
- **Readout cursor** — right-click the image → *Readout cursor*. It reports
  the value under it and, on a cut, draws the **EDC** (right) and **MDC**
  (below) through that point. **EDC ±** / **MDC ±** integrate over a window,
  shown as a shaded band on the image. Type a position into the **Cursor**
  row and press Enter to put it exactly there.
- **EDC → list** / **MDC → list** put the curve on screen in the main list as
  a dataset of its own (see [§23](#23-curves-edcs-mdcs-and-spin-resolved-data)).
- **Right-click the image → Export** saves what the panel shows — data
  (`.nxs`), an image, or vector axes (`.eps`); see [§28](#28-the-figure-composer).
- **Functions → Slice → Save slice to the main list / Freeze slice in a new
  window** keeps a copy of the slice on screen. The frozen window (and *Open
  this position in a new window* on a spatial map) opens with the panel's
  min/max and gamma, so the two can be compared directly. It has its own
  **Plot tools**, and a **Figure composer...** button to make a figure of it
  ([§28](#28-the-figure-composer)).

The **Functions** menu now holds only what works on the window itself —
the tools that take a dataset are in the launcher's tool panel
([§3](#3-the-main-window)):

| Viewer | Functions |
|---|---|
| map / kz map / k-map (contour) | Arbitrary cut; Slice: save / freeze the slice |
| cut | Straighten bands (FS correction); Slice: save / freeze |
| a map's slit cut | Straighten bands (FS correction) for the whole cube; Show the map |
| curve | Curve fit, Spin analysis; Crop, Bin, Normalise, Background, Shift axis, Add √N errors |

An entry that applies but cannot run yet is greyed, with the reason in its
tooltip.

**Brillouin zone tick box.** A k-map or converted kz map on which a zone has
been plotted ([§15](#15-brillouin-zone)) keeps it: its viewer has a
**Brillouin zone** tick box on the top row that shows or hides it.

### Plot tools: guide lines and labels on the plot

**Plot tools** on the top row of every viewer (cuts, maps, spatial scans,
curves, frozen slices) opens a small window that draws **on the plot in
front of you**:

![Guide lines on a cut](images/plot_tools_viewer.png)

| ![Plot tools](images/plot_tools.png) | |
|---|---|

- **Type**: *Vertical line*, *Horizontal line*, *Line at an angle* (through a
  point; the angle is in the plot's own units, i.e. the slope is tan(angle)
  y-units per x-unit), or *Text label*.
- **Position** (x, y), **Angle**, an optional **Label** (shown on the line,
  e.g. `E_F`, `Γ`), a **Colour** (nine one-click colours, or **More...**),
  **Width**, **Style** (solid, dashed, dotted, dash-dot) and **Text size**.
- **Add** draws it at the position typed (**Centre of view** fills that in);
  **Place by clicking** draws one wherever you click on the plot.
- **Drag** a line on the plot to move it; the window follows. Select a guide
  in the list to **edit** it — every change shows at once — or **Remove** it;
  **New guide** stops editing; **Remove all** clears the plot.
- With two panels (a spatial scan), **Draw on** chooses which.
- Guides stay while the slice moves (an E_F line stays at E_F as you scroll
  through a map) and while the Plot tools window is closed. The panel's
  **Export → Image** draws them too (*Include the Plot tools guides*).

*Plot tools* are for looking at the data. To make a figure for a paper, use
the **figure composer** and its **Figure tools** ([§28](#28-the-figure-composer)).

## 6. Maps: the contour and its two cuts

A map opens on its **constant-energy contour**. The slider at the top picks
the energy (**Pos** in eV, or **Ind** as an index; **±** integrates over a
window, **mean** or **sum**). **Scan-axis cut** and **Slit cut** open the two
orthogonal cuts, each in its own window. (The scan axis is whatever was
scanned to make the map — usually the deflector angle, but a manipulator
rotation or, for a kz map, the photon energy — so the button is named for
that.)

![Constant-energy contour of a map](images/map_contour.png)

| Slit cut (slit angle vs energy) | Scan-axis cut (deflector angle vs energy) |
|---|---|
| ![Slit cut](images/map_slit_cut.png) | ![Scan-axis cut](images/map_scan_axis_cut.png) |

**One cursor for all three windows.** Switching the readout cursor on in any
of them switches it on in the others, and moving it anywhere re-slices the
others through the new point. The axes keep one colour each everywhere:
**deflector red, slit green, energy blue**. A cursor line has the colour of
the axis it holds constant, and the EDC or MDC it feeds is drawn in the same
colour. Integration widths (the ± boxes) belong to each window and are never
shared.

From the list you can also open a map's cut straight away (right-click →
**Open slit cut** / **Open scan-axis cut**); the contour then stays hidden
until you ask for it (*Functions → Show the map*).

## 7. Arbitrary cuts and Fermi-surface correction

### 7.1 Arbitrary cut

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

### 7.2 Fermi-surface (FS) correction

*Functions → Straighten bands (FS correction)...* on a cut (or a map's slit cut) straightens a
feature that should be flat — usually the Fermi edge, which curves across the
analyser angle. Click along the feature; a polynomial (order 2 by default) is
fitted and every angle column is shifted in energy. **From a reference
(Au)...** takes the curve from a gold spectrum in the list instead. On a
map's slit cut the same shift is applied to the whole cube.

![FS correction](images/fs_correction_dialog.png)

# Part II — The Basic tab

The **Basic** tab holds what nearly every dataset goes through, grouped by
the data each tool takes: *2D & 3D* — Fermi edge fit, detector-grid
removal, detector corrections; *2D only* — background subtraction, the cut's
k conversion, the core-level check; *3D only* — the map's k conversion,
aligning and converting a photon-energy scan, the Brillouin zone.

---

## 8. Fermi edge fit

The usual order for a single cut is: **Fermi edge fit → offset the energy axis →
convert to momentum → MDC/EDC fit**. Each step is a button of the tool panel:
select the cut in the list and press it.

Tool panel → Basic → **Fermi edge fit...** fits a Fermi edge (Fermi–Dirac ×
linear DOS, convolved with a Gaussian resolution) to the EDC summed over a
box. It opens the cut's viewer together with the fit window; the two belong
together, and closing either closes the other.

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

**On a map.** Select a map (angle map or k-map) and press the same button:
the map's **slit cut** opens with the fit window, and the slit cut's slider
picks which cut is fitted. After a fit, **Offset energy axis** asks whether
to shift **the whole map by this cut's E_F** (*Yes*) or to fit **one E_F per
cut** (*No*). With one per cut, the **Fitted cuts** list on the right fills
as you go: move the slider, press **Fit**, repeat. A second fit at the same
place replaces the first, and every E_F in the list can be edited by hand.
Cuts between two fitted ones are shifted by a straight line between their
neighbours (beyond the first and the last, by the nearest fit), shown in the
plot under the list. **Offset energy axis by these fits** lists the
corrected map.

![Fermi fit on a map, one E_F per cut](images/fermi_map_cuts.png)

## 9. Removing the detector grid and detector corrections

### 9.1 Removing the detector grid

MCP (hexagonal) and mesh (square) patterns print on every image. **Do it
first**: only data still on the detector's pixels can be de-gridded, so the
button refuses anything already converted, corrected or smoothed, and says
why.

Select a map (or kz map) or a cut and press tool panel → Basic → **Remove
detector grid...**. For a map, no reference is needed — the grid stays on the same pixels in
every slice while the photoemission moves. **Find and remove the grid**:

![De-grid](images/degrid_dialog.png)

- **Before / After** a chosen slice on one colour scale (tick **Show
  −∂²I/∂E²** to make any grid obvious).
- **The grid found**, its **k-space regions** (in red), and the **contrast
  and shift** fitted for every slice.
- The power left in the grid regions (1 = no grid left).
- **Also list the grid pattern** keeps a `[grid]` dataset, which a cut taken
  later with the same lens mode and pass energy will find and use.

### 9.2 Detector corrections

On a cut, map or kz map as measured:

1. **Keystone** -- tick *Correct a trapezoidal image* and drag the four
   yellow corners onto the lit area; the right picture shows it made
   rectangular, each line's counts kept.
2. **Bad regions** -- **Add a region** and drag it over a dead or hot spot,
   or **Find bad pixels**; fill *smooth* (best for a spot) or *linear
   along* an axis (for a stripe).
3. **Apply** -- **Apply to list**. Do the grid removal, Fermi level and k
   conversion on the result.

## 10. Background subtraction

Tool panel → Basic → **Background subtraction...** (cuts) removes the
intensity that is not the bands: a **Shirley** or **Tougaard** inelastic
background, a floor common to every angle (a percentile of each energy
row), a smooth **polynomial** — or it divides out the **Fermi cut-off**
(Fermi–Dirac at the temperature and resolution the file records, starting
at the E_F it records or else at the steepest fall of the angle-integrated
EDC). The window is the same as the Advanced tab's processing windows —
source and result side by side, live preview, **Apply** lists the result —
described in [§16](#16-smoothing-derivatives-curvature-and-symmetrising).

## 11. Converting a cut to momentum

Tool panel → Basic → **k conversion (cut)...** converts the angle axis to k∥
(Å⁻¹). A **vertical line** is drawn on the cut: drag it to the slit angle of
Γ. The cut is then taken to run through Γ — not exact if it misses Γ, but
the only choice when only this cut was measured.

![The Γ line on the cut](images/cut_k_line.png)

![Cut k conversion](images/cut_k_dialog.png)

- **Advanced: Γ off this cut** is for a cut that does not run through Γ: type
  the deflector and slit angles of Γ, or **From a converted map...** to
  inherit them from a map of the same sample already converted.
- If the cut is already on E − E_F, **Energy offset** is filled in with the
  E_F that was subtracted (the note under it says so): the conversion uses
  true kinetic energies, and the result stays on E − E_F.
- **Remembered settings.** The dialog starts from the last conversion of a
  cut of the same size (same number of points on each axis) — the usual case
  for a series of cuts of one sample. The map and kz conversions do the same.

![The cut in momentum](images/cut_in_k.png)

## 12. Core-level check

A wide scan of the core levels is listed as a cut too. Select it and press
tool panel → Basic → **Core-level check...** to find out what its peaks are.

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

**Weak lines and peaks by hand.** Tick **Log counts** to see weak lines next
to a strong one; **Min height** (default 0.2 % of the strongest) sets how
small a found peak may be. To add a peak the search missed, press **Click to
add peak** and click the plot (a double-click always works), or type its
binding energy after *or at* and press **Add**: its position is refined by a
local fit, it is marked ◆, and it is named with every tabulated level within
±2 eV. **Remove selected** drops peaks; **Reset peaks** brings back the found
ones and forgets those added by hand. Assignments you typed survive a
**Re-analyse**.

**Peak table...** lists every tabulated level, filtered by element and
binding-energy range: look a peak up, show a level on the plot, or add it as
a peak with that level as its assignment. **Select range** shows a band you
drag over a region; **Zoom window...** opens it (or the current view) in a
window of its own with a cursor, a Gaussian fit at the cursor, and the
negative second derivative divided by the intensity — which shows a shoulder
or a weak line beside a strong one as clearly as the strong one.

If every peak is off by the same amount, the notes say so and **Apply
suggested shift** corrects it (a charged sample, or a wrong E_F). **Add EDC
to list** lists the EDC on the binding-energy scale; **Copy report** copies
the whole result as text. Measuring at a second photon energy settles what
stays doubtful: core levels stay put in binding energy, Auger lines and
higher-order light move.

## 13. Converting a map to momentum

Select an angle map and press tool panel → Basic → **k conversion (map)...**.
The map opens on a bare constant-energy contour — no cuts, no other
functions, just the picture — next to the conversion window; the two close
together. Drag the readout cursor to the point that should become
**k = (0, 0)** — the dialog follows it — then **Convert**. The result is a
new `k-map` row in the list.

| The contour for picking Γ | The conversion |
|---|---|
| ![Map conversion contour](images/map_k_contour.png) | ![Map k conversion](images/map_k_dialog.png) |

- **Γ deflector angle / Γ slit angle**: the angles of the point that becomes
  k = (0, 0), filled in from the cursor.
- **Sample rotation** turns the map about the surface normal; **Set rotation
  from contour...** reads it off a direction you click on the contour (two
  points along a mirror line, say).
- **Output grid**: number of points, or resolution (Å⁻¹, meV). The
  numbers of points start at the measured ones — as many kx points as
  deflector angles, ky points as slit angles, E points as energies.
- **Energy offset** is added to the energy axis before converting, for an
  axis that is not a true kinetic energy (filled in automatically for a map on
  E − E_F).
- The settings start from the last conversion of a map of the same size, so
  a second map of a series needs no retyping.

![The converted k-map](images/k_map_contour.png)

## 14. Photon-energy scans (kz)

A photon-energy scan is listed as a `kz map`; its contour is photon energy
against angle at one energy.

![A kz map](images/kz_map_contour.png)

### 14.1 Aligning the Fermi levels of a kz map

Nothing in the files says where each spectrum's Fermi level really is.
Tool panel → Basic → **Align Fermi levels (kz)...** measures it from the
spectra themselves; the kz map's slit cut opens with it, and the two close
together.

1. The slit cut opens with a box on it. Drag the box over a Fermi edge — the
   angle range sets the channels summed, the energy range the fit window.

   ![The box on the slit cut](images/kz_processing_cut.png)

2. **Fit and assemble**: the same box (by detector index) is used for every
   spectrum; each is shifted to put its own edge at zero; the cube is cropped
   to the range they all cover; and, by default, each spectrum is
   **normalised by its total intensity**, which divides out the photon flux.

![Align Fermi levels (kz map)](images/kz_processing_dialog.png)

E_F against photon energy is plotted. A smooth curve means the fits are
sound; points whose fit is not believed are marked and filled in from their
neighbours.

### 14.2 Converting a kz map to momentum

Tool panel → Basic → **kz conversion...** converts the scan to (k_z, k∥). It
starts from the last conversion of a scan of the same size.

![Convert to momentum (kz, k∥)](images/kz_conversion_dialog.png)

- **Inner potential V0** is the unknown: the preview reconverts one energy
  slice as you change it, in milliseconds.
- With the **lattice** (space group, a, b, c), the **surface normal** and the
  in-plane direction, the zone boundaries are drawn over the preview.
- **Scan V0** converts at a range of V0 values and reports where the k_z
  period matches the lattice — and how well that pins V0.
- **Measure a k_z period...**: click two points that are the same feature one
  zone apart; the matching reciprocal-lattice planes are listed.
- The yellow banner warns when the Fermi edge still bends across the analyser
  angle — correct it first (Straighten bands (FS correction) on the slit cut).
- **Convert the whole cube** lists a `kz map (k)`.

![The converted kz map](images/kz_map_in_k.png)

## 15. Brillouin zone

Select a k-map (or a converted kz map) and press tool panel → Basic →
**Brillouin zone...**. Like the map conversion, it opens a contour of its
own with the zone window, and the two close together.

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

The window is tall; it scrolls, and can be resized like every other window.

**The zone stays with the dataset.** Once plotted, it is a property of the
k-map: it is saved with it into a `.nxs` file, and its viewer shows it with
the **Brillouin zone** tick box on the top row. **Clear overlay** removes it.

![k-map with the zone shown by its tick box](images/k_map_with_bz.png)

Its exports can draw it too: **Export → Image...** has **Include the
Brillouin zone**, and **Axes only (.eps)...** has **Draw the Brillouin zone in
the box**, so the zone comes out as vector lines over the image in
Illustrator.

![Export with the zone](images/export_dialog_bz.png)

# Part III — The Advanced tab

The **Advanced** tab holds the analysis that follows: the figure composer
(*1D & 2D*, see [§28](#28-the-figure-composer)); smoothing, the band-sharpening
filters, symmetrising and ridge detection (*2D & 3D*); stack plots, MDC and
EDC fits and cut arithmetic (*2D only*); the slice series figure and the
four 3-D views (*3D only*). Curves and spatial scans are looked at in
viewers of their own.

---

## 16. Smoothing, derivatives, curvature and symmetrising

Each filter has a window of its own, opened from the tool panel with one or
more cuts (or maps) selected:

| Tab | Button | What it does |
|---|---|---|
| Basic | **Background subtraction...** | removes the intensity that is not the bands: Shirley or Tougaard inelastic background, a floor common to every angle, a smooth polynomial — or divides out the Fermi cut-off |
| Advanced | **Smooth...** | Gaussian, Savitzky–Golay or moving average, with optional spike removal |
| Advanced | **Second derivative...** | −d²I/dE² or −d²I/dk² after a pre-smooth, to sharpen the bands |
| Advanced | **Curvature...** | 1-D or 2-D curvature (Zhang et al.), which follows the band better than a second derivative |
| Advanced | **Symmetrise...** | mirror a cut about a line; for a map, mirror or rotate its constant-energy planes |

Every result is a new dataset carrying the steps that made it.

![The curvature window](images/process_dialog.png)

- **Source** and **Result** are shown side by side on the same scale, so you
  can see when a setting has eaten the band. They pan and zoom together.
- **Smooth first.** The second-derivative and curvature windows open with a
  tip: both amplify noise, so they usually work better on data smoothed
  first — run **Smooth...**, select its result in the list, and open the
  filter on that, or raise the pre-smooth widths in the window itself.
  A smoothed result carries the tag **‹smoothed›** in the list (and so does
  anything made from it); on such data the tip is replaced by a short
  "already smoothed" line and the pre-smooth widths start at 0.
- **Shape**: *Square pixels* (one data point is square — undistorted for an
  angle-vs-energy cut), *Equal axes* (one unit is the same length on both
  axes — right for angle-vs-angle or k-vs-k; chosen automatically when both
  axes have the same unit), or *Fill the box*.
- **Live preview** re-computes as you change a value; **Apply** lists the
  result.
- Starting values: **Symmetrise** mirrors about a vertical line through
  angle (or k) = 0, i.e. k → −k; the **Fermi cut-off** under *Background
  subtraction* starts at the E_F the file records, or else at the steepest
  fall of the angle-integrated EDC, with the file's temperature and
  resolution.

For a map the window works plane by plane, and asks **which planes** first —
constant-energy maps or dispersions — since the curvature of a contour and
the curvature of a cut are different quantities. **Preview at** picks the
plane shown; **Apply to the whole cube** runs it on every plane. The
operations and widths are named after the chosen plane's own axes
(*Second derivative along Kinetic energy*…), and the smoothing widths start
at 1/100 of each axis.

![Smoothing a map](images/volume_process_dialog.png)

**Truncate...**, **Compress...** and **Self-normalise...** (under the list)
change the grid or the scale of one or several datasets at once — as long as
they have the same kind and the same number of points on every axis:
**Truncate** crops to a range, **Compress** bins points together,
**Self-normalise** divides each spectrum by its own total.

**Self-normalise** divides each piece of the data by its own total (or peak):

- **Normalise** picks the piece in ARPES terms, and starts on the usual
  one: *each EDC* for a cut, *each cut (slit × energy at one scan step)*
  for a map, k-map or kz map — this takes a flux or beam-current drift
  along the scan out. Other choices: each MDC, each EDC of a map, each
  constant-energy map, or *custom*.
- **Summed over** shows the axes this sums over; the ticks can be changed by
  hand, and the sentence under them says what the current ticks do ("Each
  plane over Angle, slit × Kinetic energy is divided by its own total — one
  factor for every Angle, deflector").
- **Sum only over** restricts the sum to a window on one of the summed axes
  — the energy by default — for example the energies of a core level.

![Self-normalise](images/data_operations.png)

## 17. Ridge detection

Tool panel → Advanced → **Ridge detection...** (cuts or maps) traces bands
as ridges of the intensity, each point to a fraction of a pixel, which is
sharper than a second derivative and draws no false double lines.

- **Smoothing** along each axis, in its own unit — about the band's width.
- **Threshold**: raise it until only the bands are left (yellow points on
  the source).
- **Result to list** — the ridge image (a map plane by plane, on the planes
  chosen); **Save ridge points...** — a CSV that the tight-binding fit
  ([§26](#26-tight-binding-fit)) reads.

## 18. Stack plots

**Stack plot (EDC / MDC)...** (tool panel → Advanced, with a cut selected): a
waterfall of EDCs or MDCs, with the number of curves, binning,
normalisation, offset, shear and colouring controls; **Send to figure
composer** makes it a figure panel.

![Stack plot](images/stack_plot.png)

## 19. MDC fit and EDC fit

Two tools of the tool panel's Advanced tab fit peaks line by line:

- **MDC fit (k-cut only)...** fits the momentum distribution curves and reads
  the band off the fitted centres: Fermi velocity, effective mass,
  self-energy. It only takes a cut whose axis is already in Å⁻¹ — a velocity
  in eV·deg would be meaningless — and its button is greyed otherwise.
- **EDC fit...** fits the energy distribution curves, on any cut, optionally
  with the Fermi cut-off in the model. v_F and m\* from the fitted band need a
  cut in momentum.

![MDC fit panel](images/mdc_fit_panel.png)

**Guided fit (step by step)**, at the top of the panel, walks through the
whole fit beside the controls, one step at a time, saying what each step is
for and outlining the part of the panel it uses: what an MDC (EDC) fit
does, the line shape, the range (*Start here* / *End here* at the blue
line), the bands, the seeds (the guide moves the line to five places in the
range — click the band's peak on the curve at each), one test line, the
series, and the result. *Next* waits until a step has what the fit needs.

The pictures and the controls share the window with a handle between them,
so the settings are always in sight; the mouse wheel over a picture scrolls
the window when it is too short for everything (Ctrl + wheel zooms the
picture).

By hand:

1. **Bands** tab: **Add band** for each branch; select a band, press **Place
   a seed** and click the peak on the curve (right). Four or five seeds down a
   band are enough; **Line at** moves the line being looked at.
2. **Try this line** shows the fit of the current line (dashed: the guess).
3. **Fit settings** tab: peak shape (Lorentzian, Gaussian, Voigt),
   background, step, then **Fit the series**. Fitted positions are drawn on
   the image, sized by amplitude.
4. **Dispersion...** fits the band: **v_F** (straight line), **m\*** (parabola)
   — and checks how far out the fit is stable. **Self-energy...** gives Re Σ
   and Im Σ.

![Dispersion window](images/dispersion_window.png)

**Save table...** writes every fitted line; **Band to list** puts a band in
the main list.

**Where two straight branches cross.** In the **Dispersion...** window, with
*v_F* (a straight line) chosen, each other fitted band is fitted the same way
and the point where the two lines meet is given with its uncertainty
(propagated from both fits): the Dirac point of a cone fitted as its two
sides, or any crossing of two bands. A single fitted line whose error
collapsed (a peak fit that went to zero width) no longer takes the whole
weight of the dispersion fit: no point counts for more than ten typical
ones.

## 20. Cut arithmetic: dichroism and reference division

Two cuts combined into a third: linear dichroism (LH − LV), circular
dichroism, dividing by a gold reference, or plain differences, ratios and
sums. Select two cuts in the list and press tool panel → Advanced → **Cut
arithmetic (two cuts)...** (the first selected is A; the window can swap
them). The button is greyed until exactly two cuts are selected.

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

## 21. Slice series figure

Select a map → tool panel → Advanced → **Slice series figure...** lays out a
page of flat 2-D slices through the cube as one figure for a paper —
constant-energy contours at a series of energies, or cuts at a series of
angles or momenta — evenly spaced or at chosen values, with integration,
one shared colour scale, and labels on the outer edges only. (Not to be
confused with **3D slices**, which draws three planes in a 3-D view.)

| Settings | Result |
|---|---|
| ![Slice series settings](images/slice_figure_dialog.png) | ![Slice series figure](images/slice_figure.png) |

## 22. The 3-D views

Four tools of the Advanced tab show a map in 3-D, each in a window of its
own — drag the picture to turn it:

- **3D slices...** — the three orthogonal planes through the cube, each
  placed where you choose (tick box, slider or typed position per axis).
- **Notched cube...** — the cube with one corner cut away: a constant-energy
  map and two dispersions in one picture.
- **Volume projection...** — the whole volume projected along the view:
  maximum intensity, sum, alpha composite, or a shaded isosurface.
- **Stacked slices...** — parallel slices stacked in 3-D (below).

| 3D slices | Notched cube |
|---|---|
| ![3D slices](images/volume_slices.png) | ![Notched cube](images/volume_view.png) |

- **Sampling.** Every frame resamples the cube, so each window opens at
  *Balanced* sampling, which keeps up with the mouse. The current sampling is
  written under **Save the view...**; **Change sampling...** picks the region
  and *Fast*, *Balanced* or *Fine*, and redraws from the full cube.

  ![Change sampling](images/volume_range_dialog.png)

- **Colour.** *Apply to* sets the colours of **all slices** at once or of
  one orientation (the energy plane in one colormap, the dispersions in
  another): colormap (and *reversed*), levels (**Auto**, or **Each its own**
  for one scale per slice), gamma, **Opacity** and **Fade low values**, which
  makes the weak intensity see-through so the planes behind show.
- **Axes.** **Draw the axes, with ticks and values**. *Placement* is
  automatic — on the outer edges facing you, following the view as it turns —
  or chosen: for each axis, the edge of the cube it runs along and whether it
  counts **from min** or **from max**.
- **Camera**: azimuth, elevation, depth aspect; **Look into it** faces a
  chosen plane. **Outline the faces** and an **Orientation tripod** are
  optional.
- **Save the view...** exports the picture through the same export window as
  any panel (PNG, TIFF, JPEG at a chosen size and dpi). **To figure
  composer...** adds the view as a panel of a figure.

### Stacked slices

**Stacked slices...** (tool panel → Advanced → 3D only) draws parallel
slices of the cube stacked in 3-D: constant-energy maps one above another
up the energy axis, or cuts one behind another along k_z (or the photon
energy) — the usual picture of how a band evolves.

![Stacked slices](images/stacked_slices.png)

- **Layers** — the axis they are stacked along; **How many** spreads them
  evenly, or type the positions in **At** (comma-separated, in the axis's
  units); **Integrate** averages each layer over that width.
- **Frames** — an outline round every layer (grey, black or none).
- **Spread** stretches the box along the stacking axis, so that layers seen
  from an oblique view stand apart instead of covering one another: for a
  few cuts stacked along k_z, *Spread* 3 with the azimuth at about 40° shows
  each cut almost whole.
- **Fade low values** (*Colour*) starts at 35 % here: each layer's background
  is transparent, so the layers behind show through.
- **Mark a curve through the layers** — points in the data's own units (one
  per layer, say the crossing of two bands), joined by a red, black, blue…
  curve drawn over the stack, with an optional label at the first point.
- The axes of all 3-D views now fit the picture *and* its labels to the
  window or figure panel (no fixed margin), thin their ticks to the room
  each axis has, and write sub- and superscripts in the axis titles
  (`E - E_F`, `k_x`).

## 23. Curves: EDCs, MDCs and spin-resolved data

**EDC → list** / **MDC → list** under any image, or a spin-resolved file,
give a curve dataset, which opens in the **curve viewer**:

![Curve viewer](images/curve_viewer.png)

Channels can be hidden, offset into a waterfall (**Offset**), shown with
**Error bars** (from a σ channel), and **Range** marks an x window used by
the operations. *Functions → Transform* gives **Crop**, **Bin**,
**Normalise**, **Subtract a background** (constant, linear, Shirley),
**Shift the axis** and **Add counting errors (√N)**; each makes a new dataset
and propagates σ correctly.

### 23.1 Curve fit

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

### 23.2 Spin analysis

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
- **P to list**, **I↑ / I↓ to list**, **As a figure** (opens the figure composer).

## 24. Spatially resolved scans (nano-ARPES)

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

# Part IV — The Calculation tab

The **Calculation** tab holds what is worked out rather than measured: the
crystal in real space and its Brillouin zone, band models fitted to the
data, and band structures calculated by others. Its windows open as
numbered steps.

![The Calculation tab of the tool panel](images/main_window_calculation.png)

---

## 25. Crystal structure

Tool panel → Calculation → **Crystal structure...** (needs nothing
selected; selected k-maps can be drawn on the Brillouin zone's planes).
Five steps:

![The crystal structure](images/crystal_structure.png)

1. **The structure** — type it in (cell, space group for its centering,
   atoms), take the Brillouin-zone tool's last crystal, or **Read a CIF
   file** (its symmetry operations are applied, its space group read).
2. **The look**
   - cells along a, b, c, the atoms' overall size, cell edges and bonds
     (in the atoms' colours, or *grey, thin*);
   - per element its **colour**, **size**, **transparency** and **gloss**
     (0 %: a matt ball, 100 %: a shiny one with a sharp highlight);
   - **Polyhedra around** an element: its coordination polyhedra (its
     bonded neighbours as corners) in a colour and transparency of their
     own — the TiSe₆ octahedra, the CuS₄ tetrahedra;
   - **Key of the elements** (a ball and its name for each), **Lattice
     constants on the cell edges** (`a = 3.52 Å`, `c = 8 Å`).

   Drag the picture to turn it, wheel to zoom; *View along* a, b, c or the
   plane normal.
3. **Cleavage plane** — tick *Draw a plane*, give (hkl), slide it, or
   double-click an **atomic layer** in the list to put it there; *Hide the
   atoms above* shows the termination; atoms in the plane are ringed red.
4. **Brillouin zone** — the crystal's zone, its labelled high-symmetry
   points and a **path** through them (drawn in red, e.g. `Γ-M-K-Γ-A-L-H-A`
   or one segment, `K-H`), the k_x, k_y, k_z axes as arrows, and up to two
   **planes** normal to (hkl) — through Γ and at the zone top by default,
   each in a colour of its own.
   - **map** on a plane: a selected k-map's constant-energy map at the
     energy given (averaged over *integrate*), drawn on that plane and cut
     to the zone — the measured Fermi surface in the zone where it was
     measured. *azimuth* turns the maps against the zone's k_x.
   - **Label on the left**: the points named here (`L, M`) are labelled at
     their leftmost symmetry-equivalent point as the picture turns, clear of
     the path and the axes.

   ![The zone with measured maps on its planes](images/crystal_zone.png)
5. **Export** — the crystal or the zone: **Save picture...** (PNG at any
   size, drawn again at that size) or **To figure composer...**. A picture
   from a CIF names the CIF as its source, so the figure's data folder keeps
   it ([§28](#28-the-figure-composer)).

## 26. Tight-binding fit

With a k-map selected (or nothing, to use a file of points):

1. **The lattice** -- a 2-D layer (a, b, γ) or a 3-D crystal, or **Take the
   crystal from the Brillouin-zone tool**; the azimuth of a against kx; how
   many shells of neighbours.
2. **Band points** -- move to an energy and press **Trace this energy**
   (ridge detection), or tick *Click* and click on the band, or **Read
   points from a file**. Two or three energies are better than one: one
   contour fixes the band's shape but not its energy scale.
3. **Fit** -- the hoppings with their errors and the rms misfit; the model's
   contour is drawn in red on the plane of step 2. With one energy, hold
   t1; *Refine the azimuth* if the orientation is known only roughly. An
   orange note means the points cannot tell the shells apart: use fewer.
4. **Use the model** -- **Model to list** (its spectral function, as a map
   beside the data), **Save parameters...**, **Copy results**.

## 27. Calculated bands

Tool panel → Calculation → **Calculated bands...** puts a band-structure
calculation beside the measurement. Select what to compare it with — k-maps,
kz maps in momentum and cuts in momentum — or nothing, to look at the
calculation alone. Four steps.

Like every step-by-step window (Crystal structure, Tight-binding fit, Sum
repeated scans, Detector corrections), it **asks before closing** — the last
step's Close sits where Next was, so a quick double click no longer loses the
work; the title bar's close and Esc ask too. **Enter in a field takes its
value** (and redraws what depends on it) — in every window of the program it
never presses Next, Close or any other button, whatever the value set going
(even a computation with a progress bar); only a button that has the focus
answers Enter. Where a window redraws only on a button — the Brillouin-zone
tool's **Plot**, a viewer's selection-box **Refresh**, the core-level zoom's
**Show** — Enter in its fields does that redraw, and nothing else. What you set is **kept**: the window opens as you left
it (**Restore defaults** at the bottom left starts afresh).

**1. The calculation.** First say **what you want to see**:

- **E–k dispersion** — the bands along a path, or along your measured cut.
  Any calculation works.
- **Constant-energy contours** — to compare with a k-map (Fermi-surface
  map).
- **3D Fermi surface.**

The last two need the bands **everywhere in the zone** (a band grid or a
Hamiltonian); a band structure along a path cannot give them, and the window
says so rather than going on. The choice also sorts the file list of
**Open calculation...**: its first filter shows only the kinds of file that
can give what you want. Selected k-maps choose contours for you, selected
cuts the dispersion. A file can also be dropped on the window; **Reopen …**
opens last time's file again.

| File | Written by | What it holds | Can show |
|---|---|---|---|
| `.bxsf` (XCrySDen band grid) | Wien2k, Quantum ESPRESSO, VASP, Wannier90 | the bands everywhere in the zone, interpolated | all three |
| `_hr.dat` (+ `.win` beside it) | Wannier90, any tight-binding model | a Hamiltonian, diagonalised exactly at each k; orbital weights | all three |
| `wannier90_hr_fplo.dat` (+ `KPOINTS`, `wannier_name_center.txt`) | FPLO's Wannier Hamiltonian | as `_hr.dat`; weights per Wannier function, per shell ("Cr 3d") and per spin | all three |
| `band.dat` with rows `k1 k2 k3 E1 E2 …` (+ `KPOINTS`) | FPLO + Wannier tools | a 3-D band grid | all three |
| MATLAB `.mat` holding `Ek_3D` (+ `KPOINTS`) | the theory side's MATLAB, from that 3-D `band.dat` | a 3-D band grid | all three |
| `EIGENVAL` (+ `POSCAR`, `PROCAR`) | VASP | bands along a path or on a grid; atomic and orbital weights | a grid: all three |
| `_band.dat` (+ `labelinfo`) | Wannier90 | bands along a path | dispersion |
| `.spaghetti_ene` (+ `.klist_band`) | Wien2k | bands along a path | dispersion |
| `band.dat`, `Eigenvalue.dat` (+ `KPOINTS`, `orbital.dat`) | FPLO + Wannier tools (`k_path.exe`, `H.exe`, `plot.exe`) | bands along a path (a third column is a weight) | dispersion |

Files are recognised by their **content**, not their name: a standard file
renamed on the way (`hamiltonian.txt` that is really an `_hr.dat`) opens all
the same. The FPLO tools' `KPOINTS` (first line "segments
points-per-segment", then the direct and reciprocal vectors, then each
segment's ends with `! label`) gives the cell and the labels — keep it beside
the band file; without one, the cell is taken from the Brillouin-zone tool's
crystal. Long reads (a 3-D `band.dat` of 101³ points and 136 bands is about
2 GB; a large `_hr.dat`) show a **progress bar** with Cancel; a big grid
keeps only the bands near E_F (*Reading options ▸ Big grids*).

**Which file? ...** helps when you do not know which of a dozen files is the
calculation: choose the folder the theory side sent, and every file in it is
read and said in plain words — what it is, what it can show, what reading it
lacks (the cell, the labels), whether it is a standard file under another
name, or one file **cut into pieces** (`band.dat.1, band.dat.2`,
`band_part1.dat …`, `xaa, xab` — as files are split to pass an upload
limit). Green rows open for what you want to see; **Join the pieces, then
open** writes the pieces back into one file beside them.

*Reading options ▸* (usually nothing to change): *Set E_F myself* overrides
the file's Fermi level; *bxsf k units* is for a grid written without 2π; *Big
grids: keep bands from … to …* of E_F.

**The crystal.** The cell (a, b, c, α, β, γ) is read from the calculation
and the **space group** that names the zone's points is guessed from its
shape. Once read, both are **locked** — as are *From the Brillouin-zone
tool* — so a stray click cannot change the calculation's lattice. Tick
**Change the lattice parameters** to change them on purpose: typing a new
cell (or taking the Brillouin-zone tool's) scales the calculation to it —
its bands, held in fractional coordinates, stretch with the zone.

**A calculation in another format.** When the theory side's files are none of
these, the safest route is to ask for (or convert to) one of the standard
ones: a `.bxsf` for bands everywhere, `_hr.dat` + `.win` for a Hamiltonian,
or a two-column `_band.dat` (distance, energy; a blank line between bands)
with its `_band.labelinfo.dat` for a path. What to ask them: which columns,
which order of k, the units of k (with or without 2π, fractional or
Cartesian) and whether energies are already relative to E_F.

**2. What it looks like** — by what you chose in step 1:

- **Along a path** (dispersion). The calculation through the high-symmetry
  points (`Γ-M-K-Γ-A-L-H-A`; `G` for Γ, `|` for a jump), over the
  **measurement along the same path**, assembled from the selected maps: each
  in-plane segment from the map measured on its k_z plane (Γ-M-K-Γ from the
  Γ-plane map, A-L-H-A from the zone-top one), Γ-A from a kz map. Each map's
  k_z comes from its photon energy and the inner potential V0; correct it in
  the table if needed, untick a map to leave it out, and give the azimuth of
  the crystal's a axis on the map. **Highlight bands** draws chosen bands
  thick and in colour, the others grey — their numbers as the file counts
  them (`41-44`), or `EF` for every band crossing E_F (step 1 lists them).
- **Constant-energy contours.** The calculated contours at one energy on any
  plane of k-space, set as in the Brillouin-zone tool: its **normal (hkl)**
  (conventional reciprocal axes; (001) gives k_x-k_y), and **k⊥**, how far it
  sits from Γ along that normal — **Through point** puts it through a
  high-symmetry point (Γ K M at 0, A H L at π/c for (001)); **In-plane
  rotation** turns the two drawn axes about the normal (to match a measured
  map's azimuth, say). The **k⊥
  integration window** sums the bands over a range of k⊥ about the plane —
  along the (hkl) normal, not the absolute k_z — as a map with a finite k⊥
  resolution does; the integrated intensity is drawn grey under the
  contours of the middle plane. **Show the plane in the Brillouin zone...**
  opens a turnable 3-D view of the zone with the plane in it (red), the
  faces of its k⊥ window (light red) and the named points it holds; it
  follows the plane as you change it. The zones' outlines on the plane are
  dashed grey.
- **3D Fermi surface.** The surfaces where the bands cross the energy,
  inside the first zone, with the zone's edges and named points. Drag to
  turn, wheel to zoom; Enter in *Energy* redraws it; *Detail* trades speed
  for smoothness (a Hamiltonian is diagonalised on a grid first, with a
  progress bar). **Open in its own window** shows it large, with *Save
  picture*; the energy, detail and bands are still set here. With this goal
  there is nothing to compare, so the Compare step is left out.

Each band's **colour** (from the list, or *Other colour...*), and for a
surface its **opacity**, and whether it is **shown**, are set in the band
table. Several bands may share one colour — to group the electron pockets
in one colour and the hole pockets in another, say — and a band keeps its
colour from the contours to the 3-D view.

**Colour by orbital weight** draws fat bands along a path, or splits the
contours where the chosen orbital's weight exceeds one half.

![Along the path](images/bands_path.png)

**3. Compare.** One measurement at a time, with the calculation drawn where
it was really taken: a **cut** gets the bands along its direction — on a
fixed k_z plane, or on the curved surface its photon energy reaches (k_z
changes along the cut and with the energy); a **map** gets the calculated
contours at the energy chosen. **Fit to the cut's bands** traces the
measured bands (ridge detection) and finds **Z** — the band-width
renormalisation, measured E = Z × calculated E + shift — and the **shift**,
each point matched to its nearest band, robustly. Z and the shift apply in
every step and every output. On the demo cut it finds Z = 0.798 and a shift
of −49.8 meV, the values the demo was made with (0.8 and −50 meV).

![Compare with a cut](images/bands_compare.png)

**4. Output.** With no measurement selected, step 3 says so: press Next.

- **To figure composer...** makes the panels ticked — ticked for you by
  what you chose in step 1: the measurement along the path with the
  calculation over it, the calculation alone (the points' names on its axis),
  the compared cut or map with its calculation, the calculated contours alone
  at the energies listed (no map needed), the 3D Fermi surface as drawn in
  step 2, and — with **Constant-energy series of the compared map** — the
  measured maps at each energy with the zones dashed (top row) and the
  calculated contours (bottom row).
- **Save … to main panel** (beside it) puts the calculation in the list as
  data the main panel's viewers open, and what it saves follows step 1:
  - *dispersion* — a cut: the bands as a spectrum (Lorentzians of the width
    set) along the compared cut or along the path. A cut along a path keeps
    its high-symmetry points: the cut viewer shows their names along the
    top, with a dashed line at each (**High-symmetry points**, ticked; the k
    values stay on the bottom axis), and the figure composer puts them on
    the k axis.
  - *contours* — **Save calculated map to main panel**: a k-map of the
    calculated intensity on the plane of step 2 (its (hkl), k⊥, in-plane
    rotation and k⊥ window) at every energy of *Energies of the saved map*, to slice and
    compare like a measured map.
  - *3D Fermi surface* — greyed out: the viewers open cuts and maps, not a
    surface; the figure composer (or *Save picture*) takes it.

![Output](images/bands_output.png)

# Part V — Figures and output

How a figure is put together, saved with its data and exported — and, at
the end, two figures built step by step.

---

## 28. The figure composer

Tool panel → Advanced → **Figure composer...** with **nothing selected**
opens a blank page, like a figure of a manuscript: four empty panels (a),
(b), (c), (d) at a journal's column width. With cuts (or curves) selected,
each becomes a panel (so does a frozen slice's **Figure composer...**
button; a stack plot or a 3-D view has its own **To figure composer**
button).

![Figure composer](images/figure_composer.png)

- **Lay out the page.** Click a panel to select it; drag it to move it and
  drag its handles to resize it (Shift on a corner keeps its shape). Panels
  snap to each other's edges and to a 0.5 mm grid (**Snap**); the arrow keys
  nudge by 0.5 mm (Shift: 5 mm). *Arrange as grid* (figure tools, Panels
  tab) puts them back on rows × columns. **Re-letter** names them (a), (b)…
  in reading order.
- **Fill a panel.** Double-click an empty panel, or use **Load from
  list...** / **Load picture...** on the left: a dataset from the main list,
  or a picture from a file (PNG, JPEG, TIFF… — a sketch or a figure made
  elsewhere; it keeps its shape and has no axes).
- **Insets.** **Zoom-in inset**: drag out a region of a panel — a small
  panel magnifying it appears in the corner and the region is outlined.
  **Panel** drawn inside another panel makes an inset to fill with another
  dataset. Insets take no letter and move with their panel.
- **Draw.** **Arrow**, **Line**, **Box**, **Ellipse**, **Text**: drag (or
  click, for text) on the page. A selected drawing can be moved, its ends
  or corners dragged, and its colour, width, line style, fill and text set
  above the figure tools. Double-click a text to edit it; Del removes the
  selection; right-click for more (bring to front, delete…).
- **Axes and fonts.** The figure tools beside the page act on the selected
  panel (or all of them): ranges, ticks, labels, which edges carry numbers,
  the panel's own font size, colour map, colour bar, reference lines.
- **Undo / Redo** (Ctrl+Z / Ctrl+Y) cover every change. **PDF**, **PNG**,
  **SVG** export, **Copy** puts it on the clipboard; empty panels are never
  written. The composer's **Figure tools** are not a viewer's **Plot
  tools**, which draw guides on the live plot.

### Arranging panels

Click a panel, then **Ctrl+click** (or Shift+click) the others to select
several. The **Arrange** buttons on the left then:

| Button | What it does |
|---|---|
| ⇤ Left, Right ⇥, ⤒ Top, Bottom ⤓ | line up their left / right / top / bottom edges (at the outermost one) |
| ↔ Centres, ↕ Centres | put their centres on one vertical / horizontal line (the first-selected panel's) |
| Space ↔, Space ↕ | equal gaps between them, left to right / top to bottom (three or more) |
| = Width, = Height | the first-selected panel's width / height for all |

![Several panels selected, and the Arrange buttons](images/composer_arrange.png)

### Letters and colour bars

- **Letters in line** (*Text* tab, on by default): the letters of panels
  side by side sit at one height and those of panels one above another at
  one left edge, whatever titles or axes each panel has. **Nudge (mm)**
  moves one panel's letter (clear of a long axis title, say);
  **Re-letter** names them (a), (b)… in reading order.
- **Colour bar** (*Colour* tab): **Ends** writes the levels, the words
  *High* / *Low*, or nothing; **Length (%)** of the panel's height or width;
  **Position** *right* (upright beside the panel) or *below* (lying under
  it, centred). A 3-D view's bar is sized to the cube as drawn, not to the
  whole frame round it.

![The colour bar settings](images/composer_colour_bar.png)

### Curves: lines, symbols and error bars

*Overlays* tab, **Curves drawn on the panel**: select curves in the list (or
none, for all), then set their **colour**, **width**, **style** (solid,
dashed, dotted…), **symbols** (circle, square, triangles, diamond, with a
size) and whether they are **joined by a line**; **± σ as error bars**
turns a curve's two dotted ±σ curves into error bars at its points. **Apply
to the curves**. A curve with gaps (NaN) is drawn in pieces.

### Saving a figure to edit later, with its data

**Save** / **Save as...** write the figure in the program's own format,
`fig1.afig`: the page, its drawings and every panel with its data and
settings — open it again with **Open...** (or pick it in **Load data...**
in the main window) and go on editing exactly where you stopped. A 3-D view
is kept as the picture it shows.

Beside it, **`data for fig1/`** gets one folder per panel (`a`, `b`, …):

```
data for fig1/
  a/  processed/  what the panel shows: the image as a .nxs cut (opens in the
                  viewer again), its lines (calculated bands, contours,
                  curves) as lines.csv, a picture as .png
      original/   what it was made from: the first measurement file (a
                  k-converted, symmetrised, cut-out plane traces back to the
                  file it was first read from), the calculation with the
                  files it is read with (.win, POSCAR...), the picture file
                  or the CIF
      sources.txt which entry of the file, and the processing in between
  b/  ...
```

The original file is copied even when the panel shows something computed
from it many steps later: every dataset remembers the file it began as. A
figure opened again whose original files have moved reuses the copies in
its own earlier data folder.

## 29. Exporting a panel and saving data

Right-click any image → **Export** → **Data (.nxs, Igor)...**, **Image...** or
**Axes only (.eps)...**. It always exports exactly what the panel shows
(range and all): the data as a dataset this program reads back, or for
Igor Pro (`.itx` text or `.ibw` binary wave, with the axes scaled, the units
set and the metadata in the wave note); an image
with the colormap applied (PNG, TIFF or JPEG, at a chosen size and dpi, with
or without the axes, with the Plot tools guides if there are any, and with
the Brillouin zone if the dataset has one); or the frame, ticks and labels as
vector PostScript to place the image in, in Illustrator (optionally with the
Brillouin zone drawn in the box — see [§15](#15-brillouin-zone)).

![Export](images/export_dialog.png)

**Save** (tool panel, Datasets box) writes the selected datasets to a
`.nxs` file of your own -- several selected rows go into one file -- or, by
choosing the file type, for Igor Pro: one `.itx` holding them all, or one
`.ibw` per dataset. An Igor file is a copy; the rows stay as they were.

## 30. Worked example: two figures in the figure composer

The two figures below were made from the demo data
([§2](#2-trying-it-with-demo-data)) with the program's windows and put
together in the figure composer. `devtools/figure_tests.py` runs exactly
these steps, so it is also a script to read alongside this section.

### Figure 1: constant-energy contours

![Example figure 1](images/example_figure_1.png)

**(a) The notched cube.** Convert the Γ-plane map to momentum
([§13](#13-converting-a-map-to-momentum)); select the k-map → Advanced →
**Notched cube...**. *Change sampling...* to k_x, k_y within ±1.25 Å⁻¹ and E
up to E_F; put the notch at k_x = k_y = 0 down to the lowest energy, corner
**+ − +**; colormap *blue*, reversed; azimuth −52°, elevation 24°. **To
figure composer...**: a new figure. In the composer, *Colour* tab: **Colour
bar** on, Ends *high-low*, Position *below*, Length 60 %.

**(b) Measured and calculated contours.** Select the Γ-plane k-map →
Calculation → **Calculated bands...**; open `demo_calc_hr.dat`; in step 2
choose the **Weight** `Se-pxy`; in step 3 choose the k-map, *Integrate*
0.050 eV, then a cut to **Fit to the cut's bands** (Z and the shift apply to
the contours too); in step 4 tick **Constant-energy series at** `0, -0.3,
-0.6, -0.9`, untick the others, **To figure composer...** → the figure of
(a). Eight panels arrive: the measured maps with the zones dashed, and the
calculated contours alone, red where the weight is below one half and green
above. In each, set the colormap (*blue*, reversed), the ranges ±1.2 and
drop the titles of the bottom row (*Text* tab); add the (i)…(viii) and Γ, M,
K notes with the **Text** tool.

**(c) Stacked contours.** Select the k-map → Advanced → **Stacked
slices...**; *Layers* along E − E_F at `0, -0.15, -0.3, -0.45, -0.6, -0.75,
-0.9`, *Integrate* 0.04, *Fade low values* 45 %, z aspect 2.2, azimuth −45°,
elevation 14°. **To figure composer...** → the same figure.

**Arranging.** Drag the panels roughly into place (the 2 × 4 contours in
two rows, (a) left, (c) right), then select (a), the top row of (b) and (c)
and press **⤒ Top**; select a row of contours and **Space ↔**, and **=
Height** for each row. The letters line up by themselves; **Nudge** moved
(b)'s 6 mm left, clear of its k_y title.

![Selecting panels and lining them up with the Arrange buttons](images/composer_arrange.png)
 Only (a), (b) and (c) carry a
letter: clear the others' *Letter* box. **Save as...** `example_figure_1.afig`
(its data folder is written beside it), and **PNG** / **PDF** for the paper.

### Figure 2: a band crossing followed through k_z

![Example figure 2](images/example_figure_2.png)

**(a) The crystal structure.** Calculation → **Crystal structure...** →
*Read a CIF file* `demo_tetragonal_crystal.cif`. Step 2: Hf teal, Si dark
blue, S yellow, gloss 85 % each, atom size 0.52; *Grey, thin bonds*;
*Lattice constants on the cell edges* (no key of the elements). Step 3:
*Draw a plane* (001) at the middle of the cell, between the S layers, in
grey. Turn it to azimuth −68°, elevation 14°. Step 5: **To figure
composer...** → a new figure.

**(b) The zone with the measured maps.** Select the Γ-plane and A-plane
k-maps, then Calculation → **Crystal structure...** → *Read a CIF file*
`demo_crystal.cif`. Step 4: *Path* `K-H` (the line where the two bands
cross, in red); *Plane through Γ* → map *Γ-plane map* at −0.65 eV; *Plane
at the zone top* → map *A-plane map* at −0.55 eV; *Label on the left* `L,
M`. Step 5: *Picture of* **The Brillouin zone** → **To figure composer...**
→ the figure of (a).

**(c) The cuts stacked along k_z.** Convert the photon-energy scan to
momentum ([§14](#14-photon-energy-scans-kz)), select it → Advanced →
**Stacked slices...**: *Layers* along k_z at the four k_z of the cuts
(8.0, 8.3, 8.6 and 9.0 π/c), *Fade* 30 %, **Spread** 3, azimuth 40°,
elevation 12°. In **Mark a curve through the layers**, add the crossing on
each cut (k_z, k_∥, E) — from (d) — in red. **To figure composer...**.

**(d) The cuts, their MDCs and the crossing against k_z.**

- (i)–(iv): in the k_z map's viewer, **Functions → Slice → Save slice to the
  main list** at each k_z; open each cut in the figure composer (**Add
  panel...**). The dashed guides along the two crossing bands are the
  calculation (Calculated bands, step 3, on the plane k_z = …), shifted by
  the fitted Z and shift: in step 4, **To figure composer...** with only
  *The compared cut* ticked gives the cut with them; set their style to
  *dashed* in the *Overlays* tab. Draw the dashed box with the **Box** tool.
- (v): select the first cut → Advanced → **Stack plot (EDC / MDC)...** →
  MDCs, 40 curves of 3 rows each, over the box's energies; **Send to figure
  composer** → this figure. The guides through the MDCs are the same
  calculated bands, at each MDC's energy.
- (vi): the crossing energy at each k_z — from **MDC fit** on each cut and
  **Dispersion...**, where *Crosses band … at k = …, E = … ± …* gives it
  with its error — as points with error bars (*Overlays*: symbols *circle*,
  not joined, **± σ as error bars**), and the calculated crossing through
  k_z as the red line.

**Arranging.** (a), (b), (c) in the first row: select the three, **⤒ Top**
and **= Height**; the six panels of (d) in the second row, **⤒ Top** and
**= Height**. Letters: (a), (b), (c), (d) only. **Save as...**
`example_figure_2.afig`.

Each panel of either figure can be changed afterwards — open the `.afig`,
select the panel, and use the figure tools; the data folder keeps both the
data as shown and the measurement it came from.

# Part VI — Housekeeping

---

## 31. Autosave, memory and recovery

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
- **Memory box** (bottom of the tool panel). This program's memory, the system's, how many viewers and
  HDF5 files are open, and how much computed data is held. **Free RAM**
  drops cached copies of computed data (they are re-read from their autosaved
  file when needed); **Close all** closes every window and then frees
  memory.
- **Autosaves** (in the Memory box) shows the session folders
  in the system's file browser, to delete old working copies by hand. The
  folder of the session now running (named in the status bar) is in use.
- **Long operations** (reading a big map, converting, de-gridding, assembling
  a folder) run off the GUI thread with a progress bar and a working Cancel.
- **Closing** with unsaved computed data asks first, and names the datasets:
  **Save to a file...**, **Close anyway** (the autosaved copies are deleted)
  or **Cancel** (nothing happens).

## 32. Tips and troubleshooting

- **A greyed tool button** — hover over it: the tooltip says why (e.g. the
  MDC fit needs a cut already in Å⁻¹, cut arithmetic two cuts). Select the
  right rows in the list and it lights up.
- **A tool is not in a viewer's Functions menu** — it is in the tool panel
  (see *Where things are now* in [§3](#3-the-main-window)); check all three tabs.
- **A window is taller than the screen** — it opens no taller than the
  screen, and its contents scroll; the Apply / Close buttons stay at the
  bottom. Every window can be resized.
- **A k conversion is not offered** — the first axis of the map is not an
  angle (a photon energy, temperature…). Check *First axis of a map is* in
  the loader.
- **Cut k conversion of an E − E_F cut** — the energy offset is filled in for
  you from the Fermi fit; if you shifted the axis some other way, type the
  kinetic energy of E_F (photon energy − work function) into *Energy offset*.
- **A file is misread** — choose the reader explicitly in the loader, and see
  *Assumptions worth knowing about* in the
  [reading-data reference](../arpes_viewer/loader/README.md#assumptions-worth-knowing-about)
  (axis order in SOLEIL files).
- **Dragging the cursor on a big map is sluggish** — a map read lazily from
  its file re-slices in about 0.1 s per window, so the other windows wait for
  the mouse to pause. Once the map has been read into memory it follows live.
- **Typing numbers** — every number box takes a decimal point (`0.05`),
  whatever the system language; a comma is not accepted.
- **The machine is paging** — watch the Memory bar (amber past a third of the
  RAM, red past a half); use **Free RAM** or **Close all**.
