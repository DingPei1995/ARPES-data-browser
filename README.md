# ARPES viewer

**Version 1.2.0** (last modified 28 October 2026) — see
[CHANGELOG.md](CHANGELOG.md) for what changed. The version is in
`version.py` and is shown in the main window's title bar and at its foot,
with a link to this repository for bug reports and the latest version
(<https://github.com/DingPei1995/ARPES-data-browser>). Please quote it when
reporting a bug.

> **Terms of use.** This program is for **non-commercial use only**
> (academic research and teaching); any commercial use is not permitted.
> New and updated loaders cover data from **public large-scale facilities**
> (synchrotrons and similar user facilities) only -- formats of individual
> laboratories' own set-ups are not supported. For collaborations, or
> questions about the details, contact
> [ding.pei@synchrotron-soleil.fr](mailto:ding.pei@synchrotron-soleil.fr).

A Python/PyQt5 program for ARPES data: browse the measurements in a folder,
open them in their own windows, convert to momentum, cut, fit, process, and
lay the results out as a publication figure.

It began as a rewrite of the lab's `h5Loader*.py` for SOLEIL ANTARES `*.nxs`
files, and reading those is now one loader among several rather than the
whole program. It currently reads **SOLEIL ANTARES** (NeXus/HDF5), **SOLEIL
CASSIOPEE** (Scienta SES text, single spectra and whole folders), **ALS
MAESTRO** (BL 7.0.2 HDF5), **ALBA LOREA** (BL20 NeXus), **Diamond I05**
(HR-ARPES) and **I05-1** (nano-ARPES) NeXus, **MAX IV Bloch**, **Elettra
APE**, **BESSY One Cube** and **One Square** and **SOLARIS URANOS** (Scienta
SES `.ibw`, `.zip` and `.txt`), **ASTRID2 SGM4** (HDF5), and its own saved
format. Adding a beamline is a new file in `loader/beamlines/`.

The round-by-round history, and the reasoning behind the choices below, is in
[CHANGELOG.md](CHANGELOG.md).

**User guide with screenshots:** [English](docs/USER_GUIDE_en.md) ·
[中文](docs/USER_GUIDE_zh.md) · [Français](docs/USER_GUIDE_fr.md), and as
PDF: [English](docs/pdf/ARPES_viewer_user_guide_en.pdf) ·
[中文](docs/pdf/ARPES_viewer_user_guide_zh.pdf) ·
[Français](docs/pdf/ARPES_viewer_user_guide_fr.pdf). To try the program
without a beamtime, `python devtools/demo_data.py` writes synthetic data to
`demo_data/`.

```
python ARPES_viewer.py
```

Needs Python 3.9+, PyQt5, pyqtgraph, numpy, h5py, scipy.

**On a Mac, double-click `ARPES viewer.app`** (in this folder) instead. It
finds a Python that has the packages -- the one `python3` is in Terminal,
conda environments included, or the usual Anaconda / Homebrew / python.org
places -- remembers it (`~/.arpes_viewer/python_path`) and opens the main
panel. If none has them, it says what to install (the command is copied to
the clipboard), or lets you choose a Python. Drag it to the Dock to keep it
at hand; a copy moved elsewhere asks once where the program's folder is.
The first time, if macOS says the app is from an unidentified developer
(the repository was downloaded as a ZIP rather than cloned): **System
Settings → Privacy & Security → Open Anyway** (on older macOS: right-click
the app → Open). If the program stops with an error, a dialog shows the end
of `~/.arpes_viewer/launcher.log`.

---

## Contents

- [Getting data in](#getting-data-in)
- [What the main panel is](#what-the-main-panel-is)
- [Data kinds](#data-kinds)
- [The viewer windows](#the-viewer-windows)
- [Analysis](#analysis)
- [Figures and export](#figures-and-export)
- [The experiment log](#the-experiment-log)
- [Nothing is lost: autosave, the log, and recovery](#nothing-is-lost-autosave-the-log-and-recovery)
- [Long operations](#long-operations)
- [Adding a beamline](#adding-a-beamline)
- [Layout of the code](#layout-of-the-code)
- [Assumptions worth knowing about](#assumptions-worth-knowing-about)
- [Tests](#tests)

---

## Getting data in

**Load data...** opens the Loader window, which is where everything about
reading a file is decided before it reaches the list:

- **Files** -- pick one or many. The file-type filter is built from what the
  loaders say they read, so a beamline whose files are not `.nxs` shows up.
- **Reader** -- detection is automatic and says what it found
  (`Detected: 3 × SOLEIL CASSIOPEE`). Override it when a file from an
  unfamiliar layout is being misread. Files this program saved are recognised
  outright and never need a beamline chosen.
- **Datasets found** -- what is inside, before anything is read. A `.nxs` file
  routinely holds several complete measurements; each gets its own row.
- **Axis order** -- reorder the array's dimensions, once, at the door. The
  axis vectors and their labels come along.
- **First axis of a map is** -- a map is angle vs angle vs energy by default,
  but the same acquisition records a photon-energy, temperature or
  gate-voltage series. Saying so here keeps the wrong unit out of every later
  step and stops a k conversion being offered for an axis that is not an
  angle. Leave it at **As the file says** to keep what the reader worked out.

Only *structure* is read here -- entry names, kinds, shapes -- which is cheap
even for a multi-GB scan. A file that is listed and never opened is never
read.

### CASSIOPEE folders

A numbered folder from CASSIOPEE (`<name>_1_ROI1_.txt`, `<name>_2_ROI1_.txt`,
...) is **one measurement**, not a pile of spectra. Select any file in it --
or several, or all of them -- and the folder is listed once, as the assembled
cube. Its row is named after the folder's base name with the photon energy
and polarisation from the `_i` files (`FineMap 117eV LH`; a photon-energy
scan shows its range, `KzScan 40-120eV LH`), so repeats of one scan at
different settings can be told apart.

What the folder becomes depends on what was stepped between the spectra,
which is read from the small `_i` parameter files that sit beside them:

| what varies across the folder | first axis | kind |
|---|---|---|
| the sample's polar angle | theta (deg) | `map` |
| the monochromator | photon energy (eV) | `kz_map` |
| neither | cut index | `map` |

The stepped axis is straightened to an even grid; the motor readbacks are
kept in `info` as `cassiopee.axis0_measured`.

A **kz map** stacks the members on the first spectrum's energy axis, referred
to the Fermi level (`E - E_F = E_kin - hv + phi`) exactly as the lab's MATLAB
loader does. Nothing is interpolated or trimmed. The per-member photon
energy, work function and energy-axis start are recorded in `info`, which is
what a proper per-photon-energy Fermi-edge calibration will need later.
`load_series(..., energy_reference=...)` offers `"common"` (the default,
above), `"per_member"` (place each member by its own hv and keep only the
range they share) and `"kinetic"` (no referencing at all).

### ANTARES photon-energy scans (PhD)

An ANTARES photon-energy scan is saved as one top-level entry **per photon
energy** -- `Phd_0002`, `Phd_0003`, ... -- each holding a single MBS frame
(`MBS_0_<n>/Image32_1`) with its own energy and angle scales and its `hv`.
The browser lists the whole file as **one `kz map` row** (entry `Phd_*`);
the first entry, which holds only the beamline state, is skipped.

- Frames are sorted by `hv`; the photon-energy axis is the readback,
  straightened when it is plainly an even scan that jitters (the readbacks
  stay in `antares.axis0_measured`).
- **Energy is referred to E_F per frame**: `E - E_F = E_kin - hv + W`, with
  `W = 4.5 eV` assumed (the files record none) and written to
  `work_function_eV`, which **kz -> momentum** reads back. If the referred
  windows coincide, frames are stacked untouched; if not, they are resampled
  linearly along energy onto the range all of them cover. Run **kz map
  processing** to fit the real Fermi level per photon energy.
- A single frame can still be opened on its own, by its entry name, as a
  cut.

### ALS MAESTRO (BL 7.0.2)

MAESTRO's `.h5` files are read by **ALS MAESTRO** (`loader/beamlines/maestro.py`),
ported from the MAESTRO reader in
[TensorSpec](https://github.com/fawkesdx/TensorSpec) (MIT). A file is
recognised by its three top-level groups (`0D_Data`, `2D_Data`, `Headers`).
The scan plan in `Headers/Low_Level_Scan` says which motors were stepped,
and over what range:

| scan | opens as |
|---|---|
| one spectrum (`Process_000`) | `cut` |
| deflector scan | `map` (deflector × slit angle × energy), convertible to k |
| XY mesh | `spem_4d` |
| deflector + X in one loop (two-motor line) | one `map` per X position |
| focus XY (XY around a deflector loop) | one `spem_4d` per deflector step |
| swept `Process_000` with motors in `0D_Data` | a map or `spem_4d` if the motors form a regular grid, else a map over the scan points |

- The spectra stay in the file and are read as they are looked at, as for
  an ANTARES spatial scan.
- A fixed-mode slit axis given in pixels becomes degrees at 0.048 °/px
  (R4000 Angular30), centred on the detector; `info` says so.
- An aborted scan is read as far as it went: an XY mesh keeps its complete
  rows, anything else opens as a list of the points recorded
  (`info["maestro.aborted"]`).
- The energy axis is the analyser's own scale (kinetic energy). Set E_F
  from the data, as for CASSIOPEE spin.

### ALBA LOREA, Diamond I05 and I05-1 (NeXus)

Three NeXus readers share one module, `loader/formats/nexus_arpes.py`, and differ
only in where the spectra are and how the scanned axes are named. Each was
written against `h5info` dumps of a cut, a deflector map, a photon-energy
scan and a spatial scan from that beamline.

| loader | file | spectra | scanned axes come from |
|---|---|---|---|
| **ALBA LOREA** (`loader/beamlines/alba.py`) | `NNN_<sample>.nxs` | `entry1/data/data` | the group's `@axes` list |
| **Diamond I05 HR-ARPES** (`loader/beamlines/diamond.py`) | `i05-NNNNNN.nxs` | `entry1/analyser/data` | `@axis` on each vector |
| **Diamond I05-1 nano-ARPES** (`loader/beamlines/diamond.py`) | `i05-1-NNNNN.nxs` | `entry1/analyser/data` | `@axis` on each vector |

All three are recognised by structure. ALBA has an `NXdata` group with
`@signal`/`@axes` and no `entry1/analyser`. Diamond has `entry1/analyser`.
The two Diamond branches are told apart by `instrument/name`, then the file
name, then the I05-1 zone plate and `sm*` motors.

| scanned | ALBA | I05 | I05-1 | opens as |
|---|---|---|---|---|
| nothing | `angles, energies` | `(1, angles, energies)` | same | `cut` |
| deflector | `defl_angles` | `deflector_x` | `deflector_x` | `map`, convertible to k |
| photon energy | `photon_energies` (+ `binding_energies`) | `energy` | `value` (`scan_group`) | `kz map` |
| sample position | `saz, sax` | `sax, saz` on the grid | `smx, smy` on the grid | `spem_4d` |

- Any other scanned motor is read the same way. A polar or tilt angle opens
  as a convertible map and anything else as a map of its own name. A
  two-axis scan that is not a spatial pair opens as one map per outer step.
- **Spatial scans at Diamond** store both motors on the whole
  `(slow, fast)` grid (`@axis = '1,2'`). Which motor steps along which
  dimension is read off the values, not assumed.
- **Energy** of a cut, map or spatial scan is the analyser's kinetic
  energy, as for MAESTRO. The monochromator energy is in `info` as
  `photon_energy_eV`, and the sample temperature as `SampleTemperature_K`.
- **kz maps** are read into memory and put on one `E - E_F` axis. At ALBA
  the file's own `binding_energies` is used. It is negated if it is a
  positive binding energy, and the work function is worked out from it and
  the per-hv kinetic axes. At Diamond each spectrum's own kinetic window is
  referred by its own hv (`E_kin - hv + W`, `W = 4.5 eV` assumed, as for
  ANTARES PhD). The spectra are stacked untouched when the referred windows
  coincide, else resampled onto their common range. The hv axis is
  straightened when it is plainly a jittered even scan. The readbacks are
  kept in `info["axis0_measured"]`.
- Everything except a kz map stays in the file and is read as it is looked
  at. Integer counts are read as float32.

### MAX IV Bloch (Scienta SES `.ibw` and `.zip`)

`loader/formats/ses.py` reads what the Scienta SES software saves, and
`loader/beamlines/maxiv.py` puts the Bloch A branch's (DA30-L) name on it.
It was written against the headers of a cut, a
photon-energy scan, a spatial scan and a deflector-map zip
(`for loader/MaxIV`), and checked against the SES reader in
[peaks](https://github.com/phrgab/peaks). Nothing in it is specific to Bloch
except the manipulator axis names (`A P T X Y Z`), so other SES files of the
same kind open the same way.

| file | SES dimensions | opens as |
|---|---|---|
| `.ibw`, 2-D | energy, `Y-Scale [deg]` | `cut` |
| `.ibw`, 3-D | + `Photon Energy [eV]` (CIS) | `kz map` |
| `.ibw`, 3-D | + `Point` (one manipulator axis) | `map` of that motor |
| `.ibw`, 4-D | + `Point`, `Position` (manipulator scan) | `spem_4d` |
| `.zip` | `Energy`, `Thetax`, `Thetay` per region | `map`, convertible to k |

- **`.ibw`** is a version-5 Igor wave, little- or big-endian, read with a
  small reader of our own (no `igor2` needed). Igor stores it
  column-major, so the file is already `(*scanned, angle, energy)` and is
  memory-mapped, not read. Each axis is the wave's own `start + i * delta`.
  The wave note is the SES header. `[SES]` goes into `info` as `ses.*`,
  `[Manipulator]` as `manipulator.*`, and the photon energy as
  `photon_energy_eV`.
- **Manipulator scans.** `Point` and `Position` are named from the
  run-mode table, which lists every point's `Position, Point, A, P, T, X, Y,
  Z`. The motor whose column repeats `Point` is the one the layers step. At
  Bloch, `Point` is X (fast) and `Position` is Y (slow).
- **Binding energy.** SES's `Energy Scale=Binding` axis (positive below
  E_F) is negated into `E - E_F`, with no work function added. A CIS on a
  kinetic scale is referred by each hv (`E_kin - hv + W`).
- **`.zip`** is SES's map format. `viewer.ini` lists the regions, and each
  region has `Spectrum_<name>.bin` (float32, energy fastest, then the slit
  angle `Thetax`, then the deflector `Thetay`), `Spectrum_<name>.ini` (the
  axes) and `<name>.ini` (the SES header). A zip with several regions lists
  one dataset per region. A `.bin` stored uncompressed is read in place,
  and a compressed one is read into memory.
- The energy of a cut, map or spatial scan is the analyser's kinetic
  energy. Set E_F from the data, as for MAESTRO.

### Elettra APE and ASTRID2 SGM4 (ported from peaks, untested)

These two were ported from the readers in
[peaks](https://github.com/phrgab/peaks) (`loaders/ape.py`, `loaders/sgm4.py`)
without files from either beamline. The tests build files the way those
readers expect them. Open a real file and check it before relying on it.

- **Elettra APE** (`loader/beamlines/ape.py`) writes the same SES files as
  Bloch. `.ibw`, `.zip` and `.txt` all go through the shared SES reader
  (`loader/formats/ses.py`); a `.txt` is parsed by the CASSIOPEE Scienta
  text parser. The run-mode table
  knows APE's motor names (`polar tilt X Y Z`). Because the format is the
  same, the facility is told by the `Location` in the SES header: a file
  whose `Location` says APE or Elettra opens as APE, and anything else stays
  with Bloch (or CASSIOPEE for `.txt`). Choose the reader in the Loader
  window when `Location` was left blank. Spin data is not read, as in peaks.
- **ASTRID2 SGM4** (`loader/beamlines/sgm4.py`): `Entry/Data/ScanDetails`
  gives each slow and fast axis as start, step and length. The spectra in
  `Entry/Data/TransformedData` are one row per point, with the first slow
  axis fastest. They are read lazily through a grid view, and a scan
  stopped early is padded with NaN (`info["sgm4.incomplete"]`).
  - `ShiftX`/`ShiftY` (or `SALX`/`SALY`) are the deflectors, and open as a
    map that converts to k.
  - `FSamX/Y/Z` are the stage in µm; two of them open as a spatial scan.
  - `LEG`/`MEG`/`HEG` are the photon energy, which opens as a kz map with
    `W = 4.416 eV`. Each hv keeps one binding window, as in peaks.
  - `CapZ` (focus) and `Loop` (repeats) are listed one dataset per step.
  - A knife-edge scan holds a summed intensity per point, not spectra. It
    is listed but does not open.

### BESSY One Cube / One Square and SOLARIS URANOS (Scienta SES)

All three write SES files and go through `loader/formats/ses.py`; the
beamline modules (`loader/beamlines/bessy.py`, `loader/beamlines/solaris.py`)
only say which files are theirs. Written against the headers in
`for loader/BESSY` and `for loader/SOLARIS`.

| end station | file | opens as |
|---|---|---|
| One Cube (1^3) | `.ibw` cut | `cut` |
| | `.ibw` photon-energy scan (binding energy) | `kz map` |
| | `.zip` deflector map (`Thetay`) | `map`, convertible to k |
| One Square (1^2) | `.ibw` cut | `cut` |
| | `.ibw` map, `Point [deg]` stepping `P-Axis` | `map` over the polar angle, convertible to k |
| | `.ibw` photon-energy scan (binding energy) | `kz map` |
| | `.ibw` manipulator scan (`Position` = Z, `Point` = polar) | `spem_4d` |
| URANOS | `.ibw` cut | `cut` |
| | `.ibw` manipulator scan (`Point` = `X:`, `Position` = `Y:`) | `spem_4d`, on the motor readbacks |
| | `.zip` deflector map (`Thetay`) | `map`, convertible to k |
| | `.txt` photon-energy scan, one `[Region n]` per hv | `kz map` |

- **Which beamline.** One Square writes `Location=ARPES 1^2`, URANOS
  `Location=URANOS` (`Solaris` in the text export). One Cube leaves SES's
  default `Location=Scienta`, so it is recognised by its analyser too
  (`Instrument=DA30L_8000-11623`). A file that says neither stays with the
  Bloch reader; choose the right one in the Loader window.
- **Manipulator scans.** The motor a dimension steps is found in the
  run-mode table: the column that repeats `Point`/`Position` (Bloch), the
  one that repeats the wave's own axis (One Square, whose table starts each
  row with a running number), or the one that moves along that dimension
  only (URANOS, which numbers its points `1, 2, 3 ...` -- the `X:`/`Y:`
  readbacks then become the axes). Two manipulator dimensions always make a
  spatial scan, even when one of them is an angle.
- **URANOS photon-energy scans** are one text file with a region per hv.
  Each region is referred to E_F by its own hv (`E_kin - hv + W`, W
  assumed 4.5 eV), and resampled onto the range they share when the
  referred windows differ. Run **kz map processing** to fit the real Fermi
  level. A text file with one region, or with all regions at one hv, opens
  as cuts.

### CASSIOPEE spin end station (MBS)

The spin-resolved end station's MBS A-1 analyser writes `.krx` (binary)
and `.txt` (its text export). Both are read by **SOLEIL CASSIOPEE spin
(MBS)**, and the header -- not the file name -- decides what a file is:

| in the file | opens as |
|---|---|
| one image, main detector | `cut` (Y angle × kinetic energy) |
| one image per deflector step (`MapNoXSteps` > 1) | `map` (deflector × Y angle × kinetic energy) |
| spin system on, main detector off: one spectrum per `SpinComp#n` | `spin_edc` |

- **Energy stays kinetic.** These files record no photon energy, sample
  angles, temperature or work function. Find E_F from the data: the cut
  viewer's **Fermi level**, or the curve fit's Fermi edge (which lists a copy
  with E_F = 0).
- The `.krx` and `.txt` of the same cut read identically, to the count.
- A map with fewer images than planned (interrupted) is read as far as it
  goes and says so in `info["mbs.map_note"]`; a two-direction (X and Y)
  deflector map is refused, since it is 4-D.
- Every header field is kept in `info` as `mbs.<field>`; the spin channels'
  labels as `spin.component.<n>`.

---

## What the main panel is

A browser, deliberately. One row per *dataset*, not per file. Click a row to
read its metadata, **double-click** to open it in its own window.

Right-click a row (or a selection) for the three ways of opening it:
**Open**, **Open slit cut** and **Open scan-axis cut** (a map's cut in a
window of its own, without its contour on screen). Everything else done to
the rows -- **Experiment info**, **Log format...** / **Auto Exp. log**,
**Rename** | **Save** | **Remove** (one row), **Sum repeated scans...**,
**Session log...** -- is the
*Datasets* box of the tool panel, between the tools and the Memory box;
combining datasets is the **Combine...** button under the list.

Under the list, in a box titled *With the selected datasets*, three buttons
named for what they do (they were *Combine data*, *Data operations* and
*Process*, names that did not tell the three apart):

- **Combine datasets...** -- stack cuts into a map along an axis of your
  own, or join the pieces of a measurement taken in segments (a map's
  deflector range, a kz map's photon energies, a cut's energy window) end
  to end. See *Combining datasets* below.
- **Crop / bin / normalise...** -- change the grid or the scale: crop to a
  range, bin points together, normalise each spectrum to itself.
- **Smooth / derivative...** -- the filters that bring out bands: smoothing,
  derivatives, curvature, backgrounds, symmetrisation.

At the foot of the panel: the version and its date, and the link to GitHub
for bug reports and updates.

The menu only offers what fits the selection. An entry is greyed out when
the selected dataset is the wrong kind (a stack plot needs a cut, the 3-D
view a map) or the number of rows is wrong (arithmetic takes exactly two);
with several rows selected, one row that does not fit is enough. Hovering
over an entry shows what it does -- or, for a greyed one, which row does
not fit and why -- at once, as the pointer moves.

| entry | rows | kinds |
|---|---|---|
| Open, Experiment info, Rename | one | any |
| Log format | none needed | -- |
| Auto Exp. log | one or more | any |
| Open slit cut / deflector cut | one or more | maps (angle, kz, k) |
| Open in figure composer | one or more | cuts, EDC / MDC / spin EDC |
| Combine datasets | two or more | cuts, or maps of one kind |
| Stack plot, MDC / EDC fit, Core-level check | one | cut |
| 3D view | one | map |
| Cut arithmetic on the two | two | cuts |
| Save, Remove | one or more | any readable / any |
| Sum repeated scans | two or more | cuts, or maps of one kind |
| Ridge detection | one or more | cuts, or maps |
| Detector corrections | one | cut, map or kz map as measured |
| Tight-binding fit | none or one | a map in momentum (or points from a file) |
| Crystal structure | none needed | -- |

A cut opened from the list keeps its map's contour off screen; the cut's
**Functions → Show the map** brings it up, double-clicking the map's row
does the same, and an unseen contour closes itself when its last cut window
is closed.

**Clear list** asks first. If some datasets were computed in this session
and not saved to a file of your own, it lists them and offers **Save to a
file...** (then clears), **Clear without saving** (their auto-saved copies
go too) or **Cancel**; otherwise it is a plain yes/no. Files on disk are
never deleted, and viewers already open keep their data.

One window per thing you are looking at means several files -- or a contour
and both its cuts -- can be on screen at once and arranged freely. The
colormap chosen in the main panel applies to every open window; metadata
stays in the main panel because it describes the file rather than any one
view, and each viewer exports what it itself shows.

---

## Data kinds

| kind | shape | what it is |
|---|---|---|
| `cut` | (angle, E) | one spectrum at one deflector position |
| `map` | (angle, angle/k, E) | a deflector or polar-angle scan |
| `k_map` | (kx, ky, E) | a map converted to momentum |
| `kz_map` | (hv, angle/k, E) | a photon-energy scan |
| `kz_map_k` | (k_z, k_par, E) | that scan converted to momentum |
| `spem_1d` | (x, k, E) | a real-space line scan |
| `spem_4d` | (x, y, k, E) | a real-space raster scan |
| `edc` | (E, channel) | one or more curves against energy |
| `mdc` | (angle/k, channel) | one or more curves against angle or momentum |
| `spin_edc` | (E, spin channel) | a spin-resolved EDC: one curve per spin channel |

`map`, `k_map`, `kz_map` and `kz_map_k` are the three-axis cubes
(`CUBE_KINDS`): they share a viewer, the 3-D view and the processing panels,
and differ only in what their first two axes mean. `k_map` and `kz_map_k` are
the ones already in momentum (`MOMENTUM_KINDS`), so nothing that converts
angles is offered for them.

The three curve kinds (`CURVE_KINDS`) are tables: the physical axis, and
one or more **channels**. The channel names travel in
`info["curve.channels"]`; a channel called `σ <name>` is the standard
deviation of `<name>`, so a polarisation and its error bar, or a spin-up
spectrum and its counting error, stay together through saving and reloading.
An `mdc` is the one kind with no energy axis.

A `kz_map` gets its own conversion rather than the in-plane one: turning a
photon energy into k_z needs the inner potential, which is a property of the
sample and not of the measurement.

Which axes a kind has lives in **one** table, `loader/nxs_file.AXIS_SLOTS`.
Everything else -- the save format, the memory data wrapper, the processing
panels, which axis is the energy -- derives from it.

---

## The viewer windows

- **A cut** opens as the E-vs-k spectrum, with EDC/MDC curves and a
  draggable slice.
- **A cube** opens on its constant-energy contour. **Deflector cut** and
  **Slit cut** open the orthogonal cuts, each in a further window. They belong to the map:
  closing the map closes them too (and any dialog a viewer opened closes
  with that viewer). Snapshots popped out of a viewer are copies and stay
  open.
- **A spatial scan** opens the real-space map together with the E-vs-k
  spectrum at the cursor.
- **A curve** (EDC, MDC, spin EDC) opens in the **curve viewer** -- see
  below.

Every panel with an EDC/MDC readout has **EDC → list** and **MDC → list**
under it: the curve on screen, summed over its ± window, becomes a dataset
in the main list, named with where it was taken (and, from a map's cut
window, the slice). Its kind follows its own axis: a profile along an angle
is listed as an MDC even if it was the "EDC" of a constant-energy contour.

### The top row, and the Functions menu

Every viewer has **one row** along its top: the colormap (and Flip) and the
view controls (axis ranges, Auto, Inv, Reset, Grid) side by side, a
**Functions** menu, and **Plot tools**. A map's window also keeps its **Scan-axis cut** and
**Slit cut** buttons there, since the cuts are how a map is navigated, and,
for a dataset with a Brillouin zone, its **Brillouin zone** tick box.

Since 1.1.0 the tools that take a dataset -- Fermi edge fit, the k and kz
conversions, Align Fermi levels (kz), Remove detector grid, Brillouin zone,
the filters, MDC / EDC fit, cut arithmetic, Slice series figure, the 3-D
views, the figure composer -- are buttons of the launcher's tool panel
(Basic, Advanced and Calculation tabs), not Functions entries. What is left in
**Functions** works on the window itself:

| Viewer | Functions |
|---|---|
| map / kz map / k-map | Arbitrary cut; Slice: Save slice to the main list, Freeze slice in a new window |
| cut | Straighten bands (FS correction); Slice: same |
| slit cut of a map | Straighten bands (FS correction) of the whole cube; Show the map |
| curve (EDC, MDC, spin EDC) | Curve fit, Spin analysis (spin EDCs); Crop, Bin, Normalise, Subtract a background, Shift the axis, Add counting errors (√N) |

Several entries were renamed over time to say what they do -- *Fermi level*
→ **Fermi edge fit**, *FS correction* → **Straighten bands (FS
correction)**, *De-grid* → **Remove detector grid**, *Slice figure* /
*Multi-slice figure* → **Slice series figure**, *Open slice in a new panel* →
**Freeze slice in a new window**, *Deflector cut* → **Scan-axis cut**. The
sections below keep the old names as headings where the history is easier
to follow with them.

**Plot tools** (every viewer and every frozen slice, on the top row) draws
on the live plot: vertical, horizontal or angled guide lines and text
labels, each in its own colour, width and style, draggable on the plot and
drawn into the panel's image export too (`ui/viewers/plot_tools.py`). It is not the
figure composer's **Figure tools...**, which dresses up a figure for a paper;
the two used to share the name *Plot tools*. A frozen slice's
**Figure composer...** button is the way from one to the other.

In the curve viewer the top row keeps the display controls (error bars,
waterfall offset, the Range); an operation chosen from Functions opens the
**Operations** strip under the plot on that operation, to set it up and
**Apply** (**Close** puts the strip away).

An entry that does not apply to the data is not listed (Brillouin zone on
an angle map); one that applies but cannot run yet is greyed, and its
tooltip says why (MDC / EDC fit on a cut still in degrees). New tools go in
this menu too, unless they need to be always in sight.

Every image panel has, under its own title, the level bar, gamma and
interpolation.

### The readout cursor

Right-click an image for **Readout cursor**. It reads the data point under
it; on a cut it also draws the EDC and MDC through that point, summed over
their **±** windows. Under the image a **Cursor** row shows its position:
type a value and press **Enter** to move it there. A value outside the data
leaves the cursor where it is and says so, in red.

**On a map the cursor is one point in the cube, shared by the contour and
its cut windows.** Switching it on (or off) in any of them does the same in
the others. Moving it anywhere -- dragging, typing, or moving a window's
own slider -- moves it everywhere and re-slices the other windows through
the new point: the slit cut is taken at the cursor's deflector angle, the
deflector cut at its slit angle, and the contour at its energy.

The three axes keep one colour each in every window: **deflector red, slit
green, energy blue**. A cursor line is the colour of the axis it holds
constant, and the EDC or MDC it feeds is drawn in the same colour.

**Integration widths are each window's own and are never passed between
windows.** A cut window's EDC ± and MDC ± set only that window's EDC and MDC,
and only that window shades them around its cursor (translucent, dashed
edges, in the axis colours). Each cut's "Integrate over ±" and the contour's
energy ± are likewise its own. The contour's cursor has no width. Only the
*point* is shared.

On a large map read lazily from its file, re-slicing takes about 0.1 s per
window, so while a cursor is being dragged the other windows wait for the
mouse to pause; a map held in memory follows it live.

---

## Analysis

- **k conversion** (cube and single cut) -- the free-electron formula, with
  the inputs shown and checkable. Refused, with a reason, for an axis that is
  not an emission angle. **Set rotation from contour** reads the sample
  rotation off a direction picked on the contour: two points mark a line, so
  the order you click them in and which side of the origin they sit on make
  no difference, and what comes back is the smaller of the two turns that
  stand it vertical. A converted Map is stored in float16 by default (half
  the memory and file size; its 0.05 % steps are far below any point's
  counting noise), with float32 one click away. Values outside float16's
  range are kept in float32 automatically, and every slice is widened to
  float32 before it is summed, cut or fitted.
- **Fermi level** -- fit the edge, offset the energy axis, or divide it out.
- **MDC / EDC fitting** -- fit the peaks across a cut and read the band off
  them: `v_F`, `m*`, and the self-energy. Only accepts a cut whose in-plane
  axis is already in Å⁻¹, because every quantity built on the fitted band
  would otherwise carry meaningless units.
- **Core-level check** -- which element and orbital each peak of a cut's
  angle-integrated EDC is; see below.
- **Arbitrary cuts** -- any straight line through a contour.
- **Fermi-surface correction**, **stack plots**.
- **Cut arithmetic** -- linear and circular dichroism, dividing by a
  reference, differences, ratios and sums between two cuts; see below.
- **Processing** (2-D and over a cube) -- smoothing, derivatives, curvature,
  backgrounds, symmetrisation, despiking.
- **Crop / bin / normalise** -- truncate, self-normalise, compress, on
  several datasets at once when they really are the same format.
- **Combining datasets** -- *Stack* cuts into a map along a new axis (a
  photon-energy series of angle cuts becomes a kz map), or *Join* cuts,
  maps, kz maps or k-maps of one kind end to end along one of their own
  axes (`tools/processing/combine.py`). Pieces on one grid are copied as measured;
  a different step or offset is interpolated onto the first piece's step;
  other axes that differ are resampled onto the range all pieces cover;
  overlaps are averaged (or the earlier / later piece kept); gaps are left
  at zero and reported.
- **Brillouin zone** -- overlay a conventional or irreducible zone, or a
  moiré zone, on a converted contour. Space group first, then the lattice
  parameters it constrains; a 3-D preview shows how the cut plane sits,
  with the named high-symmetry points (Γ, M, K, Z, L, ...) marked, and a
  **Through point** list offers the offsets that put the cut through one of
  them. The points the plane holds are labelled on the contour too. The
  dialog starts from the last crystal used (`~/.arpes_viewer/brillouin_zone.json`).
- **kz map processing** and **kz -> momentum** -- see below.
- **The curve viewer**, its **curve fit** and **spin analysis** -- see below.
- **De-grid** -- remove the detector's grid (MCP pattern or mesh) from a map
  or a cut; see below.
- **Detector corrections**, **sum repeated scans**, **ridge detection**,
  the **tight-binding fit** and the **crystal structure** -- ported from
  the lab's MATLAB programs, each improved where the MATLAB version could
  be; see *Tools ported from the MATLAB programs* below.

### Core-level check

Right-click a cut in the list → **Core-level check...**. Some "cuts" are
really a wide core-level scan, and the first question about one is what
the peaks are. Three steps, each correctable before the next:

1. **Elements.** Type what the sample is made of: symbols (`Nb, S, Cl`) or
   a formula (`NbSCl`, `Bi2Se3`). The match tolerance (0.5 eV) and the
   "uncertain" window (2 eV) are set here too; the last values are kept in
   `~/.arpes_viewer/corelevel.json`.
2. **Energy reference.** The whole cut is summed over angle and the Fermi
   edge looked for where it should be (hν − W on a kinetic axis, zero on
   `E - E_F`). If the window holds it, it is fitted -- the nominal photon
   energy is often off by more than the tolerance -- and hν − E_F is shown
   as an implied work function, which should be 4-5 eV. If not (the usual
   core-level window), E_F is put at hν − W. The axis type, hν, W and E_F
   are all editable, and the E_F line can be dragged.
3. **Results.** The EDC on a binding-energy axis (deeper to the left), each
   peak marked and listed:
   - **sample** -- within the tolerance of a level of the sample's
     elements, with the orbital named (and its spin-orbit partner, when
     that peak is there too);
   - **sample (orbital ambiguous)** -- two levels of the sample fit;
   - **uncertain** -- no level within the tolerance but one within the
     uncertain window: a chemical shift, charging, or not the sample; the
     nearest levels are listed and the call is yours;
   - **not the sample** -- the lines of *common* elements near it (neither
     radioactive nor highly toxic, the usual UHV contaminants first: C, O,
     N, F, Cl, Na, K, S, Si, Ca, Ar, Ta, Mo, W, Cu, Ag, In, Au, ...) with
     where each typically comes from, plus second- and third-order light
     and Auger lines.

   The assignment column is editable (and kept through a re-analysis);
   **Copy report** includes what was changed by hand. **Add EDC to list**
   lists the integrated EDC on the binding-energy scale.

   **Weak lines next to a strong one.** *Log counts* puts the count axis on
   a log scale. *Min height* (default 0.2 % of the strongest peak) is the
   only cut besides the noise threshold; it used to be a fixed 1 %, which
   dropped real weak lines beside a strong one. A peak the search still
   misses -- a shoulder with no maximum of its own -- can be added by hand:
   press *Click to add peak* and click the plot (a double-click always
   adds), or type its binding energy and press *Add*. Its position is
   refined by a local Gaussian fit, it is marked ◆, assigned like any
   other peak, and named with every level within ±2 eV -- the sample's
   first, then all other elements. *Remove selected* drops a peak (a found
   one stays out until *Reset peaks*).

   **Peak table...** lists every tabulated level (Table 1-1, plus the
   typical chemical states), filtered by element (symbols or a formula;
   *Sample* / *All*) and by binding-energy range (*Whole window*, *Plot
   view*, or typed), sortable by any column, the sample's levels in bold.
   Picking a row marks that level on the plots (magenta dashed line);
   *Add as peak* or a double-click adds a peak there -- refined by the
   local fit -- with the level written in as its assignment.

   **Select range** shows a band on the plot to drag over a range;
   **Zoom window...** opens that range (or, with no band, what the plot
   shows) in a window of its own: from/to editable, its own log axis, a
   smoothed curve over the data, and below it the negative second
   derivative (divided by the intensity, *÷ I*, so a weak peak shows as
   clearly as a strong one) -- every peak, a shoulder with no maximum of
   its own included, is a maximum there. A cursor (drag it, or click)
   reads out its energy; *Fit at cursor* fits a Gaussian there and moves
   the cursor to its centre; *Add peak at cursor* (or a double-click) adds
   it to the results. The zoom window follows the results window's peaks
   and the level picked in the peak table.

The table is Table 1-1 of the X-ray Data Booklet (`Table_1-1.pdf` at the
top of the repository), turned once into `tools/spectra/corelevel_table.py` by
`devtools/gen_corelevels.py` -- the PDF is never read at run time. What
the check adds beyond a lookup:

- **Gas-phase values.** The table gives H, N, O, F, Cl and the rare gases
  against the vacuum level of the free gas (O 1s 543.1 eV, where an oxide
  is at ~530). Typical solid-state values (oxides, hydroxides, chlorides,
  sulfides, sulfates, adventitious carbon and its oxidised forms, PTFE,
  SiO2, implanted Ar, the common metal oxides) are compared as well; the
  checkbox turns them off.
- **Spin-orbit doublets** confirm an assignment when both components are
  found at the tabulated splitting.
- **Higher-order light.** A level deeper than hν − W can appear at E_B − hν
  (2hν) or E_B − 2hν (3hν).
- **Auger lines** (C, N, O, F, Na KLL) are named; they sit at a fixed
  kinetic energy, so measuring at a second photon energy tells them apart:
  a core level stays put in binding energy, an Auger line and higher-order
  light move.
- **A rigid shift.** If shifting every peak by the same amount would make
  clearly more of them fit (charging of an insulating sample, or a wrong
  E_F or hν), the shift is reported; **Apply suggested shift** applies it.
- **What is missing.** Levels of the sample inside the window with no peak,
  and elements with no level in the window at all (with the nearest levels
  reachable at this photon energy).

Peaks are maxima of the lightly smoothed EDC that stand 8 noise sigmas
(peak to valley) above their surroundings, the noise measured from the
curve itself and scaled as √counts for counting data; each is refined by a
Gaussian on a sloping background. **Peak threshold** lowers it for weak
lines.

### kz map processing

A photon-energy scan arrives stacked as measured, because nothing in the
files says where each spectrum's Fermi level actually is. **Functions → kz
map processing...** on a kz map's viewer measures it from the spectra
themselves:

1. The slit cut opens alongside with a selection box on it. Drag the box
   over a Fermi edge -- the angle range says which channels are summed into
   the EDC, the energy range is what the fit runs over. The box starts
   around wherever the edge appears to be in the scan as a whole.
2. That same box, **by detector index**, is used on every spectrum. An index
   range is the same channels and the same analyser window everywhere; an
   energy range would mean something different in each, since these are
   exactly the spectra whose energy scales do not yet agree.
3. Each spectrum is shifted to put its own edge at zero, and the cube is
   cropped to the energy range they all still cover.
4. **Normalise each spectrum by its total intensity** (ticked by default)
   divides the beamline out: the photon flux and analyser transmission vary
   by a large factor across a wide scan, so raw counts in a kz map are
   largely a picture of the beamline. Done after cropping, so every total is
   a sum over the same range.

E_F against photon energy is plotted as it goes. A smooth curve means the
fits are sound; scatter means some of them found something other than the
edge. Spectra whose fit is not believed -- the level pinned to the end of
the window, or a step height no larger than its own error bar -- are marked
and filled in from their neighbours rather than dropped.

The fitted levels are stored on the result as `kz.fermi_level_eV`, so the
calibration can be checked or reused later.

### kz -> momentum

**Functions → Convert to momentum (kz, k∥)...** on a kz map's viewer converts it to `(k_z, k_par, E)`
in Å⁻¹, with k_z along the first axis.

The conversion is a coordinate change:

```
k_par = sqrt(A E_kin) sin(alpha)
k_z   = sqrt(A m* (E_kin + V0) - k_par^2)      A = 0.262466 Å⁻² eV⁻¹
```

It is done by **inverting** the map rather than projecting it: for each
target `(k_z, k_par)` there is exactly one photon energy and one emission
angle that produced it, so the conversion is a resampling of the original
regular grid — one `map_coordinates` call per energy slice. Exact for every
slice, and about 120× faster than interpolating a scattered cloud (2.8 s for
a real 41×664×602 cube). Points the measurement never reached come back as
NaN rather than as the nearest sample smeared outwards.

**Flatten the Fermi surface first.** The conversion reads the energy axis
literally, so an edge that bends across the analyser angle becomes a bend in
k_z that looks like dispersion. Open the slit cut, use Fermi-surface
correction, then convert. The dialog measures the bend and says so; it is a
reminder, not a refusal.

**Choosing V₀.** The inner potential is not measured by the experiment, so
the window is built around choosing it:

- the preview reconverts one energy slice as you drag, in milliseconds;
- the zone boundaries are drawn over it from the space group, lattice
  constants, the surface normal and the in-plane direction the slit ran
  along: the real cut of the 3-D zone by the (k_z, k_par) plane, for any
  surface, not only (001). The converted map's **Brillouin zone...** entry
  draws the same cut on its contour;
- **Scan V0** converts one slice at each of a range of inner potentials,
  measures the k_z period each gives, and reports where it equals the
  lattice's. It also reports how well that pins V₀ — for a scan covering a
  couple of zones, only to a few eV, which is the measurement's limit rather
  than the method's;
- **Measure a k_z period...** lets you click two points that are the same
  feature one zone apart. Their separation is matched against every
  reciprocal lattice vector within 15% and the matching planes are listed —
  which says how the crystal cleaved. Centring is handled: in a body-centred
  lattice the period along [001] is 4π/a, not 2π/a, and the answer says so.

### The curve viewer

Every channel of a curve dataset, with error bars where it has a σ channel.
Channels can be hidden, offset into a waterfall (**Offset**), and read off
with the cursor. **Range** is a draggable x window used by the operations
and the panels below.

**Operations** (**Functions → Transform**) each make a new dataset in
the list, with the step recorded:

| operation | what it does | uncertainties |
|---|---|---|
| Crop to the range | keeps the Range | kept |
| Bin | merges *n* neighbours, by sum (counts stay counts) or mean | combined in quadrature: √Σσ², or /n for a mean |
| Normalise | to the maximum, the area, or the mean over the Range; one factor for all channels (default) or each its own | divided by the same factor |
| Subtract a background | constant (mean over the Range), linear (a line through the Range), or Shirley | left as they were |
| Shift the axis | subtracts a value; relabels as E − E_F | kept |
| Add counting errors | σ = √N for each raw-count channel without one | new |

The generic data operations refuse curves, because averaging a σ channel
would shrink it by the wrong factor.

### Curve fit

On one channel, over a fit window (typed, or taken from the viewer's Range):

- **Peaks on a background** -- Lorentzian, Gaussian or Voigt peaks
  (placed by clicking the plot, or **Find peaks**), a none / constant /
  linear / Shirley background, an optional instrumental resolution, and an
  optional Fermi-edge cut-off for EDC peaks near E_F. The solver is the
  MDC/EDC fitter's, so the same peak fitted here and on a cut gives the same
  number. Weighting is the channel's own σ when it has one, √N for raw
  counts, or uniform.
- **Fermi edge** -- the Fermi-level dialog's model. The start is the
  steepest drop in the window, so a band just below E_F does not pull the
  edge to the end of the window; a fit that converged without finding an
  edge (E_F pinned to the window's end, or no significant step) is
  **reported as such** and cannot be used to shift the data. **List a copy
  with E_F = 0** shifts the whole dataset.

Results: every parameter with its error (scaled by √χ²), the residual in
units of σ, **Fit to list** (data, fit, residual and each component as one
curve dataset, parameters in `info`), and **Copy results** (tab-separated).

### Spin analysis

From a spin EDC's viewer. The channel table shows what each `SpinComp` label
says -- manipulator setting, target magnetisation, the axis and sign it
counts as "up" -- and the axis and sign can be overridden if a label is
wrong.

- **Cross ratio** (default when the channels allow it): the geometric mean
  of the +axis channels against that of the −axis channels. With the usual
  four channels -- +Z counted at (<0,0>, +X) and (<90,180>, −X), −Z at
  (<0,0>, −X) and (<90,180>, +X) -- every setting's transmission and every
  magnetisation's target reflectivity appear once on each side and
  **cancel exactly**. The panel says which of the two cancel for the
  channels actually present.
- **One pair** (one manipulator setting): available, and drawn as dashed
  lines beside the cross ratio, but it carries the instrumental asymmetry.
- **Instrumental asymmetry ε** is measured from the two pairs
  (ε = (A₁ − A₂)/2) and reported with its error. On the uploaded file
  `Cut20260528222027_5S.krx` it is ε = +0.0186 ± 0.0005, which is a bias of
  ≈ ε/S ≈ 0.09 on a single-pair P at S = 0.2 -- larger than the polarisation
  there.
- **S_eff** defaults to 0.2, the value published for this end station's
  FERRUM VLEED detector (Sci. Rep. 2023, doi:10.1038/s41598-023-40145-1).
  Use the end station's own calibration where there is one: P and its error
  bar both scale as 1/S.
- **Bin** sums neighbouring points first (they stay counts). **Zero
  reference** optionally takes P = 0 over a window you know is unpolarised;
  it removes real polarisation too if the window has any, and is recorded.
- Counting statistics throughout: σ_A = ½(1 − A²)·√(Σ 1/Iᵢ)/k for k
  channels a side (the textbook √((1 − A²)/N) for one pair), σ_P = σ_A/S,
  and I↑,↓ = I(1 ± P)/2 with their errors. Checked against Monte Carlo:
  the predicted σ_P is 2.3 % below the actual scatter at ~400 counts per
  channel,
  and P is recovered unbiased with a planted ε = 0.05.
- **P to list** (P and σ P), **I↑ / I↓ to list** (both, with σ), and
  **As a figure**.

### De-grid

Removes the periodic pattern the detector prints on every image -- the
hexagonal MCP structure (ANTARES, the CASSIOPEE MBS end station; ~5 px) or
the square mesh (CASSIOPEE Scienta; ~11 px, 10 % contrast). **Do it first**:
only data still on the detector's pixels can be de-gridded, so the button
refuses anything k-converted, Fermi-surface corrected, kz-aligned,
interpolated along a path, smoothed or differentiated, and says why.

**A map or kz map** -- **Functions → Remove detector grid...** on the contour window
or on its slit-cut window; the whole map is done at once. No reference measurement is
needed: the map is its own reference. The grid stays on the same pixels in
every slice while the photoemission moves, so it survives an average over
the slices that the photoemission does not. From that average:

1. what depends on energy only is divided out (the Fermi edge and flat
   bands sit at the same kinetic energy in every slice and are not grid);
2. only the grid's own peaks in k-space are kept -- found automatically,
   shown in red, typically 1.5 % of k-space -- so nothing else in the image
   can be touched;
3. every slice gets its own grid contrast and sub-pixel shift, fitted to it,
   refined locally over 5 × 5 tiles; the grid is re-estimated with the
   slices shifted back into register; each slice is divided by it.

Measured on held-out slices (grid from the other half of the map), the
grid's power over that of the same k-space regions without a grid:

| map | before | after |
|---|---|---|
| ANTARES WSe2 (MBS, hexagonal) | 60 | 1.4 |
| CASSIOPEE Map80eV (Scienta θ map, 28 counts/px) | 9.7 | 1.01 |
| CASSIOPEE LHhv (Scienta hν map) | 448 | 1.5 |

1 is no grid left. The per-slice contrast matters most on a kz map, where it
falls as the count rate rises (their correlation was −0.99 on LHhv). A
113-slice ANTARES map takes about a minute, with a progress bar.

**A cut** -- **Functions → Remove detector grid...** on the cut viewer. Best is a grid from a map
taken with the same lens mode and pass energy: tick **Also list the grid
pattern** when de-gridding the map and a `[grid]` dataset goes to the list;
the cut's window finds it (same detector frame) and chooses it when the
settings match. With no such grid, the cut's own high-frequency grid peaks
are notched out of I / smooth(I): the grid goes completely, and so does the
photoemission in the same k-space regions -- the window says so when the
power left drops below 1. A grid from other settings is offered but not
chosen (pass energy 20 → 50 eV: 88 % of the grid's power removed).

The window shows the slice before and after on one colour scale (and as
−∂²I/∂E², where any grid is obvious), the grid found, its k-space regions,
the contrast and shift of every slice, and the numbers above. **Advanced**
holds the parameters; the defaults are the tested best and are what the
button uses.

### Tools ported from the MATLAB programs

Five tools from `DataAnalysisPrograms_PD` (yulinARPES, Quick_Panel,
PD_program, TB_zph_zmd, Lattice_plot_zqs), each checked against the MATLAB
original and changed where the original could be bettered. The windows that
take several decisions lead through them as numbered steps (*Step 2 of 4:
...*, with Back / Next, and Next greyed with the reason until the step has
what the next one needs; `ui/common/stepper.py`).

**Sum repeated scans** (Datasets box; `add_two_value_arrays.m`). A fine cut
or map taken in batches, so that an ageing sample can be seen, added back
together in three steps: *are these repeats?* (grids and recorded
conditions against the first scan; a slightly different grid is
interpolated onto the first), *has the sample changed?* (each scan against
the sum of the **others**, with counting statistics: χ²/point ≈ 1 for noise
only, flagged when it stands well above the set's median; the count rate;
the energy shift of its EDC, measured by cross-correlating derivatives
against the others aligned by their own shifts), and *sum* (the changed
scans start unticked; shifts undone if asked, the energy axis then cropped
to what every scan covers; sum or mean). The MATLAB tool added two arrays
with no check at all.

**Ridge detection** (Advanced; `ridge_detection.m`). Bands as ridges of the
smoothed intensity: where the Hessian curves down across the band and the
slope across it changes sign, each point placed to a fraction of a pixel.
Changed from the .m file: smoothing per axis in that axis' unit (it used
one radius in pixels of mixed units); derivatives in units of the smoothing
width, so the result does not depend on whether energy is in eV or meV;
Gaussian-derivative filters; a sign-change test instead of an absolute
slope tolerance; vectorised. On a synthetic parabolic band the points sit
within 0.2 meV (median). Results: the ridge image to the list (a map plane
by plane), the points to a CSV.

**Detector corrections** (Basic, under Remove detector grid;
`StretchDitherImage.m`, `InterpRegionAlongX/Y.m`, `inpaint_nans.m`). On the
detector's own pixels, a cut or every slice of a map: a **keystone** (the
lit trapezoid) stretched to a rectangle, as an exact inverse map with each
line's counts kept (the .m file re-gridded a Delaunay triangulation of every
pixel and did not conserve counts); and **bad regions** patched -- linear
along x or y as in MATLAB, or a harmonic (Laplace) fill, factorised once
for a whole map. *Find bad pixels* compares the mean detector image with
its own median filter and its pixel-to-pixel scatter.

**Tight-binding fit** (Calculation; `TB_zph_zmd`). E(k) = e0 + Σ t_n Σ cos(k·R)
over shells of neighbours of any 2-D or 3-D lattice (typed, or the
Brillouin-zone tool's last crystal), fitted to band points traced by ridge
detection on constant-energy planes of a k-map, clicked, or read from a
file. A **linear** least-squares fit with standard errors (the .m file ran
a nonlinear fit from all-ones starting values); robust to stray points
(Tukey reweighting); the azimuth can be refined. One contour alone fixes the
shape and not the energy scale, so t1 is then held (the .m file's fit was
undetermined there); a condition-number check says when the points cannot
tell the shells apart. The model's contours are drawn over the data, and
its spectral function goes to the list as a map.

**Crystal structure** (Calculation; `Lattice_plot_zqs`). The crystal in 3-D
-- typed in (cell, space-group centering, basis), from the Brillouin-zone
tool's last crystal, or from a **CIF** (symmetry operations applied).
Element by element colour, size, transparency and gloss; cell edges and
bonds; a **cleavage plane** (hkl) slid along its normal or put on one of the
listed atomic layers, with the atoms above faded or hidden to show the
termination. Drawn with QPainter (no OpenGL), so it renders at any size;
*To figure composer* puts it on a figure page.
The Brillouin zone in the same view, with k axes as arrows, a coloured path
and a measured k-map painted on the plane through Γ or the zone top.

**Calculated bands** (Calculation). A band structure calculated elsewhere
(Wannier90 `_hr.dat` and the other formats in the guide), with orbital
weights: bands along a path, constant-energy contours, and a comparison
with a measured cut or map -- a renormalisation Z and an energy shift fitted
to the cut's bands. Bands can be followed through crossings by eigenvector
overlap; the calculated bands go to the figure composer as dashed guides.

**Stacked slices** (Advanced). Constant-energy or constant-k_z slices of a
volume as separated, faded layers in 3-D, with curves marked through them.

### Cut arithmetic

Two cuts combined into a third: LH − LV linear dichroism, circular
dichroism, dividing by a reference spectrum (gold), and plain differences,
ratios and sums. Cuts only -- a map has to be sliced first.

Two ways in, the same window either way:

- In a cut's viewer, choose **Functions → Cut arithmetic...**. That cut is **A**. A small
  prompt asks you to click the other cut in the main list; the window opens
  on that click. Clicking a map, or A's own row, is refused in the prompt
  and it keeps waiting; closing the prompt stops the wait.
- In the main list, select two cuts, right-click, **Cut arithmetic on the
  two...**. The first selected is A.

The window shows:

- **What the two measurements recorded**, side by side: polarisation,
  photon energy, theta/tilt/phi, sample X/Y/Z, temperature, pass energy,
  lens mode. Rows in orange differ beyond a tolerance (0.01 eV, 0.05°,
  5 µm, 1 K). A dichroism map between two sample spots or two temperatures
  is a map of *that* difference, so it is worth a look before trusting the
  result. Polarisation is shown but never flagged -- it is meant to differ.
- **Warnings** from the recorded polarisations: A = LV and B = LH under
  the linear-dichroism preset (the sign is reversed -- **Swap A ↔ B**),
  both the same polarisation, or linear cuts under the circular preset.
- **Purpose** presets, each only a starting point:

  | Purpose | Result | Scale B to A |
  |---|---|---|
  | Linear dichroism (LH − LV) | A − B | same total intensity |
  | Circular dichroism | (A − B)/(A + B) | same total intensity |
  | Divide by a reference | A / B, B's angular profile, mean 1 | -- |
  | Custom | any | any |

- **Scale B to A** by total intensity, or by the intensity in a region
  (x0 y0 x1 y1, or **Take the box from A's viewer**). The two polarisations
  leave the undulator with different flux; unscaled, an LD map is mostly
  that ratio. A region above E_F, or a band known not to be dichroic, is
  the better reference when the dichroism itself changes the total.
- **Divide by** the whole reference, its angular profile (summed over
  energy: the detector's channel sensitivity, without gold's own Fermi edge
  and noise), or its energy profile. The reference is normalised to a mean
  of one, so the result keeps A's scale. **Reference floor** hides the
  ratio where the reference is below 5% of its maximum.
- **Hide where A+B below** (2% by default) blanks the background, where an
  asymmetry is a ratio of two noises.
- A and B on **one** colour scale, B as it entered the arithmetic (scaled,
  on A's grid), and the result -- signed results on a red-white-blue scale
  symmetric about zero.

B is interpolated onto A's grid when the two grids differ, and blank where
B does not reach. Axes in different units (degrees against Å⁻¹) are
refused rather than lined up.

**Uncertainty.** When both cuts are raw counts on the same grid, the
Poisson standard deviation is propagated -- `√(A + s²B)` for a difference,
`2s√(AB(A+B))/(A+sB)²` for an asymmetry -- and the report gives its median
and the fraction of pixels beyond 2σ. Checked against Monte Carlo: exact
for the difference; for the asymmetry within 5% from about 20 counts per
pixel, but 12% low at 10 counts and 32% low at 4, which the report says
when it applies (bin the cuts first). There is no σ once either cut has
been processed or interpolated, nor for a ratio.

**Result to list**, **Uncertainty to list** and **All three as a figure**
leave the window. The result is a cut named after A with `LD`, `CD`,
`norm`, `diff`, `asym`, `ratio` or `sum`, and records A, B, the operation,
the scale applied to B and the settings in its processing history.

---

## Figures and export

The **figure composer** is laid out like a manuscript page. With nothing
selected it opens a blank page of empty panels (a)–(d); each panel can be
dragged and resized on the page (snapping to the others' edges), filled
from the main list or with a picture from a file, and given an inset -- a
magnified region of itself (*Zoom-in inset*) or another panel drawn inside
it. Arrows, lines, boxes, ellipses and text are drawn straight on the page
with the tools on its left, and edited (colour, width, fill, text) when
selected; Undo/Redo covers it all. The figure tools beside the page set
the selected panel's axes, fonts, colours and the export; shared colour
scales, labels only on the outer edges and a slice-series page (a row of
constant-energy contours, or of cuts) from one cube work as before. It is
reached from the list (tool panel → Advanced → **Figure composer...**),
from a frozen slice's **Figure composer...** button, or from a map's
**Slice series figure...**. Empty slots are never written to a file.
Several panels can be selected (Ctrl/Shift-click) and lined up, spaced
or matched in size; the letters line up and run (a), (b), (c). **Save
as...** writes the composer's own `.afig` file, which reopens for editing,
and beside it a folder *data for <name>* with each panel's data as shown and
the original measurement it came from.

Export gives what the panel actually shows, in the formats that make sense
for it -- image, data, or both. Data go to a `.nxs` this program reopens, or
**for Igor Pro**: an Igor text file (`.itx`) or a binary wave (`.ibw`,
version 5), with the wave scaling, units, dimension labels and the
metadata in the wave note (`tools/output/igor.py`); an unevenly spaced axis
is written as a wave of its own beside the data. The Datasets box's
**Save** offers the same two formats for whole datasets (cuts, maps,
spatial scans, curves); an Igor export is a copy, so the rows keep pointing
where they did.

---

## The experiment log

Two buttons side by side in the tool panel's *Datasets* box write the
beamtime's log -- one CSV row per measurement, like the hand-kept DLS and
ANTARES logs -- from the files' own metadata.

- **Log format...** sets the log up, once:
  - **Log file** -- a folder and a file name. Pointed at a log that already
    exists, the window shows *that log's* format, to adjust.
  - **Format** -- the separator and the decimal mark. By default they are
    Excel's on this computer: Excel splits a CSV at the system's list
    separator, `;` wherever the decimal mark is a comma (much of Europe), so
    a comma file opens there with every row in the first column. On
    Windows they are read from the system settings; on macOS (Excel for
    Mac, Numbers) and Linux (LibreOffice) from the decimal mark. Numbers
    are written with the decimal mark, so a spreadsheet reads them as
    numbers rather than text or dates. Ranges are written as
    `start/stop/step` or `start:step:stop`.
  - **Columns** -- pick a dataset of the main list as the reference: its
    fields are listed with their values -- the file name, data type, start
    time and the scanned ranges worked out from the data, then every
    metadata field and motor position the file has. Tick a field to add it
    as a column; rename the columns, and put them in any order (Top / Up /
    Down / Bottom). *Add empty column* gives a column to fill in by hand
    (comments). *Save format as...* / *Load format...* keep one per
    beamline.
  - **Save format** keeps it (`~/.arpes_viewer/explog_settings.json`). If
    the log already holds records and its format changed, the log is
    rewritten under the new header, every record moved with its columns
    (followed through renames and moves); columns taken out of the format
    are named first, since what they hold goes. The header is only ever
    the first row.
- **Auto Exp. log** -- one click: the selected datasets are added to the log
  in its format, after its last record (blank rows a spreadsheet left at
  the end are dropped), skipping files already in its file-name column. The
  first click, before any format was saved, opens *Log format...*. *Keep
  recording* (in the format window) logs every file added to the list from
  then on, until the program closes: add each file as it is measured and
  the log keeps itself. A log continued on another computer keeps its own
  separator.
- **What a position column says** depends on the measurement, through each
  column's *If scanned* rule:
  - a cut: the sample stood still -- the motor's value;
  - a SPEM scan: the scanned stage axis as `start/stop/step` (*SPEM:
    scanned X / Y*, or *auto* when the axis carries the motor's name);
  - a map turned with the manipulator rather than the deflector: the angle
    as `start/stop/step` (*map: scanned angle*, or *auto* by name), and a
    sample position the beamline corrected along with the angle as
    `first..last`. Diamond files keep these moving positions (`<motor>
    range` among the motors); a value stored per point is summarised the
    same way;
  - a kz map: the photon energy range (*kz map: photon energy*).

New logs start with a byte-order mark, so that Excel reads hν and Δ. A log
open in Excel cannot be written on Windows: close it there first.

The work is in `tools/output/explog.py` (Qt-free), the window in
`ui/panels/explog.py`.

---

## Nothing is lost: autosave, the log, and recovery

Every derived dataset is written into this run's session folder the moment it
exists (`~/.arpes_viewer/sessions/`). That does two things: losing the
program stops meaning losing the session's work, and a listed dataset can be
a *path* rather than a resident cube, so the memory budget can drop one it is
not using and re-read it on demand.

- Everything except 4-D and 1-D spatial scans is auto-saved (those are the
  measurement itself, gigabytes of it).
- The copies are written **uncompressed** (they are working copies, written
  on the GUI thread after every operation; gzip cost seconds per large cube).
  A dataset's copy is deleted as soon as it is saved properly.
- **A normal close deletes the session folder.** Once you have confirmed
  closing (everything you want is saved, or you chose to close without
  saving), the open viewers are closed and the folder -- with any copies
  recovered from a crashed session that were in the list -- is removed at
  once.
- **Only a crash keeps them.** A folder left by a run that did not close
  normally is kept for three days, offered back at the next start, and then
  deleted.
- **Operations log** -- `~/.arpes_viewer/sessions/operations.log`, plain text,
  one line per store/save/discard. After a crash this is how you find what was
  where.
- **Recovery** -- on startup, a leftover session offers its datasets back.
- **Closing** the main panel with unsaved derived data asks first, and says
  which datasets it means: **Save to a file...**, **Close anyway** (the
  working copies are deleted) or **Cancel**.

---

## Long operations

Reading a map, converting to k-space, resampling a cube, assembling a
CASSIOPEE folder: all of it runs on a worker thread with a progress bar that
moves and a Cancel that works. One job at a time, because HDF5 here is not
thread-safe. A job that finishes within 400 ms never shows a dialog.

---

## Adding a beamline

Readers are added for public large-scale facilities only (see *Terms of
use* at the top); for one of those, or to collaborate, contact
[ding.pei@synchrotron-soleil.fr](mailto:ding.pei@synchrotron-soleil.fr).

Copy `loader/beamlines/soleil.py` (HDF5) or `loader/beamlines/cassiopee.py`
(text), implement four things, and add the module's name to
`_install_default_loaders()` in `loader/registry.py`. Nothing else changes.

```python
class MyBeamlineLoader(Loader):
    name = "MY BEAMLINE"
    patterns = ("*.dat",)
    priority = 20                      # detection order; higher goes first

    def can_open(self, path) -> bool: ...
    def list_entries(self, path) -> list: ...
    def load(self, path, entry=None, progress=None) -> NxsScan: ...

register(MyBeamlineLoader())
```

A beamline that writes Scienta SES files needs less: subclass
`loader.formats.ses.SESLoader`, set `name`, `facility`, `prefix` (and
`reads_text = True` for `.txt`), and say which headers are yours in
`claims(header)` -- see `loader/beamlines/solaris.py`.

`progress` is optional and is `callable(done, total, label)`. Accepting it
declares that the reader expects to be slow, which is what routes it through
the worker thread; leave it off for anything that opens in one call.

Return an `NxsScan` with the axis slots its `kind` names in `AXIS_SLOTS`, and
set `info["axis0.role"]` if the reader can tell what its first axis is.

---

## Layout of the code

```
ARPES_viewer.py      the launcher: the dataset list, and everything it opens
version.py           the version number and date shown in the launcher
ARPES viewer.app/    the macOS launcher: double-click to start   (bash)
loader/              reading data in, and keeping it              (no Qt)
  beamlines/         one reader per beamline
  formats/           what several readers share: NeXus layouts, Scienta SES
tools/               the algorithms                               (no Qt)
  momentum/          angle and photon energy to momentum
  crystal/           lattices, space groups, Brillouin zones
  processing/        2-D / 3-D operations, combining, de-gridding
  spectra/           EDCs/MDCs, spin, fits, core levels
  output/            export, the figure model, colormaps
  runtime/           memory, system stand-ins
ui/                  the windows                                  (Qt)
  launcher/          the main window, the Loader window, the list's menu
  viewers/           the per-dataset viewers and their widgets
  panels/            the tool windows a viewer opens
  common/            long jobs, messages, the GC guard
devtools/            one-off generators, run by hand
docs/                the user guides (en, zh, fr), their screenshots and PDFs
test/                the tests: loaders/, tools/, ui/ and the code checks
```

The top-level split is by *dependency*, which is why it is worth having:
`tools/` and `loader/` are importable from a plain script with no display,
and are tested that way; `ui/` is the only package that imports Qt. A
module's folder tells you what it may depend on. Inside each, the
subpackages group modules by what they are for.

| File | Role |
|---|---|
| `ARPES_viewer.py` | The launcher: the dataset list, the session folder, and the wiring between them. Run this. |
| **`loader/`** | |
| `registry.py` | Which beamline wrote a file and who reads it, plus the load-time axis order and scanned-axis role. |
| `nxs_file.py` | The SOLEIL parser, the saved format's reader and writer, `NxsScan`, the axis-slot table, and the lazy arrays that keep a measurement on disk. Start here. |
| `native.py` | This program's own saved format -- detected, never chosen. |
| `cubecache.py` | `CubeCache`: a map cube read into memory once it is being browsed, so the energy slider and the cut sliders slice an array rather than the file. |
| `session.py` | Where a computed dataset lives: the folder it is written to as soon as it exists, the memory budget, and the operations log. |
| `formats/nexus_arpes.py` | What the NeXus readers share: the layout, the kind it opens as, kz-map energy referencing, the lazy `NexusView`, the metadata. |
| `formats/ses.py` | Scienta SES files: Igor `.ibw`, the map `.zip`, multi-region `.txt`, the run-mode table, and `SESLoader`, which each SES beamline subclasses. |
| `beamlines/soleil.py` | SOLEIL ANTARES. **Copy this file to add an HDF5 beamline.** |
| `beamlines/cassiopee.py` | SOLEIL CASSIOPEE: Scienta text spectra, and folders of them assembled into a map or a kz map. |
| `beamlines/cassiopee_spin.py` | CASSIOPEE's spin end station: MBS A-1 `.krx` and `.txt` files, as a cut, a map or a spin EDC. |
| `beamlines/maestro.py` | ALS MAESTRO (BL 7.0.2) `.h5` files: cuts, deflector maps, XY spatial scans, two-motor lines and focus scans, read lazily. |
| `beamlines/alba.py` | ALBA LOREA (BL20) `.nxs` files. |
| `beamlines/diamond.py` | Diamond I05 (HR-ARPES) and I05-1 (nano-ARPES) `.nxs` files, as two loaders. |
| `beamlines/maxiv.py`, `beamlines/ape.py` | MAX IV Bloch and Elettra APE (SES). Bloch is also the fallback for any SES file no other beamline claims. |
| `beamlines/bessy.py`, `beamlines/solaris.py` | BESSY One Cube / One Square and SOLARIS URANOS (SES), told apart by the header's `Location` and `Instrument`. |
| `beamlines/sgm4.py` | ASTRID2 SGM4 `.h5` files, read lazily. Ported from peaks. |
| **`tools/`** | |
| `momentum/kspace.py`, `momentum/cutk.py` | Angle-to-momentum conversion, for a map and for a single cut, and the rotation read off a picked direction. |
| `momentum/kzmap.py` | Calibrating a photon-energy scan against its own Fermi edges: fit per spectrum, align, crop, normalise. |
| `momentum/kzconv.py` | The k_z conversion: the analytic forward and inverse maps, the resampling, and the inner-potential scan. |
| `crystal/lattice.py`, `crystal/spacegroups.py`, `crystal/bz3d.py`, `crystal/bz2d.py`, `crystal/moire.py`, `crystal/highsym.py` | The Brillouin-zone geometry: lattices and point groups, the 3-D and 2-D zones and their irreducible wedges, the moiré zone, and the named high-symmetry points (a static table per Bravais lattice, after Setyawan & Curtarolo 2010). |
| `crystal/cleavage.py` | From a measured k_z period to which lattice planes could have produced it. |
| `processing/analysis.py` | Arbitrary-direction cut and Fermi-surface correction. |
| `processing/dataops.py` | Truncate / self-normalise / compress, and the axis tables derived from `loader/nxs_file.py`. |
| `processing/stackcuts.py`, `processing/combine.py` | Combining datasets: stacking cuts along a new axis, and joining pieces end to end along one of their own axes. |
| `processing/process.py`, `processing/volume.py` | Smoothing, derivatives, curvature, backgrounds, symmetrisation -- in 2-D and over a cube. |
| `processing/degrid.py` | Removing the detector grid: the map's own grid estimate, the per-slice contrast / shift / local fit, the single-cut paths, and the pixel-lock check. |
| `processing/repeats.py` | Summing repeated scans: grid check, χ² of each scan against the others, energy shifts, the sum. |
| `processing/ridge.py` | Ridge detection, on a plane or plane by plane through a cube. |
| `processing/detector.py` | Keystone correction, bad-region patches, finding bad pixels, and the on-the-detector check. |
| `crystal/tightbinding.py` | Neighbour shells, the linear (and robust) tight-binding fit, the azimuth scan, reading band points. |
| `crystal/structure.py`, `crystal/last_crystal.py` | Crystal structures: CIF reading and symmetry expansion, cells of atoms, Miller planes and atomic layers; and the Brillouin-zone tool's last crystal, for the tools that reuse it. |
| `processing/cutops.py` | Arithmetic between two cuts: scaling, resampling, dichroism, reference division, Poisson uncertainty, and which recorded conditions differ. |
| `spectra/curves.py` | One-dimensional data as tables: channel names, the σ convention, and the σ-aware crop / bin / normalise / background operations. |
| `spectra/spin.py` | Spin channels from their labels, the cross ratio, the instrumental asymmetry, polarisation and spin-resolved spectra with counting errors. |
| `spectra/fermi.py`, `spectra/peaks.py`, `spectra/dispersion.py` | The Fermi-edge model, MDC/EDC peak fitting, and what the fitted band says (`v_F`, `m*`, self-energy). |
| `spectra/corelevel.py`, `spectra/corelevel_table.py` | The core-level check: the integrated EDC, the energy reference, peak finding, assignment against Table 1-1 and typical chemical states, contaminants, higher-order light, the rigid-shift test. The table is generated. |
| `output/figure.py`, `output/export.py` | What a publication figure *is*, and what leaves the program. |
| `output/colormaps.py` | Generated colormap tables. |
| `output/igor.py` | Igor Pro export: `.itx` text and `.ibw` binary waves with scaling, units and notes. |
| `output/explog.py` | The experiment log: its format, each dataset's row (scanned ranges included), continuing a CSV log after its last record, and reformatting it in place. |
| `runtime/memory.py`, `runtime/system.py` | Memory use, and the stand-ins for the lab's `tools_packages`. |
| **`ui/`** | |
| `launcher/main_window.py` | The launcher window's layout. |
| `launcher/loader_dialog.py` | The Load-data window: files, reader, axis options. |
| `launcher/list_actions.py` | Which right-click entries and tool-panel buttons of the main list fit a selection (by kind and row count), Qt-free. |
| `viewers/windows.py`, `viewers/widgets.py` | The per-file viewers, and the image panels they are built from. |
| `viewers/curves.py` | The curve viewer, its curve fit panel and the spin analysis panel. |
| `viewers/cursorlink.py` | The readout cursor shared by a map's contour and its cut windows (positions only; widths stay in each window). |
| `viewers/functions_menu.py` | The Functions menu shared by every viewer (`add_function`). |
| `viewers/plot_tools.py` | Plot tools: guide lines and labels on a live viewer, and drawing them into an image export. |
| `panels/figure.py`, `panels/fit.py`, `panels/process.py`, `panels/volume.py` | The figure composer, the MDC/EDC fit panel, and the 2-D and 3-D processing panels. |
| `panels/kzmap.py`, `panels/kzconv.py` | The Align Fermi levels (kz map) window, and the kz-to-momentum window: V0, the zone lines, the period tool. |
| `panels/stackcuts.py`, `panels/cutops.py` | The Combine datasets window (stack or join), and the cut arithmetic window. |
| `panels/corelevel.py`, `panels/degrid.py` | The core-level check (elements, energy reference, results window), and the Remove detector grid window. |
| `panels/explog.py` | The experiment log's format window, and logging the selection in that format. |
| `panels/repeats.py`, `panels/ridge.py`, `panels/detector.py` | Sum repeated scans, ridge detection, and the detector corrections. |
| `panels/tightbinding.py`, `panels/crystal.py` | The tight-binding fit and the crystal-structure viewer (the Calculation tab). |
| `common/stepper.py` | The numbered steps with Back / Next that the multi-step tool windows are built on. |
| `common/jobs.py` | Running one long operation at a time off the GUI thread. |
| `common/gcguard.py` | Cyclic garbage collection from the event loop only, never inside a Qt call (prevents random crashes on opening and closing windows). |
| `common/notify.py` | The message boxes. |
| **`devtools/`, `test/`** | |
| `devtools/gen_colormaps.py`, `devtools/gen_spacegroups.py`, `devtools/gen_corelevels.py` | Generate the tables in `tools/`; need scipy + the lab's `.mat`, `spglib`, and poppler's `pdftotext` with `Table_1-1.pdf`, none of which the program itself does. |
| `devtools/demo_data.py`, `devtools/guide_screenshots.py`, `devtools/build_guide_pdfs.py` | Synthetic data to try the program on, the script that re-takes the user guides' screenshots from it (offscreen), and the one that prints the guides to PDF (pandoc + Playwright). |
| `test/loaders/`, `test/tools/`, `test/ui/` | The tests of each package. |
| `test/test_imports.py` | Checks on the shape of the code: that no local shadows an imported module, that no package is referenced by its bare name, that everything imports, that every package has a docstring, and that `loader/` and `tools/` stay Qt-free. |

---

## Assumptions worth knowing about

### Axis order in SOLEIL `.nxs` files

MATLAB's `h5read` reports HDF5 array dimensions **in reverse order** compared
to `h5py` (documented MATLAB behaviour, not a bug), so
`load_soleil_nxs.m`'s `permute(value, [3 4 2 1])` cannot be copied
index-for-index into Python.

Rather than guessing a fixed transpose, `loader.nxs_file.align_and_transpose`
matches each raw array axis to a physical axis (x/y/k/E/deflx) **by length**
against the independently-read calibration arrays. Robust unless two axes
have the same number of points, in which case it warns and picks one.

Before trusting the output on a new set of files:

```python
from loader.nxs_file import inspect_nxs
inspect_nxs("your_file.nxs")
```

and compare the printed dataset shapes to the axis lengths of that file's
`NxsScan`. If `align_and_transpose` ever warns about an ambiguous match, that
file is where to look first.

### k conversion

For a `map`, the analyser's own slit-direction axis (`scan.k`) is already
delivered in Å⁻¹ by the MBS acquisition software. The deflector axis
(`scan.x`) is a raw angle and is converted on demand:

```
k_parallel [A^-1] = 0.5123 * sqrt(KE[eV]) * sin(angle[rad])
```

using `center_ke` from the file as the reference kinetic energy. This ignores
the inner-potential/refraction correction and treats KE as constant across
the energy window -- fine for a narrow window around a feature, worth
checking for a wide one (pass an explicit `kinetic_energy_eV` array, one
value per energy channel, to `to_kspace_cube()` for per-channel accuracy).

### CASSIOPEE energy axes

The axis scales in a Scienta text file are printed to a fixed number of
*significant figures*, so the rounding grid gets coarser as an axis crosses a
power of ten (`99.998` to three decimals, `100.01` to two). Axes are
straightened to the line through them; the check that reports a genuinely
curved scale allows for that printing grid, or it would fire on every file
whose energy window crosses 100 eV.

A kz map's absolute binding energies are only as good as the nominal photon
energies and the tabulated work function. On real scans these do not always
agree -- fit a Fermi edge per photon energy before reading absolute binding
energies off one.

### Unrecognised layouts

An unknown top-level group name raises with a clear message instead of
silently mis-parsing. Add the name to `GROUP_NAME_TO_CASE` once you know
which case it matches; new beamline-software versions add new names, and this
has already happened once ("AuK 2022" vs "2021").

---

## Tests

```
python -m pytest
```

from this folder. 859 of them, no display needed -- `conftest.py` puts the
project root on `sys.path` and pins `QT_QPA_PLATFORM=offscreen`, and
`test/conftest.py` runs garbage collection only between tests, the rule
`ui/common/gcguard.py` applies in the program.
