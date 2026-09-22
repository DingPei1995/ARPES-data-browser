# ARPES viewer

A Python/PyQt5 program for ARPES data: browse the measurements in a folder,
open them in their own windows, convert to momentum, cut, fit, process, and
lay the results out as a publication figure.

It began as a rewrite of the lab's `h5Loader*.py` for SOLEIL ANTARES `*.nxs`
files, and reading those is now one loader among several rather than the
whole program. It currently reads **SOLEIL ANTARES** (NeXus/HDF5), **SOLEIL
CASSIOPEE** (Scienta SES text, single spectra and whole folders), and its own
saved format. Adding a beamline is a new file in `loader/`.

The round-by-round history, and the reasoning behind the choices below, is in
[CHANGELOG.md](CHANGELOG.md).

```
python ARPES_viewer.py
```

Needs Python 3.9+, PyQt5, pyqtgraph, numpy, h5py, scipy.

---

## Contents

- [Getting data in](#getting-data-in)
- [What the main panel is](#what-the-main-panel-is)
- [Data kinds](#data-kinds)
- [The viewer windows](#the-viewer-windows)
- [Analysis](#analysis)
- [Figures and export](#figures-and-export)
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
cube, with the individual spectra available below it.

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

---

## What the main panel is

A browser, deliberately. One row per *dataset*, not per file. Click a row to
read its metadata, **double-click** to open it in its own window.

Right-click a row (or a selection) for everything that operates on data:
open, save, rename, remove, data operations, processing, the 3-D view, the
fit panels, export, and comparison.

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
| `spem_1d` | (x, k, E) | a real-space line scan |
| `spem_4d` | (x, y, k, E) | a real-space raster scan |

`map`, `k_map` and `kz_map` are the three-axis cubes (`CUBE_KINDS`): they
share a viewer, the 3-D view and the processing panels, and differ only in
what their first axis means. A `kz_map` is not offered a k conversion --
turning a photon energy into k_z needs the inner potential, which is a
property of the sample and not of the measurement.

Which axes a kind has lives in **one** table, `loader/nxs_file.AXIS_SLOTS`.
Everything else -- the save format, the memory data wrapper, the processing
panels, which axis is the energy -- derives from it.

---

## The viewer windows

- **A cut** opens as the E-vs-k spectrum, with EDC/MDC curves and a
  draggable slice.
- **A cube** opens on its constant-energy contour. Two buttons open the
  orthogonal cuts, each in a further window.
- **A spatial scan** opens the real-space map together with the E-vs-k
  spectrum at the cursor.

Every image panel has the same view controls in a strip rather than buried in
a right-click menu: colormap, flip, gamma, level bar, interpolation, and the
axis ranges.

---

## Analysis

- **k conversion** (cube and single cut) -- the free-electron formula, with
  the inputs shown and checkable. Refused, with a reason, for an axis that is
  not an emission angle.
- **Fermi level** -- fit the edge, offset the energy axis, or divide it out.
- **MDC / EDC fitting** -- fit the peaks across a cut and read the band off
  them: `v_F`, `m*`, and the self-energy. Only accepts a cut whose in-plane
  axis is already in Å⁻¹, because every quantity built on the fitted band
  would otherwise carry meaningless units.
- **Arbitrary cuts** -- any straight line through a contour.
- **Fermi-surface correction**, **stack plots**, **comparison windows**.
- **Processing** (2-D and over a cube) -- smoothing, derivatives, curvature,
  backgrounds, symmetrisation, despiking.
- **Data operations** -- truncate, self-normalise, compress, on several
  datasets at once when they really are the same format.
- **Brillouin zone** -- overlay a conventional or irreducible zone, or a
  moiré zone, on a converted contour. Space group first, then the lattice
  parameters it constrains; a 3-D preview shows how the cut plane sits.
- **kz map processing** -- see below.

### kz map processing

A photon-energy scan arrives stacked as measured, because nothing in the
files says where each spectrum's Fermi level actually is. **kz map
processing...** on a kz map's viewer measures it from the spectra
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

---

## Figures and export

The **figure composer** lays panels out for a paper: shared colour scales,
labels only on the outer edges, and a slice-series page (a row of
constant-energy contours, or of cuts) from one cube.

Export gives what the panel actually shows, in the formats that make sense
for it -- image, data, or both.

---

## Nothing is lost: autosave, the log, and recovery

Every derived dataset is written into this run's session folder the moment it
exists (`~/.arpes_viewer/sessions/`). That does two things: losing the
program stops meaning losing the session's work, and a listed dataset can be
a *path* rather than a resident cube, so the memory budget can drop one it is
not using and re-read it on demand.

- Everything except 4-D and 1-D spatial scans is auto-saved (those are the
  measurement itself, gigabytes of it).
- Session folders are deleted three days after they are written, or as soon
  as the dataset in them has been saved properly.
- **Operations log** -- `~/.arpes_viewer/sessions/operations.log`, plain text,
  one line per store/save/discard. After a crash this is how you find what was
  where.
- **Recovery** -- on startup, a leftover session offers its datasets back.
- **Closing** the main panel with unsaved derived data asks first, and says
  which datasets it means.

---

## Long operations

Reading a map, converting to k-space, resampling a cube, assembling a
CASSIOPEE folder: all of it runs on a worker thread with a progress bar that
moves and a Cancel that works. One job at a time, because HDF5 here is not
thread-safe. A job that finishes within 400 ms never shows a dialog.

---

## Adding a beamline

Copy `loader/soleil.py` (HDF5) or `loader/cassiopee.py` (text), implement
four things, and add the module to `_install_default_loaders()` in
`loader/registry.py`. Nothing else changes.

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

`progress` is optional and is `callable(done, total, label)`. Accepting it
declares that the reader expects to be slow, which is what routes it through
the worker thread; leave it off for anything that opens in one call.

Return an `NxsScan` with the axis slots its `kind` names in `AXIS_SLOTS`, and
set `info["axis0.role"]` if the reader can tell what its first axis is.

---

## Layout of the code

```
ARPES_viewer.py     the launcher: the dataset list, and everything it opens
loader/             reading data in, and keeping it            (no Qt)
tools/              the algorithms                             (no Qt)
ui/                 the windows                                (Qt)
devtools/           one-off generators, run by hand
test/               the tests
```

The split is by *dependency*, which is why it is worth having: `tools/` and
`loader/` are importable from a plain script with no display, and are tested
that way; `ui/` is the only package that imports Qt. A module's folder tells
you what it may depend on.

| File | Role |
|---|---|
| `ARPES_viewer.py` | The launcher: the dataset list, the session folder, and the wiring between them. Run this. |
| **`loader/`** | |
| `loader/registry.py` | Which beamline wrote a file and who reads it, plus the load-time axis order and scanned-axis role. |
| `loader/soleil.py` | The SOLEIL/ANTARES reader. **Copy this file to add an HDF5 beamline.** |
| `loader/cassiopee.py` | The SOLEIL/CASSIOPEE reader: Scienta text spectra, and folders of them assembled into a map or a kz map. |
| `loader/native.py` | This program's own saved format -- detected, never chosen. |
| `loader/nxs_file.py` | The SOLEIL parser, the saved format's reader and writer, `NxsScan`, the axis-slot table, and the lazy arrays that keep a measurement on disk. Start here. |
| `loader/session.py` | Where a computed dataset lives: the folder it is written to as soon as it exists, the memory budget, and the operations log. |
| **`tools/`** | |
| `tools/kspace.py`, `tools/cutk.py` | Angle-to-momentum conversion, for a map and for a single cut. |
| `tools/analysis.py` | Arbitrary-direction cut and Fermi-surface correction. |
| `tools/dataops.py` | Truncate / self-normalise / compress, and the axis tables derived from `loader/nxs_file.py`. |
| `tools/process.py`, `tools/volume.py` | Smoothing, derivatives, curvature, backgrounds, symmetrisation -- in 2-D and over a cube. |
| `tools/kzmap.py` | Calibrating a photon-energy scan against its own Fermi edges: fit per spectrum, align, crop, normalise. |
| `tools/fermi.py`, `tools/peaks.py`, `tools/dispersion.py` | The Fermi-edge model, MDC/EDC peak fitting, and what the fitted band says (`v_F`, `m*`, self-energy). |
| `tools/lattice.py`, `tools/spacegroups.py`, `tools/bz3d.py`, `tools/bz2d.py`, `tools/moire.py` | The Brillouin-zone geometry: lattices and point groups, the 3-D and 2-D zones and their irreducible wedges, and the moiré zone. |
| `tools/figure.py`, `tools/export.py` | What a publication figure *is*, and what leaves the program. |
| `tools/colormaps.py`, `tools/system.py` | Generated colormap tables; the stand-ins for the lab's `tools_packages`. |
| **`ui/`** | |
| `ui/main_window.py` | The launcher window's layout. |
| `ui/windows.py`, `ui/widgets.py` | The per-file viewers, and the image panels they are built from. |
| `ui/loader_dialog.py` | The Load-data window: files, reader, axis options. |
| `ui/kzmap.py` | The kz map processing window. |
| `ui/jobs.py` | Running one long operation at a time off the GUI thread. |
| `ui/figure.py`, `ui/fit.py`, `ui/process.py`, `ui/volume.py` | The figure composer, the MDC/EDC fit panel, and the 2-D and 3-D processing panels. |
| **`devtools/`, `test/`** | |
| `devtools/gen_colormaps.py`, `devtools/gen_spacegroups.py` | Generate the tables in `tools/`; need scipy + the lab's `.mat`, and `spglib`, neither of which the program itself does. |
| `test/test_imports.py` | Checks on the shape of the code: that no local shadows an imported module, that no package is referenced by its bare name, that everything imports, and that `loader/` and `tools/` stay Qt-free. |

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

from this folder. 302 of them, no display needed -- `conftest.py` puts the
project root on `sys.path` and pins `QT_QPA_PLATFORM=offscreen`.
