# Reading data: technical reference

How the loaders read each facility's files, the kinds of dataset they make, how to add a beamline, and the assumptions the file formats force on them. The code is in this folder (`beamlines/` holds one reader per facility). For how to *use* the program, see the [user guide](../docs/USER_GUIDE_en.md).

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
