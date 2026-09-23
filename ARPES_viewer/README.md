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

Right-click a row (or a selection) for everything that operates on data:
open, save, rename, remove, data operations, processing, the 3-D view, the
fit panels, export, and **Cut arithmetic on the two...** (see below).

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
view controls (axis ranges, Auto, Inv, Reset, Grid) side by side, and a
**Functions** menu. A map's window also keeps its **Deflector cut** and
**Slit cut** buttons there, since the cuts are how a map is navigated.

Everything else a viewer can do is in **Functions**, in sections:

| Viewer | Analysis | Data operations | Visualization | Slice |
|---|---|---|---|---|
| map / kz map / k-map | | Arbitrary cut, Map k conversion (angle maps), kz map processing and kz -> momentum (kz maps), De-grid map (maps as measured) | Brillouin zone (k-maps), Slice figure | Save slice to the main list, Open slice in a new panel |
| cut | Fermi level, MDC / EDC fit | FS correction, Cut k conversion, Cut arithmetic, De-grid | | same |
| slit cut of a map | | FS correction, De-grid map | | same |
| curve (EDC, MDC, spin EDC) | Curve fit, Spin analysis (spin EDCs) | Crop, Bin, Normalise, Subtract a background, Shift the axis, Add counting errors (√N) | As a figure | |

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
  stand it vertical.
- **Fermi level** -- fit the edge, offset the energy axis, or divide it out.
- **MDC / EDC fitting** -- fit the peaks across a cut and read the band off
  them: `v_F`, `m*`, and the self-energy. Only accepts a cut whose in-plane
  axis is already in Å⁻¹, because every quantity built on the fitted band
  would otherwise carry meaningless units.
- **Arbitrary cuts** -- any straight line through a contour.
- **Fermi-surface correction**, **stack plots**.
- **Cut arithmetic** -- linear and circular dichroism, dividing by a
  reference, differences, ratios and sums between two cuts; see below.
- **Processing** (2-D and over a cube) -- smoothing, derivatives, curvature,
  backgrounds, symmetrisation, despiking.
- **Data operations** -- truncate, self-normalise, compress, on several
  datasets at once when they really are the same format.
- **Brillouin zone** -- overlay a conventional or irreducible zone, or a
  moiré zone, on a converted contour. Space group first, then the lattice
  parameters it constrains; a 3-D preview shows how the cut plane sits.
- **kz map processing** and **kz -> momentum** -- see below.
- **The curve viewer**, its **curve fit** and **spin analysis** -- see below.
- **De-grid** -- remove the detector's grid (MCP pattern or mesh) from a map
  or a cut; see below.

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

**Functions → kz -> momentum...** on a kz map's viewer converts it to `(k_z, k_par, E)`
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
  constants and the surface normal;
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

**Operations** (**Functions → Data operations**) each make a new dataset in
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

**A map or kz map** -- **Functions → De-grid map...** on the contour window
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

**A cut** -- **Functions → De-grid...** on the cut viewer. Best is a grid from a map
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
| `loader/cassiopee_spin.py` | CASSIOPEE's spin end station: MBS A-1 `.krx` and `.txt` files, as a cut, a map or a spin EDC. |
| `loader/native.py` | This program's own saved format -- detected, never chosen. |
| `loader/nxs_file.py` | The SOLEIL parser, the saved format's reader and writer, `NxsScan`, the axis-slot table, and the lazy arrays that keep a measurement on disk. Start here. |
| `loader/session.py` | Where a computed dataset lives: the folder it is written to as soon as it exists, the memory budget, and the operations log. |
| **`tools/`** | |
| `tools/kspace.py`, `tools/cutk.py` | Angle-to-momentum conversion, for a map and for a single cut, and the rotation read off a picked direction. |
| `tools/analysis.py` | Arbitrary-direction cut and Fermi-surface correction. |
| `tools/dataops.py` | Truncate / self-normalise / compress, and the axis tables derived from `loader/nxs_file.py`. |
| `tools/process.py`, `tools/volume.py` | Smoothing, derivatives, curvature, backgrounds, symmetrisation -- in 2-D and over a cube. |
| `tools/kzmap.py` | Calibrating a photon-energy scan against its own Fermi edges: fit per spectrum, align, crop, normalise. |
| `tools/kzconv.py` | The k_z conversion: the analytic forward and inverse maps, the resampling, and the inner-potential scan. |
| `tools/curves.py` | One-dimensional data as tables: channel names, the σ convention, and the σ-aware crop / bin / normalise / background operations. |
| `tools/spin.py` | Spin channels from their labels, the cross ratio, the instrumental asymmetry, polarisation and spin-resolved spectra with counting errors. |
| `tools/degrid.py` | Removing the detector grid: the map's own grid estimate, the per-slice contrast / shift / local fit, the single-cut paths, and the pixel-lock check. |
| `tools/cutops.py` | Arithmetic between two cuts: scaling, resampling, dichroism, reference division, Poisson uncertainty, and which recorded conditions differ. |
| `tools/cleavage.py` | From a measured k_z period to which lattice planes could have produced it. |
| `tools/fermi.py`, `tools/peaks.py`, `tools/dispersion.py` | The Fermi-edge model, MDC/EDC peak fitting, and what the fitted band says (`v_F`, `m*`, self-energy). |
| `tools/lattice.py`, `tools/spacegroups.py`, `tools/bz3d.py`, `tools/bz2d.py`, `tools/moire.py` | The Brillouin-zone geometry: lattices and point groups, the 3-D and 2-D zones and their irreducible wedges, and the moiré zone. |
| `tools/figure.py`, `tools/export.py` | What a publication figure *is*, and what leaves the program. |
| `tools/colormaps.py`, `tools/system.py` | Generated colormap tables; the stand-ins for the lab's `tools_packages`. |
| **`ui/`** | |
| `ui/main_window.py` | The launcher window's layout. |
| `ui/windows.py`, `ui/widgets.py` | The per-file viewers, and the image panels they are built from. |
| `ui/functions_menu.py` | The Functions menu shared by every viewer (`add_function`). |
| `ui/gcguard.py` | Cyclic garbage collection from the event loop only, never inside a Qt call (prevents random crashes on opening and closing windows). |
| `ui/cursorlink.py` | The readout cursor shared by a map's contour and its cut windows (positions only; widths stay in each window). |
| `ui/loader_dialog.py` | The Load-data window: files, reader, axis options. |
| `ui/kzmap.py` | The kz map processing window. |
| `ui/kzconv.py` | The kz-to-momentum window: V0, the zone lines, the period tool. |
| `ui/cutops.py` | The cut arithmetic window. |
| `ui/degrid.py` | The De-grid window. |
| `ui/curves.py` | The curve viewer, its curve fit panel and the spin analysis panel. |
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

from this folder. 509 of them, no display needed -- `conftest.py` puts the
project root on `sys.path` and pins `QT_QPA_PLATFORM=offscreen`, and
`test/conftest.py` runs garbage collection only between tests, the rule
`ui/gcguard.py` applies in the program.
