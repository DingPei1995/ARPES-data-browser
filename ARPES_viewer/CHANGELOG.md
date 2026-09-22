# Changelog

Round-by-round history of ARPES viewer, newest first. Each entry says what
changed and, more usefully, *why* -- the measurement that settled a format
choice, the bug that a check now prevents, the assumption in the MATLAB
original that turned out not to hold on real data.

What the program does and how to use it is in [README.md](README.md); this
file is the record of how it got that way.

---

## Thirty-fourth round: arithmetic between two cuts, and one name for it

A window for combining two cuts -- LH − LV linear dichroism, circular
dichroism, dividing by a reference spectrum -- reached from a **Cut
arithmetic...** button on a cut's viewer, or from the main list's
right-click menu with two cuts selected. Cuts only, as asked.

### "Compare the two..." was already arithmetic

The main list's right-click **Compare the two...** opened a dialog that
showed the two side by side -- and also computed A − B, A / B and
(A − B)/(A + B), and exported the result to the list. The name promised a
look; the dialog did the arithmetic. Next to a new arithmetic window the two
would have been two answers to one question, and sooner or later two
different scale factors for it.

So there is one engine, `tools/cutops.py`, and one window, `ui/cutops.py`.
The right-click entry is now **Cut arithmetic on the two...** and opens the
same window the viewer button does. What Compare did well is kept: A and B
on a common colour scale, the result beside them, all three as a figure.
`CompareDialog` is gone from `ui/process.py`.

### Choosing B by clicking the list

The viewer's cut is A. Rather than a second picker listing the same rows,
a small prompt asks for a click in the main list, and the window opens on
it. The list's own record is checked *before* anything is read, so a stray
click on a sixty-spectrum map costs nothing and is refused in the prompt,
which keeps waiting. Closing the prompt, or the viewer, ends the wait.

Refusing A's own row needed the row, not the object: a cut read from a
file is re-read on every load, so the same row hands back a different
object each time, and an identity check let "A − A" through. Viewers now
carry the list key they were opened from.

### What the window checks

- Recorded conditions side by side -- photon energy, angles, sample
  position to 5 µm, temperature, pass energy, lens mode -- with the
  differing ones in orange. A dichroism map between two spots is a map of
  the two spots.
- The polarisations against the purpose: LV − LH (swap for the usual sign),
  the same polarisation twice, linear cuts under the circular preset. `LCP`
  is circular, although it starts with an L; circular is tested first.
- Axis units: degrees against Å⁻¹ is refused rather than lined up.
- B is scaled to A by total intensity or by a region, because the two
  polarisations leave the undulator with different flux and an unscaled LD
  map is mostly that ratio. A reference is normalised to a mean of one, so
  dividing by gold does not rescale the data by gold's count rate.

### The uncertainty, and where it stops being right

For raw counts on a common grid the Poisson σ is propagated to first order.
Against 400-draw Monte Carlo, the difference's σ agrees exactly (0.999).
The asymmetry's is within 5% from about 20 counts per pixel, then 12% low
at 10 counts and 32% low at 4: the estimate is built from the noisy counts,
and the ratio stops being Gaussian. The report says so below 20 counts, since
that is also where the "fraction beyond 2σ" is too generous. Processed or
interpolated cuts get no σ at all -- Poisson statistics through a smoothed
map would understate the noise.

428 tests (25 new); a headless run drives the viewer button → prompt → list
click (map and own row refused) → window → presets, swap, export, figure,
and the right-click route, on a real CASSIOPEE LH cut paired with a
synthetic LV partner (the beamtime data are all LH).

---

## Thirty-third round: a picked direction is a line, not an arrow

**Set rotation from contour** returned an angle that depended on which of
the two points was clicked first, and on which side of the origin they sat.
Both give the same *line*, and the sample does not know the difference -- so
the two orders differed by exactly 180 degrees, which is a plausible-looking
number that converts the map upside down.

`points_to_azimuth` folded its answer to (-180, 180], which is the range of
a *direction*. A line has a 180-degree period, so the range is (-90, 90]:
the smaller of the two turns, clockwise or anticlockwise, that stand it
along ky. Standing a direction along +ky and along -ky are the same
alignment of the same axis.

### The tie needed more than the interval

A horizontal line wants exactly a quarter turn, and a quarter turn either
way is equally small -- a real tie. `arctan2` returns that case as 90 plus
or minus a part in 1e14, depending on the sign of a sine that should have
been zero, and those land on opposite ends of a half-open interval. So the
same horizontal line could come back as +90 or -90 depending on floating
point. Results within a tolerance of the boundary are now pulled to +90, and
`-0.0` -- which formats as "-0.000" and reads like a real offset -- is
normalised away.

### Why it shipped

There was no `test/test_kspace.py`. There is now, with 57 tests, and the
ones that matter are properties rather than values: the answer does not
depend on the order of the two points, on which side of the origin they sit,
on where the pair is translated to, or on how far apart they are. The last
one goes end to end, turning the picked direction with the rotation matrix
the conversion actually applies and checking it lands on ky.

The Brillouin-zone dialog has a similar picker. It only *reports* an angle
rather than filling a box, and it is explicit that it measures
counter-clockwise from kx, so it keeps doing that -- but it is the number
someone then types into the zone azimuth, so it now shows the smallest
upright rotation beside it, from the same function, and says which of the
two depends on click order.

## Thirty-second round: k_z, the inner potential, and how the crystal cleaved

Ported from the lab's `plotKZ.m`, `kzconversionV4.m`, `kzconversion.m`,
`Eslice_conversion.m` and `kz_plot.m`. The formulas in those are standard and
correct; what is rebuilt here is how they are applied, and three things in
them that are not.

### The map is inverted, not projected

The MATLAB walks the measured grid forwards, works out where each sample
lands in `(k_z, k_par)`, and interpolates that scattered cloud onto a regular
grid. On a real cube that is **199 s** measured, which is why `kzconversionV4`
then approximates: it triangulates one energy slice, reuses the barycentric
weights for all the others, and corrects with a rigid shift along k_z. The
shift is right only if k_par does not move, and k_par goes as `sqrt(E_kin)`,
so every slice but the first comes out with the wrong k_par axis.

None of it is necessary, because the map inverts in closed form. Given a
target `(k_z, k_par)` and a binding energy there is exactly one photon energy
and one emission angle that produced it:

```
E_kin = [k_z² + k_par² cos²θp - A m* V0] / [A (m* - sin²θp)]
hv    = E_kin + W - E
sin α = k_par / sqrt(A E_kin)
```

So the conversion is a resampling of the original regular grid: one
`map_coordinates` call per energy slice. **0.003 s a slice against 0.32 s,
about 120×**, exact for every slice rather than for the first, and points the
measurement never reached are NaN by construction instead of whatever the
nearest triangle held. The round trip is exact to 1e-13 over 20,000 random
parameter sets, manipulator tilt included.

### Three errors in the original

**k_z is missing an in-plane component.** `Eslice_conversion.m` computes the
second in-plane momentum `kx1` that a tilted manipulator produces, comments
that it is there "to check the deviation", and then subtracts only `kpar`
inside the square root. At a manipulator angle of 20° that inflates k_z by
**0.257 Å⁻¹** — a quarter of a zone for a 6 Å repeat. `kz_plot.m:94` has the
same omission in the arc it draws.

**`kzconversionV4` ignores two of its own arguments.** `theta_offset` and
`theta_position` are in the signature and appear nowhere in the body. So the
GUI's two buttons run different physics from the same input boxes: one
subtracts an offset from the slit angle and knows nothing of the manipulator,
the other applies a 3-D rotation and does. They agree exactly when
`theta_position = 0` (the rotation reduces to `sin(α − φ)`) and diverge
silently the moment the manipulator leaves normal.

**One imaginary point abandons the whole cube.** The `isreal` check calls
`errordlg` and `return`s with the output variable never assigned, so the
caller gets an error rather than a result. Those points are physically
forbidden emission directions and belong as NaN. Worth knowing when it can
fire at all: with a free-electron final state (`m* = 1`) and `V0 > 0`, k_z is
**never** imaginary — `A E sin²α > A(E + V0)` cannot hold for `sin² ≤ 1`. It
takes `m* < 1` to reach that region, which is now a test.

### Choosing V₀, and what the scan can actually tell you

The inner potential is not measured, so the window is built around choosing
it: a preview that reconverts one slice in milliseconds, zone boundaries
drawn over it, and a scan of V₀ against the k_z period the data shows.

The honest part is the uncertainty. Two candidate criteria were measured
before either was implemented:

* **Maximising the periodicity amplitude** — peaks at the right V₀ but its
  half-width covers the whole 2–30 eV search range. Useless, and it was the
  first thing I tried.
* **Matching the measured period to 2π/d** — recovers 11.1 eV against a true
  12.0, with `d(period)/d(V0) = −0.0062 Å⁻¹/eV`. Measuring the period to 2%
  pins V₀ to **±3.4 eV**.

That ±3.4 eV is this scan's photon-energy range talking, not the fit: across
5–25 eV of V₀ the number of zones crossed only moves from 2.26 to 1.99. The
dialog reports the sensitivity and the uncertainty next to the answer rather
than presenting a single number as settled.

### Which plane the crystal cleaved along

A new tool: click two points that are the same feature one zone apart, and
their k_z separation is matched against every reciprocal lattice vector
within 15%. The matches are the candidate surface normals.

Two things make that less trivial than it sounds, and both are handled by
enumerating the **primitive** reciprocal lattice rather than the conventional
cell. Centring: in a body-centred lattice (001) is a forbidden reflection, so
the first one along that normal is (002) and the period is `4π/a`, not
`2π/a` — a factor of two, and a plausible-looking one. And which zone: a user
picking "the same feature again" may well be two zones apart, so each
candidate is tested at one, two and three orders and says which it matched at.

### Fermi-surface correction on a kz map

Checked rather than assumed, since a kz map's first axis is a photon energy
and correcting along *that* would be meaningless. It is already right: the
correction is offered from the slit cut, whose axes are (analyser angle,
energy), and applies the fit to the whole cube along dimension 1. The
deflector cut correctly offers nothing.

What the check did turn up is that the correction grows the energy axis and
pads the ends with NaN — 8% of the cube in the test. The conversion has to
leave those as holes rather than smear them, which bilinear sampling does,
and which is now a test.

The conversion reads the energy axis literally, so an edge still bending
across the analyser angle becomes a bend in k_z that looks like dispersion —
the artefact nobody catches afterwards. The dialog measures the bend on open
and says so. On the real LHhv scan it is 0.27 eV, so that reminder earns its
place. It is a reminder and not a refusal: a scan taken well away from E_F
has no edge to flatten.

### New files

| | |
|---|---|
| `tools/kzconv.py` | `forward`, `inverse`, `to_kz_cube`, `photon_arc`, `scan_inner_potential`, `edge_flatness`. No Qt. |
| `tools/cleavage.py` | `reciprocal_lengths`, `candidates` — a measured period to a list of planes. |
| `ui/kzconv.py` | The window: V0, m*, work function (read from the file, editable), geometry, the preview with zone lines, the V0 scan, the two-point period tool. |
| `test/test_kzconv.py` | 24 tests. |
| `test/test_cleavage.py` | 19 tests. |

Data kinds gained `kz_map_k` — a converted scan, `(k_z, k_par, E)` in Å⁻¹ —
and `MOMENTUM_KINDS`, the set whose first two axes are already momentum, so
that nothing offering an angle conversion has to name kinds one at a time.

## Thirty-first round: calibrating a kz map against its own Fermi edges

The loader stacks a photon-energy scan as measured and says plainly that the
absolute binding energies are only as good as the nominal photon energies and
a work-function table (thirtieth round). This is the step that stops relying
on either. Every spectrum has a Fermi edge in it, and where that edge sits
*is* the energy reference -- so fit it, shift each spectrum onto its own, and
the stack is calibrated by the measurement rather than by a table.

**kz map processing...** appears on a kz map's viewer and nowhere else.

### One box, used by index

The slit cut opens alongside with a selection box on it. The box says two
things at once: which angle channels are summed into the EDC, and what energy
range the fit runs over. It is then applied to every spectrum **by detector
index**.

By index deliberately. An index range is the same channels and the same
analyser window in every spectrum; an energy range would mean something
different in each, because these are precisely the spectra whose energy
scales do not yet agree. Picking in energy would make the fit depend on the
misalignment it exists to measure.

### Where the box starts, and why it is not the whole cut

It starts around wherever the edge appears to be, guessed from the whole scan
summed into one EDC. That was not the first version, which used the full
range, and the difference is larger than expected: on the real 41-spectrum
scan, fitting over all 629 energy points takes **37 s**; over the ±0.35 eV
the edge actually needs, **3.8 s**. The narrow window is also slightly more
accurate -- a window stretching far past the edge lets the tails pull the
background terms around.

### Fits that converge without finding anything

`fit_fermi_edge` reports success on data with no edge in it. It is a
least-squares fit, and a model with a step in it can always describe
featureless data by putting the step somewhere it does no harm. Two cases
turned up while testing, both returning numbers that look like answers:

* On flat data the step slides to the end of the fitted window and stops
  there, so E_F comes back as exactly the window bound.
* Seeded from a good neighbour -- which the fitting does, because
  neighbouring photon energies have nearly the same edge -- a dead spectrum
  converges obediently at the neighbour's level with a step height of
  essentially zero.

And on real data, a fit over the full 2.4 eV range ran away into the band
structure below the edge and returned E_F at 1.486 eV, outside the axis
entirely. That one appeared while checking the *result* of a successful run,
which is a good argument for verifying a calibration the same way it was
made.

So a fit is believed only if the level lands inside the window with at least
one energy step of room at each end, and the step height is positive and
larger than twice the error bar the fit returns on it. No hand-tuned
threshold, and nothing that assumes how bright a spectrum is -- which matters
when the flux across a scan varies ninefold. A spectrum that fails is marked
and filled in from its neighbours rather than dropped; one bad spectrum in
the middle of a scan should not take the cube with it.

### Cropping, and the order of the normalisation

Each spectrum shifts by its own amount, so afterwards they no longer cover
the same range: the cube is cropped to what all of them still reach, which
costs exactly the spread of the fitted levels. Keeping more would mean
inventing data at the edges, and a ragged top is invisible in a slice -- it
looks like an intensity drop.

The tick box normalises each spectrum by its own total intensity, on by
default. The photon flux and the analyser transmission both vary by a large
factor across a wide scan, so raw counts in a kz map are mostly a picture of
the beamline. It runs **after** aligning and cropping, so every total is a
sum over the same energy range; normalising first would divide each spectrum
by a sum over a slightly different window and bias the result by exactly the
misalignment being corrected.

### On the real scan

All 41 spectra fitted. E_F drifts smoothly from -0.14 to -0.26 eV and back,
largest step between neighbours 19 meV, total spread **0.107 eV**. After
alignment the edges sit at zero to **0.3 meV**.

That spread is worth noting against the thirtieth round: the nominal photon
energies and the work-function table put the members' windows 0.60 eV apart,
where their own Fermi edges say 0.11 eV. The loader's default -- stack as
measured, assume the operator set each window from a measured Fermi level --
was much closer to right than referring each member by its nominal hv would
have been, and this is the measurement that says so.

### New files

| | |
|---|---|
| `tools/kzmap.py` | `fit_levels`, `align`, `normalise_totals`, `process_kz_map`. No Qt, so a whole directory of scans can be calibrated from a script. |
| `ui/kzmap.py` | The window: the box, the options, and E_F plotted against photon energy as it goes. |
| `test/test_kzmap.py` | 22 tests, all built on scans whose Fermi level is known spectrum by spectrum -- the failure mode here is not a crash but a convincing cube aligned to the wrong place. |

## Thirtieth round: CASSIOPEE, folders as measurements, and a kz map

The second beamline. CASSIOPEE writes Scienta SES text, not NeXus, and it
writes a Fermi-surface map as a *folder of numbered spectra* rather than as
a cube -- so this round is as much about what a loader is allowed to decide
as about parsing.

Ported from the lab's `load_Soleil_Cassiopee_struct.m` and
`load_Soleil_Cassiopee_folder_struct.m`. Written first against one sample
spectrum, then corrected against two real folders (a 61-cut polar map and a
41-cut photon-energy scan), which is where most of what follows came from.

### Two writers, and only one of them has a `[Data]` section

The first version read the format documented by the sample file: SES 1.3.1,
which puts the counts in a `[Data n]` section. Both real folders are written
by SES **1.2.5**, which has no such section -- the counts simply follow the
`[Run Mode Information n]` block after a stray `inputA=` line and a blank
one. No `[Data n]`, no rows, no entries, nothing listed. It also uses lone
CR line endings, `key=`*tab*`value`, and -- the part that would have been
nastiest to find by eye -- a `Dimension 1 size` of **264.64** for a scale of
629 points. Not merely wrong: not a whole number.

`load_Soleil_Cassiopee_struct.m` handles the missing section by counting
three lines on from `[Run Mode Information` and trusting the offset. That
works for the files it was written against and breaks silently on the first
file with one more trailing field -- it would read counts as an axis, or an
axis as counts, and still produce a picture. The block is found here by what
it *is*: the run of lines that parse as rows of numbers. One rule, both
layouts, and no dependence on a field count. The point counts come from the
scales and never from the `size` fields.

### A false alarm that was really a printing artefact

The linearity check fired on one real file, `LHhv_34_ROI1_.txt`, claiming its
energy scale was off by 1.27 of a step. It was not. That file's window
crosses 100 eV, and the writer prints five *significant figures*: `99.998`
below the boundary, `100.01` above it. The rounding grid gets ten times
coarser part-way through one axis.

The check's allowance had been written in decimal places, which takes the
fine half of the axis for the whole of it and then reports the coarse half as
curvature. It now measures significant figures and takes the widest count in
the array -- widest, not narrowest, because a value landing on a round number
(`100.00`) looks like it was printed to fewer figures than it was, and
believing that would make the allowance enormous.

A check that fires on correct data is worse than no check, because nobody
reads the one that matters.

### A folder is one measurement

Selecting any file in a numbered folder -- or five of them, or all sixty-one
-- lists the folder **once**, as the assembled cube, with the individual
spectra below it. The series entry carries a path of its own (the folder's
first member) so that its identity is the measurement rather than whichever
file happened to be clicked; without that, selecting five files added five
identical copies of the same map.

Which quantity was stepped is read from the small `_i` parameter files, so a
folder can be characterised -- and labelled in the browser -- for about the
cost of one spectrum, rather than by reading sixty megabytes of text to find
out what to call it.

The stepped axis is straightened to an even grid. A real 61-point polar scan
asks for 0.5 deg a step and reads back 0.504, 0.504, 0.495, 0.504: that
jitter is the mechanism, not the sample, and an image places its rows evenly
whatever the axis says -- so an uneven axis is not honoured, merely
mislabelled. The readbacks are kept in `info`.

### `kz_map`

A photon-energy scan is now its own kind rather than a `map` with a label on
it. `CUBE_KINDS` is the set of three-axis cubes (`map`, `k_map`, `kz_map`)
that share a viewer, the 3-D view and the processing panels; a `kz_map`
differs in being refused a k conversion, since turning a photon energy into
k_z needs the inner potential, which is a property of the sample.

Adding it turned up two hand-written copies of "which axis is the energy"
(`{"cut": "y", "map": "z", ...}`) in two panels, both of which would have
refused to shift a `kz_map`'s energy axis with "has no energy axis to
shift". Energy is always the last constructor slot, for every kind there is,
so `energy_slot()` derives it and the copies are gone. Same for the two
copies of the kind-to-display-name map.

### What a kz map's energy axis is, and what it is not

The members of the real 41-point hv scan were each measured over a
*different* kinetic-energy window: across 40-120 eV the window drifts by
1.05 eV, a quarter of the 2.5 eV window and 250 energy channels. One axis
taken from the first member -- what the `.m` does -- puts member forty-one's
data an electronvolt from where its own photon energy says it belongs.

Referring each member by its own hv and phi(hv) and keeping only the range
they all share was tried, and is available as
`energy_reference="per_member"`. On the real scan it costs 1.2 eV of a
2.5 eV window, and the part it trims off the top is the Fermi level -- the
whole point of the measurement. Worse, it trims by exactly the amount the
calibration is wrong: the nominal photon energies and the tabulated work
function do not agree with each other here, so trimming to their verdict
compounds the error.

So the default stacks the members as measured, on the first spectrum's
referred axis, the way the lab's MATLAB already does. Nothing is
interpolated, nothing is trimmed, every count stays where the file put it.
The assumption is named rather than hidden -- that each window was set from
the Fermi level measured at that photon energy, which is how an hv scan is
normally taken -- and the per-member photon energy, work function and
energy-axis start go into `info`, which is what a proper per-photon-energy
Fermi-edge calibration will need. That calibration is a real step with its
own tool, not something to guess at during a load.

### Three places the port deliberately differs from the `.m`

* **The median filter is off.** The `.m` ran a two-point median filter along
  both axes of every spectrum, unconditionally, inside the load. On raw
  counts that is a processing step, and doing it silently means nobody
  downstream can tell it happened. It is a parameter here, off by default;
  the Process panel's despiking is the same operation, visible. Its edge is
  also fixed: `medfilt1` zero-pads, so the original's first row and column
  come out at half value -- a dark line down two sides of every spectrum,
  averaged into a map sixty times over.
* **A file name may contain underscores.** The `.m` splits names on `_` and
  takes the first piece as the base, which breaks for a folder called
  `Map80eV_Phi18` -- there is a commented-out hand-set `name_base` in it that
  says so. Membership is matched here on the `_<n>_ROI<n>_` tail.
* **The parameter file is split on the colon, not at a fixed offset.** The
  files are hand-formatted and the spacing varies; a field that moved by a
  character came back as a silent NaN. The colon has to be the first one
  *outside brackets*, because the polarisation field's own label lists its
  codes as `[0:LV, 1:LH, 2:AV, 3:AH, 4:CR]`. The first test written against
  that line found the bug. The real files then turned out to write the
  polarisation as the word `LH` rather than as a code, which is handled too.

### Reading a folder without freezing the window

A 61-cut folder is 240 MB of text and takes 7 s. `load_dataset` is called
from eight places that all expect a value back rather than a callback, so
`ui/jobs.py` gained `run_blocking`: the work runs on the worker thread as
usual, but the caller waits for the result while the event loop is pumped, so
the window keeps painting and the bar moves. What makes pumping the loop safe
here -- and not the `processEvents`-inside-a-computation pattern this module
was written to replace -- is that the progress dialog is **application**-modal
for the duration: the only thing that can be clicked is Cancel.

Two things had to be right. The runner reports itself busy during a blocking
job, or a timer firing in the pumped loop could start a second one and break
the one-job-at-a-time rule by way of the mechanism keeping the window alive.
And the wait is on the *result*, not on the thread -- a `QThread` runs an
event loop and stays "running" until `quit()`, which only comes after the
wait. Waiting on the thread deadlocked; the first run of the new tests hung
for two minutes and said so.

Which loads take that route is the reader's own declaration: a loader whose
`load` accepts a `progress` argument is one that expected to need it.

### Two things this turned up elsewhere

* The Load window's file filter was a hard-coded `*.nxs *.h5 *.hdf5`, so a
  CASSIOPEE spectrum -- a `.txt` -- could not be selected at all. It is now
  built from what the loaders say they read, with an entry per beamline.
* `load_soleil_nxs` refused any file not named `.nxs`, so a dataset this
  program saved as `.h5` could not be reopened. Latent, since the Save dialog
  forces `.nxs`, but the check was the wrong guard: detection is by a marker
  inside the file, and the extension test is there to catch a file of the
  wrong *type*. It now accepts the HDF5 family.

### The README split in two

The README had become 2,260 lines, almost all of it history. It is now a
guide to what the program does and how to use it; this file is the history.

### New files

| | |
|---|---|
| `loader/cassiopee.py` | The reader: `parse_scienta`, the parameter file, the work-function table, `series_members`, `series_kind`, `load_cut`, `load_series`. |
| `test/test_cassiopee.py` | 59 tests. The files are built in the test rather than shipped -- a real spectrum is a megabyte and what has to be checked is the format. |
| `CHANGELOG.md` | This file. |

## Twenty-ninth round: a name, a shape, and a paper trail

### It is called ARPES viewer now, and it is laid out like one program

`NxsLoader` named the file format. The program stopped being about one
format two rounds ago, when reading a beamline file became one loader among
several, so it is now `ARPES_viewer.py` and the thirty-five modules that
used to sit in one flat folder are in five:

```
loader/    reading data in, and keeping it      (no Qt)
tools/     the algorithms                       (no Qt)
ui/        the windows                          (Qt)
devtools/  one-off generators, run by hand
test/      the tests
```

The split is by **dependency**, not by taste. `tools/` and `loader/` import
no Qt and are usable from a plain script on a machine with no display --
which is how they are tested -- and `ui/` is the only package that does. A
module's folder now says what it is allowed to depend on, which a flat
folder of `nxs_*.py` could only say in a docstring.

Module names lost the `nxs_` prefix, which the folder now provides:
`nxs_kspace.py` is `tools/kspace.py`, `nxs_loaders.py` is
`loader/registry.py`, `ProcessNxsFile.py` is `loader/nxs_file.py`. Imports
are absolute from the project root (`from tools.kspace import convert_map`),
`conftest.py` puts that root on the path for pytest, and running
`python ARPES_viewer.py` from the folder works because Python puts the
script's own directory there.

Two things the move turned up, both worth naming because they are the sort
of thing a rename introduces silently. The launcher's own global `ui` --
the built main window, used on nearly every line of it -- shadowed the new
`ui` *package* the moment it was assigned, so `ui.main_window` started
looking for an attribute on a window; the launcher now imports those under
short names instead. And the marker attributes *inside* saved files
(`nxsloader_kind`, `nxsloader_format`) deliberately did **not** change:
renaming them would make every file anyone has already saved
unrecognisable. There is a test that says so.

`~/.nxsloader/` became `~/.arpes_viewer/`, and the recovery offer still
looks in the old location, because an upgrade must not be the way somebody
loses a crashed run's work.

### A bug the rename shipped, and the check that now prevents it

The Load-data window would not open:

```
for loader in loader.registry.loaders():
    UnboundLocalError: cannot access local variable 'loader'
```

`import loader.registry` binds the name **`loader`**, so every use reads
`loader.registry...` -- and `loader` is also the obvious name for a variable
holding one. Assigning it anywhere in a function makes it local for the
whole function, so the lookup on the right-hand side fails. Before the
rename the module was called `nxs_loaders` and nothing would ever collide
with it; the new short package names made the collision available, and the
loop variable took it.

Two fixes, because one of them is only the symptom. The tree no longer
refers to a package by its bare name anywhere: `from tools import kspace`
binds `kspace`, not `tools`, so a collision is a plain `NameError` at import
rather than an `UnboundLocalError` down whichever branch happens to run.
That removed about a hundred latent instances of the same trap across
seventeen files.

And `test/test_imports.py` checks the shape of the code directly, because
nothing else would have caught this: no import fails, and no test of
behaviour notices unless it happens to open that one menu. It walks the
ASTs and asserts that no local shadows an imported module, that no package
is imported by its bare name, that every module imports cleanly, and that
`loader/` and `tools/` import no Qt.

That last one immediately found a placement mistake of my own:
`tools/system.py` had the two Qt message boxes in it, in the package whose
whole claim is that it needs no display. They moved to `ui/notify.py`.

### The session folder keeps a log, and tidies up after itself

Three changes to how working copies are kept:

* **Deleted when they are no longer working copies.** Saving a dataset to a
  file of your own now re-points its row at *that* file -- the row's
  identity becomes `(your file, the entry in it)` -- and only then removes
  the session copy. The folder used to accumulate a duplicate of everything
  already saved, and the duplicate, not the saved file, stayed the thing the
  program would re-read.
* **Three days, not seven.** A copy that matters is saved; one still sitting
  there after three days is one nobody came back for.
* **An operations log.** Every write, save, deletion, prune and recovery is
  appended to `~/.arpes_viewer/sessions/operations.log` with its path and
  size -- one plain-text file, in the order things happened, readable in any
  editor on a beamline machine with nothing installed. It is in the list's
  right-click menu as "Session log...", and it is what to open when the
  question is "where did my work go".

And the program now **offers to reopen what a crash left behind**. A session
that ended properly leaves nothing -- saving removes a copy, so does
removing its row -- so anything still there is unsaved work from a run that
stopped, and it is listed at startup under its own name rather than its
filename (the native loader reports the name the dataset was saved under,
which is why a recovered row says "gold Fermi surface" and not
"007-gold Fermi surface.nxs").

### The processing panels moved off the GUI thread

The 2-D and 3-D processing panels now run through `ui/jobs.py` like the k
conversion does, with a real progress bar and a Cancel that takes effect
between planes. Their arrays are read first, on the GUI thread, because
h5py here is not built for two threads; and every widget they consult is
read there too -- reading a spin box from a worker mostly appears to work,
which is the worst kind of unsafe, and a value that changed halfway through
a run would process half the datasets differently.

Two loops stay on the GUI thread on purpose: the Fermi-edge channel fit and
the MDC/EDC line fit hand their answer straight back to the line that called
them, and turning that into a callback would restructure a working
measurement workflow for a fit that is already interruptible. What they no
longer do is leave the door open: their progress dialogs are now
**application** modal, so pumping events cannot dispatch a click on another
window into a second copy of the same code -- which was the actual hazard in
that pattern, rather than the pumping.

### One table for what axes a dataset has

"Which axis slots does this kind use" was answered in five places:
`tools/dataops.py`'s two tables, the loader's `_axis_slots`,
`session.scan_to_dict`, the launcher's `_dataset_for_saving` and
`_MemScan.__init__`. Five copies is four chances for a new kind to be
half-supported -- and one of them was already wrong: the launcher wrote a
4-D spatial scan with three of its four axes, dropping energy, so a saved
spatial scan came back with the wrong axes.

There is now one table, `loader.nxs_file.AXIS_SLOTS`, with `scan.axes()`
and `scan.array()` on `NxsScan` for callers that want to be written once for
every kind. The other five read it. This is the piece of the "make the data
model general" work that was worth doing on its own: it removes the actual
hazard (five tables disagreeing) without rewriting the 59 places that choose
a *window* by kind, which is a display decision and fine as it is.


## Twenty-eighth round: nothing is lost, and loading is a place

Four things were true of every earlier round, and all four are now false:

* every computation ran on the GUI thread, so a long one looked exactly
  like a crash;
* a dataset computed in the session existed **only** as a numpy array, so
  killing the frozen-looking program lost the morning's work;
* nothing could be unloaded, so twenty derived cubes meant twenty resident
  cubes and no ceiling;
* there was one reader, named after the beamline it reads.

### Loading has its own window

**"Load data..."** opens the loader (`ui/loader_dialog.py`): choose files,
see which reader recognised them, override that if the guess is wrong, look
at the datasets found inside, and set two things the file cannot say for
itself. Only structure is read there -- entry names, kinds, shapes -- so
pointing it at a folder of multi-GB scans is still instant, and a file that
is listed but never opened is never read.

The two axis questions are asked at the door because correcting them later
means correcting them in five views:

* **Axis order.** Two beamlines record the same measurement with the
  array's dimensions in a different order. Reorder them once, here, and the
  axis vectors and their labels come along (which is the part that is easy
  to get wrong by hand).
* **What the scanned axis actually is.** A map is angle vs angle vs energy
  by default -- but the same acquisition records a photon-energy,
  temperature or gate-voltage series. Same shape, different physics. Saying
  so keeps the wrong unit out of every later step, and makes the k
  conversion *refuse* an axis that is not an emission angle rather than
  produce a plausible, meaningless picture.

### One reader per beamline

`loader/registry.py` is a registry; each reader is its own file. Adding a
beamline is a new `nxs_loader_<name>.py` with four methods and one
`register()` call, and nothing that already worked is edited --
`loader/soleil.py` is deliberately thin and is the file to copy.

`loader/native.py` reads the program's own format and outranks every
beamline guess, because that format marks its own entries: reopening
yesterday's work never asks which synchrotron it came from.

### Derived datasets are written the moment they exist

`loader/session.py` gives each run a session folder, and every computed
dataset -- k-map, arbitrary cut, processed result -- is written to it
immediately. Two things follow. The obvious one is that stopping the
program stops losing work. The less obvious one is that a listed dataset is
now a *path*, exactly as a measured one always was, so the in-memory copy
became a cache: `MemoryBudget` keeps a bounded amount resident (1.5 GB by
default) and drops the least recently used, which is re-read from its file
when it is next opened. A dataset open in a window is held by that window
too, so eviction never pulls the rug from under anything being looked at.

Spatial scans (`spem_4d`, `spem_1d`) are **not** auto-saved: they are the
measurement itself, gigabytes of it, and a copy after every operation would
cost more than it protects.

Session folders older than a week are pruned at startup. Closing the main
panel with datasets that have not been saved anywhere of your own now asks,
with the choice of saving them, closing anyway, or cancelling -- it used to
close on the spot.

### Long operations run off the GUI thread

`ui/jobs.py` runs one job at a time on a worker thread with a progress bar
that appears only if the job lasts longer than 400 ms, and a Cancel that
works (cooperatively: it takes effect at the job's next progress report).
Loading several files and the k conversion go through it. One at a time is
not a simplification -- HDF5 here is not built thread-safe, so two readers
at once is a crash rather than a slowdown -- and the cube is read before
the worker starts for the same reason.

This replaces the four `QApplication.processEvents()` calls that used to
paper over the freeze. That trick keeps a window painting but re-enters the
event loop from inside the computation, so clicking another button ran that
handler on top of the half-finished one.

### Maps and cuts stay on disk

`LazyCube` (which the 4-D spatial scan already used) now also backs maps and
cuts, and `LazyArray` does the same for saved datasets of any rank. Opening
a map no longer reads it: the contour asks for one energy slice at a time,
and a map that is browsed and not used costs nothing. Both grew
`__array__`, so anything that genuinely needs every point -- an operation,
a fit, a save -- still works by saying `np.asarray`, and `full_cube()` on
the contour window is where that is said out loud.

Two bugs surfaced doing this, both pre-existing and both invisible until
maps were read this way. `LazyCube._normalise` tested `Ellipsis in key`,
which compares element-wise against an array index and raises "the truth
value of an array ... is ambiguous" -- and every one of the viewer's slab
reads passes an array index. And the handle-retention test in
`load_soleil_nxs` only knew about `LazyCube`, so a lazily-read value array
had its file closed underneath it.

### The saved format

Same HDF5 container, so old files still open, with two changes measured
rather than assumed (on a 120x941x96 map of smooth bands plus Poisson
noise):

| | file | write |
|---|---|---|
| raw | 86.7 MB | 0.04 s |
| gzip 4 (the old format) | 11.1 MB | 0.70 s |
| gzip 1 + shuffle | 9.5 MB | 0.47 s |
| **gzip 4 + shuffle (now)** | **8.8 MB** | **0.70 s** |
| gzip 9 + shuffle | 8.4 MB | 9.87 s |
| lzf + shuffle | 13.5 MB | 0.33 s |

The shuffle filter re-orders the bytes of each value so that the high bytes
of neighbouring numbers -- nearly always equal -- end up adjacent, which is
what gzip is good at. It is 21% smaller at the same write time, and
lossless. Level 9 costs 14x the write time for 5% and was rejected; lzf
writes faster but is half as effective.

The array is also stored in **its own dtype** now. The old writer took
whatever the caller had promoted to float64, which doubled the file for data
that arrived as float32 or as integer counts and gained nothing.

Chunking is left to h5py on purpose: a hand-tuned "whole plane" chunk reads
constant-energy slices 4x faster but the orthogonal cuts 20x slower, and the
viewer does both.

Neither change needs a new dependency -- shuffle and gzip ship with h5py.

### New files

`loader/session.py` (the session folder and the memory budget), `ui/jobs.py`
(the worker thread and its progress dialog), `loader/registry.py` (the
registry and the load-time axis options), `loader/soleil.py` and
`loader/native.py` (the two readers), and `ui/loader_dialog.py` (the
loader window). All but the last are Qt-free except `ui.jobs`, which is
Qt by nature, and all are covered by `test/test_session.py`, `test/test_jobs.py`,
`test/test_loaders.py` and `test/test_lazy.py`. Run them all with `python -m pytest` from the project root.


## Twenty-seventh round: Brillouin-zone overlay

A **"Brillouin zone..."** button appears on a map's constant-energy contour
once (and only once) it has been converted to k-space -- a Brillouin zone is
a statement about momentum-space periodicity, so it means nothing on the raw
angle-space contour, the mirror image of "Map k conversion" hiding itself
once a map already is a k-map. It opens a modeless dialog
(`BrillouinZoneDialog` in `ui/windows.py`), which works the way the
other contour tools do: the contour behind it stays live, so a direction can
be picked off it.

**Nothing is drawn until Plot is pressed.** The settings here describe a
crystal, and a half-typed lattice constant is not one, so redrawing on every
keystroke would spend its time showing zones nobody asked for. What does
follow the typing is the reference text -- which Bravais lattice the space
group implies, how far along its normal the cut plane has been slid -- and
the 3-D preview window, if it is open.

### The crystal, and what it is allowed to be

The **space-group number comes first**, and it then decides which of the six
cell parameters are the user's to choose: a cubic group offers one length
and no angles, a hexagonal group two lengths and no angles, P1 all six. The
rest are greyed out and follow along, so a cell that contradicts its own
space group cannot be typed in the first place. (The previous version
accepted anything and printed a warning underneath, which left the user
having to know the constraint anyway.) `nxs_lattice.free_parameters` is
where that table lives, and `test/test_lattice.py` checks that a cell built to it
is one `validate_lattice_parameters` has nothing to say about -- the two
halves cannot drift apart.

A rhombohedral (R) group is offered in **hexagonal axes only**. The
rhombohedral-axes setting (a = b = c, alpha = beta = gamma) describes the
same lattice equally well, but supporting both would mean one set of six
numbers meaning two different cells depending on a mode nobody would
remember setting; the label says which setting is in force.

### Conventional and irreducible zones

Either mode draws the **conventional** zone -- the first Brillouin zone,
i.e. the Wigner-Seitz cell of the reciprocal lattice -- or the
**irreducible** one, the wedge of it that the crystal's own symmetry repeats
into the whole zone.

The wedge is *constructed*, not looked up: pick a point that no symmetry
operation leaves fixed, and keep the part of the zone at least as close to
it as to every one of its images (the Dirichlet cell of that point's orbit,
a standard fundamental domain). The usual alternative -- a hard-coded list
of high-symmetry points per Bravais lattice -- is a table that is silently
wrong for any setting it did not anticipate, while this construction knows
nothing about crystal systems and is checkable: the wedge's volume is the
zone's divided by the number of operations, exactly, which is what
`test/test_bz3d.py` and `test/test_bz2d.py` assert for the cubic, hexagonal,
tetragonal, monoclinic and triclinic cases. For FCC it reproduces the
textbook six-corner Gamma-X-W-K-U-L solid, corner by corner.

In 3-D the symmetry used is the space group's point group **plus inversion**
-- the Laue class -- because time reversal makes `E(k) = E(-k)` even in a
crystal that is not itself centrosymmetric. The operations come from
`tools/spacegroups.py`, which now carries each space group's rotation matrices
as well as its Bravais lattice (still generated offline from spglib, still
no runtime dependency: 230 space groups need only 37 distinct matrix sets).

In 2-D there is no space group to ask, so the **plane lattice's own
holohedry** is used -- found by looking for the maps that preserve the Gram
matrix, so a centred-rectangular cell typed as "a = b with an odd angle" is
recognised rather than called oblique. The dialog says so, because a layer
whose basis is less symmetric than its lattice has a larger wedge than the
one drawn.

### The cut plane, and seeing it

The cut plane's normal is given as Miller indices `(h k l)`, read against
the **conventional** cell's reciprocal vectors -- `(0 0 1)` on an FCC
crystal is the cubic face everyone calls (001). (Read against the primitive
vectors, which are what the zone itself is built from, those same three
numbers would point along a <111>.) The zone is still the primitive cell's,
and still tiled by primitive translations; only the *indices* are
conventional.

**Offset along normal** slides the plane off Gamma, in A^-1, with a readout
saying what fraction of the distance from Gamma to the zone edge along that
normal it represents -- so "0.573 of the 1.571 A^-1 to the edge" rather than
a bare number. This is a geometric offset: turning a photon energy into the
`kz` that belongs here needs an inner potential, which is still the separate
kz feature's job.

**"Show cut plane in 3D..."** opens a turnable view of the zone with the
plane in it: the zone as a wireframe (faces towards the camera solid, away
dashed), the plane as a dashed square, the polygon where it actually cuts --
which is the polygon that goes onto the contour, from the same call, so the
two cannot disagree -- the normal as an arrow from Gamma, and the
conventional reciprocal axes for orientation. Drag to turn it; it re-frames
itself about 55 degrees off the normal whenever the normal changes, since a
viewpoint chosen for one plane is often exactly the wrong one for the next.

It is drawn by projecting the face loops by hand rather than through
`pyqtgraph.opengl`, for the same reason the volume viewer is: an
orthographic camera is a dozen lines, and adding PyOpenGL to a program that
installs with pip and no compiler is a poor trade for one preview window.

### Moire, and what is left out

The **moire** mode lives in the 2-D layer mode only. A moire pattern is two
stacked *layers*, so offering it for a 3-D crystal would have meant guessing
which 2-D sublattice the chosen cut plane exposes -- a different calculation
that nobody asked for, and the earlier version's one real approximation.
The moire zone comes from the two layers' nearby reciprocal-lattice
differences (`nxs_moire.moire_reciprocal_vectors`), and the closed-form
hexagonal formula ported from `moire_lattice.m` is shown beside it as a
cross-check whenever both layers are hexagonal.

Only the first zone is ever drawn, repeated by translation across the field
of view if asked -- never a second- or higher-order zone. Settings persist
into the k-map dataset's own `scan.info` as flat `bz.*` keys (matching the
existing `kconv.*`/`arbcut.*` convention), so reopening the tool on a
saved-and-reloaded k-map restores them.

### New files

`tools/spacegroups.py` (auto-generated by `devtools/gen_spacegroups.py`, which needs
`spglib` -- the app itself does not) is the space-group -> Bravais-lattice
and point-group lookup; `tools/lattice.py` builds real/reciprocal lattice
vectors, the cell constraints and the Cartesian point-group operations from
cell parameters (fixing a bug found in the MATLAB reference's
`gen_realspace_basevector.m`: no `alpha` input at all, `cos(beta)` reused
where the general triclinic formula needs `cos(alpha)` -- invisible for the
crystal systems where both happen to be zero, wrong for monoclinic,
triclinic and trigonal/rhombohedral cells). `tools/bz3d.py` and `tools/bz2d.py`
are the 3-D and 2-D zone geometry: Wigner-Seitz cell, irreducible wedge,
plane cut and tiling. The 3-D one grows its neighbour shell until the volume
converges, rather than the fixed +-1 shell `gen_brillouin.m` always uses,
which is not enough for a sufficiently oblique cell. `tools/moire.py` holds
both the general reciprocal-vector-difference method and the ported
closed-form fast path. All five are Qt-free and covered by
`test/test_lattice.py`, `test/test_bz3d.py`, `test/test_bz2d.py` and `test/test_moire.py`
against closed-form crystallography: known cell shapes and volumes, wedge
volumes equal to the zone over the group order, textbook high-symmetry
distances, the two-comb beat-frequency limit, and the hexagonal point
group's 60-degree twist folding. Run them all with `python -m pytest` from the project root.


## Twenty-sixth round: MDC/EDC fitting, the band, and what it says

*(`ui/fit.py`, `tools/peaks.py`, `tools/dispersion.py`)*

"MDC / EDC fit" in the cut viewer, or on the launcher's right-click menu.
Fit every line of a cut, then read the band off the fitted centres.

**Only k-converted cuts.** The button is there but disabled on a cut still
in degrees, and says why. Nothing in the fitting would object to an angle
axis -- the peaks would come out perfectly well -- but every quantity built
on the band would then be wrong: `dE/dk` would be in eV·degree, `m*` would
be meaningless, and neither number carries its units along to say so. Run
`Cut k conversion` first, or slice a converted k-map.

The work goes: add a band, click its top on two or three curves down the
cut, press **Fit**. Seeds are interpolated between with PCHIP, so a band
seeded at three energies has a starting centre at every energy in between,
and a monotone interpolant cannot overshoot the way a spline does. The
nonlinear fit carries only centres and widths; heights and the background
are solved exactly, by bounded linear least squares, at every iteration
(variable projection). The bands' fitted positions come back onto the
picture, sized by amplitude.

| | |
|---|---|
| Shapes | Lorentzian, Gaussian, Voigt -- widths are FWHM everywhere, so the number reported is the number a paper quotes |
| Background | none, constant, linear, Shirley (the same curve the 2-D background tab subtracts) |
| Constraints | share a width between bands (the two branches of a symmetric dispersion), fix a width, bound every centre to ±1.5 widths of its seed |
| Resolution | deconvolved as a Gaussian, giving the intrinsic width rather than the measured one |
| **EDC Fermi factor** | `I(E) = [A(E)·f(E,T)] ⊗ R(E) + background` -- the occupation multiplies **before** the resolution smearing, and the background is not multiplied by it at all |

That last row is the one that matters and it is worth stating why. An EDC
peak within a resolution width of E<sub>F</sub> is cut into by the
Fermi function, and a fit that ignores it moves the centre to compensate. On
the test curves in `test_peaks.py`, with the factor every centre comes back
right to well under 1 meV with χ² ≈ 1.1; without it a peak 10 meV *above*
E<sub>F</sub> comes back 24 meV low with χ² = 38. E<sub>F</sub>, the
temperature and the resolution are read off the dataset when the Fermi round
wrote them, so a cut that has already been offset needs nothing typed.

**Dispersion...** fits E(k) through the fitted centres, weighted properly:
an MDC series has its error on k, so the weight is
σ<sub>eff</sub> = |dE/dk|·σ<sub>k</sub>, which needs the slope being fitted
and is therefore iterated. It reports `v_F` (and in m/s) or `m*`, with the
band bottom for a parabola.

Beside it is the **window scan**, which is the honest answer to "what is the
Fermi velocity". A velocity depends on how far from E<sub>F</sub> it is
measured, so the panel fits over a series of windows and plots the answer
against the window. A flat run means a real velocity; no flat run means
there is no window-independent number to quote, which is a result about the
band rather than a failure. The plateau is anchored at the **narrowest**
window and leaves it only when the deviation exceeds both a tolerance and
two of its own error bars -- a "longest flat run" rule instead picks the
wide-window asymptote, which is the average over the kink and belongs to no
part of the band.

**Self-energy...** takes Re Σ from the difference between the band and a
bare band (fitted, or pinned through anchor energies you choose) and Im Σ
from the MDC widths, and shows the Kramers-Kronig transform of Im Σ beside
Re Σ -- the check that says whether the bare band was a defensible choice.

### Fixed in this round: how pictures are drawn (`test_orientation.py`)

Three separate faults, all of them invisible to every existing test,
because they were in the **drawing** rather than in the arithmetic. The
fitting, the smoothing, the conversion and the export were all producing
correct numbers; only what reached the screen was wrong.

1. **The axes were exchanged** in the fit panel, the 2-D processing
   previews and the comparison's third panel.
   `pg.setConfigOptions(imageAxisOrder="row-major")` is global state and
   was being applied from inside a widget's constructor, so until the
   first viewer existed images were drawn with pyqtgraph's opposite
   default -- which of the two an image got depended on the order the user
   happened to open windows in. It is now set once, at import.
2. **Upside down.** `pg.ImageView` inverts its y axis, because its default
   subject is a photograph whose first row is the top one. The three
   panels above each built a raw `pg.ImageView()` and none of them undid
   it, so E<sub>F</sub> sat at the bottom.
3. **Squashed.** `pg.ImageView` also locks the aspect ratio, which
   flattens any plot whose axes carry unrelated units (eV against Å⁻¹).

The cause of all three was the same: four places each constructed their
own `pg.ImageView` and each got a different subset of the setup right.
There is now one `plain_image_view()` and one `show_frame()` in
`ui.widgets`, and every panel goes through them, so the convention
is stated once. The arrays themselves are untouched -- they stay indexed
`[ik, iE]`, and the transpose happens on the way to the screen and nowhere
else. Transposing the stored array instead would have moved every cursor,
seed and fitted centre onto the wrong axis while the picture looked right.

The processing previews and the comparison panel also gained axis labels
in the process; they had none before, which is part of why a transposed
picture was easy to miss.

`test_orientation.py` puts one bright pixel, off-centre on both axes, on a
non-square grid through all seven drawing paths and checks where it lands.
It reads the **rendered** image rather than the numpy array: pyqtgraph's
`imageAxisOrder` decides both how an array is read and how its indices map
to the picture, so a probe that indexes the array flips with it and always
agrees with itself. The first version of the test did exactly that and
passed with the global setting removed entirely.

Separately, the fit panel's opening view range is now set explicitly
rather than by `autoRange`: its cursor is an infinite line, which an
automatic range has to treat as an item without bounds, and a view that
has not been laid out yet has no size to compute a padding against.


---

## Twenty-fifth round: data processing (2-D and 3-D)

Two new panels, reached from the launcher's **Process...** button or from
the list's right-click menu, plus the displays that go with them. The split
is by dimensionality, as asked: a cut goes to the 2-D panel, a map or
k-map to the 3-D one.

The rule everything here follows: **an operation never modifies a dataset,
it produces a new one in the list**, carrying the whole chain of steps that
made it. `Show information` displays that chain, and it is written into the
`.nxs` when the dataset is saved, so a figure can still say how it was made
six months later. Nothing is computed until "Apply" -- the preview runs on
a decimated copy so a slider stays interactive on a large cut.

### The 2-D panel (`ui.process.ProcessDialog`)

Source and result side by side, always, because curvature, second
derivatives and background subtraction all produce plausible-looking
pictures whatever the settings, and the only way to see that a setting has
eaten the band is to have the original beside it.

| Tab | What it does |
| --- | --- |
| Smooth | Gaussian (width in eV / A^-1), Savitzky-Golay (keeps peak height -- use before a fit), moving average; optional spike removal first |
| Derivative | `-d2/dE2`, `-d2/dk2`, both normalised and combined, first derivative |
| Curvature | 2-D isotropic, 1-D along either axis, and the `I/|grad I|` enhancement |
| Background | Shirley, Tougaard, angle-independent percentile, 2-D polynomial, and dividing out the Fermi cut-off |
| Symmetrise | mirror about one line, placed by clicking the image |

Above the tabs, **Shape** controls how the previews are drawn: *Square
pixels* (the default, one array element square -- the only shape that shows
the data undistorted), *Equal axes* (one data unit the same length on both
axes, right for a k-vs-k map and meaningless for k against energy), or
*Fill the box*. The two previews are panned and zoomed together, so they
can never show different parts of the data, and the cursor's position is
read out in the data's own units beside the control.

Those coordinates are what make **Pick on the image** work: arm it and
click the Source preview to put the mirror line (or, in the 3-D panel, the
rotation centre) through that point, which is then drawn on both previews
as a cross and a dashed line. This replaces an automatic search for the
symmetry centre that was written first and then removed -- whenever it was
slightly wrong it produced a beautifully symmetric picture, and nothing
about the result showed that it had happened. Clicking is as quick and is
visible. The status line reports how far the two halves actually differed,
which is the number that says whether to believe the result.

Four things are deliberately **not** here. Normalising each line is already
`Data operations -> Self-normalise`, and having it in two places invites
the two drifting apart. Rotational symmetrisation belongs to a
constant-energy map, which is a plane of a cube, so it lives in the 3-D
panel -- on a cut of momentum against energy a rotation is not a symmetry
of anything, and offering the order, the inversion and a sector wedge only
invited a meaningless choice. Shifting an axis was not worth a tab.
Peak fitting has its own panel (the round above), because fitting is a
different kind of work from filtering: it is iterative, it is judged by a
residual, and its output is numbers rather than a picture.

### The 3-D panel and viewer (`ui.volume`)

The processing panel applies any of the 2-D operations plane by plane -- you
choose which planes, because the curvature of a constant-energy map and the
curvature of a dispersion are different quantities that happen to share a
formula -- and symmetrises a whole cube.

The viewer has four modes: **orthogonal slices**, the **notched cube**,
**volume projections** (maximum intensity, sum, alpha composite) and a
**shaded isosurface**. Drag the picture to turn it.

Before it opens, a dialog asks for the momentum and energy window and the
sampling, in points per axis, and shows the resulting array size as those
are changed. This is not optional politeness: the views resample the whole
cube for every frame, so on a measured cube most of the wait is spent on
the part of the detector nobody is looking at. Binning is a NaN-aware
*mean*, not a sum, so trading resolution for speed does not move the
colour scale.

Plane positions are typed as well as dragged -- a slice at exactly E_F or
k = 0 is the one anyone actually wants -- and the sliders report only when
a drag **ends**, so dragging through two hundred positions asks for one
rebuild rather than two hundred.

It is drawn with a plain `QPainter`, not OpenGL. Every surface in these views
is planar, and the parallel projection of a rectangle is a parallelogram --
exactly an affine image transform -- so a painter's algorithm over
depth-sorted faces is pixel-exact, needs no graphics driver or PyOpenGL, and
paints at any size. That last part is the point: "Save the view..." writes
PNG, PDF or SVG at whatever resolution the journal wants, where
`volume_3d_plot_w_notch.m`'s output could only ever be as big as its window.

### Six things these do that the MATLAB originals did not

1. **Curvature is dimensionless.** `Curvature.m` feeds `I_x` (counts per
   A^-1) and `I_y` (counts per eV) into the same `(a0 + I_x^2 + I_y^2)`.
   Those differ by one to two orders of magnitude, so its "isotropic"
   curvature is nothing of the sort, and the best `a0` moves with every
   change of axis range or intensity scale. Here both axes and the intensity
   are mapped to [0,1] first, so `a0` is a pure number: the test checks that
   the result is *identical* whether the energy axis is in eV or meV, and
   whether the intensity is multiplied by 1000. `Suggest` sets `a0` to the
   scale where the denominator starts to matter.
2. **No histogram clipping.** `Curvature.m` ends by multiplying its result
   by a mask built from a 100-bin histogram, which silently sets every
   negative curvature to zero -- half the structure, and the half that says
   where the band is *not*. Clipping is a display decision, so it now lives
   in the level bar (see below).
3. **Derivatives are Savitzky-Golay**, so smoothing and differentiating
   happen in one centred step and the border is *fitted* rather than forced
   to zero. `smooth_derivation_v1.m`'s zero frame is a feature that is not
   in the data and it drags the auto-levels with it. (Its three smoothing
   buttons are also hard-coded, and the one labelled "smoother" passes
   `p = 1` to `csaps`, which is no smoothing at all.)
4. **Errors are carried through.** Every fitted parameter comes back with an
   uncertainty and the band fit uses those as weights instead of treating
   peak positions as exact; on the test data the reduced chi-squared comes
   out at 1.02, which is what says the weighting is right. The fit also
   reads the **array**, not a plotted line's `XData`/`YData` as
   `EDC_MDC_plot.m` does -- that version fits whatever offset and spline the
   display happened to apply. (The fitting itself now lives in `nxs_peaks`
   and `nxs_dispersion`; see the round on the MDC/EDC panel.)
5. **Symmetrisation says how much it invented.** The 2-D mirror reports how
   far the two halves actually differed before it averaged them; the 3-D
   symmetrisation carries a coverage map (how many copies contributed to
   each pixel), because a region built from one copy is a measurement and a
   region built from six is an average. Mirrors and inversion are available,
   not just rotations, and in 3-D the missing-data mask is taken **per
   plane** -- `data_symmetrization_pro.m` takes it from the first z-plane
   only, with a comment admitting the problem, and for a converted k-map or
   an hv scan the valid region differs in every plane. Automatic
   centre-finding is offered as a button that *fills in* the boxes; the
   typed value is still what gets used, because a slightly wrong automatic
   origin produces a beautifully symmetric picture of nothing.
6. **The notch is one construction, not eight.** `volume_3d_plot_w_notch.m`
   writes its eight orientations out longhand, about 500 lines of
   near-identical code, and computes face sizes as `round((xe-xb)/xstep)`,
   which goes to zero or negative when the notch reaches an edge and then
   fails inside `slice`. Here the orientation is three sign flips over one
   construction, the notch vertex can sit anywhere, and a face squeezed to
   nothing is simply not built.

   The faces are cut up further than the shape needs -- an L into three
   rectangles, a slice plane into quadrants -- because a rectangle is what
   can be texture-mapped and disjoint pieces are what make depth sorting
   exact. Those cuts must not be drawn: outlining each rectangle puts a
   line across the top of the cube and another down its side that are
   nowhere in the geometry. An edge shared by two faces in the same plane
   is a cut, an edge only one face has is a boundary, and that single test
   leaves the 18 visible edges of the 21 a notched box really has.

### Displays

- **Stack plot**: a waterfall of EDCs or MDCs. The offset can follow each
  curve's **real position** rather than its index, the colour runs along a
  colormap keyed to that position (readable past the eight-colour cycle
  where `stack_plot.m` starts repeating), each curve can be normalised and
  trimmed, and the whole stack can be sent to the figure composer.

  Two ways in. From the launcher, right-click a dataset -> "Stack plot".
  Or, from inside a viewer, drag the selection box over the feature and
  right-click -> "Stack plot of the selection": the box already says which
  energies and which momenta, so the window opens with the range filled in
  and the direction worked out from the box's shape -- a tall narrow box
  means EDCs, a wide flat one means MDCs. The curves are still taken from
  the full data with the box only setting the range, because a curve
  cropped before it is normalised has the wrong area.

  Curve colours are darkened when the colormap would otherwise make them
  too pale to read: `gray` and `hot` end at pure white, so on the usual
  default the topmost curve was white on a white background and simply
  absent from the plot.

  Switching between EDCs and MDCs swaps the two range pairs with it. They
  describe different axes in the two directions, and left alone they keep
  the other axis's numbers -- a stack of EDCs taken over k = 0.1 .. 0.45,
  switched to MDCs, asked for energies of 0.1 and 0.45 eV, both past the
  top of the energy axis, so both clamped to the same index and the whole
  stack collapsed to one curve. A range holding fewer points than the
  number of curves asked for now says so rather than quietly returning
  fewer lines.
- **Curve fitting** now has its own panel (see the round below): one line at
  a time with its peaks separated, its background and **its residual
  below**. The residual is the point -- a fit of two Lorentzians to a
  three-peak MDC looks convincing in the overlay and obviously wrong in the
  residual, and the MATLAB tool drew only the overlay.
- **Compare two datasets**: side by side on shared levels, with their
  difference, ratio or asymmetry on a diverging colormap. A difference
  between two different grids is refused rather than interpolated, since
  that would be a picture of the interpolation.
- **Line profile** (right-click any image): a draggable line at any angle
  with a settable integration width. EDC and MDC are its two axis-aligned
  special cases; a band running at an angle has no axis-aligned cut that
  follows it, and the usual workaround is reading positions off the picture
  by eye.
- **Level histogram**: the intensity distribution is now drawn inside the
  level bar, on a square-root scale (on a linear one an ARPES image is one
  spike at the background level and nothing else). The **Clip** button snaps
  the handles to the 1st/99th percentile -- right-click it for other pairs --
  which stops a handful of hot pixels holding the whole top of the scale.

### New files

`tools/process.py`, `tools/volume.py`, `tools/peaks.py` and `tools/dispersion.py`
are Qt-free, so every algorithm is usable from a script and is tested
directly. `ui/process.py`, `ui/volume.py` and
`ui/fit.py` hold the panels. `test_process.py`, `test_volume.py` and
`test_peaks.py` check them against data whose answer is known independently:
a derivative of an analytic function, a symmetric Fermi surface whose centre
was chosen when it was built, a Lorentzian series whose velocity and width
were put in by hand.

### Not included, on purpose

The Brillouin-zone tool (the round at the top of this file) draws on the
live k-map contour only; it does not yet hand its polygons to the figure
composer's `nxs_figure.Overlay` list (the "curves in data coordinates" hook
described under **Overlays** above), so a zone drawn on screen is not yet
exported into a saved figure -- that wiring is left for later.
High-symmetry path stitching is already covered by arbitrary cuts. `kz` /
inner-potential calibration is deliberately not part of the Brillouin-zone
tool either: it offers a geometric offset along the cut normal, in A^-1,
and leaves turning a photon energy into that offset to its own "append to
Kz map" feature.

## Twenty-fourth round: the figure composer

"Open in a new panel" now carries a **"Plot tools"** button, and it opens a
window whose job is one thing: producing the figure that goes in the paper.
It is seeded with that panel, or with a page of them from a map's new
**"Slice figure..."** button, or from the launcher's **"Plot as a figure..."**.

### One model, one painter

`tools/figure.py` holds what a figure *is* — panels, style, annotations — laid
out in **millimetres**, with no Qt in it. `ui.figure.FigurePainter`
draws that model onto any Qt paint device: the preview on screen, the PNG,
the TIFF, the vector PDF and the SVG are the same call with a different
device. A figure therefore cannot look right in the window and wrong in the
file, and "89 mm wide" means 89 mm in the PDF.

This is the opposite of the lab's `plot_tools2_demo.m`, which is forty
buttons that each call `set(gca, ...)` on whichever MATLAB axes was touched
last. Three things follow from having a model instead:

- **Panels can agree with one another.** Same colour scale, same limits,
  tick labels only on the outer edges: those are properties of the figure,
  not of one axes. "Apply to all panels" is one checkbox rather than a loop
  re-implemented inside every callback (the original's `Apply2All`).
- **What you see is what you get**, because one layout calculation feeds
  both the preview and the file.
- **A paper's figures can match.** The style — and only the style, never the
  data — goes to JSON and comes back onto the next figure.

### What the tools do

Six tabs, acting on the selected panel or on all of them:

- **Panels** — journal presets (Nature, APS, Science, Elsevier column
  widths, each with a font size that suits it), width/height/dpi, the grid,
  gaps, shared colour scale, shared ranges, edge-labels-only, the data's
  aspect, and adding/removing/reordering/lettering panels.
- **Axes** — ranges, reversal, tick direction/length/minor count/sides,
  explicit tick steps, frame width and colour, background, and a **second
  scale** on the top or right edge (`value₂ = scale·value₁ + offset`, which
  covers binding against kinetic energy, Å⁻¹ against degrees, and k in units
  of π/a).
- **Text** — title and axis labels, editable, with `_` and `^` markup:
  `k_{||} (Å^-1)`, `E - E_F`. A bare `^` takes the sign and the word after
  it, because what people type is `Å^-1`. Font, sizes for labels/ticks/
  annotations, and the panel letter's text, position, size, colour and
  whether it sits inside or outside the frame (outside by default, since a
  letter inside lands on the data).
- **Colour** — colormap and reverse, the level window, gamma, interpolation,
  a **colour bar** with its own label, and a **scale bar** for a real-space
  map, whose length is a round number chosen from the span unless typed.
- **Overlays** — **reference lines** at typed values with colour, style,
  width and a label (and "From the fit", which takes E_F out of the
  dataset's own `fitEF.ef`); **text and arrows** placed by clicking on the
  figure; an **inset** magnifying a region picked by clicking two corners,
  with the region outlined on the main panel; and a list of **curves in data
  coordinates**. That last one is the hook a Brillouin-zone tool will fill:
  it hands over `nxs_figure.Overlay` polylines and the panel draws them, so
  the BZ panel can be written separately without touching any of this. "Copy
  to all panels" is the MATLAB tool's copy-lines-to-all-axes.
- **Export** — PDF (vector, with the data embedded as an image at the dpi
  set above), SVG, PNG, TIFF, the two-layer EPS-axes-plus-PNG route from the
  twenty-first round, the clipboard, and saving/loading the style.

An ARPES panel *is* an image, so the PDF embeds a raster rather than
pretending otherwise and producing 40 MB of little rectangles; the text and
the axes around it stay vector. The raster is sized from the panel's size in
millimetres times the figure's dpi, not from the device it is drawn on,
which is what lets a 600-dpi image sit inside a PDF measured in points and a
preview stay cheap.

### Slicing a cube into a page

A map or k-map viewer's **"Slice figure..."** is the MATLAB tool's
`massplott_3D`: choose which axis to slice, where (evenly spaced, or a typed
list), how far to integrate around each position, and how many columns. Its
guesswork is gone — the panel titles are written from the sliced axis's own
label and unit rather than a hard-coded list of four units, and slicing
along an angle gives the other angle against energy without a separate code
path. A constant-energy series opens at 1:1 and on one colour scale, because
both of its axes are the same kind of quantity; a series of cuts does not.

The launcher's right-click **"Plot as a figure..."** is `massplot_2D`: one
panel per selected dataset, with anything that is not a plain 2-D dataset
reported rather than silently squeezed in.

### Deliberately not ported

`hold on/off`, `axis square`, MATLAB figure numbers and the
`shading faceted` button are artefacts of MATLAB's state machine (Smooth
unticked already *is* faceted). `exclude/include from legend` needs a legend
concept these figures do not have yet. The interpreter choice is replaced by
the `_`/`^` markup, which is what the two cases an ARPES label needs.

## Twenty-third round: k conversion for a single Cut

*(Two fixes landed after this round shipped, both about the aspect ratio:
a cut of a converted k-map no longer opens locked at 1:1 — that panel is
momentum against **energy**, and an equal scale between Å⁻¹ and eV is not a
shape, it is a coincidence of the numbers; and locking the ratio on any cut
viewer no longer lets the axes drift off the data when the window is
resized, which needed the plot box of the twenty-second round to cover
panels with EDC/MDC curves as well. Both are described where they belong,
under that round.)*

A **"Cut k conversion"** button in the cut viewer opens a dialog that turns
one E-vs-angle spectrum into E-vs-Å⁻¹. The physics lives in `tools/cutk.py`
(no Qt), in the same conventions as the Map conversion: θ is the deflector
angle, φ the angle along the slit, and a direction is
`(sinθ cosφ, sinφ, cosθ cosφ)`.

### What makes a cut different from a map

A map is a surface in the (kx, ky) plane, so its conversion has somewhere to
put both momentum components. A cut is a **line** in that plane, and the
line generally does not pass through Γ — it is a chord at some perpendicular
distance from it. That distance cannot be recovered from the cut's own
intensities (a line carries no information about where the line is), so it
has to be supplied.

The lab's `k_space_conversion_cut.m` asks instead for a "centre line"
clicked on the cut and feeds it in as `sin(α + offset)`. That conflates two
different things:

- **where normal emission is** — an instrument property, which enters
  *inside* the sine and therefore sets the k axis's **scale**;
- **where Γ is** — a sample property, a rigid translation in the (kx, ky)
  plane applied *after* the sine.

Clicking the middle of a band and calling it the offset gets the scale
wrong, not just the origin. Over ±15° at 95 eV kinetic the true kx range is
±1.292 Å⁻¹; with a centre line picked 10° off it comes out as
−2.110 … +0.435 Å⁻¹ — stretched on one side, squashed on the other. So the
two are kept separate here: the geometry fixes the whole trajectory, and
where the axis reads zero is a display choice made afterwards. There is no
centre-line picker.

### Where Γ comes from

Three routes, in descending order of trustworthiness:

1. **"From a converted map…"** — the Map conversion's θ/φ offsets *are* Γ's
   angles (the user picked the point that became k = (0,0)), so they
   transfer to any cut taken in the same alignment. The button lists the
   converted maps in the launcher, copies their `kconv.*` settings, and
   compares the two files' `SRn`/`SRz` — if the manipulator moved between
   them it says so rather than pretending the alignment still holds.
2. **The file's own geometry** — the deflector position a cut was taken at
   is in the file (`actuator_1_1`, one point; previously dropped with the
   axis, now kept as `DeflectorAngle_deg`) and fills the "Deflector angle"
   box. It is half the geometry; where normal emission sits is a beamline
   convention, so Γ still needs one of the other two routes.
3. **Typed** — Γ's deflector angle, slit angle and the sample rotation.

### The axis

Default: the momentum **along the cut**, zeroed where the cut passes
closest to Γ. Since momenta are measured from Γ, that zero falls out of the
projection with no extra bookkeeping, and it is the coordinate a band
disperses in.

**"Measure from Γ itself"** switches to the distance from Γ,
`sign(k∥)·√(k∥² + k⊥²)`. That axis is honest about the cut never reaching
Γ: it has a real gap of ±k⊥ around zero, which comes back as missing rather
than interpolated across (`np.interp` would happily draw a line over it).
It also leaves more than half the picture empty, which is why it is not the
default.

### The rest

- **Resampling** is a one-dimensional interpolation per energy row — k along
  the cut is monotonic in the slit angle (checked, and it survives a ±90°
  span; the guard fires rather than producing nonsense if it ever is not).
  The MATLAB file builds an `ii*jj` point list and calls `griddata`, whose
  Delaunay triangles mix neighbouring energies into every output point and
  smear the Fermi edge diagonally.
- **The output grid** takes the map dialog's two input styles: element
  counts, or a resolution (Å⁻¹ for k, meV for E), mutually exclusive.
- **The edges**: k scales with √E, so the lowest-energy rows do not reach as
  far as the highest. Off, the ends are missing there (NaN, as the MATLAB
  version silently produces); **"Trim to the range every energy covers"**
  cuts the output back to what all the rows share.
- **The readout** gives what the axes cannot: how far the cut runs from Γ,
  how much that varies along it and across the energy window, and the
  output's own step for scale. Phrased as how much "a line at constant k⊥"
  is an idealisation — the conversion follows the exact arc, so the number
  is not an error in the result. All of it is written into the converted
  dataset as `kcut.*`.
- **Not offered where it would be meaningless**: the button is disabled for
  a cut already in momentum, a cut taken along an arbitrary path through a
  map, or a slice through a map along the deflector axis, and says which.
- Units are Å⁻¹ only (the MATLAB file's output was in π/a, despite the
  label).

## Twenty-second round: the plot box takes the data's shape

An equal-aspect image in a panel that is not the image's shape can be
handled two ways: stretch the picture, or widen the view range past the
data. pyqtgraph does the second, and that one decision was behind all three
complaints this round. A square-ish deflector map in a wide window ended up
with a deflector axis running to ±40° when the data stopped at 12; the axes
no longer touched the data; and an export taken from that view was either
distorted (rendered at the window's proportions) or correct-but-padded
(rendered over the view's range, empty margins included). Neither is a
figure anyone would publish.

There is a third way, and it is the right one for a measurement: leave the
view range exactly on the data and make the **plot box** the shape the data
is. That is `AspectBox` in `ui/widgets.py`. It owns the image widget
and sizes it -- allowing for whatever the axes and their labels consume
inside it -- until the plotting area has the pixel aspect the data asks for,
then centres it. The empty space lands outside the frame, where it belongs.
With the range no longer being padded, pyqtgraph's own aspect lock is
switched off: the geometry already guarantees the scale.

So, at every window size and for every panel: **both axes sit exactly on the
data**, and an equal-ratio image keeps its true proportions. Zooming still
works -- the box re-shapes itself to the zoomed range instead of padding it
back out. Nothing had to be made non-interactive.

Two details, because they took some finding:

- pyqtgraph re-lays its axes out on the next turn of the event loop, not
  inside `setGeometry`. So solving the geometry twice in one pass compares a
  widget that has already been resized against a ViewBox that has not, and
  lands *further* from the answer than it started. `AspectBox.relayout` is
  therefore one pass, and `_verify` does the iterating -- one correction per
  event-loop turn, out of a budget refilled by each real resize. A good
  result deliberately does not spend the budget, because the axis widths can
  still change afterwards (shorter tick labels, a wider plot area) and knock
  the shape out again; `sigResized` then triggers a re-check.
- `set_axis_range` now marks the range as deliberate *before* moving it.
  `setXRange` emits while it runs, the box reshapes, and the resize that
  follows would otherwise see a view nobody had touched yet and fit it
  straight back to the data -- so typing a range did nothing at all.

Panels with EDC/MDC curves were briefly an exception — their image's edges
are pinned to the curves beside them and cannot move on their own — and it
showed the moment anyone ticked "Lock x/y" on a cut: the padding came
straight back. So the box now holds the image *and* its curves as one
block. The layout divides that block in fixed proportions (4:1 each way),
so resizing the block gets the image's share right, and everything stays
pinned. `AspectBox` therefore takes two arguments: the widget it resizes
and the view whose plotting area must come out the right shape — the same
widget on a plain panel, the block on a curves panel, identical arithmetic
either way.

One guard rail: a ratio can ask for a shape the panel cannot hold. 1:1
between 28° and 0.8 eV wants a strip 35 times wider than tall, which on a
curves panel leaves the image nothing at all. The wanted shape is capped to
what still leaves a 48-px plot, and the cap — not the original — is what the
geometry solves for, so it converges instead of chasing an impossible
target. The axes still sit exactly on the data; what is given up is the last
of the scale, on a setting that was already asking for something unreadable.

Consequently:

- **The export has one behaviour, not two.** The "Region" choice is gone;
  every export covers what the panel shows, and the image's proportions come
  from `export_aspect()` -- the data's own shape when a ratio is locked, the
  panel's shape when it is not, and never the window's. The same export now
  comes out of a wide window and a tall one.
- **pyqtgraph's stock "Export..." is gone from the right-click menu**
  (`scene.contextMenu = []` in `strip_stock_menu`), so the Export submenu is
  the only way out and there is no second door to a worse answer.

## Twenty-first round: one export, in the right-click menu

Exporting used to be two buttons on every viewer's toolbar ("Data (.txt)",
"View (.png)") writing whatever the window happened to hold. Both are gone.
Export now belongs to the **panel**, not the window -- a spatial scan shows
two pictures, and right-clicking the one you want is the only unambiguous
way to say which -- so it lives in each panel's right-click menu, as an
`Export` submenu with exactly three entries. They open one dialog
(`ExportDialog` in `ui/widgets.py`) on the matching tab.

All three cover the same region -- whatever the panel is showing -- which is
the reason they are one dialog rather than three menu entries: an axes-only
EPS is worthless unless it covers exactly the range of the image that goes
inside it. (This started as a "Whole data / Current view" choice at the top
of the dialog; the twenty-second round removed the choice by removing the
reason for it.)

1. **Data (.nxs).** The panel's array, both axes with their labels, and the
   file's metadata, written through the same `save_dataset` the launcher
   uses -- so the file reopens here as a dataset of its own. Three extra
   fields record where it came from (`export.source`, `export.panel`,
   `export.region`). The plain-text export it replaces wrote three separate
   `.txt` files per panel with the axes in files of their own; nothing read
   them back.

2. **Image.** Rendered by `nxs_export.render_rgba`, not screen-grabbed: the
   data is resampled onto the requested pixel grid, so 3000 px of export is
   3000 px of data rather than a blown-up copy of a 600-px window. Format
   (PNG / TIFF / JPEG), width and height in pixels with the on-screen aspect
   kept by default, and a dpi that goes into the file as metadata for
   whatever lays the figure out. The panel's colormap, level window and
   gamma are all applied exactly as on screen, and "Interpolate between data
   points" starts wherever the panel's own Smooth box is.

   **"Include the axes, ticks and labels"** is the switch that matters. On,
   the frame is drawn at the export's own resolution (`render_framed`, using
   the same ticks as the EPS) and the file grows by the margins the labels
   need. Off, the file is exactly the pixel size asked for, with everything
   outside the data -- NaNs, and the margins an equal aspect leaves --
   **transparent**, which is what drops into the EPS box.

3. **Axes only (.eps).** The frame, the ticks, the tick labels and the axis
   titles as vector PostScript with nothing inside the box, written by
   `nxs_export.axes_eps`. Box size in mm, font size, line width, minor ticks
   and whether ticks point inward. In Illustrator the text is text and the
   ticks are paths, so a figure's annotation stays editable while its data
   stays a raster. After saving, the status line reports the empty box in
   PostScript points -- position and size -- which is what to type into the
   transform palette when placing the image.

   Ticks are computed here (`nice_step` / `tick_values`, the 1/2/2.5/5
   ladder, nearest rung rather than always rounding up) rather than taken
   from pyqtgraph, because the EPS is a different physical size from the
   panel and wants its own tick density. Non-Latin-1 characters in a label
   (Å, a Unicode minus) are transliterated, since a broken EPS opens as
   nothing at all.

Also this round:

- **"→ list" is now "Save to the main list"** and **"Pop out" is "Open in a
  new panel"**.
- **The Fermi dialog's EDC is drawn over the fit window only**, not the
  whole energy axis. With the full axis the edge sat in a corner and the
  intensity scale was set by data nobody was fitting; restricted to the
  window it fills the panel and the residual is visible. Typing in the
  energy boxes or dragging the selection box redraws it immediately.
- **"Offset energy axis" moved up next to "Fit"**, since it is what you
  press straight afterwards.
- One thing worth knowing about, because it cost an afternoon: handing a
  panel a bound method of its window (`view.export_position =
  self.slice_label`) turned the whole viewer into a reference cycle, which
  meant the windows were freed by the *cyclic* collector at whatever moment
  it chose instead of at `close()`. Qt objects destroyed in the middle of a
  pyqtgraph call segfault the interpreter, reproducibly. The panels now hold
  plain values only and ask their window for the slice position through
  `self.window()` at export time; `ExportDialog` holds its view through a
  `weakref` for the same reason.

## Twentieth round: Fermi-level fitting

A **"Fermi level"** button in the cut viewer opens a modeless dialog that
fits the edge of whatever the cut shows. Rebuilt from the lab's
`FitFermiSurface` class rather than translated -- see the list of what was
wrong with the original at the end of this section.

### The model

```
I(E) = [ (a0 + a1*(E-EF)) * f(E; EF, T) + (b0 + b1*(E-EF)) ] (*) G(FWHM)
```

Seven parameters, each with a **Hold** box and an **uncertainty** column:
E_F, temperature, resolution (as FWHM, stated), the density of states at
E_F and its slope, and the background and *its* slope. The MATLAB model has
no background slope and never fits the background, the DOS offset or the
DOS slope at all -- they are estimated once from the first and last 10 % of
the spectrum and then frozen, which is exactly what trades off against E_F.

**The temperature is held by default**, at the value the file recorded
(`SampleTemperature_K`, newly kept out of the hidden `TC1.*` rows for this).
On a single edge the temperature and the resolution are not separable --
the width measured is `sqrt((3.53kT)^2 + FWHM^2)` -- so fitting both gives
two numbers neither of which means anything alone. Unticking the hold is
allowed, and the summary then prints their correlation (it comes out at
-1.00 on synthetic data, which is the point).

### What it is fitted to

The **selection box on the image is the region**: its horizontal extent is
the angle range the EDC is summed over, its vertical extent the energy
window the fit runs in. Drag it there or type the four numbers in the
dialog -- they follow each other. "Estimate" reads starting values off the
data (edge position from the steepest slope, levels from either side,
resolution from the 16-84 % width with the thermal part removed in
quadrature), so in practice nothing needs typing.

### What comes out

- **Offset energy axis**: a copy of the dataset with the energy axis shifted
  so E_F is zero, listed as `[E-Ef]` and relabelled **E - E_F (eV)**.
- The fit, its error bars, its reduced chi-squared, the held parameters and
  the window all travel with anything it produces, as `fitEF.*` in the
  metadata.

### The advanced corner (rarely needed)

- **Divide the Fermi cut-off out** (`dFD`), for superconducting gaps. It
  divides by the resolution-broadened Fermi function *only* -- not by the
  fitted DOS, which the MATLAB version divides out along with the cut-off,
  removing real spectral weight -- and takes the background out before
  dividing and puts it back after, so dark counts are not inflated near
  E_F. Everything above `E_F + n kT` (n adjustable, default 4) becomes NaN,
  because the divisor there is a small number known only to the noise.
- **Fit E_F channel by channel**, with a ± neighbours setting for a
  reference too thin to fit alone and a "fit every Nth channel" step. The
  result is drawn over the cut.

### Fermi-surface correction from a gold reference

The FS-correction dialog grew **"From a reference (Au)..."**: pick a
reference measurement (a Cut, or a Map which is summed over the deflector),
fit its edge channel by channel, and the measured points land in the same
table hand-picked points go into -- where they are smoothed by the same
polynomial and can be pruned. A progress dialog with Cancel, since hundreds
of fits take tens of seconds.

Tested end to end: a synthetic gold cut whose edge bends by
1.2e-3 eV/deg^2 comes back corrected with 5e-5 eV/deg^2 of curvature left,
i.e. flat to the per-channel noise.

### Bugs found in the MATLAB original (worth knowing if you still use it)

1. `fit.m` writes `2 * resolution` back into the table to show FWHM while
   the model reads it as HWHM -- so **the resolution doubles every time you
   press Fit**, and "Corr FS" divides with a resolution twice the fitted
   one.
2. `getFromFigure`/`getFromWorkspace` set E_F's **starting value to the
   energy span** (`max - min`) rather than a position.
3. `corr.m` divides by the whole model, background and DOS slope included,
   with no cut-off above E_F, then renormalises by the resulting (often
   divergent) maximum.
4. `fdDistribution.m` assumes an **ascending** energy axis and computes its
   step as `range/N` instead of `range/(N-1)`, which can slip the returned
   window by one bin.
5. The parameter table's default temperature (70 K) is outside its own
   bounds (5-50 K); the resolution's lower bound of 0 divides by zero in
   the kernel; `a0`'s bounds are not normalised with `a0` itself.
6. The tooltip's sign for `a0` is the opposite of the code's.
7. Statistically: unweighted least squares on Poisson counts, no
   uncertainties at all, and `fmincon` where a bound-constrained
   least-squares solver belongs.

## Nineteenth round: X/Y spatial axes, a separate information window

1. **The spatial axes are called X and Y.** The beamline uses two stage
   pairs -- ST/SZ (coarse, mm) and PIX/PIY (piezo, µm) -- and both mean the
   same two directions, so `_spatial_labels` now returns "X (mm)" / "Y (µm)"
   instead of the stage's own name. Data from another set-up gets X/Y too
   rather than a name only this beamline would recognise. The stage names
   are not lost: they are recorded as `Spatial.X_actuator` /
   `Spatial.Y_actuator` in the metadata. The first actuator the file lists
   is taken as X (the horizontal axis), which is this beamline's order --
   ST and PIX are listed first; a file that lists them the other way round
   is **flagged** (a warning, and a note in the metadata) rather than
   silently transposed, because which array axis is which is decided by
   length matching, and relabelling alone would move the label away from the
   data. If a real file ever trips that warning, say so and the parser can
   swap the axes properly.
2. **"Show information" splits in two.** The **Data Information** grid stays
   in the launcher, under the list, since it is what you compare between
   files while browsing -- its group box is titled with the dataset it
   describes. The **metadata** (a hundred-odd fields), its filter and the
   CSV export moved to a window of their own
   (`ui.main_window.InfoWindow`), made the first time it is asked for: under
   the list it left neither the list nor itself much room, and in its own
   window it can sit beside a viewer or stay open while other files are
   browsed. Both are still filled only on request.
3. **"Read cursor" is gone** from the conversion dialog. The offset boxes
   already follow the readout cursor for as long as the user has not typed a
   value of their own, which is every moment the button would have been
   useful in.
4. **Two fixes to the rotation picker**: a single picked point is now
   measured from the **chosen k-space origin** (the theta/phi offsets), not
   from the contour's (0, 0), so the direction it defines is the one that
   will actually run through the origin of the converted map; and pressing
   Convert clears the picked points, the rotation line and the readout
   cursor off the contour, instead of leaving the set-up overlays behind.

### Swap X and Y in the spatial view

A lab whose own X/Y convention is the other way round would see every map
transposed, so the spatial panel's right-click menu has **Swap X and Y**.
It exchanges the two axes *of the display only*: the labels follow, the
cursor stays on the same measured pixel (the spectrum beside it does not
even flicker), a selection box drawn on the swapped map integrates the
region it encloses, and the readout still names the position in the
measurement's own X and Y. Every export and every saved slice keeps the
orientation the measurement had, so how the map was being looked at never
leaks into the data.

Internally the window keeps the map in the file's own (y, x) order
(`_spatial_map`) and only the drawing is transposed; `_data_pixel()`
translates what the display reports back to measured indices, which is the
one place the swap has to be remembered.

## Eighteenth round: the k-conversion inputs, and a launcher that operates on data

### Sample rotation picked off the contour

"Set rotation from contour..." under the Sample rotation box turns on point
picking: click **two** points along the direction that should end up
vertical, and the rotation that stands it upright is filled in; click
**one** and the line from the origin to it is used, since that is the only
line one point defines. Rotation is about the origin, so only the direction
matters, not where the pair sits; a third click starts a new pair.

`nxs_kspace.points_to_azimuth` is `90 - atan2(dy, dx)` in degrees, folded
into (-180, 180]. The sign was **not** taken on faith: `test_kspace.py`
converts a synthetic band running at 30 degrees with the rotation this
returns and measures where the band ends up -- 90.1 degrees from kx, i.e.
vertical.

### The output grid, by resolution or by count

The conversion dialog now asks for the grid either way round: **by element
number** (as before, the default) or **by resolution** -- kx and ky in Å⁻¹
per point, E in meV per point. Picking one disables the other, so there is
never a pair of numbers disagreeing about the same axis, and a note under
the boxes always states the other reading ("→ 0.0024 / 0.0024 Å⁻¹ and 1000
meV per point").

Turning a resolution into a count needs the k box before the conversion
runs, so `nxs_kspace.k_extent()` is split out of `convert_map` and used by
both -- tested to give exactly the range the conversion itself produces.

### The launcher: information on request

Selecting a row no longer reads the dataset. Reading one to describe it is
real work -- for a GB spatial scan it is seconds -- so the information panel
is filled by **right-click → Show information** on a single row. Double-click
still opens a viewer.

The panel itself is the **Data Information** grid: min, max, num and step per
axis, with each column headed by the axis's own label so "Z" never has to be
guessed at. A cut fills X (angle) and Y (energy); a map or k_map fills X, Y
(the two in-plane axes, angles or momenta) and Z (energy); a spatial scan
fills X, Y (sample coordinates), Z (analyser angle) and a fourth column E
(energy), which only it shows. The metadata table is kept below it -- it
answers different questions (photon energy, temperature, slit) and is now
filled by the same on-demand action.

### Data operations (push button in the launcher)

One dialog for the three MATLAB data tools, acting on **every selected
dataset at once**:

- **Truncate** (`td_demo.m`): per-axis from/to, by value or by 1-based index,
  blank for "as far as the data goes". Suffix `_tk`.
- **Self-normalise** (`self_normalization_demo.m`): tick one axis to divide
  every line along it by that line's own total, or two to divide every plane
  by the plane's total; an optional window restricts what is summed, and
  "normalise to the peak" divides by the maximum instead. NaNs are excluded
  from the sums and put back afterwards. Suffix `_s_nor`.
- **Compress** (`data_comb_resamp_demo.m`, combine only): a bin factor per
  axis; values are summed over each block and the axis takes the block's
  mean, so total counts are preserved. A remainder too short for a whole
  block is dropped. Suffix `_comb`.

**Resample is deliberately not included** (you asked for it to go): it keeps
every Nth point and discards the rest, which throws away most of the counts
in a photon-starved scan, where combining keeps them and improves the noise.

A batch only means something if the datasets really are the same shape --
one set of bounds cannot describe two different axes -- so the launcher
checks kind and every axis length before opening the dialog and lists what
disagrees instead of producing quietly wrong results. Results are listed as
computed datasets and behave like any other: reopen, save, operate again.

A spatial scan is included: its 4D cube is read into memory to be operated
on, and the dialog says so and asks first when that is more than ~2 GB.
`MemoryData` grew the SPEM accessors (`spatial_overview`, `frame_at`,
`frame_over_region`, `spatial_over_region`) so a compressed scan opens in
the spatial viewer exactly like a measured one.

## Seventeenth round: gamma, arbitrary cuts, Fermi-surface correction

### Intensity scaling: Gamma and Min%/Max%, Alpha removed

Each panel's level bar now carries the MATLAB tool's three contrast
controls: **Min%**, **Max%** (the same two levels as a percentage of the
image's own range, so a setting carries from one scan to the next) and
**Gamma**. Gamma follows `PlotSlices.m`'s `plot()` exactly -- clip to the
window, normalise, raise to the exponent, map back -- and is applied in
`_for_display` against the level window the bar is showing, so the numbers
on the bar keep meaning what they say and `last_frame` (every export) stays
raw. The **Alpha** control is gone: on an image view the stock version did
nothing, and dimming the picture was never worth a control.

### Arbitrary-direction cut

"Arbitrary cut" on a map window opens a modeless dialog with six numbered
point rows, as in `arbi_cut_plot_demo.m`. Click the contour to fill the
next free row (the path is drawn as you go and the numbers stay editable),
then "Plot cut": the cube is sampled along each segment at the grid's own
spacing, the segments are laid end to end, and the result is listed as a
`cut` dataset with **distance along the path** as its x axis. It opens
straight away, with the corners marked by green dashed lines -- so a kink
in a band is not mistaken for physics when it is only a corner in the path.
"One dataset per segment" lists each segment separately instead.

Two deliberate differences from the `.m` file, both in `tools/analysis.py`:
sampling runs along the segment rather than over x (no special case for
near-vertical lines, and evenly spaced samples in every direction), and
segments are joined by dropping the duplicated corner sample rather than
writing a NaN there and calling `inpaint_nans` to repair it.

### Fermi-surface correction

"FS correction" opens a dialog that collects points clicked along a feature
that should be flat, lists them in an editable table, fits a polynomial
(order 1-4, default 2 as in `correction.m`) and draws it over the data as
you go. "Correct" shifts every angle column in energy so the feature comes
out flat, and lists the result.

- On a **cut** it corrects that cut.
- On a **map** or **k_map** the button is only on the **slit** cut window,
  because the curvature being corrected is the analyser's, along its exit
  slit -- a deflector-vs-E cut is a different axis entirely. The fit is made
  on the slit image and applied to the *whole* cube, so one fit straightens
  every deflector position, exactly as `correction.m` corrects its whole 3D
  array from a single slice. The result is listed as a `map` / `k_map`.

The energy axis grows by the span of the fit so nothing is cropped, and the
bins no column reaches are left NaN rather than filled -- there is no
measurement there. The bin width comes from the axis itself
(`mean(diff(E))`) rather than the `.m` file's `range/N`, which is off by one
bin in N. Both branches of the original (positive and negative quadratic
term) reduce to the same shift, so there is one code path here.

### Also in this round (the "B" list)

- **Slice -> list** puts the displayed slice in the launcher as a dataset of
  its own (MATLAB's "Save Slice", which put it back in the workspace), and
  **Pop out** opens it in its own frozen window ("Plot in new figure"). Both
  are on every viewer; a two-panel window asks which panel.
- **The integration window stays inside the data**: a `+/-` half-width pulls
  the slider's limits in by that much, so a slice near the end is never
  quietly integrated over fewer points than one in the middle. A width wider
  than half the axis is clamped rather than refused.
- **mean or sum** over that window, as a combo next to it. **mean is now the
  default** (it was sum): widening the window then leaves the brightness --
  and so the levels -- where they were, and the numbers agree with the
  MATLAB tool.

Computed datasets (k-maps, arbitrary cuts, saved slices, corrected maps) are
all `MemoryData`, which quacks like `NxsData`: they open in the same
viewers, export through the same code and save to `.nxs` in the same format,
with the settings that produced them (`arbcut.*`, `fscorr.*`, `slice.*`,
`kconv.*`) attached.

## Sixteenth round: drawing like `shading interp`, and typed slice positions

Two things taken from `PlotSlices.m`.

**1. "Smooth" now interpolates instead of only blurring.** The old Smooth
was a Gaussian on the *data* grid, which is why it never removed the
staircase: a deflector map is a few tens of points across, so one data
point covers a big block of screen, and blurring on that grid just gives
softer blocks. MATLAB's plots look continuous because `pcolor` is drawn
with `shading interp`, i.e. bilinear interpolation across every cell. So
`_for_display` now does:

1. optional Gaussian blur on the data grid (`smooth2d`, the "blur" box,
   default **0** -- no blur at all, which reproduces MATLAB's look);
2. bilinear resampling onto a finer grid (`upsample2d`), by
   `interp_factor(n)` per axis: enough samples to reach ~500 along each
   axis, capped at 16x, and 1x (nothing to do) for data that is already
   dense.

The blur comes first because it is far cheaper on the coarse grid, and a
smoothed redraw costs a couple of milliseconds. The resampled image covers
exactly the same extent -- only the `scale` passed to `setImage` changes,
divided by the per-axis factor -- so nothing moves under the cursor, the
readout still reports measured values, and `last_frame` (what every export
writes) stays raw. The interpolation is NaN-aware: a missing point is
filled from its neighbours instead of spreading a hole.

**2. Pos and Ind boxes on every slice slider** (`_SliceControl`), as in the
MATLAB tool: type an energy or an angle in **Pos**, or a 1-based index in
**Ind**, press Enter, and the image follows; `<` and `>` step one point.
Both boxes also read out where the slider is. A typed position snaps to the
nearest measured point and the box is rewritten to that point's value --
the data only exists there, so showing anything else would misdescribe
which slice is on screen. Indices are 1-based deliberately: that is what
the MATLAB tool shows, so the same number means the same slice in both
programs.

## Fifteenth round: the view controls come out of the right-click menu

pyqtgraph's stock context menu is gone from every plot in the app, and what
was worth keeping from it is now a visible row -- the `ViewOptionsBar` --
directly under each window's colormap:

```
View  [panel v]  X: [min][max] [x]Auto []Inv | Y: [min][max] [x]Auto []Inv | [Reset] | []Grid  Alpha [1.00]
```

- **X / Y min and max**: type a range, press Enter. The boxes also follow the
  view, so they double as a readout of where the panel is looking.
- **Auto**: this app's tight fit to the data (`fit_to_data`), not
  pyqtgraph's padded auto-range -- the panels here sit against their data.
  Typing a range, or panning/zooming by hand, switches it off, so a later
  resize cannot undo what the user just did. The two axes move together
  because an aspect lock ties them together.
- **Inv**: reverses that axis.
- **Reset**: fit back to the data and put Auto back on.
- **Grid** and **Alpha** apply to *every* panel in the window, like the
  colormap above them. Alpha is the image's own opacity: the stock "Plot
  Options -> Alpha" only ever touched curves, so on an image view it did
  nothing at all.
- A window with two panels (a spatial scan) gets a **panel chooser** in
  front of the axis boxes; a one-panel window has nothing to choose and
  hides it.

**Removed rather than moved**, because none of them mean anything here:
"visible data only" and "auto pan only" (they matter for live-updating
curves), **axis linking** (these panels have different axes -- the EDC/MDC
plots are already linked in code, which is the only linking that is ever
right here), and **mouse enabled**, which is now always on: `strip_stock_menu`
calls `setMouseEnabled(True, True)` so dragging can never be switched off by
accident. The **"Plot Options"** menu goes with them
(`PlotItem.setMenuEnabled(False, None)` -- the `None` leaves the ViewBox menu,
and therefore this app's own actions on it, alone).

What is left on the right-click menu is only this app's own entries: readout
cursor, selection box, integrate selection, open this position in a new
window.

**Layout tidying** in the same round:

- Both control rows are grouped with vertical rules and bold captions
  (`separator()`, `section_label()`): window buttons | Colormap | ... |
  Export on the first row, the view controls on the second. Exports sit at
  the far right, away from the controls that change the picture.
- Both rows live in a `scroll_strip()`: a fixed-height scroll area that
  looks like an ordinary row when there is room and grows a small horizontal
  scrollbar when there is not. Without it the sum of a row's widgets becomes
  the window's minimum width -- which would have defeated opening a
  tall-narrow map at its own aspect ratio (the map window's floor was
  1244 px; it is 520 px again).
- The map window's cut buttons are now "Deflector cut" / "Slit cut" (the
  tooltips still spell out what each one integrates over), the axis boxes
  dropped their up/down arrows (nobody nudges a view range one step at a
  time, and the arrows cost more width than the digits), and the EDC/MDC
  controls are fenced off with their own rule that appears and disappears
  with them.
- `size_to_data_aspect`'s `extra_height` grew to 260 to account for the new
  row, so a map still opens at its data's aspect.

## Fourteenth round: a usable conversion dialog, and a list you can manage

Four fixes, all reported against the k-conversion round.

1. **The conversion dialog no longer blocks the map window.** It was run
   with `exec_()`, which is modal: the prompt said "pick a point on the
   contour" while making the contour impossible to click. It is now created
   with `setModal(False)` and `show()`, kept on the window as
   `_kconv_dialog`, and wired to the contour's `readoutMoved` signal, so
   dragging the cursor updates the theta/phi boxes live while the dialog
   stays open. Clicking "map k conversion" again raises the dialog that is
   already open rather than making a second one.
2. **Converted maps can be named, and cannot collide.** The dialog has a
   name box, pre-filled with `<file> [k]` and bumped to `<file> [k] 2`,
   `... 3` if that is taken. Whatever name comes out is still passed
   through `unique_name()`, so even typing the same name twice by hand
   gives two distinct rows.
3. **Right-click menu in the launcher list** (`show_list_menu`), with
   multi-selection enabled (`ExtendedSelection`):
   - *Rename* -- single row only; refuses a name already in the list.
     Renames the record too, so a later save uses it. The file on disk is
     never touched.
   - *Remove from list* -- one or many, with a confirmation. It only drops
     rows; **nothing on disk is deleted**. The confirmation calls out that
     an unsaved converted k-map is lost when its row goes.
   - *Save* -- one or many into a single `.nxs` file, one top-level entry
     per dataset, named as they appear in the list.
4. **A cut through a k_map opens at 1:1**, like the contour it comes from
   -- both axes are Å⁻¹, so an equal aspect is the meaningful default.
   A cut through an *angle* map still does not: degrees against eV has no
   natural 1:1.

**The saved format is this program's own**, written by
`loader.nxs_file.save_dataset()`: an HDF5 container whose entries carry a
`nxsloader_kind` attribute. `_classify_entry()` checks for it first (case
5), so a saved file reloads through the ordinary "Add files" path and
comes back as ordinary rows -- values, all axes, axis labels, metadata and,
for a k-map, the `kconv.*` settings that produced it. Two things to be
aware of: saving a SPEM writes its **full 4D cube**, which for a real scan
can be GB-sized and is read into memory to do it; and a k-map lives only in
memory until it is saved.

`test_browser.py` covers all of this end to end, including the round trip
(values, axes, labels and settings identical after save + reload).

## Thirteenth round: angle -> momentum conversion for Maps

`tools/kspace.py` ports `k_space_conversion_demo.m`. The contour window grew
a **"map k conversion"** button: pick a point on the constant-energy
contour to set the theta/phi offsets (editable afterwards), choose how many
kx/ky/E elements the result should have (100/100/native by default), and
confirm -- the converted cube appears in the launcher list under the same
name with kind `k_map`, double-clickable and viewed exactly like a map,
with axes in Å⁻¹.

Implementation notes:

- The conversion is done as an **exact inverse mapping**: for each point of
  the target (kx, ky, E) grid, work out which angle it came from and
  interpolate the original cube there (`RegularGridInterpolator`). The
  `.m` file goes the other way -- forward-map every angle point, then
  Delaunay-interpolate the scattered result -- which is slower and blurs
  the data. `convert_map_forward_reference()` implements the `.m` route and
  is used by `test_kspace.py` to cross-check the fast one: 0.07 % rms,
  0.33 % max difference.
- **Sign convention:** the `.m` file's rotation does *not* put the picked
  point at k = (0, 0) -- feeding it a (5°, -3°) pick lands that point at
  k = (+0.173, -0.104) Å⁻¹. Since the whole point of the UI is "this
  feature is my origin", `origin_rotation()` negates both offsets so the
  picked point really is the origin. **The numbers in our boxes are
  therefore the negatives of the ones in MATLAB's.**
- No `cut` conversion and no kz yet, as requested.

## Twelfth round: several Cut/Map datasets in one file

A real `.nxs` file can hold up to ~6 separate Cut or Map measurements. The
loader previously showed one row per file and opened whichever entry it
picked. Now `list_datasets()` returns every measurement in a file and the
launcher lists each as its own row (keyed on `(path, entry)`); double-click
opens only that one. Classification is per-entry and content-based
(`_classify_entry`), so a navigation image saved next to a map no longer
confuses the kind label.

## Eleventh round: opening a large SPEM scan is now instant

Reported: double-clicking a SPEM file was very slow on some machines, or
appeared to do nothing until clicked several times.

**Diagnosis.** `_parse_case2` read the whole cube with `data_12[()]`. A real
scan is 46x91x96x941 = 378 million points — 0.76 GB as uint16, 1.5 GB as
int32 — and the spatial overview then summed it again on a transposed
(non-contiguous) array. The GUI thread blocked for seconds with no feedback,
and Qt discards clicks arriving while it is blocked: hence "needs several
double-clicks". On a network share (the lab's `Y:\`) it is far worse.

**On the offered trade-off: binning the energy axis to 50 meV would not have
helped.** To bin you must still read every byte; binning saves RAM, not I/O,
and the bottleneck is the read. So nothing is compressed and no resolution
is lost anywhere.

What was done instead:

1. **The cube stays on disk** (`LazyCube`). It is presented in (y, x, k, E)
   order and indexed exactly like the numpy array it replaces, but each
   index reads only the selection. A cursor move reads one (k, E) frame —
   about 0.2 ms.
2. **The spatial overview comes from the reduced preview cubes the file
   already carries.** `data_11` is the cube summed over energy and `data_01`
   over the slit angle, so summing either one's last axis gives the identical
   overview at a thousandth of the data volume. It is *validated* against the
   real cube on three pixels before being trusted — if a file's preview turns
   out not to be a plain projection, the code warns and streams the cube
   instead, row by row, so the answer is always exactly right.
3. **Region integrations stream** a spatial row at a time, so a large
   selection never needs to fit in memory.

Measured on a 194 MB benchmark (`make_bench_file.py`; the real file is ~4x
larger): **window ready in 11 ms**, against ~1.8 s for the old
read-everything path on the same data when gzip-compressed — and every lazy
result is numerically identical to the full-cube computation.

One subtlety this exposed: HDF5 gives a process **one** underlying handle per
path, so two `h5py.File` objects for the same file are not independent and
closing either invalidates both. `NxsData.acquire()` therefore hands out one
reference-counted instance per file, and the file closes only when the last
window using it does — in any order.

## Tenth round: the frame really does hug the data now

Reported against real files: blank margins around every image, the spatial
map sitting off-centre rather than in the middle of its frame, and the
contour floating in a wide empty box. Three separate causes:

1. **Aspect-lock expansion was uncentred and excessive.** pyqtgraph applies
   an aspect lock by stretching a range around the *previous* view centre,
   which left the data pushed to one side and could inflate both axes well
   past what the lock needs. `fit_to_data()` now solves the lock's own
   constraint directly — `(x_span / width_px) == aspect * (y_span / height_px)` —
   so exactly one axis ends up on the data extent and only the other grows,
   symmetrically about the data's centre. Unlocked panels (E-vs-k, cuts) come
   out tight on both axes.
2. **The fit ran before the layout settled.** With an aspect lock the correct
   range depends on the viewbox's pixel shape, which is only final after the
   last layout pass — later than any `showEvent`. The fit is now re-applied
   from the viewbox's own `sigResized`, so it stays tight through every
   resize. A deliberate zoom or pan is remembered and left alone
   (`sigRangeChangedManually`), so re-fitting never fights the user.
3. **Overlays took part in the ranging.** Cursors, crosshair lines, the
   floating readout label and the selection ROI are now added with
   `ignoreBounds=True`: they are annotations, not data, and must not drag the
   range around.

Also: the **spatial cursor is just a marker now** — its full-width crosshair
lines were noise, and the position is already in the status bar. (The EDC/MDC
crosshair on the spectrum panels stays: there the lines carry information,
showing which row and column the curves are taken from.)

**A trade-off worth knowing.** An equal aspect and a frame tight on *both*
axes are mutually exclusive unless the panel happens to have the data's
shape, so a tall-narrow contour in a wide window must have margins
left and right. Two things help: the contour window now opens sized to its
own data aspect (`size_to_data_aspect`), and unticking **Equal ratio** /
**Lock x/y** under any panel gives a fully tight fill immediately. If the
margins you are seeing are *inside* the image — empty detector rows and
columns with no counts — no fit can remove those; say the word and I will add
a "trim empty edges" option.

## Ninth round: kind labels, box axes, per-window colormaps

1. **The file list shows the data kind instead of a thumbnail.** Each row
   reads `filename   [SPEM]` / `[Cut]` / `[Map]`, or `[unknown]` for a layout
   that is not recognised (including unreadable files — `probe_kind` never
   raises). The label comes from `loader.nxs_file.probe_kind()`, which reads
   only the group structure and a couple of dataset *shapes*, never the
   detector cube, so adding a folder of multi-GB scans is instant; the old
   thumbnail column had to fully load every file. The display text carries
   the label, so the real filename travels in the item's `UserRole` data
   rather than being parsed back out.
2. **Box axes, tight to the data.** All four axes are drawn with values on
   the bottom and left only, and the view padding is zero so the frame sits
   against the data instead of floating in blank space. The spatial map and
   the constant-E contour now default to a 1:1 aspect (both of their axes
   carry the same unit; a stretched spatial map misrepresents the shape of
   what is on the sample).

   Removing the padding exposed a half-pixel error worth knowing about:
   pyqtgraph puts the *corner* of the first pixel at the `pos` you give it,
   so passing `axis[0]` shifted every image half a pixel against its tick
   labels and made it extend a full step past the last data point.
   `_image_origin()` now offsets by half a step so pixels are **centred** on
   their axis values, and a tight view spans the data range plus half a pixel
   at each end. This was invisible while the view had padding.
3. **Colormap moved into each window.** The launcher no longer has a global
   control; every viewer (and every pop-out snapshot) carries its own
   dropdown and Flip box, so two windows can be tinted differently — which is
   what makes two positions easy to tell apart when comparing them. The
   launcher still remembers a starting value between sessions, taken from the
   most recently opened window.

## Eighth round: launcher + separate viewer windows

The single window with a stacked page per data kind is gone. `ARPES_viewer.py`
is now just a **launcher** — folder, file list, colormap, metadata — and
**double-clicking a file opens it in its own window** (`ui/windows.py`):

| double-clicked | opens |
|---|---|
| real-space scan | `SpatialScanWindow`: the spatial map *and* the E-vs-k spectrum at the cursor, together |
| single cut | `CutWindow` |
| deflector map | `ContourWindow`: the constant-energy contour, with two buttons that open `MapCutWindow`s for the two orthogonal cuts |

Single click loads metadata only; double-click (or "Open selected") opens the
viewer. Re-clicking a cut button raises the window already open rather than
making a duplicate. Several files — or a contour plus both its cuts — can be
on screen at once and arranged freely, which the fixed layout could not do.

Consequences worth knowing:
- The **colormap is global**: the launcher's dropdown retints every open
  window, including cut and snapshot windows.
- **Export moved into the viewers.** Each window exports what *it* shows, so
  an export can never disagree with the window you were looking at; a map's
  cut window writes its own cut under its own name. Metadata CSV stays in the
  launcher, since metadata describes the file, not a view.

Also in this round: the readout cursor now draws **crosshair lines** like the
spatial map's, coloured to match the curves they feed (the vertical line
feeds the EDC, the horizontal one the MDC), plus **shaded bands** showing the
"+/-" integration windows, so a widened EDC/MDC shows on the image exactly
how much it is summing. The bands re-centre as the cursor moves and keep
their width.

## Seventh round: ARPES-standard curve layout, slim level bar

1. **EDC/MDC moved into the standard ARPES arrangement.** They no longer sit
   side by side under the image where they crowded it out. The **EDC** is now
   to the **right** of the image with counts horizontal and energy vertical,
   its energy axis linked to the image's; the **MDC** is **below** it with
   its angle axis linked to the image's. Panning or zooming the image moves
   all three together, and features line up across panels by construction.
   The image keeps ~80% of both dimensions; the curves take a quarter.

   Two details this needed:
   - `CurveReadout` is now a plain `QObject`, not a widget. Its two plots and
     its "+/-" boxes are separate widgets that `ImagePanel` places into its
     own grid — a single composite widget could not have put one plot beside
     the image and the other below it.
   - Linking views makes the *data* ranges match but not the on-screen pixels,
     because pyqtgraph sizes each axis to its own tick labels. `_pin_axis_geometry`
     fixes the left-axis width and bottom-axis height across the three panels,
     and the per-curve positions are reported in a label in the controls row
     instead of as plot titles — a title adds height to the top of the EDC
     only, which pushed its plot area 30 px below the image's. The smoke test
     now asserts pixel alignment, not just that the link exists.
   - The EDC hides its own energy axis and the MDC its own angle axis, since
     the image they touch already carries that axis.

2. **Level bar replaces the histogram widget.** pyqtgraph's HistogramLUTWidget
   (the intensity distribution plus a gradient strip down the side of every
   panel) is hidden. In its place, directly above each image, is a slim 28 px
   bar with two draggable handles and a numeric readout of the two values —
   no distribution drawn. It re-points itself at the new value range whenever
   the displayed data changes, so it never shows stale numbers.

## Sixth round: EDC/MDC curves and pop-out comparison windows

1. **EDC / MDC merged into the readout cursor.** Turning on the readout
   cursor on any E-vs-k figure (the spatial page's spectrum panel, a Cut,
   and the map view's two angle-vs-E panels) now also draws, under that
   panel, the **EDC** (counts vs energy through the cursor) and the **MDC**
   (counts vs angle through it). Each has its own "+/-" half-width in
   physical units: 0 is the single row/column under the cursor, anything
   larger sums every channel inside the window, which is what makes a noisy
   single-pixel line usable. The curves appear and disappear with the cursor
   and re-point themselves whenever the underlying frame changes (cursor
   move, slice change, region integration).
   The constant-E map is angle-vs-angle, not an E-vs-k figure, so it has no
   curves — EDC/MDC would be meaningless there. `CurveReadout` does fall
   back to generic "profile" labels if it is ever attached to a non-energy
   axis.
2. **Pop-out comparison windows.** Right-click the spatial panel ->
   "Open this position in a new window" snapshots the current spectrum into
   its own window. Any number can be open at once, which is the point:
   comparing several sample positions side by side. Each window holds its
   own **copy** of the data, so it stays frozen while the main cursor moves
   on, and each carries the full ImagePanel controls (Smooth, x/y ratio,
   readout cursor with EDC/MDC). Open windows follow the colormap dropdown
   and deregister themselves when closed.

   Window titles describe *what is actually displayed*, taken from a label
   recorded when the frame was built rather than re-read from the cursor
   position — those two disagree after a region integration, and a snapshot
   of an integrated region is titled `x[..], y[..] integrated` rather than
   being mislabelled as a single pixel.

## Fifth round: panel controls, smoothing, axis units, leaner metadata

1. **Selection box is easier to drive.** Thicker border (4 px) and much
   bigger drag handles (14 px vs pyqtgraph's 5). Whenever a box exists, a
   row appears under that panel with the two draggable corners as
   `x0/y0/x1/y1` boxes plus a **Refresh** button. The boxes are two-way:
   they follow the rectangle as you drag it, and typing values then hitting
   Refresh moves the rectangle there and re-runs the integration — so an
   exact region (x 110..120, y 110..120) is reproducible, which dragging by
   eye is not. Refresh and the right-click menu entry both raise the same
   signal.
2. **Smooth checkbox** under each panel, with a sigma box (pixels). The
   smoothing is **NaN-aware**: gaps are filled by normalized convolution
   rather than spreading NaN over their neighbours, which is the actual fix
   for the staircase edges of a sparsely-sampled scan. It is display-only —
   `last_frame` and therefore every export stays raw.
3. **x/y ratio control** under each panel: a lock checkbox plus a ratio box.
   The constant-E map, whose two axes are both angles, labels its checkbox
   **Equal ratio** and defaults the ratio to 1; the E cuts mix degrees with
   eV, where a 1:1 lock would be meaningless, so they get the generic
   "Lock x/y" instead.
4. **Axis labels with units**, read from the file rather than hardcoded:
   angles in `°`, energies in `eV`, and the real-space axes named after
   their actuators — `µm` when those names identify the PIX/PIY piezo stage,
   `mm` otherwise (`spatial_unit_for()`). Detection is name-based because
   the files carry no unit attribute on these datasets.
5. **Leaner metadata table**: `title`, `experiment_identifier`, `run_cycle`,
   `Machine.name/probe/type`, and the whole `User.*`, `TC1.*`, `TC2.*` and
   `ScanCfg.*` groups are no longer shown (74 rows instead of 95). The lists
   are `INFO_FIELDS_HIDDEN` / `INFO_PREFIXES_HIDDEN` in `loader/nxs_file.py`
   — delete an entry to bring a field back, nothing else needs changing.
   Note `User.mail` was not on the removal list so it is still shown; say
   the word and it goes too.

## Fourth round: export

*(Superseded by the twenty-first round: the text and PNG buttons described
here no longer exist -- exporting moved into each panel's right-click menu
and now writes `.nxs`, an image at a chosen size, or vector axes. Kept for
the record of what the original did.)*

The original `h5Loader.py` had save buttons (2D scan → txt, image →
clipboard, EDC → txt) that this port had not yet reproduced — a
`save_current_frame_txt()` existed but was wired to nothing, i.e. dead code.
Replaced with three working buttons in the left panel:

- **Data (.txt)** — writes *every* panel currently displayed, each as the
  2D array plus its two axes in separate files, following the original's
  `ScanData`/`RowRange`/`ColRange` convention so everything loads with a
  bare `np.loadtxt`. Per kind: `EvsK` + `SpatialMap` for a spatial scan,
  `Cut` for a cut, and `ConstE_map` + `Deflector_vs_E` + `Slit_vs_E` for a
  map. Arrays are written in `(x, y)` order matching their axis files.
- **View (.png)** — saves the main panel using the colormap *currently
  selected in the GUI* (flip included), so the file matches the screen.
  `nxs_tools.save_img` now takes an optional LUT instead of always using
  hardcoded `inferno`.
- **Metadata (.csv)** — the full metadata table as `field,value`.

Export always reflects what's on screen, so after a ROI integration it
writes the integrated frame, not the pre-selection one (`test_export.py`
checks this explicitly). With nothing loaded the buttons warn rather than
raise.

## Third round: colormaps, interaction, metadata

1. **32 colormaps**, in `tools/colormaps.py`: MATLAB's `gray`, `jet`, `hsv`,
   `hot`, `copper`, plus every (N,3) table from the lab's `Colormap.mat`
   (`blueblackred`, `rainbow`, `spectrum`, `terrain`, ...). A **Flip**
   checkbox next to the dropdown reverses any of them. The tables are
   embedded in the generated module as compressed base64, so the app needs
   neither scipy nor the `.mat` file at runtime — re-run `devtools/gen_colormaps.py`
   only if you want to add more. (`cm_copper` from the `.mat` and MATLAB's
   own `copper` are genuinely different tables, so both are kept, the
   latter as `copper_lab`; `cm_ring_area` is a scalar, not a colormap, and
   is skipped.)
2. **Map panel ratio fixed, controls moved above each panel.** The three
   slices are now three equal columns, each with its integration control
   directly above its own plot, so control and plot share a width by
   construction. The residual squashing came from `ImageView`'s autoRange
   re-enabling aspect locking, which `set_frame` now re-asserts off.
3. **Readout cursor** on every panel (cut, map, and both spatial-page
   panels): right-click → "Readout cursor" for a draggable crosshair that
   reports the data coordinates *and* the value under it.
4. **Rectangle selection with two-way integration** on the spatial page:
   right-click → "Selection box" on either panel, drag it, then right-click
   → "Integrate selection into the other panel". A box on the spatial map
   gives the E-vs-k spectrum summed over that patch of sample; a box on the
   E-vs-k panel gives the spatial map of that spectral feature. The update
   is explicitly triggered from the menu (not live-on-drag), and applying an
   E-vs-k selection deliberately leaves the navigation cursor where it was.
5. **Explicit file loading.** "Folder..." picks a directory and
   "Select file(s)..." opens a multi-select dialog; only the chosen files
   are loaded and thumbnailed. Nothing scans or loads a whole directory any
   more, which matters when the folder holds multi-GB scans. "Clear list"
   empties the picker.
6. **Full metadata** (~95 rows on a real file): NXentry provenance
   (`title`, `start_time`, `end_time`, `duration`, `run_cycle`), the `User`
   group, all MBS analyzer settings, monochromator, both undulators,
   machine status, both thermocouples (°C and K), scan configuration and
   trajectory, every sample/OSA/pinhole/zone-plate motor (offset-corrected),
   and the fine piezo stage — whose real path (`i12-m-cx1-ex-pi/{x,y,z}`)
   the `.m` file guessed wrong and left commented out. Device groups are
   read field-by-field rather than from a hard-coded list, so fields added
   by future beamline software show up automatically; the long
   `controller_record` Tango blobs are excluded. Keys are namespaced
   (`MBS.center_ke`, `Mono.energy`, ...) and a filter box above the table
   narrows it by substring.
7. **"Mouse Mode" removed** from every right-click menu. Note pyqtgraph
   0.14 dropped the `ViewBoxMenu.leftMenu` attribute that older versions
   exposed, so `strip_mouse_mode()` handles both layouts.

Also fixed along the way: `_SliceControl.width` would have shadowed
`QWidget.width()`, so the spin box is `width_spin`.

## Second round of real-file fixes (entry classification + UI)

1. **Entry classification is now content-based, not name-based.** The real
   ``Spatial_Scan.nxs`` bug: it has two top-level entries sharing the
   ``img1_000N`` naming pattern, but the second one doesn't actually contain
   deflector-scan data (missing `scan_data/actuator_1_1`). The old
   name-table lookup guessed "case 4" from the name and crashed deep inside
   the parser; `_classify_entry()` now inspects which datasets actually
   exist and classifies entries by content, falling through to
   "unrecognized" (never a crash) for entries that don't match anything --
   which is exactly what should happen to a broken/incomplete entry.
   `GROUP_NAME_TO_CASE` is now only a fallback for case 1 (which is
   unimplemented anyway, like the `.m` file).
2. **`map` thumbnails now show the (deflector, slit) angle-angle plane**
   (summed over E), not an E-vs-k plane.
3. **A colormap dropdown** (top of the file browser panel) applies to every
   view at once (`apply_colormap_everywhere`); pick from `inferno`,
   `viridis`, `plasma`, `magma`, `cividis`, `grey`.
4. **`KCubeExplorer` (the `map` 3D view) redesigned:**
   - Axes are now labeled/plotted as raw **angles** ("deflector", "along
     slit"), not kx/ky -- converting the deflector axis to momentum is a
     modeling choice (reference KE + small-angle approximation) that's kept
     as an opt-in (`NxsData.kcube()`) rather than baked into the default
     view.
   - Layout is a 2x2 grid for the three plot panels (equal size, proper
     aspect for each) with the slider/control column as a separate,
     fixed-width column -- it no longer competes for space with (or
     distorts) the plots.
   - Each slider now has a "+/-" half-width spin box (physical units) next
     to it: e.g. energy centered at 88 eV with a 0.25 half-width sums every
     channel in [87.75, 88.25] eV instead of showing a single-pixel slice --
     the usual way to reduce noise in a constant-energy ARPES map. Works the
     same way for both angle axes.

## Real-file quirks fixed since the first version (worth knowing about)

These were found by running against real SOLEIL files and are not visible
from the `.m` file alone:

1. **A "single deflector position" Cut can still be stored as a 3D array**
   (`(1, k, E)`) rather than the bare 2D `(k, E)` the `.m` file's shape
   implies. `_parse_case34` now accepts either and squeezes the size-1 axis
   away when present.
2. **Calibration scalars (`escalemin/mult/max`, `xscalemin/mult/max`) can be
   logged as 2D arrays** (e.g. one row per frame) instead of flat 1D
   vectors. `scalar0()` extracts a representative value the way MATLAB's
   `arr(1)` linear indexing does (first element, any shape) and warns if
   the values aren't actually uniform (which `arr(1)` implicitly assumes).
3. **A single `.nxs` file can contain more than one top-level NeXus entry**
   -- e.g. a navigation SPEM image recorded right before the actual
   deflector-angle map, both saved into the same file. `load_soleil_nxs()`
   no longer assumes exactly one; see `list_entries()` and the `entry=`
   parameter. The automatic choice (highest numeric suffix) is a heuristic,
   not a guarantee -- if a file's default choice looks wrong, pass
   `entry=` explicitly (or check the `_other_entries` field in `scan.info`,
   which lists what else was in the file).
4. `nxs_tools.save_img` used `matplotlib.cm.get_cmap`, removed in
   matplotlib >= 3.9; now tries `matplotlib.colormaps[...]` first.

`make_test_fixtures.py`/`test_process_nxs.py` have regression fixtures for
all four of these.

## Deliberate deviations from a literal line-by-line port

- `_find_actuator_group` looks up the SPEM1D case's scan actuator by name
  (`actuator_*_1`) instead of the `.m` file's positional
  `hinfo.Groups.Groups(3).Groups.Groups.Name` traversal, which only worked
  by accident for whatever one file it was written against.
- Unknown top-level group names raise immediately with a clear message
  instead of leaving `caseN` undefined (which is what the `.m` file does in
  the same situation).
- `tools/system.py` replaces `tools_packages`/`Win10Notif`/`win32*` with
  cross-platform equivalents so this runs without the lab's internal
  packages installed; the "Sync with running scan" button (which talked to
  `pyNanoScanSystem.SharedMemoryGlobals`, live-acquisition-specific) was
  dropped rather than ported, since it isn't part of the file-format change.
  Swap these back in on the lab PC if you want the original notification
  style / clipboard copy / live-sync back.
