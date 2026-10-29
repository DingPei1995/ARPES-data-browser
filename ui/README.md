# The windows: technical reference

How the main panel, the viewer windows, autosave and recovery, and long operations work. The code is in this folder (`launcher/`, `viewers/`, `panels/`, `common/`). For how to *use* them, see the [user guide](../docs/USER_GUIDE_en.md).

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

At the foot of the panel: the version and its date, the link to GitHub for
bug reports and updates, and the error logs -- **Open the error logs** (the
folder) and **Copy the newest error log** (its text to the clipboard): upload
the file with a bug report, or paste the text.

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

## Long operations

Reading a map, converting to k-space, resampling a cube, assembling a
CASSIOPEE folder: all of it runs on a worker thread with a progress bar that
moves and a Cancel that works. One job at a time, because HDF5 here is not
thread-safe. A job that finishes within 400 ms never shows a dialog.

Every operation that can take long on a big dataset or a slow computer
shows one, so the program is never silently busy: reading a whole map into
memory (browsing, a slice figure, the 3-D views), map and cut k conversion,
the FS correction of a map, an arbitrary cut, data operations (truncate,
normalise, compress), stacking or joining datasets, saving (`.nxs`, Igor),
ridge detection, the tight-binding fit, the fits along a cut and of the
Fermi edge, smoothing and the other processing, k_z conversion, de-gridding,
detector corrections, summing repeats, and in Calculated bands reading a
calculation, the bands on a plane, the calculated map and the 3-D Fermi
surface. A step that cannot say how far it is (one big read, one solver
call) shows a moving bar and the seconds it has taken. A value entered
while such a step runs (the progress window takes the focus) is not lost and
starts nothing on top of it: the step runs again once it is done, with the
latest values.
