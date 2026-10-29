# For developers: the layout of the code, and the tests

Where everything is and how it is tested. The technical references sit beside the code: [loader](../loader/README.md) · [ui](../ui/README.md) · [tools](../tools/README.md).

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

## Tests

```
python -m pytest
```

from this folder. 859 of them, no display needed -- `conftest.py` puts the
project root on `sys.path` and pins `QT_QPA_PLATFORM=offscreen`, and
`test/conftest.py` runs garbage collection only between tests, the rule
`ui/common/gcguard.py` applies in the program.
