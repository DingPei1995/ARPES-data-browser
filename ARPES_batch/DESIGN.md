# Batch processing ARPES maps with an agent -- design notes

For anyone who wants a beamtime's worth of maps processed the way the
viewer would process them one by one, without sitting at the viewer: by a
script, overnight, or by an AI agent that can run commands and look at
pictures.

## 1. Why this is possible without rewriting the viewer

The viewer keeps everything that touches data in `loader/` (reading) and
`tools/` (algorithms), and **neither imports Qt** -- `test/test_imports.py`
checks it. Reading an ANTARES `.nxs`, de-gridding, fitting a Fermi edge,
the Fermi-surface correction, the k conversion, cropping, normalising and
saving in the program's own format can all be called from a plain Python
process with no display. (The viewer's own non-Qt tests pass headless, with
only numpy, scipy and h5py installed.)

In the viewer, only two decisions need a person:

1. **E_F and the curvature of the edge along the slit** -- points clicked
   on a cut, or a reference chosen in the gold-reference dialog;
2. **where Gamma is** -- the k-space origin, picked on the contour.

A batch run has to make both without a person: E_F comes from a gold
reference fitted channel by channel, and Gamma comes from a symmetry
estimate reported with a confidence score. Previews and per-file
overrides let an agent (or a person) check each result and correct it.

## 2. Structure

```
ARPES-data-browser/
  ARPES_viewer/        unchanged
  ARPES_batch/         new; imports ARPES_viewer's loader/ and tools/ only
    arpes_batch/
      inventory.py     every entry: kind, shape, hv, pass energy, lens mode, T -> manifest.csv/json
      reference.py     gold reference: E_F per slit channel -> polynomial (curvature), E_F(kin), work function
      center.py        Gamma suggestion: inversion-symmetry centre of a constant-energy map, score 0..1
      steps.py         degrid / calibrate_energy / crop / bin / normalise / kconvert
      dataset.py       Qt-free dataset; provenance written exactly as the viewer's MemoryData writes it
      recipe.py        JSON recipe + per-file overrides (Gamma, azimuth, energy window, ...)
      runner.py        one job per entry, checkpoints and resume, processes, failure isolation, report.md
      preview.py       QC pictures (constant-energy maps + two cuts; the reference fit)
      cli.py           command line; --json for an agent to parse
      mcp_server.py    the same operations as MCP tools
    recipes/antares_maps_example.json
    CLAUDE.md          the rules and the loop an agent follows
    tests/             synthetic ANTARES maps with planted answers; end-to-end tests
```

### The chain for one map (the example recipe)

| step | viewer equivalent | function | note |
|---|---|---|---|
| `degrid` | Functions → De-grid map | `tools.degrid.degrid_map` | first: only data still on the detector pixels can be de-gridded |
| `calibrate_energy` | gold reference + FS correction + Fermi level offset | `tools.fermi.fit_channels`, `tools.analysis.fs_correction` | energy becomes E − E_F; the edge is flat along the slit |
| `crop` | Truncate | `tools.dataops.truncate` | e.g. E − E_F in [−3, 0.3] eV; also cuts memory |
| `normalise` | Self-normalise | `tools.dataops.self_normalize` | each deflector slice by its own total (beam decay) |
| `save "angle"` | Save | `loader.nxs_file.save_dataset` | checkpoint and preview |
| `kconvert` | Map k conversion | `tools.kspace.convert_map` | Gamma `"auto"` or explicit; E_kin = (E − E_F) + E_F,kin |
| `save "k"` | Save | `loader.nxs_file.save_dataset` | opens in the viewer as a k-map |

Every result is the viewer's native format, with the `proc.step.N` history
and the parameter prefixes the dialogs use (`degrid.*`, `fscorr.*`,
`fitEF.*`, `kconv.*`). The viewer's own checks therefore still apply: a
corrected map, for instance, refuses a second de-grid.

### Physics decisions

- **E_F comes from the reference, not from the map.** A semiconductor has no
  edge to fit. The reference is the one with the same lens mode and pass
  energy (the curvature depends on both), then the nearest photon energy,
  then the nearest in time.
- **Across photon energies**, `E_F,kin(hv) = E_F,kin(hv_ref) + (hv − hv_ref)`,
  i.e. the analyser work function is taken as constant. This is only as
  good as the monochromator's nominal energy, and the shift is recorded on
  the result. For a photon-energy series, a reference every few photon
  energies is worth the beamtime.
- **No matching reference stops the job** (`"required": true`) instead of
  writing an uncalibrated map that looks calibrated.
- **Gamma "auto" is a suggestion.** Below `min_centre_score` the k
  conversion is skipped; the offsets then go into `overrides`, read off the
  preview. A map that does not contain Gamma needs them given explicitly.

## 3. How an agent uses it

**A. Claude Code on the machine that holds the data (recommended).**
Sync the data folder locally (e.g. Google Drive for Desktop; a folder shared
with you has to be added to My Drive as a shortcut before it syncs), point
`ARPES_RAW` at it, and start `claude` in `ARPES_batch/` (or
`claude remote-control` to follow from a phone or the web). The agent reads
`CLAUDE.md` and works the loop:

1. `inventory`;
2. `references`, checking the fits;
3. `run --dry-run`, checking every map has a reference;
4. `run --only` on one file;
5. look at the PNGs;
6. overrides, then run everything (resumes from checkpoints);
7. summarise from `report.md`.

The agent can read the QC images directly, the files are read locally at
full speed, and the output lives outside the synced raw folder.

**B. Claude Desktop (or any MCP client) + `python -m arpes_batch.mcp_server`.**
The server exposes seven tools: `inventory`, `fit_references`, `plan`,
`run_recipe`, `suggest_centre`, `preview` and `latest_results`. This suits
anyone who prefers not to use a terminal.

**C. A cloud session.** A cloud container cannot see local files, and
connector downloads are size-limited (10 MB per file for the Google Drive
connector at the time of writing). The files have to come through the
storage provider's API instead, with a read-only credential set in the
environment. Download them to a scratch folder, process them, then delete
them.

## 4. What would change in `ARPES_viewer/` to fold this in

None of this is needed for `ARPES_batch` to work: it re-implements the
little glue it needs. But folding it into the viewer would mean one
implementation shared by the dialogs and the batch path. In order of value:

**Move the Qt-free logic out of the GUI modules**

1. Gold reference and channel-by-channel fit -- `ui/windows.py`:
   `reference_frame` (4078), `metadata_temperature` (4099),
   `run_channel_fit` (4995), `AuReferenceDialog` (5045).
   - Move the sum → fit → robust polynomial → curve part to a new
     `tools/reference.py`. `arpes_batch/reference.py` is a first draft of it.
   - Leave only the widgets in the dialog.
2. Computed-dataset container and provenance -- `ui/widgets.py`:
   `_MemScan` (871), `MemoryData` (924), `KMapData` (1012).
   - Importing `ui.widgets` imports pyqtgraph (line 45), so a script cannot
     use them.
   - Move them to `loader/derived.py` and let the GUI subclass them.
     `arpes_batch/dataset.py` then goes away.
3. Energy-axis offset -- `ui/windows.py:403 offset_energy_axis`: the shift,
   and the relabelling to `E - E_F`, belong in `tools/`.
4. k-conversion wrapper -- `ui/windows.py:3434 run_k_conversion`,
   `:3495 _finish_k_conversion`: the settings and the output labels as a
   `tools/kspace.convert_scan()`.
5. De-grid settings matching -- `ui/degrid.py:84 _SETTINGS_KEYS`,
   `:89 settings_match` live in a module that imports pyqtgraph (line 30);
   move them to `tools/degrid.py`.

**Memory, for large maps**

6. `tools/kspace.py:295` (`convert_map`) promotes the whole cube to float64,
   and its progress callback cannot cancel. Offer float32 and energy
   chunking, and allow cancelling.
7. `tools/analysis.py:175` (`fs_correction`) also promotes to float64 and
   lengthens the energy axis. Offer a dtype.
8. `tools/degrid.py:500` (`degrid_map`) reads the whole cube into memory by
   design, since it reads every slice three times.
   - A map of several hundred MB can peak at several GB.
   - Run such files one process at a time, or `bin` them first.

**Metadata (confirm the field names on a real file with `inspect_nxs`)**

9. `loader/nxs_file.py:876` (`aliases`) covers photon energy, grating, exit
   slit and resolution. Worth adding:
   - **polarisation** (from the HU60/HU256 groups);
   - **manipulator angles** (theta/tilt/phi), both to group maps by
     conditions and as a physical starting value for Gamma.
10. `loader/nxs_file.py:1403` (`list_datasets`) could return the shape and
    axis ranges, which the HDF5 shapes give cheaply. The inventory would
    then not need to build an `NxsScan` per entry.
11. **Slit-axis unit.** `README.md:757` says the MBS delivers the slit axis
    in Å⁻¹, but the loader labels it "Angle along slit (°)" and
    `convert_map` treats it as an angle.
    - One real file settles it: roughly ±15 is degrees, roughly ±1 is Å⁻¹.
    - The inventory flags it per file (`slit_axis_looks_like`).

**Small**

12. `tools/process.py:886` (`record_step`) counts the `proc.step.N.from`
    keys when numbering, so steps are numbered 1, 3, 5, ...
    - `history_of` sorts numerically, so nothing reads them wrongly today.
    - Anything that sorts the keys as text scrambles the order after the
      tenth step.
    - Count only the keys without `.from`.
13. `test/` has no synthetic ANTARES file. `ARPES_batch/tests/synthetic.py`
    could move there.
14. Optional: an "Open batch results" entry in the launcher that lists a
    run's outputs. They are native `.nxs` files, so Load data already
    opens them.

## 5. What has been checked, and what has not

**Checked, on synthetic ANTARES maps.** The files use the layout the
loader reads as its "case 4", with a planted E_F, slit curvature, Gamma, a
hexagonal detector grid and Poisson counts.

- Each planted quantity is recovered:
  - E_F within 1 meV;
  - curvature 29.7 meV (planted: 30 meV);
  - Gamma within 0.02°;
  - the grid power goes down.
- E_F moves correctly for a map taken at another photon energy.
- A map with no matching reference stops.
- The raw folder is byte-for-byte and mtime-for-mtime unchanged.
- A re-run resumes from its checkpoint (changing only the k conversion
  redoes only the k conversion).
- The outputs open in the viewer as k-maps with their history.
- Three worker processes run correctly, and the seven MCP tools respond.

**Not yet checked: real files.** On a first real beamtime, look at:

1. each entry's kind in the inventory (navigation SPEM images mixed in?),
   its axis lengths, and any "ambiguous axis" warnings;
2. the slit-axis unit (item 11);
3. whether the gold references share the maps' lens mode, pass energy and
   photon energy;
4. how well the real MCP grid is removed (the viewer's README reports
   60 → 1.4 on an ANTARES WSe₂ map);
5. memory on the largest files.
