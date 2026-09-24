# Batch processing ARPES data with an agent -- design notes

These notes are for anyone who wants a beamtime's worth of maps, kz scans
and cuts processed the way the viewer would process them one by one, but
without sitting at the viewer: by a script, overnight, or by an AI agent that
can run commands and look at pictures.

## 1. Why this is possible without rewriting the viewer

The viewer keeps everything that touches data in `loader/` (reading) and
`tools/` (algorithms), and **neither imports Qt** -- `test/test_imports.py`
checks it. Reading ANTARES, CASSIOPEE and native files, de-gridding, Fermi
edges, the Fermi-surface correction, kz alignment, the k and k_z conversions,
cropping, normalising and saving can all be called from a plain Python
process with no display.

In the viewer, some decisions need a person:

| decision | in the viewer | here |
|---|---|---|
| E_F and the edge's curvature along the slit | points clicked on a cut, or the gold-reference dialog | gold reference fitted channel by channel, matched by lens mode and pass energy |
| Gamma on a map or a cut | picked on the contour | symmetry centre with a score; explicit per-file overrides win |
| the box for kz alignment | dragged on the slit cut | central detector channels, energy window around the summed edge |
| lattice, cleavage plane, V0 | typed into the kz window | **asked of the user**, never guessed (section 4) |
| which photon energies are high-symmetry planes | read off the zone lines | also matched to a band calculation, if the user has one (section 5) |

Previews and per-file overrides let an agent (or a person) check each
result and correct it.

## 2. Structure

```
ARPES-data-browser/
  ARPES_viewer/        unchanged
  ARPES_batch/         new; imports ARPES_viewer's loader/ and tools/ only
    arpes_batch/
      inventory.py     every entry: kind, shape, hv, pass energy, lens mode, T; flags hv scans in map layout
      reference.py     gold reference: E_F per slit channel -> polynomial, E_F(kin), work function
      center.py        Gamma: inversion centre of a map, mirror centre along a slit
      sample.py        lattice, surface normal, V0, calculation: what is missing, prompts, k_z period
      calcbands.py     band files; labels from the lattice; matching a kz scan to a calculation
      steps.py         degrid, calibrate_energy, kz_calibrate, crop, bin, normalise,
                       kconvert, kz_match_calculation, kz_convert, kconvert_cut, curvature
      dataset.py       Qt-free dataset; provenance written exactly as the viewer's MemoryData writes it
      recipe.py        JSON recipe: one chain per kind, sample block, per-file overrides
      runner.py        jobs (maps before cuts), checkpoints and resume, processes, report.md
      preview.py       QC pictures; zone lines on k_z maps; the calculation-match figure
      cli.py           command line; --json for an agent; needs_input with exit status 2
      mcp_server.py    the same operations as MCP tools
    recipes/antares_maps_example.json, recipes/kz_example.json
    CLAUDE.md          the rules and the loop an agent follows
    tests/             synthetic data with planted answers; end-to-end tests
```

## 3. The three chains

**Deflector maps.** These steps are unchanged from the first version:

1. `degrid` -- the map is its own grid reference, and the grid is kept for
   the cuts;
2. `calibrate_energy` -- gold reference: slit curvature and E_F, moved by
   the photon-energy difference;
3. `crop`, then `normalise`;
4. `save "angle"`;
5. `kconvert` (Gamma `"auto"` or explicit);
6. `save "k"`.

**Cuts.**

1. `degrid` with the grid of a map taken on the same lens mode, pass energy
   and detector window.
   - The runner processes maps before cuts so that the grid exists.
   - With no such map, the cut is left as it is. `cut_fallback: notch` opts
     into the viewer's notch filter, which also removes photoemission in the
     grid's k-space regions.
2. `calibrate_energy` (same as for maps).
3. `crop`.
4. `kconvert_cut`: normal emission `"auto"` is the mirror centre along the
   slit; the deflector angle comes from the file.
5. `curvature`, with `a0` from the viewer's `suggest_a0`.

**Photon-energy (kz) scans.** CASSIOPEE folders load as kz maps. ANTARES
writes an hv scan in the deflector-map layout; the inventory flags it
(`first_axis_looks_like: photon_energy`), and listing it under
`photon_energy_scans` loads it the way the viewer's loader window does with
"First axis is photon energy".

1. `degrid`.
2. `kz_calibrate` -- the viewer's `tools.kzmap`. It has three sources:
   - `self` fits each spectrum's own edge over the same detector channels
     and window (metals);
   - `reference` sets E_F,kin = hv − W with W from gold (semiconductors,
     which have no edge; needs a kinetic energy axis);
   - `file` trusts the loader's E − E_F.

   The gold's slit curvature is applied first, since it is an analyser
   property and the same at every photon energy. Each spectrum is then
   normalised by its total.
3. `crop`, then `save "hv"`.
4. `kz_match_calculation`, if there is a calculation (section 5).
5. `kz_convert` -- the viewer's `tools.kzconv.to_kz_cube`, with V0 taken as
   a number, from `"scan"`, or from `"calculation"`.
   - Normal emission `"auto"` is the mirror centre of the scan.
   - W comes from gold, from the scan's own edges, or from the file.
6. `save "kz"`. The preview draws the Gamma planes (solid) and zone
   boundaries (dashed).

**V0 from the data.** V0 `"scan"` uses a period-overlap estimator added
here. The k_z profile at normal emission is converted at each trial V0 and
correlated with itself shifted by exactly one period.
- At the right V0 it repeats. At a wrong V0 the √(E_kin + V0) mapping
  stretches it unevenly, so no single shift lines it up.
- Several binding energies vote, and the spread of the votes is the
  uncertainty reported.
- On synthetic scans covering 2.5 zones it recovered V0 = 8, 14, 20 and
  25 eV to within 1 eV.
- The viewer's `scan_inner_potential` measures the profile's dominant period
  instead. On the same scans it found no crossing, or landed 2-5 eV high
  within its own stated ±3 eV. It is still run and reported for comparison.
- With fewer than about 1.5 zones covered, V0 cannot come from the data,
  and the step asks for a value.

Every result is the viewer's native format, with the `proc.step.N` history
and the parameter prefixes the dialogs use (`degrid.*`, `fscorr.*`,
`fitEF.*`, `kconv.*`, `kz_align.*`, `kz_to_k.*`, `kcut.*`). The viewer's
own checks therefore still apply: a corrected map, for instance, refuses a
second de-grid.

## 4. Asking the user

A kz conversion needs three things the data cannot supply:

- **the lattice** (a, b, c, α, β, γ, space group). The space group fixes
  the centering, and the centering fixes the period: along (001) of a
  body-centred cell the first reciprocal-lattice vector is (002), so the
  period is 4π/c;
- **the surface normal** (the cleavage plane);
- **V0**, or `"scan"`, or `"calculation"`.

It also asks, once, whether there is a band calculation to compare with.

These questions live in `arpes_batch/sample.py` as a table, each with its
reason and an example. `run` checks the recipe before anything is processed:

- **at a terminal**, it asks each question and offers to write the answers
  into the recipe;
- **with `--json`** (an agent) it prints
  `{"needs_input": [...], "message": "... Ask the user; do not guess them."}`
  and exits with status 2. The MCP tools return the same.

The calculation question is optional. Answering it with `null` records that
there is none, so it is not asked again. `python -m arpes_batch surfaces`
lists a lattice's low-index surfaces and their k_z periods, and, given a
measured period, the planes that would produce it -- the viewer's
`tools.cleavage`, for a user unsure which plane they cleaved.

## 5. Matching a kz scan to a band calculation

This answers which photon energies reach the high-symmetry planes, from a
calculation rather than from the lattice alone.

1. **Read the band file.**
   - Formats: two columns per band in blocks (Wannier90 `*_band.dat`,
     Quantum ESPRESSO `bands.out.gnu`, VASPKIT `BAND.dat`), or k plus one
     column per band. Give `fermi_energy_eV` if the energies are absolute.
   - The path's vertices are the repeated k points.
   - The vertices are named, in order of preference, from a labels file
     (Wannier90 `labelinfo`, VASPKIT `KLABELS`), from the `path` given, or
     -- for primitive hexagonal, tetragonal, cubic and orthorhombic lattices
     -- by matching the segment lengths to distances between the lattice's
     high-symmetry points.
   - On the example file supplied with this request this gave
     Γ–M–K–Γ–A–L–H–A with a = 6.20 Å and c = 5.77 Å, to 0.0%.
2. **Find the segment along the surface normal.** This is the one from Γ
   whose length is half the surface period (Γ–A here). Both k conventions
   (with and without 2π) are tried.
3. **Simulate the scan at normal emission.**
   - For a trial V0 and rigid energy shift δ, each (hν, E) maps to
     k_z = √(A (hν − W + E + V0)).
   - That is folded into the zone and located on the Γ–A segment. The
     calculated bands there are drawn as Gaussians, occupied states only.
4. **Compare.** The simulation is correlated with the measured −∂²I/∂E² at
   normal emission, over V0 and δ (and optionally a band renormalisation).
   The best pair gives V0, the calculation's offset, and **the photon
   energies of every Γ and boundary plane inside the scan**.
5. **Check it by eye.** `__calc_match.png` has four panels:
   - the measured normal-emission image with the calculated bands drawn
     through it and the plane photon energies marked;
   - the agreement map over V0 and δ;
   - the cut at the photon energy nearest a Γ plane, with the calculated
     Γ–M / Γ–K bands overlaid;
   - the cut nearest a boundary plane, with A–L / A–H overlaid. Both
     directions are drawn, since the slit's azimuth relative to the crystal
     is not known here.

`kz_convert` uses the matched V0 when the sample says `"calculation"`, and
reports it next to its own estimate otherwise.

On a synthetic hexagonal scan generated from its own "calculation", offset
by +50 meV at V0 = 11 eV, the match returned:
- V0 = 11.0 eV;
- δ = +0.050 eV;
- all six plane photon energies of the planted truth (21.7, 34.1, 48.8,
  65.8, 85.0 and 106.5 eV).

## 6. How an agent uses it

**A. Claude Code on the machine that holds the data (recommended).**

1. Sync the data folder locally, e.g. with Google Drive for Desktop. A
   folder shared with you has to be added to My Drive as a shortcut before
   it syncs.
2. Point `ARPES_RAW` at it.
3. Start `claude` in `ARPES_batch/`, or `claude remote-control` to follow
   from a phone or the web.

The agent reads `CLAUDE.md` and works the loop:

1. `inventory`;
2. `references`;
3. `run --dry-run`, answering `needs_input` by asking the user;
4. `run --only` on one file;
5. check the PNGs;
6. set overrides;
7. run everything;
8. summarise from `report.md`.

**B. Claude Desktop (or any MCP client) + `python -m arpes_batch.mcp_server`.**
The server exposes eight tools: `inventory`, `fit_references`, `plan`,
`run_recipe`, `surfaces`, `suggest_centre`, `preview` and `latest_results`.

**C. A cloud session.** A cloud container cannot see local files, and
connector downloads are size-limited (10 MB per file for the Google Drive
connector at the time of writing). Files have to come through the storage
provider's API with a read-only credential set in the environment. Download
them to a scratch folder, process them, then delete them.

## 7. What would change in `ARPES_viewer/` to fold this in

None of this is needed for `ARPES_batch` to work. But folding it into the
viewer would give one implementation shared by the dialogs and the batch
path, and two items below are bugs in the viewer as it stands.

**Behaviour in the viewer (worth fixing regardless)**

1. **An ANTARES hv scan cannot use the kz tools.**
   - `loader/registry.py:233` (`apply_options`) records
     `axis0.role = "photon_energy"` but leaves the kind as `map`.
   - `ui/windows.py:2965/2974` offer "kz map processing" and
     "kz -> momentum" only for kind `kz_map`.
   - Promote a cube whose first axis is the photon energy to `kz_map` there,
     as `arpes_batch.dataset.Dataset.load(as_kz=True)` does.
2. **The kz window builds a hexagonal lattice with γ = 90°.**
   - `ui/kzconv.py:331` makes `LatticeParams(a, a, c, space_group)` with the
     angles left at 90°, for every space group, the default 194 included.
   - `validate_lattice_parameters` then says "hexagonal requires gamma =
     120 deg".
   - The surface-normal list is built from that wrong reciprocal lattice.
     (001) is unaffected; any other surface is not.
   - Set the angles from the crystal system, or offer all six parameters as
     the batch `sample` block does.

**Move the Qt-free logic out of the GUI modules**

3. Gold reference and channel-by-channel fit -- `ui/windows.py`:
   `reference_frame` (4078), `metadata_temperature` (4099),
   `run_channel_fit` (4995), `AuReferenceDialog` (5045).
   - Move them to `tools/reference.py`; `arpes_batch/reference.py` is a
     first draft of it.
   - `ui/kzmap.py:351` imports `metadata_temperature` from `ui.windows` for
     the same reason.
4. Computed-dataset container and provenance -- `ui/widgets.py`:
   `_MemScan` (871), `MemoryData` (924), `KMapData` (1012).
   - Importing `ui.widgets` imports pyqtgraph (line 45).
   - Move them to `loader/derived.py`; `arpes_batch/dataset.py` then goes
     away.
5. The energy-axis offset (`ui/windows.py:403`), the k-conversion wrapper
   (`ui/windows.py:3434`, `:3495`) and the de-grid settings matching
   (`ui/degrid.py:84`, `:89`, in a module that imports pyqtgraph) belong in
   `tools/`.
6. **kz additions that belong in `tools/`:**
   - the period-overlap V0 estimator (`arpes_batch.steps.v0_by_period_overlap`),
     beside `tools.kzconv.scan_inner_potential`;
   - the calculation match (`arpes_batch/calcbands.py`), as
     `tools/calcbands.py`, with a "Compare with a calculation..." panel in
     the kz window.

**Memory, for large maps**

7. `tools/kspace.py:295` (`convert_map`) and `tools/analysis.py:175`
   (`fs_correction`) promote the whole cube to float64. Offer float32 and
   energy chunking, and let `convert_map`'s progress callback cancel.
8. `tools/degrid.py:500` (`degrid_map`) reads the whole cube into memory by
   design. Run the largest files one process at a time, or `bin` them
   first.

**Metadata (confirm the field names on a real file with `inspect_nxs`)**

9. `loader/nxs_file.py:876` (`aliases`): add the polarisation and the
   manipulator angles (theta/tilt/phi). The latter are the physical starting
   value for Gamma and for `kz_convert`'s `theta_position`.
10. `loader/nxs_file.py:1403` (`list_datasets`) could return the shapes and
    axis ranges, so the inventory does not build an `NxsScan` per entry.
11. **Slit-axis unit.**
    - `README.md:757` says the MBS delivers the slit axis in Å⁻¹, but the
      loader labels it in degrees and every conversion treats it as an
      angle.
    - One real file settles it; the inventory flags it per file.

**Small**

12. `tools/process.py:886` (`record_step`) counts the `.from` keys when
    numbering, so steps are numbered 1, 3, 5, ...
    - `history_of` sorts numerically, so nothing reads them wrongly today.
    - Anything that sorts the keys as text scrambles the order after the
      tenth step.
13. `test/` has no synthetic ANTARES or kz files.
    `ARPES_batch/tests/synthetic.py` could move there.

## 8. What has been checked, and what has not

**Checked, on synthetic data with planted answers (16 tests).**

- ANTARES maps:
  - E_F within 1 meV;
  - curvature 29.7 meV (planted: 30 meV);
  - Gamma within 0.02°;
  - the grid removed.
- ANTARES cut:
  - de-gridded with the map's grid;
  - the slit centre within 0.2°.
- ANTARES hv scan, calibrated from gold:
  - E_F = hν − W;
  - W within 5 meV.
- Body-centred kz scan:
  - per-spectrum E_F drift recovered to about 6 meV rms;
  - (002) period used;
  - V0 = 14.5 ± 0.5 eV (planted: 14).
- Hexagonal kz scan against its calculation: V0, δ and the six plane photon
  energies exact.
- A kz recipe without the sample block asks and processes nothing.
- The raw folder is byte-for-byte unchanged.
- Re-runs resume from checkpoints, and results open in the viewer.

**Not yet checked: real files.** On a first real beamtime, look at:

1. each entry's kind in the inventory (navigation SPEM images? hv scans
   in map layout?), its axis lengths, and any "ambiguous axis" warnings;
2. the slit-axis unit (item 11);
3. whether the gold references share the maps' lens mode, pass energy and
   photon energy;
4. how well the real MCP grid is removed;
5. memory on the largest files;
6. for a real calculation:
   - whether its lattice matches the sample's (a relaxed DFT cell differs
     by a few per cent; the labels are inferred with a 3 % tolerance);
   - whether a rigid shift is enough, or a renormalisation (`renormalisations`,
     e.g. `[1.0, 1.3, 1.6]`) is needed.
