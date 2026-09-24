# ARPES_batch

Batch processing of ARPES data with the ARPES viewer's own algorithms --
deflector maps, photon-energy (kz) scans and single cuts -- from a recipe, a
command line or an MCP server. It lets an agent, a script or a person process
a whole beamtime without the windows.

It **adds** a folder beside `ARPES_viewer/` and changes nothing in it. The
viewer's `loader/` and `tools/` are Qt-free by design, and this package
imports them as they are. What it writes is the viewer's own saved format,
with the processing history in it, so every result opens in the viewer.

```
python -m arpes_batch inventory  RAW_FOLDER --out OUT --json
python -m arpes_batch references "RAW/Au ref" --out OUT --json
python -m arpes_batch surfaces   --a 3.96 --c 13.02 --space-group 139 [--normal 0 0 1]
python -m arpes_batch init       --inputs RAW --references "RAW/Au ref" --output OUT [--kz] > recipe.json
python -m arpes_batch run        recipe.json --dry-run --json
python -m arpes_batch run        recipe.json [--only "map_*"] [--workers 2]
python -m arpes_batch centre     FILE.nxs            # where Gamma probably is
python -m arpes_batch preview    FILE.nxs --png OUT/quicklook.png
python -m arpes_batch.mcp_server                     # the same, as MCP tools
```

It needs what the viewer's non-Qt half needs (numpy, scipy, h5py), plus
matplotlib for the previews. The MCP server also needs `mcp`. PyQt is not
needed.

[DESIGN.md](DESIGN.md) explains why it is built this way, how an agent uses
it, and what would change in `ARPES_viewer/` to fold it in. The working rules
an agent follows here are in [CLAUDE.md](CLAUDE.md).

## The steps

A recipe holds one chain of steps per kind of data, under
`steps.map`, `steps.kz_map` and `steps.cut`.

| step | kinds | viewer equivalent | function |
|---|---|---|---|
| `degrid` | map, kz map, cut | De-grid | `tools.degrid` -- a map's grid is kept and used for the cuts taken with the same detector settings |
| `calibrate_energy` | map, cut | gold reference, FS correction, Fermi level → offset | `tools.fermi.fit_channels`, `tools.analysis.fs_correction` |
| `kz_calibrate` | kz map | kz map processing | `tools.kzmap` -- each spectrum's own edge (`self`), `hv - W` from gold (`reference`), or as loaded (`file`) |
| `crop` / `bin` / `normalise` | all | Data operations | `tools.dataops` |
| `kconvert` | map | Map k conversion | `tools.kspace.convert_map` |
| `kz_match_calculation` | kz map | *(new)* | `arpes_batch.calcbands` -- match to a band calculation |
| `kz_convert` | kz map | kz → momentum | `tools.kzconv.to_kz_cube` |
| `kconvert_cut` | cut | Cut k conversion | `tools.cutk.convert_cut` |
| `curvature` | cut | Process → curvature | `tools.process.curvature` |
| `save` | all | Save | `loader.nxs_file.save_dataset` (a checkpoint, and a preview PNG) |

## What only the user knows

The viewer's dialogs leave several decisions to a person. Here they are made
as follows.

- **E_F and the slit curvature.** They come from the gold reference with the
  same lens mode and pass energy, nearest in photon energy, then nearest in
  time. A map at another photon energy is placed at
  `E_F(hv) = E_F(hv_ref) + (hv − hv_ref)`, and this is recorded on the
  result. A semiconductor has no edge of its own, which is why this is the
  default. With no matching reference, the job **stops** rather than
  writing an uncalibrated map that looks calibrated.
- **Gamma.** It is suggested as the inversion-symmetry centre of a
  constant-energy map (for a cut, the mirror centre along the slit), with a
  score. The preview marks it. Explicit per-file offsets in `overrides` win.
- **The sample, for a kz conversion.** The lattice and space group (the
  centering sets the k_z period: body-centred (001) repeats every 4π/c),
  the surface normal, and V0 (a number, `"scan"`, or `"calculation"`).
  These are never guessed. Until they are in the recipe's `sample` block,
  `run` refuses before processing anything:
  - at a terminal, it asks for them and offers to save the answers;
  - with `--json`, it returns `needs_input` with exit status 2, for an
    agent to put the questions to the user.
- **A band calculation.** The user is also asked whether one exists; answer
  with `null` if not. If it does, `kz_match_calculation` compares the scan
  at normal emission with the calculated dispersion along Gamma → zone
  boundary on the surface normal. It gives:
  - V0 and the calculation's rigid energy offset;
  - **the photon energies at which the Gamma and boundary planes are
    reached**;
  - a figure with the cuts at those photon energies, the calculated in-plane
    bands drawn over them.

  Band files are read as two columns per band (Wannier90, Quantum ESPRESSO,
  VASPKIT) or as k plus one column per band. The high-symmetry points are
  taken from a labels file, from the path given, or -- for primitive
  lattices -- from the segment lengths themselves.

## Safety

- The output folder is refused if it is an input or reference folder, lies
  inside one, or contains one (symlinks are resolved).
- Files are opened read-only through the viewer's loader.
- Results are written atomically (to `.partial`, then renamed).
- Every result carries the `proc.step.N` history and the dialogs' own
  parameter keys (`degrid.*`, `fscorr.*`, `fitEF.*`, `kconv.*`, `kz_align.*`,
  `kz_to_k.*`, `kcut.*`). The viewer's checks therefore still apply; for
  example, it refuses to de-grid a corrected map.

## Output layout

```
OUT/
  manifest.json, manifest.csv          every entry and its conditions
  references/<ref>.json, <ref>.png     the fitted edge per reference
  grids/<lens>__<PE>__<frame>.npz      detector grids from maps, for the cuts
  <folder>/<file>__<entry>/
     ...__angle.nxs/.png, ...__k.nxs/.png           maps and cuts
     ...__hv.nxs/.png, ...__kz.nxs/.png              kz maps
     ...__kz_qc.png, ...__calc_match.png             E_F vs hv, V0; calculation match
     ...__result.json                                per-step status, timings, QC numbers
  runs/<time>/recipe.json, results.json, report.md
  latest_results.json
```

## Tests

```
python -m pytest tests
```

The tests use synthetic data with planted answers (`tests/synthetic.py`):
ANTARES maps and cuts; a CASSIOPEE-style kz scan of a body-centred metal
and one of a semiconductor on a kinetic axis; and a hexagonal crystal
together with its "calculated" band file. The planted quantities are E_F, the slit curvature, Gamma, the
detector grid, the per-spectrum E_F drift, V0 and the calculation's energy
offset.

The 16 tests check that:
- each planted quantity is recovered;
- the raw folder is byte-for-byte untouched;
- a re-run resumes from its checkpoint;
- a kz recipe without the sample block asks before running;
- the results open as the viewer's own datasets.
