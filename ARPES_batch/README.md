# ARPES_batch

Batch processing of ARPES maps with the ARPES viewer's own algorithms, from
a recipe, a command line or an MCP server -- so that an agent (or a script,
or a person at 2 a.m.) can process a beamtime without the windows.

It **adds** a folder beside `ARPES_viewer/` and changes nothing in it: the
viewer's `loader/` and `tools/` are Qt-free by design, and this package
imports them as they are. What it writes is the viewer's own saved format,
with the processing history in it, so every result opens in the viewer.

```
python -m arpes_batch inventory  RAW_FOLDER --out OUT --json
python -m arpes_batch references "RAW/Au ref" --out OUT --json
python -m arpes_batch run recipes/antares_maps_example.json --dry-run --json
python -m arpes_batch run recipes/antares_maps_example.json [--only "map_*"] [--workers 2]
python -m arpes_batch centre FILE.nxs            # where Gamma probably is
python -m arpes_batch preview FILE.nxs --png OUT/quicklook.png
python -m arpes_batch.mcp_server                 # the same, as MCP tools
```

Needs what the viewer's non-Qt half needs (numpy, scipy, h5py) plus
matplotlib for the previews; `mcp` only for the MCP server. No PyQt.

Why it is built this way, how an agent uses it, and what would change in
`ARPES_viewer/` to fold it in: [DESIGN.md](DESIGN.md). The working rules an
agent follows here: [CLAUDE.md](CLAUDE.md).

## What a run does, per map

| step | viewer equivalent | `tools/` function |
|---|---|---|
| `degrid` | Functions → De-grid map... | `tools.degrid.degrid_map` |
| `calibrate_energy` | Fermi surface from a reference, FS correction, Fermi level → offset | `tools.fermi.fit_channels`, `tools.analysis.fs_correction` |
| `crop` / `bin` / `normalise` | Data operations | `tools.dataops.truncate` / `compress` / `self_normalize` |
| `kconvert` | Functions → Map k conversion | `tools.kspace.convert_map` |
| `save` | Save... | `loader.nxs_file.save_dataset` |

The two decisions the dialogs leave to a person are made here as follows:

- **E_F and the slit curvature** come from the gold reference with the same
  lens mode and pass energy, nearest in photon energy, then in time. A map
  at another photon energy is placed at `E_F(hv) = E_F(hv_ref) + (hv - hv_ref)`
  (recorded on the result). A semiconductor has no edge of its own, which is
  why this is the default rather than fitting the map. With no matching
  reference the job **stops** (`"required": true`) rather than writing an
  uncalibrated map that looks calibrated.
- **Gamma** is suggested as the inversion-symmetry centre of a
  constant-energy map, with a score; the angle preview marks it with a cyan
  cross. Explicit per-file offsets in the recipe's `overrides` win.

## Safety

- The output folder is refused if it is, lies inside, or contains any input
  or reference folder (symlinks resolved).
- Files are opened read-only through the viewer's loader; results are
  written atomically (`.partial`, then renamed).
- Every result carries `proc.step.N` history and the dialog's own parameter
  keys (`degrid.*`, `fscorr.*`, `fitEF.*`, `kconv.*`), so the viewer's
  checks still apply -- e.g. it refuses to de-grid a corrected map.

## Output layout

```
OUT/
  manifest.json, manifest.csv          every entry and its conditions
  references/<ref>.json, <ref>.png     the fitted edge per reference
  <folder>/<file>__<entry>/
     ...__angle.nxs  ...__angle.png    de-gridded, calibrated, in angles
     ...__k.nxs      ...__k.png        in momentum
     ...__result.json                  per-step status, timings, QC numbers
  runs/<time>/recipe.json, results.json, report.md
  latest_results.json
```

## Tests

```
python -m pytest tests
```

Synthetic ANTARES maps (`tests/synthetic.py`) with a planted E_F, slit
curvature, Gamma and detector grid; the tests check each is recovered, that
the raw folder is byte-for-byte untouched, that a re-run resumes from its
checkpoint, and that the results open as the viewer's own datasets.
