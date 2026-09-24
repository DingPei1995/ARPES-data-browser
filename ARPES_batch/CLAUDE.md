# Working with ARPES_batch (instructions for an agent)

This folder drives the ARPES viewer's algorithms (`../ARPES_viewer/loader`,
`../ARPES_viewer/tools`) without its windows. Run everything from this
folder: `python -m arpes_batch ... --json`.

## Rules

- **Never write inside a raw-data folder.** Every command checks this and
  refuses; do not work around it (no symlinks, no copying outputs back).
  Outputs go to the recipe's `output`, which must be outside `inputs` and
  `references`.
- **Do not edit `../ARPES_viewer/`.** If a step needs something the viewer
  does not offer, add a step in `arpes_batch/steps.py` and a test.
- Raw files are opened read-only; never rename, move or delete them.
- A map is typically hundreds of MB on disk and several times that in
  float64 while it is processed. Keep `workers` at 1 for the largest files
  (roughly 500 MB and up); 2-3 for smaller ones on a machine with 16 GB or more.

## The loop

1. `python -m arpes_batch inventory <raw folder> --out <output>/inventory --json`
   Read the entries: kind, photon energy, pass energy, lens mode,
   temperature, `slit_axis_looks_like`, `notes`. Group the maps by the
   conditions they were taken at (a photon-energy series, polarisation
   pairs, positions on the sample).
   If any row has `slit_axis_looks_like: momentum`, stop and ask the user.
2. `python -m arpes_batch references "<Au ref folder>" --out <output> --json`
   Open `<output>/references/*.png`. A good fit: smooth points along the
   curve, rms of a few meV or less, most channels used. Report E_F, the
   work function (hv - E_F) and the resolution of each reference.
3. Copy `recipes/antares_maps_example.json`, set `inputs`, `references`,
   `output`, then `python -m arpes_batch run <recipe> --dry-run --json`.
   Every job should have a reference; if a job says "no reference with lens
   mode ...", tell the user rather than switching calibration off.
4. `python -m arpes_batch run <recipe> --only "<one file>" --json` first, open
   its `__angle.png` and `__k.png`, then run the rest.
5. For each result, check on the PNGs and in `qc`:
   - the Fermi cut-off in the cuts is flat and at 0 (`calibrate_energy`);
   - `degrid` grid power went down, and is near 1 or below after;
   - the cyan cross on `__angle.png` sits on Gamma. `centre_suggestion.score`
     below ~0.8, or a cross that is visibly off, means: put explicit
     `theta_offset_deg` / `phi_offset_deg` (deflector, slit angles of Gamma,
     read off the preview axes) into `overrides` for that file and re-run it.
     If the map does not contain Gamma at all, say so; do not invent one.
   - `azimuth_deg` so that a high-symmetry direction is along kx, if the
     user wants that; it cannot be inferred reliably without the lattice.
6. Re-running is cheap: the chain resumes from the last `save` whose inputs
   are unchanged, so changing only `kconvert` does not repeat the de-grid.
7. Summarise for the user from `<output>/runs/<latest>/report.md`: what ran,
   what was skipped or stopped and why, and which files need their eyes.
