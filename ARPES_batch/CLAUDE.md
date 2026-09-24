# Working with ARPES_batch (instructions for an agent)

This folder drives the ARPES viewer's algorithms (`../ARPES_viewer/loader`,
`../ARPES_viewer/tools`) without its windows, for deflector maps,
photon-energy (kz) scans and single cuts. Run everything from this folder:
`python -m arpes_batch ... --json`.

## Rules

- **Never write inside a raw-data folder.** Every command checks this and
  refuses; do not work around it (no symlinks, no copying outputs back).
  Outputs go to the recipe's `output`, which must be outside `inputs` and
  `references`.
- **Do not edit `../ARPES_viewer/`.** If a step needs something the viewer
  does not offer, add a step in `arpes_batch/steps.py` and a test.
- Raw files are opened read-only; never rename, move or delete them.
- **Never invent sample information.** The lattice (a, b, c, angles, space
  group), the surface normal (cleavage plane) and the inner potential V0
  come from the user. When a command answers `needs_input`, put each
  question to the user as written, with its example, and write the answers
  into the recipe's `sample` block. The band-calculation question is
  optional: if the user has none, write `"calculation": null` so it is not
  asked again.
- A map is typically hundreds of MB on disk and several times that in
  float64 while it is processed. Keep `workers` at 1 for the largest files
  (roughly 500 MB and up); 2-3 for smaller ones on a machine with 16 GB or more.

## The loop

1. `python -m arpes_batch inventory <raw folder> --out <output>/inventory --json`
   Read the entries: kind (map, cut, kz_map), photon energy, pass energy,
   lens mode, temperature, `slit_axis_looks_like`, `first_axis_looks_like`,
   `notes`.
   - Group the entries by the conditions they were taken at (a photon-energy
     series, polarisation pairs, positions on the sample).
   - If any row has `slit_axis_looks_like: momentum`, stop and ask the user.
   - A map with `first_axis_looks_like: photon_energy` is probably an hv scan
     in ANTARES's map layout. Confirm with the user, then list it under
     `photon_energy_scans`.
2. `python -m arpes_batch references "<Au ref folder>" --out <output> --json`
   Open `<output>/references/*.png`. A good fit has smooth points along the
   curve, an rms of a few meV or less, and most channels used. Report E_F,
   the work function (hv - E_F) and the resolution of each reference.
3. Copy `recipes/antares_maps_example.json` (maps and cuts) or
   `recipes/kz_example.json` (kz scans), and set `inputs`, `references` and
   `output`. Then run `python -m arpes_batch run <recipe> --dry-run --json`.
   - Every job should have a reference. If a job says "no reference with
     lens mode ...", tell the user rather than switching calibration off.
   - If `needs_input` is present, ask the user now (see Rules).
     `python -m arpes_batch surfaces --a .. --c .. --space-group ..` lists
     the candidate cleavage planes and their k_z periods, to help them answer.
4. First run a single file:
   `python -m arpes_batch run <recipe> --only "<one file>" --json`. Open its
   PNGs, then run the rest.
5. For each result, check the PNGs and `qc`:
   - **maps and cuts**
     - the Fermi cut-off is flat and at 0 (`calibrate_energy`);
     - the `degrid` grid power went down;
     - the cyan cross on `__angle.png` sits on Gamma. If
       `centre_suggestion.score` is below about 0.8, or the cross is visibly
       off, put explicit offsets (deflector and slit angles, read off the
       preview axes) into `overrides` for that file and re-run it. If the
       map does not contain Gamma, say so; do not invent one.
   - **kz maps**
     - `__kz_qc.png`: E_F against photon energy should be a smooth curve
       (`source: self`). A semiconductor has no edge of its own; use
       `source: reference`.
     - `v0_scan` gives V0 with its uncertainty; `zones_covered` below about
       1.5 means V0 cannot come from the data, so ask the user for a value.
     - `__kz.png` shows the Gamma planes (solid) and zone boundaries
       (dashed); the pattern should repeat between them.
     - With a calculation, `__calc_match.png` lists which photon energies
       reach the Gamma and boundary planes, and shows the cuts there with
       the calculated in-plane bands. Report those photon energies, the
       matched V0 and the calculation's energy shift to the user.
6. Re-running is cheap: the chain resumes from the last `save` whose inputs
   are unchanged.
7. Summarise for the user from `<output>/runs/<latest>/report.md`: what ran,
   what was skipped or stopped and why, and which files need their eyes.
