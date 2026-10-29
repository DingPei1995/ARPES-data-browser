# Spectra: technical reference

Core levels, curves, fits and spin: `corelevel.py`, `curves.py`, `peaks.py`, `spin.py` and their neighbours in this folder.

## Core-level check

Right-click a cut in the list → **Core-level check...**. Some "cuts" are
really a wide core-level scan, and the first question about one is what
the peaks are. Three steps, each correctable before the next:

1. **Elements.** Type what the sample is made of: symbols (`Nb, S, Cl`) or
   a formula (`NbSCl`, `Bi2Se3`). The match tolerance (0.5 eV) and the
   "uncertain" window (2 eV) are set here too; the last values are kept in
   `~/.arpes_viewer/corelevel.json`.
2. **Energy reference.** The whole cut is summed over angle and the Fermi
   edge looked for where it should be (hν − W on a kinetic axis, zero on
   `E - E_F`). If the window holds it, it is fitted -- the nominal photon
   energy is often off by more than the tolerance -- and hν − E_F is shown
   as an implied work function, which should be 4-5 eV. If not (the usual
   core-level window), E_F is put at hν − W. The axis type, hν, W and E_F
   are all editable, and the E_F line can be dragged.
3. **Results.** The EDC on a binding-energy axis (deeper to the left), each
   peak marked and listed:
   - **sample** -- within the tolerance of a level of the sample's
     elements, with the orbital named (and its spin-orbit partner, when
     that peak is there too);
   - **sample (orbital ambiguous)** -- two levels of the sample fit;
   - **uncertain** -- no level within the tolerance but one within the
     uncertain window: a chemical shift, charging, or not the sample; the
     nearest levels are listed and the call is yours;
   - **not the sample** -- the lines of *common* elements near it (neither
     radioactive nor highly toxic, the usual UHV contaminants first: C, O,
     N, F, Cl, Na, K, S, Si, Ca, Ar, Ta, Mo, W, Cu, Ag, In, Au, ...) with
     where each typically comes from, plus second- and third-order light
     and Auger lines.

   The assignment column is editable (and kept through a re-analysis);
   **Copy report** includes what was changed by hand. **Add EDC to list**
   lists the integrated EDC on the binding-energy scale.

   **Weak lines next to a strong one.** *Log counts* puts the count axis on
   a log scale. *Min height* (default 0.2 % of the strongest peak) is the
   only cut besides the noise threshold; it used to be a fixed 1 %, which
   dropped real weak lines beside a strong one. A peak the search still
   misses -- a shoulder with no maximum of its own -- can be added by hand:
   press *Click to add peak* and click the plot (a double-click always
   adds), or type its binding energy and press *Add*. Its position is
   refined by a local Gaussian fit, it is marked ◆, assigned like any
   other peak, and named with every level within ±2 eV -- the sample's
   first, then all other elements. *Remove selected* drops a peak (a found
   one stays out until *Reset peaks*).

   **Peak table...** lists every tabulated level (Table 1-1, plus the
   typical chemical states), filtered by element (symbols or a formula;
   *Sample* / *All*) and by binding-energy range (*Whole window*, *Plot
   view*, or typed), sortable by any column, the sample's levels in bold.
   Picking a row marks that level on the plots (magenta dashed line);
   *Add as peak* or a double-click adds a peak there -- refined by the
   local fit -- with the level written in as its assignment.

   **Select range** shows a band on the plot to drag over a range;
   **Zoom window...** opens that range (or, with no band, what the plot
   shows) in a window of its own: from/to editable, its own log axis, a
   smoothed curve over the data, and below it the negative second
   derivative (divided by the intensity, *÷ I*, so a weak peak shows as
   clearly as a strong one) -- every peak, a shoulder with no maximum of
   its own included, is a maximum there. A cursor (drag it, or click)
   reads out its energy; *Fit at cursor* fits a Gaussian there and moves
   the cursor to its centre; *Add peak at cursor* (or a double-click) adds
   it to the results. The zoom window follows the results window's peaks
   and the level picked in the peak table.

The table is Table 1-1 of the X-ray Data Booklet (`Table_1-1.pdf` at the
top of the repository), turned once into `tools/spectra/corelevel_table.py` by
`devtools/gen_corelevels.py` -- the PDF is never read at run time. What
the check adds beyond a lookup:

- **Gas-phase values.** The table gives H, N, O, F, Cl and the rare gases
  against the vacuum level of the free gas (O 1s 543.1 eV, where an oxide
  is at ~530). Typical solid-state values (oxides, hydroxides, chlorides,
  sulfides, sulfates, adventitious carbon and its oxidised forms, PTFE,
  SiO2, implanted Ar, the common metal oxides) are compared as well; the
  checkbox turns them off.
- **Spin-orbit doublets** confirm an assignment when both components are
  found at the tabulated splitting.
- **Higher-order light.** A level deeper than hν − W can appear at E_B − hν
  (2hν) or E_B − 2hν (3hν).
- **Auger lines** (C, N, O, F, Na KLL) are named; they sit at a fixed
  kinetic energy, so measuring at a second photon energy tells them apart:
  a core level stays put in binding energy, an Auger line and higher-order
  light move.
- **A rigid shift.** If shifting every peak by the same amount would make
  clearly more of them fit (charging of an insulating sample, or a wrong
  E_F or hν), the shift is reported; **Apply suggested shift** applies it.
- **What is missing.** Levels of the sample inside the window with no peak,
  and elements with no level in the window at all (with the nearest levels
  reachable at this photon energy).

Peaks are maxima of the lightly smoothed EDC that stand 8 noise sigmas
(peak to valley) above their surroundings, the noise measured from the
curve itself and scaled as √counts for counting data; each is refined by a
Gaussian on a sloping background. **Peak threshold** lowers it for weak
lines.

## The curve viewer

Every channel of a curve dataset, with error bars where it has a σ channel.
Channels can be hidden, offset into a waterfall (**Offset**), and read off
with the cursor. **Range** is a draggable x window used by the operations
and the panels below.

**Operations** (**Functions → Transform**) each make a new dataset in
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

## Curve fit

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

## Spin analysis

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
