# Processing: technical reference

De-gridding and cut arithmetic: `degrid.py`, `cutops.py` and their neighbours in this folder.

## De-grid

Removes the periodic pattern the detector prints on every image -- the
hexagonal MCP structure (ANTARES, the CASSIOPEE MBS end station; ~5 px) or
the square mesh (CASSIOPEE Scienta; ~11 px, 10 % contrast). **Do it first**:
only data still on the detector's pixels can be de-gridded, so the button
refuses anything k-converted, Fermi-surface corrected, kz-aligned,
interpolated along a path, smoothed or differentiated, and says why.

**A map or kz map** -- **Functions → Remove detector grid...** on the contour window
or on its slit-cut window; the whole map is done at once. No reference measurement is
needed: the map is its own reference. The grid stays on the same pixels in
every slice while the photoemission moves, so it survives an average over
the slices that the photoemission does not. From that average:

1. what depends on energy only is divided out (the Fermi edge and flat
   bands sit at the same kinetic energy in every slice and are not grid);
2. only the grid's own peaks in k-space are kept -- found automatically,
   shown in red, typically 1.5 % of k-space -- so nothing else in the image
   can be touched;
3. every slice gets its own grid contrast and sub-pixel shift, fitted to it,
   refined locally over 5 × 5 tiles; the grid is re-estimated with the
   slices shifted back into register; each slice is divided by it.

Measured on held-out slices (grid from the other half of the map), the
grid's power over that of the same k-space regions without a grid:

| map | before | after |
|---|---|---|
| ANTARES WSe2 (MBS, hexagonal) | 60 | 1.4 |
| CASSIOPEE Map80eV (Scienta θ map, 28 counts/px) | 9.7 | 1.01 |
| CASSIOPEE LHhv (Scienta hν map) | 448 | 1.5 |

1 is no grid left. The per-slice contrast matters most on a kz map, where it
falls as the count rate rises (their correlation was −0.99 on LHhv). A
113-slice ANTARES map takes about a minute, with a progress bar.

**A cut** -- **Functions → Remove detector grid...** on the cut viewer. Best is a grid from a map
taken with the same lens mode and pass energy: tick **Also list the grid
pattern** when de-gridding the map and a `[grid]` dataset goes to the list;
the cut's window finds it (same detector frame) and chooses it when the
settings match. With no such grid, the cut's own high-frequency grid peaks
are notched out of I / smooth(I): the grid goes completely, and so does the
photoemission in the same k-space regions -- the window says so when the
power left drops below 1. A grid from other settings is offered but not
chosen (pass energy 20 → 50 eV: 88 % of the grid's power removed).

The window shows the slice before and after on one colour scale (and as
−∂²I/∂E², where any grid is obvious), the grid found, its k-space regions,
the contrast and shift of every slice, and the numbers above. **Advanced**
holds the parameters; the defaults are the tested best and are what the
button uses.

## Cut arithmetic

Two cuts combined into a third: LH − LV linear dichroism, circular
dichroism, dividing by a reference spectrum (gold), and plain differences,
ratios and sums. Cuts only -- a map has to be sliced first.

Two ways in, the same window either way:

- In a cut's viewer, choose **Functions → Cut arithmetic...**. That cut is **A**. A small
  prompt asks you to click the other cut in the main list; the window opens
  on that click. Clicking a map, or A's own row, is refused in the prompt
  and it keeps waiting; closing the prompt stops the wait.
- In the main list, select two cuts, right-click, **Cut arithmetic on the
  two...**. The first selected is A.

The window shows:

- **What the two measurements recorded**, side by side: polarisation,
  photon energy, theta/tilt/phi, sample X/Y/Z, temperature, pass energy,
  lens mode. Rows in orange differ beyond a tolerance (0.01 eV, 0.05°,
  5 µm, 1 K). A dichroism map between two sample spots or two temperatures
  is a map of *that* difference, so it is worth a look before trusting the
  result. Polarisation is shown but never flagged -- it is meant to differ.
- **Warnings** from the recorded polarisations: A = LV and B = LH under
  the linear-dichroism preset (the sign is reversed -- **Swap A ↔ B**),
  both the same polarisation, or linear cuts under the circular preset.
- **Purpose** presets, each only a starting point:

  | Purpose | Result | Scale B to A |
  |---|---|---|
  | Linear dichroism (LH − LV) | A − B | same total intensity |
  | Circular dichroism | (A − B)/(A + B) | same total intensity |
  | Divide by a reference | A / B, B's angular profile, mean 1 | -- |
  | Custom | any | any |

- **Scale B to A** by total intensity, or by the intensity in a region
  (x0 y0 x1 y1, or **Take the box from A's viewer**). The two polarisations
  leave the undulator with different flux; unscaled, an LD map is mostly
  that ratio. A region above E_F, or a band known not to be dichroic, is
  the better reference when the dichroism itself changes the total.
- **Divide by** the whole reference, its angular profile (summed over
  energy: the detector's channel sensitivity, without gold's own Fermi edge
  and noise), or its energy profile. The reference is normalised to a mean
  of one, so the result keeps A's scale. **Reference floor** hides the
  ratio where the reference is below 5% of its maximum.
- **Hide where A+B below** (2% by default) blanks the background, where an
  asymmetry is a ratio of two noises.
- A and B on **one** colour scale, B as it entered the arithmetic (scaled,
  on A's grid), and the result -- signed results on a red-white-blue scale
  symmetric about zero.

B is interpolated onto A's grid when the two grids differ, and blank where
B does not reach. Axes in different units (degrees against Å⁻¹) are
refused rather than lined up.

**Uncertainty.** When both cuts are raw counts on the same grid, the
Poisson standard deviation is propagated -- `√(A + s²B)` for a difference,
`2s√(AB(A+B))/(A+sB)²` for an asymmetry -- and the report gives its median
and the fraction of pixels beyond 2σ. Checked against Monte Carlo: exact
for the difference; for the asymmetry within 5% from about 20 counts per
pixel, but 12% low at 10 counts and 32% low at 4, which the report says
when it applies (bin the cuts first). There is no σ once either cut has
been processed or interpolated, nor for a ratio.

**Result to list**, **Uncertainty to list** and **All three as a figure**
leave the window. The result is a cut named after A with `LD`, `CD`,
`norm`, `diff`, `asym`, `ratio` or `sum`, and records A, B, the operation,
the scale applied to B and the settings in its processing history.
