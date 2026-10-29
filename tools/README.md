# Analysis: technical reference

What each analysis does and what it assumes. The computations are in this folder, Qt-free: `momentum/`, `spectra/`, `processing/`, `crystal/`, `output/`, each with its own reference where there is one: [momentum](momentum/README.md) · [spectra](spectra/README.md) · [processing](processing/README.md) · [output](output/README.md).

## Analysis

- **k conversion** (cube and single cut) -- the free-electron formula, with
  the inputs shown and checkable. Refused, with a reason, for an axis that is
  not an emission angle. **Set rotation from contour** reads the sample
  rotation off a direction picked on the contour: two points mark a line, so
  the order you click them in and which side of the origin they sit on make
  no difference, and what comes back is the smaller of the two turns that
  stand it vertical. A converted Map is stored in float16 by default (half
  the memory and file size; its 0.05 % steps are far below any point's
  counting noise), with float32 one click away. Values outside float16's
  range are kept in float32 automatically, and every slice is widened to
  float32 before it is summed, cut or fitted.
- **Fermi level** -- fit the edge, offset the energy axis, or divide it out.
- **MDC / EDC fitting** -- fit the peaks across a cut and read the band off
  them: `v_F`, `m*`, and the self-energy. Only accepts a cut whose in-plane
  axis is already in Å⁻¹, because every quantity built on the fitted band
  would otherwise carry meaningless units.
- **Core-level check** -- which element and orbital each peak of a cut's
  angle-integrated EDC is; see below.
- **Arbitrary cuts** -- any straight line through a contour.
- **Fermi-surface correction**, **stack plots**.
- **Cut arithmetic** -- linear and circular dichroism, dividing by a
  reference, differences, ratios and sums between two cuts; see below.
- **Processing** (2-D and over a cube) -- smoothing, derivatives, curvature,
  backgrounds, symmetrisation, despiking.
- **Crop / bin / normalise** -- truncate, self-normalise, compress, on
  several datasets at once when they really are the same format.
- **Combining datasets** -- *Stack* cuts into a map along a new axis (a
  photon-energy series of angle cuts becomes a kz map), or *Join* cuts,
  maps, kz maps or k-maps of one kind end to end along one of their own
  axes (`tools/processing/combine.py`). Pieces on one grid are copied as measured;
  a different step or offset is interpolated onto the first piece's step;
  other axes that differ are resampled onto the range all pieces cover;
  overlaps are averaged (or the earlier / later piece kept); gaps are left
  at zero and reported.
- **Brillouin zone** -- overlay a conventional or irreducible zone, or a
  moiré zone, on a converted contour. Space group first, then the lattice
  parameters it constrains; a 3-D preview shows how the cut plane sits,
  with the named high-symmetry points (Γ, M, K, Z, L, ...) marked, and a
  **Through point** list offers the offsets that put the cut through one of
  them. The points the plane holds are labelled on the contour too. The
  dialog starts from the last crystal used (`~/.arpes_viewer/brillouin_zone.json`).
- **kz map processing** and **kz -> momentum** -- see below.
- **The curve viewer**, its **curve fit** and **spin analysis** -- see below.
- **De-grid** -- remove the detector's grid (MCP pattern or mesh) from a map
  or a cut; see below.
- **Detector corrections**, **sum repeated scans**, **ridge detection**,
  the **tight-binding fit** and the **crystal structure** -- ported from
  the lab's MATLAB programs, each improved where the MATLAB version could
  be; see *Tools ported from the MATLAB programs* below.

## Tools ported from the MATLAB programs

Five tools from `DataAnalysisPrograms_PD` (yulinARPES, Quick_Panel,
PD_program, TB_zph_zmd, Lattice_plot_zqs), each checked against the MATLAB
original and changed where the original could be bettered. The windows that
take several decisions lead through them as numbered steps (*Step 2 of 4:
...*, with Back / Next, and Next greyed with the reason until the step has
what the next one needs; `ui/common/stepper.py`).

**Sum repeated scans** (Datasets box; `add_two_value_arrays.m`). A fine cut
or map taken in batches, so that an ageing sample can be seen, added back
together in three steps: *are these repeats?* (grids and recorded
conditions against the first scan; a slightly different grid is
interpolated onto the first), *has the sample changed?* (each scan against
the sum of the **others**, with counting statistics: χ²/point ≈ 1 for noise
only, flagged when it stands well above the set's median; the count rate;
the energy shift of its EDC, measured by cross-correlating derivatives
against the others aligned by their own shifts), and *sum* (the changed
scans start unticked; shifts undone if asked, the energy axis then cropped
to what every scan covers; sum or mean). The MATLAB tool added two arrays
with no check at all.

**Ridge detection** (Advanced; `ridge_detection.m`). Bands as ridges of the
smoothed intensity: where the Hessian curves down across the band and the
slope across it changes sign, each point placed to a fraction of a pixel.
Changed from the .m file: smoothing per axis in that axis' unit (it used
one radius in pixels of mixed units); derivatives in units of the smoothing
width, so the result does not depend on whether energy is in eV or meV;
Gaussian-derivative filters; a sign-change test instead of an absolute
slope tolerance; vectorised. On a synthetic parabolic band the points sit
within 0.2 meV (median). Results: the ridge image to the list (a map plane
by plane), the points to a CSV.

**Detector corrections** (Basic, under Remove detector grid;
`StretchDitherImage.m`, `InterpRegionAlongX/Y.m`, `inpaint_nans.m`). On the
detector's own pixels, a cut or every slice of a map: a **keystone** (the
lit trapezoid) stretched to a rectangle, as an exact inverse map with each
line's counts kept (the .m file re-gridded a Delaunay triangulation of every
pixel and did not conserve counts); and **bad regions** patched -- linear
along x or y as in MATLAB, or a harmonic (Laplace) fill, factorised once
for a whole map. *Find bad pixels* compares the mean detector image with
its own median filter and its pixel-to-pixel scatter.

**Tight-binding fit** (Calculation; `TB_zph_zmd`). E(k) = e0 + Σ t_n Σ cos(k·R)
over shells of neighbours of any 2-D or 3-D lattice (typed, or the
Brillouin-zone tool's last crystal), fitted to band points traced by ridge
detection on constant-energy planes of a k-map, clicked, or read from a
file. A **linear** least-squares fit with standard errors (the .m file ran
a nonlinear fit from all-ones starting values); robust to stray points
(Tukey reweighting); the azimuth can be refined. One contour alone fixes the
shape and not the energy scale, so t1 is then held (the .m file's fit was
undetermined there); a condition-number check says when the points cannot
tell the shells apart. The model's contours are drawn over the data, and
its spectral function goes to the list as a map.

**Crystal structure** (Calculation; `Lattice_plot_zqs`). The crystal in 3-D
-- typed in (cell, space-group centering, basis), from the Brillouin-zone
tool's last crystal, or from a **CIF** (symmetry operations applied).
Element by element colour, size, transparency and gloss; cell edges and
bonds; a **cleavage plane** (hkl) slid along its normal or put on one of the
listed atomic layers, with the atoms above faded or hidden to show the
termination. Drawn with QPainter (no OpenGL), so it renders at any size;
*To figure composer* puts it on a figure page.
The Brillouin zone in the same view, with k axes as arrows, a coloured path
and a measured k-map painted on the plane through Γ or the zone top.

**Calculated bands** (Calculation). A band structure calculated elsewhere
(Wannier90 `_hr.dat` and the other formats in the guide), with orbital
weights: bands along a path, constant-energy contours, and a comparison
with a measured cut or map -- a renormalisation Z and an energy shift fitted
to the cut's bands. Bands can be followed through crossings by eigenvector
overlap; the calculated bands go to the figure composer as dashed guides.

**Stacked slices** (Advanced). Constant-energy or constant-k_z slices of a
volume as separated, faded layers in 3-D, with curves marked through them.
