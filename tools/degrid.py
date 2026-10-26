"""
tools/degrid.py
===============
Removing the detector grid from ARPES images, with no Qt.

What the grid is
----------------
A periodic pattern fixed on the detector: the hexagonal channel structure
of an MCP (ANTARES, the CASSIOPEE MBS end station; period ~5 px) or a
square mesh (the CASSIOPEE Scienta analyser; period ~11 px). Measured on
real maps, it is

* **multiplicative** -- ``I = S · G`` -- with an rms of 2.7 % (ANTARES) to
  10–11 % (Scienta);
* **locked to the detector pixels**: the same pattern in every slice of a
  map (grids from the odd and even slices correlate at 0.97–1.00), while
  the photoemission moves with the deflector, the polar angle or the photon
  energy;
* **not quite rigid**: from slice to slice its contrast changes by a few
  per cent (and falls with count rate -- in a photon-energy scan the
  contrast and the counts correlate at −0.99) and it moves by tenths of a
  pixel, rigidly in the ANTARES map (one displacement fits all three
  fundamental phases to 2°), with a smaller non-rigid part on top.

The method
----------
**A map is its own reference.** Everything that moves between slices
averages out of ``Σ I / Σ smooth(I)``; the grid does not. Two things in
that average are not grid and are removed before it is used:

1. what depends on energy only -- the Fermi edge and non-dispersive states
   sit at the same kinetic energy in every slice of an angle map -- by
   dividing out the average over the slit angle;
2. everything outside the grid's own Fourier peaks, by keeping only those
   regions of k-space (found automatically, ~1–2 % of it). That also takes
   the averaging noise out: the grid estimate carries ~0.05 %.

Each slice is then divided by ``1 + β·g(r − d)``, with its own contrast ``β``
and sub-pixel displacement ``d`` fitted to it (closed form, from the phases
of the fundamental peaks), and a smooth local correction of both over a
grid of tiles (a first-order Taylor term, solved by linear least squares).
The map grid is re-estimated with every slice moved back into register.

On held-out slices (grid estimated from the other half of the map) the
grid's power over that of the same k-space regions without a grid goes

* ANTARES WSe2 map (MBS, hexagonal):   60  → 1.4
* CASSIOPEE Scienta θ map (Map80eV):   9.7 → 1.01  (the counting-noise floor)
* CASSIOPEE Scienta hν map (LHhv):     448 → 1.5

-- 0.1–0.7 % of the grid's power left. The local refinement does most of
the last factor of two; a second registration pass gained nothing and is
off. Nothing outside the grid's k-space regions is touched.

**A single cut** is divided the same way by a grid pattern from a map taken
on the same detector settings (listed as ``[grid]`` when a map is
de-gridded), fitted to the cut. Without one, the fallback is a notch filter
on ``I / smooth(I)`` at the cut's own high-frequency grid peaks: it removes
the grid there completely, but also whatever photoemission lies in the same
k-space regions, and it cannot see a grid below the noise.

When not to
-----------
Only data still on the detector's pixels: nothing resampled, shifted per
column or per slice, or processed (see :func:`not_pixel_locked`). Converting
to k, a Fermi-surface correction or the kz map alignment moves the grid off
its pixels; after them there is no pattern left to find.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import fft as sfft
from scipy.ndimage import (binary_dilation, gaussian_filter, label,
                           map_coordinates, maximum_filter, median_filter)

#: The defaults: the configuration that did best on the three real maps
#: this was developed on. A user who does not know what they mean should
#: never need to touch them.
DEFAULTS = {
    "smooth": 8.0,            # px; the "smooth(I)" the grid is measured against
    "seed_threshold": 25.0,   # a grid peak stands this far above its surroundings
    "grow_threshold": 4.0,    # ...and its region extends to where it is 4× above
    "min_frequency": 0.04,    # cycles/px (period 25 px); nothing slower is grid
    "bright_fraction": 0.2,   # fit only where the intensity is above this
    "tiles": 5,               # local refinement: tiles per side (0 = off)
    "iterations": 1,          # registration passes over the map
    "notch_min_frequency": 0.08,   # single-cut fallback: leave slower peaks
}

#: Provenance that means the pixels have moved: the grid is no longer where
#: the detector put it. Parameter prefixes, as MemoryData records them.
_MOVED = {
    "fscorr.": "a Fermi-surface correction (every angle column shifted in energy)",
    "kconv.": "a k conversion (resampled onto a momentum grid)",
    "kcut.": "a cut k conversion (resampled onto a momentum grid)",
    "arbcut.": "an arbitrary cut (interpolated along a path)",
    "kz_align.": "kz map processing (every photon energy shifted in energy)",
    "kz_to_k.": "a kz conversion (resampled onto momentum)",
    "cutops.": "cut arithmetic (no longer a single measured image)",
    "degrid.": "de-gridding already",
}
#: Processing steps after which the image is no longer I = S·G.
_PROCESSED = ("smooth", "derivative", "curvature", "gradient", "background",
              "symmetry", "cut_arithmetic", "kz_align", "kz_to_k", "degrid")


class GridNotFound(ValueError):
    """No periodic pattern stands out of the noise."""


def not_pixel_locked(info: dict):
    """Why this dataset cannot be de-gridded, or None if it can."""
    info = info or {}
    for prefix, why in _MOVED.items():
        if any(str(k).startswith(prefix) for k in info):
            return f"it has been through {why}"
    if str(info.get("cassiopee.energy_reference", "")) == "per_member":
        return ("its photon energies were stacked each on its own energy "
                "axis, which shifts every slice by a different number of "
                "pixels")
    for key, value in info.items():
        if str(key).startswith("proc.step.") and not str(key).endswith(".from"):
            name = str(value).split("(")[0].strip()
            if name in _PROCESSED:
                return (f"it has been processed ({name}); de-grid the data "
                        f"as measured, then process")
    return None


# --------------------------------------------------------------------------
# Small pieces
# --------------------------------------------------------------------------
def illuminated_box(total, margin: int = 10, angle_fraction: float = 0.2,
                    energy_fraction: float = 0.05):
    """The part of the detector frame that saw electrons, less a margin.

    The dark borders are cut away because their hard edges ring across the
    whole of k-space and would drown the grid's peaks.
    """
    total = np.asarray(total, dtype=float)
    rows, cols = total.sum(1), total.sum(0)
    good_r = np.nonzero(rows > angle_fraction * rows.max())[0]
    good_c = np.nonzero(cols > energy_fraction * cols.max())[0]
    if good_r.size < 3 * margin or good_c.size < 3 * margin:
        return slice(0, total.shape[0]), slice(0, total.shape[1])
    a0, a1 = good_r[0] + margin, good_r[-1] + 1 - margin
    e0, e1 = good_c[0] + margin, good_c[-1] + 1 - margin
    if a1 - a0 < 32 or e1 - e0 < 32:
        return slice(0, total.shape[0]), slice(0, total.shape[1])
    return slice(int(a0), int(a1)), slice(int(e0), int(e1))


def _freqs(shape):
    f0 = sfft.fftfreq(shape[0])
    f1 = sfft.fftfreq(shape[1])
    return np.meshgrid(f0, f1, indexing="ij")


def _smooth(image, sigma):
    return gaussian_filter(np.asarray(image, dtype=np.float32), sigma,
                           truncate=3.0)


def find_regions(G, *, seed_threshold=DEFAULTS["seed_threshold"],
                 grow_threshold=DEFAULTS["grow_threshold"],
                 min_frequency=DEFAULTS["min_frequency"]):
    """The grid's regions of k-space in a ratio image ``G`` (≈1 on average).

    A peak is a local maximum of the (Hann-windowed, lightly smoothed) power
    that stands ``seed_threshold`` times above the median of its 41×41
    neighbourhood -- a local background, because the photoemission's own
    spectrum falls by orders of magnitude from the centre outwards and a
    global threshold would either miss the outer peaks or take the centre.
    Each region is the connected area around a peak still ``grow_threshold``
    above that background: a grid whose period varies across the detector
    gives smeared peaks, and a fixed-size notch would miss their tails.

    Returns ``(region, peaks)``; ``peaks`` is ``[(strength, f_angle, f_E)]``
    in the half plane, strongest first. Raises :class:`GridNotFound`.
    """
    G = np.asarray(G, dtype=float)
    n0, n1 = G.shape
    window = np.outer(np.hanning(n0), np.hanning(n1))
    power = np.abs(sfft.fft2((G - np.nanmean(G)) * window)) ** 2
    power = gaussian_filter(power, 1.0, mode="wrap")
    background = median_filter(power, size=41, mode="wrap")
    relative = power / np.maximum(background, 1e-30)
    F0, F1 = _freqs(G.shape)
    radius = np.hypot(F0, F1)
    seeds = ((relative == maximum_filter(relative, 15, mode="wrap"))
             & (relative > seed_threshold) & (radius > min_frequency))
    if not seeds.any():
        raise GridNotFound("no periodic pattern stands out of the noise")
    labels, _ = label(relative > grow_threshold)
    region = np.zeros(G.shape, dtype=bool)
    for i, j in np.argwhere(seeds):
        region |= labels == labels[i, j]
    region = binary_dilation(region, iterations=1) & (radius > min_frequency / 2)
    half = seeds & ((F1 > 0) | ((F1 == 0) & (F0 > 0)))
    peaks = sorted(((float(relative[i, j]), float(F0[i, j]), float(F1[i, j]))
                    for i, j in np.argwhere(half)), reverse=True)
    return region, peaks


def _fundamentals(peaks, n: int = 3, min_frequency: float = 0.05):
    """The strongest peaks, used for the phase (shift) fit. At least two
    non-parallel ones are needed for a 2-D shift; with one, the shift along
    it is still found. Only peaks of period 20 px or less: the grids seen
    are 5 and 11 px, and a slow "peak" is the likelier to be photoemission
    that did not average out -- whose phase would then drag the shift."""
    chosen = []
    for _s, a, b in peaks:
        q = np.array([a, b])
        if np.hypot(a, b) < min_frequency:
            continue
        if any(abs(q[0] * c[1] - q[1] * c[0]) < 1e-9 for c in chosen):
            continue
        chosen.append(q)
        if len(chosen) == n:
            break
    return np.array(chosen)


def _rotated(mask, degrees):
    F0, F1 = _freqs(mask.shape)
    t = np.radians(degrees)
    a = np.cos(t) * F0 - np.sin(t) * F1
    b = np.sin(t) * F0 + np.cos(t) * F1
    n0, n1 = mask.shape
    return map_coordinates(mask.astype(float), [(a * n0) % n0, (b * n1) % n1],
                           order=0, mode="wrap") > 0.5


def grid_contrast(image, region, bright=None, smooth=DEFAULTS["smooth"]):
    """How far the grid stands above the photoemission and the noise.

    Mean power of ``I/smooth(I)`` in the grid's k-space regions, over the
    mean in the same regions rotated by ±20° and 45° -- the same radii, so
    the same share of photoemission and noise, but no grid. 1 means no grid
    is left; a raw Scienta map is ~500, an ANTARES one ~50.
    """
    image = np.asarray(image, dtype=float)
    ratio = image / np.maximum(_smooth(image, smooth), 1e-9) - 1.0
    if bright is not None:
        ratio = ratio * bright
    n0, n1 = image.shape
    power = np.abs(sfft.fft2(ratio * np.outer(np.hanning(n0), np.hanning(n1)))) ** 2
    shams = [_rotated(region, d) & ~region for d in (20, -20, 45)]
    sham = np.mean([power[m].mean() for m in shams if m.any()])
    return float(power[region].mean() / max(sham, 1e-30))


# --------------------------------------------------------------------------
# The grid model and its fit to one image
# --------------------------------------------------------------------------
@dataclass
class GridModel:
    """A band-limited grid pattern on a box of the detector frame."""
    frame_shape: tuple
    box: tuple                       # (slice, slice) into the frame
    spectrum: np.ndarray             # FFT of g = G − 1 on the box, band-limited
    region: np.ndarray               # bool, k-space regions of the grid
    fundamentals: np.ndarray         # (n, 2) cycles/px, for the shift fit
    bright: np.ndarray               # bool, where the fit is weighted
    peaks: list = field(default_factory=list)

    @property
    def box_shape(self):
        return self.spectrum.shape

    def pattern(self, d=(0.0, 0.0)) -> np.ndarray:
        """g(r − d) on the box."""
        F0, F1 = _freqs(self.box_shape)
        phase = np.exp(-2j * np.pi * (F0 * d[0] + F1 * d[1]))
        return np.real(sfft.ifft2(self.spectrum * phase))

    def full(self, d=(0.0, 0.0)) -> np.ndarray:
        """``1 + g(r − d)`` over the whole frame, 1 outside the box."""
        out = np.ones(self.frame_shape)
        out[self.box] += self.pattern(d)
        return out

    @property
    def rms(self) -> float:
        g = self.pattern()
        return float(np.std(g[self.bright])) if self.bright.any() else float(np.std(g))


@dataclass
class SliceFit:
    beta: float
    shift: tuple
    local: bool = False
    empty: bool = False


def _phase_sums(r, w, qs):
    n0, n1 = r.shape
    i0 = np.arange(n0)[:, None]
    i1 = np.arange(n1)[None, :]
    out = []
    for a, b in qs:
        # separable phase factor: cheaper than a full complex image each time
        out.append(np.sum((w * r) * np.exp(-2j * np.pi * a * i0)
                          * np.exp(-2j * np.pi * b * i1)))
    return np.array(out)


def _tents(n: int, tiles: int):
    """Partition-of-unity tent functions: ``(n, tiles)``."""
    centres = (np.arange(tiles) + 0.5) * n / tiles
    width = n / tiles
    t = np.clip(1 - np.abs(np.arange(n)[:, None] - centres[None]) / width, 0, 1)
    t[: int(centres[0])] = np.eye(tiles)[0]
    t[int(centres[-1]):] = np.eye(tiles)[-1]
    return t / t.sum(1, keepdims=True)


def fit_to(model: GridModel, image_box, smooth_box, *, tiles=DEFAULTS["tiles"],
           g0=None):
    """Fit the grid to one image (the box part) and return
    ``(1 + fitted g, SliceFit)``. ``g0`` is ``model.pattern()``, passed in
    to save recomputing it for every slice of a map."""
    image_box = np.asarray(image_box, dtype=float)
    r = image_box / np.maximum(smooth_box, 1e-9) - 1.0
    w = smooth_box * model.bright
    if w.sum() <= 0 or not np.isfinite(r[model.bright]).all():
        return np.ones_like(image_box), SliceFit(0.0, (0.0, 0.0), empty=True)
    r = np.where(np.isfinite(r), r, 0.0)
    qs = model.fundamentals
    if g0 is None:
        g0 = model.pattern()
    # -- global: displacement from the phases, contrast by regression -------
    mine = _phase_sums(r, w, qs)
    theirs = _phase_sums(g0, w, qs)
    dphi = np.angle(mine * np.conj(theirs))
    shift, *_ = np.linalg.lstsq(-2 * np.pi * qs, dphi, rcond=None)
    gs = model.pattern(shift)
    denom = float(np.sum(w * gs * gs))
    beta = float(np.sum(w * r * gs) / denom) if denom > 0 else 0.0
    beta = float(np.clip(beta, 0.0, 3.0))
    fitted = beta * gs
    local = False
    # -- local: β and a small extra shift per tile, first order in the shift -
    # A tile must hold many periods for its fit to measure the grid rather
    # than the noise: at least ~100 px a side, whatever was asked for.
    tiles = min(int(tiles or 0), max(1, min(gs.shape) // 100))
    if tiles > 1 and beta > 0:
        F0, F1 = _freqs(gs.shape)
        Fgs = sfft.fft2(gs)
        gx = np.real(sfft.ifft2(Fgs * (2j * np.pi * F0)))
        gy = np.real(sfft.ifft2(Fgs * (2j * np.pi * F1)))
        t0 = _tents(gs.shape[0], tiles)
        t1 = _tents(gs.shape[1], tiles)
        basis = (gs, gx, gy)

        def tile_sums(a):
            return t0.T @ (w * a) @ t1              # (tiles, tiles)
        A = np.empty((tiles, tiles, 3, 3))
        b = np.empty((tiles, tiles, 3))
        for i in range(3):
            b[:, :, i] = tile_sums(r * basis[i])
            for j in range(i, 3):
                A[:, :, i, j] = A[:, :, j, i] = tile_sums(basis[i] * basis[j])
        total = tile_sums(np.ones_like(r))
        coef = np.zeros((tiles, tiles, 3))
        coef[:, :, 0] = beta
        ok = total > 0.01 * total.sum()
        for i, j in np.argwhere(ok):
            try:
                c = np.linalg.solve(A[i, j], b[i, j])
            except np.linalg.LinAlgError:
                continue
            # the Taylor term is only valid for a small extra shift
            if 0 < c[0] < 3 and abs(c[1] / c[0]) < 1.0 and abs(c[2] / c[0]) < 1.0:
                coef[i, j] = c
        fields = [t0 @ coef[:, :, k] @ t1.T for k in range(3)]
        fitted = fields[0] * gs + fields[1] * gx + fields[2] * gy
        local = True
    return 1.0 + fitted, SliceFit(beta, (float(shift[0]), float(shift[1])), local)


# --------------------------------------------------------------------------
# A whole map
# --------------------------------------------------------------------------
@dataclass
class DegridResult:
    values: np.ndarray               # corrected, same shape as the input
    model: GridModel
    fits: list                       # SliceFit per slice (one for a cut)
    method: str
    contrast_before: float = float("nan")    # on the preview slice / the cut
    contrast_after: float = float("nan")
    preview_index: int = 0
    notes: list = field(default_factory=list)

    def summary(self) -> str:
        m = self.model
        lines = [f"method: {self.method}",
                 f"grid: {len(m.peaks)} peaks, {100 * m.region.mean():.2f}% of "
                 f"k-space; rms {100 * m.rms:.2f}% of the intensity"]
        fits = [f for f in self.fits if not f.empty]
        if len(fits) > 1:
            betas = np.array([f.beta for f in fits])
            shifts = np.array([f.shift for f in fits])
            lines.append(f"per slice: contrast β {np.median(betas):.3f} "
                         f"(range {betas.min():.2f}–{betas.max():.2f}), shift "
                         f"std {shifts[:, 0].std():.2f} / {shifts[:, 1].std():.2f} px")
        if np.isfinite(self.contrast_before):
            lines.append(f"grid-to-background power in the grid regions: "
                         f"{self.contrast_before:.1f} → {self.contrast_after:.2f} "
                         f"(1 = no grid left)")
        return "\n".join(lines + list(self.notes))


def _slice(cube, index):
    return np.asarray(cube[index], dtype=np.float32)


def estimate_map_grid(cube, *, settings=None, progress=None):
    """The map's own grid (see the module docstring). ``cube`` is
    ``(slice, angle, energy)``, possibly lazy. Returns ``(model, fits)``."""
    s = dict(DEFAULTS, **(settings or {}))
    n = cube.shape[0]
    frame = tuple(cube.shape[1:])
    steps = 1 + int(s["iterations"])
    done = [0]

    def tick(label):
        done[0] += 1
        if progress is not None:
            progress(done[0], steps * n + n, label)

    # -- pass 1: where the detector is lit, and the plain ratio of sums -----
    total = np.zeros(frame)
    for i in range(n):
        total += _slice(cube, i)
    box = illuminated_box(total)
    shape = (box[0].stop - box[0].start, box[1].stop - box[1].start)
    numerator = np.zeros(shape, dtype=complex)
    denominator = np.zeros(shape)
    counts = np.zeros(n)
    for i in range(n):
        image = _slice(cube, i)[box]
        sm = _smooth(image, s["smooth"])
        numerator += sfft.fft2(image - sm)
        denominator += sm
        counts[i] = float(image.mean())
        tick("measuring the grid")
    if denominator.max() <= 0:
        raise GridNotFound("the map is empty")
    bright = denominator > s["bright_fraction"] * np.percentile(denominator, 99)

    def ratio_of_sums(num):
        G = 1.0 + np.real(sfft.ifft2(num)) / np.maximum(denominator, 1e-9)
        # what depends on energy only is the Fermi edge and friends, not grid
        profile = np.sum(G * bright, 0) / np.maximum(bright.sum(0), 1)
        profile = np.where(bright.sum(0) > 0, profile, 1.0)
        return G / np.where(profile > 0, profile, 1.0)[None, :]

    G = ratio_of_sums(numerator)
    region, peaks = find_regions(G, seed_threshold=s["seed_threshold"],
                                 grow_threshold=s["grow_threshold"],
                                 min_frequency=s["min_frequency"])
    taper = gaussian_filter(region.astype(float), 1.0, mode="wrap")
    qs = _fundamentals(peaks)

    def build(Gx):
        return GridModel(frame, box, sfft.fft2(Gx - 1.0) * taper, region, qs,
                         bright, peaks)

    model = build(G)
    fits = [SliceFit(1.0, (0.0, 0.0)) for _ in range(n)]
    # -- registration: fit each slice's shift, average them back into line ---
    F0, F1 = _freqs(shape)
    for _it in range(int(s["iterations"])):
        g0 = model.pattern()
        numerator = np.zeros(shape, dtype=complex)
        for i in range(n):
            image = _slice(cube, i)[box]
            sm = _smooth(image, s["smooth"])
            _, fit = fit_to(model, image, sm, tiles=0, g0=g0)
            fits[i] = fit
            d = fit.shift if not fit.empty else (0.0, 0.0)
            numerator += sfft.fft2(image - sm) * np.exp(
                2j * np.pi * (F0 * d[0] + F1 * d[1]))
            tick("registering the slices")
        model = build(ratio_of_sums(numerator))
    return model, fits


def degrid_map(cube, *, settings=None, progress=None, preview_index=None):
    """De-grid every slice of a map with the map's own grid."""
    s = dict(DEFAULTS, **(settings or {}))
    # Every slice is read three times. From a compressed file that is three
    # decompressions of the whole map; once into memory, in the file's own
    # dtype, is cheaper than that by far.
    cube = np.asarray(cube)
    model, _ = estimate_map_grid(cube, settings=s, progress=progress)
    n = cube.shape[0]
    out = np.empty(cube.shape, dtype=np.float32)
    g0 = model.pattern()
    fits = []
    counts = []
    for i in range(n):
        frame = _slice(cube, i)
        image = frame[model.box]
        sm = _smooth(image, s["smooth"])
        G, fit = fit_to(model, image, sm, tiles=s["tiles"], g0=g0)
        corrected = frame.copy()
        corrected[model.box] = image / np.maximum(G, 0.2)
        out[i] = corrected
        fits.append(fit)
        counts.append(float(image[model.bright].mean()) if model.bright.any() else 0.0)
        if progress is not None:
            progress(n * (1 + int(s["iterations"])) + i + 1,
                     n * (2 + int(s["iterations"])), "dividing out the grid")
    if preview_index is None:
        preview_index = int(np.argmax(counts))
    before = np.asarray(cube[preview_index], dtype=float)[model.box]
    after = out[preview_index][model.box]
    result = DegridResult(out, model, fits, "map grid, fitted per slice",
                          grid_contrast(before, model.region, model.bright),
                          grid_contrast(after, model.region, model.bright),
                          preview_index)
    empty = sum(f.empty for f in fits)
    if empty:
        result.notes.append(f"{empty} slice(s) with no counts in the lit "
                            f"area were left as they were")
    return result


# --------------------------------------------------------------------------
# A single cut
# --------------------------------------------------------------------------
def model_from_pattern(pattern, box=None, *, settings=None):
    """A :class:`GridModel` from a stored ``1 + g`` frame (a ``[grid]``
    dataset). ``box`` is where it is not 1; found from the pattern if not
    given. ``settings`` is accepted for symmetry with the other entry
    points; a stored pattern needs none."""
    pattern = np.asarray(pattern, dtype=float)
    if box is None:
        rows = np.nonzero(np.any(pattern != 1.0, axis=1))[0]
        cols = np.nonzero(np.any(pattern != 1.0, axis=0))[0]
        if rows.size == 0:
            raise GridNotFound("the stored grid is flat")
        box = (slice(int(rows[0]), int(rows[-1]) + 1),
               slice(int(cols[0]), int(cols[-1]) + 1))
    G = pattern[box]
    # A stored pattern is already band-limited: outside its regions the
    # spectrum is zero to rounding, so the regions are simply where it is
    # not -- the local-background test of find_regions would see nothing
    # but "peaks" against a background of zero.
    spectrum = sfft.fft2(G - 1.0)
    power = np.abs(spectrum) ** 2
    if power.max() <= 0:
        raise GridNotFound("the stored grid is flat")
    region = power > 1e-8 * power.max()
    F0, F1 = _freqs(G.shape)
    half = ((power == maximum_filter(power, 15, mode="wrap")) & region
            & ((F1 > 0) | ((F1 == 0) & (F0 > 0))))
    peaks = sorted(((float(power[i, j] / power.max()), float(F0[i, j]),
                     float(F1[i, j])) for i, j in np.argwhere(half)), reverse=True)
    peaks = [p for p in peaks if p[0] > 1e-3]
    return GridModel(pattern.shape, box, spectrum, region,
                     _fundamentals(peaks), np.ones(G.shape, bool), peaks)


def degrid_cut_with_grid(frame, pattern, *, box=None, settings=None):
    """De-grid one cut with a grid pattern measured on a map."""
    s = dict(DEFAULTS, **(settings or {}))
    frame = np.asarray(frame, dtype=float)
    pattern = np.asarray(pattern, dtype=float)
    if frame.shape != pattern.shape:
        raise ValueError(
            f"the grid is {pattern.shape[0]} × {pattern.shape[1]} detector "
            f"pixels and this cut is {frame.shape[0]} × {frame.shape[1]}: "
            f"they were not taken with the same detector window")
    model = model_from_pattern(pattern, box, settings=s)
    image = frame[model.box]
    sm = _smooth(image, s["smooth"])
    model.bright = sm > s["bright_fraction"] * np.percentile(sm, 99)
    G, fit = fit_to(model, image, sm, tiles=s["tiles"])
    out = frame.copy()
    out[model.box] = image / np.maximum(G, 0.2)
    return DegridResult(out, model, [fit], "grid from a map, fitted to this cut",
                        grid_contrast(image, model.region, model.bright),
                        grid_contrast(out[model.box], model.region, model.bright))


def degrid_cut_notch(frame, *, settings=None):
    """The fallback for a cut with no map grid: notch its own grid peaks
    out of ``I / smooth(I)``, above ``notch_min_frequency`` only."""
    s = dict(DEFAULTS, **(settings or {}))
    frame = np.asarray(frame, dtype=float)
    box = illuminated_box(frame)
    image = frame[box]
    sm = _smooth(image, s["smooth"])
    ratio = image / np.maximum(sm, 1e-9)
    region, peaks = find_regions(ratio, seed_threshold=s["seed_threshold"],
                                 grow_threshold=s["grow_threshold"],
                                 min_frequency=s["notch_min_frequency"])
    taper = gaussian_filter(region.astype(float), 1.0, mode="wrap")
    notched = np.real(sfft.ifft2(sfft.fft2(ratio) * (1.0 - taper)))
    out = frame.copy()
    out[box] = notched * sm
    bright = sm > s["bright_fraction"] * np.percentile(sm, 99)
    # what the notch took out of I/smooth(I): the grid, as far as it can tell
    model = GridModel(frame.shape, box, sfft.fft2(ratio) * taper, region,
                      _fundamentals(peaks), bright, peaks)
    result = DegridResult(out, model, [SliceFit(1.0, (0.0, 0.0))],
                          "FFT notch on this cut alone",
                          grid_contrast(image, region, bright),
                          grid_contrast(out[box], region, bright))
    if result.contrast_after < 0.5:
        result.notes.append(
            f"The power left in the grid regions ({result.contrast_after:.2f}) "
            f"is below what photoemission and noise put there elsewhere (1): "
            f"the notch has removed some of them along with the grid.")
    result.notes.append(
        "The notch also removes whatever photoemission lies in the grid's "
        "k-space regions. A grid from a map with the same detector settings "
        "does not: de-grid such a map first, then use its [grid] here.")
    return result
