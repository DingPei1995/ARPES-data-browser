"""
tools/crystal/bands.py
======================
Calculated band structures, and the arithmetic of comparing them with ARPES.

Two kinds of calculation, because they can do different things:

* **Bands everywhere** -- a band structure that can be evaluated at any
  momentum. Either eigenvalues on a 3-D grid over the zone
  (:class:`BandGrid`; an XCrySDen ``.bxsf`` from Wien2k, Quantum
  ESPRESSO, VASP or Wannier90, or a full-zone VASP run), interpolated
  periodically, or a Wannier / tight-binding Hamiltonian
  (:class:`WannierModel`; Wannier90's ``_hr.dat``), diagonalised exactly
  at each k. Either can be cut along any path, on any k_z, on the curved
  k_z surface a photon energy really measures, and contoured at any
  energy.
* **Bands along a path** (:class:`BandPath`) -- what most band-structure
  runs write: energies along a list of segments between labelled points.
  These can be drawn and compared along that path only.

Both carry an optional set of **weights** -- orbital or atomic character
(VASP's PROCAR, a Wannier function's share of each eigenvector) -- for
colouring the bands.

The comparison with a measurement (:func:`kz_on_sphere`,
:func:`bands_along_cut`, :func:`contours`, :func:`fit_alignment`,
:func:`spectral_function`) is here too, Qt-free; the files are read and
written by :mod:`tools.crystal.bandfiles`, the window is
:mod:`ui.panels.bands`.

Units: momenta in Å⁻¹ with 2π included (the program's convention), real
space in Å, energies in eV, and every energy relative to the calculation's
Fermi level once :attr:`fermi` is subtracted (the classes do that on
evaluation, so "E = 0" is E_F throughout).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import re

import numpy as np
from scipy import ndimage

__all__ = ["BandGrid", "WannierModel", "BandPath", "sample_path", "kz_on_sphere",
           "kz_of", "bands_along_cut", "contours", "contour_lines",
           "fit_alignment", "reciprocal_of", "parse_path", "named_path",
           "spectral_function", "Alignment", "K_PER_SQRT_EV", "KZ_FACTOR",
           "node_line", "follow_band"]

#: k (Å⁻¹) = K_PER_SQRT_EV * sqrt(E_kin (eV)) for a free electron
K_PER_SQRT_EV = 0.5123167
#: k² (Å⁻²) = KZ_FACTOR * E (eV)
KZ_FACTOR = K_PER_SQRT_EV ** 2


def reciprocal_of(lattice) -> np.ndarray:
    """Reciprocal vectors (rows, 2π included) of real-space ``lattice``
    (rows, Å)."""
    lattice = np.asarray(lattice, dtype=float)
    return 2.0 * np.pi * np.linalg.inv(lattice).T


# --------------------------------------------------------------------------
# The calculations
# --------------------------------------------------------------------------
class _Everywhere:
    """What :class:`BandGrid` and :class:`WannierModel` share."""

    reciprocal: np.ndarray
    fermi: float = 0.0
    weight_names: list = []
    name: str = ""

    @property
    def lattice(self) -> np.ndarray:
        return 2.0 * np.pi * np.linalg.inv(self.reciprocal).T

    def fractional(self, k_points) -> np.ndarray:
        """Cartesian k (rows, Å⁻¹) in units of the reciprocal vectors."""
        k_points = np.atleast_2d(np.asarray(k_points, dtype=float))
        if k_points.shape[1] == 2:
            k_points = np.column_stack([k_points, np.zeros(len(k_points))])
        return k_points @ np.linalg.inv(self.reciprocal)

    def energies(self, k_points) -> np.ndarray:            # pragma: no cover
        raise NotImplementedError

    def band_count(self) -> int:                           # pragma: no cover
        raise NotImplementedError


class BandGrid(_Everywhere):
    """Eigenvalues on a regular grid over the reciprocal cell.

    ``values`` is ``(bands, n1, n2, n3)``, point ``(i, j, l)`` at fractional
    ``origin + (i/n1, j/n2, l/n3)`` -- a *periodic* grid; a "general grid"
    that repeats its first plane at the end (as ``.bxsf`` files do) is
    trimmed by :meth:`from_general_grid`. Interpolated with periodic cubic
    splines, so a band is smooth between grid points and continuous across
    the zone boundary.
    """

    def __init__(self, reciprocal, values, *, origin=(0.0, 0.0, 0.0),
                 fermi: float = 0.0, weights=None, name: str = ""):
        self.reciprocal = np.asarray(reciprocal, dtype=float)
        self.values = np.asarray(values, dtype=float)
        if self.values.ndim != 4:
            raise ValueError("band values need the shape (bands, n1, n2, n3)")
        self.origin = np.asarray(origin, dtype=float)
        self.fermi = float(fermi)
        self.weights = dict(weights or {})
        self.weight_names = list(self.weights)
        self.name = name
        self._coeffs = None

    @classmethod
    def from_general_grid(cls, reciprocal, values, **options):
        """From a grid whose last point along each axis repeats its first."""
        values = np.asarray(values, dtype=float)[:, :-1, :-1, :-1]
        weights = {k: np.asarray(v)[:, :-1, :-1, :-1]
                   for k, v in (options.pop("weights", None) or {}).items()}
        return cls(reciprocal, values, weights=weights, **options)

    def band_count(self) -> int:
        return self.values.shape[0]

    def _grid_coords(self, k_points):
        frac = self.fractional(k_points) - self.origin
        n = np.array(self.values.shape[1:], dtype=float)
        return (frac * n).T                   # (3, N), periodic: wrap

    def _spline(self, array):
        return ndimage.spline_filter(array, order=3, mode="grid-wrap")

    def energies(self, k_points) -> np.ndarray:
        """``(bands, N)`` energies relative to E_F at the k points."""
        coords = self._grid_coords(k_points)
        if self._coeffs is None:
            self._coeffs = [self._spline(band) for band in self.values]
        return np.array([ndimage.map_coordinates(c, coords, order=3,
                                                 mode="grid-wrap", prefilter=False)
                         for c in self._coeffs]) - self.fermi

    def weights_at(self, k_points, name: str) -> np.ndarray:
        """``(bands, N)`` weight ``name`` at the k points (linear)."""
        coords = self._grid_coords(k_points)
        return np.array([ndimage.map_coordinates(w, coords, order=1,
                                                 mode="grid-wrap")
                         for w in self.weights[name]])


class WannierModel(_Everywhere):
    """A Wannier (or tight-binding) Hamiltonian: ``H(k) = sum_R H(R)
    exp(2 pi i f.R) / deg(R)``, with ``f`` the fractional k and ``R`` in
    lattice vectors -- Wannier90's ``_hr.dat`` convention. Diagonalised at
    every k, so there is no interpolation error at all; the weights are
    each Wannier function's share of each eigenvector (``wf1``, ``wf2``...,
    or the names given), and of each group of them in ``groups`` (name ->
    indices, e.g. every Cr 3d function), listed first."""

    def __init__(self, lattice, R, H, *, degeneracy=None, fermi: float = 0.0,
                 orbital_names=None, groups=None, name: str = ""):
        self.reciprocal = reciprocal_of(lattice)
        self.R = np.asarray(R, dtype=float).reshape(-1, 3)
        self.H = np.asarray(H, dtype=complex)
        self.degeneracy = (np.ones(len(self.R)) if degeneracy is None
                           else np.asarray(degeneracy, dtype=float))
        self.fermi = float(fermi)
        size = self.H.shape[1]
        self.orbital_names = list(orbital_names or
                                  [f"wf{i + 1}" for i in range(size)])
        self.groups = {k: list(v) for k, v in (groups or {}).items()}
        self.weight_names = list(self.groups) + list(self.orbital_names)
        self.name = name

    def band_count(self) -> int:
        return self.H.shape[1]

    def _hamiltonians(self, k_points):
        frac = self.fractional(k_points)
        phase = np.exp(2j * np.pi * frac @ self.R.T) / self.degeneracy   # (N, M)
        return np.einsum("nm,mab->nab", phase, self.H)

    def diagonalise(self, k_points):
        h = self._hamiltonians(k_points)
        h = 0.5 * (h + np.conj(np.transpose(h, (0, 2, 1))))
        return np.linalg.eigh(h)

    def energies(self, k_points) -> np.ndarray:
        values, _vectors = self.diagonalise(k_points)
        return values.T - self.fermi

    def weights_at(self, k_points, name: str) -> np.ndarray:
        index = self.groups.get(name) or [self.orbital_names.index(name)]
        _values, vectors = self.diagonalise(k_points)
        return (np.abs(vectors[:, index, :]) ** 2).sum(axis=1).T


@dataclass
class BandPath:
    """Energies along a path, as a band-structure run writes them:
    ``distance`` (N,) in Å⁻¹, ``energies`` (bands, N), the labelled points
    as ``[(distance, label), ...]``, optional ``k`` (N, 3) Cartesian, and
    ``weights`` name -> (bands, N)."""
    distance: np.ndarray
    energies: np.ndarray
    labels: list = field(default_factory=list)
    k: np.ndarray = None
    fermi: float = 0.0
    weights: dict = field(default_factory=dict)
    name: str = ""

    @property
    def weight_names(self):
        return list(self.weights)

    def relative(self) -> np.ndarray:
        """The energies relative to E_F."""
        return np.asarray(self.energies, dtype=float) - self.fermi


# --------------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------------
def sample_path(points, labels, step: float = 0.01):
    """Points along the segments between ``points`` (Cartesian, Å⁻¹), about
    ``step`` Å⁻¹ apart. ``labels`` names each point; a label "|" joins two
    segments without a line between them (X|U: a jump). Returns ``(k (N,
    3), distance (N,), ticks [(distance, label)])``."""
    points = [np.asarray(p, dtype=float) for p in points]
    points = [np.concatenate([p, [0.0]]) if p.size == 2 else p for p in points]
    ks, distance, ticks = [], [], []
    position = 0.0
    for i in range(len(points) - 1):
        a, b = points[i], points[i + 1]
        length = float(np.linalg.norm(b - a))
        n = max(2, int(np.ceil(length / step)) + 1)
        t = np.linspace(0.0, 1.0, n)
        segment = a[None, :] + t[:, None] * (b - a)[None, :]
        if i:
            segment, t = segment[1:], t[1:]
        ks.append(segment)
        distance.append(position + t * length)
        if i == 0:
            ticks.append((position, labels[0]))
        position += length
        ticks.append((position, labels[i + 1]))
    return np.vstack(ks), np.concatenate(distance), ticks


def parse_path(text: str):
    """Labels of a path written as "Γ-M-K-Γ-A" (or "G M K G A"; "G" is Γ;
    "|" marks a jump, as in "X|U")."""
    tokens = re.findall(r"\||[A-Za-zΓΣΔΛ][A-Za-z0-9_']*", str(text or ""))
    out = []
    for token in tokens:
        out.append("Γ" if token.upper() in ("G", "GAMMA") else token)
    return out


def named_path(params, text: str, step: float = 0.01):
    """``sample_path`` through the named high-symmetry points of the
    crystal ``params`` (:class:`tools.crystal.lattice.LatticeParams`), e.g.
    "Γ-M-K-Γ-A-L-H-A". Raises ``ValueError`` for an unknown label."""
    from tools.crystal.highsym import high_symmetry_points

    table = {p.label: np.asarray(p.k, dtype=float)
             for p in high_symmetry_points(params)}
    labels = parse_path(text)
    if len(labels) < 2:
        raise ValueError("a path needs at least two points")
    unknown = [l for l in labels if l != "|" and l not in table]
    if unknown:
        raise ValueError(f"not a high-symmetry point of this zone: "
                         f"{', '.join(unknown)} (it has {', '.join(table)})")
    # a jump "X|U": both points kept, no segment between them
    points, names, jumps = [], [], []
    for label in labels:
        if label == "|":
            jumps.append(len(points) - 1)
            continue
        points.append(table[label])
        names.append(label)
    k, distance, ticks = sample_path(points, names, step)
    if jumps:
        # remove the segments that were jumps: they have no length
        keep = np.ones(len(k), dtype=bool)
        edges = np.searchsorted(distance, [t[0] for t in ticks])
        for j in jumps:
            lo, hi = edges[j], edges[j + 1]
            keep[lo + 1:hi] = False
            shift = ticks[j + 1][0] - ticks[j][0]
            distance[hi:] -= shift
            ticks = (ticks[:j + 1] + [(ticks[j + 1][0] - shift,
                                       f"{ticks[j][1]}|{ticks[j + 1][1]}")]
                     + [(d - shift, t) for d, t in ticks[j + 2:]])
            ticks.pop(j)
        k, distance = k[keep], distance[keep]
    return k, distance, ticks


# --------------------------------------------------------------------------
# The k_z a photon energy measures
# --------------------------------------------------------------------------
def kz_of(k_parallel, energy, hv: float, inner_potential: float,
          work_function: float = 4.4):
    """k_z (Å⁻¹) of a photoelectron with in-plane momentum ``k_parallel``
    (Å⁻¹) and ``energy`` (eV, relative to E_F), free-electron final state:
    ``k_z = sqrt(0.2625 (E_kin + V0) - k_par²)``."""
    e_kin = hv - work_function + np.asarray(energy, dtype=float)
    return np.sqrt(np.clip(KZ_FACTOR * (e_kin + inner_potential)
                           - np.asarray(k_parallel, dtype=float) ** 2, 0.0, None))


def kz_on_sphere(k_inplane, energy, hv, inner_potential, work_function=4.4):
    """Full k vectors (N, 3) for in-plane points (N, 2) measured at ``hv``:
    each point's k_z from the free-electron final state."""
    k_inplane = np.atleast_2d(np.asarray(k_inplane, dtype=float))
    kpar = np.hypot(k_inplane[:, 0], k_inplane[:, 1])
    kz = kz_of(kpar, energy, hv, inner_potential, work_function)
    return np.column_stack([k_inplane[:, :2], kz])


def bands_along_cut(model, k_parallel, direction_deg: float = 0.0, *,
                    kz: float = None, hv: float = None, inner_potential=12.0,
                    work_function: float = 4.4, rounds: int = 3):
    """The calculated bands along a measured cut: ``(bands, N)`` energies
    at in-plane momenta ``k_parallel`` (Å⁻¹) along ``direction_deg`` from
    kx.

    With ``kz`` the cut lies on that plane; with ``hv`` (and the inner
    potential and work function) on the curved surface that photon energy
    measures. There k_z also depends on the energy itself, so each band is
    found self-consistently: k_z from the band's own energy, a few rounds.
    """
    k_parallel = np.asarray(k_parallel, dtype=float)
    phi = np.deg2rad(direction_deg)
    inplane = np.column_stack([k_parallel * np.cos(phi), k_parallel * np.sin(phi)])
    if hv is None:
        k = np.column_stack([inplane, np.full(len(inplane), float(kz or 0.0))])
        return model.energies(k)
    energy = np.zeros((model.band_count(), len(k_parallel)))
    for _round in range(max(1, rounds)):
        out = np.empty_like(energy)
        for band in range(model.band_count()):
            k = kz_on_sphere(inplane, energy[band], hv, inner_potential,
                             work_function)
            out[band] = model.energies(k)[band]
        energy = out
    return energy


def follow_band(model, k_points, start: int, band: int,
                diabatic: bool = True) -> np.ndarray:
    """One band's energies along ``k_points`` (N, 3), followed by its
    character rather than its place in the energy order: from point
    ``start`` (where it is band ``band`` counted from the bottom) each next
    point takes the eigenvector that overlaps most with the band's
    character at ``start`` (``diabatic``, the default: through crossings
    and small anticrossings alike, as a band is seen in a measurement) or
    with the previous point's (adiabatic: round an anticrossing). Either
    way a band goes straight through a true crossing (a Dirac point),
    where sorting the energies would bend it onto the other band. For a
    Hamiltonian (:class:`WannierModel`); a band grid has no eigenvectors
    and is followed by the energies' own continuity instead."""
    k_points = np.asarray(k_points, dtype=float)
    n = len(k_points)
    out = np.full(n, np.nan)
    if isinstance(model, WannierModel):
        values, vectors = model.diagonalise(k_points)
        values = values - model.fermi
        for direction in (1, -1):
            current = band
            out[start] = values[start, current]
            previous = vectors[start, :, current]
            reference = previous
            for i in range(start + direction, n if direction == 1 else -1, direction):
                overlap = np.abs(np.conj(reference) @ vectors[i]) ** 2
                current = int(np.argmax(overlap))
                previous = vectors[i, :, current]
                if not diabatic:
                    reference = previous
                out[i] = values[i, current]
        return out
    energies = model.energies(k_points)
    for direction in (1, -1):
        current, last, slope = band, energies[band, start], 0.0
        out[start] = last
        for i in range(start + direction, n if direction == 1 else -1, direction):
            current = int(np.argmin(np.abs(energies[:, i] - (last + slope))))
            slope, last = energies[current, i] - last, energies[current, i]
            out[i] = last
    return out


def node_line(model, k_parallel, kz_values, direction_deg: float = 0.0,
              bands=None):
    """Where two bands touch along a cut, followed through k_z: a Dirac
    line node. For each k_z the cut (in-plane momenta ``k_parallel`` along
    ``direction_deg``) is searched for the smallest gap between ``bands``
    (two band indices; default: the pair closest anywhere on the first
    cut). Returns ``(k (n,), E (n,), gap (n,))`` -- the node's momentum,
    energy (the two bands' mean there) and the gap left (0 for a true
    crossing, more where it opens)."""
    k_parallel = np.asarray(k_parallel, dtype=float)
    out_k, out_e, out_gap = [], [], []
    for kz in np.atleast_1d(np.asarray(kz_values, dtype=float)):
        energies = bands_along_cut(model, k_parallel, direction_deg, kz=float(kz))
        if bands is None:
            gaps = np.abs(np.diff(energies, axis=0))
            row = int(np.unravel_index(np.argmin(gaps), gaps.shape)[0])
            bands = (row, row + 1)
        i, j = bands
        gap = np.abs(energies[j] - energies[i])
        at = int(np.argmin(gap))
        # sub-sample: the gap is |linear| near a crossing, so interpolate the
        # signed difference to its zero where it changes sign
        signed = energies[j] - energies[i]
        k_at, e_at = k_parallel[at], 0.5 * (energies[i][at] + energies[j][at])
        for lo in (at - 1, at):
            if 0 <= lo < k_parallel.size - 1 and signed[lo] * signed[lo + 1] < 0:
                t = signed[lo] / (signed[lo] - signed[lo + 1])
                k_at = k_parallel[lo] + t * (k_parallel[lo + 1] - k_parallel[lo])
                e_at = float(np.interp(k_at, k_parallel, 0.5 * (energies[i] + energies[j])))
                break
        out_k.append(float(k_at))
        out_e.append(float(e_at))
        out_gap.append(float(gap[at]))
    return np.array(out_k), np.array(out_e), np.array(out_gap)


# --------------------------------------------------------------------------
# Contours
# --------------------------------------------------------------------------
_EDGES = {
    # marching squares: corner bits (bl=1, br=2, tr=4, tl=8) -> edge pairs;
    # edges: 0 bottom, 1 right, 2 top, 3 left
    1: [(3, 0)], 2: [(0, 1)], 3: [(3, 1)], 4: [(1, 2)], 5: [(3, 2), (0, 1)],
    6: [(0, 2)], 7: [(3, 2)], 8: [(2, 3)], 9: [(2, 0)], 10: [(2, 1), (0, 3)],
    11: [(2, 1)], 12: [(1, 3)], 13: [(1, 0)], 14: [(0, 3)],
}


def contour_lines(values, x, y, level: float):
    """Iso-lines of ``values[ix, iy]`` at ``level``: a list of ``(xs, ys)``
    polylines (marching squares, the segments joined end to end)."""
    values = np.asarray(values, dtype=float)
    x, y = np.asarray(x, float), np.asarray(y, float)
    above = values > level
    segments = []
    nx, ny = values.shape

    def point(i, j, edge):
        if edge == 0:            # bottom: (i, j) - (i+1, j)
            a, b, p, q = values[i, j], values[i + 1, j], (x[i], y[j]), (x[i + 1], y[j])
        elif edge == 1:          # right: (i+1, j) - (i+1, j+1)
            a, b, p, q = values[i + 1, j], values[i + 1, j + 1], (x[i + 1], y[j]), (x[i + 1], y[j + 1])
        elif edge == 2:          # top: (i+1, j+1) - (i, j+1)
            a, b, p, q = values[i + 1, j + 1], values[i, j + 1], (x[i + 1], y[j + 1]), (x[i], y[j + 1])
        else:                    # left: (i, j+1) - (i, j)
            a, b, p, q = values[i, j + 1], values[i, j], (x[i], y[j + 1]), (x[i], y[j])
        t = 0.5 if a == b else (level - a) / (b - a)
        return (p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1]))

    code = (above[:-1, :-1] * 1 + above[1:, :-1] * 2 + above[1:, 1:] * 4
            + above[:-1, 1:] * 8)
    for i, j in zip(*np.nonzero((code > 0) & (code < 15))):
        for e0, e1 in _EDGES[int(code[i, j])]:
            segments.append((point(i, j, e0), point(i, j, e1)))
    return _join(segments)


def _join(segments, places: int = 9):
    """Segments joined into polylines wherever they share an end."""
    def key(p):
        return (round(p[0], places), round(p[1], places))

    ends = {}
    for index, (a, b) in enumerate(segments):
        ends.setdefault(key(a), []).append((index, 0))
        ends.setdefault(key(b), []).append((index, 1))
    used = [False] * len(segments)
    lines = []
    for start in range(len(segments)):
        if used[start]:
            continue
        used[start] = True
        line = [segments[start][0], segments[start][1]]
        for forward in (True, False):
            while True:
                tip = line[-1] if forward else line[0]
                nxt = None
                for index, side in ends.get(key(tip), ()):
                    if not used[index]:
                        nxt = (index, side)
                        break
                if nxt is None:
                    break
                index, side = nxt
                used[index] = True
                other = segments[index][1 - side]
                if forward:
                    line.append(other)
                else:
                    line.insert(0, other)
        lines.append((np.array([p[0] for p in line]), np.array([p[1] for p in line])))
    return lines


def contours(model, kx, ky, energy: float, *, kz: float = None, hv=None,
             inner_potential: float = 12.0, work_function: float = 4.4,
             rotation_deg: float = 0.0):
    """Every band's constant-energy contour at ``energy`` on the grid
    ``kx`` x ``ky`` (Å⁻¹): ``[(band, xs, ys), ...]``. The plane is
    ``k_z = kz``, or the curved surface ``hv`` measures. ``rotation_deg``
    turns the calculation against the map (the sample azimuth)."""
    gx, gy = np.meshgrid(np.asarray(kx, float), np.asarray(ky, float),
                         indexing="ij")
    c, s = np.cos(np.deg2rad(-rotation_deg)), np.sin(np.deg2rad(-rotation_deg))
    inplane = np.column_stack([c * gx.ravel() - s * gy.ravel(),
                               s * gx.ravel() + c * gy.ravel()])
    if hv is None:
        k = np.column_stack([inplane, np.full(len(inplane), float(kz or 0.0))])
    else:
        k = kz_on_sphere(inplane, energy, hv, inner_potential, work_function)
    values = model.energies(k).reshape(model.band_count(), *gx.shape)
    out = []
    for band, grid in enumerate(values):
        if grid.min() <= energy <= grid.max():
            for xs, ys in contour_lines(grid, kx, ky, energy):
                out.append((band, xs, ys))
    return out


# --------------------------------------------------------------------------
# Aligning the calculation with the measurement
# --------------------------------------------------------------------------
@dataclass
class Alignment:
    """``E_measured = Z * E_calculated + shift``: ``Z`` the band
    renormalisation (1/Z the mass enhancement), ``shift`` the rigid
    offset; with their standard errors and the points' rms distance."""
    z: float = 1.0
    shift: float = 0.0
    z_err: float = float("nan")
    shift_err: float = float("nan")
    rms: float = float("nan")
    points: int = 0

    def apply(self, energies):
        return self.z * np.asarray(energies, dtype=float) + self.shift


def fit_alignment(k_points, energies, calculated, *, fit_z: bool = True,
                  window: float = 0.3) -> Alignment:
    """``Z`` and the shift that best put the calculated bands on measured
    band points.

    ``k_points`` are positions along the cut (any coordinate the
    ``calculated`` function takes), ``energies`` the measured energies (eV,
    E_F = 0), ``calculated(k_points) -> (bands, N)`` the calculation there.
    Each point is matched to the nearest band; points farther than
    ``window`` eV from every band do not count (a feature the calculation
    does not have). Robust least squares (soft-L1).
    """
    from scipy.optimize import least_squares

    energies = np.asarray(energies, dtype=float)
    calc = np.asarray(calculated(k_points), dtype=float)

    def residual(params):
        z = params[0] if fit_z else 1.0
        shift = params[-1]
        distance = z * calc + shift - energies[None, :]
        nearest = distance[np.argmin(np.abs(distance), axis=0),
                           np.arange(len(energies))]
        return np.clip(nearest, -window, window)

    start = [1.0, 0.0] if fit_z else [0.0]
    # a coarse search over the shift first: the nearest-band matching has
    # many local minima a shift apart
    best = min(np.linspace(-0.4, 0.4, 41),
               key=lambda s: np.sum(np.abs(residual(start[:-1] + [s]))))
    start[-1] = float(best)
    result = least_squares(residual, start, loss="soft_l1", f_scale=0.03)
    params = result.x
    res = result.fun
    kept = np.abs(res) < window
    dof = max(1, int(kept.sum()) - len(params))
    try:
        cov = np.linalg.inv(result.jac.T @ result.jac) * (np.sum(res[kept] ** 2) / dof)
        errors = np.sqrt(np.clip(np.diag(cov), 0, None))
    except np.linalg.LinAlgError:
        errors = np.full(len(params), np.nan)
    rms = float(np.sqrt(np.mean(res[kept] ** 2))) if kept.any() else float("nan")
    if fit_z:
        return Alignment(float(params[0]), float(params[1]), float(errors[0]),
                         float(errors[1]), rms, int(kept.sum()))
    return Alignment(1.0, float(params[0]), 0.0, float(errors[0]), rms,
                     int(kept.sum()))


def spectral_function(band_energies, energy_axis, width: float = 0.03,
                      weights=None) -> np.ndarray:
    """Lorentzian spectral weight ``(N, n_E)`` of bands ``(bands, N)`` on
    ``energy_axis``; ``weights`` (bands, N) scales each band."""
    band_energies = np.atleast_2d(np.asarray(band_energies, dtype=float))
    e = np.asarray(energy_axis, dtype=float)
    w = np.ones_like(band_energies) if weights is None else np.asarray(weights, float)
    out = np.zeros((band_energies.shape[1], e.size))
    for band, weight in zip(band_energies, w):
        out += weight[:, None] * (width / np.pi) / ((e[None, :] - band[:, None]) ** 2
                                                    + width ** 2)
    return out
