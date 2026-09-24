"""
A band-structure calculation, read in, and a photon-energy scan matched to it.

Which photon energies reach the high-symmetry planes is the question a kz
scan is taken to answer, and the lattice alone answers it only as well as
V0 is known. A calculation answers it better: along the surface normal
(Gamma-A in a hexagonal crystal, Gamma-Z in a tetragonal one) the bands
disperse with k_z, and the same dispersion is in the scan at normal
emission, as intensity against photon energy and binding energy. Map the
calculated bands through ``k_z(hv, E) = sqrt(A (hv - W + E + V0))`` for a
trial V0 and energy shift, compare with the measurement, keep the best --
and the photon energies of the Gamma and zone-boundary planes follow.

Reading the file
----------------
Band files come in a few shapes, all read here:

* two columns, one block per band, blocks separated by blank lines
  (Wannier90 ``*_band.dat``, Quantum ESPRESSO ``bands.out.gnu``, VASPKIT
  ``BAND.dat``; ``#`` lines are comments);
* one column of k and one column per band (VASPKIT ``REFORMATTED_BAND.dat``
  and most plotting exports).

The high-symmetry vertices come, in order of preference, from a labels file
(Wannier90 ``*_band.labelinfo.dat``, VASPKIT ``KLABELS``), from the path
given by the user, or -- for a primitive hexagonal, tetragonal, cubic or
orthorhombic lattice -- from the segment lengths themselves: a path is
usually written with each vertex repeated, and the lengths between
vertices, compared with the distances between the lattice's
high-symmetry points, name them.

Energies are taken relative to E_F; give ``fermi_energy_eV`` if the file
holds absolute energies.
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field

import numpy as np

from tools.lattice import conventional_vectors, reciprocal_vectors

#: Spellings of Gamma in band files and plots.
GAMMA = {"G", "Γ", "GAMMA", "\\GAMMA", "$\\GAMMA$", "GM", "GAMMA_POINT"}

#: High-symmetry points in fractional coordinates of the *conventional*
#: reciprocal cell, for the primitive lattices only (a centred lattice's
#: points depend on its axis ratios; give the path for those).
POINTS = {
    "hexagonal": {"G": (0, 0, 0), "M": (0.5, 0, 0), "K": (1 / 3, 1 / 3, 0),
                  "A": (0, 0, 0.5), "L": (0.5, 0, 0.5), "H": (1 / 3, 1 / 3, 0.5)},
    "tetragonal": {"G": (0, 0, 0), "X": (0, 0.5, 0), "M": (0.5, 0.5, 0),
                   "Z": (0, 0, 0.5), "R": (0, 0.5, 0.5), "A": (0.5, 0.5, 0.5)},
    "cubic": {"G": (0, 0, 0), "X": (0, 0.5, 0), "M": (0.5, 0.5, 0), "R": (0.5, 0.5, 0.5)},
    "orthorhombic": {"G": (0, 0, 0), "X": (0.5, 0, 0), "Y": (0, 0.5, 0), "Z": (0, 0, 0.5),
                     "S": (0.5, 0.5, 0), "U": (0.5, 0, 0.5), "T": (0, 0.5, 0.5),
                     "R": (0.5, 0.5, 0.5)},
}


def canonical(label: str) -> str:
    text = str(label).strip().strip("$").replace("\\", "").upper()
    return "G" if text in {g.replace("\\", "").upper() for g in GAMMA} else text


@dataclass
class BandPath:
    k: np.ndarray                   # distance along the path, as the file has it
    energies: np.ndarray            # (n_bands, n_k), E - E_F
    vertices: list                  # [(k position, label)], labels canonical ("G" for Gamma)
    source: str = ""
    notes: list = field(default_factory=list)

    def segment(self, start: str, end: str):
        """``(k_from, k_to)`` of the first segment joining the two labels,
        in either direction, or None."""
        for (k0, a), (k1, b) in zip(self.vertices, self.vertices[1:]):
            if {canonical(a), canonical(b)} == {canonical(start), canonical(end)} and k1 > k0:
                return (k0, k1) if canonical(a) == canonical(start) else (k1, k0)
        return None

    def bands_between(self, k_from: float, k_to: float, fraction):
        """Every band at ``fraction`` (0 at ``k_from``, 1 at ``k_to``) of
        a segment, shape ``(n_bands,) + fraction.shape``."""
        lo, hi = sorted((k_from, k_to))
        inside = (self.k >= lo - 1e-9) & (self.k <= hi + 1e-9)
        ks, order = self.k[inside], np.argsort(self.k[inside], kind="stable")
        ks = ks[order]
        # a repeated vertex point would make the interpolation ambiguous;
        # keep the first of each
        keep = np.concatenate([[True], np.diff(ks) > 0])
        ks = ks[keep]
        position = k_from + (k_to - k_from) * np.asarray(fraction, dtype=float)
        out = np.empty((self.energies.shape[0],) + np.shape(position))
        for n, band in enumerate(self.energies):
            values = band[inside][order][keep]
            out[n] = np.interp(position, ks, values)
        return out


# -- reading -----------------------------------------------------------------------
def _read_numbers(path: str):
    blocks, current, widths = [], [], set()
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            text = line.split("#", 1)[0].strip()
            if not text:
                if current:
                    blocks.append(np.array(current, dtype=float))
                    current = []
                continue
            try:
                row = [float(v) for v in text.replace(",", " ").split()]
            except ValueError:
                continue                        # a header line
            widths.add(len(row))
            current.append(row)
    if current:
        blocks.append(np.array(current, dtype=float))
    return blocks, widths


def _read_labels_file(path: str) -> list:
    """Wannier90 ``labelinfo`` (label, index, k, x, y, z) or VASPKIT
    ``KLABELS`` (label, k). Returns ``[(k, label)]``."""
    out = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 2:
                continue
            numbers = []
            for value in parts[1:]:
                try:
                    numbers.append(float(value))
                except ValueError:
                    break
            if not numbers:
                continue
            k = numbers[1] if len(numbers) >= 5 else numbers[0]
            out.append((k, canonical(parts[0])))
    return out


def read_bands(path: str, *, labels=None, labels_file: str = None,
               fermi_energy_eV: float = 0.0, lattice=None) -> BandPath:
    blocks, widths = _read_numbers(path)
    if not blocks:
        raise ValueError(f"{os.path.basename(path)}: no numbers found")
    notes = []
    if widths == {2}:
        k = blocks[0][:, 0]
        if any(len(b) != len(k) for b in blocks):
            raise ValueError(f"{os.path.basename(path)}: bands of different lengths")
        energies = np.array([b[:, 1] for b in blocks])
    else:
        table = np.vstack(blocks)
        k, energies = table[:, 0], table[:, 1:].T
    energies = energies - float(fermi_energy_eV or 0.0)

    vertices = []
    if labels_file is None:
        for candidate in (path.replace("band.dat", "band.labelinfo.dat"),
                          os.path.join(os.path.dirname(path), "KLABELS")):
            if candidate != path and os.path.exists(candidate):
                labels_file = candidate
    if labels_file:
        vertices = _read_labels_file(labels_file)
        notes.append(f"vertices from {os.path.basename(labels_file)}")
    if not vertices:
        repeated = np.flatnonzero(np.diff(k) == 0)
        positions = [float(k[0])] + [float(k[i]) for i in repeated] + [float(k[-1])]
        if labels:
            names = [canonical(x) for x in labels]
            if len(names) != len(positions):
                raise ValueError(
                    f"the path has {len(names)} labels but the file has {len(positions) - 1} "
                    f"segments ({len(positions)} vertices, found from repeated k points); "
                    f"give a labels file instead")
            vertices = list(zip(positions, names))
            notes.append("vertices from repeated k points, labels as given")
        elif lattice is not None and len(positions) > 2:
            names, scale, error = label_by_lengths(np.diff(positions), lattice)
            vertices = list(zip(positions, names))
            notes.append(f"labels inferred from the segment lengths and the lattice "
                         f"({'-'.join(names)}, k units x{scale:.4g}, worst mismatch "
                         f"{100 * error:.1f}%) -- check them against the plot")
        else:
            raise ValueError("cannot tell where the high-symmetry points are: give the "
                             "path labels (e.g. [\"G\",\"M\",\"K\",\"G\",\"A\",\"L\",\"H\","
                             "\"A\"]) or a labels file")
    return BandPath(np.asarray(k, dtype=float), np.asarray(energies, dtype=float),
                    [(float(p), canonical(n)) for p, n in vertices], os.path.abspath(path),
                    notes)


def _system(lattice) -> str:
    if lattice.space_group is not None:
        from tools.spacegroups import spacegroup_info
        info = spacegroup_info(lattice.space_group)
        system, centering = info.crystal_system, info.centering
    else:
        system, centering = None, (lattice.centering or "P").upper()
        if abs(lattice.gamma - 120) < 1e-3:
            system = "hexagonal"
    if centering != "P":
        raise ValueError(f"labels cannot be inferred for a {centering}-centred lattice; "
                         f"give the path labels")
    return {"trigonal": "hexagonal"}.get(system, system)


def label_by_lengths(lengths, lattice, tolerance: float = 0.03):
    """Name the vertices of a path from its segment lengths. Tries both
    k conventions (with and without 2*pi) and every walk through the
    lattice's high-symmetry points that starts at Gamma. Returns
    ``(labels, scale, worst relative error)``."""
    system = _system(lattice)
    if system not in POINTS:
        raise ValueError(f"no table of high-symmetry points for a {system} lattice; "
                         f"give the path labels")
    b = reciprocal_vectors(conventional_vectors(lattice))
    points = {name: np.asarray(f) @ b for name, f in POINTS[system].items()}
    names = list(points)
    lengths = np.asarray(lengths, dtype=float)
    best = None
    for scale in (1.0, 2 * np.pi):
        target = lengths * scale
        # depth-first: each segment must match a distance between points
        stack = [(["G"], 0.0)]
        while stack:
            path, worst = stack.pop()
            depth = len(path) - 1
            if depth == len(target):
                if best is None or worst < best[2]:
                    best = (path, scale, worst)
                continue
            here = points[path[-1]]
            for name in names:
                if name == path[-1]:
                    continue
                d = float(np.linalg.norm(points[name] - here))
                err = abs(d - target[depth]) / max(target[depth], 1e-9)
                if err <= tolerance:
                    stack.append((path + [name], max(worst, err)))
    if best is None:
        raise ValueError("the segment lengths match no path through the lattice's "
                         "high-symmetry points: the lattice given differs from the "
                         "calculation's, or the path is unusual -- give the labels")
    return best


def normal_segment(bands: BandPath, surface: dict, lattice) -> dict:
    """The segment from Gamma along the surface normal to the zone boundary:
    its two ends in the file's k, which end is Gamma, and the factor from
    the file's k units to 1/A (checked against half the surface period)."""
    half = surface["period_invA"] / 2.0
    best = None
    for (k0, a), (k1, b) in zip(bands.vertices, bands.vertices[1:]):
        if "G" not in (a, b) or k1 <= k0:
            continue
        for scale in (1.0, 2 * np.pi):
            err = abs((k1 - k0) * scale - half) / half
            if best is None or err < best[0]:
                best = (err, (k0, k1) if a == "G" else (k1, k0), scale, b if a == "G" else a)
    if best is None or best[0] > 0.05:
        found = f" (closest: {100 * best[0]:.0f}% off)" if best else ""
        raise ValueError(f"no segment from Gamma is half the surface period "
                         f"({half:.4f} 1/A){found}: the calculation's path does not run "
                         f"along this surface normal, or its lattice differs")
    err, (k_gamma, k_edge), scale, edge = best
    return {"k_gamma": k_gamma, "k_edge": k_edge, "scale": scale, "edge_label": edge,
            "length_error": err}


# -- matching ------------------------------------------------------------------------
A_CONST = 0.2624682843         # 2 m_e / hbar^2 in 1/(A^2 eV)


def normal_emission_image(cube, slit, energy, angle_offset=0.0, halfwidth_deg=1.0):
    """(hv, E) intensity at normal emission, as -d2I/dE2 (peaks, not
    backgrounds), clipped at zero and scaled per photon energy."""
    from scipy.ndimage import gaussian_filter1d
    near = np.abs(np.asarray(slit) - angle_offset) <= halfwidth_deg
    if not near.any():
        near[int(np.argmin(np.abs(np.asarray(slit) - angle_offset)))] = True
    with np.errstate(invalid="ignore"):
        image = np.nanmean(np.asarray(cube, dtype=float)[:, near, :], axis=1)
    image = np.nan_to_num(image)
    step = abs(float(energy[1] - energy[0])) if len(energy) > 1 else 0.01
    image = gaussian_filter1d(image, max(0.03 / step, 1.0), axis=1)
    curvature = -np.gradient(np.gradient(image, axis=1), axis=1)
    curvature = np.clip(curvature, 0, None)
    scale = curvature.max(axis=1, keepdims=True)
    return np.where(scale > 0, curvature / np.where(scale > 0, scale, 1), 0.0)


def kz_of(hv, energy, *, v0, work_function, effective_mass=1.0):
    kinetic = np.asarray(hv)[:, None] - work_function + np.asarray(energy)[None, :]
    return np.sqrt(np.clip(A_CONST * effective_mass * (kinetic + v0), 0, None))


def simulated_image(bands: BandPath, segment: dict, period: float, hv, energy, *,
                    v0, shift, work_function, renormalisation=1.0, width=0.06,
                    effective_mass=1.0):
    """What the calculation predicts at normal emission: Gaussians at the
    calculated energies (shifted, renormalised), occupied states only."""
    kz = kz_of(hv, energy, v0=v0, work_function=work_function,
               effective_mass=effective_mass)
    reduced = np.mod(kz, period)
    reduced = np.minimum(reduced, period - reduced)          # E(kz) = E(-kz)
    fraction = reduced / (period / 2.0)
    values = bands.bands_between(segment["k_gamma"], segment["k_edge"], fraction)
    values = values / float(renormalisation) + float(shift)
    grid = np.asarray(energy)[None, None, :]
    near = np.abs(values - grid) < 4 * width
    image = np.where(near, np.exp(-0.5 * ((values - grid) / width) ** 2), 0.0).sum(axis=0)
    occupied = 1.0 / (np.exp(np.asarray(energy) / 0.01) + 1.0)
    return image * occupied[None, :]


def _correlation(a, b):
    a = a - a.mean()
    b = b - b.mean()
    denominator = np.sqrt(np.sum(a * a) * np.sum(b * b))
    return float(np.sum(a * b) / denominator) if denominator else float("nan")


def match(cube, hv, slit, energy, bands: BandPath, segment: dict, period: float, *,
          work_function, angle_offset=0.0, v0_values=None, shifts=None,
          renormalisations=(1.0,), energy_window=(-2.0, 0.05), width=0.06):
    """Scan V0, the calculation's rigid energy shift and (optionally) its
    renormalisation for the best agreement with the scan at normal emission.
    Returns a dict with the best values, the score map and the photon
    energies of the Gamma and zone-boundary planes."""
    hv = np.asarray(hv, dtype=float)
    energy = np.asarray(energy, dtype=float)
    keep = (energy >= energy_window[0]) & (energy <= energy_window[1])
    e = energy[keep]
    measured = normal_emission_image(cube, slit, energy, angle_offset)[:, keep]
    v0_values = np.arange(2.0, 30.01, 0.5) if v0_values is None else np.asarray(v0_values)
    shifts = np.arange(-0.3, 0.301, 0.025) if shifts is None else np.asarray(shifts)
    scores = np.full((len(renormalisations), v0_values.size, shifts.size), np.nan)
    for r, z in enumerate(renormalisations):
        for i, v0 in enumerate(v0_values):
            for j, shift in enumerate(shifts):
                sim = simulated_image(bands, segment, period, hv, e, v0=v0, shift=shift,
                                      work_function=work_function, renormalisation=z,
                                      width=width)
                scores[r, i, j] = _correlation(measured, sim)
    r, i, j = np.unravel_index(int(np.nanargmax(scores)), scores.shape)
    v0, shift, z = float(v0_values[i]), float(shifts[j]), float(renormalisations[r])
    # How sharply the score picks V0: the range within 0.02 of the best
    # (at the best shift) -- a practical, not a statistical, uncertainty.
    profile = scores[r, :, j]
    good = v0_values[profile >= np.nanmax(profile) - 0.02]
    return {"inner_potential_eV": v0, "energy_shift_eV": shift, "renormalisation": z,
            "score": float(scores[r, i, j]),
            "v0_range_eV": [float(good.min()), float(good.max())],
            "planes": planes(hv, period, v0=v0, work_function=work_function,
                             edge_label=segment.get("edge_label", "boundary")),
            "v0_values": v0_values.tolist(), "shifts": shifts.tolist(),
            "scores": scores[r].tolist()}


def planes(hv, period, *, v0, work_function, binding_energy=0.0, edge_label="boundary"):
    """The photon energies inside the scan at which normal emission reaches
    a Gamma plane (k_z = n G) or a zone-boundary plane (k_z = (n + 1/2) G)."""
    hv = np.asarray(hv, dtype=float)
    lo, hi = float(hv.min()), float(hv.max())
    out = []
    n_max = int(np.sqrt(A_CONST * (hi - work_function + v0 + binding_energy)) / period) + 2
    for n in range(0, n_max + 1):
        for plane, kz in (("G", n * period), (edge_label, (n + 0.5) * period)):
            energy = kz ** 2 / A_CONST - v0 + work_function - binding_energy
            if lo <= energy <= hi:
                out.append({"plane": plane, "order": n, "k_z_invA": float(kz),
                            "hv_eV": float(energy)})
    return sorted(out, key=lambda p: p["hv_eV"])


def in_plane_directions(bands: BandPath, centre: str) -> list:
    """The segments that leave ``centre`` (e.g. Gamma-M, Gamma-K), for
    overlaying on a cut taken in that plane."""
    out = []
    for (k0, a), (k1, b) in zip(bands.vertices, bands.vertices[1:]):
        if k1 <= k0:
            continue
        if canonical(a) == canonical(centre):
            out.append((f"{a}-{b}", k0, k1))
        elif canonical(b) == canonical(centre):
            out.append((f"{b}-{a}", k1, k0))
    seen, unique = set(), []
    for name, k0, k1 in out:
        if name not in seen:
            seen.add(name)
            unique.append((name, k0, k1))
    return unique

