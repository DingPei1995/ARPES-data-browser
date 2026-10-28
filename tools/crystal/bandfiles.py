"""
tools/crystal/bandfiles.py
==========================
Calculated band structures read from (and written to) the files the
calculation codes leave behind. Three families, as agreed:

* **XCrySDen** ``.bxsf`` -- eigenvalues on a 3-D grid over the reciprocal
  cell. Written by Wien2k (through XCrySDen), Quantum ESPRESSO, VASP
  (vaspkit) and Wannier90; read as a :class:`~tools.crystal.bands.BandGrid`.
* **Wien2k** and **Wannier90** -- Wien2k's band structure along a path
  (``case.spaghetti_ene``, with the labels of ``case.klist_band``);
  Wannier90's Hamiltonian (``seedname_hr.dat``, with the cell from
  ``seedname.win``), evaluated exactly anywhere, and its band structure
  along a path (``seedname_band.dat`` and ``..._band.labelinfo.dat``).
* **VASP** -- ``EIGENVAL`` with the cell from ``POSCAR`` (and E_F from
  ``OUTCAR`` or ``DOSCAR``): along a path when ``KPOINTS`` is in line mode
  (labels from it), on a grid when the k points fill the whole zone; with
  the orbital / atomic character of ``PROCAR`` as weights.

:func:`read` looks at the file name and content and calls the right reader;
the loaders were written against the formats as their codes document them
(and the tests check each against a file written to that layout).

Units: a ``.bxsf`` gives its reciprocal vectors as the writing code chose --
Wannier90 and VASP in Å⁻¹ with 2π included, Quantum ESPRESSO's ``fs.x`` in
units of 2π/alat. ``k_scale`` multiplies them (e.g. 2π, or 2π/a); the
window compares the result with the crystal's own zone and says when they
disagree.

No Qt.
"""
from __future__ import annotations

import os
import re

import numpy as np

from tools.crystal.bands import BandGrid, BandPath, WannierModel, reciprocal_of

__all__ = ["companion_files", "read", "read_bxsf", "write_bxsf", "read_hr", "write_hr",
           "read_win_cell", "read_wannier_band", "read_spaghetti",
           "read_eigenval", "read_poscar", "read_procar", "FORMATS"]

BOHR = 0.529177210903

#: what :func:`read` recognises, for the file dialog: (label, patterns)
FORMATS = (
    ("XCrySDen band grid (*.bxsf)", ("*.bxsf",)),
    ("Wannier90 Hamiltonian (*_hr.dat)", ("*_hr.dat",)),
    ("Wannier90 band structure (*_band.dat)", ("*_band.dat",)),
    ("Wien2k band structure (*.spaghetti_ene)", ("*.spaghetti_ene",)),
    ("VASP EIGENVAL", ("EIGENVAL*",)),
)


def companion_files(path: str) -> list:
    """The files beside ``path`` that reading it also uses (the cell, E_F,
    labels, weights): a ``_hr.dat``'s ``.win``, a ``_band.dat``'s
    ``labelinfo``, an ``EIGENVAL``'s POSCAR, PROCAR, OUTCAR, DOSCAR and
    KPOINTS, a ``.spaghetti_ene``'s ``.klist_band``. Only those that exist."""
    folder, name = os.path.split(os.path.abspath(path))
    lower = name.lower()
    wanted = []
    if lower.endswith("_hr.dat"):
        wanted = [name[:-len("_hr.dat")] + ".win"]
    elif lower.endswith("_band.dat"):
        wanted = [name[:-len(".dat")] + ".labelinfo.dat",
                  name[:-len("_band.dat")] + ".win"]
    elif lower.endswith("spaghetti_ene"):
        stem = name[:-len(".spaghetti_ene")] if lower.endswith(".spaghetti_ene") else name
        wanted = [stem + ".klist_band"]
    elif lower.startswith("eigenval"):
        wanted = ["POSCAR", "CONTCAR", "PROCAR", "OUTCAR", "DOSCAR", "KPOINTS"]
    return [os.path.join(folder, w) for w in wanted
            if os.path.isfile(os.path.join(folder, w))]


def read(path: str, **options):
    """Whatever calculation ``path`` holds."""
    name = os.path.basename(path).lower()
    if name.endswith(".bxsf"):
        return read_bxsf(path, **options)
    if name.endswith("_hr.dat"):
        return read_hr(path, **options)
    if name.endswith("_band.dat"):
        return read_wannier_band(path, **options)
    if name.endswith(".spaghetti_ene") or name.endswith("spaghetti_ene"):
        return read_spaghetti(path, **options)
    if name.startswith("eigenval"):
        return read_eigenval(path, **options)
    raise ValueError(f"not a band-structure file this program reads: "
                     f"{os.path.basename(path)}")


# --------------------------------------------------------------------------
# XCrySDen .bxsf
# --------------------------------------------------------------------------
def read_bxsf(path: str, k_scale: float = 1.0, fermi=None):
    """A ``.bxsf`` band grid. Its grid is "general" (the last point along
    each axis repeats the first), and the values are listed with the third
    index fastest, as XCrySDen documents and Wannier90 writes them."""
    with open(path, encoding="utf-8", errors="replace") as handle:
        text = handle.read()
    match = re.search(r"Fermi\s+Energy\s*:\s*([-+0-9.eEdD]+)", text)
    file_fermi = float(match.group(1).replace("D", "E")) if match else 0.0
    block = re.search(r"BEGIN_BANDGRID_3D\w*(.*?)END_BANDGRID_3D", text, re.S | re.I)
    if block is None:
        raise ValueError("no BANDGRID_3D block in the file")
    lines = [line.strip() for line in block.group(1).splitlines() if line.strip()]
    nbands = int(lines[0].split()[0])
    n = [int(v) for v in lines[1].split()[:3]]
    origin = np.array([float(v) for v in lines[2].split()[:3]])
    vectors = np.array([[float(v) for v in lines[3 + i].split()[:3]] for i in range(3)])
    rest = lines[6:]
    values = []
    position = 0
    for _band in range(nbands):
        while position < len(rest) and not rest[position].upper().startswith("BAND"):
            position += 1
        position += 1
        numbers = []
        while position < len(rest) and not rest[position].upper().startswith(("BAND", "END")):
            numbers.extend(float(v) for v in rest[position].split())
            position += 1
        values.append(np.asarray(numbers[:n[0] * n[1] * n[2]]).reshape(n))
    reciprocal = vectors * float(k_scale)
    origin_frac = origin @ np.linalg.inv(vectors) if np.any(origin) else np.zeros(3)
    return BandGrid.from_general_grid(
        reciprocal, np.array(values), origin=origin_frac,
        fermi=file_fermi if fermi is None else float(fermi),
        name=os.path.basename(path))


def write_bxsf(path: str, model, n=(24, 24, 12), fermi: float = 0.0):
    """``model`` (anything with ``reciprocal`` and ``energies``) sampled on
    a general grid and written as a ``.bxsf`` -- energies absolute, so
    ``fermi`` is added back."""
    b = np.asarray(model.reciprocal, dtype=float)
    axes = [np.linspace(0.0, 1.0, m) for m in n]
    grid = np.stack(np.meshgrid(*axes, indexing="ij"), -1).reshape(-1, 3)
    values = model.energies(grid @ b) + fermi
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("BEGIN_INFO\n  # written by ARPES viewer\n"
                     f"  Fermi Energy: {fermi:.6f}\nEND_INFO\n")
        handle.write("BEGIN_BLOCK_BANDGRID_3D\n band_energies\n BEGIN_BANDGRID_3D\n")
        handle.write(f" {values.shape[0]}\n {n[0]} {n[1]} {n[2]}\n 0.0 0.0 0.0\n")
        for row in b:
            handle.write(" " + " ".join(f"{v:.8f}" for v in row) + "\n")
        for band, energies in enumerate(values, start=1):
            handle.write(f" BAND: {band}\n")
            for start in range(0, energies.size, 6):
                handle.write(" " + " ".join(f"{v:.6f}" for v in energies[start:start + 6]) + "\n")
        handle.write(" END_BANDGRID_3D\nEND_BLOCK_BANDGRID_3D\n")
    return path


# --------------------------------------------------------------------------
# Wannier90
# --------------------------------------------------------------------------
def read_win_cell(path: str):
    """The cell (rows, Å) and E_F of a Wannier90 ``.win`` file, or
    ``(None, None)``."""
    try:
        with open(path, encoding="utf-8", errors="replace") as handle:
            text = handle.read()
    except OSError:
        return None, None
    fermi = re.search(r"^\s*fermi_energy\s*[=:\s]\s*([-+0-9.eE]+)", text, re.M | re.I)
    fermi = float(fermi.group(1)) if fermi else None
    block = re.search(r"begin\s+unit_cell_cart(.*?)end\s+unit_cell_cart", text,
                      re.S | re.I)
    if block is None:
        return None, fermi
    lines = [l.split() for l in block.group(1).splitlines() if l.strip()]
    unit = 1.0
    if lines and len(lines[0]) == 1:
        unit = BOHR if lines[0][0].lower().startswith("bohr") else 1.0
        lines = lines[1:]
    cell = np.array([[float(v) for v in l[:3]] for l in lines[:3]]) * unit
    return cell, fermi


def read_hr(path: str, lattice=None, fermi=None):
    """A Wannier90 ``seedname_hr.dat``. The cell comes from ``lattice``
    (rows, Å), else from ``seedname.win`` beside it; E_F likewise from the
    ``.win`` (``fermi_energy``) unless given."""
    with open(path, encoding="utf-8", errors="replace") as handle:
        handle.readline()
        nw = int(handle.readline().split()[0])
        nr = int(handle.readline().split()[0])
        degeneracy = []
        while len(degeneracy) < nr:
            degeneracy.extend(int(v) for v in handle.readline().split())
        rows = np.loadtxt(handle, ndmin=2)
    R_all = rows[:, :3].astype(int)
    m = rows[:, 3].astype(int) - 1
    n = rows[:, 4].astype(int) - 1
    values = rows[:, 5] + 1j * rows[:, 6]
    keys, index = np.unique(R_all, axis=0, return_inverse=True)
    H = np.zeros((len(keys), nw, nw), dtype=complex)
    H[index.ravel(), m, n] = values
    # degeneracy is listed in the file's order of R; map it to the keys
    order = []
    seen = {}
    for r in map(tuple, R_all):
        if r not in seen:
            seen[r] = len(order)
            order.append(r)
    deg_by_r = {r: degeneracy[i] for i, r in enumerate(order)}
    deg = np.array([deg_by_r[tuple(r)] for r in keys], dtype=float)
    seed = path[:-len("_hr.dat")]
    win_cell, win_fermi = read_win_cell(seed + ".win")
    cell = lattice if lattice is not None else win_cell
    if cell is None:
        raise ValueError("the cell is needed: give it (or the Brillouin-zone "
                         "tool's crystal), or put seedname.win beside the file")
    return WannierModel(cell, keys, H, degeneracy=deg,
                        fermi=(fermi if fermi is not None else (win_fermi or 0.0)),
                        orbital_names=_win_orbitals(seed + ".win", nw),
                        name=os.path.basename(path))


def _win_orbitals(path: str, count: int):
    """The Wannier functions' names from a ``.win``'s "! orbitals: ..."
    line (what :func:`write_hr` writes), if it names ``count`` of them."""
    try:
        with open(path, encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.strip().lower().startswith("! orbitals:"):
                    names = line.split(":", 1)[1].split()
                    return names if len(names) == count else None
    except OSError:
        pass
    return None


def write_hr(path: str, lattice, hoppings, n_orbitals: int, fermi=None,
             orbital_names=None):
    """A Wannier90 ``_hr.dat`` from ``hoppings`` = ``{(R1, R2, R3): (nw, nw)
    complex}`` (degeneracy 1), and a ``.win`` beside it with the cell (and
    ``fermi_energy`` when given)."""
    keys = sorted(hoppings)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("written by ARPES viewer\n")
        handle.write(f"{n_orbitals:12d}\n{len(keys):12d}\n")
        degeneracy = ["    1"] * len(keys)
        for start in range(0, len(keys), 15):
            handle.write("".join(degeneracy[start:start + 15]) + "\n")
        for key in keys:
            block = np.asarray(hoppings[key], dtype=complex)
            for n in range(n_orbitals):
                for m in range(n_orbitals):
                    v = block[m, n]
                    handle.write(f"{key[0]:5d}{key[1]:5d}{key[2]:5d}{m + 1:5d}"
                                 f"{n + 1:5d}{v.real:12.6f}{v.imag:12.6f}\n")
    seed = path[:-len("_hr.dat")] if path.endswith("_hr.dat") else os.path.splitext(path)[0]
    with open(seed + ".win", "w", encoding="utf-8") as handle:
        handle.write(f"num_wann = {n_orbitals}\n")
        if fermi is not None:
            handle.write(f"fermi_energy = {fermi:.6f}\n")
        handle.write("begin unit_cell_cart\nang\n")
        for row in np.asarray(lattice, dtype=float):
            handle.write(" " + " ".join(f"{v:.8f}" for v in row) + "\n")
        handle.write("end unit_cell_cart\n")
        if orbital_names:
            handle.write("! orbitals: " + " ".join(orbital_names) + "\n")
    return path


def read_wannier_band(path: str, fermi=None):
    """Wannier90's ``seedname_band.dat`` (distance and energy, one band per
    block), with the labels of ``seedname_band.labelinfo.dat`` if present."""
    blocks, current = [], []
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = line.split()
            if len(parts) >= 2:
                current.append((float(parts[0]), float(parts[1])))
            elif current:
                blocks.append(current)
                current = []
    if current:
        blocks.append(current)
    distance = np.array([d for d, _ in blocks[0]])
    energies = np.array([[e for _, e in block] for block in blocks
                         if len(block) == len(distance)])
    labels = []
    info = path[:-len(".dat")] + ".labelinfo.dat"
    if os.path.exists(info):
        with open(info, encoding="utf-8") as handle:
            for line in handle:
                parts = line.split()
                if len(parts) >= 3:
                    labels.append((float(parts[2]), _greek(parts[0])))
    seed = path[:-len("_band.dat")]
    _cell, win_fermi = read_win_cell(seed + ".win")
    return BandPath(distance, energies, labels,
                    fermi=fermi if fermi is not None else (win_fermi or 0.0),
                    name=os.path.basename(path))


def _greek(label: str) -> str:
    return {"G": "Γ", "GAMMA": "Γ", "GM": "Γ", "\\GAMMA": "Γ", "SIGMA": "Σ",
            "DELTA": "Δ", "LAMBDA": "Λ"}.get(label.strip().upper(), label.strip())


# --------------------------------------------------------------------------
# Wien2k
# --------------------------------------------------------------------------
def read_spaghetti(path: str, fermi: float = 0.0):
    """Wien2k's ``case.spaghetti_ene``: per band a ``bandindex:`` line and
    rows of kx ky kz, the path distance and the energy (eV, already
    relative to E_F as ``spaghetti`` writes it). The labels come from
    ``case.klist_band`` beside it, if present."""
    bands, current = [], []
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if "bandindex" in line.lower():
                if current:
                    bands.append(current)
                current = []
                continue
            parts = line.split()
            if len(parts) >= 5:
                try:
                    current.append([float(v) for v in parts[:5]])
                except ValueError:
                    continue
    if current:
        bands.append(current)
    if not bands:
        raise ValueError("no bands in the file")
    first = np.array(bands[0])
    distance = first[:, 3]
    energies = np.array([np.array(b)[:, 4] for b in bands if len(b) == len(first)])
    labels = []
    klist = os.path.splitext(path)[0] + ".klist_band"
    if os.path.exists(klist):
        with open(klist, encoding="utf-8", errors="replace") as handle:
            for index, line in enumerate(handle):
                if line.startswith("END"):
                    break
                name = line[:10].strip()
                if name and index < len(distance):
                    labels.append((float(distance[index]), _greek(name)))
    return BandPath(distance, energies, labels, k=first[:, :3], fermi=fermi,
                    name=os.path.basename(path))


# --------------------------------------------------------------------------
# VASP
# --------------------------------------------------------------------------
def read_poscar(path: str):
    """``(lattice rows Å, species list, counts list)`` of a POSCAR/CONTCAR."""
    with open(path, encoding="utf-8", errors="replace") as handle:
        lines = handle.read().splitlines()
    scale = float(lines[1].split()[0])
    cell = np.array([[float(v) for v in lines[2 + i].split()[:3]] for i in range(3)])
    if scale < 0:
        cell *= (abs(scale) / abs(np.linalg.det(cell))) ** (1 / 3)
    else:
        cell *= scale
    species, counts = [], []
    tokens = lines[5].split()
    if tokens and not tokens[0].isdigit():
        species = tokens
        counts = [int(v) for v in lines[6].split()]
    else:
        counts = [int(v) for v in tokens]
    return cell, species, counts


def _vasp_fermi(folder: str):
    try:
        with open(os.path.join(folder, "OUTCAR"), encoding="utf-8",
                  errors="replace") as handle:
            found = None
            for line in handle:
                if "E-fermi" in line:
                    found = float(line.split(":")[1].split()[0])
            if found is not None:
                return found
    except (OSError, ValueError, IndexError):
        pass
    try:
        with open(os.path.join(folder, "DOSCAR"), encoding="utf-8") as handle:
            return float(handle.read().splitlines()[5].split()[3])
    except (OSError, ValueError, IndexError):
        return None


def _kpoints_labels(folder: str):
    """Line-mode labels of a KPOINTS file: ``(divisions, [labels...])`` with
    each segment's two end labels, or None."""
    try:
        with open(os.path.join(folder, "KPOINTS"), encoding="utf-8") as handle:
            lines = handle.read().splitlines()
    except OSError:
        return None
    if len(lines) < 4 or not lines[2].strip().lower().startswith("l"):
        return None
    divisions = int(lines[1].split()[0])
    labels = []
    for line in lines[4:]:
        if not line.strip():
            continue
        label = line.split("!")[1].strip() if "!" in line else ""
        labels.append(_greek(label))
    return divisions, labels


def read_eigenval(path: str, lattice=None, fermi=None, procar: bool = True):
    """VASP's ``EIGENVAL``. Along a path (``KPOINTS`` in line mode, labels
    taken from it) or on a grid (the k points fill the zone); E_F from
    OUTCAR/DOSCAR; the cell from POSCAR beside it; ``PROCAR`` weights if
    the file is there."""
    folder = os.path.dirname(os.path.abspath(path))
    with open(path, encoding="utf-8", errors="replace") as handle:
        lines = handle.read().splitlines()
    _nelect, nk, nb = (int(v) for v in lines[5].split()[:3])
    k_frac, energies = [], []
    position = 6
    for _k in range(nk):
        while not lines[position].strip():
            position += 1
        k_frac.append([float(v) for v in lines[position].split()[:3]])
        position += 1
        row = []
        for _b in range(nb):
            row.append(float(lines[position].split()[1]))
            position += 1
        energies.append(row)
    k_frac = np.array(k_frac)
    energies = np.array(energies).T                      # (bands, nk)
    if lattice is None:
        try:
            lattice, species, counts = read_poscar(os.path.join(folder, "POSCAR"))
        except OSError:
            raise ValueError("the cell is needed: put POSCAR beside EIGENVAL") from None
    else:
        species, counts = [], []
    b = reciprocal_of(lattice)
    k_cart = k_frac @ b
    file_fermi = _vasp_fermi(folder)
    fermi = fermi if fermi is not None else (file_fermi or 0.0)
    weights = {}
    if procar and os.path.exists(os.path.join(folder, "PROCAR")):
        try:
            weights = read_procar(os.path.join(folder, "PROCAR"), species, counts)
        except (ValueError, IndexError):
            weights = {}
    line_mode = _kpoints_labels(folder)
    if line_mode is not None:
        divisions, names = line_mode
        steps = np.r_[0.0, np.linalg.norm(np.diff(k_cart, axis=0), axis=1)]
        # a jump between segments (the end of one is not the start of the
        # next) adds no distance
        for start in range(divisions, nk, divisions):
            steps[start] = 0.0
        distance = np.cumsum(steps)
        labels = []
        for segment in range(len(names) // 2):
            a, b_ = names[2 * segment], names[2 * segment + 1]
            d0 = distance[segment * divisions]
            d1 = distance[min(nk - 1, (segment + 1) * divisions - 1)]
            if labels and abs(labels[-1][0] - d0) < 1e-9 and labels[-1][1] != a:
                labels[-1] = (d0, f"{labels[-1][1]}|{a}")
            elif not labels or abs(labels[-1][0] - d0) > 1e-9:
                labels.append((d0, a))
            labels.append((d1, b_))
        return BandPath(distance, energies, labels, k=k_cart, fermi=fermi,
                        weights={k: v for k, v in weights.items()},
                        name=os.path.basename(path))
    # a grid: the fractional coordinates must fill a regular mesh
    axes = [np.unique(np.round(np.mod(k_frac[:, i], 1.0), 6)) for i in range(3)]
    n = [len(a) for a in axes]
    if int(np.prod(n)) != nk:
        raise ValueError("EIGENVAL holds a k mesh reduced by symmetry; run "
                         "with ISYM = 0 for the full zone, or export a .bxsf")
    index = [np.searchsorted(axes[i], np.round(np.mod(k_frac[:, i], 1.0), 6))
             for i in range(3)]
    values = np.zeros((nb, *n))
    values[:, index[0], index[1], index[2]] = energies
    grid_weights = {}
    for key, w in weights.items():
        arr = np.zeros((nb, *n))
        arr[:, index[0], index[1], index[2]] = w
        grid_weights[key] = arr
    return BandGrid(b, values, origin=[a[0] for a in axes], fermi=fermi,
                    weights=grid_weights, name=os.path.basename(path))


def read_procar(path: str, species=(), counts=()):
    """Orbital and atomic character from a ``PROCAR``: weights named
    ``"<element> <s/p/d/f>"`` (or ``"ion<i> ..."`` without POSCAR species),
    each ``(bands, k points)``, normalised by the band's total. A
    spin-orbit file's magnetisation blocks are skipped."""
    with open(path, encoding="utf-8", errors="replace") as handle:
        text = handle.read()
    header = re.search(r"# of k-points:\s*(\d+)\s+# of bands:\s*(\d+)\s+# of ions:\s*(\d+)", text)
    if header is None:
        raise ValueError("not a PROCAR")
    nk, nb, nions = (int(v) for v in header.groups())
    owner = []
    for element, count in zip(species, counts):
        owner.extend([element] * count)
    while len(owner) < nions:
        owner.append(f"ion{len(owner) + 1}")
    sums = {}
    k_index = -1
    band_index = -1
    columns = None
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.strip().startswith("k-point"):
            k_index += 1
        elif line.strip().startswith("band"):
            band_index = int(line.split()[1]) - 1
        elif line.strip().startswith("ion") and band_index >= 0:
            columns = line.split()[1:-1]
            rows = []
            for j in range(1, nions + 1):
                rows.append([float(v) for v in lines[i + j].split()[1:-1]])
            i += nions + 1
            # the first table after "band" is the total; skip the rest
            totals = np.array(rows)
            for ion, values in enumerate(totals):
                for name, value in zip(columns, values):
                    orbital = name[0]
                    key = f"{owner[ion]} {orbital}"
                    sums.setdefault(key, np.zeros((nb, nk)))[band_index, k_index] += value
            band_index = -1
            continue
        i += 1
    total = sum(sums.values())
    total = np.where(total > 0, total, 1.0)
    return {key: value / total for key, value in sorted(sums.items())}
