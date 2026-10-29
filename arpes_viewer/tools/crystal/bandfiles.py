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
* **FPLO + Wannier tools** -- what the post-processing programs that
  evaluate an FPLO Wannier Hamiltonian (``wannier90_hr_fplo.dat``) leave
  behind (``k_path.exe``, ``H.exe``, ``plot.exe``; MATLAB scripts plot
  them): ``band.dat`` along a path (distance and energy, a blank line
  between bands; a third column -- ``orbital.dat`` -- is a weight) or on a
  3-D grid (``k1 k2 k3 E1 E2 ...`` per row, fractional k from -1/2 to 1/2),
  ``Eigenvalue.dat`` (per k a line of fractional k, then one energy per
  line), the path and the cell in their own ``KPOINTS`` (first line
  "segments points-per-segment", then the direct and reciprocal vectors,
  then the segments' end points with "! label"), and the Wannier
  functions' names in ``wannier_name_center.txt``. The names are generic,
  so :func:`read` tells these apart by their content.

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

from arpes_viewer.tools.crystal.bands import BandGrid, BandPath, WannierModel, reciprocal_of

__all__ = ["companion_files", "read", "read_bxsf", "write_bxsf", "read_hr", "write_hr",
           "read_win_cell", "read_wannier_band", "read_spaghetti",
           "read_eigenval", "read_poscar", "read_procar", "read_fplo_kpoints",
           "read_fplo_band", "read_fplo_eigenvalue", "read_fplo_grid",
           "read_matlab_grid", "detect", "survey", "split_groups", "join_pieces",
           "formats_for", "GOALS", "FORMATS"]

BOHR = 0.529177210903

#: what the user wants of a calculation, and what it takes:
#: "dispersion" (E against k: along a path, or along a measured cut),
#: "contours" (constant-energy contours) and "fs" (the 3-D Fermi surface)
#: -- the last two need the bands everywhere in the zone
GOALS = ("dispersion", "contours", "fs")
_EVERYWHERE = frozenset(GOALS)
_PATH_ONLY = frozenset({"dispersion"})

#: what :func:`read` recognises, for the file dialog: (label, patterns,
#: the goals it serves)
FORMATS = (
    ("XCrySDen band grid (*.bxsf)", ("*.bxsf",), _EVERYWHERE),
    ("Wannier90 Hamiltonian (*_hr*.dat)", ("*_hr*.dat",), _EVERYWHERE),
    ("3-D band grid of the FPLO tools (band.dat)", ("*.dat",), _EVERYWHERE),
    ("MATLAB band grid (*.mat)", ("*.mat",), _EVERYWHERE),
    ("VASP EIGENVAL", ("EIGENVAL*",), _EVERYWHERE),
    ("Wannier90 band structure (*_band.dat)", ("*_band.dat",), _PATH_ONLY),
    ("Wien2k band structure (*.spaghetti_ene)", ("*.spaghetti_ene",), _PATH_ONLY),
    ("FPLO tools along a path (band.dat, Eigenvalue.dat)", ("*.dat",), _PATH_ONLY),
)


def formats_for(goal: str = None):
    """The file-dialog filters for ``goal`` (all formats when None): first
    one that shows every file that can serve it, then each format."""
    chosen = [(label, patterns) for label, patterns, goals in FORMATS
              if goal is None or goal in goals]
    patterns = []
    for _label, pats in chosen:
        patterns += [p for p in pats if p not in patterns]
    return [("Calculations for this (" + " ".join(patterns) + ")", tuple(patterns))] + chosen


#: each kind :func:`detect` tells apart: (in plain words, the goals it
#: serves -- empty for a file that only helps another one)
KINDS = {
    "bxsf": ("bands everywhere in the zone (XCrySDen band grid)", _EVERYWHERE),
    "hr": ("a Wannier / tight-binding Hamiltonian (bands anywhere, exactly)",
           _EVERYWHERE),
    "fplo_grid": ("bands everywhere in the zone (3-D grid: k1 k2 k3 E1 E2 ...)",
                  _EVERYWHERE),
    "matlab": ("bands everywhere in the zone (MATLAB array Ek_3D)", _EVERYWHERE),
    "eigenval": ("VASP bands, along a path or on a grid", _EVERYWHERE),
    "wannier_band": ("bands along a path (Wannier90)", _PATH_ONLY),
    "spaghetti": ("bands along a path (Wien2k)", _PATH_ONLY),
    "fplo_path": ("bands along a path (distance, energy[, weight])", _PATH_ONLY),
    "fplo_eigenvalue": ("bands along a path (per k point: k, then energies)",
                        _PATH_ONLY),
    "win": ("helper: the cell and E_F of a Wannier90 run (.win)", frozenset()),
    "fplo_kpoints": ("helper: the cell and the path's labels (KPOINTS of the "
                     "FPLO tools)", frozenset()),
    "kpoints": ("helper: the k path of a VASP run (KPOINTS)", frozenset()),
    "poscar": ("helper: the cell of a VASP run (POSCAR)", frozenset()),
    "procar": ("helper: orbital weights of a VASP run (PROCAR)", frozenset()),
    "labelinfo": ("helper: the path's labels (Wannier90 labelinfo)", frozenset()),
    "klist": ("helper: the path's labels (Wien2k klist_band)", frozenset()),
    "wannier_names": ("helper: the names of the Wannier functions", frozenset()),
    "kfile": ("helper: the k points (or distances) of a path", frozenset()),
}

#: the name each kind usually has (for "renamed" notes): a test on the
#: lower-case name
_USUAL_NAME = {
    "bxsf": lambda n: n.endswith(".bxsf"),
    "hr": lambda n: bool(_HR.search(n)),
    "wannier_band": lambda n: n.endswith("_band.dat"),
    "spaghetti": lambda n: n.endswith("spaghetti_ene"),
    "eigenval": lambda n: n.startswith("eigenval"),
    "matlab": lambda n: n.endswith(".mat"),
    "fplo_grid": lambda n: n.endswith(".dat"),
    "fplo_path": lambda n: n.endswith(".dat"),
    "fplo_eigenvalue": lambda n: n.endswith(".dat"),
}

#: ``wannier90_hr.dat``, ``wannier90_hr_fplo.dat``, ``x_hr.dat``...
_HR = re.compile(r"_hr(_[^./\\]*)?\.dat$", re.I)



def companion_files(path: str) -> list:
    """The files beside ``path`` that reading it also uses (the cell, E_F,
    labels, weights): a ``_hr.dat``'s ``.win``, a ``_band.dat``'s
    ``labelinfo``, an ``EIGENVAL``'s POSCAR, PROCAR, OUTCAR, DOSCAR and
    KPOINTS, a ``.spaghetti_ene``'s ``.klist_band``. Only those that exist."""
    folder, name = os.path.split(os.path.abspath(path))
    lower = name.lower()
    wanted = []
    if _HR.search(name):
        wanted = [name[:_HR.search(name).start()] + ".win", "KPOINTS",
                  "wannier_name_center.txt"]
    elif lower.endswith("_band.dat"):
        wanted = [name[:-len(".dat")] + ".labelinfo.dat",
                  name[:-len("_band.dat")] + ".win"]
    elif lower.endswith("spaghetti_ene"):
        stem = name[:-len(".spaghetti_ene")] if lower.endswith(".spaghetti_ene") else name
        wanted = [stem + ".klist_band"]
    elif lower.startswith("eigenval"):
        wanted = ["POSCAR", "CONTCAR", "PROCAR", "OUTCAR", "DOSCAR", "KPOINTS"]
    elif lower.endswith(".dat"):
        wanted = ["KPOINTS", "orbital.dat", "kcal.dat"]
        wanted = [w for w in wanted if w.lower() != lower]
    return [os.path.join(folder, w) for w in wanted
            if os.path.isfile(os.path.join(folder, w))]


def read(path: str, **options):
    """Whatever calculation ``path`` holds, told by its content (so a file
    renamed on the way still reads). ``options`` the chosen reader does not
    take are dropped, so a caller can offer them all -- ``progress``, a
    ``report(fraction, message)`` callable, among them."""
    reader = {"bxsf": read_bxsf, "hr": read_hr,
              "wannier_band": read_wannier_band, "spaghetti": read_spaghetti,
              "eigenval": read_eigenval, "fplo_path": read_fplo_band,
              "fplo_grid": read_fplo_grid, "fplo_eigenvalue": read_fplo_eigenvalue,
              "matlab": read_matlab_grid}.get(detect(path))
    if reader is None and detect(path) in KINDS:
        raise ValueError(f"{os.path.basename(path)} is not a calculation but "
                         f"{KINDS[detect(path)][0]}: open the band file beside it")
    if reader is None:
        raise ValueError(
            f"not a band-structure file this program reads: "
            f"{os.path.basename(path)}. Read are .bxsf, *_hr.dat, *_band.dat, "
            ".spaghetti_ene, EIGENVAL, and the FPLO + Wannier tools' band.dat "
            "and Eigenvalue.dat; another layout can be written as one of "
            "them (see the user guide, 'A calculation in another format')")
    import inspect

    accepted = inspect.signature(reader).parameters
    return reader(path, **{k: v for k, v in options.items() if k in accepted})


# --------------------------------------------------------------------------
# XCrySDen .bxsf
# --------------------------------------------------------------------------
def read_bxsf(path: str, k_scale: float = 1.0, fermi=None, progress=None):
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
    for band in range(nbands):
        while position < len(rest) and not rest[position].upper().startswith("BAND"):
            position += 1
        position += 1
        start = position
        while position < len(rest) and not rest[position].upper().startswith(("BAND", "END")):
            position += 1
        numbers = np.array(" ".join(rest[start:position]).split(), dtype=float)
        values.append(numbers[:n[0] * n[1] * n[2]].reshape(n))
        if progress:
            progress((band + 1) / nbands, "Reading the band grid")
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


def read_hr(path: str, lattice=None, fermi=None, progress=None):
    """A Wannier90 ``seedname_hr.dat``. The cell comes from ``lattice``
    (rows, Å), else from ``seedname.win`` beside it; E_F likewise from the
    ``.win`` (``fermi_energy``) unless given. Read in pieces: a Hamiltonian
    of a hundred orbitals is millions of lines."""
    size = max(os.path.getsize(path), 1)
    with open(path, encoding="utf-8", errors="replace") as handle:
        handle.readline()
        nw = int(handle.readline().split()[0])
        nr = int(handle.readline().split()[0])
        degeneracy = []
        while len(degeneracy) < nr:
            degeneracy.extend(int(v) for v in handle.readline().split())
        H = np.zeros((nr, nw, nw), dtype=complex)
        order = {}                   # R -> its place in the file (= degeneracy's)
        done = 0
        for rows, read_bytes in _line_chunks(handle, 200000):
            done += read_bytes
            R = rows[:, :3].astype(int)
            unique, inverse = np.unique(R, axis=0, return_inverse=True)
            places = np.array([order.setdefault(tuple(r), len(order)) for r in unique])
            H[places[inverse.ravel()], rows[:, 3].astype(int) - 1,
              rows[:, 4].astype(int) - 1] = rows[:, 5] + 1j * rows[:, 6]
            if progress:
                progress(min(done / size, 1.0), "Reading the Hamiltonian")
    if progress:
        progress(1.0, "Reading the Hamiltonian")
    keys = np.zeros((len(order), 3))
    for r, place in order.items():
        keys[place] = r
    H = H[:len(order)]
    deg = np.array(degeneracy[:len(order)], dtype=float)
    match = _HR.search(path)
    seed = path[:match.start()] if match else os.path.splitext(path)[0]
    win_cell, win_fermi = read_win_cell(seed + ".win")
    folder = os.path.dirname(os.path.abspath(path))
    if win_cell is None:
        kpoints = read_fplo_kpoints(os.path.join(folder, "KPOINTS"))
        win_cell = kpoints["lattice"] if kpoints else None
    cell = lattice if lattice is not None else win_cell
    if cell is None:
        raise ValueError("the cell is needed: give it (or the Brillouin-zone "
                         "tool's crystal), or put seedname.win (or the FPLO "
                         "tools' KPOINTS) beside the file")
    names = _win_orbitals(seed + ".win", nw)
    groups = None
    if names is None:
        names, groups = _wannier_names(os.path.join(folder,
                                                    "wannier_name_center.txt"), nw)
    return WannierModel(cell, keys, H, degeneracy=deg,
                        fermi=(fermi if fermi is not None else (win_fermi or 0.0)),
                        orbital_names=names, groups=groups,
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


# --------------------------------------------------------------------------
# FPLO + Wannier tools (k_path.exe, H.exe, plot.exe)
# --------------------------------------------------------------------------
def _numbers(line: str):
    """The numbers of a line before any "!" comment, or None if a token
    is not one."""
    try:
        return [float(v.replace("D", "E").replace("d", "e"))
                for v in line.split("!")[0].split()]
    except ValueError:
        return None


def _fplo_layout(path: str):
    """Which FPLO-tools file ``path`` is, from its first lines: "path"
    (distance, energy[, weight]), "grid" (k1 k2 k3 E1 E2 ...),
    "eigenvalue" (a k line then one energy per line), or None."""
    rows = []
    try:
        with open(path, encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.strip():
                    rows.append(_numbers(line))
                    if len(rows) == 3:
                        break
    except OSError:
        return None
    if len(rows) < 2 or any(r is None or not r for r in rows):
        return None
    if len(rows[0]) == 3 and len(rows[1]) == 1:
        return "eigenvalue"
    if len(rows[0]) in (2, 3) and len(rows[1]) == len(rows[0]):
        return "path"
    if len(rows[0]) >= 4 and len(rows[1]) == len(rows[0]):
        return "grid"
    return None


def read_fplo_kpoints(path: str):
    """The FPLO tools' ``KPOINTS``: ``{"segments", "per_segment",
    "lattice" (rows, Å, or None), "points": [(fractional k, label), ...]}``
    -- each segment's two ends -- or None if ``path`` is not one (a VASP
    KPOINTS is not)."""
    try:
        with open(path, encoding="utf-8", errors="replace") as handle:
            lines = handle.read().splitlines()
    except OSError:
        return None
    if not lines:
        return None
    head = lines[0].split()
    try:
        segments, per_segment = int(head[0]), int(head[1])
    except (IndexError, ValueError):
        return None
    lattice, start = None, 1
    rows = [_numbers(l) for l in lines[1:4]]
    if len(rows) == 3 and all(r is not None and len(r) >= 6 for r in rows):
        lattice, start = np.array([r[:3] for r in rows]), 4
    points = []
    for line in lines[start:]:
        values = _numbers(line)
        if values is None or len(values) < 3:
            continue
        label = line.split("!", 1)[1].strip() if "!" in line else ""
        points.append((np.array(values[:3]), _greek(label)))
    if lattice is None and not points:
        return None
    return {"segments": segments, "per_segment": per_segment,
            "lattice": lattice, "points": points}


def _fplo_ticks(distance, kpoints):
    """The labelled points of a path of ``kpoints["segments"]`` segments of
    ``per_segment`` points each (both ends included, as k_path.exe writes
    them): ``[(distance, label), ...]``, a jump written "X|Y"."""
    if not kpoints:
        return []
    nseg, per = kpoints["segments"], kpoints["per_segment"]
    names = [label for _k, label in kpoints["points"]]
    if nseg * per != len(distance) or len(names) < 2 * nseg:
        return []
    ticks = [(float(distance[0]), names[0])]
    for s in range(nseg):
        label = names[2 * s + 1]
        if s + 1 < nseg and names[2 * s + 2] != label:
            label = f"{label}|{names[2 * s + 2]}"
        ticks.append((float(distance[(s + 1) * per - 1]), label))
    return ticks


def _fplo_blocks(path: str):
    """The rows of a path file, one list per band: bands are separated by
    blank lines, or -- without them -- start where the distance falls."""
    blocks, current = [], []
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            values = _numbers(line) if line.strip() else None
            if values:
                if current and values[0] < current[-1][0] - 1e-9:
                    blocks.append(current)
                    current = []
                current.append(values)
            elif current:
                blocks.append(current)
                current = []
    if current:
        blocks.append(current)
    if not blocks:
        raise ValueError("no bands in the file")
    size = len(blocks[0])
    return [np.array(b) for b in blocks if len(b) == size]


def read_fplo_band(path: str, lattice=None, fermi=None):
    """The FPLO tools' ``band.dat`` (or ``orbital.dat``) along a path:
    rows of distance (Å⁻¹, 2π included) and energy (eV, relative to E_F as
    FPLO gives it), a blank line between bands. The labels come from the
    ``KPOINTS`` beside it; a third column, or an ``orbital.dat`` of the
    same shape beside it, is a weight."""
    folder = os.path.dirname(os.path.abspath(path))
    blocks = _fplo_blocks(path)
    distance = blocks[0][:, 0]
    energies = np.array([b[:, 1] for b in blocks])
    weights = {}
    if blocks[0].shape[1] >= 3:
        weights[os.path.splitext(os.path.basename(path))[0]] = \
            np.array([b[:, 2] for b in blocks])
    other = os.path.join(folder, "orbital.dat")
    if not weights and os.path.isfile(other) and \
            os.path.abspath(other) != os.path.abspath(path):
        try:
            extra = _fplo_blocks(other)
            if len(extra) == len(blocks) and extra[0].shape == (len(distance), 3):
                weights["orbital"] = np.array([b[:, 2] for b in extra])
        except (OSError, ValueError):
            pass
    kpoints = read_fplo_kpoints(os.path.join(folder, "KPOINTS"))
    cell = lattice if lattice is not None else (kpoints or {}).get("lattice")
    k = None
    kcal = os.path.join(folder, "kcal.dat")
    if cell is not None and os.path.isfile(kcal):
        try:
            frac = np.loadtxt(kcal, ndmin=2)
            if frac.shape == (len(distance), 3):
                k = frac @ reciprocal_of(cell)
        except ValueError:
            pass
    return BandPath(distance, energies, _fplo_ticks(distance, kpoints), k=k,
                    fermi=0.0 if fermi is None else float(fermi),
                    weights=weights, name=os.path.basename(path))


def read_fplo_eigenvalue(path: str, lattice=None, fermi=None):
    """The FPLO tools' ``Eigenvalue.dat``: per k point a line of its
    coordinates -- 2π times fractional, the phase k.R that ``H.exe`` uses
    (plain fractional ones, all within ±1, are taken as they are) -- then
    one energy per line (eV, relative to E_F). The cell and labels come
    from the ``KPOINTS`` beside it."""
    folder = os.path.dirname(os.path.abspath(path))
    k_frac, energies, current = [], [], None
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            values = _numbers(line) if line.strip() else None
            if not values:
                continue
            if len(values) >= 3:
                k_frac.append(values[:3])
                current = []
                energies.append(current)
            elif current is not None:
                current.append(values[0])
    if not energies:
        raise ValueError("no k points in the file")
    count = len(energies[0])
    if any(len(e) != count for e in energies):
        raise ValueError("the k points do not all have the same number of bands")
    kpoints = read_fplo_kpoints(os.path.join(folder, "KPOINTS"))
    cell = lattice if lattice is not None else (kpoints or {}).get("lattice")
    if cell is None:
        raise ValueError("the cell is needed: put the FPLO tools' KPOINTS "
                         "beside Eigenvalue.dat (or use band.dat)")
    k_frac = np.array(k_frac)
    if np.abs(k_frac).max() > 1.0 + 1e-6:
        k_frac = k_frac / (2.0 * np.pi)
    k_cart = k_frac @ reciprocal_of(cell)
    steps = np.r_[0.0, np.linalg.norm(np.diff(k_cart, axis=0), axis=1)]
    if kpoints and kpoints["segments"] * kpoints["per_segment"] == len(k_cart):
        # a segment starts where the last ended, or after a jump: no distance
        steps[kpoints["per_segment"]::kpoints["per_segment"]] = 0.0
    distance = np.cumsum(steps)
    return BandPath(distance, np.array(energies).T, _fplo_ticks(distance, kpoints),
                    k=k_cart, fermi=0.0 if fermi is None else float(fermi),
                    name=os.path.basename(path))


def _line_chunks(handle, rows: int, usecols=None):
    """``(numbers, bytes read)`` of an open text table, ``rows`` lines at
    a time (blank lines skipped)."""
    from itertools import islice

    while True:
        lines = list(islice(handle, rows))
        if not lines:
            return
        size = sum(len(l) for l in lines)
        lines = [l for l in lines if l.strip()]
        if lines:
            yield np.loadtxt(lines, ndmin=2, usecols=usecols), size


def _row_chunks(path: str, rows: int, usecols=None):
    """``(numbers, bytes read)`` of a big text table, ``rows`` lines at a
    time."""
    with open(path, encoding="utf-8", errors="replace") as handle:
        yield from _line_chunks(handle, rows, usecols)


def read_fplo_grid(path: str, lattice=None, fermi=None, window=(-2.0, 1.0),
                   chunk: int = 20000, progress=None):
    """The FPLO tools' 3-D ``band.dat``: per row the fractional k
    (``k1 k2 k3``, each from -1/2 to 1/2, ends included) and every band's
    energy (eV, relative to E_F). Such a file is large (gigabytes for a
    101³ grid of a hundred bands), so it is read in pieces, twice, and
    only the bands that come within ``window`` (eV about E_F) are kept --
    the others cannot reach a measured cut. The cell comes from
    ``lattice`` or the FPLO tools' ``KPOINTS`` beside the file."""
    folder = os.path.dirname(os.path.abspath(path))
    if lattice is None:
        lattice = (read_fplo_kpoints(os.path.join(folder, "KPOINTS")) or {}).get("lattice")
    if lattice is None:
        raise ValueError("the cell is needed: put the KPOINTS of the band "
                         "calculation along the path beside band.dat, or "
                         "set the crystal in the Brillouin-zone tool")
    fermi = 0.0 if fermi is None else float(fermi)
    size = max(os.path.getsize(path), 1)
    lo = hi = None
    axes = [np.empty(0)] * 3
    total = done = 0
    for block, read_bytes in _row_chunks(path, chunk):
        done += read_bytes
        if progress:
            progress(0.5 * done / size, "Reading the band grid (1 of 2): which "
                                        "bands come near E_F")
        e = block[:, 3:]
        lo = e.min(0) if lo is None else np.minimum(lo, e.min(0))
        hi = e.max(0) if hi is None else np.maximum(hi, e.max(0))
        axes = [np.union1d(axes[i], np.round(block[:, i], 6)) for i in range(3)]
        total += len(block)
    if lo is None:
        raise ValueError("no k points in the file")
    if max(np.abs(a).max() for a in axes) > 1.0 + 1e-6:
        raise ValueError("the first three columns are not fractional k "
                         "(in units of the reciprocal vectors)")
    n = [len(a) for a in axes]
    if int(np.prod(n)) != total:
        raise ValueError(f"the {total} k points do not fill a regular "
                         f"{n[0]}x{n[1]}x{n[2]} grid")
    keep = np.flatnonzero((hi - fermi >= window[0]) & (lo - fermi <= window[1]))
    if keep.size == 0:
        raise ValueError(f"no band comes within {window[0]:g} to {window[1]:g} "
                         "eV of E_F")
    values = np.empty((keep.size, *n), dtype=np.float32)
    done = 0
    for block, read_bytes in _row_chunks(path, chunk, usecols=[0, 1, 2, *(3 + keep)]):
        done += read_bytes
        if progress:
            progress(0.5 + 0.5 * done / size, f"Reading the band grid (2 of 2): "
                                              f"{keep.size} bands near E_F")
        index = tuple(np.searchsorted(axes[i], np.round(block[:, i], 6))
                      for i in range(3))
        values[(slice(None), *index)] = block[:, 3:].T
    general = all(abs(a[-1] - a[0] - 1.0) < 1e-4 for a in axes)
    return _grid_model(lattice, values, [float(a[0]) for a in axes], general,
                       fermi, keep, int(lo.size), os.path.basename(path))


def _grid_model(lattice, values, origin, general, fermi, keep, total, name):
    reciprocal = reciprocal_of(lattice)
    options = dict(origin=origin, fermi=fermi, name=name)
    grid = (BandGrid.from_general_grid(reciprocal, values, **options) if general
            else BandGrid(reciprocal, values, **options))
    #: the kept bands' numbers in the file (1-based), as the scripts count them
    grid.band_numbers = np.asarray(keep) + 1
    grid.total_bands = total
    return grid


def read_matlab_grid(path: str, lattice=None, fermi=None, window=(-2.0, 1.0),
                     progress=None):
    """A MATLAB ``.mat`` holding the bands on a 3-D grid -- the ``Ek_3D``
    array (n1 x n2 x n3 x bands) the FPLO-tools scripts save from a 3-D
    ``band.dat``, on their grid: fractional k from -1/2 to 1/2 along each
    axis, ends included. Version 7.3 files (HDF5) and older ones both
    read. The cell comes from ``lattice`` or the FPLO tools' ``KPOINTS``
    beside the file."""
    folder = os.path.dirname(os.path.abspath(path))
    if lattice is None:
        lattice = (read_fplo_kpoints(os.path.join(folder, "KPOINTS")) or {}).get("lattice")
    if lattice is None:
        raise ValueError("the cell is needed: put the KPOINTS of the band "
                         "calculation along the path beside the .mat file, or "
                         "set the crystal in the Brillouin-zone tool")
    if progress:
        progress(-1, "Reading the MATLAB file")
    with open(path, "rb") as handle:
        head = handle.read(128)
    if b"7.3" in head[:116] or head[:4] == b"\x89HDF":
        import h5py

        with h5py.File(path, "r") as handle:
            arrays = {k: v for k, v in handle.items()
                      if hasattr(v, "shape") and len(v.shape) == 4}
            if not arrays:
                raise ValueError("no 4-D array (k1 x k2 x k3 x bands) in the file")
            key = next((k for k in arrays if "ek" in k.lower()), next(iter(arrays)))
            # HDF5 stores MATLAB's column-major array with its axes reversed
            values = np.asarray(arrays[key][()], dtype=np.float32).transpose(0, 3, 2, 1)
    else:
        from scipy.io import loadmat

        contents = {k: v for k, v in loadmat(path).items()
                    if isinstance(v, np.ndarray) and v.ndim == 4}
        if not contents:
            raise ValueError("no 4-D array (k1 x k2 x k3 x bands) in the file")
        key = next((k for k in contents if "ek" in k.lower()), next(iter(contents)))
        values = np.asarray(contents[key], dtype=np.float32).transpose(3, 0, 1, 2)
    fermi = 0.0 if fermi is None else float(fermi)
    lo = values.reshape(len(values), -1).min(1)
    hi = values.reshape(len(values), -1).max(1)
    keep = np.flatnonzero((hi - fermi >= window[0]) & (lo - fermi <= window[1]))
    if keep.size == 0:
        raise ValueError(f"no band comes within {window[0]:g} to {window[1]:g} "
                         "eV of E_F")
    return _grid_model(lattice, values[keep], [-0.5, -0.5, -0.5], True, fermi,
                       keep, len(values), os.path.basename(path))


def _wannier_names(path: str, count: int):
    """The Wannier functions' names in ``wannier_name_center.txt`` (rows
    of index, atom, orbital, centre) and their groups -- each element's
    shell ("Cr 3d"), and each spin -- or ``(None, None)``."""
    names, groups = [], {}
    try:
        with open(path, encoding="utf-8", errors="replace") as handle:
            for line in handle:
                parts = line.split()
                if len(parts) >= 3:
                    names.append(f"{parts[1]} {parts[2]}")
    except OSError:
        return None, None
    if len(names) != count:
        return None, None
    spins = {"UP": "spin up", "DN": "spin down", "DOWN": "spin down"}
    for i, name in enumerate(names):
        atom, orbital = name.split(" ", 1)
        element = atom.rstrip("0123456789")
        shell = re.match(r"\d*[spdf]", orbital)
        if shell:
            groups.setdefault(f"{element} {shell.group(0)}", []).append(i)
        spin = re.search(r"(UP|DN|DOWN)$", orbital.upper())
        if spin:
            groups.setdefault(spins[spin.group(1)], []).append(i)
    return names, groups


# --------------------------------------------------------------------------
# Which file is which: by content, for a folder the theory side sent
# --------------------------------------------------------------------------
def _head_lines(path: str, size: int = 65536):
    """The first lines of a text file (the last, maybe cut, dropped when
    the file is longer), or None for a binary file."""
    try:
        with open(path, "rb") as handle:
            head = handle.read(size)
            more = bool(handle.read(1))
    except OSError:
        return None
    if b"\x00" in head[:4096]:
        return None
    lines = head.decode("utf-8", errors="replace").splitlines()
    return lines[:-1] if more and len(lines) > 1 else lines


def _ints(line: str):
    try:
        return [int(v) for v in line.split()]
    except ValueError:
        return None


def _is_path_table(lines) -> bool:
    """Rows of distance, energy[, weight]: the distance never falls within
    a band, and bands are told apart by blank lines (or by the distance
    starting again)."""
    blocks, current = [], []
    for line in lines:
        values = _numbers(line) if line.strip() else None
        if values:
            if len(values) not in (2, 3):
                return False
            current.append(values[0])
        elif current:
            blocks.append(current)
            current = []
    if current:
        blocks.append(current)
    if not blocks or len(blocks[0]) < 2:
        return False
    first = np.asarray(blocks[0])
    steps = np.diff(first)
    if np.all(steps >= -1e-9):
        return True
    # no blank lines: the distance must fall back to where it started
    drops = np.flatnonzero(steps < -1e-9)
    return bool(np.allclose(first[drops + 1], first[0], atol=1e-6))


def detect(path: str):
    """What ``path`` holds, from its content (and only then its name): a
    key of :data:`KINDS`, or None."""
    name = os.path.basename(path).lower()
    try:
        with open(path, "rb") as handle:
            magic = handle.read(128)
    except OSError:
        return None
    if magic.startswith(b"MATLAB") or (name.endswith(".mat") and magic[:4] == b"\x89HDF"):
        return "matlab"
    lines = _head_lines(path)
    if lines is None:
        return None
    text = "\n".join(lines)
    upper = text.upper()
    nonblank = [l for l in lines if l.strip()]
    if not nonblank:
        return None
    if "BANDGRID_3D" in upper:
        return "bxsf"
    if "BANDINDEX" in upper:
        return "spaghetti"
    if "# OF K-POINTS" in upper or nonblank[0].strip().upper().startswith("PROCAR"):
        return "procar"
    if re.search(r"begin\s+unit_cell_cart|^\s*num_wann", text, re.I | re.M):
        return "win"
    if len(nonblank) > 4 and len(_ints(nonblank[1]) or []) == 1 \
            and len(_ints(nonblank[2]) or []) == 1 and _ints(nonblank[3]):
        rows = [r.split() for r in nonblank[3:40]]
        if any(len(r) == 7 and _ints(" ".join(r[:5])) for r in rows):
            return "hr"
    if len(lines) > 6 and len(_ints(lines[0]) or []) == 4 \
            and len(_ints(lines[5]) or []) == 3:
        return "eigenval"
    if name.endswith("klist_band"):
        return "klist"
    if read_fplo_kpoints(path) is not None and (read_fplo_kpoints(path) or {}).get("lattice") is not None:
        return "fplo_kpoints"
    if len(nonblank) >= 3 and len(_ints(nonblank[1]) or []) == 1 \
            and nonblank[2].strip()[:1].lower() in "lgmac":
        return "kpoints"
    if len(nonblank) >= 7 and _numbers(nonblank[1]) and len(_numbers(nonblank[1])) == 1 \
            and all(len(_numbers(l) or []) == 3 for l in nonblank[2:5]):
        return "poscar"
    parts = [l.split() for l in nonblank[:20]]
    if all(len(p) >= 6 and _ints(p[0]) and not _numbers(p[1]) for p in parts):
        return "wannier_names"
    if all(len(p) == 6 and not _numbers(p[0]) and _numbers(" ".join(p[1:]))
           for p in parts):
        return "labelinfo"
    layout = _fplo_layout(path)
    if layout == "eigenvalue":
        return "fplo_eigenvalue"
    if layout == "grid":
        return "fplo_grid"
    if layout == "path" or (layout is None and _numbers(nonblank[0]) and
                            len(_numbers(nonblank[0])) in (2, 3)):
        if _is_path_table(lines):
            return "wannier_band" if name.endswith("_band.dat") else "fplo_path"
        if _numbers(nonblank[0]) and len(_numbers(nonblank[0])) == 3:
            return "kfile"
    if all(len(_numbers(l) or []) == 1 for l in nonblank[:50]):
        return "kfile"
    return None



def _needs(path: str, kind: str):
    """What reading ``path`` (of ``kind``) still lacks, in plain words, and
    what is only nice to have: ``(missing, optional)`` lists."""
    folder, name = os.path.split(os.path.abspath(path))
    present = {f.lower() for f in os.listdir(folder)}
    fplo_cell = (read_fplo_kpoints(os.path.join(folder, "KPOINTS")) or {}).get("lattice") is not None
    missing, optional = [], []
    if kind == "hr":
        match = _HR.search(name)
        seed = name[:match.start()] if match else os.path.splitext(name)[0]
        if (seed + ".win").lower() not in present and not fplo_cell:
            missing.append(f"the cell: {seed}.win (or the FPLO tools' KPOINTS) "
                           "beside it -- or set the crystal in the "
                           "Brillouin-zone tool")
        if "wannier_name_center.txt" not in present and (seed + ".win").lower() not in present:
            optional.append("the orbitals' names (wannier_name_center.txt or the "
                            ".win) for orbital weights")
    elif kind in ("fplo_grid", "matlab"):
        if not fplo_cell:
            missing.append("the cell: the FPLO tools' KPOINTS beside it -- or "
                           "set the crystal in the Brillouin-zone tool")
    elif kind == "fplo_eigenvalue":
        if not fplo_cell:
            missing.append("the FPLO tools' KPOINTS beside it (the cell)")
    elif kind == "fplo_path":
        if not fplo_cell:
            optional.append("the FPLO tools' KPOINTS, for the points' names")
        if "orbital.dat" not in present:
            optional.append("orbital.dat, for an orbital weight")
    elif kind == "eigenval":
        if "poscar" not in present:
            missing.append("POSCAR beside it (the cell)")
        if "procar" not in present:
            optional.append("PROCAR, for orbital weights")
    elif kind == "wannier_band":
        if (name[:-len(".dat")] + ".labelinfo.dat").lower() not in present:
            optional.append("its _band.labelinfo.dat, for the points' names")
    elif kind == "spaghetti":
        stem = name.lower().split(".spaghetti_ene")[0]
        if stem + ".klist_band" not in present:
            optional.append("case.klist_band, for the points' names")
    return missing, optional


def split_groups(names):
    """Files that look like one file cut into pieces (to get it past an
    upload limit): ``band.dat.1, band.dat.2``, ``band_part1.dat, ...``,
    ``xaa, xab`` (what ``split`` writes). ``[(joined name, [pieces in
    order]), ...]``."""
    numbered = re.compile(r"^(?P<stem>.+?)(?P<sep>[._-]?(?:part|p)?)(?P<num>\d+)"
                          r"(?P<ext>\.[A-Za-z]\w*)?$", re.I)
    lettered = re.compile(r"^(?P<stem>.*?)(?P<sep>\.?)(?P<tag>[a-z]{2})$")
    groups = {}
    for name in names:
        m = numbered.match(name)
        if m:
            key = ("n", m["stem"], m["sep"], m["ext"] or "")
            groups.setdefault(key, []).append((int(m["num"]), name))
            continue
        m = lettered.match(name)
        if m and (m["stem"] == "x" or m["sep"]):
            key = ("a", m["stem"], m["sep"], "")
            tag = m["tag"]
            groups.setdefault(key, []).append(((ord(tag[0]) - 97) * 26 + ord(tag[1]) - 97,
                                               name))
    out = []
    for (style, stem, sep, ext), members in groups.items():
        members.sort()
        numbers = [n for n, _ in members]
        if len(members) < 2 or numbers != list(range(numbers[0], numbers[0] + len(members))):
            continue
        if style == "a" and numbers[0] != 0:
            continue
        joined = (stem.rstrip("._-") or "joined") + ext
        out.append((joined, [n for _, n in members]))
    return out


def join_pieces(paths, target: str) -> str:
    """The pieces written one after the other into ``target`` (a line
    break added where a piece does not end in one)."""
    with open(target, "wb") as out:
        for path in paths:
            with open(path, "rb") as handle:
                last = b""
                while True:
                    block = handle.read(1 << 22)
                    if not block:
                        break
                    out.write(block)
                    last = block
            if last and not last.endswith(b"\n"):
                out.write(b"\n")
    return target


def survey(folder: str):
    """Every file in ``folder``, said in plain words: a list of dicts with
    ``path``, ``name``, ``kind`` (a :data:`KINDS` key or None), ``what``,
    ``goals`` (what it can show), ``missing`` and ``optional`` (what reading
    it lacks / would gain), ``notes`` (renamed, large, in pieces) and
    ``pieces`` (for a file in pieces: their paths, else None)."""
    names = sorted(n for n in os.listdir(folder)
                   if os.path.isfile(os.path.join(folder, n)) and not n.startswith("."))
    rows = []
    in_pieces = {}
    for joined, pieces in split_groups(names):
        first = detect(os.path.join(folder, pieces[0]))
        columns = []
        for piece in pieces:
            lines = [l for l in (_head_lines(os.path.join(folder, piece), 8192) or [])
                     if l.strip()]
            values = _numbers(lines[0]) if lines else None
            columns.append(len(values) if values else None)
        if first is None and (None in columns or len(set(columns)) != 1):
            continue
        for piece in pieces:
            in_pieces[piece] = (joined, pieces)
    for name in names:
        path = os.path.join(folder, name)
        kind = detect(path)
        what, goals = KINDS.get(kind, ("not a file this program reads", frozenset()))
        missing, optional = _needs(path, kind) if kind else ([], [])
        notes = []
        usual = _USUAL_NAME.get(kind)
        if usual is not None and not usual(name.lower()):
            notes.append("an unusual name for this kind of file -- it opens all "
                         "the same")
        size = os.path.getsize(path)
        if size > 200e6:
            notes.append(f"{size / 1e9:.1f} GB: reading takes about "
                         f"{max(1, round(size / 40e6 / 60))} min")
        pieces = None
        if name in in_pieces:
            joined, members = in_pieces[name]
            pieces = [os.path.join(folder, m) for m in members]
            notes.append(f"looks like piece {members.index(name) + 1} of "
                         f"{len(members)} of one file ({joined}): join them first")
            if kind is None:
                what = "a piece of a cut-up file"
        rows.append({"path": path, "name": name, "kind": kind, "what": what,
                     "goals": set(goals), "missing": missing,
                     "optional": optional, "notes": notes, "pieces": pieces})
    return rows
