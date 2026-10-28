"""
tools/crystal/structure.py
==========================
A crystal in real space: the cell, the atoms in it, and planes through it.

The lab's ``Lattice_plot.m`` (zqs) drew the fourteen Bravais cells with
balls on the lattice points and a Miller plane cut through them -- one
picture per lattice type, written out by hand (``lattice_bcc.m``,
``lattice_fcc.m``, ``cell_hex_plot.m``... one file each), with no atoms
of the actual crystal and no way to read a structure from a file. This
module is the general version:

* **Any cell, any atoms.** The six cell parameters, and atoms at
  fractional positions -- typed in, put on the lattice points of the
  Brillouin-zone tool's last crystal (its centering included), or read
  from a **CIF file** with its symmetry operations applied (so the
  asymmetric unit a CIF lists becomes the full cell).
* **Planes as a calculation, not a drawing.** A plane (hkl) is described by
  its normal (the reciprocal vector G_hkl) and its distance from the origin
  along it; its outline is the plane clipped to the drawn block of cells.
  The atomic layers parallel to it -- the possible terminations of a
  cleave along (hkl) -- are listed with the atoms in each, so the plane can
  be put exactly on one.

Qt-free, so it can be tested on its own. The drawing is
:mod:`ui.panels.crystal`.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field

import numpy as np

__all__ = ["Atom", "Structure", "parse_cif", "read_cif", "parse_symop",
           "expand", "lattice_point_atoms", "cartesian_atoms", "plane_normal",
           "plane_outline", "atomic_layers", "element_of", "ELEMENTS",
           "covalent_radius", "element_colour"]

#: Element -> (covalent radius in Å, colour) -- radii after Cordero et al.,
#: Dalton Trans. 2008; colours the usual Jmol/CPK scheme.
ELEMENTS = {
    "H": (0.31, "#ffffff"), "He": (0.28, "#d9ffff"), "Li": (1.28, "#cc80ff"),
    "Be": (0.96, "#c2ff00"), "B": (0.84, "#ffb5b5"), "C": (0.76, "#909090"),
    "N": (0.71, "#3050f8"), "O": (0.66, "#ff0d0d"), "F": (0.57, "#90e050"),
    "Ne": (0.58, "#b3e3f5"), "Na": (1.66, "#ab5cf2"), "Mg": (1.41, "#8aff00"),
    "Al": (1.21, "#bfa6a6"), "Si": (1.11, "#f0c8a0"), "P": (1.07, "#ff8000"),
    "S": (1.05, "#ffff30"), "Cl": (1.02, "#1ff01f"), "Ar": (1.06, "#80d1e3"),
    "K": (2.03, "#8f40d4"), "Ca": (1.76, "#3dff00"), "Sc": (1.70, "#e6e6e6"),
    "Ti": (1.60, "#bfc2c7"), "V": (1.53, "#a6a6ab"), "Cr": (1.39, "#8a99c7"),
    "Mn": (1.39, "#9c7ac7"), "Fe": (1.32, "#e06633"), "Co": (1.26, "#f090a0"),
    "Ni": (1.24, "#50d050"), "Cu": (1.32, "#c88033"), "Zn": (1.22, "#7d80b0"),
    "Ga": (1.22, "#c28f8f"), "Ge": (1.20, "#668f8f"), "As": (1.19, "#bd80e3"),
    "Se": (1.20, "#ffa100"), "Br": (1.20, "#a62929"), "Kr": (1.16, "#5cb8d1"),
    "Rb": (2.20, "#702eb0"), "Sr": (1.95, "#00ff00"), "Y": (1.90, "#94ffff"),
    "Zr": (1.75, "#94e0e0"), "Nb": (1.64, "#73c2c9"), "Mo": (1.54, "#54b5b5"),
    "Tc": (1.47, "#3b9e9e"), "Ru": (1.46, "#248f8f"), "Rh": (1.42, "#0a7d8c"),
    "Pd": (1.39, "#006985"), "Ag": (1.45, "#c0c0c0"), "Cd": (1.44, "#ffd98f"),
    "In": (1.42, "#a67573"), "Sn": (1.39, "#668080"), "Sb": (1.39, "#9e63b5"),
    "Te": (1.38, "#d47a00"), "I": (1.39, "#940094"), "Xe": (1.40, "#429eb0"),
    "Cs": (2.44, "#57178f"), "Ba": (2.15, "#00c900"), "La": (2.07, "#70d4ff"),
    "Ce": (2.04, "#ffffc7"), "Pr": (2.03, "#d9ffc7"), "Nd": (2.01, "#c7ffc7"),
    "Sm": (1.98, "#8fffc7"), "Eu": (1.98, "#61ffc7"), "Gd": (1.96, "#45ffc7"),
    "Tb": (1.94, "#30ffc7"), "Dy": (1.92, "#1fffc7"), "Ho": (1.92, "#00ff9c"),
    "Er": (1.89, "#00e675"), "Tm": (1.90, "#00d452"), "Yb": (1.87, "#00bf38"),
    "Lu": (1.87, "#00ab24"), "Hf": (1.75, "#4dc2ff"), "Ta": (1.70, "#4da6ff"),
    "W": (1.62, "#2194d6"), "Re": (1.51, "#267dab"), "Os": (1.44, "#266696"),
    "Ir": (1.41, "#175487"), "Pt": (1.36, "#d0d0e0"), "Au": (1.36, "#ffd123"),
    "Hg": (1.32, "#b8b8d0"), "Tl": (1.45, "#a6544d"), "Pb": (1.46, "#575961"),
    "Bi": (1.48, "#9e4fb5"), "Th": (2.06, "#00baff"), "U": (1.96, "#008fff"),
}

#: what an atom of an unknown element is drawn as
UNKNOWN = (1.0, "#ff1493")


def element_of(text: str) -> str:
    """The element symbol in a CIF label or type: "Bi1" -> "Bi",
    "O2-" -> "O", "SE" -> "Se"."""
    match = re.match(r"([A-Za-z]{1,2})", str(text or "").strip())
    if not match:
        return "X"
    letters = match.group(1)
    two = letters[:1].upper() + letters[1:2].lower()
    if len(letters) == 2 and two in ELEMENTS:
        return two
    return letters[:1].upper()


def covalent_radius(element: str) -> float:
    return ELEMENTS.get(element, UNKNOWN)[0]


def element_colour(element: str) -> str:
    return ELEMENTS.get(element, UNKNOWN)[1]


@dataclass
class Atom:
    label: str
    element: str
    frac: tuple                    # fractional (x, y, z) in the cell
    occupancy: float = 1.0


@dataclass
class Structure:
    """The cell (a, b, c in Å; alpha, beta, gamma in degrees) and every
    atom in it, at fractional positions in [0, 1)."""
    a: float
    b: float
    c: float
    alpha: float = 90.0
    beta: float = 90.0
    gamma: float = 90.0
    atoms: list = field(default_factory=list)
    name: str = ""
    space_group: str = ""

    def vectors(self) -> np.ndarray:
        """The cell vectors as rows (Å), a along x and b in the xy plane --
        the Brillouin-zone tool's convention."""
        from tools.crystal.lattice import LatticeParams, conventional_vectors
        return conventional_vectors(LatticeParams(self.a, self.b, self.c,
                                                  self.alpha, self.beta,
                                                  self.gamma))

    def elements(self):
        seen = []
        for atom in self.atoms:
            if atom.element not in seen:
                seen.append(atom.element)
        return seen

    def formula(self) -> str:
        counts = {}
        for atom in self.atoms:
            counts[atom.element] = counts.get(atom.element, 0) + atom.occupancy
        if not counts:
            return ""
        values = np.array(list(counts.values()))
        divisor = 1
        for d in range(int(min(values)), 0, -1):
            if np.allclose(values / d, np.round(values / d), atol=0.02):
                divisor = d
                break
        parts = []
        for element, count in counts.items():
            n = count / divisor
            parts.append(element + ("" if abs(n - 1) < 0.02 else
                                    (f"{n:.0f}" if abs(n - round(n)) < 0.02
                                     else f"{n:.2g}")))
        return "".join(parts)


# --------------------------------------------------------------------------
# Symmetry operations
# --------------------------------------------------------------------------
def parse_symop(text: str):
    """``(R, t)`` from a symmetry operation written as in a CIF:
    ``"-y, x-y, z+1/3"`` -> a 3x3 integer-ish matrix and a translation."""
    parts = [p.strip().lower() for p in str(text).strip().strip("'\"").split(",")]
    if len(parts) != 3:
        raise ValueError(f"not a symmetry operation: {text!r}")
    rotation = np.zeros((3, 3))
    translation = np.zeros(3)
    for row, part in enumerate(parts):
        expression = part.replace(" ", "")
        for sign, number, axis in re.findall(
                r"([+-]?)(\d*\.?\d*(?:/\d+)?)\*?([xyz]?)", expression):
            if not number and not axis:
                continue
            value = 1.0
            if number:
                if "/" in number:
                    top, bottom = number.split("/")
                    value = float(top or 1) / float(bottom)
                else:
                    value = float(number)
            if sign == "-":
                value = -value
            if axis:
                rotation[row, "xyz".index(axis)] += value
            else:
                translation[row] += value
    return rotation, translation


def _unique_positions(positions, tolerance: float = 1e-3):
    """Fractional positions wrapped into [0, 1) with duplicates (also
    across the cell boundary) removed."""
    kept = []
    for p in np.mod(np.asarray(positions, dtype=float), 1.0):
        p = np.where(np.abs(p - 1.0) < tolerance, 0.0, p)
        if not any(np.all(np.abs(((p - q) + 0.5) % 1.0 - 0.5) < tolerance)
                   for q in kept):
            kept.append(p)
    return kept


def expand(atoms, operations):
    """Every atom of ``atoms`` (the asymmetric unit) moved by every
    operation, wrapped into the cell and with duplicates removed."""
    out = []
    for atom in atoms:
        frac = np.asarray(atom.frac, dtype=float)
        images = [rotation @ frac + translation for rotation, translation in operations]
        for i, position in enumerate(_unique_positions(images)):
            out.append(Atom(atom.label if i == 0 else f"{atom.label}",
                            atom.element, tuple(float(v) for v in position),
                            atom.occupancy))
    # A site generated from two different labels is the same atom.
    unique, seen = [], []
    for atom in out:
        p = np.asarray(atom.frac)
        if any(np.all(np.abs(((p - q) + 0.5) % 1.0 - 0.5) < 1e-3) and e == atom.element
               for q, e in seen):
            continue
        seen.append((p, atom.element))
        unique.append(atom)
    return unique


#: The lattice-point translations of each centering (conventional cell).
CENTERING_TRANSLATIONS = {
    "P": [(0, 0, 0)],
    "A": [(0, 0, 0), (0, 0.5, 0.5)],
    "B": [(0, 0, 0), (0.5, 0, 0.5)],
    "C": [(0, 0, 0), (0.5, 0.5, 0)],
    "I": [(0, 0, 0), (0.5, 0.5, 0.5)],
    "F": [(0, 0, 0), (0, 0.5, 0.5), (0.5, 0, 0.5), (0.5, 0.5, 0)],
    "R": [(0, 0, 0), (2 / 3, 1 / 3, 1 / 3), (1 / 3, 2 / 3, 2 / 3)],
}


def lattice_point_atoms(centering: str = "P", element: str = "X"):
    """One atom on every lattice point of the conventional cell: the
    Bravais lattice itself, for a crystal known only by its cell."""
    return [Atom(f"{element}{i + 1}", element, tuple(float(v) for v in t))
            for i, t in enumerate(CENTERING_TRANSLATIONS.get(centering.upper(),
                                                              [(0, 0, 0)]))]


# --------------------------------------------------------------------------
# CIF
# --------------------------------------------------------------------------
def _cif_number(text) -> float:
    """``"5.4307(2)"`` -> 5.4307; ``"?"`` / ``"."`` -> NaN."""
    text = str(text).strip()
    if text in ("?", "."):
        return float("nan")
    return float(re.sub(r"\(\d+\)$", "", text))


def _tokens(line: str):
    return [t[1:-1] if len(t) > 1 and t[0] == t[-1] and t[0] in "'\"" else t
            for t in re.findall(r"'[^']*'|\"[^\"]*\"|\S+", line)]


def parse_cif(text: str) -> Structure:
    """The first data block of a CIF: cell, symmetry operations and atom
    sites, the sites expanded by the operations into the full cell."""
    values, loops = {}, []
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if not line or line.startswith("#"):
            i += 1
            continue
        if line.lower().startswith("data_") and values:
            break                                  # the next block
        if line.lower() == "loop_":
            names, rows = [], []
            i += 1
            while i < len(lines) and lines[i].strip().startswith("_"):
                names.append(lines[i].strip().split()[0].lower())
                i += 1
            tokens = []
            while i < len(lines):
                stripped = lines[i].strip()
                if (not stripped or stripped.startswith("_")
                        or stripped.lower() == "loop_"
                        or stripped.lower().startswith("data_")):
                    break
                if stripped.startswith(";"):
                    i += 1
                    while i < len(lines) and not lines[i].startswith(";"):
                        i += 1
                    i += 1
                    tokens.append("")
                    continue
                if not stripped.startswith("#"):
                    tokens.extend(_tokens(stripped))
                i += 1
            for start in range(0, len(tokens) - len(names) + 1, max(1, len(names))):
                rows.append(dict(zip(names, tokens[start:start + len(names)])))
            loops.append((names, rows))
            continue
        if line.startswith("_"):
            parts = _tokens(line)
            key = parts[0].lower()
            if len(parts) > 1:
                values[key] = " ".join(parts[1:])
            elif i + 1 < len(lines) and lines[i + 1].startswith(";"):
                i += 2
                block = []
                while i < len(lines) and not lines[i].startswith(";"):
                    block.append(lines[i])
                    i += 1
                values[key] = "\n".join(block)
            elif i + 1 < len(lines):
                i += 1
                values[key] = " ".join(_tokens(lines[i].strip()))
        i += 1

    def number(key, default=None):
        if key not in values:
            if default is None:
                raise ValueError(f"the CIF has no {key}")
            return default
        return _cif_number(values[key])

    structure = Structure(number("_cell_length_a"), number("_cell_length_b"),
                          number("_cell_length_c"),
                          number("_cell_angle_alpha", 90.0),
                          number("_cell_angle_beta", 90.0),
                          number("_cell_angle_gamma", 90.0))
    structure.name = (values.get("_chemical_formula_sum")
                      or values.get("_chemical_name_systematic") or "").strip()
    structure.space_group = (values.get("_symmetry_space_group_name_h-m")
                             or values.get("_space_group_name_h-m_alt") or "").strip()

    operations = []
    for names, rows in loops:
        for key in ("_symmetry_equiv_pos_as_xyz", "_space_group_symop_operation_xyz"):
            if key in names:
                operations = [parse_symop(row[key]) for row in rows if row.get(key)]
    if not operations:
        operations = [(np.eye(3), np.zeros(3))]

    sites = []
    for names, rows in loops:
        if "_atom_site_fract_x" not in names:
            continue
        for row in rows:
            label = row.get("_atom_site_label", "X")
            element = element_of(row.get("_atom_site_type_symbol") or label)
            try:
                frac = tuple(_cif_number(row[f"_atom_site_fract_{c}"]) for c in "xyz")
            except (KeyError, ValueError):
                continue
            if any(np.isnan(frac)):
                continue
            occupancy = row.get("_atom_site_occupancy", "1")
            try:
                occupancy = _cif_number(occupancy)
            except ValueError:
                occupancy = 1.0
            sites.append(Atom(label, element, frac,
                              1.0 if np.isnan(occupancy) else occupancy))
    if not sites:
        raise ValueError("the CIF lists no atom sites with fractional "
                         "coordinates")
    structure.atoms = expand(sites, operations)
    if not structure.name:
        structure.name = structure.formula()
    return structure


def read_cif(path: str) -> Structure:
    with open(path, encoding="utf-8", errors="replace") as handle:
        return parse_cif(handle.read())


# --------------------------------------------------------------------------
# Geometry
# --------------------------------------------------------------------------
def cartesian_atoms(structure: Structure, repeats=(1, 1, 1),
                    include_boundary: bool = True):
    """Every atom of a block of ``repeats`` cells, as ``(atoms, xyz)`` --
    the atom record and its Cartesian position (Å). Atoms on the block's
    faces, edges and corners are repeated on the far side too (when
    ``include_boundary``), so the block is drawn whole."""
    vectors = structure.vectors()
    n = np.asarray(repeats, dtype=int)
    out_atoms, out_xyz = [], []
    eps = 1e-6
    for atom in structure.atoms:
        frac = np.asarray(atom.frac, dtype=float)
        ranges = []
        for d in range(3):
            shifts = list(range(n[d]))
            if include_boundary and abs(frac[d]) < eps:
                shifts.append(n[d])
            ranges.append(shifts)
        for i in ranges[0]:
            for j in ranges[1]:
                for k in ranges[2]:
                    out_atoms.append(atom)
                    out_xyz.append((frac + (i, j, k)) @ vectors)
    return out_atoms, np.asarray(out_xyz, dtype=float).reshape(-1, 3)


def plane_normal(structure: Structure, hkl) -> tuple:
    """``(unit normal, d_hkl)`` of the plane (hkl): the normal is along
    G_hkl = h b1 + k b2 + l b3, and ``d = 2 pi / |G|``."""
    from tools.crystal.lattice import reciprocal_vectors

    hkl = np.asarray(hkl, dtype=float)
    if not np.any(hkl):
        raise ValueError("(000) is not a plane")
    g = hkl @ reciprocal_vectors(structure.vectors())
    length = float(np.linalg.norm(g))
    return g / length, 2.0 * np.pi / length


def _box_edges(vectors, repeats):
    corners = []
    for i in (0, 1):
        for j in (0, 1):
            for k in (0, 1):
                corners.append(np.array([i * repeats[0], j * repeats[1],
                                         k * repeats[2]]) @ vectors)
    edges = []
    for a in range(8):
        for b in range(a + 1, 8):
            if bin(a ^ b).count("1") == 1:
                edges.append((corners[a], corners[b]))
    return corners, edges


def box_edges(structure: Structure, repeats=(1, 1, 1)):
    """The twelve edges of the drawn block, as pairs of points (Å)."""
    return _box_edges(structure.vectors(), repeats)[1]


def plane_outline(structure: Structure, hkl, offset: float, repeats=(1, 1, 1)):
    """The plane ``n . r = offset`` (Å, n the unit normal of (hkl))
    clipped to the block of cells: its corners in order, or an empty array
    when it misses the block."""
    normal, _d = plane_normal(structure, hkl)
    _corners, edges = _box_edges(structure.vectors(), repeats)
    points = []
    for p, q in edges:
        a, b = float(normal @ p) - offset, float(normal @ q) - offset
        if a == b:
            continue
        if a * b <= 0:
            t = a / (a - b)
            points.append(p + t * (q - p))
    if len(points) < 3:
        return np.zeros((0, 3))
    points = np.unique(np.round(np.asarray(points), 9), axis=0)
    if len(points) < 3:
        return np.zeros((0, 3))
    centre = points.mean(axis=0)
    u = points[0] - centre
    u /= np.linalg.norm(u) or 1.0
    v = np.cross(normal, u)
    angles = np.arctan2((points - centre) @ v, (points - centre) @ u)
    return points[np.argsort(angles)]


def offset_range(structure: Structure, hkl, repeats=(1, 1, 1)):
    """The smallest and largest ``offset`` at which the plane meets the
    block."""
    normal, _d = plane_normal(structure, hkl)
    corners, _edges = _box_edges(structure.vectors(), repeats)
    heights = [float(normal @ c) for c in corners]
    return min(heights), max(heights)


def atomic_layers(structure: Structure, hkl, repeats=(1, 1, 1),
                  tolerance: float = 0.05):
    """The atomic layers parallel to (hkl) in the block: ``[(height in Å,
    {element: count}), ...]`` from the lowest up -- the planes a cleave can
    leave on top."""
    normal, _d = plane_normal(structure, hkl)
    atoms, xyz = cartesian_atoms(structure, repeats)
    if not len(atoms):
        return []
    heights = xyz @ normal
    order = np.argsort(heights)
    layers = []
    for index in order:
        h = float(heights[index])
        if layers and h - layers[-1][0] <= tolerance:
            counts = layers[-1][1]
        else:
            layers.append([h, {}])
            counts = layers[-1][1]
        element = atoms[index].element
        counts[element] = counts.get(element, 0) + 1
    return [(float(h), counts) for h, counts in layers]
