"""One-off generator: converts spglib's built-in space-group database into a
self-contained tools/spacegroups.py, so the running app needs neither spglib
nor any symmetry-operation machinery at runtime -- the Brillouin-zone tool
only ever needs two pieces of static reference data per International Tables
(ITA) space-group number (1-230), which the user types in directly:

* which of the 14 Bravais lattices it belongs to (and so which cell
  parameters are free and which are fixed by symmetry), and
* its point group's rotation matrices, which is what the *irreducible*
  Brillouin zone is the fundamental domain of.

Neither depends on an atomic basis (the user gives lattice parameters, not
atom positions), so there is nothing for spglib's actual symmetry search to
do here -- it is used purely as an authoritative lookup table. Re-run this
only if spglib ships a correction; the app itself never calls it.

The rotation matrices are taken from each group's standard setting, in
*conventional-cell fractional* coordinates (which is why they are integers,
and why every entry is -1, 0 or 1); tools.lattice.point_group_operations turns
them into Cartesian rotations with the conventional cell's own vectors. Only
the rotation parts are kept: screw/glide translations and the centering
translations move a point within one lattice site's orbit and have no effect
on the point group, which is all reciprocal space sees. Identical matrix
sets are stored once and shared -- 230 space groups use only 37 distinct
sets, a few more than the 32 crystallographic point groups because some
point groups occur in two orientations relative to the conventional axes
(P321 vs P312 and friends), which are genuinely different matrix sets and
must not be merged.

Usage:  python devtools/gen_spacegroups.py   (needs `pip install spglib`)
"""
import sys


def crystal_system(number: int) -> str:
    """ITA number -> crystal system, from the standard number ranges (the
    same ranges spglib's own docs and every crystallography textbook use)."""
    if number <= 2:
        return "triclinic"
    if number <= 15:
        return "monoclinic"
    if number <= 74:
        return "orthorhombic"
    if number <= 142:
        return "tetragonal"
    if number <= 167:
        return "trigonal"
    if number <= 194:
        return "hexagonal"
    if number <= 230:
        return "cubic"
    raise ValueError(f"not a valid space group number: {number}")


# Single-letter Bravais/crystal-family prefix per crystal system, as used in
# the two-letter Bravais symbols (aP, mP, mC, oP, oC, oI, oF, tP, tI, hP, hR,
# cP, cI, cF).
_FAMILY_LETTER = {
    "triclinic": "a", "monoclinic": "m", "orthorhombic": "o",
    "tetragonal": "t", "trigonal": "h", "hexagonal": "h", "cubic": "c",
}

#: order of each of the 32 crystallographic point groups, used only to check
#: that the matrix set pulled out of spglib really is the whole point group
#: (a silently truncated set would give a silently wrong irreducible zone).
_POINT_GROUP_ORDER = {
    "1": 1, "-1": 2, "2": 2, "m": 2, "2/m": 4, "222": 4, "mm2": 4, "mmm": 8,
    "4": 4, "-4": 4, "4/m": 8, "422": 8, "4mm": 8, "-42m": 8, "4/mmm": 16,
    "3": 3, "-3": 6, "32": 6, "3m": 6, "-3m": 12, "6": 6, "-6": 6, "6/m": 12,
    "622": 12, "6mm": 12, "-6m2": 12, "6/mmm": 24, "23": 12, "m-3": 24,
    "432": 24, "-43m": 24, "m-3m": 48,
}


def _check_is_a_group(matrices, label):
    """Closure and inverses, on the integer matrices themselves.

    Cheap, and it is the property everything downstream relies on: the
    irreducible zone is a fundamental domain of this set acting on the
    Brillouin zone, which is only well defined if the set is a group.
    """
    import numpy as np

    mats = [np.array(m, dtype=int).reshape(3, 3) for m in matrices]
    known = {tuple(m.ravel().tolist()) for m in mats}
    for x in mats:
        if abs(round(float(np.linalg.det(x)))) != 1:
            raise RuntimeError(f"{label}: matrix with |det| != 1")
        for y in mats:
            if tuple((x @ y).ravel().tolist()) not in known:
                raise RuntimeError(f"{label}: not closed under multiplication")


def main():
    import spglib

    # spglib enumerates 530 Hall settings (alternate origins/axis choices)
    # covering the 230 ITA numbers; several Hall numbers can share one ITA
    # number (different settings of the *same* space group -- same Bravais
    # lattice, same crystal system). Keep the first (lowest Hall number,
    # spglib's own "standard" entry) seen for each ITA number.
    by_number = {}
    for hall in range(1, 531):
        t = spglib.get_spacegroup_type(hall)
        if t.number not in by_number:
            by_number[t.number] = (hall, t)

    if len(by_number) != 230:
        raise RuntimeError(f"expected 230 space groups, spglib gave {len(by_number)}")

    # -- point-group rotations, deduplicated across space groups -----------
    operation_sets = []          # list of unique sets, in first-seen order
    operation_index = {}         # the set itself -> its index in the list
    for number in range(1, 231):
        hall, t = by_number[number]
        rotations = spglib.get_symmetry_from_database(hall)["rotations"]
        # Unique rotation parts = the point group. Sorted so that two groups
        # with the same operations always produce the same tuple and share
        # one entry, whatever order spglib happened to list them in.
        unique = tuple(sorted({tuple(int(v) for v in r.ravel())
                               for r in rotations}))
        pg = t.pointgroup_international
        expected = _POINT_GROUP_ORDER[pg]
        if len(unique) != expected:
            raise RuntimeError(
                f"space group {number}: point group {pg} has order {expected}, "
                f"but {len(unique)} distinct rotations came back")
        if unique not in operation_index:
            _check_is_a_group(unique, f"point group {pg} (space group {number})")
            operation_index[unique] = len(operation_sets)
            operation_sets.append(unique)

    rows = []
    symbols_seen = set()
    for number in range(1, 231):
        hall, t = by_number[number]
        system = crystal_system(number)
        # international_short's first letter is the centering symbol (P, A,
        # B, C, I, F or R). A/B/C are the same physical Bravais lattice
        # (single-face centered) under different axis labelling -- ITA lists
        # some space groups only in a non-standard A/B setting -- so they are
        # normalised to C, the conventional choice. Trigonal-P and
        # hexagonal-P space groups also share one Bravais lattice (hP), which
        # is why both crystal systems map to the family letter 'h' below --
        # together this leaves exactly 14 distinct Bravais symbols across all
        # 230 groups.
        letter = t.international_short[0]
        if letter in "AB":
            letter = "C"
        bravais_symbol = _FAMILY_LETTER[system] + letter
        symbols_seen.add(bravais_symbol)
        rotations = spglib.get_symmetry_from_database(hall)["rotations"]
        unique = tuple(sorted({tuple(int(v) for v in r.ravel())
                               for r in rotations}))
        rows.append((number, t.international_short, system, letter,
                     bravais_symbol, t.pointgroup_international,
                     operation_index[unique]))

    if len(symbols_seen) != 14:
        raise RuntimeError(
            f"expected 14 Bravais lattices, derived {len(symbols_seen)}: "
            f"{sorted(symbols_seen)}")

    lines = [
        '"""',
        "tools/spacegroups.py -- AUTO-GENERATED by devtools/gen_spacegroups.py. Do not edit",
        "by hand.",
        "",
        "Static International Tables (ITA) space-group reference data for the",
        "Brillouin-zone tool: the Bravais lattice (which drives tools/lattice.py's",
        "'type a space-group number, get the cell constraints' feature) and the",
        "point group's rotation matrices (which the irreducible Brillouin zone is",
        "the fundamental domain of). Generated once from spglib's own database",
        "(see devtools/gen_spacegroups.py); the running app imports only this module, so",
        "it needs no symmetry-search dependency -- there is nothing to search,",
        "since the user supplies the space group directly rather than an atomic",
        "basis for spglib to derive it from.",
        "",
        "The rotations are integer matrices in *conventional-cell fractional*",
        "coordinates, in each group's standard setting (unique axis b for",
        "monoclinic, hexagonal axes for rhombohedral). Turn them into Cartesian",
        "rotations with tools.lattice.point_group_operations rather than applying",
        "them to Cartesian vectors directly.",
        '"""',
        "from __future__ import annotations",
        "",
        "from dataclasses import dataclass",
        "",
        "",
        "@dataclass(frozen=True)",
        "class SpacegroupInfo:",
        "    number: int",
        "    symbol: str              # ITA short Hermann-Mauguin symbol",
        "    crystal_system: str",
        "    centering: str           # P, C, I, F or R (A/B normalised to C)",
        "    bravais_symbol: str      # e.g. 'cF', 'hR', 'mC'",
        "    point_group: str",
        "    operations: int          # index into POINT_GROUP_OPERATIONS",
        "",
        "",
        "SPACEGROUPS = {",
    ]
    for number, symbol, system, letter, bravais_symbol, pg, ops in rows:
        lines.append(
            f"    {number}: SpacegroupInfo({number}, {symbol!r}, {system!r}, "
            f"{letter!r}, {bravais_symbol!r}, {pg!r}, {ops}),")
    lines.append("}")
    lines.append("")
    lines.append("#: the 14 Bravais lattices that occur, for validation/UI listing.")
    lines.append(f"BRAVAIS_SYMBOLS = {sorted({r[4] for r in rows})!r}")
    lines.append("")
    lines.append("#: the distinct point-group rotation sets the 230 space groups use,")
    lines.append("#: each a tuple of 3x3 integer matrices in conventional-cell")
    lines.append("#: fractional coordinates. Indexed by SpacegroupInfo.operations.")
    lines.append("POINT_GROUP_OPERATIONS = (")
    for ops in operation_sets:
        lines.append("    (")
        for flat in ops:
            rows3 = tuple(tuple(flat[i * 3:i * 3 + 3]) for i in range(3))
            lines.append(f"        {rows3!r},")
        lines.append("    ),")
    lines.append(")")
    lines.append("")
    lines.append("")
    lines.append("def spacegroup_info(number: int) -> SpacegroupInfo:")
    lines.append('    """Look up a space group by its ITA number (1-230)."""')
    lines.append("    try:")
    lines.append("        return SPACEGROUPS[int(number)]")
    lines.append("    except KeyError:")
    lines.append("        raise ValueError(")
    lines.append('            f"space group number must be 1-230, got {number}") from None')
    lines.append("")
    lines.append("")
    lines.append("def point_group_rotations(number: int):")
    lines.append('    """The space group\'s point-group rotations, as a tuple of 3x3')
    lines.append("    integer matrices (nested tuples) in conventional-cell fractional")
    lines.append("    coordinates. See tools.lattice.point_group_operations for the")
    lines.append('    Cartesian versions, which is what reciprocal space needs."""')
    lines.append("    return POINT_GROUP_OPERATIONS[spacegroup_info(number).operations]")
    lines.append("")

    out = "\n".join(lines)
    with open("tools/spacegroups.py", "w", encoding="utf-8") as fh:
        fh.write(out)
    print(f"wrote tools/spacegroups.py with {len(rows)} space groups, "
          f"{len(symbols_seen)} Bravais lattices, "
          f"{len(operation_sets)} distinct point-group operation sets")


if __name__ == "__main__":
    sys.exit(main())
