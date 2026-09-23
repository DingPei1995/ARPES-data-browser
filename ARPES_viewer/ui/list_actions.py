"""Which entries of the main list's right-click menu apply to a selection.

Decided from the list's own records -- the kind each row was listed as --
without reading any data, so the menu opens instantly however large the
selected files are. An entry that does not fit is greyed out, and its
tooltip says why; with several rows selected, one row that does not fit is
enough to grey it out.

Qt-free on purpose, so the rules can be tested on their own.
"""
from __future__ import annotations

#: The label a row is listed with -> the dataset kind it stands for.
LABEL_KINDS = {
    "cut": "cut",
    "map": "map",
    "k-map": "k_map",
    "kz map": "kz_map",
    "kz map (k)": "kz_map_k",
    "spem": "spem",
    "edc": "edc",
    "mdc": "mdc",
    "spin edc": "spin_edc",
}

CUBES = frozenset({"map", "k_map", "kz_map", "kz_map_k"})
CURVES = frozenset({"edc", "mdc", "spin_edc"})
KNOWN = frozenset(LABEL_KINDS.values())

#: action -> (how many rows: "one", "two" or "any", the kinds it takes
#: (None = every known kind), what to call those kinds in a tooltip)
RULES = {
    "open": ("one", None, "a dataset"),
    "info": ("one", None, "a dataset"),
    "rename": ("one", None, "a dataset"),
    "slit_cut": ("any", CUBES, "maps"),
    "deflector_cut": ("any", CUBES, "maps"),
    "figure": ("any", frozenset({"cut"}) | CURVES, "cuts or curves"),
    "stack": ("one", frozenset({"cut"}), "a cut"),
    "fit": ("one", frozenset({"cut"}), "a cut"),
    "view3d": ("one", CUBES, "a map"),
    "arithmetic": ("two", frozenset({"cut"}), "two cuts"),
    "save": ("any", None, "datasets"),
    "remove": ("any", "all", "datasets"),
    "log": ("none", "all", ""),
}


def kind_of(label) -> str:
    """The dataset kind for a row's listed kind label ("Map" -> "map");
    "unknown" for anything this program cannot open."""
    return LABEL_KINDS.get(str(label or "").strip().lower(), "unknown")


def availability(action: str, rows):
    """``(enabled, reason)`` for one menu entry.

    ``rows`` is ``[(name, kind label), ...]`` for the selected rows, in list
    order. ``reason`` is empty when the entry is enabled, and otherwise says
    what is wrong in terms of the selection -- it becomes the tooltip.
    """
    count, kinds, wanted = RULES[action]
    n = len(rows)
    if count == "none":
        return True, ""
    if n == 0:
        return False, "Select a dataset in the list first."
    if count == "one" and n != 1:
        return False, f"Works on {wanted}: select one row only."
    if count == "two" and n != 2:
        return False, f"Works on {wanted}: select exactly two rows."
    if kinds == "all":
        return True, ""
    allowed = KNOWN if kinds is None else kinds
    for name, label in rows:
        if kind_of(label) not in allowed:
            what = label or "not readable"
            return False, (f"Works on {wanted} only; “{name}” is "
                           f"{'a ' if kind_of(label) != 'unknown' else ''}"
                           f"{what}.")
    return True, ""
