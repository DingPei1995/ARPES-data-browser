"""Which of the main list's actions apply to a selection: the entries of
its right-click menu, the buttons under it, and the buttons of the tool
panel beside it.

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
PLANES = frozenset({"cut"}) | CUBES

#: The filters that work on a cut and on a map alike, but through two
#: different panels -- so a selection has to be one or the other.
PROCESSING = frozenset({"smooth", "derivative", "curvature", "symmetry"})

#: action -> (how many rows: "one", "two", "many" (two or more) or "any", the kinds it takes
#: (None = every known kind), what to call those kinds in a tooltip)
RULES = {
    "open": ("one", None, "a dataset"),
    "info": ("one", None, "a dataset"),
    "rename": ("one", None, "a dataset"),
    "slit_cut": ("any", CUBES, "maps"),
    "deflector_cut": ("any", CUBES, "maps"),
    "figure": ("any", frozenset({"cut"}) | CURVES, "cuts or curves"),
    "stack": ("one", frozenset({"cut"}), "a cut"),
    "mdcfit": ("one", frozenset({"cut"}), "a cut in momentum (k)"),
    "edcfit": ("one", frozenset({"cut"}), "a cut"),
    "corelevel": ("one", frozenset({"cut"}), "a cut"),
    "slices3d": ("one", CUBES, "a map"),
    "notch3d": ("one", CUBES, "a map"),
    "projection3d": ("one", CUBES, "a map"),
    "arithmetic": ("two", frozenset({"cut"}), "two cuts"),
    "combine": ("many", frozenset({"cut"}) | CUBES,
                "two or more cuts or maps of one kind"),
    "truncate": ("any", PLANES | {"spem"}, "cuts, maps or spatial scans"),
    "normalise": ("any", PLANES | {"spem"}, "cuts, maps or spatial scans"),
    "compress": ("any", PLANES | {"spem"}, "cuts, maps or spatial scans"),
    "smooth": ("any", PLANES, "cuts or maps"),
    "derivative": ("any", PLANES, "cuts or maps"),
    "curvature": ("any", PLANES, "cuts or maps"),
    "symmetry": ("any", PLANES, "cuts or maps"),
    "background": ("any", frozenset({"cut"}), "cuts"),
    # Only on the detector's own pixels: a map in momentum has been
    # resampled, and the grid is no longer anywhere fixed to find.
    "fermi": ("one", PLANES, "a cut or a map"),
    "kmap": ("one", frozenset({"map"}),
             "an angle map (one in k or a kz map is not converted this way)"),
    "kcut": ("one", frozenset({"cut"}), "a cut in angle"),
    "kzconv": ("one", frozenset({"kz_map"}),
               "a photon-energy scan (kz map) not yet in momentum"),
    # A Brillouin zone is a statement about momentum-space periodicity, so
    # only a map whose axes are momenta takes one.
    "bz": ("one", frozenset({"k_map", "kz_map_k"}),
           "a map in momentum (a k-map, or a converted kz map)"),
    "kzalign": ("one", frozenset({"kz_map"}),
                "a photon-energy scan (kz map) before its k conversion"),
    "degrid": ("one", frozenset({"cut", "map", "kz_map"}),
               "a cut, or a map or kz map before its k conversion"),
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
    order; a row may carry a third item, whether it is a cut in momentum
    (True / False / None for not known), which the MDC fit needs. ``reason`` is empty when the entry is enabled, and otherwise says
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
    if count == "many" and n < 2:
        return False, f"Works on {wanted}: select at least two rows."
    if kinds == "all":
        return True, ""
    allowed = KNOWN if kinds is None else kinds
    for name, label, *_ in rows:
        if kind_of(label) not in allowed:
            what = label or "not readable"
            return False, (f"Works on {wanted} only; “{name}” is "
                           f"{'a ' if kind_of(label) != 'unknown' else ''}"
                           f"{what}.")
    if action == "combine" and len({kind_of(row[1]) for row in rows}) > 1:
        return False, ("Combines datasets of one kind: cuts with cuts, maps "
                       "with maps, kz maps with kz maps.")
    if action in PROCESSING and len({kind_of(row[1]) == "cut"
                                     for row in rows}) > 1:
        return False, ("Cuts and maps are processed in different panels: "
                       "select cuts only, or maps only.")
    if action == "kcut" and any(len(row) > 2 and row[2] is True for row in rows):
        return False, (f"“{rows[0][0]}” is already in momentum.")
    if action == "mdcfit" and any(len(row) > 2 and row[2] is False
                                  for row in rows):
        return False, ("The MDC fit works on a cut in momentum only; "
                       f"“{rows[0][0]}” is still in angle. Convert it to k "
                       "first (cut viewer: Convert to momentum), or use the "
                       "EDC fit.")
    return True, ""
