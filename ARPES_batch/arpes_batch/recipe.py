"""
A recipe: which files, which reference, which steps, where the output goes.

A JSON file, so that a person can read it, an agent can write it, and a
result can be traced back to exactly what produced it (the recipe is copied
into every run folder). Paths may use ``~`` and ``$VARIABLES`` and are taken
relative to the recipe's own folder, so one recipe works on the laptop
where Google Drive is ``G:\\My Drive`` and on the one where it is
``~/Library/CloudStorage/GoogleDrive-...``::

    {
      "name": "one beamtime's deflector maps",
      "inputs": ["$ARPES_RAW/maps"],
      "references": ["$ARPES_RAW/Au ref"],
      "output": "~/ARPES_processed/my_beamtime",
      "kinds": ["map"],
      "steps": [
        {"op": "degrid"},
        {"op": "calibrate_energy", "source": "reference", "fs_correction": true},
        {"op": "crop", "energy": [-2.5, 0.3]},
        {"op": "normalise"},
        {"op": "save", "as": "angle"},
        {"op": "kconvert", "theta_offset_deg": "auto", "phi_offset_deg": "auto"},
        {"op": "save", "as": "k"}
      ],
      "overrides": {
        "offcentre_map_*": {"kconvert": {"theta_offset_deg": 1.2,
                                          "phi_offset_deg": -0.4,
                                          "azimuth_deg": 30}}
      }
    }

``overrides`` are matched against ``"<file> <entry>"`` with shell-style
patterns; a value keyed by a step name is merged into that step's
parameters, ``"skip": true`` leaves the entry out, ``"reference"`` names a
reference file to use instead of the one chosen automatically, and
``"steps"`` replaces the whole chain for that entry.

``steps`` holds one chain per kind -- ``map``, ``kz_map``, ``cut`` (see
:data:`DEFAULT_CHAINS`); a plain list is the chain for maps. A kz map needs
the ``sample`` block (lattice, surface normal, V0; see :mod:`.sample`).
ANTARES writes a photon-energy scan in the same layout as a deflector map;
name such entries in ``photon_energy_scans`` (patterns) to process them as
kz maps.
"""
from __future__ import annotations

import copy
import fnmatch
import json
import os

from .steps import STEPS

#: One chain per kind of dataset. A recipe may give a plain list, which is
#: taken as the chain for maps (what recipes looked like before kz maps and
#: cuts were handled).
DEFAULT_CHAINS = {
    "map": [
        {"op": "degrid"},
        {"op": "calibrate_energy", "source": "reference", "fs_correction": True,
         "required": True},
        {"op": "crop", "energy": [-2.5, 0.3]},
        {"op": "normalise"},
        {"op": "save", "as": "angle"},
        {"op": "kconvert", "theta_offset_deg": "auto", "phi_offset_deg": "auto",
         "n_kx": 250, "n_ky": 250, "min_centre_score": 0.5},
        {"op": "save", "as": "k"},
    ],
    "kz_map": [
        {"op": "degrid"},
        {"op": "kz_calibrate", "source": "self", "fs_correction": True,
         "required": True},
        {"op": "crop", "energy": [-2.5, 0.2]},
        {"op": "save", "as": "hv"},
        {"op": "kz_match_calculation"},
        {"op": "kz_convert", "angle_offset": "auto", "n_kz": 300, "n_kpar": 300},
        {"op": "save", "as": "kz"},
    ],
    "cut": [
        {"op": "degrid"},
        {"op": "calibrate_energy", "source": "reference", "fs_correction": True,
         "required": True},
        {"op": "crop", "energy": [-2.5, 0.3]},
        {"op": "save", "as": "angle"},
        {"op": "kconvert_cut", "gamma_slit_deg": "auto"},
        {"op": "save", "as": "k"},
        {"op": "curvature"},
        {"op": "save", "as": "curvature"},
    ],
}
DEFAULT_STEPS = DEFAULT_CHAINS["map"]

DEFAULTS = {
    "name": "batch",
    "patterns": ["*.nxs"],
    "kinds": None,                     # None: every kind that has a chain
    "references": [],
    "reference_fit": {"order": 2, "half_width": 3, "step": 4, "window_eV": 0.25},
    "steps": DEFAULT_CHAINS,
    "sample": {},
    "photon_energy_scans": [],         # map entries whose first axis is hv
    "preview": {"energies": [0.0, -0.2, -0.5, -1.0], "width": 0.03},
    "overrides": {},
    "include": [],
    "exclude": [],
    "workers": 1,
}


class RecipeError(ValueError):
    pass


def _expand(path: str, base: str) -> str:
    path = os.path.expandvars(os.path.expanduser(str(path)))
    if "$" in path:
        raise RecipeError(f"unset environment variable in path {path!r}")
    return os.path.normpath(path if os.path.isabs(path) else os.path.join(base, path))


def load(path: str) -> dict:
    with open(path, encoding="utf-8") as fh:
        raw = json.load(fh)
    return normalise(raw, base=os.path.dirname(os.path.abspath(path)), source=path)


def normalise(raw: dict, base: str = ".", source: str = "") -> dict:
    recipe = copy.deepcopy(DEFAULTS)
    recipe.update(copy.deepcopy(raw))
    recipe["_source"] = os.path.abspath(source) if source else ""
    if not recipe.get("inputs"):
        raise RecipeError("a recipe needs 'inputs': the folders or files to process")
    if not recipe.get("output"):
        raise RecipeError("a recipe needs 'output': a folder outside the raw data")
    recipe["inputs"] = [_expand(p, base) for p in recipe["inputs"]]
    recipe["references"] = [_expand(p, base) for p in recipe.get("references") or []]
    recipe["output"] = _expand(recipe["output"], base)
    if isinstance(recipe["steps"], list):
        recipe["steps"] = {"map": recipe["steps"]}
    for kind, chain in recipe["steps"].items():
        if kind not in ("map", "kz_map", "cut"):
            raise RecipeError(f"steps: no chain can be given for {kind!r}; "
                              f"the kinds processed are map, kz_map and cut")
        for index, step in enumerate(chain):
            validate_step(step, f"steps.{kind}[{index}]")
    from .sample import Sample
    try:
        Sample.from_recipe(recipe.get("sample"))
    except (ValueError, TypeError) as exc:
        raise RecipeError(str(exc)) from None
    calc = (recipe.get("sample") or {}).get("calculation")
    if isinstance(calc, str):
        calc = recipe["sample"]["calculation"] = {"file": calc}
    if isinstance(calc, dict):
        for key in ("file", "labels_file"):
            if calc.get(key):
                calc[key] = _expand(calc[key], base)
                if not os.path.exists(calc[key]):
                    raise RecipeError(f"sample.calculation.{key}: {calc[key]} does not exist")
    if recipe.get("kinds") is None:
        recipe["kinds"] = list(recipe["steps"])
    missing = [k for k in recipe["kinds"] if k not in recipe["steps"]]
    if missing:
        raise RecipeError(f"kinds {missing} have no chain under 'steps'")
    for pattern, override in (recipe.get("overrides") or {}).items():
        for key, value in override.items():
            if key == "skip":
                continue
            if key == "reference":
                override[key] = _expand(value, base)
            elif key == "steps":
                for index, step in enumerate(value):
                    validate_step(step, f"overrides[{pattern!r}].steps[{index}]")
            elif key not in STEPS:
                raise RecipeError(f"overrides[{pattern!r}]: {key!r} is not a step "
                                  f"(steps are {sorted(STEPS)})")
    return recipe


def validate_step(step: dict, where: str):
    op = step.get("op")
    if op == "save":
        if not step.get("as"):
            raise RecipeError(f"{where}: a save needs 'as', the name of the checkpoint")
        return
    if op not in STEPS:
        raise RecipeError(f"{where}: unknown op {op!r}; known: save, {', '.join(sorted(STEPS))}")


def select(row: dict, recipe: dict, only=None) -> bool:
    label = f"{row.get('file')} {row.get('entry') or ''}".strip()
    if only and not any(fnmatch.fnmatch(label, p) or fnmatch.fnmatch(row.get("file", ""), p)
                        for p in only):
        return False
    include, exclude = recipe.get("include") or [], recipe.get("exclude") or []
    if include and not any(fnmatch.fnmatch(label, p) for p in include):
        return False
    return not any(fnmatch.fnmatch(label, p) for p in exclude)


def steps_for(row: dict, recipe: dict) -> tuple:
    """``(steps, override)`` for one entry, overrides applied in the order
    they are written (a later, more specific pattern wins)."""
    label = f"{row.get('file')} {row.get('entry') or ''}".strip()
    steps = copy.deepcopy(recipe["steps"].get(row.get("kind"), []))
    merged = {}
    for pattern, override in (recipe.get("overrides") or {}).items():
        if not (fnmatch.fnmatch(label, pattern) or fnmatch.fnmatch(row.get("file", ""), pattern)):
            continue
        merged.update(copy.deepcopy(override))
        if "steps" in override:
            steps = copy.deepcopy(override["steps"])
        for op, params in override.items():
            if op in STEPS:
                for step in steps:
                    if step.get("op") == op:
                        step.update(params)
    return steps, merged
