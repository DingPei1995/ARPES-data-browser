"""
What only the user knows about the sample, and asking for it.

A photon-energy scan cannot be put into k_z from the data alone. Three
things have to come from the person who grew or chose the crystal:

* **the lattice** -- a, b, c, the angles and the space group (the space
  group fixes the centering, and the centering decides the k_z period: a
  body-centred lattice repeats every 4*pi/c along [001], not 2*pi/c);
* **the surface normal** -- which planes the crystal cleaved along, as
  Miller indices against the conventional cell;
* **the inner potential V0** -- a number; ``"scan"`` to have it chosen so
  that the data repeats with the lattice's own period; or ``"calculation"``
  to take it from a band-structure calculation matched to the scan.

And one thing that is asked for but may be declined:

* **a band-structure calculation** (``calculation``) -- a band file along a
  path through Gamma and the zone boundary on the surface normal (Gamma-A,
  Gamma-Z). If there is one, the scan at normal emission is matched to it
  (:mod:`arpes_batch.calcbands`), which says which photon energies reach
  the high-symmetry planes, and gives V0 and the calculation's energy
  offset. ``"calculation": null`` records that there is none.

None of these is guessed. A recipe that converts a kz map without them is
refused before anything runs, with a list of what is missing, what each one
is for, and an example -- printed and prompted for at a terminal, returned
as ``needs_input`` to an agent, which is to ask the user rather than invent
a crystal.

The recipe's ``sample`` block::

    "sample": {
      "name": "BaFe2As2",
      "lattice": {"a": 3.96, "b": 3.96, "c": 13.02,
                  "alpha": 90, "beta": 90, "gamma": 90, "space_group": 139},
      "surface_normal": [0, 0, 1],
      "inner_potential_eV": "scan",
      "calculation": {"file": "band.dat", "path": ["G", "M", "K", "G", "A", "L", "H", "A"],
                      "fermi_energy_eV": 0.0},
      "work_function_eV": null,
      "effective_mass": 1.0
    }
"""
from __future__ import annotations

import sys
from dataclasses import dataclass

import numpy as np

from tools import cleavage
from tools.lattice import (LatticeParams, conventional_vectors, reciprocal_vectors,
                           validate_lattice_parameters)

#: What a kz conversion needs from the user, in the order it is asked for.
#: (field, question, example, why)
KZ_INPUTS = (
    ("sample.lattice",
     "Lattice constants and space group: a, b, c in angstrom, alpha, beta, "
     "gamma in degrees, space group 1-230",
     {"a": 3.96, "b": 3.96, "c": 13.02, "alpha": 90, "beta": 90, "gamma": 90,
      "space_group": 139},
     "fixes the reciprocal lattice, hence the k_z period and where the zone "
     "boundaries are drawn"),
    ("sample.surface_normal",
     "Surface normal (cleavage plane) as Miller indices [h, k, l] of the "
     "conventional cell",
     [0, 0, 1],
     "k_z is measured along it; the period is the shortest reciprocal-"
     "lattice vector in that direction (for a centred lattice not always "
     "(001) but (002))"),
    ("sample.inner_potential_eV",
     "Inner potential V0 in eV; or \"scan\" to estimate it from the data's "
     "k_z period; or \"calculation\" to take it from a band calculation "
     "matched to the scan",
     "scan",
     "sets k_z = sqrt(A (E_kin + V0) - k_par^2); a few eV off shifts the "
     "whole k_z axis"),
)

#: Asked once, and may be answered with null.
CALCULATION_INPUT = (
    "sample.calculation",
    "Do you have a band-structure calculation (e.g. DFT) to compare with? It "
    "decides which photon energies reach the high-symmetry planes. Give the band "
    "file (two columns per band, or k plus one column per band), the path labels "
    "if there is no labels file next to it, and E_F if its energies are absolute "
    "-- or null if there is none.",
    {"file": "band.dat", "path": ["G", "M", "K", "G", "A", "L", "H", "A"],
     "fermi_energy_eV": 0.0},
    "matching the scan at normal emission to the calculated dispersion along "
    "Gamma - zone boundary gives V0, the calculation's energy offset, and the "
    "photon energies of the Gamma and boundary planes")


class NeedsInput(Exception):
    """The recipe cannot run until the user supplies these."""

    def __init__(self, missing: list):
        self.missing = missing
        names = ", ".join(m["field"] for m in missing)
        super().__init__(f"the recipe needs input from the user: {names}")

    def to_json(self) -> dict:
        return {"needs_input": self.missing,
                "message": "These describe the sample and cannot be taken from the "
                           "data. Ask the user; do not guess them. Items marked "
                           "optional may be answered with null (e.g. no band "
                           "calculation), which is then recorded in the recipe."}


_UNASKED = "__not asked__"


@dataclass
class Sample:
    name: str = ""
    calculation: object = _UNASKED      # dict, None (there is none), or not asked yet
    lattice: LatticeParams = None
    surface_normal: tuple = None
    inner_potential: object = None        # float, or "scan"
    work_function: float = None
    effective_mass: float = 1.0

    @classmethod
    def from_recipe(cls, block: dict) -> "Sample":
        block = dict(block or {})
        lattice = block.get("lattice") or {}
        if not isinstance(lattice, dict):
            raise ValueError(f"sample.lattice must be an object such as "
                             f"{KZ_INPUTS[0][2]}, not {lattice!r}")
        normal = block.get("surface_normal")
        if normal is not None and (not isinstance(normal, (list, tuple)) or len(normal) != 3):
            raise ValueError(f"sample.surface_normal must be three Miller indices such as "
                             f"[0, 0, 1], not {normal!r}")
        v0 = block.get("inner_potential_eV")
        if v0 is not None and v0 not in ("scan", "calculation") and not isinstance(
                v0, (int, float)):
            raise ValueError(f"sample.inner_potential_eV must be a number, \"scan\" or "
                             f"\"calculation\", not {v0!r}")
        params = None
        if lattice.get("a") and lattice.get("c"):
            params = LatticeParams(
                a=float(lattice["a"]), b=float(lattice.get("b") or lattice["a"]),
                c=float(lattice["c"]), alpha=float(lattice.get("alpha", 90)),
                beta=float(lattice.get("beta", 90)), gamma=float(lattice.get("gamma", 90)),
                space_group=(int(lattice["space_group"]) if lattice.get("space_group")
                             else None),
                centering=lattice.get("centering"))
        calculation = block.get("calculation", _UNASKED)
        if calculation not in (None, _UNASKED):
            if isinstance(calculation, str):
                calculation = {"file": calculation}
            if not isinstance(calculation, dict) or not calculation.get("file"):
                raise ValueError(f"sample.calculation must be null or an object with a "
                                 f"'file', such as {CALCULATION_INPUT[2]}")
        return cls(name=str(block.get("name") or ""), lattice=params,
                   calculation=calculation,
                   surface_normal=tuple(int(v) for v in normal) if normal else None,
                   inner_potential=(v0 if v0 in (None, "scan", "calculation") else float(v0)),
                   work_function=(float(block["work_function_eV"])
                                  if block.get("work_function_eV") is not None else None),
                   effective_mass=float(block.get("effective_mass") or 1.0))

    def problems(self) -> list:
        """Warnings about the numbers given (not refusals: which of the
        lengths or the space group is wrong is the user's call)."""
        if self.lattice is None:
            return []
        try:
            return list(validate_lattice_parameters(self.lattice))
        except Exception as exc:                          # noqa: BLE001
            return [str(exc)]

    def surface(self) -> dict:
        """The k_z period along the surface normal, and what it is."""
        return surface_period(self.lattice, self.surface_normal)


def missing_inputs(sample: Sample, needs_work_function: bool = False) -> list:
    """What a kz conversion still needs, as ``[{field, question, example,
    why}]``; empty when it can run."""
    have = {"sample.lattice": sample.lattice is not None
            and (sample.lattice.space_group or sample.lattice.centering),
            "sample.surface_normal": sample.surface_normal is not None,
            "sample.inner_potential_eV": sample.inner_potential is not None}
    out = [{"field": f, "question": q, "example": e, "why": w}
           for f, q, e, w in KZ_INPUTS if not have[f]]
    if sample.calculation == _UNASKED:
        f, q, e, w = CALCULATION_INPUT
        out.append({"field": f, "question": q, "example": e, "why": w, "optional": True,
                    "decline_with": None})
    elif sample.inner_potential == "calculation" and not sample.calculation:
        out.append({"field": "sample.calculation",
                    "question": "inner_potential_eV is \"calculation\" but no calculation "
                                "is given: " + CALCULATION_INPUT[1],
                    "example": CALCULATION_INPUT[2], "why": CALCULATION_INPUT[3]})
    if needs_work_function and sample.work_function is None:
        out.append({"field": "sample.work_function_eV",
                    "question": "Analyser work function in eV (E_F,kin = hv - W)",
                    "example": 4.4,
                    "why": "no gold reference matches these analyser settings and "
                           "the file does not record one"})
    if sample.lattice is not None and sample.surface_normal is not None:
        try:
            surface_period(sample.lattice, sample.surface_normal)
        except ValueError as exc:
            out.append({"field": "sample.surface_normal", "question": str(exc),
                        "example": [0, 0, 1], "why": "not a lattice plane"})
    return out


def surface_period(lattice: LatticeParams, hkl) -> dict:
    """The shortest reciprocal-lattice vector along the normal of ``(hkl)``.

    ``hkl`` is against the conventional cell, which is how planes are
    named. The answer is taken from the *primitive* reciprocal lattice
    (``tools.cleavage.reciprocal_lengths``), so centering is right by
    construction: (001) of a body-centred cell gives (002).
    """
    if lattice is None or hkl is None:
        raise ValueError("lattice and surface normal are both needed")
    b_conv = reciprocal_vectors(conventional_vectors(lattice))
    g = np.asarray(hkl, dtype=float) @ b_conv
    if not np.linalg.norm(g):
        raise ValueError(f"({' '.join(map(str, hkl))}) is not a direction")
    unit = g / np.linalg.norm(g)
    for found, length, direction in cleavage.reciprocal_lengths(lattice, max_index=6):
        if abs(abs(float(np.dot(direction, unit))) - 1.0) < 1e-6:
            note = ""
            reduced = tuple(int(v) for v in hkl)
            if tuple(abs(v) for v in found) != tuple(abs(v) for v in reduced):
                note = (f"centred lattice: the first reciprocal-lattice vector along "
                        f"({' '.join(map(str, reduced))}) is "
                        f"({' '.join(map(str, found))}), so the k_z period is "
                        f"{length:.4f} A^-1")
            return {"hkl": list(reduced), "first_reflection": list(found),
                    "period_invA": float(length), "spacing_A": float(2 * np.pi / length),
                    "note": note}
    raise ValueError(f"no reciprocal-lattice vector found along "
                     f"({' '.join(map(str, hkl))}) within |index| <= 6")


def list_surfaces(lattice: LatticeParams, limit: int = 12) -> list:
    """The low-index surfaces and their k_z periods, shortest first -- the
    choice a user is making when they name a cleavage plane."""
    out = []
    for h, L, _u in cleavage.reciprocal_lengths(lattice, max_index=3)[:limit]:
        divisor = int(np.gcd.reduce([abs(int(v)) for v in h])) or 1
        out.append({"surface": [int(v) // divisor for v in h], "first_reflection": list(h),
                    "period_invA": float(L), "spacing_A": float(2 * np.pi / L)})
    return out


def prompt(missing: list, stream_in=None, stream_out=None) -> dict:
    """Ask for each missing value at a terminal. Returns ``{field: value}``;
    an empty answer leaves that field missing."""
    import json
    stream_in = stream_in or sys.stdin
    stream_out = stream_out or sys.stderr
    answers = {}
    print("\nThis recipe converts a photon-energy scan to k_z, which needs "
          "information about the sample that is not in the data.", file=stream_out)
    for item in missing:
        optional = ("\n  (optional: press Enter or type null if there is none)"
                    if item.get("optional") else "")
        print(f"\n{item['field']}: {item['question']}\n  (why: {item['why']})\n"
              f"  example: {json.dumps(item['example'])}{optional}", file=stream_out)
        stream_out.write("  > ")
        stream_out.flush()
        text = stream_in.readline().strip()
        if not text:
            if item.get("optional"):
                answers[item["field"]] = None       # asked, and there is none
            continue
        try:
            answers[item["field"]] = json.loads(text)
        except json.JSONDecodeError:
            answers[item["field"]] = text
    return answers


def apply_answers(recipe_sample: dict, answers: dict) -> dict:
    block = dict(recipe_sample or {})
    for field, value in answers.items():
        key = field.split(".", 1)[1]
        block[key] = value
    return block
