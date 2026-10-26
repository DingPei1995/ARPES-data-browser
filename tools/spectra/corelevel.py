"""
tools/spectra/corelevel.py
==========================
A quick core-level check of a cut, as plain numpy/scipy with no Qt: which
element and orbital each peak of the angle-integrated EDC comes from.

The steps, each its own function so the window (``ui/panels/corelevel.py``) can let
the user confirm or change what came out of the one before:

1. :func:`integrated_edc` -- the whole cut summed over its angle axis.
2. :func:`detect_scale` and :func:`locate_fermi_edge` -- what the energy
   axis is (kinetic, ``E - E_F``, or binding) and where the Fermi level sits
   on it. On a kinetic axis E_F is ``hv - W``; when the window holds the
   edge it is fitted (``tools.spectra.fermi``) rather than trusted, because the
   monochromator's nominal photon energy is routinely off by more than the
   0.5 eV the comparison allows. :class:`EnergyReference` turns the axis
   into binding energies once the user has agreed the numbers.
3. :func:`find_peaks` -- the peaks of the EDC on the binding-energy scale,
   each refined by a local Gaussian fit.
4. :func:`assign` -- compared with the tabulated binding energies
   (:mod:`tools.spectra.corelevel_table`, Table 1-1 of the X-ray Data Booklet) of
   the elements the sample is made of. Within the tolerance (0.5 eV) a
   peak is the sample's and gets its orbital; further away but within a
   chemical shift it is left to the user, with the nearest lines named; a
   peak nothing in the sample explains gets the lines of common elements
   (non-radioactive, not highly toxic) near it, as possible contamination.

What a check by table alone would miss, and so is checked as well:

* **Gas-phase references.** The table gives H, N, O, F, Cl and the rare
  gases relative to the *vacuum* level of the free molecule -- O 1s 543.1,
  where every oxide is at ~530 eV. :data:`CHEMICAL_STATES` adds the usual
  solid-state values for those (and the common oxides, carbon species and
  salts), so an oxygen peak is recognised as one.
* **Spin-orbit doublets.** A p, d or f level is two peaks, split by a known
  amount with intensities near the (2j+1) ratio (1:2, 2:3, 3:4). Finding the
  partner confirms an assignment; a partner that should be in the window
  and is not is reported.
* **Higher-order light.** A grating also passes 2hv and 3hv. A level deeper
  than hv - W shows up at ``E_B - hv`` (second order) or ``E_B - 2hv``
  (third), looking like a core level that is not there.
* **Auger lines** sit at a fixed *kinetic* energy, so they move on the
  binding scale when hv changes; the commonest (C, N, O KLL...) are named.
* **A rigid shift.** Charging of an insulating sample, or a wrong E_F or hv,
  moves every peak by the same amount; :func:`suggest_shift` finds the
  shift that would make the most peaks fit and says so instead of quietly
  applying it.
* **What is not there.** Lines of the sample inside the window with no peak
  are listed, and so are elements with nothing at all in the window (with
  their nearest reachable levels at this photon energy).
"""
from __future__ import annotations

import re
import warnings
from dataclasses import dataclass, field

import numpy as np

from tools.spectra.corelevel_table import BINDING_ENERGIES, ELEMENTS

#: A peak within this distance of a tabulated level (eV) is the sample's.
TOLERANCE_EV = 0.5
#: Further than the tolerance but within this, a peak may still be the
#: sample's, chemically shifted; the user decides.
UNCERTAIN_EV = 2.0
#: How far a common element's level may be from an unexplained peak and
#: still be offered as a possible contaminant -- wider than the tolerance,
#: because a contaminant's chemical state is not known.
CONTAMINANT_EV = 1.0
#: How far (in noise sigmas, peak to surrounding valley) a maximum of the
#: smoothed EDC must stand to count as a peak. Peak-to-valley, so noise
#: alone reaches 5-6 over a long scan.
SENSITIVITY = 8.0
#: The smallest peak kept, as a fraction of the strongest one's prominence.
#: Small, because the noise test above is what rejects noise: a 1 % cut
#: threw away real weak lines next to a strong one (a 3d at 10^6 counts
#: hid a 2000-count 4p standing 40 sigma clear).
MIN_RELATIVE = 0.002
#: How far either side of a position the hints for it look (eV): a peak
#: added by hand, or the cursor, is named with every level this close.
HINT_EV = 2.0
#: The work function assumed when the file records none.
DEFAULT_WORK_FUNCTION = 4.5

#: Elements the table gives against the vacuum level (gas phase).
VACUUM_REFERENCED = frozenset({"H", "He", "N", "O", "F", "Ne", "Cl", "Ar",
                               "Kr", "Xe", "Rn"})
#: Not offered as contaminants: radioactive ...
RADIOACTIVE = frozenset({"Tc", "Pm", "Po", "At", "Rn", "Fr", "Ra", "Ac",
                         "Th", "Pa", "U"})
#: ... or highly toxic -- nobody puts these in a UHV chamber by accident.
HIGHLY_TOXIC = frozenset({"Be", "As", "Cd", "Hg", "Tl", "Pb"})
#: Hydrogen and helium have no core level an XPS spectrum shows.
NO_CORE_LEVEL = frozenset({"H", "He"})

#: Contaminants seen in practice, most likely first, with where they come
#: from. Every other element that is neither radioactive nor highly toxic is
#: still a candidate, ranked after these.
LIKELY_CONTAMINANTS = {
    "C": "adventitious carbon, carbon tape, pump oil",
    "O": "oxidation, adsorbed water / CO",
    "N": "air exposure, organic residue",
    "F": "PTFE, fluorinated grease, etching",
    "Cl": "handling, salts, etchants",
    "Na": "handling (salt), glass",
    "K": "handling, glass, alkali dosing",
    "S": "segregation from the bulk, vacuum grease",
    "Si": "silicone grease, glass, SiO2 substrates",
    "Ca": "segregation, handling",
    "Ar": "implanted by sputtering",
    "Ta": "sample plates and clamps",
    "Mo": "sample plates and clamps",
    "W": "filaments, sample plates",
    "Cu": "sample holders, Cu tape",
    "Ag": "silver epoxy / paint",
    "In": "indium foil or solder",
    "Au": "gold plates, reference foil",
    "Fe": "stainless steel",
    "Ni": "stainless steel",
    "Cr": "stainless steel",
    "Ti": "Ti sublimation pump",
    "Al": "sample plates, Al2O3 / sapphire",
    "Mg": "segregation, MgO substrates",
    "Zn": "brass parts",
    "Sn": "solder",
    "Ga": "Ga ion beam (FIB), GaAs",
    "B": "boron nitride (heaters, crucibles)",
    "P": "phosphates, segregation",
    "Pt": "Pt foil, thermocouples",
}

#: Typical solid-state binding energies (eV, Fermi level) for chemical
#: states the table does not cover: (element, orbital, eV, what). The usual
#: handbook values (Moulder et al., *Handbook of XPS*; NIST SRD 20), to
#: within the few tenths their spread is anyway.
CHEMICAL_STATES = (
    ("C", "1s", 284.8, "adventitious C (C-C/C-H)"),
    ("C", "1s", 286.3, "C-O"),
    ("C", "1s", 288.6, "C=O / O-C=O"),
    ("C", "1s", 289.8, "carbonate"),
    ("C", "1s", 292.5, "CF2 (PTFE)"),
    ("O", "1s", 529.8, "metal oxide (O2-)"),
    ("O", "1s", 531.4, "hydroxide / defective oxide"),
    ("O", "1s", 532.6, "adsorbed water, C-O, SiO2"),
    ("N", "1s", 397.0, "nitride"),
    ("N", "1s", 399.8, "organic N (C-N)"),
    ("F", "1s", 684.8, "metal fluoride"),
    ("F", "1s", 689.0, "C-F (PTFE)"),
    ("Cl", "2p3/2", 198.8, "chloride"),
    ("Cl", "2p1/2", 200.4, "chloride"),
    ("S", "2p3/2", 161.8, "sulfide"),
    ("S", "2p1/2", 163.0, "sulfide"),
    ("S", "2p3/2", 168.8, "sulfate"),
    ("Si", "2p", 103.3, "SiO2"),
    ("P", "2p", 133.5, "phosphate"),
    ("Na", "1s", 1071.8, "Na+ salt"),
    ("K", "2p3/2", 292.9, "K+ salt"),
    ("K", "2p1/2", 295.7, "K+ salt"),
    ("Ca", "2p3/2", 347.2, "Ca2+"),
    ("Ar", "2p3/2", 241.8, "implanted Ar"),
    ("Ar", "2p1/2", 243.9, "implanted Ar"),
    ("Al", "2p", 74.6, "Al2O3"),
    ("Ti", "2p3/2", 458.8, "TiO2"),
    ("Cr", "2p3/2", 576.8, "Cr2O3"),
    ("Fe", "2p3/2", 710.9, "Fe oxide"),
    ("Ni", "2p3/2", 854.0, "NiO"),
    ("Cu", "2p3/2", 933.6, "CuO"),
    ("Mo", "3d5/2", 232.6, "MoO3"),
    ("Ta", "4f7/2", 26.6, "Ta2O5"),
    ("W", "4f7/2", 35.8, "WO3"),
    ("In", "3d5/2", 444.6, "In2O3"),
    ("Sn", "3d5/2", 486.6, "SnO2"),
)

#: The strongest Auger lines of the common elements, at their kinetic
#: energy (eV), and how far that moves with the chemical state.
AUGER_LINES = (
    ("C", "KLL", 263.0), ("N", "KLL", 379.0), ("O", "KLL", 510.0),
    ("F", "KLL", 656.0), ("Na", "KLL", 990.0),
)
AUGER_EV = 6.0

#: The (2j+1) intensity ratio of the weaker to the stronger component of a
#: spin-orbit doublet, by angular momentum letter.
DOUBLET_RATIO = {"p": 0.5, "d": 2.0 / 3.0, "f": 0.75}


# --------------------------------------------------------------------------
# Elements and their lines
# --------------------------------------------------------------------------
SYMBOLS = frozenset(BINDING_ENERGIES)
_SYMBOL = re.compile(r"[A-Z][a-z]?")


def parse_elements(text: str):
    """``(elements, unknown)`` from what the user typed.

    Accepts symbols separated by anything (``Nb, S, Cl``), a formula
    (``NbSCl``, ``Bi2Se3``, ``(Ba,K)Fe2As2``) or lower-case symbols
    separated by spaces or commas (``nb s cl``). Each element appears once,
    in the order typed; anything that is not a symbol is returned in
    ``unknown`` so the dialog can say so instead of dropping it.
    """
    text = str(text or "")
    elements, unknown = [], []
    for word in re.split(r"[\s,;/+]+", text.strip()):
        if not word:
            continue
        if word[0].islower() and word.isalpha():
            # a lone lower-case symbol: "nb", "s"
            candidate = word.capitalize()
            found = [candidate] if candidate in SYMBOLS else None
            leftover = "" if found else word
        else:
            found = _SYMBOL.findall(word)
            leftover = re.sub(r"[A-Z][a-z]?|[\d.()\[\]{}x\-δ]", "", word)
        if leftover:
            unknown.append(word)
            continue
        for symbol in found:
            if symbol not in SYMBOLS:
                unknown.append(symbol)
            elif symbol not in elements:
                elements.append(symbol)
    return elements, unknown


@dataclass(frozen=True)
class Line:
    """One reference level: ``element`` ``orbital`` at ``energy`` eV."""
    element: str
    orbital: str
    energy: float
    #: "table" (Table 1-1) or the chemical state it stands for
    source: str = "table"

    @property
    def name(self) -> str:
        return f"{self.element} {self.orbital}"

    @property
    def label(self) -> str:
        if self.source == "table":
            return self.name
        return f"{self.name} ({self.source})"

    @property
    def shell(self) -> str:
        """``2p3/2`` -> ``2p``."""
        return self.orbital.split("/")[0][:2] if "/" in self.orbital \
            else self.orbital

    def partner_orbital(self):
        """The other component of the doublet, or None for an s level or a
        level the table gives only once."""
        if "/" not in self.orbital:
            return None
        n, letter = self.orbital[0], self.orbital[1]
        l_value = "spdf".index(letter)
        j2 = int(self.orbital[2:].split("/")[0])
        other = 2 * l_value - 1 if j2 == 2 * l_value + 1 else 2 * l_value + 1
        return f"{n}{letter}{other}/2"

    def is_major(self) -> bool:
        """The stronger component (j = l + 1/2) or a single level."""
        if "/" not in self.orbital:
            return True
        l_value = "spdf".index(self.orbital[1])
        return int(self.orbital[2:].split("/")[0]) == 2 * l_value + 1


def element_lines(symbol: str, *, chemical_states: bool = True) -> list:
    """Every level of one element: the table's, with the two components of
    a doublet the table gives at the same energy merged into one (``3p``),
    and -- for ``chemical_states`` -- the typical solid-state values."""
    levels = BINDING_ENERGIES.get(symbol, {})
    lines, done = [], set()
    for orbital, energy in levels.items():
        if orbital in done:
            continue
        line = Line(symbol, orbital, float(energy))
        partner = line.partner_orbital()
        if partner in levels and abs(levels[partner] - energy) < 0.05:
            line = Line(symbol, line.shell, float(energy))
            done.add(partner)
        lines.append(line)
        done.add(orbital)
    if chemical_states:
        lines += [Line(e, o, v, what) for e, o, v, what in CHEMICAL_STATES
                  if e == symbol]
    return lines


def reference_lines(elements, *, chemical_states: bool = True,
                    window=None) -> list:
    """The lines of ``elements``, by energy; only those inside ``window``
    (lo, hi in eV binding) when given."""
    lines = []
    for symbol in elements:
        lines += element_lines(symbol, chemical_states=chemical_states)
    if window is not None:
        lo, hi = window
        lines = [line for line in lines if lo <= line.energy <= hi]
    return sorted(lines, key=lambda line: line.energy)


def common_elements(exclude=()) -> list:
    """Elements offered as possible contaminants: neither radioactive nor
    highly toxic, likeliest first. ``exclude`` (the sample's own) is left
    out."""
    exclude = set(exclude)
    banned = RADIOACTIVE | HIGHLY_TOXIC | NO_CORE_LEVEL | exclude
    likely = [s for s in LIKELY_CONTAMINANTS if s not in banned]
    rest = [s for z, s in sorted(ELEMENTS.items())
            if s not in banned and s not in LIKELY_CONTAMINANTS]
    return likely + rest


# --------------------------------------------------------------------------
# The EDC and its energy scale
# --------------------------------------------------------------------------
def integrated_edc(frame, energy, energy_dim: int = 1):
    """The cut summed over its angle axis: ``(energy, intensity)``, with
    the energy ascending. NaNs (outside a converted cut's footprint) are
    skipped rather than poisoning the sum."""
    frame = np.asarray(frame, dtype=float)
    energy = np.asarray(energy, dtype=float)
    if frame.ndim != 2:
        raise ValueError("a cut is two-dimensional")
    intensity = np.nansum(frame, axis=1 - energy_dim)
    if intensity.size != energy.size:
        raise ValueError("the energy axis does not match the data")
    order = np.argsort(energy)
    return energy[order], intensity[order]


#: The three energy scales an axis can be on.
SCALES = ("kinetic", "relative", "binding")
SCALE_LABELS = {"kinetic": "Kinetic energy",
                "relative": "E - E_F (negative below E_F)",
                "binding": "Binding energy (positive below E_F)"}


def detect_scale(energy, label: str = "", photon_energy=None) -> str:
    """Which scale ``energy`` is on, from its label first and its numbers
    second: ``"kinetic"``, ``"relative"`` (E - E_F) or ``"binding"``."""
    text = str(label or "").lower().replace(" ", "")
    energy = np.asarray(energy, dtype=float)
    median = float(np.nanmedian(energy)) if energy.size else 0.0
    if "bind" in text or text.startswith("be(") or text.startswith("e_b"):
        return "relative" if median < 0 else "binding"
    if "e-e_f" in text or "e-ef" in text or "e_f" in text:
        return "relative"
    if "kin" in text:
        return "kinetic"
    if photon_energy and np.isfinite(photon_energy):
        return "kinetic" if median > 0.5 * float(photon_energy) else "relative"
    return "relative" if float(np.nanmin(energy)) < 0 else "kinetic"


def info_number(info: dict, keys, default=None):
    """The first of ``keys`` in ``info`` that reads as a finite number."""
    for key in keys:
        if key not in info:
            continue
        try:
            value = float(np.asarray(info[key], dtype=float).reshape(-1)[0])
        except (TypeError, ValueError, IndexError):
            continue
        if np.isfinite(value):
            return value
    return default


def photon_energy_of(info: dict):
    """The photon energy a dataset recorded (eV), or None."""
    value = info_number(info, ("photon_energy_eV", "PhotonEnergy",
                               "cassiopee.photon_energy_eV", "hv", "mono.energy"))
    return value if value is not None and value > 0 else None


def work_function_of(info: dict) -> float:
    return info_number(info, ("work_function_eV", "cassiopee.work_function_eV"),
                       DEFAULT_WORK_FUNCTION)


@dataclass
class FermiEstimate:
    """Where E_F is on the axis, and how that was decided."""
    level: float
    #: True when the edge was found in the data, False when E_F was put at
    #: hv - W (or 0 on a referred axis) because the window does not hold it
    found: bool
    message: str
    error: float = float("nan")
    resolution: float = float("nan")
    fit: object = None


def expected_fermi_level(scale: str, photon_energy=None,
                         work_function: float = DEFAULT_WORK_FUNCTION):
    """Where E_F should be on an axis of ``scale``: hv - W on a kinetic one,
    zero on the other two; None for a kinetic axis without hv."""
    if scale != "kinetic":
        return 0.0
    if photon_energy is None:
        return None
    return float(photon_energy) - float(work_function)


def locate_fermi_edge(energy, intensity, scale: str, *, photon_energy=None,
                      work_function: float = DEFAULT_WORK_FUNCTION,
                      temperature: float = 30.0,
                      search: float = 2.0) -> FermiEstimate:
    """Find the Fermi edge of the integrated EDC.

    The edge is looked for within ``search`` eV of where it should be (hv -
    W, or zero on a referred axis); without a photon energy, in the
    highest-energy fifth of a kinetic window, where the edge is the last
    fall. A candidate counts as the edge only if the intensity really drops
    across it, and is then fitted with the broadened Fermi function. When
    the window does not reach E_F -- the usual case for a core-level scan --
    the expected position is returned with ``found`` False.
    """
    from tools.spectra.fermi import (fit_fermi_edge, initial_guess,
                                     steepest_drop)

    e = np.asarray(energy, dtype=float)
    i = np.asarray(intensity, dtype=float)
    keep = np.isfinite(e) & np.isfinite(i)
    e, i = e[keep], i[keep]
    order = np.argsort(e)
    e, i = e[order], i[order]
    expected = expected_fermi_level(scale, photon_energy, work_function)
    step = float(np.median(np.abs(np.diff(e)))) if e.size > 1 else 0.0
    if e.size < 8 or step <= 0:
        level = expected if expected is not None else float("nan")
        return FermiEstimate(level, False, "too few points to look for an edge")

    sign = 1.0 if scale != "binding" else -1.0     # binding: E_F is a rise
    if expected is not None:
        lo, hi = expected - search, expected + search
        if hi < e[0] or lo > e[-1]:
            what = ("hν − W" if scale == "kinetic" else "zero")
            return FermiEstimate(
                expected, False,
                f"E_F is not inside this window ({e[0]:.2f} to {e[-1]:.2f}); "
                f"it is put at {what} = {expected:.3f}.")
    else:
        lo, hi = e[-1] - 0.2 * (e[-1] - e[0]), e[-1]
    inside = (e >= lo) & (e <= hi)
    if inside.sum() < 8:
        level = expected if expected is not None else float(e[-1])
        return FermiEstimate(level, False,
                             "the window holds too few points near E_F to "
                             "fit an edge")
    candidate = steepest_drop(e[inside], sign * i[inside])

    half = max(0.4, 6 * step)
    below = (e >= candidate - 3 * half) & (e < candidate - 0.5 * half)
    above = (e > candidate + 0.5 * half) & (e <= candidate + 3 * half)
    if scale == "binding":
        below, above = above, below
    if below.sum() < 2 or above.sum() < 2:
        drop_ok = False
    else:
        occupied, empty = float(np.median(i[below])), float(np.median(i[above]))
        noise = np.sqrt(max(abs(occupied), 1e-12)) if _looks_like_counts(i) \
            else _noise(i)
        drop_ok = occupied - empty > max(3 * noise, 0.15 * abs(occupied))
    if not drop_ok:
        level = expected if expected is not None else candidate
        return FermiEstimate(level, False,
                             "no clear Fermi edge where one was expected; "
                             + ("E_F is put at hν − W." if scale == "kinetic"
                                and expected is not None else
                                "check E_F by hand."))

    window = (candidate - half, candidate + half)
    near = (e >= window[0]) & (e <= window[1])
    try:
        if scale == "binding":
            # the fit expects occupied states below E_F: mirror the axis
            fit = fit_fermi_edge(-e[near][::-1], i[near][::-1],
                                 start=initial_guess(-e[near][::-1],
                                                     i[near][::-1], temperature),
                                 temperature=temperature)
            level = -fit.values["ef"]
        else:
            fit = fit_fermi_edge(e[near], i[near],
                                 start=initial_guess(e[near], i[near], temperature),
                                 temperature=temperature)
            level = fit.values["ef"]
        if not (window[0] <= level <= window[1]) or not fit.success:
            raise ValueError("the fit left the window")
        message = (f"Fermi edge fitted at {level:.3f} ± "
                   f"{fit.errors.get('ef', float('nan')):.3f} eV "
                   f"(resolution {fit.values['resolution'] * 1000:.0f} meV).")
        if expected is not None and scale == "kinetic":
            message += (f" hν − W would put it at {expected:.3f} "
                        f"({level - expected:+.3f} eV).")
        return FermiEstimate(level, True, message,
                             error=float(fit.errors.get("ef", float("nan"))),
                             resolution=float(fit.values["resolution"]), fit=fit)
    except Exception:                                    # noqa: BLE001
        return FermiEstimate(candidate, True,
                             f"Fermi edge at {candidate:.3f} eV (steepest "
                             f"drop; the edge fit did not converge).")


@dataclass
class EnergyReference:
    """What turns the axis into binding energies (positive below E_F)."""
    scale: str
    #: where E_F is on the dataset's own axis
    fermi_level: float
    photon_energy: object = None
    work_function: float = DEFAULT_WORK_FUNCTION

    def binding(self, energy):
        energy = np.asarray(energy, dtype=float)
        if self.scale == "binding":
            return energy - self.fermi_level
        return self.fermi_level - energy

    def axis(self, binding):
        """The inverse of :meth:`binding`: where a binding energy is on the
        dataset's own axis."""
        binding = np.asarray(binding, dtype=float)
        if self.scale == "binding":
            return binding + self.fermi_level
        return self.fermi_level - binding

    @property
    def max_binding(self):
        """The deepest level first-order light can reach: hv - W."""
        if self.photon_energy is None:
            return None
        return float(self.photon_energy) - float(self.work_function)

    def implied_work_function(self):
        """hv - E_F_kin, for a kinetic axis whose edge was fitted: a number
        far from 4-5 eV says the photon energy is not what the file says."""
        if self.scale != "kinetic" or self.photon_energy is None:
            return None
        return float(self.photon_energy) - float(self.fermi_level)


# --------------------------------------------------------------------------
# Peaks
# --------------------------------------------------------------------------
def _looks_like_counts(values) -> bool:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0 or values.min() < 0:
        return False
    return bool(np.all(np.abs(values - np.round(values)) < 1e-6))


def _noise(values) -> float:
    """The point-to-point noise of a curve, from the median absolute
    second difference (which a smooth peak or background barely touches)."""
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size < 4:
        return 0.0
    d2 = values[2:] - 2 * values[1:-1] + values[:-2]
    return float(np.median(np.abs(d2 - np.median(d2))) * 1.4826 / np.sqrt(6.0))


@dataclass
class Peak:
    """One peak of the EDC, on the binding-energy scale."""
    binding: float
    height: float
    fwhm: float
    prominence: float
    error: float = float("nan")
    #: how many noise sigmas the peak stands above its surroundings
    significance: float = float("nan")
    #: added by the user rather than found
    manual: bool = False


def _gaussian_line(x, amplitude, centre, sigma, offset, slope):
    return amplitude * np.exp(-0.5 * ((x - centre) / sigma) ** 2) + \
        offset + slope * (x - centre)


def find_peaks(binding, intensity, *, sensitivity: float = SENSITIVITY,
               min_relative: float = MIN_RELATIVE, smoothing: float = None) -> list:
    """The peaks of the EDC, deepest last.

    The curve is put on a uniform grid, smoothed (``smoothing`` eV, by
    default the larger of 0.08 eV and one step) and searched for maxima that
    stand ``sensitivity`` noise sigmas above their surroundings -- the noise
    measured from the curve itself, and scaled as sqrt(counts) for counting
    data so a strong background does not throw up false peaks -- and at
    least ``min_relative`` of the strongest peak. Each is then refined by a
    Gaussian on a sloping background over its own half-width.
    """
    from scipy.ndimage import gaussian_filter1d
    from scipy.optimize import curve_fit
    from scipy.signal import find_peaks as scipy_peaks, peak_widths

    b = np.asarray(binding, dtype=float)
    y = np.asarray(intensity, dtype=float)
    keep = np.isfinite(b) & np.isfinite(y)
    b, y = b[keep], y[keep]
    order = np.argsort(b)
    b, y = b[order], y[order]
    if b.size < 7:
        return []
    step = float(np.median(np.diff(b)))
    if step <= 0:
        return []
    grid = np.arange(b[0], b[-1] + 0.5 * step, step)
    raw = np.interp(grid, b, y)
    counts = _looks_like_counts(y)

    width = max(float(smoothing) if smoothing else 0.08, step)
    sigma_samples = width / step
    smooth = gaussian_filter1d(raw, sigma_samples, mode="nearest") \
        if sigma_samples > 0.3 else raw
    # noise of the smoothed curve: the raw noise over the kernel's length
    reduction = np.sqrt(max(2 * np.sqrt(np.pi) * sigma_samples, 1.0))
    if counts:
        local = np.sqrt(np.clip(smooth, 1.0, None)) / reduction
    else:
        local = np.full(smooth.shape, max(_noise(raw), 1e-12) / reduction)

    index, props = scipy_peaks(smooth, prominence=0.0, width=0.5)
    if index.size == 0:
        return []
    prominence = props["prominences"]
    significance = prominence / np.maximum(local[index], 1e-12)
    strongest = float(prominence.max())
    chosen = (significance >= sensitivity) & \
        (prominence >= min_relative * strongest)
    index, prominence, significance = \
        index[chosen], prominence[chosen], significance[chosen]
    if index.size == 0:
        return []
    widths = peak_widths(smooth, index, rel_height=0.5)[0] * step
    # nothing real is narrower than the smoothing kernel itself: a maximum
    # that is, is noise the smoothing did not quite flatten
    broad = widths >= 0.8 * 2.3548 * max(sigma_samples, 0.5) * step
    index, prominence, significance, widths = \
        index[broad], prominence[broad], significance[broad], widths[broad]

    peaks = []
    for k, p, s, w in zip(index, prominence, significance, widths):
        centre, error = float(grid[k]), float("nan")
        fwhm = float(w)
        half = max(1.2 * fwhm, 3 * step)
        near = (grid >= grid[k] - half) & (grid <= grid[k] + half)
        if near.sum() >= 6:
            x, v = grid[near], raw[near]
            start = (float(p), centre, max(fwhm / 2.355, step),
                     float(np.min(v)), 0.0)
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    popt, pcov = curve_fit(_gaussian_line, x, v, p0=start,
                                           maxfev=4000)
                if abs(popt[1] - centre) < half and popt[0] > 0 \
                        and 0 < abs(popt[2]) < 3 * half:
                    centre = float(popt[1])
                    fwhm = float(abs(popt[2]) * 2.3548)
                    error = float(np.sqrt(abs(pcov[1, 1])))
            except Exception:                            # noqa: BLE001
                pass
        if fwhm < 2 * step:
            continue            # a spike on one or two points: noise
        peaks.append(Peak(binding=centre, height=float(smooth[k]),
                          fwhm=fwhm, prominence=float(p), error=error,
                          significance=float(s)))
    return peaks


def peak_at(binding, intensity, energy: float, *, half: float = 1.0) -> Peak:
    """A peak the user put at ``energy`` (eV binding), marked ``manual``.

    The position is refined by a Gaussian on a sloping background over
    ``half`` eV either side when that fit settles near the click -- a
    shoulder on a strong line usually does -- and kept as given when it
    does not; the FWHM is then NaN. The significance is the Gaussian's
    amplitude over the noise, as for a found peak.
    """
    from scipy.optimize import curve_fit

    b = np.asarray(binding, dtype=float)
    y = np.asarray(intensity, dtype=float)
    keep = np.isfinite(b) & np.isfinite(y)
    b, y = b[keep], y[keep]
    order = np.argsort(b)
    b, y = b[order], y[order]
    energy = float(energy)
    if b.size == 0:
        return Peak(energy, float("nan"), float("nan"), float("nan"),
                    manual=True)
    height = float(np.interp(energy, b, y))
    peak = Peak(energy, height, float("nan"), float("nan"), manual=True)
    near = (b >= energy - half) & (b <= energy + half)
    if near.sum() < 6:
        return peak
    x, v = b[near], y[near]
    step = float(np.median(np.diff(x)))
    base = float(np.min(v))
    start = (max(height - base, 1e-12), energy, max(0.25 * half, step), base, 0.0)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            popt, pcov = curve_fit(_gaussian_line, x, v, p0=start, maxfev=4000)
    except Exception:                                    # noqa: BLE001
        return peak
    amplitude, centre, sigma = float(popt[0]), float(popt[1]), abs(float(popt[2]))
    if amplitude <= 0 or abs(centre - energy) > 0.5 * half or \
            not step <= sigma <= 2 * half:
        return peak
    if _looks_like_counts(y):
        noise = np.sqrt(max(float(np.interp(centre, b, y)), 1.0))
    else:
        noise = max(_noise(y), 1e-12)
    return Peak(centre, float(np.interp(centre, b, y)), sigma * 2.3548,
                amplitude, error=float(np.sqrt(abs(pcov[1, 1]))),
                significance=amplitude / noise, manual=True)


def candidates_near(energy: float, elements=(), *, within: float = HINT_EV,
                    chemical_states: bool = True, limit: int = 8):
    """What a feature at ``energy`` (eV binding) might be: ``(sample,
    others)``, the levels of the sample's ``elements`` and of every other
    element (radioactive ones and H, He aside) within ``within`` eV, each
    nearest first, one entry per element and shell, at most ``limit``
    of the others."""
    energy = float(energy)
    banned = RADIOACTIVE | NO_CORE_LEVEL

    def near(symbols):
        best = {}
        for line in reference_lines(symbols, chemical_states=chemical_states):
            delta = energy - line.energy
            if abs(delta) > within:
                continue
            key = (line.element, line.shell)
            if key not in best or abs(delta) < abs(best[key].delta):
                best[key] = Match(line, float(delta),
                                  note=LIKELY_CONTAMINANTS.get(line.element, ""))
        return sorted(best.values(), key=lambda m: abs(m.delta))

    sample = near([s for s in elements if s in SYMBOLS])
    for match in sample:
        match.note = ""
    others = near([s for z, s in sorted(ELEMENTS.items())
                   if s not in banned and s not in set(elements)])
    return sample, others[:limit]


def describe_candidates(energy: float, elements=(), **kwargs) -> str:
    """:func:`candidates_near` as one line of text."""
    sample, others = candidates_near(energy, elements, **kwargs)
    within = kwargs.get("within", HINT_EV)

    def short(m):
        return f"{m.line.label} {m.line.energy:g} (Δ {m.delta:+.2f})"
    parts = []
    if sample:
        parts.append("sample: " + ", ".join(short(m) for m in sample))
    if others:
        parts.append("other elements: " + ", ".join(short(m) for m in others))
    if not parts:
        return f"no tabulated level within ±{within:g} eV"
    return "; ".join(parts)


# --------------------------------------------------------------------------
# Assignment
# --------------------------------------------------------------------------
@dataclass
class Match:
    """A reference line that may explain a peak."""
    line: Line
    #: peak - line (eV), on the binding scale
    delta: float
    #: the order of the light: 2 and 3 for a level seen through 2hv, 3hv
    order: int = 1
    #: the index of the peak that is this line's spin-orbit partner, if found
    partner: object = None
    note: str = ""

    def describe(self) -> str:
        text = self.line.label
        if self.order > 1:
            text += f" via {self.order}hν (true E_B {self.line.energy:g})"
        text += f", Δ {self.delta:+.2f} eV"
        if self.partner is not None:
            text += ", doublet partner found"
        if self.note:
            text += f" — {self.note}"
        return text


#: What a peak was decided to be.
STATUSES = ("sample", "ambiguous", "uncertain", "unexplained")
STATUS_LABELS = {
    "sample": "sample",
    "ambiguous": "sample (orbital ambiguous)",
    "uncertain": "uncertain — check",
    "unexplained": "not the sample",
}


@dataclass
class Assignment:
    peak: Peak
    status: str
    #: sample lines within the tolerance, best first
    matches: list = field(default_factory=list)
    #: sample lines further away but within a chemical shift
    nearby: list = field(default_factory=list)
    #: lines of common elements that could explain it
    contaminants: list = field(default_factory=list)
    #: levels seen through second- or third-order light
    higher_order: list = field(default_factory=list)
    notes: list = field(default_factory=list)
    #: for a peak added by hand: (sample, others) within HINT_EV of it
    suggestions: tuple = ((), ())

    @property
    def label(self) -> str:
        """The short name a peak is marked with on the plot."""
        if self.matches:
            if self.status != "ambiguous":
                return self.matches[0].line.name
            names = []
            for match in self.matches:
                if match.line.name not in names:
                    names.append(match.line.name)
            return " / ".join(names[:2]) + "?"
        if self.status == "uncertain" and self.nearby:
            return self.nearby[0].line.name + "?"
        if self.contaminants:
            return self.contaminants[0].line.name + "?"
        if self.higher_order:
            return f"{self.higher_order[0].line.name} ({self.higher_order[0].order}hν)?"
        return "?"


@dataclass
class CoreLevelReport:
    """Everything the check found."""
    elements: list
    reference: EnergyReference
    window: tuple
    peaks: list
    assignments: list
    #: (Line, why) for the sample's lines in the window that have no peak
    missing: list = field(default_factory=list)
    #: element -> text, for sample elements with nothing in the window
    out_of_window: dict = field(default_factory=dict)
    #: (shift, peaks matched with it, peaks matched without), or None
    shift: object = None
    notes: list = field(default_factory=list)
    tolerance: float = TOLERANCE_EV

    def text(self) -> str:
        """The report as plain text, for copying out."""
        ref = self.reference
        rows = [f"Core-level check — elements: {', '.join(self.elements)}",
                f"Energy scale: {SCALE_LABELS[ref.scale]}; E_F at "
                f"{ref.fermi_level:.3f} eV on the axis"
                + (f"; hν = {ref.photon_energy:g} eV" if ref.photon_energy
                   else "")
                + f"; W = {ref.work_function:g} eV",
                f"Binding-energy window: {self.window[0]:.2f} to "
                f"{self.window[1]:.2f} eV; tolerance {self.tolerance:g} eV",
                ""]
        for n, a in enumerate(self.assignments, 1):
            p = a.peak
            rows.append(f"Peak {n}: {p.binding:.2f} eV (FWHM {p.fwhm:.2f}, "
                        f"{p.significance:.0f}σ) — {STATUS_LABELS[a.status]}"
                        + (" [added by hand]" if p.manual else ""))
            for title, group in (("matches", a.matches),
                                 ("nearby (chemical shift?)", a.nearby),
                                 ("possible contamination", a.contaminants),
                                 ("higher-order light", a.higher_order),
                                 (f"sample within ±{HINT_EV:g} eV",
                                  a.suggestions[0]),
                                 (f"other elements within ±{HINT_EV:g} eV",
                                  a.suggestions[1])):
                for match in group:
                    rows.append(f"    {title}: {match.describe()}")
            for note in a.notes:
                rows.append(f"    note: {note}")
        if self.missing:
            rows.append("")
            rows.append("Sample lines in the window with no peak:")
            for line, why in self.missing:
                rows.append(f"    {line.label} at {line.energy:g} eV{why}")
        if self.out_of_window:
            rows.append("")
            for symbol, text in self.out_of_window.items():
                rows.append(f"{symbol}: {text}")
        if self.notes:
            rows.append("")
            rows += self.notes
        return "\n".join(rows)


def _matches(peak_be, lines, tolerance, order=1, shift=0.0):
    out = []
    for line in lines:
        delta = peak_be - (line.energy - (order - 1) * shift)
        if abs(delta) <= tolerance:
            out.append(Match(line, float(delta), order=order))
    return out


def _rank(match: Match):
    """Best first: a table value before a chemical-state one at the same
    distance, the major component of a doublet before the minor."""
    return (abs(match.delta) + (0.05 if match.line.source != "table" else 0.0)
            + (0.0 if match.line.is_major() else 0.02))


def _find_partner(match: Match, index: int, peaks, tolerance):
    """The index of the peak where ``match``'s doublet partner should be."""
    partner_orbital = match.line.partner_orbital()
    if partner_orbital is None:
        return None
    if match.line.source == "table":
        partner_energy = BINDING_ENERGIES.get(match.line.element, {}).get(
            partner_orbital)
    else:
        partner_energy = next((v for e, o, v, s in CHEMICAL_STATES
                               if e == match.line.element and o == partner_orbital
                               and s == match.line.source), None)
    if partner_energy is None:
        return None
    expected = peaks[index].binding + (partner_energy - match.line.energy)
    best, best_gap = None, None
    for j, other in enumerate(peaks):
        if j == index:
            continue
        gap = abs(other.binding - expected)
        if gap <= min(tolerance, 0.5) and (best is None or gap < best_gap):
            best, best_gap = j, gap
    return best


def assign(peaks, elements, reference: EnergyReference, *,
           tolerance: float = TOLERANCE_EV, uncertain: float = UNCERTAIN_EV,
           contaminant: float = CONTAMINANT_EV, chemical_states: bool = True,
           window=None) -> CoreLevelReport:
    """Compare every peak with the sample's lines; see the module notes."""
    peaks = sorted(peaks, key=lambda p: p.binding)
    if window is None:
        bs = [p.binding for p in peaks] or [0.0, 0.0]
        window = (min(bs), max(bs))
    sample = reference_lines(elements, chemical_states=chemical_states)
    common = common_elements(exclude=elements)
    common_lines = reference_lines(common, chemical_states=chemical_states)
    priority = {symbol: n for n, symbol in enumerate(common)}
    # what second- and third-order light might show: the sample's levels,
    # and those of the likeliest contaminants
    higher_order_lines = reference_lines(
        list(elements) + [s for s in common if s in LIKELY_CONTAMINANTS],
        chemical_states=False)
    hv = reference.photon_energy
    max_be = reference.max_binding

    assignments = []
    for i, peak in enumerate(peaks):
        be = peak.binding
        matches = sorted(_matches(be, sample, tolerance), key=_rank)
        for match in matches:
            match.partner = _find_partner(match, i, peaks, tolerance)
            if match.line.source != "table":
                match.note = "typical value for this chemical state"
        # a line whose doublet partner is there beats one whose is not
        if any(m.partner is not None for m in matches):
            matches.sort(key=lambda m: (m.partner is None, _rank(m)))
        nearby = []
        if not matches:
            nearby = sorted(_matches(be, sample, uncertain), key=_rank)[:4]
        notes = []

        if matches:
            elements_hit = {m.line.element for m in matches}
            shells_hit = {(m.line.element, m.line.shell) for m in matches}
            status = "sample"
            if len(elements_hit) > 1 or len(shells_hit) > 1:
                best, second = matches[0], matches[1]
                clear = (best.partner is not None and second.partner is None) \
                    or abs(second.delta) - abs(best.delta) > 0.5 * tolerance
                if not clear:
                    status = "ambiguous"
                    notes.append("more than one level of the sample fits; "
                                 "the doublet partner or another photon "
                                 "energy decides")
        elif nearby:
            # further than the tolerance: a chemical shift, charging, or not
            # the sample at all -- the nearby lines are listed, the call is
            # the user's
            status = "uncertain"
        else:
            status = "unexplained"

        contaminants, higher = [], []
        if status in ("uncertain", "unexplained"):
            found = _matches(be, common_lines, contaminant)
            found.sort(key=lambda m: (m.line.element not in LIKELY_CONTAMINANTS,
                                      abs(m.delta) > tolerance,
                                      priority.get(m.line.element, 999),
                                      abs(m.delta)))
            seen = set()
            for match in found:
                key = (match.line.element, match.line.shell)
                if key in seen:
                    continue
                seen.add(key)
                match.partner = _find_partner(match, i, peaks, contaminant)
                source = LIKELY_CONTAMINANTS.get(match.line.element)
                if source:
                    match.note = source
                contaminants.append(match)
                if len(contaminants) >= 5:
                    break
            explained = any(m.line.element in LIKELY_CONTAMINANTS
                            and abs(m.delta) <= tolerance for m in contaminants)
            if hv is not None and not explained:
                for order in (2, 3):
                    higher += _matches(be, higher_order_lines, tolerance,
                                       order=order, shift=hv)
                higher.sort(key=lambda m: (m.line.element not in elements,
                                           m.order, abs(m.delta)))
                higher = higher[:4]
            if hv is not None:
                for element, transition, kinetic in AUGER_LINES:
                    apparent = hv - reference.work_function - kinetic
                    if abs(be - apparent) <= AUGER_EV:
                        notes.append(
                            f"possibly the {element} {transition} Auger line "
                            f"(kinetic ≈ {kinetic:g} eV, here ≈ "
                            f"{apparent:.1f} eV binding); it moves with hν")
        if max_be is not None and be > max_be:
            notes.append(f"deeper than hν − W = {max_be:.1f} eV: first-order "
                         f"light cannot reach it")
        if be < 0:
            notes.append("above E_F — check the energy reference")
        elif be < 3.0:
            notes.append("valence-band region: may be band structure")
        suggestions = candidates_near(be, elements,
                                      chemical_states=chemical_states) \
            if peak.manual else ((), ())
        assignments.append(Assignment(peak, status, matches, nearby,
                                      contaminants, higher, notes, suggestions))

    report = CoreLevelReport(list(elements), reference, tuple(window), peaks,
                             assignments, tolerance=tolerance)
    _missing_lines(report, sample, tolerance)
    _out_of_window(report, elements)
    report.shift = suggest_shift(peaks, elements, tolerance=tolerance,
                                 chemical_states=False)
    if report.shift is not None:
        shift, with_shift, without = report.shift
        report.notes.append(
            f"A rigid shift of {shift:+.2f} eV (subtract it from every binding "
            f"energy) would make {with_shift} peak(s) match the table instead "
            f"of {without}: charging of the sample, or a wrong E_F / hν.")
    if any(e in VACUUM_REFERENCED for e in elements):
        gas = ", ".join(e for e in elements if e in VACUUM_REFERENCED)
        report.notes.append(
            f"The table gives {gas} for the free gas against the vacuum "
            f"level; in a solid they sit several eV lower, so the typical "
            f"solid-state values are compared too.")
    return report


def _missing_lines(report: CoreLevelReport, sample, tolerance):
    lo, hi = report.window
    margin = 0.5
    matched = {(m.line.element, m.line.orbital, m.line.source)
               for a in report.assignments for m in a.matches}
    nearby = {(m.line.element, m.line.orbital, m.line.source)
              for a in report.assignments for m in a.nearby}
    for line in sample:
        if line.source != "table" or not (lo + margin <= line.energy <= hi - margin):
            continue
        key = (line.element, line.orbital, line.source)
        if key in matched:
            continue
        # a chemical-state line of the same level that did match counts
        if any(m.line.element == line.element and m.line.shell == line.shell
               for a in report.assignments for m in a.matches):
            continue
        why = (" (a peak nearby, shifted — see the table)" if key in nearby
               else " (weak at this hν, buried, or absent at the surface)")
        report.missing.append((line, why))


def _out_of_window(report: CoreLevelReport, elements):
    lo, hi = report.window
    max_be = report.reference.max_binding
    for symbol in elements:
        lines = [l for l in element_lines(symbol, chemical_states=False)]
        if any(lo <= l.energy <= hi for l in lines):
            continue
        reachable = [l for l in lines if max_be is None or l.energy < max_be - 5]
        if not reachable:
            report.out_of_window[symbol] = (
                "no level reachable at this photon energy"
                + (f" (hν − W = {max_be:.0f} eV)" if max_be else ""))
            continue
        nearest = sorted(reachable, key=lambda l: min(abs(l.energy - lo),
                                                      abs(l.energy - hi)))[:3]
        report.out_of_window[symbol] = (
            "no level in this window; the nearest reachable: "
            + ", ".join(f"{l.orbital} {l.energy:g} eV" for l in
                        sorted(nearest, key=lambda l: l.energy)))


def suggest_shift(peaks, elements, *, tolerance: float = TOLERANCE_EV,
                  span: float = 5.0, step: float = 0.02,
                  chemical_states: bool = False):
    """``(shift, n_with, n_without)`` when moving every peak by the same
    amount would make clearly more of them fit the sample's table values --
    at least two peaks, and at least one more than without the shift --
    else None. The shift is ``measured - tabulated``: subtract it."""
    lines = reference_lines(elements, chemical_states=chemical_states)
    if len(peaks) < 2 or not lines:
        return None
    energies = np.asarray([line.energy for line in lines])
    centres = np.asarray([p.binding for p in peaks])
    shifts = np.arange(-span, span + 0.5 * step, step)
    # distance from every shifted peak to its nearest line
    gaps = np.abs((centres[None, :] - shifts[:, None])[:, :, None]
                  - energies[None, None, :]).min(axis=2)
    counts = (gaps <= tolerance).sum(axis=1)
    spread = np.where(gaps <= tolerance, gaps, 0.0).sum(axis=1)
    zero = int(np.argmin(np.abs(shifts)))
    best_count = counts.max()
    if best_count < 2 or best_count <= counts[zero]:
        return None
    candidates = np.flatnonzero(counts == best_count)
    best = candidates[np.argmin(spread[candidates] + 1e-3 * np.abs(shifts[candidates]))]
    return float(round(shifts[best], 2)), int(best_count), int(counts[zero])


def check(binding, intensity, elements, reference: EnergyReference, *,
          sensitivity: float = SENSITIVITY, tolerance: float = TOLERANCE_EV,
          uncertain: float = UNCERTAIN_EV, chemical_states: bool = True,
          min_relative: float = MIN_RELATIVE, manual=(), ignored=(),
          merge: float = 0.3) -> CoreLevelReport:
    """Steps 3 and 4 in one call: find the peaks and assign them.

    ``manual`` are binding energies the user put peaks at (see
    :func:`peak_at`); ``ignored`` are found peaks the user removed. A
    found peak within ``merge`` eV of either is dropped -- the user's
    peak replaces it, or it is gone.
    """
    b = np.asarray(binding, dtype=float)
    finite = b[np.isfinite(b)]
    window = (float(finite.min()), float(finite.max())) if finite.size else (0, 0)
    peaks = find_peaks(b, intensity, sensitivity=sensitivity,
                       min_relative=min_relative)
    added = [peak_at(b, intensity, energy) for energy in manual]
    drop = list(ignored) + [p.binding for p in added]
    peaks = [p for p in peaks
             if all(abs(p.binding - e) > merge for e in drop)] + added
    return assign(peaks, elements, reference, tolerance=tolerance,
                  uncertain=uncertain, chemical_states=chemical_states,
                  window=window)
