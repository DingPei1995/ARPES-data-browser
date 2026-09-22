"""
tools/cutops.py
===============
Arithmetic between two cuts: dichroism, asymmetries, and dividing by a
reference.

Three measurements need it, and they are the same four operations with
different care taken around them:

* **Linear dichroism**: ``LH - LV``. The two polarisations come off the
  undulator with different flux, so B has to be scaled to A before the
  difference means anything -- otherwise the result is mostly the flux
  ratio, painted onto the bands.
* **Circular dichroism**: ``(C+ - C-) / (C+ + C-)``. Dimensionless, so it
  can be compared between bands of very different intensity -- and for the
  same reason it explodes wherever ``A + B`` is small, which is why a
  low-intensity mask comes with it.
* **Dividing by a reference**: ``A / R``, typically a polycrystalline gold
  spectrum for the detector's channel-to-channel sensitivity. What is
  divided out is usually a *profile* of the reference -- its angular
  sensitivity, summed over energy -- rather than the reference itself, whose
  own Fermi edge and noise would otherwise be printed into the result.

What this module refuses, rather than quietly doing:

* **Two cuts whose axes are in different units** (degrees against A^-1).
  Resampling one onto the other would line up numbers that do not describe
  the same thing.
* **Two cuts with no overlap at all.**

What it does, and says it did:

* **Two cuts on different grids** -- a reference taken with a different
  energy window, say -- are handled by resampling B onto A's grid over the
  region they share. The fraction of A that is covered, and whether any
  resampling happened at all, are reported, because a result built on an
  interpolation should say so.
* **Counting statistics.** When both cuts are raw counts, the Poisson
  uncertainty of the result is propagated pixel by pixel and summarised --
  the median uncertainty, and how much of the map differs from zero by more
  than twice its own error. A dichroism map is exactly the kind of picture
  in which noise at the band edges looks like signal, and that number is
  what tells them apart.

Nothing here imports Qt.
"""
from __future__ import annotations

import re
import warnings
from dataclasses import dataclass, field

import numpy as np

__all__ = ["OPERATIONS", "NORMALISATIONS", "REFERENCE_SHAPES", "PRESETS",
           "CombineResult", "combine", "metadata_differences",
           "polarisation_warnings", "axis_unit", "looks_like_counts"]

#: The operations, as the dialog names them.
OPERATIONS = {
    "difference": "A − B",
    "asymmetry": "(A − B) / (A + B)",
    "ratio": "A / B",
    "sum": "A + B",
}

#: How B is scaled to A before combining. ``region`` scales so the two have
#: the same total inside a box -- a stretch of background above E_F, or a
#: band known not to be dichroic.
NORMALISATIONS = {
    "none": "None",
    "total": "Same total intensity",
    "region": "Same intensity in a region",
}

#: What of B is divided by, for a ratio.
REFERENCE_SHAPES = {
    "full": "The whole reference",
    "angle_profile": "Its angular profile (summed over energy)",
    "energy_profile": "Its energy profile (summed over angle)",
}

#: Starting points for the three measurements this exists for. Every field
#: stays editable afterwards; a preset is a default, not a mode.
PRESETS = {
    "linear_dichroism": {
        "label": "Linear dichroism (LH − LV)",
        "operation": "difference", "normalise": "total",
        "polarisations": "linear",
    },
    "circular_dichroism": {
        "label": "Circular dichroism (C+ − C−)/(C+ + C−)",
        "operation": "asymmetry", "normalise": "total",
        "polarisations": "circular",
    },
    "reference": {
        "label": "Divide by a reference (e.g. gold)",
        "operation": "ratio", "normalise": "none",
        "reference_shape": "angle_profile",
        "polarisations": None,
    },
    "custom": {
        "label": "Custom",
        "polarisations": None,
    },
}


# --------------------------------------------------------------------------
# Units and grids
# --------------------------------------------------------------------------
def axis_unit(label: str) -> str:
    """The unit an axis label ends with: ``"Angle (deg)"`` -> ``"deg"``.

    Normalised so that the spellings in use here compare equal -- ``A^-1``,
    ``Å⁻¹`` and ``1/A`` are all one unit. An unlabelled axis
    returns ``""``, which compares equal to anything: there is nothing to
    object to.
    """
    if not label:
        return ""
    found = re.findall(r"\(([^()]*)\)\s*$", str(label))
    unit = found[-1] if found else ""
    unit = unit.strip().lower().replace(" ", "")
    for spelling in ("å⁻¹", "a^-1", "1/a", "å^-1", "1/å",
                     "inv.a", "a-1"):
        unit = unit.replace(spelling, "invA")
    for spelling in ("degrees", "degree", "deg", "°"):
        if unit == spelling:
            unit = "deg"
    return unit


def _ascending(axis, values, dim):
    """``axis`` made ascending, with ``values`` flipped along ``dim`` to match."""
    axis = np.asarray(axis, dtype=float)
    if axis.size > 1 and axis[0] > axis[-1]:
        return axis[::-1], np.flip(values, axis=dim)
    return axis, values


def _same_grid(a_axes, b_axes, tolerance: float = 0.05) -> bool:
    """Same size, and no sample more than ``tolerance`` of a step apart."""
    for a, b in zip(a_axes, b_axes):
        a = np.asarray(a, dtype=float)
        b = np.asarray(b, dtype=float)
        if a.shape != b.shape:
            return False
        step = abs(float(a[1] - a[0])) if a.size > 1 else 1.0
        if np.max(np.abs(a - b)) > tolerance * (step or 1.0):
            return False
    return True


def _resample(values, from_axes, to_axes):
    """``values`` on ``from_axes``, bilinearly onto ``to_axes``; NaN outside."""
    from scipy.interpolate import RegularGridInterpolator

    x, values = _ascending(from_axes[0], np.asarray(values, dtype=float), 0)
    y, values = _ascending(from_axes[1], values, 1)
    interpolate = RegularGridInterpolator((x, y), values, method="linear",
                                          bounds_error=False, fill_value=np.nan)
    gx, gy = np.meshgrid(np.asarray(to_axes[0], float),
                         np.asarray(to_axes[1], float), indexing="ij")
    return interpolate(np.column_stack([gx.ravel(), gy.ravel()])).reshape(gx.shape)


def looks_like_counts(values) -> bool:
    """Whether these are raw detector counts, for which Poisson statistics
    hold.

    Non-negative and whole-numbered. Anything that has been normalised,
    smoothed or interpolated fails, which is the point: propagating Poisson
    errors through a smoothed map would understate them.
    """
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0 or finite.max() <= 0:
        return False
    if finite.min() < 0:
        return False
    sample = finite if finite.size <= 200_000 else finite[:: finite.size // 200_000]
    return bool(np.mean(np.abs(sample - np.round(sample)) < 1e-6) > 0.99)


# --------------------------------------------------------------------------
# The operation
# --------------------------------------------------------------------------
@dataclass
class CombineResult:
    """The combined cut, on A's grid, and everything worth reporting."""

    values: np.ndarray
    #: Poisson standard deviation per pixel, or None when the inputs are not
    #: raw counts (or the operation is a ratio, where the reference's own
    #: noise is usually negligible and its statistics unknown).
    sigma: "np.ndarray | None"
    #: The factor B was multiplied by before combining.
    scale: float
    #: Whether B had to be interpolated onto A's grid.
    resampled: bool
    #: Fraction of A's grid that B covers.
    overlap: float
    #: B as it entered the arithmetic: on A's grid and times ``scale``.
    #: Kept so a display of the inputs shows exactly what was combined.
    b_used: np.ndarray
    #: Fraction of the covered grid hidden by the low-intensity mask or the
    #: reference floor.
    masked: float
    #: Median of ``sigma`` over what is shown, or None.
    median_sigma: "float | None" = None
    #: Fraction of the shown pixels where |result| > 2 sigma, or None.
    significant: "float | None" = None
    notes: list = field(default_factory=list)

    def summary(self) -> str:
        lines = []
        if self.scale != 1.0:
            lines.append(f"B scaled by {self.scale:.5g} to match A")
        if self.resampled:
            lines.append(f"B resampled onto A's grid; it covers "
                         f"{100 * self.overlap:.1f}% of A")
        if self.masked:
            lines.append(f"{100 * self.masked:.1f}% of the covered area masked")
        if self.median_sigma is not None:
            lines.append(f"median uncertainty {self.median_sigma:.3g} "
                         f"(Poisson, from the raw counts)")
        if self.significant is not None:
            lines.append(f"{100 * self.significant:.1f}% of the shown pixels "
                         f"differ from zero by more than 2σ")
        return "\n".join(lines + list(self.notes))


def _box_mask(axes, region):
    """Boolean mask of the grid points inside ``(x0, y0, x1, y1)``."""
    x0, y0, x1, y1 = (float(v) for v in region)
    gx, gy = np.meshgrid(np.asarray(axes[0], float), np.asarray(axes[1], float),
                         indexing="ij")
    return ((gx >= min(x0, x1)) & (gx <= max(x0, x1))
            & (gy >= min(y0, y1)) & (gy <= max(y0, y1)))


def _reference(b, shape: str):
    """The divisor for a ratio, normalised to a mean of one.

    Normalised so that dividing by it removes a *sensitivity* and leaves A's
    own scale alone: a gold reference with ten times A's counts should not
    make the result ten times smaller.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        if shape == "angle_profile":
            profile = np.nanmean(b, axis=1)[:, None] * np.ones((1, b.shape[1]))
        elif shape == "energy_profile":
            profile = np.nanmean(b, axis=0)[None, :] * np.ones((b.shape[0], 1))
        else:
            profile = np.array(b, dtype=float)
        mean = np.nanmean(profile)
    if not np.isfinite(mean) or mean == 0:
        raise ValueError("the reference is empty or zero everywhere, so there "
                         "is nothing to divide by")
    return profile / mean


def combine(a_values, a_axes, b_values, b_axes, *, operation: str,
            normalise: str = "none", region=None,
            reference_shape: str = "full", reference_floor: float = 0.05,
            intensity_floor: float = 0.0, a_labels=("", ""),
            b_labels=("", "")) -> CombineResult:
    """Combine two cuts on A's grid.

    ``a_values`` and ``b_values`` are ``(x, y)`` -- angle or momentum, then
    energy. ``region`` is ``(x0, y0, x1, y1)`` in axis units, used by
    ``normalise="region"``. ``reference_floor`` hides a ratio wherever the
    reference is below that fraction of its maximum; ``intensity_floor``
    hides a difference, asymmetry or sum wherever ``A + B`` is below that
    fraction of its maximum.
    """
    if operation not in OPERATIONS:
        raise ValueError(f"operation must be one of {sorted(OPERATIONS)}")
    if normalise not in NORMALISATIONS:
        raise ValueError(f"normalise must be one of {sorted(NORMALISATIONS)}")
    if reference_shape not in REFERENCE_SHAPES:
        raise ValueError(f"reference_shape must be one of {sorted(REFERENCE_SHAPES)}")

    for axis, (la, lb) in enumerate(zip(a_labels, b_labels)):
        ua, ub = axis_unit(la), axis_unit(lb)
        if ua and ub and ua != ub:
            raise ValueError(
                f"A's {'first' if axis == 0 else 'second'} axis is in {ua} and "
                f"B's is in {ub}. They do not describe the same thing, so "
                f"lining them up would be meaningless -- convert one of them "
                f"first.")

    a = np.asarray(a_values, dtype=float)
    b_raw = np.asarray(b_values, dtype=float)
    a_axes = [np.asarray(v, dtype=float) for v in a_axes]
    b_axes = [np.asarray(v, dtype=float) for v in b_axes]
    if a.shape != (a_axes[0].size, a_axes[1].size):
        raise ValueError(f"A is {a.shape} but its axes are "
                         f"({a_axes[0].size}, {a_axes[1].size})")
    if b_raw.shape != (b_axes[0].size, b_axes[1].size):
        raise ValueError(f"B is {b_raw.shape} but its axes are "
                         f"({b_axes[0].size}, {b_axes[1].size})")

    counts = looks_like_counts(a) and looks_like_counts(b_raw)
    if _same_grid(a_axes, b_axes):
        b = b_raw
        resampled = False
    else:
        b = _resample(b_raw, b_axes, a_axes)
        resampled = True
        counts = False     # interpolated counts are not Poisson any more
    covered = np.isfinite(b)
    overlap = float(covered.mean())
    if overlap == 0:
        raise ValueError(
            "A and B do not overlap at all: B's axes fall entirely outside A's.")

    notes = []
    if resampled and overlap < 0.5:
        notes.append(f"B covers only {100 * overlap:.0f}% of A; the rest is blank")

    # -- scale B to A --------------------------------------------------------
    scale = 1.0
    if operation != "ratio" and normalise != "none":
        usable = np.isfinite(a) & covered
        if normalise == "region":
            if region is None:
                raise ValueError("normalising to a region needs a region")
            usable &= _box_mask(a_axes, region)
            if not usable.any():
                raise ValueError("the normalisation region contains no data "
                                 "that A and B share")
        total_a = float(np.nansum(a[usable]))
        total_b = float(np.nansum(b[usable]))
        if total_b == 0 or not np.isfinite(total_b):
            raise ValueError("B has no intensity where it is being normalised, "
                             "so it cannot be scaled to A")
        scale = total_a / total_b
    b_scaled = b * scale

    # -- combine --------------------------------------------------------------
    masked = np.zeros(a.shape, dtype=bool)
    sigma = None
    with np.errstate(invalid="ignore", divide="ignore"):
        if operation == "ratio":
            divisor = _reference(b, reference_shape)
            ceiling = np.nanmax(divisor)
            masked = ~(divisor >= reference_floor * ceiling)
            out = np.where(masked, np.nan, a / divisor)
        else:
            total = a + b_scaled
            if intensity_floor > 0:
                ceiling = np.nanmax(total)
                masked = covered & ~(total >= intensity_floor * ceiling)
            if operation == "difference":
                out = a - b_scaled
            elif operation == "sum":
                out = total
            else:
                out = np.where(total > 0, (a - b_scaled) / total, np.nan)
            out = np.where(masked, np.nan, out)
            if counts:
                if operation in ("difference", "sum"):
                    sigma = np.sqrt(a + scale ** 2 * b)
                else:
                    # d/dA and d/dB' of (A - B')/(A + B'), B' = sB, with
                    # var A = A and var B' = s^2 B for raw counts.
                    sigma = np.sqrt(4.0 * scale ** 2 * a * b * (a + b)
                                    / np.where(total > 0, total, np.nan) ** 4)
                sigma = np.where(np.isfinite(out), sigma, np.nan)
    out = np.where(covered, out, np.nan)

    shown = np.isfinite(out)
    masked_fraction = (float((masked & covered).sum()) / float(covered.sum())
                       if covered.any() else 0.0)
    median_sigma = significant = None
    if sigma is not None and shown.any():
        good = shown & np.isfinite(sigma) & (sigma > 0)
        if good.any():
            median_sigma = float(np.median(sigma[good]))
            significant = float(np.mean(np.abs(out[good]) > 2.0 * sigma[good]))
            # First-order propagation, measured against Monte Carlo: for the
            # asymmetry it is within 5% from ~20 counts per pixel up, but
            # 12% low at 10 counts and 32% low at 4 -- the estimate is built
            # from the noisy counts themselves, and the ratio stops being
            # Gaussian. Thin data therefore overstate the 2-sigma fraction.
            typical = float(np.median((a + b_scaled)[good]))
            if operation == "asymmetry" and typical < 20:
                notes.append(
                    f"only ~{typical:.0f} counts per pixel: the propagated "
                    f"uncertainty is an underestimate here (by ~12% at 10 "
                    f"counts, ~32% at 4), so the 2σ fraction is too "
                    f"high. Bin the cuts first (Data operations → "
                    f"compress) for a trustworthy number.")
    return CombineResult(values=out, sigma=sigma, scale=float(scale),
                         resampled=resampled, overlap=overlap, b_used=b_scaled,
                         masked=masked_fraction, median_sigma=median_sigma,
                         significant=significant, notes=notes)


# --------------------------------------------------------------------------
# Are the two measurements comparable?
# --------------------------------------------------------------------------
#: What should be the same between two cuts that are combined, where to find
#: it in each loader's metadata, and how far apart counts as different.
#: ``None`` as a tolerance means "compare as text". Keys are tried in order;
#: patterns (strings starting with ``~``) are matched against every key.
#:
#: The positions matter more than they look: on a nano-ARPES line a 5 um
#: move is a different spot on the sample, and a dichroism map between two
#: spots is a map of the sample's inhomogeneity.
COMPARED = [
    ("Polarisation", ("cassiopee.polarisation", "~polari[sz]ation$"), None),
    ("Photon energy (eV)", ("photon_energy_eV", "cassiopee.photon_energy_eV",
                            "~^mono\\.energy$"), 0.01),
    ("Theta (deg)", ("cassiopee.sample_theta_deg", "~(^|\\.)theta$"), 0.05),
    ("Tilt (deg)", ("cassiopee.sample_tilt_deg", "~(^|\\.)tilt$"), 0.05),
    ("Phi (deg)", ("cassiopee.sample_phi_deg", "~(^|\\.)phi$"), 0.05),
    ("X (mm)", ("cassiopee.sample_X_mm",), 0.005),
    ("Y (mm)", ("cassiopee.sample_Y_mm",), 0.005),
    ("Z (mm)", ("cassiopee.sample_Z_mm",), 0.005),
    ("Temperature (K)", ("cassiopee.temperature_K", "SampleTemperature_K"), 1.0),
    ("Pass energy (eV)", ("pass_energy_eV",), 1e-6),
    ("Lens mode", ("lens_mode",), None),
]


def _lookup(info: dict, keys):
    for key in keys:
        if key.startswith("~"):
            pattern = re.compile(key[1:], re.IGNORECASE)
            for name in sorted(info):
                if pattern.search(name):
                    return info[name]
        elif key in info:
            return info[key]
    return None


def metadata_differences(info_a: dict, info_b: dict):
    """What the two measurements recorded, side by side.

    Returns ``[(quantity, a, b, differs), ...]`` for every quantity at least
    one of them recorded. ``differs`` is True, False, or None when only one
    side has a value.
    """
    rows = []
    for name, keys, tolerance in COMPARED:
        a = _lookup(info_a or {}, keys)
        b = _lookup(info_b or {}, keys)
        if a is None and b is None:
            continue
        if a is None or b is None:
            rows.append((name, a, b, None))
            continue
        if tolerance is None:
            differs = str(a).strip().upper() != str(b).strip().upper()
        else:
            try:
                differs = abs(float(a) - float(b)) > tolerance
            except (TypeError, ValueError):
                differs = str(a).strip() != str(b).strip()
        rows.append((name, a, b, differs))
    return rows


def _kind_of_polarisation(value) -> "str | None":
    if value is None:
        return None
    text = str(value).strip().upper()
    # Circular first: "LCP" starts with an L and is not linear.
    if text in ("CR", "CL", "RCP", "LCP", "C+", "C-") or text.startswith("C"):
        return "circular"
    if text in ("LH", "LV", "AH", "AV") or text.startswith("L"):
        return "linear"
    return None


def polarisation_warnings(info_a: dict, info_b: dict, preset: str):
    """What the recorded polarisations say about the chosen measurement.

    Only checks what can be checked: a preset with no expectation, or a
    measurement with no recorded polarisation, returns nothing.
    """
    expected = PRESETS.get(preset, {}).get("polarisations")
    pa = _lookup(info_a or {}, COMPARED[0][1])
    pb = _lookup(info_b or {}, COMPARED[0][1])
    if not expected or pa is None or pb is None:
        return []
    out = []
    ka, kb = _kind_of_polarisation(pa), _kind_of_polarisation(pb)
    if str(pa).upper() == str(pb).upper():
        out.append(f"A and B were both measured in {pa}: this is not a "
                   f"dichroism measurement, and the result will be noise.")
        return out
    if expected == "linear":
        if ka == "circular" or kb == "circular":
            out.append(f"A is {pa} and B is {pb}; linear dichroism wants two "
                       f"linear polarisations.")
        elif str(pa).upper() == "LV" and str(pb).upper() == "LH":
            out.append("A is LV and B is LH, so this is LV − LH. Swap them "
                       "for the usual LH − LV sign.")
    if expected == "circular" and (ka == "linear" or kb == "linear"):
        out.append(f"A is {pa} and B is {pb}; circular dichroism wants two "
                   f"circular polarisations.")
    return out
