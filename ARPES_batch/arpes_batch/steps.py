"""
The processing steps a recipe can name, each a thin call into ``tools/``.

Every step takes the dataset and the job's context and returns the next
dataset and a few numbers worth checking (``qc``). The parameters and the
record they leave (prefix, step name) are the ones the matching dialog
uses, so a batch result and a hand-made one say the same thing about how
they were made.

    step           viewer equivalent                         tools/ function
    ------------   ---------------------------------------   ------------------------------
    degrid         Functions -> De-grid map...               tools.degrid.degrid_map
    calibrate_energy  Fermi surface from a reference + FS
                   correction + Fermi level "offset"         tools.fermi / tools.analysis.fs_correction
    crop           Data operations -> Truncate               tools.dataops.truncate
    bin            Data operations -> Compress               tools.dataops.compress
    normalise      Data operations -> Self-normalise         tools.dataops.self_normalize
    kconvert       Functions -> Map k conversion             tools.kspace.convert_map

Adding a step is a function here and one line in :data:`STEPS`.
"""
from __future__ import annotations

import numpy as np

from loader import registry
from tools import dataops
from tools.analysis import fs_correction as _fs_correction

from .center import suggest_centre
from .dataset import Dataset

ANGSTROM = "Å⁻¹"


class Skip(Exception):
    """The step does not apply to this dataset; the chain goes on without it."""


def _need(ds: Dataset, kinds, step: str):
    if ds.kind not in kinds:
        raise Skip(f"{step} applies to {', '.join(kinds)}, not a {ds.kind}")


def _energy_slot(ds):
    return {"map": "z", "cut": "y"}[ds.kind]


def _slit_slot(ds):
    return {"map": "k", "cut": "x"}[ds.kind]


# -- de-grid -----------------------------------------------------------------
def degrid(ds: Dataset, ctx: dict, **settings):
    from tools import degrid as DG
    _need(ds, ("map",), "degrid")
    why = DG.not_pixel_locked(ds.info)
    if why:
        raise Skip(f"cannot de-grid: {why}")
    try:
        result = DG.degrid_map(ds.array(), settings=settings or None)
    except DG.GridNotFound as exc:
        raise Skip(f"no detector grid found ({exc})")
    m = result.model
    params = {"method": result.method, **settings,
              "grid_peaks": len(m.peaks), "grid_kspace_fraction": float(m.region.mean()),
              "grid_rms": float(m.rms), "contrast_before": float(result.contrast_before),
              "contrast_after": float(result.contrast_after)}
    qc = {"grid_rms_percent": 100 * float(m.rms),
          "grid_power_before": float(result.contrast_before),
          "grid_power_after": float(result.contrast_after)}
    return ds.derive(step="degrid", params=params, value=result.values), qc


# -- energy: E_F and the slit curvature ----------------------------------------
def calibrate_energy(ds: Dataset, ctx: dict, *, source: str = "reference",
                     fs_correction: bool = True, order: int = 2,
                     window_eV: float = 0.25):
    """Put E_F at zero and, optionally, flatten the edge along the slit.

    ``source``: ``"reference"`` -- the curvature and E_F of the reference the
    runner chose for this map (moved by the photon-energy difference);
    ``"self"`` -- the map's own edge, summed over the deflector (metals
    only: a semiconductor has no edge to fit).
    """
    from .reference import fit_reference_on
    _need(ds, ("map", "cut"), "calibrate_energy")
    e_slot, s_slot = _energy_slot(ds), _slit_slot(ds)
    energy, slit = ds.axis(e_slot), ds.axis(s_slot)
    hv = (ctx.get("row") or {}).get("photon_energy_eV")
    notes = []

    if source == "reference":
        ref = ctx.get("reference")
        if ref is None:
            raise Skip("no matching reference: " + "; ".join(ctx.get("reference_notes") or []))
        coeffs = np.array(ref.coeffs, dtype=float)
        shift = (float(hv) - float(ref.photon_energy_eV)
                 if hv is not None and ref.photon_energy_eV is not None else 0.0)
        coeffs[-1] += shift
        origin = {"reference": ref.path, "reference_entry": ref.entry,
                  "reference_hv_eV": ref.photon_energy_eV, "hv_shift_eV": shift,
                  "reference_residual_meV": ref.residual_meV}
        lo, hi = ref.slit_range
        if slit.min() < lo - 0.5 or slit.max() > hi + 0.5:
            notes.append(f"slit axis {slit.min():.3g}..{slit.max():.3g} runs beyond the "
                         f"reference's {lo:.3g}..{hi:.3g}; the curve is extrapolated there")
    elif source == "self":
        ref = fit_reference_on(ds, order=order, window_eV=window_eV)
        coeffs = np.array(ref.coeffs, dtype=float)
        origin = {"reference": "self", "reference_residual_meV": ref.residual_meV}
    else:
        raise ValueError(f"calibrate_energy: unknown source {source!r}")

    fitted = np.polyval(coeffs, slit)
    if not (energy.min() < np.median(fitted) < energy.max()):
        raise Skip(f"E_F from the reference ({np.median(fitted):.4g} eV) lies outside "
                   f"this dataset's energy window {energy.min():.4g}..{energy.max():.4g} "
                   f"-- wrong reference, or the photon energy is not what the file says")
    out = ds
    if fs_correction:
        angle_dim = ds.slots.index(s_slot)
        energy_dim = ds.slots.index(e_slot)
        corrected, new_energy = _fs_correction(ds.array(float), slit, energy, coeffs,
                                              angle_dim=angle_dim, energy_dim=energy_dim)
        axes = list(ds.axes)
        axes[energy_dim] = new_energy
        out = ds.derive(step="fs_correction", prefix="fscorr", value=corrected, axes=axes,
                        params={"order": len(coeffs) - 1, "coefficients": list(coeffs),
                                "angle_axis": ds.label(s_slot), **origin},
                        suffix="FS corr")
        ef = float(np.nanmax(fitted))       # fs_correction flattens onto the top
        energy = new_energy
    else:
        ef = float(np.polyval(coeffs, np.mean([slit.min(), slit.max()])))

    axes = list(out.axes)
    e_dim = out.slots.index(e_slot)
    axes[e_dim] = np.asarray(energy, dtype=float) - ef
    labels = dict(out.labels)
    labels[e_slot] = "E - E_F (eV)"
    out = out.derive(step="fermi_offset", prefix="fitEF", axes=axes, labels=labels,
                     params={"ef_subtracted": ef, **origin}, suffix="E-Ef",
                     extra_info={"batch.ef_kinetic_eV": ef})
    qc = {"ef_kinetic_eV": ef, "edge_curvature_span_meV": 1000 * float(np.ptp(fitted)),
          **({"notes": notes} if notes else {})}
    return out, qc


# -- data operations ------------------------------------------------------------
_AXIS_NAMES = {"map": {"deflector": "x", "slit": "k", "energy": "z"},
               "k_map": {"kx": "x", "ky": "k", "energy": "z"},
               "cut": {"slit": "x", "energy": "y"}}


def crop(ds: Dataset, ctx: dict, **bounds):
    """``energy=[lo, hi]``, ``deflector=[..]``, ``slit=[..]`` in axis units
    (``None`` for an open end)."""
    names = _AXIS_NAMES.get(ds.kind, {})
    unknown = set(bounds) - set(names)
    if unknown:
        raise ValueError(f"crop: {sorted(unknown)} are not axes of a {ds.kind} "
                         f"(it has {sorted(names)})")
    per_slot = {names[k]: tuple(v) if v is not None else (None, None) for k, v in bounds.items()}
    for key, (lo, hi) in ((k, per_slot[names[k]]) for k in bounds):
        axis = ds.axis(names[key])
        a, b = float(np.nanmin(axis)), float(np.nanmax(axis))
        if (lo is not None and lo > b) or (hi is not None and hi < a):
            hint = (" -- these look relative to E_F, but the energy axis is still kinetic "
                    "(calibrate_energy did not run)" if key == "energy" and "E_F" not in
                    ds.label(names[key]) else "")
            raise Skip(f"crop {key}={[lo, hi]} lies outside the axis {a:.4g}..{b:.4g}{hint}")
    limits = [per_slot.get(slot, (None, None)) for slot in ds.slots]
    values, axes = dataops.truncate(ds.array(), [ds.axis(s) for s in ds.slots], limits)
    return ds.derive(step="truncate", params={k: list(v) for k, v in bounds.items()},
                     value=values, axes=axes), {"shape": list(values.shape)}


def bin_(ds: Dataset, ctx: dict, *, factors):
    values, axes = dataops.compress(ds.array(), [ds.axis(s) for s in ds.slots], factors)
    return ds.derive(step="compress", params={"factors": list(factors)},
                     value=values, axes=axes), {"shape": list(values.shape)}


def normalise(ds: Dataset, ctx: dict, *, dims=(1, 2), window=None, to_peak: bool = False):
    """Default: every deflector slice divided by its own total -- takes the
    beam-current decay out of a map."""
    values = dataops.self_normalize(ds.array(), [ds.axis(s) for s in ds.slots],
                                    tuple(dims), window=tuple(window) if window else None,
                                    to_peak=to_peak)
    return ds.derive(step="self_normalise", params={"dims": list(dims), "to_peak": to_peak},
                     value=values), {}


# -- momentum ----------------------------------------------------------------------
def kconvert(ds: Dataset, ctx: dict, *, theta_offset_deg=0.0, phi_offset_deg=0.0,
             azimuth_deg: float = 0.0, n_kx: int = 200, n_ky: int = 200,
             n_energy: int = None, centre_energy: float = None,
             centre_width: float = 0.05, min_centre_score: float = 0.0):
    """Map k conversion. ``theta_offset_deg`` / ``phi_offset_deg`` are the
    deflector and slit angles that become k = 0; ``"auto"`` takes them from
    :func:`center.suggest_centre` (recorded, with its score)."""
    from tools.kspace import convert_map
    _need(ds, ("map",), "kconvert")
    scan_like = type("S", (), {"info": ds.info})()
    if not registry.role_is_angle(scan_like):
        raise Skip(f"first axis is a {ds.info.get('axis0.role')}, not an angle")
    row = ctx.get("row") or {}
    if row.get("slit_axis_looks_like") == "momentum":
        raise Skip("the slit axis looks like it is already in 1/A; converting it "
                   "as an angle would be wrong -- check the file first")

    qc = {}
    if "auto" in (theta_offset_deg, phi_offset_deg):
        found = suggest_centre(ds, energy=centre_energy, width=centre_width)
        qc["centre_suggestion"] = found
        if found["score"] < float(min_centre_score):
            raise Skip(f"Gamma not found with confidence (score {found['score']:.2f} < "
                       f"{min_centre_score}); give the offsets explicitly")
        if theta_offset_deg == "auto":
            theta_offset_deg = found["theta_offset_deg"]
        if phi_offset_deg == "auto":
            phi_offset_deg = found["phi_offset_deg"]

    ef_kin = ds.info.get("batch.ef_kinetic_eV")
    energy = ds.axis("z")
    if ef_kin is None and "E_F" in ds.label("z"):
        raise Skip("the energy axis is relative to E_F but the kinetic E_F is not "
                   "recorded; run calibrate_energy in the same recipe")
    offset = float(ef_kin or 0.0)
    settings = {"theta_offset_deg": float(theta_offset_deg),
                "phi_offset_deg": float(phi_offset_deg),
                "azimuth_deg": float(azimuth_deg), "energy_offset_eV": offset,
                "n_kx": int(n_kx), "n_ky": int(n_ky), "n_energy": n_energy}
    kx, ky, e_out, cube = convert_map(ds.axis("x"), ds.axis("k"), energy, ds.array(float),
                                      **settings)
    e_out = e_out - offset          # back onto the axis the map had
    labels = {"x": f"kx ({ANGSTROM})", "k": f"ky ({ANGSTROM})", "z": ds.label("z")}
    qc["inside_light_cone_percent"] = float(np.isfinite(cube).mean() * 100)
    out = ds.derive(step="k_conversion", prefix="kconv", kind="k_map",
                    value=cube.astype(np.float32), axes=(kx, ky, e_out), labels=labels,
                    params=settings, suffix="k")
    return out, qc


STEPS = {
    "degrid": degrid,
    "calibrate_energy": calibrate_energy,
    "crop": crop,
    "bin": bin_,
    "normalise": normalise,
    "kconvert": kconvert,
}
