"""
The processing steps a recipe can name, each a thin call into ``tools/``.

Every step takes the dataset and the job's context and returns the next
dataset and a few numbers worth checking (``qc``). The parameters and the
record they leave (prefix, step name) are the ones the matching dialog
uses, so a batch result and a hand-made one say the same thing about how
they were made.

    step           viewer equivalent                         tools/ function
    ------------   ---------------------------------------   ------------------------------
    degrid         Functions -> De-grid (map/kz map/cut)     tools.degrid.degrid_map / degrid_cut_with_grid
    calibrate_energy  Fermi surface from a reference + FS
                   correction + Fermi level "offset"         tools.fermi / tools.analysis.fs_correction
    crop           Data operations -> Truncate               tools.dataops.truncate
    bin            Data operations -> Compress               tools.dataops.compress
    normalise      Data operations -> Self-normalise         tools.dataops.self_normalize
    kconvert       Functions -> Map k conversion             tools.kspace.convert_map
    kz_calibrate   kz map -> kz map processing               tools.kzmap (fit_levels, align, normalise_totals)
    kz_match_calculation  (new) scan vs a band calculation   arpes_batch.calcbands
    kz_convert     kz map -> kz -> momentum                  tools.kzconv (to_kz_cube, scan_inner_potential)
    kconvert_cut   Cut -> Cut k conversion                   tools.cutk.convert_cut
    curvature      Process -> Curvature                      tools.process.curvature

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
    return {"map": "z", "kz_map": "z", "cut": "y"}[ds.kind]


def _slit_slot(ds):
    return {"map": "k", "kz_map": "k", "cut": "x"}[ds.kind]


def _relative(ds) -> bool:
    """Whether the energy axis is already E - E_F."""
    return "E_F" in ds.label(_energy_slot(ds)) or "E-E_F" in ds.label(_energy_slot(ds))


# -- de-grid -----------------------------------------------------------------
def _grid_key(ds: Dataset, frame_shape) -> str:
    """Which detector settings a grid pattern belongs to: the grid is the
    same for the same lens mode, pass energy and detector window."""
    import re
    from .inventory import META_KEYS, meta
    parts = [str(meta(ds.info, META_KEYS["lens_mode"]) or "lens?"),
             str(meta(ds.info, META_KEYS["pass_energy"]) or "PE?"),
             "x".join(str(n) for n in frame_shape)]
    return re.sub(r"[^A-Za-z0-9._-]+", "_", "__".join(parts))


def _store_grid(ctx, ds, result):
    """Keep a map's grid for the cuts taken with the same settings (the
    viewer's "Also list the grid pattern" as a file)."""
    import os
    folder = ctx.get("grid_dir")
    if not folder:
        return None
    os.makedirs(folder, exist_ok=True)
    pattern = result.model.full()
    path = os.path.join(folder, _grid_key(ds, pattern.shape) + ".npz")
    tmp = path + ".partial.npz"
    np.savez_compressed(tmp, pattern=pattern, source=ds.name)
    os.replace(tmp, path)
    return path


def degrid(ds: Dataset, ctx: dict, *, cut_fallback: str = "skip", **settings):
    """De-grid a map or kz map with its own grid (and keep that grid), or a
    cut with the grid of a map taken on the same detector settings.

    ``cut_fallback`` is what a cut does with no such grid: ``"skip"``
    (default) leaves it as it is; ``"notch"`` notches its own grid peaks,
    which also removes whatever photoemission lies in those k-space regions
    -- so it is a choice, not a default.
    """
    import os
    from tools import degrid as DG
    _need(ds, ("map", "kz_map", "cut"), "degrid")
    why = DG.not_pixel_locked(ds.info)
    if why:
        raise Skip(f"cannot de-grid: {why}")
    qc = {}
    try:
        if ds.kind == "cut":
            frame = ds.array(float)
            grid_file = os.path.join(ctx.get("grid_dir") or "", _grid_key(ds, frame.shape) + ".npz")
            if ctx.get("grid_dir") and os.path.exists(grid_file):
                with np.load(grid_file, allow_pickle=False) as stored:
                    result = DG.degrid_cut_with_grid(frame, stored["pattern"], settings=settings or None)
                    qc["grid_from"] = str(stored["source"])
            elif cut_fallback == "notch":
                result = DG.degrid_cut_notch(frame, settings=settings or None)
            else:
                raise Skip("no map grid with this cut's lens mode, pass energy and "
                           "detector window (de-grid such a map in the same run, or "
                           "set cut_fallback: notch)")
        else:
            result = DG.degrid_map(ds.array(), settings=settings or None)
            qc["grid_saved"] = _store_grid(ctx, ds, result)
    except DG.GridNotFound as exc:
        raise Skip(f"no detector grid found ({exc})")
    m = result.model
    params = {"method": result.method, **settings,
              "grid_peaks": len(m.peaks), "grid_kspace_fraction": float(m.region.mean()),
              "grid_rms": float(m.rms), "contrast_before": float(result.contrast_before),
              "contrast_after": float(result.contrast_after)}
    if "grid_from" in qc:
        params["grid_source"] = qc["grid_from"]
    qc.update({"method": result.method, "grid_rms_percent": 100 * float(m.rms),
               "grid_power_before": float(result.contrast_before),
               "grid_power_after": float(result.contrast_after)})
    if result.notes:
        qc["notes"] = list(result.notes)
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
               "kz_map": {"photon_energy": "x", "slit": "k", "energy": "z"},
               "kz_map_k": {"kz": "x", "kpar": "k", "energy": "z"},
               "cut": {"slit": "x", "k": "x", "energy": "y"}}


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


# -- photon-energy scans ------------------------------------------------------------
def _temperature(ds, given=None) -> float:
    from .inventory import META_KEYS, _number, meta
    return float(given or _number(meta(ds.info, META_KEYS["temperature_K"])) or 30.0)


def kz_calibrate(ds: Dataset, ctx: dict, *, source: str = "self",
                 fs_correction: bool = True, edge_window_eV: float = 0.3,
                 angle_fraction: float = 0.8, temperature: float = None,
                 normalise: bool = True):
    """Put every spectrum of a photon-energy scan on one Fermi level
    (the viewer's "kz map processing", ``tools.kzmap``).

    ``source``:
      ``"self"``      each spectrum's own edge, fitted over the same detector
                      channels and energy window in every spectrum (metals);
      ``"reference"`` E_F,kin(hv) = hv - W with W from the gold reference --
                      for a semiconductor, which has no edge; needs a kinetic
                      energy axis;
      ``"file"``      the loader already referred the axis to E_F; only the
                      curvature and the normalisation are applied.

    ``fs_correction`` flattens the edge along the slit with the reference's
    curvature first -- an analyser property, the same at every photon energy.
    """
    from tools import kzconv, kzmap
    from tools.fermi import initial_guess
    _need(ds, ("kz_map",), "kz_calibrate")
    hv, slit, energy = ds.axis("x"), ds.axis("k"), ds.axis("z")
    cube = ds.array(float)
    relative = _relative(ds)
    ref = ctx.get("reference")
    notes, params = [], {"source": source}
    T = _temperature(ds, temperature)
    out = ds

    if fs_correction:
        if ref is None:
            notes.append("no reference with these analyser settings: the edge's "
                         "curvature along the slit is not corrected")
        else:
            cube, energy = _fs_correction(cube, slit, energy, ref.coeffs,
                                          angle_dim=1, energy_dim=2)
            out = ds.derive(step="fs_correction", prefix="fscorr", value=cube,
                            axes=(hv, slit, energy), suffix="FS corr",
                            params={"order": ref.order, "coefficients": list(ref.coeffs),
                                    "angle_axis": ds.label("k"), "reference": ref.path})
    level = (float(np.nanmax(np.polyval(ref.coeffs, slit))) if (fs_correction and ref)
             else (ref.ef_kinetic() if ref else None))

    ok = np.ones(hv.size, dtype=bool)
    if source == "self":
        edc = np.nansum(cube, axis=(0, 1))
        try:
            centre = float(initial_guess(energy, edc, T)["ef"])
        except ValueError:
            centre = float(np.mean(energy))
        e0 = int(np.argmin(np.abs(energy - (centre - edge_window_eV))))
        e1 = int(np.argmin(np.abs(energy - (centre + edge_window_eV))))
        margin = int(round(slit.size * (1 - float(angle_fraction)) / 2))
        region = ((margin, slit.size - 1 - margin), (min(e0, e1), max(e0, e1)))
        try:
            ef, ok, _fits = kzmap.fit_levels(cube, energy, region, temperature=T)
        except ValueError as exc:
            raise Skip(f"{exc} -- a semiconductor has no edge of its own; use "
                       f"source: reference (gold) or file")
        params.update(angle_index_from=region[0][0], angle_index_to=region[0][1],
                      energy_index_from=region[1][0], energy_index_to=region[1][1],
                      temperature_K=T)
        if not relative:
            params["work_function_eV"] = float(np.mean(hv[ok] - ef[ok]))
    elif source == "reference":
        if relative:
            raise Skip("the energy axis is already relative to E_F (the loader "
                       "referred it); use source: self or file")
        if ref is None:
            raise Skip("no matching reference: " + "; ".join(ctx.get("reference_notes") or []))
        ef = level + (hv - float(ref.photon_energy_eV))
        params.update(reference=ref.path, work_function_eV=float(ref.work_function_eV))
    elif source == "file":
        if not relative:
            raise Skip("the energy axis is kinetic; source: file needs one the loader "
                       "already referred to E_F")
        ef = np.zeros(hv.size)
    else:
        raise ValueError(f"kz_calibrate: unknown source {source!r}")

    aligned, axis, trimmed = kzmap.align(cube, energy, ef)
    if normalise:
        aligned = kzmap.normalise_totals(aligned)
    finite = np.asarray(ef)[ok]
    spread = float(np.ptp(finite)) if finite.size else 0.0
    flat, _pos, _ang = kzconv.edge_flatness(aligned, slit, axis)
    labels = dict(ds.labels)
    labels["x"] = labels.get("x") if "hoton" in str(labels.get("x", "")) else "Photon energy (eV)"
    labels["z"] = "E - E_F (eV)"
    params.update(normalised_by_total=bool(normalise), fermi_level_spread_eV=spread,
                  energy_trimmed_eV=float(trimmed), spectra_fitted=int(ok.sum()),
                  spectra_interpolated=int((~ok).sum()))
    extra = {"kz.fermi_level_eV": ", ".join(f"{v:.5f}" for v in ef),
             "kz.fermi_level_fitted": ", ".join(str(int(v)) for v in ok)}
    if "work_function_eV" in params:
        extra["batch.work_function_eV"] = params["work_function_eV"]
    result = out.derive(step="kz_align", prefix="kz_align", value=aligned,
                        axes=(hv, slit, axis), labels=labels, params=params,
                        suffix="E-Ef", extra_info=extra)
    qc = {"source": source, "energy_trimmed_eV": float(trimmed),
          "ef_per_hv": [[float(h), float(e), bool(k)] for h, e, k in zip(hv, ef, ok)]}
    if source == "self":
        # How far the fitted edges scatter is the quality of the calibration;
        # the edge flatness after it only means something where there is an
        # edge, i.e. here.
        qc.update(fermi_level_spread_eV=spread, spectra_fitted=int(ok.sum()),
                  spectra_interpolated=int((~ok).sum()), edge_flatness_eV=float(flat))
    if "work_function_eV" in params:
        qc["work_function_eV"] = params["work_function_eV"]
    if notes:
        qc["notes"] = notes
    return result, qc


def _busy_energies(cube, slit, energy, angle_offset, count=4, separation=0.15,
                   window=1.0):
    """Up to ``count`` binding energies, at least ``separation`` apart, where
    normal emission changes most with photon energy."""
    near = np.abs(slit - angle_offset) <= window
    if not near.any():
        near[int(np.argmin(np.abs(slit - angle_offset)))] = True
    with np.errstate(invalid="ignore", divide="ignore"):
        profile = np.nanmean(cube[:, near, :], axis=1)
        contrast = np.nanstd(profile, axis=0) * np.sqrt(np.clip(np.nanmean(profile, 0), 0, None))
    contrast = np.where(energy <= 0.0, np.nan_to_num(contrast), 0.0)
    chosen = []
    for index in np.argsort(contrast)[::-1]:
        if contrast[index] <= 0:
            break
        if all(abs(energy[index] - e) >= separation for e in chosen):
            chosen.append(float(energy[index]))
        if len(chosen) == count:
            break
    return sorted(chosen)


def v0_by_period_overlap(hv, slit, energy, cube, *, period, energies, inner_potentials,
                         kpar_halfwidth=0.3, n_kz=400, **geometry):
    """The inner potential at which the k_z profile at normal emission
    overlaps itself best when shifted by exactly one period.

    At the right V0 the profile repeats with the lattice's period; at a wrong
    one the sqrt(E_kin + V0) mapping stretches it unevenly, so no single
    shift lines it up. Unlike measuring the dominant period of the profile
    (``tools.kzconv.scan_inner_potential``), this needs no sinusoidal shape
    and no Fourier peak from two or three cycles. On synthetic scans
    covering 2.5 zones it recovered V0 = 8, 14, 20, 25 eV to within 1 eV,
    where the period scan found no crossing. Each energy votes; the spread
    of the votes is the uncertainty reported.
    """
    from tools import kzconv
    curves = []
    for binding in energies:
        index = int(np.argmin(np.abs(energy - binding)))
        row = []
        for v0 in inner_potentials:
            try:
                kz, _kp, _e, out = kzconv.to_kz_cube(
                    hv, slit, energy[index:index + 1], cube[:, :, index:index + 1],
                    inner_potential=float(v0), n_kz=n_kz, n_kpar=64,
                    kpar_range=(-kpar_halfwidth, kpar_halfwidth), **geometry)
            except ValueError:
                row.append(np.nan)
                continue
            with np.errstate(invalid="ignore"):
                profile = np.nanmean(out[:, :, 0], axis=1)
            lag = int(round(period / (kz[1] - kz[0])))
            a, b = profile[:-lag], profile[lag:]
            good = np.isfinite(a) & np.isfinite(b)
            if lag < 1 or good.sum() < 20:
                row.append(np.nan)
                continue
            a, b = a[good] - a[good].mean(), b[good] - b[good].mean()
            denominator = np.sqrt(np.sum(a * a) * np.sum(b * b))
            row.append(float(np.sum(a * b) / denominator) if denominator else np.nan)
        curves.append(row)
    curves = np.array(curves, dtype=float)
    if not np.isfinite(curves).any():
        return None
    votes = [float(inner_potentials[int(np.nanargmax(c))]) for c in curves if np.isfinite(c).any()]
    combined = np.nanmean(curves, axis=0)
    best = float(inner_potentials[int(np.nanargmax(combined))])
    spread = float(np.std(votes)) if len(votes) > 1 else float("nan")
    step = float(np.mean(np.diff(inner_potentials))) if len(inner_potentials) > 1 else 0.0
    return {"best": best, "uncertainty_eV": max(spread, step), "votes": votes,
            "energies": list(map(float, energies)), "score": float(np.nanmax(combined)),
            "tried": [float(v) for v in inner_potentials],
            "curve": [None if not np.isfinite(v) else float(v) for v in combined],
            "method": "period overlap"}


def _work_function(ds, ctx, sample, given=None):
    """W for E_kin = hv - W + E: as given, from the sample block, from
    kz_calibrate (gold or the scan's own edges), from the file, from the
    reference -- or None."""
    W = given if given is not None else getattr(sample, "work_function", None)
    for key in ("batch.work_function_eV", "cassiopee.work_function_eV"):
        if W is None and ds.info.get(key) is not None:
            try:
                W = float(np.asarray(ds.info[key], dtype=float).reshape(-1)[0])
            except (TypeError, ValueError):
                pass
    ref = ctx.get("reference")
    if W is None and ref is not None:
        W = float(ref.work_function_eV)
    return W


def _normal_emission(cube, slit, energy):
    """(angle, score): the slit angle about which the scan, summed over
    photon energy near E_F, is mirror-symmetric."""
    from .center import mirror_centre
    near = (energy >= -0.6) & (energy <= 0.0)
    image = np.nansum(cube[:, :, near if near.any() else slice(None)], axis=0)
    return mirror_centre(image, slit, along=0)


def kz_match_calculation(ds: Dataset, ctx: dict, *, angle_offset="auto",
                         v0_range=(2.0, 30.0), v0_step: float = 0.5,
                         shift_range=(-0.3, 0.3), shift_step: float = 0.025,
                         renormalisations=(1.0,), energy_window=(-2.0, 0.05),
                         width: float = 0.06):
    """Match the scan at normal emission to a band calculation
    (:mod:`arpes_batch.calcbands`): V0, the calculation's rigid energy
    shift, and the photon energies at which the Gamma and zone-boundary
    planes are reached. The data is not changed; the result is recorded
    (``batch.v0_from_calculation``) for ``kz_convert`` to use when the
    sample's inner potential is ``"calculation"``, and a figure is written."""
    from . import calcbands
    _need(ds, ("kz_map",), "kz_match_calculation")
    if not _relative(ds):
        raise Skip("the energy axis is not relative to E_F yet; run kz_calibrate first")
    sample = ctx.get("sample")
    calc = getattr(sample, "calculation", None)
    if not isinstance(calc, dict):
        raise Skip("no band calculation given (sample.calculation)")
    if sample.lattice is None or sample.surface_normal is None:
        raise Skip("needs the lattice and the surface normal")
    surface = sample.surface()
    bands = calcbands.read_bands(calc["file"], labels=calc.get("path"),
                                 labels_file=calc.get("labels_file"),
                                 fermi_energy_eV=calc.get("fermi_energy_eV", 0.0),
                                 lattice=sample.lattice)
    try:
        segment = calcbands.normal_segment(bands, surface, sample.lattice)
    except ValueError as exc:
        raise Skip(str(exc))
    hv, slit, energy = ds.axis("x"), ds.axis("k"), ds.axis("z")
    cube = ds.array(float)
    W = _work_function(ds, ctx, sample)
    if W is None:
        raise Skip("needs input from the user: sample.work_function_eV")
    qc = {"calculation": bands.source, "vertices": [f"{n}@{k:.4f}" for k, n in bands.vertices],
          "normal_segment": f"G-{segment['edge_label']}", "work_function_eV": W,
          "notes": list(bands.notes)}
    if angle_offset == "auto":
        angle_offset, score = _normal_emission(cube, slit, energy)
        qc["normal_emission"] = {"angle_offset_deg": angle_offset, "score": score}
    result = calcbands.match(
        cube, hv, slit, energy, bands, segment, surface["period_invA"], work_function=W,
        angle_offset=float(angle_offset),
        v0_values=np.arange(float(v0_range[0]), float(v0_range[1]) + 1e-9, float(v0_step)),
        shifts=np.arange(float(shift_range[0]), float(shift_range[1]) + 1e-9, float(shift_step)),
        renormalisations=tuple(renormalisations), energy_window=tuple(energy_window),
        width=float(width))
    qc.update({k: v for k, v in result.items() if k not in ("scores",)})
    if ctx.get("figure_prefix"):
        from .preview import plot_calc_match
        qc["figure"] = plot_calc_match(ds, bands, segment, result, surface,
                                       ctx["figure_prefix"] + "__calc_match.png",
                                       work_function=W, angle_offset=float(angle_offset))
    params = {"calculation": bands.source, "inner_potential": result["inner_potential_eV"],
              "energy_shift": result["energy_shift_eV"],
              "renormalisation": result["renormalisation"], "score": result["score"],
              "work_function": W, "angle_offset": float(angle_offset),
              "planes": "; ".join(f"{p['plane']}{p['order']} at {p['hv_eV']:.1f} eV"
                                  for p in result["planes"])}
    out = ds.derive(step="kz_match_calculation", prefix="kzcalc", params=params,
                    suffix="matched",
                    extra_info={"batch.v0_from_calculation": result["inner_potential_eV"],
                                "batch.calc_energy_shift_eV": result["energy_shift_eV"]})
    out.name = ds.name                     # nothing was done to the data itself
    return out, qc


def _measured_period(kz_axis, kpar_axis, plane, target, halfwidth=0.3):
    band = np.abs(kpar_axis) <= halfwidth
    with np.errstate(invalid="ignore"):
        profile = np.nanmean(plane[:, band], axis=1)
    good = np.isfinite(profile)
    if good.sum() < 16:
        return float("nan")
    # (the same measurement as tools.kzconv._dominant_period)
    a, p = kz_axis[good], profile[good] - np.nanmean(profile[good])
    trial = np.linspace(0.5 * target, 2.0 * target, 600)
    phase = 2 * np.pi * a[None, :] / trial[:, None]
    amplitude = np.hypot((p * np.cos(phase)).sum(1), (p * np.sin(phase)).sum(1))
    return float(trial[int(np.argmax(amplitude))])


def kz_convert(ds: Dataset, ctx: dict, *, inner_potential=None, work_function=None,
               angle_offset=0.0, theta_position: float = None, effective_mass: float = None,
               n_kz: int = 256, n_kpar: int = 256, scan_energy: float = None,
               v0_range=(2.0, 30.0), v0_step: float = 0.5, kpar_halfwidth: float = 0.3):
    """kz map -> (k_z, k_par, E) (the viewer's "kz -> momentum").

    Needs the sample block (lattice, surface normal, V0 -- see
    :mod:`arpes_batch.sample`); the runner refuses a recipe without it
    before anything runs. ``inner_potential`` here overrides the sample's;
    ``"scan"`` picks V0 so that the data repeats with the surface's period.
    ``angle_offset`` is normal emission on the slit (``"auto"``: the
    mirror centre of the scan summed over photon energy).
    """
    from tools import cleavage, kzconv
    from .sample import missing_inputs
    _need(ds, ("kz_map",), "kz_convert")
    if not _relative(ds):
        raise Skip("the energy axis is not relative to E_F yet; run kz_calibrate first")
    sample = ctx.get("sample")
    missing = ([m for m in missing_inputs(sample) if not m.get("optional")]
               if sample is not None else [{"field": "sample"}])
    if missing:
        raise Skip("needs input from the user: " + ", ".join(m["field"] for m in missing))
    surface = sample.surface()
    period = surface["period_invA"]
    hv, slit, energy = ds.axis("x"), ds.axis("k"), ds.axis("z")
    cube = ds.array(float)
    qc, notes = {"surface": surface}, []

    W = _work_function(ds, ctx, sample, work_function)
    if W is None:
        raise Skip("needs input from the user: sample.work_function_eV (no reference, "
                   "none in the file)")
    if theta_position is None:
        theta_position = float(ds.info.get("cassiopee.sample_theta_deg", 0.0) or 0.0)
        if not theta_position:
            notes.append("manipulator theta taken as 0 (not recorded); set theta_position "
                         "if the sample was tilted")
    m_eff = float(effective_mass or sample.effective_mass or 1.0)

    if angle_offset == "auto":
        angle_offset, score = _normal_emission(cube, slit, energy)
        qc["normal_emission_suggestion"] = {"angle_offset_deg": angle_offset, "score": score}
    angle_offset = float(angle_offset)
    geometry = dict(work_function=float(W), effective_mass=m_eff,
                    angle_offset=angle_offset, theta_position=float(theta_position))

    energies = _busy_energies(cube, slit, energy, angle_offset)
    if scan_energy is None:
        scan_energy = energies[0] if energies else float(np.median(energy))
    V0 = inner_potential if inner_potential is not None else sample.inner_potential
    from_calc = ds.info.get("batch.v0_from_calculation")
    if from_calc is not None:
        qc["v0_from_calculation"] = float(np.asarray(from_calc, dtype=float).reshape(-1)[0])
    if V0 == "calculation":
        if from_calc is None:
            raise Skip("inner_potential_eV is \"calculation\" but no match to a calculation "
                       "was made; put kz_match_calculation before kz_convert")
        V0 = qc["v0_from_calculation"]
    if V0 == "scan":
        span = float(np.ptp(kzconv.kz_bounds(hv, slit, energy, inner_potential=15.0,
                                              **{k: geometry[k] for k in
                                                 ("work_function", "effective_mass")})[0]))
        if span < 1.5 * period:
            raise Skip(f"this scan covers about {span / period:.1f} zones along k_z; "
                       f"V0 cannot be estimated from fewer than ~1.5. Give "
                       f"inner_potential_eV (a value from the literature, or from a "
                       f"wider scan of the same crystal)")
        trials = np.arange(float(v0_range[0]), float(v0_range[1]) + 1e-9, float(v0_step))
        found = v0_by_period_overlap(hv, slit, energy, cube, period=period,
                                     energies=energies or [scan_energy],
                                     inner_potentials=trials,
                                     kpar_halfwidth=kpar_halfwidth, **geometry)
        # The viewer's own estimate, for comparison (it measures the
        # dominant period rather than the overlap; see v0_by_period_overlap).
        viewer = kzconv.scan_inner_potential(
            hv, slit, energy, cube, spacing=2 * np.pi / period, binding_energy=scan_energy,
            kpar_halfwidth=kpar_halfwidth, inner_potentials=trials, **geometry)
        qc["v0_scan"] = dict(found or {"best": None}, binding_energy=scan_energy,
                             period_scan={"best": viewer.best,
                                          "uncertainty_eV": viewer.uncertainty(),
                                          "periods": [float(p) for p in viewer.periods],
                                          "target": viewer.target},
                             target=period)
        if not found:
            raise Skip("the V0 scan could not compare the data with the surface's "
                       "period: check the lattice and the surface normal, or give "
                       "inner_potential_eV")
        V0 = found["best"]
        if found["best"] in (trials[0], trials[-1]):
            notes.append(f"V0 = {V0:g} eV is at the end of the range tried "
                         f"({trials[0]:g}-{trials[-1]:g} eV); widen v0_range")
    V0 = float(V0)
    kz_axis, kpar_axis, e_out, out = kzconv.to_kz_cube(
        hv, slit, energy, cube, inner_potential=V0, n_kz=int(n_kz), n_kpar=int(n_kpar),
        **geometry)
    zones = float(np.ptp(kz_axis) / period)
    qc["zones_covered"] = zones
    measured = float("nan")
    if zones >= 1.5:
        index = int(np.argmin(np.abs(e_out - scan_energy)))
        measured = _measured_period(kz_axis, kpar_axis, out[:, :, index], period,
                                    kpar_halfwidth)
        qc["period_check"] = {"expected_invA": period, "measured_invA": measured,
                              "binding_energy": scan_energy,
                              "relative_difference": (measured - period) / period
                              if np.isfinite(measured) else None}
    else:
        notes.append(f"the scan covers {zones:.2f} zones along k_z: too few to check "
                     f"the period (or to scan V0)")
    if np.isfinite(measured) and abs(measured - period) / period > 0.1 and sample.lattice:
        found = cleavage.candidates(measured, sample.lattice, tolerance=0.1)
        notes.append(f"the converted map repeats every {measured:.3f} A^-1, not the "
                     f"{period:.3f} A^-1 of the surface given; planes that would match: "
                     + "; ".join(c.describe() for c in found[:3]) if found else
                     f"the converted map repeats every {measured:.3f} A^-1, not {period:.3f}")
    qc["inside_percent"] = float(np.isfinite(out).mean() * 100)
    lattice = sample.lattice
    params = {**geometry, "inner_potential": V0, "n_kz": int(n_kz), "n_kpar": int(n_kpar),
              "space_group": lattice.space_group or lattice.centering,
              "lattice_a": lattice.a, "lattice_b": lattice.b, "lattice_c": lattice.c,
              "lattice_alpha": lattice.alpha, "lattice_beta": lattice.beta,
              "lattice_gamma": lattice.gamma, "surface_hkl": list(surface["hkl"]),
              "surface_period_invA": period}
    labels = {"x": f"k_z ({ANGSTROM})", "k": f"k_par ({ANGSTROM})", "z": ds.label("z")}
    result = ds.derive(step="kz_to_k", prefix="kz_to_k", kind="kz_map_k",
                       value=out.astype(np.float32), axes=(kz_axis, kpar_axis, e_out),
                       labels=labels, params=params, suffix="kz",
                       extra_info={"kz.inner_potential_eV": V0,
                                   "batch.kz_period_invA": period})
    if notes:
        qc["notes"] = notes
    return result, qc


# -- cuts ------------------------------------------------------------------------------
def kconvert_cut(ds: Dataset, ctx: dict, *, gamma_slit_deg=0.0, gamma_deflector_deg=0.0,
                 azimuth_deg: float = 0.0, n_k: int = 300, radial: bool = False,
                 centre_window=(-0.6, 0.0)):
    """Cut k conversion (``tools.cutk.convert_cut``). ``gamma_slit_deg``
    ``"auto"`` takes the slit angle about which the cut is mirror-symmetric
    -- right for a cut through Gamma along a mirror line, meaningless
    otherwise, so its score is reported."""
    from tools.cutk import K_LABEL, K_RADIAL_LABEL, convert_cut
    from .center import mirror_centre
    _need(ds, ("cut",), "kconvert_cut")
    if (ctx.get("row") or {}).get("slit_axis_looks_like") == "momentum":
        raise Skip("the slit axis looks like it is already in 1/A")
    slit, energy = ds.axis("x"), ds.axis("y")
    frame = ds.array(float)
    qc = {}
    offset = float(ds.info.get("batch.ef_kinetic_eV", 0.0) or 0.0)
    if not offset and _relative(ds):
        raise Skip("energy is relative to E_F but the kinetic E_F is not recorded; "
                   "run calibrate_energy in the same recipe")
    if gamma_slit_deg == "auto":
        lo, hi = centre_window if _relative(ds) else (energy.max() - 0.6, energy.max())
        rows = (energy >= lo) & (energy <= hi)
        gamma_slit_deg, score = mirror_centre(frame[:, rows if rows.any() else slice(None)],
                                              slit, along=0)
        qc["centre_suggestion"] = {"gamma_slit_deg": gamma_slit_deg, "score": score}
    deflector = float(ds.info.get("DeflectorAngle_deg", 0.0) or 0.0)
    settings = {"deflector_deg": deflector, "gamma_deflector_deg": float(gamma_deflector_deg),
                "gamma_slit_deg": float(gamma_slit_deg), "azimuth_deg": float(azimuth_deg),
                "energy_offset_eV": offset, "n_k": int(n_k), "radial": bool(radial)}
    k, e_out, values, report = convert_cut(slit, energy, frame, **settings)
    qc.update({k_: v for k_, v in report.items() if k_ != "direction"})
    labels = {"x": K_RADIAL_LABEL if radial else K_LABEL, "y": ds.label("y")}
    return ds.derive(step="cut_k_conversion", prefix="kcut", value=values,
                     axes=(k, e_out - offset), labels=labels, params=settings,
                     suffix="k"), qc


def curvature(ds: Dataset, ctx: dict, *, mode: str = "2d", a0: float = None,
              factor: float = 1.0):
    """Curvature of a cut (``tools.process.curvature``), with ``a0`` from
    ``suggest_a0`` unless given. Replaces the intensity: save before it."""
    from tools import process as P
    _need(ds, ("cut",), "curvature")
    values = np.nan_to_num(ds.array(float))
    axes = [ds.axis("x"), ds.axis("y")]
    a0 = float(a0) if a0 is not None else float(P.suggest_a0(values, axes, mode=mode,
                                                              factor=factor))
    out = P.curvature(values, axes, mode=mode, a0=a0)
    return ds.derive(step="curvature", params={"mode": mode, "a0": a0},
                     value=out), {"a0": a0}


STEPS = {
    "degrid": degrid,
    "calibrate_energy": calibrate_energy,
    "crop": crop,
    "bin": bin_,
    "normalise": normalise,
    "kconvert": kconvert,
    "kz_calibrate": kz_calibrate,
    "kz_match_calculation": kz_match_calculation,
    "kz_convert": kz_convert,
    "kconvert_cut": kconvert_cut,
    "curvature": curvature,
}
