"""
The analyser's Fermi edge, from a gold reference, without clicking.

In the viewer this is the "Fermi surface from a reference" dialog
(``ui/windows.AuReferenceDialog``): the reference is summed over the
deflector, its edge fitted channel by channel along the slit
(``tools.fermi.fit_channels``), and a polynomial through the fitted levels
becomes the Fermi-surface correction. That dialog keeps its logic in the Qt
module, so the same steps are written out here on the Qt-free pieces.

What a reference gives a batch run
----------------------------------
* the **curvature** of the edge along the slit -- the polynomial -- which
  depends on the lens mode and pass energy, not on the sample;
* the **position** of the edge in kinetic energy at the reference's photon
  energy. At another photon energy the edge moves with it
  (``E_F,kin = hv - phi``), so a reference taken at ``hv_ref`` places a map
  taken at ``hv`` at ``E_F(hv) = E_F(hv_ref) + (hv - hv_ref)``. That holds
  only as far as the monochromator's nominal energy is right; the offset
  and the photon energies it was applied across are recorded on the result.

A semiconductor (WS2, WSe2 ...) has no edge of its own to fit, which is why
the reference, and not the map, is where E_F comes from by default.
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field

import numpy as np

from tools.fermi import fit_channels, fit_fermi_edge, steepest_drop

from .dataset import Dataset
from .inventory import META_KEYS, _number, meta


@dataclass
class Reference:
    path: str
    entry: str
    photon_energy_eV: float
    pass_energy: str
    lens_mode: str
    temperature_K: float
    start_time: str
    order: int
    coeffs: list                  # E_F,kin(slit angle), highest power first
    slit_range: tuple
    ef_center_eV: float           # at the middle of the slit
    ef_flat_eV: float             # the level fs_correction flattens onto
    work_function_eV: float       # hv - E_F,kin at the middle of the slit
    resolution_eV: float
    residual_meV: float
    n_ok: int
    n_channels: int
    channels: dict = field(default_factory=dict)   # angle, ef, err, ok

    def ef_kinetic(self, angle=None, photon_energy_eV=None):
        """E_F in kinetic energy at ``angle`` (the middle of the slit by
        default), moved to another photon energy if one is given."""
        angle = np.mean(self.slit_range) if angle is None else angle
        value = np.polyval(self.coeffs, angle)
        if photon_energy_eV is not None and self.photon_energy_eV is not None:
            value = value + (float(photon_energy_eV) - float(self.photon_energy_eV))
        return value

    def to_json(self) -> dict:
        d = asdict(self)
        d["channels"] = {k: [None if not np.isfinite(v) else float(v) for v in vals]
                         if k != "ok" else [bool(v) for v in vals]
                         for k, vals in self.channels.items()}
        return d

    @classmethod
    def from_json(cls, d: dict) -> "Reference":
        d = dict(d)
        d["slit_range"] = tuple(d["slit_range"])
        d["channels"] = {k: np.array([np.nan if v is None else v for v in vals],
                                     dtype=bool if k == "ok" else float)
                         for k, vals in (d.get("channels") or {}).items()}
        return cls(**d)

    def save(self, path: str) -> str:
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(self.to_json(), fh, indent=1)
        return path


def reference_frame(ds: Dataset):
    """(slit angle, energy, frame): a cut as it is, a map summed over the
    deflector -- every deflector position sees the same slit curvature."""
    if ds.kind == "cut":
        return ds.axis("x"), ds.axis("y"), ds.array(float)
    if ds.kind == "map":
        return ds.axis("k"), ds.axis("z"), np.nansum(ds.array(np.float64), axis=0)
    raise ValueError(f"a {ds.kind} cannot be used as a Fermi-edge reference")


def _robust_polyfit(x, y, w, order, clip=3.0, rounds=5):
    keep = np.isfinite(x) & np.isfinite(y) & np.isfinite(w)
    coeffs = None
    for _ in range(rounds):
        if keep.sum() < order + 2:
            break
        coeffs = np.polyfit(x[keep], y[keep], order, w=w[keep])
        resid = y - np.polyval(coeffs, x)
        scale = 1.4826 * np.nanmedian(np.abs(resid[keep])) or np.nanstd(resid[keep])
        new = keep & (np.abs(resid) <= clip * max(scale, 1e-6))
        if new.sum() == keep.sum():
            break
        keep = new
    if coeffs is None:
        raise ValueError("too few channels with a usable edge to fit a curve through")
    return coeffs, keep


def fit_reference(path: str, entry: str = None, **options) -> Reference:
    """Fit the gold reference in ``path`` (see :func:`fit_reference_on`)."""
    ds = Dataset.load(path, entry)
    try:
        return fit_reference_on(ds, **options)
    finally:
        ds.close()


def fit_reference_on(ds: Dataset, *, order: int = 2, half_width: int = 3,
                     step: int = 4, window_eV: float = 0.25,
                     temperature: float = None) -> Reference:
    """The edge of ``ds`` along its slit: a polynomial through the
    channel-by-channel E_F, with outlying channels rejected."""
    info = ds.info
    slit, energy, frame = reference_frame(ds)
    temperature = float(temperature or _number(meta(info, META_KEYS["temperature_K"])) or 30.0)

    # Angle-integrated edge first: it places the window for every channel
    # and gives the resolution, which a single channel is too noisy for.
    edc = np.nansum(frame, axis=0)
    guess = steepest_drop(energy, edc)
    window = (guess - window_eV, guess + window_eV)
    whole = fit_fermi_edge(energy, edc, temperature=temperature, window=window)
    ef0 = whole.values["ef"]
    window = (ef0 - window_eV, ef0 + window_eV)

    ef, err, ok = fit_channels(slit, energy, frame, half_width=half_width,
                               step=step, temperature=temperature,
                               window=window, start=dict(whole.values))
    weights = np.where(ok & np.isfinite(err) & (err > 0), 1.0 / np.maximum(err, 1e-6), 0.0)
    coeffs, kept = _robust_polyfit(slit, ef, weights, order)
    ok = ok & kept
    resid = (ef - np.polyval(coeffs, slit))[ok]

    hv = _number(meta(info, META_KEYS["photon_energy_eV"]))
    mid = float(np.mean([slit.min(), slit.max()]))
    ef_mid = float(np.polyval(coeffs, mid))
    fitted = np.polyval(coeffs, slit)
    return Reference(
        path=ds.source_path, entry=str(info.get("_group", "")).lstrip("/"),
        photon_energy_eV=hv,
        pass_energy=str(meta(info, META_KEYS["pass_energy"]) or ""),
        lens_mode=str(meta(info, META_KEYS["lens_mode"]) or ""),
        temperature_K=temperature,
        start_time=str(meta(info, META_KEYS["start_time"]) or ""),
        order=int(order), coeffs=[float(c) for c in coeffs],
        slit_range=(float(slit.min()), float(slit.max())),
        ef_center_eV=ef_mid, ef_flat_eV=float(np.nanmax(fitted)),
        work_function_eV=(float(hv) - ef_mid) if hv is not None else float("nan"),
        resolution_eV=float(whole.values["resolution"]),
        residual_meV=float(np.sqrt(np.mean(resid ** 2)) * 1000) if resid.size else float("nan"),
        n_ok=int(ok.sum()), n_channels=int(np.isfinite(ef).sum()),
        channels={"angle": slit, "ef": ef, "err": err, "ok": ok})


def _same(a, b) -> bool:
    return str(a or "").strip().upper() == str(b or "").strip().upper()


def choose_reference(row: dict, references: list):
    """The reference for one map: same lens mode and pass energy (the
    curvature depends on both), then the nearest photon energy, then the
    nearest in time. Returns ``(reference, notes)``; ``None`` if nothing
    matches the analyser settings."""
    notes = []
    if not references:
        return None, ["no reference available"]
    matching = [r for r in references
                if _same(r.lens_mode, row.get("lens_mode"))
                and _same(r.pass_energy, row.get("pass_energy"))]
    if not matching:
        notes.append(f"no reference with lens mode {row.get('lens_mode')!r} and "
                     f"pass energy {row.get('pass_energy')!r}")
        return None, notes
    hv = row.get("photon_energy_eV")
    start = str(row.get("start_time") or "")

    def key(r):
        dhv = abs(float(hv) - float(r.photon_energy_eV)) if (
            hv is not None and r.photon_energy_eV is not None) else 0.0
        return (round(dhv, 2), abs(_seconds(start) - _seconds(r.start_time)))

    best = sorted(matching, key=key)[0]
    if hv is not None and best.photon_energy_eV is not None and abs(
            float(hv) - float(best.photon_energy_eV)) > 0.01:
        notes.append(f"reference taken at hv={best.photon_energy_eV:g} eV, map at "
                     f"{float(hv):g} eV: E_F moved by the difference")
    return best, notes


def _seconds(text: str) -> float:
    import datetime
    text = str(text or "").strip()
    for fmt in ("%Y-%m-%dT%H:%M:%S.%f", "%Y-%m-%dT%H:%M:%S", "%Y-%m-%d %H:%M:%S",
                "%Y-%m-%dT%H:%M:%S%z", "%Y-%m-%dT%H:%M:%S.%f%z"):
        try:
            return datetime.datetime.strptime(text[:26] if "%z" not in fmt else text, fmt).timestamp()
        except ValueError:
            continue
    return 0.0
