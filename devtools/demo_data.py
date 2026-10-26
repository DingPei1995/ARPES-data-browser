"""Synthetic ARPES data, for trying the program without a beamtime.

Writes a handful of files in the program's own saved format -- a Fermi
surface map, single cuts (including an LH / LV pair), a photon-energy scan,
a spin-resolved EDC and a small real-space scan -- generated from a simple
band model with a Fermi-Dirac cut-off, the analyser resolution and Poisson
counting noise. They open like measured data, and are what the pictures in
``docs/`` were taken from (``devtools/guide_screenshots.py``).

The model is a toy, not a material: a hexagonal electron pocket and a hole
band at Gamma on a triangular lattice (a = 3.3 A), plus a k_z-dispersing
band for the photon-energy scan (c = 6 A, V0 = 12 eV).

Usage:  python devtools/demo_data.py [output folder]     (default: demo_data/)
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from loader.nxs_file import save_dataset            # noqa: E402

A_K = 0.512317          # k [A^-1] = A_K sqrt(E_kin [eV]) sin(angle)
A_KZ = 0.262466         # k^2 [A^-2] = A_KZ E [eV]
HV = 60.0
WORK_FUNCTION = 4.4
EF_KIN = HV - WORK_FUNCTION
LATTICE_A = 3.3
KB = 8.617e-5


def _hexagonal(kx, ky, a=LATTICE_A):
    """6 - f(k) for the triangular lattice: 0 at Gamma, 9 at K."""
    return 6.0 - (2.0 * np.cos(kx * a)
                  + 4.0 * np.cos(kx * a / 2) * np.cos(np.sqrt(3) * ky * a / 2))


def _bands(kx, ky):
    """(energy, weight) of each band at (kx, ky), energies relative to E_F."""
    h = _hexagonal(kx, ky)
    return [(-0.75 + 0.120 * h, 1.0),                  # hexagonal electron pocket
            (-0.25 - 0.100 * h, 0.75),                 # hole band
            (0.25 - 0.150 * (9.0 - h), 0.8),           # hole pockets at K
            (-1.20 + 0.012 * h, 0.45)]                 # nearly flat, deeper


def pocket_k(energy):
    """|k| of the electron pocket along ky (Gamma-M) at ``energy``: where
    a cut through Gamma along the slit crosses it."""
    x = np.arccos(1.0 - (np.asarray(energy, dtype=float) + 0.75) / 0.48)
    return x / (np.sqrt(3) * LATTICE_A / 2)


def _spectral(energy, bands, temperature=30.0, resolution=0.015, ef=0.0):
    """I(..., E): Lorentzians with a lifetime growing away from E_F, a small
    background, the Fermi cut-off at ``ef`` (a number, or one per spectrum),
    then the resolution along the last axis."""
    from scipy.ndimage import gaussian_filter1d

    e = energy
    out = 0.04 + 0.03 * np.clip(-e, 0, None)          # inelastic background
    for centre, weight in bands:
        centre = np.asarray(centre, dtype=float)
        weight = np.broadcast_to(np.asarray(weight, dtype=float), centre.shape)
        width = 0.025 + 0.12 * (centre[..., None] ** 2)
        out = out + weight[..., None] * (width / np.pi) / (
            (e - centre[..., None]) ** 2 + width ** 2) * 0.05
    ef = np.asarray(ef, dtype=float)
    shift = ef[..., None] if ef.ndim else ef
    out = out / (np.exp(np.clip((e - shift) / (KB * temperature), -60, 60)) + 1.0)
    step = abs(float(energy[1] - energy[0]))
    return gaussian_filter1d(out, resolution / 2.355 / step, axis=-1, mode="nearest")


def _counts(intensity, peak, rng):
    return rng.poisson(intensity * (peak / intensity.max())).astype(np.float32)


def _info(**extra):
    info = {"photon_energy_eV": HV, "SampleTemperature_K": 30.0,
            "pass_energy_eV": 20.0, "lens_mode": "WideAngle",
            "work_function_eV": WORK_FUNCTION, "demo": "synthetic data"}
    info.update(extra)
    return info


def fermi_surface_map(rng, *, grid=0.0):
    """(deflector, slit, kinetic energy) map. ``grid`` > 0 prints a square
    detector mesh of that contrast on every slice, for the De-grid tool."""
    defl = np.linspace(-24, 24, 97)
    slit = np.linspace(-24, 24, 241)
    ekin = np.linspace(EF_KIN - 1.5, EF_KIN + 0.15, 166)
    e = ekin - EF_KIN
    D, S = np.meshgrid(np.radians(defl), np.radians(slit), indexing="ij")
    k0 = A_K * np.sqrt(EF_KIN)
    kx, ky = k0 * np.sin(D), k0 * np.cos(D) * np.sin(S)
    bands = []
    for centre, weight in _bands(kx, ky):
        # a matrix element that favours one side, as real data does
        bands.append((centre, weight * (1.0 + 0.35 * np.sin(S + 0.3))))
    cube = _counts(_spectral(e, bands), 160.0, rng)
    if grid:
        i = np.arange(slit.size)[:, None]
        j = np.arange(ekin.size)[None, :]
        mesh = 1.0 - grid * (0.5 + 0.5 * np.cos(2 * np.pi * i / 10.7)
                             * np.cos(2 * np.pi * j / 11.3))
        cube = rng.poisson(cube * mesh[None] + 0.5).astype(np.float32)
    labels = {"x": "Angle, deflector (deg)", "k": "Angle, slit (deg)",
              "z": "Kinetic energy (eV)"}
    return {"kind": "map", "axes": (defl, slit, ekin), "value": cube,
            "labels": labels, "info": _info(**{"axis0.role": "angle"})}


def single_cut(rng, *, polarisation="LH", peak=1500.0):
    """(slit, kinetic energy) cut through Gamma."""
    slit = np.linspace(-15, 15, 401)
    ekin = np.linspace(EF_KIN - 1.0, EF_KIN + 0.12, 225)
    e = ekin - EF_KIN
    S = np.radians(slit)
    k = A_K * np.sqrt(EF_KIN) * np.sin(S)
    bands = _bands(np.zeros_like(k), k)
    if polarisation == "LH":
        scale = [1.0 + 0.5 * np.sin(S + 0.2), 0.45 + 0 * S, 1.0 + 0 * S, 1.0 + 0 * S]
    else:                                   # LV: the other orbital character
        scale = [0.55 + 0 * S, 1.1 - 0.4 * np.sin(S), 1.0 + 0 * S, 1.0 + 0 * S]
    bands = [(c, w * s) for (c, w), s in zip(bands, scale)]
    cut = _counts(_spectral(e, bands, resolution=0.010), peak, rng)
    return {"kind": "cut", "axes": (slit, ekin), "value": cut,
            "labels": {"x": "Angle, slit (deg)", "y": "Kinetic energy (eV)"},
            "info": _info(**{"cassiopee.polarisation": polarisation})}


def photon_energy_scan(rng):
    """(hv, slit, E - E_F) scan with a band dispersing along k_z, each
    spectrum's Fermi level slightly off -- which kz map processing fixes."""
    hv = np.arange(30.0, 90.01, 1.5)
    slit = np.linspace(-12, 12, 121)
    e = np.linspace(-1.0, 0.15, 116)
    c, v0 = 6.0, 12.0
    H, S = np.meshgrid(hv, np.radians(slit), indexing="ij")
    ekin = H - WORK_FUNCTION
    kpar = A_K * np.sqrt(ekin) * np.sin(S)
    kz = np.sqrt(np.clip(A_KZ * (ekin + v0) - kpar ** 2, 0, None))
    centre = -0.42 - 0.28 * np.cos(kz * c) + 1.8 * kpar ** 2
    deeper = -0.85 + 0.6 * kpar ** 2 + 0.05 * np.cos(kz * c)
    # each spectrum's E_F a little off zero, as nominal photon energies leave it
    offset = (0.025 * np.sin(hv / 6.0))[:, None] + 0 * kpar
    flux = (1.0 + 2.5 * np.exp(-((hv - 55) / 25) ** 2))[:, None, None]
    intensity = _spectral(e, [(centre + offset, np.ones_like(centre)),
                              (deeper + offset, 0.6 * np.ones_like(centre))],
                          ef=offset)
    cube = _counts(intensity * flux, 120.0, rng)
    labels = {"x": "Photon energy (eV)", "k": "Angle, slit (deg)",
              "z": "E - E_F (eV)"}
    return {"kind": "kz_map", "axes": (hv, slit, e), "value": cube,
            "labels": labels,
            "info": _info(**{"axis0.role": "photon_energy",
                             "photon_energy_eV": float(hv.mean())})}


SPIN_LABELS = ["<0,0> +X  (GUI+Z)", "<90,180> +X  (GUI-Z)",
               "<0,0> -X  (GUI-Z)", "<90,180> -X  (GUI+Z)"]


def spin_edc(rng, sherman=0.2, eps=0.02):
    """Four spin channels of one EDC, with a polarisation that changes sign
    between two bands and a small instrumental asymmetry."""
    e = np.linspace(-1.0, 0.2, 121)
    spectrum = _spectral(e, [(np.array(-0.25), np.array(1.0)),
                             (np.array(-0.60), np.array(0.8))])
    spectrum = 2500.0 * spectrum / spectrum.max() + 40.0
    pol = (0.35 * np.exp(-((e + 0.25) / 0.08) ** 2)
           - 0.30 * np.exp(-((e + 0.60) / 0.10) ** 2))
    up, down = 1 + sherman * pol, 1 - sherman * pol
    n = (1.0, 0.8)
    mean = np.stack([n[0] * (1 + eps) * up * spectrum,
                     n[1] * (1 + eps) * down * spectrum,
                     n[0] * (1 - eps) * down * spectrum,
                     n[1] * (1 - eps) * up * spectrum], axis=1)
    counts = rng.poisson(mean).astype(float)
    info = _info(**{"curve.channels": "|".join(f"C{i} {l}" for i, l in enumerate(SPIN_LABELS)),
                    "curve.value_label": "Counts"})
    for i, label in enumerate(SPIN_LABELS):
        info[f"spin.component.{i}"] = label
    return {"kind": "spin_edc", "axes": (e, np.arange(4.0)), "value": counts,
            "labels": {"x": "E - E_F (eV)", "y": "Channel"}, "info": info}


def spatial_scan(rng):
    """(x, y, slit, E) raster: a flake with two domains on a substrate."""
    x = np.linspace(0, 58, 30)                       # um
    y = np.linspace(0, 48, 25)
    slit = np.linspace(-12, 12, 64)
    e = np.linspace(-1.0, 0.1, 56)
    X, Y = np.meshgrid(x, y)                          # (y, x)
    flake = ((X - 30) ** 2 / 22 ** 2 + (Y - 24) ** 2 / 17 ** 2) < 1
    domain = flake & (X + 0.6 * Y > 42)
    S = np.radians(slit)
    k = A_K * np.sqrt(EF_KIN) * np.sin(S)
    a = _spectral(e, [(c, w) for c, w in _bands(np.zeros_like(k), k)])
    b = 0.7 * _spectral(e, [(c + 0.15, w) for c, w in _bands(np.zeros_like(k), k)])
    sub = _spectral(e, [(np.full_like(k, -0.9), np.full_like(k, 0.3))])
    cube = np.where(domain[..., None, None], b[None, None],
                    np.where(flake[..., None, None], a[None, None], sub[None, None]))
    cube = _counts(cube, 60.0, rng)
    return {"kind": "spem_4d", "axes": (x, y, slit, e), "value": cube,
            "labels": {"x": "x (µm)", "y": "y (µm)",
                       "k": "Angle, slit (deg)", "z": "E - E_F (eV)"},
            "info": _info()}


def write_all(folder: str) -> "list[str]":
    os.makedirs(folder, exist_ok=True)
    rng = np.random.default_rng(2026)
    files = {
        "demo_fermi_surface_map.nxs": [("FS_map_hv60", fermi_surface_map(rng))],
        "demo_cuts.nxs": [("Gamma_cut_LH", single_cut(rng, polarisation="LH")),
                          ("Gamma_cut_LV", single_cut(rng, polarisation="LV"))],
        "demo_kz_scan.nxs": [("hv_scan_30-90eV", photon_energy_scan(rng))],
        "demo_spin_edc.nxs": [("spin_EDC_normal_emission", spin_edc(rng))],
        "demo_spatial_scan.nxs": [("SPEM_flake", spatial_scan(rng))],
        "demo_map_with_grid.nxs": [("FS_map_mesh", fermi_surface_map(rng, grid=0.12))],
    }
    written = []
    for filename, entries in files.items():
        path = os.path.join(folder, filename)
        save_dataset(path, [dict(item, name=name) for name, item in entries])
        written.append(path)
    return written


if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "demo_data")
    for path in write_all(out):
        print(path)
