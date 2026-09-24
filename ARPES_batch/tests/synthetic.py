"""
Synthetic ANTARES deflector maps, laid out the way ``loader/nxs_file.py``
reads a real one (its "case 4": ``scan_data/actuator_1_1`` for the deflector,
``data_04..06`` and ``data_07..09`` for the energy and slit scales,
``data_11`` for the cube, the metadata under ``ANTARES/...``).

Everything the pipeline is meant to recover is planted and known:

* E_F in kinetic energy at ``hv - phi``, with a parabolic bend along the
  slit of ``curvature_eV`` at its ends;
* Gamma at ``gamma_deg`` = (deflector, slit);
* a hexagonal, multiplicative detector grid locked to the (slit, energy)
  pixels;
* Poisson counts.
"""
from __future__ import annotations

import h5py
import numpy as np

K0 = 0.5123


def _grid(n_slit, n_e, contrast=0.04, period=5.0):
    i, j = np.meshgrid(np.arange(n_slit), np.arange(n_e), indexing="ij")
    g = np.zeros((n_slit, n_e))
    for angle in (0.0, 60.0, 120.0):
        a = np.radians(angle)
        g += np.cos(2 * np.pi * (i * np.cos(a) + j * np.sin(a)) / period)
    return 1.0 + contrast * g / 3.0


def ef_curve(slit, ef0, curvature_eV):
    half = max(abs(slit.min()), abs(slit.max()))
    return ef0 - curvature_eV * (slit / half) ** 2


def write_map(path, *, sample="band", n_defl=41, n_slit=96, n_e=130,
              defl_range=(-12.0, 12.0), slit_range=(-15.0, 15.0),
              hv=60.0, phi=4.35, curvature_eV=0.03, gamma_deg=(2.0, -1.5),
              counts=40.0, grid=0.04, temperature_K=20.0, lens="L4", pe="PE50",
              start_time="2026-09-18T15:46:23", entry="DeflX_0001", seed=0,
              raw_order=(0, 2, 1)):
    rng = np.random.default_rng(seed)
    ef0 = hv - phi
    e_lo, e_hi = ef0 - 1.6, ef0 + 0.25
    e_step = (e_hi - e_lo) / (n_e - 1)
    s_step = (slit_range[1] - slit_range[0]) / (n_slit - 1)
    defl = np.linspace(*defl_range, n_defl)
    slit = slit_range[0] + np.arange(n_slit) * s_step
    energy = e_lo + np.arange(n_e) * e_step

    D, S, E = np.meshgrid(defl, slit, energy, indexing="ij")
    ef = ef_curve(slit, ef0, curvature_eV)[None, :, None]
    kT = 8.617e-5 * max(temperature_K, 30.0) + 0.01          # thermal + resolution
    fermi = 1.0 / (np.exp((E - ef) / kT) + 1.0)
    if sample == "gold":
        spectrum = 1.0 + 0.1 * (ef - E)
    else:
        k0 = K0 * np.sqrt(E)
        kx = k0 * np.sin(np.radians(D - gamma_deg[0]))
        ky = k0 * np.sin(np.radians(S - gamma_deg[1]))
        k2 = kx ** 2 + ky ** 2
        eb = E - ef
        width = 0.06
        electron = -0.35 + 4.0 * k2              # a pocket crossing E_F, k_F ~ 0.3 1/A
        hole = -0.45 - 3.0 * k2                  # a hole band below it
        spectrum = (width ** 2 / ((eb - electron) ** 2 + width ** 2)
                    + 0.8 * width ** 2 / ((eb - hole) ** 2 + width ** 2) + 0.08)
    expected = counts * spectrum * fermi + 0.5
    if grid:
        expected = expected * _grid(n_slit, n_e, grid)[None]
    cube = rng.poisson(expected).astype(np.float32)
    raw = np.transpose(cube, raw_order)

    with h5py.File(path, "w") as f:
        g = f.create_group(entry)
        g["start_time"] = np.bytes_(start_time)
        g["title"] = np.bytes_(f"synthetic {sample} map")
        sd = g.create_group("scan_data")
        sd["actuator_1_1"] = defl
        for name, value in (("data_04", e_lo), ("data_05", e_step), ("data_06", energy[-1]),
                            ("data_07", slit[0]), ("data_08", s_step), ("data_09", slit[-1])):
            sd[name] = np.array([value])
        sd.create_dataset("data_11", data=raw, compression="gzip", compression_opts=1)
        an = g.create_group("ANTARES")
        mono = an.create_group("i12-m-c04-op-mono1")
        mono["energy"] = np.array([hv])
        mbs = an.create_group("MBSAcquisition_1")
        mbs["passenergy"] = np.bytes_(pe)
        mbs["lens_mode"] = np.bytes_(lens)
        mbs["center_ke"] = np.array([ef0 - 0.6])
        tc = an.create_group("i12-m-cx1-ex-tc.1")
        tc["temperature"] = np.array([temperature_K - 273.15])
    return {"ef0": ef0, "curvature_eV": curvature_eV, "gamma_deg": gamma_deg,
            "slit": slit, "defl": defl, "energy": energy}
