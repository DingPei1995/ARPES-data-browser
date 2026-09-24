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


# -- photon-energy scans --------------------------------------------------------
#: The planted crystal: body-centred tetragonal (I4/mmm, 139), cleaved (001).
#: Its k_z period is 4*pi/c (the (002) reflection), not 2*pi/c -- which is
#: the centering trap a batch conversion must not fall into.
BCT = {"a": 3.96, "b": 3.96, "c": 13.02, "alpha": 90, "beta": 90, "gamma": 90,
       "space_group": 139}


def _kz_band(k_par, k_z, c=BCT["c"]):
    """A band that disperses along k_z with the (002) period of the bct cell."""
    return -0.45 + 3.8 * k_par ** 2 - 0.25 * np.cos(k_z * c / 2.0)


def kz_intensity(hv, slit, energy, *, v0, work_function, angle_offset=0.0,
                 ef_offsets=None, curvature_eV=0.0, width=0.07, temperature_K=20.0,
                 bands=None):
    """(hv, slit, E) cube of a metal with the planted k_z dispersion.
    ``energy`` is the axis as stored; ``ef_offsets[i]`` is where spectrum i's
    true E_F sits on it (the misalignment kz_calibrate has to remove)."""
    from tools.kzconv import forward
    hv = np.asarray(hv, dtype=float)
    ef_offsets = np.zeros(hv.size) if ef_offsets is None else np.asarray(ef_offsets)
    alpha = np.radians(slit - angle_offset)
    bend = ef_curve(slit, 0.0, curvature_eV)                 # 0 at the middle
    out = np.empty((hv.size, slit.size, energy.size))
    kT = 8.617e-5 * max(temperature_K, 30.0) + 0.01
    for i, h in enumerate(hv):
        binding = energy[None, :] - ef_offsets[i] - bend[:, None]      # E - E_F
        kinetic = h - work_function + binding
        k_par, k_z = forward(kinetic, alpha[:, None] * np.ones_like(binding),
                             inner_potential=v0)
        spectral = 0.0
        for eps in (bands(k_par, k_z) if bands else [_kz_band(k_par, k_z)]):
            spectral = spectral + width ** 2 / ((binding - eps) ** 2 + width ** 2)
        fermi = 1.0 / (np.exp(binding / kT) + 1.0)
        out[i] = (spectral + 0.1) * fermi
    return out


def write_kz_native(path, *, hv=np.arange(30.0, 111.0, 1.0), n_slit=121, n_e=150,
                    slit_range=(-14.0, 14.0), v0=14.0, work_function=4.4,
                    curvature_eV=0.02, counts=60.0, drift_eV=0.04, seed=0,
                    lens="L4", pe="PE50", angle_offset=0.0, bands=None,
                    energy_range=(-1.3, 0.25)):
    """A kz map as the CASSIOPEE folder loader produces one: (hv, slit, E-E_F)
    with E_F placed from a tabulated work function -- so each spectrum's
    real edge is off by a smooth, unknown drift -- written in the viewer's own
    format (kind kz_map)."""
    from loader.nxs_file import save_dataset
    rng = np.random.default_rng(seed)
    slit = np.linspace(*slit_range, n_slit)
    energy = np.linspace(*energy_range, n_e)
    phase = rng.uniform(0, 2 * np.pi)
    drift = drift_eV * np.sin(np.linspace(0, 3, hv.size) + phase)
    cube = kz_intensity(hv, slit, energy, v0=v0, work_function=work_function,
                        angle_offset=angle_offset, ef_offsets=drift,
                        curvature_eV=curvature_eV, bands=bands)
    flux = np.linspace(1.5, 0.6, hv.size)[:, None, None]     # the beamline's share
    cube = rng.poisson(counts * cube * flux + 0.5).astype(np.float32)
    save_dataset(str(path), [{
        "name": "synthetic kz scan", "kind": "kz_map", "axes": (hv, slit, energy),
        "labels": {"x": "Photon energy (eV)", "k": "Angle along slit (°)",
                   "z": "E - E_F (eV)"},
        "value": cube,
        "info": {"axis0.role": "photon_energy", "cassiopee.work_function_eV": work_function,
                 "MBS.lens_mode": lens, "MBS.passenergy": pe,
                 "SampleTemperature_K": 20.0, "PhotonEnergy": float(hv[0])},
        "motors": {}}])
    return {"hv": hv, "slit": slit, "energy": energy, "drift": drift, "v0": v0,
            "work_function": work_function}


def write_kz_kinetic(path, *, hv=np.arange(60.0, 64.01, 0.5), n_slit=96,
                     phi=4.35, gamma_slit=0.0, counts=50.0, seed=0,
                     lens="L4", pe="PE50"):
    """A photon-energy scan of a semiconductor with its energy axis left in
    kinetic energy (CASSIOPEE's ``energy_reference="kinetic"``), written in
    the viewer's own format: bands below E_F only, no edge of its own -- so
    E_F has to come from the gold reference, as hv - W."""
    from loader.nxs_file import save_dataset
    rng = np.random.default_rng(seed)
    e_lo = hv.min() - phi - 1.2
    e_hi = hv.max() - phi + 0.3
    n_e = int(round((e_hi - e_lo) / 0.01)) + 1
    energy = np.linspace(e_lo, e_hi, n_e)
    slit = np.linspace(-15.0, 15.0, n_slit)
    cube = np.empty((hv.size, n_slit, n_e))
    for i, h in enumerate(hv):
        binding = energy[None, :] - (h - phi)
        k = K0 * np.sqrt(energy[None, :]) * np.sin(np.radians(slit[:, None] - gamma_slit))
        vbm = -0.8 - 4.0 * k ** 2 - 0.05 * np.cos(2 * np.pi * (h - hv[0]) / 3.0)
        cube[i] = 0.06 ** 2 / ((binding - vbm) ** 2 + 0.06 ** 2) + 0.05 * (binding < -0.9)
    cube = rng.poisson(counts * cube + 0.3).astype(np.float32)
    save_dataset(str(path), [{
        "name": "semiconductor kz scan", "kind": "kz_map", "axes": (hv, slit, energy),
        "labels": {"x": "Photon energy (eV)", "k": "Angle along slit (\u00b0)",
                   "z": "Kinetic energy (eV)"},
        "value": cube,
        "info": {"axis0.role": "photon_energy", "MBS.lens_mode": lens,
                 "MBS.passenergy": pe, "SampleTemperature_K": 20.0,
                 "PhotonEnergy": float(hv[0])},
        "motors": {}}])
    return {"hv": hv, "phi": phi, "energy": energy}


# -- a band calculation, and a scan of the same crystal ---------------------------
#: Primitive hexagonal, like the calculation in the user's example
#: (Gamma-M-K-Gamma-A-L-H-A, a = 6.20 A, c = 5.77 A).
HEX = {"a": 6.20, "b": 6.20, "c": 5.766, "alpha": 90, "beta": 90, "gamma": 120,
       "space_group": 191}


def hex_bands(k_par, k_z, c=HEX["c"]):
    """Three bands of a made-up hexagonal metal, with k_z dispersion of
    period 2*pi/c and different amplitudes, so Gamma and A differ."""
    phase = np.cos(k_z * c)
    return [-0.25 + 3.0 * k_par ** 2 - 0.18 * phase,
            -0.80 - 2.0 * k_par ** 2 + 0.12 * phase,
            -1.40 - 1.0 * k_par ** 2 + 0.30 * phase]


def write_band_file(path, *, bands=hex_bands, lattice=HEX,
                    labels=("G", "M", "K", "G", "A", "L", "H", "A"), points=51,
                    extra_bands=True):
    """A band file as Wannier90 writes one: k distance (1/A, 2*pi included)
    and energy, one block per band, the vertices repeated. With
    ``extra_bands``, bands far from E_F pad it out as a real one would."""
    from arpes_batch.calcbands import POINTS
    from tools.lattice import LatticeParams, conventional_vectors, reciprocal_vectors
    b = reciprocal_vectors(conventional_vectors(LatticeParams(**lattice)))
    coords = [np.asarray(POINTS["hexagonal"][name]) @ b for name in labels]
    ks, dist, total = [], [], 0.0
    for start, end in zip(coords, coords[1:]):
        for t in np.linspace(0, 1, points):
            k = start + t * (end - start)
            if ks:
                total += float(np.linalg.norm(k - ks[-1]))
            ks.append(k)
            dist.append(total)
    ks = np.array(ks)
    k_par = np.hypot(ks[:, 0], ks[:, 1])
    values = list(bands(k_par, ks[:, 2]))
    if extra_bands:
        values = [np.full(k_par.shape, -8.0) + 0.5 * np.cos(ks[:, 2] * 3)] + values + [
            2.0 + 1.5 * k_par ** 2]
    with open(path, "w") as fh:
        for band in values:
            for d, e in zip(dist, band):
                fh.write(f"  {d:.8f}  {e:.7f}\n")
            fh.write("\n")
    return {"labels": labels, "n_bands": len(values)}
