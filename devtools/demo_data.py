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


# --------------------------------------------------------------------------
# A layered crystal in three dimensions, and a "calculation" of it
# --------------------------------------------------------------------------
# The bands above, made a Hamiltonian: four orbitals on the triangular
# lattice (a = 3.3 A) stacked with c = 6 A, each with its own in-plane and
# out-of-plane hopping, and two of them mixed so their bands anticross.
# The measurement is this Hamiltonian; the "calculation" written next to it
# is the same one before renormalisation -- its bands are 1/CALC_Z wider and
# sit CALC_SHIFT higher -- which is what fitting it to the data recovers.
LATTICE_C = 6.0
INNER_POTENTIAL = 12.0
CALC_Z = 0.8                 # measured = CALC_Z * calculated + CALC_SHIFT
CALC_SHIFT = -0.05
ORBITALS = ("Ti-dz2", "Se-pxy", "Ti-dxy", "Se-pz")
#: (onsite at Gamma, t (eV, per unit of h), k_z amplitude) of each orbital:
#: E = e0 + t h(k) + kz_amp cos(k_z c), h = 6 - f(k) as in _hexagonal
#: orbitals 1 and 2 are not mixed, so their bands cross: a line of Dirac
#: points through the zone, whose energy follows orbital 2's k_z term
_ORBITAL_BANDS = ((-0.75, 0.120, 0.04), (-0.25, -0.100, -0.06),
                  (0.25 - 0.15 * 9.0, 0.150, -0.08), (-1.20, 0.012, 0.03))
_MIXING = 0.03               # orbital 0 <-> 1, eV, on site


def lattice_vectors():
    a, c = LATTICE_A, LATTICE_C
    return np.array([[a, 0.0, 0.0], [-a / 2, a * np.sqrt(3) / 2, 0.0],
                     [0.0, 0.0, c]])


def hoppings(scale: float = 1.0, shift: float = 0.0):
    """``{R: H(R)}`` (eV, R in lattice units) of the model, with every
    energy divided by ``scale`` after subtracting ``shift`` -- (1, 0) is
    the measured crystal, (CALC_Z, CALC_SHIFT) the calculation."""
    n = len(_ORBITAL_BANDS)
    out = {}

    def add(R, i, j, value):
        out.setdefault(R, np.zeros((n, n)))[i, j] += value / scale

    for o, (e0, t, kz_amp) in enumerate(_ORBITAL_BANDS):
        add((0, 0, 0), o, o, e0 + 6.0 * t - shift)
        for R in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (1, 1, 0),
                  (-1, -1, 0)):
            add(R, o, o, -t)
        if kz_amp:
            add((0, 0, 1), o, o, kz_amp / 2)
            add((0, 0, -1), o, o, kz_amp / 2)
    add((0, 0, 0), 0, 1, _MIXING)
    add((0, 0, 0), 1, 0, _MIXING)
    return out


def bulk_model(scale: float = 1.0, shift: float = 0.0):
    from tools.crystal.bands import WannierModel

    table = hoppings(scale, shift)
    R = np.array(list(table), dtype=float)
    H = np.array([table[key] for key in table], dtype=complex)
    return WannierModel(lattice_vectors(), R, H, orbital_names=list(ORBITALS),
                        name="demo crystal")


def hv_for_kz(kz: float) -> float:
    """The photon energy whose final state at normal emission and E_F has
    this k_z (free electron, V0 = INNER_POTENTIAL)."""
    return kz ** 2 / A_KZ - INNER_POTENTIAL + WORK_FUNCTION


#: the Gamma plane (k_z = 8 pi / c) and the zone-top plane (9 pi / c)
HV_GAMMA = round(hv_for_kz(8 * np.pi / LATTICE_C), 1)
HV_TOP = round(hv_for_kz(9 * np.pi / LATTICE_C), 1)


def _bulk_spectral(model, k, e, polarisation, resolution=0.012):
    """I(..., E) of the Hamiltonian at k (..., 3), each band weighted by
    its orbital character times a polarisation's matrix element."""
    shape = k.shape[:-1]
    energies, vectors = model.diagonalise(k.reshape(-1, 3))
    character = np.abs(vectors) ** 2                  # (N, orbital, band)
    element = {"LH": np.array([1.0, 0.35, 0.9, 0.5]),
               "LV": np.array([0.3, 1.0, 0.5, 0.45])}[polarisation]
    weight = np.einsum("o,nob->bn", element, character)
    bands = [(energies[:, b].reshape(shape), weight[b].reshape(shape))
             for b in range(energies.shape[1])]
    return _spectral(e, bands, resolution=resolution)


def _sphere_k(kx, ky, hv):
    e_kin = hv - WORK_FUNCTION
    kz = np.sqrt(np.clip(A_KZ * (e_kin + INNER_POTENTIAL) - kx ** 2 - ky ** 2,
                         0, None))
    return np.stack([kx, ky, kz], axis=-1)


def bulk_map(rng, hv, polarisation="LH"):
    """A deflector map at photon energy ``hv`` on the 3-D crystal: the k_z
    it reaches follows the free-electron sphere."""
    model = bulk_model()
    e_f = hv - WORK_FUNCTION
    defl = np.linspace(-22, 22, 89)
    slit = np.linspace(-22, 22, 177)
    ekin = np.linspace(e_f - 1.4, e_f + 0.12, 153)
    D, S = np.meshgrid(np.radians(defl), np.radians(slit), indexing="ij")
    k0 = A_K * np.sqrt(e_f)
    k = _sphere_k(k0 * np.sin(D), k0 * np.cos(D) * np.sin(S), hv)
    cube = _counts(_bulk_spectral(model, k, ekin - e_f, polarisation), 150.0, rng)
    return {"kind": "map", "axes": (defl, slit, ekin), "value": cube,
            "labels": {"x": "Angle, deflector (deg)", "k": "Angle, slit (deg)",
                       "z": "Kinetic energy (eV)"},
            "info": _info(photon_energy_eV=hv,
                          **{"axis0.role": "angle",
                             "cassiopee.polarisation": polarisation})}


def bulk_cut(rng, hv, direction_deg, polarisation="LH"):
    """A cut through Gamma at ``hv``, the slit along ``direction_deg``
    from k_x (0: Gamma-K, 90: Gamma-M)."""
    model = bulk_model()
    e_f = hv - WORK_FUNCTION
    slit = np.linspace(-20, 20, 361)
    ekin = np.linspace(e_f - 1.4, e_f + 0.12, 229)
    kpar = A_K * np.sqrt(e_f) * np.sin(np.radians(slit))
    phi = np.radians(direction_deg)
    k = _sphere_k(kpar * np.cos(phi), kpar * np.sin(phi), hv)
    cut = _counts(_bulk_spectral(model, k, ekin - e_f, polarisation,
                                 resolution=0.010), 1500.0, rng)
    return {"kind": "cut", "axes": (slit, ekin), "value": cut,
            "labels": {"x": "Angle, slit (deg)", "y": "Kinetic energy (eV)"},
            "info": _info(photon_energy_eV=hv,
                          **{"cassiopee.polarisation": polarisation,
                             "demo.cut_direction_deg": direction_deg})}


def bulk_hv_scan(rng):
    """A photon-energy scan along Gamma-M through a whole period of k_z."""
    model = bulk_model()
    hv = np.arange(50.0, 92.01, 1.0)
    slit = np.linspace(-12, 12, 97)
    e = np.linspace(-1.4, 0.12, 153)
    H, S = np.meshgrid(hv, np.radians(slit), indexing="ij")
    kpar = A_K * np.sqrt(H - WORK_FUNCTION) * np.sin(S)
    kz = np.sqrt(np.clip(A_KZ * (H - WORK_FUNCTION + INNER_POTENTIAL) - kpar ** 2,
                         0, None))
    k = np.stack([0 * kpar, kpar, kz], axis=-1)
    flux = (1.0 + 1.5 * np.exp(-((hv - 65) / 25) ** 2))[:, None, None]
    # LV: both branches of the line node bright (in LH one is faint)
    cube = _counts(_bulk_spectral(model, k, e, "LV") * flux, 110.0, rng)
    return {"kind": "kz_map", "axes": (hv, slit, e), "value": cube,
            "labels": {"x": "Photon energy (eV)", "k": "Angle, slit (deg)",
                       "z": "E - E_F (eV)"},
            "info": _info(**{"axis0.role": "photon_energy",
                             "photon_energy_eV": float(hv.mean()),
                             "cassiopee.polarisation": "LV"})}


#: X-ray lines (keV) of an EDS spectrum of the crystal, and their heights
EDS_LINES = (("C Kα", 0.277, 0.05), ("O Kα", 0.525, 0.08),
             ("Se Lα", 1.379, 1.0), ("Ti Kα", 4.511, 0.62),
             ("Ti Kβ", 4.932, 0.09), ("Se Kα", 11.222, 0.30),
             ("Se Kβ", 12.496, 0.05))


def eds_spectrum(rng):
    """An energy-dispersive X-ray spectrum: the lines on a bremsstrahlung
    background, the detector's width growing with energy."""
    energy = np.linspace(0.1, 14.0, 1391)
    background = 30.0 * energy * np.exp(-energy / 3.0) * (1 - np.exp(-energy / 0.4))
    counts = background.copy()
    for _name, centre, height in EDS_LINES:
        width = np.sqrt(0.055 ** 2 + 0.0025 * centre) / 2.355    # FWHM -> sigma
        counts += 2400.0 * height * np.exp(-0.5 * ((energy - centre) / width) ** 2)
    counts = rng.poisson(counts).astype(float)
    info = _info(**{"curve.value_label": "Counts",
                    "curve.channels": "EDS", "demo": "synthetic EDS spectrum"})
    return {"kind": "edc", "axes": (energy, np.zeros(1)), "value": counts[:, None],
            "labels": {"x": "Energy (keV)", "y": "Channel"}, "info": info}


CIF = """data_demo_1T
_chemical_name_common 'demo 1T-TiSe2-like crystal'
_cell_length_a {a}
_cell_length_b {a}
_cell_length_c {c}
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 120
_symmetry_space_group_name_H-M 'P -3 m 1'
_space_group_IT_number 164
loop_
_symmetry_equiv_pos_as_xyz
'x, y, z'
'-y, x-y, z'
'-x+y, -x, z'
'y, x, -z'
'x-y, -y, -z'
'-x, -x+y, -z'
'-x, -y, -z'
'y, -x+y, -z'
'x-y, x, -z'
'-y, -x, z'
'-x+y, y, z'
'x, x-y, z'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Ti1 Ti 0 0 0
Se1 Se 0.33333 0.66667 0.255
"""


#: A layered tetragonal crystal of the ZrSiS type (P4/nmm, origin choice
#: 2), for the structure figure: square Si nets between Hf-S slabs, which
#: cleave between the S layers.
TETRAGONAL_CIF = """data_demo_HfSiS_type
_chemical_name_common 'demo HfSiS-type crystal'
_cell_length_a 3.52
_cell_length_b 3.52
_cell_length_c 8.00
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_symmetry_space_group_name_H-M 'P 4/n m m'
_space_group_IT_number 129
loop_
_symmetry_equiv_pos_as_xyz
'x, y, z'
'-x+1/2, -y+1/2, z'
'-y+1/2, x, z'
'y, -x+1/2, z'
'-x, y+1/2, -z'
'x+1/2, -y, -z'
'y+1/2, x+1/2, -z'
'-y, -x, -z'
'-x, -y, -z'
'x+1/2, y+1/2, -z'
'y+1/2, -x, -z'
'-y, x+1/2, -z'
'x, -y+1/2, z'
'-x+1/2, y, z'
'-y+1/2, -x+1/2, z'
'y, x, z'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Hf1 Hf 0.25 0.25 0.255
Si1 Si 0.75 0.25 0.0
S1 S 0.25 0.25 0.62
"""


def write_calculation(folder: str) -> "list[str]":
    """The "calculation" of the crystal, as a theory group would send it:
    a Wannier90 Hamiltonian (``demo_calc_hr.dat`` + ``demo_calc.win``) and
    the same bands on a grid (``demo_calc.bxsf``); and the crystal's CIF."""
    from tools.crystal import bandfiles

    os.makedirs(folder, exist_ok=True)
    hr = os.path.join(folder, "demo_calc_hr.dat")
    bandfiles.write_hr(hr, lattice_vectors(), hoppings(CALC_Z, CALC_SHIFT),
                       len(ORBITALS), fermi=0.0, orbital_names=list(ORBITALS))
    bxsf = os.path.join(folder, "demo_calc.bxsf")
    bandfiles.write_bxsf(bxsf, bulk_model(CALC_Z, CALC_SHIFT), n=(36, 36, 16),
                         fermi=0.0)
    cif = os.path.join(folder, "demo_crystal.cif")
    with open(cif, "w", encoding="utf-8") as handle:
        handle.write(CIF.format(a=LATTICE_A, c=LATTICE_C))
    tetragonal = os.path.join(folder, "demo_tetragonal_crystal.cif")
    with open(tetragonal, "w", encoding="utf-8") as handle:
        handle.write(TETRAGONAL_CIF)
    return [hr, os.path.join(folder, "demo_calc.win"), bxsf, cif, tetragonal]


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
        "demo_bulk_bands.nxs": [
            (f"FS_map_hv{HV_GAMMA:g}_Gamma_plane", bulk_map(rng, HV_GAMMA)),
            (f"FS_map_hv{HV_TOP:g}_A_plane", bulk_map(rng, HV_TOP)),
            ("cut_GM_LH", bulk_cut(rng, HV_GAMMA, 90.0, "LH")),
            ("cut_GM_LV", bulk_cut(rng, HV_GAMMA, 90.0, "LV")),
            ("cut_GK_LH", bulk_cut(rng, HV_GAMMA, 0.0, "LH")),
            ("hv_scan_GM_50-92eV", bulk_hv_scan(rng))],
        "demo_eds.nxs": [("EDS_spectrum", eds_spectrum(rng))],
    }
    written = []
    for filename, entries in files.items():
        path = os.path.join(folder, filename)
        save_dataset(path, [dict(item, name=name) for name, item in entries])
        written.append(path)
    written += write_calculation(folder)
    return written


if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "demo_data")
    for path in write_all(out):
        print(path)
