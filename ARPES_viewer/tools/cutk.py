"""
tools/cutk.py
===========
Angle-to-momentum conversion for a single **Cut** -- one E-vs-angle
spectrum taken at one deflector position -- as the companion to
``tools.kspace``'s Map conversion, and in the same conventions.

Why this is not just the Map conversion with one row
----------------------------------------------------
A Map is a surface in the (kx, ky) plane, so its conversion has somewhere
to put both momentum components. A Cut is a *line* in that plane, and the
line generally does not pass through Gamma: it is a chord at some
perpendicular distance from it. That distance cannot be recovered from the
cut's own intensities -- a line carries no information about where the line
is -- so it has to come from the geometry: the angles at which Gamma sits.

The lab's MATLAB ``k_space_conversion_cut.m`` asks instead for a "centre
line" clicked on the cut, and feeds it in as ``sin(alpha + offset)``. That
conflates two different things:

* **where normal emission is** -- an instrument property, which enters
  *inside* the sine and therefore sets the k axis's **scale**;
* **where Gamma is** -- a sample property, which is a rigid translation in
  the (kx, ky) plane applied *after* the sine.

Clicking the middle of a band and calling it the offset gets the scale
wrong, not just the origin. With a cut over +-15 deg at 95 eV kinetic the
true kx range is +-1.292 A^-1; with a centre line picked 10 deg off it
comes out as -2.110 .. +0.435 A^-1 -- stretched on one side, squashed on
the other. So here the two are separate: the geometry (deflector angle,
Gamma's angles, azimuth) fixes the whole trajectory, and where the axis is
zeroed is a display choice made afterwards.

What the cut looks like in k space
----------------------------------
Following ``tools.kspace``'s convention -- theta is the deflector angle, phi
the angle along the slit, and a direction is
``(sin(theta) cos(phi), sin(phi), cos(theta) cos(phi))`` -- a cut at fixed
deflector traces an arc as the slit angle runs. Projected on the cut's own
direction it gives ``k_parallel``; the component across it is ``k_perp``,
the distance from Gamma. ``k_perp`` is not quite constant: the arc bows by
about 0.02 A^-1 over +-15 deg at 95 eV (comparable to a 0.2 deg angular
resolution), and it also scales with sqrt(E_kin), so it drifts by
``k_perp * dE / (2 E_kin)`` across the energy window. Both are reported by
:func:`convert_cut` rather than hidden, since they are what tells you
whether "the cut is a straight chord" is good enough for the cut in hand.

Resampling
----------
``k_parallel`` is a monotonic function of the slit angle, so the conversion
is a **one-dimensional interpolation per energy row** -- exact, cheap, and
free of the diagonal smearing a scattered/Delaunay interpolation puts into
a Fermi edge. (The MATLAB file builds an ``ii*jj`` point list and calls
``griddata``, which mixes neighbouring energies into every output point.)
"""
from __future__ import annotations

import numpy as np

from tools.kspace import (K0_PER_SQRT_EV, angle_grid_to_directions,
                        k0_of_energy, origin_rotation)

__all__ = ["cut_directions", "cut_momenta", "cut_extent", "convert_cut",
           "K_LABEL", "K_RADIAL_LABEL"]

K_LABEL = "k∥ (Å⁻¹)"
K_RADIAL_LABEL = "k from Γ (Å⁻¹)"


# --------------------------------------------------------------------------
# Geometry
# --------------------------------------------------------------------------
def cut_directions(slit_deg, *, deflector_deg: float,
                   gamma_deflector_deg: float = 0.0,
                   gamma_slit_deg: float = 0.0, azimuth_deg: float = 0.0):
    """Unit momentum directions along the cut, with Gamma at the origin.

    Returns ``(along, across, unit)``: the component of each slit angle's
    direction along the cut, the component across it, and the 2-vector the
    cut runs along. All dimensionless -- multiply by ``k0_of_energy(E)`` for
    that energy's momenta.

    ``gamma_deflector_deg`` / ``gamma_slit_deg`` are the analyser angles at
    which Gamma was observed. They are exactly the numbers the Map
    conversion calls its theta/phi offsets, which is what makes a converted
    Map's settings transferable to a cut taken in the same alignment.
    """
    slit = np.asarray(slit_deg, dtype=float)
    if slit.size < 2:
        raise ValueError("a cut needs at least two slit angles")
    R = origin_rotation(gamma_deflector_deg, gamma_slit_deg, azimuth_deg)
    # theta is the (fixed) deflector angle, phi the slit axis being scanned
    kx1, ky1 = angle_grid_to_directions([float(deflector_deg)], slit, R)
    kx1, ky1 = kx1[0], ky1[0]

    step = np.array([kx1[-1] - kx1[0], ky1[-1] - ky1[0]], dtype=float)
    length = float(np.hypot(*step))
    if length <= 0:
        raise ValueError("the cut's angles do not span a direction in k space")
    unit = step / length
    across_unit = np.array([-unit[1], unit[0]])

    along = kx1 * unit[0] + ky1 * unit[1]
    across = kx1 * across_unit[0] + ky1 * across_unit[1]
    return along, across, unit


def cut_momenta(slit_deg, energy_eV, *, deflector_deg: float,
                gamma_deflector_deg: float = 0.0, gamma_slit_deg: float = 0.0,
                azimuth_deg: float = 0.0, energy_offset_eV: float = 0.0):
    """``(k_par, k_perp)``, each ``(n_slit, n_energy)`` in A^-1.

    ``k_par`` is measured from Gamma's projection onto the cut -- it is the
    component of the momentum along the cut, and momenta are measured from
    Gamma, so its zero is the point of the cut closest to Gamma with no
    extra bookkeeping.
    """
    along, across, _unit = cut_directions(
        slit_deg, deflector_deg=deflector_deg,
        gamma_deflector_deg=gamma_deflector_deg,
        gamma_slit_deg=gamma_slit_deg, azimuth_deg=azimuth_deg)
    energy = np.asarray(energy_eV, dtype=float) + float(energy_offset_eV)
    k0 = k0_of_energy(energy)
    return along[:, None] * k0[None, :], across[:, None] * k0[None, :]


def _signed_radius(k_par, k_perp):
    """Distance from Gamma, signed by which side of Gamma's projection the
    point is on. The cut never comes closer to Gamma than ``|k_perp|``, so
    this axis has a real gap of that size around zero -- which is the point
    of offering it."""
    return np.sign(k_par) * np.hypot(k_par, k_perp)


def cut_extent(slit_deg, energy_eV, *, radial: bool = False, trim: bool = False,
               **geometry):
    """The k range :func:`convert_cut` would produce, so a resolution in
    A^-1 can be turned into a point count before anything is converted."""
    k_par, k_perp = cut_momenta(slit_deg, energy_eV, **geometry)
    values = _signed_radius(k_par, k_perp) if radial else k_par
    if trim:
        # the span every energy row actually covers
        lo = float(np.max(np.min(values, axis=0)))
        hi = float(np.min(np.max(values, axis=0)))
        if hi > lo:
            return lo, hi
    return float(values.min()), float(values.max())


# --------------------------------------------------------------------------
# Conversion
# --------------------------------------------------------------------------
def _resample_energy(values: np.ndarray, energy: np.ndarray, n_energy):
    """Linearly resample ``(n_slit, n_E)`` along its energy axis."""
    if n_energy is None or int(n_energy) == energy.size:
        return values, energy
    if int(n_energy) < 2:
        raise ValueError("n_energy must be at least 2")
    target = np.linspace(float(energy[0]), float(energy[-1]), int(n_energy))
    order = np.argsort(energy)
    src_e, src = energy[order], values[:, order]
    out = np.empty((values.shape[0], target.size), dtype=float)
    for i in range(values.shape[0]):
        out[i] = np.interp(target, src_e, src[i])
    return out, target


def convert_cut(slit_deg, energy_eV, value, *, deflector_deg: float,
                gamma_deflector_deg: float = 0.0, gamma_slit_deg: float = 0.0,
                azimuth_deg: float = 0.0, energy_offset_eV: float = 0.0,
                n_k: int = 200, n_energy: int = None, radial: bool = False,
                trim: bool = False):
    """Convert one Cut from (angle, E) to (k, E).

    ``value`` is ``(n_slit, n_energy)``, matching this program's Cut layout.
    Returns ``(k_axis, energy_axis, values, report)``.

    ``radial=False`` (the default) gives the momentum **along the cut**,
    zeroed at Gamma's projection onto it -- the coordinate a band disperses
    in, and the one to compare with a calculation along that line.
    ``radial=True`` gives the distance **from Gamma** instead, which is
    honest about the cut never reaching it: the axis has a gap of
    ``+-|k_perp|`` that comes back as missing rather than interpolated
    across.

    ``trim=True`` cuts the output down to the k range every energy row
    covers, instead of the union (where the rows that do not reach the ends
    leave NaN wedges, as the MATLAB version silently does).

    ``report`` carries what the conversion cannot put on the axes: the
    cut's distance from Gamma, how much that varies along the cut and
    across the energy window, and how much of the output was actually
    measured.
    """
    slit = np.asarray(slit_deg, dtype=float)
    energy = np.asarray(energy_eV, dtype=float) + float(energy_offset_eV)
    values = np.asarray(value, dtype=float)
    if values.shape != (slit.size, energy.size):
        raise ValueError(
            f"cut shape {values.shape} does not match its axes "
            f"(slit={slit.size}, E={energy.size})")
    if int(n_k) < 2:
        raise ValueError("n_k must be at least 2")

    values, energy = _resample_energy(values, energy, n_energy)

    geometry = dict(deflector_deg=deflector_deg,
                    gamma_deflector_deg=gamma_deflector_deg,
                    gamma_slit_deg=gamma_slit_deg, azimuth_deg=azimuth_deg)
    along, across, unit = cut_directions(slit, **geometry)
    k0 = k0_of_energy(energy)
    if not np.any(k0 > 0):
        raise ValueError("no positive kinetic energies to convert")

    k_par = along[:, None] * k0[None, :]
    k_perp = across[:, None] * k0[None, :]
    source = _signed_radius(k_par, k_perp) if radial else k_par

    # k_parallel is monotonic in the slit angle for any sane angular range;
    # if it ever is not, a 1-D interpolation per energy row is the wrong
    # tool and silently producing nonsense would be worse than stopping.
    direction = np.sign(source[-1, 0] - source[0, 0])
    if direction == 0 or not np.all(np.diff(source[:, 0]) * direction > 0):
        raise ValueError("the cut's k axis is not monotonic in the slit "
                         "angle; check the deflector and Gamma angles")

    lo, hi = cut_extent(slit, energy, radial=radial, trim=trim,
                        energy_offset_eV=0.0, **geometry)
    if not (hi > lo):
        raise ValueError("the cut covers no k range with these settings")
    k_axis = np.linspace(lo, hi, int(n_k))

    # A target point sitting exactly on the source's last sample must not be
    # lost to a rounding error -- with the grid built from that very value,
    # a bare left/right=NaN throws away the edge column of every export.
    tolerance = 1e-9 * max(abs(lo), abs(hi), 1.0)

    out = np.full((k_axis.size, energy.size), np.nan, dtype=float)
    for j in range(energy.size):
        column = source[:, j]
        row = values[:, j]
        if column[0] > column[-1]:
            column, row = column[::-1], row[::-1]
        out[:, j] = np.interp(k_axis, column, row)
        out[(k_axis < column[0] - tolerance)
            | (k_axis > column[-1] + tolerance), j] = np.nan
        if radial:
            # Inside |k| < the row's closest approach there is no data at
            # all, and np.interp would happily draw a line across the gap.
            gap = float(np.min(np.abs(k_perp[:, j])))
            out[np.abs(k_axis) < gap, j] = np.nan

    perp = np.abs(k_perp)
    # How far the arc departs from the straight chord through its ends --
    # the error in calling the cut "a line at constant k_perp" -- taken at
    # the energy where it is largest.
    span = along[-1] - along[0]
    chord = across[0] + (across[-1] - across[0]) * (along - along[0]) / span
    report = {
        "k_perp": float(np.mean(k_perp)),
        "k_perp_min": float(perp.min()),
        "k_perp_max": float(perp.max()),
        "bow": float(np.max(np.abs(across - chord)) * float(k0.max())),
        # the sqrt(E) scaling of the perpendicular offset across the window,
        # measured at the middle of the cut so the bow does not enter it
        "energy_drift": float(abs(across[across.size // 2])
                              * (float(k0.max()) - float(k0.min()))),
        "direction": (float(unit[0]), float(unit[1])),
        "k_step": float(k_axis[1] - k_axis[0]),
        "measured_fraction": float(np.isfinite(out).mean()),
    }
    return k_axis, energy, out, report
