"""
tools/kspace.py
=============
Angle-to-momentum conversion for ARPES Maps, ported from the beamline's
``k_space_conversion_demo.m`` (its 3-D branch -- the one that handles a
deflector Map; k_z is not covered here, and the single-Cut case
lives in ``tools.cutk``, which builds on the conventions below).

Physics, straight from the .m file
---------------------------------
A photoelectron leaving the sample along the angles (theta, phi) has an
in-plane momentum whose *direction* is the unit vector::

    kx0 = sin(theta) * cos(phi)
    ky0 = sin(phi)
    kz0 = cos(theta) * cos(phi)

(that really is a unit vector: cos^2(phi) * 1 + sin^2(phi) = 1), and whose
*length* depends only on the kinetic energy::

    |k| = sqrt(2 m_e E) / hbar

Because direction and length separate like that, a whole Map converts as one
direction field scaled per energy slice -- which is what makes the inversion
in :func:`convert_map` cheap.

Sample orientation enters as a rotation ``R = Rz @ Rx @ Ry`` built from the
theta offset, the phi offset and the sample azimuth, transcribed from the
.m file including its sign conventions (theta and the azimuth are negated
"just for user's habit"; phi is not).

Units: this module works in **inverse angstroms**. The .m file divided by
``k_unit = pi / lattice_constant`` to express momenta in units of pi/a, so
its lattice-constant box has no counterpart here -- ``K0_PER_SQRT_EV`` below
is the same constant with ``k_unit = 1e10 m^-1`` instead, i.e. the familiar
0.5123 * sqrt(E[eV]) A^-1.

How the resampling differs from the .m file, and why
---------------------------------------------------
The .m file goes forwards: it maps the measured angle grid into k, builds a
Delaunay triangulation of those scattered points, interpolates onto a
padded unit grid (``ceil(N*1.1)``), then ``interp2``-s once per energy.

Here the mapping is inverted instead. For a target ``(kx, ky)`` at energy
``E`` the emission angles that produced it are recovered in closed form::

    u = (kx, ky, sqrt(1 - kx^2 - ky^2)) / |k(E)|     # unit direction
    v = R^T u                                        # undo the sample rotation
    phi   = arcsin(v_y)
    theta = arctan2(v_x, v_z)

so each output point is read straight from the measured, *regular* angle
grid by ordinary bilinear interpolation. No triangulation, no padding
factor, and the result is the exact inverse of the forward map rather than a
piecewise-linear approximation of it. :func:`convert_map_forward_reference`
implements the .m file's forward-and-scatter route so the tests can check
the two against each other.

Points with ``kx^2 + ky^2 > |k|^2`` are outside the light cone and points
whose angles fall outside the measured range were never recorded; both come
back as NaN rather than as extrapolated values.
"""
from __future__ import annotations

import numpy as np

# sqrt(2 m_e e) * 2 pi / h, divided by 1e10 m^-1 per A^-1.
# The .m file's constants: ce=1.6021892e-19, me=9.109534e-31, h=6.626176e-34.
_CE = 1.6021892e-19
_ME = 9.109534e-31
_H = 6.626176e-34
K0_PER_SQRT_EV = np.sqrt(2.0 * _ME * _CE) * 2.0 * np.pi / _H / 1e10  # ~0.5123


def rotation_matrix(theta_offset_deg: float, phi_offset_deg: float,
                    azimuth_deg: float = 0.0) -> np.ndarray:
    """``R = Rz @ Rx @ Ry`` exactly as the .m file builds it, including its
    sign conventions: theta and the azimuth are negated ("negative sign is
    just for user's habit"), phi is not.

    These are the numbers as typed into the MATLAB dialog's boxes. They are
    **not** the angles of the point that lands at the k-space origin -- see
    :func:`origin_rotation`, which is what this program's GUI uses.

    Note MATLAB's ``[a', b', c']`` builds a matrix from *columns*, so e.g.
    its ``rx`` has ``-sin`` at (3,2) and ``+sin`` at (2,3) -- the transpose
    of the textbook form. That is reproduced here.
    """
    theta = -np.deg2rad(theta_offset_deg)
    phi = np.deg2rad(phi_offset_deg)
    azimuth = -np.deg2rad(azimuth_deg)

    rx = np.array([[1.0, 0.0, 0.0],
                   [0.0, np.cos(phi), np.sin(phi)],
                   [0.0, -np.sin(phi), np.cos(phi)]])
    ry = np.array([[np.cos(theta), 0.0, -np.sin(theta)],
                   [0.0, 1.0, 0.0],
                   [np.sin(theta), 0.0, np.cos(theta)]])
    rz = np.array([[np.cos(azimuth), np.sin(azimuth), 0.0],
                   [-np.sin(azimuth), np.cos(azimuth), 0.0],
                   [0.0, 0.0, 1.0]])
    return rz @ rx @ ry


def origin_rotation(theta_origin_deg: float, phi_origin_deg: float,
                    azimuth_deg: float = 0.0) -> np.ndarray:
    """Rotation that sends the emission direction ``(theta, phi)`` to normal
    emission, i.e. puts that point at the k-space origin.

    This is the sense the GUI works in: the user picks the point that should
    become k = (0, 0) on the constant-energy contour, and its own angles are
    what the boxes show.

    It differs from :func:`rotation_matrix` by a sign on both angles. The .m
    file's text boxes take the opposite convention, so its ``theta offset``
    and ``phi offset`` are the negatives of the numbers used here -- worth
    knowing when comparing against a MATLAB conversion of the same scan.
    Verified in the tests: the picked angles come back as k = (0, 0) to
    machine precision, which is not true of the raw .m convention.
    """
    return rotation_matrix(-theta_origin_deg, -phi_origin_deg, azimuth_deg)


def k0_of_energy(energy_eV) -> np.ndarray:
    """Momentum magnitude |k| in A^-1 for a kinetic energy in eV."""
    energy = np.asarray(energy_eV, dtype=float)
    return K0_PER_SQRT_EV * np.sqrt(np.clip(energy, 0.0, None))


def angle_grid_to_directions(theta_deg, phi_deg, R: np.ndarray):
    """Unit momentum directions for every point of the angle grid.

    Returns ``(kx1, ky1)``, each shaped ``(len(theta), len(phi))`` and
    dimensionless: multiply by ``k0_of_energy(E)`` for that energy's plane.
    """
    theta = np.deg2rad(np.asarray(theta_deg, dtype=float))[:, None]
    phi = np.deg2rad(np.asarray(phi_deg, dtype=float))[None, :]

    kx0 = np.sin(theta) * np.cos(phi)
    ky0 = np.broadcast_to(np.sin(phi), (theta.size, phi.size))
    kz0 = np.cos(theta) * np.cos(phi)

    kx1 = R[0, 0] * kx0 + R[0, 1] * ky0 + R[0, 2] * kz0
    ky1 = R[1, 0] * kx0 + R[1, 1] * ky0 + R[1, 2] * kz0
    return kx1, ky1


def _scaled_range(values: np.ndarray, k0_min: float, k0_max: float):
    """Range of ``values * k0`` over the whole energy span.

    The .m file takes its lower bound only from the lowest-energy plane and
    its upper bound only from the highest, which is wrong whenever the
    direction field spans zero: the most negative kx comes from the *highest*
    energy. Both endpoints are scanned here instead.
    """
    lo_candidates = (values.min() * k0_min, values.min() * k0_max)
    hi_candidates = (values.max() * k0_min, values.max() * k0_max)
    return float(min(lo_candidates)), float(max(hi_candidates))


def k_extent(theta_deg, phi_deg, energy_eV, *, theta_offset_deg: float = 0.0,
             phi_offset_deg: float = 0.0, azimuth_deg: float = 0.0,
             energy_offset_eV: float = 0.0):
    """The (kx, ky) box :func:`convert_map` would produce for these settings.

    Split out of the conversion so the dialog can answer "how many points is
    that?" before anything is converted -- asking for a resolution in A^-1
    only means something once the extent is known.
    """
    energy = np.asarray(energy_eV, dtype=float) + float(energy_offset_eV)
    R = origin_rotation(theta_offset_deg, phi_offset_deg, azimuth_deg)
    kx1, ky1 = angle_grid_to_directions(theta_deg, phi_deg, R)
    k0 = k0_of_energy(energy)
    positive = k0[k0 > 0]
    if positive.size == 0:
        raise ValueError("no positive kinetic energies to convert")
    k0_min, k0_max = float(positive.min()), float(positive.max())
    kx_lo, kx_hi = _scaled_range(kx1, k0_min, k0_max)
    ky_lo, ky_hi = _scaled_range(ky1, k0_min, k0_max)
    return kx_lo, kx_hi, ky_lo, ky_hi


def points_to_azimuth(points) -> float:
    """The smallest sample rotation that stands a picked direction upright.

    Give it two points on the constant-energy contour and it returns the
    rotation that makes the line through them vertical (along ky); give it
    one and the line is taken from the origin to that point, since that is
    the only line a single point defines.

    **Two points define a line, not an arrow.** Which of them was clicked
    first, and which side of the origin they sit on, are not things the
    sample knows about -- so the answer must not depend on either, and the
    result is folded to (-90, 90]. Standing a direction along +ky and along
    -ky are the same alignment of the same axis, and they differ by exactly
    the 180 degrees that folding removes; what is left is the smaller of the
    two turns, clockwise or anticlockwise.

    Before this was folded that far, clicking the same pair in the other
    order gave a rotation 180 degrees away -- a perfectly plausible number,
    and one that converts the map upside down.

    A high-symmetry direction picked off a Fermi surface is exactly what
    this is for: pick two points along it, and the converted map comes out
    with that direction up the ky axis.
    """
    pts = [(float(x), float(y)) for x, y in points]
    if not pts:
        raise ValueError("pick at least one point")
    if len(pts) == 1:
        dx, dy = pts[0]
    else:
        dx, dy = pts[1][0] - pts[0][0], pts[1][1] - pts[0][1]
    if abs(dx) < 1e-12 and abs(dy) < 1e-12:
        raise ValueError("the two points are the same -- they define no direction")
    # The sample rotation this program applies turns the data by -azimuth in
    # the (kx, ky) plane, so standing a direction at `angle` upright takes
    # 90 - angle. (Verified against a real conversion in test_kspace.)
    angle = np.degrees(np.arctan2(dy, dx))
    rotation = 90.0 - angle
    return _fold_to_right_angle(rotation)


def _fold_to_right_angle(rotation: float, tol: float = 1e-9) -> float:
    """Fold a rotation into (-90, 90], the range of an undirected line.

    Written the way it is so that a tie lands on +90 rather than -90: a
    quarter turn one way and a quarter turn the other are equally small, and
    picking the same one every time beats picking whichever the modulo
    happened to give.

    The tie needs the tolerance, not just the choice of interval. A
    horizontal line wants exactly a quarter turn, and ``arctan2`` returns it
    as 90 plus or minus a part in 1e14 depending on the sign of a sine that
    should have been zero -- which lands on opposite ends of a half-open
    interval. So a result within ``tol`` of the boundary is pulled to +90,
    and the same picked line gives the same number every time.
    """
    folded = -((-float(rotation) + 90.0) % 180.0 - 90.0)
    if abs(folded + 90.0) < tol or abs(folded - 90.0) < tol:
        return 90.0
    return float(folded) + 0.0      # normalise -0.0, which prints as "-0.000"


def _resample_energy(cube: np.ndarray, energy: np.ndarray, n_energy: int):
    """Linearly resample the cube along its energy axis."""
    if n_energy is None or n_energy == energy.size:
        return cube, energy
    if n_energy < 2:
        raise ValueError("n_energy must be at least 2")
    target = np.linspace(float(energy[0]), float(energy[-1]), int(n_energy))
    # interp needs an ascending sample axis
    order = np.argsort(energy)
    src_e = energy[order]
    src = cube[:, :, order]
    flat = src.reshape(-1, src.shape[2])
    out = np.empty((flat.shape[0], target.size), dtype=float)
    for i, row in enumerate(flat):
        out[i] = np.interp(target, src_e, row)
    return out.reshape(cube.shape[0], cube.shape[1], target.size), target


def convert_map(theta_deg, phi_deg, energy_eV, cube, *,
                theta_offset_deg: float = 0.0,
                phi_offset_deg: float = 0.0,
                azimuth_deg: float = 0.0,
                energy_offset_eV: float = 0.0,
                n_kx: int = 100, n_ky: int = 100,
                n_energy: int = None,
                progress=None):
    """Convert an angle-space Map cube to a regular (kx, ky, E) grid.

    ``cube`` is ``(theta, phi, E)`` -- for a deflector Map that is
    (deflector angle, angle along slit, energy), matching both the .m file's
    ``data.value`` layout and this program's own.

    ``theta_offset_deg`` / ``phi_offset_deg`` place the k-space origin: give
    them the angles of the point that should become (0, 0), which is what
    picking a point on the constant-energy contour does. (The .m file's text
    boxes use the opposite sign for these two -- see :func:`origin_rotation`.)

    ``n_kx`` / ``n_ky`` / ``n_energy`` set the output grid size, so the
    result can be coarser or finer than the measurement. ``n_energy=None``
    keeps the original energy sampling.

    Returns ``(kx, ky, energy, cube_k)`` with ``cube_k`` shaped
    ``(n_kx, n_ky, n_energy)`` and NaN wherever the target point was not
    measured.
    """
    from scipy.interpolate import RegularGridInterpolator

    theta_deg = np.asarray(theta_deg, dtype=float)
    phi_deg = np.asarray(phi_deg, dtype=float)
    energy = np.asarray(energy_eV, dtype=float) + float(energy_offset_eV)
    cube = np.asarray(cube, dtype=float)

    if cube.shape != (theta_deg.size, phi_deg.size, energy.size):
        raise ValueError(
            f"cube shape {cube.shape} does not match axes "
            f"(theta={theta_deg.size}, phi={phi_deg.size}, E={energy.size})")
    if n_kx < 2 or n_ky < 2:
        raise ValueError("n_kx and n_ky must be at least 2")

    cube, energy = _resample_energy(cube, energy, n_energy)

    R = origin_rotation(theta_offset_deg, phi_offset_deg, azimuth_deg)
    kx1, ky1 = angle_grid_to_directions(theta_deg, phi_deg, R)

    k0 = k0_of_energy(energy)
    positive = k0[k0 > 0]
    if positive.size == 0:
        raise ValueError("no positive kinetic energies to convert")
    k0_min, k0_max = float(positive.min()), float(positive.max())

    kx_lo, kx_hi = _scaled_range(kx1, k0_min, k0_max)
    ky_lo, ky_hi = _scaled_range(ky1, k0_min, k0_max)
    kx_grid = np.linspace(kx_lo, kx_hi, int(n_kx))
    ky_grid = np.linspace(ky_lo, ky_hi, int(n_ky))

    KX, KY = np.meshgrid(kx_grid, ky_grid, indexing="ij")
    Rt = R.T

    # The angle axes must ascend for the interpolator; remember the order so
    # the cube can be presented the same way.
    t_order = np.argsort(theta_deg)
    p_order = np.argsort(phi_deg)
    theta_sorted = theta_deg[t_order]
    phi_sorted = phi_deg[p_order]

    out = np.full((int(n_kx), int(n_ky), energy.size), np.nan, dtype=float)
    for m in range(energy.size):
        if k0[m] <= 0:
            continue                      # no momentum to speak of at E<=0
        kx_u = KX / k0[m]
        ky_u = KY / k0[m]
        radial = kx_u ** 2 + ky_u ** 2
        inside = radial <= 1.0            # outside is beyond the light cone
        if not inside.any():
            continue

        kz_u = np.sqrt(np.clip(1.0 - radial, 0.0, None))
        vx = Rt[0, 0] * kx_u + Rt[0, 1] * ky_u + Rt[0, 2] * kz_u
        vy = Rt[1, 0] * kx_u + Rt[1, 1] * ky_u + Rt[1, 2] * kz_u
        vz = Rt[2, 0] * kx_u + Rt[2, 1] * ky_u + Rt[2, 2] * kz_u

        phi_q = np.degrees(np.arcsin(np.clip(vy, -1.0, 1.0)))
        theta_q = np.degrees(np.arctan2(vx, vz))

        interp = RegularGridInterpolator(
            (theta_sorted, phi_sorted),
            cube[np.ix_(t_order, p_order, [m])][:, :, 0],
            bounds_error=False, fill_value=np.nan)
        points = np.stack([theta_q[inside], phi_q[inside]], axis=-1)
        plane = np.full(KX.shape, np.nan)
        plane[inside] = interp(points)
        out[:, :, m] = plane

        if progress is not None:
            progress(m + 1, energy.size)

    return kx_grid, ky_grid, energy, out


def convert_map_forward_reference(theta_deg, phi_deg, energy_eV, cube, *,
                                  theta_offset_deg: float = 0.0,
                                  phi_offset_deg: float = 0.0,
                                  azimuth_deg: float = 0.0,
                                  energy_offset_eV: float = 0.0,
                                  n_kx: int = 100, n_ky: int = 100):
    """The .m file's own route -- forward-map the angle grid, then scatter
    interpolate -- kept as a cross-check for :func:`convert_map`.

    Slower and blurrier (it is piecewise-linear over a triangulation of the
    mapped points rather than an exact inverse), so it is used by the tests
    and not by the application.
    """
    from scipy.interpolate import LinearNDInterpolator

    theta_deg = np.asarray(theta_deg, dtype=float)
    phi_deg = np.asarray(phi_deg, dtype=float)
    energy = np.asarray(energy_eV, dtype=float) + float(energy_offset_eV)
    cube = np.asarray(cube, dtype=float)

    R = origin_rotation(theta_offset_deg, phi_offset_deg, azimuth_deg)
    kx1, ky1 = angle_grid_to_directions(theta_deg, phi_deg, R)
    k0 = k0_of_energy(energy)
    k0_min, k0_max = float(k0.min()), float(k0.max())

    kx_lo, kx_hi = _scaled_range(kx1, k0_min, k0_max)
    ky_lo, ky_hi = _scaled_range(ky1, k0_min, k0_max)
    kx_grid = np.linspace(kx_lo, kx_hi, int(n_kx))
    ky_grid = np.linspace(ky_lo, ky_hi, int(n_ky))
    KX, KY = np.meshgrid(kx_grid, ky_grid, indexing="ij")

    pts = np.stack([kx1.ravel(), ky1.ravel()], axis=-1)
    out = np.full((int(n_kx), int(n_ky), energy.size), np.nan)
    for m in range(energy.size):
        if k0[m] <= 0:
            continue
        interp = LinearNDInterpolator(pts * k0[m], cube[:, :, m].ravel())
        out[:, :, m] = interp(KX, KY)
    return kx_grid, ky_grid, energy, out
