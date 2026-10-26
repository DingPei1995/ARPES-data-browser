"""
tools/moire.py
============
The moire Brillouin zone of two stacked 2-D lattices (independent lattice
constants, independent orientations), for the Brillouin-zone tool's moire
mode.

Two routes, kept side by side on purpose
------------------------------------------
:func:`moire_reciprocal_vectors` is the general method: the moire
superlattice's own reciprocal vectors are (to a very good approximation,
exact in the small-mismatch limit that is the entire point of a moire
pattern) the *differences* between nearby reciprocal-lattice points of the
two layers -- because a moire fringe is precisely where the two lattices'
Fourier components nearly cancel, i.e. where ``G_top - G_bottom`` is small.
Taking the two shortest independent such differences gives the moire
lattice's own primitive reciprocal vectors, whatever the two layers'
symmetries, lattice constants or relative rotation happen to be -- nothing
here assumes hexagonal layers or a shared lattice constant.

``moire_lattice.m`` only ever covers one case of this: two *hexagonal*
layers (its ``plotBZ_demo_old`` call has ``gamma`` hard-coded to 120
degrees), related by a closed-form formula for the moire lattice constant
and orientation. :func:`hex_moire_lattice_fast` is that formula, ported
unchanged, kept as a fast path for exactly that case (typed bilayers of the
same 2-D lattice type -- twisted bilayer graphene and its relatives are the
common real case) and, in ``test_moire.py``, as a closed-form cross-check of
the general method: for two equal hexagonal layers the two agree exactly
for a twist in [0, 30] degrees (see :func:`moire_reciprocal_vectors` for why
the closed form is only valid up to 30 degrees, and the general method is
not).
"""
from __future__ import annotations

import numpy as np

from tools.bz2d import reciprocal_vectors_2d, wigner_seitz_cell_2d

__all__ = ["hex_moire_lattice_fast", "moire_reciprocal_vectors", "moire_bz"]


def hex_moire_lattice_fast(a_top: float, a_bot: float, twist_deg: float):
    """Moire real-space lattice constant and orientation for two hexagonal
    layers, ported from ``moire_lattice.m``/``MiniBZ_plotter.m``'s
    ``Plot_MiniBZ_Callback``.

    Returns ``(a_moire, rotation_deg)``: the moire supercell's lattice
    constant (angstrom) and the angle (degrees) ``plotBZ_demo_old``'s own
    convention needs to draw it aligned with the two layers -- callers using
    the general method (:func:`moire_reciprocal_vectors`) do not need the
    second number, since that method returns actual vectors rather than an
    angle in one particular drawing convention.

    Restricted to two hexagonal (gamma = 120 deg) layers -- see the module
    docstring for why this is a fast path, not the general answer. It is
    also only exact for a twist in [0, 30] degrees: past 30 degrees the
    hexagonal point group makes a smaller moire cell available (twist theta
    gives the same pattern as 60 - theta) that this closed form does not
    check for, unlike :func:`moire_reciprocal_vectors` -- see that
    function's docstring. Equal-lattice callers past 30 degrees should fold
    the twist angle themselves (``min(twist, 60 - twist)``) before calling
    this, or just use the general method.

    A caveat specific to the ``a_top != a_bot`` branch: unlike the equal-
    lattice branch above (checked against :func:`moire_reciprocal_vectors`
    for every twist in ``test_moire.py``), this one has **not** been
    independently confirmed -- and does not, checked numerically while
    writing this, reduce to the textbook zero-twist limit ``a_top * a_bot /
    abs(a_top - a_bot)`` (the ordinary two-comb beat period), the way the
    general method does. It is ported unchanged from ``moire_lattice.m``
    for continuity with figures already made with that tool; for a new
    calculation, prefer :func:`moire_reciprocal_vectors`, which is verified
    in both limits.
    """
    top_rotation = 0.0
    bottom_rotation = float(twist_deg)
    if np.isclose(a_top, a_bot):
        theta = abs(top_rotation - bottom_rotation) / 180.0 * np.pi
        if abs(np.sin(theta / 2.0)) < 1e-12:
            raise ValueError("zero twist between two equal layers: no moire pattern")
        a_moire = a_top / (2.0 * np.sin(theta / 2.0))
        rotation = _mini_bz_rotation(top_rotation, bottom_rotation, a_top) + 60.0
    else:
        d = (a_bot - a_top) / a_bot
        theta = (bottom_rotation - top_rotation) / 180.0 * np.pi
        a_moire = (1.0 + d) * a_bot / np.sqrt(2.0 * (1.0 + d) * (1.0 - np.cos(theta) + d * d))
        rotation = (bottom_rotation
                   + np.degrees(np.arctan(np.sin(theta) / (1.0 + d - np.cos(theta))))
                   + 60.0)
    return float(a_moire), float(rotation)


def _mini_bz_rotation(top_rotation_deg, bottom_rotation_deg, lattice):
    top = np.deg2rad(top_rotation_deg)
    bot = np.deg2rad(bottom_rotation_deg)
    v_top = np.array([lattice * np.cos(top), lattice * np.sin(top)])
    v_bot = np.array([lattice * np.cos(bot), lattice * np.sin(bot)])
    v = v_bot - v_top
    return float(np.degrees(np.arccos(np.dot(v, [1.0, 0.0]) / np.linalg.norm(v))))


def moire_reciprocal_vectors(g1_top, g2_top, g1_bot, g2_bot, search: int = 1):
    """The moire superlattice's primitive reciprocal vectors, for two
    layers of *any* 2-D lattice type, lattice constant and relative
    rotation -- the two shortest linearly independent vectors among
    ``m1*g1_top + n1*g2_top - (m2*g1_bot + n2*g2_bot)`` for small integers.

    This is what a moire pattern *is*: wherever a top-layer reciprocal
    vector nearly coincides with a bottom-layer one, their difference is
    small, i.e. a long-wavelength beat -- the moire fringe. The shortest
    two independent differences among the *low-order* reciprocal vectors
    (``search`` shells; the default, 1, is each layer's own nearest star --
    +-g1, +-g2, +-(g1-g2) for a hexagonal layer) are the moire lattice's own
    reciprocal primitive vectors.

    Restricting to low order on purpose: a real periodic potential's Fourier
    weight falls off with |G|, so the moire pattern anyone would actually
    see is set by matching each layer's *strongest* components, not by
    whichever pair of high-order points happens to nearly coincide by
    numerical accident. A larger ``search`` does find shorter differences at
    some twist angles -- true near-commensurations between higher-order
    points -- but those describe a different, much weaker effect than "the
    moire pattern", so they are not what this function looks for by default.

    Confirmed against :func:`hex_moire_lattice_fast` (see ``test_moire.py``):
    for two equal hexagonal layers this matches the closed-form period
    exactly for a twist in [0, 30] degrees, and folds correctly onto it for
    (60 - twist) beyond that -- a real property of the hexagonal point group
    (rotating by theta is the same pattern as rotating by 60 - theta), which
    the closed-form formula itself does not know about and so only holds
    up to 30 degrees. This function needs no such per-symmetry knowledge:
    the fold falls out automatically from searching every low-order pairing
    rather than only matching each vector to its own index.
    """
    idx = np.arange(-search, search + 1)
    tops = np.array([m * np.asarray(g1_top) + n * np.asarray(g2_top)
                     for m in idx for n in idx])
    bots = np.array([m * np.asarray(g1_bot) + n * np.asarray(g2_bot)
                     for m in idx for n in idx])
    diffs = (tops[:, None, :] - bots[None, :, :]).reshape(-1, 2)
    norms = np.linalg.norm(diffs, axis=1)
    order = np.argsort(norms)
    diffs, norms = diffs[order], norms[order]

    nonzero = norms > 1e-9
    diffs, norms = diffs[nonzero], norms[nonzero]
    if diffs.size == 0:
        raise ValueError("the two layers are identical (no moire pattern)")

    # A real moire vector is the whole premise short: much shorter than the
    # parent lattices' own reciprocal vectors (that is what "long-wavelength
    # beat" means). If even the shortest candidate is not -- e.g. because
    # the two layers are identical or so close to it that every same-index
    # difference is exactly zero and only unrelated cross-index pairs
    # remain -- there is no well-defined moire pattern to report, and
    # returning a short-looking answer built from those unrelated pairs
    # would be worse than saying so.
    # For two hexagonal layers the legitimate ratio tops out at ~0.52 (at
    # the 30-degree fold boundary, see test_moire.py); a fully degenerate
    # (identical-layer) input's shortest surviving cross-term is a full
    # parent-vector length (ratio 1). 0.9 sits well inside the gap between
    # the two without ever tripping on a real, if extreme, twist angle.
    shortest_parent = min(np.linalg.norm(g1_top), np.linalg.norm(g2_top),
                          np.linalg.norm(g1_bot), np.linalg.norm(g2_bot))
    if norms[0] > 0.9 * shortest_parent:
        raise ValueError(
            "no small-difference (moire-forming) pair of reciprocal vectors "
            "found between the two layers -- they may be identical (no "
            "twist, no mismatch), which has no moire pattern")

    gm1 = diffs[0]
    gm1_hat = gm1 / norms[0]
    for d, n in zip(diffs[1:], norms[1:]):
        # independent: not (anti-)parallel to gm1, judged by the sine of the
        # angle between them, which is scale-free unlike a raw cross product
        cross = abs(gm1_hat[0] * d[1] - gm1_hat[1] * d[0])
        if cross > 1e-6 * n:
            return gm1, d
    raise RuntimeError(
        "could not find two independent moire reciprocal vectors within "
        f"the +-{search} search range -- the two layers may be exactly "
        "commensurate along one direction; try a larger `search`.")


def moire_bz(g1_top, g2_top, g1_bot, g2_bot, search: int = 1):
    """The moire Brillouin zone: :func:`moire_reciprocal_vectors` followed
    by :func:`tools.bz2d.wigner_seitz_cell_2d`. Returns ``(gm1, gm2,
    polygon)``.
    """
    gm1, gm2 = moire_reciprocal_vectors(g1_top, g2_top, g1_bot, g2_bot, search)
    polygon = wigner_seitz_cell_2d(gm1, gm2)
    return gm1, gm2, polygon
