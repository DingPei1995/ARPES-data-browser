"""
tools/cleavage.py
=================
From a measured k_z period to which plane the crystal cleaved along.

A photon-energy scan repeats along ``k_z`` with the period of the reciprocal
lattice *in the direction of the surface normal*. Measure that period off the
converted map -- pick two points that look like the same feature one zone
apart -- and it names the cleavage plane, because for a surface ``(hkl)`` the
period is

    |G_hkl| = 2*pi / d_hkl

and no two low-index planes of a given crystal usually have the same
spacing. So: enumerate the reciprocal lattice, and report which vectors have
the length that was measured.

Two things make this less trivial than it sounds, and both are handled here:

* **Centering.** The period belongs to the *Bravais* lattice, not to the
  conventional cell. In a body-centred cubic crystal the (001) planes of the
  conventional cell are not all lattice planes -- the shortest reciprocal
  vector along [001] is (002), and the period is ``4*pi/a``, not ``2*pi/a``.
  Enumerating the primitive reciprocal lattice gets this right by
  construction; taking ``2*pi/c`` from the conventional cell does not.
* **Which zone the two points were.** A user picking "the same feature
  again" may well have picked two zones apart, especially on a scan that
  only covers a couple. Each candidate is therefore also tested against
  half and a third of the measured distance, and says so when that is how
  it matched.

The answer is a shortlist, not a verdict: a 15% window over a few unit cells
of reciprocal space normally leaves more than one plane standing, and which
one is right is settled by knowing the material.

Nothing here imports Qt.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from tools.lattice import LatticeParams, primitive_vectors, reciprocal_vectors

__all__ = ["Candidate", "candidates", "reciprocal_lengths"]


@dataclass
class Candidate:
    """One plane whose reciprocal-lattice period matches what was measured."""

    #: Miller indices against the **conventional** cell, which is how planes
    #: are named. Integer for any real lattice vector.
    hkl: tuple
    #: The period along that normal, A^-1 -- the length of the shortest
    #: reciprocal lattice vector in that direction.
    length: float
    #: The real-space repeat along the normal, A. ``2*pi/length``.
    spacing: float
    #: How many zones the picked points were apart for this to match.
    orders: int
    #: Signed relative error, ``(measured/orders - length) / length``.
    error: float

    @property
    def family(self) -> tuple:
        """``(hkl)`` reduced to coprime integers -- the plane's usual name.

        Different from :attr:`hkl` exactly when the first allowed reflection
        along that normal is not the first index, which is centering showing
        itself: ``(002)`` reduces to ``(001)``, and the distinction is the
        difference between the right period and half of it.
        """
        values = [int(v) for v in self.hkl]
        divisor = 0
        for value in values:
            divisor = np.gcd(divisor, abs(value))
        if divisor <= 1:
            return tuple(values)
        return tuple(v // divisor for v in values)

    def describe(self) -> str:
        name = "(" + " ".join(str(int(v)) for v in self.hkl) + ")"
        text = (f"{name}  period {self.length:.4f} A^-1  "
                f"d = {self.spacing:.4f} A  off by {100 * self.error:+.1f}%")
        if self.orders > 1:
            text += f"  [picked {self.orders} zones apart]"
        if self.family != tuple(int(v) for v in self.hkl):
            text += ("  [" + " ".join(str(v) for v in self.family)
                     + " family; the first allowed reflection is not the first index]")
        return text


def reciprocal_lengths(params: LatticeParams, max_index: int = 4,
                       direction_tol: float = 1e-6):
    """Every distinct normal direction of the lattice, with its period.

    Returns ``[(hkl_conventional, length, unit_direction), ...]`` sorted by
    length, one entry per direction: the *shortest* reciprocal lattice vector
    along each, which is what a k_z period actually measures. Longer vectors
    in the same direction are higher orders of the same plane and would only
    be duplicates in the shortlist.
    """
    primitive = primitive_vectors(params)
    b = reciprocal_vectors(primitive)
    conventional_a = _conventional_rows(params)

    shortest = {}
    span = range(-int(max_index), int(max_index) + 1)
    for n1 in span:
        for n2 in span:
            for n3 in span:
                if n1 == n2 == n3 == 0:
                    continue
                g = n1 * b[0] + n2 * b[1] + n3 * b[2]
                length = float(np.linalg.norm(g))
                if length < direction_tol:
                    continue
                # (hkl) and (-h-k-l) are the same set of planes, so the sign
                # is normalised before anything else -- otherwise every
                # direction appears twice and half of them print with a
                # leading minus for no reason.
                hkl = _round_indices(conventional_a @ g / (2 * np.pi))
                if _is_negative(hkl):
                    g, hkl = -g, tuple(-v for v in hkl)
                direction = g / length
                key = tuple(np.round(direction, 6) + 0.0)
                previous = shortest.get(key)
                if previous is None or length < previous[1] - 1e-9:
                    shortest[key] = (hkl, length, direction)

    found = sorted(shortest.values(), key=lambda item: item[1])
    return found


def _conventional_rows(params: LatticeParams) -> np.ndarray:
    from tools.lattice import conventional_vectors
    return conventional_vectors(params)


def _round_indices(values, tol: float = 1e-4) -> tuple:
    """Miller indices as plain integers.

    A reciprocal lattice vector of a centred lattice has integer indices
    against the conventional cell by construction -- that is what makes the
    conventional cell conventional -- so anything that comes back
    non-integer means the two bases have been mixed up, and saying so beats
    printing ``(0.5 0 1)``.
    """
    values = np.asarray(values, dtype=float)
    rounded = np.round(values)
    if np.any(np.abs(values - rounded) > tol):
        raise ValueError(
            f"reciprocal vector {values} is not an integer combination of "
            f"the conventional cell, which should not happen: the primitive "
            f"and conventional bases do not belong to the same lattice")
    return tuple(int(v) + 0 for v in rounded)


def _is_negative(hkl) -> bool:
    """True when the first non-zero index is negative."""
    for value in hkl:
        if value:
            return value < 0
    return False


def candidates(distance: float, params: LatticeParams, *,
               tolerance: float = 0.15, max_index: int = 4,
               max_orders: int = 3, limit: int = 12):
    """Which planes could give a k_z period of ``distance`` A^-1.

    ``tolerance`` is the fractional window, 0.15 by default. ``max_orders``
    allows for the two picked points being more than one zone apart.

    Sorted by how well they match, best first.
    """
    distance = float(distance)
    if not np.isfinite(distance) or distance <= 0:
        raise ValueError(
            f"a k_z period of {distance} is not something to match; pick two "
            f"points that are actually apart in k_z")

    found = reciprocal_lengths(params, max_index=max_index)
    if not found:
        raise ValueError("this lattice produced no reciprocal vectors to "
                         "match against")

    out = []
    for orders in range(1, int(max_orders) + 1):
        measured = distance / orders
        for hkl, length, _direction in found:
            error = (measured - length) / length
            if abs(error) <= tolerance:
                out.append(Candidate(hkl=hkl, length=length,
                                     spacing=2 * np.pi / length,
                                     orders=orders, error=float(error)))
    # A plane that matches at one order and again at another is one answer,
    # not two: keep whichever fitted best.
    best = {}
    for candidate in out:
        key = tuple(int(v) for v in candidate.hkl)
        if key not in best or abs(candidate.error) < abs(best[key].error):
            best[key] = candidate
    return sorted(best.values(), key=lambda c: abs(c.error))[:int(limit)]
