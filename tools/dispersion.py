"""
tools/dispersion.py
=================
What comes out of a fitted band: the dispersion, the Fermi velocity, the
effective mass, and the self-energy. No Qt.

Three things here that the lab's MATLAB tools get wrong
-------------------------------------------------------
**The errors are on the fitted axis.** An MDC fit measures a *momentum* at a
known energy, so its uncertainty is on k; an EDC fit measures an energy at a
known momentum. ``EmassFitting.m`` and ``fit_MDC_tools.m`` both call
``fit(x, y, ...)`` unweighted and in whichever order the arrays happened to
be stored. Fitting ``E(k)`` when the error is on k is an errors-in-variables
problem, and the fix is the effective-variance method: weight by
``sigma_E**2 + (dE/dk * sigma_k)**2``, iterating because the slope appears
in its own weight. On a steep band the difference is not cosmetic.

**A velocity and a mass are local.** ``v_F`` is the slope *at E_F* and ``m*``
is the curvature *at the band extremum*; a straight line or a parabola fitted
over the whole measured range returns an average over whatever else the band
was doing, and in a material with a kink those differ by tens of percent.
:func:`window_scan` fits over a family of windows and returns the value and
its error against window width, so the plateau -- if there is one -- can be
seen and quoted, and its absence can be seen too.

**The constant.** ``EmassFitting.m`` uses ``m*/m_e = 7.631 / (2a)``, i.e.
3.8155 eV*A^2 for hbar^2/2m_e. The value is 3.80998 eV*A^2, so every mass it
has ever reported is 0.15% low. Small, but free to fix.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

#: hbar^2 / 2 m_e, in eV * A^2. E = HBAR2_OVER_2M * k^2 / (m*/m_e).
HBAR2_OVER_2M = 3.80998212
#: hbar in eV*s, for turning dE/dk in eV*A into a velocity in m/s.
HBAR_EVS = 6.582119569e-16
#: dE/dk [eV*A] * VELOCITY_FACTOR = v [m/s].
VELOCITY_FACTOR = 1e-10 / HBAR_EVS


@dataclass
class Dispersion:
    """A polynomial through a band, with what it implies."""
    coefficients: np.ndarray             # highest power first, E(k)
    errors: np.ndarray
    order: int
    k: np.ndarray
    energy: np.ndarray
    sigma: np.ndarray                    # the effective energy error used
    model: np.ndarray
    chi2: float = float("nan")
    error_axis: str = "y"
    covariance: np.ndarray = None

    def value(self, k):
        return np.polyval(self.coefficients, np.asarray(k, dtype=float))

    def slope(self, k):
        return np.polyval(np.polyder(self.coefficients),
                          np.asarray(k, dtype=float))

    def _slope_error(self, k: float) -> float:
        """Error on dE/dk at one k, from the full covariance.

        The derivative is a linear combination of the coefficients, so its
        variance needs the off-diagonal terms: for a parabola the linear and
        quadratic coefficients are strongly anticorrelated, and adding their
        variances in quadrature overestimates the error by a large factor.
        """
        if self.covariance is None or self.order < 1:
            return float("nan")
        powers = np.array([(self.order - i) * float(k) ** max(self.order - i - 1, 0)
                           for i in range(self.order)] + [0.0])
        variance = float(powers @ self.covariance @ powers)
        return float(np.sqrt(max(variance, 0.0)))

    def fermi_velocity(self, k_f: float = None):
        """``(dE/dk, error)`` at ``k_f``, in eV*A."""
        if k_f is None:
            k_f = float(np.mean(self.k))
        return float(self.slope(k_f)), self._slope_error(k_f)

    def velocity_ms(self, k_f: float = None):
        value, error = self.fermi_velocity(k_f)
        return value * VELOCITY_FACTOR, error * VELOCITY_FACTOR

    def effective_mass(self):
        """``(m*/m_e, error)`` from the quadratic term.

        Only defined for a parabola: a straight line has no curvature, and
        returning a number for it would be inventing one.
        """
        if self.order != 2:
            return float("nan"), float("nan")
        a = float(self.coefficients[0])
        if abs(a) < 1e-30:
            return float("nan"), float("nan")
        mass = HBAR2_OVER_2M / a
        error = abs(mass) * abs(float(self.errors[0]) / a)
        return mass, error

    def band_extremum(self):
        """Where a parabolic band turns over, as ``(k, E)``.

        The energy is the half that gets quoted -- the band bottom of an
        electron pocket, the top of a hole band -- so it comes back with
        the position rather than leaving the caller to evaluate the
        polynomial again and risk doing it on a different one.
        """
        if self.order != 2 or abs(self.coefficients[0]) < 1e-30:
            return float("nan"), float("nan")
        k = float(-self.coefficients[1] / (2.0 * self.coefficients[0]))
        return k, float(self.value(k))

    def crossing(self, level: float = 0.0):
        """Where the band crosses ``level`` (k_F for level = E_F)."""
        roots = np.roots(np.asarray(self.coefficients, dtype=float)
                         - np.append(np.zeros(self.order), float(level)))
        real = [float(r.real) for r in roots if abs(r.imag) < 1e-9]
        inside = [r for r in real if self.k.min() - 1e-9 <= r <= self.k.max() + 1e-9]
        return sorted(inside) if inside else sorted(real)

    def summary(self) -> str:
        lines = []
        value, error = self.fermi_velocity()
        speed, speed_error = self.velocity_ms()
        lines.append(f"dE/dk = {value:.4g} ± {error:.2g} eV·Å")
        lines.append(f"      = {speed:.4g} ± {speed_error:.2g} m/s")
        if self.order == 2:
            mass, mass_error = self.effective_mass()
            lines.append(f"m*/mₑ = {mass:.4g} ± {mass_error:.2g}")
            turn_k, turn_e = self.band_extremum()
            if np.isfinite(turn_k):
                edge = "bottom" if mass > 0 else "top"
                lines.append(f"band {edge} at k = {turn_k:.4g}, "
                             f"E = {turn_e:.4g}")
        lines.append(f"reduced χ² = {self.chi2:.3g}   "
                     f"({self.k.size} points)")
        return "\n".join(lines)


def fit_dispersion(k, energy, *, k_error=None, energy_error=None,
                   order: int = 1, iterations: int = 3) -> Dispersion:
    """Fit ``E(k)`` as a polynomial, weighted by whichever error is real.

    Give ``k_error`` for a band from MDC fits and ``energy_error`` for one
    from EDC fits. With ``k_error`` the weight is
    ``sigma_E_eff = |dE/dk| * sigma_k``, which needs the slope that is being
    fitted, so it is iterated -- three passes is far more than enough, the
    slope moves in the fourth decimal after the second.
    """
    k = np.asarray(k, dtype=float)
    energy = np.asarray(energy, dtype=float)
    good = np.isfinite(k) & np.isfinite(energy)
    if k_error is not None:
        k_error = np.asarray(k_error, dtype=float)
        good &= np.isfinite(k_error) & (k_error > 0)
    if energy_error is not None:
        energy_error = np.asarray(energy_error, dtype=float)
        good &= np.isfinite(energy_error) & (energy_error > 0)
    if good.sum() < order + 1:
        raise ValueError(f"need at least {order + 1} points for that order, "
                         f"have {int(good.sum())}")

    kf, ef = k[good], energy[good]
    design = np.vander(kf, order + 1)

    sigma = (energy_error[good] if energy_error is not None
             else np.ones(kf.shape))
    coefficients = None
    for _ in range(max(1, int(iterations))):
        scaled = design / sigma[:, None]
        coefficients, *_ = np.linalg.lstsq(scaled, ef / sigma, rcond=None)
        if k_error is None:
            break
        slope = np.polyval(np.polyder(coefficients), kf)
        extra = np.abs(slope) * k_error[good]
        base = (energy_error[good] if energy_error is not None
                else np.zeros(kf.shape))
        sigma = np.hypot(base, extra)
        sigma = np.where(sigma > 0, sigma, np.nanmedian(sigma[sigma > 0])
                         if np.any(sigma > 0) else 1.0)

    model = design @ coefficients
    dof = max(int(kf.size - coefficients.size), 1)
    chi2 = float(np.sum(((ef - model) / sigma) ** 2) / dof)
    try:
        normal = (design / sigma[:, None]).T @ (design / sigma[:, None])
        covariance = np.linalg.inv(normal) * chi2
        errors = np.sqrt(np.maximum(np.diag(covariance), 0.0))
    except np.linalg.LinAlgError:
        covariance = None
        errors = np.full(coefficients.size, np.nan)
    return Dispersion(coefficients, errors, int(order), kf, ef, sigma, model,
                      chi2, "x" if k_error is not None else "y", covariance)


@dataclass
class WindowScan:
    """How a velocity or a mass moves with the fitting window."""
    widths: np.ndarray
    values: np.ndarray
    errors: np.ndarray
    counts: np.ndarray
    quantity: str = "velocity"
    axis: str = "energy"
    centre: float = 0.0

    def plateau(self, tolerance: float = 0.05):
        """How far the window can be opened before the answer moves, as
        ``(first, last)`` indices into :attr:`widths`.

        Anchored at the **narrowest** window and extended outwards, not
        taken as the longest flat run anywhere in the scan. Both a velocity
        and a mass are local quantities -- defined at E_F and at the band
        extremum -- so the window that matters is the one closest to that
        point. Searching for the longest flat run instead finds the
        asymptote at large windows, which on a kinked band is a different
        velocity entirely and is not the one being asked for.

        A window leaves the plateau when it differs from the running median
        by more than ``tolerance`` **and** by more than ``n_sigma`` of its
        own error. Both halves are needed. The narrowest windows hold only a
        dozen points, so their velocities scatter by ten percent on noise
        alone and a bare relative test declares the plateau over
        immediately, even on a perfectly straight band; the widest hold
        hundreds, so their errors shrink to nothing and a bare statistical
        test finds a plateau nowhere.

        Returns None when even the two narrowest windows disagree, which is
        a real answer: it means the scan says nothing and a number should
        not be quoted from it.
        """
        return self._plateau(tolerance)

    def _plateau(self, tolerance: float, n_sigma: float = 2.0):
        finite = np.flatnonzero(np.isfinite(self.values))
        if finite.size < 2:
            return None
        start = int(finite[0])
        last = start
        for index in finite[1:]:
            index = int(index)
            block = self.values[start:index + 1]
            block = block[np.isfinite(block)]
            middle = float(np.median(block))
            if middle == 0:
                break
            deviation = abs(self.values[index] - middle)
            error = self.errors[index]
            allowed = tolerance * abs(middle)
            if np.isfinite(error):
                allowed = max(allowed, n_sigma * float(error))
            if deviation > allowed:
                break
            last = index
        return (start, last) if last > start else None

    def plateau_value(self, tolerance: float = 0.05):
        """``(value, error, widest_window)`` over the anchored plateau, or
        None. The error is the smallest on the plateau -- the widest window
        that is still honest -- rather than a combination that would pretend
        the windows were independent measurements."""
        found = self.plateau(tolerance)
        if found is None:
            return None
        first, last = found
        block = slice(first, last + 1)
        values = self.values[block]
        errors = self.errors[block]
        best = int(np.nanargmin(errors)) if np.isfinite(errors).any() else 0
        return float(values[best]), float(errors[best]), float(self.widths[last])


def window_scan(k, energy, *, k_error=None, energy_error=None,
                quantity: str = "velocity", centre: float = 0.0,
                widths=None, order: int = None, minimum: int = 5) -> WindowScan:
    """Refit over a family of windows and report how the answer moves.

    ``quantity`` is ``"velocity"`` -- a straight line over an energy window
    around ``centre`` (E_F) -- or ``"mass"``, a parabola over a momentum
    window around ``centre`` (the band extremum). The point is not to
    automate the choice of window but to show what it costs: a band with a
    kink gives a velocity that walks steadily with the window and no
    plateau, and that is a result about the band, not a failure of the fit.
    """
    k = np.asarray(k, dtype=float)
    energy = np.asarray(energy, dtype=float)
    if quantity == "velocity":
        axis, order = energy, 1 if order is None else order
        axis_name = "energy"
    elif quantity == "mass":
        axis, order = k, 2 if order is None else order
        axis_name = "momentum"
    else:
        raise ValueError(f"quantity {quantity!r} is not velocity or mass")

    if widths is None:
        reach = float(np.nanmax(np.abs(axis - centre)))
        widths = np.linspace(reach / 8.0, reach, 16)
    widths = np.asarray(widths, dtype=float)

    values, errors, counts = [], [], []
    for width in widths:
        inside = np.abs(axis - float(centre)) <= width
        if inside.sum() < max(order + 1, minimum):
            values.append(np.nan)
            errors.append(np.nan)
            counts.append(int(inside.sum()))
            continue
        try:
            fit = fit_dispersion(
                k[inside], energy[inside],
                k_error=None if k_error is None else np.asarray(k_error)[inside],
                energy_error=None if energy_error is None
                else np.asarray(energy_error)[inside],
                order=order)
            if quantity == "velocity":
                value, error = fit.fermi_velocity(
                    float(np.interp(centre, fit.energy[np.argsort(fit.energy)],
                                    fit.k[np.argsort(fit.energy)]))
                    if fit.energy.size > 1 else None)
            else:
                value, error = fit.effective_mass()
        except Exception:
            value, error = np.nan, np.nan
        values.append(value)
        errors.append(error)
        counts.append(int(inside.sum()))
    return WindowScan(widths, np.array(values), np.array(errors),
                      np.array(counts), quantity, axis_name, float(centre))


# ==========================================================================
# Self-energy
# ==========================================================================
@dataclass
class SelfEnergy:
    """``Re Sigma`` and ``Im Sigma`` from an MDC series and a bare band."""
    energy: np.ndarray
    real: np.ndarray
    real_error: np.ndarray
    imaginary: np.ndarray
    imaginary_error: np.ndarray
    bare: Dispersion = None
    kk_real: np.ndarray = None

    def consistency(self) -> float:
        """RMS difference between ``Re Sigma`` and the Kramers-Kronig
        transform of ``Im Sigma``, relative to the size of ``Re Sigma``.

        The two are not independent, so this is the check that decides
        whether the bare band was a defensible choice: one picked to make a
        kink look bigger fails it.
        """
        if self.kk_real is None:
            return float("nan")
        good = np.isfinite(self.kk_real) & np.isfinite(self.real)
        if not good.any():
            return float("nan")
        scale = float(np.sqrt(np.mean(self.real[good] ** 2))) or 1.0
        return float(np.sqrt(np.mean((self.kk_real[good] - self.real[good]) ** 2))
                     / scale)


def kramers_kronig(energy, imaginary):
    """``Re Sigma`` from ``Im Sigma`` by the principal-value integral."""
    energy = np.asarray(energy, dtype=float)
    imaginary = np.asarray(imaginary, dtype=float)
    order = np.argsort(energy)
    e, im = energy[order], imaginary[order]
    good = np.isfinite(im)
    out = np.full(e.shape, np.nan)
    for i in range(e.size):
        take = good.copy()
        take[i] = False
        if take.sum() < 2:
            continue
        out[i] = np.trapezoid(im[take] / (e[take] - e[i]), e[take]) / np.pi
    result = np.full(energy.shape, np.nan)
    result[order] = out
    return result


def self_energy(series, bare: Dispersion, *, band: int = 0,
                resolution: float = 0.0) -> SelfEnergy:
    """Extract the self-energy from a fitted MDC series and a bare band.

    ``Re Sigma(E) = E - eps_bare(k(E))`` and
    ``Im Sigma(E) = (1/2) dk(E) |v_bare|``, with ``dk`` the fitted FWHM.
    The resolution is removed in quadrature when it was not already folded
    into the peak shape; when it was, pass 0 and do not remove it twice.
    """
    if getattr(series, "direction", "mdc") != "mdc":
        raise ValueError("a self-energy comes from MDC widths: refit as MDCs")
    centres = np.asarray(series.centres, dtype=float)
    good = np.isfinite(centres[:, band]) & np.isfinite(series.positions)
    k = centres[good, band]
    e = np.asarray(series.positions, dtype=float)[good]
    k_error = np.asarray(series.centre_errors, dtype=float)[good, band]
    widths = np.asarray(series.widths, dtype=float)[good, band]
    width_error = np.asarray(series.width_errors, dtype=float)[good, band]

    if resolution > 0:
        widths = np.sqrt(np.maximum(widths ** 2 - float(resolution) ** 2, 0.0))

    velocity = np.abs(bare.slope(k))
    real = e - bare.value(k)
    real_error = velocity * k_error
    imaginary = 0.5 * widths * velocity
    imaginary_error = 0.5 * velocity * width_error

    kk = kramers_kronig(e, imaginary)
    if np.isfinite(kk).any() and np.isfinite(real).any():
        kk = kk - np.nanmean(kk) + np.nanmean(real)
    return SelfEnergy(e, real, real_error, imaginary, imaginary_error, bare, kk)


def bare_band_from_anchors(k, energy, anchors, *, order: int = 1) -> Dispersion:
    """A bare band through chosen anchor points.

    The bare band is an assumption, not a measurement, and the usual way to
    state it is "a straight line through the dispersion far from the kink
    and through E_F". ``anchors`` is a list of energies; the measured k at
    each is taken from the data and a polynomial put through those points.
    Keeping it explicit means the assumption is recorded rather than buried
    in whichever range happened to be fitted.
    """
    k = np.asarray(k, dtype=float)
    energy = np.asarray(energy, dtype=float)
    order_index = np.argsort(energy)
    picked_k, picked_e = [], []
    for anchor in anchors:
        picked_k.append(float(np.interp(float(anchor), energy[order_index],
                                        k[order_index])))
        picked_e.append(float(anchor))
    if len(picked_k) < order + 1:
        raise ValueError(f"need at least {order + 1} anchors for that order")
    return fit_dispersion(np.array(picked_k), np.array(picked_e), order=order)
