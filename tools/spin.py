"""
tools/spin.py
=============
Spin polarisation from the channels of a spin-resolved EDC, with no Qt.

What a channel is
-----------------
Each spin channel is one spectrum, counted with the spin manipulator in one
setting and the detector target magnetised one way. The MBS software
records all of it in the channel's label::

    <0,0> +X  (GUI+Z)        setting <0,0>, target at +X, counts +Z as "up"
    <90,180> +X  (GUI-Z)     setting <90,180>, target at +X, counts -Z as "up"
    <0,0> -X  (GUI-Z)
    <90,180> -X  (GUI+Z)

:func:`parse_channel` reads that into a :class:`SpinChannel`. The sensitive
axis and its sign are what the analysis needs; the setting and the target
magnetisation are what decide which instrumental effects cancel.

The model
---------
A channel sensitive to ``+a`` counts

    I = N_g · η_m · (1 + S·P_a)

with ``N_g`` the transmission of manipulator setting ``g``, ``η_m`` the
target's reflectivity when magnetised ``m`` (not quite equal for the two
directions -- that difference is the *instrumental asymmetry*), ``S`` the
effective Sherman function and ``P_a`` the polarisation along ``a``.

**The cross ratio.** With ``G±`` the geometric mean of the channels
sensitive to ``±a``,

    A = (G+ − G−) / (G+ + G−),     P = A / S.

When every setting and every magnetisation appears once on each side --
the design in the file above: +Z is counted at (<0,0>, +X) and
(<90,180>, −X), −Z at (<0,0>, −X) and (<90,180>, +X) -- every ``N_g`` and
every ``η_m`` appears once in each product and cancels exactly, which a
single pair of channels cannot do. :func:`design_report` says which of
the two cancel for the channels actually present.

**Pairs, and the instrumental asymmetry.** Each setting on its own gives
``A_g = (I+ − I−)/(I+ + I−) ≈ S·P ± ε``, the sign of ``ε`` depending on
which magnetisation counted ``+a``. Two settings with opposite
arrangements therefore measure ``ε = (A_1 − A_2)/2`` directly -- reported,
because an ``ε`` that is large, or that changes across the spectrum, says
the target needs attention.

**Uncertainties.** Counting statistics throughout. With
``L = ln G+ − ln G−`` and ``k`` channels on each side,
``var L = (1/k²) Σ 1/I_i`` and ``A = tanh(L/2)``, so
``σ_A = ½(1 − A²) σ_L``. For one pair this is exactly the textbook
``sqrt((1 − A²)/N)``. ``σ_P = σ_A / S``: the Sherman function divides
the error as much as the signal, which is why spin-resolved data need so
many more counts than spin-integrated ones.

**Spin-resolved spectra.** ``I↑ = I(1 + P)/2`` and ``I↓ = I(1 − P)/2``,
with ``I`` the sum of the channels used. ``I`` and ``A`` are uncorrelated
to first order for Poisson counts, so
``σ↑² = ((1 + P)/2)² I + (I/2)² σ_P²``.

**A zero reference.** Optionally, the asymmetry averaged over an energy
window believed unpolarised (background well above E_F, a known
spin-degenerate feature) is taken as the instrumental offset and removed:
``L → L − L₀``. This is a *choice*, recorded in the result, and it will
remove real polarisation if the window has any.

The effective Sherman function is a property of the detector and its target
on the day. The default of 0.2 is what has been published for this
end station's FERRUM VLEED detector (Sci. Rep. 2023,
doi:10.1038/s41598-023-40145-1); use the
end station's own calibration where there is one.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field

import numpy as np

DEFAULT_SHERMAN = 0.2

_LABEL = re.compile(
    r"<\s*(?P<setting>[^>]*)>\s*(?P<mag>[+-]\s*[XYZ])\s*"
    r"\(\s*(?:GUI)?\s*(?P<sign>[+-])\s*(?P<axis>[XYZ])\s*\)", re.I)


@dataclass
class SpinChannel:
    index: int
    label: str
    setting: str = ""
    magnetisation: str = ""          # e.g. "+X"
    axis: str = ""                   # "X", "Y" or "Z"; "" if unknown
    sign: int = 0                    # +1 / -1; 0 if unknown

    @property
    def usable(self) -> bool:
        return bool(self.axis) and self.sign in (-1, 1)

    @property
    def magnet_sign(self) -> int:
        return -1 if self.magnetisation.startswith("-") else (
            1 if self.magnetisation.startswith("+") else 0)


def parse_channel(index: int, label: str) -> SpinChannel:
    """One ``SpinComp`` label as a :class:`SpinChannel`; unrecognised
    labels come back with no axis rather than a guessed one."""
    text = " ".join(str(label).split())
    # The loader prefixes "C<n> "; drop it if present.
    text = re.sub(r"^C\d+\s+", "", text)
    match = _LABEL.search(text)
    if not match:
        return SpinChannel(index, text)
    return SpinChannel(index, text,
                       setting=" ".join(match.group("setting").split()),
                       magnetisation=match.group("mag").replace(" ", "").upper(),
                       axis=match.group("axis").upper(),
                       sign=1 if match.group("sign") == "+" else -1)


def parse_channels(labels) -> list:
    return [parse_channel(i, label) for i, label in enumerate(labels)]


def channels_from_info(info: dict, n_channels: int) -> list:
    """The channels of a spin dataset, from its ``spin.component.<n>``
    entries (or the channel names if those are missing)."""
    labels = []
    for i in range(n_channels):
        label = (info or {}).get(f"spin.component.{i}")
        if label is None:
            names = str((info or {}).get("curve.channels", "")).split("|")
            label = names[i] if i < len(names) else f"channel {i}"
        labels.append(label)
    return parse_channels(labels)


def axes_available(channels) -> list:
    """Every axis with at least one channel on each side, in X, Y, Z order."""
    out = []
    for axis in "XYZ":
        signs = {c.sign for c in channels if c.usable and c.axis == axis}
        if signs == {-1, 1}:
            out.append(axis)
    return out


def pairs(channels, axis: str) -> list:
    """``[(setting, plus_channel, minus_channel), ...]``: channels counted
    at the same manipulator setting with opposite sensitivity along
    ``axis``."""
    out = []
    settings = sorted({c.setting for c in channels
                       if c.usable and c.axis == axis})
    for setting in settings:
        plus = [c for c in channels if c.usable and c.axis == axis
                and c.setting == setting and c.sign == 1]
        minus = [c for c in channels if c.usable and c.axis == axis
                 and c.setting == setting and c.sign == -1]
        if len(plus) == 1 and len(minus) == 1:
            out.append((setting, plus[0], minus[0]))
    return out


def design_report(channels, axis: str) -> dict:
    """What cancels in the cross ratio for these channels.

    ``transmission`` -- each manipulator setting appears equally often on
    the + and − side; ``reflectivity`` -- each target magnetisation does.
    """
    used = [c for c in channels if c.usable and c.axis == axis]
    plus = [c for c in used if c.sign == 1]
    minus = [c for c in used if c.sign == -1]

    def balanced(key):
        values = {key(c) for c in used}
        return all(sum(key(c) == v for c in plus) ==
                   sum(key(c) == v for c in minus) for v in values)

    return {"n_plus": len(plus), "n_minus": len(minus),
            "transmission": bool(used) and balanced(lambda c: c.setting),
            "reflectivity": bool(used) and balanced(lambda c: c.magnetisation)}


# --------------------------------------------------------------------------
# The analysis
# --------------------------------------------------------------------------
@dataclass
class SpinResult:
    x: np.ndarray
    axis: str
    method: str
    sherman: float
    asymmetry: np.ndarray
    sigma_asymmetry: np.ndarray
    polarisation: np.ndarray
    sigma_polarisation: np.ndarray
    intensity: np.ndarray
    up: np.ndarray
    down: np.ndarray
    sigma_up: np.ndarray
    sigma_down: np.ndarray
    #: per manipulator setting: (A, σA), before any zero reference
    pair_asymmetries: dict = field(default_factory=dict)
    #: weighted mean instrumental asymmetry and its error, if measurable
    instrumental: tuple = None
    #: the log-ratio removed by the zero reference, 0 if none
    zero_offset: float = 0.0
    channels_used: list = field(default_factory=list)
    design: dict = field(default_factory=dict)
    notes: list = field(default_factory=list)

    def summary(self) -> str:
        method = "cross ratio" if self.method == "cross" else self.method
        lines = [f"P_{self.axis} = A / S with S = {self.sherman:.3g}  "
                 f"({method})"]
        d = self.design
        if d:
            cancelled = [name for name, key in (("transmission", "transmission"),
                                                ("target reflectivity",
                                                 "reflectivity")) if d.get(key)]
            lines.append(
                f"{d['n_plus']} channel(s) +{self.axis}, {d['n_minus']} "
                f"−{self.axis}; cancels: "
                + (", ".join(cancelled) if cancelled else "nothing instrumental"))
        if self.instrumental is not None:
            eps, err = self.instrumental
            lines.append(f"instrumental asymmetry ε = {eps:+.4f} ± {err:.4f}"
                         + ("  (cancelled by the cross ratio)"
                            if self.method == "cross" and d.get("reflectivity")
                            else ""))
        if self.zero_offset:
            lines.append(f"zero reference removed an asymmetry of "
                         f"{np.tanh(self.zero_offset / 2):+.4f}")
        good = np.isfinite(self.polarisation) & (self.sigma_polarisation > 0)
        if good.any():
            w = 1.0 / self.sigma_polarisation[good] ** 2
            mean = float(np.sum(w * self.polarisation[good]) / np.sum(w))
            lines.append(f"weighted mean P = {mean:+.4f} ± "
                         f"{float(np.sqrt(1 / np.sum(w))):.4f}; median σ_P = "
                         f"{float(np.median(self.sigma_polarisation[good])):.3f}")
        return "\n".join(lines + list(self.notes))


def _bin(x, counts, factor):
    factor = max(int(factor), 1)
    if factor == 1:
        return np.asarray(x, float), np.asarray(counts, float)
    n = (len(x) // factor) * factor
    if n < factor:
        raise ValueError(f"only {len(x)} points; cannot bin by {factor}")
    xb = np.asarray(x, float)[:n].reshape(-1, factor).mean(axis=1)
    cb = np.asarray(counts, float)[:n].reshape(-1, factor,
                                               counts.shape[1]).sum(axis=1)
    return xb, cb


def _log_ratio(plus_counts, minus_counts):
    """``L`` and ``σ_L`` for the geometric means of the two sides."""
    k = plus_counts.shape[1]
    floor = 0.5              # an empty channel: half a count, not log(0)
    p = np.maximum(plus_counts, floor)
    m = np.maximum(minus_counts, floor)
    L = np.mean(np.log(p), axis=1) - np.mean(np.log(m), axis=1)
    var = (np.sum(1.0 / p, axis=1) + np.sum(1.0 / m, axis=1)) / k ** 2
    return L, np.sqrt(var)


def analyse(x, counts, channels, axis: str, *, sherman: float = DEFAULT_SHERMAN,
            method: str = "cross", zero_region=None,
            bin_factor: int = 1) -> SpinResult:
    """Polarisation along ``axis`` from a ``(n_points, n_channels)`` table
    of raw counts.

    ``method`` is ``"cross"`` (every usable channel for the axis) or
    ``"pair:<setting>"`` (the two channels of one manipulator setting).
    ``zero_region`` is an ``(x0, x1)`` window taken to be unpolarised.
    ``bin_factor`` sums neighbouring points first, which keeps them counts.
    """
    if not sherman or sherman <= 0 or sherman > 1:
        raise ValueError("the Sherman function must be in (0, 1]")
    counts = np.asarray(counts, dtype=float)
    if counts.ndim != 2 or counts.shape[1] != len(channels):
        raise ValueError(f"{counts.shape} counts for {len(channels)} channels")
    if np.nanmin(counts) < 0:
        raise ValueError("negative counts: these are not raw spin channels "
                         "(background-subtracted or normalised data cannot "
                         "be analysed with counting statistics)")
    x, counts = _bin(x, counts, bin_factor)
    notes = []

    usable = [c for c in channels if c.usable and c.axis == axis]
    if method == "cross":
        plus = [c for c in usable if c.sign == 1]
        minus = [c for c in usable if c.sign == -1]
        if not plus or not minus:
            raise ValueError(f"no channel pair along {axis}")
        if len(plus) != len(minus):
            raise ValueError(f"{len(plus)} channels count +{axis} and "
                             f"{len(minus)} count −{axis}; the cross ratio "
                             f"needs the same number on each side -- use "
                             f"one pair instead")
    elif method.startswith("pair:"):
        wanted = method[5:]
        found = [p for p in pairs(channels, axis) if p[0] == wanted]
        if not found:
            raise ValueError(f"no pair at setting <{wanted}> along {axis}")
        plus, minus = [found[0][1]], [found[0][2]]
    else:
        raise ValueError(f"method {method!r} is not 'cross' or 'pair:<setting>'")
    used = plus + minus
    design = design_report(used, axis)

    plus_counts = counts[:, [c.index for c in plus]]
    minus_counts = counts[:, [c.index for c in minus]]
    L, sigma_L = _log_ratio(plus_counts, minus_counts)

    zero = 0.0
    if zero_region is not None:
        lo, hi = sorted(float(v) for v in zero_region)
        inside = (x >= lo) & (x <= hi)
        if not inside.any():
            raise ValueError("the zero-reference window contains no points")
        L0, sigma_L0 = _log_ratio(plus_counts[inside].sum(axis=0, keepdims=True),
                                  minus_counts[inside].sum(axis=0, keepdims=True))
        zero = float(L0[0])
        L = L - zero
        sigma_L = np.sqrt(sigma_L ** 2 + float(sigma_L0[0]) ** 2)

    A = np.tanh(L / 2.0)
    sigma_A = 0.5 * (1.0 - A ** 2) * sigma_L
    P = A / sherman
    sigma_P = sigma_A / sherman
    if np.any(np.abs(P) - 2 * sigma_P > 1.0):
        notes.append("|P| exceeds 1 by more than 2σ somewhere: the Sherman "
                     "function is probably too small, or the channels are "
                     "not what their labels say")

    intensity = plus_counts.sum(axis=1) + minus_counts.sum(axis=1)
    up = intensity * (1.0 + P) / 2.0
    down = intensity * (1.0 - P) / 2.0
    sigma_up = np.sqrt(((1 + P) / 2) ** 2 * intensity + (intensity / 2) ** 2
                       * sigma_P ** 2)
    sigma_down = np.sqrt(((1 - P) / 2) ** 2 * intensity + (intensity / 2) ** 2
                         * sigma_P ** 2)

    # -- per-setting asymmetries and the instrumental one -------------------
    pair_list = pairs(channels, axis)
    pair_asym = {}
    signs = {}
    for setting, p, m in pair_list:
        Lp, sLp = _log_ratio(counts[:, [p.index]], counts[:, [m.index]])
        Ap = np.tanh(Lp / 2.0)
        pair_asym[setting] = (Ap, 0.5 * (1 - Ap ** 2) * sLp)
        signs[setting] = p.magnet_sign
    instrumental = None
    plus_arr = [s for s in pair_asym if signs[s] == 1]
    minus_arr = [s for s in pair_asym if signs[s] == -1]
    if plus_arr and minus_arr:
        a1, s1 = pair_asym[plus_arr[0]]
        a2, s2 = pair_asym[minus_arr[0]]
        eps = (a1 - a2) / 2.0
        s_eps = 0.5 * np.sqrt(s1 ** 2 + s2 ** 2)
        good = np.isfinite(eps) & (s_eps > 0)
        if good.any():
            w = 1.0 / s_eps[good] ** 2
            instrumental = (float(np.sum(w * eps[good]) / np.sum(w)),
                            float(np.sqrt(1.0 / np.sum(w))))

    return SpinResult(x=x, axis=axis, method=method, sherman=float(sherman),
                      asymmetry=A, sigma_asymmetry=sigma_A, polarisation=P,
                      sigma_polarisation=sigma_P, intensity=intensity, up=up,
                      down=down, sigma_up=sigma_up, sigma_down=sigma_down,
                      pair_asymmetries=pair_asym, instrumental=instrumental,
                      zero_offset=zero, channels_used=[c.index for c in used],
                      design=design, notes=notes)
