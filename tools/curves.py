"""
tools/curves.py
===============
One-dimensional data: EDCs, MDCs, spin EDCs -- as tables, with no Qt.

A curve dataset is a table: ``x`` (energy, or angle/momentum for an MDC)
against one or more **channels**, stored ``(n_points, n_channels)``. Its
channel names travel in the metadata, ``info["curve.channels"]``, joined
with ``|``, and what the numbers are in ``info["curve.value_label"]``.

**Uncertainties are channels too.** A channel called ``σ <name>`` is the
standard deviation of the channel called ``<name>``. That keeps a
polarisation and its error bar, or a spin-up spectrum and its counting
error, in one dataset that saves and reloads like any other -- and it is why
the operations here are not the generic ones in :mod:`tools.dataops`:
binning by averaging, applied to an uncertainty column, would make it
smaller by the wrong factor. Everything here knows which columns are
errors and propagates them.
"""
from __future__ import annotations

import numpy as np

SIGMA_PREFIX = "σ "
CHANNEL_SEPARATOR = "|"

NORMALISATIONS = ("max", "area", "region")

#: numpy 2 renamed trapz; either spelling, whichever this numpy has.
_trapezoid = getattr(np, "trapezoid", None) or getattr(np, "trapz")
BACKGROUNDS = ("constant", "linear", "shirley")


# --------------------------------------------------------------------------
# Channel bookkeeping
# --------------------------------------------------------------------------
def channel_names(info: dict, n_channels: int, kind: str = "") -> list:
    """The channel names recorded with a dataset, or numbered defaults.

    A plain EDC or MDC with no names is one channel called after its kind;
    anything that does not add up (a name list of the wrong length) falls
    back to numbering rather than mislabelling.
    """
    text = str((info or {}).get("curve.channels", "") or "")
    names = [n.strip() for n in text.split(CHANNEL_SEPARATOR)] if text else []
    if len(names) == n_channels and all(names):
        return names
    if n_channels == 1:
        kind = str(kind or (info or {}).get("_kind", "")).upper()
        return [kind if kind in ("EDC", "MDC") else "intensity"]
    return [f"channel {i}" for i in range(n_channels)]


def value_label(info: dict) -> str:
    return str((info or {}).get("curve.value_label", "") or "Intensity")


def is_sigma(name: str) -> bool:
    return str(name).startswith(SIGMA_PREFIX)


def sigma_name(name: str) -> str:
    return SIGMA_PREFIX + str(name)


def data_channels(names) -> list:
    """Indices of the channels that are data rather than uncertainties."""
    return [i for i, name in enumerate(names) if not is_sigma(name)]


def sigma_of(names, index: int):
    """Index of the uncertainty channel belonging to channel ``index``, or
    None."""
    wanted = sigma_name(names[index])
    for i, name in enumerate(names):
        if name == wanted:
            return i
    return None


def looks_like_counts(values) -> bool:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    return bool(finite.size and finite.min() >= 0
                and np.all(np.abs(finite - np.round(finite)) < 1e-6))


def curve_info(names, value_text: str) -> dict:
    """The metadata entries that describe a table's channels."""
    return {"curve.channels": CHANNEL_SEPARATOR.join(str(n) for n in names),
            "curve.value_label": str(value_text)}


def table(x, columns) -> "tuple[np.ndarray, np.ndarray, list]":
    """``(x, values, names)`` from ``[(name, array), ...]``, with the
    table sorted so ``x`` ascends."""
    x = np.asarray(x, dtype=float).reshape(-1)
    names = [str(name) for name, _ in columns]
    values = np.column_stack([np.asarray(col, dtype=float).reshape(-1)
                              for _, col in columns])
    if values.shape[0] != x.size:
        raise ValueError(f"{values.shape[0]} values against {x.size} x points")
    order = np.argsort(x, kind="stable")
    return x[order], values[order], names


# --------------------------------------------------------------------------
# Operations. Each takes and returns (x, values, names).
# --------------------------------------------------------------------------
def crop(x, values, names, lo: float, hi: float):
    x = np.asarray(x, dtype=float)
    keep = (x >= min(lo, hi)) & (x <= max(lo, hi))
    if keep.sum() < 2:
        raise ValueError("fewer than two points inside that range")
    return x[keep], np.asarray(values, dtype=float)[keep], list(names)


def rebin(x, values, names, factor: int, how: str = "sum"):
    """Merge every ``factor`` neighbouring points.

    ``how="sum"`` keeps counts as counts, so their Poisson statistics still
    hold afterwards; ``"mean"`` keeps the scale. Uncertainty channels are
    combined in quadrature either way -- ``sqrt(Σσ²)`` for a sum,
    ``sqrt(Σσ²)/n`` for a mean -- which is the whole reason this is not the
    generic compress. Points left over at the end are dropped rather than
    merged into a shorter, noisier last bin.
    """
    factor = int(factor)
    if factor < 1:
        raise ValueError("the bin factor must be at least 1")
    if how not in ("sum", "mean"):
        raise ValueError("how must be 'sum' or 'mean'")
    x = np.asarray(x, dtype=float)
    values = np.asarray(values, dtype=float)
    if factor == 1:
        return x.copy(), values.copy(), list(names)
    n = (x.size // factor) * factor
    if n < factor:
        raise ValueError(f"only {x.size} points; cannot bin by {factor}")
    new_x = x[:n].reshape(-1, factor).mean(axis=1)
    blocks = values[:n].reshape(-1, factor, values.shape[1])
    out = np.empty((blocks.shape[0], values.shape[1]))
    for column, name in enumerate(names):
        block = blocks[:, :, column]
        if is_sigma(name):
            combined = np.sqrt(np.nansum(block ** 2, axis=1))
            out[:, column] = combined / factor if how == "mean" else combined
        else:
            total = np.nansum(block, axis=1)
            out[:, column] = total / factor if how == "mean" else total
    return new_x, out, list(names)


def _region_mask(x, region):
    lo, hi = sorted(float(v) for v in region)
    mask = (np.asarray(x) >= lo) & (np.asarray(x) <= hi)
    if not mask.any():
        raise ValueError("the region contains no points")
    return mask


def normalise(x, values, names, how: str = "max", region=None,
              together: bool = True):
    """Divide by the maximum, the area, or the mean over a region.

    ``together=True`` uses **one** factor for every data channel -- the one
    worked out from their sum -- which is what keeps a spin-up and a
    spin-down spectrum comparable. ``False`` normalises each channel on its
    own. An uncertainty channel is divided by its data channel's factor.
    Returns ``(x, values, names, factors)``.
    """
    if how not in NORMALISATIONS:
        raise ValueError(f"how must be one of {NORMALISATIONS}")
    x = np.asarray(x, dtype=float)
    values = np.asarray(values, dtype=float).copy()
    data = data_channels(names)

    def factor_of(curve):
        if how == "max":
            return float(np.nanmax(curve))
        if how == "area":
            return float(abs(_trapezoid(np.nan_to_num(curve), x)))
        return float(np.nanmean(curve[_region_mask(x, region)]))

    factors = {}
    if together:
        common = factor_of(np.nansum(values[:, data], axis=1))
        factors = {i: common for i in data}
    else:
        factors = {i: factor_of(values[:, i]) for i in data}
    for i, factor in factors.items():
        if not np.isfinite(factor) or factor == 0:
            raise ValueError(f"cannot normalise {names[i]!r}: its "
                             f"{how} is {factor}")
        values[:, i] /= factor
        j = sigma_of(names, i)
        if j is not None:
            values[:, j] /= abs(factor)
    return x, values, list(names), factors


def subtract_background(x, values, names, how: str = "constant",
                        region=None):
    """Subtract a background from every data channel.

    ``constant`` is the mean over ``region``; ``linear`` a straight line
    fitted over ``region`` (typically a stretch well above E_F and one
    below the features); ``shirley`` the iterative Shirley background over
    the whole curve. Uncertainties are left as they were: a background
    estimated from many points adds little to any one point's error, and
    pretending otherwise with a guessed number would be worse.
    """
    from tools.peaks import shirley_background

    if how not in BACKGROUNDS:
        raise ValueError(f"how must be one of {BACKGROUNDS}")
    x = np.asarray(x, dtype=float)
    values = np.asarray(values, dtype=float).copy()
    backgrounds = {}
    for i in data_channels(names):
        curve = values[:, i]
        if how == "shirley":
            background = shirley_background(x, curve)
        else:
            if region is None:
                raise ValueError(f"a {how} background needs a region")
            mask = _region_mask(x, region) & np.isfinite(curve)
            if how == "constant" or mask.sum() < 2:
                background = np.full_like(curve, float(np.nanmean(curve[mask])))
            else:
                slope, intercept = np.polyfit(x[mask], curve[mask], 1)
                background = slope * x + intercept
        values[:, i] = curve - background
        backgrounds[names[i]] = background
    return x, values, list(names), backgrounds


def shift_x(x, delta: float):
    return np.asarray(x, dtype=float) - float(delta)


def with_poisson_sigma(values, names):
    """Add ``σ <name>`` = sqrt(counts) for every data channel that is raw
    counts and has no uncertainty yet. Returns ``(values, names, added)``."""
    values = np.asarray(values, dtype=float)
    names = list(names)
    extra, added = [], []
    for i in data_channels(names):
        if sigma_of(names, i) is None and looks_like_counts(values[:, i]):
            extra.append(np.sqrt(np.maximum(values[:, i], 0.0)))
            added.append(sigma_name(names[i]))
    if not extra:
        return values, names, []
    return (np.column_stack([values] + extra), names + added, added)
