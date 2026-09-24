"""
Pictures to check a batch result by -- for a person, or for an agent that
can read images.

One PNG per dataset: a row of constant-energy maps and the two orthogonal
cuts through the middle (or through Gamma, once the map is in k). Contrast
is set per panel from percentiles, so a weak feature is visible and one hot
pixel does not blank the rest. These are quick looks, not figures: the
figure composer in the viewer is where a figure is made.
"""
from __future__ import annotations

import os

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                    # noqa: E402

from .center import energy_slice                                   # noqa: E402


def _levels(image, lo=1.0, hi=99.5):
    finite = np.asarray(image)[np.isfinite(image)]
    if finite.size == 0:
        return 0.0, 1.0
    a, b = np.percentile(finite, [lo, hi])
    return float(a), float(b if b > a else a + 1.0)


def _show(ax, image, x, y, xlabel, ylabel, title, cmap):
    image = np.asarray(image, dtype=float)
    vmin, vmax = _levels(image)
    ax.imshow(image.T, origin="lower", aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax,
              extent=[float(x[0]), float(x[-1]), float(y[0]), float(y[-1])],
              interpolation="nearest")
    ax.set_xlabel(xlabel, fontsize=8)
    ax.set_ylabel(ylabel, fontsize=8)
    ax.set_title(title, fontsize=8)
    ax.tick_params(labelsize=7)


def preview(ds, path: str, *, energies=(0.0, -0.2, -0.5, -1.0), width: float = 0.03,
            centre=None, cmap: str = "magma", title: str = None) -> str:
    """Write the quick look for ``ds`` to ``path`` and return it."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    if ds.kind == "cut":
        fig, ax = plt.subplots(figsize=(4, 4), dpi=110)
        _show(ax, ds.array(float), ds.axis("x"), ds.axis("y"),
              ds.label("x"), ds.label("y"), title or ds.name, cmap)
    elif ds.kind in ("map", "k_map", "kz_map", "kz_map_k"):
        cube = ds.array(np.float32)
        x, k, z = ds.axis("x"), ds.axis("k"), ds.axis("z")
        relative = "E_F" in ds.label("z")
        wanted = [e for e in energies if z.min() <= e <= z.max()] if relative else []
        if not wanted:              # kinetic axis, or nothing asked for is inside
            wanted = list(np.linspace(z.min(), z.max(), 6)[1:-1][::-1])
        n = len(wanted) + 2
        fig, axes = plt.subplots(1, n, figsize=(2.6 * n, 2.8), dpi=110)
        for ax, e in zip(axes, wanted):
            _show(ax, energy_slice(cube, z, e, width), x, k, ds.label("x"),
                  ds.label("k"), f"E = {e:+.2f} ± {width:g}", cmap)
            if centre is not None:
                ax.plot([centre[0]], [centre[1]], "c+", ms=10, mew=1.5)
        ix = int(np.argmin(np.abs(x - (centre[0] if centre else np.mean(x[[0, -1]])))))
        ik = int(np.argmin(np.abs(k - (centre[1] if centre else np.mean(k[[0, -1]])))))
        _show(axes[-2], cube[ix], k, z, ds.label("k"), ds.label("z"),
              f"cut at {ds.label('x').split(' (')[0]} = {x[ix]:.3g}", cmap)
        _show(axes[-1], cube[:, ik], x, z, ds.label("x"), ds.label("z"),
              f"cut at {ds.label('k').split(' (')[0]} = {k[ik]:.3g}", cmap)
        fig.suptitle(title or ds.name, fontsize=9)
    else:
        raise ValueError(f"no preview for a {ds.kind}")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    return path


def plot_reference(ref, path: str) -> str:
    """The edge per channel and the curve through it: the check that a
    reference fit is sound (smooth points, few rejected, small residual)."""
    ch = ref.channels
    angle, ef, err, ok = (np.asarray(ch[k]) for k in ("angle", "ef", "err", "ok"))
    fig, ax = plt.subplots(figsize=(5, 3.2), dpi=110)
    finite = np.isfinite(ef)
    ax.errorbar(angle[finite & ok], ef[finite & ok], yerr=err[finite & ok], fmt=".",
                ms=3, lw=0.5, label="channels used")
    ax.plot(angle[finite & ~ok], ef[finite & ~ok], "x", color="0.6", ms=4,
            label="rejected")
    grid = np.linspace(*ref.slit_range, 200)
    ax.plot(grid, np.polyval(ref.coeffs, grid), "r-", lw=1,
            label=f"order {ref.order}, rms {ref.residual_meV:.1f} meV")
    ax.set_xlabel("angle along slit")
    ax.set_ylabel("E_F (kinetic, eV)")
    ax.set_title(f"{os.path.basename(ref.path)}  hv={ref.photon_energy_eV}  "
                 f"{ref.lens_mode} {ref.pass_energy}", fontsize=8)
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    return path
