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
        period = _number(ds.info.get("batch.kz_period_invA")) if ds.kind == "kz_map_k" else None
        for ax, e in zip(axes, wanted):
            _show(ax, energy_slice(cube, z, e, width), x, k, ds.label("x"),
                  ds.label("k"), f"E = {e:+.2f} ± {width:g}", cmap)
            if centre is not None:
                ax.plot([centre[0]], [centre[1]], "c+", ms=10, mew=1.5)
            if period:
                _zone_lines(ax, x, period)
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


def _number(value):
    try:
        return float(np.asarray(value, dtype=float).reshape(-1)[0])
    except (TypeError, ValueError):
        return None


def _zone_lines(ax, kz, period):
    """Gamma planes (solid) and zone boundaries (dashed) along k_z, at
    multiples of the surface's period -- what V0 is judged by."""
    lo, hi = float(np.nanmin(kz)), float(np.nanmax(kz))
    for n in range(int(np.floor(lo / period)) - 1, int(np.ceil(hi / period)) + 2):
        for position, style in ((n * period, "-"), ((n + 0.5) * period, "--")):
            if lo <= position <= hi:
                ax.axvline(position, color="w", lw=0.6, ls=style, alpha=0.8)


def plot_kz_qc(qc: dict, path: str):
    """E_F against photon energy (a smooth curve means the fits are sound),
    and the V0 scan (how well the data overlaps itself one period along)."""
    calib, conv = qc.get("kz_calibrate") or {}, qc.get("kz_convert") or {}
    panels = [p for p in (calib.get("ef_per_hv"), (conv.get("v0_scan") or {}).get("curve")) if p]
    if not panels:
        return None
    fig, axes = plt.subplots(1, len(panels), figsize=(4.2 * len(panels), 3.2), dpi=110)
    axes = np.atleast_1d(axes)
    index = 0
    if calib.get("ef_per_hv"):
        rows = np.array(calib["ef_per_hv"], dtype=float)
        ok = rows[:, 2] > 0
        ax = axes[index]
        ax.plot(rows[ok, 0], rows[ok, 1], ".", label="fitted")
        ax.plot(rows[~ok, 0], rows[~ok, 1], "rx", label="interpolated")
        ax.set_xlabel("photon energy (eV)")
        ax.set_ylabel("E_F (energy axis units)")
        ax.set_title(f"spread {1000 * calib.get('fermi_level_spread_eV', 0):.0f} meV", fontsize=8)
        ax.legend(fontsize=7)
        index += 1
    scan = conv.get("v0_scan")
    if scan and scan.get("curve"):
        ax = axes[index]
        tried = np.array(scan["tried"], dtype=float)
        curve = np.array([np.nan if v is None else v for v in scan["curve"]], dtype=float)
        ax.plot(tried, curve, ".-", ms=3, label="overlap after one period")
        for vote in scan.get("votes", []):
            ax.axvline(vote, color="0.7", lw=0.6)
        if scan.get("best") is not None:
            ax.axvline(scan["best"], color="k", ls="--", lw=0.8,
                       label=f"V0 = {scan['best']:.1f} ± {scan['uncertainty_eV']:.1f} eV")
        viewer = (scan.get("period_scan") or {}).get("best")
        if viewer is not None:
            ax.axvline(viewer, color="r", ls=":", lw=0.8,
                       label=f"period scan: {viewer:.1f} eV")
        ax.set_xlabel("inner potential V0 (eV)")
        ax.set_ylabel("self-overlap")
        energies = ", ".join(f"{e:+.2f}" for e in scan.get("energies", []))
        ax.set_title(f"at E - E_F = {energies} eV (grey: each energy's vote)", fontsize=7)
        ax.legend(fontsize=7)
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


def plot_calc_match(ds, bands, segment, result, surface, path, *, work_function,
                    angle_offset=0.0):
    """The check on a match to a calculation, in four panels:

    1. the scan at normal emission (-d2I/dE2) with the calculated bands
       drawn through it at the matched V0 and shift, and the photon energies
       of the Gamma (solid) and boundary (dashed) planes;
    2. the agreement over V0 and the energy shift;
    3./4. the cut at the photon energy nearest a Gamma plane and nearest a
       boundary plane, in k_par, with the calculated in-plane bands from that
       point overlaid -- both directions, since the slit's azimuth relative
       to the crystal is not known here.
    """
    from . import calcbands
    hv, slit, energy = ds.axis("x"), ds.axis("k"), ds.axis("z")
    cube = ds.array(np.float32)
    v0, shift, z = (result["inner_potential_eV"], result["energy_shift_eV"],
                    result["renormalisation"])
    period = surface["period_invA"]
    planes = result["planes"]
    fig, axes = plt.subplots(1, 4, figsize=(17, 4.0), dpi=110)

    ax = axes[0]
    image = calcbands.normal_emission_image(cube, slit, energy, angle_offset)
    ax.imshow(image.T, origin="lower", aspect="auto", cmap="Greys",
              extent=[hv[0], hv[-1], energy[0], energy[-1]])
    grid_e = np.linspace(energy.min(), min(energy.max(), 0.05), 200)
    kz = calcbands.kz_of(hv, grid_e, v0=v0, work_function=work_function)
    reduced = np.mod(kz, period)
    reduced = np.minimum(reduced, period - reduced)
    values = bands.bands_between(segment["k_gamma"], segment["k_edge"],
                                 reduced / (period / 2)) / z + shift
    for band in values:                    # (hv, E): where E == band(hv, E)
        crossing = np.abs(band - grid_e[None, :])
        best = np.argmin(crossing, axis=1)
        e_line = grid_e[best]
        ok = crossing[np.arange(hv.size), best] < 0.02
        if ok.sum() > 2:
            ax.plot(np.where(ok, hv, np.nan), e_line, "-", color="tab:blue", lw=0.8, alpha=0.8)
    for p in planes:
        ax.axvline(p["hv_eV"], color="r", lw=0.8, ls="-" if p["plane"] == "G" else "--")
        ax.text(p["hv_eV"], energy.max(), "Γ" if p["plane"] == "G" else p["plane"],
                color="r", fontsize=7, ha="center", va="bottom")
    ax.set_ylim(energy.min(), energy.max())
    ax.set_xlabel("photon energy (eV)")
    ax.set_ylabel("E - E_F (eV)")
    ax.set_title(f"normal emission vs calculation (V0 {v0:.1f} eV, shift "
                 f"{1000 * shift:+.0f} meV)", fontsize=8)

    ax = axes[1]
    scores = np.array(result["scores"], dtype=float)
    ax.imshow(scores.T, origin="lower", aspect="auto", cmap="viridis",
              extent=[result["v0_values"][0], result["v0_values"][-1],
                      result["shifts"][0], result["shifts"][-1]])
    ax.plot([v0], [shift], "r+", ms=12, mew=2)
    ax.set_xlabel("inner potential V0 (eV)")
    ax.set_ylabel("calculation shift (eV)")
    ax.set_title(f"agreement (best {result['score']:.2f}; V0 "
                 f"{result['v0_range_eV'][0]:.1f}-{result['v0_range_eV'][1]:.1f} eV "
                 f"within 0.02)", fontsize=8)

    colours = ("tab:cyan", "tab:orange", "tab:green", "tab:pink")
    for ax, centre in zip(axes[2:], ("G", segment.get("edge_label", "A"))):
        chosen = [p for p in planes if p["plane"] == centre]
        if not chosen:
            ax.set_axis_off()
            ax.set_title(f"no {centre} plane inside the scan", fontsize=8)
            continue
        target = min(chosen, key=lambda p: abs(p["hv_eV"] - np.mean(hv)))["hv_eV"]
        i = int(np.argmin(np.abs(hv - target)))
        frame = cube[i]                                         # (slit, E)
        kinetic = hv[i] - work_function + energy
        k_par = 0.5123 * np.sqrt(np.clip(kinetic, 0, None))[None, :] * np.sin(
            np.radians(slit - angle_offset))[:, None]
        vmin, vmax = _levels(frame)
        ax.pcolormesh(k_par, np.broadcast_to(energy, frame.shape), frame, cmap="Greys",
                      vmin=vmin, vmax=vmax, shading="gouraud")
        for colour, (name, k0, k1) in zip(colours, [d for d in
                                                   calcbands.in_plane_directions(bands, centre)
                                                   if not {d[1], d[2]} == {segment["k_gamma"],
                                                                           segment["k_edge"]}]):
            fraction = np.linspace(0, 1, 120)
            length = abs(k1 - k0) * segment["scale"]
            along = fraction * length
            for n, band in enumerate(bands.bands_between(k0, k1, fraction) / z + shift):
                label = name.replace("G", "Γ") if n == 0 else None
                ax.plot(along, band, "-", color=colour, lw=0.7, label=label)
                ax.plot(-along, band, "-", color=colour, lw=0.7)
        ax.set_xlim(np.nanmin(k_par), np.nanmax(k_par))
        ax.set_ylim(energy.min(), energy.max())
        ax.set_xlabel("k_par along the slit (1/Å)")
        ax.set_ylabel("E - E_F (eV)")
        name = "Γ" if centre == "G" else centre
        ax.set_title(f"hv = {hv[i]:.1f} eV (nearest {name} plane, {target:.1f} eV)", fontsize=8)
        ax.legend(fontsize=6, loc="lower right")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    return path
