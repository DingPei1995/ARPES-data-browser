# Recipe: two-Lorentzian + linear-background MDC fit at fixed energy on a series of k-converted cuts,
# with the viewer's own fit engine (arpes_viewer.tools.spectra.peaks.fit_line).
# Writes agentic/out/mdc_fit.json (positions, errors, widths, curves) for plotting/reporting.
import json, re, os, numpy as np
from arpes_viewer.tools.spectra import peaks as PK
ROWS = r"Cut_LH_(\d+)K \[k G15\] \[E-Ef\]_s_nor$"       # normalised, k-converted cuts to fit
E0, HW = -0.1, 0.010                                    # MDC at E - E_F = E0, integrated +-HW (eV)
KLO, KHI = -1.0, -0.4                                   # fit range (1/A)
SEEDS = {"A": (-0.85, 0.10), "B": (-0.66, 0.06)}        # start centre, width (ASK the user / inspect first)

by = {r["name"]: k for k, r in A.loaded_items.items()}
rows = sorted((int(m.group(1)), n) for n in by for m in [re.match(ROWS, n)] if m)
bands = [PK.Band(n, [], "lorentzian") for n in SEEDS]
st = PK.FitSettings(background="linear", weighting="none")   # unweighted; errors from covariance x reduced chi2
res, guess = [], None
for T, name in rows:
    s = A.load_dataset(by[name]).scan
    k, E, V = np.asarray(s.x, float), np.asarray(s.y, float), np.asarray(s.value, float)
    mE = np.abs(E - E0) <= HW + 1e-9; mk = (k >= KLO) & (k <= KHI)
    y = np.nanmean(V[np.ix_(mk, mE)], axis=1); x = k[mk]; y = y / np.nanmax(y); h = np.nanmax(y) - np.nanmin(y)
    g = guess or {n: (c, w, 0.4 * h) for n, (c, w) in SEEDS.items()}
    f = PK.fit_line(x, y, bands, st, guesses=g, position=E0)
    order = np.argsort(f.centres); c = [float(f.centres[i]) for i in order]; e = [float(f.centre_errors[i]) for i in order]
    w = [float(f.widths[i]) for i in order]
    guess = {n: (c[i], w[i], 0.4 * h) for i, n in enumerate(SEEDS)}   # seed the next temperature
    dk = abs(c[1] - c[0]); dke = float(np.hypot(e[0], e[1]))
    res.append(dict(T=T, c1=c[0], c2=c[1], e1=e[0], e2=e[1], w1=w[0], w2=w[1], dk=dk, dk_err=dke, chi2=float(f.chi2),
                    x=x.tolist(), y=y.tolist(), model=np.asarray(f.model).tolist()))
    print(f"{T:5d} K  k1={c[0]:.4f}+-{e[0]:.4f}  k2={c[1]:.4f}+-{e[1]:.4f}  dk={dk:.4f}+-{dke:.4f}  ok={f.success}", flush=True)
json.dump(res, open(os.path.join(HERE, "out", "mdc_fit.json"), "w"))
