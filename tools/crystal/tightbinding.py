"""
tools/crystal/tightbinding.py
=============================
A one-band tight-binding model fitted to measured band points.

    E(k) = e0 + sum_n t_n S_n(k),      S_n(k) = sum_{R in shell n} cos(k . R)

where shell ``n`` is the set of lattice vectors at the ``n``-th distance
from the origin. One hopping per shell; in 2-D or 3-D, any lattice.

Ported from the lab's ``TB_zph_zmd`` tools (``TB_2D.m``, ``TB_3D.m`` and
their GUIs), with these changes:

* **A linear fit.** The model is linear in ``e0`` and every ``t_n``, so the
  least-squares answer is one linear solve, exact and immediate, with the
  standard errors from its covariance. The .m files built the formula as a
  string and handed it to MATLAB's nonlinear ``fit`` starting from every
  parameter = 1 with bounds of +-1000 -- which can stop away from the
  minimum, and gave no uncertainties.
* **The Fermi-surface-only case made explicit.** Points that all sit at
  one energy (a Fermi-surface contour) fix the *shape* of the band but not
  its energy scale: multiply every ``t_n`` and ``e0 - E`` by the same
  number and the contour does not move. The .m file subtracted E_F and fit
  anyway, so the answer was whatever the start point drifted to (and the
  trivial ``t = 0`` fits a single energy perfectly). Here a single-energy
  fit holds the first hopping at a value you choose (``t1``, by default -1
  eV, i.e. answers "in units of t1") and fits the rest -- the honest
  statement of what such data determine. Points at two or more energies
  (a few constant-energy contours, or a dispersion) fit everything. Even
  then a single contour fixes only about two shells: the cosines of more
  distant shells can be traded for one another along one curve, which the
  fit detects (the condition number of the scaled problem) and says.
* **Shells from a proper search.** All lattice vectors within the needed
  radius, grouped by length, for any 2-D or 3-D lattice; the .m file walked
  half a plane of a 2-D lattice by hand and had no 3-D equivalent of the
  same care.
* **Robust to stray points.** Band points traced automatically include a
  few that lie on something else (a flat region near a saddle point, a
  second band). The fit reweights iteratively (Tukey's biweight on the
  residual's robust spread), so such points are found and given no say,
  and how many were set aside is reported.
* **The in-plane orientation can be fitted.** The sample azimuth (how the
  lattice is turned against the map's kx axis) is often known only to a
  degree or two; :func:`fit_azimuth` scans it and keeps the best.

Qt-free, so it can be tested on its own.
"""
from __future__ import annotations

import numpy as np

#: Above this (scaled) condition number the fitted hoppings are not
#: determined by the points, whatever their residual.
ILL_CONDITIONED = 300.0

__all__ = ["ILL_CONDITIONED", "shells", "design_matrix", "TBFit", "fit", "evaluate",
           "fit_azimuth", "rotate", "lattice_2d", "read_points"]


def lattice_2d(a: float, b: float, gamma_deg: float) -> np.ndarray:
    """2-D real-space primitive vectors (rows) from a, b and the angle."""
    g = np.deg2rad(gamma_deg)
    return np.array([[a, 0.0], [b * np.cos(g), b * np.sin(g)]])


def rotate(vectors: np.ndarray, azimuth_deg: float) -> np.ndarray:
    """``vectors`` (rows) turned by ``azimuth_deg`` about the z axis."""
    vectors = np.asarray(vectors, dtype=float)
    c, s = np.cos(np.deg2rad(azimuth_deg)), np.sin(np.deg2rad(azimuth_deg))
    out = vectors.copy()
    out[:, 0] = c * vectors[:, 0] - s * vectors[:, 1]
    out[:, 1] = s * vectors[:, 0] + c * vectors[:, 1]
    return out


def shells(vectors: np.ndarray, count: int, tolerance: float = 1e-6):
    """The first ``count`` neighbour shells of the lattice spanned by
    ``vectors`` (rows, 2-D or 3-D).

    Returns ``[(distance, R array (m, d)), ...]``, nearest first.
    """
    vectors = np.asarray(vectors, dtype=float)
    dim = vectors.shape[0]
    if vectors.shape != (dim, dim) or dim not in (2, 3):
        raise ValueError("give 2 vectors in 2-D or 3 in 3-D")
    if abs(np.linalg.det(vectors)) < 1e-12:
        raise ValueError("the lattice vectors are not independent")
    reach = 2
    while True:
        ranges = [np.arange(-reach, reach + 1)] * dim
        grid = np.stack(np.meshgrid(*ranges, indexing="ij"), -1).reshape(-1, dim)
        grid = grid[np.any(grid != 0, axis=1)]
        points = grid @ vectors
        lengths = np.linalg.norm(points, axis=1)
        order = np.argsort(lengths)
        points, lengths = points[order], lengths[order]
        found = []
        start = 0
        while start < len(lengths) and len(found) < count:
            stop = start
            while stop < len(lengths) and lengths[stop] - lengths[start] <= \
                    tolerance * max(1.0, lengths[start]):
                stop += 1
            found.append((float(lengths[start]), points[start:stop]))
            start = stop
        # The search box must hold every vector up to the last shell found:
        # its inscribed radius is reach times the shortest height.
        heights = abs(np.linalg.det(vectors)) / np.array(
            [np.linalg.norm(np.cross(*np.delete(vectors, i, axis=0)))
             if dim == 3 else np.linalg.norm(vectors[1 - i])
             for i in range(dim)])
        if len(found) == count and found[-1][0] < reach * heights.min() - 1e-9:
            return found
        reach += 1
        if reach > 12:
            return found


def design_matrix(k_points: np.ndarray, shell_list) -> np.ndarray:
    """Columns: 1, S_1(k), S_2(k), ... for each k point (rows)."""
    k_points = np.atleast_2d(np.asarray(k_points, dtype=float))
    columns = [np.ones(len(k_points))]
    for _distance, vectors in shell_list:
        columns.append(np.cos(k_points @ np.asarray(vectors).T).sum(axis=1))
    return np.column_stack(columns)


class TBFit:
    """A fitted model: ``e0`` and ``t`` (one per shell), their standard
    errors (NaN for one held fixed), the rms residual, and whether only the
    shape was determined (one energy, ``t1`` held)."""

    def __init__(self, vectors, shell_list, e0, t, e0_err, t_err, rms,
                 n_points, shape_only, azimuth_deg=0.0):
        self.vectors = np.asarray(vectors, dtype=float)
        self.shells = shell_list
        self.e0 = float(e0)
        self.t = np.asarray(t, dtype=float)
        self.e0_err = float(e0_err)
        self.t_err = np.asarray(t_err, dtype=float)
        self.rms = float(rms)
        self.n_points = int(n_points)
        self.shape_only = bool(shape_only)
        self.azimuth_deg = float(azimuth_deg)
        self.condition = 1.0
        self.warning = ""
        self.outliers = 0
        self.weights = None
        self.score = 0.0

    def energy(self, k_points) -> np.ndarray:
        """E(k) for ``k_points`` (rows)."""
        return design_matrix(k_points, self.shells) @ np.concatenate(
            [[self.e0], self.t])

    def table(self):
        """``[(name, value, error, unit), ...]`` for display."""
        rows = [("e0", self.e0, self.e0_err, "eV")]
        for i, ((distance, vectors), value, error) in enumerate(
                zip(self.shells, self.t, self.t_err), start=1):
            rows.append((f"t{i}  (|R| = {distance:.4g} Å, {len(vectors)} neighbours)",
                         value, error, "eV"))
        return rows


def fit(k_points, energies, vectors, n_shells: int, *, weights=None,
        t1=None, energy_tolerance: float = 1e-6, robust: bool = True) -> TBFit:
    """Least-squares tight-binding fit.

    ``k_points``: (N, d) in Å⁻¹; ``energies``: (N,) in eV; ``vectors``: the
    real-space primitive vectors (rows, Å). ``t1`` holds the first hopping
    fixed: needed when every point has the same energy (see the module
    note); if None it is held at -1 eV in that case and left free
    otherwise.
    """
    k_points = np.atleast_2d(np.asarray(k_points, dtype=float))
    energies = np.asarray(energies, dtype=float).ravel()
    if len(k_points) != len(energies):
        raise ValueError("as many energies as k points are needed")
    shell_list = shells(vectors, n_shells)
    matrix = design_matrix(k_points, shell_list)
    single = float(np.ptp(energies)) <= energy_tolerance if energies.size else True
    hold = t1 if t1 is not None else (-1.0 if single else None)
    free = list(range(matrix.shape[1]))
    target = energies.copy()
    if hold is not None:
        target = target - hold * matrix[:, 1]
        free.remove(1)
    if len(energies) < len(free) + 1:
        raise ValueError(f"{len(free)} parameters need at least "
                         f"{len(free) + 1} points; there are {len(energies)}")
    base = np.ones_like(energies) if weights is None else np.asarray(weights, float)
    w = base.copy()
    for _round in range(30 if robust else 1):
        sw = np.sqrt(np.clip(w, 0, None))
        a = matrix[:, free] * sw[:, None]
        solution, *_ = np.linalg.lstsq(a, target * sw, rcond=None)
        residual = target - matrix[:, free] @ solution
        if not robust:
            break
        # Tukey's biweight on the residual in units of its robust spread:
        # a point more than ~4.7 sigma off the band (a ridge traced on
        # something else) gets no say, the rest keep nearly full weight.
        scale = 1.4826 * float(np.median(np.abs(residual - np.median(residual))))
        if scale <= 0:
            break
        u = residual / (4.685 * scale)
        new = base * np.where(np.abs(u) < 1, (1 - u ** 2) ** 2, 0.0)
        if np.allclose(new, w, atol=1e-4):
            break
        w = new
    outliers = int(np.sum(w <= 0)) if robust else 0
    dof = max(1, len(energies) - len(free))
    sigma2 = float(np.sum(w * residual ** 2) / dof)
    try:
        covariance = sigma2 * np.linalg.inv(a.T @ a)
        errors = np.sqrt(np.clip(np.diag(covariance), 0, None))
    except np.linalg.LinAlgError:
        errors = np.full(len(free), np.nan)
    params = np.zeros(matrix.shape[1])
    errs = np.full(matrix.shape[1], np.nan)
    params[free] = solution
    errs[free] = errors
    if hold is not None:
        params[1] = hold
    rms = float(np.sqrt(np.mean(residual ** 2))) if residual.size else 0.0
    if robust:
        kept = w > 0
        rms = float(np.sqrt(np.mean(residual[kept] ** 2))) if kept.any() else rms
    model = TBFit(vectors, shell_list, params[0], params[1:], errs[0], errs[1:],
                  rms, len(energies), single and t1 is None)
    model.outliers = outliers
    model.weights = w
    # What fit_azimuth ranks orientations by: the typical misfit over every
    # point, so an orientation cannot win by declaring points outliers.
    model.score = float(np.median(np.abs(residual))) if residual.size else 0.0
    # Columns scaled to unit length, so the condition number measures how
    # nearly the shells' cosines can stand in for one another on these
    # points -- not their sizes.
    scaled = a / np.maximum(np.linalg.norm(a, axis=0), 1e-300)
    model.condition = float(np.linalg.cond(scaled))
    if model.condition > ILL_CONDITIONED:
        model.warning = (
            f"The points do not pin down {n_shells} shells: on them, the "
            f"shells' cosines nearly stand in for one another (condition "
            f"number {model.condition:.3g}), so other hoppings describe the "
            f"same points about as well. Use fewer shells, or points at more "
            f"energies.")
    return model


def fit_azimuth(k_points, energies, vectors, n_shells: int, *, around=0.0,
                span=10.0, steps=81, **options) -> TBFit:
    """Fit with the lattice turned by every azimuth in ``around +- span``
    degrees (in-plane, about z) and keep the best; then refine to a tenth
    of the step."""
    def trial(angle):
        model = fit(k_points, energies, rotate(vectors, angle), n_shells,
                    **options)
        model.azimuth_deg = angle
        return model

    angles = np.linspace(around - span, around + span, steps)
    models = [trial(a) for a in angles]
    best = int(np.argmin([m.score for m in models]))
    step = angles[1] - angles[0] if steps > 1 else 1.0
    fine = np.linspace(angles[best] - step, angles[best] + step, 21)
    models = [trial(a) for a in fine]
    return models[int(np.argmin([m.score for m in models]))]


def evaluate(model: TBFit, kx, ky, kz: float = 0.0) -> np.ndarray:
    """E on the grid ``kx`` x ``ky`` (and the plane ``kz`` for a 3-D
    lattice), indexed ``[ix, iy]``."""
    kx, ky = np.asarray(kx, float), np.asarray(ky, float)
    gx, gy = np.meshgrid(kx, ky, indexing="ij")
    columns = [gx.ravel(), gy.ravel()]
    if model.vectors.shape[1] == 3:
        columns.append(np.full(gx.size, float(kz)))
    return model.energy(np.column_stack(columns)).reshape(gx.shape)


def read_points(path: str):
    """Band points from a CSV (or whitespace) table: ``(k (N, d), E (N,))``.

    The columns are found by their header -- ``kx``, ``ky``, ``kz`` and an
    energy column ("E", "Energy...", "E - E_F...") -- so the ridge points
    written by the ridge-detection window (energy, kx, ky, strength) are
    read as they are. Without a header the order is kx, ky[, kz], E.
    """
    import csv
    import re

    with open(path, encoding="utf-8-sig") as handle:
        text = handle.read()
    rows = [r for r in csv.reader(text.splitlines(), skipinitialspace=True) if r]
    if rows and len(rows[0]) == 1:
        rows = [line.split() for line in text.splitlines() if line.strip()]

    def numeric(row):
        try:
            [float(v) for v in row]
            return True
        except ValueError:
            return False

    header = None if numeric(rows[0]) else [h.strip().lower() for h in rows[0]]
    data = np.array([[float(v) for v in row] for row in rows[header is not None:]
                     if numeric(row)], dtype=float)
    if header is None:
        if data.shape[1] < 3:
            raise ValueError("need at least kx, ky and E columns")
        return data[:, :-1], data[:, -1]

    def find(pattern):
        for i, name in enumerate(header):
            if re.search(pattern, name):
                return i
        return None

    ix, iy, iz = find(r"^kx|k_x"), find(r"^ky|k_y"), find(r"^kz|k_z")
    ie = find(r"^e\b|^e[ _-]|energy|^e$|binding")
    if ix is None or iy is None or ie is None:
        raise ValueError("could not find kx, ky and energy columns in: "
                         + ", ".join(header))
    columns = [ix, iy] + ([iz] if iz is not None else [])
    return data[:, columns], data[:, ie]
