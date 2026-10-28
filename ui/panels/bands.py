"""
ui/panels/bands.py
==================
Calculated band structures, compared with the measurement, step by step.

Opened from the tool panel's Calculation tab (Calculated bands...), with
the measurements to compare selected -- maps converted to momentum
(k-maps), converted photon-energy scans (kz maps in k) and cuts converted
to k -- or with nothing selected, to look at a calculation alone. Four
steps:

1. **The calculation** -- a file from the theory side: an XCrySDen
   ``.bxsf`` band grid (Wien2k, Quantum ESPRESSO, VASP, Wannier90), a
   Wannier90 Hamiltonian (``_hr.dat``), or bands along a path (Wannier90
   ``_band.dat``, Wien2k ``.spaghetti_ene``, VASP ``EIGENVAL`` with its
   ``PROCAR``). The crystal's space group names the high-symmetry points.
2. **Along a path** -- the bands through the high-symmetry points
   (Γ-M-K-Γ-A-L-H-A), optionally coloured by orbital character, over the
   measurement along the same path, assembled from the selected maps (the
   Γ-plane and zone-top maps, a k_z scan for Γ-A).
3. **Compare** -- one measured cut or constant-energy map at a time, the
   calculation drawn over it where that measurement really was (its k_z
   plane, or the curved k_z surface its photon energy reaches), and the
   renormalisation Z and rigid shift that fit it, from the cut's own
   bands.
4. **Output** -- the calculation's spectral function to the list, and
   figure-composer panels: measurement with calculation, calculation
   alone (with the points' names on the axis), constant-energy map with
   calculated contours.

The arithmetic is :mod:`tools.crystal.bands`, :mod:`tools.crystal.compare`
and :mod:`tools.crystal.bandfiles`, none of which imports Qt.
"""
from __future__ import annotations

import os
import warnings

import numpy as np
import pyqtgraph as pg
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtWidgets import (QAbstractItemView, QCheckBox, QComboBox, QDialog,
                             QDoubleSpinBox, QFileDialog, QFormLayout,
                             QGroupBox, QHBoxLayout, QHeaderView,
                             QInputDialog, QLabel, QLineEdit, QPushButton,
                             QRadioButton, QSlider, QSpinBox, QTableWidget,
                             QTableWidgetItem, QVBoxLayout, QWidget)

from tools.crystal import bandfiles
from tools.crystal import bands as B
from tools.crystal import compare as C
from tools.crystal import last_crystal
from tools.crystal.lattice import LatticeParams
from ui.common.screen import fits_screen
from ui.common.stepper import Step, Stepper
from ui.viewers.widgets import (MemoryData, apply_colormap, fit_frame_view,
                                plain_image_view, show_frame)

__all__ = ["BandsDialog", "DEFAULT_PATHS", "params_of_lattice"]

#: a path to start from, by crystal family
DEFAULT_PATHS = {"hexagonal": "Γ-M-K-Γ-A-L-H-A",
                 "tetragonal": "Γ-X-M-Γ-Z-R-A-Z",
                 "cubic": "Γ-X-M-Γ-R-X",
                 "other": "Γ-X-S-Y-Γ-Z"}

#: colours of the weight (orbital character) classes
WEIGHT_COLOURS = ("#d62728", "#1f77b4", "#2ca02c", "#ff7f0e", "#9467bd",
                  "#8c564b")

ACCEPTED = ("k_map", "kz_map_k", "cut")


def _spin(value, lo, hi, decimals=3, suffix="", step=None):
    box = QDoubleSpinBox()
    box.setDecimals(decimals)
    box.setRange(lo, hi)
    box.setValue(value)
    box.setSuffix(suffix)
    box.setKeyboardTracking(False)
    if step:
        box.setSingleStep(step)
    return box


def params_of_lattice(lattice, space_group: int = None) -> LatticeParams:
    """The cell parameters of real-space ``lattice`` (rows, Å), with a
    guessed space group when none is given: P6/mmm for a hexagonal cell,
    P4/mmm for a tetragonal one, Pm-3m cubic, else Pmmm."""
    v = np.asarray(lattice, dtype=float)
    a, b, c = (float(np.linalg.norm(r)) for r in v)

    def angle(p, q):
        return float(np.degrees(np.arccos(np.clip(p @ q / (np.linalg.norm(p)
                                                        * np.linalg.norm(q)),
                                                  -1, 1))))

    alpha, beta, gamma = angle(v[1], v[2]), angle(v[0], v[2]), angle(v[0], v[1])
    if space_group is None:
        right = all(abs(x - 90) < 0.5 for x in (alpha, beta))
        if right and abs(a - b) < 1e-3 * a and abs(gamma - 120) < 0.5:
            space_group = 191
        elif right and abs(gamma - 60) < 0.5 and abs(a - b) < 1e-3 * a:
            # a 60-degree hexagonal setting: the same zone
            space_group, gamma = 191, 120.0
        elif right and abs(gamma - 90) < 0.5:
            if abs(a - b) < 1e-3 * a:
                space_group = 221 if abs(a - c) < 1e-3 * a else 123
            else:
                space_group = 47
        else:
            space_group = 1
    return LatticeParams(a, b, c, alpha, beta, gamma, space_group=int(space_group))


def _family(space_group: int) -> str:
    n = int(space_group or 1)
    if 143 <= n <= 194:
        return "hexagonal"
    if 75 <= n <= 142:
        return "tetragonal"
    if n >= 195:
        return "cubic"
    return "other"


def _relative_energy(axis, info) -> np.ndarray:
    """An energy axis as E - E_F: a kinetic-energy axis (all of it well
    above zero) has hv - W taken off."""
    from tools.spectra.corelevel import photon_energy_of, work_function_of

    axis = np.asarray(axis, dtype=float)
    if axis.size and float(np.min(axis)) > 3.0:
        hv = photon_energy_of(info or {})
        if hv:
            return axis - (hv - work_function_of(info or {}))
    return axis


@fits_screen
class BandsDialog(QDialog):
    """Calculated bands against the measurement. ``datasets`` is a list
    of ``(name, data)``: k-maps, kz maps in k and cuts in k."""

    datasetsCreated = pyqtSignal(list)

    def __init__(self, datasets=None, parent=None, colormap="blue",
                 existing_names=None):
        super().__init__(parent)
        self.setWindowTitle("Calculated bands")
        self.resize(1180, 760)
        self.datasets = [(n, d) for n, d in (datasets or [])
                         if getattr(d, "kind", None) in ACCEPTED]
        self.colormap = colormap
        self._existing = existing_names or (lambda: [])
        self.model = None
        self.model_path = ""
        self.params = None
        self.alignment = B.Alignment()
        self._items = []
        self._kz_typed = set()

        layout = QVBoxLayout(self)
        self.stepper = Stepper([self._step_calculation(), self._step_path(),
                                self._step_compare(), self._step_output()],
                               self, finish_text="Close")
        self.stepper.finished.connect(self.accept)
        layout.addWidget(self.stepper)

    # ==================================================================
    # Step 1: the calculation
    # ==================================================================
    def _step_calculation(self):
        page = QWidget()
        column = QVBoxLayout(page)
        row = QHBoxLayout()
        self.open_button = QPushButton("Open calculation...")
        self.open_button.clicked.connect(lambda: self.open_file())
        row.addWidget(self.open_button)
        self.file_label = QLabel("No calculation yet.")
        self.file_label.setWordWrap(True)
        row.addWidget(self.file_label, 1)
        column.addLayout(row)

        options = QGroupBox("Reading it")
        form = QFormLayout(options)
        self.fermi_box = QCheckBox("Set E_F myself")
        self.fermi_value = _spin(0.0, -1000.0, 1000.0, 4, " eV")
        fermi = QHBoxLayout()
        fermi.addWidget(self.fermi_box)
        fermi.addWidget(self.fermi_value)
        fermi.addStretch(1)
        form.addRow("Fermi level", fermi)
        self.k_unit = QComboBox()
        self.k_unit.addItems(["Å⁻¹, 2π included", "Å⁻¹, without 2π",
                              "bohr⁻¹, 2π included"])
        self.k_unit.setToolTip("Units of a .bxsf file's reciprocal vectors. "
                               "Programs differ; if the zone comes out 2π "
                               "too small or large, change this.")
        form.addRow("bxsf k units", self.k_unit)
        column.addWidget(options)

        crystal = QGroupBox("The crystal (names the high-symmetry points)")
        cform = QFormLayout(crystal)
        self.space_group = QSpinBox()
        self.space_group.setRange(1, 230)
        self.space_group.setValue(191)
        self.space_group.setToolTip("Sets the zone's shape and its labelled "
                                    "points; the cell itself comes from the "
                                    "calculation.")
        self.space_group.valueChanged.connect(lambda *_: self._crystal_changed())
        sg_row = QHBoxLayout()
        sg_row.addWidget(self.space_group)
        take = QPushButton("From the Brillouin-zone tool")
        take.setToolTip(last_crystal.describe())
        take.clicked.connect(lambda: self.take_space_group())
        sg_row.addWidget(take)
        sg_row.addStretch(1)
        cform.addRow("Space group", sg_row)
        self.crystal_label = QLabel("")
        self.crystal_label.setWordWrap(True)
        cform.addRow(self.crystal_label)
        column.addWidget(crystal)
        column.addStretch(1)
        text = ("The calculation, as the theory side sends it. Two kinds: "
                "<b>bands everywhere</b> -- a band grid over the zone "
                "(XCrySDen <code>.bxsf</code>, written by Wien2k, Quantum "
                "ESPRESSO, VASP and Wannier90) or a Wannier Hamiltonian "
                "(<code>_hr.dat</code>, with its <code>.win</code> beside it "
                "for the cell) -- which can be cut anywhere: along any "
                "measured cut, on the k_z a photon energy reaches, at any "
                "energy as contours; and <b>bands along a path</b> "
                "(<code>_band.dat</code>, <code>.spaghetti_ene</code>, VASP "
                "<code>EIGENVAL</code>), which can be drawn along that path "
                "only. Orbital weights come with a Hamiltonian (its Wannier "
                "functions) or a VASP <code>PROCAR</code>.<br>The space "
                "group names the zone's points; it is guessed from the "
                "cell's shape -- set the crystal's own for its labels.")
        return Step("The calculation", text, page,
                    ready=lambda: "" if self.model is not None
                    else "Open a calculation first.")

    def open_file(self, path=None, **options):
        if path is None:
            filters = ";;".join([label for label, _p in bandfiles.FORMATS]
                                + ["All files (*)"])
            path, _ = QFileDialog.getOpenFileName(self, "Calculated bands", "",
                                                  filters)
            if not path:
                return None
        kwargs = dict(options)
        if self.fermi_box.isChecked():
            kwargs.setdefault("fermi", self.fermi_value.value())
        if path.lower().endswith(".bxsf"):
            kwargs.setdefault("k_scale", (1.0, 2 * np.pi, 1.0 / bandfiles.BOHR)[
                self.k_unit.currentIndex()])
        try:
            model = bandfiles.read(path, **kwargs)
        except (OSError, ValueError, IndexError) as exc:
            self.file_label.setText(f"<b>Could not read it:</b> {exc}")
            return None
        self.model_path = path
        self.set_model(model)
        return model

    def set_model(self, model, space_group: int = None):
        """Use ``model`` (a :class:`tools.crystal.bands.BandGrid`,
        ``WannierModel`` or ``BandPath``)."""
        self.model = model
        if isinstance(model, B.BandPath):
            self.params = None
            kind = (f"bands along a path ({len(model.labels)} labelled "
                    f"points)")
            count = model.energies.shape[0]
        else:
            self.params = params_of_lattice(model.lattice, space_group)
            self.space_group.blockSignals(True)
            self.space_group.setValue(int(self.params.space_group))
            self.space_group.blockSignals(False)
            kind = ("a Hamiltonian, diagonalised at every k"
                    if isinstance(model, B.WannierModel)
                    else "a band grid, interpolated")
            count = model.band_count()
        weights = list(getattr(model, "weight_names", []) or [])
        name = os.path.basename(self.model_path) or getattr(model, "name", "") \
            or "calculation"
        self.file_label.setText(
            f"<b>{name}</b>: {kind}; {count} bands; E_F = "
            f"{getattr(model, 'fermi', 0.0):.4g} eV"
            + (f"; weights: {', '.join(weights)}" if weights else "") + ".")
        self.weight_combo.blockSignals(True)
        self.weight_combo.clear()
        self.weight_combo.addItem("(one colour)")
        self.weight_combo.addItems(weights)
        self.weight_combo.blockSignals(False)
        self._crystal_changed()
        self._fill_compare_list()
        self.stepper.refresh()

    def take_space_group(self):
        params = last_crystal.lattice_params()
        if params is not None and params.space_group:
            self.space_group.setValue(int(params.space_group))

    def _crystal_changed(self):
        if self.model is None or isinstance(self.model, B.BandPath):
            self.crystal_label.setText(
                "A band structure along a path carries its own labels."
                if self.model is not None else "")
            return
        self.params = params_of_lattice(self.model.lattice,
                                        self.space_group.value())
        family = _family(self.params.space_group)
        try:
            from tools.crystal.highsym import high_symmetry_points
            names = ", ".join(p.label for p in high_symmetry_points(self.params))
        except Exception as exc:                              # noqa: BLE001
            names = f"(none: {exc})"
        p = self.params
        self.crystal_label.setText(
            f"a = {p.a:.4g} Å, b = {p.b:.4g} Å, c = {p.c:.4g} Å, "
            f"α = {p.alpha:.4g}°, β = {p.beta:.4g}°, γ = {p.gamma:.4g}° "
            f"({family}); labelled points: {names}.")
        if not self.path_edit.text().strip() or self.path_edit.text() in DEFAULT_PATHS.values():
            self.path_edit.setText(DEFAULT_PATHS[family])

    # ==================================================================
    # Step 2: along a path
    # ==================================================================
    def _step_path(self):
        page = QWidget()
        column = QHBoxLayout(page)
        left = QVBoxLayout()
        form = QFormLayout()
        self.path_edit = QLineEdit("")
        self.path_edit.setPlaceholderText("Γ-M-K-Γ-A-L-H-A   (G for Γ, | for a jump)")
        self.path_edit.editingFinished.connect(lambda: self.draw_path())
        form.addRow("Path", self.path_edit)
        window = QHBoxLayout()
        self.e_lo = _spin(-1.5, -50.0, 50.0, 3, " eV", 0.1)
        self.e_hi = _spin(0.3, -50.0, 50.0, 3, " eV", 0.1)
        for box in (self.e_lo, self.e_hi):
            box.valueChanged.connect(lambda *_: self.draw_path())
            window.addWidget(box)
        form.addRow("Energy window", window)
        self.weight_combo = QComboBox()
        self.weight_combo.addItem("(one colour)")
        self.weight_combo.setToolTip("Draw the bands as fat bands: the "
                                     "thickness is the chosen orbital's "
                                     "weight in each state (optional).")
        self.weight_combo.currentIndexChanged.connect(lambda *_: self.draw_path())
        form.addRow("Weight", self.weight_combo)
        self.line_colour = QComboBox()
        self.line_colour.addItems(["black", "white", "red", "blue"])
        self.line_colour.currentIndexChanged.connect(lambda *_: self.draw_path())
        form.addRow("Line colour", self.line_colour)
        left.addLayout(form)

        group = QGroupBox("The measurement along the path")
        g = QVBoxLayout(group)
        self.path_measured = QCheckBox("Show it, from the maps below")
        self.path_measured.setChecked(True)
        self.path_measured.toggled.connect(lambda *_: self.draw_path())
        g.addWidget(self.path_measured)
        self.source_table = QTableWidget(0, 4)
        self.source_table.setHorizontalHeaderLabels(
            ["Map", "Use", "k_z (Å⁻¹) / direction (°)", "Azimuth (°)"])
        self.source_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.source_table.setSelectionMode(QAbstractItemView.NoSelection)
        self.source_table.itemChanged.connect(self._source_edited)
        g.addWidget(self.source_table)
        inner = QFormLayout()
        self.v0_box = _spin(12.0, 0.0, 100.0, 2, " eV")
        self.v0_box.setToolTip("Inner potential: gives each map's k_z from "
                               "its photon energy (free-electron final state).")
        self.v0_box.valueChanged.connect(lambda *_: self._fill_sources())
        inner.addRow("Inner potential V0", self.v0_box)
        self.fold_box = QSpinBox()
        self.fold_box.setRange(0, 6)
        self.fold_box.setSpecialValueText("from the crystal")
        self.fold_box.setToolTip("Rotation symmetry used to find a path point "
                                 "the map does not cover at an equivalent "
                                 "point it does (with the mirrors).")
        self.fold_box.valueChanged.connect(lambda *_: self.draw_path())
        inner.addRow("Symmetry fold", self.fold_box)
        g.addLayout(inner)
        left.addWidget(group)
        self.path_label = QLabel("")
        self.path_label.setWordWrap(True)
        left.addWidget(self.path_label)
        left.addStretch(1)
        column.addLayout(left, 2)
        self.path_view = plain_image_view("", "E - E_F (eV)")
        column.addWidget(self.path_view, 3)
        text = ("The calculation along a path of high-symmetry points -- "
                "the usual band-structure picture -- with the measurement "
                "under it along the <b>same</b> path. That picture is "
                "assembled from the selected maps: each in-plane segment "
                "from the map measured on its k_z plane (Γ-M-K-Γ from the "
                "zone-centre map, A-L-H-A from the zone-top one), Γ-A from a "
                "k_z scan. Each map's k_z comes from its photon energy and "
                "V0; correct it in the table if needed, untick a map to leave "
                "it out, and give the azimuth of the crystal's a axis on the "
                "map. Where a map does not reach a point, its symmetric "
                "equivalent is used.<br>The calculation is drawn as Z x E + "
                "shift (step 3 fits them); choose a weight to draw fat bands.")
        return Step("Along a path", text, page, enter=self._enter_path)

    def _enter_path(self):
        self._fill_sources()
        self.draw_path()

    def _maps(self):
        return [(n, d) for n, d in self.datasets if d.kind in ("k_map", "kz_map_k")]

    def _default_kz(self, data):
        from tools.spectra.corelevel import photon_energy_of, work_function_of

        hv = photon_energy_of(data.scan.info or {})
        if not hv:
            return 0.0
        return float(B.kz_of(0.0, 0.0, hv, self.v0_box.value(),
                             work_function_of(data.scan.info or {})))

    def _fill_sources(self):
        """One row per selected map; a k-map's k_z from its photon energy
        and V0, unless typed over."""
        maps = self._maps()
        self.source_table.blockSignals(True)
        if self.source_table.rowCount() != len(maps):
            self.source_table.setRowCount(len(maps))
            for row, (name, data) in enumerate(maps):
                item = QTableWidgetItem(name)
                item.setFlags(item.flags() & ~Qt.ItemIsEditable)
                self.source_table.setItem(row, 0, item)
                use = QTableWidgetItem("")
                use.setFlags((use.flags() | Qt.ItemIsUserCheckable)
                             & ~Qt.ItemIsEditable)
                use.setCheckState(Qt.Checked)
                self.source_table.setItem(row, 1, use)
                self.source_table.setItem(row, 2, QTableWidgetItem("0"))
                self.source_table.setItem(row, 3, QTableWidgetItem("0"))
        for row, (name, data) in enumerate(maps):
            if data.kind == "k_map" and name not in self._kz_typed:
                self.source_table.item(row, 2).setText(f"{self._default_kz(data):.4f}")
        self.source_table.blockSignals(False)

    def _source_edited(self, item):
        if item.column() == 2:
            self._kz_typed.add(self.source_table.item(item.row(), 0).text())
        self.draw_path()

    def set_source(self, name, *, use=None, value=None, azimuth=None):
        """Change one row of the map table (for scripts and tests)."""
        self._fill_sources()
        for row in range(self.source_table.rowCount()):
            if self.source_table.item(row, 0).text() != name:
                continue
            self.source_table.blockSignals(True)
            if use is not None:
                self.source_table.item(row, 1).setCheckState(
                    Qt.Checked if use else Qt.Unchecked)
            if value is not None:
                self._kz_typed.add(name)
                self.source_table.item(row, 2).setText(f"{value:.6g}")
            if azimuth is not None:
                self.source_table.item(row, 3).setText(f"{azimuth:.6g}")
            self.source_table.blockSignals(False)
        self.draw_path()

    def _sources(self):
        out = []
        maps = dict(self._maps())
        for row in range(self.source_table.rowCount()):
            name = self.source_table.item(row, 0).text()
            if (self.source_table.item(row, 1).checkState() != Qt.Checked
                    or name not in maps):
                continue
            data = maps[name]
            try:
                value = float(self.source_table.item(row, 2).text())
                azimuth = float(self.source_table.item(row, 3).text())
            except (TypeError, ValueError, AttributeError):
                continue
            scan = data.scan
            energy = _relative_energy(scan.z, scan.info)
            if data.kind == "k_map":
                out.append(C.PlaneMap(np.asarray(scan.value, dtype=float),
                                      np.asarray(scan.x, float),
                                      np.asarray(scan.k, float), energy,
                                      kz=value, rotation_deg=azimuth, name=name))
            else:
                out.append(C.KzMap(np.asarray(scan.value, dtype=float),
                                   np.asarray(scan.x, float),
                                   np.asarray(scan.k, float), energy,
                                   direction_deg=value + azimuth, name=name))
        return out

    # -- the bands along the path ----------------------------------------
    def path_k(self):
        """``(k (N, 3) in the calculation's frame, k in the zone's
        standard frame, distance, ticks)`` along the path typed."""
        from tools.crystal.lattice import primitive_vectors, reciprocal_vectors

        k_std, distance, ticks = B.named_path(self.params, self.path_edit.text(),
                                              step=0.008)
        b_std = reciprocal_vectors(primitive_vectors(self.params))
        fractional = k_std @ np.linalg.inv(b_std)
        return fractional @ self.model.reciprocal, k_std, distance, ticks

    def path_bands(self):
        """``(distance, energies (bands, N) aligned, ticks, weights or
        None)`` of the calculation along the path."""
        model = self.model
        weight = self.weight_combo.currentText() if self.weight_combo.currentIndex() > 0 else None
        if isinstance(model, B.BandPath):
            energies = self.alignment.apply(model.relative())
            weights = model.weights.get(weight) if weight else None
            return np.asarray(model.distance, float), energies, list(model.labels), weights
        k, _k_std, distance, ticks = self.path_k()
        energies = self.alignment.apply(model.energies(k))
        weights = model.weights_at(k, weight) if weight else None
        return distance, energies, ticks, weights

    def path_measurement(self, energy=None):
        """The measured intensity along the path ``(N, n_E)`` and its
        energy axis, or ``(None, None)``."""
        sources = self._sources()
        if not sources or isinstance(self.model, B.BandPath):
            return None, None, None
        if energy is None:
            steps = [abs(float(np.median(np.diff(s.energy)))) for s in sources
                     if len(s.energy) > 1]
            step = min(steps) if steps else 0.01
            lo = max(self.e_lo.value(), min(float(np.min(s.energy)) for s in sources))
            hi = min(self.e_hi.value(), max(float(np.max(s.energy)) for s in sources))
            energy = np.arange(lo, hi + 0.5 * step, step)
        _k, k_std, distance, _ticks = self.path_k()
        fold = self.fold_box.value() or C.symmetry_fold(self.params.space_group)
        from tools.crystal.lattice import primitive_vectors, reciprocal_vectors
        b = reciprocal_vectors(primitive_vectors(self.params))
        period = float(np.linalg.norm(b[2])) if abs(b[2][2]) > 1e-9 else None
        image, owner = C.experiment_along_path(k_std, energy, sources, fold=fold,
                                               kz_period=period)
        return distance, energy, image if (owner >= 0).any() else None

    def draw_path(self):
        if self.model is None:
            return
        view = self.path_view
        for item in self._items:
            view.view.removeItem(item)
        self._items = []
        try:
            distance, energies, ticks, weights = self.path_bands()
        except ValueError as exc:
            self.path_label.setText(f"<b>Path:</b> {exc}")
            return
        measured = None
        if self.path_measured.isChecked():
            try:
                _d, energy, measured = self.path_measurement()
            except Exception as exc:                          # noqa: BLE001
                self.path_label.setText(f"<b>Measurement:</b> {exc}")
        if measured is not None:
            show_frame(view, np.nan_to_num(measured, nan=np.nanmin(measured)),
                       distance, energy, auto_levels=True)
            apply_colormap(view, self.colormap, False)
            view.getImageItem().setVisible(True)
        elif view.image is not None:
            view.getImageItem().setVisible(False)
        colour = self.line_colour.currentText()
        window = (self.e_lo.value(), self.e_hi.value())
        for band, xs, ys in C.band_segments(distance, energies, ticks, window):
            item = pg.PlotDataItem(xs, ys, pen=pg.mkPen(colour, width=1.5))
            view.view.addItem(item)
            self._items.append(item)
        if weights is not None:
            w = np.clip(np.asarray(weights, float), 0, 1)
            for band in range(energies.shape[0]):
                keep = (w[band] > 0.05) & (energies[band] >= window[0]) & (energies[band] <= window[1])
                if keep.any():
                    dots = pg.ScatterPlotItem(distance[keep], energies[band][keep],
                                              size=2 + 10 * w[band][keep],
                                              pen=None,
                                              brush=pg.mkBrush(WEIGHT_COLOURS[0] + "99"))
                    view.view.addItem(dots)
                    self._items.append(dots)
        for d, _label in ticks:
            line = pg.InfiniteLine(d, angle=90, pen=pg.mkPen("#888888", width=0.8))
            view.view.addItem(line)
            self._items.append(line)
        view.view.getAxis("bottom").setTicks([[(d, l) for d, l in ticks], []])
        view.view.setXRange(float(distance[0]), float(distance[-1]), padding=0)
        view.view.setYRange(*window, padding=0)
        text = f"{energies.shape[0]} bands along {len(distance)} points"
        if measured is not None:
            covered = np.isfinite(measured).any(axis=1).mean() * 100
            text += f"; the maps cover {covered:.0f} % of the path"
        if self.alignment.z != 1.0 or self.alignment.shift != 0.0:
            text += (f"; drawn as {self.alignment.z:.3g} x E + "
                     f"{self.alignment.shift * 1000:.0f} meV")
        self.path_label.setText(text + ".")

    # ==================================================================
    # Step 3: compare
    # ==================================================================
    def _step_compare(self):
        page = QWidget()
        column = QHBoxLayout(page)
        left = QVBoxLayout()
        form = QFormLayout()
        self.compare_combo = QComboBox()
        self.compare_combo.currentIndexChanged.connect(lambda *_: self._compare_changed())
        form.addRow("Measurement", self.compare_combo)
        self.direction_box = _spin(0.0, -360.0, 360.0, 2, " °")
        self.direction_box.setToolTip("Direction of a cut from the crystal's "
                                      "k_x (its a* axis): 0 along Γ-K of a "
                                      "hexagonal zone, 90 along Γ-M. For a "
                                      "map, the crystal's azimuth on it.")
        self.direction_box.valueChanged.connect(lambda *_: self.draw_compare())
        form.addRow("Direction / azimuth", self.direction_box)
        self.offset_box = _spin(0.0, -10.0, 10.0, 4, " Å⁻¹")
        self.offset_box.setToolTip("Where the cut's k = 0 lies from Γ along "
                                   "the cut (if it was not zeroed at Γ).")
        self.offset_box.valueChanged.connect(lambda *_: self.draw_compare())
        form.addRow("k offset", self.offset_box)
        kz_row = QHBoxLayout()
        self.kz_sphere = QRadioButton("on the photon energy's k_z sphere")
        self.kz_plane_mode = QRadioButton("on the plane k_z =")
        self.kz_sphere.setChecked(True)
        self.kz_value = _spin(0.0, -50.0, 50.0, 4, " Å⁻¹")
        for widget in (self.kz_sphere, self.kz_plane_mode):
            widget.toggled.connect(lambda *_: self.draw_compare())
        self.kz_value.valueChanged.connect(lambda *_: self.draw_compare())
        kz_row.addWidget(self.kz_sphere)
        kz_row.addWidget(self.kz_plane_mode)
        kz_row.addWidget(self.kz_value)
        form.addRow("Calculation", kz_row)
        sphere = QHBoxLayout()
        self.hv_box = _spin(60.0, 1.0, 10000.0, 2, " eV")
        self.w_box = _spin(4.4, 0.0, 10.0, 3, " eV")
        for label, box in (("hν", self.hv_box), ("W", self.w_box)):
            sphere.addWidget(QLabel(label))
            sphere.addWidget(box)
            box.valueChanged.connect(lambda *_: self.draw_compare())
        form.addRow("Photon energy", sphere)
        self.energy_slider = QSlider(Qt.Horizontal)
        self.energy_slider.valueChanged.connect(lambda *_: self.draw_compare())
        self.energy_text = QLabel("")
        energy_row = QHBoxLayout()
        energy_row.addWidget(self.energy_slider, 1)
        energy_row.addWidget(self.energy_text)
        form.addRow("Map energy", energy_row)
        self.integrate_box = _spin(0.05, 0.0, 2.0, 3, " eV", 0.01)
        self.integrate_box.setToolTip("Energy window summed for a map's "
                                      "constant-energy picture.")
        self.integrate_box.valueChanged.connect(lambda *_: self.draw_compare())
        form.addRow("Integrate", self.integrate_box)
        left.addLayout(form)

        fit = QGroupBox("Z and shift: measured E = Z x calculated E + shift")
        f = QFormLayout(fit)
        self.z_box = _spin(1.0, 0.05, 5.0, 4, "", 0.01)
        self.shift_box = _spin(0.0, -5.0, 5.0, 4, " eV", 0.005)
        for box in (self.z_box, self.shift_box):
            box.valueChanged.connect(lambda *_: self._alignment_typed())
        f.addRow("Z (1/mass enhancement)", self.z_box)
        f.addRow("Shift", self.shift_box)
        self.fit_z = QCheckBox("Fit Z (otherwise only the shift)")
        self.fit_z.setChecked(True)
        f.addRow(self.fit_z)
        fit_button = QPushButton("Fit to the cut's bands")
        fit_button.setToolTip("Trace the measured bands (ridge detection) and "
                              "find the Z and shift that put the calculation "
                              "on them, each point matched to its nearest "
                              "band.")
        fit_button.clicked.connect(lambda: self.fit_alignment())
        f.addRow(fit_button)
        left.addWidget(fit)
        self.compare_label = QLabel("")
        self.compare_label.setWordWrap(True)
        left.addWidget(self.compare_label)
        left.addStretch(1)
        column.addLayout(left, 2)
        self.compare_view = plain_image_view("k (Å⁻¹)", "E - E_F (eV)")
        column.addWidget(self.compare_view, 3)
        self._compare_items = []
        self._fit_points = None
        text = ("One measurement at a time, with the calculation drawn where "
                "it was really taken. A <b>cut</b> (converted to k) gets the "
                "bands along its direction -- on a fixed k_z plane, or on the "
                "curved surface its photon energy reaches (k_z changes along "
                "the cut and with the energy). A <b>map</b> gets the "
                "calculated contours at the energy chosen.<br><b>Fit to the "
                "cut's bands</b> traces the measured bands and finds Z (the "
                "band-width renormalisation: 1/Z is the mass enhancement) and "
                "a rigid shift; they apply everywhere, in every step and "
                "output.")
        return Step("Compare", text, page, enter=self._enter_compare)

    def _fill_compare_list(self):
        current = self.compare_combo.currentText()
        self.compare_combo.blockSignals(True)
        self.compare_combo.clear()
        self.compare_combo.addItems([n for n, d in self.datasets
                                     if d.kind in ("cut", "k_map")])
        if current:
            self.compare_combo.setCurrentText(current)
        self.compare_combo.blockSignals(False)

    def _enter_compare(self):
        self._fill_compare_list()
        self._compare_changed()

    def _compare_data(self):
        name = self.compare_combo.currentText()
        return next((d for n, d in self.datasets if n == name), None)

    def _compare_changed(self):
        from tools.spectra.corelevel import photon_energy_of, work_function_of

        data = self._compare_data()
        if data is None:
            return
        info = data.scan.info or {}
        hv = photon_energy_of(info)
        for box, value in ((self.hv_box, hv), (self.w_box, work_function_of(info))):
            if value:
                box.blockSignals(True)
                box.setValue(float(value))
                box.blockSignals(False)
        is_map = data.kind == "k_map"
        self.energy_slider.setEnabled(is_map)
        self.integrate_box.setEnabled(is_map)
        self.offset_box.setEnabled(not is_map)
        if is_map:
            energy = _relative_energy(data.scan.z, info)
            self.energy_slider.blockSignals(True)
            self.energy_slider.setRange(0, len(energy) - 1)
            self.energy_slider.setValue(int(np.argmin(np.abs(energy))))
            self.energy_slider.blockSignals(False)
        self._fit_points = None
        self.draw_compare()

    def _alignment_typed(self):
        self.alignment = B.Alignment(self.z_box.value(), self.shift_box.value())
        self.draw_compare()

    def set_alignment(self, alignment):
        self.alignment = alignment
        for box, value in ((self.z_box, alignment.z), (self.shift_box, alignment.shift)):
            box.blockSignals(True)
            box.setValue(value)
            box.blockSignals(False)
        self.draw_compare()

    def _kz_options(self):
        if self.kz_plane_mode.isChecked():
            return {"kz": self.kz_value.value()}
        return {"hv": self.hv_box.value(), "inner_potential": self.v0_box.value(),
                "work_function": self.w_box.value()}

    def cut_arrays(self, data=None):
        """``(k, energy, values[k, E])`` of a cut, E relative to E_F."""
        data = data or self._compare_data()
        scan = data.scan
        label = str((scan.labels or {}).get("x", ""))
        if "deg" in label.lower() or "angle" in label.lower():
            raise ValueError("the cut is in angle: convert it to k first "
                             "(the cut viewer's k conversion)")
        return (np.asarray(scan.x, float),
                _relative_energy(scan.y, scan.info),
                np.asarray(scan.value, float))

    def cut_bands(self, k_axis, raw: bool = False):
        """The calculation along the compared cut at its momenta
        ``k_axis``: ``(bands, N)`` energies, aligned unless ``raw``."""
        if isinstance(self.model, B.BandPath):
            raise ValueError("a band structure along a path cannot be cut "
                             "along a measured direction")
        k = np.asarray(k_axis, float) + self.offset_box.value()
        energies = B.bands_along_cut(self.model, k, self.direction_box.value(),
                                     **self._kz_options())
        return energies if raw else self.alignment.apply(energies)

    def map_plane(self, data=None, energy=None):
        """``(kx, ky, plane, energy)`` of a k-map at ``energy`` (default:
        the slider's), summed over the integration window."""
        data = data or self._compare_data()
        scan = data.scan
        axis = _relative_energy(scan.z, scan.info)
        if energy is None:
            energy = float(axis[self.energy_slider.value()])
        half = 0.5 * self.integrate_box.value()
        keep = np.abs(axis - energy) <= max(half, 0.5 * abs(axis[1] - axis[0]))
        with np.errstate(all="ignore"), warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)      # all-NaN corners
            plane = np.nanmean(np.asarray(scan.value, float)[:, :, keep], axis=2)
        return np.asarray(scan.x, float), np.asarray(scan.k, float), plane, energy

    def map_contours(self, kx, ky, energy):
        """The calculation's contours at measured ``energy`` on the map's
        grid: ``[(band, xs, ys)]`` (the inverse alignment taken first)."""
        if isinstance(self.model, B.BandPath):
            raise ValueError("a band structure along a path has no contours")
        calc_energy = (energy - self.alignment.shift) / self.alignment.z
        options = self._kz_options()
        return B.contours(self.model, np.linspace(kx.min(), kx.max(), 161),
                          np.linspace(ky.min(), ky.max(), 161), calc_energy,
                          rotation_deg=self.direction_box.value(), **options)

    def contour_points_k(self, xs, ys, energy):
        """3-D k (in the calculation's frame) of contour points drawn on a
        map's (k_x, k_y) at measured ``energy`` -- for their weights."""
        c, s_ = (np.cos(np.deg2rad(-self.direction_box.value())),
                 np.sin(np.deg2rad(-self.direction_box.value())))
        xs, ys = np.asarray(xs, float), np.asarray(ys, float)
        inplane = np.column_stack([c * xs - s_ * ys, s_ * xs + c * ys])
        calc_energy = (energy - self.alignment.shift) / self.alignment.z
        options = self._kz_options()
        if "kz" in options:
            return np.column_stack([inplane, np.full(len(xs), options["kz"])])
        return B.kz_on_sphere(inplane, calc_energy, options["hv"],
                              options["inner_potential"], options["work_function"])

    def zone_overlays(self, x_range, y_range, colour="#d62728", style="dashed",
                      width=0.6):
        """The in-plane Brillouin zones (the calculation's lattice, turned
        by the azimuth) over a (k_x, k_y) range, as overlays."""
        from tools.crystal.bz2d import tile_2d, wigner_seitz_cell_2d
        from tools.output import figure as F

        if isinstance(self.model, B.BandPath):
            return []
        b = np.asarray(self.model.reciprocal, float)[:2, :2]
        angle = np.deg2rad(self.direction_box.value())
        turn = np.array([[np.cos(angle), -np.sin(angle)],
                         [np.sin(angle), np.cos(angle)]])
        g1, g2 = turn @ b[0], turn @ b[1]
        cell = wigner_seitz_cell_2d(g1, g2)
        # each edge once: neighbouring zones share their edges, and the same
        # dashed line drawn twice, offset by its own dash, comes out solid
        edges = {}
        for tile in tile_2d(cell, g1, g2, x_range, y_range, margin=1):
            for p, q in zip(tile[:-1], tile[1:]):
                key = tuple(sorted((tuple(np.round(p, 6)), tuple(np.round(q, 6)))))
                edges.setdefault(key, (p, q))
        return [F.Overlay(xs=[float(p[0]), float(q[0])], ys=[float(p[1]), float(q[1])],
                          color=colour, width_pt=width, style=style,
                          name="Brillouin zone")
                for p, q in edges.values()]

    def contour_overlays(self, kx, ky, energy, colours=("#d62728", "#2ca02c"),
                         width=0.8, style="solid", threshold=0.5):
        """The calculated contours at ``energy`` as overlays: in the first
        colour, or -- with a weight chosen in step 2 -- split where the
        chosen orbital's weight crosses ``threshold``, the parts above it in
        the second colour (a surface or orbital character told apart)."""
        from tools.output import figure as F

        weight = (self.weight_combo.currentText()
                  if self.weight_combo.currentIndex() > 0 else None)
        out = []
        for band, xs, ys in self.map_contours(kx, ky, energy):
            xs, ys = np.asarray(xs, float), np.asarray(ys, float)
            if weight is None:
                pieces = [(0, np.arange(xs.size))]
            else:
                w = self.model.weights_at(self.contour_points_k(xs, ys, energy),
                                          weight)[band]
                high = (w >= threshold).astype(int)
                cuts = np.flatnonzero(np.diff(high)) + 1
                pieces = []
                for run in np.split(np.arange(xs.size), cuts):
                    if run.size:
                        # one point into the next run, so the pieces join
                        end = min(run[-1] + 2, xs.size)
                        pieces.append((int(high[run[0]]), np.arange(run[0], end)))
            for which, idx in pieces:
                if idx.size < 2:
                    continue
                out.append(F.Overlay(xs=list(map(float, xs[idx])),
                                     ys=list(map(float, ys[idx])),
                                     color=colours[which], width_pt=width,
                                     style=style,
                                     name=f"band {band + 1}"
                                          + (f" ({weight})" if which else "")))
        return out

    def panel_contours(self, data=None, energy=None, *, zone=True,
                       colours=("#d62728", "#2ca02c")):
        """The calculated contours alone at ``energy``, over the compared
        map's (k_x, k_y) range, with the zones dashed: the calculation's
        half of a constant-energy comparison."""
        from ui.viewers.curves import curve_panel

        data = data or self._compare_data()
        kx, ky, _plane, energy = self.map_plane(data, energy)
        panel = curve_panel([float(kx.min()), float(kx.max())], [],
                            x_label="k_x (Å⁻¹)", y_label="k_y (Å⁻¹)",
                            name=f"calculation at {energy:+.2f} eV")
        panel.data.x = np.array([float(kx.min()), float(kx.max())])
        panel.data.y = np.array([float(ky.min()), float(ky.max())])
        panel.x_range = (float(kx.min()), float(kx.max()))
        panel.y_range = (float(ky.min()), float(ky.max()))
        if zone:
            panel.overlays += self.zone_overlays(panel.x_range, panel.y_range,
                                                 colour="#d62728")
        panel.overlays += self.contour_overlays(kx, ky, energy, colours)
        panel.data.sources = self.provenance([])
        return panel

    def draw_compare(self):
        view = self.compare_view
        for item in self._compare_items:
            view.view.removeItem(item)
        self._compare_items = []
        data = self._compare_data()
        if data is None or self.model is None:
            return
        try:
            if data.kind == "cut":
                k, energy, values = self.cut_arrays(data)
                show_frame(view, values, k, energy, auto_levels=True)
                view.view.setLabel("bottom", "k (Å⁻¹)")
                view.view.setLabel("left", "E - E_F (eV)")
                for band in self.cut_bands(k):
                    item = pg.PlotDataItem(k, band, pen=pg.mkPen("#e0302a", width=1.5))
                    view.view.addItem(item)
                    self._compare_items.append(item)
                if self._fit_points is not None:
                    dots = pg.ScatterPlotItem(*self._fit_points, size=3, pen=None,
                                              brush=pg.mkBrush("#2ca02c"))
                    view.view.addItem(dots)
                    self._compare_items.append(dots)
                fit_frame_view(view, k, energy)
            else:
                kx, ky, plane, energy = self.map_plane(data)
                self.energy_text.setText(f"{energy:+.3f} eV")
                show_frame(view, plane, kx, ky, auto_levels=True)
                view.view.setLabel("bottom", "k_x (Å⁻¹)")
                view.view.setLabel("left", "k_y (Å⁻¹)")
                for _band, xs, ys in self.map_contours(kx, ky, energy):
                    item = pg.PlotDataItem(xs, ys, pen=pg.mkPen("#e0302a", width=1.5))
                    view.view.addItem(item)
                    self._compare_items.append(item)
                fit_frame_view(view, kx, ky)
            apply_colormap(view, self.colormap, False)
            self.compare_label.setText(
                f"Z = {self.alignment.z:.3f}, shift = "
                f"{self.alignment.shift * 1000:+.1f} meV"
                + (f" (± {self.alignment.z_err:.3f}, ± "
                   f"{self.alignment.shift_err * 1000:.1f} meV; rms "
                   f"{self.alignment.rms * 1000:.1f} meV over "
                   f"{self.alignment.points} points)"
                   if np.isfinite(self.alignment.rms) else "") + ".")
        except ValueError as exc:
            self.compare_label.setText(f"<b>Cannot draw it:</b> {exc}")

    def fit_alignment(self):
        data = self._compare_data()
        if data is None or data.kind != "cut" or self.model is None:
            self.compare_label.setText("Choose a cut to fit on.")
            return None
        k, energy, values = self.cut_arrays(data)
        ks, es = C.band_points(values, k, energy, energy_window=(-np.inf, -0.03))
        if ks.size < 5:
            self.compare_label.setText("Too few band points were traced.")
            return None
        order = np.argsort(ks)
        ks, es = ks[order], es[order]
        try:
            result = B.fit_alignment(ks, es, lambda kk: self.cut_bands(kk, raw=True),
                                     fit_z=self.fit_z.isChecked())
        except (ValueError, np.linalg.LinAlgError) as exc:
            self.compare_label.setText(f"<b>Could not fit:</b> {exc}")
            return None
        self._fit_points = (ks, es)
        self.set_alignment(result)
        return result

    # ==================================================================
    # Step 4: output
    # ==================================================================
    def _step_output(self):
        page = QWidget()
        column = QVBoxLayout(page)
        form = QFormLayout()
        self.width_box = _spin(0.03, 0.001, 1.0, 3, " eV", 0.005)
        self.width_box.setToolTip("Lorentzian half-width given to the "
                                  "calculated bands to make a spectrum of them.")
        form.addRow("Band width (HWHM)", self.width_box)
        column.addLayout(form)
        group = QGroupBox("Figure-composer panels")
        g = QVBoxLayout(group)
        self.panel_boxes = {}
        for key, text in (("path", "Measurement along the path, calculation over it"),
                          ("path_calc", "Calculation alone along the path"),
                          ("compare", "The compared cut or map, calculation over it")):
            box = QCheckBox(text)
            box.setChecked(True)
            g.addWidget(box)
            self.panel_boxes[key] = box
        series = QHBoxLayout()
        self.series_box = QCheckBox("Constant-energy series at")
        self.series_box.setToolTip(
            "For a compared map: at each energy, the measured map with the "
            "zones dashed (top row) and the calculated contours alone "
            "(bottom row) -- coloured by the weight chosen in step 2.")
        self.series_energies = QLineEdit("0, -0.3, -0.6, -0.9")
        self.series_energies.setToolTip("Energies (eV, E - E_F), separated "
                                        "by commas.")
        series.addWidget(self.series_box)
        series.addWidget(self.series_energies, 1)
        g.addLayout(series)
        compose = QPushButton("To figure composer...")
        compose.clicked.connect(lambda: self.to_figure_composer())
        g.addWidget(compose)
        column.addWidget(group)
        buttons = QHBoxLayout()
        spectrum = QPushButton("Calculation to list")
        spectrum.setToolTip("The calculated bands as a spectrum (Lorentzians "
                            "of the width set here), along the path -- or "
                            "along the compared cut, on its own k grid -- as "
                            "a new cut in the list.")
        spectrum.clicked.connect(lambda: self.spectrum_to_list())
        buttons.addWidget(spectrum)
        buttons.addStretch(1)
        column.addLayout(buttons)
        self.output_label = QLabel("")
        self.output_label.setWordWrap(True)
        column.addWidget(self.output_label)
        column.addStretch(1)
        text = ("The comparison as figures and data. <b>To figure "
                "composer</b> makes the panels ticked: the measurement along "
                "the path with the calculation over it, the calculation alone "
                "(the points' names on its axis), and the compared cut or "
                "map with its calculation -- arranged, labelled and exported "
                "there like any panel. <b>Calculation to list</b> makes the "
                "calculation a spectrum, so it opens in the viewers beside "
                "the data or is subtracted from it.")
        return Step("Output", text, page)

    # -- panels -----------------------------------------------------------
    def _overlays(self, distance, energies, ticks, window, colour, weights=None,
                  width=0.9):
        from tools.output import figure as F

        out = []
        for band, xs, ys in C.band_segments(distance, energies, ticks, window):
            out.append(F.Overlay(xs=list(map(float, xs)), ys=list(map(float, ys)),
                                 color=colour, width_pt=width,
                                 name=f"band {band + 1}"))
        if weights is not None:
            w = np.clip(np.asarray(weights, float), 0, 1)
            for band, row in enumerate(energies):
                for start in range(0, len(distance) - 1, 6):
                    sl = slice(start, start + 7)
                    mean = float(np.mean(w[band, sl]))
                    if mean < 0.08 or not np.all(np.isfinite(row[sl])):
                        continue
                    out.append(F.Overlay(xs=list(map(float, distance[sl])),
                                         ys=list(map(float, row[sl])),
                                         color=WEIGHT_COLOURS[0],
                                         width_pt=0.5 + 3.5 * mean,
                                         name=f"weight {band + 1}"))
        return out

    def _ticks_and_lines(self, panel, ticks, colour="#000000"):
        from tools.output import figure as F

        panel.x_tick_labels = [(float(d), str(l)) for d, l in ticks]
        for d, _l in ticks[1:-1]:
            panel.lines.append(F.Line(value=float(d), axis="x", color=colour,
                                      width_pt=0.5, style="solid"))

    def panel_path(self, colormap=None):
        """The measurement along the path with the calculation over it."""
        from ui.panels.figure import panel_from_arrays

        distance, energies, ticks, weights = self.path_bands()
        _d, energy, measured = self.path_measurement()
        if measured is None:
            raise ValueError("no map covers the path: select k-maps (and a "
                             "k_z map) and tick them in step 2")
        used = [d for n, d in self._maps() if n in {src.name for src in self._sources()}]
        panel = panel_from_arrays(measured, distance, energy, "",
                                  "E - E_F (eV)", "measured along the path",
                                  colormap=colormap or self.colormap,
                                  sources=self.provenance(used))
        window = (self.e_lo.value(), self.e_hi.value())
        panel.y_range = window
        panel.overlays = self._overlays(distance, energies, ticks, window,
                                        self.line_colour.currentText(), weights)
        self._ticks_and_lines(panel, ticks, "#ffffff" if self.line_colour.currentText() == "white" else "#000000")
        return panel

    def panel_path_calc(self):
        """The calculation alone, the points' names on its axis."""
        from ui.viewers.curves import curve_panel

        distance, energies, ticks, weights = self.path_bands()
        window = (self.e_lo.value(), self.e_hi.value())
        panel = curve_panel(distance, [], x_label="", y_label="E - E_F (eV)",
                            name="calculation")
        panel.data.y = np.array(window, dtype=float)
        panel.data.x = np.array([float(distance[0]), float(distance[-1])])
        panel.data.sources = self.provenance([])
        panel.y_range = window
        panel.overlays = self._overlays(distance, energies, ticks, window,
                                        "#000000", weights)
        self._ticks_and_lines(panel, ticks, "#888888")
        return panel

    def panel_compare(self, data=None, colormap=None, colour="#ff2020",
                      energy=None, *, contours=True, zone=False):
        """The compared cut (calculation over it) or map (contours)."""
        from tools.output import figure as F
        from ui.panels.figure import panel_from_arrays

        data = data or self._compare_data()
        name = next((n for n, d in self.datasets if d is data), "measurement")
        if data.kind == "cut":
            k, e_axis, values = self.cut_arrays(data)
            panel = panel_from_arrays(values, k, e_axis, "k (Å⁻¹)", "E - E_F (eV)",
                                      name, colormap=colormap or self.colormap,
                                      sources=self.provenance([data]))
            window = (float(e_axis.min()), float(e_axis.max()))
            for band, row in enumerate(self.cut_bands(k)):
                ys = np.where((row > window[0]) & (row < window[1]), row, np.nan)
                pieces = np.split(np.arange(k.size), np.flatnonzero(np.isnan(ys)))
                for piece in pieces:
                    piece = piece[np.isfinite(ys[piece])]
                    if piece.size >= 2:
                        panel.overlays.append(F.Overlay(
                            xs=list(map(float, k[piece])),
                            ys=list(map(float, ys[piece])), color=colour,
                            width_pt=0.8, style="dashed", name=f"band {band + 1}"))
            return panel
        kx, ky, plane, energy = self.map_plane(data, energy)
        panel = panel_from_arrays(plane, kx, ky, "k_x (Å⁻¹)", "k_y (Å⁻¹)",
                                  f"{name} at {energy:+.2f} eV",
                                  colormap=colormap or self.colormap,
                                  sources=self.provenance([data]))
        if zone:
            panel.overlays += self.zone_overlays(
                (float(kx.min()), float(kx.max())),
                (float(ky.min()), float(ky.max())), colour=colour)
        if contours:
            panel.overlays += self.contour_overlays(kx, ky, energy,
                                                    (colour, "#2ca02c"))
        return panel

    def provenance(self, datasets) -> list:
        """What a panel was made from: the measurements' original files,
        and the calculation's."""
        from tools.output.figfile import source_record
        from ui.panels.figure import _dataset_sources

        out = []
        for data in datasets:
            out.extend(_dataset_sources(data))
        if self.model_path:
            out.append(source_record(self.model_path, kind="calculation",
                                     companions=bandfiles.companion_files(self.model_path),
                                     note=f"drawn as {self.alignment.z:.4g} x E "
                                          f"+ {self.alignment.shift:.4g} eV"))
        return out

    def panels(self):
        out, problems = [], []
        makers = (("path", self.panel_path), ("path_calc", self.panel_path_calc),
                  ("compare", self.panel_compare))
        for key, maker in makers:
            if not self.panel_boxes[key].isChecked():
                continue
            try:
                if key == "compare" and self._compare_data() is None:
                    continue
                out.append(maker())
            except ValueError as exc:
                problems.append(str(exc))
        if self.series_box.isChecked():
            try:
                out += self.series_panels()
            except ValueError as exc:
                problems.append(str(exc))
        return out, problems

    def series_energies_list(self):
        values = []
        for part in self.series_energies.text().replace(";", ",").split(","):
            try:
                values.append(float(part))
            except ValueError:
                continue
        return values

    def series_panels(self, data=None, energies=None):
        """The measured maps at each energy (zones dashed), then the
        calculated contours alone at the same energies: a 2-row series."""
        data = data or self._compare_data()
        if data is None or data.kind != "k_map":
            raise ValueError("a constant-energy series needs a k-map chosen in "
                             "step 3")
        energies = self.series_energies_list() if energies is None else energies
        measured = [self.panel_compare(data, energy=e, contours=False, zone=True)
                    for e in energies]
        calculated = [self.panel_contours(data, energy=e) for e in energies]
        for panel, e in zip(measured, energies):
            panel.title = "E_F" if abs(e) < 1e-9 else f"E - E_F = {e:g} eV"
        return measured + calculated

    def to_figure_composer(self, target=None, ask=True):
        from ui.panels.figure import FigureWindow, open_figure_windows
        from ui.viewers.windows import keep_alive

        if self.model is None:
            return None
        panels, problems = self.panels()
        if not panels:
            self.output_label.setText("Nothing to add" + (
                ": " + "; ".join(problems) if problems else "."))
            return None
        if target is None and ask:
            windows = open_figure_windows()
            if windows:
                names = ["A new figure"] + [w.windowTitle() for w in windows]
                choice, ok = QInputDialog.getItem(
                    self, "To figure composer", "Add the panels to:", names, 0,
                    False)
                if not ok:
                    return None
                index = names.index(choice)
                target = windows[index - 1] if index > 0 else None
        if target is None:
            target = FigureWindow(panels, "Figure — calculated bands", parent=None)
            keep_alive(target)
            target.show()
        else:
            target.add_panels(panels)
            target.raise_()
        text = f"Added {len(panels)} panel(s) to “{target.windowTitle()}”."
        if problems:
            text += " Left out: " + "; ".join(problems)
        self.output_label.setText(text)
        return target

    # -- the spectrum -------------------------------------------------------
    def build_spectrum(self):
        """The calculation as a cut: along the compared cut's k if one is
        chosen, else along the path."""
        width = self.width_box.value()
        data = self._compare_data()
        parameters = {"calculation": os.path.basename(self.model_path),
                      "z": self.alignment.z, "shift_eV": self.alignment.shift,
                      "width_eV": width}
        if data is not None and data.kind == "cut" and not isinstance(self.model, B.BandPath):
            k, energy, _values = self.cut_arrays(data)
            energies = self.cut_bands(k)
            x_label = "k (Å⁻¹)"
            base = self.compare_combo.currentText()
            parameters["direction_deg"] = self.direction_box.value()
        else:
            k, energies, ticks, _w = self.path_bands()
            energy = np.linspace(self.e_lo.value(), self.e_hi.value(), 301)
            x_label = "Path (Å⁻¹)"
            base = "path"
            parameters["path"] = " ".join(f"{l}@{d:.4f}" for d, l in ticks)
        values = B.spectral_function(energies, energy, width)
        taken = set(self._existing())
        name, n = f"{base} calculation", 2
        while name in taken:
            name, n = f"{base} calculation {n}", n + 1
        return MemoryData("cut", (np.asarray(k, float), np.asarray(energy, float)),
                          values.astype(np.float32),
                          {"x": x_label, "y": "E - E_F (eV)"},
                          source_label=name, parameters=parameters, prefix="calc",
                          source_path=self.model_path)

    def spectrum_to_list(self):
        if self.model is None:
            return None
        try:
            result = self.build_spectrum()
        except ValueError as exc:
            self.output_label.setText(f"<b>Could not make it:</b> {exc}")
            return None
        self.datasetsCreated.emit([result])
        self.output_label.setText(f"Added {result.source_label} to the list.")
        return result
