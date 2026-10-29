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
   Wannier90 Hamiltonian (``_hr.dat``), bands along a path (Wannier90
   ``_band.dat``, Wien2k ``.spaghetti_ene``, VASP ``EIGENVAL`` with its
   ``PROCAR``), or what the FPLO + Wannier tools write (``band.dat`` along
   a path or on a 3-D grid, ``Eigenvalue.dat``, with their ``KPOINTS``).
   The crystal's space group names the high-symmetry points.
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
   calculated contours, and calculated contours alone (no map needed).

With no measurement selected, steps 2 and 3 can simply be paged through:
the calculation alone goes to the output all the same.

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
                             QInputDialog, QLabel, QLineEdit, QListWidget,
                             QListWidgetItem, QPushButton, QRadioButton,
                             QSlider, QSpinBox, QStackedWidget, QTableWidget,
                             QTableWidgetItem, QVBoxLayout, QWidget)

from arpes_viewer.tools.crystal import bandfiles
from arpes_viewer.tools.crystal import bands as B
from arpes_viewer.tools.crystal import compare as C
from arpes_viewer.tools.crystal import last_crystal
from arpes_viewer.tools.crystal.lattice import LatticeParams
from arpes_viewer.ui.common.screen import fits_screen
from arpes_viewer.ui.common.jobs import one_at_a_time
from arpes_viewer.ui.common.stepper import Step, Stepper, guard_window
from arpes_viewer.ui.viewers.widgets import (MemoryData, apply_colormap, fit_frame_view,
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

#: under the goal's button in step 1: which files serve it, in plain words
GOAL_HINTS = {
    "dispersion": (
        "Any calculation works. Usually it is the file the theory side plots "
        "their band structure from: <code>band.dat</code>, "
        "<code>*_band.dat</code>, <code>Eigenvalue.dat</code>, "
        "<code>case.spaghetti_ene</code>, VASP <code>EIGENVAL</code>. A file "
        "of the bands everywhere (below) works too, and can also be cut "
        "along your own measured direction."),
    "contours": (
        "Needs the bands <b>everywhere in the zone</b>: a band grid "
        "(<code>.bxsf</code>; a <code>band.dat</code> with columns k1 k2 k3 "
        "E1 E2 ...; a MATLAB <code>Ek_3D</code> <code>.mat</code>) or a "
        "Wannier Hamiltonian (<code>*_hr.dat</code>). A band structure along "
        "a path cannot give contours."),
    "fs": (
        "Needs the bands <b>everywhere in the zone</b>, as for contours: a "
        "band grid (<code>.bxsf</code>, a 3-D <code>band.dat</code>, a MATLAB "
        "<code>Ek_3D</code> <code>.mat</code>) or a Hamiltonian "
        "(<code>*_hr.dat</code>)."),
}

#: step 2's title and text, by the goal
VIEW_STEPS = {
    "dispersion": ("Along a path", (
        "The calculation along a path of high-symmetry points -- the usual "
        "band-structure picture -- with the measurement under it along the "
        "<b>same</b> path when maps are selected (each in-plane segment from "
        "the map measured on its k_z plane, Γ-A from a k_z scan; correct a "
        "map's k_z in the table, untick it to leave it out). "
        "<b>Highlight bands</b> draws chosen bands in colour, the rest grey. "
        "The calculation is drawn as Z x E + shift (step 3 fits them).")),
    "contours": ("Constant-energy contours", (
        "The calculation's contours at one energy on one plane of k-space "
        "-- a calculated Fermi-surface map. The plane: its normal (hkl), as "
        "in the Brillouin-zone tool, and where it sits along that normal (k⊥; "
        "<b>Through point</b> puts it through a high-symmetry point). A "
        "<b>k⊥ integration window</b> sums the bands over a range of k⊥ "
        "about it, as a map's finite k⊥ resolution does. <b>Show the plane "
        "in the Brillouin zone</b> says where the cut lies. Each band's "
        "colour is yours to choose (the same for several, to group pockets). "
        "Step 3 puts the contours over a measured k-map, step 4 makes "
        "figures and a calculated map of them.")),
    "fs": ("3D Fermi surface", (
        "The surfaces where the bands cross the energy chosen, inside the "
        "first Brillouin zone. Each band's colour (the same for several, to "
        "group electron and hole pockets), opacity and visibility are yours "
        "to set. Drag to turn it; Enter in the energy redraws it. Step 4 puts "
        "the picture in the figure composer (there is nothing to compare, so "
        "the Compare step is left out).")),
}


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
    from arpes_viewer.tools.spectra.corelevel import photon_energy_of, work_function_of

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
        self._last_file = ""
        self._folder = ""
        self._which = None
        self._fs_grid = None
        self._fs_drawn = None
        self._fs_energy_asked = None
        self.params = None
        self.alignment = B.Alignment()
        self._items = []
        self._kz_typed = set()

        layout = QVBoxLayout(self)
        self.stepper = Stepper([self._step_calculation(), self._step_view(),
                                self._step_compare(), self._step_output()],
                               self, finish_text="Close")
        layout.addWidget(self.stepper)
        # the cell and space group are the calculation's, never last time's
        self._not_remembered = ("compare_combo", "weight_combo", "fs_view_combo",
                                "space_group", "change_lattice")
        guard_window(self, self.stepper, "calculated_bands", "Calculated bands",
                     close_on_finish=True,
                     extra=lambda: {"last_file": self.model_path or self._last_file,
                                    "folder": self._folder})
        from arpes_viewer.ui.common import remember

        kept = remember.load("calculated_bands").get("_extra") or {}
        self._folder = kept.get("folder") or ""
        last = kept.get("last_file") or ""
        if last and os.path.isfile(last):
            self._last_file = last
            self.reopen_button.setText(f"Reopen {os.path.basename(last)}")
            self.reopen_button.setToolTip(f"The calculation opened last time:\n{last}")
            self.reopen_button.setVisible(True)
        self._lock_lattice()
        # what was selected says what is wanted
        kinds = {d.kind for _n, d in self.datasets}
        if "k_map" in kinds and "cut" not in kinds:
            self.set_goal("contours")
        elif kinds:
            self.set_goal("dispersion")
        else:
            self._goal_changed()

    # ==================================================================
    # Step 1: the calculation
    # ==================================================================
    def _step_calculation(self):
        page = QWidget()
        column = QVBoxLayout(page)

        goal = QGroupBox("What do you want to see?")
        g = QVBoxLayout(goal)
        self.goal_dispersion = QRadioButton(
            "E–k dispersion — the bands along a path, or along my measured cut")
        self.goal_contours = QRadioButton(
            "Constant-energy contours — to compare with a k-map (Fermi surface map)")
        self.goal_fs = QRadioButton("3D Fermi surface")
        self.goal_dispersion.setChecked(True)
        for button in (self.goal_dispersion, self.goal_contours, self.goal_fs):
            g.addWidget(button)
            button.toggled.connect(lambda on: on and self._goal_changed())
        self.goal_hint = QLabel("")
        self.goal_hint.setWordWrap(True)
        self.goal_hint.setTextFormat(Qt.RichText)
        self.goal_hint.setStyleSheet("color: #333; padding-left: 4px;")
        g.addWidget(self.goal_hint)
        column.addWidget(goal)

        row = QHBoxLayout()
        self.open_button = QPushButton("Open calculation...")
        self.open_button.setToolTip("Choose the file. The list shows first the "
                                    "kinds of file that can give what you want "
                                    "to see.")
        self.open_button.clicked.connect(lambda: self.open_file())
        row.addWidget(self.open_button)
        self.which_button = QPushButton("Which file? ...")
        self.which_button.setToolTip(
            "Not sure which file of the folder is the calculation? This reads "
            "every file in it and says in plain words what each is, what it "
            "can show, what is missing, and whether a standard file was only "
            "renamed or cut into pieces.")
        self.which_button.clicked.connect(lambda: self.which_file())
        row.addWidget(self.which_button)
        self.reopen_button = QPushButton("")
        self.reopen_button.setVisible(False)
        self.reopen_button.clicked.connect(lambda: self.open_file(self._last_file))
        row.addWidget(self.reopen_button)
        row.addStretch(1)
        column.addLayout(row)
        self.file_label = QLabel("No calculation yet.")
        self.file_label.setWordWrap(True)
        self.file_label.setTextFormat(Qt.RichText)
        column.addWidget(self.file_label)

        self.options_toggle = QPushButton("Reading options ▸")
        self.options_toggle.setCheckable(True)
        self.options_toggle.setAutoDefault(False)
        self.options_toggle.setFlat(True)
        self.options_toggle.setStyleSheet("text-align: left; color: #1f6f8b;")
        self.options_toggle.setToolTip("E_F, the k units of a .bxsf, which "
                                       "bands a big grid keeps -- usually "
                                       "nothing to change.")
        column.addWidget(self.options_toggle)
        options = QGroupBox("Reading it (usually nothing to change)")
        options.setVisible(False)
        self.options_toggle.toggled.connect(lambda on: (
            options.setVisible(on),
            self.options_toggle.setText("Reading options ▾" if on else
                                        "Reading options ▸")))
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
        keep = QHBoxLayout()
        self.keep_lo = _spin(-2.0, -50.0, 0.0, 2, " eV", 0.5)
        self.keep_hi = _spin(1.0, 0.0, 50.0, 2, " eV", 0.5)
        for box in (self.keep_lo, self.keep_hi):
            box.setToolTip("A 3-D band grid (band.dat of the FPLO tools, a "
                           "MATLAB Ek_3D) holds every band on a dense grid -- "
                           "gigabytes; only the bands that come within this "
                           "window of E_F are kept.")
        keep.addWidget(self.keep_lo)
        keep.addWidget(QLabel("to"))
        keep.addWidget(self.keep_hi)
        keep.addStretch(1)
        form.addRow("Big grids: keep bands from", keep)
        column.addWidget(options)

        crystal = QGroupBox("The crystal (names the high-symmetry points)")
        cform = QFormLayout(crystal)
        self.change_lattice = QCheckBox("Change the lattice parameters")
        self.change_lattice.setToolTip(
            "The cell comes with the calculation, so it is locked once read. "
            "Tick this only to change it on purpose (to scale the calculation "
            "to another lattice constant, say, or to name its points by "
            "another space group).")
        self.change_lattice.toggled.connect(lambda *_: self._lock_lattice())
        cform.addRow(self.change_lattice)
        cell_row = QHBoxLayout()
        self.cell_boxes = []
        for name, value, top, suffix in (("a", 3.0, 1000.0, " Å"), ("b", 3.0, 1000.0, " Å"),
                                         ("c", 3.0, 1000.0, " Å"), ("α", 90.0, 179.0, "°"),
                                         ("β", 90.0, 179.0, "°"), ("γ", 90.0, 179.0, "°")):
            box = _spin(value, 0.1 if suffix == " Å" else 1.0, top, 4, suffix)
            box.editingFinished.connect(lambda: self._cell_typed())
            cell_row.addWidget(QLabel(name))
            cell_row.addWidget(box)
            self.cell_boxes.append(box)
        cform.addRow("Cell", cell_row)
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
        self.take_button = take
        sg_row.addStretch(1)
        cform.addRow("Space group", sg_row)
        self.crystal_label = QLabel("")
        self.crystal_label.setWordWrap(True)
        cform.addRow(self.crystal_label)
        column.addWidget(crystal)
        column.addStretch(1)
        self.setAcceptDrops(True)
        text = ("First say what you want to see; then open the calculation "
                "-- the file list shows first the files that can give it. Not "
                "sure which file it is? <b>Which file?</b> reads the folder "
                "and tells. (A file can also be dropped on this window.) With "
                "no measurement selected, the next steps simply show the "
                "calculation, and step 4 makes its figures.")
        return Step("The calculation", text, page, ready=self._calculation_ready)

    def _settings_applied(self):
        """Settings from last time were put back: show what they say."""
        if hasattr(self, "view_pages"):
            self._goal_changed()

    # -- what the user wants ------------------------------------------------
    def goal(self) -> str:
        if self.goal_contours.isChecked():
            return "contours"
        if self.goal_fs.isChecked():
            return "fs"
        return "dispersion"

    def set_goal(self, goal: str):
        {"dispersion": self.goal_dispersion, "contours": self.goal_contours,
         "fs": self.goal_fs}[goal].setChecked(True)
        self._goal_changed()

    def _goal_changed(self):
        goal = self.goal()
        self.goal_hint.setText(GOAL_HINTS[goal])
        title, text = VIEW_STEPS[goal]
        if hasattr(self, "view_pages"):
            self.view_pages.setCurrentIndex(bandfiles.GOALS.index(goal))
            self.stepper.retitle(1, title, text)
            # a 3-D Fermi surface has nothing to be compared with
            self.stepper.set_skipped(2, goal == "fs")
            self._output_defaults()
            self.stepper.refresh()

    def _calculation_ready(self) -> str:
        if self.model is None:
            return "Open a calculation first."
        if isinstance(self.model, B.BandPath) and self.goal() != "dispersion":
            what = ("constant-energy contours" if self.goal() == "contours"
                    else "a 3D Fermi surface")
            return (f"This file holds the bands along a path only: it cannot "
                    f"give {what}. Open a file of the bands everywhere (a band "
                    "grid or a Hamiltonian -- “Which file?” helps), or choose "
                    "E–k dispersion.")
        return ""

    # -- opening a file -----------------------------------------------------
    def dragEnterEvent(self, event):
        if event.mimeData().hasUrls():
            event.acceptProposedAction()

    def dropEvent(self, event):
        urls = [u.toLocalFile() for u in event.mimeData().urls() if u.isLocalFile()]
        if urls:
            self.open_file(urls[0])

    def which_file(self, folder=None):
        """The "Which file?" helper on ``folder`` (the last one used)."""
        from arpes_viewer.ui.panels.bandfinder import WhichFileDialog

        dialog = WhichFileDialog(folder or self._folder, self.goal(), self)
        dialog.chosen.connect(lambda path: self.open_file(path))
        self._which = dialog
        dialog.show()
        if not dialog.folder:
            dialog.choose_folder()
        return dialog

    def open_file(self, path=None, **options):
        if path is None:
            filters = [f"{label} ({' '.join(patterns)})"
                       for label, patterns in bandfiles.formats_for(self.goal())]
            path, _ = QFileDialog.getOpenFileName(
                self, "Calculated bands", self._folder or "",
                ";;".join(filters + ["All files (*)"]))
            if not path:
                return None
        if not os.path.isfile(path):
            self.file_label.setText(f"<b>Not found:</b> {path}")
            return None
        self._folder = os.path.dirname(os.path.abspath(path))
        kwargs = dict(options)
        if self.fermi_box.isChecked():
            kwargs.setdefault("fermi", self.fermi_value.value())
        if path.lower().endswith(".bxsf"):
            kwargs.setdefault("k_scale", (1.0, 2 * np.pi, 1.0 / bandfiles.BOHR)[
                self.k_unit.currentIndex()])
        kwargs.setdefault("window", (self.keep_lo.value(), self.keep_hi.value()))
        crystal = last_crystal.lattice_params()

        def work(report):
            try:
                return bandfiles.read(path, progress=report, **kwargs), False
            except ValueError as exc:
                # no cell beside the file: the Brillouin-zone tool's crystal
                if "cell is needed" not in str(exc) or crystal is None \
                        or "lattice" in kwargs:
                    raise
                from arpes_viewer.tools.crystal.lattice import primitive_vectors

                return bandfiles.read(path, progress=report,
                                      lattice=primitive_vectors(crystal),
                                      **kwargs), True

        from arpes_viewer.ui.common.jobs import JobCancelled, run_blocking

        try:
            model, borrowed = run_blocking(
                self, f"Reading {os.path.basename(path)}", work)
        except JobCancelled:
            self.file_label.setText("Reading was cancelled.")
            return None
        except (OSError, ValueError, IndexError, MemoryError) as exc:
            self.file_label.setText(
                f"<b>Could not read {os.path.basename(path)}:</b> {exc}<br>"
                "<b>Which file?</b> says what each file of the folder is.")
            return None
        self.model_path = path
        self.set_model(model)
        if borrowed:
            self.file_label.setText(self.file_label.text()
                                    + " <i>The cell is the Brillouin-zone tool's "
                                      "crystal (none came with the file).</i>")
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
            self._show_cell(self.params)
            kind = ("a Hamiltonian, diagonalised at every k"
                    if isinstance(model, B.WannierModel)
                    else "a band grid, interpolated")
            count = model.band_count()
            numbers = getattr(model, "band_numbers", None)
            if numbers is not None and len(numbers):
                count = (f"{count} (bands {numbers[0]}-{numbers[-1]} of "
                         f"{model.total_bands}, those near E_F)")
        weights = list(getattr(model, "weight_names", []) or [])
        name = os.path.basename(self.model_path) or getattr(model, "name", "") \
            or "calculation"
        crossing = self.crossing_numbers()
        shows = ("E–k dispersion along its path only" if isinstance(model, B.BandPath)
                 else "E–k dispersion, constant-energy contours and the 3D "
                      "Fermi surface")
        self.file_label.setText(
            f"<b>{name}</b>: {kind}; {count} bands; E_F = "
            f"{getattr(model, 'fermi', 0.0):.4g} eV"
            + (f"; bands crossing E_F: {', '.join(map(str, crossing))}"
               if crossing else "")
            + (f"; weights: {', '.join(weights[:12])}"
               + (f" and {len(weights) - 12} more" if len(weights) > 12 else "")
               if weights else "")
            + f".<br>It can show: <b>{shows}</b>.")
        self._fs_grid = None
        self._fs_drawn = None
        self._fs_energy_asked = None
        self.weight_combo.blockSignals(True)
        self.weight_combo.clear()
        self.weight_combo.addItem("(one colour)")
        self.weight_combo.addItems(weights)
        self.weight_combo.blockSignals(False)
        self.change_lattice.blockSignals(True)
        self.change_lattice.setChecked(False)
        self.change_lattice.blockSignals(False)
        self._lock_lattice()
        self._plane_cache = None
        if hasattr(self, "cont_through"):
            self.cont_through.clear()
        self._crystal_changed()
        self._fill_compare_list()
        self._output_defaults()
        self.stepper.refresh()

    # -- the cell: read from the calculation, locked unless asked -----------
    def _show_cell(self, params):
        for box, value in zip(self.cell_boxes, (params.a, params.b, params.c,
                                                params.alpha, params.beta,
                                                params.gamma)):
            box.blockSignals(True)
            box.setValue(float(value))
            box.blockSignals(False)

    def _lock_lattice(self):
        """The cell and space group are the calculation's: editable only
        with "Change the lattice parameters" ticked (or with no cell yet)."""
        has_cell = self.model is not None and not isinstance(self.model, B.BandPath)
        free = has_cell and self.change_lattice.isChecked()
        for widget in self.cell_boxes + [self.space_group, self.take_button]:
            widget.setEnabled(free)
        self.change_lattice.setEnabled(has_cell)
        tip = ("" if not has_cell else "Read from the calculation; tick “Change "
               "the lattice parameters” to change it." if not free else "")
        for widget in self.cell_boxes + [self.space_group]:
            widget.setToolTip(tip or widget.toolTip())

    def _cell_typed(self):
        """A cell typed with the box ticked: the calculation scaled to it."""
        if not self.change_lattice.isChecked() or self.model is None \
                or isinstance(self.model, B.BandPath):
            return
        a, b, c, alpha, beta, gamma = (box.value() for box in self.cell_boxes)
        self.set_lattice(LatticeParams(a, b, c, alpha, beta, gamma,
                                       space_group=self.space_group.value()))

    def set_lattice(self, params):
        """Give the calculation the cell ``params``: its bands, held in
        fractional coordinates, keep their place in the zone and stretch
        with it. Same angles: the calculation's own axes are kept."""
        from arpes_viewer.tools.crystal.lattice import primitive_vectors

        old = params_of_lattice(self.model.lattice, self.space_group.value())
        lattice = np.asarray(self.model.lattice, float)
        same_angles = all(abs(x - y) < 1e-6 for x, y in
                          ((old.alpha, params.alpha), (old.beta, params.beta),
                           (old.gamma, params.gamma)))
        if same_angles:
            lattice = lattice / np.linalg.norm(lattice, axis=1)[:, None] * \
                np.array([params.a, params.b, params.c])[:, None]
        else:
            lattice = primitive_vectors(params)
        self.model.reciprocal = B.reciprocal_of(lattice)
        self._fs_grid = self._fs_drawn = self._plane_cache = None
        self._fs_energy_asked = None
        self._show_cell(params_of_lattice(lattice, params.space_group))
        self._crystal_changed()

    def band_number(self, index: int) -> int:
        """Band ``index`` of the calculation as its file counts it."""
        numbers = getattr(self.model, "band_numbers", None)
        return int(numbers[index]) if numbers is not None else index + 1

    def crossing_numbers(self):
        """The bands (as the file counts them) that cross E_F, when that is
        cheap to know (a path or a grid; a Hamiltonian says it when drawn)."""
        model = self.model
        if isinstance(model, B.BandPath):
            rel = model.relative()
        elif isinstance(model, B.BandGrid):
            rel = (model.values - model.fermi).reshape(model.band_count(), -1)
        else:
            return []
        return [self.band_number(i) for i in
                np.flatnonzero((rel.min(1) <= 0) & (rel.max(1) >= 0))]

    def take_space_group(self):
        """The Brillouin-zone tool's crystal: its space group, and -- the box
        ticked -- its cell."""
        params = last_crystal.lattice_params()
        if params is None:
            return
        if params.space_group:
            self.space_group.setValue(int(params.space_group))
        if self.change_lattice.isChecked() and self.model is not None \
                and not isinstance(self.model, B.BandPath):
            self.set_lattice(params)

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
            from arpes_viewer.tools.crystal.highsym import high_symmetry_points
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
    def _path_page(self):
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
        self.weight_combo.currentIndexChanged.connect(lambda *_: self._redraw_view())
        self.highlight_edit = QLineEdit("")
        self.highlight_edit.setPlaceholderText("e.g. 41-44, or EF for those crossing E_F")
        self.highlight_edit.setToolTip(
            "Bands to draw thick and in colour, the others grey (as the "
            "theory side's scripts highlight them): their numbers as the file "
            "counts them, or EF for every band crossing E_F. Empty: all alike.")
        self.highlight_edit.editingFinished.connect(lambda: self.draw_path())
        form.addRow("Highlight bands", self.highlight_edit)
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
        return page

    # ==================================================================
    # Step 2: what the calculation looks like -- by the goal of step 1
    # ==================================================================
    def _step_view(self):
        page = QWidget()
        column = QVBoxLayout(page)
        self.view_pages = QStackedWidget()
        self.view_pages.addWidget(self._path_page())
        self.view_pages.addWidget(self._contour_page())
        self.view_pages.addWidget(self._fs_page())
        top = QHBoxLayout()
        top.addWidget(QLabel("Colour by orbital weight"))
        top.addWidget(self.weight_combo)
        top.addStretch(1)
        column.addLayout(top)
        column.addWidget(self.view_pages, 1)
        title, text = VIEW_STEPS["dispersion"]
        return Step(title, text, page, enter=self._enter_view)

    def _enter_view(self):
        goal = self.goal()
        if goal == "dispersion":
            self._enter_path()
        elif goal == "contours":
            self.draw_contours()
        else:
            if self._fs_drawn is None:
                self.draw_fermi_surface()

    def _redraw_view(self):
        goal = self.goal()
        if goal == "dispersion":
            self.draw_path()
        elif goal == "contours":
            self.draw_contours()

    def highlighted(self, count: int, energies=None):
        """The indices (0-based, of ``count`` bands) the highlight box names."""
        text = self.highlight_edit.text().strip()
        if not text:
            return set()
        if text.upper().replace("_", "") in ("EF", "E F", "FERMI"):
            if energies is None:
                return set()
            e = np.asarray(energies, float)
            return set(np.flatnonzero((np.nanmin(e, 1) <= 0) & (np.nanmax(e, 1) >= 0)))
        wanted = set()
        for part in text.replace(";", ",").split(","):
            part = part.strip()
            try:
                if "-" in part[1:]:
                    a, b = part.split("-", 1)
                    wanted.update(range(int(a), int(b) + 1))
                elif part:
                    wanted.add(int(part))
            except ValueError:
                continue
        return {i for i in range(count) if self.band_number(i) in wanted}

    def band_styles(self, distance, energies, ticks, window, colour):
        """``(band, xs, ys, colour, highlighted)`` of every drawn piece: with
        bands highlighted, those in colour (one each) and the rest grey."""
        from arpes_viewer.ui.panels.fsview import BAND_COLOURS

        chosen = sorted(self.highlighted(np.asarray(energies).shape[0], energies))
        for band, xs, ys in C.band_segments(distance, energies, ticks, window):
            if not chosen:
                yield band, xs, ys, colour, False
            elif band in chosen:
                yield (band, xs, ys,
                       BAND_COLOURS[chosen.index(band) % len(BAND_COLOURS)], True)
            else:
                yield band, xs, ys, "#9a9a9a", False

    # -- constant-energy contours -------------------------------------------
    def _contour_page(self):
        from arpes_viewer.ui.panels.fsview import BandStyleTable

        page = QWidget()
        column = QHBoxLayout(page)
        left = QVBoxLayout()
        form = QFormLayout()
        self.cont_energy = _spin(0.0, -20.0, 20.0, 3, " eV", 0.02)
        self.cont_energy.setToolTip("E - E_F of the contours (0: the Fermi surface).")
        self.cont_energy.valueChanged.connect(lambda *_: self.draw_contours())
        form.addRow("Energy E - E_F", self.cont_energy)
        hkl_row = QHBoxLayout()
        self.cont_h, self.cont_k, self.cont_l = (QSpinBox(), QSpinBox(), QSpinBox())
        for name, box, value in (("h", self.cont_h, 0), ("k", self.cont_k, 0),
                                 ("l", self.cont_l, 1)):
            box.setRange(-9, 9)
            box.setValue(value)
            box.setKeyboardTracking(False)
            box.valueChanged.connect(lambda *_: self._plane_normal_changed())
            hkl_row.addWidget(QLabel(name))
            hkl_row.addWidget(box)
        hkl_row.addStretch(1)
        form.addRow("Plane normal (hkl)", hkl_row)
        self.cont_offset = _spin(0.0, -50.0, 50.0, 4, " Å⁻¹", 0.02)
        self.cont_offset.setToolTip("How far the plane sits from Γ along its own "
                                    "normal: k⊥ relative to the (hkl) direction.")
        self.cont_offset.valueChanged.connect(lambda *_: self.draw_contours())
        form.addRow("k⊥ along the normal", self.cont_offset)
        self.cont_through = QComboBox()
        self.cont_through.setToolTip(
            "The offsets at which the plane holds a high-symmetry point (as in "
            "the Brillouin-zone tool); choosing one moves the plane there.")
        self.cont_through.activated.connect(lambda i: self._through_chosen(i))
        form.addRow("Through point", self.cont_through)
        self.cont_rotation = _spin(0.0, -360.0, 360.0, 2, " °", 5.0)
        self.cont_rotation.setToolTip(
            "Turn the plane's two drawn axes about its normal (counter-"
            "clockwise as drawn) -- to match a measured map's azimuth, say.")
        self.cont_rotation.valueChanged.connect(lambda *_: self.draw_contours())
        form.addRow("In-plane rotation", self.cont_rotation)
        self.cont_window = _spin(0.0, 0.0, 5.0, 3, " Å⁻¹", 0.02)
        self.cont_window.setToolTip(
            "Integrate over this range of k⊥ about the plane -- along the "
            "(hkl) normal, not the absolute k_z -- as a map with a finite k⊥ "
            "resolution does. 0: the plane alone. The integrated intensity is "
            "drawn grey under the contours of the middle plane.")
        self.cont_window.valueChanged.connect(lambda *_: self.draw_contours())
        form.addRow("k⊥ integration window", self.cont_window)
        self.cont_width = _spin(0.02, 0.002, 0.5, 3, " eV", 0.005)
        self.cont_width.setToolTip("Energy half-width given to each band for the "
                                   "integrated intensity.")
        self.cont_width.valueChanged.connect(lambda *_: self.draw_contours())
        form.addRow("Band width (HWHM)", self.cont_width)
        self.cont_reach = _spin(1.0, 0.5, 3.0, 2, " zone(s)", 0.25)
        self.cont_reach.setToolTip("How far out to draw: 1 shows the first zone "
                                   "and a little of its neighbours.")
        self.cont_reach.valueChanged.connect(lambda *_: self.draw_contours())
        form.addRow("Show", self.cont_reach)
        left.addLayout(form)
        zone = QPushButton("Show the plane in the Brillouin zone...")
        zone.setAutoDefault(False)
        zone.setToolTip("A turnable 3-D view of the zone with this plane (and "
                        "its k⊥ window) in it, and the points it holds.")
        zone.clicked.connect(lambda: self.open_zone_window())
        left.addWidget(zone)
        left.addWidget(QLabel("Bands (colour, shown):"))
        self.band_looks = {}
        self.cont_table = BandStyleTable(self.band_looks, opacity=False)
        self.cont_table.changed.connect(lambda: self.draw_contours(recompute=False))
        left.addWidget(self.cont_table)
        self.cont_label = QLabel("")
        self.cont_label.setWordWrap(True)
        self.cont_label.setTextFormat(Qt.RichText)
        left.addWidget(self.cont_label)
        left.addStretch(1)
        column.addLayout(left, 2)
        self.cont_view = pg.PlotWidget()
        self.cont_view.setBackground("w")
        self.cont_view.setAspectLocked(True)
        self.cont_view.setLabel("bottom", "k (Å⁻¹), first in-plane axis")
        self.cont_view.setLabel("left", "k (Å⁻¹), second in-plane axis")
        column.addWidget(self.cont_view, 3)
        self._plane_cache = None
        self._zone_window = None
        return page

    def current_plane(self):
        """The cut of the contour page as a :class:`planecut.Plane`."""
        from arpes_viewer.tools.crystal import planecut

        if self.model is None or isinstance(self.model, B.BandPath) or self.params is None:
            raise ValueError("contours need the bands everywhere: a band grid or "
                             "a Hamiltonian")
        b = np.asarray(self.model.reciprocal, float)
        reach = 0.8 * float(np.linalg.norm(b, axis=1).max()) * self.cont_reach.value()
        return planecut.plane_of(self.params, b, self.hkl(), self.cont_offset.value(),
                                 reach, self.cont_rotation.value())

    def hkl(self):
        return (self.cont_h.value(), self.cont_k.value(), self.cont_l.value())

    def named_points(self):
        from arpes_viewer.tools.crystal import planecut

        key = (id(self.model), self.params)
        if getattr(self, "_points_cache", (None,))[0] != key:
            try:
                points = planecut.named_points(self.params, self.model.reciprocal)
            except Exception:                                  # noqa: BLE001
                points = []
            self._points_cache = (key, points)
        return self._points_cache[1]

    def _plane_normal_changed(self):
        self._fill_through()
        self.draw_contours()

    def _fill_through(self):
        from arpes_viewer.tools.crystal import planecut

        combo = self.cont_through
        combo.clear()
        try:
            normal = self.current_plane().normal
            offsets = planecut.through_points(self.named_points(), normal)
        except ValueError:
            offsets = []
        combo.addItem("Choose a high-symmetry point...", None)
        for offset, labels in offsets[:24]:
            combo.addItem(f"{', '.join(labels)}   ({offset:.4f} Å⁻¹)", offset)
            combo.addItem(f"{', '.join(labels)}   ({-offset:.4f} Å⁻¹)", -offset) \
                if offset > 1e-9 else None

    def _through_chosen(self, index: int):
        offset = self.cont_through.itemData(index)
        if offset is not None:
            self.cont_offset.setValue(float(offset))

    def plane_energies(self, plane, n: int = 161):
        """The bands on ``plane`` (and its k⊥ window): ``(axis, energies)``,
        kept until the plane changes -- a new energy only re-contours."""
        from arpes_viewer.tools.crystal import planecut
        from arpes_viewer.ui.common.jobs import run_blocking

        window = self.cont_window.value()
        key = (id(self.model), tuple(np.round(plane.normal, 9)),
               tuple(np.round(plane.e1, 9)), round(plane.offset, 9),
               round(plane.reach, 9), round(window, 9), n)
        if self._plane_cache is not None and self._plane_cache[0] == key:
            return self._plane_cache[1], self._plane_cache[2]
        model = self.model
        axis, energies = run_blocking(
            self, "The bands on the plane",
            lambda report: planecut.sample(model, plane, n=n, kz_window=window,
                                           progress=report))
        energies = self.alignment.apply(energies)
        self._plane_cache = (key, axis, energies)
        return axis, energies

    def plane_overlays(self, energy: float, *, zone_colour="#888888"):
        """The zone outlines and the contours at ``energy`` on the current
        plane, as figure overlays (colours from the band table, or split by
        the chosen weight), plus ``(axis, image or None)`` -- the integrated
        intensity when there is a k⊥ window."""
        from arpes_viewer.tools.crystal import planecut
        from arpes_viewer.tools.output import figure as F

        plane = self.current_plane()
        axis, energies = self.plane_energies(plane)
        out = []
        for xs, ys in planecut.zone_outlines(self.model.reciprocal, plane):
            out.append(F.Overlay(xs=list(map(float, xs)), ys=list(map(float, ys)),
                                 color=zone_colour, width_pt=0.6, style="dashed",
                                 name="Brillouin zone"))
        weight = (self.weight_combo.currentText()
                  if self.weight_combo.currentIndex() > 0 else None)
        found = planecut.contours(energies, axis, energy)
        bands = sorted({band for band, _xs, _ys in found})
        self.cont_table.set_bands([(band, f"band {self.band_number(band)}")
                                   for band in bands]) \
            if bands != self.cont_table.bands else None
        for band, xs, ys in found:
            style = self.cont_table.style(band, bands.index(band))
            if not style["visible"]:
                continue
            pieces = [(style["colour"], np.arange(len(xs)), "")]
            if weight is not None:
                w = self.model.weights_at(plane.points(xs, ys), weight)[band]
                high = (w >= 0.5).astype(int)
                pieces = []
                for run in np.split(np.arange(len(xs)), np.flatnonzero(np.diff(high)) + 1):
                    if run.size:
                        end = min(run[-1] + 2, len(xs))
                        pieces.append(("#2ca02c" if high[run[0]] else style["colour"],
                                       np.arange(run[0], end),
                                       f" ({weight})" if high[run[0]] else ""))
            for colour, idx, tag in pieces:
                if idx.size >= 2:
                    out.append(F.Overlay(xs=list(map(float, xs[idx])),
                                         ys=list(map(float, ys[idx])), color=colour,
                                         width_pt=0.9,
                                         name=f"band {self.band_number(band)}{tag}"))
        image = (planecut.spectral_map(energies, energy, self.cont_width.value())
                 if len(energies) > 1 else None)
        return out, axis, image

    @one_at_a_time
    def draw_contours(self, recompute: bool = True):
        from arpes_viewer.ui.common.jobs import JobCancelled

        view = self.cont_view
        view.clear()
        if self.model is None or isinstance(self.model, B.BandPath):
            self.cont_label.setText("Contours need the bands everywhere: a band "
                                    "grid or a Hamiltonian.")
            return
        if self.cont_through.count() == 0:
            self._fill_through()
        energy = self.cont_energy.value()
        try:
            overlays, axis, image = self.plane_overlays(energy)
        except JobCancelled:
            self.cont_label.setText("Cancelled.")
            return
        except ValueError as exc:
            self.cont_label.setText(f"<b>Cannot draw them:</b> {exc}")
            return
        if image is not None:
            item = pg.ImageItem(image)
            step = axis[1] - axis[0]
            item.setRect(pg.QtCore.QRectF(axis[0] - step / 2, axis[0] - step / 2,
                                          axis[-1] - axis[0] + step,
                                          axis[-1] - axis[0] + step))
            item.setLookupTable(np.array([[255 - i, 255 - i, 255 - i]
                                          for i in range(0, 200)], dtype=np.uint8))
            view.addItem(item)
        for o in overlays:
            style = Qt.DashLine if o.style == "dashed" else Qt.SolidLine
            view.addItem(pg.PlotDataItem(o.xs, o.ys, pen=pg.mkPen(
                o.color, width=1.6 if o.name.startswith("band") else 0.8, style=style)))
        reach = axis[-1]
        view.setXRange(-reach, reach, padding=0.02)
        view.setYRange(-reach, reach, padding=0.02)
        plane = self.current_plane()
        window = self.cont_window.value()
        on = [p.label for p in self.named_points()
              if abs(float(np.dot(p.k, plane.normal)) - plane.offset) < 1e-4]
        self.cont_label.setText(
            f"({' '.join(map(str, self.hkl()))}) plane at k⊥ = {plane.offset:.4f} Å⁻¹"
            + (f", turned {self.cont_rotation.value():g}°"
               if abs(self.cont_rotation.value()) > 1e-9 else "")
            + (f", integrated over ±{window / 2:.3f} Å⁻¹ of k⊥" if window > 0 else "")
            + (f"; it holds {', '.join(dict.fromkeys(on))}" if on else "")
            + f". E - E_F = {energy:+.3f} eV."
            + ("" if self.cont_table.bands else " No band reaches this energy here.")
            + "<br><i>The zones are dashed grey. Step 3 puts contours over a "
              "measured k-map; step 4 makes figures of them.</i>")
        if self._zone_window is not None and self._zone_window.isVisible():
            self._zone_window.set_plane(plane.normal, plane.offset, window, self.hkl())

    def open_zone_window(self):
        """The zone, in 3-D, with the contour page's plane in it."""
        from arpes_viewer.tools.crystal.bz3d import wigner_seitz_cell
        from arpes_viewer.ui.panels.fsview import ZoneWindow

        try:
            plane = self.current_plane()
        except ValueError as exc:
            self.cont_label.setText(str(exc))
            return None
        if self._zone_window is None or not self._zone_window.isVisible():
            zone = wigner_seitz_cell(np.asarray(self.model.reciprocal, float))
            self._zone_window = ZoneWindow(zone, self.named_points(), self)
        self._zone_window.set_plane(plane.normal, plane.offset,
                                    self.cont_window.value(), self.hkl())
        self._zone_window.show()
        self._zone_window.raise_()
        return self._zone_window

    # -- the 3-D Fermi surface ------------------------------------------------
    def _fs_page(self):
        from arpes_viewer.ui.panels.fsview import VIEWS, BandStyleTable, FermiSurfaceView

        page = QWidget()
        column = QHBoxLayout(page)
        left = QVBoxLayout()
        form = QFormLayout()
        self.fs_energy = _spin(0.0, -20.0, 20.0, 3, " eV", 0.02)
        self.fs_energy.setToolTip("E - E_F of the surface (0: the Fermi surface). "
                                  "Enter draws it.")
        self.fs_energy.editingFinished.connect(lambda: self._fs_energy_entered())
        form.addRow("Energy E - E_F", self.fs_energy)
        self.fs_detail = QComboBox()
        self.fs_detail.addItems(["Coarse (fast)", "Medium", "Fine (slow)"])
        self.fs_detail.setCurrentIndex(1)
        self.fs_detail.setToolTip("How fine the surface's mesh is. A Hamiltonian "
                                  "is also diagonalised on a finer grid -- "
                                  "slower for many orbitals.")
        self.fs_detail.activated.connect(
            lambda *_: self.draw_fermi_surface() if self._fs_drawn is not None else None)
        form.addRow("Detail", self.fs_detail)
        draw = QPushButton("Draw")
        draw.setAutoDefault(False)
        draw.clicked.connect(lambda: self.draw_fermi_surface())
        form.addRow(draw)
        left.addLayout(form)
        left.addWidget(QLabel("Bands (colour -- the same for several to group "
                              "pockets -- opacity, shown):"))
        self.fs_table = BandStyleTable(self.band_looks, opacity=True)
        self.fs_table.changed.connect(lambda: self.fs_view.update())
        left.addWidget(self.fs_table)
        look = QFormLayout()
        self.fs_view_combo = QComboBox()
        self.fs_view_combo.addItems(list(VIEWS))
        self.fs_view_combo.activated.connect(
            lambda *_: self.fs_view.set_view(self.fs_view_combo.currentText()))
        look.addRow("View", self.fs_view_combo)
        self.fs_zone = QCheckBox("The Brillouin zone")
        self.fs_zone.setChecked(True)
        self.fs_labels = QCheckBox("Names of its points")
        self.fs_labels.setChecked(True)
        for box in (self.fs_zone, self.fs_labels):
            box.toggled.connect(lambda *_: self._fs_look())
            look.addRow(box)
        left.addLayout(look)
        self.fs_label = QLabel("")
        self.fs_label.setWordWrap(True)
        left.addWidget(self.fs_label)
        left.addStretch(1)
        column.addLayout(left, 2)
        right = QVBoxLayout()
        self.fs_view = FermiSurfaceView()
        self.fs_view.styles = self.band_looks
        right.addWidget(self.fs_view, 1)
        under = QHBoxLayout()
        under.addWidget(QLabel("<i>drag to turn, wheel to zoom, double-click "
                               "for the first view</i>"))
        under.addStretch(1)
        bigger = QPushButton("Open in its own window")
        bigger.setAutoDefault(False)
        bigger.setToolTip("The same surface, large, in a window of its own "
                          "(with Save picture). Change the energy, detail or "
                          "bands here, then open it again.")
        bigger.clicked.connect(lambda: self.open_fs_window())
        under.addWidget(bigger)
        right.addLayout(under)
        column.addLayout(right, 3)
        return page

    def _fs_energy_entered(self):
        # drawn, or being drawn, at this energy already: nothing to do
        # (editingFinished comes again when the progress window takes the focus)
        if self._fs_energy_asked is None or \
                abs(self._fs_energy_asked - self.fs_energy.value()) > 1e-9:
            self.draw_fermi_surface()

    def open_fs_window(self):
        """The Fermi surface drawn here, in a large window of its own."""
        from arpes_viewer.ui.panels.fsview import FermiSurfaceWindow
        from arpes_viewer.ui.viewers.windows import keep_alive

        if self._fs_drawn is None and self.draw_fermi_surface() is None:
            return None
        name = os.path.basename(self.model_path) or "calculation"
        window = FermiSurfaceWindow(
            self.fs_view, f"3D Fermi surface — {name}, E - E_F = "
                          f"{self._fs_drawn.energy:+.3f} eV")
        keep_alive(window)
        window.show()
        return window

    def _fs_look(self):
        self.fs_view.show_zone = self.fs_zone.isChecked()
        self.fs_view.show_labels = self.fs_labels.isChecked()
        self.fs_view.update()

    def zone_points(self):
        """The zone's named points (one of each) in the calculation's frame:
        ``[(label, k)]``."""
        if self.params is None:
            return []
        try:
            from arpes_viewer.tools.crystal.highsym import high_symmetry_points
            from arpes_viewer.tools.crystal.planecut import to_model_frame

            convert = to_model_frame(self.params, self.model.reciprocal)
            return [(p.label, np.asarray(p.k) @ convert)
                    for p in high_symmetry_points(self.params)]
        except Exception:                                    # noqa: BLE001
            return []

    @one_at_a_time
    def draw_fermi_surface(self):
        """Mesh the surfaces at the energy chosen (with a progress bar: a
        Hamiltonian of many orbitals takes a while) and show them."""
        from arpes_viewer.tools.crystal import fermisurface as FS
        from arpes_viewer.ui.common.jobs import JobCancelled, run_blocking

        if self.model is None or isinstance(self.model, B.BandPath):
            self.fs_view.message = ("The 3D Fermi surface needs the bands "
                                    "everywhere: a band grid or a Hamiltonian.")
            self.fs_view.set_surface(None)
            return None
        n, sample = ((32, 20), (48, 28), (72, 36))[self.fs_detail.currentIndex()]
        energy = self.fs_energy.value()
        self._fs_energy_asked = energy
        model = self.model
        cached = self._fs_grid if (self._fs_grid is not None
                                   and self._fs_grid[0] == sample) else None

        def work(report):
            grid = cached[1] if cached else FS.sampled(model, sample, report)
            fine = FS.fermi_surface(grid, energy, n=n, progress=report)
            coarse = FS.fermi_surface(grid, energy, n=24,
                                      bands=[s.band for s in fine.surfaces])
            return grid, fine, coarse

        try:
            grid, fine, coarse = run_blocking(self, "The 3D Fermi surface", work)
        except JobCancelled:
            self.fs_label.setText("Cancelled.")
            self._fs_energy_asked = None
            return None
        except (ValueError, MemoryError) as exc:
            self.fs_label.setText(f"Could not draw it: {exc}")
            self._fs_energy_asked = None
            return None
        self._fs_grid = (sample, grid)
        self._fs_drawn = fine
        self.fs_table.set_bands([(s.band, f"band {self.band_number(s.band)}")
                                 for s in fine.surfaces])
        self.fs_view.message = "No band reaches this energy."
        self.fs_view.set_surface(fine, coarse, self.zone_points())
        self._fs_look()
        count = sum(len(s.triangles) for s in fine.surfaces)
        self.fs_label.setText(
            f"E - E_F = {energy:+.3f} eV: {len(fine.surfaces)} band(s) "
            f"cross it ({count} triangles)." if fine.surfaces else
            f"No band reaches E - E_F = {energy:+.3f} eV.")
        return fine

    def _enter_path(self):
        self._fill_sources()
        self.draw_path()

    def _maps(self):
        return [(n, d) for n, d in self.datasets if d.kind in ("k_map", "kz_map_k")]

    def _default_kz(self, data):
        from arpes_viewer.tools.spectra.corelevel import photon_energy_of, work_function_of

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
        from arpes_viewer.tools.crystal.lattice import primitive_vectors, reciprocal_vectors

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
        from arpes_viewer.tools.crystal.lattice import primitive_vectors, reciprocal_vectors
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
        for band, xs, ys, pen_colour, wide in self.band_styles(
                distance, energies, ticks, window, colour):
            item = pg.PlotDataItem(xs, ys, pen=pg.mkPen(pen_colour,
                                                        width=2.6 if wide else 1.5))
            item.setZValue(2 if wide else 1)
            view.view.addItem(item)
            self._items.append(item)
        chosen = sorted(self.highlighted(energies.shape[0], energies))
        if chosen:
            self.path_label.setText("Highlighted: " + ", ".join(
                f"band {self.band_number(i)}" for i in chosen) + ".")
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
        if not self.datasets:
            self.compare_label.setText(
                "<b>No measurement was selected</b>, so there is nothing to "
                "compare: press Next for the figures of the calculation alone. "
                "(To compare, select cuts or k-maps in the main window before "
                "opening Calculated bands.)")
            return
        self._compare_changed()

    def _compare_data(self):
        name = self.compare_combo.currentText()
        return next((d for n, d in self.datasets if n == name), None)

    def _compare_changed(self):
        from arpes_viewer.tools.spectra.corelevel import photon_energy_of, work_function_of

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

    def map_contours(self, kx, ky, energy, options=None, align=True):
        """The calculation's contours at measured ``energy`` on the map's
        grid: ``[(band, xs, ys)]`` (the inverse alignment taken first);
        on the k_z of step 3 unless ``options`` say (``{"kz": ...}``)."""
        if isinstance(self.model, B.BandPath):
            raise ValueError("a band structure along a path has no contours")
        calc_energy = (energy - self.alignment.shift) / self.alignment.z
        options = options or self._kz_options()
        return B.contours(self.model, np.linspace(kx.min(), kx.max(), 161),
                          np.linspace(ky.min(), ky.max(), 161), calc_energy,
                          rotation_deg=self.direction_box.value() if align else 0.0,
                          **options)

    def contour_points_k(self, xs, ys, energy, options=None, align=True):
        """3-D k (in the calculation's frame) of contour points drawn on a
        map's (k_x, k_y) at measured ``energy`` -- for their weights."""
        angle = self.direction_box.value() if align else 0.0
        c, s_ = (np.cos(np.deg2rad(-angle)), np.sin(np.deg2rad(-angle)))
        xs, ys = np.asarray(xs, float), np.asarray(ys, float)
        inplane = np.column_stack([c * xs - s_ * ys, s_ * xs + c * ys])
        calc_energy = (energy - self.alignment.shift) / self.alignment.z
        options = options or self._kz_options()
        if "kz" in options:
            return np.column_stack([inplane, np.full(len(xs), options["kz"])])
        return B.kz_on_sphere(inplane, calc_energy, options["hv"],
                              options["inner_potential"], options["work_function"])

    def zone_overlays(self, x_range, y_range, colour="#d62728", style="dashed",
                      width=0.6, align=True):
        """The in-plane Brillouin zones (the calculation's lattice, turned
        by the azimuth) over a (k_x, k_y) range, as overlays."""
        from arpes_viewer.tools.crystal.bz2d import tile_2d, wigner_seitz_cell_2d
        from arpes_viewer.tools.output import figure as F

        if isinstance(self.model, B.BandPath):
            return []
        b = np.asarray(self.model.reciprocal, float)[:2, :2]
        angle = np.deg2rad(self.direction_box.value()) if align else 0.0
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
                         width=0.8, style="solid", threshold=0.5, options=None,
                         by_band=False):
        """The calculated contours at ``energy`` as overlays: in the first
        colour, or -- with a weight chosen in step 2 -- split where the
        chosen orbital's weight crosses ``threshold``, the parts above it in
        the second colour (a surface or orbital character told apart)."""
        from arpes_viewer.tools.output import figure as F

        weight = (self.weight_combo.currentText()
                  if self.weight_combo.currentIndex() > 0 else None)
        out = []
        from arpes_viewer.ui.panels.fsview import BAND_COLOURS

        # the calculation's own frame (no measured map's azimuth) when the
        # kz is given: the contours alone
        align = options is None
        order = {}
        for band, xs, ys in self.map_contours(kx, ky, energy, options, align):
            if by_band and weight is None:
                colour = BAND_COLOURS[order.setdefault(band, len(order)) % len(BAND_COLOURS)]
                out.append(F.Overlay(xs=list(map(float, xs)), ys=list(map(float, ys)),
                                     color=colour, width_pt=width, style=style,
                                     name=f"band {self.band_number(band)}"))
                continue
            xs, ys = np.asarray(xs, float), np.asarray(ys, float)
            if weight is None:
                pieces = [(0, np.arange(xs.size))]
            else:
                w = self.model.weights_at(self.contour_points_k(xs, ys, energy,
                                                                options, align),
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
                                     name=f"band {self.band_number(band)}"
                                          + (f" ({weight})" if which else "")))
        return out

    def panel_contours(self, data=None, energy=None, *, zone=True,
                       colours=("#d62728", "#2ca02c")):
        """The calculated contours alone at ``energy``, over the compared
        map's (k_x, k_y) range, with the zones dashed: the calculation's
        half of a constant-energy comparison."""
        from arpes_viewer.ui.viewers.curves import curve_panel

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

    def panel_contours_alone(self, energy=0.0):
        """The calculated contours at ``energy`` on the contour page's plane
        -- (hkl), k⊥, the band colours chosen there -- the zones dashed; with
        a k⊥ window, the integrated intensity under them. No measurement
        needed: the calculation's own constant-energy map."""
        from arpes_viewer.ui.panels.figure import panel_from_arrays
        from arpes_viewer.ui.viewers.curves import curve_panel

        if isinstance(self.model, B.BandPath):
            raise ValueError("a band structure along a path has no contours: "
                             "open a band grid or a Hamiltonian")
        overlays, axis, image = self.plane_overlays(energy)
        reach = float(axis[-1])
        name = (f"calculation ({' '.join(map(str, self.hkl()))}) plane, k⊥ = "
                f"{self.cont_offset.value():.3g} Å⁻¹, {energy:+.2f} eV")
        labels = ("k (Å⁻¹)", "k (Å⁻¹)") if self.hkl() != (0, 0, 1) else \
            ("k_x (Å⁻¹)", "k_y (Å⁻¹)")
        if image is not None:
            panel = panel_from_arrays(image, axis, axis, *labels, name,
                                      colormap=self.colormap)
        else:
            panel = curve_panel([-reach, reach], [], x_label=labels[0],
                                y_label=labels[1], name=name)
            panel.data.x = np.array([-reach, reach])
            panel.data.y = np.array([-reach, reach])
        panel.x_range = panel.y_range = (-reach, reach)
        panel.overlays += overlays
        panel.title = "E_F" if abs(energy) < 1e-9 else f"E - E_F = {energy:g} eV"
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
        group = QGroupBox("Figure-composer panels")
        g = QVBoxLayout(group)
        self.panel_boxes = {}
        for key, text in (("path", "Measurement along the path, calculation over it"),
                          ("path_calc", "Calculation alone along the path"),
                          ("compare", "The compared cut or map, calculation over it"),
                          ("contours", "Calculated contours alone at the energies "
                                       "below (k_z plane of step 2)"),
                          ("fs", "The 3D Fermi surface, as drawn in step 2")):
            box = QCheckBox(text)
            g.addWidget(box)
            self.panel_boxes[key] = box
        series = QHBoxLayout()
        self.series_box = QCheckBox("Constant-energy series of the compared map at")
        self.series_box.setToolTip(
            "For a compared map: at each energy, the measured map with the "
            "zones dashed (top row) and the calculated contours alone "
            "(bottom row) -- coloured by the weight chosen in step 2.")
        self.series_energies = QLineEdit("0, -0.3, -0.6, -0.9")
        self.series_energies.setToolTip("Energies (eV, E - E_F), separated "
                                        "by commas -- for the series and the "
                                        "contours alone.")
        series.addWidget(self.series_box)
        series.addWidget(self.series_energies, 1)
        g.addLayout(series)
        column.addWidget(group)
        form = QFormLayout()
        self.width_box = _spin(0.03, 0.001, 1.0, 3, " eV", 0.005)
        self.width_box.setToolTip("Lorentzian half-width given to the "
                                  "calculated bands to make a spectrum of them.")
        form.addRow("Band width (HWHM) of the saved spectrum", self.width_box)
        map_range = QHBoxLayout()
        self.map_lo = _spin(-1.0, -20.0, 20.0, 3, " eV", 0.1)
        self.map_hi = _spin(0.2, -20.0, 20.0, 3, " eV", 0.1)
        map_range.addWidget(self.map_lo)
        map_range.addWidget(QLabel("to"))
        map_range.addWidget(self.map_hi)
        map_range.addStretch(1)
        self.map_range_label = QLabel("Energies of the saved map")
        form.addRow(self.map_range_label, map_range)
        self._map_range_widgets = (self.map_lo, self.map_hi, self.map_range_label)
        column.addLayout(form)
        buttons = QHBoxLayout()
        compose = QPushButton("To figure composer...")
        compose.setAutoDefault(False)
        compose.setToolTip("The panels ticked above, in a figure to arrange, "
                           "label and export (PDF, PNG...).")
        compose.clicked.connect(lambda: self.to_figure_composer())
        buttons.addWidget(compose)
        spectrum = QPushButton("Save calculation to main panel")
        spectrum.setAutoDefault(False)
        spectrum.clicked.connect(lambda: self.spectrum_to_list())
        buttons.addWidget(spectrum)
        self.save_button = spectrum
        buttons.addStretch(1)
        column.addLayout(buttons)
        self.output_label = QLabel("")
        self.output_label.setWordWrap(True)
        column.addWidget(self.output_label)
        column.addStretch(1)
        text = ("The figures and the data. Tick the panels to make and press "
                "<b>To figure composer</b> -- arranged, labelled and exported "
                "there like any panel. <b>Save calculation to main panel</b> "
                "makes the calculation a spectrum in the list, so it opens in "
                "the viewers beside the data or is subtracted from it. "
                "<b>Close</b> (bottom right) closes this window -- it asks "
                "first, and keeps your settings for next time.")
        return Step("Output", text, page)

    def _output_defaults(self):
        """Tick the panels that suit the goal and what was selected."""
        if not hasattr(self, "panel_boxes"):
            return
        goal = self.goal()
        everywhere = self.model is not None and not isinstance(self.model, B.BandPath)
        maps = bool(self._maps())
        wanted = {"path": goal == "dispersion" and maps,
                  "path_calc": goal == "dispersion",
                  "compare": bool(self.datasets) and goal != "fs",
                  "contours": goal == "contours" and everywhere,
                  "fs": goal == "fs" and everywhere}
        for key, on in wanted.items():
            self.panel_boxes[key].setChecked(on)
        # what "save to the main panel" makes is what the main panel's
        # viewers can open: a cut, or a constant-energy map -- not a surface
        button = self.save_button
        for widget in self._map_range_widgets:
            widget.setVisible(goal == "contours")
        if goal == "fs":
            button.setText("Save calculation to main panel")
            button.setEnabled(False)
            button.setToolTip("The main panel's viewers open cuts and maps, not "
                              "a 3-D surface: use To figure composer (or, in "
                              "step 2, Open in its own window -> Save picture).")
        elif goal == "contours":
            button.setText("Save calculated map to main panel")
            button.setEnabled(everywhere)
            button.setToolTip("The calculated intensity on the plane of step 2 "
                              "-- its (hkl), k⊥ and k⊥ window -- at every "
                              "energy of the range set here (Lorentzians of the "
                              "width set), as a k-map in the main panel's "
                              "list: open it, slice it, compare it like a "
                              "measured map.")
        else:
            button.setText("Save calculation to main panel")
            button.setEnabled(self.model is not None)
            button.setToolTip("The calculated bands as a spectrum (Lorentzians "
                              "of the width set here), along the path -- or "
                              "along the compared cut, on its own k grid -- as "
                              "a new cut in the main panel's list, with the "
                              "high-symmetry points' names kept.")

    # -- panels -----------------------------------------------------------
    def _overlays(self, distance, energies, ticks, window, colour, weights=None,
                  width=0.9):
        from arpes_viewer.tools.output import figure as F

        out = []
        for band, xs, ys, line_colour, wide in self.band_styles(
                distance, energies, ticks, window, colour):
            out.append(F.Overlay(xs=list(map(float, xs)), ys=list(map(float, ys)),
                                 color=line_colour, width_pt=width * (1.8 if wide else 1),
                                 name=f"band {self.band_number(band)}"))
        # highlighted bands on top
        out.sort(key=lambda o: o.width_pt)
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
        from arpes_viewer.tools.output import figure as F

        panel.x_tick_labels = [(float(d), str(l)) for d, l in ticks]
        for d, _l in ticks[1:-1]:
            panel.lines.append(F.Line(value=float(d), axis="x", color=colour,
                                      width_pt=0.5, style="solid"))

    def panel_path(self, colormap=None):
        """The measurement along the path with the calculation over it."""
        from arpes_viewer.ui.panels.figure import panel_from_arrays

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
        from arpes_viewer.ui.viewers.curves import curve_panel

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
        from arpes_viewer.tools.output import figure as F
        from arpes_viewer.ui.panels.figure import panel_from_arrays

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
        from arpes_viewer.tools.output.figfile import source_record
        from arpes_viewer.ui.panels.figure import _dataset_sources

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
        if self.panel_boxes["contours"].isChecked():
            try:
                out += [self.panel_contours_alone(e)
                        for e in self.series_energies_list() or [0.0]]
            except ValueError as exc:
                problems.append(str(exc))
        if self.panel_boxes["fs"].isChecked():
            try:
                out.append(self.panel_fermi_surface())
            except ValueError as exc:
                problems.append(str(exc))
        return out, problems

    def panel_fermi_surface(self):
        """The 3-D Fermi surface of step 2 as a picture panel."""
        from arpes_viewer.tools.output import figure as nf
        from arpes_viewer.ui.panels.figure import Panel, _qimage_array

        if self._fs_drawn is None:
            if self.draw_fermi_surface() is None:
                raise ValueError("no 3D Fermi surface: draw it in step 2")
        rgba = _qimage_array(self.fs_view.render_image(1400, 1200))
        energy = self._fs_drawn.energy
        panel = Panel(data=nf.picture_data(
            rgba, f"Fermi surface at {energy:+.2f} eV" if abs(energy) > 1e-9
            else "Fermi surface"))
        panel.data.sources = self.provenance([])
        panel.x_labels = panel.y_labels = False
        return panel

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
        from arpes_viewer.ui.panels.figure import FigureWindow, open_figure_windows
        from arpes_viewer.ui.viewers.windows import keep_alive

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

    def build_map(self):
        """The calculation as a constant-energy map on the contour page's
        plane: ``(k1, k2, E)``, Lorentzians of the band width about every
        band, averaged over the k⊥ window."""
        from arpes_viewer.ui.common.jobs import run_blocking

        plane = self.current_plane()
        axis, energies = self.plane_energies(plane)
        lo, hi = sorted((self.map_lo.value(), self.map_hi.value()))
        energy = np.arange(lo, hi + 0.005, 0.01)
        width = self.width_box.value()
        flat = energies.reshape(energies.shape[0], energies.shape[1], -1)
        near = np.flatnonzero((flat.max((0, 2)) >= lo - 10 * width)
                              & (flat.min((0, 2)) <= hi + 10 * width))
        chosen = energies[:, near]

        def work(report):
            cube = np.empty((axis.size, axis.size, energy.size), dtype=np.float32)
            for i, e in enumerate(energy):
                weight = width / ((chosen - e) ** 2 + width ** 2) / np.pi
                cube[:, :, i] = weight.sum(axis=1).mean(axis=0)
                report((i + 1) / energy.size, "The calculated map")
            return cube

        cube = run_blocking(self, "The calculated map", work)
        hkl = " ".join(map(str, self.hkl()))
        parameters = {"calculation": os.path.basename(self.model_path),
                      "plane_hkl": hkl, "k_perp": plane.offset,
                      "in_plane_rotation_deg": self.cont_rotation.value(),
                      "k_perp_window": self.cont_window.value(),
                      "z": self.alignment.z, "shift_eV": self.alignment.shift,
                      "width_eV": width}
        taken = set(self._existing())
        name, n = f"({hkl}) calculated map", 2
        while name in taken:
            name, n = f"({hkl}) calculated map {n}", n + 1
        labels = ({"x": "kx (Å⁻¹)", "k": "ky (Å⁻¹)"} if self.hkl() == (0, 0, 1)
                  else {"x": "k1 (Å⁻¹)", "k": "k2 (Å⁻¹)"})
        labels["z"] = "E - E_F (eV)"
        return MemoryData("k_map", (axis.copy(), axis.copy(), energy), cube, labels,
                          source_label=name, parameters=parameters, prefix="calc",
                          source_path=self.model_path)

    def spectrum_to_list(self):
        from arpes_viewer.ui.common.jobs import JobCancelled

        if self.model is None or self.goal() == "fs":
            return None
        try:
            result = (self.build_map() if self.goal() == "contours"
                      else self.build_spectrum())
        except JobCancelled:
            self.output_label.setText("Cancelled.")
            return None
        except ValueError as exc:
            self.output_label.setText(f"<b>Could not make it:</b> {exc}")
            return None
        self.datasetsCreated.emit([result])
        self.output_label.setText(f"Added {result.source_label} to the list.")
        return result
