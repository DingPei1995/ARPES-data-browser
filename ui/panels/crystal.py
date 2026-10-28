"""
ui/panels/crystal.py
====================
The crystal structure in real space, in 3-D, step by step.

Opened from the tool panel's Calculation tab (Crystal structure...). Four
steps, the picture beside them all the time:

1. **The structure** -- typed in (cell and atoms), taken from the
   Brillouin-zone tool's last crystal (atoms on its lattice points, to
   which the basis can be added), or read from a CIF file.
2. **The look** -- how many cells, and for each element its colour, size,
   transparency and gloss; the cell edges and bonds.
3. **Cleavage plane** -- (hkl), slid along its normal or put on one of the
   atomic layers listed, with the atoms above it faded or hidden to show
   the termination a cleave would leave.
4. **Export** -- a picture at any size, or straight into the figure
   composer.

Drawn with a plain ``QPainter`` -- shaded spheres, depth-sorted, the same
no-OpenGL approach as the 3-D data views (:mod:`ui.panels.volume`) -- so
the picture can be rendered at any size for a paper. The crystallography
is :mod:`tools.crystal.structure`, which imports no Qt.
"""
from __future__ import annotations

import os

import numpy as np
from PyQt5.QtCore import QPointF, Qt, pyqtSignal
from PyQt5.QtGui import (QBrush, QColor, QImage, QPainter, QPen, QPolygonF,
                         QRadialGradient)
from PyQt5.QtWidgets import (QAbstractItemView, QCheckBox, QColorDialog,
                             QComboBox, QDialog, QDoubleSpinBox, QFileDialog,
                             QFormLayout, QGroupBox, QHBoxLayout, QHeaderView,
                             QInputDialog, QLabel, QListWidget, QPushButton,
                             QRadioButton, QSlider, QSpinBox, QTableWidget,
                             QTableWidgetItem, QVBoxLayout, QWidget)

from tools.crystal import last_crystal
from tools.crystal import structure as S
from tools.processing.volume import Camera
from ui.common.screen import fits_screen
from ui.common.stepper import Step, Stepper

__all__ = ["CrystalDialog", "CrystalView"]


def _blend(colour: QColor, other: QColor, f: float) -> QColor:
    f = float(np.clip(f, 0.0, 1.0))
    return QColor(int(colour.red() * (1 - f) + other.red() * f),
                  int(colour.green() * (1 - f) + other.green() * f),
                  int(colour.blue() * (1 - f) + other.blue() * f),
                  colour.alpha())


class CrystalView(QWidget):
    """The structure drawn in 3-D: drag to turn, wheel to zoom."""

    changed = pyqtSignal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumSize(420, 420)
        self.structure = None
        self.repeats = (1, 1, 1)
        self.styles = {}
        self.atom_scale = 0.45
        self.show_cell = True
        self.show_bonds = True
        self.bond_factor = 1.15
        self.plane = None          # {"hkl", "offset", "colour", "alpha"}
        self.above = "show"        # "show", "fade" or "hide"
        self.mark_plane_atoms = True
        self.background = QColor("white")
        self.camera = Camera(azimuth=-35.0, elevation=28.0)
        self.zoom = 1.0
        self._drag = None

    # -- state -------------------------------------------------------------
    def set_structure(self, structure):
        self.structure = structure
        for element in (structure.elements() if structure else []):
            self.styles.setdefault(element, {
                "colour": S.element_colour(element), "size": 1.0,
                "alpha": 1.0, "gloss": 0.6})
        self.update()

    def view_along(self, direction):
        """Turn the camera to look along ``direction`` (Cartesian)."""
        d = np.asarray(direction, dtype=float)
        d = d / (np.linalg.norm(d) or 1.0)
        self.camera.elevation = float(np.degrees(np.arcsin(np.clip(d[2], -1, 1))))
        self.camera.azimuth = float(np.degrees(np.arctan2(d[1], d[0])))
        self.update()

    # -- the scene ---------------------------------------------------------
    def _items(self):
        """Every primitive, in world coordinates: ("atom", centre, radius,
        style, marked), ("line", p, q, colour, width), ("tri", a, b, c,
        colour)."""
        items = []
        st = self.structure
        if st is None:
            return items, np.zeros((1, 3))
        atoms, xyz = S.cartesian_atoms(st, self.repeats)
        corners = np.array([p for edge in S.box_edges(st, self.repeats)
                            for p in edge])
        normal = offset = None
        if self.plane is not None:
            normal, _d = S.plane_normal(st, self.plane["hkl"])
            offset = self.plane["offset"]
        radii = np.array([S.covalent_radius(a.element) for a in atoms])
        shown = np.ones(len(atoms), dtype=bool)
        faded = np.ones(len(atoms))
        for index, (atom, centre, r) in enumerate(zip(atoms, xyz, radii)):
            style = dict(self.styles.get(atom.element, {}))
            marked = False
            if normal is not None:
                height = float(normal @ centre) - offset
                if height > 0.3:
                    if self.above == "hide":
                        shown[index] = False
                        continue
                    if self.above == "fade":
                        faded[index] = 0.15
                        style["alpha"] = style.get("alpha", 1.0) * 0.15
                marked = self.mark_plane_atoms and abs(height) <= 0.3
            items.append(("atom", centre, r * self.atom_scale * style.get("size", 1.0),
                          style, marked))
        if self.show_bonds and len(xyz) > 1:
            from scipy.spatial import cKDTree
            tree = cKDTree(xyz)
            reach = 2 * radii.max() * self.bond_factor
            for i, j in tree.query_pairs(reach):
                if not (shown[i] and shown[j]):
                    continue
                if (np.linalg.norm(xyz[i] - xyz[j])
                        <= (radii[i] + radii[j]) * self.bond_factor):
                    mid = 0.5 * (xyz[i] + xyz[j])
                    for end, k in ((xyz[i], i), (xyz[j], j)):
                        colour = QColor(self.styles.get(atoms[k].element, {}).get(
                            "colour", "#888888"))
                        colour.setAlphaF(min(faded[i], faded[j]))
                        items.append(("line", end, mid, colour, 4.0))
        if self.show_cell:
            for p, q in S.box_edges(st, self.repeats):
                for t0, t1 in zip(np.linspace(0, 1, 9)[:-1], np.linspace(0, 1, 9)[1:]):
                    items.append(("line", p + t0 * (q - p), p + t1 * (q - p),
                                  "#444444", 1.2))
        if self.plane is not None:
            outline = S.plane_outline(st, self.plane["hkl"], offset, self.repeats)
            if len(outline) >= 3:
                centre = outline.mean(axis=0)
                for k in range(len(outline)):
                    a, b = outline[k], outline[(k + 1) % len(outline)]
                    for tri in _subdivide((centre, a, b), 3):
                        items.append(("tri", *tri, self.plane["colour"],
                                      self.plane["alpha"]))
                    items.append(("line", a, b, self.plane["colour"], 2.0))
        return items, corners

    def paint_into(self, painter: QPainter, width: int, height: int):
        painter.fillRect(0, 0, width, height, self.background)
        items, corners = self._items()
        if not items:
            painter.setPen(QColor("#888"))
            painter.drawText(0, 0, width, height, Qt.AlignCenter,
                             "No structure yet")
            return
        centre = corners.mean(axis=0)
        basis = self.camera.basis()
        projected = (corners - centre) @ basis[:2].T
        extent = float(np.max(np.abs(projected))) or 1.0
        scale = self.zoom * 0.42 * min(width, height) / extent

        def screen(point):
            p = (np.asarray(point) - centre)
            return QPointF(width / 2 + scale * float(p @ basis[0]),
                           height / 2 - scale * float(p @ basis[1]))

        def depth(point):
            return float((np.asarray(point) - centre) @ basis[2])

        keyed = []
        for item in items:
            if item[0] == "atom":
                keyed.append((depth(item[1]), item))
            elif item[0] == "line":
                keyed.append((depth(0.5 * (item[1] + item[2])), item))
            else:
                keyed.append((depth((item[1] + item[2] + item[3]) / 3), item))
        keyed.sort(key=lambda pair: pair[0])
        painter.setRenderHint(QPainter.Antialiasing, True)
        for _d, item in keyed:
            if item[0] == "atom":
                _k, c, r, style, marked = item
                self._paint_atom(painter, screen(c), r * scale, style, marked)
            elif item[0] == "line":
                _k, p, q, colour, w = item
                pen = QPen(QColor(colour), w * max(1.0, scale / 60.0))
                pen.setCapStyle(Qt.RoundCap)
                painter.setPen(pen)
                painter.drawLine(screen(p), screen(q))
            else:
                _k, a, b, c, colour, alpha = item
                fill = QColor(colour)
                fill.setAlphaF(float(alpha))
                painter.setPen(Qt.NoPen)
                painter.setBrush(QBrush(fill))
                painter.drawPolygon(QPolygonF([screen(a), screen(b), screen(c)]))

    @staticmethod
    def _paint_atom(painter, centre: QPointF, radius: float, style, marked):
        if radius <= 0.3:
            return
        base = QColor(style.get("colour", "#999999"))
        alpha = float(np.clip(style.get("alpha", 1.0), 0.0, 1.0))
        gloss = float(np.clip(style.get("gloss", 0.6), 0.0, 1.0))
        focal = QPointF(centre.x() - 0.38 * radius, centre.y() - 0.38 * radius)
        gradient = QRadialGradient(centre, radius, focal)
        highlight = _blend(base, QColor("white"), 0.25 + 0.7 * gloss)
        shadow = base.darker(170 + int(80 * gloss))
        for colour in (highlight, base, shadow):
            colour.setAlphaF(alpha)
        gradient.setColorAt(0.0, highlight)
        gradient.setColorAt(0.12 + 0.35 * (1 - gloss), base)
        gradient.setColorAt(1.0, shadow)
        painter.setBrush(QBrush(gradient))
        if marked:
            painter.setPen(QPen(QColor("#e0302a"), max(1.5, radius * 0.12)))
        else:
            edge = QColor(0, 0, 0, int(60 * alpha))
            painter.setPen(QPen(edge, 0.8))
        painter.drawEllipse(centre, radius, radius)

    def render_image(self, width: int, height: int, background=None) -> QImage:
        image = QImage(int(width), int(height), QImage.Format_ARGB32)
        image.fill(background or self.background)
        painter = QPainter(image)
        try:
            self.paint_into(painter, int(width), int(height))
        finally:
            painter.end()
        return image

    # -- Qt ----------------------------------------------------------------
    def paintEvent(self, _event):
        painter = QPainter(self)
        try:
            self.paint_into(painter, self.width(), self.height())
        finally:
            painter.end()

    def mousePressEvent(self, event):
        self._drag = (event.pos(), self.camera.azimuth, self.camera.elevation)

    def mouseMoveEvent(self, event):
        if self._drag is None:
            return
        start, az, el = self._drag
        delta = event.pos() - start
        self.camera.azimuth = az - 0.5 * delta.x()
        self.camera.elevation = float(np.clip(el + 0.5 * delta.y(), -89.0, 89.0))
        self.update()

    def mouseReleaseEvent(self, _event):
        self._drag = None

    def wheelEvent(self, event):
        self.zoom *= 1.1 ** (event.angleDelta().y() / 120.0)
        self.zoom = float(np.clip(self.zoom, 0.2, 8.0))
        self.update()


def _subdivide(triangle, levels: int):
    """A triangle split into 4**levels smaller ones (for depth sorting)."""
    out = [tuple(np.asarray(p, dtype=float) for p in triangle)]
    for _ in range(levels):
        nxt = []
        for a, b, c in out:
            ab, bc, ca = (a + b) / 2, (b + c) / 2, (c + a) / 2
            nxt += [(a, ab, ca), (ab, b, bc), (ca, bc, c), (ab, bc, ca)]
        out = nxt
    return out


class _ColourButton(QPushButton):
    colourChanged = pyqtSignal(str)

    def __init__(self, colour: str):
        super().__init__("")
        self.colour = colour
        self.setFixedWidth(46)
        self._show()
        self.clicked.connect(self._pick)

    def _show(self):
        self.setStyleSheet(f"background-color: {self.colour}; border: 1px "
                           f"solid #888;")

    def _pick(self):
        chosen = QColorDialog.getColor(QColor(self.colour), self, "Colour")
        if chosen.isValid():
            self.set_colour(chosen.name())

    def set_colour(self, colour: str):
        self.colour = colour
        self._show()
        self.colourChanged.emit(colour)


@fits_screen
class CrystalDialog(QDialog):
    """The crystal structure viewer."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Crystal structure")
        self.resize(1200, 760)
        self.structure = None
        self.view = CrystalView(self)
        layout = QHBoxLayout(self)
        self.stepper = Stepper([self._step_structure(), self._step_look(),
                                self._step_plane(), self._step_export()], self,
                               finish_text="Close", show_index=False)
        self.stepper.setMinimumWidth(470)
        self.stepper.finished.connect(self.accept)
        layout.addWidget(self.stepper, 2)
        right = QVBoxLayout()
        right.addWidget(self.view, 1)
        views = QHBoxLayout()
        views.addWidget(QLabel("<i>drag to turn, wheel to zoom</i>  View along:"))
        for text, which in (("a", 0), ("b", 1), ("c", 2), ("plane normal", None)):
            button = QPushButton(text)
            button.clicked.connect(lambda _=False, w=which: self.look_along(w))
            views.addWidget(button)
        views.addStretch(1)
        right.addLayout(views)
        layout.addLayout(right, 3)
        self.manual_changed()

    # ==================================================================
    # Step 1: the structure
    # ==================================================================
    def _step_structure(self):
        page = QWidget()
        column = QVBoxLayout(page)
        source = QHBoxLayout()
        self.src_manual = QRadioButton("Type it in")
        self.src_bz = QRadioButton("From the Brillouin-zone tool")
        self.src_cif = QRadioButton("From a CIF file")
        self.src_manual.setChecked(True)
        for button in (self.src_manual, self.src_bz, self.src_cif):
            source.addWidget(button)
            button.toggled.connect(self._source_changed)
        column.addLayout(source)

        cell = QGroupBox("Cell")
        form = QFormLayout(cell)
        self.cell_boxes = {}
        row = QHBoxLayout()
        for name, value, suffix in (("a", 4.0, " Å"), ("b", 4.0, " Å"),
                                    ("c", 4.0, " Å")):
            box = QDoubleSpinBox()
            box.setDecimals(4)
            box.setRange(0.1, 200.0)
            box.setValue(value)
            box.setSuffix(suffix)
            box.setKeyboardTracking(False)
            box.valueChanged.connect(self.manual_changed)
            row.addWidget(QLabel(name))
            row.addWidget(box)
            self.cell_boxes[name] = box
        form.addRow(row)
        row = QHBoxLayout()
        for name in ("alpha", "beta", "gamma"):
            box = QDoubleSpinBox()
            box.setDecimals(3)
            box.setRange(1.0, 179.0)
            box.setValue(90.0)
            box.setSuffix(" °")
            box.setKeyboardTracking(False)
            box.valueChanged.connect(self.manual_changed)
            row.addWidget(QLabel(name))
            row.addWidget(box)
            self.cell_boxes[name] = box
        form.addRow(row)
        self.space_group = QSpinBox()
        self.space_group.setRange(1, 230)
        self.space_group.setToolTip("Sets the centering, whose translations "
                                    "are applied to the atoms below (P, A, B, "
                                    "C, I, F, R). The full symmetry "
                                    "operations come with a CIF file.")
        self.space_group.valueChanged.connect(self.manual_changed)
        form.addRow("Space group", self.space_group)
        self.centering_label = QLabel("")
        form.addRow("", self.centering_label)
        column.addWidget(cell)

        atoms = QGroupBox("Atoms (fractional coordinates; the centering "
                          "copies are added)")
        atoms_column = QVBoxLayout(atoms)
        self.atom_table = QTableWidget(0, 4)
        self.atom_table.setHorizontalHeaderLabels(["Element", "x", "y", "z"])
        self.atom_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.atom_table.itemChanged.connect(self.manual_changed)
        atoms_column.addWidget(self.atom_table)
        buttons = QHBoxLayout()
        add = QPushButton("Add atom")
        add.clicked.connect(lambda: self.add_atom("X", 0, 0, 0))
        remove = QPushButton("Remove atom")
        remove.clicked.connect(self.remove_atom)
        buttons.addWidget(add)
        buttons.addWidget(remove)
        atoms_column.addLayout(buttons)
        column.addWidget(atoms, 1)

        self.bz_button = QPushButton("Take the cell from the Brillouin-zone tool")
        self.bz_button.setToolTip(last_crystal.describe())
        self.bz_button.clicked.connect(self.take_from_bz)
        column.addWidget(self.bz_button)
        self.cif_button = QPushButton("Read a CIF file...")
        self.cif_button.clicked.connect(lambda: self.read_cif())
        column.addWidget(self.cif_button)
        self.structure_label = QLabel("")
        self.structure_label.setWordWrap(True)
        column.addWidget(self.structure_label)
        self.add_atom("X", 0, 0, 0)
        self._source_changed()
        text = ("Where the crystal comes from. <b>Type it in</b>: the cell, "
                "the space group (for its centering) and the atoms of the "
                "basis. <b>From the Brillouin-zone tool</b>: the crystal you "
                "last drew a zone for, with an atom on each lattice point -- "
                "edit them into the real basis. <b>From a CIF file</b> "
                "(ICSD, COD, Materials Project...): the full cell, its "
                "symmetry applied.")
        return Step("The structure", text, page,
                    ready=lambda: "" if self.structure is not None and
                    self.structure.atoms else "No atoms yet.")

    def _source_changed(self, *_):
        manual = self.src_manual.isChecked() or self.src_bz.isChecked()
        self.bz_button.setVisible(self.src_bz.isChecked())
        self.cif_button.setVisible(self.src_cif.isChecked())
        for box in self.cell_boxes.values():
            box.setEnabled(manual)
        self.space_group.setEnabled(manual)
        self.atom_table.setEnabled(manual)

    def add_atom(self, element, x, y, z):
        self.atom_table.blockSignals(True)
        row = self.atom_table.rowCount()
        self.atom_table.insertRow(row)
        for column, value in enumerate((element, x, y, z)):
            self.atom_table.setItem(row, column, QTableWidgetItem(str(value)))
        self.atom_table.blockSignals(False)
        self.manual_changed()

    def remove_atom(self):
        row = self.atom_table.currentRow()
        if row >= 0:
            self.atom_table.removeRow(row)
            self.manual_changed()

    def _table_atoms(self):
        atoms = []
        for row in range(self.atom_table.rowCount()):
            try:
                element = S.element_of(self.atom_table.item(row, 0).text())
                frac = tuple(float(self.atom_table.item(row, c).text())
                             for c in (1, 2, 3))
            except (AttributeError, ValueError):
                continue
            atoms.append(S.Atom(f"{element}{row + 1}", element, frac))
        return atoms

    def manual_changed(self, *_):
        if self.src_cif.isChecked() and self.structure is not None:
            return
        from tools.crystal.spacegroups import spacegroup_info
        info = spacegroup_info(self.space_group.value())
        centering = info.centering
        self.centering_label.setText(f"{info.symbol}, {info.crystal_system}: "
                                     f"centering {centering}")
        operations = [(np.eye(3), np.asarray(t, dtype=float))
                      for t in S.CENTERING_TRANSLATIONS.get(centering, [(0, 0, 0)])]
        b = self.cell_boxes
        self.set_structure(S.Structure(
            b["a"].value(), b["b"].value(), b["c"].value(), b["alpha"].value(),
            b["beta"].value(), b["gamma"].value(),
            atoms=S.expand(self._table_atoms(), operations),
            name="", space_group=str(self.space_group.value())))

    def take_from_bz(self):
        values = last_crystal.load()
        params = last_crystal.lattice_params(values)
        if params is None:
            two = last_crystal.vectors_2d(values)
            if two is None:
                self.structure_label.setText(last_crystal.describe())
                return False
            a, b_, gamma, _az = two
            params_values = (a, b_, 10.0, 90.0, 90.0, gamma, 1)
        else:
            params_values = (params.a, params.b, params.c, params.alpha,
                             params.beta, params.gamma, params.space_group or 1)
        for name, value in zip(("a", "b", "c", "alpha", "beta", "gamma"),
                               params_values[:6]):
            self.cell_boxes[name].blockSignals(True)
            self.cell_boxes[name].setValue(value)
            self.cell_boxes[name].blockSignals(False)
        self.space_group.blockSignals(True)
        self.space_group.setValue(int(params_values[6]))
        self.space_group.blockSignals(False)
        self.atom_table.setRowCount(0)
        self.add_atom("X", 0, 0, 0)
        self.structure_label.setText(
            "Taken: " + last_crystal.describe(values) + ". One atom per lattice "
            "point -- set its element and add the rest of the basis.")
        return True

    def read_cif(self, path=None):
        if path is None:
            path, _ = QFileDialog.getOpenFileName(self, "CIF file", "",
                                                  "CIF (*.cif);;All files (*)")
            if not path:
                return None
        try:
            structure = S.read_cif(path)
        except (OSError, ValueError) as exc:
            self.structure_label.setText(f"<b>Could not read it:</b> {exc}")
            return None
        self.src_cif.setChecked(True)
        self.set_structure(structure)
        self.structure_label.setText(
            f"{structure.name or os.path.basename(path)} "
            f"({structure.formula()}), {structure.space_group}: "
            f"a = {structure.a:.4g}, b = {structure.b:.4g}, c = {structure.c:.4g} Å, "
            f"{structure.alpha:.4g}/{structure.beta:.4g}/{structure.gamma:.4g}°; "
            f"{len(structure.atoms)} atoms in the cell.")
        return structure

    def set_structure(self, structure):
        self.structure = structure
        self.view.set_structure(structure)
        if hasattr(self, "style_table"):
            self._fill_styles()
            self._plane_changed()
        if hasattr(self, "stepper"):
            self.stepper.refresh()

    # ==================================================================
    # Step 2: the look
    # ==================================================================
    def _step_look(self):
        page = QWidget()
        column = QVBoxLayout(page)
        form = QFormLayout()
        repeat_row = QHBoxLayout()
        self.repeat_boxes = []
        for name in ("a", "b", "c"):
            box = QSpinBox()
            box.setRange(1, 8)
            box.setValue(1)
            box.valueChanged.connect(self._look_changed)
            repeat_row.addWidget(QLabel(f"× {name}"))
            repeat_row.addWidget(box)
            self.repeat_boxes.append(box)
        form.addRow("Cells", repeat_row)
        self.scale_box = QDoubleSpinBox()
        self.scale_box.setRange(0.05, 2.0)
        self.scale_box.setSingleStep(0.05)
        self.scale_box.setValue(0.45)
        self.scale_box.setToolTip("All atoms' radii, as a fraction of the "
                                  "covalent radius (1 = touching neighbours).")
        self.scale_box.valueChanged.connect(self._look_changed)
        form.addRow("Atom size", self.scale_box)
        self.cell_box = QCheckBox("Cell edges")
        self.cell_box.setChecked(True)
        self.bond_box = QCheckBox("Bonds (atoms closer than 1.15 x the sum of "
                                  "their covalent radii)")
        self.bond_box.setChecked(True)
        for box in (self.cell_box, self.bond_box):
            box.toggled.connect(self._look_changed)
            form.addRow(box)
        self.background_combo = QComboBox()
        self.background_combo.addItems(["White", "Black", "Transparent"])
        self.background_combo.currentIndexChanged.connect(self._look_changed)
        form.addRow("Background", self.background_combo)
        column.addLayout(form)
        self.style_table = QTableWidget(0, 5)
        self.style_table.setHorizontalHeaderLabels(
            ["Element", "Colour", "Size ×", "Transparency %", "Gloss %"])
        self.style_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.style_table.setSelectionMode(QAbstractItemView.NoSelection)
        column.addWidget(self.style_table, 1)
        text = ("How it looks. Repeat the cell along a, b and c, set the "
                "atoms' overall size, and for each element its "
                "<b>colour</b>, <b>size</b>, <b>transparency</b> and "
                "<b>gloss</b> (0 % a matt ball, 100 % a shiny one). Drag the "
                "picture to turn it, scroll to zoom; the buttons under it "
                "look straight along a, b, c or a plane's normal.")
        return Step("The look", text, page, enter=self._fill_styles)

    def _fill_styles(self):
        table = self.style_table
        elements = self.structure.elements() if self.structure else []
        table.setRowCount(len(elements))
        for row, element in enumerate(elements):
            style = self.view.styles.setdefault(element, {
                "colour": S.element_colour(element), "size": 1.0,
                "alpha": 1.0, "gloss": 0.6})
            table.setItem(row, 0, QTableWidgetItem(element))
            colour = _ColourButton(style["colour"])
            colour.colourChanged.connect(
                lambda c, e=element: self.set_style(e, colour=c))
            table.setCellWidget(row, 1, colour)
            for column, key, lo, hi, factor in ((2, "size", 0.1, 5.0, 1.0),
                                                (3, "alpha", 0, 100, None),
                                                (4, "gloss", 0, 100, 100.0)):
                box = QDoubleSpinBox()
                box.setRange(lo, hi)
                if key == "alpha":
                    box.setValue(100.0 * (1.0 - style["alpha"]))
                    box.valueChanged.connect(
                        lambda v, e=element: self.set_style(e, alpha=1 - v / 100.0))
                elif key == "gloss":
                    box.setValue(100.0 * style["gloss"])
                    box.valueChanged.connect(
                        lambda v, e=element: self.set_style(e, gloss=v / 100.0))
                else:
                    box.setSingleStep(0.1)
                    box.setValue(style["size"])
                    box.valueChanged.connect(
                        lambda v, e=element: self.set_style(e, size=v))
                table.setCellWidget(row, column, box)

    def set_style(self, element, **values):
        self.view.styles.setdefault(element, {}).update(values)
        self.view.update()

    def _look_changed(self, *_):
        self.view.repeats = tuple(box.value() for box in self.repeat_boxes)
        self.view.atom_scale = self.scale_box.value()
        self.view.show_cell = self.cell_box.isChecked()
        self.view.show_bonds = self.bond_box.isChecked()
        self.view.background = QColor(
            ("white", "black", "transparent")[self.background_combo.currentIndex()])
        if self.background_combo.currentIndex() == 2:
            self.view.background = QColor(0, 0, 0, 0)
        if hasattr(self, "plane_box"):
            self._plane_changed()
        self.view.update()

    # ==================================================================
    # Step 3: the cleavage plane
    # ==================================================================
    def _step_plane(self):
        page = QWidget()
        column = QVBoxLayout(page)
        form = QFormLayout()
        self.plane_box = QCheckBox("Draw a plane")
        self.plane_box.toggled.connect(self._plane_changed)
        form.addRow(self.plane_box)
        hkl = QHBoxLayout()
        self.hkl_boxes = []
        for name, value in (("h", 0), ("k", 0), ("l", 1)):
            box = QSpinBox()
            box.setRange(-9, 9)
            box.setValue(value)
            box.valueChanged.connect(self._plane_changed)
            hkl.addWidget(QLabel(name))
            hkl.addWidget(box)
            self.hkl_boxes.append(box)
        form.addRow("Miller indices", hkl)
        self.offset_slider = QSlider(Qt.Horizontal)
        self.offset_slider.setRange(0, 1000)
        self.offset_slider.setValue(500)
        self.offset_slider.valueChanged.connect(self._offset_moved)
        self.offset_label = QLabel("")
        offset_row = QHBoxLayout()
        offset_row.addWidget(self.offset_slider, 1)
        offset_row.addWidget(self.offset_label)
        form.addRow("Position", offset_row)
        self.above_combo = QComboBox()
        self.above_combo.addItems(["Show the atoms above it", "Fade them",
                                   "Hide them (the cleaved surface)"])
        self.above_combo.currentIndexChanged.connect(self._plane_changed)
        form.addRow("Above the plane", self.above_combo)
        look = QHBoxLayout()
        self.plane_colour = _ColourButton("#3a7bd5")
        self.plane_colour.colourChanged.connect(self._plane_changed)
        look.addWidget(self.plane_colour)
        self.plane_alpha = QSpinBox()
        self.plane_alpha.setRange(0, 100)
        self.plane_alpha.setValue(65)
        self.plane_alpha.setSuffix(" % transparent")
        self.plane_alpha.valueChanged.connect(self._plane_changed)
        look.addWidget(self.plane_alpha)
        form.addRow("Plane colour", look)
        column.addLayout(form)
        column.addWidget(QLabel("Atomic layers parallel to the plane "
                                "(double-click to put the plane on one):"))
        self.layer_list = QListWidget()
        self.layer_list.itemDoubleClicked.connect(
            lambda item: self.set_plane_height(item.data(Qt.UserRole)))
        column.addWidget(self.layer_list, 1)
        self.plane_label = QLabel("")
        self.plane_label.setWordWrap(True)
        column.addWidget(self.plane_label)
        text = ("A <b>cleavage plane</b> (hkl): tick <b>Draw</b> and give its "
                "Miller indices (in the cell above). Slide it along its "
                "normal, or double-click one of the <b>atomic layers</b> to "
                "put it exactly there -- each line says which atoms that "
                "layer holds, so the possible terminations of a cleave read "
                "off the list. <b>Hide the atoms above</b> to see the surface "
                "the cleave would leave; atoms in the plane are ringed red.")
        return Step("Cleavage plane", text, page, enter=self._plane_changed)

    def hkl(self):
        return tuple(box.value() for box in self.hkl_boxes)

    def _plane_range(self):
        return S.offset_range(self.structure, self.hkl(), self.view.repeats)

    def plane_height(self):
        lo, hi = self._plane_range()
        return lo + (hi - lo) * self.offset_slider.value() / 1000.0

    def set_plane_height(self, height):
        lo, hi = self._plane_range()
        self.offset_slider.setValue(int(round(1000 * (height - lo) / ((hi - lo) or 1))))
        self._plane_changed()

    def _offset_moved(self, *_):
        self._plane_changed(rebuild_layers=False)

    def _plane_changed(self, *_, rebuild_layers=True):
        if self.structure is None or not any(self.hkl()):
            self.view.plane = None
            self.view.update()
            return
        on = self.plane_box.isChecked()
        try:
            _normal, d = S.plane_normal(self.structure, self.hkl())
        except ValueError as exc:
            self.plane_label.setText(str(exc))
            return
        height = self.plane_height()
        self.offset_label.setText(f"{height:.3f} Å")
        self.view.plane = ({"hkl": self.hkl(), "offset": height,
                            "colour": self.plane_colour.colour,
                            "alpha": 1.0 - self.plane_alpha.value() / 100.0}
                           if on else None)
        self.view.above = ("show", "fade", "hide")[self.above_combo.currentIndex()]
        if rebuild_layers:
            self.layer_list.clear()
            for h, counts in S.atomic_layers(self.structure, self.hkl(),
                                             self.view.repeats):
                text = ", ".join(f"{e}{'' if n == 1 else f' x{n}'}"
                                 for e, n in counts.items())
                self.layer_list.addItem(f"{h:8.3f} Å   {text}")
                self.layer_list.item(self.layer_list.count() - 1).setData(
                    Qt.UserRole, h)
        self.plane_label.setText(f"d({' '.join(map(str, self.hkl()))}) = {d:.4f} Å")
        self.view.update()

    def look_along(self, which):
        if self.structure is None:
            return
        if which is None:
            if any(self.hkl()):
                normal, _d = S.plane_normal(self.structure, self.hkl())
                self.view.view_along(normal)
            return
        self.view.view_along(self.structure.vectors()[which])

    # ==================================================================
    # Step 4: export
    # ==================================================================
    def _step_export(self):
        page = QWidget()
        form = QFormLayout(page)
        size = QHBoxLayout()
        self.width_box = QSpinBox()
        self.height_box = QSpinBox()
        for box, value in ((self.width_box, 1600), (self.height_box, 1600)):
            box.setRange(100, 12000)
            box.setValue(value)
            box.setSingleStep(100)
            size.addWidget(box)
        form.addRow("Size (px)", size)
        save = QPushButton("Save picture...")
        save.clicked.connect(lambda: self.save_picture())
        form.addRow(save)
        compose = QPushButton("To figure composer...")
        compose.setToolTip("Add the picture as a panel of a figure composer "
                           "that is open, or of a new one.")
        compose.clicked.connect(lambda: self.to_figure_composer())
        form.addRow(compose)
        self.export_label = QLabel("")
        self.export_label.setWordWrap(True)
        form.addRow(self.export_label)
        text = ("The picture as it is shown, at the size you choose -- it is "
                "drawn again at that size, not enlarged. <b>Save picture</b> "
                "writes a PNG (transparent if the background is). <b>To "
                "figure composer</b> puts it on a figure page beside your "
                "data, where it is labelled and arranged like any panel.")
        return Step("Export", text, page)

    def picture(self):
        return self.view.render_image(self.width_box.value(),
                                      self.height_box.value())

    def save_picture(self, path=None):
        if path is None:
            name = (self.structure.name if self.structure else "") or "crystal"
            path, _ = QFileDialog.getSaveFileName(
                self, "Save picture", f"{name}.png", "PNG (*.png)")
            if not path:
                return None
        if not path.lower().endswith(".png"):
            path += ".png"
        self.picture().save(path, "PNG")
        self.export_label.setText(f"Wrote {path}.")
        return path

    def to_figure_composer(self, target=None, ask=True):
        from tools.output import figure as nf
        from ui.panels.figure import (FigureWindow, Panel, _qimage_array,
                                      open_figure_windows)
        from ui.viewers.windows import keep_alive

        rgba = _qimage_array(self.picture())
        name = (self.structure.name if self.structure else "") or "crystal"
        panel = Panel(data=nf.picture_data(rgba, name))
        panel.x_labels = panel.y_labels = False
        if target is None and ask:
            windows = open_figure_windows()
            if windows:
                names = ["A new figure"] + [w.windowTitle() for w in windows]
                choice, ok = QInputDialog.getItem(
                    self, "To figure composer", "Add the picture to:", names,
                    0, False)
                if not ok:
                    return None
                index = names.index(choice)
                target = windows[index - 1] if index > 0 else None
        if target is None:
            target = FigureWindow([panel], f"Figure — {name}", parent=None)
            keep_alive(target)
            target.show()
        else:
            target.add_panels([panel])
            target.raise_()
        self.export_label.setText(f"Added to “{target.windowTitle()}”.")
        return target
