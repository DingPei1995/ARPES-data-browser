"""
ui/panels/stackcuts.py
======================
Combining datasets: stacking cuts into a map, or joining the pieces of a
measurement taken in segments.

Reached from **Combine datasets...** under the main list (or its
right-click menu) with several datasets selected. Two ways to combine:

- **Stack** (cuts only): each cut becomes one slice of a new map, along an
  axis the user names -- a role (temperature, gate voltage, photon
  energy...), a name, a unit, a start and a step -- with the value each cut
  gets shown beside its name, so the order and the values are checked
  together. For a stack the launcher first asks whether to keep the
  list's order; if not, the cuts can be dragged into the order wanted.
- **Join** (cuts, maps, kz maps, k-maps): the pieces are put end to end
  along one of their own axes -- a map's deflector range taken in two
  scans, a kz map's photon-energy range split over two runs, a cut's
  energy window in two halves. The dialog shows the range of each piece
  and what the result will be (see :func:`tools.processing.combine.plan`).

The arithmetic is in :mod:`tools.processing.stackcuts` and :mod:`tools.processing.combine`,
which import no Qt.
"""
from __future__ import annotations

import numpy as np
from PyQt5.QtCore import Qt, QTimer
from PyQt5.QtWidgets import (QAbstractItemView, QButtonGroup, QComboBox,
                             QDialog, QDialogButtonBox, QDoubleSpinBox,
                             QFormLayout, QGroupBox, QHBoxLayout, QLabel,
                             QLineEdit, QListWidget, QListWidgetItem,
                             QPushButton, QRadioButton, QVBoxLayout, QWidget)

from loader.registry import AXIS_ROLES
from tools.processing import process as P
from tools.processing import combine as C
from tools.processing import stackcuts
from ui.viewers.widgets import MemoryData

__all__ = ["CombineCutsDialog", "combined_dataset", "joined_dataset"]

#: The roles offered, in the order offered. "angle" is left out: a stack of
#: separate cuts is not a deflector scan the k conversion could trust.
ROLES = ("temperature", "gate_voltage", "photon_energy", "delay",
         "position", "other")


class CombineCutsDialog(QDialog):
    """Ask how to combine, and for the order.

    ``names`` are the datasets' names in the list's order. :meth:`mode` is
    ``"stack"`` or ``"join"``; :meth:`order` gives the chosen order as
    indices into ``names``; :meth:`axis` the new axis of a stack as
    ``(role, name, unit, values)``; :meth:`join_axis` and :meth:`overlap`
    what a join is along and what it does where pieces overlap.

    ``pieces`` -- ``[(axes, labels), ...]`` per dataset, ``axes`` a tuple
    of 1-D arrays and ``labels`` their names -- enables joining; without
    it the dialog only stacks (as it did when it combined cuts only).
    ``stack_problem`` is why stacking is not possible ("" when it is);
    ``can_stack=False`` hides the choice altogether (maps cannot be
    stacked).
    """

    def __init__(self, names, reorder: bool = False, parent=None, *,
                 pieces=None, can_stack: bool = True, stack_problem: str = "",
                 mode: str = None):
        super().__init__(parent)
        self.setWindowTitle("Combine datasets")
        self.names = list(names)
        self.pieces = list(pieces) if pieces else None
        self.resize(500, 560)
        layout = QVBoxLayout(self)

        # -- how ---------------------------------------------------------------
        self.stack_radio = QRadioButton(
            "Stack cuts into a new map, along an axis you name")
        self.stack_radio.setToolTip(
            "Each cut becomes one slice of the map: a temperature series, a "
            "gate-voltage series, one cut per photon energy (a kz map)...")
        self.join_radio = QRadioButton(
            "Join pieces end to end (data taken in segments)")
        self.join_radio.setToolTip(
            "Put the pieces of one measurement back together along one of "
            "their own axes: a map's deflector range scanned in parts, a kz "
            "map's photon energies split over runs, a cut's energy window "
            "in two halves.")
        for radio in (self.stack_radio, self.join_radio):
            radio.setStyleSheet("QRadioButton { padding: 2px; }")
        self.mode_group = QButtonGroup(self)
        self.mode_group.addButton(self.stack_radio)
        self.mode_group.addButton(self.join_radio)
        how = QGroupBox("How to combine")
        how_layout = QVBoxLayout(how)
        how_layout.addWidget(self.stack_radio)
        how_layout.addWidget(self.join_radio)
        self.stack_note = QLabel(stack_problem)
        self.stack_note.setWordWrap(True)
        self.stack_note.setVisible(bool(stack_problem))
        how_layout.addWidget(self.stack_note)
        layout.addWidget(how)
        self.stack_radio.setEnabled(can_stack and not stack_problem)
        self.join_radio.setEnabled(self.pieces is not None)
        how.setVisible(can_stack and self.pieces is not None)
        if mode is None:
            mode = "stack" if self.stack_radio.isEnabled() else "join"
        (self.stack_radio if mode == "stack" else self.join_radio).setChecked(True)

        # -- the new axis of a stack -------------------------------------------
        self.stack_page = QWidget()
        form = QFormLayout(self.stack_page)
        form.setContentsMargins(0, 0, 0, 0)
        self.role = QComboBox()
        for role in ROLES:
            label, unit, _ = AXIS_ROLES[role]
            self.role.addItem(f"{label} ({unit})" if unit else label, role)
        self.name = QLineEdit()
        self.unit = QLineEdit()
        self.start = QDoubleSpinBox()
        self.step = QDoubleSpinBox()
        for box in (self.start, self.step):
            box.setDecimals(4)
            box.setRange(-1e6, 1e6)
        self.start.setValue(0.0)
        self.step.setValue(1.0)
        form.addRow("New axis is", self.role)
        form.addRow("Name", self.name)
        form.addRow("Unit", self.unit)
        form.addRow("Start", self.start)
        form.addRow("Step", self.step)
        layout.addWidget(self.stack_page)

        # -- what a join is along ----------------------------------------------
        self.join_page = QWidget()
        join_form = QFormLayout(self.join_page)
        join_form.setContentsMargins(0, 0, 0, 0)
        self.join_combo = QComboBox()
        labels = self.pieces[0][1] if self.pieces else []
        for index, label in enumerate(labels):
            self.join_combo.addItem(label or f"axis {index + 1}", index)
        if self.pieces:
            # Start on the axis the pieces actually differ along, so a cut
            # taken in two energy windows is not offered "join along angle".
            self.join_combo.setCurrentIndex(self.join_combo.findData(
                C.guess_join_axis([axes for axes, _ in self.pieces])))
        self.join_combo.setToolTip(
            "The axis the pieces were split along: a map's scanned "
            "(deflector) axis, a kz map's photon energy, or a cut's energy "
            "window. The other axes must describe the same thing.")
        self.overlap_combo = QComboBox()
        for key, text in (("average", "Average them"),
                          ("first", "Keep the earlier piece (in the order below)"),
                          ("last", "Keep the later piece (in the order below)")):
            self.overlap_combo.addItem(text, key)
        join_form.addRow("Join along", self.join_combo)
        join_form.addRow("Where pieces overlap", self.overlap_combo)
        self.join_summary = QLabel("")
        self.join_summary.setWordWrap(True)
        join_form.addRow(self.join_summary)
        layout.addWidget(self.join_page)

        # -- the order ---------------------------------------------------------
        layout.addWidget(QLabel(
            "Order (drag to rearrange):" if reorder else
            "Order (as in the list):"))
        self.list = QListWidget()
        for index, name in enumerate(self.names):
            item = QListWidgetItem(name)
            item.setData(Qt.UserRole, index)
            self.list.addItem(item)
        if reorder:
            self.list.setDragDropMode(QAbstractItemView.InternalMove)
            self.list.setDefaultDropAction(Qt.MoveAction)
            self.list.setSelectionMode(QAbstractItemView.SingleSelection)
            # A drop may arrive as a move or as an insert, depending on the
            # Qt version; relabel once it has settled either way.
            model = self.list.model()
            for signal in (model.rowsMoved, model.rowsInserted):
                signal.connect(lambda *_: QTimer.singleShot(
                    0, self._show_values))
        else:
            self.list.setSelectionMode(QAbstractItemView.NoSelection)
        layout.addWidget(self.list, 1)
        if reorder:
            buttons_row = QHBoxLayout()
            for text, delta in (("Move up", -1), ("Move down", 1)):
                button = QPushButton(text)
                button.clicked.connect(lambda _=False, d=delta: self.move(d))
                buttons_row.addWidget(button)
            reverse = QPushButton("Reverse")
            reverse.clicked.connect(self.reverse)
            buttons_row.addWidget(reverse)
            buttons_row.addStretch(1)
            layout.addLayout(buttons_row)

        self.buttons = QDialogButtonBox(QDialogButtonBox.Ok
                                        | QDialogButtonBox.Cancel)
        self.buttons.button(QDialogButtonBox.Ok).setText("Combine")
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        layout.addWidget(self.buttons)

        self.role.currentIndexChanged.connect(self._role_changed)
        for box in (self.start, self.step):
            box.valueChanged.connect(self._show_values)
        self.name.textChanged.connect(self._show_values)
        self.unit.textChanged.connect(self._show_values)
        self.join_combo.currentIndexChanged.connect(self._show_values)
        self.overlap_combo.currentIndexChanged.connect(self._show_values)
        self.mode_group.buttonToggled.connect(self._mode_changed)
        self._role_changed()
        self._mode_changed()

    # -- how -------------------------------------------------------------------
    def mode(self) -> str:
        return "stack" if self.stack_radio.isChecked() else "join"

    def set_mode(self, mode: str):
        (self.stack_radio if mode == "stack" else self.join_radio).setChecked(True)
        self._mode_changed()

    def _mode_changed(self, *_):
        stacking = self.mode() == "stack"
        self.stack_page.setVisible(stacking)
        self.join_page.setVisible(not stacking)
        self._show_values()

    def join_axis(self) -> int:
        return int(self.join_combo.currentData() or 0)

    def overlap(self) -> str:
        return self.overlap_combo.currentData()

    def join_plan(self):
        """:func:`tools.processing.combine.plan` for the current choices, or the
        message saying why they cannot be joined."""
        axes = [self.pieces[i][0] for i in self.order()]
        names = [self.names[i] for i in self.order()]
        try:
            return C.plan(axes, self.join_axis(), names)
        except ValueError as exc:
            return str(exc)

    # -- the order -------------------------------------------------------------
    def order(self) -> list:
        """Indices into ``names``, in the order the cuts will be stacked."""
        return [self.list.item(row).data(Qt.UserRole)
                for row in range(self.list.count())]

    def set_order(self, order):
        """Put the rows in ``order`` (indices into ``names``)."""
        self.list.clear()
        for index in order:
            item = QListWidgetItem(self.names[index])
            item.setData(Qt.UserRole, index)
            self.list.addItem(item)
        self._show_values()

    def move(self, delta: int):
        row = self.list.currentRow()
        target = row + delta
        if row < 0 or not 0 <= target < self.list.count():
            return
        item = self.list.takeItem(row)
        self.list.insertItem(target, item)
        self.list.setCurrentRow(target)
        self._show_values()

    def reverse(self):
        self.set_order(self.order()[::-1])

    # -- the axis --------------------------------------------------------------
    def _role_changed(self, *_):
        label, unit, _ = AXIS_ROLES[self.role.currentData()]
        self.name.setText(label)
        self.unit.setText(unit)
        self._show_values()

    def values(self) -> np.ndarray:
        return stackcuts.stack_axis(self.list.count(), self.start.value(),
                                    self.step.value())

    def axis(self):
        """``(role, name, unit, values)`` of the new axis."""
        return (self.role.currentData(), self.name.text().strip(),
                self.unit.text().strip(), self.values())

    def _show_values(self, *_):
        if not hasattr(self, "buttons"):
            return
        if self.mode() == "join":
            self._show_join()
            return
        unit = self.unit.text().strip()
        ok = self.step.value() != 0
        self.buttons.button(QDialogButtonBox.Ok).setEnabled(ok)
        if not ok:
            return
        for row, value in enumerate(self.values()):
            item = self.list.item(row)
            name = self.names[item.data(Qt.UserRole)]
            item.setText(f"{value:g} {unit}".rstrip() + f"   ←  {name}")

    def _show_join(self):
        axis = self.join_axis()
        for row in range(self.list.count()):
            item = self.list.item(row)
            index = item.data(Qt.UserRole)
            values = np.asarray(self.pieces[index][0][axis], dtype=float)
            item.setText(f"{values.min():.4g} to {values.max():.4g} "
                         f"({values.size} pts)   ←  {self.names[index]}")
        result = self.join_plan()
        ok = isinstance(result, dict)
        self.buttons.button(QDialogButtonBox.Ok).setEnabled(ok)
        if not ok:
            self.join_summary.setText(f"<b>Cannot join:</b> {result}")
            return
        joined = result["axes"][axis]
        shape = " × ".join(str(a.size) for a in result["axes"])
        text = (f"<b>Result:</b> {joined[0]:.4g} to {joined[-1]:.4g} along "
                f"the joined axis, {shape} points"
                + (", values copied as measured." if result["exact"]
                   and not result["resampled"] else "."))
        if result["notes"]:
            text += "<br>" + "<br>".join(f"• {n}" for n in result["notes"])
        self.join_summary.setText(text)


def combined_dataset(cuts, order, role, name, unit, values, label):
    """The stacked cube as a dataset for the list.

    ``cuts`` is ``[(name, data), ...]`` in the list's order, already checked
    with :func:`tools.processing.stackcuts.compatibility`; ``order`` indexes it.
    """
    ordered = [cuts[i] for i in order]
    first = ordered[0][1]
    scan = first.scan
    labels = dict(scan.labels or {})
    angle_label = labels.get("x", "Angle (deg)")
    kind = stackcuts.kind_for_role(role, angle_label)
    new_labels = {"x": stackcuts.axis_label(name, unit), "k": angle_label,
                  "z": labels.get("y", "Energy (eV)")}
    cube = stackcuts.stack([data.cut_frame for _, data in ordered])
    members = [n for n, _ in ordered]
    parameters = {"axis_role": role, "axis_name": name, "axis_unit": unit,
                  "axis_start": float(values[0]),
                  "axis_step": float(values[1] - values[0])
                  if len(values) > 1 else 0.0,
                  "members": members}
    info = P.record_step(dict(scan.info or {}),
                         P.Step("combine_cuts", {"n": len(members),
                                                 "axis": name},
                                source=", ".join(members)))
    info["axis0.role"] = role
    info["combine.axis_values"] = np.asarray(values, dtype=float)
    return MemoryData(kind, (np.asarray(values, dtype=float),
                             np.asarray(scan.x, dtype=float),
                             np.asarray(scan.y, dtype=float)),
                      cube, new_labels, source_label=label,
                      parameters=parameters, prefix="combine",
                      source_path=getattr(first, "path", ""),
                      source_info=info,
                      source_motors=dict(getattr(scan, "fourd_info", {})
                                         or {}))


def joined_dataset(datasets, order, axis, overlap, label):
    """The pieces joined end to end along ``axis``, as a dataset for the
    list.

    ``datasets`` is ``[(name, data), ...]`` in the list's order, all of one
    kind (cuts, or one kind of cube); ``order`` indexes it and decides
    which piece wins an overlap when ``overlap`` is "first" or "last".
    """
    from loader.nxs_file import axis_slots

    ordered = [datasets[i] for i in order]
    first = ordered[0][1]
    kind = first.kind
    slots = axis_slots(kind, "constructor")
    pieces = []
    for _name, data in ordered:
        axes = tuple(np.asarray(getattr(data.scan, slot), dtype=float)
                     for slot in slots)
        values = data.cut_frame if kind == "cut" else data.scan.value
        pieces.append((axes, np.asarray(values)))
    names = [n for n, _ in ordered]
    axes, values, plan = C.join(pieces, axis, overlap, names)

    scan = first.scan
    labels = dict(scan.labels or {})
    n_first = np.asarray(getattr(scan, slots[axis])).size
    info = {}
    for key, value in dict(scan.info or {}).items():
        # Per-slice records (a kz map's fitted Fermi levels, one per photon
        # energy) describe the first piece's slices only.
        if (isinstance(value, np.ndarray) and value.ndim >= 1
                and value.shape[0] == n_first and axis == 0):
            continue
        info[key] = value
    info = P.record_step(info, P.Step(
        "join", {"n": len(names), "axis": labels.get(slots[axis], slots[axis]),
                 "overlap": overlap},
        source=", ".join(names)))
    parameters = {"axis": slots[axis], "overlap": overlap, "members": names,
                  "exact": bool(plan["exact"]),
                  "notes": " ".join(plan["notes"])}
    return MemoryData(kind, tuple(axes), values, labels, source_label=label,
                      parameters=parameters, prefix="join",
                      source_path=getattr(first, "path", ""),
                      source_info=info,
                      source_motors=dict(getattr(scan, "fourd_info", {}) or {}))
