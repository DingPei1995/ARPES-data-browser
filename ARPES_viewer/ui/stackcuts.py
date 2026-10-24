"""
ui/stackcuts.py
===============
Combining cuts into a map: the order, and what the new axis is.

Reached from the main list's right-click menu with several cuts selected.
The launcher first asks whether to combine them in the list's order; if
not, this dialog lets the cuts be dragged into the order wanted. Either
way the dialog asks what the new axis is -- a role (temperature, gate
voltage, photon energy...), a name, a unit, a start and a step -- and shows
the value each cut gets beside its name, so the order and the values are
checked together.

The stacking itself is :mod:`tools.stackcuts`, which imports no Qt.
"""
from __future__ import annotations

import numpy as np
from PyQt5.QtCore import Qt, QTimer
from PyQt5.QtWidgets import (QAbstractItemView, QComboBox, QDialog,
                             QDialogButtonBox, QDoubleSpinBox, QFormLayout,
                             QHBoxLayout, QLabel, QLineEdit, QListWidget,
                             QListWidgetItem, QPushButton, QVBoxLayout)

from loader.registry import AXIS_ROLES
from tools import process as P
from tools import stackcuts
from ui.widgets import MemoryData

__all__ = ["CombineCutsDialog", "combined_dataset"]

#: The roles offered, in the order offered. "angle" is left out: a stack of
#: separate cuts is not a deflector scan the k conversion could trust.
ROLES = ("temperature", "gate_voltage", "photon_energy", "delay",
         "position", "other")


class CombineCutsDialog(QDialog):
    """Ask for the new axis and, when ``reorder`` is set, the order.

    ``names`` are the cuts' names in the list's order. :meth:`order` gives
    the chosen order as indices into ``names``; :meth:`axis` the new axis
    as ``(role, name, unit, values)``.
    """

    def __init__(self, names, reorder: bool = False, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Combine cuts into a map")
        self.names = list(names)
        self.resize(460, 480)
        layout = QVBoxLayout(self)

        # -- the new axis ------------------------------------------------------
        form = QFormLayout()
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
        form.addRow("Axis is", self.role)
        form.addRow("Name", self.name)
        form.addRow("Unit", self.unit)
        form.addRow("Start", self.start)
        form.addRow("Step", self.step)
        layout.addLayout(form)

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
        self._role_changed()

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
        unit = self.unit.text().strip()
        ok = self.step.value() != 0
        self.buttons.button(QDialogButtonBox.Ok).setEnabled(ok)
        if not ok:
            return
        for row, value in enumerate(self.values()):
            item = self.list.item(row)
            name = self.names[item.data(Qt.UserRole)]
            item.setText(f"{value:g} {unit}".rstrip() + f"   ←  {name}")


def combined_dataset(cuts, order, role, name, unit, values, label):
    """The stacked cube as a dataset for the list.

    ``cuts`` is ``[(name, data), ...]`` in the list's order, already checked
    with :func:`tools.stackcuts.compatibility`; ``order`` indexes it.
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
