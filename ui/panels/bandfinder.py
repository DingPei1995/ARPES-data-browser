"""
ui/panels/bandfinder.py
=======================
"Which file?" -- a folder the theory side sent, every file in it said in
plain words: what it is, what it can show (E-k dispersion, constant-energy
contours, the 3-D Fermi surface), what reading it still lacks, and whether
it is a standard file under another name or one file cut into pieces.

For the experimentalist who knows what they want to compare but not which
of a dozen ``.dat`` files gives it. The judging is
:func:`tools.crystal.bandfiles.survey`, by content, not by name; this is
its window. A file chosen here is opened by the Calculated bands window;
pieces are joined into one file first.
"""
from __future__ import annotations

import os

from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QColor
from PyQt5.QtWidgets import (QAbstractItemView, QDialog, QFileDialog,
                             QHBoxLayout, QHeaderView, QLabel, QMessageBox,
                             QPushButton, QTableWidget, QTableWidgetItem,
                             QVBoxLayout)

from tools.crystal import bandfiles

__all__ = ["WhichFileDialog", "GOAL_NAMES"]

#: the goals as the window says them
GOAL_NAMES = {"dispersion": "E–k dispersion",
              "contours": "constant-energy contours",
              "fs": "3D Fermi surface"}

_GREEN, _AMBER, _GREY = QColor("#e3f4e1"), QColor("#fdf1d8"), QColor("#f0f0f0")


class WhichFileDialog(QDialog):
    """The files of a folder, judged for ``goal``; :attr:`chosen` carries
    the path to open."""

    chosen = pyqtSignal(str)

    def __init__(self, folder: str = "", goal: str = "dispersion", parent=None):
        super().__init__(parent)
        self.setWindowTitle("Which file is the calculation?")
        self.resize(1000, 560)
        self.goal = goal
        self.folder = ""
        self.rows = []
        layout = QVBoxLayout(self)
        top = QHBoxLayout()
        self.folder_label = QLabel("")
        self.folder_label.setWordWrap(True)
        top.addWidget(self.folder_label, 1)
        choose = QPushButton("Choose folder...")
        choose.clicked.connect(lambda: self.choose_folder())
        top.addWidget(choose)
        layout.addLayout(top)
        self.intro = QLabel("")
        self.intro.setWordWrap(True)
        self.intro.setTextFormat(Qt.RichText)
        self.intro.setStyleSheet("background: #eef4f8; border: 1px solid #c9dbe6;"
                                 " border-radius: 4px; padding: 6px;")
        layout.addWidget(self.intro)
        self.table = QTableWidget(0, 4)
        self.table.setHorizontalHeaderLabels(["File", "What it is", "Can show",
                                              "Notes"])
        header = self.table.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(1, QHeaderView.Stretch)
        header.setSectionResizeMode(2, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(3, QHeaderView.Stretch)
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SingleSelection)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table.setWordWrap(True)
        self.table.itemSelectionChanged.connect(self._selected)
        self.table.cellDoubleClicked.connect(lambda *_: self.open_selected())
        layout.addWidget(self.table, 1)
        self.missing_label = QLabel("")
        self.missing_label.setWordWrap(True)
        self.missing_label.setTextFormat(Qt.RichText)
        layout.addWidget(self.missing_label)
        buttons = QHBoxLayout()
        self.open_button = QPushButton("Open the selected file")
        self.open_button.clicked.connect(lambda: self.open_selected())
        self.join_button = QPushButton("Join the pieces, then open")
        self.join_button.setToolTip("Writes the pieces one after the other into "
                                    "one file beside them, and opens that.")
        self.join_button.clicked.connect(lambda: self.join_selected())
        close = QPushButton("Close")
        close.setAutoDefault(False)
        close.clicked.connect(self.reject)
        buttons.addWidget(self.open_button)
        buttons.addWidget(self.join_button)
        buttons.addStretch(1)
        buttons.addWidget(close)
        layout.addLayout(buttons)
        self.open_button.setEnabled(False)
        self.join_button.setEnabled(False)
        if folder and os.path.isdir(folder):
            self.set_folder(folder)
        else:
            self.folder_label.setText("<b>No folder yet.</b> Choose the folder "
                                      "the theory side sent.")
            self._say_intro()

    # -- the folder -----------------------------------------------------------
    def choose_folder(self):
        folder = QFileDialog.getExistingDirectory(self, "The calculation's folder",
                                                  self.folder or "")
        if folder:
            self.set_folder(folder)

    def _say_intro(self):
        self.intro.setText(
            f"You want to see: <b>{GOAL_NAMES.get(self.goal, self.goal)}</b>. "
            "Every file of the folder is read (its beginning, not its name) "
            "and said in plain words. <span style='background:#e3f4e1'>Green"
            "</span>: opens for what you want. <span style='background:"
            "#fdf1d8'>Amber</span>: could, but something is missing or it is "
            "in pieces. <span style='background:#f0f0f0'>Grey</span>: a helper "
            "file the calculation reads by itself, or not data. Select one "
            "and press <b>Open</b> (or double-click it).")

    def set_folder(self, folder: str):
        self.folder = folder
        self.folder_label.setText(f"<b>Folder:</b> {folder}")
        self._say_intro()
        try:
            self.rows = bandfiles.survey(folder)
        except OSError as exc:
            self.rows = []
            self.missing_label.setText(f"<b>Cannot read the folder:</b> {exc}")
        self.table.setRowCount(len(self.rows))
        best = None
        for i, row in enumerate(self.rows):
            serves = self.goal in row["goals"]
            usable = serves and not row["missing"] and not row["pieces"]
            shows = ", ".join(GOAL_NAMES[g] for g in bandfiles.GOALS if g in row["goals"])
            notes = ([f"needs {m}" for m in row["missing"]]
                     + row["notes"]
                     + [f"(better with {o})" for o in row["optional"]])
            if row["goals"] and not serves:
                notes.insert(0, f"cannot show {GOAL_NAMES.get(self.goal, self.goal)}")
            cells = [row["name"], row["what"], shows or "—", "; ".join(notes)]
            colour = (_GREEN if usable else _AMBER if (serves or row["pieces"])
                      else _GREY)
            for column, text in enumerate(cells):
                item = QTableWidgetItem(text)
                item.setBackground(colour)
                item.setToolTip(text)
                self.table.setItem(i, column, item)
            if usable and best is None:
                best = i
        self.table.resizeRowsToContents()
        lines = []
        if not any(self.goal in r["goals"] for r in self.rows):
            lines.append(f"<b>No file here can show "
                         f"{GOAL_NAMES.get(self.goal, self.goal)}.</b>"
                         + (" Contours and the 3D Fermi surface need the bands "
                            "everywhere in the zone: a band grid (.bxsf, a 3-D "
                            "band.dat, a MATLAB Ek_3D .mat) or a Hamiltonian "
                            "(*_hr.dat)." if self.goal != "dispersion" else ""))
        self.missing_label.setText("<br>".join(lines))
        if best is not None:
            self.table.selectRow(best)
        self._selected()

    # -- choosing -------------------------------------------------------------
    def _current(self):
        rows = self.table.selectionModel().selectedRows() if self.table.selectionModel() else []
        return self.rows[rows[0].row()] if rows else None

    def _selected(self):
        row = self._current()
        self.open_button.setEnabled(bool(row and row["goals"] and not row["pieces"]))
        self.join_button.setEnabled(bool(row and row["pieces"]))

    def open_selected(self):
        row = self._current()
        if row is None or not row["goals"]:
            return None
        if row["pieces"]:
            return self.join_selected()
        self.chosen.emit(row["path"])
        self.accept()
        return row["path"]

    def join_selected(self, target: str = None):
        row = self._current()
        if row is None or not row["pieces"]:
            return None
        pieces = row["pieces"]
        name = next(j for j, members in bandfiles.split_groups(
            [os.path.basename(p) for p in pieces]))
        if target is None:
            target = os.path.join(self.folder, name)
            if os.path.exists(target):
                target = os.path.join(self.folder, "joined_" + name)
        try:
            bandfiles.join_pieces(pieces, target)
        except OSError:
            target, _ = QFileDialog.getSaveFileName(
                self, "Save the joined file", name)
            if not target:
                return None
            try:
                bandfiles.join_pieces(pieces, target)
            except OSError as exc:
                QMessageBox.warning(self, "Join the pieces", str(exc))
                return None
        self.chosen.emit(target)
        self.accept()
        return target
