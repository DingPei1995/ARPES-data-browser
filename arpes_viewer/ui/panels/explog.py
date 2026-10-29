"""
ui/panels/explog.py
===================
The experiment log, from the tool panel's Datasets box: two buttons side by
side.

* **Exp. log format...** -- the window that sets the log up, once: the
  folder and file name, the separator and decimal mark, and the columns.
  A dataset of the main list is picked as the reference and its fields are
  listed: the file name, the kind and the scanned ranges worked out from the
  data, then every metadata field and motor position the file has, each
  with its value there. Ticking a field adds a column; the columns can be
  renamed, put in any order, and told when a scanned axis's range stands in
  for the field's value (a SPEM scan's stage, a manipulator map's angle --
  see :data:`tools.output.explog.SCAN_RULES`). Pointing it at a log that
  already exists shows *that log's* format, to adjust; saving a changed
  format rewrites the log under its new header, so the header stays the
  only first row. *Keep recording* logs every file added to the main list
  from then on.
* **Auto Exp. log** -- one click: the selected datasets' records are added
  to the log, in its format, after its last record. The first click opens
  the format window when no log has been set up yet.

:mod:`tools.output.explog` does the work and is Qt-free; this module adds
the windows, and what Excel uses as a separator on this computer
(:func:`excel_conventions`).
"""
from __future__ import annotations

import os
import sys
import traceback

from PyQt5.QtCore import QLocale, Qt
from PyQt5.QtWidgets import (QAbstractItemView, QCheckBox, QComboBox, QDialog,
                             QDialogButtonBox, QFileDialog, QGridLayout,
                             QGroupBox, QHBoxLayout, QHeaderView, QLabel,
                             QLineEdit, QListWidget, QListWidgetItem,
                             QMessageBox, QPushButton, QSplitter, QTableWidget,
                             QTableWidgetItem, QVBoxLayout, QWidget)

from arpes_viewer.tools.output import explog as EL
from arpes_viewer.ui.common.screen import fits_screen

__all__ = ["ExpLogFormatDialog", "excel_conventions", "record_datasets"]

TITLE = "Experiment log"


def excel_conventions():
    """``(separator, decimal mark)`` Excel uses for a CSV on this computer.

    Excel splits a CSV at the system's *list separator*, not at commas: on
    a Windows set to French, German, Italian... that is ``;``, and a comma
    file opens with every row in the first column. Windows keeps both in
    the registry; Excel for Mac (and LibreOffice) pairs them with the
    system's decimal mark -- ``;`` where the decimal mark is a comma."""
    if sys.platform == "win32":
        try:
            import winreg
            with winreg.OpenKey(winreg.HKEY_CURRENT_USER,
                                r"Control Panel\International") as key:
                separator = winreg.QueryValueEx(key, "sList")[0]
                decimal = winreg.QueryValueEx(key, "sDecimal")[0]
            if separator in (",", ";", "\t") and decimal in (".", ","):
                return separator, decimal
        except OSError:
            pass
    point = QLocale.system().decimalPoint()
    return (";", ",") if point == "," else (",", ".")


def _describe(delimiter: str) -> str:
    return {",": "comma", ";": "semicolon", "\t": "tab"}.get(delimiter,
                                                              delimiter)


@fits_screen
class ExpLogFormatDialog(QDialog):
    """``entries()`` lists the main list as ``[(key, row text), ...]``;
    ``load(key)`` reads one as ``(data, record)`` (the data is closed here
    when done); ``selected()`` is the list's selection, as keys.

    Works on a copy of the settings; :meth:`accept` (Save format) is what
    writes them, and reformats an existing log whose format changed."""

    def __init__(self, parent=None, *, entries, load, selected, folder="",
                 keep_recording=False, system=None):
        super().__init__(parent, Qt.Window)
        self.setWindowTitle("Exp. log format")
        self.resize(1000, 760)
        self._entries, self._load, self._selected = entries, load, selected
        self.system = system or excel_conventions()
        self.settings = EL.load_settings()
        if not self.settings.folder:
            self.settings.folder = folder
        self.template = self.settings.template.copy()
        #: the format of the log on disk, when there is one with records
        self.file_template = None
        self.file_log = None
        self.reference = None            # (data, record) of the reference
        self._filling = False

        layout = QVBoxLayout(self)
        layout.addWidget(self._file_box())
        layout.addWidget(self._format_box())
        layout.addWidget(self._columns_box(), stretch=1)
        self.skip_box = QCheckBox("Skip datasets already in the log (by file "
                                  "name)")
        self.skip_box.setChecked(self.settings.skip_logged)
        self.auto_box = QCheckBox(
            "Keep recording: log each file added to the main list from now on")
        self.auto_box.setChecked(keep_recording)
        self.auto_box.setToolTip(
            "During a beamtime: add each new file to the list as it is "
            "measured (Load data or Quick add) and its record is written "
            "straight away, in this format. Until the program is closed.")
        layout.addWidget(self.skip_box)
        layout.addWidget(self.auto_box)
        buttons = QDialogButtonBox(QDialogButtonBox.Save
                                   | QDialogButtonBox.Cancel)
        buttons.button(QDialogButtonBox.Save).setText("Save format")
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

        self._sync_to_file()
        self.refresh()

    # -- building ---------------------------------------------------------
    def _file_box(self):
        box = QGroupBox("Log file")
        grid = QGridLayout(box)
        self.folder_edit = QLineEdit(self.settings.folder)
        self.folder_edit.setPlaceholderText("Folder the log is written to")
        browse = QPushButton("Folder...")
        browse.clicked.connect(self._browse)
        self.name_edit = QLineEdit(self.settings.name)
        self.name_edit.setPlaceholderText("Log file name (.csv)")
        self.name_edit.setToolTip(
            "If a log of this name is already in the folder, its own format "
            "is shown below to adjust, and new records go after its last "
            "one.")
        self.file_status = QLabel("")
        self.file_status.setWordWrap(True)
        grid.addWidget(QLabel("Folder:"), 0, 0)
        grid.addWidget(self.folder_edit, 0, 1)
        grid.addWidget(browse, 0, 2)
        grid.addWidget(QLabel("File name:"), 1, 0)
        grid.addWidget(self.name_edit, 1, 1, 1, 2)
        grid.addWidget(self.file_status, 2, 0, 1, 3)
        self.folder_edit.editingFinished.connect(self._sync_to_file)
        self.name_edit.editingFinished.connect(self._sync_to_file)
        return box

    def _format_box(self):
        box = QGroupBox("Format")
        row = QHBoxLayout(box)
        row.addWidget(QLabel("Separator:"))
        self.delimiter_combo = QComboBox()
        for value, text in EL.DELIMITERS.items():
            if value == "auto":
                text = (f"auto: {_describe(self.system[0])}, as Excel on "
                        f"this computer")
            self.delimiter_combo.addItem(text, value)
        self.delimiter_combo.setToolTip(
            "Excel opens a CSV with the computer's own list separator: a "
            "semicolon where the decimal mark is a comma (much of Europe), "
            "a comma elsewhere. A file with the other one opens with every "
            "row in the first column. LibreOffice and Numbers ask, or "
            "detect it.")
        row.addWidget(self.delimiter_combo)
        row.addWidget(QLabel("Decimal mark:"))
        self.decimal_combo = QComboBox()
        for value, text in EL.DECIMALS.items():
            if value == "auto":
                text = "auto: comma beside a semicolon"
            self.decimal_combo.addItem(text, value)
        self.decimal_combo.setToolTip(
            "Numbers are written with this decimal mark, so a spreadsheet "
            "using a decimal comma reads them as numbers (and does not turn "
            "1.5 into a date). Names and text are never changed.")
        row.addWidget(self.decimal_combo)
        row.addWidget(QLabel("Ranges as:"))
        self.style_combo = QComboBox()
        for style, text in EL.RANGE_STYLES.items():
            self.style_combo.addItem(text, style)
        row.addWidget(self.style_combo)
        row.addStretch(1)
        self.delimiter_combo.currentIndexChanged.connect(self._format_changed)
        self.decimal_combo.currentIndexChanged.connect(self._format_changed)
        self.style_combo.currentIndexChanged.connect(self._format_changed)
        self._show_format()
        return box

    def _columns_box(self):
        box = QGroupBox("Columns")
        column = QVBoxLayout(box)

        top = QHBoxLayout()
        top.addWidget(QLabel("Fields of:"))
        self.reference_combo = QComboBox()
        self.reference_combo.setToolTip(
            "A dataset of the main list whose fields the columns are chosen "
            "from. Pick one of the kind you will log most (a map, a cut, a "
            "SPEM scan); the columns then apply to every dataset logged.")
        self.reference_combo.setSizeAdjustPolicy(
            QComboBox.AdjustToMinimumContentsLengthWithIcon)
        self.reference_combo.setMinimumContentsLength(30)
        top.addWidget(self.reference_combo, stretch=1)
        read = QPushButton("Read fields")
        read.clicked.connect(self.read_reference)
        top.addWidget(read)
        column.addLayout(top)

        splitter = QSplitter(Qt.Horizontal)
        left = QWidget()
        left_layout = QVBoxLayout(left)
        left_layout.setContentsMargins(0, 0, 0, 0)
        left_layout.addWidget(QLabel("Tick the fields to log:"))
        self.field_filter = QLineEdit()
        self.field_filter.setPlaceholderText(
            "Filter fields (e.g. 'polar', 'pass', 'temp')")
        self.field_filter.textChanged.connect(self._apply_filter)
        left_layout.addWidget(self.field_filter)
        self.field_list = QListWidget()
        self.field_list.itemChanged.connect(self._field_ticked)
        left_layout.addWidget(self.field_list, stretch=1)
        splitter.addWidget(left)

        right = QWidget()
        right_layout = QVBoxLayout(right)
        right_layout.setContentsMargins(0, 0, 0, 0)
        right_layout.addWidget(QLabel(
            "The log's columns, in order (double-click a heading to rename "
            "it):"))
        self.table = QTableWidget(0, 4)
        self.table.setHorizontalHeaderLabels(
            ["Heading", "Field", "If scanned", "Value in the reference"])
        self.table.verticalHeader().setVisible(False)
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SingleSelection)
        header = self.table.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.Interactive)
        header.setSectionResizeMode(1, QHeaderView.Interactive)
        header.setSectionResizeMode(2, QHeaderView.ResizeToContents)
        header.setStretchLastSection(True)
        self.table.setColumnWidth(0, 110)
        self.table.setColumnWidth(1, 130)
        self.table.itemChanged.connect(self._heading_edited)
        right_layout.addWidget(self.table, stretch=1)
        buttons = QHBoxLayout()
        for text, slot, tip in (
                ("Top", lambda: self.move_column(-len(self.template.columns)),
                 "Move the column to the front."),
                ("Up", lambda: self.move_column(-1), "Move the column one left."),
                ("Down", lambda: self.move_column(1), "Move the column one right."),
                ("Bottom", lambda: self.move_column(len(self.template.columns)),
                 "Move the column to the end."),
                ("Remove", self.remove_column, "Take the column out."),
                ("Add empty column", self.add_blank,
                 "A column left empty, to fill in by hand: comments, the "
                 "sample's state, what was changed before this file.")):
            button = QPushButton(text)
            button.setToolTip(tip)
            button.clicked.connect(slot)
            buttons.addWidget(button)
        right_layout.addLayout(buttons)
        splitter.addWidget(right)
        splitter.setSizes([330, 650])
        column.addWidget(splitter, stretch=1)

        template_row = QHBoxLayout()
        for text, slot, tip in (
                ("Save format as...", self.save_template_as,
                 "Keep this format in a file of your own -- one per "
                 "beamline, say -- to load again for a new log."),
                ("Load format...", self.load_template_from,
                 "Use a format saved earlier with Save format as."),
                ("Default columns", self.reset_template,
                 "Start again from the default columns.")):
            button = QPushButton(text)
            button.setToolTip(tip)
            button.clicked.connect(slot)
            template_row.addWidget(button)
        template_row.addStretch(1)
        column.addLayout(template_row)
        return box

    # -- the log file -----------------------------------------------------
    def log_path(self) -> str:
        self.settings.folder = self.folder_edit.text().strip()
        self.settings.name = self.name_edit.text().strip()
        return self.settings.path()

    def _browse(self):
        folder = QFileDialog.getExistingDirectory(
            self, "Folder for the experiment log", self.folder_edit.text())
        if folder:
            self.folder_edit.setText(folder)
            self._sync_to_file()

    def _sync_to_file(self):
        """Show the format of the log the folder and name point at: for a
        log that already has records, its own (remembered, or read off its
        header); for a new one, the format being edited."""
        path = self.log_path()
        self.file_log, self.file_template = None, None
        if path:
            try:
                log = EL.read_log(path)
            except (OSError, UnicodeError):
                log = None
            if log is not None and log.header:
                known = self.settings.format_for(path, log.header)
                if known is None and self.template.headers() == log.header:
                    known = self.template.copy()
                if known is None:
                    known = EL.template_from_header(log.header)
                known.delimiter = log.delimiter
                self.file_log, self.file_template = log, known
                self.template = known.copy()
                self._show_format()
                self._fill_fields()
        self._fill_columns()

    def _update_file_status(self):
        path = self.settings.path()
        if not path:
            self.file_status.setText("Give the log a file name.")
            return
        if self.file_log is None:
            self.file_status.setText(
                f"New log: {path} is created at the first record, with these "
                f"columns as its header.")
            return
        changed = self.template.headers() != self.file_log.header or (
            self.template.delimiter not in ("auto", self.file_log.delimiter))
        text = (f"Existing log, {len(self.file_log.records)} record(s), "
                f"{_describe(self.file_log.delimiter)}-separated: below is "
                f"its format. New records go after the last one.")
        if changed:
            text += (" The format was changed: Save format rewrites the log "
                     "under the new header, its records moved with their "
                     "columns.")
        self.file_status.setText(text)

    # -- the reference dataset and its fields -----------------------------
    def refresh(self):
        """Fill the reference list from the main list, picking the
        selection."""
        selected = list(self._selected() or ())
        self.reference_combo.blockSignals(True)
        self.reference_combo.clear()
        for key, text in self._entries():
            self.reference_combo.addItem(text, key)
        want = selected[0] if selected else None
        index = self.reference_combo.findData(want) if want is not None else -1
        self.reference_combo.setCurrentIndex(max(0, index))
        self.reference_combo.blockSignals(False)
        if self.reference_combo.count():
            self.read_reference()
        else:
            self._fill_fields()
            self._fill_columns()

    def _close_reference(self):
        if self.reference is not None:
            try:
                self.reference[0].close()
            except Exception:                              # noqa: BLE001
                pass
            self.reference = None

    def read_reference(self):
        """Read the reference dataset and list its fields."""
        key = self.reference_combo.currentData()
        if key is None:
            return
        self._close_reference()
        try:
            self.reference = self._load(key)
        except Exception as exc:                           # noqa: BLE001
            traceback.print_exc()
            self.reference = None
            QMessageBox.warning(self, TITLE,
                                f"Could not read the dataset:\n{exc}")
        if self.reference is not None and self.reference[0] is None:
            self.reference = None
        self._fill_fields()
        self._fill_columns()

    def _fill_fields(self):
        self._filling = True
        self.field_list.clear()
        chosen = {c.field for c in self.template.columns}
        listed = set()
        fields = EL.available_fields(*self.reference) if self.reference else \
            [(fid, name, "") for fid, (name, _h) in EL.DERIVED.items()]
        for fid, name, value in fields:
            listed.add(fid)
            self._add_field_item(fid, name, value, fid in chosen)
        # Columns of the format that this dataset does not have.
        for c in self.template.columns:
            if c.field not in listed:
                listed.add(c.field)
                self._add_field_item(c.field, c.field.split(":", 1)[-1],
                                     "(not in this dataset)", True)
        self._filling = False
        self._apply_filter(self.field_filter.text())

    def _add_field_item(self, fid, name, value, checked):
        text = f"{name}  =  {value}" if value != "" else name
        item = QListWidgetItem(text)
        item.setData(Qt.UserRole, fid)
        item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
        item.setCheckState(Qt.Checked if checked else Qt.Unchecked)
        item.setToolTip(fid)
        self.field_list.addItem(item)

    def _apply_filter(self, text):
        text = (text or "").strip().lower()
        for i in range(self.field_list.count()):
            item = self.field_list.item(i)
            item.setHidden(bool(text) and text not in item.text().lower()
                           and text not in str(item.data(Qt.UserRole)).lower())

    def _field_ticked(self, item):
        if self._filling:
            return
        fid = item.data(Qt.UserRole)
        if item.checkState() == Qt.Checked:
            if fid == "@blank" or fid not in {c.field for c in self.template.columns}:
                self.template.columns.append(
                    EL.Column(EL.default_header(fid), fid))
        else:
            self.template.columns = [c for c in self.template.columns
                                     if c.field != fid]
        self._fill_columns()

    # -- the format -------------------------------------------------------
    def _show_format(self):
        for combo, value in ((self.delimiter_combo, self.template.delimiter),
                             (self.decimal_combo, self.template.decimal),
                             (self.style_combo, self.template.range_style)):
            combo.blockSignals(True)
            combo.setCurrentIndex(max(0, combo.findData(value)))
            combo.blockSignals(False)

    def _format_changed(self, *_):
        self.template.delimiter = self.delimiter_combo.currentData()
        self.template.decimal = self.decimal_combo.currentData()
        self.template.range_style = self.style_combo.currentData()
        self._fill_columns()

    def _decimal(self) -> str:
        log = self.file_log or EL.LogFile(self.settings.path())
        return EL.resolve(self.template, log, self.system)[1]

    # -- the columns ------------------------------------------------------
    def _fill_columns(self, select: int = None):
        self._filling = True
        # The old rows' rule boxes are only deleted later, by the event
        # loop; hidden now, so a shorter table shows no stray boxes.
        for index in range(self.table.rowCount()):
            widget = self.table.cellWidget(index, 2)
            if widget is not None:
                widget.hide()
                self.table.removeCellWidget(index, 2)
        self.table.setRowCount(len(self.template.columns))
        decimal = self._decimal()
        for index, column in enumerate(self.template.columns):
            self.table.setItem(index, 0, QTableWidgetItem(column.header))
            name = EL.DERIVED[column.field][0] if column.field in EL.DERIVED \
                else column.field.split(":", 1)[-1]
            source = QTableWidgetItem(name)
            source.setFlags(source.flags() & ~Qt.ItemIsEditable)
            source.setToolTip(column.field)
            self.table.setItem(index, 1, source)
            rules = QComboBox()
            for rule, text in EL.SCAN_RULES.items():
                rules.addItem(text, rule)
            rules.setCurrentIndex(max(0, rules.findData(column.scan)))
            rules.setEnabled(not column.field.startswith("@"))
            rules.setToolTip(
                "When a scanned axis's range (start/stop/step) is written "
                "instead of this field's value. auto: when this field names "
                "the scanned motor, or its value changed along the scan. The "
                "others apply to that kind of measurement only; any other "
                "dataset gets the field's own value.")
            rules.currentIndexChanged.connect(
                lambda _i, c=column, box=rules: self._rule_changed(c, box))
            self.table.setCellWidget(index, 2, rules)
            example = self._example(column, decimal)
            value = QTableWidgetItem(example)
            value.setFlags(value.flags() & ~Qt.ItemIsEditable)
            value.setToolTip(example)
            self.table.setItem(index, 3, value)
        self._filling = False
        if select is not None and 0 <= select < self.table.rowCount():
            self.table.selectRow(select)
        self._update_file_status()

    def _example(self, column, decimal=".") -> str:
        if self.reference is None:
            return ""
        try:
            return EL.cell(column, self.reference[0], self.reference[1],
                           self.template.range_style, decimal)
        except Exception as exc:                           # noqa: BLE001
            return f"<{exc}>"

    def _heading_edited(self, item):
        if self._filling or item.column() != 0:
            return
        row = item.row()
        if row < len(self.template.columns):
            self.template.columns[row].header = item.text().strip()
            self._update_file_status()

    def _rule_changed(self, column, box):
        column.scan = box.currentData()
        row = self.template.columns.index(column)
        item = self.table.item(row, 3)
        if item is not None:
            item.setText(self._example(column, self._decimal()))

    def current_row(self) -> int:
        rows = self.table.selectionModel().selectedRows()
        return rows[0].row() if rows else self.table.currentRow()

    def move_column(self, by: int):
        """Move the selected column ``by`` places (negative: to the left)."""
        row = self.current_row()
        columns = self.template.columns
        if not 0 <= row < len(columns):
            return
        target = min(max(row + by, 0), len(columns) - 1)
        columns.insert(target, columns.pop(row))
        self._fill_columns(select=target)

    def remove_column(self):
        row = self.current_row()
        if 0 <= row < len(self.template.columns):
            del self.template.columns[row]
            self._fill_fields()
            self._fill_columns(select=min(row, len(self.template.columns) - 1))

    def add_blank(self):
        self.template.columns.append(EL.Column("Comments", "@blank"))
        self._fill_fields()
        self._fill_columns(select=len(self.template.columns) - 1)

    def set_template(self, template):
        """Use ``template``'s columns and format. For an existing log, its
        columns are followed by heading where they keep one."""
        self.template = template
        self._show_format()
        self._fill_fields()
        self._fill_columns()

    def save_template_as(self):
        path, _ = QFileDialog.getSaveFileName(
            self, "Save the log's format", "explog_format.json",
            "Format (*.json)")
        if not path:
            return
        try:
            EL.save_template(self.template, path)
        except OSError as exc:
            QMessageBox.warning(self, TITLE,
                                f"Could not save the format:\n{exc}")

    def load_template_from(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Load a log format", "", "Format (*.json)")
        if path:
            self.set_template(EL.load_template(path))

    def reset_template(self):
        self.set_template(EL.default_template())

    # -- saving -----------------------------------------------------------
    def save(self) -> bool:
        """Keep the format; reformat the log on disk if it changed. False
        (with the reason shown) when it could not be done or was called
        off."""
        path = self.log_path()
        if not path:
            QMessageBox.warning(self, TITLE, "Give the log a file name.")
            return False
        if not self.template.columns:
            QMessageBox.warning(self, TITLE,
                                "The log has no columns: tick some fields.")
            return False
        if len(set(self.template.headers())) != len(self.template.columns):
            QMessageBox.warning(self, TITLE,
                                "Two columns have the same heading; rename "
                                "one of them.")
            return False
        log = self.file_log
        if log is not None and (
                self.template.headers() != log.header
                or self.template.delimiter not in ("auto", log.delimiter)):
            _sources, lost = EL.reformat_plan(log, self.file_template,
                                              self.template)
            question = (f"{os.path.basename(path)} already holds "
                        f"{len(log.records)} record(s). Rewrite it in the new "
                        f"format? Its records move with their columns; new "
                        f"columns are left empty for them.")
            if lost:
                question += ("\n\nThese columns are no longer in the format, "
                             "and what they hold is removed from the log: "
                             + ", ".join(lost))
            answer = QMessageBox.question(
                self, TITLE, question, QMessageBox.Yes | QMessageBox.No,
                QMessageBox.No)
            if answer != QMessageBox.Yes:
                return False
            try:
                EL.reformat(path, self.file_template, self.template,
                            self.system)
            except PermissionError as exc:
                QMessageBox.warning(
                    self, TITLE, f"Could not rewrite the log -- is it open in "
                    f"Excel? Close it there and save again.\n{exc}")
                return False
            except OSError as exc:
                QMessageBox.warning(self, TITLE,
                                    f"Could not rewrite the log:\n{exc}")
                return False
        if log is not None and self.template.delimiter == "auto":
            self.template.delimiter = EL.read_log(path).delimiter
        self.settings.template = self.template.copy()
        self.settings.skip_logged = self.skip_box.isChecked()
        self.settings.configured = True
        self.settings.remember_format(path, self.template)
        EL.save_settings(self.settings)
        return True

    def accept(self):
        if self.save():
            self._close_reference()
            super().accept()

    def reject(self):
        self._close_reference()
        super().reject()


def record_datasets(keys, load, *, parent=None, status=None, quiet=False,
                    system=None) -> int:
    """Add the records of ``keys`` to the experiment log, in its saved
    format (Auto Exp. log, and Keep recording). Returns how many were
    written. ``quiet`` reports in the status bar only."""
    status = status or (lambda text: None)

    def tell(text, warning=False):
        status(text.replace("\n", " "))
        if warning and not quiet:
            QMessageBox.warning(parent, TITLE, text)

    settings = EL.load_settings()
    path = settings.path()
    if not path or not settings.template.columns:
        tell("Set the experiment log up first (Exp. log format...).", True)
        return 0
    keys = list(keys or ())
    if not keys:
        tell("Select the datasets to log in the main list.", True)
        return 0
    system = system or excel_conventions()
    try:
        log = EL.read_log(path)
    except (OSError, UnicodeError) as exc:
        tell(f"Could not read the log {path}:\n{exc}", True)
        return 0
    template = settings.template
    if log.header and log.header != template.headers():
        template = settings.format_for(path, log.header)
        if template is None:
            tell(f"The columns of {os.path.basename(path)} are not the "
                 f"format's (was it edited by hand?). Open Exp. log format... "
                 f"to adjust it.", True)
            return 0
    _delimiter, decimal = EL.resolve(template, log, system)
    done = EL.logged_names(log, template) if settings.skip_logged else set()
    file_column = next((i for i, c in enumerate(template.columns)
                        if c.field == "@file"), None)
    rows, skipped, failed = [], 0, []
    for key in keys:
        try:
            data, record = load(key)
        except Exception as exc:                           # noqa: BLE001
            traceback.print_exc()
            failed.append(f"{key[0]}: {exc}")
            continue
        if data is None:
            failed.append(str(key[0]))
            continue
        try:
            cells = EL.row(template, data, record, decimal)
        finally:
            data.close()
        if file_column is not None:
            name = cells[file_column].strip()
            if name in done:
                skipped += 1
                continue
            done.add(name)
        rows.append(cells)
    written = 0
    if rows:
        try:
            written = EL.append(path, template, rows, system)
        except PermissionError as exc:
            tell(f"Could not write the log -- is it open in Excel? Close it "
                 f"there and log again.\n{exc}", True)
            return 0
        except (OSError, EL.FormatMismatch) as exc:
            tell(f"Could not write the log:\n{exc}", True)
            return 0
    message = f"Experiment log: {written} record(s) added to {path}"
    if skipped:
        message += f", {skipped} already in it"
    if failed:
        message += f", {len(failed)} could not be read"
    tell(message)
    if failed and not quiet:
        QMessageBox.warning(parent, TITLE,
                            "Could not read:\n" + "\n".join(failed[:10]))
    return written
