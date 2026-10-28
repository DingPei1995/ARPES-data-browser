"""
ui/panels/explog.py
===================
The Auto Exp. log window: write the experiment log (.csv) from the listed
datasets' own metadata.

Opened from the tool panel's Datasets box (**Auto Exp. log...**). Three
parts, top to bottom:

* **Log file** -- the folder and the file name. A file that is already
  there is continued after its last record; one that is not is created with
  the columns as its header.
* **Columns** -- what each row holds. A dataset of the main list is picked as
  the reference and its fields are listed: the file name, the kind and the
  scanned ranges worked out from the data, then every metadata field and
  motor position the file has, each with its value there. Ticking a field
  adds a column; the columns can be renamed, moved up and down, and told
  when a scanned axis's range stands in for the field's value (a SPEM
  scan's stage, a manipulator map's angle -- see
  :data:`tools.output.explog.SCAN_RULES`). The columns are remembered for
  the next session, and can be saved to a file and loaded from one, to keep
  one per beamline.
* **Record** -- which datasets (the selection, or the whole list) and
  whether to skip those already in the log. *Keep recording* logs every file
  added to the list from then on, which is the beamtime use: load each new
  file as it is measured and the log keeps itself.

:mod:`tools.output.explog` does the work and is Qt-free.
"""
from __future__ import annotations

import os
import traceback

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (QAbstractItemView, QCheckBox, QComboBox, QDialog,
                             QFileDialog, QGridLayout, QGroupBox, QHBoxLayout,
                             QHeaderView, QLabel, QLineEdit, QListWidget,
                             QListWidgetItem, QMessageBox, QPushButton,
                             QRadioButton, QSplitter, QTableWidget,
                             QTableWidgetItem, QVBoxLayout, QWidget)

from tools.output import explog as EL
from ui.common.screen import fits_screen

__all__ = ["ExpLogDialog"]

DEFAULT_NAME = "exp_log.csv"


@fits_screen
class ExpLogDialog(QDialog):
    """``entries()`` lists the main list as ``[(key, row text), ...]``;
    ``load(key)`` reads one as ``(data, record)`` (the data is closed here
    when done); ``selected()`` is the list's selection, as keys; ``status``
    shows a line in the main window's status bar."""

    def __init__(self, parent=None, *, entries, load, selected, folder="",
                 status=None):
        super().__init__(parent, Qt.Window)
        self.setWindowTitle("Auto Exp. log")
        self.resize(980, 720)
        self._entries, self._load, self._selected = entries, load, selected
        self._status = status or (lambda text: None)
        self.template = EL.load_template()
        self.reference = None            # (data, record) of the reference
        self._filling = False

        layout = QVBoxLayout(self)
        layout.addWidget(self._file_box(folder))
        layout.addWidget(self._columns_box(), stretch=1)
        layout.addWidget(self._record_box())
        self._fill_columns()
        self._update_file_status()

    # -- building ---------------------------------------------------------
    def _file_box(self, folder):
        box = QGroupBox("Log file")
        grid = QGridLayout(box)
        self.folder_edit = QLineEdit(folder)
        self.folder_edit.setPlaceholderText("Folder the log is written to")
        browse = QPushButton("Folder...")
        browse.clicked.connect(self._browse)
        self.name_edit = QLineEdit(DEFAULT_NAME)
        self.name_edit.setPlaceholderText("Log file name (.csv)")
        self.name_edit.setToolTip(
            "If a log of this name is already in the folder, the new records "
            "are written after its last one; nothing in it is changed.")
        self.file_status = QLabel("")
        self.file_status.setWordWrap(True)
        grid.addWidget(QLabel("Folder:"), 0, 0)
        grid.addWidget(self.folder_edit, 0, 1)
        grid.addWidget(browse, 0, 2)
        grid.addWidget(QLabel("File name:"), 1, 0)
        grid.addWidget(self.name_edit, 1, 1, 1, 2)
        grid.addWidget(self.file_status, 2, 0, 1, 3)
        self.folder_edit.textChanged.connect(self._update_file_status)
        self.name_edit.textChanged.connect(self._update_file_status)
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
        top.addWidget(QLabel("Ranges as:"))
        self.style_combo = QComboBox()
        for style, text in EL.RANGE_STYLES.items():
            self.style_combo.addItem(text, style)
        self.style_combo.setCurrentIndex(
            max(0, self.style_combo.findData(self.template.range_style)))
        self.style_combo.currentIndexChanged.connect(self._style_changed)
        top.addWidget(self.style_combo)
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
                ("Save columns as...", self.save_template_as,
                 "Keep these columns in a file of your own -- one per "
                 "beamline, say -- to load again later."),
                ("Load columns...", self.load_template_from,
                 "Use columns saved earlier with Save columns as."),
                ("Default columns", self.reset_template,
                 "Start again from the default columns.")):
            button = QPushButton(text)
            button.setToolTip(tip)
            button.clicked.connect(slot)
            template_row.addWidget(button)
        template_row.addStretch(1)
        column.addLayout(template_row)
        return box

    def _record_box(self):
        box = QGroupBox("Record")
        grid = QGridLayout(box)
        self.selected_radio = QRadioButton("The selected rows")
        self.all_radio = QRadioButton("Every dataset in the list")
        self.selected_radio.setChecked(True)
        self.skip_box = QCheckBox("Skip datasets already in the log")
        self.skip_box.setChecked(True)
        self.skip_box.setToolTip(
            "A dataset whose file name is already in the log's file-name "
            "column is not written again. Needs a File name column.")
        self.auto_box = QCheckBox(
            "Keep recording: log each file added to the main list from now on")
        self.auto_box.setToolTip(
            "During a beamtime: add each new file to the list as it is "
            "measured (Load data or Quick add) and its record is written "
            "straight away, with these columns, to this log.")
        self.write_button = QPushButton("Write records")
        self.write_button.setDefault(True)
        self.write_button.clicked.connect(lambda: self.record())
        close = QPushButton("Close")
        close.clicked.connect(self.close)
        grid.addWidget(self.selected_radio, 0, 0)
        grid.addWidget(self.all_radio, 0, 1)
        grid.addWidget(self.skip_box, 0, 2)
        grid.addWidget(self.auto_box, 1, 0, 1, 3)
        grid.addWidget(self.write_button, 0, 3)
        grid.addWidget(close, 1, 3)
        grid.setColumnStretch(2, 1)
        return box

    # -- the log file -----------------------------------------------------
    def log_path(self) -> str:
        """The log's full path, ``.csv`` added if the name has none; empty
        while no name is given."""
        name = self.name_edit.text().strip()
        if not name:
            return ""
        if not name.lower().endswith(".csv"):
            name += ".csv"
        folder = self.folder_edit.text().strip() or os.path.expanduser("~")
        return os.path.join(folder, name)

    def _browse(self):
        folder = QFileDialog.getExistingDirectory(
            self, "Folder for the experiment log", self.folder_edit.text())
        if folder:
            self.folder_edit.setText(folder)

    def _update_file_status(self, *_):
        path = self.log_path()
        if not path:
            self.file_status.setText("Give the log a file name.")
            return
        try:
            log = EL.read_log(path)
        except (OSError, UnicodeError) as exc:
            self.file_status.setText(f"Cannot read {path}: {exc}")
            return
        if not log.exists:
            self.file_status.setText(
                f"New log: {path} will be created, with the columns below as "
                f"its header.")
            return
        same = log.header == self.template.headers()
        self.file_status.setText(
            f"Existing log, {len(log.records)} record(s): new records go "
            f"after the last one."
            + ("" if same else " Its columns differ from the ones below, so "
               "a new header row is written before them."))

    # -- the reference dataset and its fields -----------------------------
    def refresh(self):
        """Fill the reference list from the main list again, keeping the
        current choice, or picking the selection."""
        current = self.reference_combo.currentData()
        selected = list(self._selected() or ())
        self.reference_combo.blockSignals(True)
        self.reference_combo.clear()
        for key, text in self._entries():
            self.reference_combo.addItem(text, key)
        want = selected[0] if selected else current
        index = self.reference_combo.findData(want) if want is not None else -1
        self.reference_combo.setCurrentIndex(max(0, index))
        self.reference_combo.blockSignals(False)
        if self.reference is None and self.reference_combo.count():
            self.read_reference()
        self._update_file_status()

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
            QMessageBox.warning(self, "Auto Exp. log",
                                f"Could not read the dataset:\n{exc}")
            return
        if self.reference is None or self.reference[0] is None:
            self.reference = None
            return
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
        # Columns of a remembered template that this dataset does not have.
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
        for index, column in enumerate(self.template.columns):
            heading = QTableWidgetItem(column.header)
            self.table.setItem(index, 0, heading)
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
            example = self._example(column)
            value = QTableWidgetItem(example)
            value.setFlags(value.flags() & ~Qt.ItemIsEditable)
            value.setToolTip(example)
            self.table.setItem(index, 3, value)
        self._filling = False
        if select is not None and 0 <= select < self.table.rowCount():
            self.table.selectRow(select)
        self._update_file_status()
        EL.save_template(self.template)

    def _example(self, column) -> str:
        if self.reference is None:
            return ""
        try:
            return EL.cell(column, self.reference[0], self.reference[1],
                           self.template.range_style)
        except Exception as exc:                           # noqa: BLE001
            return f"<{exc}>"

    def _heading_edited(self, item):
        if self._filling or item.column() != 0:
            return
        row = item.row()
        if row < len(self.template.columns):
            self.template.columns[row].header = item.text().strip()
            self._update_file_status()
            EL.save_template(self.template)

    def _rule_changed(self, column, box):
        column.scan = box.currentData()
        row = self.template.columns.index(column)
        item = self.table.item(row, 3)
        if item is not None:
            item.setText(self._example(column))
        EL.save_template(self.template)

    def _style_changed(self, *_):
        self.template.range_style = self.style_combo.currentData()
        self._fill_columns()

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
        self.template = template
        self.style_combo.blockSignals(True)
        self.style_combo.setCurrentIndex(
            max(0, self.style_combo.findData(template.range_style)))
        self.style_combo.blockSignals(False)
        self._fill_fields()
        self._fill_columns()

    def save_template_as(self):
        path, _ = QFileDialog.getSaveFileName(
            self, "Save the log's columns", "explog_columns.json",
            "Columns (*.json)")
        if not path:
            return
        try:
            EL.save_template(self.template, path)
        except OSError as exc:
            QMessageBox.warning(self, "Auto Exp. log",
                                f"Could not save the columns:\n{exc}")

    def load_template_from(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Load the log's columns", "", "Columns (*.json)")
        if path:
            self.set_template(EL.load_template(path))

    def reset_template(self):
        self.set_template(EL.default_template())

    # -- recording --------------------------------------------------------
    def auto_record(self) -> bool:
        """Whether files added to the list are to be logged as they come."""
        return self.auto_box.isChecked() and bool(self.log_path())

    def _keys_to_log(self):
        if self.all_radio.isChecked():
            return [key for key, _text in self._entries()]
        return list(self._selected() or ())

    def record(self, keys=None, quiet: bool = False) -> int:
        """Write the records of ``keys`` (by default the selection or the
        whole list, as chosen) to the log. Returns how many were written.
        ``quiet`` (Keep recording) reports in the status bar instead of in
        a message box."""
        def tell(text, warning=False):
            if quiet or not warning:
                self._status(text)
            if warning and not quiet:
                QMessageBox.warning(self, "Auto Exp. log", text)

        path = self.log_path()
        if not path:
            tell("Give the experiment log a file name first.", warning=True)
            return 0
        if not self.template.columns:
            tell("The experiment log has no columns: tick some fields first.",
                 warning=True)
            return 0
        keys = self._keys_to_log() if keys is None else list(keys)
        if not keys:
            tell("Select the datasets to log in the main list (or choose "
                 "Every dataset in the list).", warning=True)
            return 0
        try:
            existing = EL.read_log(path)
        except (OSError, UnicodeError) as exc:
            tell(f"Could not read the log {path}:\n{exc}", warning=True)
            return 0
        done = EL.logged_names(existing, self.template) \
            if self.skip_box.isChecked() else set()
        file_column = next((i for i, c in enumerate(self.template.columns)
                            if c.field == "@file"), None)
        rows, skipped, failed = [], 0, []
        for key in keys:
            try:
                data, record = self._load(key)
            except Exception as exc:                       # noqa: BLE001
                traceback.print_exc()
                failed.append(f"{key[0]}: {exc}")
                continue
            if data is None:
                failed.append(str(key[0]))
                continue
            try:
                cells = EL.row(self.template, data, record)
            finally:
                data.close()
            if file_column is not None and cells[file_column].strip() in done:
                skipped += 1
                continue
            if file_column is not None:
                done.add(cells[file_column].strip())
            rows.append(cells)
        written = 0
        if rows:
            try:
                written = EL.append(path, self.template, rows)
            except PermissionError as exc:
                tell(f"Could not write the log -- is it open in Excel? Close "
                     f"it there and write again.\n{exc}", warning=True)
                return 0
            except OSError as exc:
                tell(f"Could not write the log:\n{exc}", warning=True)
                return 0
        EL.save_template(self.template)
        self._update_file_status()
        message = f"Experiment log: {written} record(s) written to {path}"
        if skipped:
            message += f", {skipped} already in it"
        if failed:
            message += f", {len(failed)} could not be read"
        tell(message)
        if failed and not quiet:
            QMessageBox.warning(self, "Auto Exp. log",
                                "Could not read:\n" + "\n".join(failed[:10]))
        return written

    def closeEvent(self, event):
        EL.save_template(self.template)
        # Keep recording works with the window closed; only the reference
        # dataset is let go.
        self._close_reference()
        super().closeEvent(event)
