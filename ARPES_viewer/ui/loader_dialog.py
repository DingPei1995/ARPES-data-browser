"""
ui/loader_dialog.py
===================
The Loader window: choosing files, choosing who reads them, and saying the
two things about the axes that a reader cannot know.

Why loading has a window of its own
-----------------------------------
It used to be a button on the main panel that opened a file dialog and read
whatever came back with the one reader the program had. That is fine while
there is one beamline. It stops being fine the moment there are two, because
"which reader" and "what did it find" have nowhere to live, and because the
two things the file itself cannot say -- what order its axes are in, and
what the scanned axis physically *is* -- have to be asked before the data
goes anywhere, not corrected afterwards in five separate views.

So: pick files, see what was detected, override it if the guess is wrong,
look at the entries that came back, set the axis options, and only then
add them to the list.

What this window does not do
----------------------------
It does not read the data. It reads *structure* -- entry names, kinds,
shapes -- which is cheap even for a multi-GB scan, and hands the list back
to the main panel, which loads a dataset when one is actually opened. A
file that is listed and never opened is never read.
"""
from __future__ import annotations

import os

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (QComboBox, QDialog, QDialogButtonBox, QFileDialog,
                             QFormLayout, QGroupBox, QHBoxLayout, QLabel,
                             QLineEdit, QListWidget, QListWidgetItem,
                             QPushButton, QVBoxLayout, QWidget)

from loader import registry
from loader.registry import AXIS_ROLES, LoadOptions

#: The axis orders offered for a 3-D dataset. A free-form permutation box
#: would be more general and much easier to get wrong; these are the orders
#: anyone actually needs, named by what they do rather than by their digits.
PERMUTATIONS_3D = [
    ("Leave as the file has it", None),
    ("Swap the first two axes (1 <-> 2)", (1, 0, 2)),
    ("Swap the last two axes (2 <-> 3)", (0, 2, 1)),
    ("Swap the outer two axes (1 <-> 3)", (2, 1, 0)),
    ("Rotate forwards (3, 1, 2)", (2, 0, 1)),
    ("Rotate backwards (2, 3, 1)", (1, 2, 0)),
]

PERMUTATIONS_2D = [
    ("Leave as the file has it", None),
    ("Swap the two axes (transpose)", (1, 0)),
]


class LoaderDialog(QDialog):
    """Pick files, pick a loader, set the axis options, add them to the list.

    Modal, unlike the tools that work on a picture behind them: nothing else
    can usefully happen while this is open, and the result of it is a change
    to the main panel's list.
    """

    def __init__(self, parent=None, start_folder: str = ""):
        super().__init__(parent)
        self.setWindowTitle("Load data")
        self.resize(640, 620)
        self._paths = []
        self._entries = []

        layout = QVBoxLayout(self)
        intro = QLabel(
            "Choose the files, check that the right reader was recognised, "
            "and say what the axes are. Only structure is read here -- the "
            "measurements themselves are read when you open one.")
        intro.setWordWrap(True)
        layout.addWidget(intro)

        # -- the files ------------------------------------------------------
        files_group = QGroupBox("Files")
        files_layout = QVBoxLayout(files_group)
        row = QHBoxLayout()
        self.folder_box = QLineEdit(start_folder)
        self.folder_box.setPlaceholderText("Folder")
        self.folder_box.setReadOnly(True)
        self.pick_button = QPushButton("Select file(s)...")
        self.pick_button.clicked.connect(self.pick_files)
        row.addWidget(self.folder_box, stretch=1)
        row.addWidget(self.pick_button)
        files_layout.addLayout(row)

        self.file_list = QListWidget()
        self.file_list.setToolTip("The files that will be added to the list.")
        files_layout.addWidget(self.file_list)
        layout.addWidget(files_group)

        # -- the reader -----------------------------------------------------
        reader_group = QGroupBox("Reader")
        reader_form = QFormLayout(reader_group)
        self.loader_box = QComboBox()
        self.loader_box.addItem("Detect automatically", None)
        for loader in registry.loaders():
            self.loader_box.addItem(loader.name, loader.name)
        self.loader_box.setToolTip(
            "Which beamline's layout to read these files as. Detection is "
            "usually right, and is certain for files this program saved; "
            "override it when a file from an unfamiliar layout is being "
            "misread.")
        self.loader_box.currentIndexChanged.connect(self._refresh_entries)
        reader_form.addRow("Read as", self.loader_box)
        self.detected_label = QLabel("")
        self.detected_label.setWordWrap(True)
        self.detected_label.setMinimumHeight(34)
        reader_form.addRow("", self.detected_label)
        layout.addWidget(reader_group)

        # -- what was found ---------------------------------------------------
        found_group = QGroupBox("Datasets found")
        found_layout = QVBoxLayout(found_group)
        self.entry_list = QListWidget()
        found_layout.addWidget(self.entry_list)
        layout.addWidget(found_group)

        # -- the axes ---------------------------------------------------------
        axes_group = QGroupBox("Axes")
        axes_form = QFormLayout(axes_group)
        self.permutation_box = QComboBox()
        for label, value in PERMUTATIONS_3D:
            self.permutation_box.addItem(label, value)
        self.permutation_box.setToolTip(
            "Two beamlines can record the same measurement with the array's "
            "dimensions in a different order. Reorder them here, once, and "
            "the axis vectors and their labels come along.")
        axes_form.addRow("Axis order", self.permutation_box)

        self.role_box = QComboBox()
        self.role_box.addItem("As the file says", None)
        for key, (label, unit, _convertible) in AXIS_ROLES.items():
            self.role_box.addItem(f"{label} ({unit})" if unit else label, key)
        self.role_box.setToolTip(
            "A map is angle vs angle vs energy by default. The same scan is "
            "also how a photon-energy, temperature or gate-voltage series is "
            "recorded -- same shape, different physics. Saying so here keeps "
            "the wrong unit out of every later step, and stops a k "
            "conversion being offered for an axis that is not an angle.\n\n"
            "Some readers can work it out: a CASSIOPEE folder is a "
            "photon-energy scan exactly when the monochromator moved across "
            "it. Leave this at 'As the file says' to keep what the reader "
            "found, and set it only to correct one.")
        self.role_box.currentIndexChanged.connect(self._role_changed)
        axes_form.addRow("First axis of a map is", self.role_box)

        self.role_label_box = QLineEdit()
        self.role_label_box.setPlaceholderText("(use the default label)")
        self.role_label_box.setToolTip(
            "An axis label of your own, units included, for when none of the "
            "choices above says it.")
        axes_form.addRow("Label it", self.role_label_box)

        self.axes_note = QLabel("")
        self.axes_note.setWordWrap(True)
        self.axes_note.setMinimumHeight(34)
        axes_form.addRow("", self.axes_note)
        layout.addWidget(axes_group)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.button(QDialogButtonBox.Ok).setText("Add to list")
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        self.ok_button = buttons.button(QDialogButtonBox.Ok)
        self.ok_button.setEnabled(False)
        layout.addWidget(buttons)

        self._role_changed()

    # -- picking ---------------------------------------------------------------
    def pick_files(self):
        start = self.folder_box.text() or os.path.expanduser("~")
        if not os.path.isdir(start):
            start = os.path.expanduser("~")
        paths, _ = QFileDialog.getOpenFileNames(
            self, "Select data file(s)", start, self._file_filter())
        if not paths:
            return
        self.set_paths(paths)

    @staticmethod
    def _file_filter() -> str:
        """The file dialog's filter, built from what the loaders say they
        read rather than written out here.

        It was a fixed ``*.nxs *.h5 *.hdf5``, which is the sort of thing that
        goes stale the first time a beamline writes something else -- and did:
        a CASSIOPEE spectrum is a ``.txt``, and the dialog simply would not
        show one. Each loader also gets a filter of its own, for picking out
        one beamline's files in a folder holding several.
        """
        every = []
        per_loader = []
        for loader in registry.loaders():
            patterns = " ".join(loader.patterns)
            per_loader.append(f"{loader.name} ({patterns})")
            every.extend(loader.patterns)
        # dict.fromkeys: unique, in the order the loaders offered them.
        joined = " ".join(dict.fromkeys(every))
        return ";;".join([f"Data files ({joined})"] + per_loader
                         + ["All files (*)"])

    def set_paths(self, paths):
        """Take a list of files (also how a test drives this)."""
        self._paths = list(paths)
        self.folder_box.setText(os.path.dirname(self._paths[0]) if self._paths else "")
        self.file_list.clear()
        for path in self._paths:
            self.file_list.addItem(QListWidgetItem(os.path.basename(path)))
        self._detect()
        self._refresh_entries()

    def _detect(self):
        """Say which reader recognises the files, and select it.

        Reported per file rather than for the selection as a whole: a folder
        holding both a beamline file and a saved session is a normal thing
        to pick, and "3 files: 2 SOLEIL ANTARES, 1 saved here" is the useful
        answer.
        """
        if not self._paths:
            self.detected_label.setText("")
            return
        counts = {}
        for path in self._paths:
            loader = registry.detect(path)
            counts[loader.name if loader else "not recognised"] = \
                counts.get(loader.name if loader else "not recognised", 0) + 1
        parts = ", ".join(f"{n} × {name}" for name, n in sorted(counts.items()))
        self.detected_label.setText(f"Detected: {parts}")
        if "not recognised" in counts:
            self.detected_label.setText(
                self.detected_label.text() +
                " — choose a reader above for those, or they will be "
                "listed without a kind.")

    def _refresh_entries(self, *_):
        """List the datasets inside the chosen files."""
        self.entry_list.clear()
        self._entries = []
        options = self.options()
        for path in self._paths:
            try:
                entries = registry.list_entries(path, options)
            except Exception as exc:                     # noqa: BLE001
                entries = []
                self.entry_list.addItem(f"{os.path.basename(path)}: {exc}")
            if not entries:
                self.entry_list.addItem(
                    f"{os.path.basename(path)}: nothing recognised")
                continue
            for info in entries:
                self._entries.append((path, info))
                name = os.path.basename(path)
                entry = info.get("entry")
                label = f"{name}  ·  {entry}   [{info.get('kind')}]" if entry \
                    else f"{name}   [{info.get('kind')}]"
                self.entry_list.addItem(label)
        self.ok_button.setEnabled(bool(self._paths))

    def _role_changed(self, *_):
        role = self.role_box.currentData()
        if role is None:
            self.axes_note.setText(
                "Whatever the reader worked out is kept. For most files that "
                "is an ordinary angle-vs-angle-vs-energy map; a reader that "
                "can tell otherwise -- a CASSIOPEE folder that stepped the "
                "monochromator, say -- says so itself.")
            return
        label, unit, convertible = AXIS_ROLES.get(role, AXIS_ROLES["other"])
        if convertible:
            self.axes_note.setText(
                "Angle vs angle vs energy -- the ordinary deflector map, and "
                "the only case a k conversion applies to.")
        else:
            shown = f"{label} ({unit})" if unit else label
            self.axes_note.setText(
                f"The first axis will be labelled “{shown}”, and the "
                "map's k conversion will be refused for it, since only an "
                "emission angle converts to momentum.")

    # -- the result --------------------------------------------------------------
    def options(self) -> LoadOptions:
        return LoadOptions(
            loader=self.loader_box.currentData(),
            permutation=self.permutation_box.currentData(),
            axis0_role=self.role_box.currentData(),
            axis0_label=self.role_label_box.text().strip(),
        )

    def chosen(self):
        """``(paths, options)`` for the main panel to add."""
        return list(self._paths), self.options()
