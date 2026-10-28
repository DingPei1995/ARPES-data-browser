"""
ui/common/stepper.py
====================
A window's work laid out as numbered steps, one page each.

The new tools that take several decisions in a row -- summing repeated
scans, the detector corrections, the tight-binding fit, the crystal viewer
-- all lead the user through them the way the guided MDC/EDC fit does
(``ui.panels.fit.FitGuide``): a title "Step 2 of 4: ...", two or three
sentences saying what the step is for and what to do, the step's own
controls, and Back / Next. Next stays greyed until the step has what the
following one needs, and its tooltip says what is missing.

Unlike the fit's guide this is the window itself, not a helper beside a
panel: every control lives on exactly one step, so a user who has never seen
the tool cannot miss one, and one who knows it pages through in a few
clicks (the step list on the left jumps straight to any step already
reached).
"""
from __future__ import annotations

from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtWidgets import (QHBoxLayout, QLabel, QListWidget, QPushButton,
                             QStackedWidget, QVBoxLayout, QWidget)


class Step:
    """One step: a title, what to do, its page, and when Next may be
    pressed (``ready`` returns ``""`` when it may, otherwise what is
    missing). ``enter`` runs every time the step is shown."""

    def __init__(self, title: str, text: str, page: QWidget, *,
                 ready=None, enter=None):
        self.title = title
        self.text = text
        self.page = page
        self.ready = ready or (lambda: "")
        self.enter = enter or (lambda: None)


class Stepper(QWidget):
    """Numbered pages with Back / Next; the last step's Next reads
    ``finish_text`` and emits :attr:`finished`."""

    #: the step now shown
    stepChanged = pyqtSignal(int)
    #: Next pressed on the last step
    finished = pyqtSignal()

    def __init__(self, steps, parent=None, *, finish_text: str = "Done",
                 show_index: bool = True):
        super().__init__(parent)
        self.steps = list(steps)
        self.finish_text = finish_text
        self.index = 0
        self.reached = 0

        outer = QHBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        self.index_list = QListWidget()
        self.index_list.setMaximumWidth(190)
        self.index_list.setToolTip("Every step; click one already reached to "
                                   "go back to it.")
        for i, step in enumerate(self.steps):
            self.index_list.addItem(f"{i + 1}. {step.title}")
        self.index_list.currentRowChanged.connect(self._jump)
        self.index_list.setVisible(show_index)
        outer.addWidget(self.index_list)

        column = QVBoxLayout()
        outer.addLayout(column, 1)
        self.title = QLabel("")
        self.title.setStyleSheet("font-size: 11pt; font-weight: bold; "
                                 "color: #1f6f8b;")
        column.addWidget(self.title)
        self.text = QLabel("")
        self.text.setWordWrap(True)
        self.text.setTextFormat(Qt.RichText)
        self.text.setStyleSheet("background: #eef4f8; border: 1px solid "
                                "#c9dbe6; border-radius: 4px; padding: 6px;")
        column.addWidget(self.text)
        self.pages = QStackedWidget()
        for step in self.steps:
            self.pages.addWidget(step.page)
        column.addWidget(self.pages, 1)

        nav = QHBoxLayout()
        self.progress = QLabel("")
        self.progress.setStyleSheet("color: #666;")
        nav.addWidget(self.progress)
        self.missing = QLabel("")
        self.missing.setStyleSheet("color: #b0602a;")
        self.missing.setWordWrap(True)
        nav.addWidget(self.missing, 1)
        self.back_button = QPushButton("← Back")
        self.back_button.clicked.connect(lambda: self.go(self.index - 1))
        nav.addWidget(self.back_button)
        self.next_button = QPushButton("Next →")
        self.next_button.setDefault(True)
        self.next_button.clicked.connect(self._next)
        nav.addWidget(self.next_button)
        column.addLayout(nav)
        self.go(0)

    # -- navigation --------------------------------------------------------
    def go(self, index: int):
        index = max(0, min(index, len(self.steps) - 1))
        self.index = index
        self.reached = max(self.reached, index)
        step = self.steps[index]
        self.pages.setCurrentIndex(index)
        self.title.setText(f"Step {index + 1} of {len(self.steps)}: {step.title}")
        self.text.setText(step.text)
        self.index_list.blockSignals(True)
        self.index_list.setCurrentRow(index)
        self.index_list.blockSignals(False)
        step.enter()
        self.refresh()
        self.stepChanged.emit(index)

    def _jump(self, row: int):
        if 0 <= row <= self.reached:
            self.go(row)
        else:
            self.index_list.blockSignals(True)
            self.index_list.setCurrentRow(self.index)
            self.index_list.blockSignals(False)

    def _next(self):
        if self.refresh():
            return
        if self.index == len(self.steps) - 1:
            self.finished.emit()
        else:
            self.go(self.index + 1)

    def refresh(self) -> str:
        """Re-ask the current step whether it is ready; returns what is
        missing ("" when Next may be pressed)."""
        missing = self.steps[self.index].ready() or ""
        last = self.index == len(self.steps) - 1
        self.next_button.setText(self.finish_text if last else "Next →")
        self.next_button.setEnabled(not missing)
        self.next_button.setToolTip(missing)
        self.missing.setText(missing)
        self.back_button.setEnabled(self.index > 0)
        self.progress.setText("● " * (self.index + 1)
                              + "○ " * (len(self.steps) - self.index - 1))
        for i in range(self.index_list.count()):
            item = self.index_list.item(i)
            flags = item.flags()
            item.setFlags(flags | Qt.ItemIsEnabled if i <= self.reached
                          else flags & ~Qt.ItemIsEnabled)
        return missing
