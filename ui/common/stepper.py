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

:func:`guard_window` keeps the work from being lost by a slip: closing the
window (its last step's Close, which sits where Next was, the title bar's
close button, or Esc) asks first, Enter never presses the last step's
button, and what was set is kept and comes back the next time the window
opens (with **Restore defaults** to start afresh).
"""
from __future__ import annotations

import os

from PyQt5.QtCore import QEvent, QObject, Qt, pyqtSignal
from PyQt5.QtWidgets import (QHBoxLayout, QLabel, QListWidget, QMessageBox,
                             QPushButton, QStackedWidget, QVBoxLayout, QWidget)

__all__ = ["Step", "Stepper", "guard_window", "confirm_close"]


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
        #: steps left out (they do not apply to what the user is doing)
        self.skipped = set()

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
        self.defaults_button = QPushButton("Restore defaults")
        self.defaults_button.setAutoDefault(False)
        self.defaults_button.setToolTip("Set every control back to how the "
                                        "window first opens (what was kept "
                                        "from last time is forgotten).")
        self.defaults_button.setVisible(False)
        nav.addWidget(self.defaults_button)
        self.progress = QLabel("")
        self.progress.setStyleSheet("color: #666;")
        nav.addWidget(self.progress)
        self.missing = QLabel("")
        self.missing.setStyleSheet("color: #b0602a;")
        self.missing.setWordWrap(True)
        nav.addWidget(self.missing, 1)
        self.back_button = QPushButton("← Back")
        self.back_button.setAutoDefault(False)
        self.back_button.clicked.connect(lambda: self._back())
        nav.addWidget(self.back_button)
        self.next_button = QPushButton("Next →")
        self.next_button.setDefault(True)
        self.next_button.clicked.connect(self._next)
        nav.addWidget(self.next_button)
        column.addLayout(nav)
        self.go(0)

    # -- navigation --------------------------------------------------------
    def _shown(self):
        return [i for i in range(len(self.steps)) if i not in self.skipped]

    def _numbered(self, index: int) -> str:
        shown = self._shown()
        place = shown.index(index) + 1 if index in shown else index + 1
        return f"Step {place} of {len(shown)}: {self.steps[index].title}"

    def set_skipped(self, index: int, skipped: bool = True):
        """Leave step ``index`` out (or put it back): Next and Back pass it
        by, and it leaves the list and the count."""
        if skipped:
            self.skipped.add(index)
        else:
            self.skipped.discard(index)
        shown = self._shown()
        for i in range(self.index_list.count()):
            self.index_list.item(i).setHidden(i in self.skipped)
            if i in shown:
                self.index_list.item(i).setText(f"{shown.index(i) + 1}. {self.steps[i].title}")
        if self.index in self.skipped:
            self.go(self.index)
        else:
            self.title.setText(self._numbered(self.index))
            self.refresh()

    def go(self, index: int):
        index = max(0, min(index, len(self.steps) - 1))
        if index in self.skipped:
            later = [i for i in self._shown() if i > index]
            earlier = [i for i in self._shown() if i < index]
            index = later[0] if later else earlier[-1]
        self.index = index
        self.reached = max(self.reached, index)
        step = self.steps[index]
        self.pages.setCurrentIndex(index)
        self.title.setText(self._numbered(index))
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
        later = [i for i in self._shown() if i > self.index]
        if not later:
            self.finished.emit()
        else:
            self.go(later[0])

    def _back(self):
        earlier = [i for i in self._shown() if i < self.index]
        if earlier:
            self.go(earlier[-1])

    def refresh(self) -> str:
        """Re-ask the current step whether it is ready; returns what is
        missing ("" when Next may be pressed)."""
        missing = self.steps[self.index].ready() or ""
        shown = self._shown()
        last = self.index == shown[-1]
        self.next_button.setText(self.finish_text if last else "Next →")
        # Enter never presses the last step's button (Close sits where Next
        # was): only a click does
        self.next_button.setDefault(not last)
        self.next_button.setAutoDefault(not last)
        self.next_button.setEnabled(not missing)
        self.next_button.setToolTip(missing)
        self.missing.setText(missing)
        self.back_button.setEnabled(self.index > shown[0])
        place = shown.index(self.index) + 1 if self.index in shown else 1
        self.progress.setText("● " * place + "○ " * (len(shown) - place))
        for i in range(self.index_list.count()):
            item = self.index_list.item(i)
            flags = item.flags()
            item.setFlags(flags | Qt.ItemIsEnabled if i <= self.reached
                          else flags & ~Qt.ItemIsEnabled)
        return missing

    def retitle(self, index: int, title: str, text: str = None):
        """Change step ``index``'s title (and what it says)."""
        step = self.steps[index]
        step.title = title
        if text is not None:
            step.text = text
        shown = self._shown()
        if index in shown:
            self.index_list.item(index).setText(f"{shown.index(index) + 1}. {title}")
        if index == self.index:
            self.title.setText(self._numbered(index))
            self.text.setText(step.text)


def confirm_close(parent, title: str) -> bool:
    """Ask whether to close ``title``; "Keep working" is the default, so a
    second hasty Enter keeps it open. (``ARPES_VIEWER_NO_CONFIRM`` set:
    no question, for scripts.)"""
    if os.environ.get("ARPES_VIEWER_NO_CONFIRM"):
        return True
    box = QMessageBox(QMessageBox.Question, title, f"Close “{title}”?",
                      parent=parent)
    box.setInformativeText("What you have set is kept: it comes back the next "
                           "time you open this window.")
    close = box.addButton("Close", QMessageBox.AcceptRole)
    keep = box.addButton("Keep working", QMessageBox.RejectRole)
    box.setDefaultButton(keep)
    box.setEscapeButton(keep)
    box.exec_()
    return box.clickedButton() is close


class _CloseGuard(QObject):
    """Asks before ``window`` closes, and keeps its settings when it does."""

    def __init__(self, window, key, title, extra):
        super().__init__(window)
        self.window = window
        self.key = key
        self.title = title
        self.extra = extra
        self.ask = True

    def save(self):
        from ui.common import remember

        remember.save(self.window, self.key,
                      self.extra() if callable(self.extra) else None)

    def eventFilter(self, obj, event):
        if obj is not self.window:
            return False
        if event.type() == QEvent.KeyPress and event.key() == Qt.Key_Escape:
            self.window.close()
            return True
        if event.type() == QEvent.Close and self.window.isVisible():
            if self.ask and not confirm_close(self.window, self.title):
                event.ignore()
                return True
            self.save()
        return False


def guard_window(window, stepper: Stepper, key: str, title: str, *,
                 close_on_finish: bool = False, extra=None):
    """Guard ``window`` (whose steps are ``stepper``) against closing by a
    slip, and keep its settings under ``key`` between uses. With
    ``close_on_finish`` the last step's button closes it (asking first).
    ``extra()`` returns more to keep (a dict); the window reads it back
    with ``remember.load(key).get("_extra")``. Returns True when settings
    from last time were put back."""
    from ui.common import enterkeys, remember

    enterkeys.install()
    defaults = remember.snapshot(window)

    def reset():
        remember.apply(window, defaults)
        remember.forget(key)
        stepper.go(stepper.index)

    restored = remember.restore(window, key)
    stepper.defaults_button.setVisible(True)
    stepper.defaults_button.clicked.connect(reset)
    if restored:
        stepper.defaults_button.setToolTip(
            "These are your settings from last time. Click to set every "
            "control back to how the window first opens.")
    guard = _CloseGuard(window, key, title, extra)
    window.installEventFilter(guard)
    window._close_guard = guard
    if close_on_finish:
        stepper.finished.connect(window.close)
    return restored
