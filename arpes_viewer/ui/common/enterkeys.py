"""
ui/common/enterkeys.py
======================
Enter in a field means "take this value" -- and redraw what it shows --
never "press a button".

The rule, for every window of the program, now and for any added later
(the filter is on the application, so a new dialog has it without asking):

* Enter never presses Next, Close, OK, Save or any other button that moves
  on, closes, saves or makes something. Only a button that itself has the
  focus answers Enter (Space too).
* A field updates what depends on it when its value is entered: connect
  its ``valueChanged`` (with ``setKeyboardTracking(False)``) or
  ``editingFinished`` to the redraw.
* A window that redraws only on an explicit action (the Brillouin-zone
  tool's Plot) names it ``enter_redraws``; Enter then does that, and only
  that.

A Qt dialog presses its default button when Enter reaches it -- and Enter
typed in a spin box or a line edit does reach it: the field takes the value
and passes the key on. In a step-by-step window the default button is Next
(or, on the last step, Close), so typing an energy and pressing Enter moved
on a step instead of redrawing. This filter stops the key at the window once
the field has had it: the value is taken (the field's ``editingFinished`` /
``valueChanged`` still fire, so the window updates), and no button is
pressed. A button that itself has the focus still answers Enter, and Qt's
own small dialogs (file, message, input, colour, progress) keep their usual
Enter = OK.
"""
from __future__ import annotations

from PyQt5.QtCore import QEvent, QObject, Qt
from PyQt5.QtWidgets import (QApplication, QColorDialog, QDialog, QFileDialog,
                             QFontDialog, QInputDialog, QMessageBox,
                             QProgressDialog)

__all__ = ["install"]

_QT_DIALOGS = (QFileDialog, QMessageBox, QInputDialog, QColorDialog,
               QFontDialog, QProgressDialog)


class _EnterTakesTheValue(QObject):
    """Stops Return/Enter at every dialog of the program.

    The key reaches the window only when nothing in it wanted it: a field
    (spin box, line edit, combo box) has already taken its value and passed
    it on, a check box or a list ignores it, while a button with the focus
    clicks itself and stops it there. So stopping it at the window loses
    nothing -- and keeps the window's default button (Next, Close, OK) from
    being pressed by an Enter meant for a field. It does not look at where
    the focus is: the value just entered may have started work that opened
    a progress window and moved the focus away before the key gets here,
    which is how Enter in a field once pressed Next all the same."""

    def eventFilter(self, obj, event):
        if event.type() != QEvent.KeyPress or event.key() not in (Qt.Key_Return,
                                                                    Qt.Key_Enter):
            return False
        # on its way up from the field: the nearest widget that names its
        # redraw does it (a viewer's selection box, a zoom window's range)
        redraw = getattr(obj, "enter_redraws", None)
        if callable(redraw) and not isinstance(obj, _QT_DIALOGS):
            event.accept()
            redraw()
            return True
        if isinstance(obj, QDialog) and not isinstance(obj, _QT_DIALOGS):
            event.accept()
            return True
        return False


def install(app=None):
    """Install the filter on the application, once."""
    app = app or QApplication.instance()
    if app is None or getattr(app, "_enter_takes_the_value", None) is not None:
        return
    app._enter_takes_the_value = _EnterTakesTheValue(app)
    app.installEventFilter(app._enter_takes_the_value)
