"""
ui/common/enterkeys.py
======================
Enter in a field means "take this value", never "press a button".

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
from PyQt5.QtWidgets import (QAbstractSpinBox, QApplication, QColorDialog,
                             QComboBox, QDialog, QFileDialog, QFontDialog,
                             QInputDialog, QLineEdit, QMessageBox,
                             QProgressDialog)

__all__ = ["install"]

_QT_DIALOGS = (QFileDialog, QMessageBox, QInputDialog, QColorDialog,
               QFontDialog, QProgressDialog)
_FIELDS = (QAbstractSpinBox, QLineEdit, QComboBox)


class _EnterTakesTheValue(QObject):
    def eventFilter(self, obj, event):
        if (event.type() == QEvent.KeyPress
                and event.key() in (Qt.Key_Return, Qt.Key_Enter)
                and isinstance(obj, QDialog) and not isinstance(obj, _QT_DIALOGS)):
            focus = QApplication.focusWidget()
            field = focus
            # the line edit inside a spin box or an editable combo box
            while field is not None and not isinstance(field, _FIELDS):
                field = field.parentWidget()
                if field is obj:
                    field = None
            if field is not None and field.window() is obj:
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
