"""
ui/common/crashguard.py
=======================
An error is written down and said; it does not close the program.

PyQt5 aborts the whole program when a Python exception escapes a slot and
no ``sys.excepthook`` of the program's own is set -- a bug in one button
took every open window, and any unsaved work, with it. :func:`install`
sets one: the traceback goes to the error stream, which
:mod:`tools.runtime.errorlog` copies into ``error_logs/`` in the program's
folder, and a message says what happened and where the details are (once
every few seconds at most, so a fault repeating in a redraw does not bury
the screen in boxes). The program then carries on.
"""
from __future__ import annotations

import sys
import threading
import time

from PyQt5.QtCore import QTimer
from PyQt5.QtWidgets import QApplication, QMessageBox

from tools.runtime import errorlog

__all__ = ["install"]

_QUIET_SECONDS = 5.0
_last = [0.0]


def _say(summary: str, path: str):
    app = QApplication.instance()
    if app is None or time.monotonic() - _last[0] < _QUIET_SECONDS:
        return
    _last[0] = time.monotonic()
    where = (f"The details were saved in:\n{path}\n\nWith a bug report, "
             "upload that file (on GitHub Issues: drag it into the text box) "
             "or paste its text -- the main panel's foot has “Copy the newest "
             "error log”." if path else "")

    def show():
        box = QMessageBox(QMessageBox.Warning, "ARPES viewer -- an error",
                          f"Something went wrong:\n\n{summary}\n\nThe program "
                          "keeps running; if in doubt, save your work.")
        box.setInformativeText(where)
        box.exec_()

    # not from inside the failing call (a paint event, say): from the loop
    QTimer.singleShot(0, show)


def _hook(exc_type, exc, tb, where=""):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc, tb)
        return
    path = errorlog.report(exc_type, exc, tb, where)
    _say(f"{exc_type.__name__}: {exc}", path)


def install(program_folder: str, version: str = ""):
    """Keep the errors in ``program_folder``/error_logs, and keep the
    program running through them."""
    errorlog.install(program_folder, version)
    sys.excepthook = _hook
    threading.excepthook = lambda args: _hook(
        args.exc_type, args.exc_value, args.exc_traceback,
        f"thread {getattr(args.thread, 'name', '')}")
