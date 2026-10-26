"""
ui/notify.py
============
The message boxes that stand in for the lab's ``Win10Notif`` toasts.

These live here rather than beside the rest of the ``tools_packages``
stand-ins in ``tools/system.py`` for one reason: they need Qt, and nothing
in ``tools/`` may. That rule is what makes ``tools/`` and ``loader/``
importable from a script on a machine with no display -- which is how they
are tested -- so the two functions that break it move rather than the rule
bending. ``test/test_imports.py`` checks it.
"""
from __future__ import annotations

from PyQt5.QtWidgets import QMessageBox


def info(message: str, title: str = "Info"):
    QMessageBox.information(None, title, message)


def warning(message: str, title: str = "Warning"):
    QMessageBox.warning(None, title, message)
