"""
ui/common/remember.py
=====================
What a tool window's controls were set to, kept between uses.

A window that takes many settings (the step-by-step tools) loses them all
when it is closed -- by a click meant for Next, say. This keeps them: the
window's input widgets (its attributes that are spin boxes, check boxes,
radio buttons, combo boxes, line edits) are saved by attribute name in the
program's settings folder when it closes, and set again when it opens.
Signals are blocked while they are set, so nothing is recomputed half-way;
the window's ``_settings_applied()`` (if it has one) runs once after.

A window lists attributes it does not want kept in ``_not_remembered``.
"""
from __future__ import annotations

import json
import os

from PyQt5.QtWidgets import (QCheckBox, QComboBox, QDoubleSpinBox, QLineEdit,
                             QRadioButton, QSpinBox)

__all__ = ["path", "snapshot", "apply", "save", "load", "restore", "forget"]

FILE_NAME = "tool_windows.json"
_INPUTS = (QDoubleSpinBox, QSpinBox, QCheckBox, QRadioButton, QComboBox, QLineEdit)


def path() -> str:
    folder = os.environ.get("ARPES_VIEWER_CONFIG_DIR") or os.path.join(
        os.path.expanduser("~"), ".arpes_viewer")
    return os.path.join(folder, FILE_NAME)


def _inputs(owner):
    skip = set(getattr(owner, "_not_remembered", ()))
    for name, widget in sorted(vars(owner).items()):
        if name.startswith("__") or name in skip or not isinstance(widget, _INPUTS):
            continue
        yield name, widget


def snapshot(owner) -> dict:
    """The values of ``owner``'s input widgets, by attribute name."""
    values = {}
    for name, widget in _inputs(owner):
        if isinstance(widget, (QDoubleSpinBox, QSpinBox)):
            values[name] = widget.value()
        elif isinstance(widget, (QCheckBox, QRadioButton)):
            values[name] = widget.isChecked()
        elif isinstance(widget, QComboBox):
            values[name] = widget.currentText()
        elif isinstance(widget, QLineEdit):
            values[name] = widget.text()
    return values


def apply(owner, values: dict):
    """Set ``owner``'s input widgets to ``values`` (what does not fit --
    a choice no longer offered -- is left as it is)."""
    widgets = dict(_inputs(owner))
    blocked = []
    for name, value in (values or {}).items():
        widget = widgets.get(name)
        if widget is None:
            continue
        blocked.append((widget, widget.blockSignals(True)))
        try:
            if isinstance(widget, QDoubleSpinBox):
                widget.setValue(float(value))
            elif isinstance(widget, QSpinBox):
                widget.setValue(int(value))
            elif isinstance(widget, QRadioButton):
                if value:                   # the group unchecks the others
                    widget.setChecked(True)
            elif isinstance(widget, QCheckBox):
                widget.setChecked(bool(value))
            elif isinstance(widget, QComboBox):
                index = widget.findText(str(value))
                if index >= 0:
                    widget.setCurrentIndex(index)
                elif widget.isEditable():
                    widget.setEditText(str(value))
            elif isinstance(widget, QLineEdit):
                widget.setText(str(value))
        except (TypeError, ValueError):
            continue
    for widget, was in blocked:
        widget.blockSignals(was)
    hook = getattr(owner, "_settings_applied", None)
    if callable(hook):
        hook()


def _all() -> dict:
    try:
        with open(path(), encoding="utf-8") as handle:
            data = json.load(handle)
        return data if isinstance(data, dict) else {}
    except (OSError, ValueError):
        return {}


def load(key: str) -> dict:
    entry = _all().get(key)
    return entry if isinstance(entry, dict) else {}


def save(owner, key: str, extra: dict = None):
    """Keep ``owner``'s settings under ``key`` (and ``extra`` with them,
    under ``"_extra"``)."""
    data = _all()
    entry = snapshot(owner)
    if extra:
        entry["_extra"] = extra
    data[key] = entry
    try:
        os.makedirs(os.path.dirname(path()), exist_ok=True)
        with open(path(), "w", encoding="utf-8") as handle:
            json.dump(data, handle, indent=1, ensure_ascii=False)
    except OSError:
        pass


def restore(owner, key: str) -> bool:
    """Set ``owner`` as it was last kept under ``key``; False if never."""
    values = load(key)
    values = {k: v for k, v in values.items() if k != "_extra"}
    if not values:
        return False
    apply(owner, values)
    return True


def forget(key: str):
    data = _all()
    if data.pop(key, None) is not None:
        try:
            with open(path(), "w", encoding="utf-8") as handle:
                json.dump(data, handle, indent=1, ensure_ascii=False)
        except OSError:
            pass
