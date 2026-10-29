"""
ui/common/screen.py
===================
Windows that fit the screen they open on.

A dialog cannot be made smaller than its contents need, so one whose
settings, plots and notes add up to more than the screen's height opens with
its bottom -- often its buttons -- off screen, and cannot be shrunk to fit.
:func:`fit_to_screen` (or the :func:`fits_screen` class decorator) puts the
contents in a scroll area, keeps a trailing row of buttons in sight below it,
and opens the window no larger than the screen.
"""
from __future__ import annotations

from PyQt5.QtCore import QEvent, QObject, Qt
from PyQt5.QtWidgets import (QAbstractButton, QApplication, QDialogButtonBox,
                             QFrame, QMainWindow, QPushButton, QScrollArea,
                             QVBoxLayout, QWidget)


#: Off for the guide's screenshots (devtools/guide_screenshots.py), whose
#: virtual screen is smaller than the windows are meant to be shown at.
LIMIT_TO_SCREEN = True


def _is_button_widget(widget) -> bool:
    """A button box: Qt's, or this program's ``ButtonRow`` (by name, so this
    module needs nothing from the viewers)."""
    return (isinstance(widget, QDialogButtonBox)
            or type(widget).__name__ == "ButtonRow")


def _is_button_row(item) -> bool:
    """Whether a layout item is a row of buttons (to keep in sight)."""
    widget = item.widget()
    if widget is not None:
        return _is_button_widget(widget)
    layout = item.layout()
    if layout is None or layout.count() == 0:
        return False
    widgets = [layout.itemAt(i).widget() for i in range(layout.count())]
    widgets = [w for w in widgets if w is not None]
    return bool(widgets) and all(
        isinstance(w, QAbstractButton) or _is_button_widget(w)
        for w in widgets) and any(isinstance(w, QPushButton) for w in widgets)


def fit_to_screen(window, fraction: float = 0.88):
    """Let a tall window be resized smaller than its contents, and open it
    no taller than the screen.

    A dialog cannot be shrunk below what its contents need, so one whose
    settings, plots and notes add up to more than the screen's height opens
    with its bottom (often its buttons) off screen and cannot be made to fit.
    Here the contents go into a scroll area: the window can be any height,
    and what does not fit scrolls. A row of buttons at the bottom stays
    outside it, always in sight. The window then opens at the size it asked
    for, or the screen's, whichever is smaller. Call it at the end of the
    window's construction (after any ``resize`` of its own).
    """
    if getattr(window, "_fits_screen", False):
        return window
    window._fits_screen = True
    wanted = window.size() if window.testAttribute(Qt.WA_Resized) else None
    area = QScrollArea()
    area.setWidgetResizable(True)
    area.setFrameShape(QFrame.NoFrame)
    if isinstance(window, QMainWindow):
        inner = window.takeCentralWidget()
        if inner is None:
            return window
        area.setWidget(inner)
        window.setCentralWidget(area)
    else:
        layout = window.layout()
        if layout is None:
            return window
        trailing = []
        while layout.count() and _is_button_row(layout.itemAt(layout.count() - 1)):
            trailing.insert(0, layout.takeAt(layout.count() - 1))
        inner = QWidget()
        inner.setLayout(layout)
        area.setWidget(inner)
        outer = QVBoxLayout(window)
        margins = layout.contentsMargins()
        outer.setContentsMargins(0, 0, 0, margins.bottom())
        outer.addWidget(area, 1)
        for item in trailing:
            row = QWidget()
            holder = QVBoxLayout(row)
            holder.setContentsMargins(margins.left(), 4, margins.right(), 0)
            if item.widget() is not None:
                holder.addWidget(item.widget())
            elif item.layout() is not None:
                holder.addLayout(item.layout())
            outer.addWidget(row)
    hint = inner.sizeHint()
    screen = (window.screen() if hasattr(window, "screen") and window.screen()
              else QApplication.primaryScreen())
    available = screen.availableGeometry() if screen is not None else None
    width = max(wanted.width() if wanted else 0,
                hint.width() + area.verticalScrollBar().sizeHint().width() + 4)
    height = (wanted.height() if wanted else hint.height() + 60)
    if available is not None and LIMIT_TO_SCREEN:
        width = min(width, int(available.width() * 0.95))
        height = min(height, int(available.height() * fraction))
    window.resize(width, height)
    return window


def fits_screen(cls):
    """Class decorator: :func:`fit_to_screen` on every new window of
    ``cls``, once its own ``__init__`` has built it."""
    build = cls.__init__

    def __init__(self, *args, **kwargs):
        build(self, *args, **kwargs)
        fit_to_screen(self)

    __init__.__doc__ = build.__doc__
    cls.__init__ = __init__
    return cls




class _WheelScrollsWindow(QObject):
    """See :func:`wheel_scrolls_window`."""

    def eventFilter(self, watched, event):
        if event.type() != QEvent.Wheel or event.modifiers() & Qt.ControlModifier:
            return False
        widget = watched.parentWidget() if watched is not None else None
        while widget is not None and not isinstance(widget, QScrollArea):
            widget = widget.parentWidget()
        if widget is None:
            return False
        bar = widget.verticalScrollBar()
        if bar.maximum() <= bar.minimum():
            return False            # nothing to scroll: the plot zooms
        QApplication.sendEvent(bar, event)
        return True


def wheel_scrolls_window(*views):
    """Make the mouse wheel over these plots scroll the window they are in
    (when it scrolls at all) rather than zoom the plot; Ctrl + wheel still
    zooms.

    A pyqtgraph plot takes every wheel turn as a zoom, so in a window too
    tall for the screen, a turn of the wheel with the pointer over the
    picture -- most of the window -- zoomed the picture instead of bringing
    the settings below it into sight. ``views`` are QGraphicsViews
    (``PlotWidget``, ``ImageView.ui.graphicsView``); the filter goes on
    their viewports, which are what receive the wheel.
    """
    for view in views:
        viewport = view.viewport() if hasattr(view, "viewport") else view
        guard = _WheelScrollsWindow(viewport)
        viewport.installEventFilter(guard)
        viewport._wheel_guard = guard
