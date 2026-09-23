"""Garbage collection only at safe moments, from the event loop.

Why this exists
---------------
Python's cyclic garbage collector runs whenever enough objects have been
allocated -- which can be in the middle of a Qt call: while Qt is deleting a
window, walking a widget's children, or building a pyqtgraph plot, as soon
as some Python code runs inside it and allocates. If that collection frees
the last Python reference to *other* Qt objects, their C++ side is deleted
underneath the C++ code that is still running, and the program crashes later
at an unrelated place (a segmentation fault, no traceback): here while a
window was being closed or opened, after several windows had been opened and
closed in a row.

The standard remedy (pyqtgraph ships the same thing as
``pyqtgraph.GarbageCollector``) is to switch automatic collection off and
collect from a timer instead, so a collection only ever runs from the event
loop, between Qt calls, never inside one. Reference counting -- which frees
almost everything -- is unaffected; only reference *cycles* wait up to one
timer interval to be freed.

Measured on this program's test suite, which opens and closes a few hundred
windows: with automatic collection 6 of 8 full runs crashed (and the round-38
package 5 of 5); collecting only between tests, none did (see CHANGELOG,
round 39, for the counts). The crash surfaced in whatever Qt call came next
-- walking a closing window's children, or pyqtgraph building a plot --
which is why it looked random.
"""
from __future__ import annotations

import gc

from PyQt5.QtCore import QObject, QTimer

#: How often to look, in ms. pyqtgraph uses 1 s.
INTERVAL_MS = 1000


class SafeGarbageCollector(QObject):
    """Turn automatic collection off; collect from a timer instead.

    Each tick collects exactly the generations automatic collection would
    have collected by now, using the interpreter's own thresholds, so memory
    behaves as before -- only *when* it happens changes.
    """

    def __init__(self, parent=None, interval_ms: int = INTERVAL_MS):
        super().__init__(parent)
        self.thresholds = gc.get_threshold()
        gc.disable()
        self.timer = QTimer(self)
        self.timer.timeout.connect(self.check)
        self.timer.start(int(interval_ms))

    def check(self):
        t0, t1, t2 = self.thresholds
        c0, c1, c2 = gc.get_count()
        if c2 > t2:
            gc.collect(2)
        elif c1 > t1:
            gc.collect(1)
        elif c0 > t0:
            gc.collect(0)

    def stop(self):
        """Back to automatic collection (tests, or on the way out)."""
        self.timer.stop()
        gc.enable()


_guard = None


def install(app) -> SafeGarbageCollector:
    """Install once, owned by the application object."""
    global _guard
    if _guard is None:
        _guard = SafeGarbageCollector(app)
    return _guard
