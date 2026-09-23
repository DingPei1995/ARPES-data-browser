"""One readout cursor shared by a map's contour and its two cut windows.

A map is a cube over three axes -- deflector angle, slit angle, energy --
and each of its windows is a plane through it:

=====================  ============  ============  ======================
window                 horizontal    vertical      position along the rest
=====================  ============  ============  ======================
constant-E contour     deflector     slit          energy slider
slit cut               slit          energy        "integrate over" slider
deflector cut          deflector     energy        "integrate over" slider
=====================  ============  ============  ======================

So one point in the cube is a cursor in all three windows at once, and this
keeps them that way: moving the cursor in any window, or any window's slider,
moves the cursor in the others and re-slices their images through the new
point. Switching the cursor on or off in one switches it in all.

Only the *point* is shared. Integration widths are each window's own and
are not passed between windows: the contour's energy "+/-", each cut's
"integrate over" "+/-", and the cuts' EDC and MDC "+/-" are all set
independently, and each is drawn only where it is set -- a cut shades its
own EDC and MDC windows around its cursor; the contour's cursor shows no
width at all.

Updates to the *other* windows are coalesced on a short timer: dragging a
cursor on a large map would otherwise re-read a slice of the file for every
pixel the mouse crosses. How long it waits adapts to what a re-slice costs:
an in-memory map follows the mouse live, while one read lazily from a large
file (~0.1 s per slice for a 113 x 800 x 983 ANTARES map) is re-sliced when
the mouse pauses, so the cursor being dragged never stutters.
"""
from __future__ import annotations

import time

from PyQt5.QtCore import QObject, QTimer

#: How long to gather cursor moves before re-slicing the other windows, ms.
COALESCE_MS = 30
#: A re-slice slower than this (ms) switches to waiting for the mouse to rest...
SLOW_FLUSH_MS = 60
#: ...for this long (ms).
REST_MS = 120


class _Member:
    def __init__(self, window, panel, x_axis, y_axis, control, control_axis):
        self.window = window
        self.panel = panel
        self.view = panel.view
        self.x_axis, self.y_axis = x_axis, y_axis
        self.control, self.control_axis = control, control_axis
        self.slots = []            # (signal, slot) pairs, to disconnect on removal

    def connect(self, signal, slot):
        signal.connect(slot)
        self.slots.append((signal, slot))

    def disconnect_all(self):
        for signal, slot in self.slots:
            try:
                signal.disconnect(slot)
            except (TypeError, RuntimeError):
                pass
        self.slots.clear()


class CursorLink(QObject):
    """See the module docstring. Owned by the contour window."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.members = []
        #: axis name -> index of the shared point along it
        self.index = {}
        self.on = False
        self._busy = False          # applying our own changes: ignore echoes
        self._pending = set()       # members whose change is being passed on
        self._last_flush_ms = 0.0
        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(COALESCE_MS)
        self._timer.timeout.connect(self.flush)

    # -- membership -----------------------------------------------------------
    def add(self, window, panel, *, x_axis: str, y_axis: str, control=None,
            control_axis: str = None):
        """Link one window. ``control`` is its slider along the third axis
        (``control_axis``), whose "+/-" is the window this view is summed
        over -- shown in the other windows."""
        member = _Member(window, panel, x_axis, y_axis, control, control_axis)
        self.members.append(member)
        view = member.view
        member.connect(view.readoutToggled, lambda on, m=member: self._toggled(m, on))
        member.connect(view.readoutMoved, lambda ix, iy, m=member: self._moved(m, ix, iy))
        if control is not None:
            member.connect(control.changed, lambda m=member: self._slice_changed(m))
        closed = getattr(window, "closed", None)
        if closed is not None:
            member.connect(closed, lambda *_, m=member: self.remove(m))
        if self.on:
            self._busy = True
            try:
                view.set_readout_cursor_visible(True)
                self._apply(member)
            finally:
                self._busy = False
        return member

    def remove(self, member):
        if member not in self.members:
            return
        member.disconnect_all()
        self.members.remove(member)
        self._pending.discard(member)

    def member_for(self, window):
        for member in self.members:
            if member.window is window:
                return member
        return None

    # -- cursor on / off ------------------------------------------------------
    def _toggled(self, member, on: bool):
        if self._busy:
            return
        self._busy = True
        try:
            if on:
                if not self.on or not self._complete():
                    self._adopt_from(member)
                self.on = True
                for other in self.members:
                    if other is not member and other.view.readout_cursor is None:
                        other.view.set_readout_cursor_visible(True)
                for other in self.members:
                    self._apply(other)
            else:
                self.on = False
                self._pending.clear()
                for other in self.members:
                    if other is not member and other.view.readout_cursor is not None:
                        other.view.set_readout_cursor_visible(False)
        finally:
            self._busy = False

    def _complete(self) -> bool:
        needed = set()
        for member in self.members:
            needed.update((member.x_axis, member.y_axis))
            if member.control_axis:
                needed.add(member.control_axis)
        return needed.issubset(self.index)

    def _adopt_from(self, member):
        """Take the shared point from the window the cursor was switched on
        in: where its cursor is, and where its slider is."""
        indices = member.view.readout_indices()
        if indices is not None:
            self.index[member.x_axis], self.index[member.y_axis] = indices
        if member.control is not None:
            self.index[member.control_axis] = member.control.index
        # Axes this window does not show come from the windows that do.
        for other in self.members:
            if other.control is not None and other.control_axis not in self.index:
                self.index[other.control_axis] = other.control.index
            other_indices = other.view.readout_indices()
            for axis, value in zip((other.x_axis, other.y_axis),
                                   other_indices or (None, None)):
                if axis not in self.index and value is not None:
                    self.index[axis] = value

    # -- moves --------------------------------------------------------------
    def _moved(self, member, ix: int, iy: int):
        if self._busy or not self.on:
            return
        self.index[member.x_axis] = int(ix)
        self.index[member.y_axis] = int(iy)
        self._schedule(member)

    def _slice_changed(self, member):
        if self._busy or not self.on:
            return
        self.index[member.control_axis] = member.control.index
        self._schedule(member)

    def _schedule(self, member):
        self._pending.add(member)
        if self._last_flush_ms > SLOW_FLUSH_MS:
            self._timer.start(REST_MS)          # restarted by every move
        elif not self._timer.isActive():
            self._timer.start(COALESCE_MS)

    def flush(self):
        """Pass the latest point on to every window. Called by the timer;
        tests call it directly.

        The window the change came from is included: applying is a no-op for
        a window already at the point (a cursor is only moved if it is on a
        different data point, so one being dragged is never snapped), and a
        slider whose integration window cannot reach the point -- it keeps
        its window inside the data -- is where the others learn that."""
        self._timer.stop()
        if not self._pending:
            return
        self._pending = set()
        self._busy = True
        started = time.perf_counter()
        try:
            for member in list(self.members):
                self._apply(member)
        finally:
            self._busy = False
        self._last_flush_ms = 1000.0 * (time.perf_counter() - started)

    def _apply(self, member):
        """Put ``member`` at the shared point: its slider first (which
        re-slices its image), then its cursor."""
        if member.control is not None and member.control_axis in self.index:
            wanted = self.index[member.control_axis]
            if member.control.index != wanted:
                member.control.set_index(wanted)
        if member.view.readout_cursor is not None:
            member.view.set_readout_index(self.index.get(member.x_axis),
                                          self.index.get(member.y_axis))
