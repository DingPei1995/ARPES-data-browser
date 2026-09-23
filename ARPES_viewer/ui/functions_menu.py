"""The **Functions** menu every viewer carries (image and curve viewers).

A viewer keeps on its own top row only what is used all the time while
looking at the data -- a map's two cut buttons, the colormap, the view or
display controls -- and offers everything else here. New tools go here too
(:meth:`FunctionsMenuMixin.add_function`) unless there is a reason for them
to be always in sight.
"""
from __future__ import annotations

from PyQt5.QtCore import QEvent, QObject, Qt
from PyQt5.QtWidgets import QAction, QMenu, QToolButton, QToolTip


class _InstantTooltips(QObject):
    """Event filter behind :func:`instant_tooltips`."""

    def __init__(self, menu: QMenu):
        super().__init__(menu)
        self.menu = menu
        self.current = None

    def show_for(self, action):
        if action is self.current:
            return
        self.current = action
        text = action.toolTip() if action is not None else ""
        # Qt fills an empty tooltip with the entry's own text; that says
        # nothing the entry does not, so it is not shown.
        # (Qt's version of it also drops a trailing "...".)
        own = action.text().replace("&", "") if action is not None else ""
        if (action is None or action.isSeparator() or not text
                or text in (own, own.rstrip(".").rstrip())):
            QToolTip.hideText()
            return
        rect = self.menu.actionGeometry(action)
        QToolTip.showText(self.menu.mapToGlobal(rect.topRight()), text,
                          self.menu, rect)

    def eventFilter(self, obj, event):
        kind = event.type()
        if kind == QEvent.MouseMove:
            self.show_for(self.menu.actionAt(event.pos()))
        elif kind in (QEvent.Leave, QEvent.Hide):
            self.current = None
            QToolTip.hideText()
        elif kind == QEvent.ToolTip:
            return True                     # Qt's delayed tooltip: ours instead
        return False


def instant_tooltips(menu: QMenu) -> QMenu:
    """Show a menu entry's tooltip the moment the pointer reaches it.

    Qt's own menu tooltips wait for the pointer to rest (about 0.7 s) and,
    moving from one entry to the next, leave the previous entry's text up
    until it has rested again -- which reads as the menu lagging. Here the
    tooltip follows the pointer instead: it changes as soon as another entry
    is under it, disappears over an entry that has none, and sits beside the
    entry rather than under the pointer.

    Driven by the menu's mouse moves (and by ``hovered``, for the keyboard),
    not by ``hovered`` alone: whether a *disabled* entry is ever "hovered"
    depends on the platform's style, and a disabled entry is the one whose
    tooltip -- why it is disabled -- matters most.
    """
    menu.setToolTipsVisible(True)
    tips = _InstantTooltips(menu)
    menu.installEventFilter(tips)
    menu.hovered.connect(tips.show_for)
    menu._instant_tooltips = tips
    return menu


class FunctionsMenuMixin:
    """Mix into a ``QMainWindow`` subclass (listed first)."""

    #: Sections, in the order they are shown. A function added under any
    #: other name gets a section of its own, after these and before "Slice".
    FUNCTION_SECTIONS = ("Analysis", "Data operations", "Visualization")

    def add_function(self, text: str, slot, tooltip: str = "", *,
                     section: str = "Data operations", visible: bool = True,
                     enabled: bool = True) -> QAction:
        """Offer ``slot`` in this window's **Functions** menu.

        Returns the action: hide it (``setVisible``) where it does not apply,
        or disable it with a tooltip saying why, which the menu shows.
        """
        action = QAction(text, self)
        if tooltip:
            action.setToolTip(tooltip)
            action.setStatusTip(tooltip.split("\n")[0])
        action.triggered.connect(lambda *_: slot())
        action.setVisible(visible)
        action.setEnabled(enabled)
        self.__dict__.setdefault("_functions", []).append((section, action))
        return action

    def function_actions(self, section: str = None):
        """The menu's actions (in one section, or all), visible or not."""
        return [action for name, action in getattr(self, "_functions", [])
                if section is None or name == section]

    def make_functions_button(self) -> QToolButton:
        """The **Functions** button with its menu, for the window's top row."""
        self.functions_menu = instant_tooltips(QMenu("Functions", self))
        self.functions_menu.aboutToShow.connect(self._rebuild_functions_menu)
        self.functions_button = QToolButton()
        self.functions_button.setText("Functions")
        self.functions_button.setToolTip(
            "Everything this window can do: analysis, operations and figures.")
        self.functions_button.setPopupMode(QToolButton.InstantPopup)
        self.functions_button.setToolButtonStyle(Qt.ToolButtonTextOnly)
        self.functions_button.setMenu(self.functions_menu)
        self._rebuild_functions_menu()
        return self.functions_button

    def _rebuild_functions_menu(self):
        """Lay the menu out afresh each time it opens: sections in a fixed
        order, and none left showing a heading over nothing, since which
        tools apply depends on the data (a kz map, a spin EDC...)."""
        menu = getattr(self, "functions_menu", None)
        if menu is None:
            return
        menu.clear()
        entries = getattr(self, "_functions", [])
        named = [section for section, _ in entries]
        order = [sec for sec in self.FUNCTION_SECTIONS if sec in named]
        order += [sec for sec in dict.fromkeys(named)
                  if sec not in order and sec != "Slice"]
        order += ["Slice"] if "Slice" in named else []
        for section in order:
            actions = [a for sec, a in entries if sec == section and a.isVisible()]
            if not actions:
                continue
            menu.addSection(section)
            for action in actions:
                menu.addAction(action)
