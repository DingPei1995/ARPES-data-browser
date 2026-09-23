"""The **Functions** menu every viewer carries (image and curve viewers).

A viewer keeps on its own top row only what is used all the time while
looking at the data -- a map's two cut buttons, the colormap, the view or
display controls -- and offers everything else here. New tools go here too
(:meth:`FunctionsMenuMixin.add_function`) unless there is a reason for them
to be always in sight.
"""
from __future__ import annotations

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QAction, QMenu, QToolButton


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
        self.functions_menu = QMenu("Functions", self)
        self.functions_menu.setToolTipsVisible(True)
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
