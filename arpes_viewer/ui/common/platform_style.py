"""
ui/common/platform_style.py
===========================
What differs between macOS, Windows and Linux in how the windows look, kept
in one place.

* **The widget style.** macOS's native ("macintosh") style draws push
  buttons as fixed-height Aqua bezels with wide outer margins, ignores the
  padding the style sheets give them, and draws no text on a progress bar.
  In the launcher's tool panel -- a column of buttons two pixels apart --
  the bezels overlapped and the headings slid under them. Fusion lays a
  button out from the style sheet like the Windows style does, so on macOS
  the application uses Fusion and the panels look the same on every system.
* **The interface font.** "Lucida Sans Unicode" is a Windows font; asked for
  on a Mac, Qt spends a moment searching for aliases at start-up and then
  falls back to whatever it finds. The first installed family of a short
  list is used instead.
* **The modifier key.** Qt reads the Command key as ``Qt.ControlModifier``
  on macOS, so "Ctrl" in a hint is "⌘" there.
"""
import sys

IS_MACOS = sys.platform == "darwin"

#: The key Qt calls Control, as the user knows it.
MOD_KEY = "⌘" if IS_MACOS else "Ctrl"

#: The interface font, best first: the one the panels were designed with
#: (Windows), its macOS sibling, then each system's own.
UI_FONTS = ("Lucida Sans Unicode", "Lucida Grande", "Segoe UI",
            "Helvetica Neue", "DejaVu Sans")


def apply_style(app):
    """Give ``app`` the same widget style on every system: Fusion on macOS,
    where the native style does not fit the compact panels; the system's own
    elsewhere."""
    from PyQt5.QtCore import Qt
    from PyQt5.QtWidgets import QApplication
    # Sharp icons and pixmaps on Retina and high-DPI screens.
    QApplication.setAttribute(Qt.AA_UseHighDpiPixmaps, True)
    if IS_MACOS:
        app.setStyle("Fusion")


def ui_font_family():
    """The first of ``UI_FONTS`` installed here; Qt's default family if none
    is. Needs a ``QApplication``."""
    from PyQt5.QtGui import QFontDatabase
    from PyQt5.QtWidgets import QApplication
    installed = set(QFontDatabase().families())
    for family in UI_FONTS:
        if family in installed:
            return family
    return QApplication.font().family()
