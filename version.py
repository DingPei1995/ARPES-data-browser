"""The program's version and the date it was last updated.

Shown in the launcher's title bar and footer, so that a user (and a bug
report) can say which build is in use. The version number changes only when
the maintainer decides it does -- not with every change: later work is
recorded under the current version's CHANGELOG entry. A small change moves
only the last number (1.1.0 -> 1.1.1); new tools move the middle one
(1.1.1 -> 1.2.0). Keep UPDATED current, and the
version in ``ARPES viewer.app/Contents/Info.plist`` (the macOS launcher)
the same as VERSION -- a test checks it.
"""
VERSION = "1.2.0"
UPDATED = "2026-10-28"
REPOSITORY = "https://github.com/DingPei1995/ARPES-data-browser"


def title() -> str:
    """``"ARPES viewer v1.2.0 (updated 2026-10-28)"``."""
    return f"ARPES viewer v{VERSION} (updated {UPDATED})"
