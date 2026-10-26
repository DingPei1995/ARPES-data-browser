"""The program's version and the date it was last updated.

Shown in the launcher's title bar and footer, so that a user (and a bug
report) can say which build is in use. Bump both with every change that
reaches users, and give the CHANGELOG entry the same version.
"""
VERSION = "1.0.7"
UPDATED = "2026-10-26"
REPOSITORY = "https://github.com/DingPei1995/ARPES-data-browser"


def title() -> str:
    """``"ARPES viewer v1.0.0 (updated 2026-10-25)"``."""
    return f"ARPES viewer v{VERSION} (updated {UPDATED})"
