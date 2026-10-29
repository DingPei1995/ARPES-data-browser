"""The version written by hand in the documents matches ``version.py``.

The program reads its version from ``version.py``; the README and the user
guides state it in their first lines as plain text, so they fall behind
when only ``version.py`` is changed.
"""
import datetime
import os
import re

from arpes_viewer import version
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _head(*parts, lines=8):
    with open(os.path.join(ROOT, *parts), encoding="utf-8") as f:
        return "".join(f.readline() for _ in range(lines))


def test_the_readme_states_the_current_version_and_date():
    head = _head("README.md")
    assert f"**Version {version.VERSION}**" in head
    day = datetime.date.fromisoformat(version.UPDATED)
    assert f"{day.day} {day:%B %Y}" in head


def test_the_user_guides_state_the_current_version():
    for lang in ("en", "zh", "fr"):
        head = _head("docs", f"USER_GUIDE_{lang}.md")
        assert f"v{version.VERSION}" in head, lang


def test_the_changelog_opens_on_the_current_version():
    with open(os.path.join(ROOT, "CHANGELOG.md"), encoding="utf-8") as f:
        first = next(line for line in f if line.startswith("## "))
    assert re.match(rf"## {re.escape(version.VERSION)} — {version.UPDATED}", first)
