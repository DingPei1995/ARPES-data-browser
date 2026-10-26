"""Make ``python -m pytest`` work from this folder.

The tests live in ``test/`` and import the packages beside it (``from
tools.bz3d import ...``). pytest puts the *test file's* directory on
sys.path, not the project root, so without this every test would fail on
its first import. Putting this file at the root is the documented way to
say "the project root is here".

It also pins Qt to the offscreen platform when nothing else has chosen one,
so the suite runs over ssh and on a build machine with no display -- Qt
aborts the whole interpreter rather than raising when it cannot open one,
which would take the non-GUI tests down with it.
"""
import os
import sys

ROOT = os.path.dirname(os.path.abspath(__file__))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
