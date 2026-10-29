"""
tools/runtime/remembered.py
===========================
Settings remembered between uses of a tool, per size of data.

A beamtime measures the same sample over and over -- one map per photon
energy, per polarisation -- and every one of those maps has the same number
of points and wants the same conversion: the same Γ, the same rotation, the
same output grid. Retyping them for each is the chore this saves. The key is
the data's shape, because that is what "the same measurement again" looks
like from here: a map of another size is another measurement, and settings
carried over to it by name alone would be wrong without anyone noticing.

Kept in ``remembered.json`` beside the program's config (moved by
``ARPES_VIEWER_CONFIG_DIR``, as the tests do), so they survive a restart.
Qt-free.
"""
from __future__ import annotations

import json
import os

FILE_NAME = "remembered.json"

#: how many sizes each tool keeps, newest first
KEEP = 20


def _path() -> str:
    folder = os.environ.get("ARPES_VIEWER_CONFIG_DIR") or os.path.join(
        os.path.expanduser("~"), ".arpes_viewer")
    return os.path.join(folder, FILE_NAME)


def shape_key(shape) -> str:
    """``(97, 241, 166)`` -> ``"97x241x166"``."""
    return "x".join(str(int(n)) for n in tuple(shape))


def _load() -> dict:
    try:
        with open(_path(), encoding="utf-8") as handle:
            stored = json.load(handle)
        return stored if isinstance(stored, dict) else {}
    except (OSError, ValueError):
        return {}


def recall(tool: str, shape):
    """The settings ``tool`` was last used with on data of this ``shape``,
    or None."""
    entry = (_load().get(tool) or {}).get(shape_key(shape))
    return dict(entry) if isinstance(entry, dict) else None


def remember(tool: str, shape, settings: dict) -> None:
    """Keep ``settings`` (plain JSON values) for the next use of ``tool``
    on data of this ``shape``. A failure to write is not an error: this is
    a convenience, and the conversion it follows has already happened."""
    stored = _load()
    by_shape = stored.get(tool) if isinstance(stored.get(tool), dict) else {}
    key = shape_key(shape)
    by_shape.pop(key, None)
    by_shape = {key: dict(settings), **by_shape}
    stored[tool] = dict(list(by_shape.items())[:KEEP])
    path = _path()
    try:
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w", encoding="utf-8") as handle:
            json.dump(stored, handle, indent=1)
    except OSError:
        pass
