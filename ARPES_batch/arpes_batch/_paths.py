"""
Where the viewer's code is, and where output may go.

The viewer is imported as it stands (``loader``, ``tools``), from the
``ARPES_viewer`` folder beside this one, or from ``$ARPES_VIEWER_ROOT``.
Nothing in it is edited or copied.
"""
from __future__ import annotations

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_VIEWER_ROOT = os.path.normpath(os.path.join(HERE, "..", "..", "ARPES_viewer"))


def viewer_root() -> str:
    return os.path.abspath(os.environ.get("ARPES_VIEWER_ROOT", DEFAULT_VIEWER_ROOT))


def ensure_viewer_on_path() -> str:
    root = viewer_root()
    if not os.path.isfile(os.path.join(root, "loader", "registry.py")):
        raise ImportError(
            f"ARPES viewer not found at {root}; set ARPES_VIEWER_ROOT to the "
            f"folder that holds loader/ and tools/")
    if root not in sys.path:
        sys.path.insert(0, root)
    return root


def _real(path: str) -> str:
    return os.path.normcase(os.path.realpath(os.path.abspath(path)))


def _inside(child: str, parent: str) -> bool:
    child, parent = _real(child), _real(parent)
    return child == parent or child.startswith(parent.rstrip(os.sep) + os.sep)


class UnsafeOutput(ValueError):
    """The output folder is, or contains, or sits inside, a raw-data folder."""


def check_output_dir(output: str, inputs) -> str:
    """Refuse an output folder that would put anything next to the raw data.

    Three cases are refused: the output *is* an input folder, lies *inside*
    one, or *contains* one (a later clean-up of the output would then reach
    the raw data). Symlinks are resolved first, so a link into the raw
    folder is caught too. Returns the absolute output path.
    """
    output = os.path.abspath(output)
    for raw in inputs or ():
        raw = os.path.dirname(raw) if os.path.isfile(raw) else raw
        if _inside(output, raw):
            raise UnsafeOutput(
                f"output {output} is inside the raw-data folder {raw}; "
                f"choose a folder outside it")
        if _inside(raw, output):
            raise UnsafeOutput(
                f"output {output} contains the raw-data folder {raw}; "
                f"choose a folder that does not")
    return output
