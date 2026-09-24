"""
arpes_batch -- ARPES viewer's algorithms, driven without the windows.

The viewer keeps everything that touches data in ``loader/`` and ``tools/``,
with no Qt in either (``test/test_imports.py`` checks it). That split is what
makes this package possible: a batch run, a command line and an agent all
call the same functions the viewer's dialogs call, and what they write is
the viewer's own saved format, so every result opens in the viewer with its
processing history attached.

What this package adds is only what the dialogs used to decide by hand:

* which files and entries there are (:mod:`.inventory`);
* the analyser's Fermi-edge curvature and E_F, from the gold reference
  rather than from points clicked on a cut (:mod:`.reference`);
* where Gamma is, as a suggestion that can be checked on a preview
  (:mod:`.center`);
* the order of the steps and their parameters, as a recipe
  (:mod:`.recipe`, :mod:`.steps`, :mod:`.runner`);
* pictures to check the result by (:mod:`.preview`).

Nothing here writes next to the raw data: every output path is checked
against the input folders first (:func:`._paths.check_output_dir`).
"""
from ._paths import ensure_viewer_on_path

ensure_viewer_on_path()

__version__ = "0.1.0"
