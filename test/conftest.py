"""Shared test setup.

Cyclic garbage collection runs only between tests, never in the middle of
one -- the same rule the program applies with ui.common.gcguard (collection from
the event loop only, never inside a Qt call). With automatic collection the
window tests crash the interpreter at random, exactly as the program did.
"""
import gc
import sys

import pytest


@pytest.fixture(autouse=True)
def _collect_between_tests():
    was_enabled = gc.isenabled()
    gc.disable()
    yield
    gc.collect()
    if was_enabled:
        gc.enable()


@pytest.fixture(autouse=True)
def _private_config_dir(tmp_path, monkeypatch):
    """Settings a dialog remembers between uses (the Brillouin-zone tool's
    last crystal) go to a folder of this test's own, not the user's home --
    so no test sees what another one, or a real session, left behind."""
    monkeypatch.setenv("ARPES_VIEWER_CONFIG_DIR", str(tmp_path / "config"))


#: The QApplication, once a test has made one, held for the whole session.
_QAPP = []


@pytest.fixture(autouse=True)
def _keep_the_qapplication():
    """Never let the QApplication be destroyed between test modules.

    Each GUI module has its own ``app`` fixture, which drops its reference
    when the module ends; if that was the last one, PyQt deletes the
    QApplication and the next module makes a new one. Any widget still
    alive at that point -- a dialog held by a ``QTimer.singleShot`` that
    never ran, for want of an event loop -- outlives the application it was
    made under, and a later module crashes with a segmentation fault in
    whatever Qt call comes next. The program itself has one QApplication
    for its whole life; so, now, does the suite. Qt is not imported here
    unless a test already did.
    """
    yield
    widgets = sys.modules.get("PyQt5.QtWidgets")
    if widgets is not None and not _QAPP:
        app = widgets.QApplication.instance()
        if app is not None:
            _QAPP.append(app)
