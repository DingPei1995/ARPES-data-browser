"""Shared test setup.

Cyclic garbage collection runs only between tests, never in the middle of
one -- the same rule the program applies with ui.gcguard (collection from
the event loop only, never inside a Qt call). With automatic collection the
window tests crash the interpreter at random, exactly as the program did.
"""
import gc

import pytest


@pytest.fixture(autouse=True)
def _collect_between_tests():
    was_enabled = gc.isenabled()
    gc.disable()
    yield
    gc.collect()
    if was_enabled:
        gc.enable()
