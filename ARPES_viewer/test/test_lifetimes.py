"""Who holds a dataset, and for how long.

Regression tests for two crashes:

* "Can't synchronously read data (identifier is not of specified type)"
  when a viewer was opened soon after another one on the same file was
  closed -- a file closed under a dataset somebody was still handed out;
* a segmentation fault some time after closing a map whose slit cut was
  open -- pyqtgraph unregistering named views during garbage collection.

The windowed tests need Qt; run headless with QT_QPA_PLATFORM=offscreen.
The segfault is tested in a child interpreter, since a crash there would
otherwise take the whole test run down with it.
"""
import gc
import os
import subprocess
import sys
import textwrap

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from loader.nxs_file import HANDLES, save_dataset          # noqa: E402
from loader.session import MemoryBudget                    # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _map_file(tmp_path, name="maps.nxs", entries=("A",)):
    rng = np.random.default_rng(0)
    cube = rng.random((6, 20, 24)).astype(np.float32) + 1
    axes = (np.arange(6.0), np.arange(20.0), np.arange(24.0))
    path = str(tmp_path / name)
    save_dataset(path, [dict(name=e, kind="map", axes=axes, labels={},
                             value=cube * (i + 1), info={})
                        for i, e in enumerate(entries)])
    return path


class Counted:
    """Stands in for NxsData: counts references, dies at zero."""

    def __init__(self):
        self.refs = 1
        self.scan = type("Scan", (), {"value": np.zeros(100)})()   # 800 bytes

    def retain(self):
        self.refs += 1
        return self

    def close(self):
        self.refs -= 1

    def alive(self):
        return self.refs > 0


# -- the memory budget owns a reference of its own ---------------------------
def test_budget_takes_and_gives_back_its_own_reference():
    budget = MemoryBudget()
    data = Counted()                      # the caller's reference
    budget.put("a", data)
    assert data.refs == 2
    data.close()                          # the caller is done (a viewer closed)
    assert budget.get("a") is data and data.alive()
    budget.discard("a")
    assert data.refs == 0


def test_budget_releases_on_eviction_and_replacement():
    budget = MemoryBudget(limit_bytes=0)
    first, second = Counted(), Counted()
    budget.put("a", first)
    budget.put("b", second)               # evicts "a"
    assert first.refs == 1 and second.refs == 2
    third = Counted()
    budget.put("b", third)                # replaces "b"
    assert second.refs == 1 and third.refs == 2
    budget.clear()
    assert third.refs == 1


def test_budget_never_hands_out_a_dead_dataset():
    budget = MemoryBudget()
    data = Counted()
    budget.put("a", data)
    data.refs = 0                         # its file closed behind the budget's back
    assert budget.get("a") is None
    assert "a" not in budget._entries


# -- one handle per file, and short looks never close it ---------------------
def test_borrow_does_not_close_a_file_in_use(tmp_path):
    path = _map_file(tmp_path)
    held = HANDLES.acquire(path)
    try:
        dset = held["A/data"] if "A/data" in held else None
        with HANDLES.borrow(path) as f:
            assert f is held
        assert held.id.valid
        if dset is not None:
            assert dset.id.valid
    finally:
        HANDLES.release(path, held)
    assert not HANDLES.is_open(path)


def test_a_release_for_a_replaced_handle_is_ignored(tmp_path):
    path = _map_file(tmp_path)
    old = HANDLES.acquire(path)
    old.close()                           # invalidated (as some builds do)
    new = HANDLES.acquire(path)           # the registry opens it afresh
    assert new is not old and new.id.valid
    HANDLES.release(path, old)            # the stale holder lets go...
    assert new.id.valid                   # ...without closing the new one
    HANDLES.release(path, new)
    assert not HANDLES.is_open(path)


def test_saving_over_a_file_open_in_a_viewer_is_refused(tmp_path):
    path = _map_file(tmp_path)
    held = HANDLES.acquire(path)
    try:
        with pytest.raises(OSError, match="open"):
            _map_file(tmp_path)           # same name
    finally:
        HANDLES.release(path, held)
    _map_file(tmp_path)                   # closed: fine


# -- NxsData -------------------------------------------------------------------
@pytest.fixture(scope="module")
def app():
    pytest.importorskip("PyQt5")
    from PyQt5.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def test_acquire_replaces_a_dataset_whose_file_died(app, tmp_path):
    from ui.widgets import NxsData
    path = _map_file(tmp_path)
    first = NxsData.acquire(path, "A")
    HANDLES._open[os.path.abspath(path)][0].close()      # the file dies
    assert not first.alive()
    second = NxsData.acquire(path, "A")
    assert second is not first and second.alive()
    np.asarray(second.angle_cube[3][:, :, 2])             # and it reads
    second.close()
    first.close()
    assert not HANDLES.is_open(path)


def test_closing_a_viewer_twice_releases_once(app, tmp_path):
    from ui.widgets import NxsData
    from ui import windows as W
    path = _map_file(tmp_path)
    launcher = NxsData.acquire(path, "A")                 # the list's reference
    window = W.ContourWindow(NxsData.acquire(path, "A"), "A", "gray", False)
    window.show()
    window.close()
    window.close()                                        # a second close event
    assert launcher._refs == 1 and launcher.alive()
    np.asarray(launcher.angle_cube[3][:, :, 1])
    launcher.close()
    assert not HANDLES.is_open(path)


def test_closing_a_map_closes_its_cut_windows_and_dialogs(app, tmp_path):
    from ui.widgets import NxsData
    from ui import windows as W
    path = _map_file(tmp_path)
    data = NxsData.acquire(path, "A")
    window = W.ContourWindow(data, "A", "gray", False)
    window.show()
    window.open_cut("slit")
    window.open_cut("deflector")
    cuts = list(window.cut_windows.values())
    dialog = window.open_degrid()
    assert data._refs == 3
    window.close()
    assert not any(c.isVisible() for c in cuts)
    assert dialog is None or not dialog.isVisible()
    assert data._refs == 0 and not HANDLES.is_open(path)


def test_a_popped_out_snapshot_outlives_its_viewer(app):
    from ui.widgets import MemoryData
    from ui import windows as W
    cut = MemoryData("cut", (np.arange(20.0), np.arange(24.0)),
                     np.random.rand(20, 24), {}, source_label="c")
    window = W.CutWindow(cut, "c", "gray", False)
    window.show()
    snapshot = window.pop_out_slice()
    assert snapshot is not None and snapshot in window.popout_windows
    window.close()
    del window
    gc.collect()
    assert snapshot.isVisible() and snapshot in W._DETACHED
    snapshot.close()
    assert snapshot not in W._DETACHED


def test_no_view_is_registered_under_a_name(app):
    import pyqtgraph as pg
    from ui.widgets import MemoryData
    from ui import windows as W
    cube = MemoryData("map", (np.arange(6.0), np.arange(20.0), np.arange(24.0)),
                      np.random.rand(6, 20, 24), {}, source_label="m")
    window = W.ContourWindow(cube, "m", "gray", False)
    window.open_cut("slit")
    assert not pg.ViewBox.NamedViews
    window.close()


def test_closing_a_map_with_its_slit_cut_open_does_not_crash():
    """The original report: open a map, open its slit cut, close the map,
    carry on. The crash came at the next garbage collection."""
    script = textwrap.dedent(f"""
        import gc, os, sys
        os.environ["QT_QPA_PLATFORM"] = "offscreen"
        sys.path.insert(0, {ROOT!r})
        import numpy as np
        from PyQt5.QtWidgets import QApplication
        app = QApplication([])
        from ui.widgets import MemoryData
        from ui import windows as W
        d = MemoryData("map", (np.arange(8.), np.arange(40.), np.arange(50.)),
                       np.random.rand(8, 40, 50), {{}}, source_label="m")
        for _ in range(3):
            w = W.ContourWindow(d, "m", "gray", False)
            w.show(); w.open_cut("slit"); w.open_cut("deflector")
            app.processEvents()
            w.close(); app.processEvents()
            del w
            gc.collect(); app.processEvents()
        print("survived")
    """)
    result = subprocess.run([sys.executable, "-c", script], capture_output=True,
                            text=True, timeout=120, cwd=ROOT)
    assert result.returncode == 0 and "survived" in result.stdout, result.stderr[-2000:]
