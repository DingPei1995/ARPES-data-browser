"""The right-click menus, entry by entry, on the windows the guides show.

* Every image panel offers the readout cursor, the selection box, the line
  profile, the stack plot and the three exports; the SPEM window adds the
  two-way box integration, the snapshot window and Swap X and Y.
* SPEM: a box on the spatial map gives the E-k spectrum summed over those
  pixels, and a box on the spectrum gives the map summed over that k-E
  window -- swapped axes or not.
* Swap X and Y carries the selection box with the pixels it framed (it used
  to stay where it was on screen, over another part of the sample).
* The spectrum panel's title says which pixel, or which box, it shows.
* After an E-k box has been integrated onto the map, the map's title says
  so, and "Spatial map over the whole spectrum" brings the overview back.
* Saving from the list is written to the operations log.

Needs Qt; run headless with QT_QPA_PLATFORM=offscreen.
"""
import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PyQt5")
from PyQt5.QtWidgets import QApplication                      # noqa: E402

from devtools import demo_data                                # noqa: E402
from ui.widgets import NxsData                                # noqa: E402
from ui import windows as W                                   # noqa: E402


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


@pytest.fixture(scope="module")
def files(tmp_path_factory):
    paths = demo_data.write_all(str(tmp_path_factory.mktemp("demo")))
    return {os.path.basename(p): p for p in paths}


def pump(app, n=5):
    for _ in range(n):
        app.processEvents()


def _open(app, path, entry=None):
    window = W.open_viewer(NxsData(path, entry=entry), os.path.basename(path),
                           "gray", False)
    window.show()
    pump(app)
    return window


def _visible_entries(view):
    return [a.text() for a in view.view_box.menu.actions()
            if a.isVisible() and a.text()]


COMMON = ["Readout cursor", "Selection box", "Line profile...", "Export"]


def _check_common_entries(app, view):
    for entry in COMMON:
        assert entry in _visible_entries(view)
    view.action_readout.trigger()
    assert view.readout_value() is not None
    view.action_readout.trigger()
    assert view.readout_cursor is None
    view.action_selection.trigger()
    assert view.selection_corners() is not None
    view._sync_stack_action()
    assert view.action_stack.text().startswith("Stack plot of the selection (")
    view.action_stack.trigger()
    pump(app)
    assert view._stack_window.isVisible()
    view._stack_window.close()
    view.action_profile.trigger()
    pump(app)
    profile = view._profile_window
    distance, intensity = profile.profile()
    assert distance.size == profile.samples_box.value()
    assert np.isfinite(intensity).any()
    profile.close()
    for page, action in enumerate((view.action_export_nxs, view.action_export_image,
                                   view.action_export_eps)):
        action.trigger()
        pump(app)
        assert view._export_dialog.tabs.currentIndex() == page
        view._export_dialog.close()
    view.action_selection.trigger()
    assert view.selection_corners() is None


@pytest.mark.parametrize("name, entry", [("demo_cuts.nxs", "Gamma_cut_LH"),
                                         ("demo_fermi_surface_map.nxs", None),
                                         ("demo_kz_scan.nxs", None)])
def test_every_viewer_offers_the_common_entries(app, files, name, entry):
    window = _open(app, files[name], entry)
    view = window.view
    assert not view.action_apply.isVisible()        # SPEM only
    _check_common_entries(app, view)
    window.close()


@pytest.fixture
def spem(app, files):
    window = _open(app, files["demo_spatial_scan.nxs"])
    yield window
    window.close()


def test_spem_panels_offer_their_own_entries(app, spem):
    spatial = _visible_entries(spem.spatial_view)
    assert "Integrate selection -> E vs k panel" in spatial
    assert "Open this position in a new window" in spatial
    assert "Swap X and Y" in spatial
    assert "Spatial map over the whole spectrum" not in spatial
    assert "Integrate selection -> spatial panel" in _visible_entries(spem.frame_view)
    for view in (spem.spatial_view, spem.frame_view):
        _check_common_entries(app, view)


def test_a_box_on_the_map_sums_the_spectrum_over_its_pixels(app, spem):
    view, cube = spem.spatial_view, spem.data.scan.value4d
    assert not view.action_apply.isEnabled()        # no box yet
    view.set_selection_corners(10, 10, 30, 20)
    assert view.action_apply.isEnabled()
    (ix0, ix1), (iy0, iy1) = view.selection_index_ranges()
    view.action_apply.trigger()
    np.testing.assert_allclose(spem.frame_view.last_frame,
                               cube[iy0:iy1 + 1, ix0:ix1 + 1].sum(axis=(0, 1)))
    assert "integrated" in spem.frame_label
    assert spem.frame_panel.title_label.text().startswith("E vs k, summed over x 10 to 30")
    # typed corners + Refresh do the same
    spem.spatial_panel.corner_spins[0].setValue(0.0)
    spem.spatial_panel.refresh_button.click()
    (ix0, ix1), (iy0, iy1) = view.selection_index_ranges()
    assert ix0 == 0
    np.testing.assert_allclose(spem.frame_view.last_frame,
                               cube[iy0:iy1 + 1, ix0:ix1 + 1].sum(axis=(0, 1)))


def test_a_box_on_the_spectrum_maps_where_it_comes_from_and_back(app, spem):
    frame, cube = spem.frame_view, spem.data.scan.value4d
    frame.set_selection_corners(-5, -0.4, 5, -0.1)
    (k0, k1), (e0, e1) = frame.selection_index_ranges()
    frame.action_apply.trigger()
    expected = cube[:, :, k0:k1 + 1, e0:e1 + 1].sum(axis=(2, 3))
    np.testing.assert_allclose(spem.spatial_view._image2d, expected)
    title = spem.spatial_panel.title_label.text()
    assert title.startswith("Spatial map, integrated over Angle, slit")
    assert spem.action_full_map.isVisible()
    spem.action_full_map.trigger()
    np.testing.assert_allclose(spem.spatial_view._image2d, spem.data.spatial_overview)
    assert spem.spatial_panel.title_label.text() == "Spatial overview"
    assert not spem.action_full_map.isVisible()


def test_swapping_the_axes_keeps_the_box_on_the_same_pixels(app, spem):
    view, cube = spem.spatial_view, spem.data.scan.value4d
    view.set_selection_corners(10, 10, 30, 20)
    before = view.selection_index_ranges()
    view.action_swap.trigger()
    assert spem.swap_xy
    x0, y0, x1, y1 = view.selection_corners()
    assert (x0, x1) == pytest.approx((10, 20), abs=1.01)   # was y
    assert (y0, y1) == pytest.approx((10, 30), abs=1.01)   # was x
    assert view.selection_index_ranges() == (before[1], before[0])
    # and integrating it still sums the same measured pixels
    view.action_apply.trigger()
    (iy0, iy1), (ix0, ix1) = view.selection_index_ranges()
    np.testing.assert_allclose(spem.frame_view.last_frame,
                               cube[iy0:iy1 + 1, ix0:ix1 + 1].sum(axis=(0, 1)))
    # the k-E box lands on the swapped map transposed
    spem.frame_view.set_selection_corners(-5, -0.4, 5, -0.1)
    (k0, k1), (e0, e1) = spem.frame_view.selection_index_ranges()
    spem.frame_view.action_apply.trigger()
    np.testing.assert_allclose(spem.spatial_view._image2d,
                               cube[:, :, k0:k1 + 1, e0:e1 + 1].sum(axis=(2, 3)).T)


def test_a_snapshot_window_holds_the_current_spectrum(app, spem):
    spem.spatial_view.action_popout.trigger()
    pump(app)
    assert len(spem.popout_windows) == 1
    snapshot = spem.popout_windows[0]
    np.testing.assert_allclose(snapshot.view.last_frame, spem.frame_view.last_frame)
    snapshot.close()


def test_saving_from_the_list_is_logged(tmp_path, monkeypatch, app, files):
    import ARPES_viewer as av
    from loader.session import SessionStore, read_log
    store = SessionStore(root=str(tmp_path / "sessions"))
    monkeypatch.setattr(av, "session", store)
    target = str(tmp_path / "mine.nxs")
    key = (files["demo_cuts.nxs"], "Gamma_cut_LH")
    monkeypatch.setattr(av, "selected_keys", lambda: [key])
    monkeypatch.setattr(av, "loaded_items", {key: {
        "name": "Gamma_cut_LH", "kind": "Cut", "entry": "Gamma_cut_LH",
        "path": key[0], "memory": None, "backing": None}})
    monkeypatch.setattr(av, "_rebind_saved_rows", lambda *a: None)
    monkeypatch.setattr(av, "QFileDialog", type("D", (), {
        "getSaveFileName": staticmethod(lambda *a, **k: (target, ""))}))

    class _Status:
        def showMessage(self, *a):
            pass
    monkeypatch.setattr(av, "ui", type("U", (), {
        "FilePathLE": type("L", (), {"displayText": lambda self: ""})(),
        "statusbar": _Status()})(), raising=False)
    monkeypatch.setattr(av, "win", None, raising=False)

    def no_box(*args, **kwargs):
        raise AssertionError(f"unexpected message box: {args[1:3]}")
    monkeypatch.setattr(av, "QMessageBox", type("B", (), {
        "warning": staticmethod(no_box), "information": staticmethod(no_box)}))
    av.save_selected()
    assert os.path.exists(target)
    rows = read_log(store.root)
    assert any(row[2].strip() == "SAVE" and target in row[4] for row in rows)
