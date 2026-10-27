"""The three 3-D view windows (ui/panels/volume.py): opened straight away
at the balanced sampling, resampled on request, coloured per orientation
with transparency, with coordinate axes, and exported as figures."""
import numpy as np
import pytest
from PyQt5.QtWidgets import QApplication, QDialog

DEFL = np.linspace(-10, 10, 480)
SLIT = np.linspace(-8, 8, 41)
EKIN = np.linspace(54, 56, 31)
LABELS = ("Angle, deflector (deg)", "Angle, slit (deg)", "Kinetic energy (eV)")


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


def _cube():
    rng = np.random.default_rng(1)
    return rng.random((DEFL.size, SLIT.size, EKIN.size)) * 100


def _window(mode="slices"):
    from ui.panels.volume import VolumeWindow
    return VolumeWindow(_cube(), (DEFL, SLIT, EKIN), LABELS, "m", mode=mode)


@pytest.mark.parametrize("mode", ["slices", "notch", "projection"])
def test_each_view_opens_alone_at_the_balanced_sampling(app, mode):
    from ui.panels.volume import BALANCED_POINTS
    window = _window(mode)
    assert not hasattr(window, "tabs")
    assert window.values.shape[0] <= BALANCED_POINTS < DEFL.size
    assert window.values.shape[1:] == (SLIT.size, EKIN.size)
    assert "of 480 pts" in window.sampling_label.text()
    assert window.windowTitle().startswith(
        {"slices": "3D slices", "notch": "Notched cube",
         "projection": "Volume projection"}[mode])
    window.close()


def test_change_sampling_resamples_and_redraws(app, monkeypatch):
    from ui.panels import volume
    window = _window("slices")

    class Chosen:
        def __init__(self, *args, **kwargs):
            assert kwargs["targets"] == window.sampling[1]   # starts from now
        def exec_(self):
            return QDialog.Accepted
        def selection(self):
            return [(-5.0, 5.0), None, None], [30, 41, 31]
    monkeypatch.setattr(volume, "VolumeRangeDialog", Chosen)
    assert window.change_sampling()
    assert window.axes[0].min() >= -5 and window.values.shape[0] <= 60
    assert window.slice_at[0].axis.max() <= 5
    assert "cropped" in window.sampling_label.text()
    window.close()


def test_one_slice_gets_its_own_colours_and_all_changes_only_what_changed(app):
    window = _window("slices")
    styles = window.scene.styles
    window.style_target.setCurrentIndex(3)                 # energy slices only
    window.level_hi.setValue(42.0)
    window.colormap.setCurrentIndex(1)
    assert styles[2].levels[1] == 42.0 and styles[0].levels[1] != 42.0
    window.style_target.setCurrentIndex(0)                 # all
    window.gamma.setValue(2.0)
    assert all(styles[d].gamma == 2.0 for d in range(3))
    assert styles[2].levels[1] == 42.0                     # kept its own
    window.opacity.setValue(50.0)
    assert not styles[0].opaque()
    window.auto_levels_each()
    assert all(styles[d].levels[1] > styles[d].levels[0] for d in range(3))
    window.close()


def test_transparent_and_faded_faces(app):
    from tools.output import colormaps
    from ui.panels.volume import _colour_face
    image = np.array([[0.0, 5.0], [10.0, np.nan]])
    lut = colormaps.get_lut("gray")
    qimage = _colour_face(image, lut, (0.0, 10.0), opacity=0.5, fade=0.5)
    alpha = [qimage.pixelColor(x, y).alpha() for y in range(2) for x in range(2)]
    # (u, v) -> pixel (v, u): 0 -> transparent, 5 and 10 -> half opaque
    assert alpha[0] == 0 and alpha[3] == 0
    assert alpha[1] == alpha[2] == 128


def test_axes_placement_by_hand(app):
    window = _window("notch")
    assert window.scene.axes.edges is None                 # automatic
    assert not window.axes_edge[0].isEnabled()
    window.axes_edge[0].setCurrentIndex(3)
    window._axes_chosen()
    assert window.scene.axes.edges[0] == (1, 1)
    assert window.axes_auto.currentData() is False
    window.axes_on.setChecked(False)
    assert window.scene.axes is None
    window.close()


def test_the_view_exports_as_a_figure_and_joins_a_composer(app, tmp_path):
    window = _window("notch")
    window.resize(900, 700)
    for name in ("v.png", "v.pdf", "v.svg"):
        assert window.save_view(str(tmp_path / name))
        assert (tmp_path / name).stat().st_size > 1000
    composer = window.to_figure_composer()
    assert composer.figure.panels[0].data.kind == "scene3d"
    window.to_figure_composer(target=composer)
    assert len(composer.figure.panels) == 2
    composer.close()
    window.close()


@pytest.fixture(autouse=True)
def _no_modal_boxes(monkeypatch):
    """A warning box would wait for a click that never comes."""
    from PyQt5.QtWidgets import QMessageBox

    def fail(parent, title, text, *args, **kwargs):
        raise AssertionError(f"{title}: {text}")
    for name in ("warning", "information", "critical"):
        monkeypatch.setattr(QMessageBox, name, staticmethod(fail))
