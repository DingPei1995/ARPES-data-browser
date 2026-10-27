"""Plot tools: guide lines and labels drawn on a live viewer
(ui/viewers/plot_tools.py), kept apart from the figure composer's Figure tools."""
import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


@pytest.fixture(scope="module")
def app():
    pytest.importorskip("PyQt5")
    from PyQt5.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def _cut_window():
    from ui.viewers import windows as W
    from ui.viewers.widgets import MemoryData
    data = MemoryData("cut", (np.linspace(-10, 10, 41), np.linspace(-1, 0.2, 61)),
                      np.random.default_rng(0).random((41, 61)),
                      {"x": "Angle (deg)", "y": "E - E_F (eV)"},
                      source_label="c")
    return W.CutWindow(data, "c", "gray", False)


def test_viewer_plot_tools_draw_edit_and_remove_guides(app):
    import pyqtgraph as pg
    from ui.viewers.plot_tools import open_plot_tools, guides_of
    window = _cut_window()
    assert window.plot_tools_button.text() == "Plot tools"
    dialog = open_plot_tools(window)
    assert open_plot_tools(window) is dialog            # one per window
    plot = dialog.plot()

    dialog.kind_combo.setCurrentIndex(dialog.kind_combo.findData("horizontal"))
    dialog.y_box.setValue(0.0)
    dialog.set_colour("#ff9500")
    guide = dialog.add()
    assert isinstance(guide.item, pg.InfiniteLine) and guide.item.angle == 0
    assert guide.item.getViewBox() is plot.vb
    # editing the selected guide redraws it
    dialog.y_box.setValue(-0.3)
    assert guide.y == pytest.approx(-0.3)
    assert guide.item.pos().y() == pytest.approx(-0.3)
    # dragging on the plot updates the guide
    guide.item.setPos(pg.Point(0, -0.5))
    assert guide.y == pytest.approx(-0.5)

    dialog.list.setCurrentRow(-1)
    dialog.kind_combo.setCurrentIndex(dialog.kind_combo.findData("angled"))
    dialog.angle_box.setValue(30)
    angled = dialog.add()
    assert angled.item.angle == pytest.approx(30)
    dialog.list.setCurrentRow(-1)
    dialog.kind_combo.setCurrentIndex(dialog.kind_combo.findData("text"))
    dialog.text_edit.setText("Γ")
    label = dialog.add()
    assert isinstance(label.item, pg.TextItem)
    assert len(guides_of(window.view)) == 3

    dialog.list.setCurrentRow(0)
    dialog.remove_selected()
    assert guide.item is None and len(dialog.current_guides()) == 2
    dialog.clear()
    assert not guides_of(window.view)
    dialog.close()
    window.close()


def test_guides_are_drawn_on_the_exported_image(app):
    from ui.viewers.plot_tools import Guide, draw_guides
    blank = np.zeros((50, 100, 4), dtype=np.uint8)
    blank[..., 3] = 255
    out = draw_guides(blank, [Guide("vertical", x=0.5, color="#ff0000",
                                    style="solid", width=2)],
                      x_range=(0, 1), y_range=(0, 1))
    assert out.shape == blank.shape
    assert out[25, 50, 0] > 200 and out[25, 5, 0] == 0


def test_frozen_slice_names_figure_composer_apart_from_plot_tools(app):
    from ui.viewers.widgets import FrameWindow
    window = FrameWindow("t", np.ones((5, 6)), np.arange(5.0), np.arange(6.0))
    assert window.plot_tools_button.text() == "Plot tools"
    assert window.figure_button.text() == "Figure composer..."
    figure = window.open_figure()
    assert figure.tools_button.text() == "Figure tools"
    assert figure.open_tools().windowTitle() == "Figure tools"
    assert figure.open_tools().isVisible()      # beside the page, not a dialog
    figure.close()
    window.close()
