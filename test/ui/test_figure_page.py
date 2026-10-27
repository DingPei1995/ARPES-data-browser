"""The figure composer as a manuscript page: a blank page of slots,
panels dragged and resized, insets, drawings, undo, and clean exports."""
import sys

import numpy as np
import pytest

from PyQt5.QtCore import QEvent, QPoint, QPointF, Qt
from PyQt5.QtGui import QImage, QMouseEvent
from PyQt5.QtWidgets import QApplication

import tools.output.figure as nf


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication(sys.argv)


def _send(widget, kind, point):
    event_type = {"press": QEvent.MouseButtonPress, "move": QEvent.MouseMove,
                  "release": QEvent.MouseButtonRelease}[kind]
    button = Qt.NoButton if kind == "move" else Qt.LeftButton
    buttons = Qt.NoButton if kind == "release" else Qt.LeftButton
    QApplication.sendEvent(widget, QMouseEvent(event_type, QPointF(point),
                                               button, buttons, Qt.NoModifier))


def _drag(canvas, start_mm, end_mm):
    a = canvas.to_px(*start_mm).toPoint()
    b = canvas.to_px(*end_mm).toPoint()
    _send(canvas, "press", a)
    _send(canvas, "move", QPoint(a.x() + 6, a.y() + 6))
    _send(canvas, "move", b)
    _send(canvas, "release", b)


def _panel():
    from ui.panels.figure import panel_from_arrays
    k = np.linspace(-1, 1, 60)
    energy = np.linspace(-0.5, 0.1, 50)
    kk, ee = np.meshgrid(k, energy, indexing="ij")
    return panel_from_arrays(np.exp(-((ee + 0.4 * kk ** 2) / 0.03) ** 2), k, energy,
                             "k (Å^-1)", "E - E_F (eV)", name="cut")


@pytest.fixture
def window(app):
    from ui.panels.figure import FigureWindow
    window = FigureWindow([], "blank")
    window.resize(1300, 820)
    window.show()
    app.processEvents()
    window.canvas._ensure_image()
    yield window
    window.close()


def test_a_blank_page_has_four_lettered_slots(window):
    figure = window.figure
    assert figure.free and len(figure.panels) == 4
    assert all(p.data.kind == "empty" for p in figure.panels)
    assert [p.label for p in figure.panels] == ["(a)", "(b)", "(c)", "(d)"]


def test_slots_are_filled_and_never_exported_empty(window, tmp_path):
    figure = window.figure
    window.fill_panel(figure.panels[0], _panel())
    assert figure.panels[0].data.kind == "cut" and figure.panels[0].label == "(a)"
    path = str(tmp_path / "f.png")
    window.write(path, "png")
    image = QImage(path)
    # the empty slot (b) leaves the page white where it sits
    x, y, w, h = figure.panels[1].frame
    scale = image.width() / figure.style.width_mm
    pixel = image.pixelColor(int((x + w / 2) * scale), int((y + h / 2) * scale))
    assert (pixel.red(), pixel.green(), pixel.blue()) == (255, 255, 255)


def test_a_panel_is_dragged_snapped_and_undone(window):
    canvas, panel = window.canvas, window.figure.panels[1]
    x, y, w, h = panel.frame
    _drag(canvas, (x + w / 2, y + h / 2), (x + w / 2 + 7.2, y + h / 2 + 3.1))
    nx, ny, nw, nh = panel.frame
    assert (nw, nh) == (w, h) and nx > x and ny > y
    window.undo()
    assert window.figure.panels[1].frame == (x, y, w, h)


def test_a_handle_resizes_the_panel(window):
    canvas, panel = window.canvas, window.figure.panels[0]
    canvas.select("panel", 0)
    x, y, w, h = panel.frame
    _drag(canvas, (x + w, y + h), (x + w - 6, y + h - 4))
    assert panel.frame[2] < w and panel.frame[3] < h
    assert panel.frame[:2] == (x, y)


def test_a_zoom_inset_magnifies_a_region_and_follows_its_panel(window):
    figure, canvas = window.figure, window.canvas
    window.fill_panel(figure.panels[0], _panel())
    window.set_tool("zoom")
    x, y, w, h = figure.panels[0].frame
    _drag(canvas, (x + 0.2 * w, y + 0.2 * h), (x + 0.5 * w, y + 0.6 * h))
    inset = figure.panels[-1]
    assert inset.role == "inset" and inset.zoom_of == figure.panels[0].uid
    assert inset.label == "" and inset.x_range[0] < inset.x_range[1]
    assert canvas.tool == "select"
    old = inset.frame
    figure.move_panel(figure.panels[0], (x + 10, y, w, h))
    assert inset.frame[0] == pytest.approx(old[0] + 10)


def test_drawings_are_drawn_edited_and_deleted(window):
    canvas = window.canvas
    for tool in ("arrow", "rect", "line", "ellipse"):
        window.set_tool(tool)
        _drag(canvas, (10, 10), (30, 20))
    kinds = [s.kind for s in window.figure.shapes]
    assert kinds == ["arrow", "rect", "line", "ellipse"]
    assert canvas.selected_shape() is window.figure.shapes[-1]
    assert window.shape_editor.isVisible()
    window.shape_editor.color.set_color("#ff0000")
    window.shape_editor.color.colorChanged.emit("#ff0000")
    assert window.figure.shapes[-1].color == "#ff0000"
    canvas.delete_selection()
    assert len(window.figure.shapes) == 3


def test_a_slot_drawn_inside_a_panel_is_an_inset(window):
    figure = window.figure
    window.set_tool("slot")
    x, y, w, h = figure.panels[0].frame
    _drag(window.canvas, (x + 2, y + 2), (x + w / 2, y + h / 2))
    assert figure.panels[-1].role == "inset" and figure.panels[-1].label == ""
    window.set_tool("slot")
    _drag(window.canvas, (2, figure.style.height_mm - 1), (20, figure.style.height_mm - 15))
    assert figure.panels[-1].role == "panel"
    assert sorted(p.label for p in figure.main_panels()) == \
        ["(a)", "(b)", "(c)", "(d)", "(e)"]


def test_a_picture_fills_a_slot_with_its_shape(window, tmp_path):
    image = QImage(40, 20, QImage.Format_RGB32)
    image.fill(Qt.red)
    path = str(tmp_path / "sketch.png")
    image.save(path)
    panel = window.load_picture(window.figure.panels[2], path)
    assert panel.data.kind == "picture"
    assert panel.frame[3] == pytest.approx(panel.frame[2] / 2)


def test_old_callers_still_get_a_grid(app):
    from ui.panels.figure import FigureWindow
    window = FigureWindow([_panel(), _panel(), _panel()], "grid")
    assert not window.figure.free
    window.add_panels([_panel()])
    assert len(window.figure.panels) == 4
    window.close()
