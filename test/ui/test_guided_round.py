"""The smoothed tag, and the guided MDC/EDC fit."""
import sys

import numpy as np
import pytest

from PyQt5.QtWidgets import QApplication

from ui.viewers.widgets import MemoryData
import tools.processing.process as P


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication(sys.argv)


def _cut(info=None):
    k = np.linspace(-0.5, 0.5, 101)
    energy = np.linspace(-0.4, 0.05, 80)
    kk, ee = np.meshgrid(k, energy, indexing="ij")
    values = 1 / (1 + ((kk - (ee + 0.1) / 2.0) / 0.02) ** 2) + 0.05
    return MemoryData("cut", (k, energy), values,
                      {"x": "k (Å⁻¹)", "y": "E - E_F (eV)"},
                      source_label="k", source_info=dict(info or {}))


def test_a_smooth_step_tags_the_result_and_its_descendants():
    info = P.record_step({}, P.Step("smooth", {"method": "gaussian"}))
    assert P.is_smoothed(info) and P.tags_of(info) == ["smoothed"]
    later = P.record_step(info, P.Step("curvature", {"a0": 1.0}))
    assert P.is_smoothed(later)
    assert not P.is_smoothed(P.record_step({}, P.Step("derivative", {})))
    # read back off the history alone, as from a file saved before the tag
    assert P.is_smoothed({"proc.step.1": "smooth(method=gaussian)"})


def test_no_smooth_first_tip_on_smoothed_data(app):
    from ui.panels.process import ProcessDialog
    raw = ProcessDialog([("k", _cut())], start="curvature")
    assert not raw.smoothed and raw.curv_smooth_x.value() > 0
    smoothed = _cut(P.record_step({}, P.Step("smooth", {})))
    dialog = ProcessDialog([("k_sm", smoothed)], start="curvature")
    assert dialog.smoothed
    assert dialog.curv_smooth_x.value() == 0 and dialog.curv_smooth_y.value() == 0
    texts = [w.text() for w in dialog.findChildren(type(dialog.note))]
    assert not any("Tip:" in t for t in texts)
    assert any("already smoothed" in t for t in texts)
    raw.close()
    dialog.close()


def test_the_guided_fit_walks_to_a_fitted_series(app):
    import tools.spectra.peaks as PK
    from ui.panels.fit import FitPanel
    panel = FitPanel(_cut(), "k")
    panel.show()
    guide = panel.start_guide()
    assert guide.isVisible() and not panel.guide_button.isEnabled()
    for _ in range(4):
        guide.go(guide.step + 1)
    assert guide.step == 4 and not guide.next_button.isEnabled()
    assert panel.seed_button.isChecked()          # armed for the first click
    visited = []
    for _ in range(guide.SEED_LINES):
        position = panel.current_position()
        visited.append(position)
        centre = (position + 0.1) / 2.0
        along, line = panel.current_line()
        height, width = PK.suggest_seed(along, line, centre)
        panel.selected_band().seeds.append(PK.Seed(position, centre, width, height))
        panel.sync_bands()
        panel.seedPlaced.emit()
    assert len(set(np.round(visited, 6))) == guide.SEED_LINES
    assert guide.next_button.isEnabled()
    guide.go(6)
    assert not guide.next_button.isEnabled()
    panel.run_fit()
    assert panel.series is not None and guide.next_button.isEnabled()
    guide.go(len(guide.steps))
    assert not guide.isVisible() and panel.guide_button.isEnabled()
    panel.close()


def test_the_wheel_scrolls_a_scrolling_fit_panel(app):
    from PyQt5.QtCore import QPoint, QPointF, Qt
    from PyQt5.QtGui import QWheelEvent
    from ui.panels.fit import FitPanel
    panel = FitPanel(_cut(), "k")
    panel.resize(900, 300)                      # too short: the panel scrolls
    panel.show()
    app.processEvents()
    area = panel.centralWidget()
    bar = area.verticalScrollBar()
    assert bar.maximum() > 0
    viewport = panel.image.ui.graphicsView.viewport()
    before = panel.image.view.vb.viewRange()

    def turn(modifiers):
        event = QWheelEvent(QPointF(50, 50), QPointF(50, 50), QPoint(0, 0),
                            QPoint(0, -120), Qt.NoButton, modifiers,
                            Qt.NoScrollPhase, False)
        QApplication.sendEvent(viewport, event)

    turn(Qt.NoModifier)
    assert bar.value() > 0                         # the window scrolled...
    assert panel.image.view.vb.viewRange() == before   # ...the picture did not
    panel.close()
