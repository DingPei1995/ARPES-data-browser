"""Enter in a field takes the value (and redraws what the window names as
its redraw); it never presses Next, Close, OK or any other button --
whatever the value just entered set going (a progress window that takes
the focus away was how Enter in the contour page's rotation pressed Next)."""
import time

import pytest
from PyQt5.QtCore import Qt
from PyQt5.QtTest import QTest
from PyQt5.QtWidgets import (QApplication, QCheckBox, QDialog, QDoubleSpinBox,
                             QMainWindow, QPushButton, QVBoxLayout, QWidget)

from ui.common import enterkeys


@pytest.fixture(scope="module")
def app():
    app = QApplication.instance() or QApplication([])
    enterkeys.install(app)
    return app


def _dialog(slow=False):
    dialog = QDialog()
    layout = QVBoxLayout(dialog)
    spin = QDoubleSpinBox()
    spin.setKeyboardTracking(False)
    layout.addWidget(spin)
    box = QCheckBox("a box")
    layout.addWidget(box)
    button = QPushButton("Next")
    button.setDefault(True)
    layout.addWidget(button)
    pressed = []
    button.clicked.connect(lambda: pressed.append(True))
    changed = []
    spin.valueChanged.connect(changed.append)
    if slow:
        from ui.common.jobs import run_blocking

        spin.valueChanged.connect(lambda _v: run_blocking(
            dialog, "long", lambda report: time.sleep(0.7)))
    dialog.show()
    return dialog, spin, box, button, pressed, changed


@pytest.mark.parametrize("slow", [False, True])
def test_enter_in_a_field_takes_the_value_and_presses_no_button(app, slow):
    dialog, spin, box, button, pressed, changed = _dialog(slow)
    spin.setFocus()
    app.processEvents()
    QTest.keyClicks(spin.lineEdit(), "7")
    QTest.keyClick(spin.lineEdit(), Qt.Key_Return)
    app.processEvents()
    assert changed and not pressed        # the value taken, Next not pressed
    box.setFocus()
    QTest.keyClick(box, Qt.Key_Return)
    assert not pressed
    button.setFocus()
    QTest.keyClick(button, Qt.Key_Return)  # a focused button still answers
    assert pressed
    dialog.close()


def test_a_window_that_names_its_redraw_gets_it(app):
    dialog, spin, _box, _button, pressed, _changed = _dialog()
    redrawn = []
    dialog.enter_redraws = lambda: redrawn.append(spin.value())
    spin.setFocus()
    QTest.keyClicks(spin.lineEdit(), "3")
    QTest.keyClick(spin.lineEdit(), Qt.Key_Return)
    assert redrawn and not pressed
    dialog.close()

    window = QMainWindow()
    holder = QWidget()
    inner = QVBoxLayout(holder)
    field = QDoubleSpinBox()
    inner.addWidget(field)
    window.setCentralWidget(holder)
    holder.enter_redraws = lambda: redrawn.append("holder")
    window.show()
    field.setFocus()
    QTest.keyClick(field.lineEdit(), Qt.Key_Return)
    assert redrawn[-1] == "holder"
    window.close()


def test_the_rotation_of_the_contour_plane_does_not_press_next(app, tmp_path):
    import numpy as np

    from tools.crystal import bandfiles
    from ui.panels.bands import BandsDialog

    a = 3.3
    lattice = np.array([[a, 0, 0], [-a / 2, a * np.sqrt(3) / 2, 0], [0, 0, 6.0]])
    hop = {(0, 0, 0): np.array([[0.5]])}
    for R in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (1, 1, 0), (-1, -1, 0)):
        hop[R] = np.array([[-0.15]])
    hr = str(tmp_path / "x_hr.dat")
    bandfiles.write_hr(hr, lattice, hop, 1)
    dialog = BandsDialog([])
    dialog.show()
    dialog.open_file(hr)
    dialog.set_goal("contours")
    dialog.stepper.go(1)
    dialog.cont_rotation.setFocus()
    app.processEvents()
    QTest.keyClicks(dialog.cont_rotation.lineEdit(), "5")
    QTest.keyClick(dialog.cont_rotation.lineEdit(), Qt.Key_Return)
    app.processEvents()
    assert dialog.stepper.index == 1
    assert "turned" in dialog.cont_label.text()
    dialog.close()


def test_enter_in_the_fermi_surface_energy_draws_once_without_an_error(
        app, tmp_path, monkeypatch):
    """The progress window of a slow surface takes the focus, and the energy
    field's editingFinished fires a second time while the first is drawn:
    that started a second job and failed ("another operation is already
    running"). Now it is one drawing, and no error."""
    import sys

    import numpy as np

    from tools.crystal import bands as B
    from ui.panels.bands import BandsDialog

    errors = []
    monkeypatch.setattr(sys, "excepthook", lambda *a: errors.append(a))
    rng = np.random.default_rng(0)
    size = 24                               # enough orbitals to be slow
    a = 3.3
    lattice = np.array([[a, 0, 0], [-a / 2, a * np.sqrt(3) / 2, 0], [0, 0, 6.0]])
    R = [(0, 0, 0), (1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
    H = np.zeros((len(R), size, size), complex)
    H[0] = np.diag(np.linspace(-2, 2, size))
    for i, r in enumerate(R[1:], start=1):
        t = rng.normal(0, 0.05, (size, size))
        H[i] = t
    for i in range(1, len(R), 2):           # H(-R) = H(R)^dagger
        H[i + 1] = H[i].conj().T
    model = B.WannierModel(lattice, np.array(R, float), H)
    dialog = BandsDialog([])
    dialog.show()
    dialog.set_goal("fs")
    dialog.model_path = "slow_hr.dat"
    dialog.set_model(model)
    calls = []
    real = B.WannierModel.energies

    def slow(self, k):
        time.sleep(0.05)
        return real(self, k)

    monkeypatch.setattr(B.WannierModel, "energies", slow)
    from tools.crystal import fermisurface as FS

    real_surface = FS.fermi_surface

    def counted(grid, energy, **options):
        if options.get("n") != 24:                 # the fine mesh: one drawing
            calls.append(energy)
        return real_surface(grid, energy, **options)

    monkeypatch.setattr(FS, "fermi_surface", counted)
    dialog.fs_detail.setCurrentIndex(0)
    dialog.stepper.go(1)                     # the first drawing, at E_F
    for _ in range(30):
        app.processEvents()
        time.sleep(0.02)
    calls.clear()
    finished = []
    dialog.fs_energy.editingFinished.connect(lambda: finished.append(True))
    QApplication.setActiveWindow(dialog)
    dialog.activateWindow()
    app.processEvents()
    dialog.fs_energy.setFocus()
    app.processEvents()
    assert QApplication.focusWidget() is not None
    dialog.fs_energy.lineEdit().selectAll()
    QTest.keyClicks(dialog.fs_energy.lineEdit(), "0.1")
    QTest.keyClick(dialog.fs_energy.lineEdit(), Qt.Key_Return)
    for _ in range(150):                    # the redraw a second call queues
        app.processEvents()
        time.sleep(0.02)
    assert len(finished) == 2                # the progress window took the focus
    assert not errors, errors
    assert calls == [pytest.approx(0.1)]     # drawn once, at the new energy
    assert dialog._fs_drawn is not None and dialog._fs_drawn.energy == pytest.approx(0.1)
    assert dialog.stepper.index == 1
    dialog.close()


def test_busy_work_runs_off_the_gui_thread_and_never_refuses(app):
    import threading

    from ui.common import jobs

    gui = threading.get_ident()
    where = jobs.busy_work(None, "a long step",
                           lambda: (time.sleep(0.5), threading.get_ident())[1])
    assert where != gui                          # in the background
    inner = []

    def outer(report):
        # a step started while a job runs: done here, not refused
        inner.append(jobs.busy_work(None, "inside", lambda: threading.get_ident()))
        return True

    assert jobs.run_blocking(None, "outer", outer)
    assert inner and inner[0] != gui             # ran in the outer job's thread
