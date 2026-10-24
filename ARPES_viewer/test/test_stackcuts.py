"""Combining several cuts into a map along an axis the user names
(tools/stackcuts.py, ui/stackcuts.py and the main list's right-click entry):
only cuts on one grid are accepted, the order can be changed, and the new
axis's role decides whether the result is a map or a kz map."""
import os

import numpy as np
import pytest

from tools import stackcuts as S
from ui.list_actions import availability

ANGLE = np.linspace(-10, 10, 21)
ENERGY = np.linspace(-0.5, 0.1, 31)
LABELS = {"x": "Angle (deg)", "y": "E - E_F (eV)"}


def cut(name, offset=0.0, angle=ANGLE, energy=ENERGY, labels=LABELS):
    values = np.full((angle.size, energy.size), offset, dtype=float)
    return (name, (angle, energy), values, dict(labels))


def test_matching_cuts_are_accepted():
    assert S.compatibility([cut("a"), cut("b", 1), cut("c", 2)]) == ""


@pytest.mark.parametrize("other, words", [
    (cut("b", energy=ENERGY[:-1]), "31"),
    (cut("b", energy=ENERGY + 0.05), "energy axis runs"),
    (cut("b", labels={"x": "k (Å⁻¹)", "y": "E - E_F (eV)"}), "invA"),
])
def test_mismatched_cuts_are_refused_by_name(other, words):
    problem = S.compatibility([cut("a"), other])
    assert "“b”" in problem and words in problem


def test_jitter_within_tolerance_is_the_same_grid():
    step = ENERGY[1] - ENERGY[0]
    assert S.compatibility([cut("a"), cut("b", energy=ENERGY + 0.01 * step)]) == ""


def test_one_cut_is_not_enough():
    assert S.compatibility([cut("a")])


def test_axis_stack_and_kind():
    np.testing.assert_allclose(S.stack_axis(4, 10, 5), [10, 15, 20, 25])
    with pytest.raises(ValueError):
        S.stack_axis(3, 0, 0)
    cube = S.stack([np.zeros((2, 3)), np.ones((2, 3))])
    assert cube.shape == (2, 2, 3) and cube[1].sum() == 6
    assert S.axis_label("Temperature", "K") == "Temperature (K)"
    assert S.axis_label("Gate", "") == "Gate"
    assert S.kind_for_role("photon_energy", "Angle (deg)") == "kz_map"
    assert S.kind_for_role("photon_energy", "k (Å⁻¹)") == "map"
    assert S.kind_for_role("temperature", "Angle (deg)") == "map"


def test_list_menu_entry_needs_two_or_more_cuts():
    c, m = ("c", "Cut"), ("m", "Map")
    assert availability("combine", [c, c])[0]
    assert availability("combine", [c, c, c, c])[0]
    assert not availability("combine", [c])[0]
    assert not availability("combine", [c, m])[0]


# -- the dialog and the launcher, with Qt ---------------------------------------
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


@pytest.fixture(scope="module")
def app():
    pytest.importorskip("PyQt5")
    from PyQt5.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def test_dialog_order_and_values(app):
    from ui.stackcuts import CombineCutsDialog
    dialog = CombineCutsDialog(["a", "b", "c"], reorder=True)
    dialog.role.setCurrentIndex(dialog.role.findData("gate_voltage"))
    assert dialog.unit.text() == "V"
    dialog.start.setValue(-2)
    dialog.step.setValue(0.5)
    dialog.list.setCurrentRow(2)
    dialog.move(-1)                              # c above b
    assert dialog.order() == [0, 2, 1]
    dialog.reverse()
    assert dialog.order() == [1, 2, 0]
    role, name, unit, values = dialog.axis()
    assert (role, name, unit) == ("gate_voltage", "Gate voltage", "V")
    np.testing.assert_allclose(values, [-2, -1.5, -1])
    assert dialog.list.item(0).text().startswith("-2 V")
    assert dialog.list.item(0).text().endswith("b")
    dialog.step.setValue(0)
    assert not dialog.buttons.button(dialog.buttons.Ok).isEnabled()


def _list_cuts(launcher, monkeypatch, cuts):
    from PyQt5 import QtCore
    from PyQt5.QtWidgets import QAbstractItemView, QListWidget, QListWidgetItem
    from ui.widgets import MemoryData

    class _UI:
        FilePathListWidget = QListWidget()
        statusbar = type("S", (), {"showMessage": staticmethod(lambda *a: None)})

    _UI.FilePathListWidget.setSelectionMode(QAbstractItemView.ExtendedSelection)
    monkeypatch.setattr(launcher, "ui", _UI)
    # A message box would wait for a click that never comes.
    monkeypatch.setattr(launcher.QMessageBox, "information",
                        lambda *args: pytest.fail(args[2]))
    monkeypatch.setattr(launcher, "loaded_items", {})
    monkeypatch.setattr(launcher, "autosave", lambda data, name: None)
    for name, axes, values, labels in cuts:
        data = MemoryData("cut", axes, values, labels, source_label=name)
        key = ("<converted>", name)
        launcher.loaded_items[key] = {"name": name, "kind": "Cut", "memory": data,
                                "entry": None, "path": ""}
        item = QListWidgetItem(name)
        item.setData(QtCore.Qt.UserRole, key)
        _UI.FilePathListWidget.addItem(item)
        item.setSelected(True)


def test_launcher_combines_in_a_chosen_order(app, monkeypatch):
    import ARPES_viewer as av
    _list_cuts(av, monkeypatch, [cut("T10", 10), cut("T20", 20), cut("T30", 30)])

    def run(dialog):
        dialog.role.setCurrentIndex(dialog.role.findData("temperature"))
        dialog.start.setValue(30)
        dialog.step.setValue(-10)
        dialog.reverse()
        return True

    result = av.combine_selected_cuts(ask=lambda names: False,
                                      dialog_hook=run)
    assert result.kind == "map"
    assert result.scan.info["axis0.role"] == "temperature"
    assert result.scan.labels["x"] == "Temperature (K)"
    np.testing.assert_allclose(result.scan.x, [30, 20, 10])
    assert [float(result.scan.value[i, 0, 0]) for i in range(3)] == [30, 20, 10]
    np.testing.assert_allclose(result.scan.k, ANGLE)
    np.testing.assert_allclose(result.scan.z, ENERGY)
    assert any(r.get("memory") is result for r in av.loaded_items.values())


def test_launcher_photon_energy_series_is_a_kz_map(app, monkeypatch):
    import ARPES_viewer as av
    _list_cuts(av, monkeypatch, [cut("a"), cut("b", 1)])

    def run(dialog):
        dialog.role.setCurrentIndex(dialog.role.findData("photon_energy"))
        dialog.start.setValue(60)
        dialog.step.setValue(2)
        return True

    result = av.combine_selected_cuts(ask=lambda names: True, dialog_hook=run)
    assert result.kind == "kz_map"
    np.testing.assert_allclose(result.scan.x, [60, 62])


def test_launcher_refuses_mismatched_cuts(app, monkeypatch):
    import ARPES_viewer as av
    shown = []
    _list_cuts(av, monkeypatch, [cut("a"), cut("b", energy=ENERGY[:-3])])
    monkeypatch.setattr(av.QMessageBox, "information",
                        lambda *args: shown.append(args[2]))
    assert av.combine_selected_cuts(ask=lambda n: True,
                                    dialog_hook=lambda d: True) is None
    assert shown and "“b”" in shown[0]


def test_combined_map_survives_a_save(app, tmp_path):
    from loader import session as nxs_session
    from loader.nxs_file import load_soleil_nxs, save_dataset
    from ui.stackcuts import combined_dataset
    from ui.widgets import MemoryData
    cuts = [(n, MemoryData("cut", axes, v, labels, source_label=n))
            for n, axes, v, labels in (cut("a"), cut("b", 1))]
    data = combined_dataset(cuts, [1, 0], "gate_voltage", "Gate", "V",
                            np.array([0.0, 0.5]), "stack")
    path = str(tmp_path / "stack.nxs")
    save_dataset(path, [nxs_session.scan_to_dict(data.kind, data.scan, "stack")])
    back = load_soleil_nxs(path)
    assert back.kind == "map" and back.value.shape == (2, 21, 31)
    np.testing.assert_allclose(back.x, [0, 0.5])
    assert back.labels["x"] == "Gate (V)"
    assert back.info["axis0.role"] == "gate_voltage"
    assert float(np.asarray(back.value)[0, 0, 0]) == 1.0
    back.close()
