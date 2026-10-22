"""The core-level check (tools/corelevel.py, ui/corelevel.py): the table
read from Table 1-1, the energy reference, peak finding and assignment on
a synthetic NbSCl spectrum, and the window run from the main list.

Needs Qt for the last part; run headless with QT_QPA_PLATFORM=offscreen.
"""
import os

import numpy as np
import pytest

from tools import corelevel as CL
from tools.corelevel_table import BINDING_ENERGIES, ELEMENTS

HV, W = 600.0, 4.5
#: (E_B, height) of the synthetic lines: Nb 3d, Nb 4p, Nb 4s, S 2p,
#: chloride Cl 2p, and adventitious carbon the sample does not contain
LINES = [(202.3, 4000), (205.0, 2700), (30.8, 600), (32.6, 300), (56.4, 300),
         (162.5, 2000), (163.6, 1000), (198.8, 1500), (200.4, 750),
         (284.8, 800)]


def _spectrum(lo=-1.0, hi=320.0, step=0.05, shift=0.0, seed=1, angles=30):
    """A kinetic-energy cut of NbSCl at hv = 600 eV: Gaussians on a sloping
    background, cut off at E_F, with Poisson noise. ``shift`` moves every
    core level deeper (charging)."""
    rng = np.random.default_rng(seed)
    ekin = np.arange(HV - W - hi, HV - W - lo, step)
    be = HV - W - ekin
    intensity = 200 + 0.3 * be
    for centre, height in LINES:
        intensity = intensity + height * np.exp(
            -0.5 * ((be - centre - shift) / (0.6 / 2.355)) ** 2)
    intensity = intensity / (np.exp(-be / 0.003) + 1) + 5
    frame = rng.poisson(np.tile(intensity / angles, (angles, 1))).astype(float)
    return np.linspace(-15, 15, angles), ekin, frame


def _report(shift=0.0, **kwargs):
    _, ekin, frame = _spectrum(shift=shift, **kwargs)
    e, i = CL.integrated_edc(frame, ekin)
    fermi = CL.locate_fermi_edge(e, i, "kinetic", photon_energy=HV,
                                 work_function=W)
    reference = CL.EnergyReference("kinetic", fermi.level, HV, W)
    return CL.check(reference.binding(e), i, ["Nb", "S", "Cl"], reference)


def _at(report, energy):
    return min(report.assignments, key=lambda a: abs(a.peak.binding - energy))


# -- the table --------------------------------------------------------------
def test_the_table_has_every_element_and_reads_its_values():
    assert len(ELEMENTS) == 92 and ELEMENTS[1] == "H" and ELEMENTS[92] == "U"
    assert BINDING_ENERGIES["C"] == {"1s": 284.2}
    assert BINDING_ENERGIES["Nb"]["3d5/2"] == 202.3
    assert BINDING_ENERGIES["Au"]["4f7/2"] == 84.0
    assert BINDING_ENERGIES["Cd"]["4d5/2"] == 10.7          # "l0.7" in the PDF
    assert BINDING_ENERGIES["Ca"]["3s"] == 44.3              # "44.3 †"
    assert "3s" not in BINDING_ENERGIES["Pm"]                # "—"
    assert "4d3/2" not in BINDING_ENERGIES["Gd"]
    assert BINDING_ENERGIES["Gd"]["4d5/2"] == 142.6


def test_the_generator_reproduces_the_table():
    """The committed table is what the generator makes from the PDF."""
    pdf = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__)))), "Table_1-1.pdf")
    import shutil
    import subprocess
    if not os.path.exists(pdf) or shutil.which("pdftotext") is None:
        pytest.skip("needs Table_1-1.pdf and pdftotext")
    from devtools.gen_corelevels import parse
    text = subprocess.run(["pdftotext", "-layout", pdf, "-"], check=True,
                          capture_output=True, text=True).stdout
    table = parse(text)
    assert {table[z][0]: table[z][1] for z in table} == BINDING_ENERGIES


@pytest.mark.parametrize("text, elements", [
    ("NbSCl", ["Nb", "S", "Cl"]), ("Nb, S, Cl", ["Nb", "S", "Cl"]),
    ("nb s cl", ["Nb", "S", "Cl"]), ("Bi2Se3", ["Bi", "Se"]),
    ("(Ba,K)Fe2As2", ["Ba", "K", "Fe", "As"]), ("Fe3O4 Fe", ["Fe", "O"])])
def test_elements_are_read_from_symbols_or_a_formula(text, elements):
    assert CL.parse_elements(text) == (elements, [])


def test_what_is_not_an_element_is_reported():
    assert CL.parse_elements("Nb Xy")[1] == ["Xy"]
    assert CL.parse_elements("nb foo")[1] == ["foo"]


def test_equal_doublet_components_merge_and_partners_pair_up():
    names = [line.orbital for line in CL.element_lines("Fe", chemical_states=False)]
    assert "3p" in names and "3p3/2" not in names        # 52.7 and 52.7
    line = CL.Line("Nb", "3d5/2", 202.3)
    assert line.partner_orbital() == "3d3/2" and line.is_major()
    assert CL.Line("Au", "4f5/2", 87.6).partner_orbital() == "4f7/2"


def test_contaminants_leave_out_radioactive_and_toxic_elements():
    common = CL.common_elements(exclude=["C"])
    assert common[0] == "O" and "C" not in common
    for symbol in ("U", "Tc", "Po", "Hg", "Cd", "Pb", "As", "Be", "H"):
        assert symbol not in common
    assert "Ag" in common and "Ta" in common


# -- the energy reference -----------------------------------------------------
@pytest.mark.parametrize("label, energy, hv, scale", [
    ("Kinetic energy (eV)", [550, 590], 600, "kinetic"),
    ("E - E_F (eV)", [-10, 0.2], 600, "relative"),
    ("Binding Energy (eV)", [0, 300], None, "binding"),
    ("Energy (eV)", [550, 590], 600, "kinetic"),
    ("Energy (eV)", [-5, 0.3], None, "relative")])
def test_the_scale_is_read_from_the_label_then_the_numbers(label, energy, hv, scale):
    assert CL.detect_scale(np.asarray(energy, float), label, hv) == scale


def test_the_fermi_edge_is_fitted_not_taken_from_the_photon_energy():
    _, ekin, frame = _spectrum()
    e, i = CL.integrated_edc(frame, ekin)
    # the nominal photon energy is 0.4 eV off; the edge says where E_F is
    fermi = CL.locate_fermi_edge(e, i, "kinetic", photon_energy=HV + 0.4,
                                 work_function=W)
    assert fermi.found and abs(fermi.level - (HV - W)) < 0.02
    reference = CL.EnergyReference("kinetic", fermi.level, HV + 0.4, W)
    assert abs(reference.implied_work_function() - (W + 0.4)) < 0.02


def test_a_window_without_the_edge_falls_back_to_hv_minus_w():
    _, ekin, frame = _spectrum(lo=150, hi=320)
    e, i = CL.integrated_edc(frame, ekin)
    fermi = CL.locate_fermi_edge(e, i, "kinetic", photon_energy=HV,
                                 work_function=W)
    assert not fermi.found and fermi.level == pytest.approx(HV - W)
    assert "not inside this window" in fermi.message


def test_binding_energies_from_each_scale():
    assert CL.EnergyReference("kinetic", 595.5).binding(395.5) == 200.0
    assert CL.EnergyReference("relative", 0.1).binding(-199.9) == 200.0
    assert CL.EnergyReference("binding", -0.2).binding(199.8) == 200.0


# -- peaks and assignment -----------------------------------------------------
def test_every_line_is_found_close_to_where_it_was_put():
    report = _report()
    found = sorted(p.binding for p in report.peaks)
    assert len(found) == len(LINES)
    for centre, _ in LINES:
        assert min(abs(np.asarray(found) - centre)) < 0.05


def test_sample_peaks_get_their_orbital_and_doublet_partner():
    report = _report()
    nb = _at(report, 202.3)
    assert nb.status == "sample" and nb.matches[0].line.name == "Nb 3d5/2"
    assert nb.matches[0].partner is not None
    assert _at(report, 205.0).matches[0].line.name == "Nb 3d3/2"
    assert _at(report, 162.5).matches[0].line.name == "S 2p3/2"
    # chloride sits 1.2 eV below the gas-phase table value: the chemical
    # state table is what recognises it
    cl = _at(report, 198.8)
    assert cl.status == "sample" and cl.matches[0].line.source == "chloride"


def test_a_peak_the_sample_does_not_explain_gets_common_elements():
    carbon = _at(_report(), 284.8)
    assert carbon.status == "unexplained"
    assert carbon.contaminants[0].line.element == "C"
    assert "adventitious" in carbon.contaminants[0].note
    assert carbon.label.startswith("C 1s")


def test_a_shifted_peak_is_left_to_the_user():
    report = _report(shift=1.2)
    nb = _at(report, 203.5)
    assert nb.status == "uncertain"
    assert nb.nearby[0].line.name == "Nb 3d5/2"
    assert nb.nearby[0].delta == pytest.approx(1.2, abs=0.05)


def test_a_rigid_shift_is_suggested_for_a_charged_sample():
    report = _report(shift=1.2)
    assert report.shift is not None
    assert report.shift[0] == pytest.approx(1.2, abs=0.06)
    assert report.shift[1] > report.shift[2]
    assert "rigid shift" in report.text()


def test_missing_lines_and_elements_out_of_the_window_are_reported():
    report = _report()
    missing = {line.name for line, _ in report.missing}
    assert "S 2s" in missing and "Nb 3d5/2" not in missing
    _, ekin, frame = _spectrum(lo=150, hi=210)
    e, i = CL.integrated_edc(frame, ekin)
    reference = CL.EnergyReference("kinetic", HV - W, HV, W)
    report = CL.check(reference.binding(e), i, ["Nb", "O"], reference)
    assert "O" in report.out_of_window and "1s" in report.out_of_window["O"]


def test_second_order_light_is_recognised():
    """Nb 2p3/2 (2371 eV) cannot be reached by 600 eV light; through the
    third order (1800 eV) it can't either, but through 2hv ... a level at
    E_B shows at E_B - hv. Put one there and see it named."""
    hv = 1200.0
    be = np.arange(0.0, 1300.0, 0.1)
    intensity = 100 + 1500 * np.exp(-0.5 * ((be - (2371.0 - hv)) / 0.4) ** 2)
    reference = CL.EnergyReference("binding", 0.0, hv, W)
    report = CL.check(be, np.round(intensity), ["Nb"], reference)
    peak = _at(report, 1171.0)
    assert peak.status == "unexplained"
    assert any(m.order == 2 and m.line.name == "Nb 2p3/2"
               for m in peak.higher_order)


# -- the window ---------------------------------------------------------------
@pytest.fixture(scope="module")
def app():
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    pytest.importorskip("PyQt5")
    from PyQt5.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def _cut_data():
    from ui.widgets import MemoryData
    angles, ekin, frame = _spectrum()
    return MemoryData("cut", (angles, ekin), frame,
                      {"x": "Angle (deg)", "y": "Kinetic energy (eV)"},
                      source_label="NbSCl",
                      source_info={"photon_energy_eV": HV,
                                   "work_function_eV": W})


def test_the_dialogs_and_window_run_end_to_end(app):
    from ui import corelevel as UC

    data = _cut_data()
    energy, intensity, label = UC.cut_edc(data)
    elements = UC.ElementsDialog(elements="NbSCl", settings={})
    assert elements.buttons.button(elements.buttons.Ok).isEnabled()
    reference = UC.EnergyReferenceDialog(energy, intensity, label,
                                         data.scan.info)
    assert reference.scale() == "kinetic"
    assert abs(reference.ef_box.value() - (HV - W)) < 0.02
    assert reference.estimate.found

    window = UC.start(data, "NbSCl", elements_dialog=elements,
                      reference_dialog=reference)
    report = window.report
    assert window.table.rowCount() == len(report.assignments) == len(LINES)
    assert "Nb 3d5/2" in window.table.item(
        report.assignments.index(_at(report, 202.3)), 5).text()

    created = []
    window.datasetsCreated.connect(created.extend)
    data_edc = window.edc_to_list()
    assert created == [data_edc] and data_edc.kind == "edc"
    assert data_edc.scan.labels["x"] == "Binding energy (eV)"

    # an assignment typed by hand shows up in the report
    window.table.item(0, 5).setText("my call")
    assert "Peak 1: my call" in window.report_text()
    window.close()


def test_the_shift_button_moves_the_reference(app):
    from ui import corelevel as UC
    from ui.widgets import MemoryData

    angles, ekin, frame = _spectrum(shift=1.2)
    data = MemoryData("cut", (angles, ekin), frame,
                      {"x": "Angle (deg)", "y": "Kinetic energy (eV)"},
                      source_label="charged",
                      source_info={"photon_energy_eV": HV,
                                   "work_function_eV": W})
    energy, intensity, label = UC.cut_edc(data)
    window = UC.CoreLevelWindow(
        energy, intensity, label, ["Nb", "S", "Cl"],
        UC.CL.EnergyReference("kinetic", HV - W, HV, W), "charged")
    assert window.shift_button.isEnabled()
    report = window.apply_shift()
    assert _at(report, 202.3).status == "sample"
    window.close()


def test_the_list_menu_offers_the_check_for_one_cut():
    from ui.list_actions import availability
    assert availability("corelevel", [("c", "Cut")])[0]
    assert not availability("corelevel", [("m", "Map")])[0]
    assert not availability("corelevel", [("c", "Cut"), ("d", "Cut")])[0]
