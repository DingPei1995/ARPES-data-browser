"""The experiment log (tools.output.explog): what each kind of measurement
writes in a column, and how a log on disk is continued."""
import csv

import numpy as np
import pytest

from arpes_viewer.tools.output import explog as EL
from arpes_viewer.ui.viewers.widgets import MemoryData

SLIT = np.linspace(-15, 15, 31)
EKIN = np.linspace(84, 90, 61)


def _cut(**info):
    return MemoryData("cut", (SLIT, EKIN), np.ones((SLIT.size, EKIN.size)),
                      {"x": "Angle (deg)", "y": "Kinetic energy (eV)"},
                      source_label="cut", source_info=info,
                      source_motors={"X": 0.31, "Y": -5.58})


def _spem():
    x, y = np.linspace(38.3, 40.1, 19), np.linspace(-3.5, -1.5, 21)
    return MemoryData("spem_4d", (x, y, SLIT, EKIN),
                      np.ones((y.size, x.size, SLIT.size, EKIN.size),
                              dtype=np.float32),
                      {"x": "S-Z (mm)", "y": "S-T (mm)", "k": "Angle (deg)",
                       "z": "Energy (eV)"},
                      source_label="spem", source_info={"photon_energy_eV": 95.0},
                      source_motors={"S-Z": 38.3, "S-T": -3.5, "S-N": -0.2})


def _polar_map():
    polar = np.linspace(-14, 14, 57)
    info = {"axis0.role": "angle", "photon_energy_eV": 79.0,
            # the sample position, corrected along with the angle
            "manipulator.sax": np.array([-0.37, -0.2, 0.1, 0.69])}
    motors = {"polar range": "-14 .. 14 deg (57 steps)",
              "saz range": "-3.75 .. -3.12 mm (57 steps)"}
    return MemoryData("map", (polar, SLIT, EKIN),
                      np.ones((polar.size, SLIT.size, EKIN.size),
                              dtype=np.float32),
                      {"x": "polar (deg)", "k": "Angle (deg)",
                       "z": "Kinetic energy (eV)"},
                      source_label="map", source_info=info,
                      source_motors=motors)


def test_ranges_are_written_as_the_beamline_logs_write_them():
    assert EL.axis_range([3.0]) == "3"
    assert EL.axis_range([2.5, 2.5, 2.5]) == "2.5"
    assert EL.axis_range(np.linspace(-9, 9, 37)) == "-9/9/0.5"
    assert EL.axis_range(np.linspace(-9, 9, 37), "colon") == "-9:0.5:9"
    assert EL.axis_range([-0.37, -0.1, 0.69]) == "-0.37..0.69"
    assert EL.plain(np.float64(0.1)) == "0.1"
    assert EL.plain(b"LH") == "LH"


def test_a_cut_logs_its_fixed_position():
    template = EL.Template([EL.Column("File", "@file"),
                            EL.Column("Type", "@kind"),
                            EL.Column("X", "motor:X"),
                            EL.Column("hv", "info:photon_energy_eV")])
    row = EL.row(template, _cut(photon_energy_eV=90.0), {"name": "i05-1.nxs"})
    assert row == ["i05-1.nxs", "Cut", "0.31", "90"]


def test_a_spem_scan_logs_the_start_stop_and_step_of_its_stage():
    data = _spem()
    template = EL.Template([EL.Column("S-Z", "motor:S-Z"),
                            EL.Column("S-T", "motor:S-T", "spatial_y"),
                            EL.Column("S-N", "motor:S-N")])
    # S-Z by name (auto), S-T by its rule, S-N stood still
    assert EL.row(template, data) == ["38.3/40.1/0.1", "-3.5/-1.5/0.1", "-0.2"]
    template.columns[0].scan = "fixed"
    assert EL.row(template, data)[0] == "38.3"
    # the same columns on a cut: the motors' own values
    template.columns[1].scan = "spatial_y"
    cut = MemoryData("cut", (SLIT, EKIN), np.ones((SLIT.size, EKIN.size)),
                     {}, source_label="c",
                     source_motors={"S-Z": 39.8, "S-T": -1.2})
    assert EL.row(template, cut)[:2] == ["39.8", "-1.2"]


def test_a_manipulator_map_logs_its_angle_and_the_corrected_position():
    data = _polar_map()
    template = EL.Template([EL.Column("Type", "@kind"),
                            EL.Column("Polar", "info:manipulator.polar_deg"),
                            EL.Column("Z", "motor:saz"),
                            EL.Column("X", "info:manipulator.sax"),
                            EL.Column("Scan", "@scan1")])
    assert EL.row(template, data) == [
        "Map", "-14/14/0.5", "-3.75..-3.12", "-0.37..0.69",
        "polar (deg): -14/14/0.5"]
    # the angle rule names the scanned axis whatever the field is called
    template.columns[1] = EL.Column("D-Polar", "info:deflector", "angle_scan")
    assert EL.row(template, data)[1] == "-14/14/0.5"


def test_the_fields_offered_are_the_derived_ones_then_the_metadata():
    fields = EL.available_fields(_cut(photon_energy_eV=90.0, _group="/x"),
                                 {"name": "a.nxs"})
    ids = [f[0] for f in fields]
    assert ids[:len(EL.DERIVED)] == list(EL.DERIVED)
    assert "info:photon_energy_eV" in ids and "motor:X" in ids
    assert not any(i.startswith("info:_") for i in ids)
    assert dict((f[0], f[2]) for f in fields)["@file"] == "a.nxs"


def test_a_log_is_created_then_continued_after_its_last_record(tmp_path):
    path = str(tmp_path / "log.csv")
    template = EL.Template([EL.Column("File name", "@file"),
                            EL.Column("hν", "info:photon_energy_eV")])
    assert EL.append(path, template, [["a.nxs", "90"]]) == 1
    with open(path, "rb") as handle:
        assert handle.read(3) == b"\xef\xbb\xbf"       # for Excel
    # a spreadsheet saving it leaves empty rows after the last record
    with open(path, "a", encoding="utf-8", newline="") as handle:
        handle.write(",\r\n,\r\n")
    EL.append(path, template, [["b.nxs", "95"]])
    with open(path, encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.reader(handle))
    assert rows == [["File name", "hν"], ["a.nxs", "90"], ["b.nxs", "95"]]
    log = EL.read_log(path)
    assert log.header == ["File name", "hν"] and len(log.records) == 2
    assert EL.logged_names(log, template) == {"a.nxs", "b.nxs"}


def test_a_log_takes_records_only_under_its_own_header(tmp_path):
    path = str(tmp_path / "log.csv")
    first = EL.Template([EL.Column("File", "@file")])
    EL.append(path, first, [["a.nxs"]])
    second = EL.Template([EL.Column("File", "@file"), EL.Column("T", "info:T")])
    with pytest.raises(EL.FormatMismatch):
        EL.append(path, second, [["b.nxs", "12"]])
    assert EL.read_log(path).records == [["a.nxs"]]


def test_reformatting_rewrites_the_header_and_moves_the_records(tmp_path):
    path = str(tmp_path / "log.csv")
    old = EL.Template([EL.Column("File", "@file"), EL.Column("hv", "info:hv"),
                       EL.Column("Notes", "@blank")])
    EL.append(path, old, [["a.nxs", "90", "nice"], ["b.nxs", "95", ""]])
    new = old.copy()
    new.columns[1].header = "hν/Δ"                 # renamed: followed by uid
    new.columns.insert(0, new.columns.pop(2))     # moved to the front
    new.columns.append(EL.Column("T /K", "info:T"))
    log = EL.read_log(path)
    assert EL.reformat_plan(log, old, new) == ([2, 0, 1, None], [])
    EL.reformat(path, old, new)
    EL.append(path, new, [["", "c.nxs", "100", "12"]])
    with open(path, encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.reader(handle))
    assert rows == [["Notes", "File", "hν/Δ", "T /K"],
                    ["nice", "a.nxs", "90", ""], ["", "b.nxs", "95", ""],
                    ["", "c.nxs", "100", "12"]]
    # a column dropped from the format loses what it held: said beforehand
    fewer = EL.Template(new.columns[1:])
    assert EL.reformat_plan(EL.read_log(path), new, fewer)[1] == ["Notes"]


def test_a_log_unknown_here_is_matched_by_its_headings(tmp_path):
    path = tmp_path / "hand.csv"
    path.write_bytes("File name;Polar\r\ni05-1.nxs;-14\r\n".encode())
    log = EL.read_log(str(path))
    assert log.delimiter == ";" and log.header == ["File name", "Polar"]
    found = EL.template_from_header(log.header)
    assert [c.field for c in found.columns] == ["@file", "@blank"]
    assert EL.resolve(found, log) == (";", ",")


def test_separator_and_decimal_follow_excel_on_this_computer(tmp_path):
    path = str(tmp_path / "eu.csv")
    template = EL.Template([EL.Column("File", "@file"),
                            EL.Column("hv", "info:photon_energy_eV"),
                            EL.Column("X", "motor:X"),
                            EL.Column("Slit", "@angle")])
    log = EL.read_log(path)
    europe = (";", ",")
    delimiter, decimal = EL.resolve(template, log, europe)
    assert (delimiter, decimal) == (";", ",")
    row = EL.row(template, _cut(photon_energy_eV=90.5), {"name": "s_1.5K.nxs"},
                 decimal)
    # numbers take the decimal comma; a file name keeps its point
    assert row == ["s_1.5K.nxs", "90,5", "0,31", "-15/15/1"]
    EL.append(path, template, [row], europe)
    with open(path, encoding="utf-8-sig") as handle:
        assert handle.read() == "File;hv;X;Slit\n" \
            "s_1.5K.nxs;90,5;0,31;-15/15/1\n"
    # appending on a computer with other settings keeps the file's own
    assert EL.resolve(template, EL.read_log(path), (",", ".")) == (";", ",")
    explicit = EL.Template(template.columns, delimiter="\t", decimal=".")
    assert EL.resolve(explicit, EL.LogFile("x"), europe) == ("\t", ".")


def test_an_existing_log_without_a_bom_keeps_its_encoding(tmp_path):
    path = tmp_path / "old.csv"
    path.write_bytes(b"File,Type\r\ni05-1.nxs,quick map\r\n")
    template = EL.Template([EL.Column("File", "@file"),
                            EL.Column("Type", "@kind")])
    EL.append(str(path), template, [["i05-2.nxs", "Map"]])
    raw = path.read_bytes()
    assert not raw.startswith(b"\xef\xbb\xbf")
    assert raw.endswith(b"i05-1.nxs,quick map\r\ni05-2.nxs,Map\r\n")


def test_the_settings_and_formats_are_remembered(tmp_path):
    template = EL.Template([EL.Column("Polar", "info:p", "angle_scan")],
                           range_style="colon", delimiter=";")
    settings = EL.Settings(str(tmp_path), "beamtime", template)
    assert settings.path() == str(tmp_path / "beamtime.csv")
    settings.remember_format(settings.path(), template)
    EL.save_settings(settings)
    again = EL.load_settings()
    assert again.template.to_json() == template.to_json()
    assert again.format_for(settings.path(), ["Polar"]).columns[0].uid \
        == template.columns[0].uid
    assert again.format_for(settings.path(), ["Other"]) is None
    EL.save_template(template, str(tmp_path / "mine.json"))
    assert EL.load_template(str(tmp_path / "mine.json")).delimiter == ";"
    assert EL.load_template(str(tmp_path / "none.json")).headers() \
        == EL.default_template().headers()
