"""The MAX IV Bloch (Scienta SES) reader.

The files are built here as SES writes them, from the headers in
``for loader/MaxIV``: version-5 Igor waves whose note is the SES header
(lines separated by ``\\r``, the run-mode table's columns by ``\\x0b``),
and the deflector-map zip with its ``viewer.ini``, ``Spectrum_<name>.ini``
and ``<name>.ini``. Every buffer holds ``np.arange`` in C order, so a value
says exactly which point and pixel it came from.
"""
import struct
import zipfile

import numpy as np
import pytest

from loader.formats import ses as maxiv_loader
from loader import registry

NA, NE = 5, 6


def _note(scale="Kinetic", run_mode=(), hv=43.9998687499):
    lines = ["[SES]", "Version=1.2.5", "Region Name=fcut",
             "Lens Mode=DA30L_01", "Pass Energy=20",
             f"Excitation Energy={hv}", f"Energy Scale={scale}",
             "Acquisition Mode=Swept", "Spectrum Name=fcut004",
             "Location=Bloch", "Sample=WS", "Date=2025-10-23",
             "Time=19:54:45", "ThetaX=0", "ThetaY=2",
             "[Run Mode Information]", *run_mode, "",
             "[User Interface Information]", "[Manipulator]", "A=-24.386",
             "P=-0.031", "T=-7.292", "X=-1.096", "Y=274.634", "Z=-2.373", ""]
    return "\r".join(lines)


def _write_ibw(path, data, starts, deltas, units, note, order="<"):
    """``data`` in C order, i.e. ``(chunks, layers, columns, rows)``."""
    data = np.ascontiguousarray(data, dtype=order + "f4")
    igor_shape = list(reversed(data.shape)) + [0] * (4 - data.ndim)
    body = data.tobytes()
    eunits = [u.encode("latin-1") for u in units] + [b""] * (4 - len(units))
    note_b = note.encode("latin-1")
    bin_header = struct.pack(order + "hhiiii4i4i3i", 5, 0, 320 + len(body),
                             0, len(note_b), 0, *[len(u) for u in eunits],
                             0, 0, 0, 0, 0, 0, 0)
    wave = bytearray(320)
    struct.pack_into(order + "i", wave, 12, data.size)
    struct.pack_into(order + "h", wave, 16, 2)
    wave[28:28 + 8] = b"wave001\0"
    struct.pack_into(order + "4i", wave, 68, *igor_shape)
    pad = lambda v: list(v) + [0.0] * (4 - len(v))          # noqa: E731
    struct.pack_into(order + "4d", wave, 84, *pad(deltas))
    struct.pack_into(order + "4d", wave, 116, *pad(starts))
    with open(path, "wb") as f:
        f.write(bin_header + bytes(wave) + body + note_b + b"".join(eunits))


def _cube(shape):
    return np.arange(int(np.prod(shape)), dtype=np.float32).reshape(shape)


def test_cut(tmp_path):
    path = str(tmp_path / "WS0004.ibw")
    data = _cube((NA, NE))
    _write_ibw(path, data, (38.2, -16.9), (0.003, 0.0474),
               ("Kinetic Energy [eV]", "Y-Scale [deg]"), _note())
    assert registry.detect(path).name == "MAX IV Bloch"
    [entry] = registry.list_entries(path)
    assert entry["kind"] == "Cut" and entry["title"] == "fcut004"
    scan = registry.load(path, entry["entry"])
    assert scan.kind == "cut"
    np.testing.assert_array_equal(np.asarray(scan.value), data)
    np.testing.assert_allclose(scan.x, -16.9 + 0.0474 * np.arange(NA))
    np.testing.assert_allclose(scan.y, 38.2 + 0.003 * np.arange(NE))
    assert scan.labels["y"] == "Kinetic energy (eV)"
    assert scan.info["photon_energy_eV"] == pytest.approx(43.9998687499)
    assert scan.info["analyser.pass_energy"] == 20.0
    assert scan.info["manipulator.T"] == -7.292
    assert scan.info["facility"] == "MAX IV Bloch"


def test_big_endian_wave(tmp_path):
    path = str(tmp_path / "be.ibw")
    data = _cube((NA, NE))
    _write_ibw(path, data, (38.2, -16.9), (0.003, 0.0474),
               ("Kinetic Energy [eV]", "Y-Scale [deg]"), _note(), order=">")
    np.testing.assert_array_equal(np.asarray(registry.load(path).value), data)


@pytest.mark.parametrize("step", [-0.1, 0.1])
def test_photon_energy_scan_in_binding_energy(tmp_path, step):
    """SES's positive binding axis becomes ascending E - E_F, the spectra
    flipped with it when it has to be; no work function is added. SES
    writes it descending (``Delta`` < 0), which needs no flip."""
    path = str(tmp_path / "WS0008.ibw")
    data = _cube((3, NA, NE))
    run_mode = ["Name=CIS", "Point 1=30 eV", "Point 2=31.5 eV",
                "Point 3=33 eV"]
    _write_ibw(path, data, (7.9, -16.9, 30.0), (step, 0.0474, 1.5),
               ("Binding Energy [eV]", "Y-Scale [deg]",
                "Photon Energy [eV]"), _note("Binding", run_mode, 30.0))
    [entry] = registry.list_entries(path)
    assert entry["kind"] == "kz map"
    scan = registry.load(path)
    assert scan.kind == "kz_map"
    assert scan.info["axis0.role"] == "photon_energy"
    np.testing.assert_allclose(scan.x, [30.0, 31.5, 33.0])
    energy = -(7.9 + step * np.arange(NE))
    expected = data
    if step > 0:
        energy, expected = energy[::-1], data[..., ::-1]
    np.testing.assert_allclose(scan.z, energy)
    np.testing.assert_array_equal(scan.value, expected)
    assert scan.labels["z"].startswith("E - E_F")


def test_photon_energy_scan_in_kinetic_energy(tmp_path):
    """A kinetic CIS keeps one binding window: the spectra stack as they
    are, on E_kin - hv + W."""
    path = str(tmp_path / "cis.ibw")
    data = _cube((3, NA, NE))
    _write_ibw(path, data, (25.0, -16.9, 30.0), (0.01, 0.0474, 1.5),
               ("Kinetic Energy [eV]", "Y-Scale [deg]",
                "Photon Energy [eV]"), _note())
    scan = registry.load(path)
    assert scan.kind == "kz_map"
    np.testing.assert_array_equal(scan.value, data)
    w = scan.info["work_function_eV"]
    np.testing.assert_allclose(scan.z, 25.0 + 0.01 * np.arange(NE) - 30 + w)


def _manipulator_rows(ys, xs, separator="\x0b"):
    rows = [separator.join(["Position [unit]", "Point [unit]",
                            "A", "P", "T", "X", "Y", "Z"])]
    for y in ys:
        for x in xs:
            rows.append(separator.join(f"{v:.3f}" for v in
                                       (y, x, 0.039, -0.004, -0.042, x, y,
                                        -2.275)))
    return ["Name=Manipulator Scan"] + rows


@pytest.mark.parametrize("separator", ["\x0b", "\r"])
def test_spatial_scan(tmp_path, separator):
    """``Point`` steps X and ``Position`` steps Y, read off the run-mode
    table -- with SES's own column separator, or with every value on a line
    of its own as in a text dump of the note."""
    xs = -1.76 + 0.05 * np.arange(4)
    ys = 274.523 + 0.05 * np.arange(3)
    path = str(tmp_path / "NNS0001.ibw")
    data = _cube((3, 4, NA, NE))
    _write_ibw(path, data, (38.2, -16.9, -1.76, 274.523),
               (0.0015, 0.0474, 0.05, 0.05),
               ("Kinetic Energy [eV]", "Y-Scale [deg]", "Point [unit]",
                "Position [unit]"),
               _note(run_mode=_manipulator_rows(ys, xs, separator)))
    [entry] = registry.list_entries(path)
    assert entry["kind"] == "SPEM"
    scan = registry.load(path)
    assert scan.kind == "spem_4d"
    np.testing.assert_allclose(scan.x, xs)
    np.testing.assert_allclose(scan.y, ys)
    assert scan.labels["x"] == "X (mm)" and scan.labels["y"] == "Y (mm)"
    np.testing.assert_array_equal(scan.value4d[1, 2], data[1, 2])
    np.testing.assert_array_equal(np.asarray(scan.value4d), data)


def test_line_scan_is_a_map_of_position(tmp_path):
    xs = np.arange(4) * 0.1
    path = str(tmp_path / "line.ibw")
    data = _cube((4, NA, NE))
    rows = ["Name=Manipulator Scan",
            "\x0b".join(["Point [unit]", "A", "P", "T", "X", "Y", "Z"])]
    rows += ["\x0b".join(f"{v:.3f}" for v in (x, 0, 0, 0, x, 1, 2))
             for x in xs]
    _write_ibw(path, data, (38.2, -16.9, 0.0), (0.0015, 0.0474, 0.1),
               ("Kinetic Energy [eV]", "Y-Scale [deg]", "Point [unit]"),
               _note(run_mode=rows))
    scan = registry.load(path)
    assert scan.kind == "map"
    assert scan.info["axis0.role"] == "position"
    assert scan.labels["x"] == "X (mm)"
    np.testing.assert_array_equal(np.asarray(scan.value), data)


def _ses_zip(path, regions, compression=zipfile.ZIP_DEFLATED):
    """``regions``: ``{name: data}`` with data ``(depth, height, width)``."""
    with zipfile.ZipFile(path, "w", compression) as z:
        viewer = ["[viewer]", "name=SES 3D Viewer",
                  "region_list=" + " ".join(f"region_{i}"
                                            for i in range(len(regions)))]
        for i, (name, data) in enumerate(regions.items()):
            depth, height, width = data.shape
            viewer += [f"[viewer.region_{i}]", f"ini_path={name}.ini",
                       f"name={name}", "channel_list=channel_0",
                       f"[viewer.region_{i}.channel_0]",
                       f"path=Spectrum_{name}.bin"]
            z.writestr(f"Spectrum_{name}.ini", "\n".join([
                "[spectrum]", f"width={width}", f"height={height}",
                f"depth={depth}", "byteperpoint=4",
                "widthoffset=38.233926", "widthdelta=0.001517",
                "heightoffset=-16.907537", "heightdelta=0.047426",
                "depthoffset=-15.000000", "depthdelta=0.500000",
                "widthlabel=Energy [eV]", "heightlabel=Thetax [deg]",
                "depthlabel=Thetay [deg]", f"name={name}"]))
            z.writestr(f"{name}.ini", _note().replace("\r", "\r\n"))
            z.writestr(f"Spectrum_{name}.bin",
                       np.ascontiguousarray(data, "<f4").tobytes())
        z.writestr("viewer.ini", "\n".join(viewer))


@pytest.mark.parametrize("compression",
                         [zipfile.ZIP_DEFLATED, zipfile.ZIP_STORED])
def test_deflector_map_zip(tmp_path, compression):
    path = str(tmp_path / "WS0001.zip")
    data = _cube((4, NA, NE))
    _ses_zip(path, {"qmap": data}, compression)
    assert registry.detect(path).name == "MAX IV Bloch"
    [entry] = registry.list_entries(path)
    assert entry["kind"] == "Map" and entry["name"] == "WS0001"
    scan = registry.load(path, entry["entry"])
    assert scan.kind == "map"
    assert registry.role_is_angle(scan)
    assert scan.labels["x"].startswith("Deflector Thetay")
    np.testing.assert_allclose(scan.x, -15 + 0.5 * np.arange(4))
    np.testing.assert_allclose(scan.k, -16.907537 + 0.047426 * np.arange(NA))
    np.testing.assert_allclose(scan.z, 38.233926 + 0.001517 * np.arange(NE))
    np.testing.assert_array_equal(np.asarray(scan.value), data)
    np.testing.assert_array_equal(scan.value[2, 1:3], data[2, 1:3])
    assert scan.info["photon_energy_eV"] == pytest.approx(43.9998687499)


def test_zip_with_two_regions(tmp_path):
    path = str(tmp_path / "two.zip")
    first, second = _cube((2, NA, NE)), _cube((3, NA, NE)) + 1000
    _ses_zip(path, {"a": first, "b": second})
    entries = registry.list_entries(path)
    assert [e["entry"] for e in entries] == ["region a", "region b"]
    np.testing.assert_array_equal(
        np.asarray(registry.load(path, "region b").value), second)
    with pytest.raises(KeyError):
        registry.load(path, "region c")


def test_other_files_are_not_claimed(tmp_path):
    loader = registry.get_loader("MAX IV Bloch")
    plain = tmp_path / "plain.zip"
    with zipfile.ZipFile(plain, "w") as z:
        z.writestr("readme.txt", "hello")
    assert not loader.can_open(str(plain))
    junk = tmp_path / "junk.ibw"
    junk.write_bytes(b"\0" * 1000)
    assert not loader.can_open(str(junk))
    # An Igor wave without an SES note is someone else's.
    other = str(tmp_path / "other.ibw")
    _write_ibw(other, _cube((NA, NE)), (0, 0), (1, 1), ("", ""), "")
    assert not loader.can_open(other)


def test_header_of_the_uploaded_cut():
    """The note of the real cut parses to the values it states."""
    import os
    here = os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))))
    path = os.path.join(here, "for loader", "MaxIV", "Cut.txt")
    if not os.path.exists(path):
        pytest.skip("reference dump not present")
    with open(path, encoding="latin-1") as f:
        text = f.read()
    sections = maxiv_loader.parse_ses_note(text.split("Note:", 1)[1])
    info = maxiv_loader._info(sections)
    assert info["photon_energy_eV"] == pytest.approx(43.9998687499)
    assert info["analyser.lens_mode"] == "DA30L_01"
    assert info["manipulator.A"] == -24.386
