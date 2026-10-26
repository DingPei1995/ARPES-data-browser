"""Tests for loader/cassiopee_spin.py -- MBS A-1 .krx and .txt files.

The files are written here in the layout read off the real ones (a 64-bit
pointer table, int32 images, a length-prefixed text header after each), so
the tests need no data shipped with the program."""
import os

import numpy as np
import pytest

from loader import cassiopee_spin as M
from loader import registry as L


def header(**over):
    fields = {
        "Lines": "72", "Frames Per Step": "500", "No. Steps": "20",
        "Pass Energy": "PE020", "Lens Mode": "L4Ang0d8", "No Scans": "1",
        "RegNo": "7", "Start K.E.": "40.7500000", "Step Size": "0.0050000",
        "End K.E.": "40.8500000", "Center K.E.": "40.8000000",
        "DeflX": "-2.7000000", "DeflY": "0.0000000", "AcqMode": "Swept",
        "NoS": "11", "RegName": "Cut34eV", "STim": "29/05/2026   00:21",
        "ScaleMult": "0.5000000", "ScaleMax": "2.5000000",
        "ScaleMin": "-2.5000000", "ScaleName": "Y Angle(Degrees)",
        "EScaleName": "Kinetic Energy (eV)",
        "XMapScaleName": "X Map Angle(Degrees)", "MapNoXSteps": "1",
        "MapStartX": "-10.0000000", "MapEndX": "10.0000000",
        "MapNoYSteps": "1", "MapStartY": "0", "MapEndY": "0",
        "MapCoordinate": "NONE", "MainDetOn": "YES", "SpinSystemOn": "NO",
        "MBS A1Soft": "LV 2016 Ver.16.0f5",
    }
    fields.update(over)
    extra = fields.pop("_extra", "")
    lines = [f"{k}\t{v}" for k, v in fields.items()]
    # RegNo appears twice in real headers
    lines.insert(20, "RegNo\t7")
    return "\r\n".join(lines) + ("\r\n" + extra if extra else "") + "\r\nDATA:\t\r\n"


def write_krx(path, images, headers, bits=64):
    n = len(images)
    width = bits // 8
    table_bytes = width * (1 + 3 * n)
    pos = (table_bytes + 3) // 4 + 5          # a few words of padding
    rows, blobs = [], []
    for image, text in zip(images, headers):
        ny, nx = image.shape
        rows.append((pos, ny, nx))
        raw = text.encode("latin-1")
        blob = (np.asarray(image, "<i4").tobytes()
                + np.int32(len(raw)).tobytes() + raw)
        blob += b"\x00" * ((-len(blob)) % 4)
        blobs.append((pos, blob))
        pos += len(blob) // 4 + 3
    dtype = "<i8" if bits == 64 else "<i4"
    table = np.array([3 * n] + [v for row in rows for v in row], dtype)
    size = pos * 4 + 64
    buffer = bytearray(size)
    buffer[:table.nbytes] = table.tobytes()
    for offset, blob in blobs:
        buffer[offset * 4: offset * 4 + len(blob)] = blob
    with open(path, "wb") as f:
        f.write(bytes(buffer))
    return path


def write_txt(path, text_header, energy, image):
    rows = ["\t".join([f"{e:.7f}"] + [f"{v:.7f}" for v in image[:, i]])
            for i, e in enumerate(energy)]
    with open(path, "wb") as f:
        f.write((text_header + "\r\n".join(rows) + "\r\n").encode("latin-1"))
    return path


@pytest.fixture
def cut_image():
    rng = np.random.default_rng(0)
    return rng.poisson(50, size=(11, 20)).astype(np.int32)


# -- the formats ---------------------------------------------------------------
@pytest.mark.parametrize("bits", [64, 32])
def test_cut_krx(tmp_path, cut_image, bits):
    path = write_krx(tmp_path / "Cut1.krx", [cut_image], [header()], bits)
    loader = L.detect(str(path))
    assert loader is not None and loader.name == M.CassiopeeSpinLoader.name
    entries = loader.list_entries(str(path))
    assert entries[0]["kind"] == "Cut" and "Cut34eV" in entries[0]["name"]
    scan = L.load(str(path))
    assert scan.kind == "cut"
    np.testing.assert_array_equal(scan.value, cut_image)
    np.testing.assert_allclose(scan.x, -2.5 + 0.5 * np.arange(11))
    np.testing.assert_allclose(scan.y, 40.75 + 0.005 * np.arange(20))
    assert scan.labels == {"x": "Y angle (deg)", "y": "Kinetic energy (eV)"}
    assert scan.info["pass_energy_eV"] == 20.0
    assert scan.info["lens_mode"] == "L4Ang0d8"
    assert scan.info["mbs.deflector_x_deg"] == -2.7
    assert scan.info["energy_reference"] == "kinetic"
    assert scan.info["mbs.RegNo"] == "7"


def test_txt_matches_krx(tmp_path, cut_image):
    krx = write_krx(tmp_path / "c.krx", [cut_image], [header()])
    energy = 40.75 + 0.005 * np.arange(20)
    txt = write_txt(tmp_path / "c.txt", header(), energy, cut_image)
    assert L.detect(str(txt)).name == M.CassiopeeSpinLoader.name
    a, b = L.load(str(krx)), L.load(str(txt))
    np.testing.assert_array_equal(a.value, b.value)
    np.testing.assert_allclose(a.y, b.y)
    np.testing.assert_allclose(a.x, b.x)


def test_map(tmp_path, cut_image):
    images = [cut_image + i for i in range(5)]
    h = header(MapNoXSteps="5", MapStartX="-2.0000000", MapEndX="2.0000000",
               MapCoordinate="X Direction", AcqMode="Fixed")
    path = write_krx(tmp_path / "Map.krx", images, [h] * 5)
    assert L.list_entries(str(path))[0]["kind"] == "Map"
    scan = L.load(str(path))
    assert scan.kind == "map" and scan.value.shape == (5, 11, 20)
    np.testing.assert_allclose(scan.x, [-2, -1, 0, 1, 2])
    np.testing.assert_array_equal(scan.value[3], cut_image + 3)
    assert scan.labels["x"] == "X map angle (deg)"
    assert scan.info["axis0.role"] == "angle"


def test_interrupted_map(tmp_path, cut_image):
    h = header(MapNoXSteps="5", MapStartX="-2", MapEndX="2",
               MapCoordinate="X Direction")
    path = write_krx(tmp_path / "Map.krx", [cut_image] * 3, [h] * 3)
    scan = L.load(str(path))
    np.testing.assert_allclose(scan.x, [-2, -1, 0])
    assert "interrupted" in scan.info["mbs.map_note"]


def test_two_direction_map_is_refused(tmp_path, cut_image):
    h = header(MapNoXSteps="2", MapNoYSteps="2")
    path = write_krx(tmp_path / "Map.krx", [cut_image] * 4, [h] * 4)
    with pytest.raises(ValueError, match="4-D"):
        L.load(str(path))


SPIN_EXTRA = ("SpinComp#0\t<0,0> +X  (GUI+Z)\r\nSpinComp#1\t<90,180> +X  (GUI-Z)"
              "\r\nSpinComp#2\t<0,0> -X  (GUI-Z)\r\nSpinComp#3\t<90,180> -X  (GUI+Z)"
              "\r\nManMagDir\t-X")


def test_spin_krx(tmp_path):
    rng = np.random.default_rng(1)
    spectra = [rng.poisson(1000, size=(1, 30)).astype(np.int32) for _ in range(4)]
    h = header(**{"No. Steps": "30", "Start K.E.": "16.9650000",
                  "Step Size": "0.0150000", "MainDetOn": "NO",
                  "SpinSystemOn": "YES", "NoS": "705", "_extra": SPIN_EXTRA})
    path = write_krx(tmp_path / "Cut_5S.krx", spectra, [h] * 4)
    assert L.list_entries(str(path))[0]["kind"] == "Spin EDC"
    scan = L.load(str(path))
    assert scan.kind == "spin_edc" and scan.value.shape == (30, 4)
    for i in range(4):
        np.testing.assert_array_equal(scan.value[:, i], spectra[i][0])
    np.testing.assert_allclose(scan.x, 16.965 + 0.015 * np.arange(30))
    np.testing.assert_array_equal(scan.y, np.arange(4))
    assert scan.info["spin.component.1"] == "<90,180> +X (GUI-Z)"
    assert scan.info["curve.channels"].split("|")[3] == "C3 <90,180> -X (GUI+Z)"
    assert scan.info["spin.magnetisation_now"] == "-X"


def test_spin_txt_one_block(tmp_path):
    rng = np.random.default_rng(2)
    table = rng.poisson(500, size=(4, 12)).astype(float)
    h = header(**{"No. Steps": "12", "MainDetOn": "NO", "SpinSystemOn": "YES",
                  "_extra": SPIN_EXTRA})
    path = write_txt(tmp_path / "s.txt", h, 40.75 + 0.005 * np.arange(12), table)
    scan = L.load(str(path))
    assert scan.kind == "spin_edc"
    np.testing.assert_array_equal(scan.value, table.T)


# -- detection -------------------------------------------------------------------
def test_not_claimed(tmp_path):
    scienta = tmp_path / "scienta.txt"
    scienta.write_text("[Region 1]\nDimension 1 scale=1 2 3\n")
    assert not M.CassiopeeSpinLoader().can_open(str(scienta))
    junk = tmp_path / "junk.krx"
    junk.write_bytes(os.urandom(256))
    assert not M.CassiopeeSpinLoader().can_open(str(junk))
    empty = tmp_path / "empty.krx"
    empty.write_bytes(b"")
    assert not M.CassiopeeSpinLoader().can_open(str(empty))


def test_header_parsing():
    h = M.parse_header(header(_extra=SPIN_EXTRA))
    assert h["_spin_components"][0] == "<0,0> +X  (GUI+Z)"
    assert h["Start K.E."] == "40.7500000"
    assert M.classify(h, 1) == "cut"
    assert M.classify(dict(h, MainDetOn="NO", SpinSystemOn="YES"), 4) == "spin_edc"
    assert M.classify(h, 3) == "map"


def test_pretty_labels():
    assert M._pretty_scale("Y Angle(Degrees)", "") == "Y angle (deg)"
    assert M._pretty_scale("Kinetic Energy (eV)", "") == "Kinetic energy (eV)"
    assert M._pretty_scale("", "fallback") == "fallback"
