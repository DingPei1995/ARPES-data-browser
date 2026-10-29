"""The BESSY One Cube / One Square and SOLARIS URANOS (Scienta SES) readers.

The files are built here as SES writes them, with the headers of the
uploads in ``for loader/BESSY`` and ``for loader/SOLARIS``: Igor waves
whose note is the SES header, the deflector-map zip, and the multi-region
text export. Where the uploaded header dumps are present, waves are also
built from them directly -- their own note, their own scan axes -- with a
tiny angle and energy range, so the run-mode tables are read exactly as
they were written.
"""
import os
import re
import zipfile

import numpy as np
import pytest

from arpes_viewer.loader import registry
from arpes_viewer.loader.formats import ses
from test.loaders.test_maxiv_loader import _cube, _write_ibw

NA, NE = 5, 6
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))
UPLOADS = os.path.join(ROOT, "for loader")


def _note(location, instrument, manipulator, run_mode=(), scale="Kinetic",
          hv=100.0):
    lines = ["[SES]", "Version=1.2.5", "Region Name=scan",
             "Lens Mode=DA30L_01", "Pass Energy=20",
             f"Excitation Energy={hv}", f"Energy Scale={scale}",
             "Acquisition Mode=Fixed", "Spectrum Name=scan",
             f"Instrument={instrument}", f"Location={location}",
             "Sample=S", "Date=2025-12-19", "Time=10:47:29",
             "[Run Mode Information]", *run_mode, "",
             "[User Interface Information]", "[Manipulator]", *manipulator,
             ""]
    return "\r".join(lines)


def _cube_note(**kwargs):
    return _note("Scienta", "DA30L_8000-11623", ["Z=0.090", "Polar=0.000"],
                 **kwargs)


def _square_note(**kwargs):
    return _note("ARPES 1^2", "R8000-8ES202",
                 ["P-Axis=-6.800", "R-Axis=0.000", "Z-Axis=90.598"], **kwargs)


def _uranos_note(**kwargs):
    return _note("URANOS", "VLeed-12609",
                 ["X:=1.800", "Y:=1.300", "Z:=3.030", "R1:=90.000",
                  "R2:=0.000", "R3:=0.000"], **kwargs)


def _cut(path, note):
    _write_ibw(path, _cube((NA, NE)), (116.5, -17.6), (0.003, 0.0456),
               ("Kinetic Energy [eV]", "Y-Scale [deg]"), note)


# --------------------------------------------------------------------------
# Telling the beamlines apart
# --------------------------------------------------------------------------
@pytest.mark.parametrize("note, name", [
    (_cube_note(), "BESSY One Cube"),
    (_square_note(), "BESSY One Square"),
    (_uranos_note(), "SOLARIS URANOS"),
    (_note("Bloch", "DA30L", ["X=0"]), "MAX IV Bloch"),
])
def test_cut_is_claimed_by_its_beamline(tmp_path, note, name):
    path = str(tmp_path / "cut.ibw")
    _cut(path, note)
    assert registry.detect(path).name == name
    scan = registry.load(path)
    assert scan.kind == "cut" and scan.info["facility"] == name
    np.testing.assert_array_equal(np.asarray(scan.value), _cube((NA, NE)))


def test_one_cube_is_also_told_by_its_location(tmp_path):
    path = str(tmp_path / "cut.ibw")
    _cut(path, _note("BESSY 1^3", "DA30L", ["Z=0"]))
    assert registry.detect(path).name == "BESSY One Cube"


# --------------------------------------------------------------------------
# One Cube
# --------------------------------------------------------------------------
def test_one_cube_photon_energy_scan(tmp_path):
    """Binding energy, Delta < 0, hv in the layers: a kz map on E - E_F."""
    path = str(tmp_path / "kz.ibw")
    data = _cube((3, NA, NE))
    run_mode = ["Name=CIS", "Point 1=20 eV", "Point 2=22 eV", "Point 3=24 eV"]
    _write_ibw(path, data, (5.7628, -19.23, 20.0), (-0.001368, 0.0905, 2.0),
               ("Binding Energy [eV]", "Y-Scale [deg]",
                "Photon Energy [eV]"),
               _cube_note(run_mode=run_mode, scale="Binding", hv=20.0))
    assert registry.detect(path).name == "BESSY One Cube"
    [entry] = registry.list_entries(path)
    assert entry["kind"] == "kz map"
    scan = registry.load(path)
    assert scan.kind == "kz_map"
    np.testing.assert_allclose(scan.x, [20, 22, 24])
    np.testing.assert_allclose(scan.z, -(5.7628 - 0.001368 * np.arange(NE)))
    np.testing.assert_array_equal(scan.value, data)


def _zip(path, note, data, depth_label="Thetay [deg]"):
    depth, height, width = data.shape
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("viewer.ini", "\r\n".join([
            "[viewer.region_1]", "ini_path=qmap.ini", "name=qmap",
            "channel_list=channel_0", "[viewer]", "name=SES 3D Viewer",
            "region_list=region_1", "current_region=region_1",
            "[viewer.region_1.channel_0]", f"width={width}",
            f"height={height}", f"depth={depth}", "name=qmap",
            "path=Spectrum_qmap.bin"]))
        z.writestr("Spectrum_qmap.ini", "\r\n".join([
            "[spectrum]", f"width={width}", f"height={height}",
            f"depth={depth}", "byteperpoint=4", "widthoffset=49.236520",
            "widthdelta=0.001368", "heightoffset=-19.230769",
            "heightdelta=0.090498", "depthoffset=-15.000000",
            "depthdelta=0.500000", "widthlabel=Energy [eV]",
            "heightlabel=Thetax [deg]", f"depthlabel={depth_label}",
            "name=qmap"]))
        z.writestr("qmap.ini", note.replace("\r", "\r\n"))
        z.writestr("Spectrum_qmap.bin",
                   np.ascontiguousarray(data, "<f4").tobytes())


@pytest.mark.parametrize("note, name", [(_cube_note(), "BESSY One Cube"),
                                        (_uranos_note(), "SOLARIS URANOS")])
def test_deflector_map_zip(tmp_path, note, name):
    path = str(tmp_path / "qmap.zip")
    data = _cube((4, NA, NE))
    _zip(path, note, data)
    assert registry.detect(path).name == name
    [entry] = registry.list_entries(path)
    assert entry["kind"] == "Map"
    scan = registry.load(path, entry["entry"])
    assert scan.kind == "map" and registry.role_is_angle(scan)
    np.testing.assert_allclose(scan.x, -15 + 0.5 * np.arange(4))
    np.testing.assert_array_equal(np.asarray(scan.value), data)
    assert scan.info["facility"] == name


# --------------------------------------------------------------------------
# One Square
# --------------------------------------------------------------------------
def test_one_square_polar_map(tmp_path):
    """Layers labelled ``Point [deg]``; the run-mode table (here with every
    value on a line of its own, as in the upload) lists a running number and
    the ``P-Axis`` readback, which repeats the wave's own axis -- so the map
    is over the polar angle and can go to k."""
    polar = -15 + 0.5 * np.arange(4)
    rows = ["Point [deg]", "P-Axis", "R-Axis", "Z-Axis"]
    for i, p in enumerate(polar):
        rows += [str(i + 1 if i == 0 else i + 2), f"{p:.3f}", "0.000",
                 "0.000"]
    path = str(tmp_path / "map.ibw")
    data = _cube((4, NA, NE))
    _write_ibw(path, data, (97.97, -17.58, -15.0), (0.0015, 0.0456, 0.5),
               ("Kinetic Energy [eV]", "Y-Scale [deg]", "Point [deg]"),
               _square_note(run_mode=rows))
    assert registry.detect(path).name == "BESSY One Square"
    [entry] = registry.list_entries(path)
    assert entry["kind"] == "Map"
    scan = registry.load(path)
    assert scan.kind == "map" and registry.role_is_angle(scan)
    assert scan.labels["x"].startswith("Polar")
    np.testing.assert_allclose(scan.x, polar)


def test_one_square_manipulator_scan_is_spatial(tmp_path):
    """``Point`` steps the polar (P-Axis), ``Position`` steps Z: two motors,
    so a spatial scan, the slow one first."""
    points = -0.5 + 0.1 * np.arange(4)
    positions = -0.2 + 0.1 * np.arange(3)
    rows = ["Name=Manipulator Scan", "Position [mm]", "Point [deg]",
            "P-Axis", "R-Axis", "Z-Axis"]
    for z in positions:
        for p in points:
            rows += [f"{z:.3f}", f"{p:.3f}", f"{p:.3f}", "0.000", f"{z:.3f}"]
    path = str(tmp_path / "scan2D.ibw")
    data = _cube((3, 4, NA, NE))
    _write_ibw(path, data, (91.97, -17.58, -0.5, -0.2),
               (0.0015, 0.0456, 0.1, 0.1),
               ("Kinetic Energy [eV]", "Y-Scale [deg]", "Point [deg]",
                "Position [mm]"), _square_note(run_mode=rows))
    [entry] = registry.list_entries(path)
    assert entry["kind"] == "SPEM"
    scan = registry.load(path)
    assert scan.kind == "spem_4d"
    np.testing.assert_allclose(scan.y, positions)
    np.testing.assert_allclose(scan.x, points)
    assert scan.labels["y"] == "Z (mm)"
    assert scan.labels["x"].startswith("Polar")
    np.testing.assert_array_equal(np.asarray(scan.value4d), data)


# --------------------------------------------------------------------------
# URANOS
# --------------------------------------------------------------------------
@pytest.mark.parametrize("separator", ["\x0b", "\r"])
def test_uranos_spatial_scan_takes_its_axes_from_the_motors(tmp_path,
                                                            separator):
    """SES numbers the points ``1, 2, 3 ...`` (``[au]``); the motors that
    move along each dimension -- ``X:`` fast, ``Y:`` slow -- give the real
    positions."""
    xs = 1.2 + 0.16 * np.arange(4)
    ys = 0.6 + 0.16 * np.arange(3)
    header = ["Position [au]", "Point [au]", "X:", "Y:", "Z:", "R1:", "R2:",
              "R3:"]
    rows = [separator.join(header)]
    for j, y in enumerate(ys):
        for i, x in enumerate(xs):
            rows.append(separator.join(
                f"{v:.3f}" for v in (j + 1, i + 1, x, y, 3.03, 90, 4, 0)))
    path = str(tmp_path / "2D_scan.ibw")
    data = _cube((3, 4, NA, NE))
    _write_ibw(path, data, (43.11, -15.04, 1.0, 1.0),
               (0.0008, 0.055, 1.0, 1.0),
               ("Kinetic Energy [eV]", "Y-Scale [deg]", "Point [au]",
                "Position [au]"),
               _uranos_note(run_mode=["Name=Manipulator Scan"] + rows))
    assert registry.detect(path).name == "SOLARIS URANOS"
    scan = registry.load(path)
    assert scan.kind == "spem_4d"
    np.testing.assert_allclose(scan.x, xs)
    np.testing.assert_allclose(scan.y, ys)
    assert scan.labels["x"] == "X (mm)" and scan.labels["y"] == "Y (mm)"
    np.testing.assert_array_equal(np.asarray(scan.value4d), data)


def _write_text(path, regions, location="Solaris"):
    """SES 1.3.1 text export: ``regions`` is ``[(hv, energy, angle,
    values (angle, energy))]``."""
    lines = ["[Info]", "Number of Regions=%04d" % len(regions),
             "Version=1.3.1", ""]
    for n, (hv, e, a, v) in enumerate(regions, start=1):
        lines += [f"[Region {n}]", f"Region Name=Region_{n - 1}",
                  "Dimension 1 name=Kinetic Energy [eV]",
                  f"Dimension 1 size={len(e)}",
                  "Dimension 1 scale=" + " ".join(f"{x:.6g}" for x in e),
                  "Dimension 2 name=Y-Scale [deg]",
                  f"Dimension 2 size={len(a)}",
                  "Dimension 2 scale=" + " ".join(f"{x:.6g}" for x in a),
                  "", f"[Info {n}]", f"Region Name=Region_{n - 1}",
                  "Lens Mode=DA30L_01", "Pass Energy=20",
                  f"Excitation Energy={hv:g}", "Energy Scale=Kinetic",
                  "Instrument=DA30L-9ES211L", f"Location={location}",
                  "Date=2023-04-27", "Time=22:53:54",
                  f"[Run Mode Information {n}]", "Name=Add Dimension", "",
                  f"[User Interface Information {n}]", f"[Manipulator {n}]",
                  "X=-0.425", "", "", f"[Data {n}]"]
        for row, energy in enumerate(e):
            lines.append(f" {energy:.14E}  " + "  ".join(
                f"{x:.14E}" for x in np.asarray(v)[:, row]))
        lines.append("")
    with open(path, "w", newline="\r\n", encoding="latin-1") as f:
        f.write("\n".join(lines))


@pytest.mark.parametrize("offset", [0.0, 0.02])
def test_uranos_text_photon_energy_scan(tmp_path, offset):
    """One region per hv, each window moving with hv: a kz map on E - E_F.
    When the windows line up after referring, the spectra stack untouched;
    when they do not, they go onto the range they share."""
    angle = -16.4 + 0.048 * np.arange(NA)
    hvs = [18.0, 21.0, 24.0]
    regions = []
    for i, hv in enumerate(hvs):
        energy = 10.0 + (hv - 18.0) + i * offset + 0.04 * np.arange(NE)
        values = (i + 1) * 1000 + _cube((NA, NE))
        regions.append((hv, energy, angle, values))
    path = str(tmp_path / "Kz.txt")
    _write_text(path, regions)
    assert registry.detect(path).name == "SOLARIS URANOS"
    [entry] = registry.list_entries(path)
    assert entry["kind"] == "kz map" and entry["name"] == "Kz 18-24eV"
    scan = registry.load(path, entry["entry"])
    assert scan.kind == "kz_map" and scan.info["axis0.role"] == "photon_energy"
    np.testing.assert_allclose(scan.x, hvs)
    np.testing.assert_allclose(scan.k, angle, atol=1e-5)
    w = scan.info["work_function_eV"]
    assert scan.labels["z"].startswith("E - E_F")
    if offset == 0.0:
        np.testing.assert_allclose(scan.z, 10 - 18 + w + 0.04 * np.arange(NE))
        np.testing.assert_allclose(
            scan.value, np.stack([r[3] for r in regions]), rtol=1e-6)
    else:
        assert "energy_resampled" in scan.info
        assert scan.z[0] == pytest.approx(10 - 18 + w + 2 * offset)
        assert np.all(np.isfinite(scan.value))
    assert scan.info["facility"] == "SOLARIS URANOS"


def test_uranos_text_single_region_is_a_cut(tmp_path):
    angle = -16.4 + 0.048 * np.arange(NA)
    energy = 10.0 + 0.04 * np.arange(NE)
    path = str(tmp_path / "cut.txt")
    _write_text(path, [(18.0, energy, angle, _cube((NA, NE)))])
    [entry] = registry.list_entries(path)
    assert entry["kind"] == "Cut"
    scan = registry.load(path, entry["entry"])
    assert scan.kind == "cut" and scan.info["facility"] == "SOLARIS URANOS"
    np.testing.assert_allclose(scan.value, _cube((NA, NE)))


def test_a_text_file_from_elsewhere_stays_with_cassiopee(tmp_path):
    angle = -16.4 + 0.048 * np.arange(NA)
    energy = 10.0 + 0.04 * np.arange(NE)
    path = str(tmp_path / "cut.txt")
    _write_text(path, [(18.0, energy, angle, _cube((NA, NE)))],
                location="Soleil")
    assert registry.detect(path).name == "SOLEIL CASSIOPEE"


# --------------------------------------------------------------------------
# The uploaded headers themselves
# --------------------------------------------------------------------------
_DIM = re.compile(r"^(Rows|Columns|Layers|Chunks):\s*(\d+)\s+Start:\s*(\S+)"
                  r"\s+Delta:\s*(\S+)\s+Units:\s*(.*)$")


def _wave_from_dump(dump, path, n_angle=NA, n_energy=NE):
    """An ``.ibw`` with the dump's note and scan axes, a small detector."""
    if not os.path.exists(dump):
        pytest.skip("reference dump not present")
    with open(dump, encoding="latin-1") as f:
        text = f.read().replace("\r\n", "\n")
    head, note = text.split("Note:\n", 1)
    dims = []
    for line in head.splitlines():
        m = _DIM.match(line.strip())
        if m:
            dims.append((int(m.group(2)), float(m.group(3)),
                         float(m.group(4)), m.group(5).strip()))
    sizes = [n_energy, n_angle] + [d[0] for d in dims[2:]]
    data = _cube(tuple(reversed(sizes)))
    _write_ibw(path, data, [d[1] for d in dims], [d[2] for d in dims],
               [d[3] for d in dims], note.replace("\n", "\r"))
    return data


def test_uploaded_one_square_spem(tmp_path):
    path = str(tmp_path / "SPEM.ibw")
    data = _wave_from_dump(os.path.join(UPLOADS, "BESSY", "One Square",
                                        "SPEM.txt"), path)
    assert registry.detect(path).name == "BESSY One Square"
    scan = registry.load(path)
    assert scan.kind == "spem_4d"
    assert scan.labels["y"] == "Z (mm)" and scan.labels["x"].startswith("Polar")
    np.testing.assert_allclose(scan.y, -0.2 + 0.1 * np.arange(10))
    np.testing.assert_allclose(scan.x, -0.5 + 0.1 * np.arange(18))
    assert np.asarray(scan.value4d).shape == data.shape


def test_uploaded_one_square_map(tmp_path):
    path = str(tmp_path / "Map.ibw")
    _wave_from_dump(os.path.join(UPLOADS, "BESSY", "One Square", "Map.txt"),
                    path)
    scan = registry.load(path)
    assert scan.kind == "map" and registry.role_is_angle(scan)
    assert scan.labels["x"].startswith("Polar")
    np.testing.assert_allclose(scan.x, -15 + 0.5 * np.arange(61))


@pytest.mark.parametrize("station, hv", [("One Cube", 20 + 2 * np.arange(51)),
                                         ("One Square", 50 + 2 * np.arange(41))])
def test_uploaded_bessy_kz(tmp_path, station, hv):
    path = str(tmp_path / "kz.ibw")
    _wave_from_dump(os.path.join(UPLOADS, "BESSY", station, "KzMap.txt"),
                    path)
    assert registry.detect(path).name == f"BESSY {station}"
    scan = registry.load(path)
    assert scan.kind == "kz_map"
    assert scan.labels["z"].startswith("E - E_F")
    np.testing.assert_allclose(scan.x, hv)


def test_uploaded_uranos_spem(tmp_path):
    path = str(tmp_path / "SPEM.ibw")
    _wave_from_dump(os.path.join(UPLOADS, "SOLARIS", "SPEM.txt"), path)
    assert registry.detect(path).name == "SOLARIS URANOS"
    scan = registry.load(path)
    assert scan.kind == "spem_4d"
    np.testing.assert_allclose(scan.x, 1.2 + 0.16 * np.arange(11), atol=1e-9)
    np.testing.assert_allclose(scan.y, 0.6 + 0.16 * np.arange(11), atol=1e-9)
    assert scan.labels["x"] == "X (mm)" and scan.labels["y"] == "Y (mm)"


@pytest.mark.parametrize("station", ["One Cube", "One Square"])
def test_uploaded_bessy_cut(tmp_path, station):
    path = str(tmp_path / "cut.ibw")
    _wave_from_dump(os.path.join(UPLOADS, "BESSY", station, "Cut.txt"), path)
    assert registry.detect(path).name == f"BESSY {station}"
    assert registry.load(path).kind == "cut"


def test_uploaded_uranos_cut(tmp_path):
    path = str(tmp_path / "cut.ibw")
    _wave_from_dump(os.path.join(UPLOADS, "SOLARIS", "Cut.txt"), path)
    assert registry.detect(path).name == "SOLARIS URANOS"


def test_uploaded_uranos_text_scan():
    """The real text export, cut down to its first photon energy (and the
    header of the second): one spectrum, read as a cut."""
    path = os.path.join(UPLOADS, "SOLARIS",
                        "Kz_only first cut other cuts deleted.txt")
    if not os.path.exists(path):
        pytest.skip("reference file not present")
    assert registry.detect(path).name == "SOLARIS URANOS"
    [entry] = registry.list_entries(path)
    scan = registry.load(path, entry["entry"])
    assert scan.kind == "cut"
    assert np.asarray(scan.value).shape == (661, 96)
    assert scan.info["photon_energy_eV"] == 18.0


def test_uploaded_zip_headers():
    """The ``.ini`` files of the two map zips, zipped with a small ``.bin``
    in place of the real one."""
    for folder, name in ((("BESSY", "One Cube"), "qmap"),
                         (("SOLARIS",), "60 eV map")):
        base = os.path.join(UPLOADS, *folder)
        if not os.path.exists(os.path.join(base, f"{name}.ini")):
            pytest.skip("reference headers not present")
        with open(os.path.join(base, f"{name}.ini"), encoding="latin-1") as f:
            header = ses._keys(ses.parse_ses_note(f.read()).get("SES"))
        expected = "BESSY One Cube" if folder[0] == "BESSY" \
            else "SOLARIS URANOS"
        claimed = [loader.name for loader in registry.loaders()
                   if isinstance(loader, ses.SESLoader)
                   and loader.name != "MAX IV Bloch" and loader.claims(header)]
        assert claimed == [expected]
