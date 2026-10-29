"""The ALBA LOREA and Diamond I05 / I05-1 readers.

The files are built here in the layouts the beamlines' ``h5info`` dumps
show (``for loader/*.mat``): ALBA's ``entry1/data`` with an ``@axes`` list,
Diamond's ``entry1/analyser`` with ``@axis`` on every vector and the motors
of a two-dimensional scan stored on the whole grid. Every spectra buffer
holds ``np.arange``, so a value says exactly which point and pixel it came
from, and the tests check the opened array against it.
"""
import h5py
import numpy as np
import pytest

from arpes_viewer.loader.beamlines import alba as alba_loader, diamond as diamond_loader
from arpes_viewer.loader.formats import nexus_arpes
from arpes_viewer.loader import registry
from arpes_viewer.loader.nxs_file import HANDLES

NA, NE = 5, 6


def _cube(shape, integer=False):
    return np.arange(int(np.prod(shape)),
                     dtype=np.int32 if integer else np.float32).reshape(shape)


# --------------------------------------------------------------------------
# ALBA
# --------------------------------------------------------------------------
def _write_alba(path, axes, vectors, shape, extra=None):
    """``axes`` the ``@axes`` list, ``vectors`` ``{name: values}``."""
    data = _cube(shape)
    with h5py.File(path, "w") as f:
        f.attrs["default"] = "entry1"
        entry = f.create_group("entry1")
        entry.attrs["NX_class"] = "NXentry"
        entry.attrs["default"] = "data"
        entry["title"] = "a scan"
        group = entry.create_group("data")
        group.attrs["NX_class"] = "NXdata"
        group.attrs["signal"] = "data"
        group.attrs["axes"] = axes
        group["data"] = data
        for name, values in vectors.items():
            group[name] = np.asarray(values, dtype=float)
            group[name].attrs["units"] = "eV" if "energ" in name else "degrees"
        instrument = entry.create_group("instrument")
        analyser = instrument.create_group("analyser")
        analyser["pass_energy"] = np.array([20.0])
        analyser["data"] = group["data"]
        instrument.create_group("monochromator")["energy"] = np.array([70.0])
        manipulator = instrument.create_group("manipulator")
        manipulator["sapolar"] = np.array([1.0])
        entry.create_group("sample")["temperature"] = 18.0
        for path_, values in (extra or {}).items():
            entry[path_] = values
    return data


def _alba_cut(path):
    return _write_alba(path, ["angles", "energies"],
                       {"angles": np.linspace(-10, 10, NA),
                        "energies": np.linspace(60, 61, NE)}, (NA, NE))


def test_alba_cut(tmp_path):
    path = str(tmp_path / "042_Ge.nxs")
    data = _alba_cut(path)
    assert registry.detect(path).name == "ALBA LOREA"
    [entry] = registry.list_entries(path)
    assert entry["kind"] == "Cut" and entry["title"] == "a scan"
    with registry.load(path, entry["entry"]) as scan:
        assert scan.kind == "cut"
        np.testing.assert_array_equal(np.asarray(scan.value), data)
        np.testing.assert_allclose(scan.x, np.linspace(-10, 10, NA))
        np.testing.assert_allclose(scan.y, np.linspace(60, 61, NE))
        assert scan.labels["y"] == "Kinetic energy (eV)"
        assert scan.info["photon_energy_eV"] == 70.0
        assert scan.info["SampleTemperature_K"] == 18.0
        assert scan.info["analyser.pass_energy"] == 20.0
        assert scan.info["facility"].startswith("ALBA")


def test_alba_deflector_map(tmp_path):
    path = str(tmp_path / "018_Ge.nxs")
    data = _write_alba(path, ["defl_angles", "angles", "energies"],
                       {"defl_angles": np.linspace(-5, 5, 4),
                        "angles": np.linspace(-10, 10, NA),
                        "energies": np.linspace(60, 61, NE)}, (4, NA, NE))
    with registry.load(path) as scan:
        assert scan.kind == "map"
        assert scan.info["axis0.role"] == "angle"
        assert registry.role_is_angle(scan)
        assert scan.labels["x"].startswith("Deflector angle")
        np.testing.assert_allclose(scan.x, np.linspace(-5, 5, 4))
        np.testing.assert_array_equal(np.asarray(scan.value), data)
        np.testing.assert_array_equal(scan.value[2, 1:4, -1], data[2, 1:4, -1])


def test_alba_spatial_scan(tmp_path):
    path = str(tmp_path / "007_Ga1.nxs")
    data = _write_alba(path, ["saz", "sax", "angles", "energies"],
                       {"saz": np.linspace(0, 1, 3), "sax": np.linspace(2, 3, 4),
                        "angles": np.linspace(-10, 10, NA),
                        "energies": np.linspace(60, 61, NE)}, (3, 4, NA, NE))
    with registry.load(path) as scan:
        assert scan.kind == "spem_4d"
        np.testing.assert_allclose(scan.y, np.linspace(0, 1, 3))
        np.testing.assert_allclose(scan.x, np.linspace(2, 3, 4))
        assert scan.labels["x"].startswith("sax") \
            and scan.labels["y"].startswith("saz")
        np.testing.assert_array_equal(scan.value4d[1, 2], data[1, 2])
        np.testing.assert_array_equal(np.asarray(scan.value4d), data)


@pytest.mark.parametrize("positive_binding", [False, True])
def test_alba_photon_energy_scan(tmp_path, positive_binding):
    """The file's binding-energy axis is used, in E - E_F whichever sign
    convention it was written in, and the work function comes from
    comparing it with the kinetic axes."""
    hv = np.linspace(40, 60, 3)
    e_ef = np.linspace(-1.0, 0.2, NE)
    # A positive binding energy is written ascending, so E - E_F of each
    # pixel then runs the other way.
    binding = np.linspace(-0.2, 1.0, NE) if positive_binding else e_ef
    pixel = -binding if positive_binding else binding
    kinetic = hv[:, None] - 4.3 + pixel[None, :]
    path = str(tmp_path / "019_Ge.nxs")
    data = _write_alba(
        path, ["photon_energies", "angles", "binding_energies"],
        {"photon_energies": hv, "angles": np.linspace(-10, 10, NA),
         "binding_energies": binding}, (3, NA, NE),
        extra={"instrument/analyser/energies": kinetic})
    [entry] = registry.list_entries(path)
    assert entry["kind"] == "kz map"
    scan = registry.load(path)
    assert scan.kind == "kz_map"
    assert scan.info["axis0.role"] == "photon_energy"
    np.testing.assert_allclose(scan.x, hv)
    np.testing.assert_allclose(scan.z, e_ef)
    expected = data[..., ::-1] if positive_binding else data
    np.testing.assert_array_equal(scan.value, expected)
    assert scan.info["work_function_eV"] == pytest.approx(4.3)
    assert scan.labels["z"] == "E - E_F (eV)"
    assert "photon_energy_eV" not in scan.info


def test_alba_is_not_claimed_by_the_other_readers(tmp_path):
    path = str(tmp_path / "042_Ge.nxs")
    _alba_cut(path)
    for name in ("SOLEIL ANTARES", "Diamond I05 HR-ARPES",
                 "Diamond I05-1 nano-ARPES", "ALS MAESTRO"):
        assert not registry.get_loader(name).can_open(path), name


# --------------------------------------------------------------------------
# Diamond
# --------------------------------------------------------------------------
def _vector(group, name, values, axis, primary=False, units=None):
    group[name] = np.asarray(values, dtype=float)
    group[name].attrs["axis"] = axis
    if primary:
        group[name].attrs["primary"] = 1
    if units:
        group[name].attrs["units"] = units


def _write_diamond(path, scan_shape, vectors=(), energies=None, nano=False,
                   name=None):
    """``vectors`` are ``(name, values, axis, primary, units)``; the spectra
    are ``(*scan_shape, NA, NE)``."""
    n_scan = len(scan_shape)
    data = _cube(tuple(scan_shape) + (NA, NE), integer=True)
    with h5py.File(path, "w") as f:
        entry = f.create_group("entry1")
        entry.attrs["NX_class"] = "NXentry"
        entry["title"] = "a scan"
        group = entry.create_group("analyser")
        group.attrs["NX_class"] = "NXdata"
        group["data"] = data
        group["data"].attrs["signal"] = 1
        _vector(group, "angles", np.linspace(-15, 15, NA), str(n_scan + 1),
                True, "degree")
        energies = np.linspace(80, 81, NE) if energies is None else energies
        _vector(group, "energies", energies, str(n_scan + 2), True, "eV")
        _vector(group, "current", np.full(scan_shape, 300.0),
                ",".join(str(i + 1) for i in range(n_scan)), n_scan == 1, "mA")
        for vector in vectors:
            _vector(group, *vector)
        instrument = entry.create_group("instrument")
        if name is not None:
            instrument["name"] = name
        analyser = instrument.create_group("analyser")
        analyser["data"] = group["data"]
        analyser["pass_energy"] = np.array([10])
        instrument.create_group("monochromator")["energy"] = np.array([90.0])
        manipulator = instrument.create_group("manipulator")
        for motor in (("smpolar", "smx") if nano else ("sapolar", "sax")):
            manipulator[motor] = np.array([0.5])
        if nano:
            instrument.create_group("zone_plate")["zpx"] = np.array([1.0])
        entry.create_group("sample")["temperature"] = np.array([11.0])
    return data


@pytest.mark.parametrize("nano", [False, True])
def test_diamond_cut_is_a_scan_of_one_point(tmp_path, nano):
    path = str(tmp_path / "cut.nxs")
    data = _write_diamond(path, (1,), nano=nano)
    expected = "Diamond I05-1 nano-ARPES" if nano else "Diamond I05 HR-ARPES"
    assert registry.detect(path).name == expected
    [entry] = registry.list_entries(path)
    assert entry["kind"] == "Cut"
    with registry.load(path) as scan:
        assert scan.kind == "cut"
        assert scan.value.shape == (NA, NE)
        np.testing.assert_array_equal(np.asarray(scan.value), data[0])
        assert scan.value.dtype == np.float32
        assert scan.info["photon_energy_eV"] == 90.0
        assert scan.info["SampleTemperature_K"] == 11.0
        assert ("I05-1" in scan.info["facility"]) == nano


def test_diamond_branch_from_the_instrument_name(tmp_path):
    path = str(tmp_path / "plain.nxs")
    _write_diamond(path, (1,), name="i05-1")
    assert registry.detect(path).name == "Diamond I05-1 nano-ARPES"
    _write_diamond(path, (1,), nano=True, name="i05")
    assert registry.detect(path).name == "Diamond I05 HR-ARPES"


def test_diamond_branch_from_the_file_name(tmp_path):
    path = str(tmp_path / "i05-1-64334.nxs")
    with h5py.File(path, "w") as f:
        group = f.create_group("entry1/analyser")
        for name in ("data", "angles", "energies"):
            group[name] = np.zeros(2)
        f.create_group("entry1/instrument")
    assert diamond_loader.branch(h5py.File(path, "r"), path) == "I05-1"


def test_diamond_deflector_map(tmp_path):
    path = str(tmp_path / "map.nxs")
    data = _write_diamond(path, (4,), [
        ("deflector_x", np.linspace(-8, 8, 4), "1", True, None)])
    with registry.load(path) as scan:
        assert scan.kind == "map"
        assert scan.info["axis0.role"] == "angle"
        assert scan.labels["x"].startswith("Deflector angle")
        np.testing.assert_allclose(scan.x, np.linspace(-8, 8, 4))
        np.testing.assert_array_equal(np.asarray(scan.value), data)
        np.testing.assert_array_equal(scan.value[::-1, 0], data[::-1, 0])


def test_diamond_polar_map_keeps_the_positions_that_followed_the_angle(
        tmp_path):
    """A manipulator map: the beamline moved the sample along with the
    polar angle to keep the spot on it. Those positions are motors of the
    scan ("<name> range"), so the experiment log can say where it was."""
    from arpes_viewer.tools.output import explog as EL
    path = str(tmp_path / "polar.nxs")
    _write_diamond(path, (4,), [
        ("sapolar", np.linspace(-6, 6, 4), "1", True, "deg"),
        ("sax", np.array([0.1, 0.12, 0.17, 0.25]), "1", False, "mm"),
        ("say", np.full(4, -5.6), "1", False, "mm")])
    with registry.load(path) as scan:
        assert scan.kind == "map"
        np.testing.assert_allclose(scan.x, np.linspace(-6, 6, 4))
        assert scan.fourd_info["sax range"] == "0.1 .. 0.25 mm (4 steps)"
        assert "say range" not in scan.fourd_info         # it stood still

        class Data:
            kind = scan.kind
        Data.scan = scan
        template = EL.Template([EL.Column("Polar", "motor:sapolar"),
                                EL.Column("X", "info:manipulator.sax")])
        assert EL.row(template, Data()) == ["-6/6/4", "0.1..0.25"]


@pytest.mark.parametrize("hv_name,nano", [("energy", False), ("value", True)])
def test_diamond_photon_energy_scan(tmp_path, hv_name, nano):
    """Each spectrum's own kinetic window, referred by its own hv: when the
    analyser followed hv exactly, the spectra stack untouched."""
    hv = np.array([60.0, 50.0, 70.0])                 # not in order
    e_ef = np.linspace(-1.0, 0.25, NE)
    kinetic = hv[:, None] - 4.5 + e_ef[None, :]
    path = str(tmp_path / "kz.nxs")
    data = _write_diamond(path, (3,), [
        (hv_name, hv, "1", True, "eV"),
        ("centre_energy", hv - 4.5, "1", hv_name == "energy", None)],
        energies=kinetic, nano=nano)
    [entry] = registry.list_entries(path)
    assert entry["kind"] == "kz map"
    scan = registry.load(path)
    assert scan.kind == "kz_map"
    np.testing.assert_allclose(scan.x, np.sort(hv))
    np.testing.assert_allclose(scan.z, e_ef)
    np.testing.assert_array_equal(scan.value, data[np.argsort(hv)])
    assert scan.info["work_function_eV"] == 4.5
    assert "energy_resampled" not in scan.info


def test_diamond_photon_energy_scan_with_moving_windows_is_resampled(tmp_path):
    hv = np.array([50.0, 60.0, 70.0])
    kinetic = hv[:, None] - 4.5 + np.linspace(-1.0, 0.25, NE)[None, :] \
        + np.array([0.0, 0.1, 0.2])[:, None]
    path = str(tmp_path / "kz.nxs")
    _write_diamond(path, (3,), [("energy", hv, "1", True, "eV")],
                   energies=kinetic)
    scan = registry.load(path)
    assert scan.kind == "kz_map"
    assert scan.z[0] == pytest.approx(-0.8) and scan.z[-1] <= 0.25 + 1e-9
    assert scan.value.shape[:2] == (3, NA)
    assert "energy_resampled" in scan.info


@pytest.mark.parametrize("x_fast", [True, False])
def test_diamond_spatial_scan_reads_the_motors_off_the_grid(tmp_path, x_fast):
    """Both motors are stored on the whole (slow, fast) grid; which one
    steps along which dimension is read from their values."""
    slow, fast = np.meshgrid(np.arange(3), np.arange(4), indexing="ij")
    first = 1.0 + 0.1 * fast
    second = 5.0 + 0.2 * slow
    x, z = (first, second) if x_fast else (second, first)
    path = str(tmp_path / "spem.nxs")
    data = _write_diamond(path, (3, 4), [
        ("sax", x, "1,2", True, "mm"), ("saz", z, "1,2", False, "mm")])
    with registry.load(path) as scan:
        assert scan.kind == "spem_4d"
        assert scan.value4d.shape == (3, 4, NA, NE)
        np.testing.assert_allclose(scan.y, 5.0 + 0.2 * np.arange(3))
        np.testing.assert_allclose(scan.x, 1.0 + 0.1 * np.arange(4))
        fast_name = "sax" if x_fast else "saz"
        assert scan.labels["x"].startswith(fast_name)
        np.testing.assert_array_equal(np.asarray(scan.value4d), data)
        np.testing.assert_array_equal(scan.value4d[2, :, 1], data[2, :, 1])


def test_diamond_two_axis_scan_without_a_spatial_pair_is_one_map_per_step(
        tmp_path):
    polar, defl = np.meshgrid(np.array([1.0, 2.0]), np.linspace(-4, 4, 3),
                              indexing="ij")
    path = str(tmp_path / "polar_defl.nxs")
    data = _write_diamond(path, (2, 3), [
        ("sapolar", polar, "1,2", True, "deg"),
        ("deflector_x", defl, "1,2", True, None)])
    entries = registry.list_entries(path)
    assert [e["entry"] for e in entries] == ["point 0", "point 1"]
    assert all(e["kind"] == "Map" for e in entries)
    with registry.load(path, "point 1") as scan:
        assert scan.kind == "map"
        np.testing.assert_allclose(scan.x, np.linspace(-4, 4, 3))
        np.testing.assert_array_equal(np.asarray(scan.value), data[1])
    with pytest.raises(KeyError):
        registry.load(path, "point 5")


def test_diamond_is_not_claimed_by_the_other_readers(tmp_path):
    path = str(tmp_path / "cut.nxs")
    _write_diamond(path, (1,))
    for name in ("SOLEIL ANTARES", "ALBA LOREA", "ALS MAESTRO",
                 "Diamond I05-1 nano-ARPES"):
        assert not registry.get_loader(name).can_open(path), name


def test_the_file_closes_when_the_scan_does(tmp_path):
    path = str(tmp_path / "map.nxs")
    _write_diamond(path, (4,), [
        ("deflector_x", np.linspace(-8, 8, 4), "1", True, None)])
    scan = registry.load(path)
    assert HANDLES.is_open(path)
    scan.close()
    assert not HANDLES.is_open(path)
    # a kz map is read whole and keeps nothing open
    _write_diamond(path, (3,), [("energy", [50.0, 60, 70], "1", True, "eV")],
                   energies=np.array([50.0, 60, 70])[:, None] - 4.5
                   + np.linspace(-1, 0.2, NE)[None, :])
    registry.load(path)
    assert not HANDLES.is_open(path)


# --------------------------------------------------------------------------
# The lazy view
# --------------------------------------------------------------------------
def test_lazy_view_indexing_matches_numpy(tmp_path):
    path = str(tmp_path / "v.h5")
    data = _cube((3, 4, 5, 6))
    with h5py.File(path, "w") as f:
        f["d"] = data
    with h5py.File(path, "r") as f:
        perm = (1, 0, 3, 2)
        view = nexus_arpes.NexusView(f["d"], perm)
        expected = np.transpose(data, perm)
        assert view.shape == expected.shape
        for key in [(1,), (slice(None), 2), (Ellipsis, 3), (-1, 0, 1, 2),
                    (slice(None, None, -1), 1), (1, [0, 2], slice(1, 4)),
                    (slice(1, 3), Ellipsis, slice(None, None, 2))]:
            np.testing.assert_array_equal(view[key], expected[key], str(key))
        pinned = nexus_arpes.NexusView(f["d"], (0, 1, 2, 3), {1: 2})
        assert pinned.shape == (3, 5, 6)
        np.testing.assert_array_equal(np.asarray(pinned), data[:, 2])
        np.testing.assert_array_equal(pinned[1, :, 0], data[1, 2, :, 0])
