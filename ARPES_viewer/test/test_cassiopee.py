"""The SOLEIL CASSIOPEE reader.

The files are built here rather than shipped: a real spectrum is 6 MB of
text, and what has to be checked is the *format* -- section headers, the
axis scales written out number by number, the leading column of the data
block -- which a four-by-three file exercises exactly as well as a real one.
The writer below is the format as the Scienta software emits it, CRLF and
all, and was checked against a real ``0001.txt`` from the beamline.

The parts worth testing are the ones where this reader had to decide
something the file does not say: which quantity a folder stepped, whether a
scale that is not quite linear is rounding or a real non-linearity, and what
a file name's ``_12_ROI1_`` tail means when the base has underscores in it
too.
"""
import os
import warnings

import numpy as np
import pytest

from loader import cassiopee, registry


# --------------------------------------------------------------------------
# Building files that look like the beamline's
# --------------------------------------------------------------------------
def write_spectrum(path, energy, angle, values, region_name="fixed_cut",
                   excitation_energy=150.0, extra_regions=()):
    """A Scienta text spectrum. ``values`` is ``(angle, energy)``."""
    lines = ["[Info]", "Number of Regions=%04d" % (1 + len(extra_regions)),
             "Version=1.3.1", ""]
    regions = [(region_name, energy, angle, values)] + list(extra_regions)
    for number, (name, e, a, v) in enumerate(regions, start=1):
        lines += [
            f"[Region {number}]",
            f"Region Name={name}",
            "Dimension 1 name=Kinetic Energy [eV]",
            f"Dimension 1 size={len(e)}",
            "Dimension 1 scale=" + " ".join(f"{x:.5f}" for x in e),
            "Dimension 2 name=Y-Scale [deg]",
            f"Dimension 2 size={len(a)}",
            "Dimension 2 scale=" + " ".join(f"{x:.5f}" for x in a),
            "",
            f"[Info {number}]",
            f"Region Name={name}",
            "Lens Mode=Angular30",
            "Pass Energy=20",
            f"Excitation Energy={excitation_energy:.4f}",
            "Location=Soleil",
            "Date=3/23/2025",
            "Time=17:02:41",
            "",
            f"[Run Mode Information {number}]",
            "Name=Fixed Energies",
            "",
            f"[Data {number}]",
        ]
        for row, energy_value in enumerate(e):
            numbers = " ".join(f"{x:.5f}" for x in np.asarray(v)[:, row])
            lines.append(f"{energy_value:.5f} {numbers}")
        lines.append("")
    with open(path, "w", newline="\r\n", encoding="latin-1") as handle:
        handle.write("\n".join(lines))
    return path


def write_parameters(path, theta=0.0, photon_energy=150.0, phi=18.0,
                     spacing=" : "):
    """The ``_i`` file beside a spectrum. ``spacing`` varies the punctuation
    because the real files are hand-formatted and it does."""
    lines = [
        f"x (mm){spacing}1.2340",
        f"y (mm){spacing}-0.5000",
        f"z (mm){spacing}10.0000",
        f"theta (deg){spacing}{theta:.4f}",
        f"phi (deg){spacing}{phi:.4f}",
        f"tilt (deg){spacing}0.5000",
        f"P(mbar){spacing}3.2e-11",
        f"hv mono (eV){spacing}{photon_energy:.4f}",
        f"T_B (K){spacing}17.5",
        f"Polarisation [0:LV, 1:LH, 2:AV, 3:AH, 4:CR]{spacing}1",
    ]
    with open(path, "w", newline="\r\n", encoding="latin-1") as handle:
        handle.write("\n".join(lines))
    return path


def make_series(folder, base="Map80eV_Phi18", count=4, thetas=None,
                photon_energies=None, shape=(5, 4), roi=1):
    """A numbered folder series, plus its parameter files."""
    angle = np.linspace(-15.0, 15.0, shape[0])
    energy = np.linspace(115.6, 116.0, shape[1])
    paths = []
    for index in range(1, count + 1):
        values = np.full(shape, float(index))
        spectrum = os.path.join(folder, f"{base}_{index}_ROI{roi}_.txt")
        write_spectrum(spectrum, energy, angle, values,
                       excitation_energy=(photon_energies[index - 1]
                                          if photon_energies else 150.0))
        write_parameters(
            os.path.join(folder, f"{base}_{index}_i.txt"),
            theta=(thetas[index - 1] if thetas else 0.0),
            photon_energy=(photon_energies[index - 1]
                           if photon_energies else 150.0))
        paths.append(spectrum)
    return paths


@pytest.fixture
def spectrum(tmp_path):
    angle = np.linspace(-15.0, 15.0, 5)
    energy = np.linspace(115.6, 116.0, 4)
    values = np.arange(20, dtype=float).reshape(5, 4)
    path = write_spectrum(str(tmp_path / "0001.txt"), energy, angle, values)
    return path, energy, angle, values


# --------------------------------------------------------------------------
# Reading one spectrum
# --------------------------------------------------------------------------
def test_a_spectrum_comes_back_as_angle_by_energy(spectrum):
    path, energy, angle, values = spectrum
    regions = cassiopee.parse_scienta(path)
    assert len(regions) == 1
    region = regions[0]
    assert region.name == "fixed_cut"
    # The file is written a row per energy; the cut is (angle, energy).
    assert region.values.shape == (angle.size, energy.size)
    assert np.allclose(region.values, values)
    assert np.allclose(region.energy, energy)
    assert np.allclose(region.angle, angle)


def test_the_leading_column_of_each_data_row_is_dropped(spectrum):
    """Each data row starts with its own energy, repeating Dimension 1.
    Keeping it would put the energy scale into the counts as a bright edge
    channel -- the MATLAB loader drops it, and so must this."""
    path, energy, _angle, values = spectrum
    region = cassiopee.parse_scienta(path)[0]
    assert region.values.shape[1] == energy.size
    assert np.allclose(region.values[:, 0], values[:, 0])
    assert not np.any(np.isclose(region.values, energy[0] + 0.0)
                      & (region.values > 100))


def test_axis_labels_are_taken_from_the_file(spectrum):
    region = cassiopee.parse_scienta(spectrum[0])[0]
    assert region.energy_label == "Kinetic Energy (eV)"
    assert region.angle_label == "Y-Scale (deg)"


def test_acquisition_settings_are_kept(spectrum):
    region = cassiopee.parse_scienta(spectrum[0])[0]
    assert region.info["photon_energy_eV"] == 150.0
    assert region.info["lens_mode"] == "Angular30"
    assert region.info["beamline"] == "Soleil"
    assert region.info["scienta.Pass Energy"] == "20"


def test_every_region_of_a_multi_region_file_is_read(tmp_path):
    angle = np.linspace(-10.0, 10.0, 4)
    energy = np.linspace(20.0, 20.3, 3)
    second = ("second_cut", energy, angle, np.ones((4, 3)) * 7.0)
    path = write_spectrum(str(tmp_path / "two.txt"), energy, angle,
                          np.zeros((4, 3)), extra_regions=[second])
    regions = cassiopee.parse_scienta(path)
    assert [r.name for r in regions] == ["fixed_cut", "second_cut"]
    assert np.allclose(regions[1].values, 7.0)


def test_a_data_block_of_the_wrong_width_is_an_error_not_a_reshape(tmp_path):
    """Silently reshaping a short block would produce a picture with the
    rows walked one channel sideways, which looks like dispersion."""
    angle = np.linspace(-10.0, 10.0, 4)
    energy = np.linspace(20.0, 20.3, 3)
    path = write_spectrum(str(tmp_path / "short.txt"), energy, angle,
                          np.zeros((4, 3)))
    lines = open(path, encoding="latin-1").read().splitlines()
    index = lines.index("[Data 1]") + 1
    lines[index] = " ".join(lines[index].split()[:-1])      # one column short
    with open(path, "w", newline="\r\n", encoding="latin-1") as handle:
        handle.write("\n".join(lines))
    with pytest.raises(ValueError, match="not a whole number of rows"):
        cassiopee.parse_scienta(path)


# --------------------------------------------------------------------------
# The axis scales
# --------------------------------------------------------------------------
def test_rounding_jitter_is_straightened_without_complaint(tmp_path):
    angle = np.linspace(-15.0, 15.0, 6)
    energy = np.linspace(115.6, 116.0, 5)
    path = write_spectrum(str(tmp_path / "a.txt"), energy, angle,
                          np.zeros((6, 5)))
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        region = cassiopee.parse_scienta(path)
        steps = np.diff(region[0].energy)
    assert np.allclose(steps, steps[0])


def test_a_genuinely_non_linear_scale_is_reported(tmp_path):
    """Replacing a curved scale by a straight line moves the data onto the
    wrong energies. The MATLAB version does it silently; this says so."""
    angle = np.linspace(-15.0, 15.0, 4)
    energy = np.array([115.6, 115.7, 116.4, 116.5])      # not evenly stepped
    path = write_spectrum(str(tmp_path / "bent.txt"), energy, angle,
                          np.zeros((4, 4)))
    with pytest.warns(UserWarning, match="not linear"):
        cassiopee.parse_scienta(path)


def test_regularisation_can_be_turned_off(tmp_path):
    angle = np.linspace(-15.0, 15.0, 4)
    energy = np.array([115.6, 115.7, 116.4, 116.5])
    path = write_spectrum(str(tmp_path / "bent.txt"), energy, angle,
                          np.zeros((4, 4)))
    region = cassiopee.parse_scienta(path, regularise=False)[0]
    assert np.allclose(region.energy, energy)


# --------------------------------------------------------------------------
# The median filter
# --------------------------------------------------------------------------
def test_the_median_filter_is_off_unless_asked_for(spectrum):
    path, _energy, _angle, values = spectrum
    assert np.allclose(cassiopee.parse_scienta(path)[0].values, values)
    filtered = cassiopee.parse_scienta(path, median_filter=True)[0].values
    assert not np.allclose(filtered, values)


def test_the_median_filter_is_the_two_point_mean_along_both_axes():
    values = np.array([[0.0, 4.0], [8.0, 0.0]])
    out = cassiopee._median_filter_2(values)
    # axis 0 first: [[0, 4], [4, 2]]; then axis 1: [[0, 2], [4, 3]]
    assert np.allclose(out, [[0.0, 2.0], [4.0, 3.0]])


def test_the_first_row_is_not_halved():
    """MATLAB's medfilt1 zero-pads, so its first row and column come out at
    half value -- a dark line down two sides of every spectrum."""
    values = np.full((3, 3), 10.0)
    out = cassiopee._median_filter_2(values)
    assert np.allclose(out, 10.0)


def test_using_the_filter_is_recorded(spectrum):
    region = cassiopee.parse_scienta(spectrum[0], median_filter=True)[0]
    assert region.info["cassiopee.median_filter"] == 1


# --------------------------------------------------------------------------
# The parameter file
# --------------------------------------------------------------------------
def test_the_parameter_file_is_found_beside_the_spectrum(tmp_path):
    spectrum = str(tmp_path / "Map_3_ROI1_.txt")
    open(spectrum, "w").close()
    parameters = write_parameters(str(tmp_path / "Map_3_i.txt"))
    assert cassiopee.parameter_path(spectrum) == parameters


def test_a_spectrum_with_no_parameter_file_is_not_an_error(tmp_path):
    spectrum = str(tmp_path / "Map_3_ROI1_.txt")
    open(spectrum, "w").close()
    assert cassiopee.parameter_path(spectrum) is None
    assert cassiopee.read_parameter_file(None) == {}


@pytest.mark.parametrize("spacing", [" : ", ":", " :", ": ", "  :  "])
def test_the_parameter_file_is_split_on_the_colon_not_at_a_column(tmp_path,
                                                                 spacing):
    """The MATLAB version reads from a fixed character offset. A field that
    moved by one character would come back as a silent NaN."""
    path = write_parameters(str(tmp_path / "p.txt"), theta=12.5,
                            photon_energy=80.0, spacing=spacing)
    values = cassiopee.read_parameter_file(path)
    assert values["sample_theta_deg"] == pytest.approx(12.5)
    assert values["photon_energy_eV"] == pytest.approx(80.0)
    assert values["sample_X_mm"] == pytest.approx(1.234)
    assert values["temperature_K"] == pytest.approx(17.5)


def test_the_polarisation_code_is_named(tmp_path):
    path = write_parameters(str(tmp_path / "p.txt"))
    assert cassiopee.read_parameter_file(path)["polarisation"] == "LH"


# --------------------------------------------------------------------------
# The work-function table
# --------------------------------------------------------------------------
def test_the_work_function_passes_through_its_calibration_points():
    table = cassiopee.ANALYSER_WORK_FUNCTION
    for hv, phi in zip(table["photon_energy_eV"], table["work_function_eV"]):
        assert cassiopee.work_function(hv) == pytest.approx(phi, abs=1e-9)


def test_the_work_function_interpolates_between_them():
    low = cassiopee.work_function(45.0)
    high = cassiopee.work_function(75.0)
    middle = cassiopee.work_function(60.0)
    assert low < middle < high


def test_extrapolating_the_work_function_warns():
    """A cubic through a table that dips leaves its range fast, and a wrong
    work function shifts a whole photon-energy scan's binding energies while
    the picture still looks fine."""
    with pytest.warns(UserWarning, match="outside the calibrated"):
        cassiopee.work_function(300.0)
    with pytest.warns(UserWarning, match="outside the calibrated"):
        cassiopee.work_function(10.0)


# --------------------------------------------------------------------------
# Finding a series by name
# --------------------------------------------------------------------------
def test_a_series_is_found_even_when_the_base_has_underscores(tmp_path):
    """``Map80eV_Phi18_3_ROI1_`` -- the MATLAB version splits on ``_`` and
    takes the first piece, which gives ``Map80eV``. It has a commented-out
    hand-set base in it that says as much."""
    paths = make_series(str(tmp_path), base="Map80eV_Phi18", count=4)
    members = cassiopee.series_members(paths[2])
    assert [index for index, _ in members] == [1, 2, 3, 4]
    assert members[0][1] == paths[0]


def test_a_different_base_in_the_same_folder_is_not_swept_in(tmp_path):
    make_series(str(tmp_path), base="MapA", count=3)
    other = make_series(str(tmp_path), base="MapB", count=2)
    assert len(cassiopee.series_members(other[0])) == 2


def test_a_different_roi_is_a_different_series(tmp_path):
    make_series(str(tmp_path), base="Map", count=3, roi=1)
    second = make_series(str(tmp_path), base="Map", count=2, roi=2)
    assert len(cassiopee.series_members(second[0])) == 2


def test_the_members_come_back_in_numerical_not_alphabetical_order(tmp_path):
    paths = make_series(str(tmp_path), base="Map", count=11)
    members = cassiopee.series_members(paths[0])
    assert [index for index, _ in members] == list(range(1, 12))


def test_a_standalone_spectrum_belongs_to_no_series(tmp_path):
    angle = np.linspace(-1.0, 1.0, 3)
    energy = np.linspace(10.0, 10.2, 3)
    path = write_spectrum(str(tmp_path / "0001.txt"), energy, angle,
                          np.zeros((3, 3)))
    assert cassiopee.series_members(path) == []
    with pytest.raises(ValueError, match="not part of a numbered"):
        cassiopee.load_series(path)


# --------------------------------------------------------------------------
# Assembling a folder
# --------------------------------------------------------------------------
def test_a_stepped_polar_angle_makes_a_map(tmp_path):
    thetas = [-6.0, -3.0, 0.0, 3.0]
    paths = make_series(str(tmp_path), count=4, thetas=thetas)
    scan = cassiopee.load_series(paths[0])
    assert scan.kind == "map"
    assert scan.value.shape == (4, 5, 4)            # step, angle, energy
    assert np.allclose(scan.x, thetas)
    assert scan.info["axis0.role"] == "angle"
    assert scan.labels["x"] == "Polar angle theta (deg)"
    # Each member was written as a constant equal to its index.
    assert np.allclose(scan.value[2], 3.0)


def test_a_stepped_photon_energy_makes_an_hv_scan(tmp_path):
    energies = [80.0, 90.0, 100.0]
    paths = make_series(str(tmp_path), count=3, photon_energies=energies)
    scan = cassiopee.load_series(paths[0])
    assert np.allclose(scan.x, energies)
    assert scan.info["axis0.role"] == "photon_energy"
    assert scan.labels["x"] == "Photon energy (eV)"


def test_an_hv_scan_is_referred_to_the_fermi_level(tmp_path):
    """Each member was taken at a different hv, so the same kinetic energy
    is a different state; only E - E_F is comparable across the set."""
    energies = [80.0, 90.0, 100.0]
    paths = make_series(str(tmp_path), count=3, photon_energies=energies)
    scan = cassiopee.load_series(paths[0])
    phi = cassiopee.work_function(80.0)
    expected = np.linspace(115.6, 116.0, 4) - 80.0 + phi
    assert np.allclose(scan.z, expected)
    assert scan.labels["z"] == "E - E_F (eV)"
    assert scan.info["cassiopee.work_function_eV"] == pytest.approx(phi)


def test_a_series_that_stepped_neither_is_indexed(tmp_path):
    """Repeated acquisitions at one position. Calling the first axis an
    angle would offer a k conversion of a cut number."""
    paths = make_series(str(tmp_path), count=3)
    scan = cassiopee.load_series(paths[0])
    assert np.allclose(scan.x, [1, 2, 3])
    assert scan.info["axis0.role"] == "other"
    assert scan.labels["x"] == "Index (cut)"


def test_a_map_is_not_offered_a_k_conversion_when_it_stepped_photon_energy(
        tmp_path):
    paths = make_series(str(tmp_path), count=3,
                        photon_energies=[80.0, 90.0, 100.0])
    scan = cassiopee.load_series(paths[0])
    assert registry.role_is_angle(scan) is False


def test_a_theta_map_is_offered_one(tmp_path):
    paths = make_series(str(tmp_path), count=3, thetas=[-3.0, 0.0, 3.0])
    assert registry.role_is_angle(cassiopee.load_series(paths[0])) is True


def test_a_member_of_a_different_shape_is_an_error(tmp_path):
    """Stacking them anyway would raise deep inside numpy, or worse,
    broadcast."""
    paths = make_series(str(tmp_path), base="Map", count=2,
                        thetas=[0.0, 1.0])
    angle = np.linspace(-15.0, 15.0, 7)
    energy = np.linspace(115.6, 116.0, 4)
    write_spectrum(str(tmp_path / "Map_3_ROI1_.txt"), energy, angle,
                   np.zeros((7, 4)))
    write_parameters(str(tmp_path / "Map_3_i.txt"), theta=2.0)
    with pytest.raises(ValueError, match="cannot be one map"):
        cassiopee.load_series(paths[0])


def test_the_progress_callback_is_called_once_per_member(tmp_path):
    paths = make_series(str(tmp_path), count=5, thetas=[0, 1, 2, 3, 4])
    seen = []
    cassiopee.load_series(paths[0], progress=lambda d, t, l: seen.append((d, t)))
    assert [d for d, _ in seen] == [0, 1, 2, 3, 4]
    assert {t for _, t in seen} == {5}


# --------------------------------------------------------------------------
# Through the registry, as the program actually uses it
# --------------------------------------------------------------------------
def test_the_loader_is_registered():
    assert "SOLEIL CASSIOPEE" in [loader.name for loader in registry.loaders()]


def test_a_spectrum_is_detected(spectrum):
    found = registry.detect(spectrum[0])
    assert found is not None and found.name == "SOLEIL CASSIOPEE"


def test_a_parameter_file_is_not_claimed(tmp_path):
    """They are .txt in the same folder, and listing them would put rows in
    the browser that cannot be opened."""
    path = write_parameters(str(tmp_path / "Map_1_i.txt"))
    assert cassiopee.CassiopeeLoader().can_open(path) is False


def test_binary_is_not_claimed(tmp_path):
    """This loader is asked before the HDF5 readers, and latin-1 decodes
    arbitrary bytes without complaint."""
    path = str(tmp_path / "thing.nxs")
    with open(path, "wb") as handle:
        handle.write(b"\x89HDF\r\n\x1a\n" + b"Dimension 1 scale=\x00" * 40)
    assert cassiopee.CassiopeeLoader().can_open(path) is False


def test_a_missing_file_is_not_claimed(tmp_path):
    assert cassiopee.CassiopeeLoader().can_open(
        str(tmp_path / "nothing.txt")) is False


def test_a_lone_spectrum_lists_one_entry(spectrum):
    entries = registry.list_entries(spectrum[0])
    assert [e["entry"] for e in entries] == ["cut"]
    assert entries[0]["kind"] == "Cut"


def test_a_member_of_a_series_lists_both_itself_and_the_series(tmp_path):
    """The same file is a spectrum in its own right and a slice of a map,
    and which one is wanted depends on what is being looked at."""
    paths = make_series(str(tmp_path), base="Map80eV_Phi18", count=4,
                        thetas=[-3.0, -1.0, 1.0, 3.0])
    entries = registry.list_entries(paths[1])
    kinds = {e["entry"]: e["kind"] for e in entries}
    assert kinds == {"cut": "Cut", "series": "Map"}
    series = [e for e in entries if e["entry"] == "series"][0]
    assert "Map80eV_Phi18" in series["name"] and "4" in series["name"]


def test_loading_through_the_registry_gives_a_cut(spectrum):
    scan = registry.load(spectrum[0], "cut")
    assert scan.kind == "cut"
    assert scan.info["loader.name"] == "SOLEIL CASSIOPEE"


def test_loading_the_series_entry_gives_the_map(tmp_path):
    paths = make_series(str(tmp_path), count=4, thetas=[-3.0, -1.0, 1.0, 3.0])
    scan = registry.load(paths[0], "series")
    assert scan.kind == "map" and scan.value.shape == (4, 5, 4)


def test_the_readers_own_role_survives_a_default_load(tmp_path):
    """The dialog's "as the file says" must not relabel a photon-energy scan
    as an angle -- that would put a k conversion back on the menu for it."""
    paths = make_series(str(tmp_path), count=3,
                        photon_energies=[80.0, 90.0, 100.0])
    scan = registry.load(paths[0], "series", registry.LoadOptions())
    assert scan.info["axis0.role"] == "photon_energy"
    assert scan.labels["x"] == "Photon energy (eV)"


def test_the_user_can_still_override_the_role(tmp_path):
    """A reader's guess is a guess; a series that stepped temperature while
    theta happened to drift would be read as a map."""
    paths = make_series(str(tmp_path), count=3, thetas=[10.0, 11.0, 12.0])
    scan = registry.load(paths[0], "series",
                         registry.LoadOptions(axis0_role="temperature"))
    assert scan.info["axis0.role"] == "temperature"
    assert scan.labels["x"] == "Temperature (K)"
    assert registry.role_is_angle(scan) is False


def test_progress_reaches_the_loader_through_the_registry(tmp_path):
    paths = make_series(str(tmp_path), count=3, thetas=[0.0, 1.0, 2.0])
    seen = []
    registry.load(paths[0], "series", progress=lambda d, t, l: seen.append(d))
    assert seen == [0, 1, 2]


def test_a_loader_that_takes_no_progress_is_still_callable(tmp_path):
    """The argument was added to the protocol after the other loaders were
    written; passing it to one of those must not be a TypeError."""
    class Old(registry.Loader):
        name = "old-style"

        def can_open(self, path):
            return True

        def list_entries(self, path):
            return []

        def load(self, path, entry=None):
            return cassiopee.load_cut(path)

    old = Old()
    angle = np.linspace(-1.0, 1.0, 3)
    energy = np.linspace(10.0, 10.2, 3)
    path = write_spectrum(str(tmp_path / "x.txt"), energy, angle,
                          np.zeros((3, 3)))
    assert registry._accepts_progress(old.load) is False
    registry.register(old)
    try:
        scan = registry.load(path, None,
                             registry.LoadOptions(loader="old-style"),
                             progress=lambda *a: None)
    finally:
        registry._REGISTRY.pop("old-style", None)
    assert scan.kind == "cut"


def test_the_permutation_option_works_on_an_assembled_series(tmp_path):
    paths = make_series(str(tmp_path), count=4, thetas=[-3.0, -1.0, 1.0, 3.0])
    scan = registry.load(paths[0], "series",
                         registry.LoadOptions(permutation=(1, 0, 2)))
    assert scan.value.shape == (5, 4, 4)
