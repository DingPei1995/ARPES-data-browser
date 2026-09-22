"""Tests for the loader registry: that the program's own format is
recognised without being chosen, that a beamline loader can be added without
touching anything, and that the two load-time axis options (permutation and
what the scanned axis really is) do what they say to the axes, the labels
and the array together.
"""
import numpy as np
import pytest

import loader.registry as L
from loader.nxs_file import NxsScan
from loader.session import SessionStore


@pytest.fixture()
def store(tmp_path):
    return SessionStore(root=str(tmp_path))


def _map_scan(shape=(2, 3, 4)):
    scan = NxsScan(kind="map")
    scan.x = np.arange(float(shape[0]))
    scan.k = np.arange(float(shape[1]))
    scan.z = np.arange(float(shape[2]))
    scan.value = np.arange(int(np.prod(shape)), dtype=float).reshape(shape)
    scan.labels = {"x": "Angle (deflector) (deg)", "k": "ky (A-1)",
                   "z": "Energy (eV)"}
    return scan


def _stored_map(store, shape=(2, 3, 4)):
    scan = _map_scan(shape)
    return store.store({"name": "m", "kind": "map",
                        "axes": (scan.x, scan.k, scan.z), "labels": scan.labels,
                        "value": scan.value, "info": {}, "motors": {}})


# -- the registry -------------------------------------------------------------
def test_both_shipped_loaders_are_registered():
    names = [loader.name for loader in L.loaders()]
    assert "SOLEIL ANTARES" in names
    assert "This program's own format" in names


def test_the_native_format_outranks_the_beamline_guess():
    """Detection order matters: the native format is self-describing, while
    a beamline loader recognises a layout it might share with another."""
    ranked = L.loaders()
    assert ranked[0].name == "This program's own format"
    assert ranked[0].priority > max(l.priority for l in ranked[1:])


def test_a_saved_dataset_is_detected_without_being_told(store):
    loader = L.detect(_stored_map(store))
    assert loader is not None
    assert loader.name == "This program's own format"


def test_detect_returns_none_for_something_unrecognisable(tmp_path):
    path = tmp_path / "not-a-scan.nxs"
    path.write_bytes(b"this is not HDF5 at all")
    assert L.detect(str(path)) is None


def test_load_without_a_loader_says_so_usefully(tmp_path):
    path = tmp_path / "junk.nxs"
    path.write_bytes(b"nope")
    with pytest.raises(ValueError, match="no loader recognises"):
        L.load(str(path))


def test_a_broken_loader_does_not_stop_the_others(store, monkeypatch):
    class Exploding(L.Loader):
        name = "exploding test loader"
        priority = 1000

        def can_open(self, path):
            raise RuntimeError("this loader is broken")

    broken = Exploding()
    L.register(broken)
    try:
        assert L.detect(_stored_map(store)).name == "This program's own format"
    finally:
        L._REGISTRY.pop(broken.name, None)


def test_get_loader_names_what_is_available():
    with pytest.raises(ValueError, match="no loader called"):
        L.get_loader("Diamond I05")


def test_adding_a_beamline_is_one_class_and_one_register(store):
    """The claim the whole module exists to support, exercised: a loader
    defined here, registered here, is used by the ordinary load path."""
    class FakeBeamline(L.Loader):
        name = "test beamline"
        priority = 50

        def can_open(self, path):
            return path.endswith(".fake")

        def list_entries(self, path):
            return [{"entry": "e1", "kind": "Map", "title": None,
                     "start_time": None}]

        def load(self, path, entry=None):
            return _map_scan()

    L.register(FakeBeamline())
    try:
        scan = L.load("whatever.fake")
        assert scan.kind == "map"
        assert scan.info["loader.name"] == "test beamline"
        assert L.list_entries("whatever.fake")[0]["entry"] == "e1"
    finally:
        L._REGISTRY.pop("test beamline", None)


# -- permuting the axes ---------------------------------------------------------
def test_permutation_moves_array_axes_and_their_vectors_together():
    scan = L.permute_axes(_map_scan((2, 3, 4)), (1, 0, 2))
    assert scan.value.shape == (3, 2, 4)
    assert scan.x.shape == (3,)          # what used to be k
    assert scan.k.shape == (2,)          # what used to be x
    assert scan.z.shape == (4,)


def test_permutation_takes_the_labels_with_it():
    """The part that is easy to get wrong by hand: an axis that moves must
    keep its own name and unit, or every later plot is mislabelled."""
    scan = L.permute_axes(_map_scan(), (1, 0, 2))
    assert scan.labels["x"] == "ky (A-1)"
    assert scan.labels["k"] == "Angle (deflector) (deg)"
    assert scan.labels["z"] == "Energy (eV)"


def test_permutation_matches_numpy_transpose():
    original = _map_scan((2, 3, 4)).value.copy()
    scan = L.permute_axes(_map_scan((2, 3, 4)), (2, 0, 1))
    assert np.array_equal(np.asarray(scan.value), np.transpose(original, (2, 0, 1)))


def test_the_identity_permutation_changes_nothing():
    scan = L.permute_axes(_map_scan(), (0, 1, 2))
    assert scan.value.shape == (2, 3, 4)
    assert scan.labels["x"] == "Angle (deflector) (deg)"


def test_a_permutation_that_is_not_an_ordering_is_refused():
    for bad in [(0, 0, 1), (0, 1), (0, 1, 3)]:
        with pytest.raises(ValueError, match="not an ordering"):
            L.permute_axes(_map_scan(), bad)


def test_permutation_works_through_a_lazy_array(store):
    """A stored dataset comes back lazy; asking for a different axis order
    has to read it, and must still be correct."""
    path = _stored_map(store, (2, 3, 4))
    expected = np.transpose(np.arange(24, dtype=float).reshape(2, 3, 4), (1, 0, 2))
    scan = L.load(path, options=L.LoadOptions(permutation=(1, 0, 2)))
    assert np.allclose(np.asarray(scan.value), expected)


# -- what the scanned axis is ------------------------------------------------------
def test_the_first_axis_is_an_angle_unless_told_otherwise():
    scan = L.apply_options(_map_scan(), L.LoadOptions())
    assert scan.info["axis0.role"] == "angle"
    assert L.role_is_angle(scan) is True
    assert scan.labels["x"] == "Angle (deflector) (deg)"


@pytest.mark.parametrize("role,expected", [
    ("photon_energy", "Photon energy (eV)"),
    ("temperature", "Temperature (K)"),
    ("gate_voltage", "Gate voltage (V)"),
    ("delay", "Delay (ps)"),
])
def test_naming_the_scanned_axis_relabels_it(role, expected):
    scan = L.apply_options(_map_scan(), L.LoadOptions(axis0_role=role))
    assert scan.labels["x"] == expected
    assert scan.info["axis0.role"] == role


def test_a_non_angle_axis_is_not_offered_to_the_k_conversion():
    """Converting a temperature axis to momentum would produce a plausible
    and meaningless picture; this is the flag that lets the viewer refuse."""
    for role in ("photon_energy", "temperature", "gate_voltage", "other"):
        scan = L.apply_options(_map_scan(), L.LoadOptions(axis0_role=role))
        assert L.role_is_angle(scan) is False


def test_a_custom_label_beats_the_role_default():
    scan = L.apply_options(
        _map_scan(), L.LoadOptions(axis0_role="other", axis0_label="Bias (mV)"))
    assert scan.labels["x"] == "Bias (mV)"


def test_a_scan_read_before_roles_existed_counts_as_an_angle():
    """Every file read before this feature existed was an angle scan, and
    must keep behaving as one."""
    scan = _map_scan()
    scan.info.pop("axis0.role", None)
    assert L.role_is_angle(scan) is True
