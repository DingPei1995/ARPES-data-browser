"""Tests for loader/session.py and the native storage format it writes: that a
derived dataset survives being written and read back exactly, that reading
one back costs no memory until it is used, and that the memory budget
evicts by least-recent use.
"""
import os

import numpy as np
import pytest

from loader.nxs_file import (save_dataset, load_soleil_nxs, LazyArray,
                            NATIVE_FORMAT_VERSION)
from loader.session import SessionStore, MemoryBudget, scan_to_dict


def _arpes_cube(shape=(20, 30, 16), dtype=np.float32):
    """Smooth band plus counting noise: compressible the way real data is,
    unlike uniform random numbers."""
    nx, nk, ne = shape
    x = np.linspace(-1, 1, nx)[:, None, None]
    k = np.linspace(-1, 1, nk)[None, :, None]
    e = np.linspace(0, 1, ne)[None, None, :]
    rate = 400.0 * np.exp(-((e - 0.5 - 0.3 * (x ** 2 + k ** 2)) ** 2) / 0.01) + 20.0
    return np.random.default_rng(0).poisson(np.broadcast_to(rate, shape)).astype(dtype)


def _map_item(name="derived", cube=None):
    cube = _arpes_cube() if cube is None else cube
    nx, nk, ne = cube.shape
    return {"name": name, "kind": "map",
            "axes": (np.linspace(-1, 1, nx), np.linspace(-2, 2, nk),
                     np.linspace(20, 21, ne)),
            "labels": {"x": "kx (A-1)", "k": "ky (A-1)", "z": "Energy (eV)"},
            "value": cube, "info": {"kconv.azimuth_deg": 12.5}, "motors": {}}


@pytest.fixture()
def store(tmp_path):
    return SessionStore(root=str(tmp_path))


# -- the format ---------------------------------------------------------------
@pytest.mark.parametrize("dtype", [np.float64, np.float32, np.uint16])
def test_round_trip_is_exact_and_keeps_the_dtype(store, dtype):
    """Version 2 writes the array in its own dtype. Version 1 wrote whatever
    the caller had promoted to float64, which doubled the file for data that
    arrived as float32 and gained nothing."""
    cube = _arpes_cube(dtype=dtype)
    scan = store.load_scan(store.store(_map_item(cube=cube)))
    assert scan.value.dtype == cube.dtype
    assert np.array_equal(np.asarray(scan.value), cube)


def test_stored_file_is_smaller_than_the_raw_array(store):
    cube = _arpes_cube((40, 60, 32))
    path = store.store(_map_item(cube=cube))
    assert os.path.getsize(path) < 0.7 * cube.nbytes


def test_axes_labels_and_settings_survive(store):
    scan = store.load_scan(store.store(_map_item()))
    assert scan.kind == "map"
    assert scan.labels["x"] == "kx (A-1)"
    assert scan.labels["z"] == "Energy (eV)"
    assert scan.info["kconv.azimuth_deg"] == pytest.approx(12.5)
    assert scan.x.shape == (20,) and scan.k.shape == (30,) and scan.z.shape == (16,)


def test_value_comes_back_lazy_not_resident(store):
    """The point of storing derived data: a listed dataset costs a path, not
    a cube, until something actually asks for the numbers."""
    scan = store.load_scan(store.store(_map_item()))
    assert isinstance(scan.value, LazyArray)
    frame = scan.value[:, :, 3]                 # one slice, not the cube
    assert isinstance(frame, np.ndarray)
    assert frame.shape == (20, 30)


def test_lazy_value_survives_being_handed_to_numpy(store):
    """Anything that genuinely needs every point -- an operation, a fit, a
    save -- must keep working on a lazy array without knowing it is one."""
    cube = _arpes_cube()
    scan = store.load_scan(store.store(_map_item(cube=cube)))
    assert float(np.asarray(scan.value).sum()) == pytest.approx(float(cube.sum()))
    assert np.asarray(scan.value, dtype=np.float64).dtype == np.float64


def test_a_version_1_file_still_reads(tmp_path):
    """Files written before the format changed must keep opening: uncompressed,
    float64, no shuffle -- the reader must not depend on any of that."""
    import h5py
    from loader.nxs_file import NATIVE_ATTR, NATIVE_FORMAT_ATTR
    path = str(tmp_path / "old.nxs")
    cube = _arpes_cube(dtype=np.float64)
    with h5py.File(path, "w") as f:
        grp = f.create_group("old map")
        grp.attrs[NATIVE_ATTR] = "map"
        grp.attrs[NATIVE_FORMAT_ATTR] = 1
        grp.attrs["nxsloader_name"] = "old map"
        for name, values in zip(("axis_x", "axis_y", "axis_z"),
                                (np.arange(20.), np.arange(30.), np.arange(16.))):
            grp.create_dataset(name, data=values)
        grp.create_dataset("value", data=cube, compression="gzip", compression_opts=4)
        grp.create_group("info")
        grp.create_group("motors")
    scan = load_soleil_nxs(path)
    assert scan.kind == "map"
    assert np.array_equal(np.asarray(scan.value), cube)


def test_format_version_is_recorded(store):
    import h5py
    path = store.store(_map_item())
    with h5py.File(path, "r") as f:
        grp = f[list(f.keys())[0]]
        assert int(grp.attrs["nxsloader_format"]) == NATIVE_FORMAT_VERSION
        assert grp["value"].compression == "gzip"
        assert grp["value"].shuffle is True


# -- the store ------------------------------------------------------------------
def test_two_datasets_of_the_same_name_do_not_collide(store):
    first = store.store(_map_item("same name"))
    second = store.store(_map_item("same name"))
    assert first != second
    assert os.path.isfile(first) and os.path.isfile(second)


def test_discard_removes_the_file(store):
    path = store.store(_map_item())
    store.discard(path)
    assert not os.path.exists(path)
    store.discard(path)          # and again is harmless


def test_prune_old_leaves_this_session_alone(tmp_path):
    import time
    old = tmp_path / "20200101-000000-1"
    old.mkdir()
    os.utime(old, (time.time() - 30 * 86400,) * 2)
    store = SessionStore(root=str(tmp_path))
    store.store(_map_item())                     # so this session has a folder
    assert store.prune_old(keep_days=7) == 1
    assert not old.exists()
    assert os.path.isdir(store.folder)


def test_bytes_on_disk_counts_what_was_written(store):
    assert store.bytes_on_disk() == 0
    path = store.store(_map_item())
    assert store.bytes_on_disk() == os.path.getsize(path)


def test_scan_to_dict_round_trips_a_cut(store):
    from loader.nxs_file import NxsScan
    scan = NxsScan(kind="cut")
    scan.x, scan.y = np.arange(10.0), np.arange(12.0)
    scan.value = np.random.default_rng(1).random((10, 12))
    scan.labels = {"x": "k", "y": "E"}
    item = scan_to_dict("cut", scan, "a cut")
    back = store.load_scan(store.store(item))
    assert back.kind == "cut"
    assert np.allclose(np.asarray(back.value), scan.value)
    assert back.y.shape == (12,)


def test_scan_to_dict_refuses_a_kind_it_cannot_describe():
    from loader.nxs_file import NxsScan
    with pytest.raises(ValueError, match="plain axes"):
        scan_to_dict("unsupported", NxsScan(kind="unsupported"), "x")


# -- the memory budget -----------------------------------------------------------
class _Fake:
    def __init__(self, nbytes):
        self.scan = type("S", (), {"value": np.zeros(nbytes // 8)})()


def test_budget_evicts_least_recently_used_first():
    budget = MemoryBudget(limit_bytes=2400)
    for name in ("a", "b", "c"):
        budget.put(name, _Fake(800))
    budget.get("a")                       # a is now the most recent
    budget.put("d", _Fake(800))           # over budget -> b goes, not a
    assert budget.get("b") is None
    assert budget.get("a") is not None
    assert budget.get("d") is not None


def test_budget_never_evicts_the_only_entry():
    """A single dataset larger than the whole budget still has to be usable;
    the budget governs how much is kept on spec, not what may exist."""
    budget = MemoryBudget(limit_bytes=100)
    budget.put("huge", _Fake(80000))
    assert budget.get("huge") is not None


def test_budget_discard_and_clear():
    budget = MemoryBudget()
    budget.put("a", _Fake(80))
    budget.discard("a")
    assert budget.get("a") is None
    budget.put("b", _Fake(80))
    budget.clear()
    assert budget.total_bytes() == 0


# -- the operations log ---------------------------------------------------------
def test_every_store_is_logged(store):
    from loader.session import read_log
    path = store.store(_map_item("logged one"))
    rows = read_log(store.root)
    assert any(row[2].strip() == "STORE" and row[4] == path for row in rows), rows


def test_saving_a_dataset_logs_it_and_removes_the_working_copy(store):
    from loader.session import read_log
    path = store.store(_map_item("done with"))
    store.released(path, "/home/me/results.nxs::done with")
    assert not os.path.exists(path)
    rows = read_log(store.root)
    saved = [row for row in rows if row[2].strip() == "SAVED"]
    assert saved and "results.nxs" in saved[-1][4]


def test_discard_says_why_in_the_log(store):
    from loader.session import read_log
    path = store.store(_map_item())
    store.discard(path, reason="removed from the list")
    rows = [row for row in read_log(store.root) if row[2].strip() == "DISCARD"]
    assert rows and "removed from the list" in rows[-1][5]


def test_the_log_survives_a_corrupted_tail(store):
    """A crash mid-write leaves a half line. The log is consulted *after* a
    crash, so it has to read back anyway."""
    from loader.session import read_log
    store.store(_map_item("good"))
    with open(store.log_path, "a", encoding="utf-8") as handle:
        handle.write("this line was cut off mid-wr")
    rows = read_log(store.root)
    assert any(row[2].strip() == "STORE" for row in rows)


def test_the_log_is_plain_text_anyone_can_open(store):
    store.store(_map_item("readable"))
    text = open(store.log_path, encoding="utf-8").read()
    assert "STORE" in text and "readable" in text
    assert text.count(" | ") >= 5


def test_retention_is_three_days():
    """Long enough to come back on Monday for Friday's crash is not the
    trade any more: a copy that matters is saved, and one that is still here
    after three days is one nobody came back for."""
    from loader.session import KEEP_DAYS
    assert KEEP_DAYS == 3


# -- finding what a crash left behind -------------------------------------------
def test_leftovers_lists_another_sessions_datasets(tmp_path):
    from loader.session import leftover_sessions
    crashed = SessionStore(root=str(tmp_path))
    crashed.store(_map_item("half-finished"))
    fresh = SessionStore(root=str(tmp_path))
    found = fresh.leftovers()
    assert len(found) == 1
    assert found[0]["folder"] == crashed.folder
    assert len(found[0]["files"]) == 1
    assert found[0]["bytes"] > 0
    assert leftover_sessions(str(tmp_path)) == found


def test_a_session_that_saved_everything_leaves_nothing_to_recover(tmp_path):
    tidy = SessionStore(root=str(tmp_path))
    path = tidy.store(_map_item())
    tidy.released(path, "/home/me/results.nxs::m")
    assert SessionStore(root=str(tmp_path)).leftovers() == []


def test_leftovers_ignores_this_session(tmp_path):
    store = SessionStore(root=str(tmp_path))
    store.store(_map_item())
    assert store.leftovers() == []


def test_leftovers_are_newest_first(tmp_path):
    import time as _time
    older = SessionStore(root=str(tmp_path))
    older.store(_map_item())
    _time.sleep(0.01)
    newer = SessionStore(root=str(tmp_path))
    assert newer.folder != older.folder      # same second, still distinct
    newer.store(_map_item())
    os.utime(newer.folder, (_time.time() + 10,) * 2)
    found = SessionStore(root=str(tmp_path)).leftovers()
    assert [entry["folder"] for entry in found][0] == newer.folder


def test_pruning_records_what_it_removed(tmp_path):
    import time as _time
    from loader.session import read_log
    old = SessionStore(root=str(tmp_path))
    old.store(_map_item())
    os.utime(old.folder, (_time.time() - 10 * 86400,) * 2)
    fresh = SessionStore(root=str(tmp_path))
    assert fresh.prune_old(keep_days=3) == 1
    rows = [row for row in read_log(str(tmp_path)) if row[2].strip() == "PRUNE"]
    assert rows and old.name in rows[-1][3]


# -- one table for "what axes does this kind have" -------------------------------
def test_a_four_dimensional_spatial_scan_keeps_all_four_axes(store):
    """Regression: the launcher's save path wrote a spem_4d with only three
    of its four axes -- (x, y, k), dropping energy -- because it carried its
    own copy of the axis table and that copy was wrong. Everything now reads
    loader.nxs_file.AXIS_SLOTS, so there is one place for this to be right.
    """
    from loader.nxs_file import NxsScan
    scan = NxsScan(kind="spem_4d")
    scan.x = np.linspace(0, 1, 4)        # spatial x
    scan.y = np.linspace(0, 2, 5)        # spatial y
    scan.k = np.linspace(-1, 1, 6)       # analyser angle
    scan.z = np.linspace(20, 21, 7)      # energy
    scan.value4d = np.random.default_rng(0).random((5, 4, 6, 7))
    scan.value = scan.value4d
    scan.labels = {"x": "x (um)", "y": "y (um)", "k": "k", "z": "E"}

    item = scan_to_dict("spem_4d", scan, "a spatial scan")
    assert len(item["axes"]) == 4
    back = store.load_scan(store.store(item))
    assert back.kind == "spem_4d"
    assert back.x.shape == (4,) and back.y.shape == (5,)
    assert back.k.shape == (6,)
    assert back.z is not None and back.z.shape == (7,), "the energy axis was lost"
    assert np.allclose(back.z, scan.z)


@pytest.mark.parametrize("kind,n_axes", [
    ("cut", 2), ("map", 3), ("k_map", 3), ("kz_map", 3), ("kz_map_k", 3),
    ("spem_1d", 3), ("spem_4d", 4),
])
def test_every_kind_has_one_axis_slot_per_dimension(kind, n_axes):
    from loader.nxs_file import axis_slots
    assert len(axis_slots(kind, "array")) == n_axes
    assert len(axis_slots(kind, "constructor")) == n_axes
    # the two orders are permutations of each other, never different sets
    assert set(axis_slots(kind, "array")) == set(axis_slots(kind, "constructor"))


def test_every_kind_has_an_energy_axis():
    """``energy_slot`` is derived from the axis table rather than written
    out, so a kind added to the table gets one for free -- which is the
    point, since the two panels that used to carry their own copy would
    otherwise refuse to shift a new kind's energy axis."""
    from loader.nxs_file import AXIS_SLOTS, energy_slot
    for kind in AXIS_SLOTS:
        assert energy_slot(kind) is not None, kind
    assert energy_slot("cut") == "y"          # a cut is (angle, energy)
    assert energy_slot("kz_map") == "z"
    assert energy_slot("unsupported") is None


def test_the_cube_kinds_are_exactly_the_three_axis_cubes():
    """CUBE_KINDS is what the viewers, the 3-D view and the processing
    panels branch on. It has to stay in step with the axis table, or a kind
    gets a window and no processing, or the reverse."""
    from loader.nxs_file import AXIS_SLOTS, CUBE_KINDS
    for kind in CUBE_KINDS:
        assert AXIS_SLOTS[kind]["array"] == ("x", "k", "z"), kind
    assert set(CUBE_KINDS) == {"map", "k_map", "kz_map", "kz_map_k"}


def test_every_kind_has_a_name_a_person_can_read():
    from loader.nxs_file import AXIS_SLOTS, KIND_LABELS
    for kind in AXIS_SLOTS:
        assert KIND_LABELS.get(kind), kind


def test_the_derived_tables_match_the_authority():
    from loader.nxs_file import AXIS_SLOTS
    from tools.dataops import ARRAY_AXES, CONSTRUCTOR_AXES
    from loader.registry import _axis_slots
    assert ARRAY_AXES == {k: v["array"] for k, v in AXIS_SLOTS.items()}
    assert CONSTRUCTOR_AXES == {k: v["constructor"] for k, v in AXIS_SLOTS.items()}
    for kind in AXIS_SLOTS:
        assert _axis_slots(kind) == AXIS_SLOTS[kind]["array"]


def test_memory_data_fills_the_slots_the_table_names():
    from ui.widgets import MemoryData
    data = MemoryData("spem_4d",
                      (np.arange(4.), np.arange(5.), np.arange(6.), np.arange(7.)),
                      np.zeros((5, 4, 6, 7)), {}, source_label="s")
    assert data.scan.x.shape == (4,)
    assert data.scan.y.shape == (5,)
    assert data.scan.k.shape == (6,)
    assert data.scan.z.shape == (7,)


def test_memory_data_rejects_the_wrong_number_of_axes():
    from ui.widgets import MemoryData
    with pytest.raises(ValueError, match="takes 3 axes"):
        MemoryData("map", (np.arange(4.), np.arange(5.)), np.zeros((4, 5)), {},
                   source_label="s")


def test_a_session_folder_from_before_the_rename_is_still_found(tmp_path, monkeypatch):
    """The program used to be called NxsLoader and kept its sessions
    elsewhere. Upgrading must not be a way to lose a crashed session's work,
    so the old location is still looked in."""
    from loader import session as S
    legacy_root = tmp_path / "old" / "sessions"
    legacy_root.mkdir(parents=True)
    old = SessionStore(root=str(legacy_root))
    old.store(_map_item("work from before the rename"))

    new_root = tmp_path / "new" / "sessions"
    monkeypatch.setattr(S, "DEFAULT_ROOT", str(new_root))
    monkeypatch.setattr(S, "LEGACY_ROOTS", (str(legacy_root),))
    fresh = SessionStore(root=str(new_root))
    found = fresh.leftovers()
    assert [entry["folder"] for entry in found] == [old.folder]


def test_a_file_saved_by_the_old_program_still_opens(tmp_path):
    """The marker attributes inside a saved file keep their original names
    on purpose: renaming them would make every file anyone has already saved
    unrecognisable."""
    from loader.nxs_file import NATIVE_ATTR, NATIVE_FORMAT_ATTR
    assert NATIVE_ATTR == "nxsloader_kind"
    assert NATIVE_FORMAT_ATTR == "nxsloader_format"
    store = SessionStore(root=str(tmp_path))
    path = store.store(_map_item("written now"))
    from loader import registry
    assert registry.detect(path).name == "This program's own format"
