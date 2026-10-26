"""Tests for loader.cubecache.CubeCache: the map cube that is read into
memory once browsing starts, instead of from the file on every slice.

The contract is the lazy arrays' own -- indexing is indistinguishable from
indexing the equivalent numpy array -- whether the cube is still on disk,
held whole in memory, or (too big for the budget) served from blocks.
"""
import numpy as np
import pytest
import h5py

from loader.cubecache import CubeCache, simplify_key
from loader.nxs_file import LazyCube, LazyArray


@pytest.fixture()
def dataset(tmp_path):
    array = np.arange(6 * 7 * 9, dtype=np.float32).reshape(6, 7, 9)
    path = tmp_path / "cube.h5"
    with h5py.File(path, "w") as f:
        f.create_dataset("v", data=array)
    handle = h5py.File(path, "r")
    yield handle["v"], array
    handle.close()


class _Counting:
    """A lazy array that counts how often it is read."""

    def __init__(self, array):
        self.array = array
        self.shape, self.dtype = array.shape, array.dtype
        self.reads = 0

    def __getitem__(self, key):
        self.reads += 1
        return self.array[key].copy()


RUN = np.array([2, 3, 4])
PATTERNS = {
    "contour slab": lambda a: a[:, :, RUN],
    "deflector cut": lambda a: a[:, RUN, :],
    "slit cut": lambda a: a[RUN, :, :],
    "one plane": lambda a: a[:, :, 5],
    "negative plane": lambda a: a[:, :, -1],
    "leading index": lambda a: a[2],
    "ellipsis": lambda a: a[..., 3],
    "scattered indices": lambda a: a[:, :, np.array([0, 3, 5])],
    "two axes": lambda a: a[1:3, :, 2],
    "everything": lambda a: a[:],
}


@pytest.mark.parametrize("name", list(PATTERNS))
@pytest.mark.parametrize("budget", [None, 1])
@pytest.mark.parametrize("permutation", [(0, 1, 2), (2, 0, 1)])
def test_indexes_like_the_array_it_stands_for(dataset, name, budget, permutation):
    """On disk (first read), in memory, and -- with a budget too small for
    the cube -- from blocks: always the same numbers."""
    dset, array = dataset
    reference = np.transpose(array, permutation)
    cached = CubeCache(LazyCube(dset, permutation),
                       budget=None if budget is None else reference.nbytes // 3)
    fn = PATTERNS[name]
    for _ in range(3):
        assert np.array_equal(fn(cached), fn(reference))


def test_first_read_is_lazy_then_the_cube_is_read_once():
    source = _Counting(np.random.default_rng(0).random((5, 6, 40)))
    cached = CubeCache(source, budget=10 ** 9)
    cached[:, :, 3]
    assert source.reads == 1 and not cached.loaded       # opening stays cheap
    assert cached.load_pending()
    for e in range(40):
        np.testing.assert_array_equal(cached[:, :, e], source.array[:, :, e])
    assert source.reads == 2 and cached.loaded            # one whole read


def test_a_cube_over_budget_is_never_read_whole():
    array = np.random.default_rng(1).random((4, 5, 64))
    source = _Counting(array)
    cached = CubeCache(source, budget=array.nbytes // 2)
    assert not cached.fits and not cached.load_pending()
    for e in range(64):
        np.testing.assert_array_equal(cached[:, :, e], array[:, :, e])
        np.testing.assert_array_equal(cached[:, :, e:e + 3].mean(axis=2),
                                      array[:, :, e:e + 3].mean(axis=2))
    assert not cached.loaded
    # Blocks of neighbouring planes: far fewer reads than slices asked for.
    assert source.reads < 40
    held = sum(block.nbytes for block in cached._blocks.values())
    assert held <= cached.budget


def test_slices_handed_out_cannot_change_the_cache():
    cached = CubeCache(_Counting(np.ones((3, 4, 5))), budget=10 ** 9)
    cached[0]
    view = cached[:, :, 1]
    with pytest.raises(ValueError):
        view[0, 0] = 5.0
    whole = np.asarray(cached)                 # a copy, and writable
    whole[:] = 7.0
    assert cached[0, 0, 0] == 1.0


def test_an_array_in_memory_is_used_as_it_is():
    array = np.arange(24.0).reshape(2, 3, 4)
    cached = CubeCache(array)
    assert cached.loaded
    assert np.shares_memory(cached[:, :, 1:3], array)


def test_lazy_array_and_dtype(dataset):
    dset, array = dataset
    cached = CubeCache(LazyArray(dset))
    assert cached.shape == array.shape and cached.dtype == array.dtype
    assert len(cached) == array.shape[0]
    assert np.asarray(cached, dtype=float).dtype == np.float64
    np.testing.assert_array_equal(np.asarray(cached), array)


def test_contiguous_runs_become_slices():
    key = simplify_key((slice(None), np.array([4, 5, 6]), [1]), 3)
    assert key == (slice(None), slice(4, 7), slice(1, 2))
    kept = simplify_key((np.array([1, 3]),), 3)[0]
    assert isinstance(kept, np.ndarray)
    assert simplify_key((Ellipsis, 2), 3) == (slice(None), slice(None), 2)
