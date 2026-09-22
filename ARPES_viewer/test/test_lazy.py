"""Tests for the lazy arrays that keep a measurement on disk: LazyCube (an
axis permutation applied on the way out) and LazyArray (no permutation, any
rank).

The contract is that indexing them is indistinguishable from indexing the
equivalent numpy array -- so every case below compares against exactly that.
The GUI's own access patterns are in here by name, because a map is now read
this way and those are the reads that happen on every cursor move.
"""
import numpy as np
import pytest
import h5py

from loader.nxs_file import LazyCube, LazyArray


@pytest.fixture()
def dataset(tmp_path):
    array = np.arange(6 * 7 * 5, dtype=np.float32).reshape(6, 7, 5)
    path = tmp_path / "lazy.h5"
    with h5py.File(path, "w") as f:
        f.create_dataset("v", data=array)
    handle = h5py.File(path, "r")
    yield handle["v"], array
    handle.close()


PATTERNS = {
    # what the constant-energy contour asks for on every energy change
    "contour slab [:, :, idxs]": lambda a, i: a[:, :, i],
    # what the two orthogonal cuts ask for
    "deflector cut [:, idxs, :]": lambda a, i: a[:, i, :],
    "slit cut [idxs, :, :]": lambda a, i: a[i, :, :],
    "single plane [:, :, 7]": lambda a, i: a[:, :, 2],
    "leading index [2]": lambda a, i: a[2],
    "ellipsis at the end [..., 3]": lambda a, i: a[..., 3],
    "ellipsis in the middle": lambda a, i: a[1, ..., 2],
    "plain slice [2:5]": lambda a, i: a[2:5],
    "everything [:]": lambda a, i: a[:],
}


@pytest.mark.parametrize("name", list(PATTERNS))
@pytest.mark.parametrize("permutation", [(0, 1, 2), (1, 0, 2), (2, 0, 1)])
def test_lazy_cube_indexes_like_the_array_it_stands_for(dataset, name, permutation):
    dset, array = dataset
    lazy = LazyCube(dset, permutation)
    reference = np.transpose(array, permutation)
    index = np.array([1, 2, 3])
    fn = PATTERNS[name]
    assert np.array_equal(fn(lazy, index), fn(reference, index))


def test_an_array_index_does_not_confuse_the_ellipsis_check(dataset):
    """Regression: `Ellipsis in key` compares element-wise against an array
    index and raises "truth value of an array ... is ambiguous". Every one
    of the GUI's slab reads passes an array index, so this broke the moment
    maps started being read lazily."""
    dset, array = dataset
    lazy = LazyCube(dset, (0, 1, 2))
    index = np.array([0, 2, 4])
    assert np.array_equal(lazy[:, :, index], array[:, :, index])


def test_lazy_cube_reports_its_shape_in_its_own_axis_order(dataset):
    dset, array = dataset
    assert LazyCube(dset, (0, 1, 2)).shape == (6, 7, 5)
    assert LazyCube(dset, (1, 0, 2)).shape == (7, 6, 5)
    assert LazyCube(dset, (2, 1, 0)).shape == (5, 7, 6)


def test_lazy_cube_materialises_in_its_own_axis_order(dataset):
    dset, array = dataset
    lazy = LazyCube(dset, (1, 0, 2))
    assert np.array_equal(lazy.materialise(), np.transpose(array, (1, 0, 2)))
    assert np.array_equal(np.asarray(lazy), np.transpose(array, (1, 0, 2)))


def test_lazy_cube_reports_size_and_footprint(dataset):
    dset, array = dataset
    lazy = LazyCube(dset, (0, 1, 2))
    assert lazy.size == array.size
    assert lazy.nbytes == array.nbytes      # what it *would* cost, resident


# -- LazyArray ----------------------------------------------------------------
@pytest.mark.parametrize("shape", [(11,), (6, 7), (6, 7, 5), (3, 4, 5, 2)])
def test_lazy_array_works_at_any_rank(tmp_path, shape):
    array = np.arange(int(np.prod(shape)), dtype=np.float32).reshape(shape)
    path = tmp_path / "a.h5"
    with h5py.File(path, "w") as f:
        f.create_dataset("v", data=array)
    with h5py.File(path, "r") as f:
        lazy = LazyArray(f["v"])
        assert lazy.shape == shape
        assert lazy.ndim == len(shape)
        assert lazy.dtype == array.dtype
        assert np.array_equal(np.asarray(lazy), array)
        assert np.array_equal(lazy[0], array[0])


def test_lazy_array_passes_dtype_through_asarray(tmp_path):
    array = np.arange(10, dtype=np.uint16)
    path = tmp_path / "a.h5"
    with h5py.File(path, "w") as f:
        f.create_dataset("v", data=array)
    with h5py.File(path, "r") as f:
        lazy = LazyArray(f["v"])
        assert np.asarray(lazy).dtype == np.uint16
        assert np.asarray(lazy, dtype=np.float64).dtype == np.float64


def test_numpy_functions_work_on_a_lazy_array(tmp_path):
    """np.sum(lazy) and friends go through __array__; lazy.sum() does not
    exist and is not meant to -- callers that want the whole array say so."""
    array = np.arange(24, dtype=np.float32).reshape(4, 6)
    path = tmp_path / "a.h5"
    with h5py.File(path, "w") as f:
        f.create_dataset("v", data=array)
    with h5py.File(path, "r") as f:
        lazy = LazyArray(f["v"])
        assert float(np.sum(lazy)) == pytest.approx(float(array.sum()))
        assert np.allclose(np.mean(lazy, axis=0), array.mean(axis=0))
        assert not hasattr(lazy, "sum")
