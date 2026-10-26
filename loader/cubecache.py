"""
loader/cubecache.py
===================
A map cube that is read from its file once while it is being browsed,
instead of once per slice.

Why this exists
---------------
A map is opened lazily (``loader.nxs_file.LazyCube``/``LazyArray``): the
numbers stay in the file and each slice is read as it is asked for. That is
right for *opening* a map, and it is right for the orthogonal cuts, which are
contiguous runs on disk. It is exactly wrong for the constant-energy contour.

Energy is the fastest-varying axis in the file, so one constant-energy plane
is one value out of every ``len(E)`` -- every page of the dataset, or every
compressed chunk, has to be read (and decompressed) to collect it. One tick of
the energy slider therefore cost about as much I/O as reading the whole map,
and dragging the slider repeated that for every step: ~0.1-0.2 s a tick on a
local disk, far more from a network share. No amount of memory helps that; it
is the access pattern.

What this does
--------------
:class:`CubeCache` wraps the lazy array and indexes exactly like it. The
first read goes straight to the file, so opening a map stays instant and a
map that is only glanced at is never read whole. From the second read on --
the moment someone is actually browsing -- the cube is read into memory once
(one sequential read, the same I/O as a *single* contour used to cost) and
every later slice is a view of that array: microseconds instead of a file
read.

A cube too large for the memory budget is never read whole. It falls back to
a small LRU cache of *blocks* along whichever axis is being browsed: one read
of a block of neighbouring planes costs about what one plane did, and the
next ticks of the slider land in it.

Index arrays that are contiguous runs (what every integration window is) are
turned into slices on the way through: to h5py a slice is one hyperslab,
while an index array is a much slower point selection, and to numpy a slice
is a view rather than a copy.

Nothing here imports Qt.
"""
from __future__ import annotations

from collections import OrderedDict

import numpy as np

#: Never let the cached cube claim more than this, whatever the machine has.
MAX_BUDGET = 4 * 1024 ** 3
#: The share of *currently available* memory a cached cube may take.
AVAILABLE_SHARE = 0.5
#: Used when the available memory cannot be found out.
FALLBACK_BUDGET = 1024 ** 3
#: In the block fallback, how many blocks the budget is split between.
BLOCKS_IN_BUDGET = 4


def default_budget() -> int:
    """How many bytes a cached cube may take on this machine, now."""
    try:
        from tools.memory import system_memory
        _total, available = system_memory()
    except Exception:                                       # noqa: BLE001
        available = None
    if not available:
        return FALLBACK_BUDGET
    return int(min(MAX_BUDGET, AVAILABLE_SHARE * available))


def _as_slice(sel):
    """``sel`` as a slice if it is a 1-D run of consecutive ascending
    integers, else ``sel`` unchanged. Indexing with the slice gives the same
    shape and values as with the run."""
    if isinstance(sel, (list, tuple)):
        sel = np.asarray(sel)
    if (isinstance(sel, np.ndarray) and sel.ndim == 1 and sel.size
            and np.issubdtype(sel.dtype, np.integer)):
        first = int(sel[0])
        if first >= 0 and np.array_equal(sel, np.arange(first, first + sel.size)):
            return slice(first, first + sel.size)
    return sel


def simplify_key(key, ndim: int) -> tuple:
    """A full-length index tuple, with contiguous index runs as slices."""
    if not isinstance(key, tuple):
        key = (key,)
    positions = [i for i, sel in enumerate(key) if sel is Ellipsis]
    if positions:
        i = positions[0]
        key = key[:i] + (slice(None),) * (ndim - len(key) + 1) + key[i + 1:]
    key = key + (slice(None),) * (ndim - len(key))
    return tuple(_as_slice(sel) for sel in key)


class CubeCache:
    """A lazy array (or a plain one) that is read into memory once it is
    being browsed. See the module docstring.

    Indexes like the array it wraps; ``shape``, ``dtype``, ``ndim``,
    ``size``, ``nbytes``, ``len()`` and ``np.asarray()`` all work, so it
    can stand wherever the lazy cube stood.

    The in-memory copy is read-only: slices handed out are views of it, and
    a caller writing into one would otherwise silently change the map for
    every later slice.
    """

    def __init__(self, source, budget: int = None, lazy_reads: int = 1):
        self._source = source
        self.shape = tuple(int(n) for n in source.shape)
        self.dtype = np.dtype(source.dtype)
        self.ndim = len(self.shape)
        self.size = int(np.prod(self.shape)) if self.shape else 1
        self.nbytes = self.size * self.dtype.itemsize
        self.budget = default_budget() if budget is None else int(budget)
        self._lazy_reads = int(lazy_reads)
        self._reads = 0
        self._full = None
        if isinstance(source, np.ndarray):
            # Already in memory: nothing to read, nothing to copy.
            self._full = source
        # Block fallback: (axis, block number) -> array, least recent first.
        self._blocks = OrderedDict()
        self._block_bytes = 0

    def __len__(self):
        return self.shape[0] if self.shape else 0

    def __repr__(self):
        state = "in memory" if self.loaded else "on disk"
        return f"CubeCache(shape={self.shape}, dtype={self.dtype}, {state})"

    # -- state ---------------------------------------------------------------
    @property
    def loaded(self) -> bool:
        return self._full is not None

    @property
    def fits(self) -> bool:
        return self.nbytes <= self.budget

    def load_pending(self) -> bool:
        """True if the next read will read the whole cube into memory --
        so a GUI can put up a wait cursor first."""
        return (not self.loaded and self.fits
                and self._reads >= self._lazy_reads)

    def load(self):
        """Read the whole cube into memory now (if it fits the budget)."""
        if self.loaded or not self.fits:
            return
        full = np.asarray(self._source[(slice(None),) * self.ndim])
        if not full.flags.c_contiguous:
            # A LazyCube hands back a transposed view of what it read;
            # planes along any axis are then all equally quick to take.
            full = np.ascontiguousarray(full)
        full.setflags(write=False)
        self._full = full
        self._blocks.clear()
        self._block_bytes = 0

    def release(self):
        """Forget the in-memory copy; the next reads go to the file again."""
        if self._full is not self._source:
            self._full = None
        self._blocks.clear()
        self._block_bytes = 0
        self._reads = 0

    # -- reading -------------------------------------------------------------
    def __getitem__(self, key):
        key = simplify_key(key, self.ndim)
        if self._full is not None:
            return self._full[key]
        self._reads += 1
        if self._reads > self._lazy_reads:
            if self.fits:
                self.load()
                return self._full[key]
            blocked = self._from_blocks(key)
            if blocked is not None:
                return blocked
        return self._source[key]

    def materialise(self) -> np.ndarray:
        """The whole cube as a new, writable array."""
        if self._full is None and self.fits:
            self.load()
        if self._full is not None:
            return np.array(self._full)
        return np.asarray(self._source[(slice(None),) * self.ndim])

    def __array__(self, dtype=None, copy=None):
        """A fresh array, like the lazy arrays' own ``__array__``: the
        caller may write into it without touching the cache."""
        data = self.materialise()
        return data if dtype is None else data.astype(dtype, copy=False)

    # -- the fallback for a cube too big to hold -----------------------------
    def _from_blocks(self, key):
        """Serve a read along one axis from cached blocks of planes, or
        return None if ``key`` is not that shape of read."""
        axes = [ax for ax, sel in enumerate(key)
                if not (isinstance(sel, slice) and sel == slice(None))]
        if len(axes) != 1:
            return None
        axis = axes[0]
        sel = key[axis]
        n = self.shape[axis]
        if isinstance(sel, (int, np.integer)):
            start = int(sel) + n if sel < 0 else int(sel)
            stop, squeeze = start + 1, True
        elif isinstance(sel, slice) and sel.step in (None, 1):
            start, stop, _ = sel.indices(n)
            squeeze = False
        else:
            return None
        if not 0 <= start < stop <= n:
            return None

        plane_bytes = max(1, self.nbytes // max(n, 1))
        per_block = max(1, (self.budget // BLOCKS_IN_BUDGET) // plane_bytes)
        if stop - start > per_block * (BLOCKS_IN_BUDGET - 1):
            return None                 # wider than the cache can hold
        parts = []
        for number in range(start // per_block, (stop - 1) // per_block + 1):
            b0 = number * per_block
            block = self._block(axis, number, b0, min(n, b0 + per_block))
            lo, hi = max(start, b0) - b0, min(stop, b0 + per_block) - b0
            parts.append(np.take(block, np.arange(lo, hi), axis=axis))
        out = parts[0] if len(parts) == 1 else np.concatenate(parts, axis=axis)
        return np.take(out, 0, axis=axis) if squeeze else out

    def _block(self, axis, number, b0, b1):
        tag = (axis, number)
        block = self._blocks.get(tag)
        if block is not None:
            self._blocks.move_to_end(tag)
            return block
        key = [slice(None)] * self.ndim
        key[axis] = slice(b0, b1)
        block = np.asarray(self._source[tuple(key)])
        self._blocks[tag] = block
        self._block_bytes += block.nbytes
        while self._block_bytes > self.budget and len(self._blocks) > 1:
            _, old = self._blocks.popitem(last=False)
            self._block_bytes -= old.nbytes
        return block
