"""
loader/nxs_file.py
==================
Read SOLEIL/ANTARES nano-ARPES ``*.nxs`` (NeXus/HDF5) files in Python.

This is a line-by-line port of ``load_soleil_nxs.m`` (kept in this repo for
reference) with **one deliberate, important deviation**, documented below.

--------------------------------------------------------------------------
WHY THIS IS NOT A BLIND TRANSLATION OF THE .m FILE (read this first)
--------------------------------------------------------------------------
MATLAB's ``h5read``/``h5info`` report the dimensions of a multi-dimensional
HDF5 dataset in **reverse order** compared to ``h5py`` (this is documented
MATLAB behaviour: HDF5/C stores arrays row-major, MATLAB is column-major, so
MATLAB flips the reported shape rather than copying data). Concretely, a
dataset written with a Python/numpy shape ``(a, b, c, d)`` is reported by
``size(h5read(...))`` in MATLAB as ``(d, c, b, a)``.

That means the hard-coded ``permute(value, [3 4 2 1])`` (and similar) calls
in ``load_soleil_nxs.m`` cannot be copied verbatim into ``h5py``-based code
-- the axis indices refer to MATLAB's *reversed* dimension order. Doing so
silently would very likely scramble the (x, y, k, E) axes without raising
any error, which is worse than not converting at all.

Instead, this module determines which raw HDF5 axis is which physical
quantity by **matching axis lengths** against the independently-read
calibration arrays (the actuator/scale datasets, e.g. ``data.x``, ``data.y``,
``data.k``, ``data.z`` in the .m file) and then transposes into an explicit,
named, canonical order. See :func:`align_and_transpose`.

This is more robust than an index transcription, but it is only as good as
the assumption that no two axes happen to have the same length. Use
:func:`inspect_nxs` on a real file to sanity-check the result before trusting
it (print shapes, compare to the calibration-array lengths) -- this is
strongly recommended, since no real ``*.nxs`` file was available while
writing this module (only the ``.m`` source was provided).

--------------------------------------------------------------------------
Which entry, and which case
--------------------------------------------------------------------------
A ``*.nxs`` file can contain more than one top-level NeXus entry (e.g. a
real-space navigation image recorded right before the actual k-space scan,
both saved in the same file). Which entry is parsed, and as which "case"
(2/3/4, see below), is determined by **inspecting what datasets actually
exist inside each entry** (:func:`_classify_entry`) -- entry names like
``img1_0001``/``img1_0002`` turned out on real files to be generic,
sequentially-numbered labels, not a reliable indicator of content (the same
name pattern held a valid scan in one file and an incomplete/unrelated
entry in another). See :func:`list_entries` and the ``entry=`` parameter of
:func:`load_soleil_nxs`.

--------------------------------------------------------------------------
Data "kinds" produced
--------------------------------------------------------------------------
- ``"spem_4d"``  : real-space scan, full 4D cube (x, y, k, E) + a reduced 3D
                   preview cube. This is the nano-ARPES *imaging* mode.
- ``"spem_1d"``  : real-space scan along a single spatial line (x, k/angle, E
                   -- named ``data.y``/``data.z`` in the .m file for this
                   branch, kept as-is here).
- ``"cut"``      : a single E-vs-angle(k) spectrum (deflector fixed).
- ``"map"``      : a deflector-angle scan built into a (deflx, k, E) cube --
                   this is the k-space "Fermi surface" mode, i.e. the raw
                   material for a (kx, ky, E) visualization. See
                   :func:`deflector_angle_to_kx` to convert the raw deflector
                   angle axis to a momentum axis.
- ``"unsupported"``: case 1 (``Scan2D_MBS_vs_PIX_PIY``), left unimplemented
                   in the .m file too ("not finish for 2019 version").
"""
from __future__ import annotations

import os
import warnings
from dataclasses import dataclass, field
from typing import Optional

import h5py
import numpy as np

# --------------------------------------------------------------------------
# Group name -> case dispatch.
#
# IMPORTANT: names like "img1_0001"/"img1_0002" turned out on real files to
# be generic, sequentially-numbered entry labels, NOT a reliable indicator
# of which case an entry holds -- the same name pattern can be a valid
# real-space scan in one file and an incomplete/unrelated entry (missing
# the datasets a parser needs) in another. Case dispatch is therefore
# primarily CONTENT-based now (see `_classify_entry`), which inspects what
# datasets actually exist inside the entry. This name table is kept only
# as a fallback for entries whose content doesn't match any known
# structure (currently just case 1, which is unimplemented -- like the .m
# file -- and so has nothing to detect by content).
# --------------------------------------------------------------------------
#: Marks an entry written by :func:`save_dataset`, so it is recognised on
#: the way back in regardless of what the user named it.
NATIVE_ATTR = "nxsloader_kind"
NATIVE_FORMAT_ATTR = "nxsloader_format"
#: 1: gzip-4, arrays written as whatever was handed in (in practice float64,
#:    since every caller promoted first), axes named axis_x/y/z only.
#: 2: gzip-4 **with the shuffle filter**, the array's own dtype preserved,
#:    and a fourth axis (axis_w) for a 4-D spatial scan. Measured on a
#:    120x941x96 map: 11.1 -> 8.8 MB at the same write time (shuffle
#:    re-orders the bytes of each value so that the high bytes of
#:    neighbouring numbers, which are nearly always equal, end up adjacent,
#:    which is what gzip is good at). Version 1 files still read.
NATIVE_FORMAT_VERSION = 2

#: How arrays are written. Chosen by measurement rather than by taste, on
#: ARPES-shaped data (smooth bands plus Poisson noise):
#:
#:   raw                 86.7 MB   write 0.04 s
#:   gzip 4              11.1 MB   write 0.70 s     <- version 1
#:   gzip 1 + shuffle     9.5 MB   write 0.47 s
#:   gzip 4 + shuffle     8.8 MB   write 0.70 s     <- version 2
#:   gzip 9 + shuffle     8.4 MB   write 9.87 s
#:   lzf + shuffle       13.5 MB   write 0.33 s
#:
#: Level 9 costs 14x the write time for 5% -- no. lzf is faster to write but
#: half as effective. Chunking is left to h5py: a hand-tuned "whole plane"
#: chunk reads constant-energy slices 4x faster but the orthogonal cuts 20x
#: slower, and the viewer does both, while h5py's own guess is within a few
#: per cent of the best balanced choice.
NATIVE_COMPRESSION = dict(compression="gzip", compression_opts=4,
                          shuffle=True, chunks=True)

GROUP_NAME_TO_CASE = {
    "/Scan2D_MBS_vs_PIX_PIY": 1,   # unimplemented in the original .m file too
}


def _classify_entry(f: h5py.File, g: str) -> Optional[int]:
    """Determine which parsing case (2, 3 or 4) the entry at group path
    ``g`` actually contains, by checking which datasets exist inside it --
    robust to arbitrary/reused entry names. Returns None if nothing
    recognized (falls back to :data:`GROUP_NAME_TO_CASE` by name in that
    case, then to raising "unrecognized" if that also comes up empty)."""
    sd, traj = f"{g}/scan_data", f"{g}/scan_config/trajectory"

    # Case 5: written by this program's own "Save dataset" (see
    # save_dataset). Checked first because it is unambiguous -- the entry
    # says so in an attribute -- and because a saved k-map has axes that are
    # momenta, which no beamline layout would describe.
    try:
        if NATIVE_ATTR in f[g].attrs:
            return 5
    except (KeyError, AttributeError):
        pass

    # Case 2 (SPEM real-space scan): two independent spatial actuators
    # under scan_config/trajectory, plus the 4D-candidate cube dataset
    # (data_12 is read regardless of whether it turns out to be 4D or not
    # -- that decision happens inside _parse_case2).
    if f"{sd}/data_12" in f and traj in f:
        traj_grp = f[traj]
        actuators = [k for k in traj_grp.keys() if k.startswith("actuator_") and k.endswith("_1")]
        if len(actuators) >= 1 and f"{sd}/data_04" in f:
            return 2

    # Case 3/4 (deflector-angle scan): a single scan actuator living
    # directly under scan_data (NOT under scan_config/trajectory -- that
    # distinguishes it from case 2's spatial actuators).
    if f"{sd}/actuator_1_1" in f:
        has_case4 = all(f"{sd}/{d}" in f for d in
                         ("data_04", "data_05", "data_06", "data_07", "data_08", "data_09", "data_11"))
        has_case3 = all(f"{sd}/{d}" in f for d in
                         ("data_01", "data_02", "data_03", "data_04", "data_05", "data_06", "data_09"))
        if has_case4:
            return 4
        if has_case3:
            return 3

    return GROUP_NAME_TO_CASE.get(g)


# --------------------------------------------------------------------------
# Small helpers
# --------------------------------------------------------------------------
def _decode(x):
    """Decode bytes / numpy bytes_ scalars & arrays to str, pass through else."""
    if isinstance(x, (bytes, np.bytes_)):
        return x.decode("utf-8", errors="replace")
    if isinstance(x, np.ndarray) and x.dtype.kind == "S":
        return np.char.decode(x, "utf-8")
    return x


def _read(f: h5py.File, path: str):
    """h5read-equivalent: read a dataset, decoding strings, tolerant of a
    trailing '/' (the .m file appends one after every path)."""
    path = path.rstrip("/")
    return _decode(f[path][()])


def _read_opt(f: h5py.File, path: str, default=None):
    path = path.rstrip("/")
    if path not in f:
        return default
    return _decode(f[path][()])


def matlab_colon(start: float, step: float, stop: float) -> np.ndarray:
    """Reproduce MATLAB's ``start:step:stop`` (inclusive of ``stop`` up to
    floating point tolerance). Defensively coerces its arguments to plain
    Python floats first -- see :func:`scalar0` for why that isn't always a
    no-op on real files."""
    start = float(np.ravel(start)[0]) if np.ndim(start) else float(start)
    step = float(np.ravel(step)[0]) if np.ndim(step) else float(step)
    stop = float(np.ravel(stop)[0]) if np.ndim(stop) else float(stop)
    if step == 0:
        return np.array([start], dtype=float)
    n = int(np.floor((stop - start) / step + 1e-9)) + 1
    n = max(n, 0)
    return start + step * np.arange(n)


def scalar0(arr, name: str = "value") -> float:
    """Extract a single representative scalar from a calibration dataset,
    the way MATLAB's ``arr(1)`` linear indexing does regardless of the
    array's actual shape.

    On real files these "single-value" calibration datasets (escalemin,
    xscalemult, ...) are sometimes logged as 2D arrays (e.g. one row per
    acquisition frame) rather than the flat 1D vectors implicitly assumed
    by the .m file's ``arr(1)`` indexing -- plain numpy ``arr[0]`` then
    returns a sub-array instead of a number, which is a bug, not a
    replication of the MATLAB behaviour. This flattens first (matching
    "just give me the very first stored element") and warns if the values
    are not all (approximately) equal, since ``arr(1)`` implicitly assumes
    they are.
    """
    a = np.asarray(arr)
    flat = a.ravel()
    if flat.size == 0:
        raise ValueError(f"{name}: empty array, cannot extract a scalar (shape={a.shape}).")
    v0 = flat[0]
    if flat.size > 1 and np.issubdtype(flat.dtype, np.number):
        spread = float(np.ptp(flat.astype(float)))
        if spread > 1e-6 * (abs(float(v0)) + 1e-12):
            warnings.warn(
                f"{name}: expected a constant-valued calibration array (the "
                f".m file always reads element (1) of it) but values span "
                f"{spread:.6g} across shape {a.shape}; using the first "
                f"element ({v0!r}). If this file's calibration genuinely "
                f"drifts, treating it as constant is wrong -- inspect it by hand."
            )
    return float(v0)


def squeeze_to_1d(arr, name: str = "array") -> np.ndarray:
    """Collapse a nominally-1D scan-axis array (e.g. an actuator readback)
    to a true 1D array, tolerating harmless extra singleton dimensions
    (shape (N,1), (1,N), ...) that show up on some real files. Warns
    (rather than silently guessing) if a genuinely multi-valued extra
    dimension is found, since flattening that would silently interleave
    unrelated data."""
    a = np.asarray(arr)
    squeezed = np.squeeze(a)
    if squeezed.ndim == 0:
        return squeezed.reshape(1).astype(float)
    if squeezed.ndim != 1:
        warnings.warn(
            f"{name}: expected a 1D scan axis, got shape {a.shape} (still "
            f"{squeezed.ndim}D after squeezing singleton dims); flattening "
            f"in C order -- verify this is meaningful for your file with "
            f"inspect_nxs()."
        )
        squeezed = squeezed.reshape(-1)
    return squeezed.astype(float)


# --------------------------------------------------------------------------
# Metadata fields deliberately NOT shown in the GUI's info table.
#
# These are read-able but were judged noise for day-to-day work at the
# beamline: fixed provenance strings, the User group, the parts of the
# machine record that never change, the thermocouples, and the scan_config
# summary (whose useful content is duplicated under Traj.*). Delete an entry
# here to bring a field back; nothing else needs changing.
# --------------------------------------------------------------------------
INFO_FIELDS_HIDDEN = {
    "title", "experiment_identifier", "run_cycle",
    "Machine.name", "Machine.probe", "Machine.type",
}
INFO_PREFIXES_HIDDEN = ("User.", "TC1.", "TC2.", "ScanCfg.")


def _hide_noisy_info(info: dict) -> dict:
    return {k: v for k, v in info.items()
            if k not in INFO_FIELDS_HIDDEN and not k.startswith(INFO_PREFIXES_HIDDEN)}


# --------------------------------------------------------------------------
# Axis labels and units
# --------------------------------------------------------------------------
#: The analyser's slit-direction axis and the deflector axis are both
#: reported in degrees, and energies in eV, throughout these files.
ANGLE_UNIT = "\u00b0"
ENERGY_UNIT = "eV"


def spatial_unit_for(actuator_names) -> str:
    """Unit for the real-space scan axes, inferred from the actuator names
    recorded in ``scan_config/trajectory/actuator_*_1/name``.

    The piezo stage (PIX/PIY) moves in micrometres; the coarse sample
    stages are in millimetres, which is the default when the names don't
    identify a piezo axis. Detection is name-based because the files carry
    no unit attribute on these datasets.
    """
    joined = " ".join(str(n).lower() for n in actuator_names if n is not None)
    if any(tag in joined for tag in ("pix", "piy", "pi_x", "pi_y", "ex-pi")):
        return "\u00b5m"
    return "mm"
# --------------------------------------------------------------------------
def resolve_axis_order(shape, axis_lengths: "dict[str, int]", order: "list[str]",
                       what: str = "array") -> "tuple[int, ...]":
    """Work out which raw axis is which named axis, by matching lengths.

    Split out of :func:`align_and_transpose` so the same reasoning can be
    applied to a dataset that is never read into memory (see
    :class:`LazyCube`): only ``shape`` is needed, not the data.
    """
    if len(shape) != len(order):
        raise ValueError(
            f"resolve_axis_order: {what} has {len(shape)} dims but {len(order)} "
            f"named axes were requested ({order}). Raw shape={tuple(shape)}. "
            f"Use inspect_nxs() to look at the file and adjust the parsing."
        )
    remaining = list(range(len(shape)))
    perm = []
    for name in order:
        target = axis_lengths[name]
        candidates = [ax for ax in remaining if shape[ax] == target]
        if not candidates:
            raise ValueError(
                f"resolve_axis_order: no remaining raw axis of length {target} "
                f"for '{name}'. Raw shape={tuple(shape)}, remaining={remaining}, "
                f"expected lengths={axis_lengths}."
            )
        if len(candidates) > 1:
            warnings.warn(
                f"resolve_axis_order: axis '{name}' (length {target}) is ambiguous "
                f"-- {len(candidates)} raw axes share that length ({candidates}). "
                f"Picking axis {candidates[0]}; verify with inspect_nxs() if this "
                f"looks wrong."
            )
        perm.append(candidates[0])
        remaining.remove(candidates[0])
    return tuple(perm)


class _HandleRegistry:
    """One open HDF5 handle per file path, reference counted.

    Two reasons this is not left to `h5py.File` per caller:

    * A file often holds several measurements (six Cuts or Maps is normal),
      and the user may have several of them open at once. One handle serves
      them all instead of one per window.
    * HDF5 keeps a single underlying handle per path per process, so two
      `h5py.File` objects for the same file are not independent: on some
      builds, closing either one invalidates datasets held by the other
      (which is exactly what a LazyCube holds). Routing every open through
      here means the file is closed once, when the last user lets go.
    """

    def __init__(self):
        self._open = {}          # abspath -> [h5py.File, refcount]

    def acquire(self, path: str) -> h5py.File:
        key = os.path.abspath(path)
        entry = self._open.get(key)
        if entry is not None and entry[0].id.valid:
            entry[1] += 1
            return entry[0]
        handle = h5py.File(path, "r")
        self._open[key] = [handle, 1]
        return handle

    def release(self, path: str) -> None:
        key = os.path.abspath(path)
        entry = self._open.get(key)
        if entry is None:
            return
        entry[1] -= 1
        if entry[1] <= 0:
            self._open.pop(key, None)
            try:
                entry[0].close()
            except Exception:
                pass

    def open_count(self) -> int:
        """Number of files currently held open (used by the tests)."""
        return len(self._open)


HANDLES = _HandleRegistry()


class LazyCube:
    """A 4D cube that stays on disk, presented in (y, x, k, E) order.

    A real spatial scan is 46x91x96x941 -- around a gigabyte, and often on a
    network share. Reading it whole just to open a file froze the GUI for
    many seconds (long enough that clicks were dropped, so the window seemed
    to need several double-clicks). Almost nothing needs the whole cube: a
    cursor move needs one (k, E) frame, and the spatial overview comes from
    the reduced ``data_11``/``data_01`` datasets the file already carries.

    Indexing works exactly as on the equivalent numpy array -- ints and
    slices in (y, x, k, E) order -- and reads only the selection. The file
    handle stays open, so the owning :class:`NxsScan` must be closed when
    finished with (``NxsData`` does this).
    """

    def __init__(self, dataset, perm: "tuple[int, ...]"):
        self._dset = dataset
        self._perm = perm
        self.shape = tuple(dataset.shape[ax] for ax in perm)
        self.dtype = dataset.dtype
        self.ndim = len(self.shape)

    def __len__(self):
        return self.shape[0]

    def _normalise(self, key):
        if not isinstance(key, tuple):
            key = (key,)
        # `Ellipsis in key` and `key.index(Ellipsis)` both compare with ==,
        # which on an array index returns an array and then raises "the
        # truth value of an array ... is ambiguous". Identity is what is
        # meant anyway -- there is only one Ellipsis object.
        positions = [i for i, sel in enumerate(key) if sel is Ellipsis]
        if positions:
            i = positions[0]
            key = key[:i] + (slice(None),) * (self.ndim - len(key) + 1) + key[i + 1:]
        return key + (slice(None),) * (self.ndim - len(key))

    def __getitem__(self, key):
        key = self._normalise(key)
        # Reorder the selection into the dataset's own axis order, read, then
        # put the surviving axes back into (y, x, k, E) order.
        src_key = [None] * self.ndim
        for out_ax, sel in enumerate(key):
            src_key[self._perm[out_ax]] = sel
        data = self._dset[tuple(src_key)]

        surviving_src = [ax for ax in range(self.ndim)
                         if not isinstance(src_key[ax], (int, np.integer))]
        wanted_src = [self._perm[out_ax] for out_ax, sel in enumerate(key)
                      if not isinstance(sel, (int, np.integer))]
        if len(wanted_src) < 2:
            return data
        return np.transpose(data, [surviving_src.index(ax) for ax in wanted_src])

    @property
    def size(self) -> int:
        return int(np.prod(self.shape)) if self.shape else 1

    @property
    def nbytes(self) -> int:
        return self.size * self.dtype.itemsize

    def materialise(self) -> np.ndarray:
        """Read the whole cube into memory, in this object's axis order.
        Only for callers that genuinely need every point; everything in the
        GUI avoids this."""
        return self[(slice(None),) * self.ndim]

    def __array__(self, dtype=None):
        """So ``np.asarray(cube)`` works -- see :class:`LazyArray`."""
        data = self.materialise()
        return data if dtype is None else data.astype(dtype)


class LazyArray:
    """An array of any rank that stays in its HDF5 file until indexed.

    :class:`LazyCube` above does this for the 4-D spatial scan with an axis
    permutation; this is the same idea without the permutation and for any
    rank, which is what a map, a cut and a saved dataset need. The two are
    kept separate rather than merged because the permutation is the whole
    point of the one and dead weight in the other.

    It exposes ``shape``/``dtype``/``ndim`` and passes indexing straight
    through to h5py, so slicing code cannot tell the difference. The one
    addition that makes it safe to drop in everywhere is ``__array__``:
    ``np.asarray(lazy)`` reads the whole thing, so any routine that really
    does need every point -- an operation in ``tools.process``, a save, a
    fit -- keeps working unchanged, and only pays for what it asked for.

    The file handle belongs to the owning :class:`NxsScan`; once that is
    closed this can no longer be read.
    """

    def __init__(self, dataset):
        self._dset = dataset
        self.shape = tuple(dataset.shape)
        self.dtype = dataset.dtype
        self.ndim = len(self.shape)
        self.size = int(np.prod(self.shape)) if self.shape else 1
        self.nbytes = self.size * self.dtype.itemsize

    def __len__(self):
        return self.shape[0] if self.shape else 0

    def __getitem__(self, key):
        return self._dset[key]

    def materialise(self) -> np.ndarray:
        return self._dset[()]

    def __array__(self, dtype=None):
        data = self.materialise()
        return data if dtype is None else data.astype(dtype)

    def __repr__(self):
        return f"LazyArray(shape={self.shape}, dtype={self.dtype})"


def align_and_transpose(arr: np.ndarray, axis_lengths: "dict[str, int]", order: "list[str]") -> np.ndarray:
    """Reorder ``arr``'s axes to match ``order`` by matching each named
    axis's expected length (``axis_lengths[name]``) against ``arr.shape``.

    Raises ValueError if ``arr.ndim != len(order)``. Warns (does not raise)
    if two candidate axes share the same length -- in that case the first
    unused match is taken and the ambiguity is reported, since it cannot be
    resolved from lengths alone.
    """
    return np.transpose(arr, resolve_axis_order(arr.shape, axis_lengths, order,
                                               what="array"))


# --------------------------------------------------------------------------
# Momentum conversion (only needed for the deflector-scanned angle axis --
# the analyzer's own slit-direction scale, "xscale" below, is already
# delivered in k by the MBS acquisition software, matching the .m file's own
# "SPEM_kmin/kmax" naming for that axis).
# --------------------------------------------------------------------------
def deflector_angle_to_k(angle_deg: np.ndarray, kinetic_energy_eV) -> np.ndarray:
    """Standard free-electron final-state small-angle ARPES conversion::

        k_parallel [A^-1] = 0.5123 * sqrt(KE[eV]) * sin(angle[rad])

    ``kinetic_energy_eV`` may be a scalar (single reference KE, e.g. center
    kinetic energy) or an array broadcastable against ``angle_deg`` (e.g. one
    KE per energy channel, to convert a whole (angle, E) plane at once).

    Caveats (verify against your own calibration before publishing numbers):
    - assumes normal-incidence / zero inner-potential correction (no k_z
      refraction correction at the sample surface),
    - assumes ``angle_deg`` is the true emission angle in degrees (this
      matches the ANTARES "DeflX" deflector convention, but was not
      independently verified against a real data file),
    - the analyzer slit-direction axis is NOT passed through this function
      -- it is already in k according to the acquisition software.
    """
    angle_rad = np.deg2rad(np.asarray(angle_deg, dtype=float))
    KE = np.asarray(kinetic_energy_eV, dtype=float)
    return 0.5123 * np.sqrt(KE) * np.sin(angle_rad)


# --------------------------------------------------------------------------
# Result container
# --------------------------------------------------------------------------
#: Which of :class:`NxsScan`'s axis slots hold each kind's array dimensions,
#: **in the array's own order**, and which slot holds its constructor's.
#:
#: This is the one authority for the question "what are this dataset's
#: axes?". It used to be answered in five places -- ``tools.dataops``'s two
#: tables, the loader's ``_axis_slots``, ``loader.session.scan_to_dict``, the
#: launcher's ``_dataset_for_saving`` and ``_MemScan.__init__`` -- which is
#: four chances for a new kind to be half-supported, and exactly the kind of
#: duplication that makes adding a fifth axis (photon energy, delay, kz) a
#: survey of the whole program rather than an edit.
#:
#: ``array`` is the order the dimensions are stored in; ``constructor`` is
#: the order :class:`ui.widgets.MemoryData` takes them in. They differ
#: for a 4-D spatial scan, whose cube is stored (y, x, k, E) because that is
#: the order the file's own dataset is in, and constructed (x, y, k, E).
AXIS_SLOTS = {
    "cut": {"array": ("x", "y"), "constructor": ("x", "y")},
    "map": {"array": ("x", "k", "z"), "constructor": ("x", "k", "z")},
    "k_map": {"array": ("x", "k", "z"), "constructor": ("x", "k", "z")},
    "kz_map": {"array": ("x", "k", "z"), "constructor": ("x", "k", "z")},
    "spem_1d": {"array": ("x", "y", "z"), "constructor": ("x", "y", "z")},
    "spem_4d": {"array": ("y", "x", "k", "z"),
                "constructor": ("x", "y", "k", "z")},
}

#: The kinds that are a three-axis cube of (scanned axis, analyser angle or
#: momentum, energy). They share a viewer, a 3-D view and a processing
#: panel, and differ only in what their first axis *means* -- which is
#: exactly the distinction the program is otherwise careless about, so it is
#: worth having the set written down once:
#:
#: ``map``     deflector or polar angle, and so convertible to momentum
#: ``k_map``   already converted: in-plane momentum
#: ``kz_map``  photon energy, i.e. an out-of-plane (k_z) scan. Not
#:             convertible by the in-plane formula -- turning a photon
#:             energy into k_z needs the inner potential, which is a
#:             property of the sample and not of the measurement.
CUBE_KINDS = ("map", "k_map", "kz_map")

#: What each kind is called in the file browser's "kind" column. A kind with
#: no entry shows its internal name, which is better than showing nothing.
KIND_LABELS = {
    "cut": "Cut",
    "map": "Map",
    "k_map": "k-map",
    "kz_map": "kz map",
    "spem_4d": "SPEM",
    "spem_1d": "SPEM",
}


def axis_slots(kind: str, order: str = "array") -> tuple:
    """The axis slot names for ``kind``; ``()`` for one with no plain
    axes-and-array form (``unsupported``). ``order`` is "array" or
    "constructor" -- see :data:`AXIS_SLOTS`."""
    return AXIS_SLOTS.get(kind, {}).get(order, ())


def energy_slot(kind: str):
    """Which axis of ``kind`` is the energy, or ``None`` if it has none.

    Derived rather than tabulated. Energy is the last axis a dataset is
    constructed with, for every kind there is: a cut is (angle, energy), a
    cube is (scanned, angle, energy), a 4-D scan is (x, y, angle, energy).
    Two panels each carried their own hand-written ``{"cut": "y", "map":
    "z", ...}`` copy of that, which is two more places to forget when a kind
    is added -- and both would have silently refused to shift the energy
    axis of a ``kz_map`` with "has no energy axis to shift".
    """
    slots = axis_slots(kind, "constructor")
    return slots[-1] if slots else None


@dataclass
class NxsScan:
    kind: str                       # 'spem_4d' | 'spem_1d' | 'cut' | 'map' | 'unsupported'
    filename_prefix: str = ""
    info: dict = field(default_factory=dict)
    fourd_info: dict = field(default_factory=dict)

    # axes (only the ones relevant to `kind` are populated)
    x: Optional[np.ndarray] = None   # spatial x, or deflector angle (Map), or k/angle (Cut)
    y: Optional[np.ndarray] = None   # spatial y, or analyzer k axis (Map), or Energy (Cut)
    k: Optional[np.ndarray] = None   # analyzer momentum axis (spem_4d, map)
    z: Optional[np.ndarray] = None   # Energy axis
    kx: Optional[np.ndarray] = None  # converted deflector-angle -> k (map only, see .to_kspace())

    value: Optional[np.ndarray] = None     # main 2D/3D array for this `kind`
    #: spem_4d only: the (y, x, k, E) cube. A :class:`LazyCube` reading from
    #: the still-open file, not a numpy array -- index it as usual.
    value4d: object = None
    #: spem_4d only: {name: (dataset, perm)} for the reduced preview cubes,
    #: left unread; see NxsData.spatial_overview.
    previews: dict = field(default_factory=dict)
    #: open file handle when the scan reads lazily; released by close().
    _h5file: object = None
    #: the path that handle belongs to, for releasing it in the registry.
    _h5path: object = None

    def axis_slots(self, order: str = "array") -> tuple:
        """The names of the axis attributes this scan's dimensions live in
        (see :data:`AXIS_SLOTS`)."""
        return axis_slots(self.kind, order)

    def axes(self, order: str = "array") -> list:
        """This scan's axis vectors, in the array's order.

        Asking the scan rather than looking the kind up in a table is what
        lets a caller be written once for every kind: ``for axis in
        scan.axes()`` works for a 2-D cut and a 4-D spatial scan alike.
        """
        return [getattr(self, slot) for slot in self.axis_slots(order)]

    def array(self):
        """The dataset's own array -- the 4-D cube for a spatial scan, the
        value for everything else. May be lazy; ``np.asarray`` it if every
        point is needed."""
        if self.kind == "spem_4d" and self.value4d is not None:
            return self.value4d
        return self.value

    def close(self):
        """Release this scan's claim on the file. Safe to call more than
        once. The file itself is only closed when the last scan using it
        lets go, so other entries of the same file keep working. After this
        a LazyCube in ``value4d`` can no longer be indexed."""
        if self._h5file is not None:
            self._h5file = None
            if self._h5path is not None:
                HANDLES.release(self._h5path)
                self._h5path = None

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    #: Display labels (including units) for each axis this scan populates,
    #: keyed 'x', 'y', 'k', 'z'. Filled by the parsers; the GUI uses these
    #: verbatim for its axis titles so units live in one place.
    labels: dict = field(default_factory=dict)


# --------------------------------------------------------------------------
# Common ("FourDinfo" + acquisition/beamline) metadata, identical for every
# case in the .m file.
# --------------------------------------------------------------------------
def _scalarize(value):
    """Collapse a 1-element array to a plain scalar for display; leave real
    arrays and strings alone. Metadata datasets in these files are usually
    shape-(1,) rather than true scalars, which displays as '[3.5]' unless
    unwrapped."""
    if isinstance(value, np.ndarray):
        if value.size == 1:
            return value.reshape(-1)[0]
        if value.size == 0:
            return None
    return value


def _read_group_fields(f: h5py.File, group_path: str, prefix: str, out: dict,
                        skip=("controller_record",)):
    """Copy every scalar/string dataset directly under ``group_path`` into
    ``out`` as ``prefix.<field>``. Used instead of a hand-maintained field
    list so new fields the beamline adds show up automatically rather than
    being silently dropped.

    ``controller_record`` is skipped by default: it is a long Tango
    configuration blob present on every device group, useless in a summary
    table and long enough to swamp it."""
    if group_path not in f:
        return
    grp = f[group_path]
    if not isinstance(grp, h5py.Group):
        return
    for key in sorted(grp.keys()):
        if key in skip:
            continue
        item = grp[key]
        if not isinstance(item, h5py.Dataset):
            continue
        # Skip bulk data: metadata fields are scalars or short strings, and
        # a full detector frame in the info table would be useless anyway.
        if item.size > 16:
            continue
        try:
            out[f"{prefix}.{key}"] = _scalarize(_decode(item[()]))
        except (OSError, TypeError) as exc:
            out[f"{prefix}.{key}"] = f"<unreadable: {exc}>"


def _read_common_info(f: h5py.File, g: str) -> dict:
    """Read the experiment metadata worth showing in the GUI's info table.

    Field selection follows what actually exists in these files (verified
    against a real ``h5info`` dump of a SOLEIL ANTARES nano-ARPES file):
    NXentry-level provenance, the User group, the MBS analyzer settings,
    the monochromator/undulator/machine state, both thermocouples, and
    every sample/optics motor position. Device groups are read
    field-by-field via :func:`_read_group_fields` rather than through a
    hard-coded list, so fields added by future beamline-software versions
    appear automatically.
    """
    info = {}

    # -- NXentry-level provenance -------------------------------------
    for key in ("title", "experiment_identifier", "start_time", "end_time",
                "duration", "run_cycle"):
        value = _read_opt(f, f"{g}/{key}")
        if value is not None:
            info[key] = _scalarize(value)

    # -- who ran it ----------------------------------------------------
    _read_group_fields(f, f"{g}/User", "User", info)

    # -- analyzer (MBS) ------------------------------------------------
    _read_group_fields(f, f"{g}/ANTARES/MBSAcquisition_1", "MBS", info)

    # -- beamline: monochromator, undulators, machine ------------------
    _read_group_fields(f, f"{g}/ANTARES/i12-m-c04-op-mono1", "Mono", info)
    _read_group_fields(f, f"{g}/ANTARES/hu60", "HU60", info)
    _read_group_fields(f, f"{g}/ANTARES/hu256", "HU256", info)
    _read_group_fields(f, f"{g}/ANTARES/ans-ca-machinestatus", "Machine", info)

    # -- temperatures (both thermocouples; converted to K like the .m file)
    for tc, label in (("i12-m-cx1-ex-tc.1", "TC1"), ("i12-m-cx1-ex-tc.2", "TC2")):
        temp = _read_opt(f, f"{g}/ANTARES/{tc}/temperature")
        if temp is not None:
            temp = _scalarize(temp)
            info[f"{label}.temperature_C"] = temp
            info[f"{label}.temperature_K"] = temp + 273.15
            # The per-thermocouple rows are hidden as noise (they were asked
            # to be), but the sample temperature itself is not noise: the
            # Fermi-edge fit holds the temperature at it, so one plain row
            # survives with the first reading available.
            info.setdefault("SampleTemperature_K", temp + 273.15)
        tctype = _read_opt(f, f"{g}/ANTARES/{tc}/type")
        if tctype is not None:
            info[f"{label}.type"] = _scalarize(tctype)

    # -- scan configuration --------------------------------------------
    _read_group_fields(f, f"{g}/scan_config", "ScanCfg", info)
    traj = f"{g}/scan_config/trajectory"
    if traj in f and isinstance(f[traj], h5py.Group):
        for actuator in sorted(f[traj].keys()):
            _read_group_fields(f, f"{traj}/{actuator}", f"Traj.{actuator}", info)

    # -- convenience aliases, matching the .m file's own naming ---------
    # (kept so anything that read the old key names still works)
    aliases = {
        "PhotonEnergy": f"{g}/ANTARES/i12-m-c04-op-mono1/energy",
        "Grating": f"{g}/ANTARES/i12-m-c04-op-mono1/current_grating_name",
        "exitSlit": f"{g}/ANTARES/i12-m-c04-op-mono1/exit_slit_aperture",
        "resolution": f"{g}/ANTARES/i12-m-c04-op-mono1/resolution",
    }
    for alias, path in aliases.items():
        value = _read_opt(f, path)
        if value is not None:
            info[alias] = _scalarize(value)
    return _hide_noisy_info(info)


#: Which of the two real-space stages is the horizontal ("X") axis, by the
#: actuator name the file records. The beamline uses two pairs -- the coarse
#: sample stages ST/SZ and the piezo stage PIX/PIY -- and in both the first
#: named one is the horizontal axis. Anything else is taken in the order the
#: file lists it, so data from another set-up is labelled X/Y too rather than
#: carrying a stage name only this beamline would recognise.
SPATIAL_HORIZONTAL = ("st", "pix", "pi_x")
SPATIAL_VERTICAL = ("sz", "piy", "pi_y")


def _spatial_labels(f: h5py.File, g: str) -> "tuple[str, str, dict]":
    """Axis titles for the two real-space scan axes, plus what they were
    called in the file.

    The titles are always **X** and **Y** with the unit
    (see :func:`spatial_unit_for`), not the stage's own name: "ST"/"SZ" and
    "PIX"/"PIY" mean the same two directions, and a third set-up's names
    would mean them again. The names themselves are returned separately and
    end up in the metadata table, so nothing is lost.

    The first actuator the file lists is taken as the horizontal axis, which
    is the convention on this beamline (ST and PIX are the horizontal ones,
    and they are listed first). A file that lists them the other way round
    is flagged rather than silently transposed -- which of the two array
    axes is which is decided by length matching, not by this label, so
    quietly relabelling would move the label away from the data.
    """
    traj = f"{g}/scan_config/trajectory"
    names = []
    if traj in f and isinstance(f[traj], h5py.Group):
        for actuator in sorted(f[traj].keys()):
            names.append(_read_opt(f, f"{traj}/{actuator}/name"))
    unit = spatial_unit_for(names)
    first = _decode(names[0]).strip() if len(names) > 0 and names[0] is not None else ""
    second = _decode(names[1]).strip() if len(names) > 1 and names[1] is not None else ""

    detail = {"Spatial.X_actuator": first or "(unnamed)",
              "Spatial.Y_actuator": second or "(unnamed)"}
    if first.lower() in SPATIAL_VERTICAL or second.lower() in SPATIAL_HORIZONTAL:
        message = (f"the file lists {first!r} before {second!r}, the opposite of this "
                   f"beamline's usual order; X is still the first scanned axis")
        warnings.warn(f"{g}: {message}")
        detail["Spatial.axis_order"] = message
    return f"X ({unit})", f"Y ({unit})", detail


def _motor_minus_offset(f: h5py.File, g: str, motor: str):
    pos = _read_opt(f, f"{g}/ANTARES/{motor}/position")
    off = _read_opt(f, f"{g}/ANTARES/{motor}/offset")
    if pos is None:
        return None
    return _scalarize(pos) - (_scalarize(off) if off is not None else 0)


# Motor groups worth reporting, as (device name, friendly label). The first
# eight reproduce the .m file's "FourDinfo"; the rest cover the remaining
# sample/optics stages present in these files (OSA, pinhole) that the .m
# file never read.
_MOTOR_GROUPS = [
    ("i12-m-cx1-ex-sample-mt_sn", "focus"),
    ("i12-m-cx1-ex-sample-mt_st", "roughY"),
    ("i12-m-cx1-ex-sample-mt_sz", "roughX"),
    ("i12-m-cx1-ex-sample-mt_srn", "SRn"),
    ("i12-m-cx1-ex-sample-mt_srz", "SRz"),
    ("i12-m-cx1-op-zp-mt_zs", "ZP_ZS"),
    ("i12-m-cx1-op-zp-mt_zx", "ZP_ZX"),
    ("i12-m-cx1-op-zp-mt_zz", "ZP_ZT"),
    ("i12-m-cx1-op-osa-mt_os", "OSA_S"),
    ("i12-m-cx1-op-osa-mt_ox", "OSA_X"),
    ("i12-m-cx1-op-osa-mt_oz", "OSA_Z"),
    ("i12-m-cx1-op-pin-mt_ps", "Pinhole_S"),
    ("i12-m-cx1-op-pin-mt_px", "Pinhole_X"),
    ("i12-m-cx1-op-pin-mt_pz", "Pinhole_Z"),
]


def _read_fourd_info(f: h5py.File, g: str) -> dict:
    """Motor positions (offset-corrected, as in the .m file). Also reads the
    fine piezo stage, whose real path (``i12-m-cx1-ex-pi/{x,y,z}``) the .m
    file guessed wrong and left commented out -- it stores x/y/z directly
    with no separate offset dataset, so it is read as-is."""
    out = {}
    for device, label in _MOTOR_GROUPS:
        value = _motor_minus_offset(f, g, device)
        if value is not None:
            out[label] = _scalarize(value)

    for axis in ("x", "y", "z"):
        value = _read_opt(f, f"{g}/ANTARES/i12-m-cx1-ex-pi/{axis}")
        if value is not None:
            out[f"fine_{axis}"] = _scalarize(value)
    return out


def _find_actuator_group(f: h5py.File, trajectory_path: str) -> str:
    """Robust replacement for the .m file's fragile positional lookup
    (``hinfo.Groups.Groups(3).Groups.Groups.Name``): find an
    ``actuator_*_1`` subgroup under ``trajectory_path`` by name instead of
    by ordinal position."""
    grp = f[trajectory_path]
    candidates = [k for k in grp.keys() if k.startswith("actuator_") and k.endswith("_1")]
    if not candidates:
        raise KeyError(f"No actuator_*_1 group found under {trajectory_path}")
    if len(candidates) > 1:
        warnings.warn(
            f"Multiple actuator groups under {trajectory_path}: {candidates}; "
            f"using '{candidates[0]}'."
        )
    return f"{trajectory_path}/{candidates[0]}"


# --------------------------------------------------------------------------
# Case 2: real-space (SPEM) scan
# --------------------------------------------------------------------------
def _parse_case2(f: h5py.File, g: str, is_zpalign: bool) -> NxsScan:
    ds = f[f"{g}/scan_data/data_12"]
    if ds.ndim == 4:
        # NOT ds[()]: see LazyCube -- reading the whole cube here is what
        # made opening a spatial scan take many seconds.
        raw4d_shape = ds.shape

        traj = f"{g}/scan_config/trajectory"
        relxdelta = scalar0(_read(f, f"{traj}/actuator_1_1/delta"), "actuator_1_1/delta")
        relxfrom = scalar0(_read(f, f"{traj}/actuator_1_1/from"), "actuator_1_1/from")
        relxto = scalar0(_read(f, f"{traj}/actuator_1_1/to"), "actuator_1_1/to")

        relydelta = scalar0(_read(f, f"{traj}/actuator_2_1/delta"), "actuator_2_1/delta")
        relyfrom = scalar0(_read(f, f"{traj}/actuator_2_1/from"), "actuator_2_1/from")
        relyto = scalar0(_read(f, f"{traj}/actuator_2_1/to"), "actuator_2_1/to")

        escalemin = scalar0(f[f"{g}/scan_data/data_04"][()], "data_04 (escalemin)")
        escalemult = scalar0(f[f"{g}/scan_data/data_05"][()], "data_05 (escalemult)")
        escalemax = scalar0(f[f"{g}/scan_data/data_06"][()], "data_06 (escalemax)")
        xscalemin = scalar0(f[f"{g}/scan_data/data_07"][()], "data_07 (xscalemin)")
        xscalemult = scalar0(f[f"{g}/scan_data/data_08"][()], "data_08 (xscalemult)")
        xscalemax = scalar0(f[f"{g}/scan_data/data_09"][()], "data_09 (xscalemax)")

        x = matlab_colon(relxfrom, relxdelta, relxto)
        y = matlab_colon(relyfrom, relydelta, relyto)
        k = matlab_colon(xscalemin, xscalemult, xscalemax)
        z = matlab_colon(escalemin, escalemult, escalemax)

        perm = resolve_axis_order(raw4d_shape,
                                  {"y": len(y), "x": len(x), "k": len(k), "e": len(z)},
                                  ["y", "x", "k", "e"], what="data_12")
        value4d = LazyCube(ds, perm)

        # The reduced preview cubes the file already carries. data_11 is the
        # cube summed over energy and data_01 the cube summed over the slit
        # angle, so either one collapses to the spatial overview at a
        # thousandth of the data volume. They are kept as open datasets,
        # unread, and NxsData validates one against the real cube before
        # trusting it.
        previews = {}
        for name, axes in (("data_11", ("y", "x", "k")), ("data_01", ("y", "x", "e"))):
            path = f"{g}/scan_data/{name}"
            if path not in f:
                continue
            dset = f[path]
            try:
                lengths = {"y": len(y), "x": len(x), "k": len(k), "e": len(z)}
                previews[name] = (dset, resolve_axis_order(
                    dset.shape, {a: lengths[a] for a in axes}, list(axes), what=name))
            except (ValueError, KeyError) as exc:
                warnings.warn(f"{name}: not usable as a spatial preview ({exc})")
        value3d = None

        scan = NxsScan(kind="spem_4d",
                        filename_prefix="ZPAlign_" if is_zpalign else "SPEM_")
        scan.x, scan.y, scan.k, scan.z = x, y, k, z
        scan.value4d = value4d
        scan.value = value3d
        scan.previews = previews
        scan.info["SPEM_kmin"] = xscalemin
        scan.info["SPEM_kmax"] = xscalemax
        x_label, y_label, spatial_detail = _spatial_labels(f, g)
        scan.info.update(spatial_detail)
        scan.labels = {
            "x": x_label,
            "y": y_label,
            "k": f"Angle along slit ({ANGLE_UNIT})",
            "z": f"Energy ({ENERGY_UNIT})",
        }
        return scan

    # --- fallback: 1D spatial scan ("SPEM1D") ---
    raw3d = f[f"{g}/scan_data/data_11"][()]
    actuator_path = _find_actuator_group(f, f"{g}/scan_config/trajectory")
    relxdelta = scalar0(_read(f, f"{actuator_path}/delta"), "actuator/delta")
    relxfrom = scalar0(_read(f, f"{actuator_path}/from"), "actuator/from")
    relxto = scalar0(_read(f, f"{actuator_path}/to"), "actuator/to")

    escalemin = scalar0(f[f"{g}/scan_data/data_04"][()], "data_04 (escalemin)")
    escalemult = scalar0(f[f"{g}/scan_data/data_05"][()], "data_05 (escalemult)")
    escalemax = scalar0(f[f"{g}/scan_data/data_06"][()], "data_06 (escalemax)")
    xscalemin = scalar0(f[f"{g}/scan_data/data_07"][()], "data_07 (xscalemin)")
    xscalemult = scalar0(f[f"{g}/scan_data/data_08"][()], "data_08 (xscalemult)")
    xscalemax = scalar0(f[f"{g}/scan_data/data_09"][()], "data_09 (xscalemax)")

    x = matlab_colon(relxfrom, relxdelta, relxto)
    yaxis = matlab_colon(xscalemin, xscalemult, xscalemax)  # kept as .y, matching the .m file
    z = matlab_colon(escalemin, escalemult, escalemax)

    value = align_and_transpose(raw3d, {"x": len(x), "y": len(yaxis), "e": len(z)}, ["x", "y", "e"])

    scan = NxsScan(kind="spem_1d", filename_prefix="SPEM1D_")
    scan.x, scan.y, scan.z = x, yaxis, z
    scan.value = value
    x_label, _y_label, spatial_detail = _spatial_labels(f, g)
    scan.info.update(spatial_detail)
    scan.labels = {
        "x": x_label,
        # .y holds the analyser slit axis for this branch, matching the .m file
        "y": f"Angle along slit ({ANGLE_UNIT})",
        "z": f"Energy ({ENERGY_UNIT})",
    }
    return scan


def _best_effort_transpose_3d(arr: np.ndarray, axis_lengths: dict) -> np.ndarray:
    """The .m file's ``data.value`` 3D preview cube isn't spelled out
    axis-by-axis; try the 3 most-likely orders against the known axis
    lengths and fall back to the raw array (with a warning) if none match
    cleanly."""
    names = list(axis_lengths.keys())
    if arr.ndim != 3:
        warnings.warn(f"Unexpected ndim={arr.ndim} for the 3D preview cube; returning raw array as-is.")
        return arr
    from itertools import permutations, combinations
    for combo in combinations(names, 3):
        for order in permutations(combo):
            if all(arr.shape[i] == axis_lengths[order[i]] for i in range(3)):
                return np.transpose(arr, list(range(3)))  # already in a matching order; identity transpose
    warnings.warn(
        f"Could not match the 3D preview cube's shape {arr.shape} to any "
        f"combination of known axis lengths {axis_lengths}; returning it "
        f"unmodified. Use inspect_nxs() to check axis order by hand."
    )
    return arr


# --------------------------------------------------------------------------
# Case 3 / 4: deflector-angle scan (k-space "Cut" or "Map")
# --------------------------------------------------------------------------
def _parse_case34(f: h5py.File, g: str, e_idx: tuple, x_idx: tuple, v_idx: str) -> NxsScan:
    """``e_idx``/``x_idx`` are the ('data_NN','data_NN','data_NN') triplets
    for (escale) and (xscale) that differ between case 3 (01-03 / 04-06) and
    case 4 (04-06 / 07-09); ``v_idx`` is the value dataset name (data_09 /
    data_11)."""
    deflx = squeeze_to_1d(f[f"{g}/scan_data/actuator_1_1"][()], "actuator_1_1 (deflx)")

    escalemin = scalar0(f[f"{g}/scan_data/{e_idx[0]}"][()], f"{e_idx[0]} (escalemin)")
    escalemult = scalar0(f[f"{g}/scan_data/{e_idx[1]}"][()], f"{e_idx[1]} (escalemult)")
    escalemax = scalar0(f[f"{g}/scan_data/{e_idx[2]}"][()], f"{e_idx[2]} (escalemax)")
    xscalemin = scalar0(f[f"{g}/scan_data/{x_idx[0]}"][()], f"{x_idx[0]} (xscalemin)")
    xscalemult = scalar0(f[f"{g}/scan_data/{x_idx[1]}"][()], f"{x_idx[1]} (xscalemult)")
    xscalemax = scalar0(f[f"{g}/scan_data/{x_idx[2]}"][()], f"{x_idx[2]} (xscalemax)")

    escale = matlab_colon(escalemin, escalemult, escalemax)
    xscale = matlab_colon(xscalemin, xscalemult, xscalemax)

    # The dataset, not its contents: the axis order is worked out from the
    # *shape*, which costs nothing to ask for, and the numbers are then left
    # on disk behind a LazyCube. A deflector map is tens to hundreds of MB
    # and reading it whole here was the largest single allocation the
    # program made -- for a view that only ever shows one slice at a time.
    raw_dataset = f[f"{g}/scan_data/{v_idx}"]

    # On real files this comes both as a bare 2D (k,E) array for a true
    # single-position Cut *and* as a 3D (scan_points, k, E) array with
    # scan_points==1 even for a single position (the acquisition software
    # apparently doesn't special-case a 1-point scan) -- handle both rather
    # than assuming the .m file's implicit 2D shape.
    if raw_dataset.ndim == 2:
        perm = resolve_axis_order(raw_dataset.shape,
                                  {"k": len(xscale), "e": len(escale)},
                                  ["k", "e"], what="array")
        value = LazyCube(raw_dataset, tuple(perm))
    elif raw_dataset.ndim == 3:
        perm = resolve_axis_order(
            raw_dataset.shape,
            {"deflx": deflx.size, "k": len(xscale), "e": len(escale)},
            ["deflx", "k", "e"], what="array")
        value3d = LazyCube(raw_dataset, tuple(perm))
        # One deflector position is a Cut, and a Cut is small enough that
        # keeping it on disk buys nothing; read that one frame now.
        value = value3d[0] if deflx.size == 1 else value3d
    else:
        raise ValueError(
            f"{v_idx}: unexpected ndim={raw_dataset.ndim} (shape={raw_dataset.shape}); "
            f"expected 2 (bare Cut) or 3 (scan_points, k, E). Use inspect_nxs() to check."
        )

    if deflx.size == 1:
        scan = NxsScan(kind="cut", filename_prefix="Cut_")
        scan.x, scan.y = xscale, escale
        scan.value = value
        scan.labels = {
            "x": f"Angle along slit ({ANGLE_UNIT})",
            "y": f"Energy ({ENERGY_UNIT})",
        }
        # The one deflector position the cut was taken at. It is not an axis
        # here -- there is only one of it -- but it is half the geometry a
        # k conversion of this cut needs (the other half being where Gamma
        # sits), so it is kept rather than dropped with the axis.
        scan.info["DeflectorAngle_deg"] = float(deflx[0])
        return scan

    # multiple deflector positions -> 3D (deflx, k, E) map
    scan = NxsScan(kind="map", filename_prefix="Map_")
    scan.x, scan.k, scan.z = deflx, xscale, escale
    scan.value = value
    scan.labels = {
        "x": f"Angle, deflector ({ANGLE_UNIT})",
        "k": f"Angle along slit ({ANGLE_UNIT})",
        "z": f"Energy ({ENERGY_UNIT})",
    }
    return scan


# --------------------------------------------------------------------------
# Public API
# --------------------------------------------------------------------------
def _entry_kind(f: h5py.File, g: str) -> Optional[str]:
    """Label one entry: "SPEM", "Cut", "Map", or None if unrecognised.

    Reads only the entry's structure and the *shape* of the deflector
    actuator -- never a data array -- so a directory of multi-GB scans can
    be labelled instantly.
    """
    case = _classify_entry(f, g)
    if case == 5:
        # A dataset this program saved; it states its own kind, which may be
        # one no beamline layout produces (a converted k_map).
        kind = _decode(f[g].attrs.get(NATIVE_ATTR, ""))
        if isinstance(kind, np.ndarray):
            kind = str(kind.reshape(-1)[0])
        return KIND_LABELS.get(str(kind), str(kind) or None)
    if case == 2:
        return "SPEM"
    if case in (3, 4):
        # Cut vs Map is decided by how many deflector positions were
        # scanned; the shape alone answers that, no data needed.
        actuator = f.get(f"{g}/scan_data/actuator_1_1")
        if actuator is None:
            return None
        n_positions = int(np.prod(actuator.shape)) if actuator.shape else 1
        return "Cut" if n_positions <= 1 else "Map"
    return None


def _first_string(f: h5py.File, path: str) -> Optional[str]:
    """A short string dataset's value, or None if absent/unreadable."""
    try:
        if path not in f:
            return None
        value = _decode(f[path][()])
        if isinstance(value, np.ndarray):
            value = value.reshape(-1)[0] if value.size else None
        return None if value is None else str(value).strip() or None
    except Exception:
        return None


def _parse_native(f: h5py.File, g: str) -> NxsScan:
    """Read back an entry written by :func:`save_dataset`.

    Deliberately simple and self-describing: the axes, the cube, the axis
    labels and the metadata are all stored under plain names, so nothing has
    to be inferred the way a beamline file's ``data_09`` does. That also
    means a saved k-map keeps its momentum axes labelled as momenta instead
    of being mistaken for angles.
    """
    grp = f[g]
    kind = _decode(grp.attrs[NATIVE_ATTR])
    if isinstance(kind, np.ndarray):
        kind = str(kind.reshape(-1)[0])
    kind = str(kind)

    scan = NxsScan(kind=kind, filename_prefix="")
    # Left on disk behind a LazyArray, and in the dtype it was written in:
    # a saved dataset is opened exactly like a beamline file, so a session's
    # worth of them can sit in the list without any of them being resident.
    # Axes are small and always read.
    scan.value = LazyArray(f[f"{g}/value"])
    scan.x = np.asarray(f[f"{g}/axis_x"][()], dtype=float)

    if kind == "spem_4d":
        # (x, y, k, E): the cube belongs in value4d, where the viewer looks
        # for it, and value mirrors it as MemoryData does.
        scan.y = np.asarray(f[f"{g}/axis_y"][()], dtype=float)
        scan.k = np.asarray(f[f"{g}/axis_z"][()], dtype=float)
        if f"{g}/axis_w" in f:
            scan.z = np.asarray(f[f"{g}/axis_w"][()], dtype=float)
        scan.value4d = scan.value
    elif f"{g}/axis_z" in f:
        # 3-D: a map or k-map (x, k, E), or a 1-D spatial scan (x, y, E).
        if kind == "spem_1d":
            scan.y = np.asarray(f[f"{g}/axis_y"][()], dtype=float)
        else:
            scan.k = np.asarray(f[f"{g}/axis_y"][()], dtype=float)
        scan.z = np.asarray(f[f"{g}/axis_z"][()], dtype=float)
    else:
        scan.y = np.asarray(f[f"{g}/axis_y"][()], dtype=float)

    labels = {}
    for axis in ("x", "y", "k", "z"):
        name = f"{g}/label_{axis}"
        if name in f:
            value = _first_string(f, name)
            if value:
                labels[axis] = value
    scan.labels = labels

    info_group = f"{g}/info"
    if info_group in f:
        for key, value in f[info_group].attrs.items():
            scan.info[key] = _decode(value)
    motors_group = f"{g}/motors"
    if motors_group in f:
        for key, value in f[motors_group].attrs.items():
            scan.fourd_info[key] = _decode(value)
    return scan


def save_dataset(path: str, datasets: "list[dict]", progress=None) -> None:
    """Write datasets to a ``.nxs`` file this program can open again.

    ``datasets`` is a list of ``{"name", "kind", "axes", "labels", "value",
    "info", "motors"}`` dicts, one per entry; ``axes`` is ``(x, y)`` for a
    cut and ``(x, y, z)`` otherwise, up to four for a spatial scan. Several
    datasets go into one file as several top-level entries, exactly like a
    beamline file holding several measurements -- so saving a selection and
    reopening it gives the same rows back.

    The entries are marked with :data:`NATIVE_ATTR`, which is what
    :func:`_classify_entry` keys on, so the name of the entry is free to be
    whatever the user called the dataset, and a file written here is
    recognised on sight without anyone having to say which beamline it came
    from.

    ``progress`` is an optional ``callable(done, total, label)`` for a
    caller showing a progress bar; it is called once per dataset. Raising
    from it (which is how a Cancel button reports itself) aborts the write.

    The array is stored in **its own dtype**. Earlier versions wrote
    whatever the caller had promoted to float64, which doubled the file for
    data that arrived as float32 or as integer counts, and lost nothing by
    doing so.

    Returns the entry names it actually used, in the order given. They are
    not always the names asked for -- HDF5 group names cannot contain "/",
    and two datasets may be called the same thing -- and the caller needs
    the real ones to be able to point a list row at what was just written.
    """
    used = set()
    written = []
    total = len(datasets)
    with h5py.File(path, "w") as f:
        for index, item in enumerate(datasets):
            if progress is not None:
                progress(index, total, str(item.get("name", "")))
            # HDF5 group names cannot contain '/', and two datasets may well
            # have been given the same name by the user.
            base = str(item["name"]).replace("/", "_").strip() or "dataset"
            name, n = base, 2
            while name in used:
                name, n = f"{base}_{n}", n + 1
            used.add(name)
            written.append(name)

            grp = f.create_group(name)
            grp.attrs[NATIVE_ATTR] = str(item["kind"])
            grp.attrs[NATIVE_FORMAT_ATTR] = NATIVE_FORMAT_VERSION
            grp.attrs["nxsloader_name"] = str(item["name"])

            axes = item["axes"]
            # axis_w is the fourth, for a spatial scan; the first three keep
            # their version-1 names so version-1 files stay readable.
            for axis_name, values in zip(("axis_x", "axis_y", "axis_z", "axis_w"),
                                         axes):
                grp.create_dataset(axis_name, data=np.asarray(values, dtype=float))
            value = np.asarray(item["value"])
            # A 0-d or tiny array cannot be chunked, and compressing it would
            # cost more in filter overhead than it saves.
            options = dict(NATIVE_COMPRESSION) if value.size > 1024 else {}
            grp.create_dataset("value", data=value, **options)

            for axis, label in (item.get("labels") or {}).items():
                if label:
                    grp.create_dataset(f"label_{axis}", data=np.bytes_(str(label).encode()))

            info = grp.create_group("info")
            for key, val in (item.get("info") or {}).items():
                try:
                    info.attrs[str(key)] = val if isinstance(
                        val, (int, float, np.integer, np.floating)) else str(val)
                except (TypeError, ValueError):
                    info.attrs[str(key)] = str(val)
            motors = grp.create_group("motors")
            for key, val in (item.get("motors") or {}).items():
                try:
                    motors.attrs[str(key)] = float(val)
                except (TypeError, ValueError):
                    motors.attrs[str(key)] = str(val)
    return written


def list_datasets(path: str) -> "list[dict]":
    """Every independently-openable dataset in a ``.nxs`` file.

    A single file routinely holds several measurements -- up to six Cuts or
    Maps recorded one after another, each a complete top-level NeXus entry
    with its own ``scan_data``. They are separate measurements that happen
    to share a file, so the browser lists one row per entry and opening one
    opens only that entry.

    Returns a list of dicts, in file order, with keys:

    ``entry``       the NeXus group name, e.g. ``"DeflX_0003"``
    ``kind``        "SPEM" / "Cut" / "Map"
    ``title``       the entry's own ``title``, if it has one
    ``start_time``  when it was recorded, if recorded

    Entries whose layout is not recognised are skipped, and an unreadable
    file yields an empty list rather than raising -- a browser listing must
    not fail because of one bad file.
    """
    out = []
    try:
        with h5py.File(path, "r") as f:
            for name in f.keys():
                g = "/" + name
                try:
                    kind = _entry_kind(f, g)
                except Exception:
                    kind = None
                if kind is None:
                    continue
                out.append({
                    "entry": name,
                    "kind": kind,
                    "title": _first_string(f, f"{g}/title"),
                    "start_time": _first_string(f, f"{g}/start_time"),
                })
    except Exception:
        return []
    return out


def probe_kind(path: str) -> str:
    """The kind of a file's *default* entry (see :func:`_pick_entry`), for
    callers that want one label per file. A browser listing wants
    :func:`list_datasets` instead, which enumerates every entry."""
    try:
        with h5py.File(path, "r") as f:
            top_level = list(f.keys())
            if not top_level:
                return "unknown"
            return _entry_kind(f, _pick_entry(f, path, top_level)) or "unknown"
    except Exception:
        return "unknown"


def list_entries(path: str) -> "list[tuple[str, Optional[int]]]":
    """List every top-level NeXus entry in ``path`` as ``(group_name, case)``
    pairs (``case`` is ``None`` if the entry's content doesn't match any
    known structure). Real SOLEIL files can contain *more than one*
    top-level entry in a single .nxs file -- e.g. a real-space navigation
    image (case 2) recorded right before the actual deflector-angle map
    (case 3/4) that the file is "about", or a stale/incomplete entry left
    over from an aborted acquisition. Use this to see what's really in a
    file before calling :func:`load_soleil_nxs` with an explicit ``entry=``."""
    with h5py.File(path, "r") as f:
        return [("/" + k, _classify_entry(f, "/" + k)) for k in f.keys()]


def _pick_entry(f: h5py.File, path: str, top_level: "list[str]") -> str:
    """Choose which top-level entry to parse when a file has more than one.
    Prefers a recognized entry with the highest numeric suffix (SOLEIL's
    ``..._0001``, ``..._0002``, ... naming), on the assumption that a lower
    number is an earlier, preliminary acquisition (e.g. a navigation image
    taken before the actual scan) -- this is a heuristic, not a guarantee;
    pass ``entry=`` explicitly to load_soleil_nxs() if it picks the wrong
    one for your file. Entries whose *content* doesn't match any known case
    (see :func:`_classify_entry`) are never picked automatically, even if
    their name looks like a match."""
    import re
    candidates = [("/" + k, _classify_entry(f, "/" + k)) for k in top_level]
    recognized = [(g, c) for g, c in candidates if c is not None]
    if not recognized:
        return "/" + top_level[0]  # let the caller raise a clear "unrecognized" error

    if len(recognized) == 1:
        return recognized[0][0]

    def _suffix(name: str) -> int:
        m = re.search(r"(\d+)$", name)
        return int(m.group(1)) if m else -1

    recognized.sort(key=lambda gc: _suffix(gc[0]))
    chosen, _ = recognized[-1]
    warnings.warn(
        f"{path}: multiple usable top-level entries found: {recognized}. "
        f"Using '{chosen}' (highest numeric suffix -- typically the main scan "
        f"rather than a preliminary navigation image, but this is a heuristic). "
        f"Pass entry='{chosen.lstrip('/')}' explicitly to silence this warning, "
        f"or a different entry name from the list above if this guess is wrong "
        f"for your data."
    )
    return chosen


def load_soleil_nxs(path: str, entry: Optional[str] = None) -> NxsScan:
    """Load a SOLEIL ANTARES nano-ARPES ``*.nxs`` file. Returns an
    :class:`NxsScan`. Raises for unsupported/unrecognized files rather than
    silently returning garbage (unlike the .m file, which would error out
    inside an undefined ``switch`` variable in the same situation).

    ``entry``: which top-level NeXus entry to parse, if the file has more
    than one (see :func:`list_entries`). By default the entry is picked
    automatically among the entries whose *content* is recognized (see
    :func:`_pick_entry`); pass this explicitly to override that guess.
    """
    # Any of the HDF5 extensions, not just .nxs. The check is here to catch
    # a file of the wrong *type* being pointed at -- a text spectrum, a
    # spreadsheet -- and for that the family is the right granularity. NeXus
    # is an HDF5 layout, this program's own saved format is written with
    # h5py, and a dataset saved as .h5 or handed over as .hdf5 is readable
    # by exactly the same code; refusing it on its name alone meant a file
    # this program itself wrote could not be reopened.
    if not path.lower().endswith((".nxs", ".h5", ".hdf5", ".hdf", ".nx5")):
        raise ValueError(
            f"Not an HDF5/NeXus file (.nxs, .h5, .hdf5): {path}")

    # Opened through the registry, not with a `with`: a spatial scan keeps
    # reading from the file as the user moves around it (see LazyCube), so
    # the handle has to outlive this call, and several entries of the same
    # file must share one handle. Anything that raises, or that read
    # everything it needed, releases it before returning.
    f = HANDLES.acquire(path)
    keep_open = False
    try:
        top_level = list(f.keys())
        if entry is not None:
            g = entry if entry.startswith("/") else "/" + entry
            if g.lstrip("/") not in f:
                raise KeyError(f"Entry '{g}' not found in {path}; available: {top_level}")
        else:
            g = _pick_entry(f, path, top_level)

        case = _classify_entry(f, g)
        if case is None:
            raise KeyError(
                f"Entry '{g}' in {path} doesn't match any known case (its "
                f"content was inspected, not just its name) -- add support "
                f"for it in loader/nxs_file.py once you know its structure. "
                f"Run inspect_nxs('{path}') to look at the file. "
                f"All top-level entries in this file: {top_level}"
            )

        if case == 1:
            scan = NxsScan(kind="unsupported", filename_prefix="")
        elif case == 2:
            scan = _parse_case2(f, g, is_zpalign=(g == "/ZPAlign_0001"))
        elif case == 3:
            scan = _parse_case34(f, g, ("data_01", "data_02", "data_03"), ("data_04", "data_05", "data_06"), "data_09")
        elif case == 5:
            scan = _parse_native(f, g)
        elif case == 4:
            scan = _parse_case34(f, g, ("data_04", "data_05", "data_06"), ("data_07", "data_08", "data_09"), "data_11")
        else:
            raise AssertionError(f"unreachable case {case}")

        # A natively-saved entry carries its own metadata, already read by
        # _parse_native; the beamline readers would find nothing there.
        if scan.kind != "unsupported" and case != 5:
            scan.info.update(_read_common_info(f, g))
            scan.fourd_info.update(_read_fourd_info(f, g))
        scan.info["_group"] = g
        scan.info["_case"] = case
        if len(top_level) > 1:
            scan.info["_other_entries"] = [t for t in top_level if "/" + t != g]

        # Anything still pointing into the file -- a lazy cube, a lazy value
        # array, or an unread preview -- means the handle has to outlive this
        # call. Getting this condition wrong does not fail here; it fails
        # later, at the first slice, with an h5py "identifier is not of
        # specified type", which is why it is one expression rather than
        # being spread across the parsers.
        if (isinstance(scan.value4d, (LazyCube, LazyArray))
                or isinstance(scan.value, (LazyCube, LazyArray))
                or scan.previews):
            scan._h5file = f
            scan._h5path = path
            keep_open = True
        return scan
    finally:
        if not keep_open:
            HANDLES.release(path)


def to_kspace_cube(scan: NxsScan, kinetic_energy_eV=None, work_function_eV: float = 4.5):
    """For a ``kind == "map"`` scan, convert the raw deflector-angle axis
    (``scan.x``, in degrees) into momentum ``kx`` using
    :func:`deflector_angle_to_k`, returning ``(kx, ky, E, cube)`` where
    ``ky`` is the analyzer's own (already-in-k) axis and ``cube`` has shape
    ``(len(kx), len(ky), len(E))``.

    ``kinetic_energy_eV``: reference kinetic energy for the conversion.
        - If None, uses ``scan.info["center_ke"]`` if available.
        - May be a scalar (one KE for the whole cube -- an approximation,
          fine for a narrow energy window) or an array of length
          ``len(scan.z)`` (one KE per energy channel, for a wide window
          where KE varies non-negligibly across it).
    ``work_function_eV``: only used if you pass binding energies instead of
        kinetic energies for ``scan.z`` and want KE = hv - BE - work
        function; unused by default. Provided for convenience/consistency
        checks, not applied automatically.

    This function does not resample onto a rectilinear (kx, ky) grid -- kx
    is energy-dependent in general (kx = f(angle, KE)) and a proper
    rectilinear (kx, ky, E) volume requires per-slice interpolation. For a
    single representative slice this is usually a fine approximation; for
    quantitative work across a wide energy range, interpolate per-E-slice
    before treating this as a uniform-grid volume.
    """
    if scan.kind != "map":
        raise ValueError(f"to_kspace_cube expects kind=='map', got '{scan.kind}'")
    if kinetic_energy_eV is None:
        # 'MBS.center_ke' is the current namespaced key; 'center_ke' is the
        # older flat one, still accepted so existing scripts keep working.
        kinetic_energy_eV = scan.info.get("MBS.center_ke", scan.info.get("center_ke"))
        if kinetic_energy_eV is None:
            raise ValueError(
                "No kinetic_energy_eV given and neither scan.info['MBS.center_ke'] "
                "nor scan.info['center_ke'] is available in this file.")
    kx = deflector_angle_to_k(scan.x, kinetic_energy_eV)
    return kx, scan.k, scan.z, scan.value


# --------------------------------------------------------------------------
# Diagnostics
# --------------------------------------------------------------------------
def inspect_nxs(path: str, max_depth: int = 6) -> None:
    """Print the HDF5 tree (groups/datasets/shapes) and the detected case,
    for validating the axis-alignment assumptions above against a real
    file. Read-only, side-effect-free besides printing."""
    with h5py.File(path, "r") as f:
        top_level = list(f.keys())
        print(f"File: {path}")
        print(f"Top-level entries and their detected case (content-based, not name-based):")
        for k in top_level:
            g = "/" + k
            case = _classify_entry(f, g)
            print(f"  {g}: {'UNRECOGNIZED' if case is None else f'case {case}'}")
        print("-" * 70)

        def _walk(name, obj, depth):
            if depth > max_depth:
                return
            indent = "  " * depth
            if isinstance(obj, h5py.Dataset):
                print(f"{indent}{name.split('/')[-1]}  shape={obj.shape} dtype={obj.dtype}")
            else:
                print(f"{indent}{name.split('/')[-1]}/")

        f.visititems(lambda name, obj: _walk(name, obj, name.count("/")))


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        inspect_nxs(sys.argv[1])
    else:
        print("Usage: python loader/nxs_file.py <path-to-file.nxs>   # prints file structure")
