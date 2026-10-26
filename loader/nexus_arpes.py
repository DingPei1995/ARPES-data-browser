"""
loader/nexus_arpes.py
=====================
What the ALBA and Diamond loaders share: both beamlines write NeXus files
in which one ``NXdata`` group holds the analyser's spectra as
``(*scanned axes, slit angle, energy)`` beside the vectors of those axes.
The two differ in where that group is and in how the scanned axes are
described (an ``@axes`` list at ALBA, ``@axis``/``@primary`` on each
vector at Diamond), so each beamline's module works out a :class:`Layout`
and hands it here, and everything after that -- which kind the data opens
as, the photon-energy referencing of a kz map, the lazy array, the
metadata -- is written once.

How a layout becomes datasets
-----------------------------
==========================================  ================================
scanned axes                                opens as
==========================================  ================================
none (or only axes of length 1)             ``cut``
one photon energy                           ``kz_map``
one anything else                           ``map``; the axis's role says
                                            whether it is an angle
two spatial motors                          ``spem_4d``
anything more                               one ``map`` (or ``spem_4d``)
                                            per step of the extra axes
==========================================  ================================

A kz map is the one kind read into memory: every photon energy has its own
kinetic window, and putting them on one energy axis means referring each to
the Fermi level by its own ``hv`` (``E - E_F = E_kin - hv + W``) and, when
the referred windows do not coincide, resampling. Everything else stays in
the file behind :class:`NexusView` and is read as it is looked at.

Nothing in here imports Qt.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from itertools import product
from typing import Optional

import h5py
import numpy as np

from loader.nxs_file import (ANGLE_UNIT, ENERGY_UNIT, HANDLES, KIND_LABELS,
                             NxsScan, PHD_WORK_FUNCTION_EV, _decode,
                             _read_group_fields, _resample_energy, _scalarize)

__all__ = ["Dim", "Layout", "NexusView", "NexusArpesLoader", "text",
           "unit_text", "motor_role", "plan_entries", "build_scan"]

#: What the analyser work function is taken to be when a file records none.
#: The same figure as an ANTARES PhD scan, for the same reason: a wrong value
#: is a rigid shift of the energy axis, which the Fermi-edge tools correct.
WORK_FUNCTION_EV = PHD_WORK_FUNCTION_EV


# --------------------------------------------------------------------------
# Small helpers
# --------------------------------------------------------------------------
def text(value) -> str:
    """A dataset or attribute value as plain text."""
    value = _decode(value)
    if isinstance(value, np.ndarray):
        if value.size == 0:
            return ""
        if value.size == 1:
            return text(value.reshape(-1)[0])
        return ", ".join(text(v) for v in value.reshape(-1))
    return str(value).strip()


def attr(item, name: str, default=None):
    try:
        return item.attrs[name] if name in item.attrs else default
    except (OSError, TypeError):                         # pragma: no cover
        return default


def unit_text(unit) -> str:
    folded = text(unit).casefold()
    if folded in {"um", "micron", "microns", "µm"}:
        return "µm"
    if folded in {"deg", "degree", "degrees", "°"}:
        return ANGLE_UNIT
    return text(unit)


def number(f: h5py.File, path: str) -> Optional[float]:
    """The first value of a numeric dataset, or None."""
    if path not in f or not isinstance(f[path], h5py.Dataset):
        return None
    try:
        values = np.asarray(f[path][()], dtype=float).reshape(-1)
    except (TypeError, ValueError, OSError):
        return None
    return float(values[0]) if values.size and np.isfinite(values[0]) \
        else None


_SPATIAL = {"sax", "say", "saz", "smx", "smy", "smz", "x", "y", "z"}


def motor_role(name: str) -> tuple:
    """``(role, label)`` for a scanned motor, from its name.

    ``role`` is what :data:`loader.registry.AXIS_ROLES` keys on for a map
    (``angle`` is the one a k conversion makes sense for), or ``x``/``y``
    style ``position`` for a spatial motor."""
    folded = name.casefold()
    if "defl" in folded:
        return "angle", f"Deflector angle ({ANGLE_UNIT})"
    if "polar" in folded:
        return "angle", f"Polar angle ({ANGLE_UNIT})"
    if "tilt" in folded:
        return "angle", f"Tilt angle ({ANGLE_UNIT})"
    if "azimuth" in folded:
        return "other", f"Azimuth ({ANGLE_UNIT})"
    if folded in _SPATIAL:
        return "position", None
    if "photon" in folded or folded in {"energy", "hv", "value"}:
        return "photon_energy", f"Photon energy ({ENERGY_UNIT})"
    if "temperature" in folded or "cryostat" in folded:
        return "temperature", "Temperature (K)"
    return "other", None


# --------------------------------------------------------------------------
# The layout a beamline module works out
# --------------------------------------------------------------------------
@dataclass
class Dim:
    """One scanned dimension of the spectra."""
    name: str
    unit: str
    axis: np.ndarray
    role: str                   # an AXIS_ROLES key, or "position"/"index"
    label: str = ""

    def __post_init__(self):
        self.axis = np.asarray(self.axis, dtype=float).reshape(-1)
        if not self.label:
            unit = unit_text(self.unit)
            self.label = f"{self.name} ({unit})" if unit else self.name


@dataclass
class Layout:
    """Where the spectra are and what their axes mean.

    ``perm`` puts the dataset's dimensions in the order
    ``(*dims, angle, energy)``: ``perm[i]`` is the dataset axis that becomes
    axis ``i``. ``energies`` is set when every scan point has its own
    kinetic-energy axis (a photon-energy scan), shaped ``(points, energy)``;
    ``energy_scale`` says what ``energy`` itself is -- ``"kinetic"``, or
    ``"binding"`` for an axis the beamline already referred to E_F.
    """
    dataset: h5py.Dataset
    perm: tuple
    dims: list
    angle: np.ndarray
    energy: np.ndarray
    facility: str
    prefix: str
    energy_scale: str = "kinetic"
    energies: Optional[np.ndarray] = None
    info: dict = field(default_factory=dict)


# --------------------------------------------------------------------------
# Reading the spectra lazily
# --------------------------------------------------------------------------
def _is_int(key) -> bool:
    return isinstance(key, (int, np.integer))


class NexusView:
    """A NeXus spectra dataset as an array that stays on disk.

    Presented with its axes in ``perm`` order and with the dataset axes in
    ``fixed`` pinned to one index (that is how one step of an extra scan
    axis becomes a dataset of its own, and how a scan axis of length one
    disappears from a cut). Indexing works as on the numpy array and reads
    only the selection; ``np.asarray`` reads the lot. Integer counts are
    returned as float32, so a sum or a normalisation cannot overflow.
    """

    def __init__(self, dataset, perm, fixed: dict = None):
        self._dset = dataset
        self._fixed = {int(k): int(v) for k, v in (fixed or {}).items()}
        self._perm = [int(p) for p in perm if int(p) not in self._fixed]
        self.shape = tuple(int(dataset.shape[p]) for p in self._perm)
        self.ndim = len(self.shape)
        dtype = np.dtype(dataset.dtype)
        self.dtype = dtype if np.issubdtype(dtype, np.floating) \
            else np.dtype(np.float32)
        self.size = int(np.prod(self.shape)) if self.shape else 1
        self.nbytes = self.size * self.dtype.itemsize

    def __len__(self):
        return self.shape[0] if self.shape else 0

    def __repr__(self):
        return f"NexusView(shape={self.shape}, dtype={self.dtype})"

    def _normalise(self, key):
        if not isinstance(key, tuple):
            key = (key,)
        positions = [i for i, k in enumerate(key) if k is Ellipsis]
        if positions:
            i = positions[0]
            key = key[:i] + (slice(None),) * (self.ndim - len(key) + 1) \
                + key[i + 1:]
        if len(key) > self.ndim:
            raise IndexError(f"too many indices for a {self.ndim}-D array")
        return key + (slice(None),) * (self.ndim - len(key))

    def __getitem__(self, key):
        key = self._normalise(key)
        source = [slice(None)] * self._dset.ndim
        for axis, index in self._fixed.items():
            source[axis] = index
        later = {}                      # out axis -> what h5py was not asked
        for out_axis, sel in enumerate(key):
            src = self._perm[out_axis]
            if _is_int(sel):
                sel = int(sel)
                source[src] = sel if sel >= 0 else self.shape[out_axis] + sel
            elif isinstance(sel, slice) and (sel.step is None or sel.step > 0):
                source[src] = sel
            else:
                later[out_axis] = sel
        data = np.asarray(self._dset[tuple(source)])

        surviving = [a for a in range(self._dset.ndim) if not _is_int(source[a])]
        wanted = [self._perm[a] for a, sel in enumerate(key) if not _is_int(sel)]
        if len(wanted) > 1:
            data = np.transpose(data, [surviving.index(a) for a in wanted])
        kept = [a for a, sel in enumerate(key) if not _is_int(sel)]
        for out_axis, sel in later.items():
            index = [slice(None)] * data.ndim
            index[kept.index(out_axis)] = sel if isinstance(sel, slice) \
                else np.asarray(sel)
            data = data[tuple(index)]
        return data.astype(self.dtype, copy=False)

    def materialise(self) -> np.ndarray:
        return self[(slice(None),) * self.ndim]

    def __array__(self, dtype=None, copy=None):
        data = self.materialise()
        return data if dtype is None else data.astype(dtype)


# --------------------------------------------------------------------------
# From a layout to the datasets the browser lists
# --------------------------------------------------------------------------
def plan_entries(layout: Layout) -> tuple:
    """``(kind, shown dims, pinned dims)`` for ``layout``.

    Dimensions of length one are pinned silently -- a Diamond cut is a scan
    of one point -- and do not count towards the listing."""
    dims = layout.dims
    real = [i for i, d in enumerate(dims) if len(d.axis) > 1]
    trivial = [i for i in range(len(dims)) if i not in real]
    roles = [dims[i].role for i in real]
    if not real:
        return "cut", [], trivial
    spatial = [i for i in real if dims[i].role == "position"]
    if len(spatial) >= 2:
        shown = spatial[-2:]
        kind = "spem_4d"
    elif len(real) == 1 and roles[0] == "photon_energy":
        shown, kind = real, "kz_map"
    else:
        shown, kind = [real[-1]], "map"
        if dims[real[-1]].role == "photon_energy":
            kind = "kz_map"
    pinned = [i for i in range(len(dims)) if i not in shown]
    return kind, shown, pinned


def _listed_pins(layout: Layout, pinned: list) -> list:
    """The pinned dims that actually have more than one step."""
    return [i for i in pinned if len(layout.dims[i].axis) > 1]


def _entry_key(indices) -> str:
    return "point " + ",".join(str(i) for i in indices)


def _entry_label(layout: Layout, pins: list, indices) -> str:
    parts = []
    for dim_index, i in zip(pins, indices):
        dim = layout.dims[dim_index]
        unit = unit_text(dim.unit)
        parts.append(f"{dim.name} = {dim.axis[i]:.4g}"
                     + (f" {unit}" if unit else ""))
    return ", ".join(parts)


def list_layout(layout: Layout, base: str, title: str = None,
                start_time: str = None) -> list:
    kind, _shown, pinned = plan_entries(layout)
    pins = _listed_pins(layout, pinned)
    label = KIND_LABELS.get(kind, kind)
    if not pins:
        return [{"entry": "entry1", "kind": label, "name": base,
                 "title": title, "start_time": start_time}]
    combos = product(*[range(len(layout.dims[i].axis)) for i in pins])
    return [{"entry": _entry_key(indices), "kind": label,
             "name": f"{base} {_entry_label(layout, pins, indices)}",
             "title": title, "start_time": start_time}
            for indices in combos]


def pinned_indices(layout: Layout, entry: str) -> dict:
    """``{dim: index}`` for every pinned dim, from an entry key."""
    _kind, _shown, pinned = plan_entries(layout)
    pins = _listed_pins(layout, pinned)
    fixed = {i: 0 for i in pinned}
    if pins:
        if entry and entry.startswith("point "):
            indices = tuple(int(i) for i in entry[6:].split(","))
        else:
            indices = (0,) * len(pins)
        if len(indices) != len(pins) or any(
                not 0 <= i < len(layout.dims[d].axis)
                for d, i in zip(pins, indices)):
            raise KeyError(f"no dataset {entry!r}")
        fixed.update(zip(pins, indices))
    return fixed


# --------------------------------------------------------------------------
# Photon-energy scans
# --------------------------------------------------------------------------
def straighten(hv: np.ndarray) -> np.ndarray:
    """A monochromator readback that is plainly a jittered even scan, put on
    the even scan it was asked to be -- the viewers place rows evenly."""
    hv = np.asarray(hv, dtype=float)
    if hv.size <= 2:
        return hv
    index = np.arange(hv.size, dtype=float)
    slope, intercept = np.polyfit(index, hv, 1)
    straight = index * slope + intercept
    if np.max(np.abs(straight - hv)) < 0.25 * abs(slope):
        return straight
    return hv


def _ascending(axis: np.ndarray, cube: np.ndarray):
    if axis.size > 1 and axis[0] > axis[-1]:
        return axis[::-1].copy(), cube[..., ::-1]
    return axis, cube


def _refer_kz(layout: Layout, cube: np.ndarray, hv: np.ndarray, notes: dict):
    """``(cube, energy, label)``: the spectra of a photon-energy scan on one
    energy axis referred to E_F. ``cube`` is ``(hv, angle, energy)``."""
    energy = np.asarray(layout.energy, dtype=float)
    kinetic = layout.energies
    if layout.energy_scale == "binding":
        # Already referred by the beamline. A positive binding energy grows
        # as the kinetic energy falls; turn that into E - E_F.
        if kinetic is not None and kinetic.shape[-1] == energy.size \
                and energy.size > 1:
            same_way = np.sign(kinetic[0, -1] - kinetic[0, 0]) \
                == np.sign(energy[-1] - energy[0])
        else:
            same_way = np.median(energy) <= 0
        referred = energy if same_way else -energy
        if not same_way:
            notes["energy_reference"] = "binding energy in the file, negated"
        else:
            notes["energy_reference"] = "E - E_F as recorded in the file"
        if kinetic is not None and kinetic.shape == (hv.size, energy.size):
            w = np.median(referred[None, :] - kinetic + hv[:, None])
            if np.isfinite(w) and 0.0 < w < 10.0:
                notes["work_function_eV"] = float(w)
                notes["work_function_source"] = "from the file's two axes"
        notes.setdefault("work_function_eV", WORK_FUNCTION_EV)
        referred, cube = _ascending(referred, cube)
        return cube, referred, f"E - E_F ({ENERGY_UNIT})"

    if kinetic is None:
        kinetic = np.broadcast_to(energy, (hv.size, energy.size))
    w = WORK_FUNCTION_EV
    notes["work_function_eV"] = w
    notes["work_function_source"] = "assumed (not recorded)"
    notes["energy_reference"] = (
        f"each spectrum referred by its own hv: E_kin - hv + {w:g} eV")
    referred = [np.asarray(row, dtype=float) - h + w
                for row, h in zip(kinetic, hv)]
    first = referred[0]
    step = abs(first[1] - first[0]) if first.size > 1 else 1.0
    if all(np.max(np.abs(axis - first)) < 0.01 * step for axis in referred):
        axis, cube = _ascending(first, cube)
        return cube, axis, f"E - E_F ({ENERGY_UNIT})"
    low = max(float(axis.min()) for axis in referred)
    high = min(float(axis.max()) for axis in referred)
    if high <= low:
        raise ValueError(
            "the energy windows of this photon-energy scan do not overlap "
            "once referred to the Fermi level")
    common = low + step * np.arange(int(np.floor((high - low) / step
                                                 + 1e-9)) + 1)
    out = np.empty(cube.shape[:2] + (common.size,), dtype=np.float32)
    for index, axis in enumerate(referred):
        frame = cube[index]
        if axis[0] > axis[-1]:
            axis, frame = axis[::-1], frame[:, ::-1]
        out[index] = _resample_energy(frame, axis, common)
    notes["energy_resampled"] = (
        f"spectra resampled onto their common range, {common[0]:.4g} to "
        f"{common[-1]:.4g} eV")
    return out, common, f"E - E_F ({ENERGY_UNIT})"


# --------------------------------------------------------------------------
# Building the scan
# --------------------------------------------------------------------------
def energy_label(layout: Layout) -> str:
    if layout.energy_scale == "binding":
        return f"E - E_F ({ENERGY_UNIT})"
    return f"Kinetic energy ({ENERGY_UNIT})"


def build_scan(layout: Layout, entry: str = None) -> NxsScan:
    """The :class:`NxsScan` for one listed dataset of ``layout``."""
    kind, shown, pinned = plan_entries(layout)
    fixed_dims = pinned_indices(layout, entry)
    n_dims = len(layout.dims)
    # The order the shown dims come out of NexusView: the layout's, with the
    # pinned ones gone. A spatial scan wants (y, x), which is that order
    # when the slower motor is listed first, as both beamlines write it.
    out_perm = [layout.perm[i] for i in shown] \
        + [layout.perm[n_dims], layout.perm[n_dims + 1]]
    fixed = {layout.perm[d]: i for d, i in fixed_dims.items()}
    view = NexusView(layout.dataset, out_perm, fixed)

    angle_label = f"Angle along slit ({ANGLE_UNIT})"
    scan = NxsScan(kind=kind, filename_prefix=layout.prefix)
    scan.info.update(layout.info)
    scan.info["facility"] = layout.facility
    scan.info["nexus.dataset"] = layout.dataset.name
    for d, i in fixed_dims.items():
        dim = layout.dims[d]
        scan.fourd_info[dim.name] = \
            f"{dim.axis[i]:.6g} {unit_text(dim.unit)}".strip()

    if kind == "cut":
        scan.x, scan.y = layout.angle, layout.energy
        scan.value = view
        scan.labels = {"x": angle_label, "y": energy_label(layout)}
    elif kind == "map":
        dim = layout.dims[shown[0]]
        scan.x, scan.k, scan.z = dim.axis, layout.angle, layout.energy
        scan.value = view
        scan.labels = {"x": dim.label, "k": angle_label,
                       "z": energy_label(layout)}
        scan.info["axis0.role"] = dim.role if dim.role in (
            "angle", "photon_energy", "temperature", "position") else "other"
    elif kind == "kz_map":
        dim = layout.dims[shown[0]]
        hv = dim.axis
        point = fixed_dims
        energies = layout.energies
        if energies is not None and energies.ndim == n_dims + 1:
            # One kinetic axis per scan point: keep the ones of this slice.
            energies = energies[tuple(point.get(i, slice(None))
                                      for i in range(n_dims))]
        sliced = Layout(**{**layout.__dict__, "energies": energies})
        cube = np.asarray(view, dtype=np.float32)
        order = np.argsort(hv, kind="stable")
        hv, cube = hv[order], cube[order]
        if sliced.energies is not None:
            sliced.energies = np.asarray(sliced.energies)[order]
        notes = {}
        try:
            cube, energy, e_label = _refer_kz(sliced, cube, hv, notes)
        except ValueError as exc:
            # No common window: keep the spectra as they are, as a map
            # whose first axis is a photon energy.
            scan.kind = "map"
            scan.x, scan.k, scan.z = hv, layout.angle, layout.energy
            scan.value = cube
            scan.labels = {"x": dim.label, "k": angle_label,
                           "z": energy_label(layout)}
            scan.info["axis0.role"] = "photon_energy"
            scan.info["energy_reference"] = f"not referred: {exc}"
            return scan
        stepped = straighten(hv)
        scan.x, scan.k, scan.z = stepped, layout.angle, energy
        scan.value = cube
        scan.labels = {"x": f"Photon energy ({ENERGY_UNIT})",
                       "k": angle_label, "z": e_label}
        scan.info.update(notes)
        scan.info["axis0.role"] = "photon_energy"
        scan.info["axis0_measured"] = hv
        scan.info["photon_energy_range_eV"] = f"{hv.min():g} to {hv.max():g}"
        scan.info.pop("photon_energy_eV", None)
    elif kind == "spem_4d":
        y_dim, x_dim = layout.dims[shown[0]], layout.dims[shown[1]]
        scan.y, scan.x = y_dim.axis, x_dim.axis
        scan.k, scan.z = layout.angle, layout.energy
        scan.value4d = view
        scan.labels = {"x": x_dim.label, "y": y_dim.label, "k": angle_label,
                       "z": energy_label(layout)}
    else:                                                  # pragma: no cover
        raise AssertionError(kind)

    for dim in layout.dims:
        if len(dim.axis) > 1 and dim.role != "index":
            scan.fourd_info.setdefault(
                f"{dim.name} range",
                f"{dim.axis[0]:.6g} .. {dim.axis[-1]:.6g} "
                f"{unit_text(dim.unit)} ({len(dim.axis)} steps)".strip())
    return scan


# --------------------------------------------------------------------------
# Metadata
# --------------------------------------------------------------------------
def read_info(f: h5py.File, entry: str, groups: dict) -> dict:
    """The metadata worth showing: the entry's provenance, and every scalar
    of each group in ``groups`` (``{path under the entry: prefix}``)."""
    info = {}
    for key in ("title", "start_time", "end_time", "scan_command",
                "experiment_identifier", "entry_identifier", "notes"):
        path = f"{entry}/{key}"
        if path in f and isinstance(f[path], h5py.Dataset):
            value = text(f[path][()])
            if value:
                info[key] = value
    for path, prefix in groups.items():
        _read_group_fields(f, f"{entry}/{path}", prefix, info)
    for key, value in list(info.items()):
        if isinstance(value, np.ndarray) and value.size == 0:
            del info[key]
        elif isinstance(value, (bytes, np.bytes_)):
            info[key] = text(value)
        else:
            info[key] = _scalarize(value)
    return info


# --------------------------------------------------------------------------
# The loader, given a way to read a layout
# --------------------------------------------------------------------------
class NexusArpesLoader:
    """Mixed into a :class:`loader.registry.Loader`: ``read_layout(f)`` and
    ``recognises(f)`` are the beamline's, the rest is shared."""

    def recognises(self, f: h5py.File) -> bool:            # pragma: no cover
        raise NotImplementedError

    def read_layout(self, f: h5py.File) -> Layout:         # pragma: no cover
        raise NotImplementedError

    def can_open(self, path: str) -> bool:
        try:
            with HANDLES.borrow(path) as f:
                return bool(self.recognises(f))
        except Exception:
            return False

    def list_entries(self, path: str) -> list:
        import os
        base = os.path.splitext(os.path.basename(path))[0]
        try:
            with HANDLES.borrow(path) as f:
                layout = self.read_layout(f)
                return list_layout(layout, base, layout.info.get("title"),
                                   layout.info.get("start_time"))
        except Exception:
            return []

    def load(self, path: str, entry: str = None, progress=None):
        f = HANDLES.acquire(path)
        keep_open = False
        try:
            layout = self.read_layout(f)
            try:
                scan = build_scan(layout, entry)
            except KeyError:
                raise KeyError(f"no dataset {entry!r} in {path}") from None
            if isinstance(scan.array(), NexusView):
                scan._h5file = f
                scan._h5path = path
                keep_open = True
            return scan
        finally:
            if not keep_open:
                HANDLES.release(path, f)
