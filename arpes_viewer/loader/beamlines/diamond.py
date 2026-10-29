"""
loader/beamlines/diamond.py
===========================
Diamond Light Source's two ARPES branches -- **I05** (HR-ARPES) and
**I05-1** (nano-ARPES) -- write the same GDA NeXus dialect, so one reader
serves both; each is registered as a loader of its own so the loader
window names the right one and its metadata says which.

The layout, from ``h5info`` dumps of a cut, a deflector map, a
photon-energy scan and a spatial scan from each branch (``for
loader/DLS_HR.mat``, ``for loader/DLS_Nano.mat``):

* **The spectra** are ``/entry1/analyser/data`` (an ``NXdata`` group,
  ``@signal`` on the dataset), shaped ``(*scan, angles, energies)`` in C
  order. A single cut is a scan of one point: ``(1, angles, energies)``.
* **Each vector says which dimension it belongs to** in ``@axis``, 1-based
  and in the same order: ``angles`` is ``'2'`` in a cut and ``'3'`` in a
  spatial scan, and a motor of a two-dimensional scan is ``'1,2'`` and
  stored on the whole grid, ``(slow, fast)``. Which of the two it steps
  along is read off its values.
* **What was scanned**, by branch:

  ====================  ===================================  ===============
  scan                  I05                                  I05-1
  ====================  ===================================  ===============
  deflector map         ``deflector_x``                      ``deflector_x``
  photon energy         ``energy`` (and ``centre_energy``)   ``value`` (and
                                                             ``centre_energy``)
  spatial               ``sax``, ``saz`` (mm)                ``smx``, ``smy``
                                                             (µm)
  ====================  ===================================  ===============

  A photon-energy scan also stores ``energies`` per point, ``(hv, energy)``:
  the kinetic window moves with ``hv``. ``current`` (the ring current) is
  recorded against every scan and is never the scanned axis.
* **Which branch.** ``/entry1/instrument/name`` when it says, else the
  structure: I05-1 has a ``zone_plate`` and its manipulator's motors are
  ``sm*``; I05's are ``sa*``.
"""
from __future__ import annotations

import os

import h5py
import numpy as np

from arpes_viewer.loader.formats.nexus_arpes import (Dim, Layout, NexusArpesLoader, attr,
                                        motor_role, number, read_info, text)
from arpes_viewer.loader.registry import Loader, register

__all__ = ["DiamondI05Loader", "DiamondI051Loader", "read_layout",
           "is_diamond", "branch"]

_ROOT = "/entry1"
_DATA = "/entry1/analyser"
#: Recorded against the scan but never what was scanned.
_NOT_AXES = {"data", "angles", "energies", "current", "cps"}
#: When several vectors belong to one dimension, these are the scanned
#: photon energy (``energy`` at I05, ``value`` at I05-1); ``centre_energy``
#: is the analyser following it.
_HV_NAMES = ("energy", "value", "photon_energy", "hv")

_INFO_GROUPS = {
    "instrument/analyser": "analyser",
    "instrument/monochromator": "mono",
    "instrument/monochromator/grating": "grating",
    "instrument/insertion_device": "id",
    "instrument/insertion_device/beam": "beam",
    "instrument/manipulator": "manipulator",
    "instrument/source": "source",
    "instrument/zone_plate": "zone_plate",
    "instrument/order_sorting_aperture": "osa",
    "sample": "sample",
    "user01": "user",
}


def is_diamond(f: h5py.File) -> bool:
    group = f.get(_DATA)
    return isinstance(group, h5py.Group) and all(
        isinstance(group.get(name), h5py.Dataset)
        for name in ("data", "angles", "energies")) \
        and f"{_ROOT}/instrument" in f


def branch(f: h5py.File, path: str = "") -> str:
    """``"I05-1"`` or ``"I05"``."""
    name = ""
    if f"{_ROOT}/instrument/name" in f:
        try:
            name = text(f[f"{_ROOT}/instrument/name"][()]).casefold()
        except (OSError, TypeError):
            name = ""
    if "i05-1" in name or "i05_1" in name:
        return "I05-1"
    if name.startswith("i05"):
        return "I05"
    if os.path.basename(path).casefold().startswith("i05-1-"):
        return "I05-1"
    manipulator = f.get(f"{_ROOT}/instrument/manipulator")
    if f"{_ROOT}/instrument/zone_plate" in f or (
            isinstance(manipulator, h5py.Group)
            and any(k.startswith("sm") for k in manipulator.keys())):
        return "I05-1"
    return "I05"


def _axis_numbers(item) -> list:
    """``@axis`` as 0-based dimensions: ``'1,2'`` -> ``[0, 1]``."""
    value = attr(item, "axis")
    if value is None:
        return []
    out = []
    for part in text(value).replace(" ", "").split(","):
        try:
            out.append(int(float(part)) - 1)
        except ValueError:
            continue
    return out


def _is_primary(item) -> bool:
    try:
        return int(float(text(attr(item, "primary", "0")))) == 1
    except ValueError:
        return False


def _stepped_dim(values: np.ndarray):
    """Which axis of a motor stored on the whole grid it steps along, and
    its positions there; ``None`` if it does not move."""
    spreads = []
    for axis in range(values.ndim):
        spreads.append(float(np.mean(np.ptp(values, axis=axis))))
    best = int(np.argmax(spreads))
    scale = max(1.0, float(np.max(np.abs(values))))
    if spreads[best] <= 1e-9 * scale:
        return None
    others = tuple(a for a in range(values.ndim) if a != best)
    return best, values.mean(axis=others) if others else values


def read_layout(f: h5py.File, path: str = "") -> Layout:
    if not is_diamond(f):
        raise ValueError("not a Diamond I05 / I05-1 NeXus file")
    group = f[_DATA]
    dataset = group["data"]
    shape = tuple(int(n) for n in dataset.shape)
    if len(shape) < 2:
        raise ValueError(f"spectra of shape {shape} have no detector plane")
    angle = np.asarray(group["angles"][()], dtype=float).reshape(-1)
    energies = np.asarray(group["energies"][()], dtype=float)
    n_energy = energies.shape[-1] if energies.ndim else 1

    angle_dim = _axis_numbers(group["angles"])
    angle_dim = angle_dim[0] if angle_dim and angle_dim[0] < len(shape) \
        and shape[angle_dim[0]] == angle.size else len(shape) - 2
    energy_dim = _axis_numbers(group["energies"])
    energy_dim = energy_dim[-1] if energy_dim and energy_dim[-1] < len(shape) \
        and shape[energy_dim[-1]] == n_energy else len(shape) - 1
    if shape[angle_dim] != angle.size or shape[energy_dim] != n_energy:
        raise ValueError(
            f"spectra {shape} do not match {angle.size} angles and "
            f"{n_energy} energies")
    scan_axes = [i for i in range(len(shape)) if i not in (angle_dim,
                                                            energy_dim)]
    scan_shape = tuple(shape[i] for i in scan_axes)
    per_point = energies.ndim > 1 and energies.shape[:-1] == scan_shape
    energy = energies.reshape(-1, n_energy)[0] if per_point \
        else energies.reshape(-1)

    # Every vector that could be a scanned axis, per scan dimension.
    candidates = {i: [] for i in range(len(scan_axes))}
    constants = {}
    for name in group.keys():
        item = group[name]
        if name in _NOT_AXES or not isinstance(item, h5py.Dataset) \
                or not np.issubdtype(item.dtype, np.number):
            continue
        values = np.asarray(item[()], dtype=float)
        unit = text(attr(item, "units", ""))
        if values.shape == scan_shape and values.ndim > 1:
            stepped = _stepped_dim(values)
            if stepped is None:
                constants[name] = float(values.reshape(-1)[0])
                continue
            dim, axis = stepped
        else:
            values = values.reshape(-1)
            numbers = [scan_axes.index(a) for a in _axis_numbers(item)
                       if a in scan_axes]
            dim = numbers[0] if len(numbers) == 1 else next(
                (i for i, n in enumerate(scan_shape) if n == values.size
                 and len(scan_shape) == 1), None)
            if dim is None or scan_shape[dim] != values.size:
                continue
            axis = values
        candidates[dim].append((name, unit, axis, _is_primary(item)))

    dims, followers = [], []
    for i, n in enumerate(scan_shape):
        found = candidates[i]
        hv_scan = False
        if per_point and n > 1:
            # The kinetic window moves along this dimension: it is the
            # photon energy that was stepped.
            hv_scan = float(np.ptp(energies[..., 0], axis=i).max()) > 0
        pick = None
        if hv_scan:
            pick = next((c for name in _HV_NAMES for c in found
                         if c[0] == name), None)
        if pick is None:
            ranked = sorted(found, key=lambda c: (
                not c[3], c[0] == "centre_energy", c[0] in _HV_NAMES))
            pick = ranked[0] if ranked else None
        if pick is None:
            dims.append(Dim("Scan point", "", np.arange(n, dtype=float),
                            "index", "Scan point (index)"))
            continue
        name, unit, axis, _primary = pick
        for other, other_unit, other_axis, _ in found:
            if other != name and other_axis.size > 1 \
                    and float(np.ptp(other_axis)) > 0:
                followers.append(Dim(other, other_unit, other_axis, "position"))
        role, label = motor_role(name)
        if hv_scan:
            role, label = "photon_energy", "Photon energy (eV)"
        elif role == "photon_energy" and name == "value":
            role, label = "other", None
        dims.append(Dim(name, unit, axis, role, label or ""))

    info = read_info(f, _ROOT, _INFO_GROUPS)
    info.update({f"analyser.{k}": v for k, v in constants.items()})
    mono = f"{_ROOT}/instrument/monochromator/energy"
    if mono in f and f[mono].size == 1:
        info["photon_energy_eV"] = number(f, mono)
    temperature = number(f, f"{_ROOT}/sample/temperature")
    if temperature is not None:
        info["SampleTemperature_K"] = temperature
    which = branch(f, path)
    facility = "Diamond I05-1 (nano-ARPES)" if which == "I05-1" \
        else "Diamond I05 (HR-ARPES)"
    return Layout(dataset=dataset,
                  perm=tuple(scan_axes) + (angle_dim, energy_dim),
                  dims=dims, angle=angle, energy=energy, facility=facility,
                  prefix="I05-1_" if which == "I05-1" else "I05_",
                  energies=energies if per_point else None, info=info,
                  followers=followers)


class _DiamondLoader(NexusArpesLoader, Loader):
    patterns = ("*.nxs",)
    priority = 12
    branch_name = ""

    def recognises(self, f: h5py.File) -> bool:
        return is_diamond(f) and branch(f, f.filename) == self.branch_name

    def read_layout(self, f: h5py.File) -> Layout:
        return read_layout(f, f.filename)


class DiamondI05Loader(_DiamondLoader):
    name = "Diamond I05 HR-ARPES"
    description = ("NeXus .nxs files from the high-resolution ARPES branch "
                   "of I05 at Diamond: cuts, deflector maps, photon-energy "
                   "scans and spatial scans.")
    branch_name = "I05"


class DiamondI051Loader(_DiamondLoader):
    name = "Diamond I05-1 nano-ARPES"
    description = ("NeXus .nxs files from the nano-ARPES branch (I05-1) at "
                   "Diamond: cuts, deflector maps, photon-energy scans and "
                   "spatial scans.")
    branch_name = "I05-1"


register(DiamondI05Loader())
register(DiamondI051Loader())
