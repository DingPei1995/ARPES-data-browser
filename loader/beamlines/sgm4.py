"""
loader/beamlines/sgm4.py
========================
ASTRID2 SGM4 (Aarhus): the HDF5 files of its micro/nano-ARPES end station
(SPECS Phoibos 150 SAL).

Ported from the SGM4 reader in peaks
(<https://github.com/phrgab/peaks>, ``peaks/core/fileIO/loaders/sgm4.py``).
What is kept from it:

* **The scan plan** is ``Entry/Data/ScanDetails``: ``SlowAxis_names``,
  ``_start``, ``_step`` and ``_length`` for the scanned axes, and the same
  four ``FastAxis_*`` for the detector's (``Kinetic Energy`` and
  ``OrdinateRange``, the angle along the slit). Every axis is ``start + i *
  step``.
* **The spectra** are ``Entry/Data/TransformedData``, one row per scan
  point, ``(points, *fast axes)``. The points run with the *first* slow
  axis fastest, so the grid is ``SlowAxis_length`` reversed. A scan stopped
  early has fewer rows than the grid; the missing points are NaN.
* **A knife-edge scan** has no fast axes, only a summed intensity per point
  (``Entry/Process/SumData``). It is listed, but has no spectra to open.
* **Axis names**: ``ShiftX``/``ShiftY`` are the deflectors (``SALX``/
  ``SALY`` in older files), ``FSamX/Y/Z`` the sample stage (``SamX/Y/Z`` in
  a knife-edge scan) in µm, ``CapZ`` the capillary focus, ``LEG``/``MEG``/
  ``HEG`` the photon energy of the grating in use, and ``Loop`` a repeat.
* **A photon-energy scan** keeps one binding window: each hv's kinetic axis
  is the first one moved by the change in hv. The analyser work function is
  4.416 eV.
* **Metadata**: ``Entry/Instrument/Detector/Info`` holds the analyser
  settings as ``Key: value`` lines; the monochromator, temperatures,
  positioners and Keithleys are under ``Entry/Instrument``. Polar angles
  are also given as peaks computes them (sample polar ``/1000 - 45``, the
  analyser's ``+ 14``).

How the scanned axes become a kind is shared with the NeXus readers
(:mod:`loader.formats.nexus_arpes`): one deflector gives a map that converts to k,
the photon energy a kz map, two stage axes a spatial scan; anything else is
listed one dataset per extra step. Everything but a kz map is read lazily.

Untested on real SGM4 files: none were available when it was written.

Nothing in here imports Qt.
"""
from __future__ import annotations

import os

import h5py
import numpy as np

from loader.formats.nexus_arpes import (Dim, Layout, NexusArpesLoader,
                                        motor_role, text)
from loader.nxs_file import ANGLE_UNIT, ENERGY_UNIT, HANDLES
from loader.registry import Loader, register

__all__ = ["Sgm4Loader", "is_sgm4", "read_layout", "FlatGrid",
           "WORK_FUNCTION_EV"]

WORK_FUNCTION_EV = 4.416
_DETAILS = "Entry/Data/ScanDetails"
_SPECTRA = "Entry/Data/TransformedData"
_LEGACY = {"SALX": "ShiftX", "SALY": "ShiftY"}


def is_sgm4(f: h5py.File) -> bool:
    return isinstance(f.get(_DETAILS), h5py.Group) \
        and "SlowAxis_names" in f[_DETAILS] \
        and (_SPECTRA in f or "Entry/Process/SumData" in f)


def is_knife_edge(f: h5py.File) -> bool:
    return "FastAxis_start" not in f[_DETAILS]


def _names(values) -> list:
    names = [text(v) for v in np.atleast_1d(values)]
    return [_LEGACY.get(n, n) for n in names]


def _axes(details: h5py.Group, which: str) -> tuple:
    """``(names, lengths, axes)`` of the slow or fast axes."""
    names = _names(details[f"{which}_names"][()])
    starts = np.atleast_1d(details[f"{which}_start"][()]).astype(float)
    steps = np.atleast_1d(details[f"{which}_step"][()]).astype(float)
    lengths = [int(n) for n in np.atleast_1d(details[f"{which}_length"][()])]
    axes = [start + step * np.arange(n)
            for start, step, n in zip(starts, steps, lengths)]
    return names, lengths, axes


def _dim(name: str, axis) -> Dim:
    folded = name.casefold()
    if folded in {"shiftx", "shifty"}:
        return Dim(name, "deg", axis, "angle",
                   f"Deflector {name} ({ANGLE_UNIT})")
    if folded in {"leg", "meg", "heg", "hv"}:
        return Dim(name, "eV", axis, "photon_energy",
                   f"Photon energy ({ENERGY_UNIT})")
    if folded in {"fsamx", "fsamy", "fsamz", "samx", "samy", "samz"}:
        return Dim(name, "µm", axis, "position", f"{name[-1]} (µm)")
    if folded == "capz":
        return Dim(name, "µm", axis, "other", "Defocus, CapZ (µm)")
    if folded == "loop":
        return Dim(name, "", axis, "other", "Repeat (index)")
    if folded == "sampolar":
        return Dim(name, "deg", axis, "angle", f"Polar angle ({ANGLE_UNIT})")
    if folded == "samtilt":
        return Dim(name, "deg", axis, "angle", f"Tilt angle ({ANGLE_UNIT})")
    role, label = motor_role(name)
    return Dim(name, "", axis, role, label or "")


# --------------------------------------------------------------------------
# One row per point, presented as a grid
# --------------------------------------------------------------------------
class FlatGrid:
    """``TransformedData`` -- ``(points, *fast)`` -- as the array
    ``(*grid, angle, energy)``, read only where it is indexed.

    ``grid`` is the slow axes' lengths, slowest first; point ``p`` of the
    file is grid index ``unravel_index(p, grid)``. Points the file does not
    have (a scan stopped early) read as NaN. ``fast_perm`` puts the file's
    fast axes in ``(angle, energy)`` order."""

    def __init__(self, dataset, grid, fast_perm):
        self._dset = dataset
        self._grid = tuple(int(n) for n in grid)
        self._fast_perm = tuple(int(p) for p in fast_perm)
        self._rows = int(dataset.shape[0])
        fast = tuple(int(dataset.shape[1 + p]) for p in self._fast_perm)
        self.shape = self._grid + fast
        self.ndim = len(self.shape)
        self.dtype = np.dtype(np.float32)
        self.name = dataset.name

    def __getitem__(self, key):
        if not isinstance(key, tuple):
            key = (key,)
        key = key + (slice(None),) * (self.ndim - len(key))
        n_grid = len(self._grid)
        picks, drop = [], []
        for axis, sel in enumerate(key[:n_grid]):
            chosen = np.arange(self._grid[axis])[sel]
            if np.ndim(chosen) == 0:
                drop.append(axis)
            picks.append(np.atleast_1d(chosen))
        flat = np.ravel_multi_index(np.meshgrid(*picks, indexing="ij"),
                                    self._grid) if n_grid else np.zeros(
                                        (), dtype=int)
        flat = np.asarray(flat)
        present = flat < self._rows
        wanted = np.unique(flat[present])
        if wanted.size and wanted[-1] - wanted[0] + 1 == wanted.size:
            rows = self._dset[int(wanted[0]):int(wanted[-1]) + 1]
        elif wanted.size:
            rows = self._dset[wanted.tolist()]
        else:
            rows = np.empty((0,) + tuple(self._dset.shape[1:]))
        rows = np.transpose(np.asarray(rows, dtype=np.float32),
                            (0,) + tuple(1 + p for p in self._fast_perm))
        out = np.full(flat.shape + rows.shape[1:], np.nan, dtype=np.float32)
        out[present] = rows[np.searchsorted(wanted, flat[present])]
        out = out[(slice(None),) * flat.ndim + tuple(key[n_grid:])]
        if drop:
            out = out.reshape(tuple(n for i, n in enumerate(out.shape)
                                    if i not in drop))
        return out


# --------------------------------------------------------------------------
# Metadata
# --------------------------------------------------------------------------
def _read_info(f: h5py.File) -> dict:
    info = {}
    for path, key in (("Entry/Data/Command", "scan_command"),
                      ("Entry/Data/Timestamp", "start_time")):
        if path in f:
            value = np.atleast_1d(f[path][()])
            if value.size:
                info[key] = text(value[0])
    detector_info = f.get("Entry/Instrument/Detector/Info")
    if isinstance(detector_info, h5py.Dataset):
        for line in text(detector_info[()]).splitlines():
            if ":" in line:
                name, value = line.split(":", 1)
                key = "analyser." + name.strip().casefold().replace(" ", "_")
                try:
                    info[key] = float(value)
                except ValueError:
                    info[key] = value.strip()
    instrument = f.get("Entry/Instrument")
    if isinstance(instrument, h5py.Group):
        for group_name, group in instrument.items():
            if not isinstance(group, h5py.Group):
                continue
            prefix = group_name.casefold().replace(" ", "_")
            for name, item in group.items():
                if not isinstance(item, h5py.Dataset) or name == "Info":
                    continue
                key = f"{prefix}.{name.replace(' ', '_')}"
                try:
                    values = np.asarray(item[()], dtype=float).reshape(-1)
                except (TypeError, ValueError):
                    info[key] = text(item[()])
                    continue
                values = values[np.isfinite(values)]
                if values.size:
                    info[key] = float(values.mean())
    hv = f.get("Entry/Instrument/Monochromator/Photon Energy")
    if isinstance(hv, h5py.Dataset) and hv.size:
        info["photon_energy_eV"] = float(np.atleast_1d(hv[()])[0])
    if "temperature.Sample_Stage" in info:
        info["SampleTemperature_K"] = info["temperature.Sample_Stage"]
    # Polar angles as peaks reads them: the sample's in millidegrees about
    # a 45 deg offset, the analyser's about -14 deg.
    if "positioner.SamPolar" in info:
        info["manipulator.polar_deg"] = info["positioner.SamPolar"] / 1000 - 45
    if "detector.AnaPolar" in info:
        info["analyser.polar_deg"] = info["detector.AnaPolar"] + 14
    info["analyser.work_function_eV"] = WORK_FUNCTION_EV
    return info


# --------------------------------------------------------------------------
# The layout
# --------------------------------------------------------------------------
def read_layout(f: h5py.File) -> Layout:
    if not is_sgm4(f):
        raise ValueError("not an SGM4 file")
    if is_knife_edge(f):
        # Not a ValueError: the viewer answers that by retrying an .h5 with
        # the ANTARES reader, whose complaint would hide this one.
        raise NotImplementedError(
            "a knife-edge scan holds a summed intensity per point, not "
            "spectra, and is not read here")
    details = f[_DETAILS]
    slow_names, slow_lengths, slow_axes = _axes(details, "SlowAxis")
    fast_names, _fast_lengths, fast_axes = _axes(details, "FastAxis")
    dataset = f[_SPECTRA]
    if dataset.ndim != 1 + len(fast_names) or len(fast_names) != 2:
        raise ValueError(f"spectra of shape {dataset.shape} do not match "
                         f"the fast axes {fast_names}")
    energy_index = next((i for i, n in enumerate(fast_names)
                         if "energy" in n.casefold()), 0)
    angle_index = 1 - energy_index
    grid = tuple(reversed(slow_lengths))
    if dataset.shape[0] > int(np.prod(grid)):
        raise ValueError(f"{dataset.shape[0]} spectra for a grid of "
                         f"{int(np.prod(grid))} points")

    dims = [_dim(name, axis) for name, axis in zip(reversed(slow_names),
                                                   reversed(slow_axes))]
    info = _read_info(f)
    if dataset.shape[0] < int(np.prod(grid)):
        info["sgm4.incomplete"] = (f"{dataset.shape[0]} of "
                                   f"{int(np.prod(grid))} points recorded; "
                                   f"the rest are NaN")
    energy = fast_axes[energy_index]
    energies = None
    hv_dims = [d for d in dims if d.role == "photon_energy"
               and len(d.axis) > 1]
    if len(hv_dims) == 1:
        hv = hv_dims[0].axis
        energies = energy[None, :] + (hv - hv[0])[:, None]
    return Layout(dataset=FlatGrid(dataset, grid, (angle_index, energy_index)),
                  perm=tuple(range(len(grid) + 2)), dims=dims,
                  angle=fast_axes[angle_index], energy=energy,
                  facility="ASTRID2 SGM4", prefix="SGM4_", energies=energies,
                  info=info, work_function=WORK_FUNCTION_EV)


class Sgm4Loader(NexusArpesLoader, Loader):
    name = "ASTRID2 SGM4"
    description = ("HDF5 files from the SGM4 micro-ARPES end station at "
                   "ASTRID2 (Aarhus): cuts, deflector maps, photon-energy "
                   "scans and spatial scans.")
    patterns = ("*.h5", "*.hdf5")
    priority = 13

    def recognises(self, f: h5py.File) -> bool:
        return is_sgm4(f)

    def read_layout(self, f: h5py.File) -> Layout:
        return read_layout(f)

    def list_entries(self, path: str) -> list:
        try:
            with HANDLES.borrow(path) as f:
                knife_edge = is_sgm4(f) and is_knife_edge(f)
        except Exception:
            return []
        if knife_edge:
            base = os.path.splitext(os.path.basename(path))[0]
            return [{"entry": "entry1", "kind": "Knife edge", "name": base,
                     "title": "knife-edge scan (no spectra; not read here)",
                     "start_time": None}]
        return super().list_entries(path)


register(Sgm4Loader())
