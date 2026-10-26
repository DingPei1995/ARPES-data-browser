"""
loader/beamlines/alba.py
========================
ALBA's ARPES beamline (BL20, LOREA): the NeXus ``.nxs`` files its
acquisition writes.

The layout, from ``h5info`` dumps of a cut, a deflector map, a
photon-energy scan and a spatial scan (``for loader/ALBA.mat``):

* **One entry, ``/entry1``**, named by the root's ``@default``; its own
  ``@default`` names the ``NXdata`` group, ``data``.
* **The spectra** are ``data/data``. The group's ``@signal`` names them and
  its ``@axes`` lists the vector of each dimension, in the array's order:

  =======================  ============================================
  ``@axes``                 measurement
  =======================  ============================================
  angles, energies          one cut
  defl_angles, angles,      deflector map
  energies
  photon_energies, angles,  photon-energy scan; ``binding_energies`` is
  binding_energies          already referred to E_F, and the kinetic
                            axis of each photon energy is
                            ``instrument/analyser/energies``
                            ``(hv, energy)``
  saz, sax, angles,         spatial scan, the slower motor first
  energies
  =======================  ============================================

* **Metadata** is in ``instrument/{analyser, monochromator,
  insertion_device, manipulator, source}`` and ``sample``.

Anything else the beamline scans (a polar angle, a single motor) is read the
same way, from ``@axes``: the analyser's ``angles`` and the energy are found
by name and the rest are scanned axes. What each becomes is decided in
:mod:`loader.formats.nexus_arpes`.
"""
from __future__ import annotations

import h5py
import numpy as np

from loader.formats.nexus_arpes import (Dim, Layout, NexusArpesLoader, attr,
                                        motor_role, number, read_info, text)
from loader.registry import Loader, register

__all__ = ["AlbaLoader", "read_layout", "is_alba"]

_ENERGY_NAMES = ("energies", "binding_energies", "kinetic_energies",
                 "energy")
_INFO_GROUPS = {
    "instrument": "instrument",
    "instrument/analyser": "analyser",
    "instrument/monochromator": "mono",
    "instrument/insertion_device": "id",
    "instrument/insertion_device/beam": "beam",
    "instrument/manipulator": "manipulator",
    "instrument/source": "source",
    "sample": "sample",
}


def _entry(f: h5py.File):
    name = text(attr(f, "default", "")) or "entry1"
    if name in f and isinstance(f[name], h5py.Group):
        return "/" + name
    for key in f.keys():
        if isinstance(f[key], h5py.Group) \
                and text(attr(f[key], "NX_class", "")) == "NXentry":
            return "/" + key
    return None


def _nxdata(f: h5py.File, entry: str):
    name = text(attr(f[entry], "default", "")) or "data"
    path = f"{entry}/{name}"
    if path in f and isinstance(f[path], h5py.Group):
        return path
    return None


def is_alba(f: h5py.File) -> bool:
    """An entry whose ``NXdata`` group names its signal and axes, beside an
    ``instrument/analyser`` -- and none of Diamond's ``entry1/analyser``
    ``NXdata``, which is how the two NeXus dialects tell apart."""
    entry = _entry(f)
    if entry is None or f"{entry}/instrument/analyser" not in f:
        return False
    analyser = f.get(f"{entry}/analyser")
    if isinstance(analyser, h5py.Group):
        return False
    data = _nxdata(f, entry)
    if data is None:
        return False
    group = f[data]
    signal = text(attr(group, "signal", ""))
    return bool(signal) and signal in group \
        and attr(group, "axes") is not None


def _axis_names(group: h5py.Group, ndim: int) -> list:
    """The vector named for each dimension: ``@axes``, else each vector's
    own 1-based ``@axis``."""
    axes = attr(group, "axes")
    if axes is not None:
        names = [text(a) for a in np.atleast_1d(_as_list(axes))]
        if len(names) == ndim:
            return names
    names = [None] * ndim
    for key in group.keys():
        number_ = attr(group[key], "axis")
        try:
            index = int(text(number_)) - 1
        except (TypeError, ValueError):
            continue
        if 0 <= index < ndim and names[index] is None:
            names[index] = key
    return names


def _as_list(value):
    if isinstance(value, (bytes, str, np.bytes_)):
        return [v.strip() for v in text(value).split(",")]
    return list(np.atleast_1d(value))


def read_layout(f: h5py.File) -> Layout:
    if not is_alba(f):
        raise ValueError("not an ALBA NeXus file")
    entry = _entry(f)
    group = f[_nxdata(f, entry)]
    dataset = group[text(attr(group, "signal"))]
    shape = tuple(int(n) for n in dataset.shape)
    names = _axis_names(group, len(shape))

    def vector(i):
        name = names[i]
        if name and name in group and isinstance(group[name], h5py.Dataset):
            values = np.asarray(group[name][()], dtype=float).reshape(-1)
            if values.size == shape[i]:
                return values, text(attr(group[name], "units", ""))
        return np.arange(shape[i], dtype=float), ""

    angle_dim = next((i for i, n in enumerate(names) if n == "angles"), None)
    energy_dim = next((i for i, n in enumerate(names)
                       if n in _ENERGY_NAMES and i != angle_dim), None)
    if angle_dim is None or energy_dim is None:
        if len(shape) < 2:
            raise ValueError(f"spectra of shape {shape} have no detector plane")
        angle_dim = len(shape) - 2 if angle_dim is None else angle_dim
        energy_dim = len(shape) - 1 if energy_dim is None else energy_dim
    angle, _ = vector(angle_dim)
    energy, _ = vector(energy_dim)

    dims = []
    scan_axes = [i for i in range(len(shape)) if i not in (angle_dim,
                                                            energy_dim)]
    for i in scan_axes:
        values, unit = vector(i)
        name = names[i] or f"axis {i + 1}"
        role, label = motor_role(name)
        if names[i] is None:
            role, label = "index", "Scan point (index)"
        dims.append(Dim(name, unit, values, role, label or ""))

    info = read_info(f, entry, _INFO_GROUPS)
    hv = number(f, f"{entry}/instrument/monochromator/energy")
    if hv is not None:
        info["photon_energy_eV"] = hv
    temperature = number(f, f"{entry}/sample/temperature")
    if temperature is not None:
        info["SampleTemperature_K"] = temperature

    energies = None
    energy_scale = "binding" if "binding" in (names[energy_dim] or "") \
        else "kinetic"
    if any(d.role == "photon_energy" for d in dims):
        path = f"{entry}/instrument/analyser/energies"
        if path in f and isinstance(f[path], h5py.Dataset):
            values = np.asarray(f[path][()], dtype=float)
            n_hv = next(len(d.axis) for d in dims
                        if d.role == "photon_energy")
            if values.shape == (energy.size, n_hv):
                values = values.T
            if values.shape == (n_hv, energy.size):
                energies = values
                if len(dims) > 1:
                    energies = None         # cannot place them on the grid

    return Layout(dataset=dataset, perm=tuple(scan_axes) + (angle_dim,
                                                            energy_dim),
                  dims=dims, angle=angle, energy=energy,
                  facility="ALBA LOREA (BL20)", prefix="ALBA_",
                  energy_scale=energy_scale, energies=energies, info=info)


class AlbaLoader(NexusArpesLoader, Loader):
    name = "ALBA LOREA"
    description = ("NeXus .nxs files from the LOREA ARPES beamline (BL20) "
                   "at ALBA: cuts, deflector maps, photon-energy scans and "
                   "spatial scans.")
    patterns = ("*.nxs",)
    priority = 12

    def recognises(self, f: h5py.File) -> bool:
        return is_alba(f)

    def read_layout(self, f: h5py.File) -> Layout:
        return read_layout(f)


register(AlbaLoader())
