"""
A dataset between steps, with no Qt.

The viewer's in-memory dataset (``ui.widgets.MemoryData``) lives in the Qt
package, because it also serves the windows; its bookkeeping is simple and
is reproduced here exactly, so what a batch run writes is indistinguishable
from what the dialogs write:

* ``info["_kind"]`` / ``info["_source"]`` -- what it is and what it is called;
* ``info["<prefix>.<parameter>"]`` -- the settings of the step that made it,
  under the prefix the dialog uses (``degrid.``, ``fscorr.``, ``fitEF.``,
  ``kconv.``) -- which is also what ``tools.degrid.not_pixel_locked`` reads
  to refuse a second de-grid after a correction;
* ``info["proc.step.N"]`` -- the history, via ``tools.process.record_step``.
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field

import numpy as np

from loader import registry
from loader.nxs_file import axis_slots, save_dataset
from tools.process import Step, record_step


@dataclass
class Dataset:
    kind: str
    axes: tuple                     # constructor order (AXIS_SLOTS)
    value: object                   # numpy array, or the loader's lazy cube
    labels: dict
    info: dict = field(default_factory=dict)
    motors: dict = field(default_factory=dict)
    name: str = "dataset"
    source_path: str = ""
    _scan: object = None            # the NxsScan holding the file open, if any

    # -- in -----------------------------------------------------------------
    @classmethod
    def load(cls, path: str, entry: str = None, name: str = None) -> "Dataset":
        scan = registry.load(path, entry)
        slots = axis_slots(scan.kind, "constructor")
        if not slots:
            scan.close()
            raise ValueError(f"{os.path.basename(path)} [{entry}]: a {scan.kind!r} "
                             f"has no axes-and-array form to process")
        value = scan.value4d if scan.kind == "spem_4d" else scan.value
        label = name or f"{os.path.splitext(os.path.basename(path))[0]}" + (
            f" {entry}" if entry else "")
        return cls(scan.kind, tuple(getattr(scan, s) for s in slots), value,
                   dict(scan.labels or {}), dict(scan.info or {}),
                   dict(getattr(scan, "fourd_info", {}) or {}), label,
                   os.path.abspath(path), scan)

    def close(self):
        if self._scan is not None:
            self._scan.close()
            self._scan = None

    # -- looking at it -------------------------------------------------------
    @property
    def slots(self) -> tuple:
        return axis_slots(self.kind, "constructor")

    def axis(self, slot: str) -> np.ndarray:
        return np.asarray(self.axes[self.slots.index(slot)], dtype=float)

    def label(self, slot: str) -> str:
        return str(self.labels.get(slot, slot))

    def array(self, dtype=None) -> np.ndarray:
        """Every point, read now if it was still on disk."""
        value = self.value
        if hasattr(value, "materialise"):
            value = value.materialise()
        return np.asarray(value, dtype=dtype) if dtype else np.asarray(value)

    @property
    def scan(self):
        """What ``tools.process.history_of`` and friends read (``.info``)."""
        return self

    @property
    def shape(self) -> tuple:
        return tuple(getattr(self.value, "shape", np.shape(self.value)))

    # -- out -------------------------------------------------------------------
    def derive(self, *, step: str, params: dict = None, prefix: str = None,
               value=None, axes=None, labels=None, kind=None,
               suffix: str = None, extra_info: dict = None) -> "Dataset":
        """The next dataset in the chain, with the step recorded the way the
        viewer records it. Anything not given is carried over."""
        params = dict(params or {})
        kind = kind or self.kind
        name = f"{self.name} [{suffix or step}]"
        info = record_step(dict(self.info), Step(step, params, source=self.name))
        info["_kind"] = kind
        info["_source"] = name
        for key, val in params.items():
            info[f"{prefix or step}.{key}"] = _flat(val)
        info.update(extra_info or {})
        return Dataset(kind, tuple(axes if axes is not None else self.axes),
                       self.value if value is None else value,
                       dict(labels if labels is not None else self.labels),
                       info, dict(self.motors), name, self.source_path, None)

    def to_save_dict(self) -> dict:
        return {"name": self.name, "kind": self.kind,
                "axes": tuple(np.asarray(a, dtype=float) for a in self.axes),
                "labels": dict(self.labels), "value": self.array(),
                "info": {k: _flat(v) for k, v in self.info.items()},
                "motors": dict(self.motors)}

    def save(self, path: str) -> str:
        """Write in the viewer's own format, atomically: a crash half-way
        leaves the previous file (or nothing), never a truncated one."""
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
        tmp = path + ".partial"
        save_dataset(tmp, [self.to_save_dict()])
        os.replace(tmp, path)
        return path


def _flat(value):
    """``save_dataset`` writes scalars and strings only; the dialogs flatten
    lists the same way (see ``kconv.*`` in ui/windows.py)."""
    if isinstance(value, (list, tuple, np.ndarray)):
        return ", ".join(f"{v:.6g}" if isinstance(v, (float, np.floating)) else str(v)
                         for v in np.asarray(value, dtype=object).ravel())
    if isinstance(value, (bool, np.bool_)):
        return int(value)
    if value is None:
        return "None"
    return value
