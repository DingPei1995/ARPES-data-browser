"""
loader/maestro_loader.py
========================
ALS beamline 7.0.2 (MAESTRO): the HDF5 files its acquisition software
writes.

Ported from the MAESTRO reader in TensorSpec
(<https://github.com/fawkesdx/TensorSpec>,
``tensorspec/core/io/loaders/maestro``, MIT licence), which worked out the
layout against real files. What is kept from it, in one place:

* **The signature.** A MAESTRO file has three top-level groups,
  ``0D_Data``, ``2D_Data`` and ``Headers``; nothing else here writes that.
* **The scan plan.** ``Headers/Low_Level_Scan`` is a table of
  ``(longname, name, value, comment)`` rows. ``lwlvnm`` is the scan's name,
  ``lwlvlpn`` the number of loops, ``nmsbdv<i>`` the motors in loop ``i``,
  and ``nm_/un_/st_/en_/n_<i>_<j>`` each motor's name, unit, start, end and
  number of steps. The motors are stepped on an even grid, so the axes come
  from the plan, not from the motor readbacks.
* **Where the spectra are.** ``2D_Data/Process_000`` when present, else the
  first ``2D_Data/Fixed_Spectra*``. A fixed-mode (``Headers/DAQ_Fixed``)
  buffer is ``(angle, energy, point)`` or ``(point, angle, energy)``: the
  detector plane is always angle then energy, whatever order ``unitNames``
  lists the two scales in. The point axis is the one whose length matches
  the plan (or ``Headers/Main`` ``num_cycles``).
* **The order of the points.** Which motor steps fastest depends on the
  kind of scan, and is taken from TensorSpec kind by kind: an XY mesh has X
  fastest; a focus scan (an XY loop around a deflector loop) has the
  deflector fastest, then X, then Y; a two-motor line (deflector and X in
  one loop) has the motor listed first stepping *slowest*.
* **A detector axis in pixels** is turned into degrees at 0.048 deg per
  pixel, the Scienta R4000 Angular30 figure.
* **An aborted scan** keeps what was recorded: an XY mesh is cut to its
  complete rows; anything else becomes a list of scan points.

The kinds, and what each becomes here:

=====================  ====================================  ==============
MAESTRO scan           axes                                  opens as
=====================  ====================================  ==============
one spectrum           angle, energy                         ``cut``
deflector scan         deflector, angle, energy              ``map``
XY mesh                Y, X, angle, energy                   ``spem_4d``
two-motor line         X, deflector, angle, energy           one ``map``
                                                             per X position
focus XY               Y, X, deflector, angle, energy        one
                                                             ``spem_4d`` per
                                                             deflector step
swept ``Process_000``  0D motors (if a grid), angle, energy  as above
=====================  ====================================  ==============

The five-axis and four-axis-with-deflector scans have no single kind in this
program, so the file is listed as several datasets -- one per step of the
extra axis -- each of which is an ordinary map or spatial scan.

Nothing is read until a dataset is opened, and then only what is looked at:
the spectra stay in the file behind :class:`MaestroCube`, the way a SOLEIL
spatial scan does behind :class:`loader.nxs_file.LazyCube`.

Nothing in here imports Qt.
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field
from itertools import product

import h5py
import numpy as np

from loader.nxs_file import HANDLES, NxsScan
from loader.registry import Loader, register

__all__ = ["MaestroLoader", "MaestroCube", "ScanMotor", "ScanPlan",
           "is_maestro", "parse_low_level_scan", "read_layout",
           "PIXEL_ANGLE_STEP"]

#: Degrees per detector pixel for the Scienta R4000 in Angular30, used when
#: a fixed-mode buffer gives its slit axis in pixels.
PIXEL_ANGLE_STEP = 0.048

_REQUIRED_GROUPS = ("0D_Data", "2D_Data", "Headers")
_X_NAMES = {"scan x", "sample x"}
_Y_NAMES = {"scan y", "sample y"}


# --------------------------------------------------------------------------
# The scan plan
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class ScanMotor:
    name: str
    units: str
    start: float
    end: float
    n: int

    @property
    def axis(self) -> np.ndarray:
        return np.linspace(self.start, self.end, self.n)


@dataclass(frozen=True)
class ScanPlan:
    mode_name: str
    loops: tuple             # of tuples of ScanMotor; loop 0 is the inner one
    parallel: bool = False

    @property
    def motors(self) -> list:
        """Every motor, in the order the plan lists them."""
        return [motor for loop in self.loops for motor in loop]

    @property
    def expected_cycles(self) -> int:
        return int(np.prod([motor.n for motor in self.motors])) \
            if self.motors else 1

    def motors_slow_to_fast(self) -> list:
        """Loop 0 is the innermost (fastest) loop, and within a loop the
        motor listed first is the faster one."""
        return list(reversed(self.motors))

    def angle_motors(self) -> list:
        return [motor for motor in self.motors if _is_angle_motor(motor.name)]

    def xy_motors(self):
        """``(x, y)`` from the one loop that steps both, or None."""
        for loop in self.loops:
            x = next((m for m in loop if m.name.casefold() in _X_NAMES), None)
            y = next((m for m in loop if m.name.casefold() in _Y_NAMES), None)
            if x is not None and y is not None:
                return x, y
        return None


def _is_angle_motor(name: str) -> bool:
    label = name.casefold()
    if label.endswith(" x") or label.endswith(" y") \
            or " x" in label or " y" in label:
        return False
    return any(word in label
               for word in ("defl", "deflection", "theta", "tilt", "phi"))


def _is_spatial_motor(name: str) -> bool:
    label = name.casefold()
    return label in _X_NAMES or label in _Y_NAMES


def _text(value) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, np.ndarray) and value.size == 1:
        return _text(value.reshape(-1)[0])
    return str(value)


def _unquote(value) -> str:
    text = _text(value).strip()
    if len(text) >= 2 and text[0] == text[-1] and text[0] in "'\"":
        return text[1:-1]
    return text


def parse_low_level_scan(table) -> ScanPlan:
    """The scan plan from the rows of ``Headers/Low_Level_Scan``."""
    fields = {}
    for row in table:
        longname, name, value = row[0], row[1], row[2]
        tag = _text(longname) or _text(name)
        fields[tag.casefold().strip()] = value

    def get(tag, default="0"):
        return _unquote(fields.get(tag, default)) or default

    loops = []
    for i in range(int(float(get("lwlvlpn")))):
        motors = []
        for j in range(int(float(get(f"nmsbdv{i}")))):
            key = f"{i}_{j}"
            motors.append(ScanMotor(
                name=get(f"nm_{key}", ""), units=get(f"un_{key}", ""),
                start=float(get(f"st_{key}")), end=float(get(f"en_{key}")),
                n=int(float(get(f"n_{key}")))))
        loops.append(tuple(motors))
    parallel = get("scanpar", "").casefold() in {
        "t", "true", "p", "parallel", "y", "yes"}
    return ScanPlan(mode_name=get("lwlvnm", ""), loops=tuple(loops),
                    parallel=parallel)


# --------------------------------------------------------------------------
# Looking at a file
# --------------------------------------------------------------------------
def is_maestro(f: h5py.File) -> bool:
    return all(name in f for name in _REQUIRED_GROUPS)


def _is_fixed(f: h5py.File) -> bool:
    headers = f.get("Headers")
    return headers is not None and "DAQ_Fixed" in headers


def _spectra_dataset(f: h5py.File) -> h5py.Dataset:
    group = f["2D_Data"]
    if isinstance(group.get("Process_000"), h5py.Dataset):
        return group["Process_000"]
    for name in sorted(group.keys()):
        if name.startswith("Fixed_Spectra") \
                and isinstance(group[name], h5py.Dataset):
            return group[name]
    for name in group.keys():
        if isinstance(group[name], h5py.Dataset):
            return group[name]
    raise ValueError("no spectra dataset in 2D_Data")


def _scales(dataset):
    """``(units, offsets, deltas)`` from the dataset's attributes, or None."""
    attrs = dataset.attrs
    if not all(key in attrs for key in ("unitNames", "scaleOffset",
                                        "scaleDelta")):
        return None
    return ([_text(u) for u in attrs["unitNames"]],
            [float(v) for v in attrs["scaleOffset"]],
            [float(v) for v in attrs["scaleDelta"]])


def _is_energy_unit(unit: str) -> bool:
    unit = unit.casefold()
    return "ev" in unit or "energy" in unit


def _is_angle_unit(unit: str) -> bool:
    unit = unit.casefold()
    return any(word in unit for word in ("deg", "pixel", "angle"))


def _header_num_cycles(f: h5py.File):
    main = f.get("Headers/Main")
    if main is None:
        return None
    try:
        if isinstance(main, h5py.Group):
            for key in main.keys():
                if key.casefold() == "num_cycles":
                    return int(float(_unquote(main[key][()])))
            for key, value in main.attrs.items():
                if key.casefold() == "num_cycles":
                    return int(float(_unquote(value)))
            return None
        values = main[()]
        if getattr(values.dtype, "names", None):
            for name in values.dtype.names:
                if name.casefold() == "num_cycles":
                    return int(float(_unquote(values[name])))
        flat = np.asarray(values).reshape(-1)
        for i, value in enumerate(flat[:-1]):
            if _text(value).strip().casefold() == "num_cycles":
                return int(float(_unquote(flat[i + 1])))
    except (TypeError, ValueError):
        return None
    return None


def _mode_name(f: h5py.File) -> str:
    header = f["Headers"].get("DAQ_Fixed" if _is_fixed(f) else "DAQ_Swept")
    if isinstance(header, h5py.Dataset):
        try:
            for row in header[()]:
                if len(row) >= 3 and "mode" in _text(row[1]).casefold():
                    return _text(row[2]).strip()
        except TypeError:
            pass
    return ""


# --------------------------------------------------------------------------
# Layout: what the spectra are, and how the points are arranged
# --------------------------------------------------------------------------
@dataclass
class Dim:
    """One scanned dimension: a motor, or the plain point number."""
    name: str
    unit: str
    axis: np.ndarray
    role: str          # "x", "y", "angle", "index" or "other"


@dataclass
class Layout:
    """Everything needed to read any part of a MAESTRO file.

    The spectra are stored as a run of points along ``points_axis`` of
    ``dataset``. ``acq_shape`` is how that run folds into a grid (slowest
    first, C order), and ``dims`` are the same grid's dimensions in the
    order this program shows them, ``dims[i]`` being ``acq_shape``'s axis
    ``order[i]``.
    """
    dataset: object
    points_axis: object            # 0, 2, or None for a single spectrum
    plane: tuple                   # (angle axis, energy axis) of the plane
    angle: np.ndarray
    energy: np.ndarray
    acq_shape: tuple
    order: tuple
    dims: list
    kind_id: str
    mode: str = ""
    info: dict = field(default_factory=dict)
    angle_unit: str = "deg"
    energy_unit: str = "eV"


def _plane_axes(dataset, n1: int, n2: int, fixed: bool):
    """``(angle axis index, angle, energy, angle unit, energy unit,
    from_pixels)`` for a detector plane of ``(n1, n2)``.

    A fixed-mode plane is always ``(angle, energy)``; a swept one is in the
    order ``unitNames`` gives."""
    scales = _scales(dataset)
    if scales is None:
        if fixed:
            angle_index = 0
        else:
            angle_index = 1
        n_angle = (n1, n2)[angle_index]
        n_energy = (n1, n2)[1 - angle_index]
        angle = (np.arange(n_angle) - n_angle // 2) * PIXEL_ANGLE_STEP \
            if fixed else np.linspace(-1.0, 1.0, n_angle)
        return (angle_index, angle, np.linspace(0.0, 1.0, n_energy),
                "deg", "eV", fixed)

    units, offsets, deltas = scales
    e_scale = next((i for i, u in enumerate(units) if _is_energy_unit(u)), 0)
    a_scale = next((i for i, u in enumerate(units)
                    if i != e_scale and _is_angle_unit(u)), 1 - e_scale)
    if fixed:
        angle_index = 0
    else:
        # Swept: the scales are listed in the plane's own order.
        angle_index = 1 if e_scale == 0 else 0
    n_angle = (n1, n2)[angle_index]
    n_energy = (n1, n2)[1 - angle_index]
    energy = offsets[e_scale] + np.arange(n_energy) * deltas[e_scale]
    angle_unit = units[a_scale] if a_scale < len(units) else "deg"
    from_pixels = "pixel" in angle_unit.casefold()
    if from_pixels and fixed:
        angle = (np.arange(n_angle) - n_angle // 2) * PIXEL_ANGLE_STEP
    elif a_scale < len(offsets):
        angle = offsets[a_scale] + np.arange(n_angle) * deltas[a_scale]
    else:
        angle = np.arange(n_angle, dtype=float)
    if from_pixels or not angle_unit:
        angle_unit = "deg"
    return (angle_index, angle, energy, angle_unit,
            units[e_scale] or "eV", from_pixels and fixed)


def _points_axis(shape, expected: int, num_cycles):
    """Which end of a 3-D buffer counts the points, and how many there are.

    The plan's count is tried first, then ``num_cycles``; failing both, the
    longest end that is still shorter than the plan is taken as an aborted
    scan."""
    for wanted in (expected, num_cycles):
        if wanted is None:
            continue
        exact = [axis for axis in (0, 2) if shape[axis] == wanted]
        if len(exact) >= 1:
            return exact[0], int(shape[exact[0]])
    shorter = [axis for axis in (0, 2) if shape[axis] < expected]
    if shorter:
        longest = max(shape[axis] for axis in shorter)
        winners = [axis for axis in shorter if shape[axis] == longest]
        if len(winners) == 1:
            return winners[0], int(longest)
    raise ValueError(
        f"no end of the spectra buffer {tuple(shape)} matches the scan "
        f"plan's {expected} points")


def _motor_dim(motor: ScanMotor) -> Dim:
    label = motor.name.casefold()
    if label in _X_NAMES:
        role = "x"
    elif label in _Y_NAMES:
        role = "y"
    elif _is_angle_motor(motor.name):
        role = "angle"
    else:
        role = "other"
    return Dim(motor.name, motor.units, motor.axis, role)


def _match_kind(plan: ScanPlan, fixed: bool):
    """The TensorSpec kind for this plan, or None."""
    if not fixed:
        return None
    xy = plan.xy_motors()
    angles = plan.angle_motors()
    loops = plan.loops
    if len(loops) == 2 and xy is not None and len(angles) == 1:
        sizes = sorted(len(loop) for loop in loops)
        if sizes == [1, 2] and any(
                len(loop) == 1 and _is_angle_motor(loop[0].name)
                for loop in loops):
            return "focus_xy_fine_5d"
    if len(loops) == 1 and xy is None and len(loops[0]) == 2 \
            and len(angles) == 1 and any(
                m.name.casefold() in _X_NAMES for m in loops[0]):
        return "defl_x_line_4d"
    if len(loops) == 1 and xy is None and len(loops[0]) == 1 \
            and len(angles) == 1:
        return "fermi_defl_3d"
    if len(loops) == 1 and xy is not None and not angles:
        return "xy_fine_4d"
    return None


def _kind_grid(kind_id: str, plan: ScanPlan):
    """``(acquisition motors, shown motors)`` for a kind: the order the
    points were taken in, and the order this program shows the axes in."""
    if kind_id == "xy_fine_4d":
        x, y = plan.xy_motors()
        return [y, x], [y, x]
    if kind_id == "focus_xy_fine_5d":
        x, y = plan.xy_motors()
        return plan.motors_slow_to_fast(), [y, x, plan.angle_motors()[0]]
    if kind_id == "defl_x_line_4d":
        loop = plan.loops[0]
        x = next(m for m in loop if m.name.casefold() in _X_NAMES)
        return list(loop), [x, plan.angle_motors()[0]]
    if kind_id == "fermi_defl_3d":
        return list(plan.loops[0]), list(plan.loops[0])
    raise ValueError(kind_id)


def read_layout(f: h5py.File) -> Layout:
    """Work out how ``f``'s spectra are arranged, reading only headers."""
    if not is_maestro(f):
        raise ValueError("not a MAESTRO file: missing "
                         + ", ".join(n for n in _REQUIRED_GROUPS if n not in f))
    fixed = _is_fixed(f)
    dataset = _spectra_dataset(f)
    headers = f["Headers"]
    plan = None
    if "Low_Level_Scan" in headers:
        plan = parse_low_level_scan(headers["Low_Level_Scan"][()])
    kind_id = _match_kind(plan, fixed) if plan is not None else None

    if kind_id is not None:
        return _kind_layout(f, dataset, plan, kind_id)
    if dataset.name.endswith("/Process_000"):
        return _generic_layout(f, dataset, plan, fixed)
    if plan is not None:
        raise ValueError(
            f"no MAESTRO scan kind matches {plan.mode_name!r} "
            f"(fixed mode: {fixed}, {plan.expected_cycles} points)")
    raise ValueError("MAESTRO file has neither Headers/Low_Level_Scan nor "
                     "2D_Data/Process_000")


def _kind_layout(f, dataset, plan: ScanPlan, kind_id: str) -> Layout:
    shape = tuple(int(n) for n in dataset.shape)
    if len(shape) != 3:
        raise ValueError(f"expected a 3-D spectra buffer, got {shape}")
    expected = plan.expected_cycles
    num_cycles = _header_num_cycles(f)
    paxis, actual = _points_axis(shape, expected, num_cycles)
    n1, n2 = (shape[1], shape[2]) if paxis == 0 else (shape[0], shape[1])
    angle_index, angle, energy, angle_unit, energy_unit, from_pixels = \
        _plane_axes(dataset, n1, n2, fixed=True)

    acquisition, shown = _kind_grid(kind_id, plan)
    acq_shape = tuple(m.n for m in acquisition)
    order = tuple(acquisition.index(m) for m in shown)
    dims = [_motor_dim(m) for m in shown]
    info = {"maestro.kind": kind_id,
            "maestro.scan_name": plan.mode_name,
            "maestro.points_expected": expected,
            "maestro.points_recorded": actual}
    if num_cycles is not None:
        info["maestro.num_cycles"] = num_cycles

    if actual < expected:
        partial = None
        if kind_id == "xy_fine_4d" and actual >= acq_shape[1]:
            rows = actual // acq_shape[1]
            partial = (rows, acq_shape[1])
        if partial is not None:
            acq_shape = partial
            dims[0] = Dim(dims[0].name, dims[0].unit,
                          dims[0].axis[:partial[0]], dims[0].role)
            info["maestro.aborted"] = (
                f"scan stopped after {actual} of {expected} points; kept "
                f"the {partial[0]} complete rows")
        else:
            acq_shape, order = (actual,), (0,)
            dims = [Dim("Scan point", "", np.arange(actual, dtype=float),
                        "index")]
            info["maestro.aborted"] = (
                f"scan stopped after {actual} of {expected} points; opened "
                f"as a list of the points recorded")
    if from_pixels:
        info["maestro.angle_from_pixels"] = \
            f"{PIXEL_ANGLE_STEP} deg per pixel (R4000 Angular30)"
    return Layout(dataset=dataset, points_axis=paxis,
                  plane=(angle_index, 1 - angle_index), angle=angle,
                  energy=energy, acq_shape=acq_shape, order=order, dims=dims,
                  kind_id=kind_id, mode=plan.mode_name, info=info,
                  angle_unit=angle_unit, energy_unit=energy_unit)


def _generic_layout(f, dataset, plan, fixed: bool) -> Layout:
    """``Process_000``: points along the leading axis, if any, and the motor
    positions in ``0D_Data``, one value per point."""
    shape = tuple(int(n) for n in dataset.shape)
    if len(shape) < 2:
        raise ValueError("Process_000 has fewer than two dimensions")
    n1, n2 = shape[-2], shape[-1]
    if len(shape) == 2:
        source, paxis, actual = dataset, None, 1
    else:
        actual = int(np.prod(shape[:-2]))
        # A buffer with more than one leading axis cannot be flattened in
        # the file; such files are small enough to read whole.
        source = dataset if len(shape) == 3 else \
            np.asarray(dataset[()]).reshape(actual, n1, n2)
        paxis = 0
    # Unlike Fixed_Spectra, Process_000 lists its scales in the plane's own
    # order, fixed mode or not.
    angle_index, angle, energy, angle_unit, energy_unit, _ = \
        _plane_axes(dataset, n1, n2, fixed=False)

    mode = _mode_name(f) or (plan.mode_name if plan is not None else "")
    info = {"maestro.kind": "process000_generic",
            "maestro.points_recorded": actual}
    if mode:
        info["maestro.scan_name"] = mode

    dims, acq_shape, order = [], (), ()
    if paxis is not None:
        motors = _zero_d(f, plan, actual, info)
        grid = _motor_grid(motors, actual)
        if grid is not None:
            dims, acq_shape = grid
            order = tuple(range(len(dims)))
        else:
            dims = [Dim("Scan point", "", np.arange(actual, dtype=float),
                        "index")]
            acq_shape, order = (actual,), (0,)
    return Layout(dataset=source, points_axis=paxis,
                  plane=(angle_index, 1 - angle_index), angle=angle,
                  energy=energy, acq_shape=acq_shape, order=order, dims=dims,
                  kind_id="process000_generic", mode=mode, info=info,
                  angle_unit=angle_unit, energy_unit=energy_unit)


def _zero_d(f, plan, actual: int, info: dict) -> dict:
    """The ``0D_Data`` channels that may be motors: one value per point.

    ``0D_Data`` also records readbacks that drift (beam current, a
    temperature), which would make any grid look irregular; when the scan
    plan names its motors, only those are used. A channel that does not
    change is a setting, not an axis, and goes into ``info``."""
    names = None
    if plan is not None and plan.motors:
        wanted = {m.name.casefold() for m in plan.motors}
        names = [n for n in f["0D_Data"].keys() if n.casefold() in wanted]
    if not names:
        names = list(f["0D_Data"].keys())
    motors = {}
    for name in names:
        item = f["0D_Data"][name]
        if not isinstance(item, h5py.Dataset):
            continue
        values = np.asarray(item[()], dtype=float).reshape(-1)
        if len(values) < actual:
            continue
        values = values[:actual]
        if np.ptp(values) <= 1e-9 * max(1.0, float(np.abs(values).max())):
            info[f"maestro.{name}"] = float(values[0])
            continue
        motors[name] = values
    return motors


def _motor_grid(motors: dict, actual: int):
    """The regular grid ``motors`` step through, slowest first, or None.

    A motor's place in the grid is how often it changes: the fastest changes
    at almost every point, the slowest only once per sweep of all the
    others. The grid is accepted only if those motors, laid out that way,
    reproduce every recorded position."""
    if not motors:
        return None
    changes = {name: int(np.count_nonzero(np.abs(np.diff(v)) > _tol(v)))
               for name, v in motors.items()}
    names = sorted(motors, key=lambda name: changes[name])     # slow first
    stride = actual
    axes, shape = [], []
    for name in names:
        values = motors[name]
        # Points between changes of this motor = product of the faster ones.
        runs = np.flatnonzero(np.abs(np.diff(values)) > _tol(values))
        run = int(runs[0]) + 1 if len(runs) else actual
        count = stride // run if run else 0
        if count < 2 or stride % run:
            return None
        axis = values[:stride:run][:count]
        steps = np.diff(axis)
        if not (np.all(steps > 0) or np.all(steps < 0)):
            return None                 # an axis has to run one way
        axes.append(axis)
        shape.append(count)
        stride = run
    if stride != 1 or int(np.prod(shape)) != actual:
        return None
    mesh = np.meshgrid(*axes, indexing="ij")
    for name, grid in zip(names, mesh):
        if not np.allclose(grid.reshape(-1), motors[name],
                           atol=10 * _tol(motors[name])):
            return None
    dims = []
    for name, axis in zip(names, axes):
        dim = _motor_dim(ScanMotor(name, "", float(axis[0]), float(axis[-1]),
                                   len(axis)))
        dim.axis = np.asarray(axis, dtype=float)
        dims.append(dim)
    return dims, tuple(shape)


def _tol(values: np.ndarray) -> float:
    span = float(np.ptp(values))
    return max(span * 1e-4, 1e-9)


# --------------------------------------------------------------------------
# Reading the spectra lazily
# --------------------------------------------------------------------------
class MaestroCube:
    """Part of a MAESTRO spectra buffer, as an array that stays on disk.

    Presented as ``(*scan dims, angle, energy)``. ``fixed`` pins some of the
    layout's scan dimensions to one index -- that is how one deflector step
    of a focus scan becomes a 4-D spatial scan of its own. Indexing works as
    on the numpy array and reads only the points it needs, so a cursor move
    reads one spectrum; ``np.asarray`` reads the lot.
    """

    def __init__(self, layout: Layout, fixed: dict = None):
        self._layout = layout
        self._fixed = dict(fixed or {})
        self._free = [i for i in range(len(layout.dims)) if i not in self._fixed]
        dataset = layout.dataset
        self.dtype = np.dtype(dataset.dtype)
        if not np.issubdtype(self.dtype, np.floating):
            self.dtype = np.dtype(np.float32)
        scan = tuple(layout.acq_shape[layout.order[i]] for i in self._free)
        self.shape = scan + (len(layout.angle), len(layout.energy))
        self.ndim = len(self.shape)
        self.size = int(np.prod(self.shape))
        self.nbytes = self.size * self.dtype.itemsize

    def __len__(self):
        return self.shape[0]

    def __repr__(self):
        return f"MaestroCube(shape={self.shape}, dtype={self.dtype})"

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
        layout = self._layout
        n_scan = len(self._free)
        scan_key, (a_key, e_key) = key[:n_scan], key[n_scan:]

        # Which index of each shown dimension, then the points in the file.
        picks = []
        for dim in range(len(layout.dims)):
            if dim in self._fixed:
                picks.append(np.array([self._fixed[dim]]))
            else:
                sel = scan_key[self._free.index(dim)]
                size = layout.acq_shape[layout.order[dim]]
                picks.append(np.arange(size)[sel].reshape(-1))
        mesh = np.meshgrid(*picks, indexing="ij") if picks else []
        acq = [None] * len(layout.acq_shape)
        for dim, grid in enumerate(mesh):
            acq[layout.order[dim]] = grid
        if acq:
            flat = np.ravel_multi_index(tuple(acq), layout.acq_shape) \
                .reshape(-1)
        else:
            flat = np.zeros(1, dtype=int)
        frames = self._read(flat, a_key, e_key)

        block = frames.reshape(tuple(len(p) for p in picks)
                               + frames.shape[1:])
        # Drop what was pinned, and what an integer index asked to drop.
        drop = [i for i in range(len(layout.dims)) if i in self._fixed
                or _is_int(scan_key[self._free.index(i)])]
        index = tuple(0 if i in drop else slice(None)
                      for i in range(len(layout.dims)))
        block = block[index + (slice(None),) * (block.ndim - len(index))]
        tail = tuple(0 if _is_int(k) else slice(None) for k in (a_key, e_key))
        result = block[(Ellipsis,) + tail]
        return result.astype(self.dtype, copy=False)

    def _read(self, flat: np.ndarray, a_key, e_key) -> np.ndarray:
        """``(len(flat), angle, energy)`` -- the selected points, each
        cropped to the detector selection while it is read."""
        layout = self._layout
        dataset = layout.dataset
        a_sel = _as_slice(a_key)
        e_sel = _as_slice(e_key)
        plane_sel = [None, None]
        plane_sel[layout.plane[0]] = a_sel
        plane_sel[layout.plane[1]] = e_sel

        if layout.points_axis is None:
            frame = np.asarray(dataset[tuple(plane_sel)])
            frames = frame[None]
        else:
            order = np.argsort(flat, kind="stable")
            ordered = flat[order]
            unique, inverse = np.unique(ordered, return_inverse=True)
            sel = _point_selection(unique)
            if layout.points_axis == 0:
                raw = np.asarray(dataset[(sel,) + tuple(plane_sel)])
            else:
                raw = np.asarray(dataset[tuple(plane_sel) + (sel,)])
                raw = np.moveaxis(raw, -1, 0)
            if isinstance(sel, slice) and len(raw) != len(unique):
                raw = raw[unique - unique[0]]
            frames = np.empty((len(flat),) + raw.shape[1:], dtype=raw.dtype)
            frames[order] = raw[inverse]
        if layout.plane[0] == 1:
            frames = np.swapaxes(frames, 1, 2)
        # What h5py could not be asked for (a list, a reversed slice) is
        # applied now, to what was read.
        if a_sel is not a_key and not _is_int(a_key):
            frames = frames[:, a_key if isinstance(a_key, slice)
                            else np.asarray(a_key)]
        if e_sel is not e_key and not _is_int(e_key):
            frames = frames[:, :, e_key if isinstance(e_key, slice)
                            else np.asarray(e_key)]
        return frames

    def materialise(self) -> np.ndarray:
        return self[(slice(None),) * self.ndim]

    def __array__(self, dtype=None, copy=None):
        data = self.materialise()
        return data if dtype is None else data.astype(dtype)


def _is_int(key) -> bool:
    return isinstance(key, (int, np.integer))


def _as_slice(key):
    """What to ask h5py for along one detector axis: the key itself if it is
    an int or a slice with a positive step, everything otherwise (the rest
    is done in memory). An int is kept as a length-1 slice so that the
    plane keeps both its axes until the end."""
    if _is_int(key):
        key = int(key)
        return slice(key, key + 1) if key != -1 else slice(-1, None)
    if isinstance(key, slice) and (key.step is None or key.step > 0):
        return key
    return slice(None)


def _point_selection(points: np.ndarray):
    """A slice where the sorted, unique ``points`` are evenly spaced (the
    usual case: a row of a map, every n-th deflector step), else the list."""
    if len(points) == 1:
        return slice(int(points[0]), int(points[0]) + 1)
    steps = np.diff(points)
    if np.all(steps == steps[0]):
        return slice(int(points[0]), int(points[-1]) + 1, int(steps[0]))
    if points[-1] - points[0] + 1 <= 4 * len(points):
        return slice(int(points[0]), int(points[-1]) + 1)
    return points.tolist()


# --------------------------------------------------------------------------
# From a layout to the datasets the browser lists
# --------------------------------------------------------------------------
def _plan_entries(layout: Layout):
    """``(kind, shown dims, pinned dims)`` describing each dataset in the
    file, plus the pinned combinations: one entry when the scan fits one of
    this program's kinds, one per step of the extra axis when it does not.
    """
    dims = layout.dims
    roles = [d.role for d in dims]
    if not dims:
        return "cut", [], []
    if "x" in roles and "y" in roles:
        shown = [roles.index("y"), roles.index("x")]
        kind = "spem_4d"
    else:
        if "angle" in roles:
            first = roles.index("angle")
        elif "index" in roles:
            first = roles.index("index")
        else:
            first = len(dims) - 1          # the fastest motor
        shown = [first]
        kind = "map"
    pinned = [i for i in range(len(dims)) if i not in shown]
    return kind, shown, pinned


def _combinations(layout: Layout, pinned: list):
    sizes = [len(layout.dims[i].axis) for i in pinned]
    return list(product(*[range(n) for n in sizes]))


def _entry_key(indices) -> str:
    return "point " + ",".join(str(i) for i in indices)


def _entry_label(layout: Layout, pinned: list, indices) -> str:
    parts = []
    for dim_index, i in zip(pinned, indices):
        dim = layout.dims[dim_index]
        unit = f" {dim.unit}" if dim.unit else ""
        parts.append(f"{dim.name} = {dim.axis[i]:.4g}{unit}")
    return ", ".join(parts)


def _axis_label(dim: Dim) -> str:
    unit = _unit_text(dim.unit)
    if dim.role == "index":
        return "Scan point (index)"
    return f"{dim.name} ({unit})" if unit else dim.name


def _unit_text(unit: str) -> str:
    folded = unit.casefold()
    if folded in {"um", "micron", "microns"}:
        return "µm"
    if folded in {"deg", "degree", "degrees"}:
        return "deg"
    return unit


_KIND_LABELS = {"cut": "Cut", "map": "Map", "spem_4d": "SPEM"}


def _role_for_axis0(dim: Dim) -> str:
    if dim.role == "angle":
        return "angle"
    if dim.role in ("x", "y"):
        return "position"
    return "other"


def _build_scan(layout: Layout, kind: str, shown: list, pinned: list,
                indices) -> NxsScan:
    fixed = dict(zip(pinned, indices))
    # MaestroCube shows the free dims in the layout's order; a 4-D spatial
    # scan wants (y, x), which _plan_entries already put first in `shown`.
    free = [i for i in range(len(layout.dims)) if i not in fixed]
    cube = MaestroCube(layout, fixed)
    if free != shown:
        cube = _Transposed(cube, [free.index(i) for i in shown] + [
            len(free), len(free) + 1])

    angle_label = f"Angle along slit ({_unit_text(layout.angle_unit)})"
    energy_label = f"Energy ({layout.energy_unit})"
    scan = NxsScan(kind=kind, filename_prefix="MAESTRO_")
    if kind == "cut":
        scan.x, scan.y = layout.angle, layout.energy
        scan.value = cube
        scan.labels = {"x": angle_label, "y": energy_label}
    elif kind == "map":
        dim = layout.dims[shown[0]]
        scan.x, scan.k, scan.z = dim.axis, layout.angle, layout.energy
        scan.value = cube
        scan.labels = {"x": _axis_label(dim), "k": angle_label,
                       "z": energy_label}
        scan.info["axis0.role"] = _role_for_axis0(dim)
    elif kind == "spem_4d":
        y_dim, x_dim = layout.dims[shown[0]], layout.dims[shown[1]]
        scan.y, scan.x = y_dim.axis, x_dim.axis
        scan.k, scan.z = layout.angle, layout.energy
        scan.value4d = cube
        scan.labels = {"x": _axis_label(x_dim), "y": _axis_label(y_dim),
                       "k": angle_label, "z": energy_label}
    else:                                                  # pragma: no cover
        raise AssertionError(kind)

    scan.info.update(layout.info)
    scan.info["facility"] = "ALS MAESTRO (BL 7.0.2)"
    if layout.mode:
        scan.info["maestro.mode"] = layout.mode
    scan.info["maestro.dataset"] = layout.dataset.name \
        if hasattr(layout.dataset, "name") else "2D_Data/Process_000"
    for dim_index, i in fixed.items():
        dim = layout.dims[dim_index]
        scan.fourd_info[dim.name] = f"{dim.axis[i]:.6g} {dim.unit}".strip()
    if fixed:
        scan.info["maestro.slice"] = _entry_label(layout, pinned, indices)
    for dim in layout.dims:
        if dim.role != "index":
            scan.fourd_info.setdefault(
                f"{dim.name} range",
                f"{dim.axis[0]:.6g} .. {dim.axis[-1]:.6g} {dim.unit} "
                f"({len(dim.axis)} steps)".strip())
    return scan


class _Transposed:
    """A :class:`MaestroCube` with its axes in another order."""

    def __init__(self, cube: MaestroCube, perm: list):
        self._cube = cube
        self._perm = list(perm)
        self.shape = tuple(cube.shape[p] for p in perm)
        self.dtype = cube.dtype
        self.ndim = cube.ndim
        self.size = cube.size
        self.nbytes = cube.nbytes

    def __len__(self):
        return self.shape[0]

    def __getitem__(self, key):
        key = self._cube._normalise(key)
        source = [None] * self.ndim
        for out_axis, sel in enumerate(key):
            source[self._perm[out_axis]] = sel
        data = self._cube[tuple(source)]
        kept_source = [a for a in range(self.ndim) if not _is_int(source[a])]
        wanted = [self._perm[a] for a, sel in enumerate(key) if not _is_int(sel)]
        if len(wanted) < 2:
            return data
        return np.transpose(data, [kept_source.index(a) for a in wanted])

    def materialise(self):
        return self[(slice(None),) * self.ndim]

    def __array__(self, dtype=None, copy=None):
        data = self.materialise()
        return data if dtype is None else data.astype(dtype)


# --------------------------------------------------------------------------
# The loader itself
# --------------------------------------------------------------------------
class MaestroLoader(Loader):
    name = "ALS MAESTRO"
    description = ("HDF5 files from beamline 7.0.2 (MAESTRO) at the ALS: "
                   "single spectra, deflector maps, XY spatial scans, "
                   "two-motor lines and focus scans.")
    patterns = ("*.h5", "*.hdf5")
    priority = 15

    def can_open(self, path: str) -> bool:
        """Recognised by the three top-level groups every MAESTRO file has
        (``0D_Data``, ``2D_Data``, ``Headers``), not by its name."""
        try:
            with HANDLES.borrow(path) as f:
                return is_maestro(f)
        except Exception:
            return False

    def list_entries(self, path: str) -> list:
        """One row for the scan, or one per step of the axis that does not
        fit a kind here (see the module docstring)."""
        base = os.path.splitext(os.path.basename(path))[0]
        try:
            with HANDLES.borrow(path) as f:
                layout = read_layout(f)
        except Exception:
            return []
        kind, _shown, pinned = _plan_entries(layout)
        title = layout.mode or layout.kind_id
        combos = _combinations(layout, pinned)
        if not pinned:
            return [{"entry": "scan", "kind": _KIND_LABELS[kind],
                     "name": base, "title": title, "start_time": None}]
        return [{"entry": _entry_key(indices), "kind": _KIND_LABELS[kind],
                 "name": f"{base} {_entry_label(layout, pinned, indices)}",
                 "title": title, "start_time": None}
                for indices in combos]

    def load(self, path: str, entry: str = None, progress=None):
        f = HANDLES.acquire(path)
        keep_open = False
        try:
            layout = read_layout(f)
            kind, shown, pinned = _plan_entries(layout)
            if pinned:
                if entry and entry.startswith("point "):
                    indices = tuple(int(i) for i in entry[6:].split(","))
                else:
                    indices = (0,) * len(pinned)
                if len(indices) != len(pinned) or any(
                        not 0 <= i < len(layout.dims[d].axis)
                        for d, i in zip(pinned, indices)):
                    raise KeyError(f"no dataset {entry!r} in {path}")
            else:
                indices = ()
            scan = _build_scan(layout, kind, shown, pinned, indices)
            if isinstance(layout.dataset, h5py.Dataset):
                scan._h5file = f
                scan._h5path = path
                keep_open = True
            return scan
        finally:
            if not keep_open:
                HANDLES.release(path, f)


register(MaestroLoader())
