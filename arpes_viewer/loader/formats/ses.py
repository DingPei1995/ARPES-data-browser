"""
loader/formats/ses.py
=====================
Scienta SES files, as every SES beamline writes them: Igor binary waves
(``.ibw``) for cuts, photon-energy scans and manipulator scans, the SES
``.zip`` for a deflector map, and the multi-region ``.txt`` export.

Written first for MAX IV Bloch against the headers of a cut, a photon-energy
scan, a spatial scan and the ``.ini`` files of a deflector-map zip, and
checked against the SES reader in peaks
(<https://github.com/phrgab/peaks>, ``peaks/core/fileIO/loaders``); then
extended with the headers from BESSY One Cube / One Square and SOLARIS
URANOS (``for loader/BESSY``, ``for loader/SOLARIS``). Nothing here is
specific to a beamline: each one is a :class:`SESLoader` subclass in
``loader/beamlines`` that says which files are its own (:meth:`claims`).

The ``.ibw``
------------
One Igor wave, format version 5. Igor stores a wave column-major, so the
array read in C order is ``(chunks, layers, columns, rows)`` -- which for
SES is already ``(*scanned axes, slit angle, energy)``:

=========  ===========================================  ==================
dimension  SES writes                                   becomes
=========  ===========================================  ==================
rows       ``Kinetic Energy [eV]`` or ``Binding         energy
           Energy [eV]``
columns    ``Y-Scale [deg]``, the angle along the slit  angle
layers     ``Photon Energy [eV]`` (CIS), ``Point``      the scanned axis
           (manipulator scan), a deflector angle
chunks     ``Position`` (the slow manipulator axis)     the second one
=========  ===========================================  ==================

Each axis is ``start + i * delta`` from the wave's own scaling. The wave
note is the SES header: ``[SES]`` (region, lens, pass energy, hv, ...),
``[Run Mode Information]`` and ``[Manipulator]`` (``A P T X Y Z``).

A manipulator scan's ``Point`` and ``Position`` are named in the run-mode
table, which lists every point's ``Position, Point`` and the motors
(``A P T X Y Z`` at Bloch, ``P-Axis R-Axis Z-Axis`` at One Square,
``X: Y: Z: R1: R2: R3:`` at URANOS). The motor a dimension steps is found
in that table: the one whose column repeats ``Point`` (Bloch), the one that
repeats the wave's own axis (One Square), or the one that moves along that
dimension and nowhere else -- whose readback then replaces an axis SES only
numbered ``1, 2, 3 ...`` (URANOS). Two manipulator axes make a spatial scan
(``spem_4d``) even when one of them is an angle; one makes a map.

A binding-energy axis (SES's ``Energy Scale=Binding``, positive below E_F)
is turned into ``E - E_F`` by negating it, without a work function: SES has
already referred it. A CIS scan recorded on a kinetic scale keeps one
binding window, so each hv's kinetic axis is the first one shifted by the
change in hv.

The ``.zip``
------------
SES's own map format: ``viewer.ini`` lists the regions, and each region has
``Spectrum_<name>.bin`` (float32, energy fastest, then the slit angle, then
the deflector), ``Spectrum_<name>.ini`` (the three axes' size, offset,
step and label) and ``<name>.ini`` (the SES header). A zip with several
regions lists one dataset per region. A ``.bin`` stored uncompressed is
read in place; a compressed one is read into memory.

The ``.txt``
------------
SES's text export, parsed by :func:`loader.beamlines.cassiopee.parse_scienta`:
one ``[Region n]`` per spectrum. A file whose regions were taken at
different photon energies (URANOS writes a whole photon-energy scan into
one file this way) opens as one kz map, each spectrum referred to E_F by its
own ``hv``; otherwise every region is a cut.

Everything but a kz map stays on disk and is read as it is looked at.

Nothing in here imports Qt.
"""
from __future__ import annotations

import configparser
import os
import re
import struct
import zipfile
from dataclasses import dataclass

import numpy as np

from arpes_viewer.loader.formats.nexus_arpes import (Dim, Layout, build_scan,
                                        list_layout, motor_role, straighten)
from arpes_viewer.loader.nxs_file import ANGLE_UNIT, ENERGY_UNIT, KIND_LABELS, NxsScan
from arpes_viewer.loader.registry import Loader

__all__ = ["SESLoader", "ses_header", "ses_location", "IgorWave",
           "read_ibw", "parse_ses_note", "read_layout_ibw",
           "read_layouts_zip", "read_layout_text", "is_ses_ibw",
           "is_ses_zip", "is_ses_text", "text_header"]

FACILITY = "Scienta SES"
PREFIX = "SES_"

# --------------------------------------------------------------------------
# Igor binary wave, version 5
# --------------------------------------------------------------------------
_BIN_HEADER = 64
_WAVE_HEADER = 320
_IGOR_TYPES = {
    0x02: "f4", 0x04: "f8",
    0x08: "i1", 0x10: "i2", 0x20: "i4",
    0x48: "u1", 0x50: "u2", 0x60: "u4",
}


@dataclass
class IgorWave:
    """A wave's header, and where its numbers are in the file.

    ``shape`` is Igor's ``(rows, columns, layers, chunks)`` cut to the
    dimensions in use; ``start``/``delta``/``units`` go with it. The array
    itself is not read: :meth:`array` maps it.
    """
    path: str
    name: str
    shape: tuple
    start: tuple
    delta: tuple
    units: tuple
    dtype: np.dtype
    offset: int
    note: str
    data_units: str = ""

    def array(self) -> np.ndarray:
        """The numbers, memory-mapped, in C order: ``shape`` reversed."""
        return np.memmap(self.path, dtype=self.dtype, mode="r",
                         offset=self.offset, shape=tuple(reversed(self.shape)))

    def axis(self, i: int) -> np.ndarray:
        return self.start[i] + self.delta[i] * np.arange(self.shape[i])


def _cstr(raw: bytes) -> str:
    return raw.split(b"\0", 1)[0].decode("latin-1").strip()


def read_ibw(path: str) -> IgorWave:
    """The header of a version-5 ``.ibw``. Raises ``ValueError`` on
    anything else."""
    with open(path, "rb") as f:
        head = f.read(_BIN_HEADER + _WAVE_HEADER)
        if len(head) < _BIN_HEADER + _WAVE_HEADER:
            raise ValueError("too short for an Igor binary wave")
        for order in "<>":
            if struct.unpack(order + "h", head[:2])[0] == 5:
                break
        else:
            raise ValueError("not a version-5 Igor binary wave")
        (_v, _chk, wfm_size, formula_size, note_size, data_eunits_size
         ) = struct.unpack(order + "hhiiii", head[:20])
        dim_eunits = struct.unpack(order + "4i", head[20:36])
        w = head[_BIN_HEADER:]
        npnts, = struct.unpack(order + "i", w[12:16])
        type_code, = struct.unpack(order + "h", w[16:18])
        name = _cstr(w[28:60])
        n_dim = struct.unpack(order + "4i", w[68:84])
        sf_a = struct.unpack(order + "4d", w[84:116])
        sf_b = struct.unpack(order + "4d", w[116:148])
        data_units = _cstr(w[148:152])
        dim_units = [_cstr(w[152 + 4 * i:156 + 4 * i]) for i in range(4)]
        if type_code not in _IGOR_TYPES:
            raise ValueError(f"unsupported Igor number type {type_code:#x}")
        dtype = np.dtype(order + _IGOR_TYPES[type_code])
        ndim = next((i for i, n in enumerate(n_dim) if n == 0), 4)
        shape = tuple(int(n) for n in n_dim[:ndim]) or (int(npnts),)
        if int(np.prod(shape)) != npnts:
            raise ValueError("wave dimensions do not match its point count")

        f.seek(_BIN_HEADER + wfm_size + formula_size)
        note = f.read(max(note_size, 0)).decode("latin-1")
        extended = f.read(max(data_eunits_size, 0)).decode("latin-1")
        data_units = extended.strip("\0 ") or data_units
        units = []
        for i in range(4):
            text = f.read(max(dim_eunits[i], 0)).decode("latin-1").strip("\0 ")
            units.append(text or dim_units[i])
    return IgorWave(path=path, name=name, shape=shape,
                    start=tuple(sf_b[:len(shape)]),
                    delta=tuple(sf_a[:len(shape)]),
                    units=tuple(units[:len(shape)]), dtype=dtype,
                    offset=_BIN_HEADER + _WAVE_HEADER, note=note,
                    data_units=data_units)


# --------------------------------------------------------------------------
# The SES header
# --------------------------------------------------------------------------
def parse_ses_note(note: str) -> dict:
    """``{section: [lines]}`` of an SES header (the wave note or a region's
    ``.ini``). Lines keep their ``\\x0b`` separators: the run-mode table
    uses them between columns."""
    sections: dict = {}
    current = None
    for line in re.split(r"\r\n|\r|\n", note):
        stripped = line.strip()
        match = re.fullmatch(r"\[(.+)\]", stripped)
        if match:
            current = match.group(1).strip()
            sections.setdefault(current, [])
        elif current is not None and stripped:
            sections[current].append(line.strip("\r\n"))
    return sections


def _keys(lines) -> dict:
    out = {}
    for line in lines or ():
        if "=" in line and "\x0b" not in line:
            key, value = line.split("=", 1)
            out.setdefault(key.strip(), value.strip())
    return out


def _float(value):
    try:
        return float(str(value).split()[0])
    except (TypeError, ValueError, IndexError):
        return None


def _run_mode_table(lines) -> dict:
    """``{column: values}`` from a manipulator scan's run-mode table.

    SES separates the columns with ``\\x0b`` and the rows with line breaks;
    a dump that turned every separator into a line break is read too, by
    counting the header's columns."""
    tokens = []
    for line in lines or ():
        if "=" in line and "\x0b" not in line:
            continue                    # Name=..., Point 1=30 eV, ...
        tokens.extend(t.strip() for t in line.split("\x0b"))
    tokens = [t for t in tokens if t]
    header = []
    for token in tokens:
        if _float(token) is not None and not re.search(r"[A-Za-z]", token):
            break
        header.append(token)
    numbers = [_float(t) for t in tokens[len(header):]]
    if not header or not numbers or any(v is None for v in numbers):
        return {}
    rows = len(numbers) // len(header)
    if rows == 0:
        return {}
    table = np.asarray(numbers[:rows * len(header)]).reshape(rows,
                                                             len(header))
    return {name: table[:, i] for i, name in enumerate(header)}


#: Manipulator columns of the run-mode table, by their folded name:
#: ``(label, unit)``. Bloch writes ``A P T X Y Z``, APE ``polar tilt X Y Z``,
#: One Square ``P-Axis R-Axis Z-Axis``, URANOS ``X: Y: Z: R1: R2: R3:``.
_MOTORS = {
    "x": ("X", "mm"), "y": ("Y", "mm"), "z": ("Z", "mm"),
    "x-axis": ("X", "mm"), "y-axis": ("Y", "mm"), "z-axis": ("Z", "mm"),
    "p": ("Polar", "deg"), "polar": ("Polar", "deg"),
    "p-axis": ("Polar", "deg"),
    "t": ("Tilt", "deg"), "tilt": ("Tilt", "deg"),
    "t-axis": ("Tilt", "deg"),
    "a": ("Azimuth", "deg"), "azi": ("Azimuth", "deg"),
    "azimuth": ("Azimuth", "deg"), "r-axis": ("Rotation", "deg"),
    "r1": ("R1", "deg"), "r2": ("R2", "deg"), "r3": ("R3", "deg"),
}


def _fold(name: str) -> str:
    """A run-mode column's name as :data:`_MOTORS` keys it: URANOS ends
    its motor names with a colon (``X:``)."""
    return name.strip().rstrip(":").strip().casefold()


def _motor_for(table: dict, column: str, axis=None, index: int = 0,
               counts: tuple = ()):
    """``(label, unit, readback)`` of the manipulator axis that scanned
    dimension ``index`` of ``counts`` (slowest first) steps, or None.

    Tried in turn, spatial motors first each time: the motor whose column
    repeats the dimension's own column (``Point``, ``Position``); the one
    that repeats the wave's ``axis`` point by point; and the one that moves
    along this dimension and stays put along every other. ``readback`` is
    that motor's positions along the dimension when they are what places
    it -- the third case, where the wave's own axis is only a count -- and
    None otherwise."""
    motors = sorted((name for name in table if _fold(name) in _MOTORS),
                    key=lambda n: _MOTORS[_fold(n)][1] != "mm")
    values = table.get(column)
    if values is not None:
        for name in motors:
            other = table[name]
            if other.size == values.size and np.allclose(other, values,
                                                         atol=1e-6):
                return _MOTORS[_fold(name)] + (None,)
    if axis is None or not counts:
        return None
    rows = int(np.prod(counts))
    if not motors or table[motors[0]].size != rows:
        return None
    axis = np.asarray(axis, dtype=float)
    grids = np.meshgrid(*[np.arange(n) for n in counts], indexing="ij")
    expected = axis[grids[index].reshape(-1)]
    step = abs(axis[1] - axis[0]) if axis.size > 1 else 1.0
    for name in motors:
        if np.allclose(table[name], expected, atol=max(1e-6, 1e-3 * step)):
            return _MOTORS[_fold(name)] + (None,)
    for name in motors:
        grid = table[name].reshape(counts)
        along = np.moveaxis(grid, index, -1).reshape(-1, counts[index])
        if counts[index] < 2 or np.ptp(along[0]) <= 1e-9:
            continue
        if np.allclose(along, along[0], atol=1e-6 + 1e-3 * np.ptp(along[0])):
            return _MOTORS[_fold(name)] + (along[0].astype(float),)
    return None


def _split_label(label: str) -> tuple:
    """``("Photon Energy", "eV")`` from ``"Photon Energy [eV]"``."""
    match = re.fullmatch(r"\s*(.*?)\s*\[(.*)\]\s*", label or "")
    if match:
        return match.group(1), match.group(2)
    return (label or "").strip(), ""


def _scan_dim(label: str, axis, table: dict, index: int = 0,
              counts: tuple = ()) -> Dim:
    """What a scanned SES dimension is, from its label -- and, for a
    manipulator scan, from the run-mode table: ``index`` is this
    dimension's place in ``counts``, the scanned dimensions' sizes slowest
    first, which is the order the table lists its points in."""
    name, unit = _split_label(label)
    folded = name.casefold()
    if "photon" in folded or folded in {"hv", "excitation energy"}:
        return Dim(name or "Photon energy", unit or "eV", axis,
                   "photon_energy", f"Photon energy ({ENERGY_UNIT})")
    if folded.startswith("theta") or "deflector" in folded:
        return Dim(name, unit or "deg", axis, "angle",
                   f"Deflector {name} ({ANGLE_UNIT})")
    if folded in {"point", "position"}:
        motor = _motor_for(table, label) or _motor_for(
            table, name, axis, index, counts)
        if motor is not None:
            pretty, motor_unit, readback = motor
            if readback is not None:
                axis = straighten(readback)
            role, text = motor_role(pretty)
            if motor_unit == "mm":
                role, text = "position", f"{pretty} ({motor_unit})"
            return Dim(pretty, motor_unit, axis, role, text or "")
        return Dim(name, "" if unit == "unit" else unit, axis, "position")
    role, text = motor_role(name)
    return Dim(name or "Scan", unit, axis, role, text or "")


def _info(sections: dict) -> dict:
    """The metadata worth showing, from an SES header."""
    info = {}
    ses = _keys(sections.get("SES"))
    for key, value in ses.items():
        number = _float(value)
        info[f"ses.{key}"] = number if number is not None \
            and re.fullmatch(r"[-+0-9.eE]+", value) else value
    for key, value in _keys(sections.get("Manipulator")).items():
        number = _float(value)
        info[f"manipulator.{key}"] = number if number is not None else value
    if ses.get("Spectrum Name"):
        info["title"] = ses["Spectrum Name"]
    if ses.get("Date"):
        info["start_time"] = f"{ses['Date']} {ses.get('Time', '')}".strip()
    hv = _float(ses.get("Excitation Energy"))
    if hv is not None:
        info["photon_energy_eV"] = hv
    for key, out in (("Pass Energy", "analyser.pass_energy"),
                     ("Lens Mode", "analyser.lens_mode"),
                     ("Acquisition Mode", "analyser.acquisition_mode")):
        if key in ses:
            info[out] = _float(ses[key]) if key == "Pass Energy" \
                else ses[key]
    return info


# --------------------------------------------------------------------------
# A memory-mapped array as NexusView wants it
# --------------------------------------------------------------------------
class _Source:
    """Wraps a (memory-mapped) array so that every read is a copy in native
    float32 -- a view of a read-only map must not leak into code that
    writes to what it gets."""

    def __init__(self, array: np.ndarray, name: str):
        self._array = array
        self.name = name
        self.shape = tuple(array.shape)
        self.ndim = array.ndim
        self.dtype = np.dtype(np.float32)

    def __getitem__(self, key):
        return np.array(self._array[key], dtype=np.float32)


def _energy(axis: np.ndarray, label: str, array: np.ndarray):
    """``(energy, scale, array)``: a binding axis turned into ``E - E_F``,
    ascending, with the array's last axis flipped to match (a view)."""
    name, _unit = _split_label(label)
    axis = np.asarray(axis, dtype=float)
    scale = "kinetic"
    if "binding" in name.casefold():
        axis = -axis
        scale = "binding"
    if axis.size > 1 and axis[0] > axis[-1]:
        axis = axis[::-1].copy()
        array = array[..., ::-1]
    return axis, scale, array


def _layout(array, name, energy_axis, energy_label, angle_axis, dims,
            info) -> Layout:
    energy, scale, array = _energy(energy_axis, energy_label, array)
    energies = None
    hv_dims = [d for d in dims if d.role == "photon_energy"]
    if scale == "binding":
        info["energy_reference"] = "binding energy from SES, negated"
    elif len(dims) == 1 and hv_dims:
        # CIS on a kinetic scale: one binding window, so each hv's kinetic
        # axis is the first one moved by the change in hv.
        hv = hv_dims[0].axis
        energies = energy[None, :] + (hv - hv[0])[:, None]
    return Layout(dataset=_Source(array, name),
                  perm=tuple(range(array.ndim)), dims=dims,
                  angle=np.asarray(angle_axis, dtype=float), energy=energy,
                  facility=FACILITY, prefix=PREFIX, energy_scale=scale,
                  energies=energies, info=info)


# --------------------------------------------------------------------------
# .ibw
# --------------------------------------------------------------------------
def is_ses_ibw(path: str) -> bool:
    if not path.lower().endswith(".ibw"):
        return False
    wave = read_ibw(path)
    return "[SES]" in wave.note and len(wave.shape) >= 2


def read_layout_ibw(path: str) -> Layout:
    wave = read_ibw(path)
    if len(wave.shape) < 2:
        raise ValueError("a one-dimensional wave is not a spectrum")
    sections = parse_ses_note(wave.note)
    info = _info(sections)
    info["ibw.wave"] = wave.name
    table = _run_mode_table(sections.get("Run Mode Information"))
    # C order: (chunks, layers, columns, rows) -> scanned axes slowest first.
    scanned = list(reversed(range(2, len(wave.shape))))
    counts = tuple(wave.shape[i] for i in scanned)
    dims = [_scan_dim(wave.units[i], wave.axis(i), table, n, counts)
            for n, i in enumerate(scanned)]
    manipulator = [d for i, d in zip(scanned, dims)
                   if _split_label(wave.units[i])[0].casefold()
                   in {"point", "position"} and d.axis.size > 1]
    if len(manipulator) == 2:
        # A two-motor manipulator scan is a spatial scan, even when one of
        # the motors is an angle (One Square steps the polar against Z).
        for dim in manipulator:
            dim.role = "position"
    return _layout(wave.array(), wave.name, wave.axis(0), wave.units[0],
                   wave.axis(1), dims, info)


# --------------------------------------------------------------------------
# .zip
# --------------------------------------------------------------------------
def _ini(text: str) -> configparser.ConfigParser:
    parser = configparser.ConfigParser(strict=False, interpolation=None)
    parser.optionxform = str
    parser.read_string(text)
    return parser


def _zip_regions(z: zipfile.ZipFile) -> list:
    """``[(region name, bin member, spectrum ini, SES ini)]``."""
    names = set(z.namelist())
    regions = []
    if "viewer.ini" in names:
        viewer = _ini(z.read("viewer.ini").decode("latin-1"))
        for region in viewer.get("viewer", "region_list",
                                 fallback="").split():
            section = f"viewer.{region}"
            if not viewer.has_section(section):
                continue
            title = viewer.get(section, "name", fallback=region)
            ses_ini = viewer.get(section, "ini_path", fallback=f"{title}.ini")
            for channel in viewer.get(section, "channel_list",
                                      fallback="").split():
                path = viewer.get(f"{section}.{channel}", "path",
                                  fallback=f"Spectrum_{title}.bin")
                if path in names:
                    regions.append((title, path,
                                    path[:-4] + ".ini", ses_ini))
                    break
    if not regions:
        for member in sorted(names):
            match = re.fullmatch(r"Spectrum_(.+)\.bin", member)
            if match:
                regions.append((match.group(1), member,
                                f"Spectrum_{match.group(1)}.ini",
                                f"{match.group(1)}.ini"))
    return regions


def is_ses_zip(path: str) -> bool:
    if not path.lower().endswith(".zip") or not zipfile.is_zipfile(path):
        return False
    with zipfile.ZipFile(path) as z:
        return bool(_zip_regions(z))


def _member_array(path: str, z: zipfile.ZipFile, member: str, shape):
    """A ``.bin`` member as float32 of ``shape``: mapped in place when it is
    stored uncompressed, else read."""
    item = z.getinfo(member)
    count = int(np.prod(shape))
    if item.compress_type == zipfile.ZIP_STORED and item.file_size == 4 * count:
        with open(path, "rb") as f:
            f.seek(item.header_offset)
            local = f.read(30)
        name_len, extra_len = struct.unpack("<HH", local[26:30])
        offset = item.header_offset + 30 + name_len + extra_len
        return np.memmap(path, dtype="<f4", mode="r", offset=offset,
                         shape=shape)
    raw = z.read(member)
    if len(raw) < 4 * count:
        raise ValueError(f"{member} holds {len(raw)} bytes, "
                         f"{4 * count} expected")
    return np.frombuffer(raw, dtype="<f4", count=count).reshape(shape)


def read_layouts_zip(path: str) -> dict:
    """``{region name: Layout}`` for an SES map zip."""
    layouts = {}
    with zipfile.ZipFile(path) as z:
        names = set(z.namelist())
        for title, member, spectrum_ini, ses_ini in _zip_regions(z):
            if spectrum_ini not in names:
                continue
            spec = _ini(z.read(spectrum_ini).decode("latin-1"))
            section = "spectrum" if spec.has_section("spectrum") \
                else spec.sections()[0]
            get = lambda key, default=None: spec.get(section, key,  # noqa: E731
                                                      fallback=default)
            sizes = [int(get(k, 1)) for k in ("width", "height", "depth")]
            axes = [float(get(f"{k}offset", 0)) + float(get(f"{k}delta", 1))
                    * np.arange(n) for k, n in zip(("width", "height",
                                                    "depth"), sizes)]
            labels = [get(f"{k}label", "") for k in ("width", "height",
                                                     "depth")]
            sections = parse_ses_note(z.read(ses_ini).decode("latin-1")) \
                if ses_ini in names else {}
            info = _info(sections)
            info.setdefault("title", title)
            info["zip.region"] = title
            array = _member_array(path, z, member,
                                  (sizes[2], sizes[1], sizes[0]))
            dims = [_scan_dim(labels[2], axes[2], {})]
            layouts[title] = _layout(array, member, axes[0], labels[0],
                                     axes[1], dims, info)
    return layouts


def ses_header(path: str) -> dict:
    """The ``[SES]`` keys of an SES ``.ibw`` or ``.zip`` -- ``Location``,
    ``Instrument``, ... -- or of the first region of a ``.txt``."""
    lower = path.lower()
    if lower.endswith(".txt"):
        return text_header(path)
    if lower.endswith(".ibw"):
        note = read_ibw(path).note
    else:
        with zipfile.ZipFile(path) as z:
            names = set(z.namelist())
            note = next((z.read(ini).decode("latin-1")
                         for _t, _m, _s, ini in _zip_regions(z)
                         if ini in names), "")
    return _keys(parse_ses_note(note).get("SES"))


def ses_location(path: str) -> str:
    """The ``Location`` an SES file names in its header (``Bloch`` at MAX
    IV), or ``""``."""
    return ses_header(path).get("Location", "")


# --------------------------------------------------------------------------
# .txt
# --------------------------------------------------------------------------
def is_ses_text(path: str) -> bool:
    if not path.lower().endswith(".txt"):
        return False
    from arpes_viewer.loader.beamlines.cassiopee import CassiopeeLoader
    return CassiopeeLoader().can_open(path)


def text_header(path: str) -> dict:
    """The ``key=value`` lines of a text export's first ``[Info n]``
    section (where SES writes ``Location``, ``Instrument``, ...), read from
    the head of the file only -- a photon-energy scan is many megabytes."""
    with open(path, "rb") as handle:
        head = handle.read(65536).decode("latin-1")
    out, inside = {}, False
    for line in head.splitlines():
        stripped = line.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            if inside:
                break
            inside = bool(re.fullmatch(r"\[Info \d+\]", stripped))
            continue
        if inside and "=" in stripped:
            key, value = stripped.split("=", 1)
            out.setdefault(key.strip(), value.strip())
    return out


def _photon_energies(regions) -> np.ndarray:
    return np.array([_float(r.info.get("photon_energy_eV")) or np.nan
                     for r in regions], dtype=float)


def is_photon_energy_scan(regions) -> bool:
    """Several regions, each at its own photon energy."""
    hv = _photon_energies(regions)
    return len(regions) > 1 and np.all(np.isfinite(hv)) \
        and np.unique(np.round(hv, 4)).size == hv.size


def read_layout_text(path: str, regions=None) -> Layout:
    """A multi-region text export taken at one photon energy per region, as
    the layout of a kz map: ``(hv, angle, energy)`` with each region's own
    kinetic axis in ``energies``, which :func:`build_scan` refers to E_F
    region by region."""
    from arpes_viewer.loader.beamlines.cassiopee import parse_scienta
    regions = parse_scienta(path) if regions is None else regions
    if not is_photon_energy_scan(regions):
        raise ValueError(f"{os.path.basename(path)} is not a photon-energy "
                         f"scan: its regions do not each have their own hv")
    angle = regions[0].angle
    if any(r.angle.size != angle.size for r in regions):
        raise ValueError(f"{os.path.basename(path)}: the regions do not "
                         f"share one angle axis")
    count = min(r.energy.size for r in regions)
    info = {key: value for key, value in regions[0].info.items()
            if not key.startswith("scienta.Region")}
    if any(r.energy.size != count for r in regions):
        info["energy_trimmed"] = (f"regions cut to the {count} energy points "
                                  f"they all have")
    cube = np.stack([np.asarray(r.values[:, :count], dtype=np.float32)
                     for r in regions])
    energies = np.stack([np.asarray(r.energy[:count], dtype=float)
                         for r in regions])
    hv = _photon_energies(regions)
    info["title"] = info.get("scienta.Spectrum Name") or regions[0].name
    if info.get("date"):
        info["start_time"] = f"{info['date']} {info.get('time', '')}".strip()
    info["ses.regions"] = len(regions)
    dims = [Dim("Photon energy", "eV", hv, "photon_energy",
                f"Photon energy ({ENERGY_UNIT})")]
    return Layout(dataset=_Source(cube, os.path.basename(path)),
                  perm=(0, 1, 2), dims=dims, angle=np.asarray(angle, float),
                  energy=energies[0], facility=FACILITY, prefix=PREFIX,
                  energy_scale="kinetic", energies=energies, info=info)


# --------------------------------------------------------------------------
# The loaders
# --------------------------------------------------------------------------
class SESLoader(Loader):
    """Scienta SES ``.ibw`` and ``.zip`` -- and ``.txt`` when
    ``reads_text`` is set -- for one facility.

    A subclass sets ``name``, ``facility`` and ``prefix`` and, when it
    should only claim its own files, :meth:`claims`, which is handed the
    file's ``[SES]`` keys: SES writes the ``Location`` and ``Instrument``
    set up on the acquisition PC, which is how two beamlines writing the
    same format are told apart."""

    facility = FACILITY
    prefix = PREFIX
    patterns = ("*.ibw", "*.zip")
    reads_text = False

    def claims(self, header: dict) -> bool:
        return True

    def can_open(self, path: str) -> bool:
        try:
            lower = path.lower()
            if lower.endswith(".txt"):
                if not (self.reads_text and is_ses_text(path)):
                    return False
            elif not (is_ses_ibw(path) or is_ses_zip(path)):
                return False
            return self.claims(ses_header(path))
        except Exception:
            return False

    def _own(self, layout: Layout) -> Layout:
        layout.facility, layout.prefix = self.facility, self.prefix
        return layout

    def list_entries(self, path: str) -> list:
        base = os.path.splitext(os.path.basename(path))[0]
        try:
            lower = path.lower()
            if lower.endswith(".txt"):
                return self._list_text(path, base)
            if lower.endswith(".zip"):
                layouts = read_layouts_zip(path)
                rows = []
                for region, layout in layouts.items():
                    name = base if len(layouts) == 1 else f"{base} {region}"
                    for row in list_layout(layout, name,
                                           layout.info.get("title"),
                                           layout.info.get("start_time")):
                        row["entry"] = f"region {region}"
                        rows.append(row)
                return rows
            layout = read_layout_ibw(path)
            return list_layout(layout, base, layout.info.get("title"),
                               layout.info.get("start_time"))
        except Exception:
            return []

    def load(self, path: str, entry: str = None, progress=None):
        lower = path.lower()
        if lower.endswith(".txt"):
            return self._load_text(path, entry)
        if lower.endswith(".zip"):
            layouts = read_layouts_zip(path)
            if not layouts:
                raise ValueError(f"no SES region in {path}")
            if entry and entry.startswith("region "):
                region = entry[len("region "):]
                if region not in layouts:
                    raise KeyError(f"no dataset {entry!r} in {path}")
                layout = layouts[region]
            else:
                layout = next(iter(layouts.values()))
            return build_scan(self._own(layout), None)
        layout = self._own(read_layout_ibw(path))
        try:
            return build_scan(layout, entry)
        except KeyError:
            raise KeyError(f"no dataset {entry!r} in {path}") from None

    # -- text exports ------------------------------------------------------
    def _list_text(self, path: str, base: str) -> list:
        from arpes_viewer.loader.beamlines.cassiopee import parse_scienta
        regions = parse_scienta(path)
        if is_photon_energy_scan(regions):
            hv = _photon_energies(regions)
            return [{"entry": "hv scan", "kind": KIND_LABELS["kz_map"],
                     "name": f"{base} {hv.min():g}-{hv.max():g}eV",
                     "title": regions[0].name,
                     "start_time": regions[0].info.get("date")}]
        return [{"entry": f"region {i + 1}" if len(regions) > 1 else "cut",
                 "kind": KIND_LABELS["cut"], "name": region.name or base,
                 "title": region.name,
                 "start_time": region.info.get("date")}
                for i, region in enumerate(regions)]

    def _load_text(self, path: str, entry: str = None):
        from arpes_viewer.loader.beamlines.cassiopee import parse_scienta
        regions = parse_scienta(path)
        if not regions:
            raise ValueError(f"{os.path.basename(path)}: no spectrum in it")
        if entry in (None, "", "hv scan") and is_photon_energy_scan(regions):
            return build_scan(self._own(read_layout_text(path, regions)))
        index = 0
        if entry and entry.startswith("region "):
            index = int(entry.split()[1]) - 1
        region = regions[min(index, len(regions) - 1)]
        scan = NxsScan(kind="cut", filename_prefix=self.prefix)
        scan.x, scan.y = region.angle, region.energy
        scan.value = region.values
        scan.labels = {"x": region.angle_label, "y": region.energy_label}
        scan.info.update(region.info)
        scan.info["facility"] = self.facility
        return scan
