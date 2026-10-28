"""
tools/output/igor.py
====================
Datasets written for Igor Pro: an Igor text file (``.itx``) holding any
number of waves, or one Igor binary wave (``.ibw``, version 5) per dataset.

The MATLAB tools did this with ``Export2ScientaData.m`` /
``ExportAsScientaData.m``, which wrote a Scienta-style text block that Igor
then had to be taught to read (and which only took a 2-D cut). Writing Igor's
own formats instead means a double-click in Igor opens the data with its
axes already scaled, its units set, its dimension labels named and the
dataset's metadata in the wave note -- for a curve, a cut, a 3-D map or a
4-D spatial scan alike.

What is kept, and how
---------------------
* **Values** as double precision (``.itx``) or as the array's own float
  type (``.ibw``: float32 stays float32, everything else becomes float64).
* **Axes** as Igor's wave scaling (start and step per dimension). Igor's
  scaling is linear; an axis that is not evenly spaced (beyond 0.1 % of a
  step) is written as a wave of its own beside the data (``<name>_x``...),
  and the scaling then spans its first to last value, so the image is
  still drawn the right size.
* **Units** taken from the axis label's parentheses ("Kinetic energy (eV)"
  -> ``eV``), the rest of the label as the dimension's name.
* **Metadata** (the dataset's ``info``) as ``key=value`` lines of the wave
  note, the way the Scienta files Igor users already know carry theirs.

Array order is this program's (constructor order for the axes, the array as
it is stored), so a cut is (angle, energy) -- angle along Igor's rows --
exactly as it is drawn here.

Qt-free, so it can be tested headlessly.
"""
from __future__ import annotations

import re
import struct
import time

import numpy as np

__all__ = ["igor_name", "split_label", "wave_from_dataset", "itx_text",
           "write_itx", "ibw_bytes", "write_ibw", "IgorWaveSpec", "axis_specs",
           "write_datasets"]

#: Igor's wave names are at most 31 characters in a version-5 file.
MAX_NAME = 31
#: How uneven an axis may be, as a fraction of its mean step, and still be
#: written as plain scaling.
UNEVEN = 1e-3

# 1904-01-01, Igor's epoch, in Unix seconds.
_IGOR_EPOCH = -2082844800


def igor_name(text: str, taken=()) -> str:
    """A legal Igor wave name for ``text``: letters, digits and ``_``,
    starting with a letter, at most 31 characters, and not in ``taken``."""
    name = re.sub(r"[^A-Za-z0-9_]+", "_", str(text or "")).strip("_") or "wave"
    if not name[0].isalpha():
        name = "w" + name
    name = name[:MAX_NAME]
    base, n = name, 1
    while name in taken:
        n += 1
        suffix = f"_{n}"
        name = base[:MAX_NAME - len(suffix)] + suffix
    return name


def split_label(label: str):
    """``("Kinetic energy", "eV")`` from ``"Kinetic energy (eV)"``; the unit
    is empty when the label has no trailing parentheses."""
    label = str(label or "").strip()
    match = re.match(r"^(.*?)\s*\(([^()]*)\)\s*$", label)
    if match:
        return match.group(1).strip(), match.group(2).strip()
    return label, ""


def _ascii_unit(unit: str) -> str:
    """A unit Igor's legacy (Windows-1252) text can hold: Å stays, the
    superscript minus-one becomes ``^-1``."""
    return (unit.replace("⁻¹", "^-1").replace("⁻", "^-")
            .replace("¹", "1").replace("²", "2").replace("³", "3"))


class IgorWaveSpec:
    """One wave as Igor will see it: name, values, and per dimension its
    start, step, unit and label, plus the explicit axis values of any
    dimension that is not evenly spaced."""

    def __init__(self, name, values, axes, labels, note="", data_unit=""):
        self.name = name
        self.values = np.asarray(values)
        self.dims = []
        for axis, label in zip(axes, labels):
            axis = np.asarray(axis, dtype=float).ravel()
            title, unit = split_label(label)
            n = axis.size
            start = float(axis[0]) if n else 0.0
            step = float((axis[-1] - axis[0]) / (n - 1)) if n > 1 else 1.0
            even = True
            if n > 2:
                steps = np.diff(axis)
                even = bool(np.all(np.abs(steps - step)
                                   <= UNEVEN * max(abs(step), 1e-300)))
            self.dims.append({"start": start, "step": step or 1.0,
                              "unit": unit, "label": title,
                              "axis": None if even else axis})
        self.note = note
        self.data_unit = data_unit


def _note_text(info: dict, header: str = "") -> str:
    lines = [header] if header else []
    for key in sorted(info or {}):
        value = info[key]
        if isinstance(value, np.ndarray):
            if value.size > 64:
                continue
            value = " ".join(f"{v:.8g}" if isinstance(v, (float, np.floating))
                             else str(v) for v in value.ravel())
        elif isinstance(value, (float, np.floating)):
            value = f"{float(value):.10g}"
        text = str(value).replace("\r", " ").replace("\n", " ")
        lines.append(f"{key}={text}")
    return "\r".join(lines)


def wave_from_dataset(item: dict, taken=()) -> IgorWaveSpec:
    """An :class:`IgorWaveSpec` from the plain dataset form
    (:func:`loader.session.scan_to_dict`: name, kind, axes, labels, value,
    info).

    The axes come in constructor order; a 4-D spatial scan is stored as
    (y, x, k, E), so its first two are swapped to match the array.
    """
    from loader.nxs_file import axis_slots

    kind = item.get("kind", "")
    value = np.asarray(item["value"])
    axes = list(item["axes"])
    slots = list(axis_slots(kind, "constructor")) or ["x", "y", "z", "t"][:value.ndim]
    array_slots = list(axis_slots(kind, "array")) or slots
    by_slot = dict(zip(slots, axes))
    axes = [by_slot[s] for s in array_slots] if set(array_slots) <= set(by_slot) else axes
    labels = dict(item.get("labels") or {})
    names = [labels.get(s, s) for s in array_slots]
    name = igor_name(item.get("name", "wave"), taken)
    header = (f"Written by ARPES viewer on {time.strftime('%Y-%m-%d %H:%M')}"
              f"\rkind={kind}\rdataset={item.get('name', '')}")
    if kind in ("edc", "mdc", "spin_edc"):
        channels = (item.get("info") or {}).get("curve.channels")
        if channels is not None:
            header += "\rchannels=" + ";".join(str(c) for c in np.ravel(channels))
    return IgorWaveSpec(name, value, axes[:value.ndim], names[:value.ndim],
                        _note_text(item.get("info") or {}, header), "counts")


# --------------------------------------------------------------------------
# Igor text (.itx)
# --------------------------------------------------------------------------
def _quote(text: str) -> str:
    return str(text).replace("\\", "\\\\").replace('"', '\\"')


def _itx_wave(spec: IgorWaveSpec) -> list:
    values = np.asarray(spec.values, dtype=np.float64)
    shape = values.shape
    out = []
    dims = ",".join(str(n) for n in shape)
    out.append(f"WAVES/D/N=({dims}) '{spec.name}'")
    out.append("BEGIN")
    if values.ndim == 1:
        out.extend(_row(v) for v in values)
    else:
        # Igor reads a line as one row (fixed first index) across the
        # columns; higher dimensions come as blocks, layer after layer,
        # separated by a blank line.
        flat = values.reshape(shape[0], shape[1], -1)
        for layer in range(flat.shape[2]):
            block = flat[:, :, layer]
            out.extend("\t".join(_num(v) for v in row) for row in block)
            if flat.shape[2] > 1:
                out.append("")
    out.append("END")
    commands = []
    letters = "xyzt"
    for i, dim in enumerate(spec.dims[:4]):
        unit = _quote(dim["unit"])
        if dim["axis"] is not None:
            first, last = float(dim["axis"][0]), float(dim["axis"][-1])
            commands.append(f"SetScale/I {letters[i]} {first!r},{last!r},"
                            f"\"{unit}\", '{spec.name}'")
        else:
            commands.append(f"SetScale/P {letters[i]} {float(dim['start'])!r},"
                            f"{float(dim['step'])!r},\"{unit}\", '{spec.name}'")
        if dim["label"]:
            commands.append(f"SetDimLabel {i},-1,{igor_name(dim['label'])},"
                            f"'{spec.name}'")
    commands.append(f"SetScale d 0,0,\"{_quote(spec.data_unit)}\", '{spec.name}'")
    for command in commands:
        out.append("X " + command)
    if spec.note:
        # One Note command per line keeps each command short; Igor joins
        # them with returns ("Note" appends, separated by a CR).
        for line in spec.note.split("\r"):
            out.append(f"X Note '{spec.name}', \"{_quote(line)}\"")
    return out


def _num(value) -> str:
    value = float(value)
    if np.isnan(value):
        return "NaN"
    if np.isinf(value):
        return "INF" if value > 0 else "-INF"
    return f"{value:.9g}"


def _row(value) -> str:
    return _num(value)


def axis_specs(spec: IgorWaveSpec) -> list:
    """The explicit axis waves ``spec`` needs: one per dimension that is not
    evenly spaced, named ``<name>_x``, ``_y``, ``_z``, ``_t``."""
    out = []
    for i, dim in enumerate(spec.dims):
        if dim["axis"] is not None:
            out.append(IgorWaveSpec(
                igor_name(f"{spec.name[:28]}_{'xyzt'[i]}"), dim["axis"],
                [np.arange(dim["axis"].size)], [""],
                f"axis {i} of {spec.name}: {dim['label']}", dim["unit"]))
    return out


def itx_text(specs) -> str:
    """The text of an ``.itx`` file holding every wave in ``specs``, each
    followed by its axis waves for unevenly spaced dimensions."""
    lines = ["IGOR"]
    for spec in specs:
        lines.extend(_itx_wave(spec))
        for axis_spec in axis_specs(spec):
            lines.extend(_itx_wave(axis_spec))
    return "\r\n".join(lines) + "\r\n"


def write_itx(path: str, specs) -> str:
    """Write ``specs`` to ``path`` as an Igor text file (UTF-8, with the
    byte-order mark Igor 7 and later use to recognise UTF-8)."""
    with open(path, "w", encoding="utf-8-sig", newline="") as handle:
        handle.write(itx_text(specs))
    return path


# --------------------------------------------------------------------------
# Igor binary wave, version 5 (.ibw)
# --------------------------------------------------------------------------
_BIN_HEADER = 64
_WAVE_HEADER = 320
_NT_FP32, _NT_FP64 = 2, 4
_LABEL = 32          # MAX_DIM_LABEL_CHARS + 1 in a version-5 file


def ibw_bytes(spec: IgorWaveSpec) -> bytes:
    """The bytes of a version-5 Igor binary wave (little-endian)."""
    values = np.asarray(spec.values)
    if values.ndim > 4:
        raise ValueError("Igor waves have at most four dimensions")
    if values.dtype == np.float32:
        type_code, data = _NT_FP32, values.astype("<f4")
    else:
        type_code, data = _NT_FP64, values.astype("<f8")
    # Igor stores a wave column-major: the first index runs fastest.
    payload = np.asfortranarray(data).ravel(order="F").tobytes()
    npnts = int(values.size)

    enc = "cp1252"

    def encode(text):
        return _ascii_unit(str(text)).encode(enc, errors="replace")

    note = encode(spec.note.replace("\n", "\r"))
    data_eunits = encode(spec.data_unit) if len(spec.data_unit) > 3 else b""
    dim_eunits, dim_units, dim_labels = [], [], []
    for i in range(4):
        dim = spec.dims[i] if i < len(spec.dims) and i < values.ndim else None
        unit = encode(dim["unit"]) if dim else b""
        if len(unit) > 3:
            dim_eunits.append(unit)
            dim_units.append(b"")
        else:
            dim_eunits.append(b"")
            dim_units.append(unit)
        if dim and dim["label"]:
            block = bytearray(_LABEL * (values.shape[i] + 1))
            label = encode(dim["label"])[:_LABEL - 1]
            block[:len(label)] = label
            dim_labels.append(bytes(block))
        else:
            dim_labels.append(b"")

    wfm_size = _WAVE_HEADER + len(payload)
    bin_header = bytearray(struct.pack(
        "<hhiiii4i4iiii", 5, 0, wfm_size, 0, len(note), len(data_eunits),
        *[len(u) for u in dim_eunits], *[len(b) for b in dim_labels], 0, 0, 0))
    assert len(bin_header) == _BIN_HEADER

    now = int(time.time()) - _IGOR_EPOCH
    wave = bytearray(_WAVE_HEADER)
    struct.pack_into("<IIi", wave, 4, now & 0xFFFFFFFF, now & 0xFFFFFFFF, npnts)
    struct.pack_into("<h", wave, 16, type_code)
    struct.pack_into("<h", wave, 26, 1)                       # whVersion
    name = spec.name.encode("ascii", errors="replace")[:MAX_NAME]
    wave[28:28 + len(name)] = name
    n_dim = list(values.shape) + [0] * (4 - values.ndim)
    struct.pack_into("<4i", wave, 68, *n_dim)
    sf_a, sf_b = [1.0] * 4, [0.0] * 4
    for i, dim in enumerate(spec.dims[:values.ndim]):
        if dim["axis"] is not None and dim["axis"].size > 1:
            sf_b[i] = float(dim["axis"][0])
            sf_a[i] = float((dim["axis"][-1] - dim["axis"][0])
                            / (dim["axis"].size - 1))
        else:
            sf_a[i], sf_b[i] = dim["step"], dim["start"]
    struct.pack_into("<4d", wave, 84, *sf_a)
    struct.pack_into("<4d", wave, 116, *sf_b)
    data_unit = encode(spec.data_unit) if len(spec.data_unit) <= 3 else b""
    wave[148:148 + len(data_unit)] = data_unit
    for i, unit in enumerate(dim_units):
        wave[152 + 4 * i:152 + 4 * i + len(unit)] = unit

    # The checksum makes the 16-bit sum of both headers zero.
    head = bytes(bin_header) + bytes(wave)
    total = sum(struct.unpack(f"<{len(head) // 2}h", head)) & 0xFFFF
    struct.pack_into("<h", bin_header, 2,
                     struct.unpack("<h", struct.pack("<H", (-total) & 0xFFFF))[0])
    return (bytes(bin_header) + bytes(wave) + payload + note + data_eunits
            + b"".join(dim_eunits) + b"".join(dim_labels))


def write_ibw(path: str, spec: IgorWaveSpec) -> str:
    """Write one wave to ``path`` as a version-5 Igor binary wave."""
    with open(path, "wb") as handle:
        handle.write(ibw_bytes(spec))
    return path


def write_datasets(path: str, items, fmt: str = "itx") -> list:
    """Write datasets (the plain form of :func:`wave_from_dataset`) for
    Igor and return the files written.

    ``itx``: every dataset, and its axis waves, into the one file ``path``.
    ``ibw``: one file per wave -- ``path`` for a single dataset, otherwise
    ``<path stem>_<wave name>.ibw`` beside it -- and one more per unevenly
    spaced axis.
    """
    import os

    taken, specs = set(), []
    for item in items:
        spec = wave_from_dataset(item, taken)
        taken.add(spec.name)
        specs.append(spec)
    if fmt == "itx":
        return [write_itx(path, specs)]
    stem, _ext = os.path.splitext(path)
    written = []
    for spec in specs:
        target = path if len(specs) == 1 else f"{stem}_{spec.name}.ibw"
        written.append(write_ibw(target, spec))
        for axis in axis_specs(spec):
            written.append(write_ibw(
                f"{os.path.splitext(target)[0]}_{axis.name[-1]}axis.ibw", axis))
    return written
