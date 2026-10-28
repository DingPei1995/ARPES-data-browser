"""
tools/output/explog.py
======================
The experiment log: one CSV row per measurement, written from the data
files' own metadata instead of typed in by hand during the beamtime.

A beamline log (``DLS_20230508.xlsx``, the ANTARES template) is a table:
file name, data type, photon energy, polarisation, sample position, the
manipulator angles, pass energy, lens mode, centre energy, temperature...
Which columns, under which headings and in which order differs between
beamlines and between people, so the columns are a :class:`Template` the
user builds once from one dataset's fields -- its metadata (``info``), its
motor positions (``motors``) and a few fields worked out from the data
itself (file name, kind, the scanned axes' ranges) -- and that is then
applied to every dataset logged.

What a position column says depends on the measurement:

- a cut: the sample stood still, so one number;
- a spatial (SPEM) scan: the scanned stage axes as ``start/stop/step``;
- a map taken by turning the manipulator rather than with the deflector:
  the scanned angle as ``start/stop/step``, and the sample position, which
  the beamline may correct along with the angle to keep the spot on the
  sample, as ``first..last`` when it moved.

Each column therefore has a *scan rule* (:data:`SCAN_RULES`) saying when its
value is to be replaced by the range of a scanned axis; the default, "auto",
recognises a scanned motor by name and a moving one by its per-point values.

An existing log is appended to after its last record, never rewritten
above it, so a log kept open across a beamtime (and edited by hand in
between) keeps everything already in it. Qt-free, so it can be tested on its
own.
"""
from __future__ import annotations

import csv
import io
import json
import os
import re
from dataclasses import dataclass, field

import numpy as np

#: Fields worked out from the dataset rather than read from its metadata:
#: id -> (what the field list calls it, the default column heading).
DERIVED = {
    "@file": ("File name", "File name"),
    "@entry": ("Entry in the file", "Entry"),
    "@kind": ("Data type (cut, map, SPEM...)", "Data Type"),
    "@title": ("Title", "Title"),
    "@start": ("Start time", "Start time"),
    "@scan1": ("Scanned axis (map / SPEM X)", "Scan axis"),
    "@scan2": ("Second scanned axis (SPEM Y)", "Scan axis 2"),
    "@angle": ("Analyser angle range (slit)", "Slit angle"),
    "@energy": ("Energy window", "Energy window"),
    "@points": ("Number of points per axis", "Points"),
    "@blank": ("(empty column, filled in by hand)", "Comments"),
}

#: When a column's value is replaced by the range of a scanned axis:
#: rule -> what the rule list calls it.
SCAN_RULES = {
    "auto": "auto",
    "fixed": "never (fixed value)",
    "spatial_x": "SPEM: scanned X",
    "spatial_y": "SPEM: scanned Y",
    "angle_scan": "map: scanned angle",
    "hv_scan": "kz map: photon energy",
}

#: How a range is written: style -> (what the list calls it, the format).
RANGE_STYLES = {
    "slash": "start/stop/step",
    "colon": "start:step:stop",
}

#: Where each kind keeps its scanned axes, its analyser angle and its
#: energy, as attribute names on the scan.
SCANNED_SLOTS = {
    "map": ("x",), "k_map": ("x",), "kz_map": ("x",), "kz_map_k": ("x",),
    "spem_4d": ("x", "y"), "spem_1d": ("x",),
}
ANGLE_SLOT = {"cut": "x", "map": "k", "k_map": "k", "kz_map": "k",
              "kz_map_k": "k", "spem_4d": "k", "spem_1d": "y", "mdc": "x"}
ENERGY_SLOT = {"cut": "y", "map": "z", "k_map": "z", "kz_map": "z",
               "kz_map_k": "z", "spem_4d": "z", "spem_1d": "z", "edc": "x",
               "spin_edc": "x"}
SPEM = ("spem_4d", "spem_1d")
MAP_KINDS = ("map", "k_map")
HV_KINDS = ("kz_map", "kz_map_k")

#: What the kinds are called in the log's Data Type column.
KIND_NAMES = {
    "cut": "Cut", "map": "Map", "k_map": "k-map", "kz_map": "kz map",
    "kz_map_k": "kz map (k)", "spem_4d": "SPEM", "spem_1d": "SPEM (line)",
    "edc": "EDC", "mdc": "MDC", "spin_edc": "Spin EDC",
}

#: Words of a metadata key that are units, not the motor's name:
#: "manipulator.polar_deg" is the motor "polar".
_UNIT_WORDS = {"deg", "mm", "um", "ev", "k", "mev", "s", "ms", "a", "mbar"}

#: A loader's "<motor> range" entry: "-14 .. 14 deg (57 steps)".
_RANGE_TEXT = re.compile(
    r"^\s*(-?[\d.eE+-]+)\s*\.\.\s*(-?[\d.eE+-]+)\s*(.*?)\s*\((\d+) steps?\)\s*$")


@dataclass
class Column:
    """One column of the log: its heading, the field it is filled from,
    and when a scanned axis's range replaces that field's value."""
    header: str
    field: str
    scan: str = "auto"

    def to_json(self) -> dict:
        return {"header": self.header, "field": self.field, "scan": self.scan}

    @classmethod
    def from_json(cls, item: dict) -> "Column":
        scan = item.get("scan", "auto")
        return cls(str(item.get("header", "")), str(item.get("field", "")),
                   scan if scan in SCAN_RULES else "auto")


@dataclass
class Template:
    """The log's columns, in order, and how ranges are written."""
    columns: list = field(default_factory=list)
    range_style: str = "slash"

    def headers(self) -> list:
        return [column.header for column in self.columns]

    def to_json(self) -> dict:
        return {"columns": [c.to_json() for c in self.columns],
                "range_style": self.range_style}

    @classmethod
    def from_json(cls, stored) -> "Template":
        if not isinstance(stored, dict):
            return default_template()
        columns = [Column.from_json(item) for item in stored.get("columns", ())
                   if isinstance(item, dict) and item.get("field")]
        style = stored.get("range_style", "slash")
        return cls(columns, style if style in RANGE_STYLES else "slash")


def default_template() -> Template:
    """What a first log holds before anyone has chosen its columns: what
    every dataset has, plus a column for notes."""
    return Template([Column(DERIVED[f][1], f) for f in
                     ("@file", "@kind", "@start", "@scan1", "@scan2",
                      "@energy", "@blank")])


# --------------------------------------------------------------------------
# The template, remembered between sessions
# --------------------------------------------------------------------------
TEMPLATE_FILE = "explog_template.json"


def _config_path() -> str:
    folder = os.environ.get("ARPES_VIEWER_CONFIG_DIR") or os.path.join(
        os.path.expanduser("~"), ".arpes_viewer")
    return os.path.join(folder, TEMPLATE_FILE)


def load_template(path: str = None) -> Template:
    """The template kept at ``path`` (by default the one used last); the
    default template when there is none or it cannot be read."""
    try:
        with open(path or _config_path(), encoding="utf-8") as handle:
            return Template.from_json(json.load(handle))
    except (OSError, ValueError):
        return default_template()


def save_template(template: Template, path: str = None) -> None:
    """Keep ``template`` at ``path`` (by default for the next session).
    Raises OSError when ``path`` was given and cannot be written; the
    remembered copy is a convenience, so failing to write it is not."""
    target = path or _config_path()
    try:
        os.makedirs(os.path.dirname(target) or ".", exist_ok=True)
        with open(target, "w", encoding="utf-8") as handle:
            json.dump(template.to_json(), handle, indent=1)
    except OSError:
        if path is not None:
            raise


# --------------------------------------------------------------------------
# Values
# --------------------------------------------------------------------------
def number(value) -> str:
    """A number as a log shows it: six significant digits, no ``.0``."""
    value = float(value)
    if value == 0:
        return "0"
    return f"{value:.6g}"


def _uniform_step(values: np.ndarray):
    """The step of an evenly spaced axis, or None if it is not one."""
    if values.size < 2:
        return None
    steps = np.diff(values)
    step = (values[-1] - values[0]) / (values.size - 1)
    if step == 0:
        return None
    scale = max(abs(step), 1e-12)
    return float(step) if np.all(np.abs(steps - step) <= 0.02 * scale) else None


def axis_range(values, style: str = "slash") -> str:
    """An axis or a list of per-point values, as the log writes it.

    One value (or the same value throughout) is that value; an evenly
    spaced one is ``start/stop/step`` (or ``start:step:stop``); one that
    moved unevenly -- a sample position corrected along with an angle -- is
    ``first..last``.
    """
    values = np.asarray(values, dtype=float).ravel()
    values = values[np.isfinite(values)]
    if values.size == 0:
        return ""
    first, last = values[0], values[-1]
    if values.size == 1 or np.ptp(values) <= 1e-9 * max(1.0, abs(first)):
        return number(first)
    step = _uniform_step(values)
    if step is None:
        return f"{number(first)}..{number(last)}"
    if style == "colon":
        return f"{number(first)}:{number(step)}:{number(last)}"
    return f"{number(first)}/{number(last)}/{number(step)}"


def _range_from_text(text: str, data, style: str):
    """A loader's ``"-14 .. 14 deg (57 steps)"`` as a log range; None if
    ``text`` is not one.

    The text has only the ends and the count. When they are those of one
    of the dataset's scanned axes, it is that axis, written with its step;
    otherwise it is a motor that moved along with the scan (a sample
    position corrected as the manipulator turned), whose steps need not be
    even, so only where it started and ended is written."""
    match = _RANGE_TEXT.match(str(text))
    if not match:
        return None
    first, last, count = float(match[1]), float(match[2]), int(match[4])
    if count < 2 or first == last:
        return number(first)
    for slot in SCANNED_SLOTS.get(getattr(data, "kind", ""), ()):
        axis = _axis(data.scan, slot)
        if axis is not None and axis.size == count and np.allclose(
                [axis[0], axis[-1]], [first, last],
                rtol=1e-4, atol=1e-6 * max(1.0, abs(first), abs(last))):
            return axis_range(axis, style)
    return f"{number(first)}..{number(last)}"


def plain(value, style: str = "slash") -> str:
    """Any metadata value as a log cell."""
    if value is None:
        return ""
    if isinstance(value, bytes):
        value = value.decode("utf-8", "replace")
    if isinstance(value, (list, tuple)) and value and all(
            isinstance(v, (int, float, np.number)) for v in value):
        value = np.asarray(value, dtype=float)
    if isinstance(value, np.ndarray):
        if value.size == 0:
            return ""
        if np.issubdtype(value.dtype, np.number):
            return axis_range(value, style)
        return " ".join(str(v) for v in value.ravel())
    if isinstance(value, (bool, np.bool_)):
        return str(bool(value))
    if isinstance(value, (int, float, np.number)):
        return number(value)
    return " ".join(str(value).split())


def _words(name: str) -> str:
    """A motor or axis name, reduced to compare: ``"Polar (deg)"`` and
    ``"manipulator.polar_deg"`` both become ``"polar"``."""
    name = str(name).split("(")[0].strip()
    name = name.split(".")[-1]
    parts = [p for p in re.split(r"[^0-9a-zA-Z]+", name.lower()) if p]
    while len(parts) > 1 and parts[-1] in _UNIT_WORDS:
        parts.pop()
    return "".join(parts)


def _axis(scan, slot):
    values = getattr(scan, slot, None)
    if values is None:
        return None
    values = np.asarray(values, dtype=float).ravel()
    return values if values.size else None


def _scan_role(data) -> str:
    """What a map's first axis is: "angle", "photon_energy", "position",
    "temperature", "other"."""
    kind = getattr(data, "kind", "")
    if kind in HV_KINDS:
        return "photon_energy"
    role = (getattr(data.scan, "info", {}) or {}).get("axis0.role")
    return str(role) if role else "angle"


def _rule_range(rule: str, data, style: str):
    """The range a column's explicit scan rule gives this dataset, or None
    when the rule does not apply to it (then the field's own value is
    used)."""
    kind, scan = getattr(data, "kind", ""), data.scan
    if rule == "spatial_x" and kind in SPEM:
        axis = _axis(scan, "x")
    elif rule == "spatial_y" and kind == "spem_4d":
        axis = _axis(scan, "y")
    elif rule == "angle_scan" and kind in MAP_KINDS \
            and _scan_role(data) in ("angle", "other"):
        axis = _axis(scan, "x")
    elif rule == "hv_scan" and (kind in HV_KINDS or (
            kind == "map" and _scan_role(data) == "photon_energy")):
        measured = (scan.info or {}).get("axis0_measured")
        axis = np.asarray(measured, dtype=float).ravel() \
            if measured is not None else _axis(scan, "x")
    else:
        return None
    return axis_range(axis, style) if axis is not None else None


def _auto_range(key: str, data, style: str):
    """The range of the scanned axis that ``key`` names, if it names one:
    a loader's ``"<motor> range"`` motor entry, or a scanned axis whose
    label is the motor's name. None otherwise."""
    wanted = _words(key)
    if not wanted:
        return None
    scan = data.scan
    for name, value in (getattr(scan, "fourd_info", {}) or {}).items():
        if str(name).endswith(" range") and _words(str(name)[:-6]) == wanted:
            found = _range_from_text(value, data, style)
            if found is not None:
                return found
    labels = getattr(scan, "labels", {}) or {}
    for slot in SCANNED_SLOTS.get(getattr(data, "kind", ""), ()):
        label = labels.get(slot)
        axis = _axis(scan, slot)
        if label and axis is not None and _words(label) == wanted:
            return axis_range(axis, style)
    if getattr(data, "kind", "") in HV_KINDS and (
            "photonenergy" in wanted or wanted == "hv"):
        return _rule_range("hv_scan", data, style)
    return None


def _source(field_id: str, data):
    """``(found, value)`` of a metadata field on this dataset."""
    scan = data.scan
    if field_id.startswith("info:"):
        info = getattr(scan, "info", {}) or {}
        key = field_id[5:]
        return (key in info), info.get(key)
    if field_id.startswith("motor:"):
        motors = getattr(scan, "fourd_info", {}) or {}
        key = field_id[6:]
        return (key in motors), motors.get(key)
    return False, None


def _derived(field_id: str, data, record: dict, style: str) -> str:
    record = record or {}
    scan = data.scan
    info = getattr(scan, "info", {}) or {}
    kind = getattr(data, "kind", "")
    if field_id == "@file":
        return str(record.get("name") or os.path.basename(
            str(record.get("path") or getattr(data, "path", "") or "")))
    if field_id == "@entry":
        return str(record.get("entry") or "")
    if field_id == "@kind":
        name = KIND_NAMES.get(kind, kind)
        role = _scan_role(data) if kind == "map" else ""
        if role == "photon_energy":
            return "Map (photon energy)"
        if role in ("position", "temperature"):
            return f"Map ({role})"
        return name
    if field_id == "@title":
        return plain(info.get("title") or record.get("title"))
    if field_id == "@start":
        return plain(info.get("start_time") or record.get("start_time"))
    if field_id in ("@scan1", "@scan2"):
        slots = SCANNED_SLOTS.get(kind, ())
        index = 0 if field_id == "@scan1" else 1
        if index >= len(slots):
            return ""
        slot = slots[index]
        axis = _axis(scan, slot)
        if axis is None:
            return ""
        label = (getattr(scan, "labels", {}) or {}).get(slot, "")
        text = axis_range(axis, style)
        return f"{label}: {text}" if label else text
    if field_id == "@angle":
        axis = _axis(scan, ANGLE_SLOT.get(kind, ""))
        return axis_range(axis, style) if axis is not None else ""
    if field_id == "@energy":
        axis = _axis(scan, ENERGY_SLOT.get(kind, ""))
        if axis is None:
            return ""
        return axis_range(np.sort(axis) if axis[0] > axis[-1] else axis, style)
    if field_id == "@points":
        sizes = []
        for slot in ("x", "y", "k", "z"):
            axis = _axis(scan, slot)
            if axis is not None:
                sizes.append(str(axis.size))
        return "x".join(sizes)
    return ""


def cell(column: Column, data, record: dict = None,
         style: str = "slash") -> str:
    """What ``column`` holds for this dataset."""
    if column.field.startswith("@"):
        return _derived(column.field, data, record, style)
    if column.scan != "fixed":
        if column.scan == "auto":
            key = column.field.split(":", 1)[-1]
            found = _auto_range(key, data, style)
        else:
            found = _rule_range(column.scan, data, style)
        if found is not None:
            return found
    _present, value = _source(column.field, data)
    return plain(value, style)


def row(template: Template, data, record: dict = None) -> list:
    """The log row for one dataset, in the template's column order."""
    return [cell(column, data, record, template.range_style)
            for column in template.columns]


def available_fields(data, record: dict = None) -> list:
    """Every field a column can be filled from, for this dataset:
    ``[(field id, name to show, its value here), ...]`` -- the derived
    fields first, then the metadata, then the motors, each in the file's
    own order."""
    fields = [(fid, name, _derived(fid, data, record, "slash"))
              for fid, (name, _header) in DERIVED.items()]
    scan = data.scan
    for key, value in (getattr(scan, "info", {}) or {}).items():
        if str(key).startswith("_"):
            continue                       # the reader's own bookkeeping
        fields.append((f"info:{key}", str(key), plain(value)))
    for key, value in (getattr(scan, "fourd_info", {}) or {}).items():
        fields.append((f"motor:{key}", f"Motor.{key}", plain(value)))
    return fields


def default_header(field_id: str) -> str:
    """The heading a column starts with when its field is chosen."""
    if field_id in DERIVED:
        return DERIVED[field_id][1]
    key = field_id.split(":", 1)[-1]
    return key.split(".")[-1] if "." in key else key


# --------------------------------------------------------------------------
# The CSV file
# --------------------------------------------------------------------------
def _blank(line_cells) -> bool:
    return all(not str(c).strip() for c in line_cells)


@dataclass
class LogFile:
    """What is already in a log on disk: its last header row, the records
    under it, and whether it starts with a byte-order mark."""
    path: str
    header: list = field(default_factory=list)
    records: list = field(default_factory=list)
    exists: bool = False
    bom: bool = False

    def column(self, header: str) -> list:
        """The values already in the column headed ``header``."""
        if header not in self.header:
            return []
        index = self.header.index(header)
        return [r[index] for r in self.records if index < len(r)]


def read_log(path: str) -> LogFile:
    """Read an existing log (or describe one that does not exist yet).

    The header is the last row the program would have written as one: a log
    whose columns changed half way has a second header row, and the records
    after it are the ones that header describes."""
    log = LogFile(path)
    if not os.path.isfile(path):
        return log
    log.exists = True
    with open(path, "rb") as handle:
        raw = handle.read()
    log.bom = raw.startswith(b"\xef\xbb\xbf")
    text = _decode(raw)
    rows = [r for r in csv.reader(io.StringIO(text)) if not _blank(r)]
    if not rows:
        return log
    # A later header row is written when the template changed (append());
    # the records after the last one are the ones it describes.
    start = next((i for i in range(len(rows) - 1, -1, -1)
                  if rows[i] and rows[i][0] == _HEADER_MARK), None)
    header = rows[0] if start is None else rows[start][1:]
    while header and not header[-1].strip():
        header = header[:-1]          # a spreadsheet's trailing empty cells
    log.header = header
    log.records = rows[(0 if start is None else start) + 1:]
    return log


#: The first cell of a header row written part way down a log, when its
#: columns changed: it says what the row is, so it is not taken for a
#: measurement, and lets the next append find the header that applies.
_HEADER_MARK = "# columns"


def _decode(raw: bytes) -> str:
    for encoding in ("utf-8-sig", "cp1252", "latin-1"):
        try:
            return raw.decode(encoding)
        except UnicodeDecodeError:
            continue
    return raw.decode("utf-8", "replace")


def append(path: str, template: Template, rows: list) -> int:
    """Add ``rows`` to the log at ``path`` after its last record, creating
    it (with its header) if it does not exist. Returns how many were
    written.

    Blank lines after the last record -- a spreadsheet saving the file
    leaves rows of empty cells -- are dropped, so the new records follow the
    last one directly. When the file's header is not the template's, a new
    header row (marked ``# columns``) is written first, and everything
    already in the file stays as it was. Written to a temporary file and
    moved into place, so a failure cannot leave half a log.
    """
    headers = template.headers()
    existing = read_log(path)
    out = io.StringIO()
    writer = csv.writer(out, lineterminator="\r\n")
    if existing.exists:
        with open(path, "rb") as handle:
            text = _decode(handle.read())
        lines = text.splitlines()
        while lines and _blank(next(csv.reader([lines[-1]]), [])):
            lines.pop()
        body = "\r\n".join(lines)
        if body:
            body += "\r\n"
        out.write(body)
        if existing.header != headers:
            writer.writerow([_HEADER_MARK] + headers if lines else headers)
    else:
        writer.writerow(headers)
    for cells in rows:
        writer.writerow(cells)
    folder = os.path.dirname(os.path.abspath(path))
    os.makedirs(folder, exist_ok=True)
    temporary = os.path.join(folder, f".{os.path.basename(path)}.writing")
    # Excel reads a CSV as the system's code page unless it starts with a
    # byte-order mark, which turns "hν" and "Δ" into mojibake: new logs get
    # one, an existing log keeps whatever it had.
    bom = existing.bom or not existing.exists
    with open(temporary, "w", encoding="utf-8-sig" if bom else "utf-8",
              newline="") as handle:
        handle.write(out.getvalue())
    os.replace(temporary, path)
    return len(rows)


def logged_names(log: LogFile, template: Template) -> set:
    """The file names already logged, read from the log's file-name column
    (the template's ``@file`` column, under its heading); empty when the
    template has none."""
    for column in template.columns:
        if column.field == "@file":
            return {v.strip() for v in log.column(column.header) if v.strip()}
    return set()
