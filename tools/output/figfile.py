"""
tools/output/figfile.py
=======================
The figure composer's own file: a figure saved to be opened and edited
again, and the data it was made from kept beside it.

``fig1.afig`` is a zip archive holding

* ``figure.json`` -- the page: its style, its drawings, and every panel's
  settings (ranges, levels, colour map, lines, notes, arrows, overlays,
  insets, frames, letters) and the description of its data (labels, kind,
  where it came from);
* ``arrays.npz`` -- each panel's arrays (the image, its axes), exactly.

Nothing in it is executed when it is read: JSON and ``.npz`` with
``allow_pickle=False``.

Beside it, :func:`write_data_folder` makes ``data for fig1/`` with one
sub-folder per panel -- ``a``, ``b``... -- holding both

* ``processed/`` -- the data as the panel shows it: the image as a cut in
  the program's own .nxs (it opens in the viewer again), the lines drawn
  on it (calculated bands, contours, curves) as CSV, a picture as PNG;
* ``original/`` -- what it was made from, traced as far back as the
  datasets record it: the **original measurement file** (a k-converted,
  symmetrised, cut-out plane points back to the file it was first read
  from, through the ``origin.path`` its info carries), the calculation
  file with the files it is read with, or the picture file a panel was
  loaded from;

and ``sources.txt``, saying which entry of the file it was and the
processing between it and the panel.

Qt-free: a 3-D scene panel is turned into a picture by the caller before
saving (:func:`save` takes a ``flatten`` callback).
"""
from __future__ import annotations

import dataclasses
import io
import json
import os
import re
import shutil
import zipfile

import numpy as np

from tools.output import figure as nf

__all__ = ["EXTENSION", "FORMAT", "save", "load", "data_folder_of",
           "write_data_folder", "sources_from_info", "source_record"]

EXTENSION = ".afig"
FORMAT = "arpes-viewer-figure"
VERSION = 1

#: Panel fields that are tuples in the model (JSON hands lists back)
_TUPLE_FIELDS = ("levels", "x_range", "y_range", "frame")
#: the drawing dataclasses a panel holds, by name
_TYPES = {cls.__name__: cls for cls in (nf.Line, nf.TextNote, nf.Arrow,
                                         nf.Overlay, nf.ScaleBar, nf.Inset,
                                         nf.SecondaryAxis, nf.Shape)}


# --------------------------------------------------------------------------
# Where a panel came from
# --------------------------------------------------------------------------
def source_record(path: str, *, kind: str = "dataset", entry: str = "",
                  history=(), note: str = "", companions=()) -> dict:
    """One thing a panel was made from: a measurement file (and its
    entry), a calculation file, or a picture. ``companions`` are files
    the first needs beside it to be read again (a ``.win``, a POSCAR)."""
    return {"kind": kind, "path": os.path.abspath(path) if path else "",
            "entry": str(entry or ""), "history": [str(h) for h in history],
            "note": str(note or ""),
            "companions": [os.path.abspath(c) for c in companions]}


def _history(info: dict) -> list:
    steps = []
    for key, value in (info or {}).items():
        match = re.fullmatch(r"proc\.step\.(\d+)", str(key))
        if match:
            steps.append((int(match.group(1)), str(value)))
    return [text for _i, text in sorted(steps)]


def sources_from_info(info: dict, path: str = "") -> list:
    """The origin a dataset's info records -- ``origin.path`` and
    ``origin.entry``, set when the measurement was first read and carried
    through every processing step -- else ``path`` (the dataset's own
    file). ``[]`` when there is neither."""
    info = info or {}
    origin = str(info.get("origin.path", "") or "") or str(path or "")
    if not origin:
        return []
    return [source_record(origin, entry=str(info.get("origin.entry", "") or ""),
                          history=_history(info))]


def _json_safe(value):
    """``value`` as plain JSON: numbers, strings, lists, dicts."""
    if isinstance(value, dict):
        return {str(k): _json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(v) for v in value]
    if isinstance(value, np.ndarray):
        return _json_safe(value.tolist()) if value.size <= 64 else f"<array {value.shape}>"
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        value = float(value)
    if isinstance(value, float):
        # NaN breaks a drawn line in two; Python's json keeps it (as NaN)
        return value
    if isinstance(value, (str, int, bool)) or value is None:
        return value
    if isinstance(value, bytes):
        return value.decode("utf-8", "replace")
    return str(value)


# --------------------------------------------------------------------------
# Writing
# --------------------------------------------------------------------------
def _encode(value):
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        out = {"__type__": type(value).__name__}
        for f in dataclasses.fields(value):
            out[f.name] = _encode(getattr(value, f.name))
        return out
    if isinstance(value, (list, tuple)):
        return [_encode(v) for v in value]
    return _json_safe(value)


def _panel_dict(panel, index: int, arrays: dict, flatten=None) -> dict:
    data = panel.data
    if data is not None and data.kind == "scene3d":
        if flatten is None:
            raise ValueError("a 3-D view panel needs flattening to be saved")
        rgba = flatten(panel)
        data = nf.picture_data(rgba, data.name)
        data.info = dict(panel.data.info or {})
        data.sources = list(getattr(panel.data, "sources", []) or [])
    out = {}
    for f in dataclasses.fields(panel):
        if f.name == "data":
            continue
        out[f.name] = _encode(getattr(panel, f.name))
    if data is not None:
        key = f"p{index}"
        arrays[f"{key}_array"] = np.asarray(data.array)
        arrays[f"{key}_x"] = np.asarray(data.x, dtype=float)
        arrays[f"{key}_y"] = np.asarray(data.y, dtype=float)
        out["data"] = {"key": key, "x_label": data.x_label,
                       "y_label": data.y_label, "name": data.name,
                       "kind": data.kind, "info": _json_safe(data.info or {}),
                       "sources": _json_safe(getattr(data, "sources", []) or [])}
    return out


def save(path: str, figure, *, flatten=None, app_version: str = "") -> str:
    """Write ``figure`` (:class:`tools.output.figure.Figure`) to ``path``
    (``.afig`` added if missing)."""
    if not path.lower().endswith(EXTENSION):
        path += EXTENSION
    arrays = {}
    payload = {
        "format": FORMAT, "version": VERSION, "app_version": app_version,
        "title": figure.title, "rows": figure.rows, "cols": figure.cols,
        "style": dataclasses.asdict(figure.style),
        "shapes": [_encode(s) for s in figure.shapes],
        "panels": [_panel_dict(p, i, arrays, flatten)
                   for i, p in enumerate(figure.panels)],
    }
    buffer = io.BytesIO()
    np.savez_compressed(buffer, **arrays)
    folder = os.path.dirname(os.path.abspath(path))
    os.makedirs(folder, exist_ok=True)
    temporary = path + ".part"
    with zipfile.ZipFile(temporary, "w", zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("figure.json", json.dumps(payload, indent=1,
                                                   ensure_ascii=False))
        archive.writestr("arrays.npz", buffer.getvalue())
    os.replace(temporary, path)
    return path


# --------------------------------------------------------------------------
# Reading
# --------------------------------------------------------------------------
def _decode(value):
    if isinstance(value, dict) and "__type__" in value:
        cls = _TYPES.get(value["__type__"])
        if cls is None:
            return None
        known = {f.name for f in dataclasses.fields(cls)}
        return cls(**{k: _decode(v) for k, v in value.items() if k in known})
    if isinstance(value, list):
        return [_decode(v) for v in value]
    return value


def load(path: str):
    """The :class:`tools.output.figure.Figure` saved in ``path``."""
    with zipfile.ZipFile(path) as archive:
        payload = json.loads(archive.read("figure.json").decode("utf-8"))
        with np.load(io.BytesIO(archive.read("arrays.npz")),
                     allow_pickle=False) as npz:
            arrays = {key: npz[key] for key in npz.files}
    if payload.get("format") != FORMAT:
        raise ValueError(f"{os.path.basename(path)} is not a figure saved by "
                         f"this program")
    if int(payload.get("version", 0)) > VERSION:
        raise ValueError("this figure was saved by a newer version of the "
                         "program; update it to open the figure")
    style = nf.FigureStyle()
    for key, value in (payload.get("style") or {}).items():
        if hasattr(style, key):
            setattr(style, key, value)
    panel_fields = {f.name for f in dataclasses.fields(nf.Panel)}
    panels = []
    for record in payload.get("panels", []):
        values = {}
        for name, value in record.items():
            if name in panel_fields and name != "data":
                value = _decode(value)
                if name in _TUPLE_FIELDS and value is not None:
                    value = tuple(value)
                if name == "x_tick_labels":
                    value = [tuple(v) for v in value or []]
                values[name] = value
        panel = nf.Panel(**values)
        meta = record.get("data")
        if meta is not None:
            key = meta["key"]
            data = nf.PanelData(array=arrays[f"{key}_array"],
                                x=arrays[f"{key}_x"], y=arrays[f"{key}_y"],
                                x_label=meta.get("x_label", ""),
                                y_label=meta.get("y_label", ""),
                                name=meta.get("name", ""),
                                kind=meta.get("kind", "cut"),
                                info=dict(meta.get("info") or {}))
            data.sources = list(meta.get("sources") or [])
            panel.data = data
        else:
            panel.data = nf.empty_data()
        panels.append(panel)
    figure = nf.Figure(panels, rows=int(payload.get("rows", 1)),
                       cols=int(payload.get("cols", 1)), style=style,
                       title=str(payload.get("title", "")))
    figure.shapes = [s for s in (_decode(v) for v in payload.get("shapes", []))
                     if s is not None]
    #: where it was read from: its data folder still holds the originals
    #: if they have moved since
    figure.loaded_from = os.path.abspath(path)
    return figure


# --------------------------------------------------------------------------
# The data beside the figure
# --------------------------------------------------------------------------
def data_folder_of(path: str) -> str:
    """``data for fig1`` beside ``fig1.afig``."""
    stem = os.path.splitext(os.path.basename(path))[0]
    return os.path.join(os.path.dirname(os.path.abspath(path)), f"data for {stem}")


def _panel_folders(figure):
    """A folder name for every panel: its letter (``a``), an inset its
    parent's letter and ``inset``, an unlettered panel ``panel 3``."""
    names, used = [], set()
    for index, panel in enumerate(figure.panels):
        letter = re.sub(r"[^0-9A-Za-z]+", "", str(panel.label or ""))
        if panel.role == "inset":
            parent = next((p for p in figure.panels if p.uid == panel.zoom_of), None)
            base = (re.sub(r"[^0-9A-Za-z]+", "", str(parent.label or "")) + " inset"
                    if parent is not None and parent.label else "inset")
        else:
            base = letter or f"panel {index + 1}"
        name, n = base, 2
        while name in used:
            name, n = f"{base} {n}", n + 1
        used.add(name)
        names.append(name)
    return names


def _copy(source: str, folder: str) -> str:
    """Copy ``source`` into ``folder`` (kept as is if an identical copy is
    there already); a name taken by another file gets a number."""
    name = os.path.basename(source)
    target = os.path.join(folder, name)
    stem, ext = os.path.splitext(name)
    n = 2
    while os.path.exists(target):
        if (os.path.getsize(target) == os.path.getsize(source)
                and int(os.path.getmtime(target)) == int(os.path.getmtime(source))):
            return target
        target = os.path.join(folder, f"{stem} ({n}){ext}")
        n += 1
    shutil.copy2(source, target)
    return target


def _safe_name(text: str, default: str = "panel") -> str:
    name = re.sub(r'[\\/:*?"<>|]+', "_", str(text or "")).strip(" ._")
    return name[:80] or default


def _write_processed(panel, folder: str, picture_writer=None) -> list:
    """What the panel shows, as data: an image as a cut in the program's
    own .nxs (it opens in the viewer again), its drawn lines (calculated
    bands, contours, curves) as CSV, a picture as PNG."""
    from loader.nxs_file import save_dataset

    data = panel.data
    os.makedirs(folder, exist_ok=True)
    files = []
    name = _safe_name(data.name)
    if data.kind == "picture":
        name = _safe_name(os.path.splitext(name)[0]
                          if os.path.splitext(name)[1].lower() in
                          (".png", ".jpg", ".jpeg", ".tif", ".tiff", ".bmp",
                           ".gif", ".svg") else name)
        if picture_writer is not None:
            target = os.path.join(folder, f"{name}.png")
            picture_writer(np.asarray(data.array), target)
            files.append(target)
    elif data.kind not in ("stack", "empty", "scene3d"):
        values = np.asarray(data.array, dtype=float)
        if values.ndim == 2 and values.size > 4:
            target = os.path.join(folder, f"{name}.nxs")
            info = {k: v for k, v in _json_safe(data.info or {}).items()
                    if isinstance(v, (str, int, float, bool))}
            info["figure.panel"] = str(panel.label or "")
            save_dataset(target, [{
                "name": name, "kind": "cut",
                "axes": (np.asarray(data.x, float), np.asarray(data.y, float)),
                "labels": {"x": data.x_label or "x", "y": data.y_label or "y"},
                "value": values, "info": info}])
            files.append(target)
    if panel.overlays:
        target = os.path.join(folder, "lines.csv")
        with open(target, "w", encoding="utf-8") as handle:
            handle.write("# the lines drawn on the panel (calculated bands, "
                         "contours, curves): one block per line\n")
            handle.write(f"# x = {data.x_label or 'x'}, y = {data.y_label or 'y'}\n")
            for overlay in panel.overlays:
                handle.write(f"\n# {overlay.name or 'line'}\n")
                for x, y in zip(overlay.xs, overlay.ys):
                    handle.write(f"{float(x):.8g},{float(y):.8g}\n")
        files.append(target)
    return files


def write_data_folder(path: str, figure, *, picture_writer=None) -> dict:
    """Make ``data for <name>/`` beside the figure ``path``. Per panel
    (``a``, ``b``...): ``processed/`` -- the data the panel shows, as it
    shows it (:func:`_write_processed`) -- and ``original/`` -- the files it
    was made from, traced back (the first measurement, the calculation
    with the files it needs, the picture file) -- and ``sources.txt``
    saying which entry and what was done in between. ``picture_writer(rgba,
    target)`` writes a PNG (the caller's). Returns ``{folder: [files]}``."""
    root = data_folder_of(path)
    os.makedirs(root, exist_ok=True)
    written = {}
    for panel, name in zip(figure.panels, _panel_folders(figure)):
        data = panel.data
        if data is None or data.kind == "empty":
            continue
        folder = os.path.join(root, name)
        os.makedirs(folder, exist_ok=True)
        files = _write_processed(panel, os.path.join(folder, "processed"),
                                 picture_writer)
        lines = [f"Panel {panel.label or name}: {data.name}", "",
                 "processed/  what the panel shows: "
                 + (", ".join(os.path.basename(f) for f in files) or "(nothing)"),
                 "original/   what it was made from:", ""]
        sources = list(getattr(data, "sources", []) or [])
        if not sources:
            sources = sources_from_info(data.info)
        original = os.path.join(folder, "original")
        earlier = getattr(figure, "loaded_from", "")
        earlier = (os.path.join(data_folder_of(earlier), name, "original")
                   if earlier else "")

        def on_disk(file_path):
            """The file, or the copy the figure's earlier data folder kept."""
            if file_path and os.path.isfile(file_path):
                return file_path
            kept = os.path.join(earlier, os.path.basename(file_path or ""))
            if earlier and file_path and os.path.isfile(kept):
                return kept
            return ""

        for record in sources:
            origin = record.get("path", "")
            lines.append(f"{record.get('kind', 'dataset')}: {origin or '(no file)'}")
            if record.get("entry"):
                lines.append(f"  entry: {record['entry']}")
            for step in record.get("history", []):
                lines.append(f"  then: {step}")
            if record.get("note"):
                lines.append(f"  {record['note']}")
            found = on_disk(origin)
            if found and os.path.abspath(os.path.dirname(found)) == os.path.abspath(original):
                files.append(found)                  # saved over itself
                lines.append(f"  kept as: original/{os.path.basename(found)}")
            elif found:
                os.makedirs(original, exist_ok=True)
                copied = _copy(found, original)
                files.append(copied)
                lines.append(f"  copied as: original/{os.path.basename(copied)}")
                for extra in record.get("companions", []):
                    extra = on_disk(extra)
                    if extra:
                        files.append(_copy(extra, original))
                        lines.append(f"  with: original/{os.path.basename(extra)}")
            elif origin:
                lines.append("  (not on disk any more: not copied)")
        if not sources:
            lines.append("(no file is known: only the processed data is kept)")
        with open(os.path.join(folder, "sources.txt"), "w", encoding="utf-8") as handle:
            handle.write("\n".join(lines) + "\n")
        written[folder] = files
    return written
