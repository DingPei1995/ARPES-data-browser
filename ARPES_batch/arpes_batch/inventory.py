"""
What is in a folder of measurements, as a table.

Opening a map reads only its structure and axes -- the cube stays on disk
behind the loader's lazy array -- so listing a beamtime's worth of files
costs seconds, not the gigabytes they hold. The table is what an agent (or
a person) reads to decide what to process and with which reference.
"""
from __future__ import annotations

import csv
import fnmatch
import json
import os
import re
import warnings

import numpy as np

from loader import registry
from loader.nxs_file import axis_slots

#: Where each quantity is found in the metadata, tried in order. A string
#: starting with "~" is a pattern matched against every key. The ANTARES
#: names are the ones ``loader/nxs_file._read_common_info`` produces.
META_KEYS = {
    "photon_energy_eV": ("PhotonEnergy", "photon_energy_eV", "~^mono\\.energy$"),
    "pass_energy": ("MBS.passenergy", "pass_energy_eV", "~pass.?energy"),
    "lens_mode": ("MBS.lens_mode", "lens_mode", "~lens.?mode"),
    "center_ke_eV": ("MBS.center_ke", "center_ke"),
    "temperature_K": ("SampleTemperature_K", "~temperature_K$"),
    "polarisation": ("~polari[sz]ation",),
    "start_time": ("start_time",),
    "title": ("title",),
}

REFERENCE_PATTERNS = ("*au*", "*gold*", "*ref*")


def meta(info: dict, keys) -> object:
    for key in keys:
        if key.startswith("~"):
            pattern = re.compile(key[1:], re.IGNORECASE)
            for name in sorted(info):
                if pattern.search(str(name)):
                    return info[name]
        elif key in info:
            return info[key]
    return None


def _number(value):
    """``"PE100"`` -> 100.0, ``[21.2]`` -> 21.2, anything else -> None."""
    if value is None:
        return None
    try:
        return float(np.asarray(value, dtype=float).reshape(-1)[0])
    except (TypeError, ValueError):
        match = re.search(r"[-+]?\d+(\.\d+)?", str(value))
        return float(match.group(0)) if match else None


def _axis_summary(axis) -> dict:
    axis = np.asarray(axis, dtype=float).ravel()
    if axis.size == 0:
        return {"n": 0}
    step = float(np.mean(np.diff(axis))) if axis.size > 1 else 0.0
    return {"min": float(np.nanmin(axis)), "max": float(np.nanmax(axis)),
            "n": int(axis.size), "step": step}


def find_files(roots, patterns=("*.nxs",), recursive: bool = True) -> list:
    found = []
    for root in roots:
        if os.path.isfile(root):
            found.append(os.path.abspath(root))
            continue
        for folder, _dirs, files in os.walk(root):
            for name in sorted(files):
                if any(fnmatch.fnmatch(name.lower(), p.lower()) for p in patterns):
                    found.append(os.path.abspath(os.path.join(folder, name)))
            if not recursive:
                break
    return sorted(set(found))


def is_reference(path: str, patterns=REFERENCE_PATTERNS) -> bool:
    """A gold reference, by its file or folder name (``Au ref/Au_...nxs``)."""
    parts = [p.lower() for p in os.path.normpath(path).split(os.sep)[-2:]]
    return any(fnmatch.fnmatch(part, p) for part in parts for p in patterns)


def describe_entry(path: str, entry: str, listed: dict = None) -> dict:
    """One row: the entry's kind, axes and the conditions it was taken at."""
    row = {"file": os.path.basename(path), "path": path,
           "folder": os.path.basename(os.path.dirname(path)),
           "size_MB": round(os.path.getsize(path) / 1e6, 1),
           "entry": entry, "role": "reference" if is_reference(path) else "sample"}
    row.update({k: v for k, v in (listed or {}).items() if k in ("kind", "title", "start_time")})
    notes = []
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        try:
            scan = registry.load(path, entry)
        except Exception as exc:                          # noqa: BLE001
            row["status"] = "unreadable"
            row["error"] = f"{type(exc).__name__}: {exc}"
            return row
    notes += [str(w.message) for w in caught]
    try:
        info = dict(scan.info or {})
        row["kind"] = scan.kind
        row["loader"] = info.get("loader.name", "")
        value = scan.value4d if scan.kind == "spem_4d" else scan.value
        row["shape"] = list(getattr(value, "shape", ()))
        for slot in axis_slots(scan.kind, "constructor"):
            row[f"axis_{slot}"] = _axis_summary(getattr(scan, slot))
            row[f"label_{slot}"] = (scan.labels or {}).get(slot, "")
        for key, keys in META_KEYS.items():
            value = meta(info, keys)
            if value is not None:
                row[key] = value if isinstance(value, (int, float, str)) else str(value)
        row["pass_energy_eV"] = _number(row.get("pass_energy"))
        row["photon_energy_eV"] = _number(row.get("photon_energy_eV"))
        row["temperature_K"] = _number(row.get("temperature_K"))
        row["axis0_role"] = str(info.get("axis0.role", "angle"))
        if scan.kind in ("map", "cut"):
            slit_slot = "k" if scan.kind == "map" else "x"
            slit = row.get(f"axis_{slit_slot}", {})
            span = max(abs(slit.get("min", 0.0)), abs(slit.get("max", 0.0)))
            label = row.get(f"label_{slit_slot}", "")
            # The README says the MBS may deliver the slit axis in 1/A while
            # the loader labels it in degrees. Which one a file really has
            # decides whether a k conversion is right, so say what it looks like.
            row["slit_axis_looks_like"] = "momentum" if span < 2.0 else "angle"
            if row["slit_axis_looks_like"] == "momentum" and "°" in label:
                notes.append(f"slit axis is labelled {label!r} but spans only "
                             f"+-{span:.2f}: check whether it is already in 1/A")
        row["status"] = "ok"
    finally:
        scan.close()
    if notes:
        row["notes"] = notes
    return row


def inventory(roots, patterns=("*.nxs",), recursive: bool = True,
              progress=None) -> list:
    rows = []
    files = find_files(roots, patterns, recursive)
    for index, path in enumerate(files):
        if progress:
            progress(index, len(files), os.path.basename(path))
        loader = registry.detect(path)
        if loader is None:
            rows.append({"file": os.path.basename(path), "path": path,
                         "status": "no loader recognises it"})
            continue
        try:
            listed = list(loader.list_entries(path) or [])
        except Exception as exc:                          # noqa: BLE001
            rows.append({"file": os.path.basename(path), "path": path,
                         "status": "unreadable", "error": str(exc)})
            continue
        if not listed:
            rows.append({"file": os.path.basename(path), "path": path,
                         "status": "no recognised entries"})
            continue
        for item in listed:
            rows.append(describe_entry(path, item.get("entry"), item))
    return rows


FLAT_COLUMNS = ("file", "folder", "entry", "kind", "role", "status", "size_MB",
                "start_time", "photon_energy_eV", "pass_energy", "lens_mode",
                "temperature_K", "polarisation", "shape", "slit_axis_looks_like",
                "axis0_role", "path")


def write_manifest(rows: list, out_dir: str) -> dict:
    os.makedirs(out_dir, exist_ok=True)
    json_path = os.path.join(out_dir, "manifest.json")
    csv_path = os.path.join(out_dir, "manifest.csv")
    with open(json_path, "w", encoding="utf-8") as fh:
        json.dump(rows, fh, indent=1, ensure_ascii=False, default=str)
    axis_cols = sorted({k for r in rows for k in r if k.startswith("axis_") and k != "axis0_role"})
    with open(csv_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(list(FLAT_COLUMNS[:-1]) + axis_cols + ["notes", "path"])
        for r in rows:
            axes = []
            for col in axis_cols:
                a = r.get(col) or {}
                axes.append(f"{a.get('min', ''):.4g}..{a.get('max', ''):.4g} ({a.get('n')})"
                            if "min" in a else "")
            writer.writerow([r.get(c, "") for c in FLAT_COLUMNS[:-1]] + axes
                            + [" | ".join(r.get("notes", [])), r.get("path", "")])
    return {"json": json_path, "csv": csv_path}
