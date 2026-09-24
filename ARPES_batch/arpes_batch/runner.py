"""
Running a recipe over a folder: one job per map, checkpoints, a report.

* **The raw data is only read.** The output folder is checked against every
  input and reference folder before anything is written
  (:func:`._paths.check_output_dir`), files are opened read-only by the
  loader, and results are written under the output folder only.
* **One failure is one row.** A job that raises is recorded with its
  traceback and the run goes on; a step that does not apply (no grid, no
  matching reference) is recorded as skipped and the chain continues.
* **Re-running is cheap.** Every ``save`` in the chain writes a checkpoint
  with a hash of everything that produced it (the raw file's size and
  modification time, the entry, the steps up to there). A re-run starts
  from the last checkpoint whose hash still matches, so changing only the
  k conversion does not repeat a de-grid.
* **Separate processes, not threads.** HDF5 here is not thread-safe, and a
  map is hundreds of MB; ``workers`` processes each take whole files.
"""
from __future__ import annotations

import datetime
import hashlib
import json
import os
import shutil
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

from . import inventory as inv
from . import recipe as R
from ._paths import check_output_dir
from .dataset import Dataset
from .preview import plot_reference, preview
from .reference import Reference, choose_reference, fit_reference
from .steps import STEPS, Skip


def _job_key(row: dict) -> str:
    stem = os.path.splitext(row["file"])[0]
    entry = str(row.get("entry") or "").strip("/")
    return f"{stem}__{entry}" if entry else stem


def _chain_hash(row: dict, steps: list, upto: int, reference_path: str) -> str:
    stat = os.stat(row["path"])
    payload = {"path": os.path.basename(row["path"]), "size": stat.st_size,
               "mtime": int(stat.st_mtime), "entry": row.get("entry"),
               "reference": os.path.basename(reference_path or ""),
               "steps": steps[:upto + 1]}
    return hashlib.sha1(json.dumps(payload, sort_keys=True, default=str).encode()).hexdigest()[:16]


def _json_safe(value):
    if isinstance(value, dict):
        return {str(k): _json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(v) for v in value]
    if isinstance(value, (np.floating, float)):
        return None if not np.isfinite(value) else float(value)
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.ndarray):
        return _json_safe(value.tolist())
    return value


# -- references ---------------------------------------------------------------
def prepare_references(recipe: dict, out_dir: str, log=print) -> list:
    """Fit every reference once (cached in ``<output>/references``)."""
    if not recipe["references"]:
        return []
    ref_dir = os.path.join(out_dir, "references")
    os.makedirs(ref_dir, exist_ok=True)
    rows = inv.inventory(recipe["references"], recipe["patterns"])
    fit_options = dict(recipe.get("reference_fit") or {})
    refs = []
    for row in rows:
        if row.get("status") != "ok" or row.get("kind") not in ("map", "cut"):
            continue
        key = _job_key(row)
        stat = os.stat(row["path"])
        tag = hashlib.sha1(json.dumps([stat.st_size, int(stat.st_mtime), fit_options],
                                      sort_keys=True).encode()).hexdigest()[:10]
        cache = os.path.join(ref_dir, f"{key}.json")
        if os.path.exists(cache):
            with open(cache, encoding="utf-8") as fh:
                stored = json.load(fh)
            if stored.get("_tag") == tag:
                refs.append(Reference.from_json(stored["reference"]))
                continue
        try:
            ref = fit_reference(row["path"], row.get("entry"), **fit_options)
        except Exception as exc:                          # noqa: BLE001
            log(f"reference {key}: fit failed ({exc})")
            continue
        with open(cache, "w", encoding="utf-8") as fh:
            json.dump({"_tag": tag, "reference": _json_safe(ref.to_json())}, fh, indent=1)
        plot_reference(ref, os.path.join(ref_dir, f"{key}.png"))
        log(f"reference {key}: E_F(mid) {ref.ef_center_eV:.4f} eV at hv "
            f"{ref.photon_energy_eV}, curvature {1000 * np.ptp(np.polyval(ref.coeffs, np.linspace(*ref.slit_range, 50))):.1f} meV, "
            f"rms {ref.residual_meV:.1f} meV, {ref.n_ok}/{ref.n_channels} channels")
        refs.append(ref)
    return refs


# -- one job ---------------------------------------------------------------------
def run_job(job: dict) -> dict:
    """Process one entry. Runs in a worker process; everything it needs is
    in ``job`` and everything it reports is in the returned dict."""
    row, steps, out_dir = job["row"], job["steps"], job["job_dir"]
    reference = Reference.from_json(job["reference"]) if job.get("reference") else None
    ctx = {"row": row, "reference": reference, "reference_notes": job.get("reference_notes")}
    key = _job_key(row)
    result = {"key": key, "file": row["file"], "entry": row.get("entry"),
              "status": "ok", "steps": [], "outputs": [],
              "reference": reference.path if reference else None,
              "reference_notes": job.get("reference_notes") or []}
    os.makedirs(out_dir, exist_ok=True)
    started = time.time()
    source = None
    try:
        # Resume from the last checkpoint that is still valid.
        start_at, ds = 0, None
        if not job.get("force"):
            for index in range(len(steps) - 1, -1, -1):
                step = steps[index]
                if step.get("op") != "save":
                    continue
                path = os.path.join(out_dir, f"{key}__{step['as']}.nxs")
                meta = path + ".json"
                if os.path.exists(path) and os.path.exists(meta):
                    with open(meta, encoding="utf-8") as fh:
                        stored = json.load(fh)
                    if stored.get("chain") == _chain_hash(row, steps, index, result["reference"]):
                        ds = Dataset.load(path)
                        ds.name = stored.get("name", ds.name)
                        ds.info = dict(ds.info)
                        ds.value = ds.array()            # read it now...
                        ds.close()                       # ...and let the file go
                        start_at = index + 1
                        result["resumed_from"] = step["as"]
                        break
        if ds is None:
            source = ds = Dataset.load(row["path"], row.get("entry"), name=key)

        for index in range(start_at, len(steps)):
            step = dict(steps[index])
            op = step.pop("op")
            required = bool(step.pop("required", False))
            t0 = time.time()
            record = {"op": op, "status": "ok"}
            if op == "save":
                name = step["as"]
                path = os.path.join(out_dir, f"{key}__{name}.nxs")
                ds.save(path)
                with open(path + ".json", "w", encoding="utf-8") as fh:
                    json.dump({"chain": _chain_hash(row, steps, index, result["reference"]),
                               "name": ds.name, "kind": ds.kind, "shape": list(ds.shape),
                               "steps": steps[:index + 1]}, fh, indent=1, default=str)
                png = preview(ds, os.path.join(out_dir, f"{key}__{name}.png"),
                              **job.get("preview", {}),
                              centre=job.get("_centre") if ds.kind == "map" else
                              ((0.0, 0.0) if ds.kind == "k_map" else None))
                result["outputs"] += [path, png]
                record["path"] = path
            else:
                try:
                    ds_next, qc = STEPS[op](ds, ctx, **step)
                    record["qc"] = _json_safe(qc)
                    if op == "kconvert" and "centre_suggestion" in qc:
                        job["_centre"] = (qc["centre_suggestion"]["theta_offset_deg"],
                                          qc["centre_suggestion"]["phi_offset_deg"])
                    ds = ds_next
                except Skip as why:
                    record.update(status="skipped", note=str(why))
            record["seconds"] = round(time.time() - t0, 2)
            result["steps"].append(record)
            if required and record["status"] == "skipped":
                # Everything after a required step assumes it ran (an energy
                # axis relative to E_F, say); going on would write results
                # that look right and are not.
                result["status"] = "stopped"
                result["error"] = f"required step {op} was skipped: {record['note']}"
                break
        # A map shown before its k conversion carries the Gamma suggestion
        # only once kconvert has run; redraw the angle preview with it.
        if job.get("_centre"):
            for index, step in enumerate(steps):
                if step.get("op") == "save":
                    path = os.path.join(out_dir, f"{key}__{step['as']}.nxs")
                    if os.path.exists(path):
                        saved = Dataset.load(path)
                        try:
                            if saved.kind == "map":
                                preview(saved, path[:-4] + ".png", **job.get("preview", {}),
                                        centre=job["_centre"])
                        finally:
                            saved.close()
        if result["status"] == "ok" and any(s["status"] == "skipped" for s in result["steps"]):
            result["status"] = "partial"
    except Exception as exc:                              # noqa: BLE001
        result["status"] = "failed"
        result["error"] = f"{type(exc).__name__}: {exc}"
        result["traceback"] = traceback.format_exc()
    finally:
        if source is not None:
            source.close()
    result["seconds"] = round(time.time() - started, 1)
    with open(os.path.join(out_dir, f"{key}__result.json"), "w", encoding="utf-8") as fh:
        json.dump(_json_safe(result), fh, indent=1)
    return _json_safe(result)


# -- the whole run -------------------------------------------------------------------
def plan(recipe: dict, only=None, log=print) -> tuple:
    out_dir = check_output_dir(recipe["output"], recipe["inputs"] + recipe["references"])
    os.makedirs(out_dir, exist_ok=True)
    rows = inv.inventory(recipe["inputs"], recipe["patterns"])
    inv.write_manifest(rows, out_dir)
    refs = prepare_references(recipe, out_dir, log=log)
    by_path = {os.path.abspath(r.path): r for r in refs}
    jobs = []
    for row in rows:
        if row.get("status") != "ok" or row.get("role") == "reference":
            continue
        if row.get("kind") not in recipe["kinds"] or not R.select(row, recipe, only):
            continue
        steps, override = R.steps_for(row, recipe)
        if override.get("skip"):
            continue
        if override.get("reference"):
            chosen = by_path.get(os.path.abspath(override["reference"]))
            notes = [] if chosen else [f"override reference {override['reference']} was not fitted"]
        else:
            chosen, notes = choose_reference(row, refs)
        jobs.append({"row": row, "steps": steps,
                     "job_dir": os.path.join(out_dir, row.get("folder") or "data", _job_key(row)),
                     "reference": _json_safe(chosen.to_json()) if chosen else None,
                     "reference_notes": notes, "preview": recipe.get("preview", {})})
    return out_dir, rows, refs, jobs


def run(recipe: dict, *, only=None, workers: int = None, force: bool = False,
        dry_run: bool = False, log=print) -> dict:
    out_dir, rows, refs, jobs = plan(recipe, only, log=log)
    stamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    run_dir = os.path.join(out_dir, "runs", stamp)
    summary = {"output": out_dir, "run_dir": run_dir, "entries_found": len(rows),
               "references": [os.path.basename(r.path) for r in refs],
               "jobs": [{"key": _job_key(j["row"]), "reference": j["reference"] and
                         os.path.basename(j["reference"]["path"]),
                         "reference_notes": j["reference_notes"]} for j in jobs]}
    if dry_run:
        summary.pop("run_dir")
        return summary
    os.makedirs(run_dir, exist_ok=True)
    if recipe.get("_source"):
        shutil.copy2(recipe["_source"], os.path.join(run_dir, "recipe.json"))
    for job in jobs:
        job["force"] = force
    workers = int(workers or recipe.get("workers") or 1)
    results = []
    if workers <= 1 or len(jobs) <= 1:
        for index, job in enumerate(jobs):
            log(f"[{index + 1}/{len(jobs)}] {_job_key(job['row'])}")
            results.append(run_job(job))
            log(f"    -> {results[-1]['status']} ({results[-1]['seconds']} s)")
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = {pool.submit(run_job, job): job for job in jobs}
            for done, future in enumerate(as_completed(futures), 1):
                results.append(future.result())
                log(f"[{done}/{len(jobs)}] {results[-1]['key']} -> {results[-1]['status']}")
    results.sort(key=lambda r: r["key"])
    summary["results"] = results
    summary["counts"] = {s: sum(r["status"] == s for r in results)
                         for s in ("ok", "partial", "stopped", "failed")}
    with open(os.path.join(run_dir, "results.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=1)
    write_report(summary, os.path.join(run_dir, "report.md"))
    shutil.copy2(os.path.join(run_dir, "results.json"), os.path.join(out_dir, "latest_results.json"))
    return summary


def write_report(summary: dict, path: str) -> str:
    base = os.path.dirname(path)
    lines = [f"# Batch run {os.path.basename(base)}", "",
             f"Output: `{summary['output']}`  ", f"References: {', '.join(summary['references']) or 'none'}",
             "", "| entry | status | reference | E_F kin (eV) | grid power | Γ suggestion (score) | time (s) |",
             "|---|---|---|---|---|---|---|"]
    for r in summary.get("results", []):
        qc = {s["op"]: s.get("qc", {}) for s in r["steps"]}
        ef = qc.get("calibrate_energy", {}).get("ef_kinetic_eV")
        grid = qc.get("degrid", {})
        centre = qc.get("kconvert", {}).get("centre_suggestion") or {}
        lines.append("| {} | {} | {} | {} | {} | {} | {} |".format(
            r["key"], r["status"] + (f": {r.get('error')}" if r.get("error") else ""),
            os.path.basename(r.get("reference") or "-"),
            f"{ef:.4f}" if ef is not None else "-",
            f"{grid.get('grid_power_before', float('nan')):.1f} → {grid.get('grid_power_after', float('nan')):.2f}" if grid else "-",
            f"({centre['theta_offset_deg']:.2f}, {centre['phi_offset_deg']:.2f}) ({centre['score']:.2f})" if centre else "-",
            r["seconds"]))
    lines.append("")
    for r in summary.get("results", []):
        skipped = [f"{s['op']}: {s.get('note')}" for s in r["steps"] if s["status"] == "skipped"]
        pngs = [p for p in r["outputs"] if p.endswith(".png")]
        lines.append(f"## {r['key']}")
        lines += [f"- skipped {s}" for s in skipped]
        lines += [f"- note: {n}" for n in r.get("reference_notes") or []]
        lines += [f"![{os.path.basename(p)}]({os.path.relpath(p, base)})" for p in pngs]
        lines.append("")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))
    return path
