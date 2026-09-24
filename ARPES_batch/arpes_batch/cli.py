"""
The command line -- for a person, and for an agent.

    python -m arpes_batch inventory  FOLDER... --out DIR   [--json]
    python -m arpes_batch references FOLDER... --out DIR   [--json]
    python -m arpes_batch centre     FILE [--entry E] [--energy E]
    python -m arpes_batch preview    FILE [--entry E] --png OUT.png
    python -m arpes_batch init       --inputs .. --references .. --output .. > recipe.json
    python -m arpes_batch run        RECIPE.json [--only GLOB] [--workers N] [--force] [--dry-run] [--json]

With ``--json`` the result goes to stdout as one JSON document and the
progress to stderr, so an agent can parse the one and ignore the other.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

from . import inventory as inv
from . import recipe as R
from ._paths import check_output_dir


def _log(message):
    print(message, file=sys.stderr, flush=True)


def _emit(obj, as_json: bool):
    if as_json:
        json.dump(obj, sys.stdout, indent=1, default=str, ensure_ascii=False)
        sys.stdout.write("\n")
    else:
        print(json.dumps(obj, indent=1, default=str, ensure_ascii=False))


def cmd_inventory(args):
    out = check_output_dir(args.out, args.roots)
    rows = inv.inventory(args.roots, args.pattern, progress=lambda i, n, f: _log(f"[{i + 1}/{n}] {f}"))
    paths = inv.write_manifest(rows, out)
    brief = [{k: r.get(k) for k in ("file", "entry", "kind", "role", "status",
                                    "photon_energy_eV", "pass_energy", "lens_mode",
                                    "temperature_K", "shape", "slit_axis_looks_like", "notes")}
             for r in rows]
    _emit({"manifest": paths, "entries": brief}, args.json)


def cmd_references(args):
    from .runner import prepare_references
    out = check_output_dir(args.out, args.roots)
    recipe = R.normalise({"inputs": args.roots, "references": args.roots, "output": out})
    refs = prepare_references(recipe, out, log=_log)
    _emit([{k: v for k, v in r.to_json().items() if k != "channels"} for r in refs], args.json)


def cmd_centre(args):
    from .center import suggest_centre
    from .dataset import Dataset
    ds = Dataset.load(args.file, args.entry)
    try:
        _emit(suggest_centre(ds, energy=args.energy, width=args.width), args.json)
    finally:
        ds.close()


def cmd_preview(args):
    from .dataset import Dataset
    from .preview import preview
    check_output_dir(os.path.dirname(os.path.abspath(args.png)), [args.file])
    ds = Dataset.load(args.file, args.entry)
    try:
        _emit({"png": preview(ds, args.png)}, args.json)
    finally:
        ds.close()


def cmd_init(args):
    recipe = {"name": args.name, "inputs": args.inputs, "references": args.references or [],
              "output": args.output, "kinds": ["map"], "steps": R.DEFAULT_STEPS,
              "preview": R.DEFAULTS["preview"], "overrides": {}, "workers": 1}
    R.normalise(recipe)                     # fail now, not at run time
    json.dump(recipe, sys.stdout, indent=2)
    sys.stdout.write("\n")


def cmd_run(args):
    from .runner import run
    recipe = R.load(args.recipe)
    summary = run(recipe, only=args.only, workers=args.workers, force=args.force,
                  dry_run=args.dry_run, log=_log)
    if not args.json and "results" in summary:
        for r in summary["results"]:
            skipped = [s["op"] for s in r["steps"] if s["status"] == "skipped"]
            _log(f"{r['status']:8s} {r['key']}" + (f"  skipped: {', '.join(skipped)}" if skipped else "")
                 + (f"  {r.get('error')}" if r.get("error") else ""))
        _log(f"report: {os.path.join(summary['run_dir'], 'report.md')}")
    _emit(summary if args.json else summary.get("counts", summary), args.json)
    return 1 if summary.get("counts", {}).get("failed") else 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="arpes_batch", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--json", action="store_true", help="machine-readable output on stdout")
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("inventory", parents=[common], help="list every dataset and the conditions it was taken at")
    p.add_argument("roots", nargs="+")
    p.add_argument("--out", required=True)
    p.add_argument("--pattern", nargs="+", default=["*.nxs"])
    p.set_defaults(func=cmd_inventory)

    p = sub.add_parser("references", parents=[common], help="fit the Fermi edge of gold references")
    p.add_argument("roots", nargs="+")
    p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_references)

    p = sub.add_parser("centre", parents=[common], help="suggest where Gamma is on a map")
    p.add_argument("file")
    p.add_argument("--entry")
    p.add_argument("--energy", type=float)
    p.add_argument("--width", type=float, default=0.05)
    p.set_defaults(func=cmd_centre)

    p = sub.add_parser("preview", parents=[common], help="a quick-look PNG of one dataset")
    p.add_argument("file")
    p.add_argument("--entry")
    p.add_argument("--png", required=True)
    p.set_defaults(func=cmd_preview)

    p = sub.add_parser("init", parents=[common], help="write a starting recipe to stdout")
    p.add_argument("--name", default="batch")
    p.add_argument("--inputs", nargs="+", required=True)
    p.add_argument("--references", nargs="*")
    p.add_argument("--output", required=True)
    p.set_defaults(func=cmd_init)

    p = sub.add_parser("run", parents=[common], help="run a recipe")
    p.add_argument("recipe")
    p.add_argument("--only", nargs="+", help="file or 'file entry' patterns")
    p.add_argument("--workers", type=int)
    p.add_argument("--force", action="store_true", help="ignore checkpoints")
    p.add_argument("--dry-run", action="store_true", help="plan only: what would run, with which reference")
    p.set_defaults(func=cmd_run)

    args = parser.parse_args(argv)
    return int(args.func(args) or 0)
