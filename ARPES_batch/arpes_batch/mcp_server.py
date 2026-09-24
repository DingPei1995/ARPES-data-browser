"""
The same operations as the command line, as MCP tools.

For an agent that talks MCP rather than running shell commands (Claude
Desktop, or any MCP client). Needs the ``mcp`` package (``pip install mcp``);
nothing else in this package does. Run it on the machine that holds the
data, e.g. in Claude Desktop's ``claude_desktop_config.json``::

    "mcpServers": {
      "arpes-batch": {
        "command": "python",
        "args": ["-m", "arpes_batch.mcp_server"],
        "cwd": "C:/path/to/ARPES-data-browser/ARPES_batch"
      }
    }

Every tool that writes goes through the same raw-data guard as the command
line. A whole beamtime can take longer than a client is willing to wait for
one call; ``plan`` first, then ``run_recipe`` with ``only`` on a few files at
a time.
"""
from __future__ import annotations

import json
import os

try:                                            # mcp >= 2
    from mcp.server.mcpserver import MCPServer
except ImportError:                             # mcp 1.x
    from mcp.server.fastmcp import FastMCP as MCPServer

from . import inventory as inv
from . import recipe as R
from ._paths import check_output_dir

server = MCPServer(
    "arpes-batch",
    instructions=(
        "Batch processing of ARPES maps with the ARPES viewer's own algorithms. "
        "Typical order: inventory -> fit_references -> plan -> run_recipe (a few "
        "files at a time) -> look at the preview PNGs -> put per-file Gamma / "
        "azimuth / energy window into the recipe's overrides -> run_recipe again "
        "(it resumes from checkpoints). Never point an output inside a raw-data folder."),
)


def _brief(row: dict) -> dict:
    keys = ("file", "entry", "kind", "role", "status", "photon_energy_eV", "pass_energy",
            "lens_mode", "temperature_K", "shape", "slit_axis_looks_like", "notes", "path")
    return {k: row.get(k) for k in keys if row.get(k) is not None}


@server.tool()
def inventory(roots: list[str], out_dir: str) -> dict:
    """List every dataset under ``roots`` (folders or .nxs files) with its kind,
    shape, photon energy, pass energy, lens mode and temperature. Writes
    manifest.json / manifest.csv to ``out_dir`` (which must be outside ``roots``)."""
    out = check_output_dir(out_dir, roots)
    rows = inv.inventory(roots)
    return {"manifest": inv.write_manifest(rows, out), "entries": [_brief(r) for r in rows]}


@server.tool()
def fit_references(roots: list[str], out_dir: str) -> list:
    """Fit the Fermi edge of every gold reference under ``roots``: E_F in
    kinetic energy, its curvature along the slit, the resolution, and a PNG
    of the fit in ``out_dir``/references."""
    from .runner import prepare_references
    out = check_output_dir(out_dir, roots)
    recipe = R.normalise({"inputs": roots, "references": roots, "output": out})
    return [{k: v for k, v in r.to_json().items() if k != "channels"}
            for r in prepare_references(recipe, out, log=lambda *_: None)]


@server.tool()
def plan(recipe_path: str, only: list[str] | None = None) -> dict:
    """What a recipe would do, without doing it: the jobs, and which
    reference each would be calibrated with (and why not, where none fits)."""
    from .runner import run
    return run(R.load(recipe_path), only=only, dry_run=True, log=lambda *_: None)


@server.tool()
def run_recipe(recipe_path: str, only: list[str] | None = None, force: bool = False) -> dict:
    """Run a recipe (optionally only the files matching ``only``). Returns,
    per entry, the status, what each step measured (E_F, grid power, Gamma
    suggestion and its score) and the paths of the .nxs results and PNGs."""
    from .runner import run
    summary = run(R.load(recipe_path), only=only, force=force, log=lambda *_: None)
    for result in summary.get("results", []):
        result.pop("traceback", None)
    return summary


@server.tool()
def suggest_centre(file: str, entry: str | None = None, energy: float | None = None) -> dict:
    """Where Gamma probably is on a map (deflector, slit angles), with a
    symmetry score 0..1. A suggestion: check it on the preview."""
    from .center import suggest_centre as suggest
    from .dataset import Dataset
    ds = Dataset.load(file, entry)
    try:
        return suggest(ds, energy=energy)
    finally:
        ds.close()


@server.tool()
def preview(file: str, png: str, entry: str | None = None) -> str:
    """Write a quick-look PNG of a raw or processed dataset (constant-energy
    maps and two cuts) and return its path, to be opened and looked at."""
    from .dataset import Dataset
    from .preview import preview as draw
    check_output_dir(os.path.dirname(os.path.abspath(png)), [file])
    ds = Dataset.load(file, entry)
    try:
        return draw(ds, png)
    finally:
        ds.close()


@server.tool()
def latest_results(out_dir: str) -> dict:
    """The summary of the last run written to ``out_dir``."""
    with open(os.path.join(out_dir, "latest_results.json"), encoding="utf-8") as fh:
        return json.load(fh)


if __name__ == "__main__":
    server.run()
