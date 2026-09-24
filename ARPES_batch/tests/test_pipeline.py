"""End-to-end checks on synthetic ANTARES maps whose answers are planted."""
import hashlib
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import arpes_batch                                            # noqa: E402
from arpes_batch import recipe as R                           # noqa: E402
from arpes_batch._paths import UnsafeOutput, check_output_dir  # noqa: E402
from arpes_batch.center import suggest_centre                 # noqa: E402
from arpes_batch.dataset import Dataset                       # noqa: E402
from arpes_batch.reference import fit_reference               # noqa: E402
from arpes_batch.runner import run                            # noqa: E402
from synthetic import write_map                               # noqa: E402

arpes_batch.ensure_viewer_on_path()

FRAME = dict(n_slit=160, n_e=220, grid=0.05, counts=60)


@pytest.fixture(scope="module")
def beamtime(tmp_path_factory):
    root = tmp_path_factory.mktemp("beamtime")
    raw = root / "raw"
    (raw / "TMD").mkdir(parents=True)
    (raw / "Au ref").mkdir()
    truth = {
        "WS_a.nxs": write_map(raw / "TMD" / "WS_a.nxs", seed=1, **FRAME),
        "WS_b.nxs": write_map(raw / "TMD" / "WS_b.nxs", seed=2, hv=62.0,
                              gamma_deg=(-1.0, 0.8), **FRAME),
        "CSS_c.nxs": write_map(raw / "TMD" / "CSS_c.nxs", seed=4, lens="L2", **FRAME),
        "Au.nxs": write_map(raw / "Au ref" / "Au.nxs", sample="gold", n_defl=5,
                            counts=400, seed=3, n_slit=160, n_e=220, grid=0.05),
    }
    return root, raw, truth


def _fingerprint(folder):
    out = {}
    for dirpath, _dirs, files in os.walk(folder):
        for name in files:
            path = os.path.join(dirpath, name)
            with open(path, "rb") as fh:
                out[path] = (os.path.getmtime(path), hashlib.sha1(fh.read()).hexdigest())
    return out


def _recipe(root, raw, **extra):
    recipe = {"inputs": [str(raw / "TMD")], "references": [str(raw / "Au ref")],
              "output": str(root / "processed"),
              "steps": [{"op": "degrid"},
                        {"op": "calibrate_energy", "required": True},
                        {"op": "crop", "energy": [-1.2, 0.15]},
                        {"op": "normalise"},
                        {"op": "save", "as": "angle"},
                        {"op": "kconvert", "theta_offset_deg": "auto",
                         "phi_offset_deg": "auto", "n_kx": 80, "n_ky": 80},
                        {"op": "save", "as": "k"}]}
    recipe.update(extra)
    return R.normalise(recipe)


# -- safety --------------------------------------------------------------------
def test_output_may_not_be_inside_or_around_the_raw_data(tmp_path):
    raw = tmp_path / "raw"
    raw.mkdir()
    with pytest.raises(UnsafeOutput):
        check_output_dir(str(raw), [str(raw)])
    with pytest.raises(UnsafeOutput):
        check_output_dir(str(raw / "processed"), [str(raw)])
    with pytest.raises(UnsafeOutput):
        check_output_dir(str(tmp_path), [str(raw)])
    link = tmp_path / "link"
    link.symlink_to(raw)
    with pytest.raises(UnsafeOutput):
        check_output_dir(str(link / "out"), [str(raw)])
    assert check_output_dir(str(tmp_path / "out"), [str(raw)])


# -- the pieces -------------------------------------------------------------------
def test_reference_recovers_the_planted_edge_and_curvature(beamtime):
    _root, raw, truth = beamtime
    ref = fit_reference(str(raw / "Au ref" / "Au.nxs"))
    t = truth["Au.nxs"]
    assert abs(ref.ef_center_eV - t["ef0"]) < 0.003
    span = np.ptp(np.polyval(ref.coeffs, t["slit"]))
    assert abs(span - t["curvature_eV"]) < 0.003
    assert ref.lens_mode == "L4" and ref.pass_energy == "PE50"
    assert abs(ref.work_function_eV - 4.35) < 0.003


def test_centre_suggestion_finds_gamma(beamtime):
    _root, raw, truth = beamtime
    ds = Dataset.load(str(raw / "TMD" / "WS_b.nxs"))
    try:
        found = suggest_centre(ds)
    finally:
        ds.close()
    assert found["theta_offset_deg"] == pytest.approx(-1.0, abs=0.15)
    assert found["phi_offset_deg"] == pytest.approx(0.8, abs=0.15)
    assert found["score"] > 0.8


# -- the whole run ---------------------------------------------------------------------
def test_run_processes_every_map_and_leaves_the_raw_data_alone(beamtime):
    root, raw, truth = beamtime
    before = _fingerprint(raw)
    summary = run(_recipe(root, raw), log=lambda *_: None)
    assert _fingerprint(raw) == before                  # not one byte, not one mtime
    by_key = {r["key"]: r for r in summary["results"]}

    a = by_key["WS_a__DeflX_0001"]
    assert a["status"] == "ok", a
    steps = {s["op"]: s for s in a["steps"]}
    assert steps["calibrate_energy"]["qc"]["ef_kinetic_eV"] == pytest.approx(55.65, abs=0.005)
    assert steps["degrid"]["qc"]["grid_power_after"] < steps["degrid"]["qc"]["grid_power_before"]

    b = by_key["WS_b__DeflX_0001"]
    ef_b = {s["op"]: s for s in b["steps"]}["calibrate_energy"]["qc"]["ef_kinetic_eV"]
    assert ef_b == pytest.approx(57.65, abs=0.005)       # moved with hv, 60 -> 62 eV

    # No reference for lens mode L2: calibration is required, so it stops
    # rather than writing an uncalibrated map that looks calibrated.
    c = by_key["CSS_c__DeflX_0001"]
    assert c["status"] == "stopped"
    assert [s["op"] for s in c["steps"]] == ["degrid", "calibrate_energy"]

    # What was written opens in the viewer as its own format, with history.
    k_path = [p for p in a["outputs"] if p.endswith("__k.nxs")][0]
    assert os.path.commonpath([k_path, str(root / "processed")]) == str(root / "processed")
    k = Dataset.load(k_path)
    try:
        assert k.kind == "k_map"
        assert "E_F" in k.label("z")
        from tools.process import history_of
        history = history_of(k)
        assert [h.split("(")[0] for h in history] == [
            "degrid", "fs_correction", "fermi_offset", "truncate", "self_normalise",
            "k_conversion"]
        z = k.axis("z")
        assert -1.25 < z.min() and z.max() < 0.2
        # The planted pocket is centred on Gamma: the k-map is symmetric.
        e0 = int(np.argmin(np.abs(z - 0.0)))
        from scipy.ndimage import gaussian_filter
        image = gaussian_filter(np.nan_to_num(k.array(float)[:, :, e0]), 1.5)
        ring = image > 0.5 * image.max()          # the pocket, not the background
        KX, KY = np.meshgrid(k.axis("x"), k.axis("k"), indexing="ij")
        assert abs(KX[ring].mean()) < 0.03 and abs(KY[ring].mean()) < 0.03
    finally:
        k.close()

    report = os.path.join(summary["run_dir"], "report.md")
    assert os.path.exists(report)


def test_a_rerun_resumes_from_the_last_valid_checkpoint(beamtime):
    root, raw, _truth = beamtime
    recipe = _recipe(root, raw)
    run(recipe, only=["WS_a*"], log=lambda *_: None)
    again = run(recipe, only=["WS_a*"], log=lambda *_: None)["results"][0]
    assert again["resumed_from"] == "k" and again["steps"] == []

    recipe["overrides"] = {"WS_a*": {"kconvert": {"theta_offset_deg": 0.0,
                                                  "phi_offset_deg": 0.0}}}
    changed = run(R.normalise(recipe), only=["WS_a*"], log=lambda *_: None)["results"][0]
    assert changed["resumed_from"] == "angle"
    assert [s["op"] for s in changed["steps"]] == ["kconvert", "save"]


def test_recipe_errors_are_caught_before_anything_runs(tmp_path):
    with pytest.raises(R.RecipeError):
        R.normalise({"inputs": ["x"], "output": "y", "steps": [{"op": "smoothify"}]})
    with pytest.raises(R.RecipeError):
        R.normalise({"inputs": ["x"], "output": "y", "overrides": {"*": {"kconvrt": {}}}})
    with pytest.raises(R.RecipeError):
        R.normalise({"inputs": ["x"]})
