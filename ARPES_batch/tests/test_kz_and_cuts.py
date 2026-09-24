"""kz maps, cuts, the sample questions and the match to a band calculation,
on synthetic data whose answers are planted."""
import io
import json
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import arpes_batch                                                  # noqa: E402
from arpes_batch import calcbands, cli                              # noqa: E402
from arpes_batch import recipe as R                                 # noqa: E402
from arpes_batch.dataset import Dataset                             # noqa: E402
from arpes_batch.runner import run                                  # noqa: E402
from arpes_batch.sample import (NeedsInput, Sample, apply_answers,  # noqa: E402
                                prompt, surface_period)
from synthetic import (BCT, HEX, hex_bands, write_band_file,        # noqa: E402
                       write_hv_scan_antares, write_kz_native, write_map)

arpes_batch.ensure_viewer_on_path()
from tools.lattice import LatticeParams                             # noqa: E402

FRAME = dict(n_slit=160, n_e=220, grid=0.05, counts=60)
SAMPLE = {"name": "bct metal", "lattice": BCT, "surface_normal": [0, 0, 1],
          "inner_potential_eV": "scan", "calculation": None}


@pytest.fixture(scope="module")
def beamtime(tmp_path_factory):
    root = tmp_path_factory.mktemp("kz_beamtime")
    raw = root / "raw"
    for sub in ("maps", "cuts", "kz", "hv", "Au ref"):
        (raw / sub).mkdir(parents=True)
    truth = {
        "map": write_map(raw / "maps" / "map_a.nxs", seed=1, **FRAME),
        "cut": write_map(raw / "cuts" / "cut_a.nxs", seed=5, n_defl=1,
                         defl_range=(2.0, 2.0), **FRAME),
        "kz": write_kz_native(raw / "kz" / "kz_bct.nxs", seed=3),
        "hv": write_hv_scan_antares(raw / "hv" / "hv_scan.nxs"),
    }
    write_map(raw / "Au ref" / "Au.nxs", sample="gold", n_defl=5, counts=400, seed=3,
              n_slit=160, n_e=220, grid=0.05, curvature_eV=0.02)
    return root, raw, truth


def _recipe(root, raw, sample=SAMPLE, **extra):
    recipe = {
        "inputs": [str(raw / d) for d in ("maps", "cuts", "kz", "hv")],
        "references": [str(raw / "Au ref")], "output": str(root / "processed"),
        "photon_energy_scans": ["hv_scan*"],
        "steps": {
            "map": [{"op": "degrid"}, {"op": "calibrate_energy", "required": True},
                    {"op": "save", "as": "angle"}],
            "kz_map": [{"op": "kz_calibrate", "source": "self", "required": True},
                       {"op": "crop", "energy": [-1.2, 0.1]},
                       {"op": "save", "as": "hv"},
                       {"op": "kz_match_calculation"},
                       {"op": "kz_convert", "angle_offset": "auto", "n_kz": 160, "n_kpar": 120},
                       {"op": "save", "as": "kz"}],
            "cut": [{"op": "degrid"}, {"op": "calibrate_energy", "required": True},
                    {"op": "crop", "energy": [-1.2, 0.15]},
                    {"op": "kconvert_cut", "gamma_slit_deg": "auto", "gamma_deflector_deg": 2.0},
                    {"op": "save", "as": "k"}, {"op": "curvature"},
                    {"op": "save", "as": "curvature"}]},
        "overrides": {"hv_scan*": {"kz_calibrate": {"source": "reference"},
                                   "kz_convert": {"inner_potential": 12.0}}},
    }
    if sample is not None:
        recipe["sample"] = sample
    recipe.update(extra)
    return R.normalise(recipe)


# -- the sample, and asking for it ------------------------------------------------------
def test_the_kz_period_follows_the_centering():
    bct = surface_period(LatticeParams(**BCT), (0, 0, 1))
    assert bct["first_reflection"] == [0, 0, 2]
    assert bct["period_invA"] == pytest.approx(4 * np.pi / BCT["c"], rel=1e-6)
    assert "centred" in bct["note"]
    hexagonal = surface_period(LatticeParams(**HEX), (0, 0, 1))
    assert hexagonal["period_invA"] == pytest.approx(2 * np.pi / HEX["c"], rel=1e-6)


def test_a_kz_recipe_without_the_sample_asks_before_running(beamtime, capsys):
    root, raw, _truth = beamtime
    recipe = _recipe(root, raw, sample=None)
    plan = run(recipe, dry_run=True, log=lambda *_: None)
    asked = [m["field"] for m in plan["needs_input"]]
    assert asked == ["sample.lattice", "sample.surface_normal",
                     "sample.inner_potential_eV", "sample.calculation"]
    assert [m.get("optional", False) for m in plan["needs_input"]] == [False] * 3 + [True]
    with pytest.raises(NeedsInput):
        run(recipe, log=lambda *_: None)
    assert not os.path.exists(os.path.join(recipe["output"], "kz"))   # nothing ran

    # the same, for an agent: JSON on stdout, exit status 2
    path = root / "recipe.json"
    path.write_text(json.dumps({k: v for k, v in recipe.items() if not k.startswith("_")}))
    capsys.readouterr()
    assert cli.main(["run", str(path), "--json"]) == 2
    out = json.loads(capsys.readouterr().out)
    assert out["needs_input"][0]["field"] == "sample.lattice"
    assert "do not guess" in out["message"]


def test_the_terminal_prompt_fills_the_sample_block():
    missing = [{"field": "sample.lattice", "question": "?", "example": {}, "why": "."},
               {"field": "sample.surface_normal", "question": "?", "example": [], "why": "."},
               {"field": "sample.inner_potential_eV", "question": "?", "example": 1, "why": "."},
               {"field": "sample.calculation", "question": "?", "example": {}, "why": ".",
                "optional": True}]
    answers = prompt(missing, stream_in=io.StringIO(
        json.dumps(BCT) + "\n[0, 0, 1]\nscan\n\n"), stream_out=io.StringIO())
    block = apply_answers({}, answers)
    assert block["calculation"] is None                      # asked, and there is none
    sample = Sample.from_recipe(block)
    assert sample.surface()["period_invA"] == pytest.approx(4 * np.pi / BCT["c"])


def test_malformed_sample_values_are_refused():
    base = {"inputs": ["x"], "output": "y"}
    with pytest.raises(R.RecipeError):
        R.normalise(dict(base, sample={"surface_normal": "001"}))
    with pytest.raises(R.RecipeError):
        R.normalise(dict(base, sample={"inner_potential_eV": "about 12"}))


# -- the whole beamtime ----------------------------------------------------------------
@pytest.fixture(scope="module")
def processed(beamtime):
    root, raw, truth = beamtime
    summary = run(_recipe(root, raw), log=lambda *_: None)
    return {r["key"]: r for r in summary["results"]}, summary, truth


def _steps(result):
    return {s["op"]: s for s in result["steps"]}


def test_every_kind_is_processed(processed):
    results, summary, _truth = processed
    assert {r["kind"] for r in results.values()} == {"map", "cut", "kz_map"}
    assert all(r["status"] in ("ok", "partial") for r in results.values()), \
        {k: (r["status"], r.get("error")) for k, r in results.items()}
    assert summary["sample"]["surface"]["first_reflection"] == [0, 0, 2]


def test_a_kz_map_is_aligned_and_converted(processed):
    results, _summary, truth = processed
    r = results["kz_bct__synthetic kz scan"]
    steps = _steps(r)
    ef = np.array(steps["kz_calibrate"]["qc"]["ef_per_hv"])[:, 1]
    drift = truth["kz"]["drift"]
    assert np.sqrt(np.mean(((ef - ef.mean()) - (drift - drift.mean())) ** 2)) < 0.015
    assert steps["kz_match_calculation"]["status"] == "skipped"          # none given
    scan = steps["kz_convert"]["qc"]["v0_scan"]
    assert scan["best"] == pytest.approx(truth["kz"]["v0"], abs=1.5)
    kz = Dataset.load([p for p in r["outputs"] if p.endswith("__kz.nxs")][0])
    try:
        assert kz.kind == "kz_map_k" and kz.label("x").startswith("k_z")
        from tools.process import history_of
        names = [h.split("(")[0] for h in history_of(kz)]
        assert names[-1] == "kz_to_k" and "kz_align" in names
    finally:
        kz.close()


def test_an_antares_hv_scan_is_calibrated_from_gold(processed):
    results, _summary, truth = processed
    r = results["hv_scan__hvscan_0001"]
    assert r["kind"] == "kz_map"
    rows = np.array(_steps(r)["kz_calibrate"]["qc"]["ef_per_hv"])
    # E_F moves with the photon energy: hv - W, W from the gold
    assert np.allclose(rows[:, 1] - rows[:, 0], rows[0, 1] - rows[0, 0], atol=1e-9)
    assert rows[0, 0] - rows[0, 1] == pytest.approx(truth["hv"]["phi"], abs=0.005)
    assert _steps(r)["kz_convert"]["qc"]["zones_covered"] < 1.5


def test_a_cut_uses_the_maps_grid_and_finds_its_centre(processed):
    results, _summary, _truth = processed
    steps = _steps(results["cut_a__DeflX_0001"])
    assert steps["degrid"]["qc"]["grid_from"].startswith("map_a")
    centre = steps["kconvert_cut"]["qc"]["centre_suggestion"]
    assert centre["gamma_slit_deg"] == pytest.approx(-1.5, abs=0.2)
    assert steps["curvature"]["status"] == "ok"


# -- the match to a calculation ----------------------------------------------------------
def test_band_file_vertices_are_named_from_the_lattice(tmp_path):
    write_band_file(tmp_path / "band.dat")
    bands = calcbands.read_bands(str(tmp_path / "band.dat"), lattice=LatticeParams(**HEX))
    assert [n for _k, n in bands.vertices] == ["G", "M", "K", "G", "A", "L", "H", "A"]
    with pytest.raises(ValueError):
        calcbands.read_bands(str(tmp_path / "band.dat"), labels=["G", "M", "K"])


def test_matching_a_calculation_names_the_high_symmetry_photon_energies(tmp_path):
    raw = tmp_path / "raw"
    raw.mkdir()
    write_band_file(tmp_path / "band.dat")
    measured = lambda k_par, k_z: [e + 0.05 for e in hex_bands(k_par, k_z)]  # noqa: E731
    write_kz_native(raw / "kz_hex.nxs", hv=np.arange(20.0, 120.1, 1.0), v0=11.0,
                    bands=measured, n_e=200, energy_range=(-1.9, 0.25),
                    curvature_eV=0.0, seed=7)
    recipe = R.normalise({
        "inputs": [str(raw)], "output": str(tmp_path / "out"),
        "sample": {"lattice": HEX, "surface_normal": [0, 0, 1],
                   "inner_potential_eV": "calculation",
                   "calculation": {"file": str(tmp_path / "band.dat")}},
        "steps": {"kz_map": [{"op": "kz_calibrate", "fs_correction": False},
                             {"op": "kz_match_calculation"},
                             {"op": "kz_convert", "n_kz": 120, "n_kpar": 100},
                             {"op": "save", "as": "kz"}]}})
    result = run(recipe, log=lambda *_: None)["results"][0]
    match = _steps(result)["kz_match_calculation"]["qc"]
    assert match["inner_potential_eV"] == pytest.approx(11.0, abs=0.5)
    assert match["energy_shift_eV"] == pytest.approx(0.05, abs=0.03)
    truth = calcbands.planes(np.arange(20.0, 120.1), 2 * np.pi / HEX["c"], v0=11.0,
                             work_function=4.4, edge_label="A")
    assert [(p["plane"], p["order"]) for p in match["planes"]] == \
        [(p["plane"], p["order"]) for p in truth]
    assert np.allclose([p["hv_eV"] for p in match["planes"]],
                       [p["hv_eV"] for p in truth], atol=1.5)
    assert os.path.exists(match["figure"])
    assert _steps(result)["kz_convert"]["qc"]["v0_from_calculation"] == match["inner_potential_eV"]
