"""Which right-click entries of the main list apply to a selection
(ui/list_actions.py): by kind and by number of rows, one misfit row being
enough to grey an entry out, and a tooltip that says which row and why."""
import pytest

from ui.list_actions import availability, kind_of, RULES

CUT, MAP, KZ, KMAP, SPIN, EDC, SPEM, BAD = (
    ("c", "Cut"), ("m", "Map"), ("kz", "kz map"), ("k", "k-map"),
    ("s", "Spin EDC"), ("e", "EDC"), ("x", "SPEM"), ("u", "unknown"))


def ok(action, *rows):
    return availability(action, list(rows))[0]


def test_labels_map_to_kinds():
    assert kind_of("Map") == "map" and kind_of("kz map (k)") == "kz_map_k"
    assert kind_of("Spin EDC") == "spin_edc" and kind_of("nonsense") == "unknown"
    assert kind_of(None) == "unknown"


@pytest.mark.parametrize("rows", [(MAP,), (KZ,), (KMAP,), (MAP, KZ)])
def test_map_cuts_for_maps(rows):
    assert ok("slit_cut", *rows) and ok("deflector_cut", *rows)


@pytest.mark.parametrize("rows", [(CUT,), (SPIN,), (MAP, CUT), (SPEM,), ()])
def test_map_cuts_greyed_otherwise(rows):
    assert not ok("slit_cut", *rows)


def test_one_misfit_row_greys_the_entry_and_is_named():
    enabled, reason = availability("figure", [CUT, MAP, EDC])
    assert not enabled and "“m”" in reason and "Map" in reason


def test_counts():
    assert ok("stack", CUT) and not ok("stack", CUT, CUT)
    assert ok("arithmetic", CUT, CUT) and not ok("arithmetic", CUT)
    assert not ok("arithmetic", CUT, SPIN)
    assert ok("view3d", MAP) and not ok("view3d", CUT)
    assert not ok("open", CUT, MAP) and ok("open", SPEM)
    assert "select one row" in availability("rename", [CUT, CUT])[1]


def test_save_remove_log():
    assert ok("save", CUT, MAP, SPIN, SPEM) and not ok("save", CUT, BAD)
    assert ok("remove", CUT, BAD) and ok("log")
    assert not ok("remove")


def test_no_process_entry():
    assert "process" not in RULES
