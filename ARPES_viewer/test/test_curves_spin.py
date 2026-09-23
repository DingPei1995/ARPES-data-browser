"""Tests for tools/curves.py and tools/spin.py -- one-dimensional data and
spin polarisation -- and the pieces they rely on elsewhere."""
import numpy as np
import pytest

from tools import curves as C
from tools import spin as S


# ==========================================================================
# tools/curves.py
# ==========================================================================
def test_channel_names_and_sigma_pairing():
    names = C.channel_names({"curve.channels": "I↑|I↓|σ I↑|σ I↓"}, 4)
    assert names == ["I↑", "I↓", "σ I↑", "σ I↓"]
    assert C.data_channels(names) == [0, 1]
    assert C.sigma_of(names, 0) == 2 and C.sigma_of(names, 1) == 3
    assert C.channel_names({}, 1, "edc") == ["EDC"]
    assert C.channel_names({"curve.channels": "a|b"}, 3) == [
        "channel 0", "channel 1", "channel 2"]          # wrong length: numbered


def test_rebin_propagates_sigma():
    x = np.arange(10.0)
    values = np.column_stack([np.full(10, 4.0), np.full(10, 2.0)])
    names = ["N", "σ N"]
    xs, vs, _ = C.rebin(x, values, names, 3, "sum")
    assert xs.size == 3                         # the leftover point is dropped
    np.testing.assert_allclose(xs, [1, 4, 7])
    np.testing.assert_allclose(vs[:, 0], 12.0)
    np.testing.assert_allclose(vs[:, 1], np.sqrt(3) * 2.0)
    _, vm, _ = C.rebin(x, values, names, 3, "mean")
    np.testing.assert_allclose(vm[:, 0], 4.0)
    np.testing.assert_allclose(vm[:, 1], 2.0 / np.sqrt(3))
    with pytest.raises(ValueError):
        C.rebin(x, values, names, 20)


def test_normalise_together_keeps_ratios():
    x = np.linspace(0, 1, 11)
    up, down = 3 * np.ones(11), 1 * np.ones(11)
    values = np.column_stack([up, down, 0.3 * np.ones(11), 0.1 * np.ones(11)])
    names = ["u", "d", "σ u", "σ d"]
    _, v, _, factors = C.normalise(x, values, names, "max", together=True)
    assert v[0, 0] / v[0, 1] == pytest.approx(3.0)
    assert factors[0] == factors[1] == pytest.approx(4.0)
    np.testing.assert_allclose(v[:, 2], 0.3 / 4.0)
    _, v2, _, _ = C.normalise(x, values, names, "max", together=False)
    np.testing.assert_allclose(v2[:, 0], 1.0)
    np.testing.assert_allclose(v2[:, 1], 1.0)
    _, v3, _, _ = C.normalise(x, values, names, "region", region=(0.2, 0.4),
                              together=False)
    np.testing.assert_allclose(v3[:, 0], 1.0)
    with pytest.raises(ValueError):
        C.normalise(x, values, names, "region", region=(5, 6))


def test_background_subtraction():
    x = np.linspace(-1, 1, 201)
    peak = np.exp(-(x / 0.1) ** 2)
    line = 2.0 + 0.5 * x
    values = (peak + line)[:, None]
    _, v, _, _ = C.subtract_background(x, values, ["y"], "linear",
                                       region=(0.6, 1.0))
    # a line fitted where the peak is not recovers the line everywhere
    np.testing.assert_allclose(v[:, 0], peak, atol=1e-6)
    _, vc, _, _ = C.subtract_background(x, values, ["y"], "constant",
                                        region=(0.9, 1.0))
    # the mean over 0.9..1.0 is the line at 0.95, 0.025 below its end
    assert vc[-1, 0] == pytest.approx(0.025, abs=1e-9)
    with pytest.raises(ValueError):
        C.subtract_background(x, values, ["y"], "linear")


def test_crop_shift_and_poisson():
    x = np.linspace(0, 1, 11)
    values = np.arange(11.0)[:, None]
    xc, vc, _ = C.crop(x, values, ["n"], 0.25, 0.75)
    np.testing.assert_allclose(xc, [0.3, 0.4, 0.5, 0.6, 0.7])
    np.testing.assert_allclose(C.shift_x(x, 0.5)[0], -0.5)
    v, names, added = C.with_poisson_sigma(values, ["n"])
    assert names == ["n", "σ n"] and added == ["σ n"]
    np.testing.assert_allclose(v[:, 1], np.sqrt(np.arange(11)))
    # already has one, or is not counts: nothing added
    assert C.with_poisson_sigma(v, names)[2] == []
    assert C.with_poisson_sigma(values + 0.5, ["n"])[2] == []


def test_table_sorts_x():
    x, v, names = C.table([3, 1, 2], [("a", [30, 10, 20]), ("b", [3, 1, 2])])
    np.testing.assert_allclose(x, [1, 2, 3])
    np.testing.assert_allclose(v[:, 0], [10, 20, 30])


# ==========================================================================
# tools/spin.py
# ==========================================================================
LABELS = ["<0,0> +X  (GUI+Z)", "<90,180> +X  (GUI-Z)",
          "<0,0> -X  (GUI-Z)", "<90,180> -X  (GUI+Z)"]


def mean_counts(P, sherman, eps, n=(400.0, 320.0)):
    """The expected counts of the four channels above for polarisation P,
    Sherman function S, instrumental asymmetry eps, and transmissions n."""
    up, down = 1 + sherman * P, 1 - sherman * P
    plus_m, minus_m = 1 + eps, 1 - eps
    return np.array([n[0] * plus_m * up, n[1] * plus_m * down,
                     n[0] * minus_m * down, n[1] * minus_m * up])


def test_parse_channels():
    ch = S.parse_channels(LABELS)
    assert [c.axis for c in ch] == ["Z"] * 4
    assert [c.sign for c in ch] == [1, -1, -1, 1]
    assert [c.setting for c in ch] == ["0,0", "90,180", "0,0", "90,180"]
    assert [c.magnetisation for c in ch] == ["+X", "+X", "-X", "-X"]
    assert S.parse_channel(0, "C0 <0,0> +Y (GUI-X)").axis == "X"
    assert not S.parse_channel(0, "something else").usable
    assert S.axes_available(ch) == ["Z"]
    assert [p[0] for p in S.pairs(ch, "Z")] == ["0,0", "90,180"]
    report = S.design_report(ch, "Z")
    assert report == {"n_plus": 2, "n_minus": 2, "transmission": True,
                      "reflectivity": True}
    # one pair shares its setting on both sides, so the transmission
    # cancels; the target was magnetised differently, so the reflectivity not
    single = S.design_report([ch[0], ch[2]], "Z")
    assert single["transmission"] and not single["reflectivity"]


def test_cross_ratio_is_exact_on_expected_counts():
    ch = S.parse_channels(LABELS)
    P = np.linspace(-0.5, 0.5, 11)
    counts = np.array([mean_counts(p, 0.25, 0.07) for p in P])
    r = S.analyse(np.arange(11.0), counts, ch, "Z", sherman=0.25)
    np.testing.assert_allclose(r.polarisation, P, atol=1e-12)
    assert r.instrumental[0] == pytest.approx(0.07, abs=2e-3)
    # one pair alone is biased by the instrumental asymmetry
    single = S.analyse(np.arange(11.0), counts, ch, "Z", sherman=0.25,
                       method="pair:0,0")
    assert np.all(single.polarisation - P > 0.2)


def test_uncertainty_matches_monte_carlo():
    rng = np.random.default_rng(3)
    ch = S.parse_channels(LABELS)
    lam = mean_counts(0.3, 0.2, 0.04)
    draws = rng.poisson(lam, size=(4000, 4)).astype(float)
    r = S.analyse(np.arange(4000.0), draws, ch, "Z", sherman=0.2)
    assert np.mean(r.polarisation) == pytest.approx(0.3, abs=0.01)
    assert np.median(r.sigma_polarisation) == pytest.approx(
        np.std(r.polarisation), rel=0.06)
    # spin-resolved spectra: I↑ + I↓ = I, and their errors match the scatter
    np.testing.assert_allclose(r.up + r.down, r.intensity)
    assert np.median(r.sigma_up) == pytest.approx(np.std(r.up), rel=0.08)


def test_single_pair_sigma_is_textbook():
    ch = S.parse_channels(LABELS)
    counts = np.array([[600.0, 1.0, 400.0, 1.0]])
    r = S.analyse([0.0], counts, ch, "Z", sherman=1.0, method="pair:0,0")
    A, N = 0.2, 1000.0
    assert r.asymmetry[0] == pytest.approx(A)
    assert r.sigma_asymmetry[0] == pytest.approx(np.sqrt((1 - A ** 2) / N))


def test_zero_reference_and_binning():
    ch = S.parse_channels(LABELS)
    P = np.r_[np.zeros(10), np.full(10, 0.4)]
    counts = np.array([mean_counts(p, 0.2, 0.0) for p in P])
    counts[:, 0] *= 1.05                     # a drift the cross ratio keeps
    raw = S.analyse(np.arange(20.0), counts, ch, "Z", sherman=0.2)
    assert abs(raw.polarisation[:10].mean()) > 0.05
    fixed = S.analyse(np.arange(20.0), counts, ch, "Z", sherman=0.2,
                      zero_region=(0, 9))
    np.testing.assert_allclose(fixed.polarisation[:10], 0.0, atol=1e-9)
    assert fixed.zero_offset != 0.0
    binned = S.analyse(np.arange(20.0), counts, ch, "Z", sherman=0.2,
                       bin_factor=4)
    assert binned.x.size == 5


def test_refusals():
    ch = S.parse_channels(LABELS)
    counts = np.ones((5, 4)) * 100
    with pytest.raises(ValueError):
        S.analyse(np.arange(5.0), counts, ch, "Z", sherman=0.0)
    with pytest.raises(ValueError, match="negative"):
        S.analyse(np.arange(5.0), counts - 200, ch, "Z")
    with pytest.raises(ValueError):
        S.analyse(np.arange(5.0), counts, ch, "X")
    unbalanced = ch[:3]
    with pytest.raises(ValueError, match="same number"):
        S.analyse(np.arange(5.0), counts[:, :3], unbalanced, "Z")
    with pytest.raises(ValueError):
        S.analyse(np.arange(5.0), counts, ch, "Z", method="pair:nope")


def test_summary_mentions_what_matters():
    ch = S.parse_channels(LABELS)
    counts = np.array([mean_counts(0.1, 0.2, 0.02) for _ in range(5)])
    text = S.analyse(np.arange(5.0), counts, ch, "Z").summary()
    assert "cross ratio" in text and "S = 0.2" in text
    assert "transmission" in text and "reflectivity" in text
    assert "instrumental asymmetry" in text


# ==========================================================================
# Elsewhere
# ==========================================================================
def test_steepest_drop_finds_the_edge_beside_a_band():
    from tools.fermi import steepest_drop, fermi_edge_model
    e = np.linspace(-0.3, 0.15, 91)
    band = 0.8 * np.exp(-((e + 0.05) / 0.04) ** 2)
    i = fermi_edge_model(e, 0.0, 30.0, 0.02, 1.0, 0.0, 0.1, 0.0) * (1 + band)
    assert abs(steepest_drop(e, i)) < 0.02


def test_curve_kinds_in_the_tables():
    from loader.nxs_file import CURVE_KINDS, AXIS_SLOTS, KIND_LABELS, energy_slot
    for kind in CURVE_KINDS:
        assert AXIS_SLOTS[kind]["array"] == ("x", "y")
        assert kind in KIND_LABELS
    assert energy_slot("mdc") is None


def test_stack_panel_extent_is_the_data_range():
    import tools.figure as F
    panel = F.Panel(data=F.PanelData(array=np.full((2, 2), np.nan),
                                     x=np.array([17.0, 17.8]),
                                     y=np.array([0.0, 1.0]), kind="stack"))
    x0, x1, y0, y1 = panel.extent()
    assert (x0, x1) == (17.0, 17.8)
    assert y0 == pytest.approx(-0.04) and y1 == pytest.approx(1.04)
    image = F.Panel(data=F.PanelData(array=np.zeros((2, 2)),
                                     x=np.array([0.0, 1.0]),
                                     y=np.array([0.0, 1.0])))
    assert image.extent() == (-0.5, 1.5, -0.5, 1.5)   # pixels keep their halves
