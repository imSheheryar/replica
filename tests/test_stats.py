import numpy as np

from replica.stats import bh, by, combine_pvalues, maxt_adjusted, npc_test, permutation_pvalue


def test_combine_pvalues_orders_evidence():
    weak = combine_pvalues(np.array([0.4, 0.5, 0.6]), "liptak")
    strong = combine_pvalues(np.array([0.001, 0.002, 0.003]), "liptak")
    assert strong > weak


def test_bh_known_values():
    p = np.array([0.01, 0.04, 0.03, 0.20])
    adjusted = bh(p)
    expected = np.array([0.04, 0.04 * 4 / 3, 0.04 * 4 / 3, 0.20])
    assert np.allclose(adjusted, expected)


def test_by_is_at_least_as_conservative_as_bh():
    p = np.array([0.001, 0.01, 0.02, 0.5])
    assert np.all(by(p) >= bh(p) - 1e-12)


def test_permutation_pvalue_plus_one():
    null = np.arange(100, dtype=float)
    assert permutation_pvalue(1000, null) == 1 / 101
    assert permutation_pvalue(-1, null) == 1.0


def test_npc_detects_consistent_signal():
    rng = np.random.default_rng(1)
    null = rng.normal(0, 1, size=(2000, 4))
    observed = np.array([3.0, 2.8, 3.1, 2.9])
    result = npc_test(observed, null, method="liptak")
    assert result["combined_pvalue"] < 0.01
    assert len(result["component_pvalues"]) == 4


def test_npc_null_is_calibrated():
    rng = np.random.default_rng(2)
    null = rng.normal(0, 1, size=(2000, 4))
    observed = np.array([0.1, -0.2, 0.3, 0.0])
    result = npc_test(observed, null, method="fisher")
    assert 0.05 < result["combined_pvalue"] < 0.95


def test_maxt_adjusted_bounds():
    rng = np.random.default_rng(3)
    null = rng.normal(size=(500, 3))
    observed = np.array([0.0, 2.0, -1.0])
    adjusted = maxt_adjusted(observed, null)
    assert adjusted.shape == (3,)
    assert np.all((adjusted > 0) & (adjusted <= 1))
    assert adjusted[1] < adjusted[0]
