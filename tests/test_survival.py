import numpy as np
import pandas as pd

from replica.survival import (
    cindex,
    fit_cox,
    logrank_statistic,
    permutation_logrank,
    risk_score,
)


def _synthetic(n=240, seed=7):
    rng = np.random.default_rng(seed)
    group = rng.choice(["LumA", "LumB", "Basal", "Her2"], n)
    hazard = pd.Series({"LumA": 1.0, "LumB": 1.3, "Basal": 2.0, "Her2": 1.6}).loc[group].to_numpy()
    time = rng.exponential(1.0 / hazard)
    censor = rng.exponential(4.0, n)
    event = (time <= censor).astype(int)
    time = np.minimum(time, censor)
    return pd.DataFrame({"time": time, "event": event, "subtype4": group})


def test_logrank_detects_separation():
    frame = _synthetic()
    stat = logrank_statistic(frame.time, frame.event, frame.subtype4)
    assert stat > 0


def test_permutation_logrank_is_calibrated_under_null():
    rng = np.random.default_rng(11)
    n = 180
    time = rng.exponential(10, n)
    event = rng.integers(0, 2, n)
    group = rng.choice(["A", "B", "C"], n)
    result = permutation_logrank(time, event, group, n_perm=199, seed=3)
    assert 0.05 <= result["p_value"] <= 0.95


def test_cox_and_cindex():
    frame = _synthetic()
    fit = fit_cox(frame, "time", "event", [], subtype_col="subtype4")
    risk = risk_score(fit, frame)
    assert cindex(frame.time.to_numpy(), frame.event.to_numpy(), risk.to_numpy()) > 0.5
    assert "subtype_Basal" in fit.summary().index
