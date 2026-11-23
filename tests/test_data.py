import numpy as np
import pandas as pd

from replica.data import parse_event, parse_receptor, parse_stage_metabric, parse_stage_tcga


def test_parse_event():
    series = pd.Series(["1:DECEASED", "0:LIVING", None, "weird"])
    out = parse_event(series)
    assert out.iloc[0] == 1
    assert out.iloc[1] == 0
    assert np.isnan(out.iloc[2])
    assert np.isnan(out.iloc[3])


def test_parse_stage_tcga():
    series = pd.Series(["STAGE I", "Stage IIA", "STAGE IIIC", "STAGE IV", "STAGE X"])
    out = parse_stage_tcga(series)
    assert out.tolist()[:4] == [1, 2, 3, 4]
    assert np.isnan(out.iloc[4])


def test_parse_stage_metabric():
    series = pd.Series(["0", "1", "2", "3", "4", None])
    out = parse_stage_metabric(series)
    assert np.isnan(out.iloc[0])
    assert out.iloc[1] == 1
    assert out.iloc[4] == 4


def test_parse_receptor():
    series = pd.Series(["Positive", "Negative", "Indeterminate", None])
    out = parse_receptor(series)
    assert out.iloc[0] == 1
    assert out.iloc[1] == 0
    assert np.isnan(out.iloc[2])
    assert np.isnan(out.iloc[3])
