"""Test the stats helpers."""

from multiprocessing.sharedctypes import Value
import micom.stats as ms
import numpy as np
import pandas as pd
from pytest import approx, raises


def test_fdr_correct():
    p = [0.01, 0.1, 0.3]
    assert ms.fdr_adjust(p) == approx([0.03, 0.15, 0.3])
    p = [0.01, 0.1, 0.1]
    assert ms.fdr_adjust(p) == approx([0.03, 0.1, 0.1])
    with raises(ValueError):
        ms.fdr_adjust([0.2, 0.3, float("nan")])


single_df = lambda: pd.DataFrame(
    {
        "metabolite": [f"metabolite_{i}" for i in range(1, 5)],
        "flux": np.random.normal(100, 10, 4),
    }
)


def make_grouped_fluxes(n=10):
    np.random.seed(42)
    dfs = []
    for group in range(1, 4):
        for i in range(1, n + 1):
            df = single_df()
            df["group"] = f"group_{group}"
            df["sample_id"] = f"sample_{(group - 1) * n + i}"
            dfs.append(df)
    dfs = pd.concat(dfs)
    dfs.loc[
        dfs.group.str.contains("3") & dfs.metabolite.str.contains("1|2"), "flux"
    ] += 100.0
    return dfs


def make_correlated_fluxes(n=10):
    np.random.seed(42)
    dfs = []
    for t in range(1, 9):
        for i in range(1, n + 1):
            df = single_df()
            df["time"] = t
            df["sample_id"] = f"sample_{(t - 1) * n + i}"
            dfs.append(df)
    dfs = pd.concat(dfs)
    dfs.loc[dfs.metabolite.str.contains("1|2"), "flux"] += (
        dfs.loc[dfs.metabolite.str.contains("1|2"), "time"] * 10
    )
    return dfs


def test_comparison_binary():
    data = make_grouped_fluxes()
    tests = ms.compare_groups(
        data, "group", groups=["group_1", "group_3"], progress=False
    )
    assert "metabolite_1" in tests[tests.p < 0.01].metabolite.values
    assert "metabolite_2" in tests[tests.p < 0.01].metabolite.values
    assert "metabolite_3" in tests[tests.p > 0.01].metabolite.values
    assert "metabolite_4" in tests[tests.p > 0.01].metabolite.values
    assert "log_fold_change" in tests.columns
    assert "comparison" in tests.columns


def test_comparison_many():
    data = make_grouped_fluxes()
    tests = ms.compare_groups(data, "group", progress=False)
    assert "metabolite_1" in tests[tests.p < 0.01].metabolite.values
    assert "metabolite_2" in tests[tests.p < 0.01].metabolite.values
    assert "metabolite_3" in tests[tests.p > 0.01].metabolite.values
    assert "metabolite_4" in tests[tests.p > 0.01].metabolite.values
    assert "log_mean_std" in tests.columns
    assert "covariate" in tests.columns


def test_correlation():
    data = make_correlated_fluxes()
    tests = ms.correlate_fluxes(data, "time", progress=False)
    assert "metabolite_1" in tests[tests.p < 0.01].metabolite.values
    assert "metabolite_2" in tests[tests.p < 0.01].metabolite.values
    assert "metabolite_3" in tests[tests.p > 0.01].metabolite.values
    assert "metabolite_4" in tests[tests.p > 0.01].metabolite.values
    assert "covariate" in tests.columns
    assert (tests.covariate == "time").all()
