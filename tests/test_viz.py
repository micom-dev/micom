"""Test visualization."""

from fixtures import growth_data, tradeoff_data, check_viz
import micom.viz as viz
from os import path
import pandas as pd
import pytest
import sys


def test_plot_growth(growth_data, tmp_path):
    v = viz.plot_growth(growth_data.growth_rates, str(tmp_path))
    check_viz(v)


def test_plot_tradeoff(tradeoff_data, tmp_path):
    v = viz.plot_tradeoff(tradeoff_data, str(tmp_path))
    check_viz(v)


def test_plot_sample_exchanges(growth_data, tmp_path):
    v = viz.plot_exchanges_per_sample(growth_data.exchanges, str(tmp_path))
    check_viz(v)
    v = viz.plot_exchanges_per_sample(
        growth_data.exchanges, str(tmp_path), direction="export")
    check_viz(v)
    v = viz.plot_exchanges_per_sample(
        growth_data.exchanges, str(tmp_path), cluster=False)
    check_viz(v)
    with pytest.raises(ValueError):
        v = viz.plot_exchanges_per_sample(
            growth_data.exchanges, str(tmp_path), direction="dog")


@pytest.mark.skipif(sys.platform == "darwin",
                    reason="llvmlite problems on MacOS")
def test_plot_taxon_exchanges(growth_data, tmp_path):
    v = viz.plot_exchanges_per_taxon(growth_data.exchanges, str(tmp_path))
    check_viz(v)
    v = viz.plot_exchanges_per_taxon(
        growth_data.exchanges, str(tmp_path), direction="export")
    check_viz(v)
    with pytest.raises(ValueError):
        v = viz.plot_exchanges_per_taxon(
            growth_data.exchanges, str(tmp_path), direction="dog")
        print(v.data)


def test_fit(growth_data, tmp_path):
    meta = pd.Series([0, 0, 1, 1],
                     index=growth_data.growth_rates.sample_id.unique())
    v = viz.plot_fit(growth_data.exchanges, meta, out_folder=str(tmp_path),
                     min_coef=0)
    check_viz(v)
    v = viz.plot_fit(growth_data.exchanges, meta, variable_type="continuous",
                     out_folder=str(tmp_path), min_coef=0)
    check_viz(v)

    with pytest.raises(ValueError):
        v = viz.plot_fit(growth_data.exchanges, meta, variable_type="dog",
                         out_folder=str(tmp_path))
