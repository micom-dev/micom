"""Test visualization."""

from .fixtures import growth_data, tradeoff_data, check_viz
import micom.viz as viz
from os import path
import pandas as pd
import pytest
import random


def random_groups(res):
    """Create some random groups for samples."""
    sids = res.growth_rates.sample_id.unique()
    groups = pd.Series(
        random.choices(["a", "b"], k=len(sids)),
        index=sids,
        name="random",
    )
    return groups


def test_plot_growth(growth_data, tmp_path):
    v = viz.plot_growth(growth_data, str(tmp_path / "viz.html"))
    check_viz(v)
    v = viz.plot_growth(
        growth_data, str(tmp_path / "viz.html"), groups=random_groups(growth_data)
    )
    check_viz(v)


def test_plot_tradeoff(tradeoff_data, tmp_path):
    v = viz.plot_tradeoff(tradeoff_data, str(tmp_path / "viz.html"))
    check_viz(v)


def test_plot_sample_exchanges(growth_data, tmp_path):
    v = viz.plot_exchanges_per_sample(growth_data, str(tmp_path / "viz.html"))
    check_viz(v)
    v = viz.plot_exchanges_per_sample(
        growth_data, str(tmp_path / "viz.html"), direction="export"
    )
    check_viz(v)
    v = viz.plot_exchanges_per_sample(
        growth_data, str(tmp_path / "viz.html"), cluster=False
    )
    check_viz(v)
    with pytest.raises(ValueError):
        v = viz.plot_exchanges_per_sample(
            growth_data, str(tmp_path / "viz.html"), direction="dog"
        )


def test_plot_taxon_exchanges(growth_data, tmp_path):
    v = viz.plot_exchanges_per_taxon(growth_data, str(tmp_path / "viz.html"))
    check_viz(v)
    v = viz.plot_exchanges_per_taxon(
        growth_data, str(tmp_path / "viz.html"), direction="export"
    )
    check_viz(v)
    v = viz.plot_exchanges_per_taxon(
        growth_data, str(tmp_path / "viz.html"), groups=random_groups(growth_data)
    )
    check_viz(v)
    with pytest.raises(ValueError):
        v = viz.plot_exchanges_per_taxon(
            growth_data, str(tmp_path / "viz.html"), direction="dog"
        )


def test_association(growth_data, tmp_path):
    meta = pd.Series([0, 0, 1, 1], index=growth_data.growth_rates.sample_id.unique())
    v = viz.plot_association(
        growth_data, meta, filename=str(tmp_path / "viz.html"), fdr_threshold=0.5
    )
    check_viz(v)
    v = viz.plot_association(
        growth_data,
        meta,
        variable_type="continuous",
        filename=str(tmp_path / "viz.html"),
        fdr_threshold=0.5,
    )
    check_viz(v)

    with pytest.raises(ValueError):
        v = viz.plot_association(
            growth_data, meta, variable_type="dog", filename=str(tmp_path / "viz.html")
        )


def test_plot_focal_interactions(growth_data, tmp_path):
    v = viz.plot_focal_interactions(
        growth_data, taxon="Escherichia_coli_2", filename=str(tmp_path / "viz.html")
    )
    check_viz(v)


def test_plot_mes(growth_data, tmp_path):
    v = viz.plot_mes(growth_data, filename=str(tmp_path / "viz.html"))
    check_viz(v)
    v = viz.plot_mes(
        growth_data,
        groups=random_groups(growth_data),
        filename=str(tmp_path / "viz.html"),
    )
    check_viz(v)
