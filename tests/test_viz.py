"""Test visualization."""

from .fixtures import batch_grown, batch_tradeoff, check_viz
import micom.viz as viz
from os import path
import pandas as pd
import pytest
import random


def random_groups(res):
    """Create some random groups for samples."""
    groups = pd.Series(
        random.choices(["a", "b"], k=4),
        index=res.growth_rates.sample_id.unique(),
        name="random",
    )
    return groups


def test_plot_growth(batch_grown, tmp_path):
    v = viz.plot_growth(batch_grown.results, str(tmp_path / "viz.html"))
    check_viz(v)
    v = viz.plot_growth(
        batch_grown.results,
        str(tmp_path / "viz.html"),
        groups=random_groups(batch_grown.results),
    )
    check_viz(v)


def test_plot_tradeoff(batch_tradeoff, tmp_path):
    v = viz.plot_tradeoff(batch_tradeoff.tradeoffs, str(tmp_path / "viz.html"))
    check_viz(v)


def test_plot_sample_exchanges(batch_grown, tmp_path):
    v = viz.plot_exchanges_per_sample(batch_grown.results, str(tmp_path / "viz.html"))
    check_viz(v)
    v = viz.plot_exchanges_per_sample(
        batch_grown.results, str(tmp_path / "viz.html"), direction="export"
    )
    check_viz(v)
    v = viz.plot_exchanges_per_sample(
        batch_grown.results, str(tmp_path / "viz.html"), cluster=False
    )
    check_viz(v)
    with pytest.raises(ValueError):
        v = viz.plot_exchanges_per_sample(
            batch_grown.results, str(tmp_path / "viz.html"), direction="dog"
        )


def test_plot_taxon_exchanges(batch_grown, tmp_path):
    v = viz.plot_exchanges_per_taxon(batch_grown.results, str(tmp_path / "viz.html"))
    check_viz(v)
    v = viz.plot_exchanges_per_taxon(
        batch_grown.results, str(tmp_path / "viz.html"), direction="export"
    )
    check_viz(v)
    v = viz.plot_exchanges_per_taxon(
        batch_grown.results,
        str(tmp_path / "viz.html"),
        groups=random_groups(batch_grown.results),
    )
    check_viz(v)
    with pytest.raises(ValueError):
        v = viz.plot_exchanges_per_taxon(
            batch_grown.results, str(tmp_path / "viz.html"), direction="dog"
        )


def test_association(batch_grown, tmp_path):
    meta = pd.Series(
        [0, 0, 1, 1], index=batch_grown.results.growth_rates.sample_id.unique()
    )
    v = viz.plot_association(
        batch_grown.results,
        meta,
        filename=str(tmp_path / "viz.html"),
        fdr_threshold=0.5,
    )
    check_viz(v)
    v = viz.plot_association(
        batch_grown.results,
        meta,
        variable_type="continuous",
        filename=str(tmp_path / "viz.html"),
        fdr_threshold=0.5,
    )
    check_viz(v)

    with pytest.raises(ValueError):
        v = viz.plot_association(
            batch_grown.results,
            meta,
            variable_type="dog",
            filename=str(tmp_path / "viz.html"),
        )


def test_association_fillna(batch_grown, tmp_path):
    meta = pd.Series(
        [0, 0, 1, 1], index=batch_grown.results.growth_rates.sample_id.unique()
    )
    v = viz.plot_association(
        batch_grown.results,
        meta,
        fillna=1e-6,
        filename=str(tmp_path / "viz.html"),
        fdr_threshold=0.5,
    )
    check_viz(v)

    with pytest.raises(ValueError):
        v = viz.plot_association(
            batch_grown.results,
            meta,
            fillna=1e-6,
            variable_type="dog",
            filename=str(tmp_path / "viz.html"),
        )


def test_plot_focal_interactions(batch_grown, tmp_path):
    v = viz.plot_focal_interactions(
        batch_grown.results, taxon="strain_b", filename=str(tmp_path / "viz.html")
    )
    check_viz(v)


def test_plot_mes(batch_grown, tmp_path):
    v = viz.plot_mes(batch_grown.results, filename=str(tmp_path / "viz.html"))
    check_viz(v)
    v = viz.plot_mes(
        batch_grown.results,
        groups=random_groups(batch_grown.results),
        filename=str(tmp_path / "viz.html"),
    )
    check_viz(v)
