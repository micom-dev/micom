"""Test the interaction module."""

from .fixtures import results, batch_grown
import micom.interaction as mi
import pytest


def test_focal(results):
    """Test single focal taxon."""
    print(results.growth_rates.columns)
    ints = mi.interactions(results, taxa="s__Akkermansia_muciniphila", progress=False)
    assert all(ints.focal == "s__Akkermansia_muciniphila")
    assert ints.partner.nunique() > 10
    assert all(ints.flux > 0)


def test_taxon_correction(results):
    """Test taxon name correction."""
    ints = mi.interactions(results, taxa="s__Akkermansia muciniphila", progress=False)
    assert all(ints.focal == "s__Akkermansia_muciniphila")
    assert ints.partner.nunique() > 10
    assert all(ints.flux > 0)


def test_wrong_taxa(results):
    """Test incorrect taxa names."""
    with pytest.raises(ValueError):
        mi.interactions(results, taxa="s__Akkermansia_muciniphilos", progress=False)
    with pytest.raises(ValueError):
        mi.interactions(
            results, taxa=["s__Akkermansia_muciniphila", "blub"], progress=False
        )


def test_summary(results):
    """Test the the results summary."""
    ints = mi.interactions(results, taxa="s__Akkermansia_muciniphila", progress=False)
    summ = mi.summarize_interactions(ints)
    for col in ["mass_flux", "flux", "C_flux", "N_flux", "n_ints"]:
        assert col in summ.columns
    for cl in ["provided", "received", "co-consumed"]:
        assert cl in summ["class"].unique()
    assert all(summ.groupby(["sample_id", "focal", "partner"]).flux.count() <= 3)


def test_all_interactions(batch_grown):
    """Test all vs all."""
    ints = mi.interactions(batch_grown.results, taxa=None, progress=False)
    summ = mi.summarize_interactions(ints)
    assert ints.focal.nunique() == 4
    assert ints.partner.nunique() == 4
    assert all(ints.flux > 0)
    assert summ.focal.nunique() == 4
    assert summ.partner.nunique() == 4
    for col in ["mass_flux", "flux", "C_flux", "N_flux", "n_ints"]:
        assert col in summ.columns


def test_mes(results):
    """Test MES score."""
    mes = mi.MES(results)
    assert "MES" in mes.columns
    assert "metabolite" in mes.columns
    assert all(mes.MES >= 0)


def test_basic(results):
    """Test basic producer, consumer, and flux statistics."""
    basic = mi.basic(results)
    row = basic.iloc[0]
    exchanges = results.exchanges[
        (results.exchanges.metabolite == row.metabolite)
        & (results.exchanges.sample_id == row.sample_id)
        & (results.exchanges.taxon != "medium")
    ]
    producers = exchanges[exchanges.direction == "export"]
    consumers = exchanges[exchanges.direction == "import"]

    assert row.producers == producers.taxon.nunique()
    assert row.consumers == consumers.taxon.nunique()
    assert row.consumption_flux == pytest.approx(
        (consumers.abundance * consumers.flux.abs()).sum()
    )
    assert row.production_flux == pytest.approx(
        (producers.abundance * producers.flux.abs()).sum()
    )
