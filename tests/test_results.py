"""Test for the GrowthResults and conversions."""

from micom.batch import GrowthResults, save_results, load_results
from .fixtures import community, batch_grown
import pytest


def test_conversion(community):
    sol = community.cooperative_tradeoff(fraction=0.5, fluxes=False)
    with pytest.raises(ValueError):
        results = GrowthResults.from_solution(sol, community)
    sol = community.cooperative_tradeoff(fraction=0.5, fluxes=True)
    results = GrowthResults.from_solution(sol, community)
    assert hasattr(results, "growth_rates")
    assert results.growth_rates.shape[0] == 4
    assert hasattr(results, "exchanges")
    assert "glc__D_m" in results.exchanges.metabolite.values
    assert hasattr(results, "annotations")
    assert "glc__D_m" in results.annotations.metabolite.values


def test_add(community):
    sol = community.cooperative_tradeoff(fraction=0.5, fluxes=True)
    community.id = "sample1"
    r1 = GrowthResults.from_solution(sol, community)
    community.id = "sample2"
    r2 = GrowthResults.from_solution(sol, community)
    combined = r1 + r2
    assert combined.growth_rates.shape[0] == 8
    assert combined.exchanges.shape[0] == 2 * r1.exchanges.shape[0]
    assert combined.annotations.shape == r1.annotations.shape


def test_results_saving(batch_grown, tmp_path):
    results_file = tmp_path / "test.zip"
    batch_grown.results.save(results_file)
    assert results_file.exists()

    grown = batch_grown.results
    save_results(grown, results_file)
    loaded = load_results(results_file)
    assert isinstance(loaded, GrowthResults)
    assert loaded.growth_rates.shape == grown.growth_rates.shape
    assert loaded.exchanges.shape == grown.exchanges.shape
    assert loaded.annotations.shape == grown.annotations.shape
