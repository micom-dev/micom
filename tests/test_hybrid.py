"""Test the hybrid solver.

We only test the core functions with relaxed accuracies as hybrid is not good with LPs.
"""

from .fixtures import batch, community
import cobra.util.solver as su
import micom.data as md
from micom.batch import Batch, GrowthResults, Configuration
from micom.qiime_formats import load_qiime_medium
from micom.solution import CommunitySolution, OptimizationError
import pytest
from pytest import approx

pytestmark = pytest.mark.skipif(
    "hybrid" not in su.solvers, reason="hybrid not functional here"
)

data = md.test_data()
medium = load_qiime_medium(md.test_medium)
db = md.test_db
config = Configuration(solver="hybrid")


@pytest.mark.parametrize("strategy", ["none", "minimal imports", "pFBA"])
def test_grow(tmp_path, strategy):
    batch = Batch(data, medium=medium, model_db=db, config=config)
    batch.build(str(tmp_path))
    grown = batch.grow()
    assert batch.has_results
    assert isinstance(grown, GrowthResults)
    assert "growth_rate" in grown.growth_rates.columns
    assert "flux" in grown.exchanges.columns


def test_tradeoff(tmp_path):
    batch = Batch(data, medium=medium, model_db=db, config=config)
    batch.build(str(tmp_path))
    rates = batch.tradeoff()
    assert "growth_rate" in rates.columns
    assert "tradeoff" in rates.columns
    assert rates.dropna().shape[0] < rates.shape[0]


def test_media(tmp_path):
    batch = Batch(data, medium=medium, model_db=db, config=config)
    batch.build(str(tmp_path))
    media = batch.minimal_medium(0.5)
    assert media.shape[0] > 3
    assert "flux" in media.columns
    assert "reaction" in media.columns


def test_complete_community_medium(tmp_path):
    batch = Batch(data, medium=medium, model_db=db, config=config)
    batch.build(str(tmp_path))
    bad_medium = medium.iloc[0:2, :]
    fixed = batch.complete_medium(0.5, 0.001, medium=bad_medium)
    assert fixed.shape[0] > 3
    assert "description" in fixed.columns


def test_community_objective(community):
    community.solver = "hybrid"
    x = community.optimize()
    y = community.optimize(fluxes=True)
    assert isinstance(x, CommunitySolution)
    assert x.growth_rate == approx(0.873922, 1e-3, 1e-3)
    assert x.members.growth_rate.dropna().sum() == approx(4 * 0.873922, 1e-3, 1e-3)
    assert isinstance(y, CommunitySolution)
    assert y.fluxes.shape[0] == 5


def test_cooperative_tradeoff(community):
    community.solver = "hybrid"
    sol = community.cooperative_tradeoff(fraction=1.0)
    for g in sol.members.growth_rate.dropna():
        assert g == approx(0.874, 1e-3, 1e-3)


def test_multiple_tradeoffs(community):
    community.solver = "hybrid"
    fs = [1.0, 0.5, 0.3, 0.1]
    sol = community.cooperative_tradeoff(fraction=fs)
    assert sol.shape == (4, 2)
    for i, f in sol.iterrows():
        assert f.solution.growth_rate == approx(fs[i] * 0.873922, 1e-3, 1e-3)
