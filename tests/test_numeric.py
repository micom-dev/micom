"""Test numerical stabilization methods."""

from .fixtures import community
import micom.solution as ms
from os import path
from pytest import mark, raises, approx


def test_crossover(community):
    sol = community.cooperative_tradeoff(fraction=1.0)
    gcs = sol.members.growth_rate.dropna()
    sol = ms.crossover(community, sol)
    for i, g in enumerate(gcs):
        assert sol.members.growth_rate[i] == approx(g)
    sol = ms.crossover(community, sol, fluxes=True)
    for i, g in enumerate(gcs):
        assert sol.members.growth_rate[i] == approx(g)
        assert sol.fluxes.shape[1] > 100


def test_reset(community):
    sol = community.optimize()
    ms.reset_solver(community)
    assert community.optimize().growth_rate == approx(sol.growth_rate)


def test_retry(community):
    sol = ms.optimize_with_retry(community)
    assert sol == approx(0.874, 0.001)


def test_fraction(community):
    sol = ms.optimize_with_fraction(community, 0.5)
    assert sol.growth_rate == approx(0.874, 0.001)
