"""Test optimization functions."""

from fixtures import community
import cobra
import numpy as np


def test_community_objective(community):
    x = community.optimize()

    assert isinstance(x, cobra.core.solution.Solution)
    assert np.allclose(x.f, 0.874 * len(community.taxonomy), 1e-3, 1e-3)


def test_benchmark_community_objective(community, benchmark):
    benchmark(community.optimize)


def test_individual_objective(community):
    growth_rates = community.optimize_all()
    gc = community.optimize_single(0)
    assert all(growth_rates > 0.5)
    assert np.allclose(gc, growth_rates[0])


def test_benchmark_individual(community, benchmark):
    benchmark(community.optimize_single, 0)
