"""Test optimization functions."""

from .fixtures import community
import cobra
from cobra.util.solver import solvers
import numpy as np
import pytest
from micom.solution import CommunitySolution

stable = ["glpk", "cplex"]
solvers = [s for s in solvers.keys() if s in stable]


def test_community_objective(community):
    x = community.optimize()
    y = community.optimize(fluxes=True)
    assert isinstance(x, CommunitySolution)
    assert np.allclose(x.growth_rate, 0.873922)
    assert np.allclose(x.members.growth_rate.dropna().sum(), 4 * 0.873922)
    assert isinstance(y, CommunitySolution)
    assert y.fluxes.shape[0] == 5


def test_individual_objective(community):
    growth_rates = community.optimize_all()
    assert np.allclose(growth_rates, 4 * 0.873922)
