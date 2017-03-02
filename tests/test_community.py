"""Tests for basic construction of a community."""

from fixtures import community
from micom import Community
from micom.data import test_taxonomy
import numpy as np


def test_construction():
    tax = test_taxonomy()
    nr, mr = tax.reactions, tax.metabolites
    com = Community(tax)
    assert len(com.objectives) == 5
    assert len(com.taxonomy) == 5
    assert len(com.reactions) == tax.reactions.sum()
    assert len(com.metabolites) == tax.metabolites.sum()


def test_abundance_cutoff():
    tax = test_taxonomy(n=3)
    tax["abundance"] = [1.0, 2.0, 1e-6]
    com = Community(tax)
    assert len(com.objectives) == 2
    assert len(com.taxonomy) == 2


def test_abundances(community):
    assert np.allclose(community.abundances, np.ones(5)/5)

    ab = np.array([1.0, 2.0, 1e-6, 3.0, 1e-8])
    expected = np.array([1.0/6, 2.0/6, 1e-6, 3.0/6, 1e-6])
    community.abundances = ab
    assert np.allclose(community.abundances, expected)
