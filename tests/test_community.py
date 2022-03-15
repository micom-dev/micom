"""Tests for basic construction of a community."""

from .fixtures import community
from micom import Community, load_pickle
from micom.data import test_taxonomy
import numpy as np


def test_construction():
    tax = test_taxonomy()
    com = Community(tax)
    assert len(com.taxa) == 4
    assert len(com.taxonomy) == 4
    assert len(com.reactions) > tax.reactions.sum()
    assert len(com.metabolites) > tax.metabolites.sum()


def test_abundance_cutoff():
    tax = test_taxonomy(n=3)
    tax["abundance"] = [1.0, 2.0, 1e-6]
    com = Community(tax)
    assert len(com.taxa) == 2
    assert len(com.taxonomy) == 2


def test_abundances(community):
    assert np.allclose(community.abundances, np.ones(4) / 4)

    ab = np.array([1.0, 2.0, 1e-8, 3.0])
    expected = np.array([1.0 / 6, 2.0 / 6, 1e-6, 3.0 / 6])
    community.abundances = ab
    assert np.allclose(community.abundances, expected)


def test_exchanges(community):
    assert "glc__D_m" in community.metabolites
    assert "EX_glc__D_m" in community.reactions
    r = community.reactions.EX_glc__D_e__Escherichia_coli_1
    glc_e = community.metabolites.get_by_id("glc__D_e__Escherichia_coli_1")
    glc_m = community.metabolites.get_by_id("glc__D_m")
    assert r.metabolites[glc_e] == -1
    assert r.metabolites[glc_m] == 0.25


def test_get_taxonomy(community):
    tax = community.taxonomy
    tax["id"] = "bla"
    assert all(tax["id"] != community.taxonomy["id"])


def test_community_pickling(community, tmpdir):
    filename = str(tmpdir.join("com.pickle"))
    community.to_pickle(filename)
    loaded = load_pickle(filename)
    assert len(community.reactions) == len(loaded.reactions)
