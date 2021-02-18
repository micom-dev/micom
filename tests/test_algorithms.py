"""Test the micom algorithms."""

from .fixtures import community
import numpy as np
import micom
import micom.algorithms as algo

inclusion = np.ones((5, 95), dtype=int)


def test_jaccard():
    j = algo.jaccard(inclusion)
    assert np.allclose(j, np.zeros((5, 5)))


def test_euclidean():
    j = algo.euclidean(inclusion)
    assert np.allclose(j, np.zeros((5, 5)))


def test_reaction_matrix():
    tax = micom.data.test_taxonomy()
    m = algo.reaction_matrix([tax.loc[0, "file"]])
    assert np.allclose(m, np.ones((1, 95)))


def test_metabolic_dist(community):
    md = algo.metabolic_dist(community.reactions)
    real = np.zeros((5, 5))
    real[4, :] = 1
    real[:, 4] = 1
    real[4, 4] = 0
    assert np.allclose(md, real)
