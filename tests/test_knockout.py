"""Test knockouts."""

from .fixtures import community
import micom.solution as ms
import numpy as np
from os import path
from pytest import mark, raises, approx


def test_knockout(community):
    ko = community.knockout_taxa()
    assert ko.shape == (4, 4)
    for i in range(4):
        assert ko.iloc[i, i] == approx(-0.874, 0.001)
    print(ko)
    assert ko.sum().sum() > 0.0


def test_single_knockout(community):
    ko = community.knockout_taxa(taxa=community.taxa[0])
    assert ko.shape == (1, 4)
    assert ko.iloc[0, :].values == approx([-0.874, 0.305, 0.305, 0.305], 0.01)
