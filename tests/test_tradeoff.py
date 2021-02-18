"""Test cooperative tradeoff."""

from .fixtures import community
import micom.qiime_formats as qf
from os import path
from pytest import mark, raises, approx


def test_cooperative_tradeoff(community):
    sol = community.cooperative_tradeoff(fraction=1.0)
    for g in sol.members.growth_rate.dropna():
        assert g == approx(0.874, 0.001)


def test_multiple_tradeoffs(community):
    fs = [1.0, 0.5, 0.3, 0.1]
    sol = community.cooperative_tradeoff(fraction=fs)
    assert sol.shape == (4, 2)
    for i, f in sol.iterrows():
        assert f.solution.growth_rate == approx(fs[i] * 0.874, 0.001)
