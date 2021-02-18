"""Test interventions."""

from .fixtures import community
from micom.elasticity import elasticities
from pytest import approx


def test_elasticities(community):
    el = elasticities(community, fraction=1.0)
    assert "reaction" in el.columns
    assert "effector" in el.columns
    s = el[(el.reaction == "EX_glc__D_m") & (el.effector == "EX_glc__D_m")]
    print(s)
    assert s.elasticity.iloc[0] == approx(1.0)
