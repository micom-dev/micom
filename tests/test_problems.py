from .fixtures import community
import micom.problems as probs
from micom.optcom import optcom
import numpy as np
import pytest

# Only test linear ones
strategies = ["original", "lmoma"]


@pytest.mark.parametrize("strategy", strategies)
def test_is_consistent(community, strategy):
    community.solver = "glpk"
    sol = optcom(community, strategy=strategy)
    assert np.allclose(sol.growth_rate, 0.873922)
    assert np.allclose(sol.members.growth_rate.dropna().sum(), 4 * 0.873922)


@pytest.mark.parametrize("f", [0, 0.25, "0.5", [0, 0.2, 0.4, 0], np.ones(4) * 0.7])
def test_good_mingrowth(community, f):
    community.solver = "glpk"
    sol = optcom(community, strategy="lmoma", min_growth=f)
    assert np.allclose(sol.growth_rate, 0.873922)
    assert np.allclose(sol.members.growth_rate.dropna().sum(), 4 * 0.873922)


@pytest.mark.parametrize("f", ["a", [1, 2], np.ones(6)])
def test_bad_mingrowth(community, f):
    community.solver = "glpk"
    with pytest.raises(ValueError):
        sol = optcom(community, strategy="lmoma", min_growth=f)


@pytest.mark.parametrize("strategy", strategies)
def test_must_share(community, strategy):
    community.solver = "glpk"
    community.reactions.EX_glc__D_m.lower_bound = -5
    sol = optcom(community, strategy=strategy, min_growth=1.0)
    assert sol is None
    sol = optcom(community, strategy=strategy, min_growth=0.0)
    assert sol.growth_rate > 0.4
    assert (sol.members.growth_rate.dropna() > 0.4).any()


def test_conservation(community):
    community.solver = "glpk"
    community.reactions.EX_glc__D_m.lower_bound = -5
    sol = optcom(community, strategy="lmoma", fluxes=True)
    imports = sol.fluxes[sol.fluxes.index.str.startswith("EX_glc__D_e")][0:4].values
    total_influx = community.abundances.dot(imports)
    assert np.allclose(total_influx, -5)
