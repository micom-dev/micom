from fixtures import community
import micom.problems as probs
import numpy as np
import pytest

# Only test linear ones
strategies = ['original', 'lmoma', 'linear', 'linear lagrangian']


class TestOptcom():

    @pytest.mark.parametrize("strategy", strategies)
    def test_is_consistent(self, community, strategy):
        sol = community.optcom(strategy=strategy)
        assert np.allclose(sol.growth_rate, 0.873922)
        assert np.allclose(sol.members.growth_rate, 0.873922)

    @pytest.mark.parametrize("f", [0, 0.25, "0.5", [0, 0.2, 0.4, 0, 0],
                                   np.ones(5) * 0.7])
    def test_good_mingrowth(self, community, f):
        sol = community.optcom(strategy="linear", min_growth=f)
        assert np.allclose(sol.growth_rate, 0.873922)
        assert np.allclose(sol.members.growth_rate, 0.873922)

    @pytest.mark.parametrize("f", ["a", [1, 2], np.ones(6)])
    def test_bad_mingrowth(self, community, f):
        with pytest.raises(ValueError):
            sol = community.optcom(strategy="linear", min_growth=f)

    @pytest.mark.parametrize("strategy", strategies)
    def test_must_share(self, community, strategy):
        community.reactions.EX_glc__D_m.lower_bound = -5
        sol = community.optcom(strategy=strategy, min_growth=1.0)
        assert sol is None
        sol = community.optcom(strategy=strategy, min_growth=0.0)
        assert sol.growth_rate > 0.4
        assert (sol.members.growth_rate > 0.4).any()

    def test_conservation(self, community):
        community.reactions.EX_glc__D_m.lower_bound = -5
        sol = community.optcom(strategy="linear", fluxes=True)
        imports = sol["EX_glc__D_e"][0:5]
        total_influx = community.abundances.dot(imports)
        assert np.allclose(total_influx, -5)
