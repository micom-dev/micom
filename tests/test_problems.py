from fixtures import community
import micom.problems as probs
import numpy as np
import pytest

# Only test linear ones
strategies = ['original', 'lmoma', 'linear', 'linear lagrangian']


class TestLinearOptcom():

    @pytest.mark.parametrize("strategy", strategies)
    def test_is_consistent(self, community, strategy):
        cobj, gcs = community.optcom(strategy=strategy)
        assert np.allclose(cobj, 0.873922)
        assert np.allclose(gcs, 0.873922)

    @pytest.mark.parametrize("f", [0, 0.25, "0.5", [0, 0.2, 0.4, 0, 0],
                                   np.ones(5) * 0.7])
    def test_good_mingrowth(self, community, f):
        cobj, gcs = community.optcom(min_growth=f)
        assert np.allclose(cobj, 0.873922)
        assert np.allclose(gcs, 0.873922)

    @pytest.mark.parametrize("f", ["a", [1, 2], np.ones(6)])
    def test_bad_mingrowth(self, community, f):
        with pytest.raises(ValueError):
            cobj, gcs = community.optcom(min_growth=f)

    @pytest.mark.parametrize("strategy", strategies)
    def test_must_share(self, community, strategy):
        community.reactions.EX_glc__D_m.lower_bound = -5
        cobj, gcs = community.optcom(min_growth=1.0)
        assert cobj is None
        assert np.isnan(gcs).all()
        cobj, gcs = community.optcom(min_growth=0.0)
        assert cobj > 0.4
        assert (gcs > 0.4).any()

    def test_conservation(self, community):
        community.reactions.EX_glc__D_m.lower_bound = -5
        cobj, fluxes = community.optcom(strategy="linear", fluxes=True)
        total_influx = community.abundances.dot(fluxes["EX_glc__D_e"])
        assert np.allclose(total_influx, -5)
