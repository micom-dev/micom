from fixtures import community
import micom.problems as probs
import numpy as np
import pytest


class TestLinearOptcom():

    def test_is_consistent(self, community):
        cobj, gcs = probs.linear_optcom(community, 1.0)
        assert np.allclose(cobj, 0.873922)
        assert np.allclose(gcs, 0.873922)

    @pytest.mark.parametrize("f", [0, 1.0, "1.0", [0, 0.2, 0.4, 0, 0],
                                   np.ones(5)])
    def test_good_fractions(self, community, f):
        cobj, gcs = community.optcom(method="linear", fractions=f)
        assert np.allclose(cobj, 0.873922)
        assert np.allclose(cobj, 0.873922)

    @pytest.mark.parametrize("f", ["a", [1, 2], np.ones(6)])
    def test_bad_fractions(self, community, f):
        with pytest.raises(ValueError):
            cobj, gcs = community.optcom(method="linear", fractions=f)

    def test_must_share(self, community):
        community.reactions.EX_glc__D_m.lower_bound = -5
        cobj, gcs = community.optcom(method="linear", fractions=1.0)
        assert cobj is None
        assert np.isnan(gcs).all()
        cobj, gcs = community.optcom(method="linear", fractions=0.0)
        assert cobj > 0.4
        assert (gcs > 0.3).any()

    def test_conservation(self, community):
        community.reactions.EX_glc__D_m.lower_bound = -5
        cobj, fluxes = community.optcom(method="linear", fluxes=True)
        total_influx = community.abundances.dot(fluxes["EX_glc__D_e"])
        assert np.allclose(total_influx, -5)
