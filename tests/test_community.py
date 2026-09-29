"""Tests for basic construction of a community."""

from .fixtures import community, community_with_host
from micom import Community, load_pickle
from micom.data import test_taxonomy
import numpy as np
import pandas as pd
import pytest


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
    assert np.allclose(community.microbial_abundances, np.ones(4) / 4)

    ab = np.array([1.0, 2.0, 1e-8, 3.0])
    series = pd.Series(ab, index=community.taxa)
    expected = np.array([1.0 / 6, 2.0 / 6, 1e-6, 3.0 / 6])
    community.microbial_abundances = series
    assert np.allclose(community.microbial_abundances, expected)

    with pytest.raises(TypeError):
        community.microbial_abundances = ab


def test_exchanges(community):
    assert "glc__D_m" in community.metabolites
    assert "EX_glc__D_m" in community.reactions
    r = community.reactions.EX_glc__D_e__strain_a
    glc_e = community.metabolites.get_by_id("glc__D_e__strain_a")
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


@pytest.mark.parametrize(
    "strategy", ["resource constraint", "resource coupling", "enzyme coupling"]
)
def test_add_coupling_constraints(community, strategy):
    community.add_coupling_constraints(
        ids=[community.taxa[0]], strategy=strategy, include_exchanges=True
    )
    assert any(
        constraint.name.startswith("resource_constraint__")
        or constraint.name.startswith("coupling_")
        for constraint in community.constraints
    )


def test_microbial_abundance_validation_and_clamping(community):
    with pytest.raises(TypeError, match="must be a pandas Series"):
        community.set_microbial_abundance([1.0])
    with pytest.raises(ValueError, match="not in the community"):
        community.set_microbial_abundance(pd.Series({"unknown": 1.0}))

    taxon = community.taxa[0]
    community.set_microbial_abundance(pd.Series({taxon: 0.0}))
    assert community.microbial_abundances[taxon] >= community._rtol

    community.set_microbial_abundance(pd.Series({taxon: 0.3}), normalize=False)
    assert community.microbial_abundances[taxon] == 0.3


def test_build_metrics_require_database(community):
    with pytest.raises(ValueError, match="only available for models build"):
        community.build_metrics


def test_host_abundance_validation_and_update(community_with_host):
    with pytest.raises(TypeError, match="must be a pandas Series"):
        community_with_host.set_host_abundance([1.0])
    with pytest.raises(ValueError, match="not in the community"):
        community_with_host.set_host_abundance(pd.Series({"unknown": 1.0}))

    host = community_with_host.host[0]
    community_with_host.set_host_abundance(pd.Series({host: 0.5}))
    assert community_with_host.host_abundances[host] == 1.0
