"""Test the inclusion of host models."""

import pytest
import micom as mm


def test_host_model():
    """Test that the host model is included in the community."""
    com = mm.Community(mm.data.test_taxonomy(host=True), progress=False)
    assert len(com.taxa) == 4
    assert len(com.host) == 1
    assert "h" in com.compartments

def test_host_db():
    """Test that an error is raised if no model db is provided."""
    tax = mm.data.test_taxonomy(host=True)
    host_db = {"human": tax.loc[tax.is_host, "file"].values[0]}
    com = mm.Community(tax, model_db=mm.data.test_db, host_db=host_db, progress=False)
    assert com.host is not None
    assert com.host_abundances.shape[0] == 1

def test_missing_host_db():
    """Test that an error is raised if no host model is provided."""
    tax = mm.data.test_taxonomy(host=True)
    del tax["file"]
    with pytest.raises(ValueError):
        mm.Community(tax, model_db=mm.data.test_db, progress=False)

def test_host_medium():
    """Test that the host model is included in the medium."""
    com = mm.Community(mm.data.test_taxonomy(host=True), progress=False)
    com.host_medium = {"EX_glc__D_h": -10}
    assert com.host_medium["EX_glc__D_h"] == -10

    com.host_medium = {"EX_glc__D_h": -10, "blub": -100}
    assert len(com.host_medium) == 1

    with pytest.raises(ValueError):
        com.host_medium = {"EX_glc__D_e": -10, "EX_o2_e": -10}

def test_host_optimization():
    """Test that the host model is included in the optimization."""
    com = mm.Community(mm.data.test_taxonomy(host=True), progress=False)
    com.host_medium = {"EX_o2_h": -10}
    sols = com.optimize_host()
    assert "human" in sols
    assert all(s.host_rates.sum() > 0 for s in sols.values())

def test_host_repr():
    """Test that the host model is included in the repr."""
    com = mm.Community(mm.data.test_taxonomy(host=True), progress=False)
    com.host_medium = {"EX_o2_h": -10}
    sol = com.optimize_host()["human"]
    assert "host" in repr(sol)
    assert "host maintenance" in sol._repr_html_()