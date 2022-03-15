"""Helper fixtures for mico."""

import micom
import micom.data as md
from micom.workflows import build, grow, tradeoff
import os.path as path
import pytest
import pandas as pd

this_dir, _ = path.split(__file__)
medium = micom.qiime_formats.load_qiime_medium(md.test_medium)


@pytest.fixture
def community():
    """A simple community containing 4 species."""
    return micom.Community(micom.data.test_taxonomy(), progress=False)


@pytest.fixture
def linear_community():
    """A simple community containing 4 species."""
    return micom.Community(micom.data.test_taxonomy(), progress=False, solver="glpk")


def check_viz(v):
    """Check a visualization."""
    for d in v.data:
        assert isinstance(v.data[d], pd.DataFrame)
    assert path.exists(v.filename)


@pytest.fixture
def growth_data(tmp_path):
    """Generate some growth simulation data."""
    data = md.test_data()
    built = build(data, md.test_db, str(tmp_path), cutoff=0)
    grown = grow(built, str(tmp_path), medium, 0.5)
    return grown


@pytest.fixture
def tradeoff_data(tmp_path):
    """Generate some growth simulation data."""
    data = md.test_data()
    built = build(data, md.test_db, str(tmp_path), cutoff=0)
    rates = tradeoff(built, str(tmp_path), medium)
    return rates
