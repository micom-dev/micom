"""Helper fixtures for mico."""

import micom
import micom.data as md
from micom.batch import Batch
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
def community_with_host():
    """A simple community containing 4 species."""
    return micom.Community(micom.data.test_taxonomy(host=True), progress=False)


@pytest.fixture
def results():
    """A more complex results example."""
    res = md.test_results()
    return res


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
    batch = Batch(md.test_data(), medium=medium, model_db=md.test_db, progress=False)
    batch.build(str(tmp_path))
    batch.grow()
    return batch


@pytest.fixture
def tradeoff_data(tmp_path):
    """Generate some growth simulation data."""
    data = md.test_data()
    batch = Batch(md.test_data(), medium=medium, model_db=md.test_db, progress=False)
    batch.tradeoff()
    return batch
