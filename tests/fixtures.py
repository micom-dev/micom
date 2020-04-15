"""Helper fixtures for mico."""

import pytest
import micom
import os.path as path

this_dir, _ = path.split(__file__)

@pytest.fixture
def community():
    """A simple community containing 4 species."""
    return micom.Community(micom.data.test_taxonomy(), progress=False)
