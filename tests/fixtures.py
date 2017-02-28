"""Helper fixtures for mico."""

import pytest
import mico


@pytest.fixture
def community():
    """A simple community containing 3 species."""
    return mico.Community(mico.data.test_taxonomy())
