"""Helper fixtures for mico."""

import pytest
import micom


@pytest.fixture
def community():
    """A simple community containing 3 species."""
    return micom.Community(micom.data.test_taxonomy())
