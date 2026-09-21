"""Test the inclusion of host models."""

import pytest
import micom as mm


def test_host_model():
    """Test that the host model is included in the community."""
    com = mm.Community(mm.data.test_taxonomy(host=True), progress=False)
    assert len(com.taxa) == 4
    assert len(com.host) == 1
    assert "h" in com.compartments