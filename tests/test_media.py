"""Test growth media agorithms."""

from fixtures import community
import numpy as np
import micom.media as media


def test_medium_linear(community):
    medium = media.minimal_medium(community, 0.8, 0.1)
    assert len(medium) <= 4
    assert all(medium > 1e-9)


def test_medium_mip(community):
    medium = media.minimal_medium(community, 0.8, 0.1,
                                    minimize_components=True)
    assert len(medium) <= 4
    assert all(medium > 1e-9)

    # Anaerobic growth
    medium = media.minimal_medium(community, 0.1, 0.1,
                                    minimize_components=True)
    assert len(medium) <= 3
    assert all(medium > 1e-9)


def test_complete(community):
    m = media.minimal_medium(community, 0.85, 0.85,
                             minimize_components=True)
    medium = media.complete_medium(community, m[0:2], 0.8, max_import=20)
    assert len(medium) > 2


def test_complete_mip(community):
    m = media.minimal_medium(community, 0.85, 0.85,
                             minimize_components=True)
    medium = media.complete_medium(community, m[0:2], 0.8, max_import=20,
                                   minimize_components=True)
    assert len(medium) > 2
