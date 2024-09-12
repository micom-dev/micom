"""Test growth media agorithms."""

from .fixtures import community
import numpy as np
import micom.media as media
from pandas.testing import assert_series_equal
import pytest


def test_medium_linear(community):
    medium = media.minimal_medium(community, 0.8, 0.8)
    assert len(medium) <= 4
    assert all(medium > 1e-9)


def test_nonmatching_medium(community):
    with pytest.raises(ValueError):
        community.medium = {"EX_blabla": 100}
    assert "EX_blabla" not in community.medium


def test_partial_matching_medium(community):
    community.medium = {"EX_glc__D_m": 10, "EX_blabla_m": 2}


def test_matching_medium(community):
    community.medium = {"EX_glc__D_m": 10}


def test_medium_mip(community):
    medium = media.minimal_medium(community, 0.8, 0.8, minimize_components=True)
    assert len(medium) <= 4
    assert all(medium > 1e-9)

    # Anaerobic growth
    medium = media.minimal_medium(community, 0.1, 0.1, minimize_components=True)
    assert len(medium) <= 3
    assert all(medium > 1e-9)


def test_medium_mass(community):
    medium = media.minimal_medium(community, 0.8, 0.8, weights="mass")
    assert len(medium) <= 4
    medium["EX_o2_m"] > 20.0

    # Anaerobic growth
    medium = media.minimal_medium(community, 0.1, 0.1, weights="mass")
    assert len(medium) <= 4
    medium["EX_o2_m"] > 4.0


def test_medium_element(community):
    medium = media.minimal_medium(community, 0.8, 0.8, weights="C")
    assert len(medium) <= 4
    assert medium["EX_glc__D_m"] < 9.5

    # Anaerobic growth
    medium = media.minimal_medium(community, 0.1, 0.1, weights="C")
    assert len(medium) <= 4
    assert medium["EX_glc__D_m"] < 2.0


def test_medium_wrong_element(community):
    with pytest.raises(ValueError):
        medium = media.minimal_medium(community, 0.8, 0.8, weights="Cat")


def test_complete_strict(community):
    m = media.minimal_medium(community, 0.85, 0.85, minimize_components=False)
    medium = media.complete_medium(
        community, m[0:2], 0.8, max_import=20, strict=m.index
    )
    assert len(medium) > 2


def test_complete_weights(community):
    m = media.minimal_medium(community, 0.85, 0.85, minimize_components=False)
    medium = media.complete_medium(
        community, m[0:2], 0.8, max_import=20, strict=m.index, weights="C"
    )
    assert len(medium) == 4


def test_complete_non_strict(community):
    m = media.minimal_medium(community, 0.85, 0.85, minimize_components=False)
    # request growth rates not fulfillable with previous bounds
    medium = media.complete_medium(community, m, 0.95, 0.95, max_import=20)
    community.medium = medium
    new_min = media.minimal_medium(community, 0.95, 0.95)
    assert len(medium) > 2
    print(new_min - medium)
    assert_series_equal(new_min, medium)  # check if indeed minimal


def test_complete_mip_strict(community):
    m = media.minimal_medium(community, 0.85, 0.85, minimize_components=True)
    medium = media.complete_medium(
        community,
        m[0:2],
        0.8,
        0.8,
        max_import=20,
        strict=m.index,
        minimize_components=True,
    )
    assert len(medium) == 4


def test_complete_mip_non_strict(community):
    m = media.minimal_medium(community, 0.85, 0.85, minimize_components=True)
    medium = media.complete_medium(
        community, m, 0.95, 0.95, max_import=20, minimize_components=True
    )
    assert len(medium) == 5
