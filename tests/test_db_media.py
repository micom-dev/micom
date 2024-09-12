"""Test media helpers for model databases."""

import micom.data as md
from micom.qiime_formats import load_qiime_medium
from micom.workflows import complete_db_medium
import pytest

db = md.test_db
medium = load_qiime_medium(md.test_medium)
medium["global_id"] = medium["reaction"].replace("_m$", "_e", regex=True)


def test_complete_strict():
    pruned = medium.iloc[0:2]
    manifest, fixed = complete_db_medium(
        db, growth=0.85, medium=pruned, strict=pruned.global_id, max_added_import=20
    )
    assert fixed.shape[0] > 2
    assert manifest.can_grow.all()
    assert manifest.added.mean() > 0
    assert manifest.added_flux.mean() > 1.0


def test_complete_non_strict():
    manifest, fixed = complete_db_medium(
        db, growth=0.95, medium=medium, max_added_import=20
    )
    assert fixed.shape[0] > 2
    assert manifest.added.mean() > 0
    assert manifest.added_flux.mean() > 1.0


def test_complete_weighted():
    pruned = medium.iloc[0:2]
    manifest, fixed = complete_db_medium(
        db,
        growth=0.85,
        medium=pruned,
        strict=pruned.global_id,
        max_added_import=20,
        weights="mass",
    )
    assert fixed.shape[0] > 2
    assert manifest.can_grow.all()
    assert manifest.added.mean() > 0
    assert manifest.added_flux.mean() > 1.0
