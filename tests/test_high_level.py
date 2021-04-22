"""Test the high level API."""

import micom.data as md
from micom.workflows import (
    build_database,
    build,
    grow,
    tradeoff,
    minimal_media,
    fix_medium,
)
from micom.qiime_formats import load_qiime_medium, load_qiime_manifest
from micom.solution import OptimizationError
import pytest

medium = load_qiime_medium(md.test_medium)
db = md.test_db


def test_db(tmp_path):
    manifest = load_qiime_manifest(db)
    tax = md.test_taxonomy()
    manifest.file = tax.file[0]
    built = build_database(manifest, str(tmp_path), rank="species")
    assert built.shape[0] == 3
    for fi in built.file:
        assert (tmp_path / fi).exists()
    built = build_database(manifest, str(tmp_path / "db.zip"))
    assert (tmp_path / "db.zip").exists()


def test_build(tmp_path):
    data = md.test_data()
    built = build(data, db, str(tmp_path), cutoff=0, threads=1)
    assert built.shape[0] == 4
    assert "sample_id" in built.columns
    assert "found_fraction" in built.columns
    assert "file" in built.columns
    for fi in built.file:
        assert (tmp_path / fi).exists()


def test_build_no_db(tmp_path):
    data = md.test_data()
    built = build(data, None, str(tmp_path), cutoff=0)
    assert built.shape[0] == 4
    assert "sample_id" in built.columns
    assert "found_fraction" not in built.columns
    assert "file" in built.columns
    for fi in built.file:
        assert (tmp_path / fi).exists()


@pytest.mark.parametrize("strategy", ["none", "minimal imports", "pFBA"])
def test_grow(tmp_path, strategy):
    data = md.test_data()
    built = build(data, db, str(tmp_path), cutoff=0)
    grown = grow(built, str(tmp_path), medium, 0.5, strategy=strategy)
    assert len(grown) == 3
    assert "growth_rate" in grown.growth_rates.columns
    assert "flux" in grown.exchanges.columns
    with pytest.raises(OptimizationError):
        grow(built, str(tmp_path), medium, 1.5)


def test_grow_bad_strategy(tmp_path):
    data = md.test_data()
    built = build(data, db, str(tmp_path), cutoff=0)
    with pytest.raises(ValueError):
        grown = grow(built, str(tmp_path), medium, 0.5, strategy="blub")


def test_tradeoff(tmp_path):
    data = md.test_data()
    built = build(data, db, str(tmp_path), cutoff=0)
    rates = tradeoff(built, str(tmp_path), medium)
    assert "growth_rate" in rates.columns
    assert "tradeoff" in rates.columns
    assert rates.dropna().shape[0] < rates.shape[0]


def test_media(tmp_path):
    data = md.test_data()
    built = build(data, db, str(tmp_path), cutoff=0)
    media = minimal_media(built, str(tmp_path), 0.5)
    assert media.shape[0] > 3
    assert "flux" in media.columns
    assert "reaction" in media.columns


def test_fix_medium(tmp_path):
    data = md.test_data()
    built = build(data, db, str(tmp_path), cutoff=0)
    bad_medium = medium.iloc[0:2, :]
    fixed = fix_medium(built, str(tmp_path), bad_medium, 0.5, 0.001, 10)
    assert fixed.shape[0] > 3
    assert "description" in fixed.columns
