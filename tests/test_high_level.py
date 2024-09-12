"""Test the high level API."""

import micom.data as md
from micom.workflows import (
    build_database,
    build,
    grow,
    tradeoff,
    minimal_media,
    complete_community_medium,
    save_results,
    load_results,
    GrowthResults,
)
from micom.logger import logger
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
    assert built.shape[0] == 4
    for fi in built.file:
        assert (tmp_path / fi).exists()
    built = build_database(manifest, str(tmp_path / "db.zip"))
    assert (tmp_path / "db.zip").exists()


def test_build(tmp_path, caplog):
    data = md.test_data()
    built = build(data, db, str(tmp_path), cutoff=0, threads=1)
    assert built.shape[0] == 4
    assert "sample_id" in built.columns
    assert "found_fraction" in built.columns
    assert "file" in built.columns
    for fi in built.file:
        assert (tmp_path / fi).exists()
    second_build = build(data, db, str(tmp_path), cutoff=0, threads=1)
    print(logger.level)
    assert "Found existing models for 4 samples." in caplog.text


def test_build_no_db(tmp_path):
    data = md.test_data(uses_db=False)
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
    assert isinstance(grown, GrowthResults)
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


def test_media_no_summary(tmp_path):
    data = md.test_data()
    built = build(data, db, str(tmp_path), cutoff=0)
    media = minimal_media(built, str(tmp_path), growth=0.5, summarize=False)
    assert media.shape[0] > 3 * built.shape[0]
    assert "flux" in media.columns
    assert "reaction" in media.columns


def test_media_solution(tmp_path):
    data = md.test_data()
    built = build(data, db, str(tmp_path), cutoff=0)
    media, res = minimal_media(
        built, str(tmp_path), growth=0.5, summarize=False, solution=True
    )
    assert media.shape[0] > 3 * built.shape[0]
    assert "flux" in media.columns
    assert "reaction" in media.columns
    assert isinstance(res, GrowthResults)
    assert res.growth_rates.shape[0] == 12


@pytest.mark.parametrize("w", [None, "mass", "C"])
def test_media_weights(tmp_path, w):
    data = md.test_data()
    built = build(data, db, str(tmp_path), cutoff=0)
    media = minimal_media(
        built, str(tmp_path), community_growth=0.5, weights=w, summarize=True
    )
    assert media.shape[0] > 3
    assert "flux" in media.columns
    assert "reaction" in media.columns


def test_complete_community_medium(tmp_path):
    data = md.test_data()
    built = build(data, db, str(tmp_path), cutoff=0)
    bad_medium = medium.iloc[0:2, :]
    fixed = complete_community_medium(built, str(tmp_path), bad_medium, 0.5, 0.001, 10)
    assert fixed.shape[0] > 3
    assert "description" in fixed.columns


def test_complete_community_medium_no_summary(tmp_path):
    data = md.test_data()
    built = build(data, db, str(tmp_path), cutoff=0)
    bad_medium = medium.iloc[0:2, :]
    fixed = complete_community_medium(
        built, str(tmp_path), bad_medium, 0.5, 0.001, 10, summarize=False
    )
    assert fixed.shape[0] > 3 * built.shape[0]
    assert "description" in fixed.columns


def test_results_saving(tmp_path):
    data = md.test_data()
    built = build(data, db, str(tmp_path), cutoff=0)
    grown = grow(built, str(tmp_path), medium, 0.5)
    results_file = str(tmp_path / "test.zip")
    save_results(grown, results_file)
    assert (tmp_path / "test.zip").exists()
    loaded = load_results(results_file)
    assert isinstance(loaded, GrowthResults)
    assert loaded.growth_rates.shape == grown.growth_rates.shape
    assert loaded.exchanges.shape == grown.exchanges.shape
    assert loaded.annotations.shape == grown.annotations.shape
