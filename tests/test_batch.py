"""Test the Batch interface."""

import pytest
import pandas as pd

from .fixtures import batch, batch_built
from micom.batch import GrowthResults, Batch, Configuration
import micom.data as md
from micom.solution import OptimizationError


def _mark_built(batch, tmp_path):
    batch.build_manifest = pd.DataFrame(
        {"sample_id": ["sample"], "file": ["sample.pickle"]}
    )
    batch.out_folder = tmp_path

def test_build(batch, tmp_path, caplog):
    batch.build(str(tmp_path))
    manifest = batch.build_manifest

    assert batch.is_built
    assert manifest.shape[0] == 4
    assert "sample_id" in manifest.columns
    assert "found_fraction" in manifest.columns
    assert "file" in manifest.columns
    for fi in manifest.file:
        assert (tmp_path / fi).exists()
    batch.build(str(tmp_path))
    assert "Found existing models for 4 samples." in caplog.text


def test_build_no_db(tmp_path):
    data = md.test_data(uses_db=False)
    config = Configuration()
    config.dbs.microbial = None
    batch = Batch(data, medium=md.test_medium, model_db=None, config=config)
    built = batch.build(str(tmp_path))
    assert built.shape[0] == 4
    assert "sample_id" in built.columns
    assert "found_fraction" not in built.columns
    assert "file" in built.columns
    for fi in built.file:
        assert (tmp_path / fi).exists()


@pytest.mark.parametrize("strategy", ["none", "minimal imports", "pFBA"])
def test_grow(batch_built, strategy):
    batch_built.config.simulation.flux_method = strategy
    grown = batch_built.grow()
    assert isinstance(grown, GrowthResults)
    assert isinstance(batch_built.results, GrowthResults)
    assert "growth_rate" in grown.growth_rates.columns
    assert "flux" in grown.exchanges.columns


def test_tradeoff(batch_built):
    rates = batch_built.tradeoff()
    assert "growth_rate" in rates.columns
    assert "tradeoff" in rates.columns
    assert rates.dropna().shape[0] < rates.shape[0]


def test_media(batch_built):
    media = batch_built.minimal_medium(0.5)
    assert media.shape[0] > 3
    assert "flux" in media.columns
    assert "reaction" in media.columns


def test_media_no_summary(batch_built):
    media = batch_built.minimal_medium(0.5, summarize=False)
    assert media.shape[0] > 12
    assert "flux" in media.columns
    assert "reaction" in media.columns
    assert "sample_id" in media.columns


@pytest.mark.parametrize("w", [None, "mass", "C"])
def test_media_weights(batch_built, w):
    media = batch_built.minimal_medium(0.5, minimize=w)
    assert media.shape[0] > 3
    assert "flux" in media.columns
    assert "reaction" in media.columns


def test_complete_community_medium(batch_built):
    bad_medium = batch_built.medium.iloc[0:2, :]
    fixed = batch_built.complete_medium(community_growth=0.5, taxa_growth=0.001, medium=bad_medium)
    assert fixed.shape[0] > 3
    assert "description" in fixed.columns



def test_complete_community_medium_no_summary(batch_built):
    bad_medium = batch_built.medium.iloc[0:2, :]
    fixed = batch_built.complete_medium(community_growth=0.5, taxa_growth=0.001, medium=bad_medium, summarize=False)
    assert fixed.shape[0] > 12
    assert "description" in fixed.columns
    assert "sample_id" in fixed.columns


def test_batch_not_built_errors(batch):
    with pytest.raises(ValueError, match="has not been built"):
        batch.grow()
    with pytest.raises(ValueError, match="has not been built"):
        batch.tradeoff()
    with pytest.raises(ValueError, match="has not been built"):
        batch.minimal_medium()
    with pytest.raises(ValueError, match="has not been built"):
        batch.complete_medium()


def test_batch_state_and_representations_incomplete(batch):
    assert not batch.is_built
    assert not batch.has_results
    assert not batch.has_tradeoffs
    assert "Batch" in repr(batch)
    assert "samples" in str(batch)
    assert "<table>" in batch._repr_html_()
    with pytest.raises(AttributeError, match="cannot be changed"):
        batch.taxonomy = batch.taxonomy

def test_batch_state_and_representations_complete(batch_built):
    batch_built.grow()
    batch_built.tradeoff()
    assert batch_built.is_built
    assert batch_built.has_results
    assert batch_built.has_tradeoffs
    assert "Batch" in repr(batch_built)
    assert "samples" in str(batch_built)
    assert "μᵢ" in str(batch_built)
    assert "μᵢ" in batch_built._repr_html_()
    assert "<table>" in batch_built._repr_html_()
    assert "growing fraction" in batch_built._repr_html_()
    assert "growing fraction" in str(batch_built)
    with pytest.raises(AttributeError, match="cannot be changed"):
        batch_built.taxonomy = batch_built.taxonomy


@pytest.mark.parametrize("tradeoffs", [[-0.1], [1.1]])
def test_tradeoff_rejects_invalid_values(batch_built, tradeoffs):
    with pytest.raises(ValueError, match="tradeoff values must between 0 and 1"):
        batch_built.tradeoff(tradeoffs)


def test_grow_preconditions(batch, tmp_path):
    _mark_built(batch, tmp_path)
    batch._medium = None
    with pytest.raises(ValueError, match="No medium has been specified"):
        batch.grow()

    batch.medium = md.test_medium
    strategy = batch.config.simulation.strategy
    batch.config.simulation.strategy = "SteadyCom"
    with pytest.raises(NotImplementedError, match="SteadyCom is not yet implemented"):
        batch.grow()
    batch.config.simulation.strategy = strategy


def test_grow_all_samples_fail(batch, tmp_path, monkeypatch):
    _mark_built(batch, tmp_path)
    monkeypatch.setattr("micom.batch.batch.workflow", lambda *args, **kwargs: [None])
    with pytest.raises(OptimizationError, match="All numerical optimizations failed"):
        batch.grow()


def test_tradeoff_all_samples_fail(batch, tmp_path, monkeypatch):
    _mark_built(batch, tmp_path)
    monkeypatch.setattr("micom.batch.batch.workflow", lambda *args, **kwargs: [None])
    with pytest.raises(OptimizationError, match="All numerical optimizations failed"):
        batch.tradeoff([0.5])


def test_minimal_medium_partial_failure(batch, tmp_path, monkeypatch, caplog):
    _mark_built(batch, tmp_path)
    medium = pd.DataFrame(
        {"reaction": ["EX_test"], "flux": [1.0], "sample_id": ["sample"]}
    )
    monkeypatch.setattr(
        "micom.batch.batch.workflow",
        lambda *args, **kwargs: [None, {"medium": medium}],
    )
    result = batch.minimal_medium()
    assert result.reaction.tolist() == ["EX_test"]
    assert "some samples" in caplog.text


def test_minimal_medium_all_samples_fail(batch, tmp_path, monkeypatch):
    _mark_built(batch, tmp_path)
    monkeypatch.setattr("micom.batch.batch.workflow", lambda *args, **kwargs: [None])
    with pytest.raises(OptimizationError, match="Could not find a growth medium"):
        batch.minimal_medium()


def test_complete_medium_all_samples_fail(batch, tmp_path, monkeypatch):
    _mark_built(batch, tmp_path)
    monkeypatch.setattr("micom.batch.batch.workflow", lambda *args, **kwargs: [None])
    with pytest.raises(OptimizationError, match="All optimizations failed"):
        batch.complete_medium()