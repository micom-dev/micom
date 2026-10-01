"""Test media helpers for model databases."""

import micom.data as md
import micom.batch.db_media as db_media
from micom.qiime_formats import load_qiime_medium
from micom.batch import complete_db_medium
import pytest
import pandas as pd

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


def _manifest():
    return pd.DataFrame(
        {
            "id": ["model_a", "model_b"],
            "file": ["a.json", "b.json"],
            "summary_rank": ["genus", "genus"],
        }
    )


def test_check_db_medium_classifies_growth(monkeypatch):
    monkeypatch.setattr(db_media, "load_manifest", lambda _: _manifest())
    monkeypatch.setattr(
        db_media,
        "workflow",
        lambda *args, **kwargs: [
            {"id": "model_a", "growth_rate": 0.2},
            {"id": "model_b", "growth_rate": 0.0},
        ],
    )
    result = db_media.check_db_medium("models", medium)
    assert result.can_grow.tolist() == [True, False]


def test_complete_db_medium_collects_partial_results(monkeypatch):
    monkeypatch.setattr(db_media, "load_manifest", lambda _: _manifest())
    imports = [
        (True, 1, 0.5, pd.Series({"EX_a": 0.5})),
        (False, float("nan"), float("nan"), pd.Series(dtype=float)),
    ]
    monkeypatch.setattr(db_media, "workflow", lambda *args, **kwargs: imports)
    manifest, fixed = db_media.complete_db_medium("models", medium)
    assert manifest.can_grow.tolist() == [True, False]
    assert manifest.added.iloc[0] == 1
    assert fixed.index.tolist() == ["model_a", "model_b"]


def test_try_complete_failure(monkeypatch):
    monkeypatch.setattr(db_media, "load_model", lambda _: object())
    monkeypatch.setattr(db_media, "find_external_compartment", lambda _: "e")

    def fail(*args, **kwargs):
        raise db_media.OptimizationError("infeasible")

    monkeypatch.setattr(db_media.mm, "complete_medium", fail)
    can_grow, added, flux, fixed = db_media._try_complete(
        ("model.json", pd.Series({"EX_a": 1.0}), 0.1, 1.0, False, None, [])
    )
    assert not can_grow
    assert pd.isna(added) and pd.isna(flux)
    assert fixed.isna().all()


def test_try_complete_success(monkeypatch):
    monkeypatch.setattr(db_media, "load_model", lambda _: object())
    monkeypatch.setattr(db_media, "find_external_compartment", lambda _: "C_e")
    monkeypatch.setattr(
        db_media.mm,
        "complete_medium",
        lambda *args, **kwargs: pd.Series({"EX_a_e": 2.0, "EX_b_e": 3.0}),
    )
    can_grow, added, added_flux, fixed = db_media._try_complete(
        ("model.json", pd.Series({"EX_a_e": 1.0}), 0.1, 1.0, False, None, [])
    )
    assert can_grow
    assert added == 1
    assert added_flux == 4.0
    assert fixed.index.tolist() == ["EX_a_m", "EX_b_m"]


def test_db_annotations_deduplicates(monkeypatch):
    manifest = _manifest()
    manifest["summary_rank"] = "genus"
    monkeypatch.setattr(db_media, "load_manifest", lambda _: manifest)
    monkeypatch.setattr(
        db_media,
        "workflow",
        lambda *args, **kwargs: [
            pd.DataFrame(
                {
                    "reaction": ["EX_a", "EX_b"],
                    "metabolite": ["a", "b"],
                }
            ),
            pd.DataFrame({"reaction": ["EX_a"], "metabolite": ["a"]}),
        ],
    )
    result = db_media.db_annotations("models")
    assert result.reaction.tolist() == ["EX_a", "EX_b"]


def test_grow_helper_warns_without_matching_exchanges(monkeypatch, caplog):
    model = type(
        "Model",
        (),
        {"exchanges": [], "slim_optimize": lambda self: 0.0},
    )()
    monkeypatch.setattr(db_media, "load_model", lambda _: model)
    result = db_media._grow(("model_a", "model.json", pd.Series({"EX_a": 1.0})))
    assert result == {"id": "model_a", "growth_rate": 0.0}
    assert "Could not find any reactions" in caplog.text


@pytest.mark.parametrize("suffix", [".zip", ".qza"])
def test_check_db_medium_compressed_dispatch(monkeypatch, tmp_path, suffix):
    manifest = _manifest()
    loader = "load_zip_model_db" if suffix == ".zip" else "load_qiime_model_db"
    monkeypatch.setattr(db_media, loader, lambda *args: manifest)
    monkeypatch.setattr(
        db_media,
        "workflow",
        lambda *args, **kwargs: [
            {"id": "model_a", "growth_rate": 0.2},
            {"id": "model_b", "growth_rate": 0.0},
        ],
    )
    result = db_media.check_db_medium(f"models{suffix}", medium)
    assert result.can_grow.tolist() == [True, False]


@pytest.mark.parametrize("suffix", [".zip", ".qza"])
def test_complete_db_medium_compressed_dispatch(monkeypatch, suffix):
    manifest = _manifest()
    loader = "load_zip_model_db" if suffix == ".zip" else "load_qiime_model_db"
    monkeypatch.setattr(db_media, loader, lambda *args: manifest)
    monkeypatch.setattr(
        db_media,
        "workflow",
        lambda *args, **kwargs: [
            (True, 1, 0.5, pd.Series({"EX_a": 0.5})),
            (False, float("nan"), float("nan"), pd.Series(dtype=float)),
        ],
    )
    result, imports = db_media.complete_db_medium(f"models{suffix}", medium)
    assert result.can_grow.tolist() == [True, False]
    assert imports.index.tolist() == ["model_a", "model_b"]


@pytest.mark.parametrize("suffix", [".zip", ".qza"])
def test_db_annotations_compressed_dispatch(monkeypatch, suffix):
    manifest = _manifest()
    manifest["summary_rank"] = "genus"
    loader = "load_zip_model_db" if suffix == ".zip" else "load_qiime_model_db"
    monkeypatch.setattr(db_media, loader, lambda *args: manifest)
    monkeypatch.setattr(
        db_media,
        "workflow",
        lambda *args, **kwargs: [
            pd.DataFrame({"reaction": ["EX_a"], "metabolite": ["a"]}),
            pd.DataFrame({"reaction": ["EX_b"], "metabolite": ["b"]}),
        ],
    )
    result = db_media.db_annotations(f"models{suffix}")
    assert result.reaction.tolist() == ["EX_a", "EX_b"]
