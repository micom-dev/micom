"""Test Qiime 2 compatibility."""

from .fixtures import this_dir
import micom.qiime_formats as qf
from micom.data import test_db, test_medium
from os import path, environ
from pytest import mark, raises, approx

db, medium = test_db, test_medium
models = path.join(this_dir, "data", "build.qza")


def test_qiime_db(tmp_path):
    meta = qf.metadata(db)
    assert "uuid" in meta
    assert meta["type"] == "MetabolicModels[JSON]"
    manifest = qf.load_qiime_model_db(db, str(tmp_path))
    assert manifest.shape[0] == 4
    assert all(path.exists(f) for f in manifest.file)


@mark.parametrize("arti", [db, models])
def test_good_manifest(arti):
    manifest = qf.load_qiime_manifest(arti)
    assert "file" in manifest.columns


def test_bad_manifest():
    with raises(ValueError):
        manifest = qf.load_qiime_manifest(medium)


def test_qiime_medium():
    m = qf.load_qiime_medium(medium)
    assert "reaction" in m.columns
    assert "flux" in m.columns


def test_qiime_model():
    manifest = qf.load_qiime_manifest(models)
    assert "sample_id" in manifest.columns
    for i in manifest.sample_id:
        com = qf.load_qiime_model(models, i)
        assert len(com.abundances) == 3
        assert com.optimize().growth_rate == approx(0.874, 0.001)
    with raises(ValueError):
        qf.load_qiime_model(models, "blub")
    with raises(ValueError):
        qf.load_qiime_model(db, "blub")
