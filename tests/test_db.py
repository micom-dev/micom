"""Test model db creation."""

from .fixtures import this_dir
import micom as mm
import micom.workflows as mw
import micom.db as mdb
from os import path
from pytest import approx, mark, raises

db = mm.data.test_db


def test_qiime_community():
    tax = mm.data.test_taxonomy()
    tax["species"] = tax["id"].str.replace("_", " ")
    tax["abundance"] = [1, 2, 3, 4]
    del tax["file"]
    com = mm.Community(tax, db, progress=False)
    assert len(com.abundances) == 4
    m = com.build_metrics
    assert m[0] == 4
    assert m[1] == 4
    assert m[2] == approx(1.0)
    assert m[3] == approx(1.0)


@mark.parametrize("rank", ["genus", "species"])
def test_dir_build(tmp_path, rank):
    manifest = mm.data.test_taxonomy()
    with raises(ValueError):
        dman = mw.build_database(manifest, str(tmp_path), rank=rank, progress=False)
    for co in ["kingdom", "phylum", "class", "order", "family"]:
        manifest[co] = "fake"
    dman = mw.build_database(manifest, str(tmp_path), rank=rank, progress=False)
    assert dman.shape[0] == 1
    assert path.exists(str(tmp_path / dman.file[0]))
    assert path.exists(str(tmp_path / "manifest.csv"))
    tax = mm.data.test_taxonomy()
    com = mm.Community(tax, str(tmp_path), progress=False)
    m = com.build_metrics
    assert m[0] == 4
    assert m[1] == 4
    assert m[2] == 1.0
    assert m[3] == 1.0


@mark.parametrize("rank", ["genus", "species"])
def test_zip_build(tmp_path, rank):
    manifest = mm.data.test_taxonomy()
    with raises(ValueError):
        dman = mw.build_database(manifest, str(tmp_path), rank=rank, progress=False)
    for co in ["kingdom", "phylum", "class", "order", "family"]:
        manifest[co] = "fake"
    dman = mw.build_database(
        manifest, str(tmp_path / "test.zip"), rank=rank, progress=False
    )
    assert dman.shape[0] == 1
    assert path.exists(str(tmp_path / "test.zip"))
    tax = mm.data.test_taxonomy()
    com = mm.Community(tax, str(tmp_path / "test.zip"), progress=False)
    m = com.build_metrics
    assert m[0] == 4
    assert m[1] == 4
    assert m[2] == 1.0
    assert m[3] == 1.0
