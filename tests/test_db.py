"""Test model db creation."""

from .fixtures import this_dir
import micom as mm
import micom.batch as mb
import micom.db as mdb
from os import path, environ
from pathlib import Path
from pytest import approx, mark, raises

db = mm.data.test_db


def test_community_model_db():
    tax = mm.data.test_taxonomy()
    del tax["file"]
    com = mm.Community(tax, db, progress=False)
    assert len(com.microbial_abundances) == 4
    m = com.build_metrics
    assert m[0] == 4
    assert m[1] == 4
    assert m[2] == approx(1.0)
    assert m[3] == approx(1.0)


@mark.parametrize("rank", ["genus", "species"])
def test_dir_build(tmp_path, rank):
    manifest = mm.data.test_taxonomy()
    with raises(ValueError):
        dman = mb.build_database(manifest, str(tmp_path), rank=rank, progress=False)
    for co in ["kingdom", "phylum", "class", "order", "family"]:
        manifest[co] = "fake"
    dman = mb.build_database(manifest, str(tmp_path), rank=rank, progress=False)
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
        dman = mb.build_database(manifest, str(tmp_path), rank=rank, progress=False)
    for co in ["kingdom", "phylum", "class", "order", "family"]:
        manifest[co] = "fake"
    dman = mb.build_database(
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

@mark.xfail(condition=environ.get("GITHUB_ACTIONS") == "true", reason="Fails on GitHub Actions")
@mark.parametrize("loc", ["default://agora103_gtdb207_genus_1.qza", "https://zenodo.org/records/7739096/files/agora103_gtdb207_genus_1.qza?download=1"])
def test_model_db_download(loc, tmp_path):
    db = mdb.get_database(loc, tmp_path)
    man = mm.qiime_formats.load_qiime_model_db(db, tmp_path / "model_db")
    assert all(man.summary_rank == "genus")
    assert "file" in man.columns
    assert "genus" in man.columns
    assert "family" in man.columns

@mark.xfail(condition=environ.get("GITHUB_ACTIONS") == "true", reason="Fails on GitHub Actions")
@mark.parametrize("loc", ["default://himalaya.qza", "https://raw.githubusercontent.com/micom-dev/media/refs/heads/main/media/vmh_high_fiber_agora.qza"])
def test_media_db_download(loc, tmp_path):
    db = mdb.get_database(loc, tmp_path, what="media")
    medium = mm.qiime_formats.load_qiime_medium(db)
    assert all(medium.flux > 0)
    assert "flux" in medium.columns
    assert "reaction" in medium.columns
    assert "metabolite" in medium.columns

def test_get_database_trivial():
    db = mdb.get_database(mm.data.test_db, "whatever")
    assert path.exists(db)
    assert Path(mm.data.test_db).resolve() == Path(db).resolve()

    db = mdb.get_database(None, "whatever")
    assert db is None