"""Provides support for Qiime formats."""

from micom.util import load_pickle
import os
from os import path
import pandas as pd
from zipfile import ZipFile
from ruamel.yaml import YAML
from tempfile import TemporaryDirectory

yaml = YAML()
_has_manifest = ["CommunityModels[Pickle]", "MetabolicModels[JSON]",
                 "MetabolicModels[SBML]"]


def metadata(artifact):
    """Read metadata from a Qiime 2 artifact."""
    with ZipFile(artifact) as zf:
        files = zf.namelist()
        meta = [fi for fi in files
                if "metadata.yaml" in fi and "provenance" not in fi]
        if len(meta) == 0:
            raise ValueError(
                "%s is not a valid Qiime 2 artifact :(" % artifact)
        with zf.open(meta[0]) as mf:
            meta = yaml.load(mf)
        return meta


def load_qiime_model_db(artifact, extract_path):
    """Prepare a model database for use."""
    if not path.exists(extract_path):
        os.mkdir(extract_path)
    meta = metadata(artifact)
    if meta["type"] != "MetabolicModels[JSON]":
        raise ValueError("%s is not a q2-micom model database :(" % artifact)
    uuid = meta["uuid"]
    with ZipFile(artifact) as zf:
        zf.extractall(extract_path)
    manifest = pd.read_csv(
        path.join(extract_path, uuid, "data", "manifest.csv"))
    manifest["file"] = [
        path.join(extract_path, uuid, "data", f)
        for f in manifest.file]
    return manifest


def load_qiime_manifest(artifact):
    """Prepare community models for use."""
    meta = metadata(artifact)
    if meta["type"] not in _has_manifest:
        raise ValueError(
            "%s is not a supported q2-micom artifact :(" % artifact)
    uuid = meta["uuid"]
    with ZipFile(artifact) as zf, TemporaryDirectory(prefix="micom_") as td:
        zf.extract(uuid + "/data/manifest.csv", str(td))
        manifest = pd.read_csv(
            path.join(str(td), uuid, "data", "manifest.csv"))
    return manifest


def load_qiime_model(artifact, id):
    """Load a model from a Qiime 2 artifact."""
    meta = metadata(artifact)
    if meta["type"] != "CommunityModels[Pickle]":
        raise ValueError(
            "%s is not a q2-micom community model collection :(" % artifact)
    uuid = meta["uuid"]
    with ZipFile(artifact) as zf, TemporaryDirectory(prefix="micom_") as td:
        try:
            zf.extract("%s/data/%s.pickle" % (uuid, id), str(td))
        except Exception:
            raise ValueError(
                "Could not extract model with ID `%s` :(. "
                "Are you sure the ID is valid?" % id)
        model = load_pickle(path.join(str(td), uuid, "data", "%s.pickle" % id))
    return model


def load_qiime_medium(artifact):
    """Load a growth medium/diet from a Qiime 2 artifact."""
    meta = metadata(artifact)
    if not meta["type"].startswith("MicomMedium["):
        raise ValueError(
            "%s is not a q2-micom medium :(" % artifact)
    uuid = meta["uuid"]
    with ZipFile(artifact) as zf, TemporaryDirectory(prefix="micom_") as td:
        zf.extract(uuid + "/data/medium.csv", str(td))
        medium = pd.read_csv(
            path.join(str(td), uuid, "data", "medium.csv"))
    medium.index = medium.reaction
    return medium


def load_qiime_feature_table(artifact):
    """Load a feature table from a Qiime 2 artifact."""
    try:
        import biom
    except ImportError:
        raise ImportError(
            "Reading Qiime 2 FeatureTables requires the `biom-format` package."
            "You can install it with:\n pip install numpy Cython\n"
            "pip install biom-format")
    meta = metadata(artifact)
    if not meta["type"].startswith("FeatureTable["):
        raise ValueError(
            "%s is not a Qiime 2 FeatureTable :(" % artifact)
    uuid = meta["uuid"]
    with ZipFile(artifact) as zf, TemporaryDirectory(prefix="micom_") as td:
        zf.extract(uuid + "/data/feature-table.biom", str(td))
        table = biom.load_table(
            path.join(str(td), uuid, "data", "feature-table.biom"))
    return table


def load_qiime_taxonomy(artifact):
    """Load taxonomy feature data from a Qiime 2 artifact."""
    meta = metadata(artifact)
    if not meta["type"].startswith("FeatureData[Taxonomy]"):
        raise ValueError(
            "%s is not a Qiime 2 FeatureData object :(" % artifact)
    uuid = meta["uuid"]
    with ZipFile(artifact) as zf, TemporaryDirectory(prefix="micom_") as td:
        zf.extract(uuid + "/data/taxonomy.tsv", str(td))
        taxa = pd.read_csv(
            path.join(str(td), uuid, "data", "taxonomy.tsv"),
            sep="\t",
            index_col=0
        )["Taxon"]
    return taxa
