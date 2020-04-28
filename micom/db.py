"""Build a database of organism metabolic models."""

import pandas as pd
import os
from os import path
from zipfile import ZipFile


def load_manifest(folder):
    """Get the manifest from a model DB."""
    mpath = path.join(folder, "manifest.csv")
    if not path.exists(mpath):
        raise ValueError(
            "No manifest found. `%s` does not look like a valid "
            "model database." % folder
        )
    manifest = pd.read_csv(path.join(folder, "manifest.csv"))
    if "file" not in manifest.columns:
        raise ValueError("Invalid manifest for model database :(")
    manifest.file = [path.join(folder, f) for f in manifest.file]
    return manifest


def load_zip_model_db(artifact, extract_path):
    """Prepare a model database for use."""
    if not path.exists(extract_path):
        os.mkdir(extract_path)
    with ZipFile(artifact) as zf:
        zf.extractall(extract_path)
    manifest = load_manifest(extract_path)
    manifest["file"] = [path.join(extract_path, f) for f in manifest.file]
    return manifest
