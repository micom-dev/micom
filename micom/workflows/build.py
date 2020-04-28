"""Worflow to build models for several samples."""

from micom import Community
from micom.workflows.core import workflow
from os import path
import pandas as pd


def build_and_save(args):
    """Build a single community model."""
    s, tax, out = args
    com = Community(tax, id=s, progress=False)
    com.to_pickle(out)


def build(
    taxonomy, model_db, threads, out_folder, cutoff=0.0001,
):
    """Build the community models."""
    samples = taxonomy.sample_id.unique()
    out_path = pd.Series(
        {s: path.join(out_folder, s + ".pickle") for s in samples}
    )
    args = [
        [s, taxonomy[taxonomy.sample_id == s], out_path[s]] for s in samples
    ]

    workflow(build_and_save, args, threads)
    taxonomy["file"] = taxonomy.sample_id + ".pickle"
    taxonomy.to_csv(path.join(out_folder, "manifest.csv"), index=False)
    return taxonomy
