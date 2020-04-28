"""Worflow to build models for several samples."""

from cobra.io import read_sbml_model, save_json_model
from micom.logger import logger
from micom.util import join_models
from micom import Community
from micom.workflows.core import workflow
import os
import pandas as pd
from tempfile import TemporaryDirectory
from zipfile import ZipFile


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
        {s: os.path.join(out_folder, s + ".pickle") for s in samples}
    )
    args = [
        [s, taxonomy[taxonomy.sample_id == s], out_path[s]] for s in samples
    ]

    workflow(build_and_save, args, threads)
    taxonomy["file"] = taxonomy.sample_id + ".pickle"
    taxonomy.to_csv(os.path.join(out_folder, "manifest.csv"), index=False)
    return taxonomy


REQ_FIELDS = pd.Series(
    [
        "file",
        "kingdom",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
    ]
)


def _reduce_group(df):
    keep = df.columns[df.nunique() == 1]
    new = df.iloc[0, :][keep]
    new["file"] = "|".join(df.file.astype(str))
    return pd.DataFrame.from_records([new])


def _summarize_models(args):
    tid, row, new_path = args
    files = row["file"].split("|")
    if len(files) > 1:
        mod = join_models(files, id=tid)
    else:
        mod = read_sbml_model(files[0])
    save_json_model(mod, new_path)


def build_database(
    manifest, out_path, rank="genus", threads=1, compress=None, progress=True
):
    """Create a model database from a set of SBML files.

    Note
    ----
    A manifest for the joined models will also be written to the output folder
    as "manifest.csv". This may contain NA entries for additional columns
    that had different values within the summarized models.

    Parameters
    ----------
    manifest : pandas.DataFrame
        A manifest of SBML files containing their filepath as well as taxonomy.
        Must contain the columns "file", "kingdom", "phylum", "class",
        "order", "family", "genus", and "species". May contain additional
        columns.
    out_path : str
        The directory where the joined models will be written.
    compress : bool
        Whether to compress the output. Default is True if out_path ends with
        ".zip" otherwise no.
    progress : bool
        Whether to show a progress bar.

    Returns
    -------
    pd.DataFrame
        The manifest of the joined models. Will still contain information
        from the original metadata.
    """
    meta = manifest.copy()
    meta.columns = meta.columns.str.lower()
    compress = out_path.endswith(".zip")

    if not REQ_FIELDS.isin(meta.columns).all():
        raise ValueError(
            "Metadata File needs to have the following "
            "columns %s." % ", ".join(REQ_FIELDS)
        )
    bad = meta.file.apply(lambda x: not os.path.exists(x))
    if any(bad):
        raise ValueError(
            "The following models are in the manifest but do "
            "not exist at the specified path: %s" % meta.file[bad]
        )

    meta = meta.groupby(rank).apply(_reduce_group).reset_index(drop=True)
    logger.info("Building %d models on rank `%s`." % (meta.shape[0], rank))
    meta.index = meta[rank].str.replace("[^\\w\\_]", "_")
    meta["id"] = meta.index
    meta["summary_rank"] = rank

    if compress:
        with TemporaryDirectory(prefix="micom_") as tdir:
            args = [
                (tid, row, os.path.join(tdir, "%s.json" % tid))
                for tid, row in meta.iterrows()
            ]
            workflow(_summarize_models, args, threads)
            meta.file = meta.index + ".json"
            meta.to_csv(os.path.join(tdir, "manifest.csv"), index=False)
            with ZipFile(out_path, "w") as zf:
                [zf.write(a[2], os.path.basename(a[2])) for a in args]
                zf.write(os.path.join(tdir, "manifest.csv"), "manifest.csv")
    else:
        args = [
            (tid, row, os.path.join(out_path, "%s.json" % tid))
            for tid, row in meta.iterrows()
        ]
        workflow(_summarize_models, args, threads)
        meta.file = meta.index + ".json"
        meta.to_csv(os.path.join(out_path, "manifest.csv"), index=False)

    return meta
