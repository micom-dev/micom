"""Worflow to build models for several samples."""

from cobra.io import read_sbml_model, save_json_model
from micom.logger import logger
from micom.util import join_models
from micom.community import Community, _ranks
from micom.workflows.core import workflow
import os
import pandas as pd
from tempfile import TemporaryDirectory
from zipfile import ZipFile


def _reduce_group(df):
    keep = df.columns[df.nunique() == 1]
    new = df.iloc[0, :][keep]
    if "file" in df.columns:
        new["file"] = "|".join(df.file.astype(str))
    return pd.DataFrame.from_records([new])


def build_and_save(args):
    """Build a single community model."""
    s, tax, db, out, cutoff, solver = args
    com = Community(
        tax, model_db=db, id=s, progress=False, rel_threshold=cutoff, solver=solver
    )
    com.to_pickle(out)
    if db is None:
        metrics = pd.DataFrame({"sample_id": s}, index=[0])
    else:
        metrics = com.build_metrics.to_frame().T
        metrics["sample_id"] = s
    return metrics


def build(
    taxonomy,
    model_db,
    out_folder,
    cutoff=0.0001,
    threads=1,
    solver=None,
):
    """Build a series of community models.

    This is a best-practice implementation of building community models
    for several samples in parallel.

    Parameters
    ----------
    taxonomy : pandas.DataFrame
        The taxonomy used for building the model. Must have at least the
        columns "id" and "sample_id". This must also
        contain at least a column with the same name as the rank used in
        the model database. Thus, for a genus-level database you will need
        a column `genus`. Additional taxa ranks can also be specified and
        will be used to be more stringent in taxa matching.
        Finally, the taxonomy should contain a column `abundance`. It will
        be used to quantify each individual in the community. If absent,
        MICOM will assume all individuals are present in the same amount.
    model_db : str
        A pre-built model database. If ending in `.qza` must be a Qiime 2
        artifact of type `MetabolicModels[JSON]`. Can also be a folder,
        zip (must end in `.zip`) file or None if the taxonomy contains a
        column `file`.
    out_folder : str
        The built models and a manifest file will be written to this
        folder.
    cutoff : float in [0.0, 1.0]
        Abundance cutoff. Taxa with a relative abundance smaller than this
        will not be included in the model.
    threads : int >=1
        The number of parallel workers to use when building models. As a
        rule of thumb you will need around 1GB of RAM for each thread.
    solver : str
        Name of the solver used for the linear and quadratic problems.

    Returns
    -------
    pandas.DataFrame
        The manifest for the built models. Contains taxa abundances,
        build metrics and file basenames.

    """
    os.makedirs(out_folder, exist_ok=True)
    samples = taxonomy.sample_id.unique()
    out_path = pd.Series({s: os.path.join(out_folder, s + ".pickle") for s in samples})
    args = [
        [s, taxonomy[taxonomy.sample_id == s], model_db, out_path[s], cutoff, solver]
        for s in samples
    ]
    res = workflow(build_and_save, args, threads)
    metrics = pd.concat(res)
    taxonomy = (
        taxonomy.groupby("sample_id")
        .apply(_reduce_group)
        .dropna(axis=1)
        .reset_index(drop=True)
    )
    taxonomy = taxonomy.loc[:, ~taxonomy.columns.isin(_ranks)]
    taxonomy["file"] = taxonomy.sample_id + ".pickle"
    taxonomy = pd.merge(taxonomy, metrics, on="sample_id")
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
    threads : int >=1
        The number of parallel workers to use when building models. As a
        rule of thumb you will need around 1GB of RAM for each thread.
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
        os.makedirs(out_path, exist_ok=True)
        args = [
            (tid, row, os.path.join(out_path, "%s.json" % tid))
            for tid, row in meta.iterrows()
        ]
        workflow(_summarize_models, args, threads)
        meta.file = meta.index + ".json"
        meta.to_csv(os.path.join(out_path, "manifest.csv"), index=False)

    return meta
