"""Worflow to build models for several samples."""

from cobra.io import save_json_model
from glob import glob
from micom.logger import logger
from micom.util import join_models, load_pickle, _read_model
from micom.community import Community, _ranks
from micom.workflows.core import workflow
import os
import pandas as pd
from tempfile import TemporaryDirectory
import zipfile


def _reduce_group(df):
    keep = df.columns[df.nunique() == 1]
    new = df.iloc[0, :][keep]
    if "file" in df.columns:
        new["file"] = "|".join(df.file.astype(str))
    return pd.DataFrame.from_records([new])


def build_and_save(args):
    """Build a single community model."""
    s, tax, db, out, cutoff, solver = args
    if os.path.exists(out):
        com = load_pickle(out)
    else:
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
        folder. Will continue
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
    if os.path.exists(out_folder):
        existing = [
            s.split(".pickle")[0] for s in glob(os.path.join(out_folder, "*.pickle"))
        ]
        if len(existing) > 0:
            logger.warning(
                f"Found existing models for {len(existing)} samples. Will skip those. "
                "Delete the output folder if you would like me to rebuild them."
            )
    else:
        os.makedirs(out_folder)

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
        mod = _read_model(files[0])
    save_json_model(mod, new_path)


def build_database(
    manifest,
    out_path,
    rank="genus",
    threads=1,
    compress=None,
    compresslevel=6,
    progress=True,
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
        The directory or zip file where the joined models will be written.
    threads : int >=1
        The number of parallel workers to use when building models. As a
        rule of thumb you will need around 1GB of RAM for each thread.
    compress : str (default None)
        Compression method to use. Must be "zlib", "bz2", "lzma" or None.
        This parameter is ignored if out_path does not end with ".zip".
    compresslevel : int [1-9] (default: 6)
        Level of compression. Only used if compress is not None.
        This parameter is ignored if out_path does not end with ".zip".
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
    meta.index = meta[rank].str.replace("[^\\w\\_]", "_", regex=True)
    meta["id"] = meta.index
    meta["summary_rank"] = rank

    # compress is ignored if outpath does not end with ".zip"
    if out_path.endswith(".zip"):
        # Explicitly check compression level
        if compresslevel not in range(1, 10):
            raise ValueError("compresslevel parameter must be an int between 1 and 9")

        # Explicitly check for supported zipfile compression options
        compressdict = {
            None: zipfile.ZIP_STORED,
            "zlib": zipfile.ZIP_DEFLATED,
            "bz2": zipfile.ZIP_BZIP2,
            "lzma": zipfile.ZIP_LZMA,
        }
        if compress not in compressdict:
            raise ValueError('compress parameter must be "zlib", "bz2", "lzma" or None')
        compressopt = compressdict[compress]
        # Check if zipfile compression dependencies are installed
        # Raise RuntimeError if the module is missing
        zipfile._check_compression(compressopt)

        # Store model database as zipfile
        with TemporaryDirectory(prefix="micom_") as tdir:
            args = [
                (tid, row, os.path.join(tdir, "%s.json" % tid))
                for tid, row in meta.iterrows()
            ]
            workflow(_summarize_models, args, threads, progress=progress)
            meta.file = meta.index + ".json"
            meta.to_csv(os.path.join(tdir, "manifest.csv"), index=False)
            with zipfile.ZipFile(
                out_path,
                mode="w",
                compression=compressopt,
                compresslevel=compresslevel,
            ) as zf:
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
