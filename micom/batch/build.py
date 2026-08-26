"""Worflow to build models for several samples."""

from cobra.io import save_json_model
from ..util import join_models, load_pickle, _read_model
from ..community import Community
from .core import workflow
import logging
from pathlib import Path
import pandas as pd
from tempfile import TemporaryDirectory
import zipfile

logger = logging.getLogger(__name__)


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
