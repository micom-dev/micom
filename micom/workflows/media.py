"""Example workflows for micom."""

from os import path
import pandas as pd
from micom import load_pickle
from micom.db import load_manifest, load_zip_model_db
import micom.media as mm
from micom.util import load_model, clean_ids
from micom.workflows.core import workflow
from micom.media import complete_medium
from micom.logger import logger
from micom.qiime_formats import load_qiime_model_db
from micom.solution import OptimizationError
from tempfile import TemporaryDirectory


def process_medium(medium, samples):
    """Prepare a medium for simulation."""
    medium.index = medium.reaction
    if "sample_id" not in medium.columns:
        meds = []
        for s in samples:
            m = medium.copy()
            m["sample_id"] = s
            meds.append(m)
        medium = pd.concat(meds, axis=0)
    return medium


def _medium(args):
    """Get minimal medium for a single model."""
    s, p, min_growth = args
    com = load_pickle(p)
    # open the bounds
    for ex in com.exchanges:
        ex.bounds = (-1000.0, 1000.0)
    try:
        medium = mm.minimal_medium(com, 0.0, min_growth=min_growth).to_frame()
    except Exception:
        return None
    medium.columns = ["flux"]
    medium["sample_id"] = s
    medium.index.name = "reaction"
    return medium.reset_index()


def minimal_media(
    manifest, model_folder, summarize=True, min_growth=0.1, threads=1
):
    """Calculate the minimal medium for a set of community models."""
    samples = manifest.sample_id.unique()
    paths = [
        (
            s,
            path.join(
                model_folder, manifest[manifest.sample_id == s].file.iloc[0]
            ),
        )
        for s in samples
    ]
    args = [[s, p, min_growth] for s, p in paths]
    results = workflow(_medium, args, threads)
    if any(r is None for r in results):
        raise OptimizationError(
            "Could not find a growth medium that allows the specified "
            "growth rate for all taxa in all samples :("
        )
    results = pd.concat(results, axis=0)
    if summarize:
        medium = results.groupby("reaction").flux.max().reset_index()
    medium["metabolite"] = medium.reaction.str.replace("EX_", "")
    return medium


def _fix_medium(args):
    """Get the fixed medium for a model."""
    mid, file, medium, min_growth, max_import, min_c = args
    model = load_model(file)
    for r in model.reactions:
        r.id = clean_ids(r.id)
    try:
        fixed = complete_medium(
            model,
            medium,
            min_growth=min_growth,
            max_import=max_import,
            minimize_components=min_c,
        )
    except Exception:
        logger.warning(
            "Can't reach the specified growth rate for model %s." % mid
        )
        return None
    fixed = pd.DataFrame({"reaction": fixed.index, "flux": fixed.values})
    fixed["metabolite"] = [
        model.reactions.get_by_id(r).reactants[0].id for r in fixed.reaction
    ]
    fixed["description"] = [
        model.reactions.get_by_id(r).reactants[0].name for r in fixed.reaction
    ]
    return fixed


def fix_medium(
    model_db,
    medium,
    min_growth=0.1,
    max_import=1,
    minimize_components=True,
    n_jobs=4,
):
    """Augment a growth medium so all community members can grow in it.

    Arguments
    ---------
    tax : pandas.Dataframe
        A taxonomy specification as passed to `micom.Community`.
    medium : pandas.Series or pandas.DataFrame
        A growth medium with exchange reaction IDs as index and positive
        import fluxes as values. If a DataFrame needs columns `flux` and
        `reaction`.
    min_growth : positive float
        The minimum biomass production required for growth.
    max_import : positive float
        The maximum import rate for added imports.
    minimize_components : boolean
        Whether to minimize the number of media components rather than the
        total flux.
    n_jobs: int
        The number of processes to use.

    Returns
    -------
    pandas.Series
        A new growth medium with the smallest amount of augmentations such
        that all members of the community can grow in it.

    """
    if isinstance(medium, pd.DataFrame):
        medium = medium.copy()
        medium.index = medium.reaction
        medium = medium.flux
    elif not isinstance(medium, pd.Series):
        raise ValueError("`medium` must be a DataFrame or Series.")

    compressed = model_db.endswith(".qza") or model_db.endswith(".zip")
    if compressed:
        tdir = TemporaryDirectory(prefix="micom_")
    if model_db.endswith(".qza"):
        manifest = load_qiime_model_db(model_db, tdir.name)
    elif model_db.endswith(".zip"):
        manifest = load_zip_model_db(model_db, tdir.name)
    else:
        manifest = load_manifest(model_db)
    manifest["file"] = [path.join(tdir.name, f) for f in manifest.file]

    if medium[medium < 1e-6].any():
        medium[medium < 1e-6] = 1e-6
        logger.info(
            "Some import rates were to small and were adjusted to 1e-6."
        )
    args = [
        (row.id, row.file, medium, min_growth, max_import, minimize_components)
        for _, row in manifest.iterrows()
    ]
    res = workflow(_fix_medium, args, n_jobs=n_jobs, unit="model(s)")
    if all(r is None for r in res):
        raise OptimizationError(
            "All optimizations failed. You may need to increase `max_import` "
            "or lower the target growth rate."
        )
    res = pd.concat(res)
    final = (
        res.groupby(["reaction", "metabolite", "description"])
        .flux.max()
        .reset_index()
    )
    return final
