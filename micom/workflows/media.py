"""Example workflows for micom."""

from cobra.util.solver import OptimizationError
from optlang.interface import OPTIMAL
from os import path
import pandas as pd
from micom import load_pickle
import micom.media as mm
from micom.util import load_model, clean_ids
from micom.workflows.core import workflow
from micom.media import complete_medium
from micom.logger import logger


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
    p, min_growth = args
    com = load_pickle(p)
    # open the bounds
    for ex in com.exchanges:
        ex.bounds = (-1000.0, 1000.0)
    try:
        medium = mm.minimal_medium(com, 0.0, min_growth=min_growth).to_frame()
    except Exception:
        return None
    medium.columns = ["flux"]
    medium.index.name = "reaction"
    return medium.reset_index()


def minimal_media(
    manifest, model_folder, min_growth=0.1, threads=1
):
    """Calculate the minimal medium for a set of community models."""
    samples = manifest.sample_id.unique()
    paths = [
        path.join(model_folder, manifest[manifest.sample_id == s].file.iloc[0])
        for s in samples
    ]
    args = [[p, min_growth] for p in paths]
    results = workflow(_medium, args, threads)
    if any(r is None for r in results):
        raise OptimizationError(
            "Could not find a growth medium that allows the specified "
            "growth rate for all taxa in all samples :("
        )
    results = pd.concat(results, axis=0)
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
        fixed = medium.copy()
    if model.solver.status != OPTIMAL:
        logger.warning(
            "Can't reach the specified growth rate for model %s." % mid
        )
        fixed = medium.copy()
    fixed.name = mid
    return fixed


def fix_community_medium(
    tax,
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
    medium : pandas.Series
        A growth medium with exchange reaction IDs as index and positive
        import fluxes as values.
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
    if medium[medium < 1e-6].any():
        medium[medium < 1e-6] = 1e-6
        logger.info(
            "Some import rates were to small and were adjusted to 1e-6."
        )
    args = [
        (row.id, row.file, medium, min_growth, max_import, minimize_components)
        for _, row in tax.iterrows()
    ]
    res = workflow(_fix_medium, args, n_jobs=n_jobs, unit="model(s)")
    return pd.concat(res, axis=1).max(axis=1)
