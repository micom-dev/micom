"""Example workflows for micom."""

from os import path
import pandas as pd
from micom import load_pickle
from micom.workflows.core import workflow
import micom.media as mm
from micom.logger import logger
from micom.solution import OptimizationError


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
    return medium.drop_duplicates(subset=["reaction", "sample_id"])


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
        logger.error("Could not get a minimal medium for sample %s." % s)
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
    sid, p, growth, min_growth, max_import, mip, medium, weights = args
    com = load_pickle(p)
    try:
        fixed = mm.complete_medium(
            com,
            medium,
            growth=growth,
            min_growth=min_growth,
            max_import=max_import,
            minimize_components=mip,
            weights=weights
        )
    except Exception:
        logger.error(
            "Can't reach the specified growth rates for model %s." % sid
        )
        return None
    fixed = pd.DataFrame({"reaction": fixed.index, "flux": fixed.values})
    fixed["metabolite"] = [
        list(com.reactions.get_by_id(r).metabolites.keys())[0].id
        for r in fixed.reaction
    ]
    fixed["description"] = [
        list(com.reactions.get_by_id(r).metabolites.keys())[0].name
        for r in fixed.reaction
    ]
    fixed["sample_id"] = sid
    return fixed


def fix_medium(
    manifest,
    model_folder,
    medium,
    community_growth=0.1,
    min_growth=0.001,
    max_import=1,
    minimize_components=False,
    summarize=True,
    weights=None,
    threads=1,
):
    """Augment a growth medium so all community members can grow in it.

    Arguments
    ---------
    manifest : pandas.DataFrame
        The manifest as returned by the `build` workflow.
    model_folder : str
        The folder in which to find the files mentioned in the manifest.
    medium : pandas.Series or pandas.DataFrame
        A growth medium with exchange reaction IDs as index and positive
        import fluxes as values. If a DataFrame needs columns `flux` and
        `reaction`.
    community_growth : positive float
        The minimum community-wide growth rate that has to be achieved on the created
        medium.
    min_growth : positive float
        The minimum biomass production required for growth.
    max_import : positive float
        The maximum import rate for added imports.
    minimize_components : boolean
        Whether to minimize the number of media components rather than the
        total flux.
    summarize: boolean
        Whether to summarize the medium across all samples. If False will
        return a medium for each sample.
    weights : str
        Will scale the fluxes by a weight factor. Can either be "mass" which will
        scale by molecular mass, a single element which will scale by
        the elemental content (for instance "C" to scale by carbon content).
        If None every metabolite will receive the same weight.
        Will be ignored if `minimize_components` is True.
    threads: int
        The number of processes to use.

    Returns
    -------
    pandas.DataFrame
        A new growth medium with the smallest amount of augmentations such
        that all members of the community can grow in it.

    """
    if not isinstance(medium, pd.DataFrame):
        raise ValueError("`medium` must be a DataFrame.")

    samples = manifest.sample_id.unique()
    paths = {
        s: path.join(
            model_folder, manifest[manifest.sample_id == s].file.iloc[0])
        for s in samples
    }
    medium = process_medium(medium, samples)
    if medium.flux[medium.flux < 1e-6].any():
        medium.loc[medium.flux < 1e-6, "flux"] = 1e-6
        logger.info(
            "Some import rates were to small and were adjusted to 1e-6."
        )
    args = [
        [s, p, community_growth, min_growth, max_import, minimize_components,
         medium.flux[medium.sample_id == s], weights]
        for s, p in paths.items()
    ]
    res = workflow(_fix_medium, args, threads=threads, description="Augmenting media")
    if all(r is None for r in res):
        raise OptimizationError(
            "All optimizations failed. You may need to increase `max_import` "
            "or lower the target growth rate."
        )
    final = pd.concat(res)
    if summarize:
        final = (
            final.groupby(["reaction", "metabolite", "description"])
            .flux.max()
            .reset_index()
        )
    return final
