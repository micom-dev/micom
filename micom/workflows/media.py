"""Example workflows for micom."""

from os import path
import pandas as pd
from micom import load_pickle
from micom.workflows.core import workflow
from micom.workflows.results import GrowthResults, combine_results
import micom.media as mm
from micom.logger import logger
from micom.solution import OptimizationError

DIRECTION = pd.Series(["import", "export"], index=[0, 1])


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
    elif not all(s in medium.sample_id.unique() for s in samples):
        missing = [s for s in samples if s not in medium.sample_id.unique()]
        raise ValueError(
            f"The medium is missing samples from the manifest: {', '.join(missing)}."
        )
    return medium.drop_duplicates(subset=["reaction", "sample_id"])


def _medium(args):
    """Get minimal medium for a single model."""
    s, p, com_growth, growth, mc, weights, solution = args
    com = load_pickle(p)

    tol = com.solver.configuration.tolerances.feasibility

    res = mm.minimal_medium(
        com,
        community_growth=com_growth,
        min_growth=growth,
        minimize_components=mc,
        open_exchanges=True,
        solution=solution,
        weights=weights,
        atol=tol,
        rtol=tol,
    )
    if res is None:
        logger.info("Could not get a minimal medium for sample %s." % s)
        return None
    result = dict()
    if solution:
        medium = res["medium"].to_frame()
        result["growth"] = GrowthResults.from_solution(res["solution"], com)
    else:
        medium = res.to_frame()
    medium.columns = ["flux"]
    medium["sample_id"] = s
    medium.index.name = "reaction"
    result["medium"] = medium.reset_index()
    return result


def minimal_media(
    manifest: pd.DataFrame,
    model_folder: str,
    community_growth: float = 0.0,
    growth: float = 0.1,
    minimize_components: bool = False,
    weights: str = None,
    summarize: bool = True,
    solution: bool = False,
    threads: int = 1,
) -> pd.DataFrame:
    """Calculate the minimal medium for a set of community models.

    This requires specification of either the minimal community growth rate,
    a minimal taxon growth rate that has to be reachable by all taxa in the sample
    simultaneously, or a combination of both. All imports will be opened and the
    minimal medium allowing those growth rates will be returned. What exactly is being
    minimized (mass flux, carbon flux, number of components) can be specified through
    the `weights` and `minimize_components` options.

    Note
    ----
    A common usage example would be to request some realistic growth rate for the entire
    community and a very low growth rate for all taxa to ensure they are growing ("alive")
    in the medium. The returned solution comes from the medium minimization problem and
    does not have to correspond to the cooperative tradeoff solution with the same medium.

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
    growth : positive float, dict, or pd.Series
        The taxon-specific growth rates that have to be achieved. If a single float gives
        the growth rate for each individual taxon. If a dict or Series gives the growth
        rate for each taxon specified that way. Here keys are the IDs for the taxon.
    minimize_components : boolean
        Whether to minimize the number of media components rather than the
        total flux. This will ignore the weight argument and might be very slow.
    weights : str
        Will scale the fluxes by a weight factor. Can either be "mass" which will
        scale by molecular mass, a single element which will scale by
        the elemental content (for instance "C" to scale by carbon content).
        If None every metabolite will receive the same weight.
        Will be ignored if `minimize_components` is True.
    summarize: boolean
        Whether to summarize the medium across all samples. If False will
        return a medium for each sample.
    threads: int
        The number of processes to use.

    Returns
    -------
    pandas.DataFrame or tuple of pandas.DataFrame and GrowthResult
        Either the medium or, if `solution=True` a tuple of the medium and the
        growth results.
    """
    samples = manifest.sample_id.unique()
    args = [
        (
            s,
            path.join(model_folder, manifest[manifest.sample_id == s].file.iloc[0]),
            community_growth,
            growth,
            minimize_components,
            weights,
            solution,
        )
        for s in samples
    ]
    results = workflow(_medium, args, threads)
    if all(r is None for r in results):
        raise OptimizationError(
            "Could not find a growth medium that allows the specified "
            "growth rate for any sample :("
        )
    elif any(r is None for r in results):
        logger.error(
            "For some samples I could not find a medium that fulfills "
            "the growth rate requirements. Returning media only for the "
            "succesful samples."
        )
    medium = pd.concat(r["medium"] for r in results if r is not None)
    if summarize:
        medium = medium.groupby("reaction").flux.max().reset_index()
    medium["metabolite"] = medium.reaction.str.replace("EX_", "")

    if solution:
        results = combine_results(r["growth"] for r in results if r is not None)
        return medium, results

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
            weights=weights,
        )
    except Exception:
        logger.error("Can't reach the specified growth rates for model %s." % sid)
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


def complete_community_medium(
    manifest: pd.DataFrame,
    model_folder: str,
    medium: pd.DataFrame,
    community_growth: float = 0.1,
    min_growth: float = 0.001,
    max_import: float = 1,
    minimize_components: float = False,
    summarize: bool = True,
    weights: str = None,
    threads: int = 1,
) -> pd.DataFrame:
    """Augment a growth medium so a community or specific taxa can grow on it.

    Note
    ----
    This will complete a growth medium for a single community/sample. For building
    growth media that work for arbitrary samples/compositions of taxa see
    `complete_db_medium` In contrast to `complete_db_medium` this will account for
    taxon-taxon interactions. However, growth rates will no longer be an emergent
    property of the simulation, because one needs to specify the community growth rate
    or growth rates for individual taxa.

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
        s: path.join(model_folder, manifest[manifest.sample_id == s].file.iloc[0])
        for s in samples
    }
    medium = process_medium(medium, samples)
    if medium.flux[medium.flux < 1e-6].any():
        medium.loc[medium.flux < 1e-6, "flux"] = 1e-6
        logger.info("Some import rates were to small and were adjusted to 1e-6.")
    args = [
        [
            s,
            p,
            community_growth,
            min_growth,
            max_import,
            minimize_components,
            medium.flux[medium.sample_id == s],
            weights,
        ]
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
