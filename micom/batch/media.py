"""Example workflows for micom."""

import pandas as pd
from ..util import load_pickle
from .results import GrowthResults
from ..media import minimal_medium, complete_medium
import logging

logger = logging.getLogger(__name__)

DIRECTION = pd.Series(["import", "export"], index=[0, 1])


def process_medium(medium, samples):
    """Prepare a medium for simulation.

    This will set the index to the reaction ID and convert it to a sample-wise medium if it is not already.

    Parameters
    ----------
    medium : pd.DataFrame
        The medium to process. Should have columns `reaction`, `flux`, and optionally `sample_id`.
    samples : list of str
        The sample IDs to include in the medium.

    Returns
    -------
    pd.DataFrame
        The processed medium with columns `reaction`, `flux`, and `sample_id`.

    """
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

    res = minimal_medium(
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


def _fix_medium(args):
    """Get the fixed medium for a model."""
    sid, p, growth, min_growth, max_import, mip, medium, weights = args
    com = load_pickle(p)
    try:
        fixed = complete_medium(
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
