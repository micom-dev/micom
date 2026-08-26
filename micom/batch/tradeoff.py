"""Workflow to run cooperative tradeoff with various tradeoff values."""

from ..util import load_pickle
import numpy as np
import pandas as pd
import logging

logger = logging.getLogger(__name__)


def _tradeoff(args):
    p, tradeoffs, medium, atol, rtol, presolve = args

    com = load_pickle(p)
    ex_ids = [r.id for r in com.exchanges]
    logger.info(
        "%d/%d import reactions found in model.",
        medium.index.isin(ex_ids).sum(),
        len(medium),
    )
    com.medium = medium[medium.index.isin(ex_ids)]
    com.solver.configuration.presolve = presolve

    try:
        sol = com.optimize(rtol=rtol, atol=atol)
    except Exception:
        logger.error(
            "Sample %s could not be optimized (%s)." % (com.id, com.solver.status),
        )
        return None
    rates = sol.members
    rates["taxon"] = rates.index
    rates["tradeoff"] = np.nan
    rates["sample_id"] = com.id
    df = [rates]

    # Get growth rates
    try:
        sol = com.cooperative_tradeoff(fraction=tradeoffs)
    except Exception:
        logger.info(
            "Sample %s could not be optimized with cooperative tradeoff (%s)."
            % (com.id, com.solver.status),
        )
        return None
    for i, s in enumerate(sol.solution):
        rates = s.members
        rates["taxon"] = rates.index
        rates["tradeoff"] = sol.tradeoff[i]
        rates["sample_id"] = com.id
        df.append(rates)
    df = pd.concat(df)
    return df[df.taxon != "medium"]
