"""Workflow to run cooperative tradeoff with various tradeoff values."""

from cobra.util.solver import OptimizationError
from micom import load_pickle
from micom.logger import logger
from micom.workflows.core import workflow
from micom.workflows.media import process_medium
import numpy as np
from os import path
import pandas as pd


def _tradeoff(args):
    p, tradeoffs, medium = args
    com = load_pickle(p)
    ex_ids = [r.id for r in com.exchanges]
    logger.info(
        "%d/%d import reactions found in model.",
        medium.index.isin(ex_ids).sum(),
        len(medium),
    )
    com.medium = medium[medium.index.isin(ex_ids)]
    try:
        sol = com.optimize()
    except Exception:
        logger.warning(
            "Sample %s could not be optimized (%s)." %
            (com.id, com.solver.status),
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
        logger.warning(
            "Sample %s could not be optimized (%s)." %
            (com.id, com.solver.status),
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


def tradeoff(
    manifest,
    model_folder,
    medium,
    tradeoffs=np.arange(0.1, 1.0 + 1e-6, 0.1),
    threads=1,
):
    """Run tradeoff analysis."""
    samples = manifest.sample_id.unique()
    paths = {
        s: path.join(model_folder, manifest[manifest.sample_id == s].file[0])
        for s in samples
    }
    if any(t < 0.0 or t > 1.0 for t in tradeoffs):
        raise ValueError(
            "tradeoff values must between 0 and 1 :("
        )
    medium = process_medium(medium, samples)
    args = [
        [p, tradeoffs, medium.flux[medium.sample_id == s]]
        for s, p in paths.items()
    ]
    results = workflow(_tradeoff, args, threads)
    if all(r is None for r in results):
        raise OptimizationError(
            "All numerical optimizations failed. This indicates a problem "
            "with the solver or numerical instabilities. Check that you have "
            "CPLEX or Gurobi installed. You may also increase the abundance "
            "cutoff in `qiime micom build` to create simpler models."
        )
    results = pd.concat(results)
    return results
