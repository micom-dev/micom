"""Performs growth and exchange analysis for several models."""

from cobra.util.solver import interface_to_str, OptimizationError
from collections import namedtuple
from micom import load_pickle
from micom.logger import logger
from micom.media import minimal_medium
from micom.workflows.core import workflow
from micom.workflows.media import process_medium
import pandas as pd

DIRECTION = pd.Series(["import", "export"], index=[0, 1])
GrowthResults = namedtuple("GrowthResults", ["growth_rates", "exchanges"])


def _growth(args):
    p, tradeoff, medium = args
    com = load_pickle(p)

    if "glpk" in interface_to_str(com.solver.interface):
        logger.error(
            "Community models were not built with a QP-capable solver. "
            "This means that you did not install CPLEX or Gurobi. "
            "If you did install one of the two please file a bug report "
            "at https://github.com/micom-dev/q2-micom/issues."
        )
        return None

    ex_ids = [r.id for r in com.exchanges]
    logger.info(
        "%d/%d import reactions found in model.",
        medium.index.isin(ex_ids).sum(),
        len(medium),
    )
    com.medium = medium[medium.index.isin(ex_ids)]

    # Get growth rates
    try:
        sol = com.cooperative_tradeoff(fraction=tradeoff)
        rates = sol.members
        rates["taxon"] = rates.index
        rates["tradeoff"] = tradeoff
        rates["sample_id"] = com.id
    except Exception:
        logger.warning("Could not solve cooperative tradeoff for %s." % com.id)
        return None

    # Get the minimal medium
    med = minimal_medium(com, 0.95 * sol.growth_rate,
                         0.95 * rates.growth_rate.drop("medium"))
    # Apply medium and reoptimize
    com.medium = med[med > 0]
    sol = com.cooperative_tradeoff(fraction=tradeoff, fluxes=True, pfba=False)
    fluxes = sol.fluxes.loc[:, sol.fluxes.columns.str.startswith("EX_")].copy()
    fluxes["sample_id"] = com.id
    return {"growth": rates, "exchanges": fluxes}


def grow(
    manifest,
    model_folder,
    medium,
    tradeoff,
    threads,
):
    """Simulate growth for a set of community models."""
    samples = manifest.sample_id.unique()
    paths = {s: manifest[manifest.sample_id == s].file[0] for s in samples}
    medium = process_medium(medium, samples)
    args = [
        [p, tradeoff, medium.flux[medium.sample_id == s]]
        for s, p in paths.items()
    ]
    results = workflow(_growth, args, threads)
    if all([r is None for r in results]):
        raise OptimizationError(
            "All numerical optimizations failed. This indicates a problem "
            "with the solver or numerical instabilities. Check that you have "
            "CPLEX or Gurobi installed. You may also increase the abundance "
            "cutoff to create simpler models."
        )
    growth = pd.concat(r["growth"] for r in results if r is not None)
    growth = growth[growth.taxon != "medium"]
    exchanges = pd.concat(r["exchanges"] for r in results)
    exchanges["taxon"] = exchanges.index
    exchanges = exchanges.melt(
        id_vars=["taxon", "sample_id"], var_name="reaction", value_name="flux"
    )
    abundance = growth[["taxon", "sample_id", "abundance"]]
    exchanges = pd.merge(exchanges, abundance,
                         on=["taxon", "sample_id"], how="outer")
    exchanges["metabolite"] = exchanges.reaction.str.replace("EX_", "")
    exchanges["direction"] = DIRECTION[
        (exchanges.flux > 0.0).astype(int)
    ].values

    return GrowthResults(growth, exchanges)
