"""Performs growth and exchange analysis for several models."""

from cobra.util.solver import interface_to_str, OptimizationError
from ..util import load_pickle
from ..annotation import annotate_metabolites_from_exchanges
from ..constants import DIRECTION, ARGS
from ..media import minimal_medium
from os import path
import pandas as pd
import warnings
import logging

logger = logging.getLogger(__name__)

warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)

DIRECTION = pd.Series(["import", "export"], index=[0, 1])
ARGS = {
    "none": {"fluxes": True, "pfba": False},
    "minimal imports": {"fluxes": False, "pfba": False},
    "pFBA": {"fluxes": True, "pfba": True},
}


def _growth(args):
    p, tradeoff, medium, conf = args
    com = load_pickle(p)

    if "glpk" in interface_to_str(com.solver.interface):
        logger.error(
            "Community models were not built with a QP-capable solver. "
            "This means that you did not install CPLEX or Gurobi. "
            "If you did install one of the two please file a bug report "
            "at https://github.com/micom-dev/micom/issues."
        )
        return None

    ex_ids = [r.id for r in com.exchanges]
    logger.info(
        "%d/%d import reactions found in model.",
        medium.index.isin(ex_ids).sum(),
        len(medium),
    )
    com.medium = medium[medium.index.isin(ex_ids)]

    # Add coupling or resource constraints if requested
    if conf.coupling.enabled:
        com.add_coupling_constraints(
            ids=None,
            strategy=conf.coupling.strategy,
            include_exchanges=conf.coupling.include_exchanges,
            constraint=conf.coupling.constraint,
            lower=conf.coupling.lower,
        )

    # Get growth rates
    args = ARGS[conf.simulation.flux_method].copy()
    args["atol"] = conf.tolerance
    args["rtol"] = conf.tolerance
    args["fraction"] = tradeoff
    try:
        sol = com.cooperative_tradeoff(**args)
        rates = sol.members
        rates["taxon"] = rates.index
        rates["tradeoff"] = tradeoff
        rates["sample_id"] = com.id
    except Exception as e:
        logger.error(
            "Could not solve cooperative tradeoff for %s. "
            "This can often be fixed by enabling `presolve`, choosing more "
            "permissive atol and rtol arguments, or by checking that medium "
            "fluxes are > atol.\nAdditional context: %s" % (com.id, e)
        )
        return None

    if conf.simulation.flux_method == "minimal imports":
        # Get the minimal medium and the solution at the same time
        med = minimal_medium(
            com,
            exchanges=None,
            community_growth=sol.growth_rate,
            min_growth=rates.growth_rate.drop("medium"),
            solution=True,
            weights=conf.media.weights,
            atol=conf.tolerance,
            rtol=conf.tolerance,
        )
        if med is None:
            logger.error(
                "The minimal medium optimization failed for %s. "
                "This can often be fixed by enabling `presolve`, choosing more "
                "permissive atol and rtol arguments, or by checking that medium "
                "fluxes are > atol." % com.id
            )
            return None
        sol = med["solution"]

    exs = list({r.global_id for r in com.internal_exchanges + com.exchanges})
    fluxes = sol.fluxes.loc[:, exs].copy()
    fluxes["sample_id"] = com.id
    fluxes["tolerance"] = conf.tolerance
    anns = annotate_metabolites_from_exchanges(com)
    return {"growth": rates, "exchanges": fluxes, "annotations": anns}
