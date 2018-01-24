"""Implements tradeoff optimization between community and egoistic growth."""

from micom.util import (_format_min_growth, _apply_min_growth,
                        check_modification, get_context)
from micom.logger import logger
from micom.solution import solve
from optlang.symbolics import Zero
from collections import Sized
from functools import partial
import pandas as pd
from tqdm import tqdm


def reset_min_community_growth(com):
    """Reset the lower bound for the community growth."""
    com.variables.community_objective.lb = 0.0


def regularize_l2_norm(community, min_growth):
    """Add an objective to find the most "egoistic" solution.

    This adds an optimization objective finding a solution that maintains a
    (sub-)optimal community growth rate but is the closest solution to the
    community members individual maximal growth rates. So it basically finds
    the best possible tradeoff between maximizing community growth and
    individual (egoistic) growth. Here the objective is given as the sum of
    squared differences between the individuals current and maximal growth
    rate. In the linear case squares are substituted by absolute values
    (Manhattan distance).

    Arguments
    ---------
    community : micom.Community
        The community to modify.
    min_community_growth : positive float
        The minimal community growth rate that has to be mantained.
    linear : boolean
        Whether to use a non-linear (sum of squares) or linear version of the
        cooperativity cost. If set to False requires a QP-capable solver.
    max_gcs : None or dict
        The precomputed maximum individual growth rates.

    """
    logger.info("adding L2 norm to %s" % community.id)
    l2 = Zero
    community.variables.community_objective.lb = min_growth
    context = get_context(community)
    if context is not None:
        context(partial(reset_min_community_growth, community))

    for sp in community.species:
        species_obj = community.constraints["objective_" + sp]
        ex = sum(v for v in species_obj.variables if (v.ub - v.lb) > 1e-6)
        l2 += (ex**2).expand()
    community.objective = -l2
    community.modification = "l2 norm"
    logger.info("finished adding tradeoff objective to %s" % community.id)


def cooperative_tradeoff(community, min_growth, fraction, fluxes, pfba):
    """Find the best tradeoff between community and individual growth."""
    with community as com:
        check_modification(community)
        min_growth = _format_min_growth(min_growth, community.species)
        _apply_min_growth(community, min_growth)

        com.objetive = 1.0 * com.variables.community_objective
        min_growth = com.slim_optimize()

        if not isinstance(tradeoff, Sized):
            tradeoff = [tradeoff]

        # Add needed variables etc.
        regularize_l2_norm(com, 0.0)
        results = []
        for to in tradeoff:
            com.variables.community_objective.lb = to * min_growth
            results.append((to, solve(community, fluxes=fluxes, pfba=pfba)))
        if len(results) == 1:
            return results[0][1]
        return pd.DataFrame.from_records(results,
                                         columns=["tradeoff", "solution"])


def __is_needed(r, s):
    """Find the smallest numbers of reactions that knock-out the individual."""
    return (r.community_id == s) & \
           ((r.bounds[0] * r.bounds[1] > 0.0) |
            ("m" in {m.compartment for m in r.metabolites}))


def zero_growth(com, s):
    """Force zero_growth for a given species."""
    const = com.constraints["objective_" + s]
    bounds = (const.lb, const.ub)

    def reset():
        const.lb = bounds[0]
        const.ub = bounds[1]

    const.lb = const.ub = 0.0
    context = get_context(com)
    if context:
        context(reset)


def knockout_species(community, species, fraction, method, progress):
    """Knockout a species from the community."""
    with community as com:
        check_modification(com)
        min_growth = _format_min_growth(0.0, com.species)
        _apply_min_growth(com, min_growth)

        min_growth = com.slim_optimize()
        regularize_l2_norm(com, fraction * min_growth)
        old = com.optimize().members["growth_rate"]
        results = []

        if progress:
            species = tqdm(species, unit="knockout(s)")
        for sp in species:
            with com:
                logger.info("getting egoistic tradeoff growth rates for "
                            "%s knockout" % sp)
                com.variables.community_objective.lb = 0.0
                [r.knock_out() for r in
                 com.reactions.query(lambda ri: __is_needed(ri, sp))]
                # this should not be necessary but leaving the fixed zero
                # variables in the objective sometimes leads to problems
                # for some solver (cplex for instance)
                zero_growth(com, sp)
                with com:
                    com.objective = 1.0 * com.variables.community_objective
                    min_growth = com.slim_optimize()
                com.variables.community_objective.lb = fraction * min_growth
                sol = com.optimize()
                new = sol.members["growth_rate"]
                if "change" in method:
                    new = new - old
                if "relative" in method:
                    new /= old
                results.append(new)

        return pd.DataFrame(results, index=species).drop("medium", 1)
