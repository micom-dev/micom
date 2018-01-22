"""Implements tradeoff optimization between community and egoistic growth."""

from micom.util import (_format_min_growth, _apply_min_growth,
                        check_modification)
from micom.logger import logger
from micom.solution import solve
from optlang.symbolics import Zero
from collections import Sized
from cobra.util import get_context
import pandas as pd


def reset_min_community_growth(com):
    """Reset the lower bound for the community growth."""
    com.variables.community_objective.lb = 0.0


def add_tradeoff(community, tradeoff, linear=False,
                 max_gcs=None):
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
    logger.info("adding egoistic objective to %s" % community.id)
    if max_gcs is None:
        max_gcs = community.optimize_all(progress=False)
    egoistic_expr = Zero
    for sp in community.species:
        species_obj = community.constraints["objective_" + sp]
        ex = max_gcs[sp] - species_obj.expression
        if not linear:
            ex = (ex**2).expand()
        egoistic_expr += ex
    community.objective = (tradeoff * community.variables.community_objective -
                           (1.0 - tradeoff) * egoistic_expr)
    if linear:
        community.modification = "linear cooperative tradeoff"
    else:
        community.modification = "quadratic cooperative tradeoff"
    logger.info("finished adding tradeoff objective to %s" % community.id)


def cooperative_tradeoff(community, linear, min_growth, tradeoff, fluxes,
                         pfba):
    """Find the best tradeoff between community and individual growth.

    Finds the set of growth rates which is closest to the community members
    individual maximal growth rates and still yields a (sub-)optimal community
    objective.

    Parameters
    ----------
    community : micom.Community
        The community to optimize.
    linear : boolean
        Whether to use a non-linear (sum of squares) or linear version of the
        cooperativity cost. If set to False requires a QP-capable solver.
    min_growth : float or array-like
        The minimal growth rate required for each individual. May be a
        single value or an array-like object with the same length as there
        are individuals.
    tradeoff : float or list of floats in [0, 1]
        How much the optimization should concentrate on community growth.
        0 would mean only optimize individual growth rates and 1 only optimize
        the community growth rate.
    fluxes : boolean
        Whether to return the fluxes as well.
    pfba : boolean
        Whether to obtain fluxes by parsimonious FBA rather than
        "classical" FBA. This is highly recommended.

    Returns
    -------
    micom.CommunitySolution
        The solution of the optimization. If fluxes==False will only contain
        the objective value, community growth rate and individual growth rates.

    """
    with community as com:
        check_modification(community)
        min_growth = _format_min_growth(min_growth, community.species)
        _apply_min_growth(community, min_growth)

        if not isinstance(tradeoff, Sized):
            tradeoff = [tradeoff]

        max_gcs = community.optimize_all(progress=False)
        # Add needed variables etc.
        results = []
        for to in tradeoff:
            with com:
                add_tradeoff(community, to, linear, max_gcs=max_gcs)
                results.append((to, solve(community, fluxes=fluxes,
                                          pfba=pfba)))
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


def knockout_species(community, species, linear, tradeoff, method):
    """Knockout a species from the community."""
    with community as com:
        check_modification(com)
        min_growth = _format_min_growth(0.0, com.species)
        _apply_min_growth(com, min_growth)

        max_gcs = community.optimize_all(progress=False)
        add_tradeoff(com, tradeoff, linear, max_gcs)
        old = com.optimize().members["growth_rate"]
        results = []

        for sp in com.species:
            with com:
                logger.info("getting egoistic tradeoff growth rates for "
                            "%s knockout" % sp)
                [r.knock_out() for r in
                 com.reactions.query(lambda ri: __is_needed(ri, sp))]
                # this should not be necessary but leaving the fixed zero
                # variables in the objective sometimes leads to problems
                # for some solver (cplex for instance)
                zero_growth(com, sp)
                sol = com.optimize()
                new = sol.members["growth_rate"]
                if "change" in method:
                    new = new - old
                if "relative" in method:
                    new /= old
                results.append(new)

        return pd.DataFrame(results, index=species).drop("medium", 1)
