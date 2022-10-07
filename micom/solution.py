"""A community solution object."""

import numpy as np
import pandas as pd
from collections import Counter
from optlang.interface import (
    OPTIMAL,
    NUMERIC,
    FEASIBLE,
    SUBOPTIMAL,
    ITERATION_LIMIT,
)
from optlang.symbolics import Zero
from itertools import chain
from functools import partial
from cobra.exceptions import OptimizationError
from cobra.core import Solution
from cobra.util import interface_to_str, get_context
from micom.logger import logger
from micom.util import reset_min_community_growth, _apply_min_growth
from swiglpk import glp_adv_basis


good = [OPTIMAL, NUMERIC, FEASIBLE, SUBOPTIMAL, ITERATION_LIMIT]
"""Solver states that permit returning the solution."""


def _group_taxa(values, ids, taxa, what="reaction"):
    """Format a list of values by id and taxa."""
    df = pd.DataFrame({values.name: values, what: ids, "compartment": taxa})
    df = df.pivot(index="compartment", columns=what, values=values.name)
    df.name = values.name
    return df


class CommunitySolution(Solution):
    """An FBA solution for an entire community.

    Attributes
    ----------
    objective_value : float
        The (optimal) value for the objective function.
    members : pandas.Series
        Contains basic info about the individual compartments/members of the
        community such as id, abundance and growth rates. Will also include
        one row for the external medium (without abundance and growth rate).
    growth_rate : float
        The overall growth rate for the community normalized to 1 gDW.
    status : str
        The solver status related to the solution.
    strategy : str
        The optimization strategy used to obtain the solution (may be empty).
    fluxes : pandas.DataFrame
        Contains the reaction fluxes (primal values of variables) stratified
        by compartment. Columns denote individual fluxes and rows denote
        compartments: one for every taxon plus one for the external medium.
        Fluxes will be NA if the reaction does not exist in the organism.
    reduced_costs : pandas.Series
        Contains reaction reduced costs (dual values of variables) stratified
        by taxa. Columns denote individual fluxes and rows denote taxa.
        Reduced costs will be NA if the reaction does not exist in the
        organism.
    shadow_prices : pandas.Series
        Contains metabolite shadow prices (dual values of constraints)
        stratified by taxa. Columns denote individual metabolites and rows
        denote taxa. Shadow prices will be NA if the metabolite does not
        exist in the organism.

    """

    def __init__(self, community, slim=False, reactions=None, metabolites=None):
        """Get the solution from a community model."""
        if reactions is None:
            reactions = community.reactions
        if metabolites is None:
            metabolites = community.metabolites
        rids = np.array([(r.global_id, r.community_id) for r in reactions])
        mids = np.array([(m.global_id, m.community_id) for m in metabolites])
        if not slim:
            var_primals = community.solver.primal_values
            fluxes = pd.Series(
                [var_primals[r.id] - var_primals[r.reverse_id] for r in reactions],
                name="fluxes",
            )
            super(CommunitySolution, self).__init__(
                community.solver.objective.value,
                community.solver.status,
                _group_taxa(fluxes, rids[:, 0], rids[:, 1]),
                None,
                None,
            )
        else:
            super(CommunitySolution, self).__init__(
                community.solver.objective.value,
                community.solver.status,
                None,
                None,
                None,
            )
        gcs = pd.Series(dtype="float64")
        for sp in community.taxa:
            gcs[sp] = community.constraints["objective_" + sp].primal
        # Workaround for an optlang bug (PR #120)
        if interface_to_str(community.problem) == "gurobi":
            gcs = gcs.abs()
        self.strategy = community.modification
        self.members = pd.DataFrame(
            {
                "abundance": community.abundances,
                "growth_rate": gcs,
                "reactions": pd.Series(Counter(rids[:, 1])),
                "metabolites": pd.Series(Counter(mids[:, 1])),
            }
        )
        self.members.index.name = "compartments"
        self.growth_rate = sum(community.abundances * gcs)
        # Save estimated accuracy for osqp
        if interface_to_str(community.problem) == "osqp":
            self.primal_residual = community.solver.problem.info.pri_res
            self.dual_residual = community.solver.problem.info.dua_res

    def _repr_html_(self):
        if self.status in good:
            with pd.option_context("display.max_rows", 10):
                html = (
                    "<strong>community growth:</strong> {:.3f}"
                    "<br><strong>status:</strong> {}"
                    "<br><strong>taxa:</strong>{}".format(
                        sum(
                            self.members.abundance.dropna()
                            * self.members.growth_rate.dropna()
                        ),
                        self.status,
                        self.members._repr_html_(),
                    )
                )
        else:
            html = "<strong>{}</strong> solution :(".format(self.status)
        return html

    def __repr__(self):
        """Convert CommunitySolution instance to string representation."""
        if self.status not in good:
            return "<CommunitySolution {0:s} at 0x{1:x}>".format(self.status, id(self))
        return "<CommunitySolution {0:.3f} at 0x{1:x}>".format(
            self.growth_rate, id(self)
        )


def add_pfba_objective(community, atol=1e-6, rtol=1e-6):
    """Add pFBA objective.

    Add objective to minimize the summed flux of all reactions to the
    current objective. This one will work with any objective (even non-linear
    ones).

    See Also
    --------
    pfba

    Parameters
    ----------
    community : micom.Community
        The community to add the objective to.
    """
    # Fix all growth rates
    rates = {
        sp: community.constraints["objective_" + sp].primal for sp in community.taxa
    }
    _apply_min_growth(community, rates, atol, rtol)

    if community.solver.objective.name == "_pfba_objective":
        raise ValueError("model already has pfba objective")
    reaction_variables = (
        (rxn.forward_variable, rxn.reverse_variable) for rxn in community.reactions
    )
    variables = chain(*reaction_variables)
    community.objective = Zero
    community.objective_direction = "min"
    community.objective.set_linear_coefficients(dict.fromkeys(variables, 1.0))
    if interface_to_str(community.solver.interface) == "osqp":
        community.objective += 1e-6 * community.variables.community_objective ** 2
    if community.modification is None:
        community.modification = "pFBA"
    else:
        community.modification += " and pFBA"


def solve(community, fluxes=True, pfba=True, raise_error=False, atol=1e-6, rtol=1e-6):
    """Get all fluxes stratified by taxa."""
    solver_name = interface_to_str(community.solver.interface)
    term = None
    if solver_name == "osqp" and community.objective.is_Linear:
        # This improves OSQP by soooo much
        term = (
            1e-6
            * community.solver.problem.direction
            * community.variables.community_objective ** 2
        )
        community.objective += term
    community.solver.optimize()
    status = community.solver.status
    if status in good:
        if status != OPTIMAL:
            if raise_error:
                raise OptimizationError("solver returned the status %s." % status)
            else:
                logger.info(
                    "solver returned the status %s," % status
                    + " returning the solution anyway."
                )
        if pfba and fluxes:
            add_pfba_objective(community, atol, rtol)
            community.solver.optimize()
        if fluxes:
            sol = CommunitySolution(community)
        else:
            sol = CommunitySolution(community, slim=True)
        if term:
            correction = 1e-6 * community.variables.community_objective.primal ** 2
            sol.objective_value -= community.solver.problem.direction * correction
            community.objective -= term
        return sol
    logger.warning("solver encountered an error %s" % status)
    return None


def reset_solver(community):
    """Reset the solver."""
    interface = interface_to_str(community.solver.interface)
    logger.info("resetting solver, hoping for the best.")
    if interface == "cplex":
        community.solver.configuration.lp_method = "network"
        community.solver.configuration.lp_method = "barrier"
    elif interface == "gurobi":
        community.solver.problem.reset()
    elif interface == "glpk":
        glp_adv_basis(community.solver.problem, 0)
    elif interface == "osqp":
        community.solver.problem.reset()


def optimize_with_retry(com, message="could not get optimum."):
    """Try to reset the solver."""
    sol = com.optimize()
    if sol is None:
        reset_solver(com)
        sol = com.optimize()
    if sol is None:
        raise OptimizationError(message)
    else:
        return sol.objective_value


def crossover(community, sol, fluxes=False):
    """Get the crossover solution."""
    gcs = sol.members.growth_rate.drop("medium")
    com_growth = sol.growth_rate
    logger.info("Starting crossover...")
    with community as com:
        logger.info("constraining growth rates.")
        context = get_context(community)
        if context is not None:
            context(partial(reset_min_community_growth, com))
        reset_min_community_growth(com)
        com.variables.community_objective.lb = 0.0
        com.variables.community_objective.ub = com_growth + 1e-6
        com.objective = com.scale * com.variables.community_objective
        for sp in com.taxa:
            const = com.constraints["objective_" + sp]
            const.ub = gcs[sp]
        logger.info("finding closest feasible solution")
        s = com.optimize()
        if s is None:
            reset_solver(com)
            s = com.optimize()
        if s is not None:
            s = CommunitySolution(com, slim=not fluxes)
        for sp in com.taxa:
            com.constraints["objective_" + sp].ub = None
    if s is None:
        raise OptimizationError(
            "crossover could not converge (status = %s)." % community.solver.status
        )
    s.objective_value /= com.scale
    return s


def optimize_with_fraction(com, fraction, growth_rate=None, fluxes=False):
    """Optimize with a constrained community growth rate."""
    com.variables.community_objective.lb = 0
    com.variables.community_objective.ub = None
    if growth_rate is None:
        with com:
            com.objective = com.scale * com.variables.community_objective
            growth_rate = optimize_with_retry(
                com, message="could not get community growth rate."
            )
            growth_rate /= com.scale
    com.variables.community_objective.lb = fraction * growth_rate
    com.variables.community_objective.ub = growth_rate
    sol = com.optimize()
    if sol.status != OPTIMAL:
        sol = crossover(com, sol, fluxes=fluxes)
    else:
        sol = CommunitySolution(com, slim=not fluxes)
    return sol
