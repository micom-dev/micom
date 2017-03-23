"""Implements optimization and model problems."""

from micom.util import fluxes_from_primals
import pandas as pd
from cobra.flux_analysis.parsimonious import add_pfba
from micom.duality import fast_dual
from sympy.core.singleton import S
from functools import partial


def all_fluxes(community, pfba=False):
    """Get all fluxes stratified by species."""
    community.solver.optimize()
    if community.solver.status == "optimal":
        if pfba:
            add_pfba(community)
            community.solver.optimize()
        fluxes = (fluxes_from_primals(community, row)
                  for _, row in community.taxonomy.iterrows())
        fluxes = pd.concat(fluxes, axis=1).T
        return fluxes
    return None


def check_modification(community):
    """Check whether a community already carries a modification."""
    if community.modification is not None:
        raise ValueError("Community already carries a modification "
                         "({})!".format(community.modification))


def _format_min_growth(min_growth, species):
    """Format min_growth into a pandas series."""
    try:
        min_growth = float(min_growth)
    except (TypeError, ValueError):
        if len(min_growth) != len(species):
            raise ValueError(
                "min_growth must be single value or an array-like "
                "object with an entry for each species in the model.")
    return pd.Series(min_growth, species)


def add_linear_optcom(community, min_growth=0.1):
    """Add a linear version of optcom."""
    check_modification(community)
    species = list(community.objectives.keys())
    min_growth = _format_min_growth(min_growth, species)

    prob = community.solver.interface
    to_add = []
    for sp in species:
        obj = prob.Constraint(community.objectives[sp],
                              name="objective_" + sp,
                              lb=min_growth[sp])
        to_add.append(obj)

    community.add_cons_vars(to_add)
    community.modification = "linear optcom"


def add_lagrangian(community, tradeoff, linear=False):
    """Adds a Lagrangian optimization target to a linear OptCom model."""
    species = list(community.objectives.keys())
    max_gcs = community.optimize_all()
    prob = community.problem
    com_expr = S.Zero
    cost_expr = S.Zero
    abundances = community.abundances
    for sp in species:
        com_expr += abundances[sp] * community.objectives[sp]
        v_max_gc = prob.Variable("gc_constant_" + sp, lb=max_gcs[sp],
                                 ub=max_gcs[sp])
        community.add_cons_vars([v_max_gc])
        ex = v_max_gc - community.objectives[sp]
        if not linear:
            ex = ex**2
        cost_expr += ex.expand()
    community.objective = (S.One - tradeoff) * com_expr - tradeoff * cost_expr
    if "with lagrangian" not in community.modification:
        community.modification += " with lagrangian"


def add_dualized_optcom(community, min_growth):
    """Add dual Optcom variables and constraints to a community."""
    check_modification(community)
    species = list(community.objectives.keys())
    min_growth = _format_min_growth(min_growth, species)

    prob = community.solver.interface

    # Temporarily subtitute objective with sum of individual objectives
    # for correct dual variables
    old_obj = community.objective
    community.objective = S.Zero
    for expr in community.objectives.values():
        community.objective += expr

    to_add = []
    for sp in species:
        obj = prob.Constraint(community.objectives[sp],
                              name="objective_" + sp,
                              lb=min_growth[sp])
        to_add.append(obj)
    community.add_cons_vars(to_add)
    dual_coefs = fast_dual(community)

    for sp in species:
        primal_const = community.constraints["objective_" + sp]
        coefs = primal_const.get_linear_coefficients(primal_const.variables)
        coefs.update({dual_var: -coef for dual_var, coef in
                      dual_coefs.items() if sp in dual_var.name})
        obj_constraint = prob.Constraint(
            S.Zero, lb=0, ub=0, name="optcom_suboptimality_" + sp)
        community.add_cons_vars([obj_constraint])
        community.solver.update()
        obj_constraint.set_linear_coefficients(coefs)

    community.objective = old_obj
    community.modification = "dual optcom"


def add_moma_optcom(community, min_growth, linear=False):
    """Add a dualized MOMA version of OptCom."""
    check_modification(community)
    species = list(community.objectives.keys())
    min_growth = _format_min_growth(min_growth, species)

    prob = community.solver.interface
    old_obj = community.objective
    coefs = old_obj.get_linear_coefficients(old_obj.variables)

    # Get maximum individual growth rates
    max_gcs = community.optimize_all()

    to_add = []
    for sp in species:
        obj = prob.Constraint(community.objectives[sp],
                              name="objective_" + sp,
                              lb=min_growth[sp])
        to_add.append(obj)
    community.add_cons_vars(to_add)
    dual_coefs = fast_dual(community)
    coefs.update({
        v: -coef for v, coef in
        dual_coefs.items()})
    obj_constraint = prob.Constraint(
        S.Zero, lb=0, ub=0,
        name="optcom_suboptimality_" + sp)
    community.add_cons_vars([obj_constraint])
    community.solver.update()
    obj_constraint.set_linear_coefficients(coefs)
    obj_expr = S.Zero
    for sp in species:
        v = prob.Variable("gc_constant_" + sp, lb=max_gcs[sp], ub=max_gcs[sp])
        community.add_cons_vars([v])
        ex = v - community.objectives[sp]
        if not linear:
            ex = ex**2
        obj_expr += ex.expand()
    community.objective = prob.Objective(obj_expr, direction="min")
    community.modification = "moma optcom"


_methods = {"linear": [add_linear_optcom],
            "lagrangian": [add_linear_optcom, add_lagrangian],
            "linear lagrangian": [add_linear_optcom,
                                  partial(add_lagrangian, linear=True)],
            "original": [add_dualized_optcom],
            "moma": [add_moma_optcom],
            "lmoma": [partial(add_moma_optcom, linear=True)]}


def optcom(community, strategy, min_growth, tradeoff, fluxes, pfba):
    """Obtain optimal growth rates with a dualized version of OptCom."""
    if strategy not in _methods:
        raise ValueError("strategy most be one of {}!".format(
                         ",".join(_methods)))
    species = list(community.objectives.keys())
    funcs = _methods[strategy]

    community_obj = None
    gcs = pd.Series(None, species)
    with community as com:
        # Add needed variables etc.
        funcs[0](com, min_growth)
        if "lagrangian" in strategy:
            funcs[1](com, tradeoff)
        com.solver.optimize()
        if com.solver.status == "optimal":
            community_obj = com.variables.community_objective.primal
            for sp in species:
                gcs[sp] = com.constraints["objective_" + sp].primal
            if fluxes:
                com_fluxes = all_fluxes(com, pfba=pfba)

    if fluxes:
        return community_obj, com_fluxes
    else:
        return community_obj, gcs
