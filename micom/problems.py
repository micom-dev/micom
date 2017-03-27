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
    """Check whether a community already carries a modification.

    Arguments
    ---------
    community : micom.Community
        The community class to check.

    Raises
    ------
    ValueError
        If the community already carries a modification and adding another
        would not be safe.
    """
    if community.modification is not None:
        raise ValueError("Community already carries a modification "
                         "({})!".format(community.modification))


def _format_min_growth(min_growth, species):
    """Format min_growth into a pandas series.

    Arguments
    ---------
    min_growth : positive float or array-like object.
        The minimum growth rate for each individual in the community. Either
        a single value applied to all individuals or one value for each.
    species : array-like
        The ID for each individual model in the community.

    Returns
    -------
    pandas.Series
        A pandas Series mapping each individual to its minimum growth rate.
    """
    try:
        min_growth = float(min_growth)
    except (TypeError, ValueError):
        if len(min_growth) != len(species):
            raise ValueError(
                "min_growth must be single value or an array-like "
                "object with an entry for each species in the model.")
    return pd.Series(min_growth, species)


def add_linear_optcom(community, min_growth=0.1):
    """Add a linear version of optcom.

    Adds constraints to the community such that each individual grows with
    at least its minimal growth rate.

    Arguments
    ---------
    community : micom.Community
        The community to modify.
    min_growth : positive float or array-like object.
        The minimum growth rate for each individual in the community. Either
        a single value applied to all individuals or one value for each.
    """
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
    """Adds a Lagrangian optimization target to a linear OptCom model.

    Lagrangian forms in `micom` are usually objectives of the form
    (1 - tradeoff) * community_objective + tradeoff * cooperativity_cost.
    Here, cooperativity cost specifies the "sacrifice" in growth rate an
    individual has to make in order to maximize community growth. It is
    calculated as the sum of squared differences between the individuals
    current and maximal growth rate. In the linear case squares are substituted
    by absolute values (Manhattan distance).

    Arguments
    ---------
    community : micom.Community
        The community to modify.
    tradeoff : float in [0,1]
        The tradeoff between community growth and the individual "egoistic"
        growth target. 0 means optimize only community growth and 1 means
        optimize only the individuals own growth target.
    linear : boolean
        Whether to use a non-linear (sum of squares) or linear version of the
        cooperativity cost. If set to False requires a QP-capable solver.
    """
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
    """Add dual Optcom variables and constraints to a community.

    Uses the original formulation of OptCom and solves the following
    multi-objective problem::

        maximize community_growth
        s.t. maximize growth_rate_i for all i
             s.t. Sv_i = 0
                  lb_i >= v_i >= ub_i

    Notes
    -----
    This method will only find one arbitrary solution from the Pareto front.
    There may exist several other optimal solutions.

    Arguments
    ---------
    community : micom.Community
        The community to modify.
    min_growth : positive float or array-like object.
        The minimum growth rate for each individual in the community. Either
        a single value applied to all individuals or one value for each.
    """
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
    """Add a dualized MOMA version of OptCom.

    Solves a MOMA (minimization of metabolic adjustment) formulation of OptCom
    given by::

        minimize cooperativity_cost
        s.t. maximize community_objective
             s.t. Sv = 0
                  lb >= v >= ub
        where community_cost = sum (growth_rate - max_growth)**2
              if linear=False or
              community_cost = sum |growth_rate - max_growth|
              if linear=True

    Arguments
    ---------
    community : micom.Community
        The community to modify.
    min_growth : positive float or array-like object.
        The minimum growth rate for each individual in the community. Either
        a single value applied to all individuals or one value for each.
    linear : boolean
        Whether to use a non-linear (sum of squares) or linear version of the
        cooperativity cost. If set to False requires a QP-capable solver.
    """
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
    """Run OptCom for the community.

    OptCom methods are a group of optimization procedures to find community
    solutions that provide a tradeoff between the cooperative community
    growth and the egoistic growth of each individual [1]. `micom`
    provides several strategies that can be used to find optimal solutions:

    - "linear": Applies a lower bound for the individual growth rates and
      finds the optimal community growth rate. This is the fastest methods
      but also ignores that individuals might strive to optimize their
      individual growth instead of community growth.
    - "lagrangian": Optimizes a joint objective containing the community
      objective (maximized) as well as a cooperativity cost which
      represents the  distance to the individuals "egoistic" maximum growth
      rate (minimized). Requires the `tradeoff` parameter. This method is
      still relatively fast and does require only few additional variables.
    - "linear lagrangian": The same as "lagrangian" only with a linear
      representation of the cooperativity cost (absolute value).
    - "moma": Minimization of metabolic adjustment. Simultaneously
      optimizes the community objective (maximize) and the cooperativity
      cost (minimize). This method finds an exact maximum but doubles the
      number of required variables, thus being slow.
    - "lmoma": The same as "moma" only with a linear
      representation of the cooperativity cost (absolute value).
    - "original": Solves the multi-objective problem described in [1].
      Here, the community growth rate is maximized simultanously with all
      individual growth rates. Note that there are usually many
      Pareto-optimal solutions to this problem and the method will only
      give one solution. This is also the slowest method.

    Parameters
    ----------
    community : micom.Community
        The community to optimize.
    strategy : str
        The strategy used to solve the OptCom formulation. Defaults to
        "lagrangian" which gives a decent tradeoff between speed and
        correctness.
    min_growth : float or array-like
        The minimal growth rate required for each individual. May be a
        single value or an array-like object with the same length as there
        are individuals.
    tradeoff : float in [0, 1]
        Only used for lagrangian strategies. Must be between 0 and 1 and
        describes the strength of the cooperativity cost / egoism. 1 means
        optimization will only minimize the cooperativity cost and zero
        means optimization will only maximize the community objective.
    fluxes : boolean
        Whether to return the fluxes as well.
    pfba : boolean
        Whether to obtain fluxes by parsimonious FBA rather than
        "classical" FBA.

    Returns
    -------
    tuple
        For fluxes=False a tuple of (community_gc, gcs) containing the
        overall community growth rates and a pandas series containing the
        individual growth rates. For fluxes=True a tuple
        (community_gc, fluxes) containing the overall community growth
        rates and a pandas data frame containing the fluxes.

    References
    ----------
    .. [1] OptCom: a multi-level optimization framework for the metabolic
       modeling and analysis of microbial communities.
       Zomorrodi AR, Maranas CD. PLoS Comput Biol. 2012 Feb;8(2):e1002363.
       doi: 10.1371/journal.pcbi.1002363, PMID: 22319433
    """
    if strategy not in _methods:
        raise ValueError("strategy must be one of {}!".format(
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
