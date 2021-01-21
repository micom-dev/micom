"""Manages functions for growth media analysis and manipulation."""

from cobra.core.formula import Formula, elements_and_molecular_weights
from optlang.symbolics import Zero
import numpy as np
import pandas as pd
from micom import Community
from micom.util import (
    _format_min_growth,
    _apply_min_growth,
    check_modification,
)
from micom.logger import logger
from micom.solution import OptimizationError
import warnings


def add_linear_obj(community, exchanges, weights):
    """Add a linear version of a minimal medium to the community.

    Changes the optimization objective to finding the growth medium requiring
    the smallest total import flux::

        minimize sum w_i|r_i| for r_i in import_reactions

    Arguments
    ---------
    community : micom.Community
        The community to modify.
    exchanges : list of cobra.Reaction
        The reactions to constrain.
    weights : dict
        Maps each exchange reaction to a weight that is used in the
        minimization.
    """
    check_modification(community)
    coefs = {}
    for rxn in exchanges:
        export = len(rxn.reactants) == 1 or (
            len(rxn.metabolites) == 2 and rxn.products[0].compartment == "m"
        )
        met = list(rxn.metabolites)[0]
        if export:
            coefs[rxn.reverse_variable] = weights[met] * community.scale
        else:
            coefs[rxn.forward_variable] = weights[met] * community.scale
    community.objective.set_linear_coefficients(coefs)
    community.objective.direction = "min"
    community.modification = "minimal medium linear"


def add_mip_obj(community, exchanges):
    """Add a mixed-integer version of a minimal medium to the community.

    Changes the optimization objective to finding the medium with the least
    components::

        minimize size(R) where R part of import_reactions

    Arguments
    ---------
    community : micom.Community
        The community to modify.
    exchanges : list of cobra.Reaction
        The reactions to constrain.
    """
    check_modification(community)
    if len(community.variables) > 1e4:
        logger.warning(
            "the MIP version of minimal media is extremely slow for"
            " models that large :("
        )
    boundary_rxns = exchanges
    M = max(np.max(np.abs(r.bounds)) for r in boundary_rxns)
    prob = community.problem
    coefs = {}
    to_add = []
    for rxn in boundary_rxns:
        export = len(rxn.reactants) == 1 or (
            len(rxn.reactants) == 2 and rxn.products[0].compartment == "m"
        )
        indicator = prob.Variable("ind_" + rxn.id, lb=0, ub=1, type="binary")
        if export:
            vrv = rxn.reverse_variable
            indicator_const = prob.Constraint(
                vrv - indicator * M, ub=0, name="ind_constraint_" + rxn.id
            )
        else:
            vfw = rxn.forward_variable
            indicator_const = prob.Constraint(
                vfw - indicator * M, ub=0, name="ind_constraint_" + rxn.id
            )
        to_add.extend([indicator, indicator_const])
        coefs[indicator] = 1
    community.add_cons_vars(to_add)
    community.solver.update()
    community.objective.set_linear_coefficients(coefs)
    community.objective.direction = "min"
    community.modification = "minimal medium mixed-integer"


def safe_weight(met):
    """Get the weight of a molecule."""
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            w = max(Formula(met.formula).weight, 1.0)
    except Exception:
        w = 1.0
    return w


def weight(exchanges, what):
    """Obtain elemental weights for metabolites."""
    mets = [list(r.metabolites)[0] for r in exchanges]
    if what is None:
        weights = {m: 1.0 for m in mets}
    elif what == "mass":
        weights = {m: safe_weight(m) for m in mets}
    elif what in elements_and_molecular_weights:
        weights = {m: Formula(m.formula).elements.get(what, 1e-2) for m in mets}
    else:
        raise ValueError(
            "%s is not a valid elements. Must be one of: %s."
            % (what, ", ".join(elements_and_molecular_weights))
        )
    return weights


def minimal_medium(
    community,
    community_growth,
    min_growth=0.0,
    exports=False,
    exchanges=None,
    minimize_components=False,
    open_exchanges=False,
    solution=False,
    weights=None,
    atol=None,
    rtol=None,
):
    """Find the minimal growth medium for the community.

    Finds the minimal growth medium for the community which allows for
    community as well as individual growth. Here, a minimal medium can either
    be the medium requiring the smallest total import flux or the medium
    requiring the least components (ergo ingredients).

    Arguments
    ---------
    community : micom.Community
        The community to modify.
    community_growth : positive float
        The minimum community-wide growth rate.
    min_growth : positive float or array-like object.
        The minimum growth rate for each individual in the community. Either
        a single value applied to all individuals or one value for each.
    exports : boolean
        Whether to include export fluxes in the returned medium. Defaults to
        False which will only return import fluxes.
    exchanges : list of cobra.Reactions
        The list of exchange reactions that are penalized.
    minimize_components : boolean
        Whether to minimize the number of components instead of the total
        import flux. Might be more intuitive if set to True but may also be
        slow to calculate for large communities.
    open_exchanges : boolean or number
        Whether to ignore currently set bounds and make all exchange reactions
        in the model possible. If set to a number all exchange reactions will
        be opened with (-number, number) as bounds.
    solution : boolean
        Whether to also return the entire solution and all fluxes for the
        minimal medium.
    weights : str
        Will scale the fluxes by a weight factor. Can either be "mass" which will
        scale by molecular mass, a single element which will scale by
        the elemental content (for instance "C" to scale by carbon content).
        If None every metabolite will receive the same weight.
        Will be ignored if `minimize_components` is True.
    atol : float
        Absolute tolerance for the growth rates. If None will use the solver tolerance.
    rtol : float
        Relative tolerqance for the growth rates. If None will use the solver tolerance.


    Returns
    -------
    pandas.Series or dict
        A series {rid: flux} giving the import flux for each required import
        reaction. If `solution` is True retuns a dictionary
        {"medium": panas.Series, "solution": micom.CommunitySolution}.

    """
    logger.info("calculating minimal medium for %s" % community.id)

    if atol is None:
        atol = community.solver.configuration.tolerances.optimality
    if rtol is None:
        rtol = community.solver.configuration.tolerances.optimality

    if exchanges is None:
        boundary_rxns = community.exchanges
    else:
        boundary_rxns = community.reactions.get_by_any(exchanges)
    if isinstance(open_exchanges, bool):
        open_bound = 1000
    else:
        open_bound = open_exchanges
    min_growth = _format_min_growth(min_growth, community.taxa)
    with community as com:
        if open_exchanges:
            logger.info("opening exchanges for %d imports" % len(boundary_rxns))
            for rxn in boundary_rxns:
                rxn.bounds = (-open_bound, open_bound)
        logger.info("applying growth rate constraints")
        const = community.problem.Constraint(
            community.objective.expression,
            lb=community_growth,
            name="micom_growth_const",
        )
        community.add_cons_vars([const])
        _apply_min_growth(community, min_growth, atol, rtol)
        com.objective = Zero
        logger.info("adding new media objective")
        if minimize_components:
            add_mip_obj(com, boundary_rxns)
        else:
            scales = weight(boundary_rxns, weights)
            add_linear_obj(com, boundary_rxns, scales)
        sol = com.optimize(fluxes=True, pfba=False)
        if sol is None:
            logger.warning("minimization of medium was unsuccessful")
            return None

        logger.info("formatting medium")
        medium = pd.Series()
        set_medium = community.medium
        ex = [r for r in com.exchanges if r.id in set_medium]
        for rxn in ex:
            export = len(rxn.reactants) == 1
            flux = sol.fluxes.loc["medium", rxn.id]
            if abs(flux) < atol:
                continue
            if export:
                medium[rxn.id] = -flux
            elif not export:
                medium[rxn.id] = flux
        if not exports:
            medium = medium[medium > 0.0]

    if solution:
        return {"medium": medium, "solution": sol}
    else:
        return medium


def complete_medium(
    model,
    medium,
    growth=0.1,
    min_growth=0.001,
    max_import=1,
    minimize_components=False,
    weights=None,
):
    """Fill in missing components in a growth medium.

    Finds the minimal number of additions to make a model form biomass. In
    order to avoid bias all added reactions will have a maximum import
    rate of `max_import`.

    Note
    ----
    This function fixes the growth medium for a single cobra Model. We also
    provide a function `fix_medium` in `micom.workflows` that fixes a growth
    medium for an entire model database.

    Arguments
    ---------
    model : cobra.Model
        The model to use.
    medium : pandas.Series
        A growth medium. Must contain positive floats as elements and
        exchange reaction ids as index. Note that reactions not present in the
        model will be removed from the growth medium.
    growth : positive float
        The minimum overall growth rate that has to be achieved. For single COBRA
        model this is just the biomass flux and for community models this is the
        community biomass flux.
    min_growth : positive float or array-like object.
        The minimum growth rate for each individual in the community. Either
        a single value applied to all individuals or one value for each. Only used
        if model is a `micom.Community` model.
    minimize_components : boolean
        Whether to minimize the number of components instead of the total
        import flux. Might be more intuitive if set to True but may also be
        slow to calculate for large communities.
    max_import: positive float
        The import rate applied for the added exchanges.
    weights : str
        Will scale the fluxes by a weight factor. Can either be "mass" which will
        scale by molecular mass, a single element which will scale by
        the elemental content (for instance "C" to scale by carbon content).
        If None every metabolite will receive the same weight.
        Will be ignored if `minimize_components` is True.


    Returns
    -------
    pandas.Series or dict
        A series {rid: flux} giving the import flux for each required import
        reaction. This will include the initial `medium` as passed to the
        function as well as a minimal set of additional changes such that the
        model produces biomass with a rate >= `min_growth`.

    """
    exids = [r.id for r in model.exchanges]
    candidates = [r for r in model.exchanges if r.id not in medium.index]
    medium = medium[[i for i in medium.index if i in exids]]
    tol = model.solver.configuration.tolerances.feasibility
    with model:
        model.modification = None
        const = model.problem.Constraint(
            model.objective.expression,
            lb=growth,
            name="micom_growth_const",
        )
        if isinstance(model, Community):
            min_growth = _format_min_growth(min_growth, model.taxa)
            _apply_min_growth(model, min_growth, tol, tol)
        model.add_cons_vars([const])
        model.objective = Zero
        model.medium = medium.to_dict()
        for ex in candidates:
            export = len(ex.reactants) == 1
            if export:
                ex.lower_bound = -max_import
            else:
                ex.upper_bound = max_import
        if minimize_components:
            add_mip_obj(model, candidates)
        else:
            scales = weight(candidates, weights)
            add_linear_obj(model, candidates, scales)
        if isinstance(model, Community):
            sol = model.optimize(fluxes=True, pfba=False)
            fluxes = sol.fluxes.loc["medium", :]
        else:
            try:
                sol = model.optimize(raise_error=True)
                fluxes = sol.fluxes
            except OptimizationError:
                sol = None
    if sol is None:
        raise OptimizationError(
            "Could not find a solution that completes the medium :("
        )
    completed = pd.Series()
    for rxn in model.exchanges:
        export = len(rxn.reactants) == 1
        if rxn.id in medium.index:
            completed[rxn.id] = medium[rxn.id]
            continue
        else:
            flux = -fluxes[rxn.id] if export else fluxes[rxn.id]
        if abs(flux) < tol:
            continue
        completed[rxn.id] = flux

    return completed[completed > 0]
