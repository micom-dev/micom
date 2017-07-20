"""Manages functions for growth media analysis and manipulation."""

from numbers import Number
from sympy.core.singleton import S
from optlang.interface import OPTIMAL
import numpy as np
import pandas as pd
from micom.util import (_format_min_growth, _apply_min_growth,
                        check_modification)
from micom.logger import logger


default_excludes = ["biosynthesis", "transcription", "replication", "sink",
                    "demand", "DM_", "SN_", "SK_"]
"""A list of sub-strings in reaction IDs that usually indicate that
the reaction is *not* an exchange reaction."""


def add_linear_obj(community):
    """Add a linear version of a minimal medium to the community.

    Changes the optimization objective to finding the growth medium requiring
    the smallest total import flux::

        minimize sum |r_i| for r_i in import_reactions

    Arguments
    ---------
    community : micom.Community
        The community to modify.
    """
    check_modification(community)
    coefs = {}
    for rxn in community.exchanges:
        export = len(rxn.reactants) == 1
        if export:
            coefs[rxn.reverse_variable] = 1
        else:
            coefs[rxn.forward_variable] = 1
    community.objective.set_linear_coefficients(coefs)
    community.objective.direction = "min"
    community.modification = "minimal medium linear"


def add_mip_obj(community):
    """Add a mixed-integer version of a minimal medium to the community.

    Changes the optimization objective to finding the medium with the least
    components::

        minimize size(R) where R part of import_reactions

    Arguments
    ---------
    community : micom.Community
        The community to modify.
    """
    check_modification(community)
    if len(community.variables) > 1e4:
        logger.warning("the MIP version of minimal media is extremely slow for"
                       " models that large :(")
    boundary_rxns = community.exchanges
    M = max(np.max(np.abs(r.bounds)) for r in boundary_rxns)
    prob = community.problem
    coefs = {}
    to_add = []
    for rxn in boundary_rxns:
        export = len(rxn.reactants) == 1
        indicator = prob.Variable("ind_" + rxn.id, lb=0, ub=1, type="binary")
        if export:
            vrv = rxn.reverse_variable
            indicator_const = prob.Constraint(
                vrv - indicator * M, ub=0, name="ind_constraint_" + rxn.name)
        else:
            vfw = rxn.forward_variable
            indicator_const = prob.Constraint(
                vfw - indicator * M, ub=0, name="ind_constraint_" + rxn.name)
        to_add.extend([indicator, indicator_const])
        coefs[indicator] = 1
    community.add_cons_vars(to_add)
    community.solver.update()
    community.objective.set_linear_coefficients(coefs)
    community.objective.direction = "min"
    community.modification = "minimal medium mixed-integer"


def minimal_medium(community, community_growth, min_growth=0.1,
                   minimize_components=False, open_exchanges=False):
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
    minimize_components : boolean
        Whether to minimize the number of components instead of the total
        import flux. Might be more intuitive if set to True but may also be
        slow to calculate for large communities.
    open_exchanges : boolean or number
        Whether to ignore currently set bounds and make all exchange reactions
        in the model possible. If set to a number all exchange reactions will
        be opened with (-number, number) as bounds.

    Returns
    -------
    pandas.Series
        A series {rid: flux} giving the import flux for each required import
        reaction.

    """
    logger.info("calculating minimal medium for %s" % community.id)
    boundary_rxns = community.exchanges
    if isinstance(open_exchanges, Number):
        open_bound = open_exchanges
    else:
        open_bound = 1000
    min_growth = _format_min_growth(
        min_growth, list(community.objectives.keys()))
    with community as com:
        if open_exchanges:
            logger.info("opening exchanges for %d imports" %
                        len(boundary_rxns))
            for rxn in boundary_rxns:
                rxn.bounds = (-open_bound, open_bound)
        logger.info("applying growth rate constraints")
        obj_const = com.problem.Constraint(
            com.objective.expression, lb=community_growth,
            name="medium_obj_constraint")
        com.add_cons_vars([obj_const])
        com.solver.update()
        _apply_min_growth(community, min_growth)
        com.objective = S.Zero
        logger.info("adding new media objective")
        if minimize_components:
            add_mip_obj(com)
        else:
            add_linear_obj(com)
        com.solver.optimize()
        if com.solver.status != OPTIMAL:
            logger.warning("minimization of medium was infeasible")
            return None

        logger.info("formatting medium")
        medium = pd.Series()
        for rxn in boundary_rxns:
            export = len(rxn.reactants) == 1
            flux = rxn.flux
            if export and flux < 0:
                medium[rxn.id] = -flux
            elif not export and flux > 0:
                medium[rxn.id] = flux

    return medium
