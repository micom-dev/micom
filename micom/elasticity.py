"""Calculate elasticity coefficients.

Functions to calculate elasticity coefficients for various community
quantities.
"""

from functools import partial
import pandas as pd
import numpy as np
from tqdm.auto import tqdm
from cobra.util import get_context
from micom.util import reset_min_community_growth
from micom.problems import regularize_l2_norm
from micom.solution import optimize_with_fraction


STEP = 0.1


def _get_fluxes(sol, reactions):
    """Get the primal values for a set of variables."""
    fluxes = {
        r.id: sol.fluxes.loc[r.community_id, r.global_id] for r in reactions
    }
    return pd.Series(fluxes)


def _derivatives(before, after):
    """Get the elasticities for fluxes."""
    before_signs = np.sign(before)
    after_signs = np.sign(after)
    if any(np.abs(before_signs - after_signs) > 2):
        ValueError(
            "Some of the fluxes changed sign. " "Can't compute elasticities :("
        )
    direction = np.repeat("zero", len(before)).astype("<U8")
    direction[(before > 1e-6) | (after > 1e-6)] = "forward"
    direction[(before < -1e-6) | (after < -1e-6)] = "reverse"
    derivs = (np.log(after.abs() + 1e-6) - np.log(before.abs() + 1e-6)) / STEP
    return derivs, direction


def elasticities_by_medium(com, reactions, fraction, growth_rate, progress):
    """Get the elasticity coefficients for a set of variables.

    Arguments
    ---------
    com : micom.Community
        The community for wrhich to calculate elasticities.
    variables : list of optlang.Variable
        The variables for which to calculate the elasticities. All of these
        must have non-zero primal vaues in the previous solution.

    Returns
    -------
    pandas.Dataframe
        The long/tidy version of the elasticities. Contains columns variable,
        effector, and elasticity.
    """
    regularize_l2_norm(com, 0.0)
    sol = optimize_with_fraction(com, fraction, growth_rate, True)
    before = _get_fluxes(sol, reactions)
    import_fluxes = pd.Series()
    dfs = []

    for ex in com.exchanges:
        export = len(ex.reactants) == 1
        flux = sol.fluxes.loc[ex.community_id, ex.global_id]
        if export and (flux < -1e-6):
            import_fluxes[ex] = flux
        elif not export and (flux > 1e-6):
            import_fluxes[ex] = -flux
        else:
            continue

    fluxes = import_fluxes.index
    if progress:
        fluxes = tqdm(fluxes, unit="optimizations", desc="medium")
    for r in fluxes:
        flux = import_fluxes[r]
        with com:
            if flux < -1e-6:
                r.lower_bound *= np.exp(STEP)
            else:
                r.upper_bound *= np.exp(STEP)
            sol = optimize_with_fraction(com, fraction, growth_rate, True)
            after = _get_fluxes(sol, reactions)
        deriv, dirs = _derivatives(before, after)
        res = pd.DataFrame(
            {
                "reaction": [rx.global_id for rx in reactions],
                "taxon": [list(r.compartments)[0] for r in reactions],
                "effector": r.id,
                "direction": dirs,
                "elasticity": deriv,
            }
        )
        dfs.append(res)

    return pd.concat(dfs)


def elasticities_by_abundance(com, reactions, fraction, growth_rate, progress):
    """Get the elasticity coefficients for a set of variables.

    Arguments
    ---------
    com : micom.Community
        The community for which to calculate elasticities.
    variables : list of optlang.Variable
        The variables for which to calculate the elasticities. All of these
        must have non-zero primal vaues in the previous solution.

    Returns
    -------
    pandas.Dataframe
        The long/tidy version of the elasticities. Contains columns variable,
        effector, and elasticity.
    """
    regularize_l2_norm(com, 0.0)
    sol = optimize_with_fraction(com, fraction, growth_rate, True)
    before = _get_fluxes(sol, reactions)
    dfs = []

    abundance = com.abundances.copy()
    taxa = abundance.index

    if progress:
        taxa = tqdm(taxa, unit="optimizations", desc="taxa abundances")
    for sp in taxa:
        old = abundance[sp]
        abundance.loc[sp] *= np.exp(STEP)
        com.set_abundance(abundance, normalize=False)
        sol = optimize_with_fraction(com, fraction, growth_rate, True)
        after = _get_fluxes(sol, reactions)
        abundance.loc[sp] = old
        com.set_abundance(abundance, normalize=False)
        deriv, dirs = _derivatives(before, after)
        res = pd.DataFrame(
            {
                "reaction": [r.global_id for r in reactions],
                "taxon": [list(r.compartments)[0] for r in reactions],
                "effector": sp,
                "direction": dirs,
                "elasticity": deriv,
            }
        )
        dfs.append(res)

    return pd.concat(dfs)


def elasticities(com, fraction=0.5, reactions=None, progress=True):
    """Calculate elasticities for reactions.

    Calculates elasticity coefficients using the specified reactions as
    response and exchange bounds (diet) and taxa abundances as
    effectors/parameters. Will use an arbitrary flux distribution as base.

    Arguments
    ---------
    com : micom.Community
        The community for wrhich to calculate elasticities.
    fraction : double
        The tradeoff to use for the cooperative tradeoff method. Fraction of
        maximal community growth to enforce.
    reactions : iterable
        A list of reactions to get elasticities for. Elements can either be
        reactions from the model, strings specifying the ids of reactions
        or ints specifying the indices of reactions. Defaults to using all
        reactions.
    progress : boolean
        Whether to shwo progress bars. Will show two, one for the diet
        optimizations and another one for the taxa abundances.

    Returns
    -------
    pandas.DataFrame
        A data frame with the following columns:
        "reaction" - the exchange reaction (response),
        "taxon" - the taxon the reaction is from,
        "effector" - the parameter that was changed,
        "direction" - whether the flux runs in the forward or reverse
            direction,
        "elasticity" - the elasticity coefficient,
        "type" - the type of effector either "exchange" for diet or "abundance"
        for taxa abundances.
    """
    growth_rate = None
    if reactions is None:
        reactions = com.reactions
    reactions = com.reactions.get_by_any(reactions)
    with com:
        context = get_context(com)
        context(partial(reset_min_community_growth, com))
        by_medium = elasticities_by_medium(
            com, reactions, fraction, growth_rate, progress
        )
        by_medium["type"] = "exchanges"

        by_abundance = elasticities_by_abundance(
            com, reactions, fraction, growth_rate, progress
        )
        by_abundance["type"] = "abundance"

    both = pd.concat([by_medium, by_abundance]).reset_index(drop=True)
    both.loc[both.taxon == "m", "taxon"] = "medium"
    return both
