"""Calculate elasticity coefficients.

Functions to calculate elasticity coefficients for various community
quantities.
"""

from functools import partial
import pandas as pd
import numpy as np
from tqdm import tqdm
from cobra.util import get_context
from micom.util import reset_min_community_growth
from micom.problems import regularize_l2_norm
from micom.solution import optimize_with_fraction
from micom.media import minimal_medium


STEP = 0.1


def _get_fluxes(sol, reactions):
    """Get the primal values for a set of variables."""
    fluxes = {r.id: sol.fluxes.loc[r.community_id, r.global_id]
              for r in reactions}
    return pd.Series(fluxes)


def _derivatives(before, after):
    """Get the elasticities for fluxes."""
    before_signs = np.sign(before)
    after_signs = np.sign(after)
    if any(np.abs(before_signs - after_signs) > 2):
        ValueError("Some of the fluxes changed sign. "
                   "Can't compute elasticities :(")
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
            fluxes = tqdm(fluxes, unit="optimizations")
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
        res = pd.DataFrame({"reaction": [rx.id for rx in reactions],
                            "effector": r.id,
                            "direction": dirs,
                            "elasticity": deriv})
        dfs.append(res)

    return pd.concat(dfs)


def elasticities_by_abundance(com, reactions, fraction, growth_rate, progress):
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
    dfs = []

    abundance = com.abundances.copy()
    species = abundance.index

    if progress:
            species = tqdm(species, unit="optimizations")
    for sp in species:
        old = abundance[sp]
        abundance.loc[sp] *= np.exp(STEP)
        com.set_abundance(abundance, normalize=False)
        sol = optimize_with_fraction(com, fraction, growth_rate, True)
        after = _get_fluxes(sol, reactions)
        abundance.loc[sp] = old
        com.set_abundance(abundance, normalize=False)
        deriv, dirs = _derivatives(before, after)
        res = pd.DataFrame({"reaction": [r.id for r in reactions],
                            "effector": sp,
                            "direction": dirs,
                            "elasticity": deriv})
        dfs.append(res)

    return pd.concat(dfs)


def exchange_elasticities(com, fraction=0.5, min_medium=True,
                          progress=True):
    """Calculate elasticities for exchange reactions."""
    growth_rate = None
    if min_medium:
        sol = com.cooperative_tradeoff(fraction)
        gcs = sol.members.growth_rate.drop("medium")
        med = minimal_medium(com, 0.95 * sol.growth_rate,
                             min_growth=0.95 * gcs)
        fraction = 1.0
        growth_rate = sol.growth_rate
    with com:
        context = get_context(com)
        context(partial(reset_min_community_growth, com))
        if min_medium:
            com.medium = med
        rxns = com.exchanges
        by_medium = elasticities_by_medium(com, rxns, fraction,
                                           growth_rate, progress)
        by_medium["type"] = "exchanges"

        by_abundance = elasticities_by_abundance(com, rxns, fraction,
                                                 growth_rate, progress)
        by_abundance["type"] = "abundance"

    return pd.concat([by_medium, by_abundance])
