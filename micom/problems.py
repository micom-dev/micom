"""Implements optimization and model problems."""

from micom.util import fluxes_from_primals
import pandas as pd
from cobra.flux_analysis.parsimonious import add_pfba


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


def linear_optcom(community, fractions=0.0, fluxes=False, pfba=True):
    """Get a community solution from a linear version of optcom."""
    species = list(community.objectives.keys())
    try:
        fractions = float(fractions)
    except (TypeError, ValueError):
        if len(fractions) != len(species):
            raise ValueError(
                "fraction must be single value or an array-like "
                "object with an entry for each species in the model.")

    fractions = pd.Series(fractions, species)

    max_gcs = community.optimize_all()
    prob = community.solver.interface
    to_add = []
    for sp in species:
        obj = prob.Constraint(community.objectives[sp],
                              name="objective_" + sp,
                              lb=max_gcs[sp] * fractions[sp])
        to_add.append(obj)

    community_obj = None
    gcs = pd.Series(None, species)
    com_fluxes = None
    with community as com:
        com.add_cons_vars(to_add)
        com.solver.optimize()
        if com.solver.status == "optimal":
            for sp in species:
                gcs[sp] = com.solver.constraints["objective_" + sp].primal
            community_obj = com.solver.objective.value

            if fluxes:
                com_fluxes = all_fluxes(com, pfba=pfba)

    if fluxes:
        return community_obj, com_fluxes
    else:
        return community_obj, gcs
