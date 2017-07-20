"""A community solution object."""

import numpy as np
import pandas as pd
from optlang.interface import OPTIMAL
from cobra.core import Solution, get_solution


def _group_species(values, ids, species, what="reaction"):
    """Format a list of values by id and species."""
    df = pd.DataFrame({values.name: values, what: ids, "species": species})
    df = df.pivot(index="species", columns=what, values=values.name)
    df.name = values.name
    return df


class CommunitySolution(Solution):
    """An FBA solution for an entire community.

    Attributes
    ----------
    objective_value : float
        The (optimal) value for the objective function.
    members : pandas.Series
        Contains basic info about the individual members of the community such
        as id, abundance and growth rates.
    growth_rate : float
        The overall growth rate for the community normalized to 1 gDW.
    status : str
        The solver status related to the solution.
    fluxes : pandas.DataFrame
        Contains the reaction fluxes (primal values of variables) stratified
        by species. Columns denote individual fluxes and rows denote species.
        Fluxes will be NA if the reaction does not exist in the organism.
    reduced_costs : pandas.Series
        Contains reaction reduced costs (dual values of variables) stratified
        by species. Columns denote individual fluxes and rows denote species.
        Reduced costs will be NA if the reaction does not exist in the
        organism.
    shadow_prices : pandas.Series
        Contains metabolite shadow prices (dual values of constraints)
        stratified by species. Columns denote individual metabolites and rows
        denote species. Shadow prices will be NA if the metabolite does not
        exist in the organism.

    """

    def __init__(self, community, slim=False,
                 reactions=None, metabolites=None):
        """Get the solution from a community model."""
        if reactions is None:
            reactions = community.reactions
        if metabolites is None:
            metabolites = community.metabolites
        if not slim:
            rids = np.array([(r.global_id, r.community_id) for r in reactions])
            mids = np.array([(m.global_id, m.community_id)
                             for m in metabolites])
            sol = get_solution(community, reactions, metabolites)
            super(CommunitySolution, self).__init__(
                community.solver.objective.value, community.solver.status,
                np.unique(rids[:, 0]),
                _group_species(sol.fluxes, rids[:, 0], rids[:, 1]),
                _group_species(sol.reduced_costs, rids[:, 0], rids[:, 1]),
                np.unique(mids[:, 0]),
                _group_species(sol.shadow_prices, mids[:, 0], mids[:, 1],
                               what="metabolites"))
        else:
            super(CommunitySolution, self).__init__(
                community.solver.objective.value, community.solver.status,
                None, None, None, None, None)
        gcs = pd.Series()
        for sp in community.objectives:
            gcs[sp] = community.constraints["objective_" + sp].primal
        self.members = pd.DataFrame({"id": gcs.index,
                                     "abundance": community.abundances,
                                     "growth_rate": gcs})
        self.growth_rate = sum(community.abundances * gcs)
        del self.reactions
        del self.metabolites

    def __repr__(self):
        """Convert CommunitySolution instance to string representation."""
        if self.status != OPTIMAL:
            return "<CommunitySolution {0:s} at 0x{1:x}>".format(
                self.status, id(self))
        return "<CommunitySolution {0:.3f} at 0x{1:x}>".format(
            self.community_growth, id(self))
