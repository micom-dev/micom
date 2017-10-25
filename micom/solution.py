"""A community solution object."""

import numpy as np
import pandas as pd
from collections import Counter
from optlang.interface import OPTIMAL
from cobra.core import Solution, get_solution


def _group_species(values, ids, species, what="reaction"):
    """Format a list of values by id and species."""
    df = pd.DataFrame({values.name: values, what: ids, "compartment": species})
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
        compartments: one for every species plus one for the external medium.
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
        rids = np.array([(r.global_id, r.community_id) for r in reactions])
        mids = np.array([(m.global_id, m.community_id)
                         for m in metabolites])
        if not slim:
            sol = get_solution(community, reactions, metabolites)
            super(CommunitySolution, self).__init__(
                community.solver.objective.value, community.solver.status,
                _group_species(sol.fluxes, rids[:, 0], rids[:, 1]),
                _group_species(sol.reduced_costs, rids[:, 0], rids[:, 1]),
                _group_species(sol.shadow_prices, mids[:, 0], mids[:, 1],
                               what="metabolites"))
        else:
            super(CommunitySolution, self).__init__(
                community.solver.objective.value, community.solver.status,
                None, None, None)
        gcs = pd.Series()
        for sp in community.species:
            gcs[sp] = community.constraints["objective_" + sp].primal
        self.strategy = community.modification
        self.members = pd.DataFrame({
            "abundance": community.abundances,
            "growth_rate": gcs,
            "reactions": pd.Series(Counter(rids[:, 1])),
            "metabolites": pd.Series(Counter(mids[:, 1]))})
        self.members.index.name = "compartments"
        self.growth_rate = sum(community.abundances * gcs)

    def __repr__(self):
        """Convert CommunitySolution instance to string representation."""
        if self.status != OPTIMAL:
            return "<CommunitySolution {0:s} at 0x{1:x}>".format(
                self.status, id(self))
        return "<CommunitySolution {0:.3f} at 0x{1:x}>".format(
            self.growth_rate, id(self))
