"""A class for tracking results for growth simulations."""

from __future__ import annotations
from functools import reduce
from dataclasses import dataclass
import pandas as pd
import typing
from zipfile import ZipFile
from ..annotation import annotate_metabolites_from_exchanges

if typing.TYPE_CHECKING:
    from ..community import Community
    from ..solution import CommunitySolution

DIRECTION = pd.Series(["import", "export"], index=[0, 1])


@dataclass
class GrowthResults:
    growth_rates: pd.DataFrame
    exchanges: pd.DataFrame
    annotations: pd.DataFrame

    def save(self, path: str):
        """Save growth results to a file.

        This will write all tables as CSV into a single ZIP file.

        Arguments
        ---------
        path : str
            A filepath for the generated file. Should end in `.zip`.
        """
        with ZipFile(path, "w") as zippy:
            for attr in ["growth_rates", "exchanges", "annotations"]:
                getattr(self, attr).to_csv(zippy.open(f"{attr}.csv", "w"), index=False)

    @staticmethod
    def load(path: str) -> GrowthResults:
        """Load growth results from a file.

        Arguments
        ---------
        path : str
            Path to saved `GrowthResults`.

        Returns
        -------
        GrowthResults
            The loaded growth results.
        """
        tables = []
        with ZipFile(path, "r") as zippy:
            for attr in ["growth_rates", "exchanges", "annotations"]:
                tab = pd.read_csv(zippy.open(f"{attr}.csv", "r"))
                tables.append(tab)
        return GrowthResults(*tables)

    def __add__(self, other: GrowthResults) -> GrowthResults:
        """Combine two GrowthResults objects.

        Arguments
        ---------
        other : GrowthResults
            The other result to merge with the current one.

        Returns
        -------
        GrowthResult
            A merged growth result containing data from both.
        """
        rates = pd.concat([self.growth_rates, other.growth_rates])
        exs = pd.concat([self.exchanges, other.exchanges])
        anns = pd.concat([self.annotations, other.annotations]).drop_duplicates(
            subset=["reaction"]
        )
        anns.index = anns.reaction
        return GrowthResults(rates, exs, anns)

    @staticmethod
    def from_solution(sol: CommunitySolution, com: Community) -> GrowthResults:
        """Convert a solution to growth results."""
        if sol.fluxes is None:
            raise ValueError("Solution needs to contain fluxes to be converted.")
        tol = com.solver.configuration.tolerances.feasibility

        # Get the main objects
        rates = sol.members.drop("medium")
        rates["taxon"] = rates.index
        rates["sample_id"] = com.id
        fluxes = sol.exchange_fluxes
        fluxes["sample_id"] = com.id
        fluxes["tolerance"] = tol
        anns = annotate_metabolites_from_exchanges(com)
        anns.drop_duplicates(subset=["reaction"], inplace=True)

        anns.index = anns.reaction
        exchanges = pd.merge(fluxes, anns[["metabolite"]], on="reaction", how="left")
        exchanges = exchanges[exchanges.flux.abs() > exchanges.tolerance]
        return GrowthResults(rates, exchanges, anns)


def save_results(results: GrowthResults, path: str):
    """Save growth results to a file.

    This will write all tables as CSV into a single ZIP file.

    Arguments
    ---------
    results : GrowthResults
        The results as returned from `grow`.
    path : str
        A filepath for the generated file. Should end in `.zip`.
    """
    results.save(path)


def load_results(path):
    """Load growth results from a file.

    Arguments
    ---------
    path : str
        Path to saved `GrowthResults`.

    Returns
    -------
    GrowthResults
        The saved GrowthResults.
    """
    return GrowthResults.load(path)


def combine_results(it: typing.Iterable[GrowthResults]) -> GrowthResults:
    """Combine several GrowthResults.

    Arguments
    ---------
    it : Iterable of GrowthResults
        The growth results to combine.

    Returns
    -------
    GrowthResults
        The merged results.
    """
    return reduce(lambda x, y: x + y, it)
