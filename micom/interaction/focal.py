"""Quantify metabolic interactions between taxa."""

from ..taxonomy import taxon_id
from ..workflows import GrowthResults, workflow
import pandas as pd
from typing import List, Union


def _metabolite_interaction(
    fluxes: pd.DataFrame, taxon: str, partner: str
) -> pd.DataFrame:
    """Checks if and how taxa interact."""
    tol = fluxes.tolerance.max()
    f = fluxes[(fluxes.flux.abs() * fluxes.abundance) > tol]
    if (f.shape[0] < 2) or (f.direction == "export").all():
        return None
    if (f.direction == "import").sum() == 2:
        int_type = "co-consumed"
    elif (f.loc[f.taxon == taxon, "direction"] == "export").all():
        int_type = "provided"
    else:
        int_type = "received"

    return pd.DataFrame(
        {
            "focal": taxon,
            "partner": partner,
            "class": int_type,
            "flux": (f.flux.abs() * f.abundance).min(),
        },
        index=[0],
    )


def sample_interactions(
    fluxes: pd.DataFrame, sample_id: str, taxon: str
) -> pd.DataFrame:
    """Quantify interactions in a single sammple.

    Arguments
    ---------
    fluxes : pandas.DataFrame
        A table of exchange fluxes.
    sample_id : str
        The sample id to use.
    taxon : str
        The focal taxon to use.

    Returns
    -------
    pandas.DataFrame
        The mapped interactions between the focal taxon and all other taxa.
    """
    ex = fluxes[fluxes.sample_id == sample_id]
    partners = pd.Series(ex.taxon.unique())
    partners = partners[(partners != taxon) & (partners != "medium")]
    ints = []
    for p in partners:
        fluxes = ex[ex.taxon.isin((taxon, p))]
        ints.append(
            fluxes.groupby("metabolite")
            .apply(lambda df: _metabolite_interaction(df, taxon, p))
            .reset_index()
        )
    ints = pd.concat([i for i in ints if i is not None])
    ints["sample_id"] = sample_id
    return ints


def _interact(args: List) -> pd.DataFrame:
    """Quantify interactions of a focal taxon with other taxa."""
    results, taxon = args
    ex = results.exchanges[results.exchanges.taxon != "medium"]

    ints = (
        ex.groupby("sample_id")
        .apply(lambda df: sample_interactions(df, df.name, taxon))
        .reset_index(drop=True)
        .drop(["level_1", "index"], axis=1, errors="ignore")
        .merge(results.annotations, on="metabolite")
    )

    return ints


def interactions(
    results: GrowthResults,
    taxa: Union[None, str, List[str]],
    threads: int = 1,
    progress: bool = True,
) -> pd.DataFrame:
    """Quantify interactions of a focal/reference taxon with other taxa.

    Arguments
    ---------
    results : GrowthResults
        The growth results to use.
    taxa : str, list of str, or None
        The focal taxa to use. Can be a single taxon, a list of taxa or None in which
        case all taxa are considered.

    Returns
    -------
    pandas.DataFrame
        The mapped interactions between the focal taxon and all other taxa.
    """
    if isinstance(taxa, str):
        return _interact([results, taxon_id(taxa, results.growth_rates)])
    elif taxa is None:
        taxa = results.growth_rates.taxon.unique()

    taxa = [taxon_id(t, results.growth_rates) for t in taxa]
    ints = pd.concat(
        workflow(
            _interact, [[results, t] for t in taxa], threads=threads, progress=progress
        )
    )
    return ints
