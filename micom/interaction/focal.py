"""Quantify metabolic interactions between taxa."""

from ..taxonomy import taxon_id
from ..batch import GrowthResults, workflow
import numpy as np
import pandas as pd
import numpy as np
from typing import List, Union


def sample_interactions(fluxes: pd.DataFrame, taxon: str) -> pd.DataFrame:
    """Quantify interactions in a single sample.

    This attempts to quantify the interactions of a focal taxon with other taxa in a single sample. It does so by comparing the direction
    of fluxes of the focal taxon with those of other taxa.

    The interaction is classified as follows:

    - co-consumed: both the focal taxon and the partner taxon are consuming the same metabolite (both have import fluxes).
    - provided: the focal taxon is producing a metabolite that the partner taxon is consuming (focal has export flux, partner has import flux).
    - received: the focal taxon is consuming a metabolite that the partner taxon is producing (focal has import flux, partner has export flux).


    Parameters
    ----------
    fluxes : pd.DataFrame
        The fluxes of a single sample.
    taxon : str
        The focal taxon to quantify interactions for.

    Returns
    -------
    pd.DataFrame or None
        A dataframe with the interactions of the focal taxon with other taxa in the sample. If there are no interactions, returns None.
    """
    # Add scale column to indicate direction of flux relative to focal taxon
    df = fluxes.copy()
    df["scale"] = np.where(df["direction"] == "import", -1, 1)

    # Extract focal taxon data across all samples
    focal = df[df["taxon"] == taxon]
    if focal.empty:
        return None

    focal_side = focal[["sample_id", "metabolite", "scale"]].rename(
        columns={"scale": "focal_scale"}
    )
    focal_side["focal_flux"] = focal["flux"].abs() * focal["abundance"]

    # Extract partner data across all samples (excluding focal and medium)
    partners = df[(df["taxon"] != taxon) & (df["taxon"] != "medium")]
    if partners.empty:
        return None

    partner_side = partners[["sample_id", "metabolite", "taxon", "scale"]].rename(
        columns={"taxon": "partner"}
    )
    partner_side["partner_flux"] = partners["flux"].abs() * partners["abundance"]

    # Add the focal scale to the partner flux
    merged = pd.merge(focal_side, partner_side, on=["sample_id", "metabolite"])
    if merged.empty:
        return None

    conditions = [
        (merged["scale"] < 0) & (merged["focal_scale"] < 0),
        (merged["scale"] < 0) & (merged["focal_scale"] > 0),
        (merged["scale"] > 0) & (merged["focal_scale"] < 0),
    ]
    choices = ["co-consumed", "provided", "received"]

    merged["class"] = np.select(conditions, choices, default="none")
    merged = merged[merged["class"] != "none"]
    if merged.empty:
        return None

    merged["focal"] = taxon
    merged["flux"] = merged[["focal_flux", "partner_flux"]].min(axis=1)

    return merged[
        ["focal", "partner", "metabolite", "class", "flux", "sample_id"]
    ].reset_index(drop=True)


def _interact(args: List) -> pd.DataFrame:
    """Quantify interactions of a focal taxon with other taxa."""
    results, taxon = args
    ex = results.exchanges[results.exchanges.taxon != "medium"]

    ints = sample_interactions(ex, taxon).merge(
        results.annotations.drop_duplicates(subset="metabolite"), on="metabolite"
    )

    return ints


def interactions(
    results: GrowthResults,
    taxa: Union[None, str, List[str]],
    threads: int = 1,
    progress: bool = True,
) -> pd.DataFrame:
    """Quantify interactions of a focal/reference taxon with other taxa.

    This parallelizes across taxa. Samples are not parallelized, as sample interactions can be vectorized quite well.

    Note
    ----
    The function will attempt to resolve taxa names to taxon IDs. If a taxon name cannot be resolved, a ValueError will be raised.

    Parameters
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

    Raises
    ------
    ValueError
        If a taxon name cannot be resolved to a taxon ID.
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
