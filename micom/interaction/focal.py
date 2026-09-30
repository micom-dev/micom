"""Quantify metabolic interactions between taxa."""

from ..taxonomy import taxon_id
from ..batch import GrowthResults, workflow
import numpy as np
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


import numpy as np


def sample_interactions(fluxes: pd.DataFrame, taxon: str) -> pd.DataFrame:
    """Quantify interactions in a single sample (high-performance)."""
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

    # Extract partner data across all samples (excluding focal and medium)
    partners = df[(df["taxon"] != taxon) & (df["taxon"] != "medium")]
    if partners.empty:
        return None

    partner_side = partners[
        ["sample_id", "metabolite", "taxon", "flux", "scale"]
    ].rename(columns={"taxon": "partner"})
    partner_side.loc

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
    merged["flux"] = merged["flux"].abs()

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
