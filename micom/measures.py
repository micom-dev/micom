"""Measures calculated based on fluxes."""

import pandas as pd
import numpy as np


def production_rates(results):
    """Calculate the production rates for a set of metabolites.

    Note
    ----
    Production rates are not net fluxes but the flux of a particular metabolite
    into the external environment independent whether it is taken up by another
    taxon. The net rates can be obtained directly by inspecting the exchange fluxes
    in the "medium" compartment.
    The (transient) production rates reported here are the fluxes other surrounding
    cells/consumers would have access to. Thus, they are often more interesting for
    a phenotype than the net rates which are the excess productin rates, even in the
    absence of another consumer.

    Arguments
    ---------
    results : micom.GrowthResults
        A growth results as returned by grow.

    Returns
    -------
    pandas.DataFrame
        A new flux DataFrame containing the intermediate production rates for each
        metabolite.
    """
    fluxes = results.exchanges
    pos = fluxes[(fluxes.direction == "export") & (fluxes.taxon != "medium")]
    rates = pos.groupby(["sample_id", "metabolite"]).apply(
        lambda df: pd.Series({"flux": np.sum(df.abundance * df.flux.abs())})
    ).reset_index()
    rates = pd.merge(rates, results.annotations, on="metabolite")
    return rates


def consumption_rates(results):
    """Calculate the production rates for a set of metabolites.

    Note
    ----
    Consumption rates are not net fluxes but the total consumption of a particular
    metabolite by all taxa in the community independent whether it is secreted again.
    The net rates can be obtained directly by inspecting the exchange fluxes
    in the "medium" compartment.
    The (transient) consumption rates reported here is the total flux of a metabolite
    imported/consumed by taxa in the community.

    Arguments
    ---------
    results : micom.GrowthResults
        A growth results as returned by grow.

    Returns
    -------
    pandas.DataFrame
        A new flux DataFrame containing the intermediate production rates for each
        metabolite.
    """
    fluxes = results.exchanges
    neg = fluxes[(fluxes.direction == "import") & (fluxes.taxon != "medium")]
    rates = neg.groupby(["sample_id", "metabolite"]).apply(
        lambda df: pd.Series({"flux": np.sum(df.abundance * df.flux.abs())})
    ).reset_index()
    rates = pd.merge(rates, results.annotations, on="metabolite")
    return rates
