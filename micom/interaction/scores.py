"""Various interaction scores."""

from collections import Counter
import pandas as pd
from ..batch import GrowthResults


def _mes(df: pd.DataFrame) -> float:
    """Helper to calculate the MES score."""
    cn = Counter(df.direction)
    p, c = cn["export"], cn["import"]
    return pd.Series(2.0 * p * c / (p + c), index=["MES"])


def basic(results: GrowthResults, cutoff: float = None) -> pd.DataFrame:
    """Calculate basic exchange statistics for each metabolite and sample.

    Parameters
    ----------
    results : GrowthResults
        The growth results to use.
    cutoff : float
        The smallest flux to consider. Defaults to the solver tolerance.

    Returns
    -------
    pandas.DataFrame
        Producer and consumer counts and abundance-weighted consumption and
        production fluxes for each metabolite and sample, including annotations.
    """
    if cutoff is None:
        cutoff = results.exchanges.tolerance.max()
    fluxes = results.exchanges[
        (results.exchanges.flux.abs() > cutoff) & (results.exchanges.taxon != "medium")
    ].copy()
    groups = ["metabolite", "sample_id"]
    fluxes["weighted_flux"] = fluxes.abundance * fluxes.flux.abs()

    counts = fluxes.pivot_table(
        index=groups,
        columns="direction",
        values="taxon",
        aggfunc="nunique",
        fill_value=0,
    ).reindex(columns=["export", "import"], fill_value=0.0)
    counts.columns = ["producers", "consumers"]

    rates = fluxes.pivot_table(
        index=groups,
        columns="direction",
        values="weighted_flux",
        aggfunc="sum",
        fill_value=0,
    ).reindex(columns=["import", "export"], fill_value=0.0)
    rates.columns = ["consumption_flux", "production_flux"]

    scores = counts.join(rates).reset_index()
    annotations = results.annotations.drop_duplicates(subset=["metabolite"])
    return scores.merge(annotations, on="metabolite", how="inner")


def MES(results: GrowthResults, cutoff: float = None) -> pd.DataFrame:
    """Calculate the Metabolic Exchange Score (MES) for each metabolite.

    MES is the harmonic mean of producers and consumers for a chosen metabolite
    in one sample. High values indicate a large prevalence of cross-feeding for the
    particular metabolite. A value of zero indicates an absence of cross-feeding for
    the particular metabolite.

    Parameters
    ---------
    results : GrowthResults
        The growth results to use.
    cutoff : float
        The smallest flux to consider in the analysis. Will default to the
        solver tolerance if set to None.

    Returns
    -------
    pandas.DataFrame
        The scores for each metabolite and each sample including metabolite annotations.

    References
    ----------
    .. [1] Marcelino, V.R., et al.
           Disease-specific loss of microbial cross-feeding interactions in the human gut
           Nat Commun 14, 6546 (2023). https://doi.org/10.1038/s41467-023-42112-w
    """
    if cutoff is None:
        cutoff = results.exchanges.tolerance[0]
    fluxes = results.exchanges[
        (results.exchanges.flux.abs() > cutoff) & (results.exchanges.taxon != "medium")
    ]
    mes = fluxes.groupby(["metabolite", "sample_id"]).apply(_mes).reset_index()
    mes = mes.merge(
        results.annotations.drop_duplicates(subset=["metabolite"]),
        on="metabolite",
        how="inner",
    )
    return mes
