"""Functions to summarize interactions over all metabolites."""

import pandas as pd


def _summarize(ints: pd.DataFrame) -> pd.DataFrame:
    """Summarize the overall interactions."""
    return ints.groupby("class").apply(
        lambda df: pd.DataFrame(
            {
                "flux": df.flux.sum(),
                "mass_flux": (df.flux * df.molecular_weight).sum() * 1e-3,
                "C_flux": (df.flux * df.C_number).sum(),
                "N_flux": (df.flux * df.N_number).sum(),
                "n_ints": df.metabolite.count(),
            },
            index=[0],
        )
    )


def summarize_interactions(ints: pd.DataFrame) -> pd.DataFrame:
    """Summarize interactions to key quantities.

    Arguments
    ---------
    ints : pandas.DataFrame
        The interactions for individual metabolites calculated before.

    Returns
    -------
    pandas.DataFrame
        The summarized interactions contaning the total flux, mass flux, carbon flux,
        nitrogen flux and number of interactions between any pair of taxa in that
        sample.
    """
    return (
        ints.groupby(["sample_id", "focal", "partner"])
        .apply(_summarize)
        .reset_index()
        .drop("level_4", axis=1)
    )
