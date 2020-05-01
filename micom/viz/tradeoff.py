"""Visualizations for tradeoff analysis."""

from datetime import datetime
from micom.viz import Visualization
import numpy as np
import pandas as pd


def plot_tradeoff(
    tradeoff_rates,
    out_folder="tradeoff_%s" % datetime.now().strftime("%Y%m%d"),
):
    """Plot diagnostics for varying tradeoff values.

    Parameters
    ----------
    tradeoff_rates : pandas.DataFrame
        The growth rates returned by the `tradeoff` workflow.
    out_folder : str
        The folder where the visualization will be saved.

    Returns
    -------
    Visualization
        A MICOM visualization. Can be served with `viz.serve`.

    """
    rates = tradeoff_rates
    data = {"tradeoff": rates}
    viz = Visualization(out_folder, data, "tradeoff.html")
    growth = rates[
        ["taxon", "sample_id", "abundance", "tradeoff", "growth_rate"]
    ].copy()
    growth.tradeoff = growth.tradeoff.round(6).astype(str)
    growth.loc[growth.tradeoff == "nan", "tradeoff"] = "none"
    growth.loc[growth.growth_rate < 1e-6, "growth_rate"] = 1e-6
    growth.loc[:, "log_growth_rate"] = np.log10(growth.growth_rate)
    tradeoff = (
        growth.groupby(["tradeoff", "sample_id"])
        .apply(
            lambda df: pd.Series(
                {
                    "n_taxa": df.shape[0],
                    "n_growing": df[df.growth_rate > 1e-6].shape[0],
                    "fraction_growing": (
                        df[df.growth_rate > 1e-6].shape[0] / df.shape[0]
                    ),
                }
            )
        )
        .reset_index()
    )
    viz.save(
        growth=growth.to_json(orient="records"),
        tradeoff=tradeoff.to_json(orient="records"),
        extent=[growth.log_growth_rate.min(), growth.log_growth_rate.max()],
        width=400,
        height=300,
    )

    return viz
