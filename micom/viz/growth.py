"""Visualization for growth rates."""

from datetime import datetime
from micom.viz import Visualization
from micom.workflows.results import GrowthResults
import pandas as pd


def plot_growth(
    results: GrowthResults,
    filename: str = "growth_rates_%s.html" % datetime.now().strftime("%Y%m%d"),
    groups: pd.Series = None,
    tolerance: float = 1e-6,
) -> None:
    """Plot the taxa growth rates.

    Parameters
    ----------
    results : micom.workflows.GrowthResults
        The results returned by the `grow` workflow.
    filename : str
        The HTML file where the visualization will be saved.
    groups : pandas.Series
        Additional metadata to color points. The index must correspond to the
        `sample_id` in the results and values must be categorical. The `.name` attribute
        will be used to name the groups. If not provided will color by taxon.
    tolerance : float
        Smallest growth rate that will be considered. Everything lower will be
        considered not growing.

    Returns
    -------
    Visualization
        A MICOM visualization. Can be served with `viz.serve`.
    """
    rates = results.growth_rates
    rates = rates[rates.growth_rate > tolerance][
        ["taxon", "sample_id", "abundance", "growth_rate"]
    ]

    if groups is not None:
        if groups.dtype not in ["object", "category", "bool"]:
            raise ValueError("Groups need to be categorical.")
        name = "group" if groups.name is None else groups.name
        rates[name] = groups[rates.sample_id].values
    else:
        name = "taxon"

    data = {"growth_rates": rates}
    n_taxa = rates.taxon.nunique()
    viz = Visualization(filename, data, "growth.html")
    viz.save(data=rates.to_json(orient="records"), color=name, n_taxa=n_taxa)

    return viz
