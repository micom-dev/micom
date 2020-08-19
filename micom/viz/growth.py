"""Visualization for growth rates."""

from datetime import datetime
from micom.viz import Visualization


def plot_growth(
    results,
    filename="growth_rates_%s.html" % datetime.now().strftime("%Y%m%d"),
    tolerance=1e-6
):
    """Plot the taxa growth rates.

    Parameters
    ----------
    results : micom.workflows.GrowthResults
        The results returned by the `grow` workflow.
    filename : str
        The HTML file where the visualization will be saved.
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
    data = {"growth_rates": rates}
    viz = Visualization(filename, data, "growth.html")
    viz.save(
        data=rates.to_json(orient="records"), width=800, height=400
    )

    return viz
