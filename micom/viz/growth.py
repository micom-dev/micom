"""Visualization for growth rates."""

from datetime import datetime
from micom.viz import Visualization


def plot_growth(
    growth_rates,
    out_folder="growth_rates_%s" % datetime.now().strftime("%Y%m%d")
):
    """Plot the taxa growth rates.

    Parameters
    ----------
    growth_rates : pandas.DataFrame
        The growth rates returned by the `grow` workflow.
    out_folder : str
        The folder where the visualization will be saved.

    Returns
    -------
    Visualization
        A MICOM visualization. Can be served with `viz.serve`.
    """
    rates = growth_rates
    rates = rates[rates.growth_rate > 1e-6][
        ["taxon", "sample_id", "abundance", "growth_rate"]
    ]
    data = {"growth_rates": rates}
    viz = Visualization(out_folder, data, "growth.html")
    viz.save(
        data=rates.to_json(orient="records"), width=800, height=400
    )

    return viz
