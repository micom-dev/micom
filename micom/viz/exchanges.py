"""Visualization for exchanges."""

from datetime import datetime
from fastcluster import linkage
from scipy.cluster.hierarchy import leaves_list
from micom.viz.core import Visualization
import pandas as pd
from sklearn.manifold import TSNE


def plot_exchanges_per_sample(
    results,
    filename="sample_exchanges_%s.html" % datetime.now().strftime("%Y%m%d"),
    direction="import",
    cluster=True,
) -> None:
    """Plot the per sample exchange fluxes.

    Parameters
    ----------
    results : micom.workflows.GrowthResults
        The results returned by the `grow` workflow.
    filename : str
        The HTML file where the visualization will be saved.
    direction : str either "import" or "export"
        The direction of fluxes to plot.
    cluster : bool
        Whether to reorder samples so that samples with similar exchange
        fluxes are close to another.

    Returns
    -------
    Visualization
        A MICOM visualization. Can be served with `viz.serve`.
    """
    exchanges = results.exchanges
    anns = results.annotations
    anns.index = anns.metabolite
    tol = exchanges.tolerance.iloc[0]
    if direction not in ["import", "export"]:
        ValueError("Not a valid flux direction. Must be `import` or `export`.")
    exchanges = exchanges[
        (exchanges.taxon == "medium")
        & (exchanges.direction == direction)
    ].copy()
    exchanges.flux = exchanges.flux.abs()
    mat = exchanges.pivot_table(
        values="flux", index="metabolite", columns="sample_id", fill_value=tol
    )
    sample_order = leaves_list(linkage(mat.values.T, method="average"))
    reaction_order = leaves_list(linkage(mat.values, method="average"))
    mat = mat.iloc[reaction_order, sample_order]

    mat["metabolite"] = mat.index
    data = mat.melt(
        id_vars="metabolite", var_name="sample_id", value_name="flux"
    )
    data["description"] = anns.loc[data.metabolite, "name"].values
    data = {"exchange_fluxes": data}
    viz = Visualization(filename, data, "sample_heatmap.html")
    w = mat.shape[1] * 10
    viz.save(
        data=data["exchange_fluxes"].to_json(orient="records"),
        width=w,
        height=w * mat.shape[0] / mat.shape[1],
    )
    return viz


def plot_exchanges_per_taxon(
    results,
    filename="taxon_exchanges_%s.html" % datetime.now().strftime("%Y%m%d"),
    direction="import",
    **tsne_args
) -> None:
    """Plot the exchange fluxes per taxon.

    Parameters
    ----------
    results : micom.workflows.GrowthResults
        The exchanges returned by the `grow` workflow.
    filename : str
        The HTML file where the visualization will be saved.
    direction : str either "import" or "export"
        The direction of fluxes to plot.
    tsne_args : dict
        Additional arguments passed to TSNE.

    Returns
    -------
    Visualization
        A MICOM visualization. Can be served with `viz.serve`.

    """
    exchanges = results.exchanges
    tol = exchanges.tolerance.iloc[0]

    if direction not in ["import", "export"]:
        ValueError("Not a valid flux direction. Must be `import` or `export`.")
    exchanges = exchanges[
        (exchanges.taxon != "medium") & (exchanges.direction == direction)
    ].copy()
    exchanges["flux"] = exchanges.flux.abs() * exchanges.abundance
    mat = exchanges.pivot_table(
        values="flux",
        index=["sample_id", "taxon"],
        columns="reaction",
        fill_value=tol,
    )
    reduced = TSNE(**tsne_args).fit_transform(mat.values)
    reduced = pd.DataFrame(
        reduced, index=mat.index, columns=["TSNE 1", "TSNE 2"]
    ).reset_index()
    data = {"reduced": reduced}
    viz = Visualization(filename, data, "reduced.html")
    viz.save(data=reduced.to_json(orient="records"), width=600, height=500)

    return viz
