"""Visualization for exchanges."""

from datetime import datetime
from fastcluster import linkage
from scipy.cluster.hierarchy import leaves_list
from micom.viz.core import Visualization
import pandas as pd


def plot_exchanges_per_sample(
    results,
    out_folder="sample_exchanges_%s" % datetime.now().strftime("%Y%m%d"),
    direction="import",
    cluster=True,
) -> None:
    """Plot the per sample exchange fluxes.

    Parameters
    ----------
    results : micom.workflows.GrowthResults
        The results returned by the `grow` workflow.
    out_folder : str
        The folder where the visualization will be saved.
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
    if direction not in ["import", "export"]:
        ValueError("Not a valid flux direction. Must be `import` or `export`.")
    exchanges = exchanges[
        (exchanges.taxon == "medium")
        & (exchanges.direction == direction)
        & (exchanges.flux.abs() > 1e-6)
    ].copy()
    exchanges.flux = exchanges.flux.abs()
    mat = exchanges.pivot_table(
        values="flux", index="metabolite", columns="sample_id", fill_value=1e-6
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
    viz = Visualization(out_folder, data, "sample_heatmap.html")
    w = mat.shape[1] * 10
    viz.save(
        data=data["exchange_fluxes"].to_json(orient="records"),
        width=w,
        height=w * mat.shape[0] / mat.shape[1],
    )
    return viz


def plot_exchanges_per_taxon(
    results,
    out_folder="taxon_exchanges_%s" % datetime.now().strftime("%Y%m%d"),
    direction="import",
    **umap_args
) -> None:
    """Plot the exchange fluxes per taxon.

    Parameters
    ----------
    results : micom.workflows.GrowthResults
        The exchanges returned by the `grow` workflow.
    out_folder : str
        The folder where the visualization will be saved.
    direction : str either "import" or "export"
        The direction of fluxes to plot.

    Returns
    -------
    Visualization
        A MICOM visualization. Can be served with `viz.serve`.

    """
    from umap import UMAP  # allows skipping tests

    exchanges = results.exchanges
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
        fill_value=0,
    )
    umapped = UMAP(**umap_args).fit_transform(mat.values)
    umapped = pd.DataFrame(
        umapped, index=mat.index, columns=["UMAP 1", "UMAP 2"]
    ).reset_index()
    data = {"umap": umapped}
    viz = Visualization(out_folder, data, "umap.html")
    viz.save(data=umapped.to_json(orient="records"), width=600, height=500)

    return viz
