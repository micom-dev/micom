"""Visualization for exchanges."""

from datetime import datetime
from scipy.cluster.hierarchy import linkage, leaves_list
from micom.logger import logger
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
        & (exchanges.flux.abs() > tol)
    ].copy()
    exchanges.flux = exchanges.flux.abs()
    mat = exchanges.pivot_table(
        values="flux", index="metabolite", columns="sample_id", fill_value=tol
    )
    if cluster:
        sample_order = leaves_list(linkage(mat.values.T, method="average"))
    else:
        sample_order = range(mat.shape[1])
    reaction_order = leaves_list(linkage(mat.values, method="average"))
    mat = mat.iloc[reaction_order, sample_order]

    mat["metabolite"] = mat.index
    data = mat.melt(id_vars="metabolite", var_name="sample_id", value_name="flux")
    data["description"] = anns.loc[data.metabolite, "name"].values
    data = {"exchange_fluxes": data}
    viz = Visualization(filename, data, "sample_heatmap.html")
    long = mat.shape[0] > mat.shape[1]
    w = mat.shape[1] * 10 if long else mat.shape[0] * 10
    height = mat.shape[0] * 10 if long else mat.shape[1] * 10
    viz.save(
        data=data["exchange_fluxes"].to_json(orient="records"),
        width=w,
        height=height,
        long=long,
    )
    return viz


def plot_exchanges_per_taxon(
    results,
    filename="taxon_exchanges_%s.html" % datetime.now().strftime("%Y%m%d"),
    direction="import",
    use_total_flux=False,
    **tsne_args,
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
    use_total_fluxes : bool
        Whether to use fluxes normalized to 1gDW of bacteria or the total flux.
    tsne_args : dict
        Additional arguments passed to TSNE.

    Returns
    -------
    Visualization
        A MICOM visualization. Can be served with `viz.serve`.

    """
    exchanges = results.exchanges

    if direction not in ["import", "export"]:
        ValueError("Not a valid flux direction. Must be `import` or `export`.")
    exchanges = exchanges[
        (exchanges.taxon != "medium") & (exchanges.direction == direction)
    ].copy()
    if use_total_flux:
        exchanges["flux"] = exchanges.flux.abs() * exchanges.abundance
    else:
        exchanges["flux"] = exchanges.flux.abs()
    mat = exchanges.pivot_table(
        values="flux",
        index=["sample_id", "taxon"],
        columns="reaction",
        fill_value=0.0,
    )

    n = exchanges.sample_id.nunique()
    if "init" not in tsne_args:
        tsne_args["init"] = "pca"
    if "learning_rate" not in tsne_args:
        tsne_args["learning_rate"] = 200.0
    if "perplexity" not in tsne_args and n <= 30:
        logger.warn(f"Not enough samples. Adjusting T-SNE perplexity to {n // 2}.")
        tsne_args["perplexity"] = n // 2

    reduced = TSNE(**tsne_args).fit_transform(mat.values)
    reduced = pd.DataFrame(
        reduced, index=mat.index, columns=["TSNE 1", "TSNE 2"]
    ).reset_index()
    data = {"reduced": reduced}
    viz = Visualization(filename, data, "reduced.html")
    viz.save(data=reduced.to_json(orient="records"), width=600, height=500)

    return viz
