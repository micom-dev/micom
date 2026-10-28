"""Visualization for exchanges."""

from datetime import datetime
from scipy.cluster.hierarchy import linkage, leaves_list
from sklearn.preprocessing import StandardScaler
from ..batch import GrowthResults
from .core import Visualization
import numpy as np
import pandas as pd
from umap import UMAP
import logging

logger = logging.getLogger(__name__)


def plot_exchanges_per_sample(
    results: GrowthResults,
    filename: str = "sample_exchanges_%s.html" % datetime.now().strftime("%Y%m%d"),
    direction: str = "import",
    cluster: bool = True,
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
    anns = results.annotations.copy().drop_duplicates(subset="metabolite")
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
    results: GrowthResults,
    filename: str = "taxon_exchanges_%s.html" % datetime.now().strftime("%Y%m%d"),
    direction: str = "import",
    groups: pd.Series = None,
    standardize: bool = False,
    use_total_flux: bool = False,
    **umap_args,
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
    groups : pandas.Series
        Additional metadata to color exchanges. The index must correspond to the
        `sample_id` in the results and values must be categorical. The `.name` attribute
        will be used to name the groups. If not provided will color by taxon.
    standardize : bool
        Whether to standardize the fluxes per taxon. Default is False.
    use_total_flux : bool
        Whether to use the total flux (flux * abundance) or just the flux. Default is False.
    umap_args : dict
        Additional arguments passed to UMAP.

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
    ex = exchanges.pivot_table(
        values="flux",
        index=["sample_id", "taxon"],
        columns="reaction",
        fill_value=0.0,
    )

    n = exchanges.sample_id.nunique()
    if "init" not in umap_args:
        umap_args["init"] = "spectral"
    if "n_neighbors" not in umap_args:
        umap_args["n_neighbors"] = min(8, n // 2)
    if "min_dist" not in umap_args:
        umap_args["min_dist"] = 0.1
    if "metric" not in umap_args:
        umap_args["metric"] = "euclidean"

    if umap_args["metric"] in ["braycurtis", "canberra", "jaccard", "dice"]:
        mat = ex.values
    else:
        mat = np.log10(ex.values + 1e-6)

    if standardize:
        mat = StandardScaler().fit_transform(mat)

    reduced = UMAP(**umap_args).fit_transform(mat)
    reduced = pd.DataFrame(
        reduced, index=ex.index, columns=["UMAP 1", "UMAP 2"]
    ).reset_index()

    if groups is not None:
        if groups.dtype not in ["str", "object", "category", "bool"]:
            raise ValueError("Groups need to be categorical.")
        name = "group" if groups.name is None else groups.name
        reduced[name] = groups[reduced.sample_id].values
    else:
        name = "taxon"

    data = {"reduced": reduced}
    viz = Visualization(filename, data, "reduced.html")
    viz.save(data=reduced.to_json(orient="records"), color=name, width=600, height=500)

    return viz
