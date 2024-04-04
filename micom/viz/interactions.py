"""Visualizations for interactions."""

from datetime import datetime
from micom.interaction import interactions, summarize_interactions, MES
from micom.logger import logger
from micom.workflows import GrowthResults
from micom.viz.core import Visualization
import pandas as pd
import re

UNITS = {
    "flux": "mmol/[gDW·h]",
    "mass": "g/[gDW·h]",
    "C": "C/[gDW·h]",
    "N": "N/[gDW·h]"
}


def plot_focal_interactions(
    results : GrowthResults,
    taxon : str,
    filename : str = "focal_interactions_%s.html" % datetime.now().strftime("%Y%m%d"),
    kind : str = "mass"
) -> None:
    """Plot metabolic interactions between a focal taxa and all other taxa.

    This will visualize metabolic interaction between a taxon of interest (focal taxon)
    and all other taxa across all samples.

    Parameters
    ----------
    results : micom.workflows.GrowthResults
        The results returned by the `grow` workflow.
    taxon : str
        The focal taxon to use as a reference. Must be one of the taxa appearing
        in `results.growth_rates.taxon`.
    filename : str
        The HTML file where the visualization will be saved.
    kind : str
        Which kind of flux to use. Either
        - "flux": molar flux of a metabolite
        - "mass" (default): the mass flux (flux normalized by molecular weight)
        - "C": carbon flux
        - "N": nitrogen flux

    Returns
    -------
    Visualization
        A MICOM visualization. Can be served with `viz.serve`.
    """
    if not kind in UNITS:
        raise ValueError(
            f"Not a supported flux type. Please choose one of {','.join(UNITS)}.")
    ints = interactions(results, taxon, progress=False)
    n_taxa = ints.partner.nunique()
    n_mets = ints.metabolite.nunique()
    data = {"interactions": ints}
    viz = Visualization(filename, data, "focal_interactions.html")
    summ = summarize_interactions(ints)
    data["summary"] = summ
    summ["flux"] = summ["flux" if kind == "flux" else kind + "_flux"]
    viz.save(
        summary=summ.to_json(orient="records"),
        interactions=data["interactions"].to_json(orient="records"),
        n_taxa=n_taxa,
        n_mets=n_mets,
        taxon=re.sub(r"\w__", "", taxon).replace("_", " "),
        unit=UNITS[kind]
    )
    return viz


def plot_mes(
    results : GrowthResults,
    taxon : str,
    filename : str = "mes_%s.html" % datetime.now().strftime("%Y%m%d"),
    metadata : pd.DataFrame = None
) -> None:
    """Plot metabolic interactions between a focal taxa and all other taxa.

    This will visualize metabolic interaction between a taxon of interest (focal taxon)
    and all other taxa across all samples.

    Parameters
    ----------
    results : micom.workflows.GrowthResults
        The results returned by the `grow` workflow.
    taxon : str
        The focal taxon to use as a reference. Must be one of the taxa appearing
        in `results.growth_rates.taxon`.
    filename : str
        The HTML file where the visualization will be saved.
    metadata : pamdas.DataFrame
        Additional metadata to stratify MES score. Must contain a column called
        `sample_id` and other categorical columns.

    Returns
    -------
    Visualization
        A MICOM visualization. Can be served with `viz.serve`.
    """
    if metadata is not None:
        if "sample_id" not in metadata.columns:
            raise ValueError("Metadata must contain a column named `sample_id`.")
        categorical = metadata.drop(columns=["sample_id"]).select_dtypes(
            include=["object", "category", "bool"]).columns
        if len(categorical) == 0:
            raise ValueError("Metadata contains no categorical columns.")
    tol = results.exchanges.tolerance.max()
    scores = MES(results, tol)
    scores = pd.merge(scores, metadata[categorical], on="sample_id")
    n_mets = scores.metabolite.nunique()
    data = {"scores": scores}
    viz = Visualization(filename, data, "scores.html")
    viz.save(
        scores=scores.to_json(orient="records"),
        n_mets=n_mets,
        name="Metabolic Exchange Score (MES)",
        col_name="mes",
        categories=", ".join(categorical)
    )
    return viz

