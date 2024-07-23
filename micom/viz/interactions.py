"""Visualizations for interactions."""

from datetime import datetime
import json
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
    "N": "N/[gDW·h]",
}


def plot_focal_interactions(
    results: GrowthResults,
    taxon: str,
    filename: str = "focal_interactions_%s.html" % datetime.now().strftime("%Y%m%d"),
    kind: str = "mass",
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
            f"Not a supported flux type. Please choose one of {','.join(UNITS)}."
        )
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
        unit=UNITS[kind],
    )
    return viz


def plot_mes(
    results: GrowthResults,
    filename: str = "mes_%s.html" % datetime.now().strftime("%Y%m%d"),
    groups: pd.Series = None,
    prevalence: float = 0.5,
) -> None:
    """Plot metabolic interactions between a focal taxa and all other taxa.

    This will plot the metabolic exchange score across samples and metabolites.
    The metabolic exchange score (MES) is defined as the geometric mean of the
    number of producers and consumers for a given metabolite in a sample.

    $$
    MES = 2\cdot\frac{|p||c|}{|p| + |c|}
    $$

    Parameters
    ----------
    results : micom.workflows.GrowthResults
        The results returned by the `grow` workflow.
    filename : str
        The HTML file where the visualization will be saved.
    groups : pandas.Series
        Additional metadata to stratify MES score. The index must correspond to the
        `sample_id` in the results and values must be categorical. The `.name` attribute
        will be used to name the groups.
    prevalence : float in [0, 1]
        In what proportion of samples the metabolite has to have a non-zero MES to
        be shown on the plots. Can be used to show only very commonly exchanged
        metabolites.

    Notes
    -----
    The CSV files will always include all MES scores independent of the prevalence
    filter which is only used for visualization.

    Returns
    -------
    Visualization
        A MICOM visualization. Can be served with `viz.serve`.
    """
    tol = results.exchanges.tolerance.max()
    raw = MES(results, tol)
    scores = raw[raw.MES > 0]
    prev = scores.metabolite.value_counts() / scores.sample_id.nunique()
    prev = prev[prev > prevalence].index
    scores = scores[scores.metabolite.isin(prev)]

    if groups is not None:
        if groups.dtype not in ["object", "category", "bool"]:
            raise ValueError("Groups need to be categorical.")
        name = "group" if groups.name is None else groups.name
        scores[name] = groups[scores.sample_id].values
    else:
        scores["group"] = "all"
        name = "group"
    n_mets = scores.metabolite.nunique()
    data = {"scores": raw}
    viz = Visualization(filename, data, "scores.html")
    viz.save(
        scores=scores.to_json(orient="records"),
        n_mets=n_mets,
        name="Metabolic Exchange Score (MES)",
        col_name="MES",
        cat=name,
    )
    return viz
