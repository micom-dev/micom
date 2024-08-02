"""Helpers to convert external data to a MICOM taxonomy."""

from .constants import RANKS
from .qiime_formats import (
    load_qiime_feature_table,
    load_qiime_taxonomy,
)
import pandas as pd
import re


def build_from_qiime(
    abundance: pd.DataFrame,
    taxonomy: pd.Series,
    collapse_on: str = "genus",
    trim_rank_prefix: bool = False,
) -> pd.DataFrame:
    """Build the specification for the community models."""
    if trim_rank_prefix:
        taxa = taxonomy.str.replace("[\\w_]+__|\\[|\\]", "", regex=True)
    else:
        taxa = taxonomy
    taxa = taxa.str.split(";\\s*", expand=True).replace("", None)
    taxa.columns = RANKS[0 : taxa.shape[1]]
    taxa["taxid"] = taxonomy.index
    taxa.index == taxa.taxid

    if isinstance(collapse_on, str):
        collapse_on = [collapse_on]

    ranks = [r for r in collapse_on if r in taxa.columns]
    taxa["mapping_ranks"] = taxa[ranks].apply(
        lambda s: "|".join(s.astype("str")), axis=1
    )

    abundance = (
        abundance.collapse(
            lambda id_, x: taxa.loc[id_, "mapping_ranks"],
            axis="observation",
            norm=False,
        )
        .to_dataframe(dense=True)
        .T
    )
    abundance["sample_id"] = abundance.index

    abundance = abundance.melt(
        id_vars="sample_id", var_name="mapping_ranks", value_name="abundance"
    )
    abundance = pd.merge(
        abundance[abundance.abundance > 0.0],
        taxa[ranks + ["mapping_ranks"]].drop_duplicates(),
        on="mapping_ranks",
    )
    abundance["id"] = abundance["mapping_ranks"].replace(
        r"[^A-Za-z0-9_]+", "_", regex=True
    )
    del abundance["mapping_ranks"]
    abundance.dropna(subset=ranks, inplace=True)
    depth = abundance.groupby("sample_id").abundance.sum()
    abundance["relative"] = abundance.abundance / depth[abundance.sample_id].values

    return abundance


def qiime_to_micom(feature_table, taxonomy, collapse_on="genus"):
    """Load a micom taxonomy from Qiime 2 data.

    Parameters
    ----------
    feature_table : str
        Path to a Qiime 2 FeatureTable artifact.
    taxonomy : str
        Path to a Qiime 2 FeatureData[Taxonomy] artifact.
    collapse_on : str or List[str]
        The taxa ranks to collapse on. This will dictate how strict the database
        matching will be as well.

    Returns
    -------
    pd.DataFrame
        A micom taxonomy containing abundances and taxonomy calls in long
        format.
    """
    table = load_qiime_feature_table(feature_table)
    taxonomy = load_qiime_taxonomy(taxonomy)

    return build_from_qiime(table, taxonomy, collapse_on)


def rank_prefixes(manifest: pd.DataFrame) -> pd.Series:
    """Get the used prefixes for taxonomic ranks.

    Arguments
    ---------
    manifest : pandas.DataFrame
        A model database manifest.

    Returns
    -------
    pandas.Series
        The detected prefix for each taxonomic rank in the manifest.
    """
    ranks = [c for c in manifest.columns if c.lower() in RANKS]
    prefixes = pd.Series(
        {r: manifest[r].str.extract(r"^([a-z]__)").iloc[0, 0] for r in ranks}
    )

    return prefixes


def unify_rank_prefixes(taxonomy: pd.DataFrame, manifest: pd.DataFrame) -> pd.DataFrame:
    """Handle taxonomic rank prefixes in the taxonomy or database manifest.

    Arguments
    ---------
    taxonomy : pandas.DataFrame
        A taxonomy table.
    manifest : pandas.DataFrame
        A database manifest.

    Returns
    -------
    tuple of pandas.DataFrame
        The taxonomy with adjusted taxa names consistent with the database.
    """
    tax_prefixes = rank_prefixes(taxonomy)
    db_prefixes = rank_prefixes(manifest)
    ranks = tax_prefixes.index[tax_prefixes.index.isin(db_prefixes.index)]
    if all(tax_prefixes[ranks] == db_prefixes[ranks]):
        return taxonomy

    taxonomy = taxonomy.copy()
    ranks = [c for c in taxonomy.columns if c.lower() in RANKS]
    if db_prefixes.isna().all():
        for r in ranks:
            taxonomy[r] = taxonomy[r].str.replace(r"^[a-z]__", "", regex=True)
    else:
        for r in ranks:
            taxonomy[r] = db_prefixes[r] + taxonomy[r]

    return taxonomy


def taxon_id(term : str, rates : pd.DataFrame) -> str:
    """Find the ID for a taxon.

    Arguments
    ---------
    term : str
        The search term.
    rates : pandas.Dataframe
        A table of growth rates.

    Returns
    -------
    str
        The ID of the found taxon. Will raise a ValueError if not found.
    """
    if term in rates.taxon.values:
        return term

    term = re.sub(r"[^A-Za-z0-9_]+", "_", term)
    if term in rates.taxon.values:
        return term

    raise ValueError(
        f"The (cleaned) taxon name {term} is not in the data set. "
        f"Possible taxa are: {','.join(rates.taxon.unique())}."
    )