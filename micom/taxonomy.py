"""Helpers to convert external data to a MICOM taxonomy."""

from micom.community import _ranks
from micom.qiime_formats import (
    load_qiime_feature_table,
    load_qiime_taxonomy,
)
import pandas as pd


def build_from_qiime(
    abundance,
    taxonomy: pd.Series,
    collapse_on="genus"
) -> pd.DataFrame:
    """Build the specification for the community models."""
    taxa = taxonomy.str.replace("[\\w_]+__|\\[|\\]", "", regex=True)
    taxa = taxa.str.split(";\\s*", expand=True).replace("", None)
    taxa.columns = _ranks[0 : taxa.shape[1]]
    taxa["taxid"] = taxonomy.index
    taxa.index == taxa.taxid

    if isinstance(collapse_on, str):
        collapse_on = [collapse_on]

    ranks = [
        r
        for r in collapse_on
        if r in taxa.columns
    ]
    taxa["mapping_ranks"] = taxa[ranks].apply(
        lambda s: "|".join(s.astype("str")), axis=1
    )

    abundance = (
        abundance.collapse(
            lambda id_, x: taxa.loc[id_, "mapping_ranks"],
            axis="observation", norm=False
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
        on="mapping_ranks"
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
