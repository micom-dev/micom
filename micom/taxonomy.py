"""Helper to convert external data to a micom taxonomy."""

from micom.db import load_manifest, load_zip_model_db
from micom.community import _ranks
from micom.qiime_formats import (
    load_qiime_model_db,
    load_qiime_feature_table,
    load_qiime_taxonomy,
)
import pandas as pd
from tempfile import TemporaryDirectory


def build_from_qiime(
    abundance,
    taxonomy: pd.Series,
    manifest: pd.DataFrame,
    cutoff: float = 1e-4,
    strict: bool = True,
) -> pd.DataFrame:
    """Build the specification for the community models."""
    taxa = taxonomy.str.replace("[\\w_]+__|\\[|\\]", "")
    taxa = taxa.str.split(";\\s*", expand=True).replace("", None)
    taxa.columns = _ranks[0 : taxa.shape[1]]
    taxa["taxid"] = taxonomy.index
    taxa.index == taxa.taxid

    rank = manifest.summary_rank[0]
    if strict:
        ranks = [
            r
            for r in _ranks[0 : (_ranks.index(rank) + 1)]
            if r in taxa.columns and r in manifest.columns
        ]
    else:
        ranks = [rank]
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
    del abundance["mapping_ranks"]
    abundance.dropna(subset=ranks, inplace=True)
    depth = abundance.groupby("sample_id").abundance.sum()
    abundance["relative"] = abundance.abundance / depth[abundance.sample_id].values

    micom_taxonomy = pd.merge(manifest, abundance, on=ranks)
    micom_taxonomy = micom_taxonomy[micom_taxonomy.relative > cutoff]
    print("Taxa per sample:")
    print(micom_taxonomy.sample_id.value_counts().describe(), "\n")
    return abundance


def qiime_to_micom(feature_table, taxonomy, model_db, cutoff, strict=False):
    """Load a micom taxonomy from Qiime 2 data.

    Parameters
    ----------
    feature_table : str
        Path to a Qiime 2 FeatureTable artifact.
    taxonomy : str
        Path to a Qiime 2 FeatureData[Taxonomy] artifact.
    model_db : str
        Path to a model database folder.
    cutoff : float
        Minimum relative abundance considered present in a sample.
    strict : bool
        Whether to match taxa on all ranks or only the target taxonomy rank.

    Returns
    -------
    pd.DataFrame
        A micom taxonomy containing abundances and taxonomy calls in long
        format.
    """
    compressed = model_db.endswith(".qza") or model_db.endswith(".zip")
    if compressed:
        tdir = TemporaryDirectory(prefix="micom_")
    if model_db.endswith(".qza"):
        manifest = load_qiime_model_db(model_db, tdir.name)
    elif model_db.endswith(".zip"):
        manifest = load_zip_model_db(model_db, tdir.name)
    else:
        manifest = load_manifest(model_db)

    table = load_qiime_feature_table(feature_table)
    taxonomy = load_qiime_taxonomy(taxonomy)

    return build_from_qiime(table, taxonomy, manifest, cutoff, strict)
