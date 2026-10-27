"""Test helper for taxonomy handling."""

import pytest
import pandas as pd

from micom.data import test_taxonomy, test_data
import micom.taxonomy as mt
from micom.types import check_taxonomy


no_prefix = test_taxonomy()
with_prefix = no_prefix.copy()
with_prefix["genus"] = "g__" + with_prefix["genus"]
with_prefix["species"] = "s__" + with_prefix["species"]
with_prefix["strain"] = "t__" + with_prefix["strain"]


def test_get_prefixes():
    assert mt.rank_prefixes(no_prefix).isna().all()
    assert mt.rank_prefixes(with_prefix)["species"] == "s__"


def test_unify_all_good():
    tax = mt.unify_rank_prefixes(no_prefix, no_prefix)
    assert all(tax.species == no_prefix.species)

    tax = mt.unify_rank_prefixes(with_prefix, with_prefix)
    assert all(tax.species == with_prefix.species)


def test_remove_prefix():
    tax = mt.unify_rank_prefixes(with_prefix, no_prefix)
    assert all(tax.genus == no_prefix.genus)
    assert all(tax.species == no_prefix.species)


def test_add_prefix():
    tax = mt.unify_rank_prefixes(no_prefix, with_prefix)
    assert all(tax.genus == with_prefix.genus)
    assert all(tax.species == with_prefix.species)


def test_check_taxonomy():
    tax = test_taxonomy()
    check_taxonomy(no_prefix, samples=False)

    with pytest.raises(ValueError):
        check_taxonomy(tax, samples=True)

    del tax["abundance"]
    with pytest.raises(ValueError):
        check_taxonomy(tax, samples=False)

    tax = test_data()
    check_taxonomy(tax, samples=True)


def test_build_from_qiime_collapses_and_normalizes():
    class FeatureTable:
        def __init__(self):
            self.mapping = None
            self.data = pd.DataFrame(
                [[2.0, 0.0], [0.0, 3.0]],
                index=["feature_a", "feature_b"],
                columns=["sample_a", "sample_b"],
            )

        def collapse(self, mapping, axis, norm):
            self.mapping = mapping
            assert axis == "observation"
            assert not norm
            labels = [mapping(feature, None) for feature in self.data.index]
            self.data = self.data.groupby(labels).sum()
            return self

        def to_dataframe(self, dense):
            assert dense
            return self.data

    taxonomy = pd.Series(
        [
            "k__Bacteria;p__Firmicutes;c__Clostridia;o__Order;f__Family;"
            "g__Genus_a;s__species_a;t__strain_a",
            "k__Bacteria;p__Firmicutes;c__Clostridia;o__Order;f__Family;"
            "g__Genus_b;s__species_b;t__strain_b",
        ],
        index=["feature_a", "feature_b"],
    )
    table = FeatureTable()
    result = mt.build_from_qiime(table, taxonomy, trim_rank_prefix=True)

    assert table.mapping("feature_a", None) == "Genus_a"
    assert result[["sample_id", "genus", "id"]].values.tolist() == [
        ["sample_a", "Genus_a", "Genus_a"],
        ["sample_b", "Genus_b", "Genus_b"],
    ]
    assert result.relative.tolist() == [1.0, 1.0]

    table = FeatureTable()
    result = mt.build_from_qiime(
        table, taxonomy, collapse_on=["genus", "species"], trim_rank_prefix=True
    )
    assert table.mapping("feature_b", None) == "Genus_b|species_b"
    assert result.species.tolist() == ["species_a", "species_b"]
