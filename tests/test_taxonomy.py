"""Test helper for taxonomy handling."""

from micom.data import test_taxonomy
import micom.taxonomy as mt


no_prefix = test_taxonomy()
with_prefix = no_prefix.copy()
with_prefix["genus"] = "g__" + with_prefix["genus"]
with_prefix["species"] = "s__" + with_prefix["species"]


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
