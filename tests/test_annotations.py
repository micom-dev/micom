"""Test for pulling out annotations from models."""

from cobra.io import load_model
from .fixtures import community
from micom.annotation import annotate, annotate_metabolites_from_exchanges


def test_annotation_cobra_reactions():
    mod = load_model("textbook")
    rids = [r.id for r in mod.reactions]
    anns = annotate(rids, mod, what="reaction")
    assert "reaction" in anns.columns
    assert anns.reaction.isin(rids).all()
    assert "bigg.reaction" in anns.columns


def test_annotation_cobra_metabolites():
    mod = load_model("textbook")
    mids = [m.id for m in mod.metabolites]
    anns = annotate(mids, mod, what="metabolite")
    assert "metabolite" in anns.columns
    assert anns.metabolite.isin(mids).all()
    assert all(
        col in anns.columns for col in ["molecular_weight", "C_number", "N_number"]
    )


def test_annotation_micom_reactions(community):
    mod = community
    rids = [r.id for r in mod.reactions]
    anns = annotate(rids, mod, what="reaction")
    assert "reaction" in anns.columns
    assert "bigg.reaction" in anns.columns


def test_annotation_micom_metabolites(community):
    mod = community
    mids = [m.id for m in mod.metabolites]
    anns = annotate(mids, mod, what="metabolite")
    assert "metabolite" in anns.columns
    assert all(
        col in anns.columns for col in ["molecular_weight", "C_number", "N_number"]
    )


def test_annotations_cobra_exchanges():
    mod = load_model("textbook")
    anns = annotate_metabolites_from_exchanges(mod)
    assert anns.shape[0] == len(mod.exchanges)
    assert all(
        col in anns.columns
        for col in ["molecular_weight", "C_number", "N_number", "kegg.compound"]
    )


def test_annotations_micom_exchanges(community):
    mod = community
    anns = annotate_metabolites_from_exchanges(mod)
    assert anns.shape[0] == len(mod.exchanges) * 2
    assert all(
        col in anns.columns
        for col in ["molecular_weight", "C_number", "N_number", "kegg.compound"]
    )
