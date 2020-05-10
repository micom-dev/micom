"""Helper to annotate metabolites and species."""

import pandas as pd


def flatten(d):
    """Flatten a dictionary into strings."""
    return {k: str(v) for k, v in d.items()}


def annotate(ids, community, what="reaction"):
    """Annotate a list of entities."""
    if what == "reaction":
        elems = community.reactions
    elif what == "metabolite":
        elems = community.metabolites
    else:
        raise ValueError("Invalid value for `what` :(")

    objs = elems.get_by_any(ids)
    anns = [
        pd.Series({
            what: o.global_id,
            "name": o.name,
            **flatten(o.annotation)
        })
        for o in objs
    ]
    return pd.DataFrame.from_records(anns).drop_duplicates()


def annotate_metabolites_from_exchanges(com):
    """Annotate exchange reactions by their metabolite."""
    mets = pd.DataFrame.from_records([
        {"metabolite": r.reactants[0], "mid":r.reactants[0].global_id,
         "rid": r.global_id}
        for r in com.reactions if r.id.startswith("EX_")
    ])
    anns = annotate(mets.metabolite.tolist(), com, "metabolite")
    idmap = mets[["mid", "rid"]].drop_duplicates()
    idmap.index = idmap.mid
    anns["reaction"] = idmap.loc[anns.metabolite, "rid"].values
    return anns
