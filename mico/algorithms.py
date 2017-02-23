"""Implements additional analysis algorithms for communities."""

import pandas as pd


def jaccard(inclusion):
    """Calculate jaccard distances for a community."""
    a = inclusion.reshape(inclusion.shape[0], 1, inclusion.shape[1])
    jaccard = (a - inclusion)
    jaccard = (a & inclusion).sum(2) / (a | inclusion).sum(2)

    return 1 - jaccard


def euclidean(inclusion):
    """Calculate euclidean distances for a community."""
    a = inclusion.reshape(inclusion.shape[0], 1, inclusion.shape[1])
    euclidean = (a - inclusion) ** 2
    euclidean = euclidean.sum(2).sqrt()

    return euclidean


def metabolic_dist(model, metric=jaccard):
    """Calculate the metabolic distances between all members."""
    rids = [(r.global_id, r.community_id) for r in model.reactions]
    rlist = pd.DataFrame(rids, columns=["reaction", "id"])
    rlist["value"] = 1
    rlist = rlist.pivot_table(values="value", index="id", columns="reaction")
    inclusion = rlist.fillna(0).astype(int).as_matrix()

    dists = metric(inclusion)

    return pd.DataFrame(dists, index=rlist.index, columns=rlist.index)
