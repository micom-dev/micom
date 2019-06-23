"""Implements additional analysis algorithms for communities."""

import pandas as pd
import numpy as np
from micom.util import load_model
from micom.logger import logger


def jaccard(inclusion):
    """Calculate jaccard distances for a community."""
    logger.info(
        "calculating jaccard distance for {}x{} input matrix".format(
            *inclusion.shape
        )
    )
    jaccard = np.apply_along_axis(
        lambda a: (a & inclusion).sum(1), 1, inclusion
    )
    jaccard = jaccard / np.apply_along_axis(
        lambda a: (a | inclusion).sum(1), 1, inclusion
    )

    return 1 - jaccard


def euclidean(inclusion):
    """Calculate euclidean distances for a community."""
    logger.info(
        "calculating euclidean distance for {}x{} input matrix".format(
            *inclusion.shape
        )
    )
    euclidean = np.apply_along_axis(
        lambda a: ((a - inclusion) ** 2).sum(1), 1, inclusion
    )

    return np.sqrt(euclidean)


def reaction_matrix(files):
    """Create a matrix of reactions x models."""
    ids = []
    for f in files:
        model = load_model(f)
        ids.extend([(r.id, model.name) for r in model.reactions])
    rlist = pd.DataFrame(ids, columns=["reaction", "id"])
    rlist["value"] = 1
    rlist = rlist.pivot_table(values="value", index="id", columns="reaction")

    return rlist.fillna(0).astype(int)


def metabolic_dist(reactions, metric=jaccard):
    """Calculate the metabolic distances between all members."""
    rids = [(r.global_id, r.community_id) for r in reactions]
    rlist = pd.DataFrame(rids, columns=["reaction", "id"])
    rlist["value"] = 1
    rlist = rlist.pivot_table(values="value", index="id", columns="reaction")
    inclusion = rlist.fillna(0).astype(int).values

    dists = metric(inclusion)

    return pd.DataFrame(dists, index=rlist.index, columns=rlist.index)
