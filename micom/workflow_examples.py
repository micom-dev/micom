"""Example workflows for micom."""

from optlang.interface import OPTIMAL
import pandas as pd
from micom.util import load_model, clean_ids
from micom.workflows import workflow
from micom.media import complete_medium
from micom.logger import logger


def _fix_medium(args):
    """Get the fixed medium for a model."""
    mid, file, medium, min_growth, max_import, min_c = args
    model = load_model(file)
    for r in model.reactions:
        r.id = clean_ids(r.id)
    try:
        fixed = complete_medium(
            model,
            medium,
            min_growth=min_growth,
            max_import=max_import,
            minimize_components=min_c,
        )
    except Exception:
        fixed = medium.copy()
    if model.solver.status != OPTIMAL:
        logger.warning(
            "Can't reach the specified growth rate for model %s." % mid
        )
        fixed = medium.copy()
    fixed.name = mid
    return fixed


def fix_community_medium(
    tax,
    medium,
    min_growth=0.1,
    max_import=1,
    minimize_components=True,
    n_jobs=4,
):
    """Augment a growth medium so all community members can grow in it.

    Arguments
    ---------
    tax : pandas.Dataframe
        A taxonomy specification as passed to `micom.Community`.
    medium : pandas.Series
        A growth medium with exchange reaction IDs as index and positive
        import fluxes as values.
    min_growth : positive float
        The minimum biomass production required for growth.
    max_import : positive float
        The maximum import rate for added imports.
    minimize_components : boolean
        Whether to minimize the number of media components rather than the
        total flux.
    n_jobs: int
        The number of processes to use.

    Returns
    -------
    pandas.Series
        A new growth medium with the smallest amount of augmentations such
        that all members of the community can grow in it.

    """
    if medium[medium < 1e-6].any():
        medium[medium < 1e-6] = 1e-6
        logger.info(
            "Some import rates were to small and were adjusted to 1e-6."
        )
    args = [
        (row.id, row.file, medium, min_growth, max_import, minimize_components)
        for _, row in tax.iterrows()
    ]
    res = workflow(_fix_medium, args, n_jobs=n_jobs, unit="model(s)")
    return pd.concat(res, axis=1).max(axis=1)
