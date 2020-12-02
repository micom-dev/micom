"""Test growth media for a model database."""

import pandas as pd
from micom.db import load_zip_model_db, load_manifest
from micom.workflows.core import workflow
from micom.workflows.media import process_medium
import micom.media as mm
from micom.logger import logger
from micom.solution import OptimizationError
from micom.util import load_model
from micom.qiime_formats import load_qiime_model_db
from tempfile import TemporaryDirectory


def _try_complete(args):
    """Try to complete the medium for a model."""
    file, med, growth, max_import, mip, w = args
    mod = load_model(file)
    try:
        fixed = mm.complete_medium(
            mod, med, growth, max_import=max_import, minimize_components=mip, weights=w
        )
        added = fixed.index.apply(lambda i: i not in med.index).sum()
        can_grow = True
        logger.info("Could grow `%s` by adding %d import." % (file, added))
    except OptimizationError:
        fixed = pd.Series(float("nan"), index=med.index)
        added = float("nan")
        can_grow = False
        logger.info("Could not grow `%s`." % file)

    return (can_grow, added, fixed)


def check_db_medium(
    model_db,
    medium,
    growth=0.1,
    max_added_import=1,
    minimize_components=False,
    weights=None,
    threads=1
):
    """Check and complete a growth medium for all models in a database.

    Arguments
    ---------
    model_db : str
        A pre-built model database. If ending in `.qza` must be a Qiime 2
        artifact of type `MetabolicModels[JSON]`. Can also be a folder,
        zip (must end in `.zip`) file or None if the taxonomy contains a
        column `file`.
    medium : pd.DataFrame
        A growth medium. Must have columns "reaction" and "flux" denoting
        exchange reactions and their respective maximum flux. Can not be sample
        specific.
    growth : positive float
        The minimum growth rate the model has to achieve with the (fixed) medium.
    max_added_import : positive float
        Maximum import flux for each added additional import not included in the growth
        medium. If set to zero will simply check for growth on the given medium. If
        positive will expand the medium with additional imports in order to fulfill
        the growth objective.
    minimize_components : boolean
        Whether to minimize the number of components instead of the total
        import flux. Might be more intuitive if set to True but may also be
        slow to calculate.
    weights : str
        Will scale the fluxes by a weight factor. Can either be "mass" which will
        scale by molecular mass, a single element which will scale by
        the elemental content (for instance "C" to scale by carbon content).
        If None every metabolite will receive the same weight.
        Will be ignored if `minimize_components` is True.
    threads : int >=1
        The number of parallel workers to use when building models. As a
        rule of thumb you will need around 1GB of RAM for each thread.

    Returns
    -------
    tuple of (manifest, import fluxes)
        Returns an annotated manifest file with a column `can_grow` that tells you
        whether the model can grow on the (fixed) medium, and a column `added` that
        gives the number of added imports apart from the ones in the medium.
    """
    medium = process_medium(medium, "dummy").flux
    compressed = model_db.endswith(".qza") or model_db.endswith(".zip")
    if compressed:
        tdir = TemporaryDirectory(prefix="micom_")
    if model_db.endswith(".qza"):
        manifest = load_qiime_model_db(model_db, tdir.name)
    elif model_db.endswith(".zip"):
        manifest = load_zip_model_db(model_db, tdir.name)
    else:
        manifest = load_manifest(model_db)
    rank = manifest["summary_rank"][0]
    logger.info(
        "Checking %d %s-level models on a medium with %d components."
        % (manifest.shape[0], rank, len(medium))
    )

    args = [
        (f, medium, growth, max_added_import, minimize_components, weights)
        for f in manifest.file
    ]
    results = workflow(_try_complete, args, threads)
    manifest["can_grow"] = [r[0] for r in results]
    manifest["added"] = [r[1] for r in results]
    imports = pd.DataFrame.from_records([r[2] for r in results]).fillna(0.0)
    imports.index = manifest.id

    if compressed:
        tdir.cleanup()

    return (manifest, imports)
