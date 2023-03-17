"""Test growth media for a model database."""

import pandas as pd
from cobra.medium import find_external_compartment
from micom.annotation import annotate_metabolites_from_exchanges
from micom.db import load_zip_model_db, load_manifest
from micom.workflows.core import workflow
from micom.workflows.media import process_medium
import micom.media as mm
from micom.logger import logger
from micom.solution import OptimizationError
from micom.util import load_model
from micom.qiime_formats import load_qiime_model_db
import re
from rich import print as rprint
from tempfile import TemporaryDirectory


def _grow(args):
    """Get the maximum growth rate under a given medium."""
    file, med = args
    mod = load_model(file)
    good = med[med.index.isin([r.id for r in mod.exchanges])]
    if len(good) == 0:
        logger.warning(
            "Could not find any reactions from the medium in `%s`. "
            "Maybe a mismatch in IDs?"
        )
    mod.medium = med[med.index.isin([r.id for r in mod.exchanges])]
    rate = mod.slim_optimize()
    return rate


def _try_complete(args):
    """Try to complete the medium for a model."""
    file, med, growth, max_import, mip, w, strict = args
    mod = load_model(file)
    exc = find_external_compartment(mod)
    if exc.startswith("C_"):  # for CARVEME models
        exc = exc[2:]
    try:
        fixed = mm.complete_medium(
            mod,
            med,
            growth,
            max_import=max_import,
            strict=strict,
            minimize_components=mip,
            weights=w,
        )
        added = sum(i not in med.index for i in fixed.index)
        added_flux = fixed.sum() - med[med.index.isin(fixed.index)].sum()
        can_grow = True
        logger.info(
            "Could grow `%s` by adding %d imports"
            "and %g additional .flux" % (file, added, added_flux)
        )
    except OptimizationError:
        fixed = pd.Series(float("nan"), index=med.index)
        added = float("nan")
        added_flux = float("nan")
        can_grow = False
        logger.info("Could not grow `%s`." % file)
    fixed.index = [
        re.sub(
            "(_{}$)|([^a-zA-Z0-9 :]{}[^a-zA-Z0-9 :]$)".format(exc, exc),
            "_m",
            rid,
        )
        for rid in fixed.index
    ]

    return (can_grow, added, added_flux, fixed)


def check_db_medium(model_db, medium, threads=1):
    """Complete a growth medium for all models in a database.

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
    threads : int >=1
        The number of parallel workers to use when building models. As a
        rule of thumb you will need around 1GB of RAM for each thread.

    Returns
    -------
    pd.DataFrame
        Returns an annotated manifest file with a column `can_grow` that tells you
        whether the model can grow on the (fixed) medium, and a column `growth_rate`
        that gives the growth rate.
    """
    medium = process_medium(medium, ["dummy"])
    medium.index = medium.global_id
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
    rprint(
        "Checking %d %s-level models on a medium with %d components."
        % (manifest.shape[0], rank, len(medium))
    )

    args = [(f, medium.flux) for f in manifest.file]
    results = workflow(_grow, args, threads)
    manifest["growth_rate"] = results
    manifest["can_grow"] = manifest.growth_rate.notna() & (manifest.growth_rate > 1e-6)

    if compressed:
        tdir.cleanup()

    return manifest


def complete_db_medium(
    model_db,
    medium,
    growth=0.001,
    max_added_import=1,
    minimize_components=False,
    weights=None,
    threads=1,
    strict=list(),
):
    """Complete a growth medium for all models in a database.

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
    growth : positive float or pandas.Series
        The minimum growth rate the model has to achieve with the (fixed) medium. If
        a Series will have a minimum growth rate for each id/taxon in the model db.
    max_added_import : positive float
        Maximum import flux for each added additional import not included in the growth
        medium. If positive will expand the medium with additional imports in order to
        fulfill the growth objective.
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
    strict : list
        Whether to match the imports in the predefined medium exactly. For reactions IDs
        listed here will not allow additional import of the components in the provided
        medium. For example, if your input medium has a flux of 10 mmol/(gDW*h) defined
        and the requested growth rate can only be fulfilled by ramping this up that
        would be allowed in non-strict mode but forbidden in strict mode. To match all
        medium components to strict mode use `strict=medium.global_id`.

    Returns
    -------
    tuple of (manifest, import fluxes)
        Returns an annotated manifest file with a column `can_grow` that tells you
        whether the model can grow on the (fixed) medium, and a column `added` that
        gives the number of added imports apart from the ones in the medium.
    """
    medium = process_medium(medium, ["dummy"])
    medium.index = medium.global_id
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
    rprint(
        "Completing %d %s-level models on a medium with %d components"
        " ([red]%d strict[/red])." % (manifest.shape[0], rank, len(medium), len(strict))
    )
    if not isinstance(growth, pd.Series):
        growth = pd.Series(growth, index=manifest.id)

    manifest.index = manifest.id
    args = [
        (
            manifest.loc[i, "file"],
            medium.flux,
            growth[i],
            max_added_import,
            minimize_components,
            weights,
            strict,
        )
        for i in manifest.index
    ]
    results = workflow(_try_complete, args, threads)
    manifest["can_grow"] = [r[0] for r in results]
    manifest["added"] = [r[1] for r in results]
    manifest["added_flux"] = [r[2] for r in results]
    imports = pd.DataFrame.from_records([r[3] for r in results]).fillna(0.0)
    imports.index = manifest.id

    if compressed:
        tdir.cleanup()

    metrics = (
        manifest.can_grow.sum(),
        manifest.added_flux.dropna().mean(),
        imports.mean().sum(),
    )
    rprint(
        "Obtained growth for %d models adding additional"
        " flux of %.2f/%.2f on average." % metrics
    )

    return (manifest, imports)


def _annotate(f):
    """Get annotation for a model."""
    mod = load_model(f)
    return annotate_metabolites_from_exchanges(mod)


def db_annotations(
    model_db,
    threads=1,
):
    """Get metabolite annotations from a model DB.

    Arguments
    ---------
    model_db : str
        A pre-built model database. If ending in `.qza` must be a Qiime 2
        artifact of type `MetabolicModels[JSON]`. Can also be a folder,
        zip (must end in `.zip`) file or None if the taxonomy contains a
        column `file`.
    threads : int >=1
        The number of parallel workers to use when building models. As a
        rule of thumb you will need around 1GB of RAM for each thread.

    Returns
    -------
    pd.DataFrame
        Annotations for all exchanged metabolites.
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
    rank = manifest["summary_rank"][0]
    logger.info(
        "Getting annotations from %d %s-level models ." % (manifest.shape[0], rank)
    )

    args = manifest.file.tolist()
    results = workflow(_annotate, args, threads)
    anns = pd.concat(results).drop_duplicates(subset=["reaction"])

    if compressed:
        tdir.cleanup()

    return anns
