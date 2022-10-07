"""Makes it easier to run analyses on several samples in parallel."""

from collections import namedtuple, abc
from multiprocessing import Pool
import pandas as pd
from rich.progress import track
import warnings
from zipfile import ZipFile

GrowthResults = namedtuple(
    "GrowthResults", ["growth_rates", "exchanges", "annotations"]
)


def save_results(results, path):
    """Save growth results to a file.

    This will write all tables as CSV into a single ZIP file.

    Arguments
    ---------
    results : GrowthResults
        The results as returned from `grow`.
    path : str
        A filepath for the generated file. Should end in `.zip`.
    """
    with ZipFile(path, "w") as zippy:
        for attr in ["growth_rates", "exchanges", "annotations"]:
            getattr(results, attr).to_csv(zippy.open(f"{attr}.csv", "w"), index=False)


def load_results(path):
    """Load growth results from a file.

    Arguments
    ---------
    path : str
        Path to saved `GrowthResults`.

    Returns
    -------
    GrowthResults
        The saved GrowthResults.
    """
    tables = []
    with ZipFile(path, "r") as zippy:
        for attr in ["growth_rates", "exchanges", "annotations"]:
            tab = pd.read_csv(zippy.open(f"{attr}.csv", "r"))
            tables.append(tab)
    return GrowthResults(*tables)


def workflow(func, args, threads=4, description=None, progress=True):
    """Run analyses for several samples in parallel.

    This will analyze several samples in parallel. Includes a workaround for
    optlang memory leak.

    Arguments
    ---------
    func : function
        A function that takes a single argument (can be any object) and
        that performs your analysis for a single sample.
    args : array-like object
        An array-like object (list, tuple, numpy array, pandas Series, etc.)
        that contains the arguments for each sample.
    threads : positive int
        How many samples to analyze in parallel at once.
    description : str
        The dewscription shown in front of the progress bar.
    progress : bool
        Whether to show a progress bar.
    """
    if not isinstance(args, abc.Sized):
        ValueError("`args` must have a length.")
    if description is None:
        description = "Running"

    # Don't generate overhead if single thread
    if threads == 1:
        it = map(func, args)
        if progress:
            it = track(it, total=len(args), description=description)
        return list(it)

    # We don't use the context  manager because of
    # https://pytest-cov.readthedocs.io/en/latest/subprocess-support.html
    pool = Pool(processes=threads, maxtasksperchild=1)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            it = pool.imap_unordered(func, args)
            if progress:
                it = track(it, total=len(args), description="Running")
            results = list(it)
    finally:
        pool.close()
        pool.join()
    return results
