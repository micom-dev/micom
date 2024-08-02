"""Makes it easier to run analyses on several samples in parallel."""

import logging
from collections import abc
from multiprocessing import get_context
from rich.progress import track
import warnings

logger = logging.getLogger(__name__)

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
    logger.setLevel("ERROR")

    # Don't generate overhead if single thread
    if threads == 1:
        it = map(func, args)
        if progress:
            it = track(it, total=len(args), description=description)
        return list(it)

    # We don't use the context  manager because of
    # https://pytest-cov.readthedocs.io/en/latest/subprocess-support.html
    pool = get_context("spawn").Pool(processes=threads, maxtasksperchild=1)
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

    logger.setLevel("WARNING")
    return results