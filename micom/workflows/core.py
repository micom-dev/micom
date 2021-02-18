"""Makes it easier to run analyses on several samples in parallel."""

from collections import Sized, namedtuple
from multiprocessing import Pool
from rich.progress import track

GrowthResults = namedtuple(
    "GrowthResults",
    ["growth_rates", "exchanges", "annotations"]
)


def workflow(func, args, n_jobs=4, unit="sample(s)", progress=True):
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
    n_jobs : positive int
        How many samples to analyze in parallel at once.
    unit : str
        The unit used for the progress bar.
    progress : bool
        Whether to show a progress bar.
    """
    if not isinstance(args, Sized):
        ValueError("`args` must have a length.")

    with Pool(processes=n_jobs, maxtasksperchild=1) as pool:
        it = pool.imap_unordered(func, args)
        if progress:
            it = track(it, total=len(args), description="Running")
        results = list(it)
    return results
