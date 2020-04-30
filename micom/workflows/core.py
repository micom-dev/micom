"""Makes it easier to run analyses on several samples in parallel."""

from collections import Sized
from loky import get_reusable_executor
from tqdm import tqdm
import warnings


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

    executor = get_reusable_executor(max_workers=n_jobs, reuse=False,
                                     timeout=3600)
    it = executor.map(func, args)
    if progress:
        it = tqdm(it, total=len(args), unit=unit)
    # loky will raise a warning because processes are not reused
    # which is something we do on purpose to work around memory leaks
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        results = list(it)
    return results
