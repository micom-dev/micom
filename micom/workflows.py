"""Makes it easier to run analyses on several samples in parallel."""

from collections import Sized
from tqdm import tqdm
from multiprocessing import Process, Queue


def _process(f, args, queue):
    res = f(args)
    queue.put(res)


def _consume(processes, queue, results, max_procs):
    for p in processes:
        results.append(queue.get())
    for p in processes:
        p.join()
    processes = []
    return results, processes


def workflow(func, args, n_jobs=4, unit="sample(s)"):
    if not isinstance(args, Sized):
        ValueError("`args` must have a length.")

    results = []
    processes = []
    q = Queue()

    for arg in tqdm(args, unit=unit):
        p = Process(target=_process, args=(func, arg, q))
        p.start()
        processes.append(p)
        if len(processes) >= n_jobs:
            results, processes = _consume(processes, q, results, n_jobs)
    results, processes = _consume(processes, q, results, n_jobs)
    return results
