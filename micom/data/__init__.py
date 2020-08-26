"""Submodule including some common data sets."""

from os.path import split, join
from numpy.random import randint
import pandas as pd
import pickle

__all__ = ("agora", "test_taxonomy")
this_dir, _ = split(__file__)

agora = pd.read_csv(join(this_dir, "agora.csv"))
agora["file"] = agora["id"] + ".xml"

test_db = join(this_dir, "artifacts", "species_models.qza")
test_medium = join(this_dir, "artifacts", "medium.qza")


def test_taxonomy(n=4):
    """Create a simple test taxonomy.

    Parameters
    ----------
    n : positive int
        How many species to include.

    Returns
    -------
    pandas.DataFrame
        Taxonomy specification for a.

    """
    ecoli_file = join(this_dir, "e_coli_core.xml.gz")
    ids = ["Escherichia_coli_{}".format(i) for i in range(1, n + 1)]
    taxa = pd.DataFrame({"id": ids})
    taxa["genus"] = "Escherichia"
    taxa["species"] = "Escherichia coli"
    taxa["reactions"] = 95
    taxa["metabolites"] = 72
    taxa["file"] = ecoli_file
    return taxa


def test_data(n_samples=4):
    """Create a simple test data set.

    Parameters
    ----------
    n_samples : positive int
        How many samples to include.

    Returns
    -------
    pandas.DataFrame
        Taxonomy specification for the example data.

    """
    samples = ["sample_%d" % i for i in range(1, n_samples + 1)]
    data = [test_taxonomy() for s in samples]
    for i, d in enumerate(data):
        d["sample_id"] = samples[i]
        d["species"] += " " + d.index.astype("str")
    data = pd.concat(data)
    data["abundance"] = randint(1, 1000, data.shape[0])
    return(data)


def crc_results():
    """Return growth results for a realistic example data sets.

    This will give you access to the results as produced by the community tutorial
    at https://github.com/micom-dev/q2-micom/tree/master/docs.

    Returns
    -------
    GrowthResults
        The results for running the `grow` workflow with a tradeoff of 0.3.
    """
    return pickle.load(open(join(this_dir, "artifacts", "crc_results.pickle"), "rb"))


def crc_tradeoff():
    """Return tradeoff results for a realistic example data set.

    This will give you access to the tradeoffs as produced by the community tutorial
    at https://github.com/micom-dev/q2-micom/tree/master/docs.

    Returns
    -------
    pandas.DataFrame
        Growth rates for varying tradeoff values.
    """
    return pd.read_csv(join(this_dir, "artifacts", "crc_tradeoff.csv"), index_col=0)
