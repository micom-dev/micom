"""Submodule including some common data sets."""

from os.path import split, join
import pandas as pd

__all__ = ("agora", "test_taxonomy")
this_dir, _ = split(__file__)

agora = pd.read_csv(join(this_dir, "agora.csv"))
agora["file"] = agora["id"] + ".xml"


def test_taxonomy(n=5):
    """Create a simple test taxonomy.

    Parameters
    ----------
    n : positive int
        How many species to include, maximum 3.

    Returns
    -------
    pandas.DataFrame
        Taxonomy specification for a.

    """
    ecoli_file = join(this_dir, "e_coli_core.xml.gz")
    ids = ["Escherichia_coli_{}".format(i) for i in range(1, n+1)]
    taxa = pd.DataFrame({"id": ids})
    taxa["genus"] = "Escherichia"
    taxa["species"] = "Eschericia coli"
    taxa["reactions"] = 95
    taxa["metabolites"] = 72
    taxa["file"] = ecoli_file
    return taxa
