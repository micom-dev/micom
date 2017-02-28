"""Submodule including some common data sets."""

from os.path import split, join
import pandas as pd

__all__ = ["agora"]
this_dir, _ = split(__file__)

agora = pd.read_csv(join(this_dir, "agora.csv"))
