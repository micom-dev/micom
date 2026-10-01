"""Generate random names."""

import random
import pandas as pd
from pathlib import Path

this_dir = Path(__file__).parent.absolute()
files = (
    this_dir / "british_english_adjectives.txt",
    this_dir / "obsolete_occupations.txt",
    this_dir / "gut_microbiome_genera.txt",
)


def generate_random_name(adjective=False):
    """Generates a random name by combining an adjective, an obsolete

    profession, and a gut microbiome genus from text files using pandas.

    Parameters
    ----------
    adjective : bool, optional
        If True, include an adjective in the generated name. Default is False. This will result
        in more diverse names, but may also make them longer and more complex.

    Returns
    -------
    str
        A string in the format '[adjective]_[profession]_[genus]' where all words
        are lowercased and any non-letter characters are substituted with
        underscores.
    """
    lists = files
    if not adjective:
        lists = lists[1:]

    clean_load = lambda filepath: (
        pd.read_csv(filepath, header=None, names=["word"])["word"]
        .astype(str)
        .str.lower()
        .str.replace(r"[^a-zA-Z]", "_", regex=True)
        .str.replace(r"_+", "_", regex=True)
        .str.strip("_")
        .loc[lambda s: s != ""]
        .tolist()
    )
    words = [random.choice(clean_load(f)) for f in lists]

    return "_".join(words)
