"""Provides a class for a microbial or tsiiue community."""

import os.path
import pandas as pd
import cobra
import cobra.io as io


class Community(cobra.Model):
    """A class holding a community of models."""

    def __init__(self, taxonomy, id=None, name=None):
        """Construct a community from a taxonomy."""
        super(self).__init__(id, name)

        if not (isinstance(taxonomy, pd.DataFrame) and
                all(["genus", "species", "file"] in taxonomy.columns)):
            raise ValueError("`taxonomy` must be a pandas DataFrame with"
                             "at least columns genus, species and file.")

        
