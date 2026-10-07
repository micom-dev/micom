"""Constants used across the code base."""

import pandas as pd

RANKS = ["kingdom", "phylum", "class", "order", "family", "genus", "species", "strain"]

DB_URL = "https://zenodo.org/records/7739096/files/%s?download=1"
MEDIA_URL = "https://github.com/micom-dev/media/raw/refs/heads/main/media/%s"

DIRECTION = pd.Series(["import", "export"], index=[0, 1])
ARGS = {
    "none": {"fluxes": True, "pfba": False},
    "minimal imports": {"fluxes": False, "pfba": False},
    "pFBA": {"fluxes": True, "pfba": True},
}
