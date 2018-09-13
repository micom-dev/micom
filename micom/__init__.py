"""Simple init file for mico."""

from micom.community import Community, load_pickle
from micom import (algorithms, problems, util, data, duality, media,
                   solution, workflows)
from micom._version import get_versions


__all__ = ("Community", "algorithms", "problems", "optcom", "util", "data",
           "duality", "media", "solution", "load_pickle", "workflows")

__version__ = get_versions()['version']
del get_versions
