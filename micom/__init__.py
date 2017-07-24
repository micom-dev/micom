"""Simple init file for mico."""

from micom.community import Community, load_pickle
from micom import algorithms, problems, util, data, duality, media, solution
from micom._version import get_versions


__all__ = ("Community", "algorithms", "problems", "util", "data",
           "duality", "media", "solution", "load_pickle")

__version__ = get_versions()['version']
del get_versions
