"""Simple init file for mico."""

from micom.community import Community
from micom import algorithms, problems, util, data, duality, media, solution
from micom._version import get_versions


__all__ = ("Community", "algorithms", "problems", "util", "data",
           "duality", "media", "solution")

__version__ = get_versions()['version']
del get_versions
