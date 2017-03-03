"""Simple init file for mico."""

from micom.community import Community
from micom import algorithms, problems, util, data


__all__ = ("Community", "algorithms", "problems", "util", "data")

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
