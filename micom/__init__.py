"""Simple init file for micom."""

from importlib.metadata import PackageNotFoundError, version as _distribution_version

from .batch import Batch, Configuration
from .community import Community
from .deps import show_versions
from .util import load_pickle
from . import (
    algorithms,
    problems,
    util,
    data,
    duality,
    elasticity,
    media,
    qiime_formats,
    solution,
    batch,
    interaction,
    names,
)

__all__ = (
    "Community",
    "batch",
    "Batch",
    "Configuration",
    "algorithms",
    "db",
    "problems",
    "util",
    "data",
    "duality",
    "elasticity",
    "interaction",
    "media",
    "qiime_formats",
    "names",
    "solution",
    "load_pickle",
    "logger",
    "workflows",
    "show_versions",
)

try:
    __version__ = _distribution_version("micom")
except PackageNotFoundError:
    __version__ = "0+unknown"
