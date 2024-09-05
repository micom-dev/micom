"""Simple init file for micom."""

from micom.community import Community
from micom.deps import show_versions
from micom.logger import logger
from micom.util import load_pickle
from micom import (
    algorithms,
    problems,
    util,
    data,
    duality,
    elasticity,
    media,
    qiime_formats,
    solution,
    workflows,
    interaction,
)


__all__ = (
    "Community",
    "algorithms",
    "db",
    "problems",
    "optcom",
    "util",
    "data",
    "duality",
    "elasticity",
    "media",
    "qiime_formats",
    "solution",
    "load_pickle",
    "logger",
    "workflows",
    "show_versions",
)

__version__ = "0.36.4"
