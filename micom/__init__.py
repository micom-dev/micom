"""Simple init file for mico."""

from micom.community import Community, load_pickle
from micom import (
    algorithms,
    problems,
    util,
    data,
    duality,
    media,
    solution,
    workflows,
    representation,
)


__all__ = (
    "Community",
    "algorithms",
    "problems",
    "optcom",
    "util",
    "data",
    "duality",
    "media",
    "solution",
    "load_pickle",
    "workflows",
    "representation",
)

__version__ = "0.9.4"
