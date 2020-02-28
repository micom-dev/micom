"""Simple init file for mico."""

from micom.community import Community, load_pickle
from micom import (
    algorithms,
    problems,
    util,
    data,
    duality,
    elasticity,
    media,
    solution,
    workflows,
    workflow_examples,
)


__all__ = (
    "Community",
    "algorithms",
    "problems",
    "optcom",
    "util",
    "data",
    "duality",
    "elasticity",
    "media",
    "solution",
    "load_pickle",
    "workflows",
    "workflow_examples",
)

__version__ = "0.10.4"
