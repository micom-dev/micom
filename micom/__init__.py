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
<<<<<<< HEAD
    representation,
=======
    workflow_examples,
>>>>>>> 4821cdad888314f13791e0bfe2ebd8ee716f4c30
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
    "workflow_examples",
)

__version__ = "0.9.8"
