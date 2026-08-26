"""Manage the configuration for a set of MICOM analyses."""

from dataclasses import dataclass, field
from pathlib import Path
from ruamel.yaml import YAML
from typing import Any, Dict
from zipfile import Path
import json


@dataclass
class Configuration:
    """Manage the configuration for a set of MICOM analyses."""

    model_db : Path
    """Path to the model database."""
    tolerance : float = 1e-6
    """Tolerance for the solver."""
    solver : str = "cplex"
    """Name of the solver used for the linear and quadratic problems."""
    threads : int = 4
    """The number of parallel workers to use when running models."""
    tradeoff : float = 0.5
    """The tradeoff between community and individual growth."""

    simulation : Dict[str, Any] = field(default_factory=lambda: {
        "strategy": "ctFBA",
        "flux-method": "minimal imports",
        "host-method": "before",
        "host-relative": True,
        "host-growth": 0.5,
    })

    coupling : Dict[str, Any] = field(default_factory=lambda: {
        "enabled": True,
        "strategy": "resource coupling",
        "include-exchanges": True,
        "constraint": 400.0,
        "lower": 0.0,
    })

    build : Dict[str, Any] = field(default_factory=lambda: {
        "cutoff": 1e-4,
        "force-rebuild": False,
    })

    dbs : Dict[str, Any] = field(default_factory=lambda: {
        "microbial": "default://agora201_refseq216_species_1.qza",
        "host": None
        "download-location": Path("databases"),
    })

    media : Dict[str, Any] = field(default_factory=lambda: {
        "weights": None,
        "max_import": 10.0,
    })

    def to_json(self, path: Path) -> Dict[str, Any]:
        """Write the configuration to a JSON file.

        Parameters
        ----------
        path : Path
            Path to the JSON file.

        Returns
        -------
        Nothing.

        """
        with open(path, "w") as f:
            json.dump(self.__dict__, f, indent=4)

    @staticmethod
    def from_json(path: Path) -> "Configuration":
        """Load the configuration from a JSON file.

        Parameters
        ----------
        path : Path
            Path to the JSON file.

        Returns
        -------
        Configuration
            The configuration object.

        """
        with open(path, "r") as f:
            config_dict = json.load(f)
        return Configuration(**config_dict)

    def to_yaml(self, path: Path) -> Dict[str, Any]:
        """Write the configuration to a YAML file.

        Parameters
        ----------
        path : Path
            Path to the YAML file.

        Returns
        -------
        Nothing.

        """
        yaml = YAML()
        with open(path, "w") as f:
            yaml.dump(self.__dict__, f)

    @staticmethod
    def from_yaml(path: Path) -> "Configuration":
        """Load the configuration from a YAML file.

        Parameters
        ----------
        path : Path
            Path to the YAML file.

        Returns
        -------
        Configuration
            The configuration object.

        """
        yaml = YAML()
        with open(path, "r") as f:
            config_dict = yaml.load(f)
        return Configuration(**config_dict)






