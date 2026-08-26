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

    model_db: str = "default://agora201_refseq216_species_1.qza"
    """Path or URL to the model database."""
    tolerance: float = 1e-6
    """Tolerance for the solver."""
    solver: str = "cplex"
    """Name of the solver used for the linear and quadratic problems."""
    threads: int = 4
    """The number of parallel workers to use when running models."""
    tradeoff: float = 0.5
    """The tradeoff between community and individual growth."""

    simulation: Dict[str, Any] = field(
        default_factory=lambda: {
            "strategy": "ctFBA",
            "flux-method": "minimal imports",
            "host-method": "before",
            "host-relative": True,
            "host-growth": 0.5,
        }
    )
    """Simulation parameters for the cooperative tradeoff and related methods.

    - strategy: The strategy to use for the cooperative tradeoff. Can be one of "ctFBA" or "SteadyCom".
    - flux-method: The method to obtain the fluxes for the community under a given growth rate. Can be one of "minimal imports", "pFBA", or "none".
    - host-method: The method to include the host in the simulation. Can be one of "before" or "ctFBA".
    - host-relative: Whether to use the host growth rate relative to it's own maximum.
    - host-growth: The growth rate of the host. If host-relative is True, this must be between 0 and 1 and denotes the percentage of the host's maximum growth rate.
    """

    coupling: Dict[str, Any] = field(
        default_factory=lambda: {
            "enabled": True,
            "strategy": "resource coupling",
            "include-exchanges": True,
            "constraint": 400.0,
            "lower": 0.0,
        }
    )
    """Parameters for flux coupling to the growth rates of the individual taxa.

    - enabled: Whether to enable flux coupling.
    - strategy: The strategy to use for flux coupling. Can be one of "resource constraint", "resource coupling", or "enzyme coupling".
    - include-exchanges: Whether to include the exchange reactions in the coupling.
    - constraint: The constraint to use for the coupling. Its interpretation depends on the strategy used.
    - lower: A small lower bound for the enzyme usage. Only used for coupled enzyme usages.
    """

    build: Dict[str, Any] = field(
        default_factory=lambda: {
            "cutoff": 1e-4,
            "force-rebuild": False,
        }
    )
    """Parameters for building the community models.

    - cutoff: The cutoff for the relative abundance of taxa to include in the model.
    - force-rebuild: Whether to force rebuilding the models even if they already exist in the output folder.
    """

    dbs: Dict[str, Any] = field(
        default_factory=lambda: {
            "microbial": "default://agora201_refseq216_species_1.qza",
            "host": None,
            "download-location": "databases",
        }
    )
    """Parameters for the model databases.

    - microbial: The path to the microbial model database. Can be a local path or a URL.
    - host: The host database. Currently a dictionary of host_id -> host_model_path. If None, no host model will be used.
    - download-location: The location where the databases will be downloaded if they are not local already.
    """

    media: Dict[str, Any] = field(
        default_factory=lambda: {
            "weights": None,
            "max_import": 10.0,
        }
    )
    """Parameters for the media.

    - weights: The weights for the media components.
    - max_import: The maximum added import rate for the media components.
    """

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
