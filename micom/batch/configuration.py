"""Manage the configuration for a set of MICOM analyses."""

from pathlib import Path
from ruamel.yaml import YAML
from pydantic import BaseModel, Field
from typing import Any, Optional, Union, Self, Literal, Dict
import json

from ..util import pathify

symbols = {
    "simulation": "📊",
    "coupling": "🔗",
    "build": "🔨",
    "dbs": "📁",
    "media": "🍲",
}


class SimulationConfig(BaseModel):
    """Simulation parameters for the cooperative tradeoff and related methods."""

    strategy: Literal["ctFBA", "SteadyCom"] = "ctFBA"
    """The strategy to use for the cooperative tradeoff. Can be one of "ctFBA" or "SteadyCom"."""
    flux_method: Literal["minimal imports", "pFBA", "none"] = "minimal imports"
    """The method to obtain the fluxes for the community under a given growth rate. Can be one of "minimal imports", "pFBA", or "none"."""
    host_method: Literal["before", "ctFBA"] = "before"
    """The method to include the host in the simulation. Can be one of "before" or "ctFBA"."""
    host_relative: bool = True
    """Whether to use the host growth rate relative to it's own maximum."""
    host_growth: float = 0.5
    """The growth rate of the host. If host_relative is True, this must be between 0 and 1 and denotes the percentage of the host's maximum growth rate."""


class CouplingConfig(BaseModel):
    """Parameters for flux coupling to the growth rates of the individual taxa."""

    enabled: bool = True
    """Whether to enable flux coupling."""
    strategy: Literal["resource constraint", "resource coupling", "enzyme coupling"] = (
        "resource coupling"
    )
    """The strategy to use for flux coupling. Can be one of "resource constraint", "resource coupling", or "enzyme coupling"."""
    include_exchanges: bool = True
    """Whether to include the exchange reactions in the coupling."""
    constraint: float = 400.0
    """The constraint to use for the coupling. Its interpretation depends on the strategy used."""
    lower: float = 0.0
    """A small lower bound for the enzyme usage. Only used for coupled enzyme usages."""


class BuildConfig(BaseModel):
    """Parameters for building the community models."""

    cutoff: float = 1e-4
    """The cutoff for the relative abundance of taxa to include in the model."""
    force_rebuild: bool = False
    """Whether to force rebuilding the models even if they already exist in the output folder."""


class DBConfig(BaseModel):
    """Parameters for the model databases."""

    microbial: Path = Path("default://agora201_refseq216_species_1.qza")
    """The path to the microbial model database. Can be a local path or a URL."""
    host: Optional[Dict[str, Path]] = None
    """The host database. Currently a dictionary of host_id -> host_model_path. If None, no host model will be used."""
    download_location: Path = Path("databases")
    """The location where the databases will be downloaded if they are not local already."""


class MediaConfig(BaseModel):
    """Parameters for the media."""

    weights: str = None
    """The weights for the media components."""
    max_import: float = 10.0
    """The maximum added import rate for the media components."""


class Configuration(BaseModel):
    """Manage the configuration for a set of MICOM analyses."""

    tolerance: float = 1e-6
    """Tolerance for the solver."""
    solver: str = "cplex"
    """Name of the solver used for the linear and quadratic problems."""
    threads: int = 4
    """The number of parallel workers to use when running models."""
    tradeoff: float = 0.5
    """The tradeoff between community and individual growth."""
    progress: bool = True
    """Whether to show a progress bar when running models."""

    simulation: SimulationConfig = Field(default_factory=SimulationConfig)
    coupling: CouplingConfig = Field(default_factory=CouplingConfig)
    build: BuildConfig = Field(default_factory=BuildConfig)
    dbs: DBConfig = Field(default_factory=DBConfig)
    media: MediaConfig = Field(default_factory=MediaConfig)

    @pathify
    def to_json(self: Self, path: Union[Path, str]) -> None:
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
            json.dump(self.model_dump(mode="json"), f, indent=4)

    @staticmethod
    @pathify
    def from_json(path: Union[Path, str]) -> Self:
        """Load the configuration from a JSON file.

        Parameters
        ----------
        path : Union[Path, str]
            Path to the JSON file.

        Returns
        -------
        Configuration
            The configuration object.

        """
        with open(path, "r") as f:
            config_dict = json.load(f)
        return Configuration(**config_dict)

    @pathify
    def to_yaml(self: Self, path: Union[Path, str]) -> None:
        """Write the configuration to a YAML file.

        Parameters
        ----------
        path : Union[Path, str]
            Path to the YAML file.

        Returns
        -------
        Nothing.

        """
        yaml = YAML()
        with open(path, "w") as f:
            yaml.dump(self.model_dump(mode="json"), f)

    @staticmethod
    @pathify
    def from_yaml(path: Union[Path, str]) -> Self:
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

    def __repr__(self: Self) -> str:
        """Return a string representation of the configuration.

        Returns
        -------
        str
            String representation of the configuration.

        """
        s = "# Configuration\n"
        for key, value in self.model_dump().items():
            if isinstance(value, dict):
                s += f"  {symbols.get(key, '├')} {key}:\n"
                for k, v in value.items():
                    s += f"     ├ {k}: {v}\n"
            else:
                s += f"  {symbols.get(key, '├')} {key}: {value}\n"
        return s

    def __str__(self: Self) -> str:
        """Return a string representation of the configuration.

        Returns
        -------
        str
            String representation of the configuration.

        """
        return self.__repr__()

    def _repr_html_(self: Self) -> str:
        """Return an HTML representation of the configuration.

        Returns
        -------
        str
            HTML representation of the configuration.

        """
        s = "<h2>Configuration</h2>\n<ul>\n"
        for key, value in self.model_dump().items():
            if isinstance(value, dict):
                s += f"<li><h3>{symbols.get(key, '├')} {key}</h3>\n<ul>\n"
                for k, v in value.items():
                    s += f"  <li>{k}: {v}</li>\n"
                s += "</ul></li>\n"
            else:
                s += f"<li>{key}: {value}</li>\n"
        return s
