"""Manage a batch of community models and the corresponding steps and configuration."""

from pathlib import Path
import pandas as pd
import logging
import numpy as np
from typing import Any, Union, Self

from ..constants import RANKS, DIRECTION
from ..db import get_database
from .configuration import Configuration
from .core import workflow
from .build import build_and_save, _reduce_group
from .media import _fix_medium, _medium, process_medium
from .grow import _growth
from .results import GrowthResults
from .tradeoff import _tradeoff
from ..solution import OptimizationError
from ..types import check_medium, pathify, check_taxonomy
from ..qiime_formats import load_qiime_medium

logger = logging.getLogger(__name__)


class Batch(object):
    """Manage a batch of community models and the corresponding steps and configuration."""

    def __init__(
        self,
        taxonomy: pd.DataFrame,
        medium: Union[pd.DataFrame, str, Path] = None,
        config: Configuration = Configuration(),
    ) -> None:
        """Initialize the batch with a configuration.

        Parameters
        ----------
        taxonomy : pd.DataFrame
            The taxonomy of the community models.
        medium : Union[pd.DataFrame, str, Path]
            The growth medium to use for the batch. Can be a DataFrame, a path to a CSV file,
            a path to a QZA file or a download URL.
        config : Configuration
            The configuration object.

        """
        self.build_manifest = None
        self.results = None
        self.tradeoffs = None
        self.out_folder = None
        check_taxonomy(taxonomy)
        self._taxonomy = taxonomy
        self.medium = medium
        self.config = config

    @property
    def taxonomy(self: Self) -> pd.DataFrame:
        """Get the taxonomy of the batch.

        Returns
        -------
        pd.DataFrame
            The taxonomy of the community models.

        """
        return self._taxonomy.copy()

    @taxonomy.setter
    def taxonomy(self: Self, taxonomy: Any) -> None:
        """Set the taxonomy of the batch.

        This is not allowed after initialization. Please create a new Batch object instead.

        Parameters
        ----------
        taxonomy : Any
            The taxonomy of the community models.

        """
        raise AttributeError(
            "The taxonomy cannot be changed after initialization."
            " Please create a new Batch object instead."
        )

    @property
    def medium(self: Self) -> pd.DataFrame:
        """Get the growth medium for the batch.

        Returns
        -------
        pd.DataFrame
            The growth medium to use for the batch. Should contain columns
            `reaction`, `flux`, and `sample_id`. If `sample_id` is not present,
            the same medium will be used for all samples.

        """
        return self._medium

    @medium.setter
    def medium(self: Self, medium: Union[pd.DataFrame, str, Path]) -> None:
        """Set the growth medium for the batch.

        Parameters
        ----------
        medium : pd.DataFrame
            The growth medium to use for the batch. Should contain columns
            `reaction`, `flux`, and `sample_id`. If `sample_id` is not present,
            the same medium will be used for all samples.

        """
        if not isinstance(medium, pd.DataFrame):
            path = get_database(
                medium, Path(self.config.dbs.download_location), what="media"
            )
            if medium.suffix == ".qza":
                medium = load_qiime_medium(path)
            elif medium.suffix == ".csv":
                medium = pd.read_csv(path)
            else:
                raise ValueError(
                    "If the medium is a path, it must be a CSV or QZA file."
                )
        check_medium(medium)
        self._medium = medium

    @pathify
    def build(
        self: Self,
        out_folder: Union[Path, str],
    ) -> pd.DataFrame:
        """Build a series of community models.

        This is a best-practice implementation of building community models
        for several samples in parallel.

        Parameters
        ----------
        out_folder : str
            The built models and a manifest file will be written to this
            folder. Will skip existing models if the folder already exists, contains models,
            *and* if config.build.force_rebuild is False.

        Returns
        -------
        pandas.DataFrame
            The manifest for the built models. Contains taxa abundances,
            build metrics and file basenames.

        """
        if out_folder.exists():
            existing = [s.name.split(".pickle")[0] for s in out_folder.glob("*.pickle")]
            if (len(existing) > 0) and (not self.config.build.force_rebuild):
                logger.warning(
                    f"Found existing models for {len(existing)} samples. Will skip those. "
                    "set `config.build.force_rebuild = True` to rebuild all models."
                )
        else:
            out_folder.mkdir(parents=True)

        # Several checks if the taxonomy table makes sense
        tax = self.taxonomy.copy()
        conf = self.config
        sample_abundances = tax.groupby("sample_id").abundance.sum()
        if any(sample_abundances == 0):
            bad = sample_abundances.index[sample_abundances == 0]
            logger.warning(
                "The following samples sum to a zero abundance and will be excluded: "
                f"{', '.join(bad)}"
            )
            tax = tax[~tax.sample_id.isin(bad)]
        if self.config.dbs.microbial is not None:
            if "file" in tax.columns:
                logger.warning(
                    "The table includes a `file` column even though a model database "
                    "is used. Will ignore it and use the model database instead. "
                    "If you want to use the `file` column please set `conf.dbs.microbial = None`."
                )
                del tax["file"]
            db = get_database(conf.dbs.microbial, Path(conf.dbs.download_location))

        samples = tax.sample_id.unique()
        out_path = pd.Series({s: out_folder / (s + ".pickle") for s in samples})
        args = [[s, tax[tax.sample_id == s], db, out_path[s], conf] for s in samples]
        res = workflow(
            build_and_save,
            args,
            conf.threads,
            description="Assembling models",
            progress=conf.progress,
        )
        metrics = pd.concat(res)
        manifest = (
            tax.groupby("sample_id")
            .apply(_reduce_group, include_groups=False)
            .reset_index(level=0)
            .dropna(axis=1)
            .reset_index(drop=True)
        )
        manifest = manifest.loc[:, ~manifest.columns.isin(RANKS)]
        manifest["file"] = manifest.sample_id + ".pickle"
        manifest = pd.merge(manifest, metrics, on="sample_id")

        if db is not None:
            if any(manifest.found_taxa == 0):
                missing = manifest.sample_id[manifest.found_taxa == 0]
                logger.warning(
                    "The following samples had no taxon matches in the model "
                    "database and will be excluded. "
                    "We recommend to verify that the taxon names and "
                    "ranks match the database and version. "
                    f"Missing samples: {', '.join(missing)} ."
                )
                manifest = manifest[manifest.found_taxa > 0]
            frac = manifest.found_abundance_fraction
            if any((frac > 0) & (frac < 0.5)):
                low = manifest.sample_id[(frac > 0) & (frac < 0.5)]
                logger.warning(
                    "Less than 50%% of the abundance could be matched to the "
                    f"model database for these samples: {', '.join(low)} ."
                )

        manifest.to_csv(out_folder / "manifest.csv", index=False)
        self.build_manifest = manifest
        self.out_folder = out_folder
        return manifest

    def grow(self: Self) -> GrowthResults:
        """Simulate growth for a set of community models.

        Note
        ----
        The strategy `mimimal imports` can become unstable for common carbon sources since
        it will add in infeasible imports that are very small but import some high-C
        molecules. If you use it check that only components from your medium have been used
        and molecules that should be essential are indeed consumed.

        Returns
        -------
        GrowthResults
            A named tuple containing the growth rates and exchange fluxes for all
            samples/models.
        """
        if not self.is_built():
            raise ValueError(
                "The batch has not been built yet. Please run `Batch.build()` first."
            )

        if self.medium is None:
            raise ValueError(
                "No medium has been specified. Please set `batch.medium` to a valid medium first."
            )

        if self.config.simulation.strategy == "SteadyCom":
            raise NotImplementedError(
                "SteadyCom is not yet implemented in the batch workflow. Please use `config.simulation.strategy = 'ctFBA'` for now."
            )

        flux_method = self.config.simulation.flux_method
        weights = self.config.media.weights
        medium = self.medium
        tradeoff = self.config.tradeoff
        man = self.build_manifest
        samples = man.sample_id.unique()
        paths = {
            s: self.out_folder / man[man.sample_id == s].file.iloc[0] for s in samples
        }
        medium = process_medium(medium, samples)
        args = [
            [
                p,
                tradeoff,
                medium.flux[medium.sample_id == s],
                weights,
                flux_method,
                None,
                None,
                False,
            ]
            for s, p in paths.items()
        ]
        results = workflow(
            _growth,
            args,
            self.config.threads,
            description="Simulating growth",
            progress=self.config.progress,
        )
        if all([r is None for r in results]):
            raise OptimizationError(
                "All numerical optimizations failed. This indicates a problem "
                "with the solver or numerical instabilities. Check that you have "
                "CPLEX or Gurobi installed. You may also increase the abundance "
                "cutoff to create simpler models."
            )
        growth = pd.concat(r["growth"] for r in results if r is not None)
        growth = growth[growth.taxon != "medium"]
        exchanges = pd.concat(r["exchanges"] for r in results if r is not None)
        exchanges["taxon"] = exchanges.index.values
        exchanges = exchanges.melt(
            id_vars=["taxon", "sample_id", "tolerance"],
            var_name="reaction",
            value_name="flux",
        ).dropna(subset=["flux"])
        abundance = growth[["taxon", "sample_id", "abundance"]]
        exchanges = pd.merge(
            exchanges, abundance, on=["taxon", "sample_id"], how="outer"
        )
        anns = pd.concat(
            r["annotations"] for r in results if r is not None
        ).drop_duplicates(subset=["reaction"])
        anns.index = anns.reaction
        exchanges = pd.merge(exchanges, anns[["metabolite"]], on="reaction", how="left")
        exchanges["direction"] = DIRECTION[(exchanges.flux > 0.0).astype(int)].values
        exchanges = exchanges[exchanges.flux.abs() > exchanges.tolerance]

        self.results = GrowthResults(growth, exchanges, anns)

        return self.results

    def tradeoff(
        self: Self,
        tradeoffs: np.array = np.arange(0.1, 1.0 + 1e-6, 0.1),
    ) -> pd.DataFrame:
        """Run growth rate predictions for varying tradeoff values.

        Parameters
        ----------
        tradeoffs : array of floats in (0.0, 1.0]
            An array of tradeoff values to be tested. One simulation without
            a tradeoff (no cooperative tradeoff) will always be run additionally
            and will have a tradeoff of "NaN".

        Returns
        -------
        pandas.DataFrame
            The predicted growth rates.

        """
        if not self.is_built():
            raise ValueError(
                "The batch has not been built yet. Please run `Batch.build()` first."
            )
        man = self.build_manifest

        samples = man.sample_id.unique()
        paths = {
            s: self.out_folder / man[man.sample_id == s].file.iloc[0] for s in samples
        }
        if any(t < 0.0 or t > 1.0 for t in tradeoffs):
            raise ValueError("tradeoff values must between 0 and 1 :(")
        medium = process_medium(self.medium, samples)
        args = [
            [p, tradeoffs, medium.flux[medium.sample_id == s], None, None, False]
            for s, p in paths.items()
        ]
        results = workflow(
            _tradeoff,
            args,
            self.config.threads,
            description="Testing tradeoffs",
            progress=self.config.progress,
        )
        if all(r is None for r in results):
            raise OptimizationError(
                "All numerical optimizations failed. This indicates a problem "
                "with the solver or numerical instabilities. Check that you have "
                "CPLEX or Gurobi installed. You may also increase the abundance "
                "cutoff in `qiime micom build` to create simpler models or choose "
                "a more permissive solver tolerance."
            )
        results = pd.concat(results)
        if any(r is None for r in results):
            missing = set(samples) - set(results.sample_id)
            raise OptimizationError(
                "Some numerical optimizations failed. The following samples could"
                f"not be solved: {', '.join(missing)}."
            )
        self.tradeoffs = results

        return results

    def minimal_medium(
        self: Self,
        community_growth: float = 0.1,
        taxa_growth: float = 0.001,
        minimize: str = "mass",
        summarize: bool = True,
    ) -> pd.DataFrame:
        """Calculate the minimal medium for a set of community models.

        This requires specification of either the minimal community growth rate,
        a minimal taxon growth rate that has to be reachable by all taxa in the sample
        simultaneously, or a combination of both. All imports will be opened and the
        minimal medium allowing those growth rates will be returned. What exactly is being
        minimized (mass flux, carbon flux, number of components) can be specified through
        the `weights` and `minimize_components` options.

        Note
        ----
        A common usage example would be to request some realistic growth rate for the entire
        community and a very low growth rate for all taxa to ensure they are growing ("alive")
        in the medium. The returned solution comes from the medium minimization problem and
        does not have to correspond to the cooperative tradeoff solution with the same medium.

        Parameters
        ---------
        community_growth : positive float
            The minimum community-wide growth rate that has to be achieved on the created
            medium.
        taxa_growth : positive float
            The minimum growth rate achievable for each individual taxon in the community on the created medium.
        minimize: str
            What is minimized when completing the medium. The, default, "flux" will minimize the total added.
            "components" will minimize the number of added compounds in the medium (this may be *very* slow for large models).
            "mass" will minimize the added mass flux. Can also be set to any single element like "C" or "N" to minimize the added
            flux of that element. If None will use the criterion from the configuration.
        summarize: boolean
            Whether to summarize the medium across all samples. If False will
            return a medium for each sample (sample-specific medium).

        Returns
        -------
        pandas.DataFrame
            A new growth medium with the smallest amount of augmentations such
            that all communities can grow in it.

        """
        if not self.is_built():
            raise ValueError(
                "The batch has not been built yet. Please run `Batch.build()` first."
            )
        man = self.build_manifest

        samples = man.sample_id.unique()
        args = [
            (
                s,
                self.out_folder / man[man.sample_id == s].file.iloc[0],
                community_growth,
                taxa_growth,
                True if minimize == "components" else False,
                medium.flux[medium.sample_id == s],
                minimize if minimize not in ["components", "flux"] else None,
            )
            for s in samples
        ]
        results = workflow(
            _medium,
            args,
            self.config.threads,
            description="Finding minimal media",
            progress=self.config.progress,
        )
        if all(r is None for r in results):
            raise OptimizationError(
                "Could not find a growth medium that allows the specified "
                "growth rate for any sample :("
            )
        elif any(r is None for r in results):
            logger.error(
                "For some samples I could not find a medium that fulfills "
                "the growth rate requirements. Returning media only for the "
                "succesful samples."
            )
        medium = pd.concat(r["medium"] for r in results if r is not None)
        if summarize:
            medium = medium.groupby("reaction").flux.max().reset_index()
        medium["metabolite"] = medium.reaction.str.replace("EX_", "")

        return medium

    def complete_medium(
        self: Self,
        community_growth: float = 0.1,
        taxa_growth: float = 0.001,
        minimize: str = "mass",
        summarize: bool = True,
        medium: pd.DataFrame = None,
    ) -> pd.DataFrame:
        """Augment a growth medium so a community or specific taxa can grow on it.

        Note
        ----
        This will complete a growth medium for a single community/sample. For building
        growth media that work for arbitrary samples/compositions of taxa see
        `complete_db_medium` In contrast to `complete_db_medium` this will account for
        taxon-taxon interactions. However, growth rates will no longer be an emergent
        property of the simulation, because one needs to specify the community growth rate
        or growth rates for individual taxa.

        Parameters
        ---------
        community_growth : positive float
            The minimum community-wide growth rate that has to be achieved on the created
            medium.
        taxa_growth : positive float
            The minimum growth rate achievable for each individual taxon in the community on the created medium.
        minimize: str
            What is minimized when completing the medium. The, default, "flux" will minimize the total added.
            "components" will minimize the number of added compounds in the medium (this may be *very* slow for large models).
            "mass" will minimize the added mass flux. Can also be set to any single element like "C" or "N" to minimize the added
            flux of that element. If None will use the criterion from the configuration.
        summarize: boolean
            Whether to summarize the medium across all samples. If False will
            return a medium for each sample (sample-specific medium).
        medium: pd.DataFrame
            The medium to augment. If None will use the batch medium.

        Returns
        -------
        pandas.DataFrame
            A new growth medium with the smallest amount of augmentations such
            that all members of the community can grow in it.

        """
        if not self.is_built():
            raise ValueError(
                "The batch has not been built yet. Please run `Batch.build()` first."
            )
        if medium is None:
            medium = self.medium
        man = self.build_manifest

        samples = man.sample_id.unique()
        paths = {
            s: self.out_folder / man[man.sample_id == s].file.iloc[0] for s in samples
        }
        medium = process_medium(medium, samples)
        if medium.flux[medium.flux < 1e-6].any():
            medium.loc[medium.flux < 1e-6, "flux"] = 1e-6
            logger.info("Some import rates were to small and were adjusted to 1e-6.")
        args = [
            [
                s,
                p,
                community_growth,
                taxa_growth,
                self.config.media.max_import,
                True if minimize == "components" else False,
                medium.flux[medium.sample_id == s],
                minimize if minimize not in ["components", "flux"] else None,
            ]
            for s, p in paths.items()
        ]
        res = workflow(
            _fix_medium,
            args,
            threads=self.config.threads,
            description="Augmenting media",
            progress=self.config.progress,
        )
        if all(r is None for r in res):
            raise OptimizationError(
                "All optimizations failed. You may need to increase `max_import` "
                "or lower the target growth rate."
            )
        final = pd.concat(res)
        if summarize:
            final = (
                final.groupby(["reaction", "metabolite", "description"])
                .flux.max()
                .reset_index()
            )
        return final

    def is_built(self: Self) -> bool:
        """Check if the batch has been built.

        Returns
        -------
        bool
            True if the batch has been built, False otherwise.

        """
        check = (
            self.build_manifest is not None
            and self.out_folder is not None
            and self.out_folder.exists()
        )
        return check

    def has_results(self) -> bool:
        """Check if the batch has results.

        Returns
        -------
        bool
            True if the batch has results, False otherwise.

        """

        return self.results is not None

    def has_tradeoffs(self) -> bool:
        """Check if the batch has tradeoff results.

        Returns
        -------
        bool
            True if the batch has tradeoff results, False otherwise.

        """

        return self.tradeoffs is not None
