"""A class representing a microbial or tissue community."""

import re
import six
import six.moves.cPickle as pickle
import cobra
import pandas as pd
from optlang.symbolics import Zero
from tqdm.auto import tqdm
from micom.db import load_zip_model_db, load_manifest
from micom.util import (
    load_model,
    join_models,
    add_var_from_expression,
    adjust_solver_config,
    clean_ids,
)
from micom.logger import logger
from micom.optcom import optcom, solve
from micom.problems import cooperative_tradeoff, knockout_taxa
from micom.qiime_formats import load_qiime_model_db
from tempfile import TemporaryDirectory

_ranks = ["kingdom", "phylum", "class", "order", "family", "genus", "species", "strain"]

cobra.io.sbml.LOGGER.setLevel("ERROR")


class Community(cobra.Model):
    """A community of models.

    This class represents a community of individual models. It was designed for
    microbial communities but may also be used for multi-tissue or tissue-cell
    mixture models as long as all individuals exist within a single enclosing
    compartment.
    """

    def __init__(
        self,
        taxonomy,
        model_db=None,
        id=None,
        name=None,
        rel_threshold=1e-6,
        solver=None,
        progress=True,
        max_exchange=100,
        mass=1,
    ):
        """Create a new community object.

        `micom` builds a community from a taxonomy which may simply be a list
        of model files in its simplest form. Usually, the taxonomy will contain
        additional information such as annotations for the individuals (for
        instance phylum, organims or species) and abundances.

        The recommended way to build a micom model is to supply a
        quantification of taxa (called "taxonomy" here) which specifies the
        taxonomic ranks for a taxon and its abundance, and a model database
        for a specific rank (for instance "genus"). MICOM will match the
        ranks from your taxonomy to the model database and assemble the
        community models from that. You will also get information about the
        construction process by calling `Community.build_metrics`.

        The most customizable way only takes a single table where summarization
        and matching to the reference database has already occured. In this
        case you will also provide paths to model files for each taxon. This is
        the "old" way but may still be applicable if you want to use a custom
        database or want full control of matching your data to reference
        models.

        Notes
        -----
        `micom` will automatically add exchange fluxes and and a community
        objective maximizing the overall growth rate of the community.

        Parameters
        ----------
        taxonomy : pandas.DataFrame
            The taxonomy used for building the model. Must have at least the
            column "id". If no model database is specified in the next argument
            it furthermore requires a column "file" which specifies a filepath
            for each model. Valid file extensions are ".pickle", ".xml",
            ".xml.gz" and ".json". If a model database is specified this must
            contain at least a column with the same name as the rank used in
            the model database. Thus, for a genus-level database you will need
            a column `genus`. Additional taxa ranks can also be specified and
            will be used to be more stringent in taxa matching.
            Finally, the taxonomy should contain a column `abundance`. It will
            be used to quantify each individual in the community. If absent,
            MICOM will assume all individuals are present in the same amount.
        model_db : str
            A pre-built model database. If ending in `.qza` must be a Qiime 2
            artifact of type `MetabolicModels[JSON]`. Can also be a folder,
            zip (must end in `.zip`) file or None if the taxonomy contains a
            column `file`.
        id : str, optional
            The ID for the community. Should only contain letters and numbers,
            otherwise it will be formatted as such.
        name : str, optional
            The name for the community.
        rel_threshold : float < 1, optional
            The relative abundance threshold that will be used. Describes the
            smallest relative amount of an individual that will be considered
            non-zero. All individuals with a smaller relative amount will be
            omitted.
        solver : str, optional
            Which solver to use. Will default to cplex if available which is
            better suited for large problems.
        progress : bool, optional
            Show a progress bar.
        max_exchange : positive float, optional
            During model constructions exchange reactions are duplicated into
            internal and external exchange reactions. This specifies the new
            import flux bound for the *internal* exchange reaction. Import
            rates for the exchanges between the medium and outside are still
            mantained.
        mass : positive float, optional
            The total mass of the community in gDW. Used to adjust import
            fluxes which are assumed to be given as mmol/gDW*h for the
            entire community. As a consequence all import fluxes will be
            divided by that number.

        Attributes
        ----------
        taxa : list
            A list of taxa IDs in the community.

        """
        super(Community, self).__init__(id, name)

        logger.info("building new micom model {}.".format(id))
        if not solver:
            solver = [
                s
                for s in ["cplex", "osqp", "gurobi", "glpk"]
                if s in cobra.util.solver.solvers
            ][0]
        logger.info("using the %s solver." % solver)
        if solver == "glpk":
            logger.warning(
                "No QP solver found, will use GLPK. A lot of functionality "
                "in MICOM will require a QP solver :/"
            )
        self.solver.configuration.lp_method = "auto"
        self.solver.configuration.qp_method = "auto"
        self.solver.configuration.presolve = False
        self.solver = solver
        self._rtol = rel_threshold
        self._modification = None
        self.mass = mass
        self.__db_metrics = None
        adjust_solver_config(self.solver)
        taxonomy = taxonomy.copy()
        if "abundance" not in taxonomy.columns:
            taxonomy["abundance"] = 1
        taxonomy.abundance /= taxonomy.abundance.sum()
        logger.info(
            "{} individuals with abundances below threshold".format(
                (taxonomy.abundance <= self._rtol).sum()
            )
        )
        taxonomy = taxonomy[taxonomy.abundance > self._rtol]

        if not (isinstance(taxonomy, pd.DataFrame) and "id" in taxonomy.columns):
            raise ValueError(
                "`taxonomy` must be a pandas DataFrame with at" "least a column `id` :("
            )
        if model_db is None and "file" not in taxonomy.columns:
            raise ValueError(
                "If no model database is specified you need to pass "
                "file names for models in a `file` column as well."
            )
        compressed = False
        if model_db is not None:
            compressed = model_db.endswith(".qza") or model_db.endswith(".zip")
            if compressed:
                tdir = TemporaryDirectory(prefix="micom_")
            if "file" in taxonomy.columns:
                del taxonomy["file"]
            if model_db.endswith(".qza"):
                manifest = load_qiime_model_db(model_db, tdir.name)
            elif model_db.endswith(".zip"):
                manifest = load_zip_model_db(model_db, tdir.name)
            else:
                manifest = load_manifest(model_db)
            rank = manifest["summary_rank"][0]
            if rank not in taxonomy.columns:
                raise ValueError("Missing the column `%s` from the taxonomy." % rank)
            keep_cols = [
                r
                for r in _ranks[0 : (_ranks.index(rank) + 1)]
                if r in taxonomy.columns and r in manifest.columns
            ]
            manifest = manifest[keep_cols + ["file"]]
            merged = pd.merge(taxonomy, manifest, on=keep_cols)

            self.__db_metrics = pd.Series(
                {
                    "found_taxa": merged.shape[0],
                    "total_taxa": taxonomy.shape[0],
                    "found_fraction": merged.shape[0] / taxonomy.shape[0],
                    "found_abundance_fraction": merged.abundance.sum(),
                }
            )
            logger.info(
                "Matched %g%% of total abundance in model DB."
                % (100.0 * self.__db_metrics[3])
            )
            if self.__db_metrics["found_abundance_fraction"] < 0.5:
                logger.warning(
                    "Less than 50%% of the abundance could be matched to the "
                    "model database. Model `%s` may not be representative "
                    "of the sample" % self.id
                )
            taxonomy = merged
            taxonomy["abundance"] /= taxonomy["abundance"].sum()

        if taxonomy.id.str.contains(r"[^A-Za-z0-9_]", regex=True).any():
            logger.warning(
                "Taxa IDs contain prohibited characters and" " will be reformatted."
            )
            taxonomy.id = taxonomy.id.replace(
                [r"[^A-Za-z0-9_\s]", r"\s+"], ["", "_"], regex=True
            )

        self.__taxonomy = taxonomy
        self.__taxonomy.index = self.__taxonomy.id

        obj = Zero
        self.taxa = []
        index = self.__taxonomy.index
        index = tqdm(index, unit="models") if progress else index
        for idx in index:
            row = self.__taxonomy.loc[idx]
            if isinstance(row.file, list):
                if len(row.file) > 1:
                    model = join_models(row.file)
                    logger.info("joined {} models".format(len(row.file)))
                else:
                    model = load_model(row.file[0])
            else:
                model = load_model(row.file)
            suffix = "__" + idx.replace(" ", "_").strip()
            logger.info("converting IDs for {}".format(idx))
            external = cobra.medium.find_external_compartment(model)
            logger.info(
                "Identified %s as the external compartment for %s. "
                "If that is wrong you may be in trouble..." % (external, idx)
            )
            for r in model.reactions:
                r.global_id = clean_ids(r.id)
                r.id = r.global_id + suffix
                r.community_id = idx
                # avoids https://github.com/opencobra/cobrapy/issues/926
                r._compartments = None
                # SBO terms may not be maintained
                if "sbo" in r.annotation:
                    del r.annotation["sbo"]
            for m in model.metabolites:
                m.global_id = clean_ids(m.id)
                m.id = m.global_id + suffix
                m.compartment += suffix
                m.community_id = idx
            logger.info("adding reactions for {} to community".format(idx))
            self.add_reactions(model.reactions)
            o = self.solver.interface.Objective.clone(
                model.objective, model=self.solver
            )
            obj += o.expression * row.abundance
            self.taxa.append(idx)
            taxa_obj = self.problem.Constraint(
                o.expression, name="objective_" + idx, lb=0.0
            )
            self.add_cons_vars([taxa_obj])
            self.__add_exchanges(
                model.reactions,
                row,
                external_compartment=external,
                internal_exchange=max_exchange,
            )
            self.solver.update()  # to avoid dangling refs due to lazy add

        if compressed:
            tdir.cleanup()
        com_obj = add_var_from_expression(self, "community_objective", obj, lb=0)
        self.objective = self.problem.Objective(com_obj, direction="max")

    def __add_exchanges(
        self, reactions, info, external_compartment="e", internal_exchange=1000,
    ):
        """Add exchange reactions for a new model."""
        for r in reactions:
            # Some sanity checks for whether the reaction is an exchange
            ex = external_compartment + "__" + r.community_id
            if not cobra.medium.is_boundary_type(r, "exchange", ex):
                continue
            if not r.id.lower().startswith("ex"):
                logger.warning(
                    "Reaction %s seems to be an exchange " % r.id
                    + "reaction but its ID does not start with 'EX_'..."
                )

            export = len(r.reactants) == 1
            if export:
                lb = r.lower_bound / self.mass
                ub = r.upper_bound
            else:
                lb = -r.upper_bound / self.mass
                ub = -r.lower_bound
            if lb < 0.0 and lb > -1e-6:
                logger.info(
                    "lower bound for %r below numerical accuracy "
                    "-> adjusting to stabilize model."
                )
                lb = -1e-6
            if ub > 0.0 and ub < 1e-6:
                logger.info(
                    "upper bound for %r below numerical accuracy "
                    "-> adjusting to stabilize model."
                )
                ub = 1e-6
            met = (r.reactants + r.products)[0]
            old_compartment = met.compartment.replace("__" + r.community_id, "")
            medium_id = re.sub(
                "(_{}$)|([^a-zA-Z0-9 :]{}[^a-zA-Z0-9 :]$)".format(
                    old_compartment, old_compartment
                ),
                "",
                met.global_id,
            )
            medium_id += "_m"
            if medium_id == met.id:
                medium_id += "_medium"
            if medium_id not in self.metabolites:
                # If metabolite does not exist in medium add it to the model
                # and also add an exchange reaction for the medium
                logger.info("adding metabolite %s to external medium" % medium_id)
                medium_met = met.copy()
                medium_met.id = medium_id
                medium_met.compartment = "m"
                medium_met.global_id = medium_id
                medium_met.community_id = "medium"
                ex_medium = cobra.Reaction(
                    id="EX_" + medium_met.id,
                    name=medium_met.id + " medium exchange",
                    lower_bound=lb,
                    upper_bound=ub,
                )
                ex_medium.add_metabolites({medium_met: -1})
                ex_medium.global_id = ex_medium.id
                ex_medium.community_id = "medium"
                self.add_reactions([ex_medium])
            else:
                logger.info(
                    "updating import rate for external metabolite %s" % medium_id
                )
                medium_met = self.metabolites.get_by_id(medium_id)
                ex_medium = self.reactions.get_by_id("EX_" + medium_met.id)
                ex_medium.lower_bound = min(lb, ex_medium.lower_bound)
                ex_medium.upper_bound = max(ub, ex_medium.upper_bound)

            coef = info.abundance
            r.add_metabolites({medium_met: coef if export else -coef})
            if export:
                r.lower_bound = -internal_exchange
            else:
                r.upper_bound = internal_exchange

    def __update_exchanges(self):
        """Update exchanges."""
        logger.info("updating exchange reactions for %s" % self.id)
        for met in self.metabolites.query(lambda x: x.compartment == "m"):
            for r in met.reactions:
                if r.boundary:
                    continue
                coef = self.__taxonomy.loc[r.community_id, "abundance"]
                if met in r.products:
                    r.add_metabolites({met: coef}, combine=False)
                else:
                    r.add_metabolites({met: -coef}, combine=False)

    def __update_community_objective(self):
        """Update the community objective."""
        logger.info("updating the community objective for %s" % self.id)
        v = self.variables.community_objective
        const = self.constraints.community_objective_equality
        self.remove_cons_vars([const])
        com_obj = Zero
        for sp in self.taxa:
            ab = self.__taxonomy.loc[sp, "abundance"]
            taxa_obj = self.constraints["objective_" + sp]
            com_obj += ab * taxa_obj.expression
        const = self.problem.Constraint(
            (v - com_obj).expand(), lb=0, ub=0, name="community_objective_equality",
        )
        self.add_cons_vars([const])

    def optimize_single(self, id):
        """Optimize growth rate for one individual.

        `optimize_single` will calculate the maximal growth rate for one
        individual member of the community.

        Notes
        -----
        This might well mean that growth rates for all other individuals are
        low since the individual may use up all available resources.

        Parameters
        ----------
        id : str
            The ID of the individual to be optimized.
        fluxes : boolean, optional
            Whether to return all fluxes. Defaults to just returning the
            maximal growth rate.

        Returns
        -------
        float
            The maximal growth rate for the given taxa.

        """
        if isinstance(id, six.string_types):
            if id not in self.__taxonomy.index:
                raise ValueError(id + " not in taxonomy!")
            info = self.__taxonomy.loc[id]
        elif isinstance(id, int) and id >= 0 and id < len(self.__taxonomy):
            info = self.__taxonomy.iloc[id]
        else:
            raise ValueError("`id` must be an id or positive index!")

        logger.info("optimizing for {}".format(info.name))

        obj = self.constraints["objective_" + info.name]
        with self as m:
            m.objective = obj.expression
            m.solver.optimize()
            return m.objective.value

    def optimize_all(self, fluxes=False, progress=False):
        """Return solutions for individually optimizing each model.

        Notes
        -----
        This might well mean that growth rates for all other individuals are
        low since the individual may use up all available resources. As a
        consequence the reported growth rates may usually never be obtained
        all at once.

        Parameters
        ----------
        fluxes : boolean, optional
            Whether to return all fluxes. Defaults to just returning the
            maximal growth rate.
        progress : boolean, optional
            Whether to show a progress bar.

        Returns
        -------
        pandas.Series
            The maximal growth rate for each taxa.

        """
        index = self.__taxonomy.index
        if progress:
            index = tqdm(self.__taxonomy.index, unit="optimizations")

        individual = (self.optimize_single(id) for id in index)
        return pd.Series(individual, self.__taxonomy.index)

    def optimize(
        self, fluxes=False, pfba=True, raise_error=False, atol=1e-6, rtol=1e-6
    ):
        """Optimize the model using flux balance analysis.

        Parameters
        ----------
        slim : boolean, optional
            Whether to return a slim solution which does not contain fluxes,
            just growth rates.
        raise_error : boolean, optional
            Should an error be raised if the solution is not optimal. Defaults
            to False which will either return a solution with a non-optimal
            status or None if optimization fails.

        Returns
        -------
        micom.CommunitySolution
            The solution after optimization or None if there is no optimum.

        """
        with self:
            solution = solve(self, fluxes=fluxes, pfba=pfba, atol=atol, rtol=rtol)
        return solution

    @property
    def abundances(self):
        """pandas.Series: The normalized abundances.

        Setting this attribute will also trigger the appropriate updates in
        the exchange fluxes and the community objective.
        """
        return self.__taxonomy.abundance

    @abundances.setter
    def abundances(self, value):
        self.set_abundance(value, normalize=True)

    def set_abundance(self, value, normalize=True):
        """Change abundances for one or more taxa.

        Parameters
        ----------
        value : array-like object
            The new abundances. Must contain one value for each taxon. Can
            be a named object like a pandas Series.
        normalize : boolean, optional
            Whether to normalize the abundances to a total of 1.0. Many things
            in micom asssume that this is always the case. Only change this
            if you know what you are doing :O
        """
        try:
            self.__taxonomy.abundance = value
        except Exception:
            raise ValueError(
                "value must be an iterable with an entry for " "each taxa/tissue"
            )

        logger.info("setting new abundances for %s" % self.id)
        ab = self.__taxonomy.abundance
        if normalize:
            self.__taxonomy.abundance /= ab.sum()
            small = ab < self._rtol
            logger.info(
                "adjusting abundances for %s to %g"
                % (str(self.__taxonomy.index[small]), self._rtol)
            )
            self.__taxonomy.loc[small, "abundance"] = self._rtol
        self.__update_exchanges()
        self.__update_community_objective()

    @property
    def taxonomy(self):
        """pandas.DataFrame: The taxonomy used within the model.

        This attribute only returns a copy.
        """
        return self.__taxonomy.copy()

    @property
    def modification(self):
        """str: Denotes modifications to the model currently applied.

        Will be None if the community is unmodified.
        """
        return self._modification

    @modification.setter
    @cobra.util.context.resettable
    def modification(self, mod):
        self._modification = mod

    @property
    def exchanges(self):
        """list: Returns all exchange reactions in the model.

        Uses several heuristics based on the reaction name and compartments
        to exclude reactions that are *not* exchange reactions.
        """
        return cobra.medium.find_boundary_types(self, "exchange", "m")

    @property
    def build_metrics(self):
        """pd.Series: Returns general metrics for database matching.

        Only available when built using a model database.
        """
        if self.__db_metrics is not None:
            return self.__db_metrics
        else:
            raise ValueError(
                "Metrics are only available for models build using a model "
                "database :("
            )

    def optcom(self, strategy="lagrangian", min_growth=0.0, fluxes=False, pfba=True):
        """Run OptCom for the community.

        OptCom methods are a group of optimization procedures to find community
        solutions that provide a tradeoff between the cooperative community
        growth and the egoistic growth of each individual [#c1]_. `micom`
        provides several strategies that can be used to find optimal solutions:

        - "moma": Minimization of metabolic adjustment. Simultaneously
          optimizes the community objective (maximize) and the cooperativity
          cost (minimize). This method finds an exact maximum but doubles the
          number of required variables, thus being slow.
        - "lmoma": The same as "moma" only with a linear
          representation of the cooperativity cost (absolute value).
        - "original": Solves the multi-objective problem described in [#c1]_.
          Here, the community growth rate is maximized simultanously with all
          individual growth rates. Note that there are usually many
          Pareto-optimal solutions to this problem and the method will only
          give one solution. This is also the slowest method.

        Parameters
        ----------
        community : micom.Community
            The community to optimize.
        strategy : str
            The strategy used to solve the OptCom formulation. Defaults to
            "lagrangian" which gives a decent tradeoff between speed and
            correctness.
        min_growth : float or array-like
            The minimal growth rate required for each individual. May be a
            single value or an array-like object with the same length as there
            are individuals.
        fluxes : boolean
            Whether to return the fluxes as well.
        pfba : boolean
            Whether to obtain fluxes by parsimonious FBA rather than
            "classical" FBA.

        Returns
        -------
        micom.CommunitySolution
            The solution of the optimization. If fluxes==False will only
            contain the objective value, community growth rate and individual
            growth rates.

        References
        ----------
        .. [#c1] OptCom: a multi-level optimization framework for the metabolic
           modeling and analysis of microbial communities.
           Zomorrodi AR, Maranas CD. PLoS Comput Biol. 2012 Feb;8(2):e1002363.
           doi: 10.1371/journal.pcbi.1002363, PMID: 22319433

        """
        return optcom(self, strategy, min_growth, fluxes, pfba)

    def cooperative_tradeoff(
        self, min_growth=0.0, fraction=1.0, fluxes=False, pfba=True
    ):
        """Find the best tradeoff between community and individual growth.

        Finds the set of growth rates which maintian a particular community
        growth and spread up growth across all taxa as much as possible.
        This is done by minimizing the L2 norm of the growth rates with a
        minimal community growth.

        Parameters
        ----------
        min_growth : float or array-like, optional
            The minimal growth rate required for each individual. May be a
            single value or an array-like object with the same length as there
            are individuals.
        fraction : float or list of floats in [0, 1]
            The minum percentage of the community growth rate that has to be
            maintained. For instance 0.9 means maintain 90% of the maximal
            community growth rate. Defaults to 100%.
        fluxes : boolean, optional
            Whether to return the fluxes as well.
        pfba : boolean, optional
            Whether to obtain fluxes by parsimonious FBA rather than
            "classical" FBA. This is highly recommended.

        Returns
        -------
        micom.CommunitySolution or pd.Series of solutions
            The solution of the optimization. If fluxes==False will only
            contain the objective value, community growth rate and individual
            growth rates. If more than one fraction value is given will return
            a pandas Series of solutions with the fractions as indices.
        """
        return cooperative_tradeoff(self, min_growth, fraction, fluxes, pfba)

    def knockout_taxa(
        self, taxa=None, fraction=1.0, method="change", progress=True, diag=True,
    ):
        """Sequentially knowckout a list of taxa in the model.

        This uses cooperative tradeoff as optimization criterion in order to
        get unqiue solutions for individual growth rates. Requires a QP
        solver to work.

        Parameters
        ----------
        taxa : str or list of strs
            Names of taxa to be knocked out.
        fraction : float in [0, 1], optional
            Percentage of the maximum community growth rate that has to be
            maintained. Defaults to 100%.
        method : str, optional
            One of "raw", "change" or "relative change" that dictates whether
            to return the new growth rate (raw), the change in growth rate
            new - old or the relative change ([new - old] / old).
        progress : bool, optional
            Whether to show a progress bar. On by default.
        diag : bool, optional
            Whether the diagonal should contain values as well. If False will
            be filled with NaNs.

        Returns
        -------
        pandas.DataFrame
            A data frame with one row for each knockout and growth rates in the
            columns. Here the row name indicates which taxon has been knocked
            out and the columns contain the growth changes for all taxa in
            that knockout.

        """
        if taxa is None:
            taxa = self.taxa
        if isinstance(taxa, six.string_types):
            taxa = [taxa]
        if any(sp not in self.taxa for sp in taxa):
            raise ValueError(
                "At least one of the arguments is not a taxon " "in the community."
            )
        if method not in ["raw", "change", "relative change"]:
            raise ValueError(
                "`method` must be one of 'raw', 'change', " "or 'relative change'."
            )
        return knockout_taxa(self, taxa, fraction, method, progress, diag)

    @property
    def scale(self):
        """Get a scale to improve numerical properties."""
        scale = 1000.0
        if cobra.util.interface_to_str(self.problem) == "osqp":
            scale = 1.0
        return scale

    def to_pickle(self, filename):
        """Save a community in serialized form.

        Parameters
        ----------
        filename : str
            Where to save the pickled community.

        Returns
        -------
        Nothing

        """
        with open(filename, mode="wb") as out:
            pickle.dump(self, out, protocol=2)
