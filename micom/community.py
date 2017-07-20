"""A class representing a microbial or tissue community."""

import re
import six
import cobra
import pandas as pd
from sympy.core.singleton import S
from micom.util import load_model, join_models, add_var_from_expression
from micom.logger import logger
from micom.media import default_excludes
from micom.problems import optcom, solve

_taxonomy_cols = ["id", "file"]


class Community(cobra.Model):
    """A community of models.

    This class represents a community of individual models. It was designed for
    microbial communities but may also be used for multi-tissue or tissue-cell
    mixture models as long as all individuals exist within a single enclosing
    compartment.
    """

    def __init__(self, taxonomy, id=None, name=None, rel_threshold=1e-6,
                 solver=None):
        """Create a new community object.

        `micom` builds a community from a taxonomy which may simply be a list
        of model files in its simplest form. Usually, the taxonomy will contain
        additional information such as annotations for the individuals (for
        instance phylum, organims or species) and abundances.

        Notes
        -----
        `micom` will automatically add exchange fluxes and and a community
        objective maximizing the overall growth rate of the community.

        Parameters
        ----------
        taxonomy : pandas.DataFrame
            The taxonomy used for building the model. Must have at least the
            two columns "id" and "file" which specify an ID and the filepath
            for each model. Valid file extensions are ".pickle", ".xml",
            ".xml.gz" and ".json". If the taxonomy includes a column named
            "abundance" it will be used to quantify each individual in the
            community. If absent `micom` will assume all individuals are
            present in the same amount.
        id : str, optional
            The ID for the community.
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

        Attributes
        ----------
        objectives : dict
            A dict of {id: sympy_expression} denoting the individual growth
            objectives for each model in the community.

        """
        super(Community, self).__init__(id, name)

        logger.info("building new micom model {}.".format(id))
        if not solver:
            self.solver = ("cplex" if "cplex" in cobra.util.solver.solvers
                           else "glpk")
        else:
            self.solver = solver

        if not (isinstance(taxonomy, pd.DataFrame) and
                all(col in taxonomy.columns for col in _taxonomy_cols)):
            raise ValueError("`taxonomy` must be a pandas DataFrame with at"
                             "least columns id and file :(")

        self._rtol = rel_threshold
        self._modification = None

        taxonomy = taxonomy.copy()
        if "abundance" not in taxonomy.columns:
            taxonomy["abundance"] = 1
        taxonomy.abundance /= taxonomy.abundance.sum()
        logger.info("{} individuals with abundances below threshold".format(
                    (taxonomy.abundance <= self._rtol).sum()))
        taxonomy = taxonomy[taxonomy.abundance > self._rtol]

        self.__taxonomy = taxonomy.copy()
        self.__taxonomy.index = self.__taxonomy.id

        obj = S.Zero
        self.objectives = {}
        for idx, row in self.__taxonomy.iterrows():
            if isinstance(row.file, list):
                model = join_models(row.file)
                if len(row.file) > 1:
                    logger.info("joined {} models".format(len(row.file)))
            else:
                model = load_model(row.file)
            suffix = "__" + idx.replace(" ", "_").strip()
            logger.info("converting IDs for {}".format(idx))
            for r in model.reactions:
                r.global_id = r.id
                r.id += suffix
                r.community_id = idx
            for m in model.metabolites:
                m.global_id = m.id
                m.id += suffix
                m.compartment += suffix
                m.community_id = idx
            logger.info("adding reactions for {} to community".format(idx))
            self.add_reactions(model.reactions)
            o = self.solver.interface.Objective.clone(model.objective,
                                                      model=self.solver)
            obj += o.expression * row.abundance
            self.objectives[idx] = o.expression
            species_obj = self.problem.Constraint(
                o.expression, name="objective_" + idx, lb=0.0)
            self.add_cons_vars([species_obj])
            self.__add_exchanges(model.reactions, row)

        com_obj = add_var_from_expression(self, "community_objective",
                                          obj, lb=0)
        self.objective = self.problem.Objective(com_obj, direction="max")

    def __add_exchanges(self, reactions, info, exclude=default_excludes,
                        external_compartment="e"):
        """Add exchange reactions for a new model."""
        for r in reactions:
            # Some sanity checks for whether the reaction is an exchange
            ex = external_compartment + "__" + r.community_id
            if (not r.boundary or any(ex in r.id for ex in exclude) or
                    ex not in r.compartments):
                continue
            if not r.id.lower().startswith("ex"):
                logger.warning(
                    "Reaction %s seems to be an exchange " % r.id +
                    "reaction but its ID does not start with 'EX_'...")

            export = len(r.reactants) == 1
            lb, ub = r.bounds if export else (-r.upper_bound, -r.lower_bound)
            met = (r.reactants + r.products)[0]
            medium_id = re.sub("_{}$".format(met.compartment), "", met.id)
            if medium_id in exclude:
                continue
            medium_id += "_m"
            if medium_id == met.id:
                medium_id += "_medium"
            if medium_id not in self.metabolites:
                # If metabolite does not exist in medium add it to the model
                # and also add an exchange reaction for the medium
                medium_met = met.copy()
                medium_met.id = medium_id
                medium_met.compartment = "m"
                medium_met.global_id = medium_id
                medium_met.community_id = "medium"
                ex_medium = cobra.Reaction(
                    id="EX_" + medium_met.id,
                    name=medium_met.id + " medium exchange",
                    lower_bound=lb,
                    upper_bound=ub)
                ex_medium.add_metabolites({medium_met: -1})
                ex_medium.global_id = ex_medium.id
                ex_medium.community_id = "medium"
                self.add_reactions([ex_medium])
            else:
                medium_met = self.metabolites.get_by_id(medium_id)
                ex_medium = self.reactions.get_by_id("EX_" + medium_met.id)
                ex_medium.lower_bound = min(lb, ex_medium.lower_bound)
                ex_medium.upper_bound = max(ub, ex_medium.upper_bound)

            coef = info.abundance
            r.add_metabolites({medium_met: coef if export else -coef})

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
        com_obj = S.Zero
        for sp, expr in self.objectives.items():
            ab = self.__taxonomy.loc[sp, "abundance"]
            com_obj += ab * expr
        const = self.problem.Constraint(v - com_obj, lb=0, ub=0,
                                        name="community_objective_equality")
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
            The maximal growth rate for the given species.

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

        obj = self.objectives[info.name]
        with self as m:
            m.objective = obj
            m.solver.optimize()
            return m.objective.value

    def optimize_all(self, fluxes=False):
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

        Returns
        -------
        pandas.Series
            The maximal growth rate for each species.

        """
        individual = (self.optimize_single(id) for id in
                      self.__taxonomy.index)

        return pd.Series(individual, self.__taxonomy.index)

    def optimize(self, slim=True):
        """Optimize the model using flux balance analysis.

        Parameters
        ----------
        slim : boolean
            Whether to return a slim solution which does not contain fluxes,
            just growth rates.

        Returns
        -------
        micom.CommunitySolution
            The solution after optimization or None if there is no optimum.

        """
        self.solver.optimize()
        with self:
            solution = solve(self, fluxes=not slim)
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
        try:
            self.__taxonomy.abundance = value
        except Exception:
            raise ValueError("value must be an iterable with an entry for "
                             "each species/tissue")

        logger.info("setting new abundances for %s" % self.id)
        ab = self.__taxonomy.abundance
        self.__taxonomy.abundance /= ab.sum()
        small = ab < self._rtol
        logger.info("adjusting abundances for %s to %g" %
                    (str(self.__taxonomy.index[small]), self._rtol))
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
        return self.reactions.query(
            lambda x: x.boundary and not
            any(ex in x.id for ex in default_excludes) and
            "m" in x.compartments)

    def optcom(self, strategy="lagrangian", min_growth=0.1, tradeoff=0.5,
               fluxes=False, pfba=True):
        """Run OptCom for the community.

        OptCom methods are a group of optimization procedures to find community
        solutions that provide a tradeoff between the cooperative community
        growth and the egoistic growth of each individual [#c1]_. `micom`
        provides several strategies that can be used to find optimal solutions:

        - "linear": Applies a lower bound for the individual growth rates and
          finds the optimal community growth rate. This is the fastest methods
          but also ignores that individuals might strive to optimize their
          individual growth instead of community growth.
        - "lagrangian": Optimizes a joint objective containing the community
          objective (maximized) as well as a cooperativity cost which
          represents the  distance to the individuals "egoistic" maximum growth
          rate (minimized). Requires the `tradeoff` parameter. This method is
          still relatively fast and does require only few additional variables.
        - "linear lagrangian": The same as "lagrangian" only with a linear
          representation of the cooperativity cost (absolute value).
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
        tradeoff : float in [0, 1]
            Only used for lagrangian strategies. Must be between 0 and 1 and
            describes the strength of the cooperativity cost / egoism. 1 means
            optimization will only minimize the cooperativity cost and zero
            means optimization will only maximize the community objective.
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
        return optcom(self, strategy, min_growth, tradeoff, fluxes, pfba)
