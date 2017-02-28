"""A class representing a microbial or tissue community."""

import six
import cobra
import pandas as pd
from sympy.core.singleton import S
from mico.util import load_model, fluxes_from_primals
from mico.logger import logger


_taxonomy_cols = ["id", "file"]


class Community(cobra.Model):
    """A community of models."""

    def __init__(self, taxonomy, id=None, name=None, idmap=None,
                 rel_threshold=1e-6, solver=None):
        """Constructor for the class."""
        super(Community, self).__init__(id, name)

        logger.info("building new mico model {}.".format(id))
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

        if "abundance" not in taxonomy.columns:
            taxonomy["abundance"] = 1
        taxonomy.abundance /= taxonomy.abundance.sum()
        logger.info("{} models with abundances below threshold".format(
                    (taxonomy.abundance <= self._rtol).sum()))
        taxonomy = taxonomy[taxonomy.abundance > self._rtol]

        self.taxonomy = taxonomy.set_index("id")

        obj = S.Zero
        self.objectives = {}
        for idx, row in self.taxonomy.iterrows():
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
            obj += o.expression
            self.objectives[idx] = o.expression
            self.__add_exchanges(model)

        self.objective = self.solver.interface.Objective(obj, direction="max")

    def __add_exchanges(self, model):
        """Add exchange reactions for a new model."""
        pass

    def optimize_single(self, id, fluxes=False):
        """Optimize growth rate for a single model in the community."""
        if isinstance(id, six.string_types):
            if id not in self.taxonomy.index:
                raise ValueError(id + " not in taxonomy!")
            info = self.taxonomy.loc[id]
        elif isinstance(id, int) and id >= 0 and id < len(self.taxonomy):
            info = self.taxonomy.iloc[id]
        else:
            raise ValueError("`id` must be an id or positive index!")

        logger.info("optimizing for {}".format(info.name))

        obj = self.objectives[info.name]
        with self as m:
            m.objective = obj
            m.solver.optimize()
            if fluxes:
                res = fluxes_from_primals(m, info)
            else:
                res = m.objective.value

        return res

    def optimize_all(self, fluxes=False):
        """Return solutions for individually optimizing each model."""
        individual = (self.optimize_single(id, fluxes) for id in
                      self.taxonomy.index)

        if fluxes:
            return pd.concat(individual, axis=1).T
        else:
            return pd.Series(individual, self.taxonomy.index)

    @property
    def abundances(self):
        """The relative abundances of species/tissues in the community."""
        return self.taxonomy.abundance

    @abundances.setter
    def abundances(self, value):
        """Set new abundance levels."""
        try:
            self.taxonomy.abundance = value
        except Exception:
            raise ValueError("value must be an iterable with an entry for "
                             "each species/tissue")

        ab = self.taxonomy.abundance
        self.taxonomy.abundance /= ab.sum()
        self.taxonomy.abundance[ab < self._rtol] = self._rtol
