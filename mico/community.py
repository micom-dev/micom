"""A class representing a microbial or tissue community."""

import six
import cobra
import pandas as pd
from sympy.core.singleton import S
from mico.util import load_model, fluxes_from_primals
from mico.logger import logger
import optlang


_taxonomy_cols = ["id", "file"]


class Community(cobra.Model):
    """A community of models."""

    def __init__(self, taxonomy, id=None, name=None, idmap=None):
        """Constructor for the class."""
        super(Community, self).__init__(id, name)

        logger.debug("building new mico model.")

        if not (isinstance(taxonomy, pd.DataFrame) and
                all(col in taxonomy.columns for col in _taxonomy_cols)):
            raise ValueError("`taxonomy` must be a pandas DataFrame with at"
                             "least columns id and file :(")

        if "abundance" in taxonomy.columns:
            taxonomy.abundance /= taxonomy.abundance.sum()

        self.taxonomy = taxonomy

        obj = S.Zero
        self.objectives = {}
        for idx, row in taxonomy.iterrows():
            local_obj = {}
            model = load_model(row.file)
            suffix = "__" + row.id.replace(" ", "_").strip()
            for r in model.reactions:
                r.id += suffix
                obj_coef = r.objective_coefficient
                if obj_coef != 0:
                    local_obj[r.id] = obj_coef
            for m in model.metabolites:
                m.id += suffix
                m.compartment += suffix
            self.add_reactions(model.reactions)
            expr = model.objective
            expr = optlang.Objective._substitute_variables(
                   expr, model=self.solver)
            obj += expr
            self.objectives[row.id] = expr
            self.__add_exchanges(model)

        self.objective = self.solver.interface.Objective(obj, direction="max")

    def __add_exchanges(self, model):
        """Add exchange reactions for a new model."""
        pass

    def optimize_single(self, id, fluxes=True):
        """Optimize growth rate for a single model in the community."""
        if isinstance(id, six.string_types):
            if id not in self.taxonomy.id:
                raise ValueError(id + " not in taxonomy!")
            info = self.taxonomy[self.taxonomy.id == id]
        elif isinstance(id, int) and id >= 0:
            info = self.taxonomy.iloc[id]
        else:
            raise ValueError("`id` must be an id or positive index!")

        obj = self.objectives[info.id]
        with self as m:
            m.objective = obj
            m.solver.optimize()
            if fluxes:
                res = fluxes_from_primals(m, info)
            else:
                res = m.objective.value

        return res
