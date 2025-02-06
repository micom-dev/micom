"""Adds different coupling and resource constraints to a model."""

from itertools import chain
from optlang.symbolics import Zero
from uuid import uuid4


def add_coupling(model, reactions, objective_variable, coupling=400, lower=0.0, exclude_variables=None):
    """Adds coupling constraints to a model.

    Note that this will add a lot of constraints which can be slow for large models.

    Parameters
    ----------
    model : cobra.Model
        The model to add the constraints to.
    reactions : iterable
        The reactions to couple to the biomass reaction.
    biomass_expression : Variable
        The biomass variable to couple to.
    coupling : float, optional
        The coupling constraint to add (Default: 400).
    lower : float, optional
        The lower bound for the coupled reaction constraint (Default: 0.0).
    exclude_variables: iterable
        Variables to exclude from the constraint. For instance if they are associated
        with the biomass reaction.

    Returns
    -------
    Nothing. The constraints are added to the model in place.
    """
    for reaction in reactions:
        v = reaction.forward_variable
        if v in exclude_variables:
            continue
        const_fwd = model.problem.Constraint(
            Zero,
            lb=None,
            ub=lower,
            name="coupling_" + v.name,
        )
        model.add_cons_vars([const_fwd])
        model.constraints[const_fwd.name].set_linear_coefficients({objective_variable: -coupling, v: 1.0})
        if reaction.reversibility:
            v = reaction.reverse_variable
            if v in exclude_variables:
                continue
            const_rev = model.problem.Constraint(
                objective_variable - coupling * v,
                lb=None,
                ub=lower,
                name="coupling_" + v.name,
            )
            model.add_cons_vars([const_rev])
            model.constraints[const_rev.name].set_linear_coefficients({objective_variable: -coupling, v: 1.0})
    model.solver.update()


def add_resource_constraint(
    model,
    reactions,
    objective_variable,
    constraint=1000,
    lower=0.0,
    coupled=True,
    name=f"resource_constraint_{uuid4().hex[:8]}",
    exclude_variables=None,
):
    """Adds a resource constraint to a model.

    Parameters
    ----------
    model : cobra.Model
        The model to add the constraints to.
    reactions : iterable
        The reactions to add the resource constraint to.
    objective_variable : Variable
        The biomass variable to couple to.
    constraint : float, optional
        The resource constraint to add (Default: 35).
    coupled : bool, optional
        If True, the resource constraint is coupled to the biomass reaction (Default: True).
    name: str
        The name of the resource constraint.
    exclude_variables: iterable
        Variables to exclude from the constraint. For instance if they are associated
        with the biomass reaction.

    Returns
    -------
    Nothing. The constraints are added to the model in place.

    """
    reaction_variables = (
        (rxn.forward_variable, rxn.reverse_variable) for rxn in reactions
    )
    variables = [v for v in chain(*reaction_variables) if v not in exclude_variables]
    const = model.problem.Constraint(Zero, name=name, lb=None, ub=0)
    model.add_cons_vars([const])
    coeffs = {v: 1.0 for v in variables}
    if coupled:
        coeffs[objective_variable] = -constraint
    model.constraints[const.name].set_linear_coefficients(coeffs)
    const.ub = lower if coupled else constraint
    model.solver.update()
