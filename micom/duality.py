"""Implements a fast dual formulation."""

from sympy.core.singleton import S
from micom.logger import logger


def fast_dual(model, prefix="dual_"):
    """Add dual formulation to the problem.

    A mathematical optimization problem can be viewed as a primal and a dual
    problem. If the primal problem is a minimization problem the dual is a
    maximization problem, and the optimal value of the dual is a lower bound of
    the optimal value of the primal. For linear problems, strong duality holds,
    which means that the optimal values of the primal and dual are equal
    (duality gap = 0). This functions takes an optlang Model representing a
    primal linear problem
    and returns a new Model representing the dual optimization problem. The
    provided model must have a linear objective, linear constraints and only
    continuous variables. Furthermore, the problem must be in standard form,
    i.e. all variables should be non-negative. Both minimization and
    maximization problems are allowed.

    Attributes
    ----------
    model : cobra.Model
        The model to be dualized.
    prefix : str
        The string that will be prepended to all variable and constraint names
        in the returned dual problem.

    Returns
    -------
    dict
        The coefficients for the new dual objective.

    """
    logger.info("adding dual variables")
    if len(model.variables) > 1e5:
        logger.warning(
            "the model has a lot of variables,"
            "dual optimization will be extremely slow :O"
        )
    prob = model.problem

    maximization = model.objective.direction == "max"

    if maximization:
        sign = 1
    else:
        sign = -1

    coefficients = {}
    dual_objective = {}
    to_add = []

    # Add dual variables from primal constraints:
    for constraint in model.constraints:
        if constraint.expression == 0:
            continue  # Skip empty constraint
        if not constraint.is_Linear:
            raise ValueError(
                "Non-linear problems are not supported: " + str(constraint)
            )
        if constraint.lb is None and constraint.ub is None:
            logger.warning("skipped free constraint %s" % constraint.name)
            continue  # Skip free constraint
        if constraint.lb == constraint.ub:
            const_var = prob.Variable(
                prefix + constraint.name + "_constraint", lb=None, ub=None
            )
            to_add.append(const_var)
            if constraint.lb != 0:
                dual_objective[const_var.name] = sign * constraint.lb
            coefs = constraint.get_linear_coefficients(constraint.variables)
            for variable, coef in coefs.items():
                coefficients.setdefault(variable.name, {})[const_var.name] = (
                    sign * coef
                )
        else:
            if constraint.lb is not None:
                lb_var = prob.Variable(
                    prefix + constraint.name + "_constraint_lb", lb=0, ub=None
                )
                to_add.append(lb_var)
                if constraint.lb != 0:
                    dual_objective[lb_var.name] = -sign * constraint.lb
            if constraint.ub is not None:
                ub_var = prob.Variable(
                    prefix + constraint.name + "_constraint_ub", lb=0, ub=None
                )
                to_add.append(ub_var)
                if constraint.ub != 0:
                    dual_objective[ub_var.name] = sign * constraint.ub

            if not (
                constraint.expression.is_Add or constraint.expression.is_Mul
            ):
                raise ValueError(
                    "Invalid expression type: "
                    + str(type(constraint.expression))
                )
            if constraint.expression.is_Add:
                coefficients_dict = constraint.get_linear_coefficients(
                    constraint.variables
                )
            else:  # constraint.expression.is_Mul:
                args = constraint.expression.args
                coefficients_dict = {args[1]: args[0]}

            for variable, coef in coefficients_dict.items():
                if constraint.lb is not None:
                    coefficients.setdefault(variable.name, {})[lb_var.name] = (
                        -sign * coef
                    )
                if constraint.ub is not None:
                    coefficients.setdefault(variable.name, {})[ub_var.name] = (
                        sign * coef
                    )

    # Add dual variables from primal bounds
    for variable in model.variables:
        if not variable.type == "continuous":
            raise ValueError(
                "Integer variables are not supported: " + str(variable)
            )
        if variable.lb is not None and variable.lb < 0:
            raise ValueError(
                "Problem is not in standard form ("
                + variable.name
                + " can be negative)"
            )
        if variable.lb > 0:
            bound_var = prob.Variable(
                prefix + variable.name + "_lb", lb=0, ub=None
            )
            to_add.append(bound_var)
            coefficients.setdefault(variable.name, {})[bound_var.name] = -sign
            dual_objective[bound_var.name] = -sign * variable.lb
        if variable.ub is not None:
            bound_var = prob.Variable(
                prefix + variable.name + "_ub", lb=0, ub=None
            )
            to_add.append(bound_var)
            coefficients.setdefault(variable.name, {})[bound_var.name] = sign
            if variable.ub != 0:
                dual_objective[bound_var.name] = sign * variable.ub

    model.add_cons_vars(to_add)

    # Add dual constraints from primal objective
    primal_objective_dict = model.objective.get_linear_coefficients(
        model.objective.variables
    )
    for variable in model.objective.variables:
        obj_coef = primal_objective_dict[variable]
        if maximization:
            const = prob.Constraint(
                S.Zero, lb=obj_coef, name=prefix + variable.name
            )
        else:
            const = prob.Constraint(
                S.Zero, ub=obj_coef, name=prefix + variable.name
            )
        model.add_cons_vars([const])
        model.solver.update()
        coefs = {
            model.variables[vid]: coef
            for vid, coef in coefficients[variable.name].items()
        }
        const.set_linear_coefficients(coefs)

    # Make dual objective
    coefs = {
        model.variables[vid]: coef
        for vid, coef in dual_objective.items()
        if coef != 0
    }
    logger.info("dual model has {} terms in objective".format(len(coefs)))

    return coefs
