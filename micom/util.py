"""Holds utility functions for other modules."""

import cobra.io as io
from cobra.util.context import get_context
from cobra.util.solver import interface_to_str, linear_reaction_coefficients
from cobra import Reaction
import os.path as path
from functools import partial
import pickle
from urllib.parse import urlparse
import urllib.request as urlreq
import tempfile
from shutil import rmtree
import pandas as pd
import re
from micom.logger import logger


_read_funcs = {
    ".xml": io.read_sbml_model,
    ".gz": io.read_sbml_model,
    ".mat": io.load_matlab_model,
    ".json": io.load_json_model,
    ".pickle": lambda fn: pickle.load(open(fn, "rb")),
}

COMPARTMENT_RE = "(_{}$)|([^a-zA-Z0-9 :]{}[^a-zA-Z0-9 :]$)"


def download_model(url, folder="."):
    """Download a model."""
    dest = path.join(folder, path.basename(url))
    urlreq.urlretrieve(url, dest)

    return dest


def _read_model(file):
    """Read a model from a local file."""
    _, ext = path.splitext(file)
    read_func = _read_funcs[ext]
    return read_func(file)


def load_model(filepath):
    """Load a cobra model from several file types."""
    logger.info("reading model from {}".format(filepath))
    parsed = urlparse(filepath)
    if parsed.scheme and parsed.netloc:
        tmpdir = tempfile.mkdtemp()
        logger.info("created temporary directory {}".format(tmpdir))
        filepath = download_model(filepath, folder=tmpdir)
        model = _read_model(filepath)
        rmtree(tmpdir)
        logger.info("deleted temporary directory {}".format(tmpdir))
    else:
        model = _read_model(filepath)
    return model


def load_pickle(filename):
    """Load a community model from a pickled version.

    Parameters
    ----------
    filename : str
        The file the community is stored in.

    Returns
    -------
    micom.Community
        The loaded community model.

    """
    with open(filename, mode="rb") as infile:
        mod = pickle.load(infile)
        adjust_solver_config(mod.solver)
        return mod


def serialize_models(files, dir="."):
    """Convert several models to Python pickles."""
    for f in files:
        fname = path.basename(f).split(".")[0]
        model = load_model(f)
        logger.info("serializing {}".format(f))
        pickle.dump(
            model, open(path.join(dir, fname + ".pickle"), "wb"), protocol=2
        )  # required for Python 2 compat


def chr_or_input(m):
    """Return ascii character for the ordinal or the original string."""
    i = int(m.groups()[0])
    if i > 31 and i < 128:
        return chr(i)
    else:
        return "__%d__" % i


def clean_ids(id):
    """Clean ids up a bit."""
    return re.sub("__(\\d+)__", chr_or_input, id)


def compartment_id(micom_obj):
    """Get the compartment id for a cobra object and prune the prefix if needed."""
    comp_id = micom_obj.compartment.replace("__" + micom_obj.community_id, "")
    pruned = comp_id.replace("C_", "")
    match_original = (
        re.search(COMPARTMENT_RE.format(comp_id, comp_id), micom_obj.global_id)
        is not None
    )
    match_pruned = (
        re.search(COMPARTMENT_RE.format(pruned, pruned), micom_obj.global_id)
        is not None
    )
    if not match_original and match_pruned:
        return pruned
    return comp_id


def join_models(model_files, id=None):
    """Join several models into one.

    This requires all the models to use the same ID system.

    Arguments
    ----------
    model_files : list of strings
        The files to be joined.
    id : str
        The new ID for the model. Will be the ID of the first model if None.

    Returns
    -------
    cobra.Model
        The joined cobra Model.

    """
    model = load_model(model_files[0])
    n = len(model_files)
    if id:
        model.id = id
    biomass = Reaction(
        id="micom_combined_biomass",
        name="combined biomass reaction from model joining",
        subsystem="biomass production",
        lower_bound=0,
        upper_bound=1000,
    )
    coefs = linear_reaction_coefficients(model, model.reactions)
    for r, coef in coefs.items():
        biomass += r * (coef / n)
    rids = set(r.id for r in model.reactions)
    for filepath in model_files[1:]:
        other = load_model(filepath)
        new = [r.id for r in other.reactions if r.id not in rids]
        model.add_reactions(other.reactions.get_by_any(new))
        coefs = linear_reaction_coefficients(other, other.reactions)
        for r, coef in coefs.items():
            biomass += model.reactions.get_by_id(r.id) * (coef / n)
        rids.update(new)
    model.add_reactions([biomass])
    model.objective = biomass

    return model


def fluxes_from_primals(model, info):
    """Extract a list of fluxes from the model primals."""
    primals = model.solver.primal_values
    rxns = model.reactions.query(lambda r: info.id == r.community_id)
    rids = [r.global_id for r in rxns]

    fluxes = (
        primals[rxn.forward_variable.name] - primals[rxn.reverse_variable.name]
        for rxn in rxns
    )
    fluxes = pd.Series(fluxes, rids, name=info.id)

    return fluxes


def add_var_from_expression(model, name, expr, lb=None, ub=None):
    """Add a variable to a model equaling an expression."""
    var = model.problem.Variable(name, lb=lb, ub=ub)
    const = model.problem.Constraint(
        (var - expr).expand(), lb=0, ub=0, name=name + "_equality"
    )
    model.add_cons_vars([var, const])
    return var


def check_modification(community):
    """Check whether a community already carries a modification.

    Arguments
    ---------
    community : micom.Community
        The community class to check.

    Raises
    ------
    ValueError
        If the community already carries a modification and adding another
        would not be safe.

    """
    if community.modification is not None:
        raise ValueError(
            "Community already carries a modification "
            "({})!".format(community.modification)
        )


def _format_min_growth(min_growth, taxa):
    """Format min_growth into a pandas series.

    Arguments
    ---------
    min_growth : positive float or array-like object.
        The minimum growth rate for each individual in the community. Either
        a single value applied to all individuals or one value for each.
    taxa : array-like
        The ID for each individual model in the community.

    Returns
    -------
    pandas.Series
        A pandas Series mapping each individual to its minimum growth rate.

    """
    try:
        min_growth = float(min_growth)
    except (TypeError, ValueError):
        if len(min_growth) != len(taxa):
            raise ValueError(
                "min_growth must be single value or an array-like "
                "object with an entry for each taxon in the model."
            )
    return pd.Series(min_growth, taxa)


def _apply_min_growth(community, min_growth, atol=1e-6, rtol=1e-6):
    """Set minimum growth constraints on a model.

    Will integrate with the context.
    """
    context = get_context(community)

    def reset(taxon, lb):
        logger.info("resetting growth rate constraint for %s" % taxon)
        community.constraints["objective_" + taxon].ub = None
        community.constraints["objective_" + taxon].lb = lb

    for sp in community.taxa:
        logger.info("setting growth rate constraint for %s" % sp)
        obj = community.constraints["objective_" + sp]
        if context:
            context(partial(reset, sp, obj.lb))
        if min_growth[sp] > atol:
            obj.lb = (1.0 - rtol) * min_growth[sp] - atol
        else:
            logger.info(
                "minimal growth rate smaller than tolerance," " setting to zero."
            )
            obj.lb = 0


def adjust_solver_config(solver):
    """Adjust the optlang solver configuration for larger problems."""
    interface = interface_to_str(solver.interface)

    solver.configuration.tolerances.feasibility = 1e-6
    if interface == "cplex":
        solver.configuration.lp_method = "barrier"
        solver.configuration.qp_method = "barrier"
        solver.problem.parameters.threads.set(1)
        solver.problem.parameters.barrier.convergetol.set(1e-9)
    if interface == "gurobi":
        solver.configuration.lp_method = "barrier"
        solver.problem.Params.BarConvTol = 1e-9
        solver.problem.Params.BarIterLimit = 10000
        solver.problem.Params.Threads = 1
        solver.problem.Params.LogToConsole = 0
    if interface == "glpk":
        solver.configuration.presolve = True
    if interface == "osqp":
        # as a direct solver OSQP has trouble getting to good accuracies for LPs
        solver.configuration.tolerances.optimality = 1e-3
        solver.configuration.tolerances.feasibility = 1e-4
        solver.configuration.presolve = False


def reset_min_community_growth(com):
    """Reset the lower bound for the community growth."""
    com.variables.community_objective.lb = 0.0
    com.variables.community_objective.ub = None
