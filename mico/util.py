"""Holds utility functions for other modules."""

import cobra.io as io
import os.path as path
from six.moves.urllib.parse import urlparse
import six.moves.urllib.request as urlreq
import tempfile


_read_funcs = {".xml": io.read_sbml_model,
               ".mat": io.load_matlab_model,
               ".json": io.load_json_model}


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
    with tempfile.TemporaryDirectory() as tmpdir:
        parsed = urlparse(filepath)
        if parsed.scheme and parsed.netloc:
            filepath = download_model(filepath, folder=tmpdir)
        return _read_model(filepath)


def fluxes_from_primals(model, info):
    """Extract a list of fluxes from the model primals."""
    suffix = "__" + info.name.replace(" ", "_").strip()
    primals = model.solver.primal_values
    rxns = model.reactions.query(lambda r: suffix in r.id)

    fluxes = {primals[rxn.forward_variable.id] -
              primals[rxn.reverse_variable.id] for rxn in rxns}

    return fluxes
