"""Test utilities."""

from os.path import basename
import numpy as np
import micom
import micom.util as util
from fixtures import community

URL = "http://bigg.ucsd.edu/static/models/e_coli_core.xml.gz"
tax = micom.data.test_taxonomy()


def test_download(tmpdir):
    print(tmpdir.dirpath())
    util.download_model(URL, str(tmpdir))
    assert tmpdir.join("e_coli_core.xml.gz").check()

    model = util.load_model(URL)
    assert len(model.reactions) == 95
    assert len(model.metabolites) == 72


def test_load_model():
    row = tax.loc[0]
    model = util.load_model(row.file)
    assert len(model.reactions) == 95
    assert len(model.metabolites) == 72


def test_serialization(tmpdir):
    row = tax.loc[0]
    util.serialize_models([row.file], str(tmpdir))
    assert tmpdir.join("e_coli_core.pickle").check()


def test_fluxes_from_primals(community):
    community.solver.optimize()
    fluxes = util.fluxes_from_primals(community, tax.loc[0])
    assert len(fluxes) < len(community.reactions)
    assert len(fluxes) == 95


def test_join_models():
    single = util.load_model(tax.file[0])
    single_coefs = {
        v.name: coef for v, coef in
        single.objective.get_linear_coefficients(single.variables).items()
    }
    mod = util.join_models(tax.file, id="test_model")
    coefs = {
        v.name: coef for v, coef in
        mod.objective.get_linear_coefficients(mod.variables).items()
    }
    assert len(mod.reactions) == len(single.reactions) + 1  # added biomass
    assert len(mod.metabolites) == len(single.metabolites)
    assert np.allclose(single.slim_optimize(), mod.slim_optimize())
