"""Test utilities."""

import cobra
from cobra.io import read_sbml_model
import numpy as np
import micom
import micom.util as util
from .fixtures import community

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
        v.name: coef
        for v, coef in single.objective.get_linear_coefficients(
            single.variables
        ).items()
    }
    mod = util.join_models(tax.file, id="test_model")
    coefs = {
        v.name: coef
        for v, coef in mod.objective.get_linear_coefficients(mod.variables).items()
    }
    assert len(mod.reactions) == len(single.reactions) + 1  # added biomass
    assert len(mod.metabolites) == len(single.metabolites)
    assert np.allclose(single.slim_optimize(), mod.slim_optimize())


def test_compartment_id():
    met = cobra.core.Metabolite(id="test_met_e__taxon")
    met.community_id = "taxon"
    met.global_id = "text_met_e"
    met.compartment = "e__taxon"
    assert util.compartment_id(met) == "e"
    met.compartment = "C_e__taxon"
    assert util.compartment_id(met) == "e"
    met.global_id = "test_met[e]"
    assert util.compartment_id(met) == "e"
    met.global_id = "test_met(e)"
    assert util.compartment_id(met) == "e"


def test_fix_demands(tmp_path):
    fpath = str(tmp_path / "test.xml")
    model = read_sbml_model(micom.data.test_taxonomy().file[0])
    model.exchanges[0].lower_bound = 0.1
    cobra.io.write_sbml_model(model, fpath)
    model = util.load_model(fpath)
    assert model.exchanges[0].lower_bound == 0.0
