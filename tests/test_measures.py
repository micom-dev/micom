"""Test some measures."""

import micom.measures as mea
import micom as mm

data = mm.data.test_data()
db = mm.data.test_db
medium = mm.qiime_formats.load_qiime_medium(mm.data.test_medium)


def test_production(tmp_path):
    built = mm.workflows.build(data, db, str(tmp_path), cutoff=0)
    grown = mm.workflows.grow(built, str(tmp_path), medium, 0.5)
    rates = mea.production_rates(grown)
    assert "taxon" not in rates.columns
    assert all(rates.flux >= 0)
    assert all(rates.flux[rates.metabolite == "co2_e"] > 0)


def test_consumption(tmp_path):
    built = mm.workflows.build(data, db, str(tmp_path), cutoff=0)
    grown = mm.workflows.grow(built, str(tmp_path), medium, 0.5)
    rates = mea.consumption_rates(grown)
    assert "taxon" not in rates.columns
    assert all(rates.flux >= 0)
    assert all(rates.flux[rates.metabolite == "glc__D_m"] > 0)
