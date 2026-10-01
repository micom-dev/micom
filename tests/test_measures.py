"""Test some measures."""

from .fixtures import batch_grown

import micom.measures as mea
import micom as mm

data = mm.data.test_data()
db = mm.data.test_db
medium = mm.qiime_formats.load_qiime_medium(mm.data.test_medium)


def test_production(batch_grown):
    rates = mea.production_rates(batch_grown.results)
    assert "taxon" not in rates.columns
    assert all(rates.flux >= 0)
    assert all(rates.flux[rates.metabolite == "co2_e"] > 0)


def test_consumption(batch_grown):
    rates = mea.consumption_rates(batch_grown.results)
    assert "taxon" not in rates.columns
    assert all(rates.flux >= 0)
    assert all(rates.flux[rates.metabolite == "glc__D_m"] > 0)
