"""Test the configuration interface."""

from pathlib import Path

import pytest

from micom.batch import Configuration


def test_config():
    """Test the configuration interface."""
    conf = Configuration()
    assert conf.simulation.flux_method == "minimal imports"
    assert conf.tolerance == 1e-6


def test_simulation_config():
    """Test the simulation configuration interface."""
    conf = Configuration()
    conf.simulation.flux_method = "pFBA"
    assert conf.simulation.flux_method == "pFBA"

    with pytest.raises(ValueError):
        conf.simulation.flux_method = "invalid"


def test_dbs_config():
    """Test the dbs configuration interface."""
    conf = Configuration()
    assert conf.dbs.microbial is not None
    conf.dbs.microbial = None
    assert conf.dbs.microbial is None

    with pytest.raises(ValueError):
        conf.dbs.microbial = 12

    conf.dbs.host = {
        "enterocytes": "/home/test/test.sbml",
        "goblet_cells": "/home/test/goblet_cells.xml",
    }
    assert conf.dbs.host["enterocytes"].name == "test.sbml"

    with pytest.raises(ValueError):
        conf.dbs.host = "default://test.sbml"

    conf.dbs.download_location = Path("/home/test/downloads")
    assert conf.dbs.download_location.name == "downloads"

    with pytest.raises(ValueError):
        conf.dbs.download_location = None


def test_media_config():
    """Test the media configuration interface."""
    conf = Configuration()
    assert conf.media.max_import == 10.0
    conf.media.max_import = 5.0
    assert conf.media.max_import == 5.0

    with pytest.raises(ValueError):
        conf.media.max_import = "large"


def test_coupling_config():
    """Test the coupling configuration interface."""
    conf = Configuration()
    assert conf.coupling.enabled is True
    conf.coupling.enabled = False
    assert conf.coupling.enabled is False

    with pytest.raises(ValueError):
        conf.coupling.strategy = "invalid"

    conf.coupling.strategy = "resource constraint"
    assert conf.coupling.strategy == "resource constraint"

    with pytest.raises(ValueError):
        conf.coupling.constraint = "high"


def test_build_config():
    """Test the build configuration interface."""
    conf = Configuration()
    assert conf.build.cutoff == 1e-4
    conf.build.cutoff = 1e-3
    assert conf.build.cutoff == 1e-3

    with pytest.raises(ValueError):
        conf.build.cutoff = "low"


def test_config_json(tmp_path):
    """Test the configuration json interface."""
    conf = Configuration()
    conf.simulation.flux_method = "pFBA"
    conf.dbs.microbial = None
    conf.media.max_import = 5.0
    conf.coupling.enabled = False
    conf.build.cutoff = 1e-3

    conf_file = tmp_path / "config.json"
    conf.to_json(conf_file)
    assert conf_file.exists()

    new_conf = Configuration.from_json(conf_file)
    assert new_conf.simulation.flux_method == "pFBA"
    assert new_conf.dbs.microbial is None
    assert new_conf.media.max_import == 5.0
    assert new_conf.coupling.enabled is False
    assert new_conf.build.cutoff == 1e-3


def test_config_yaml(tmp_path):
    """Test the configuration yaml interface."""
    conf = Configuration()
    conf.simulation.flux_method = "pFBA"
    conf.dbs.microbial = None
    conf.media.max_import = 5.0
    conf.coupling.enabled = False
    conf.build.cutoff = 1e-3

    conf_file = tmp_path / "config.yaml"
    conf.to_yaml(conf_file)
    assert conf_file.exists()

    new_conf = Configuration.from_yaml(conf_file)
    assert new_conf.simulation.flux_method == "pFBA"
    assert new_conf.dbs.microbial is None
    assert new_conf.media.max_import == 5.0
    assert new_conf.coupling.enabled is False
    assert new_conf.build.cutoff == 1e-3


def test_config_str():
    """Test the configuration string representation."""
    conf = Configuration()
    conf_str = str(conf)
    assert "Configuration" in conf_str
    assert "simulation" in conf_str
    assert "coupling" in conf_str
    assert "build" in conf_str
    assert "dbs" in conf_str
    assert "media" in conf_str


def test_config_repr_html():
    """Test the configuration HTML representation."""
    conf = Configuration()
    html_str = conf._repr_html_()
    assert "<li>" in html_str
    assert "simulation" in html_str
    assert "coupling" in html_str
    assert "build" in html_str
    assert "dbs" in html_str
    assert "media" in html_str
