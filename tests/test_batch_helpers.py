"""Tests for the lower-level batch workflow helpers."""

from types import SimpleNamespace

import pandas as pd
import pytest

import micom.batch.grow as grow
import micom.batch.media as media
import micom.batch.tradeoff as tradeoff
from micom.batch.configuration import Configuration


def test_process_medium_duplicates_and_deduplicates():
    raw = pd.DataFrame({"reaction": ["EX_a"], "flux": [1.0]})
    result = media.process_medium(raw, ["sample_a", "sample_b"])
    assert result.sample_id.tolist() == ["sample_a", "sample_b"]
    assert result.index.tolist() == ["EX_a", "EX_a"]

    duplicate = pd.DataFrame(
        {
            "reaction": ["EX_a", "EX_a"],
            "flux": [1.0, 2.0],
            "sample_id": ["sample_a", "sample_a"],
        }
    )
    result = media.process_medium(duplicate, ["sample_a"])
    assert len(result) == 1


def test_process_medium_requires_each_sample():
    raw = pd.DataFrame(
        {"reaction": ["EX_a"], "flux": [1.0], "sample_id": ["sample_a"]}
    )
    with pytest.raises(ValueError, match="missing samples.*sample_b"):
        media.process_medium(raw, ["sample_a", "sample_b"])


@pytest.mark.parametrize("with_solution", [False, True])
def test_minimal_medium_helper(with_solution, monkeypatch):
    medium = pd.Series({"EX_a": 1.0})
    solution = object()
    result = {"medium": medium, "solution": solution} if with_solution else medium
    community = SimpleNamespace(solver=_fake_solver())
    monkeypatch.setattr(media, "load_pickle", lambda _: community)
    monkeypatch.setattr(media, "minimal_medium", lambda *args, **kwargs: result)
    growth_results = object()
    monkeypatch.setattr(
        media.GrowthResults,
        "from_solution",
        lambda *args, **kwargs: growth_results,
    )

    output = media._medium(
        ("sample", "model.pickle", 0.2, 0.01, False, None, with_solution)
    )
    assert output["medium"].reaction.tolist() == ["EX_a"]
    assert output["medium"].sample_id.tolist() == ["sample"]
    if with_solution:
        assert output["growth"] is growth_results


def test_minimal_medium_helper_failure(monkeypatch):
    community = SimpleNamespace(solver=_fake_solver())
    monkeypatch.setattr(media, "load_pickle", lambda _: community)
    monkeypatch.setattr(media, "minimal_medium", lambda *args, **kwargs: None)
    assert media._medium(("sample", "model", 0.2, 0.01, False, None, False)) is None


def test_fix_medium_success(monkeypatch):
    metabolite = type("Metabolite", (), {"id": "a", "name": "A"})()
    reaction = SimpleNamespace(metabolites={metabolite: -1})
    community = SimpleNamespace(
        reactions=SimpleNamespace(get_by_id=lambda _: reaction)
    )
    monkeypatch.setattr(media, "load_pickle", lambda _: community)
    monkeypatch.setattr(
        media,
        "complete_medium",
        lambda *args, **kwargs: pd.Series({"EX_a": 2.0}),
    )

    result = media._fix_medium(
        ("sample", "model", 0.2, 0.01, 10, False, pd.Series(dtype=float), None)
    )
    assert result.loc[0, "metabolite"] == "a"
    assert result.loc[0, "description"] == "A"
    assert result.loc[0, "sample_id"] == "sample"


def test_fix_medium_failure(monkeypatch):
    monkeypatch.setattr(media, "load_pickle", lambda _: object())
    monkeypatch.setattr(
        media,
        "complete_medium",
        lambda *args, **kwargs: (_ for _ in ()).throw(RuntimeError("infeasible")),
    )
    result = media._fix_medium(
        ("sample", "model", 0.2, 0.01, 10, False, pd.Series(dtype=float), None)
    )
    assert result is None


def _fake_solver(interface="cplex"):
    tolerances = SimpleNamespace(feasibility=1e-7)
    configuration = SimpleNamespace(tolerances=tolerances, presolve=None)
    return SimpleNamespace(
        interface=interface, configuration=configuration, status="optimal"
    )


def _growth_config(flux_method="none"):
    return Configuration(
        coupling={"enabled": False},
        simulation={"flux_method": flux_method},
    )


def test_growth_rejects_glpk(monkeypatch):
    community = SimpleNamespace(solver=_fake_solver("glpk"))
    monkeypatch.setattr(grow, "load_pickle", lambda _: community)
    monkeypatch.setattr(grow, "interface_to_str", lambda _: "glpk_interface")
    assert grow._growth(
        ("model", 0.5, pd.Series(dtype=float), _growth_config())
    ) is None


def test_growth_tradeoff_failure(monkeypatch):
    def fail(**kwargs):
        raise RuntimeError("infeasible")

    community = SimpleNamespace(
        id="sample",
        solver=_fake_solver(),
        exchanges=[SimpleNamespace(id="EX_a")],
        cooperative_tradeoff=fail,
    )
    monkeypatch.setattr(grow, "load_pickle", lambda _: community)
    monkeypatch.setattr(grow, "interface_to_str", lambda _: "cplex_interface")
    assert grow._growth(
        ("model", 0.5, pd.Series({"EX_a": 1.0}), _growth_config())
    ) is None


def test_growth_success(monkeypatch):
    rates = pd.DataFrame({"growth_rate": [0.2]}, index=["taxon"])
    solution = SimpleNamespace(
        members=rates,
        fluxes=pd.DataFrame({"EX_a": [0.5]}, index=["taxon"]),
    )
    exchange = SimpleNamespace(id="EX_a", global_id="EX_a")
    community = SimpleNamespace(
        id="sample",
        solver=_fake_solver(),
        exchanges=[exchange],
        internal_exchanges=[exchange],
        cooperative_tradeoff=lambda **kwargs: solution,
    )
    monkeypatch.setattr(grow, "load_pickle", lambda _: community)
    monkeypatch.setattr(grow, "interface_to_str", lambda _: "cplex_interface")
    monkeypatch.setattr(
        grow,
        "annotate_metabolites_from_exchanges",
        lambda _: pd.DataFrame({"reaction": ["EX_a"], "metabolite": ["a"]}),
    )

    result = grow._growth(
        ("model", 0.5, pd.Series({"EX_a": 1.0}), _growth_config())
    )
    assert result["growth"].sample_id.tolist() == ["sample"]
    assert result["exchanges"].loc["taxon", "sample_id"] == "sample"


def test_growth_minimal_import_failure(monkeypatch):
    rates = pd.DataFrame(
        {"growth_rate": [0.2, 0.2]}, index=["taxon", "medium"]
    )
    solution = SimpleNamespace(members=rates, growth_rate=0.2)
    community = SimpleNamespace(
        id="sample",
        solver=_fake_solver(),
        exchanges=[SimpleNamespace(id="EX_a")],
        cooperative_tradeoff=lambda **kwargs: solution,
    )
    monkeypatch.setattr(grow, "load_pickle", lambda _: community)
    monkeypatch.setattr(grow, "interface_to_str", lambda _: "cplex_interface")
    monkeypatch.setattr(grow, "minimal_medium", lambda *args, **kwargs: None)

    assert grow._growth(
        (
            "model",
            0.5,
            pd.Series({"EX_a": 1.0}),
            _growth_config("minimal imports"),
        )
    ) is None


def test_tradeoff_optimizer_failure(monkeypatch):
    def fail(**kwargs):
        raise RuntimeError("infeasible")

    community = SimpleNamespace(
        id="sample",
        solver=_fake_solver(),
        exchanges=[SimpleNamespace(id="EX_a")],
        optimize=fail,
    )
    monkeypatch.setattr(tradeoff, "load_pickle", lambda _: community)
    assert tradeoff._tradeoff(
        ("model", [0.5], pd.Series({"EX_a": 1.0}), None, None, False)
    ) is None


def test_tradeoff_cooperative_failure(monkeypatch):
    rates = pd.DataFrame({"growth_rate": [0.2]}, index=["taxon"])

    def fail(**kwargs):
        raise RuntimeError("infeasible")

    community = SimpleNamespace(
        id="sample",
        solver=_fake_solver(),
        exchanges=[SimpleNamespace(id="EX_a")],
        optimize=lambda **kwargs: SimpleNamespace(members=rates),
        cooperative_tradeoff=fail,
    )
    monkeypatch.setattr(tradeoff, "load_pickle", lambda _: community)
    assert tradeoff._tradeoff(
        ("model", [0.5], pd.Series({"EX_a": 1.0}), None, None, False)
    ) is None


def test_tradeoff_success(monkeypatch):
    initial = pd.DataFrame({"growth_rate": [0.2]}, index=["taxon"])
    cooperative = pd.DataFrame({"growth_rate": [0.3]}, index=["taxon"])
    community = SimpleNamespace(
        id="sample",
        solver=_fake_solver(),
        exchanges=[SimpleNamespace(id="EX_a")],
        optimize=lambda **kwargs: SimpleNamespace(members=initial),
        cooperative_tradeoff=lambda **kwargs: SimpleNamespace(
            solution=[SimpleNamespace(members=cooperative)], tradeoff=[0.5]
        ),
    )
    monkeypatch.setattr(tradeoff, "load_pickle", lambda _: community)
    result = tradeoff._tradeoff(
        ("model", [0.5], pd.Series({"EX_a": 1.0}), None, None, False)
    )
    assert pd.isna(result.tradeoff.iloc[0])
    assert result.tradeoff.iloc[1] == 0.5
    assert result.sample_id.tolist() == ["sample", "sample"]