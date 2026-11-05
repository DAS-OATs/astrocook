"""The tools of the recipes, generated from the registry (phase 6-ter, PR 1).

The principle of the phase (``ROADMAP.md``, phase 6-ter): a recipe added to
the registry reaches the agents with no code of the server. The tests check
it with a recipe that exists only in the test.
"""
import asyncio
from types import MappingProxyType

import numpy as np
import pytest

pytest.importorskip("mcp")

from mcp import Client  # noqa: E402

from astrocook.core.session import SessionV2, load_session_from_file  # noqa: E402
from astrocook.mcp_server.schemas import parameters_schema  # noqa: E402
from astrocook.mcp_server.server import build_server, instructions  # noqa: E402
from astrocook.mcp_server.tools import (  # noqa: E402
    EXCLUDED, SERVER_PARAMS, as_tools, tool_name)
from astrocook.mocks import MockComponent, MockTruth, generate  # noqa: E402
from astrocook.recipes import registry  # noqa: E402
from astrocook.recipes.flux import RecipeFluxV2  # noqa: E402
from astrocook.recipes.params import typed_recipe  # noqa: E402
from astrocook.recipes.registry import all_recipes, build_spec, get_recipe  # noqa: E402


def calls(server, *steps):
    async def go():
        async with Client(server) as client:
            return [await client.call_tool(tool, args) for tool, args in steps]
    return asyncio.run(go())


def list_tools(server):
    async def go():
        async with Client(server) as client:
            return {t.name: t for t in (await client.list_tools()).tools}
    return asyncio.run(go())


def error_text(result):
    assert result.is_error
    return result.content[0].text


@pytest.fixture(scope="module")
def mock_archive(tmp_path_factory):
    """A small mock (one CIV doublet, z_em = 2.0, R = 45000) as .acs2."""
    truth = MockTruth(components=[MockComponent("CIV", 1.9, 13.8, 10.0)],
                      lam_min_nm=440.0, lam_max_nm=460.0, dv_kms=2.5,
                      resol=45000.0, snr=50.0, z_em=2.0, seed=1)
    path = tmp_path_factory.mktemp("mock") / "mock.acs2"
    SessionV2("mock", None, spec=generate(truth)).save(str(path))
    return path


# --- Every recipe is a tool, or excluded with its reason -----------------

def test_every_recipe_is_a_tool_or_excluded_with_a_reason():
    tools = list_tools(build_server())
    recipes = {f"{s.category}.{s.name}": s for s in all_recipes()}
    for key, reason in EXCLUDED.items():
        assert key in recipes, f"{key} excluded but not in the registry"
        assert reason.strip()
        assert tool_name(recipes[key]) not in tools
    for key, spec in recipes.items():
        assert key in EXCLUDED or tool_name(spec) in tools, (
            f"The recipe {key} is neither a tool nor excluded with a reason")


def test_the_schemas_of_the_tools_are_those_of_the_registry():
    tools = list_tools(build_server())
    for spec in all_recipes():
        name = tool_name(spec)
        if name not in tools:
            continue
        schema = tools[name].input_schema
        expected = parameters_schema(spec)
        for param, prop in expected["properties"].items():
            assert schema["properties"][param] == prop, (name, param)
        assert set(schema["properties"]) == set(expected["properties"]) | set(SERVER_PARAMS)
        assert schema["required"] == ["session_id", *expected["required"]]
        assert schema["additionalProperties"] is False
        # The parameters not offered to agents are not in the schema
        for p in spec.params:
            assert (p.name in schema["properties"]) == p.agents, (name, p.name)
        assert tools[name].description.startswith(spec.brief)
        assert tools[name].output_schema["properties"]["recipe"]


def test_the_instructions_list_the_menus_of_the_registry():
    text = instructions()
    for spec in all_recipes():
        assert (tool_name(spec) in text) == (f"{spec.category}.{spec.name}" not in EXCLUDED)


def test_the_messages_name_the_tools():
    assert as_tools("set it with Set session properties (edit.set_properties).") == (
        "set it with Set session properties (tool edit_set_properties).")


# --- A recipe added to the registry, with no code of the server ----------

@pytest.fixture
def new_recipe(monkeypatch):
    """A recipe that exists only in this test: Flux > Scale (test), which
    writes the flux times ``factor`` in a new column, and Flux > Measure
    (test), which returns a value."""
    def scale_test(self, factor: float = 2.0, column: str = "y_scaled"):
        return self._session.edit.apply_expression(column, f"y * {factor}")

    def measure_test(self):
        return float(np.nanmedian(self._session.spec.y.value))

    for fn in (scale_test, measure_test):
        monkeypatch.setattr(RecipeFluxV2, fn.__name__, typed_recipe(fn), raising=False)
    flux = dict(registry.RECIPES["flux"])
    flux["scale_test"] = build_spec("flux", "scale_test", {
        "brief": "Scale the flux (test).",
        "details": "Writes the flux times a factor in a new column.",
        "params": [
            {"name": "factor", "type": float, "default": 2.0, "gt": 0,
             "doc": "Factor of the flux."},
            {"name": "column", "type": str, "default": "y_scaled", "kind": "column",
             "doc": "Column written."}]})
    flux["measure_test"] = build_spec("flux", "measure_test", {
        "brief": "Measure the flux (test).", "details": "The median flux."})
    recipes = dict(registry.RECIPES)
    recipes["flux"] = MappingProxyType(flux)
    monkeypatch.setattr(registry, "RECIPES", MappingProxyType(recipes))


def test_a_new_recipe_is_a_tool(new_recipe):
    tools = list_tools(build_server())
    tool = tools["flux_scale_test"]
    assert tool.description.startswith("Scale the flux (test). Writes the flux")
    assert tool.input_schema["properties"]["factor"] == {
        "type": "number", "description": "Factor of the flux.",
        "exclusiveMinimum": 0, "default": 2.0}
    assert "flux_measure_test" in tools
    assert "flux_scale_test" in instructions()


def test_a_new_recipe_runs_is_checked_and_logged(new_recipe, mock_archive, tmp_path):
    _, bad, ran, saved = calls(
        build_server(),
        ("session_open", {"path": str(mock_archive)}),
        ("flux_scale_test", {"session_id": "mock", "factor": -1}),
        ("flux_scale_test", {"session_id": "mock", "factor": 3}),
        ("session_save", {"session_id": "mock", "path": str(tmp_path / "s.acs2")}))
    # The value is checked by the registry, with the message of the recipe
    assert "factor" in error_text(bad)
    info = ran.structured_content
    assert info["recipe"] == "flux.scale_test"
    assert info["changes"]["columns_added"] == ["y_scaled"]
    assert info["summary"] == "Scale the flux (test): column 'y_scaled' added."
    assert info["log_entries"] == 1
    reloaded = load_session_from_file(str(tmp_path / "s.acs2"), "s", None, "auto")
    (entry,) = reloaded.log_manager.to_dict()["entries"]
    assert entry == {"recipe_name": "scale_test",
                     "params": {"factor": 3.0, "column": "y_scaled"}}
    y = reloaded.spec.y.value
    assert reloaded.spec.get_column("y_scaled").value == pytest.approx(3 * y)


def test_a_recipe_that_returns_a_value(new_recipe, mock_archive):
    _, result = calls(build_server(),
                      ("session_open", {"path": str(mock_archive)}),
                      ("flux_measure_test", {"session_id": "mock"}))
    info = result.structured_content
    assert info["value"] == pytest.approx(np.nanmedian(
        load_session_from_file(str(mock_archive), "m", None, "auto").spec.y.value))
    # The session does not change, nothing is logged
    assert info["changes"] is None and info["log_entries"] == 0


# --- The recipes of the registry -------------------------------------------

def test_unknown_parameters_are_refused(mock_archive):
    _, result = calls(build_server(),
                      ("session_open", {"path": str(mock_archive)}),
                      ("flux_smooth", {"session_id": "mock", "fwhm": 3}))
    message = error_text(result)
    assert "Unknown parameters for flux_smooth: fwhm" in message
    params = [p.name for p in get_recipe("flux", "smooth").params if p.agents]
    assert all(p in message for p in params)


def test_a_value_out_of_bounds_is_refused_with_the_message_of_the_recipe(mock_archive):
    _, result = calls(build_server(),
                      ("session_open", {"path": str(mock_archive)}),
                      ("continuum_estimate_auto", {"session_id": "mock", "kappa": -1}))
    assert "kappa" in error_text(result)


def test_a_recipe_reads_another_open_session_by_its_identifier(mock_archive):
    # A parameter of kind 'session' names an open session: its identifier
    _, _, result, missing = calls(
        build_server(),
        ("session_open", {"path": str(mock_archive)}),
        ("session_open", {"path": str(mock_archive)}),
        ("flux_equalize", {"session_id": "mock", "reference_session": "mock-2"}),
        ("flux_equalize", {"session_id": "mock", "reference_session": "nope"}))
    assert not result.is_error, error_text(result)
    assert result.structured_content["log_entries"] == 1
    assert "nope" in error_text(missing)


def test_a_branching_recipe_opens_a_session_with_the_log_of_its_parent(
        mock_archive, tmp_path):
    _, _, branched, saved = calls(
        build_server(),
        ("session_open", {"path": str(mock_archive), "z_em": 2.0}),
        ("flux_smooth", {"session_id": "mock"}),
        ("edit_split", {"session_id": "mock", "expression": "x > 450"}),
        ("session_save", {"session_id": "mock-1", "path": str(tmp_path / "b.acs2")}))
    info = branched.structured_content
    assert not branched.is_error, error_text(branched)
    # The parent does not change; the new session has its own identifier
    assert info["log_entries"] == 2 and info["changes"] is None
    (new,) = info["new_sessions"]
    assert new["session_id"] == "mock-1"
    assert new["spectrum"]["n_pixels"] < info["spectrum"]["n_pixels"]
    reloaded = load_session_from_file(str(tmp_path / "b.acs2"), "b", None, "auto")
    names = [e["recipe_name"] for e in reloaded.log_manager.to_dict()["entries"]]
    assert names == ["set_properties", "smooth", "split"]


def test_the_components_removed_are_listed_with_their_reason(mock_archive):
    _, first, second = calls(
        build_server(),
        ("session_open", {"path": str(mock_archive), "z_em": 2.0, "resol": 45000}),
        ("absorbers_doublets_auto", {"session_id": "mock", "provide_missing": True}),
        ("absorbers_doublets_auto", {"session_id": "mock"}))
    added = first.structured_content["fit"]["components"]
    changes = second.structured_content["changes"]
    assert changes["n_components_removed"] == len(added) > 0
    removed = changes["components_removed"]
    assert {c["uuid"] for c in removed} == {c["uuid"] for c in added}
    assert all(c["reason"] for c in removed)
