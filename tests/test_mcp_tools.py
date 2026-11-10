"""The tools of the recipes, generated from the registry (phase 6-ter, PR 1).

The principle of the phase (``ROADMAP.md``, phase 6-ter): a recipe added to
the registry reaches the agents with no code of the server. The tests check
it with a recipe that exists only in the test.
"""
import asyncio
import logging
import sys
import threading
import time
from types import MappingProxyType

import numpy as np
import pytest

pytest.importorskip("mcp")

from mcp import Client, StdioServerParameters  # noqa: E402

from astrocook.core.session import SessionV2, load_session_from_file  # noqa: E402
from astrocook.mcp_server.schemas import parameters_schema  # noqa: E402
from astrocook.mcp_server.server import build_server, instructions  # noqa: E402
from astrocook.mcp_server.tools import (  # noqa: E402
    EXCLUDED, SERVER_PARAMS, as_tools, tool_name)
from astrocook.mocks import MockComponent, MockTruth, generate  # noqa: E402
from astrocook.recipes import registry  # noqa: E402
from astrocook.recipes.flux import RecipeFluxV2  # noqa: E402
from astrocook.recipes.params import typed_recipe  # noqa: E402
from astrocook.recipes.progress import check_cancelled  # noqa: E402
from astrocook.recipes.registry import all_recipes, build_spec, get_recipe  # noqa: E402


logger = logging.getLogger("astrocook.test_mcp_tools")


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

    def progress_test(self, steps: int = 4):
        for i in range(steps + 1):
            logger.info("Step %d of %d", i, steps)
            logger.info(">> PROGRESS: %d", 100 * i // steps)
            if i == 1:
                logger.info(">> PROGRESS: %d", 100 * i // steps)  # repeated
        return float(steps)

    for fn in (scale_test, measure_test, progress_test):
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
    flux["progress_test"] = build_spec("flux", "progress_test", {
        "brief": "Log its progress (test).", "details": "Writes progress lines.",
        "params": [{"name": "steps", "type": int, "default": 4, "gt": 0,
                    "doc": "Number of steps."}]})
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


def test_export_ascii_as_a_tool(mock_archive, tmp_path):
    # F106: a tool since the recipe declares its prefix and does not
    # overwrite unless asked
    args = {"session_id": "mock", "targets": "F, dF", "path": str(tmp_path),
            "prefix": "out"}
    _, first, again, forced = calls(
        build_server(),
        ("session_open", {"path": str(mock_archive)}),
        ("file_export_ascii", args),
        ("file_export_ascii", args),
        ("file_export_ascii", {**args, "overwrite": True}))
    assert first.structured_content["value"] is True
    assert "out_spec.csv already exists" in error_text(again)
    assert forced.structured_content["value"] is True
    assert (tmp_path / "out_spec.csv").is_file()


# --- Progress notifications (phase 6-ter, PR 4) ----------------------------

def calls_with_progress(server, *steps):
    """Run the steps; return their results and, for each, the notifications
    of progress received as (progress, total, message)."""
    async def go():
        async with Client(server) as client:
            out = []
            for tool, args in steps:
                seen = []

                async def on_progress(progress, total, message, seen=seen):
                    seen.append((progress, total, message))
                out.append((await client.call_tool(tool, args, progress_callback=on_progress),
                            seen))
            return out
    return asyncio.run(go())


def test_the_context_is_not_a_parameter_of_the_tools(new_recipe):
    for name, tool in list_tools(build_server()).items():
        assert "ctx" not in tool.input_schema["properties"], name


def test_progress_lines_are_sent_as_notifications(new_recipe, mock_archive):
    _, (result, seen) = calls_with_progress(
        build_server(),
        ("session_open", {"path": str(mock_archive)}),
        ("flux_progress_test", {"session_id": "mock", "steps": 4}))
    assert not result.is_error
    values = [p for p, _, _ in seen]
    # The start, the lines (the repeated one and the 0 of the first are not
    # sent: the value must increase), the end
    assert values == [0, 25, 50, 75, 99, 100]
    assert all(t == 100 for _, t, _ in seen)
    # The message is the last status line, as in the GUI; the lines of the
    # recipe stop at 99, the end of the tool is 100
    assert [m for _, _, m in seen[:-1]] == [
        "Running flux_progress_test", "Step 1 of 4", "Step 2 of 4", "Step 3 of 4",
        "Step 4 of 4"]
    assert seen[-1][2].startswith("flux_progress_test done in ")


def test_a_recipe_without_progress_lines_sends_a_start_and_an_end(new_recipe, mock_archive):
    _, (_, seen) = calls_with_progress(
        build_server(),
        ("session_open", {"path": str(mock_archive)}),
        ("flux_scale_test", {"session_id": "mock"}))
    assert [(p, t) for p, t, _ in seen] == [(0, 100), (100, 100)]
    assert seen[0][2] == "Running flux_scale_test"


def test_a_failed_recipe_leaves_the_logger_as_it_was(new_recipe, mock_archive):
    package = logging.getLogger("astrocook")
    level, handlers = package.level, list(package.handlers)
    _, (bad, seen) = calls_with_progress(
        build_server(),
        ("session_open", {"path": str(mock_archive)}),
        ("flux_scale_test", {"session_id": "mock", "factor": -1}))
    assert bad.is_error
    assert package.level == level and package.handlers == handlers


def test_no_progress_is_sent_when_the_client_does_not_ask(new_recipe, mock_archive):
    _, ok = calls(build_server(),
                  ("session_open", {"path": str(mock_archive)}),
                  ("flux_progress_test", {"session_id": "mock"}))
    assert not ok.is_error and ok.structured_content["value"] == 4.0


def test_progress_over_stdio(new_recipe, mock_archive):
    # The notifications go through the protocol, with the lines of the log
    # written to stderr: stdout stays clean (a real recipe, no fake one: the
    # subprocess has not the recipe of the test)
    async def go():
        seen = []

        async def on_progress(progress, total, message):
            seen.append((progress, message))
        params = StdioServerParameters(
            command=sys.executable, args=["-m", "astrocook.mcp_server.server"])
        async with Client(params) as client:
            await client.call_tool("session_open", {"path": str(mock_archive)})
            result = await client.call_tool(
                "flux_smooth", {"session_id": "mock", "sigma_kms": 20},
                progress_callback=on_progress)
        return result, seen
    result, seen = asyncio.run(go())
    assert not result.is_error
    assert [p for p, _ in seen] == [0, 100]
    assert seen[0][1] == "Running flux_smooth"


# --- Read-only tools while a recipe runs, and the heartbeat (F113) ---------

RELEASE = threading.Event()


@pytest.fixture
def wait_recipe(new_recipe, monkeypatch):
    """Flux > Wait (test): waits for RELEASE (``seconds`` = 0) or sleeps
    ``seconds`` without logging, then logs ``then`` as its progress (if >= 0);
    with ``seconds`` < 0 it waits for RELEASE checking the cancellation."""
    def wait_test(self, seconds: float = 0.0, then: int = -1):
        if seconds > 0:
            time.sleep(seconds)
        elif seconds < 0:
            while not RELEASE.wait(0.01):
                check_cancelled()
        else:
            assert RELEASE.wait(30), "never released"
        if then >= 0:
            logger.info("Halfway")
            logger.info(">> PROGRESS: %d", then)
            time.sleep(seconds)
        return self._session.edit.apply_expression("y_waited", "y * 1")

    monkeypatch.setattr(RecipeFluxV2, "wait_test", typed_recipe(wait_test), raising=False)
    flux = dict(registry.RECIPES["flux"])
    flux["wait_test"] = build_spec("flux", "wait_test", {
        "brief": "Wait (test).", "details": "Waits, silent.",
        "params": [{"name": "seconds", "type": float, "default": 0.0,
                    "doc": "Seconds."},
                   {"name": "then", "type": int, "default": -1, "doc": "Progress after."}]})
    recipes = dict(registry.RECIPES)
    recipes["flux"] = MappingProxyType(flux)
    monkeypatch.setattr(registry, "RECIPES", MappingProxyType(recipes))
    RELEASE.clear()
    yield
    RELEASE.set()


def test_the_read_only_tools_answer_while_a_recipe_runs(wait_recipe, mock_archive, tmp_path):
    saved = tmp_path / "during.acs2"

    async def go():
        async with Client(build_server()) as client:
            await client.call_tool("session_open", {"path": str(mock_archive)})
            waiting = asyncio.create_task(
                client.call_tool("flux_wait_test", {"session_id": "mock"}))
            # Another recipe waits for the first one
            queued = asyncio.create_task(
                client.call_tool("flux_scale_test", {"session_id": "mock"}))
            for _ in range(200):
                listed = (await client.call_tool("session_list", {})).structured_content
                if listed["sessions"][0]["running"]:
                    break
                await asyncio.sleep(0.02)
            reads = {tool: await asyncio.wait_for(client.call_tool(tool, args), 10)
                     for tool, args in [
                         ("session_history", {"session_id": "mock"}),
                         ("session_describe", {"session_id": "mock"}),
                         ("components_list", {"session_id": "mock"}),
                         ("spectrum_plot", {"session_id": "mock"}),
                         ("session_note", {"session_id": "mock", "text": "Waiting."}),
                         ("session_save", {"session_id": "mock", "path": str(saved)})]}
            busy = not waiting.done() and not queued.done()
            RELEASE.set()
            done = [await asyncio.wait_for(t, 30) for t in (waiting, queued)]
            history = await client.call_tool("session_history", {"session_id": "mock"})
            return listed, reads, busy, done, history

    listed, reads, busy, done, history = asyncio.run(go())
    assert listed["sessions"][0]["running"].startswith("flux_wait_test, for ")
    assert busy
    for tool, result in reads.items():
        assert not result.is_error, tool
    # What is read and saved is the state before the recipe
    assert reads["session_history"].structured_content["states"][-1]["index"] == 0
    save = reads["session_save"].structured_content
    assert save["log_entries"] == 0
    assert any("A recipe is running on 'mock'" in w for w in save["warnings"])
    assert saved.is_file()
    assert all(not r.is_error for r in done)
    states = history.structured_content["states"]
    assert [s["recipe"] for s in states] == [None, "wait_test", "scale_test"]
    assert states[0]["notes"][0]["text"] == "Waiting."


def test_a_silent_recipe_sends_a_heartbeat(wait_recipe, mock_archive, monkeypatch):
    from astrocook.mcp_server import progress
    monkeypatch.setattr(progress, "HEARTBEAT_S", 0.05)
    _, (result, seen) = calls_with_progress(
        build_server(),
        ("session_open", {"path": str(mock_archive)}),
        ("flux_wait_test", {"session_id": "mock", "seconds": 0.4, "then": 50}))
    assert not result.is_error
    values = [p for p, _, _ in seen]
    # Always increasing: the heartbeats below 1 before the line of 50, below
    # 51 after it; the end is 100
    assert values == sorted(set(values))
    beats = [(p, m) for p, _, m in seen if m.startswith("Still running flux_wait_test: ")]
    assert any(0 < p < 1 for p, _ in beats)
    assert any(50 < p < 51 and m.endswith("; last: Halfway") for p, m in beats)
    assert 50 in values and values[-1] == 100


def test_no_heartbeat_when_the_client_does_not_ask(wait_recipe, mock_archive, monkeypatch):
    from astrocook.mcp_server import progress
    monkeypatch.setattr(progress, "HEARTBEAT_S", 0.05)
    _, result = calls(build_server(), ("session_open", {"path": str(mock_archive)}),
                      ("flux_wait_test", {"session_id": "mock", "seconds": 0.2}))
    assert not result.is_error


# --- Cancellation (phase 6-ter, PR 7-quater; F113) --------------------------

async def _wait_until_running(client):
    for _ in range(500):
        listed = (await client.call_tool("session_list", {})).structured_content
        if listed["sessions"][0]["running"]:
            return
        await asyncio.sleep(0.01)
    raise AssertionError("the recipe never started")


def test_a_client_that_cancels_its_call_stops_the_recipe(wait_recipe, mock_archive):
    # It went on (in trial 2 for hours, with the session out of reach)
    async def go():
        async with Client(build_server()) as client:
            await client.call_tool("session_open", {"path": str(mock_archive)})
            waiting = asyncio.create_task(
                client.call_tool("flux_wait_test", {"session_id": "mock", "seconds": -1}))
            await _wait_until_running(client)
            waiting.cancel()
            with pytest.raises(asyncio.CancelledError):
                await waiting
            for _ in range(500):
                listed = (await client.call_tool("session_list", {})).structured_content
                if not listed["sessions"][0]["running"]:
                    break
                await asyncio.sleep(0.01)
            history = await client.call_tool("session_history", {"session_id": "mock"})
            return listed, history

    listed, history = asyncio.run(go())
    assert listed["sessions"][0]["running"] is None, "still running"
    assert not RELEASE.is_set()
    assert [s["recipe"] for s in history.structured_content["states"]] == [None]


def test_session_cancel_stops_the_recipe(wait_recipe, mock_archive):
    async def go():
        async with Client(build_server()) as client:
            await client.call_tool("session_open", {"path": str(mock_archive)})
            idle = await client.call_tool("session_cancel", {"session_id": "mock"})
            waiting = asyncio.create_task(
                client.call_tool("flux_wait_test", {"session_id": "mock", "seconds": -1}))
            await _wait_until_running(client)
            cancelled = await client.call_tool("session_cancel", {"session_id": "mock"})
            result = await asyncio.wait_for(waiting, 10)
            history = await client.call_tool("session_history", {"session_id": "mock"})
            again = await client.call_tool("flux_scale_test", {"session_id": "mock"})
            return idle, cancelled, result, history, again

    idle, cancelled, result, history, again = asyncio.run(go())
    assert idle.structured_content["cancelled"] is None
    assert idle.structured_content["message"] == "No recipe is running on 'mock'."
    assert cancelled.structured_content["cancelled"].startswith("flux_wait_test, for ")
    assert result.is_error
    assert "cancelled while it ran: the session is as before it" in error_text(result)
    assert [s["recipe"] for s in history.structured_content["states"]] == [None]
    # The next recipe is not cancelled
    assert not again.is_error


def test_a_recipe_cancelled_after_its_last_check_is_not_kept(wait_recipe, mock_archive,
                                                             monkeypatch):
    from astrocook.mcp_server import tools
    real = tools.run_registry

    def then_cancelled(*args, **kwargs):
        result = real(*args, **kwargs)
        (event,) = store_cancels.values()
        event.set()
        return result

    store_cancels = {}
    real_init = tools.SessionStore.__init__

    def init(self):
        real_init(self)
        self.cancels = store_cancels

    monkeypatch.setattr(tools.SessionStore, "__init__", init)
    monkeypatch.setattr(tools, "run_registry", then_cancelled)
    _, result = calls(build_server(), ("session_open", {"path": str(mock_archive)}),
                      ("flux_scale_test", {"session_id": "mock"}))
    assert result.is_error
    assert "flux_scale_test was cancelled: the session is unchanged." in error_text(result)
