"""MCP server (pilot, ROADMAP.md section 4.1): tools called through the SDK
client, connected to the server in process."""
import asyncio
import logging
import threading

import pytest

pytest.importorskip("mcp")

from mcp import Client  # noqa: E402
from mcp.server.mcpserver.exceptions import ToolError  # noqa: E402

from astrocook import __version__  # noqa: E402
from astrocook.core.session import SessionV2, load_session_from_file  # noqa: E402
from astrocook.mcp_server.errors import run_recipe  # noqa: E402
from astrocook.mcp_server.server import build_server  # noqa: E402
from astrocook.mocks import MockComponent, MockTruth, generate  # noqa: E402

from conftest import J0103_SESSION  # noqa: E402

TOOLS = {"session_open", "session_save"}


def call(server, tool, **args):
    """Call one tool and return the result (structured content or error)."""
    async def go():
        async with Client(server) as client:
            return await client.call_tool(tool, args)
    return asyncio.run(go())


def calls(server, *steps):
    """Call several tools in sequence on the same server."""
    async def go():
        async with Client(server) as client:
            return [await client.call_tool(tool, args) for tool, args in steps]
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


def test_tools_have_descriptions_and_schemas():
    async def go():
        async with Client(build_server()) as client:
            return (await client.list_tools()).tools
    tools = {t.name: t for t in asyncio.run(go())}
    assert set(tools) == TOOLS
    for tool in tools.values():
        assert tool.description
        assert tool.output_schema["properties"]
        for name, prop in tool.input_schema["properties"].items():
            assert prop.get("description"), (tool.name, name)


def test_open(mock_archive):
    result = call(build_server(), "session_open", path=str(mock_archive))
    assert not result.is_error
    info = result.structured_content
    assert info["session_id"] == "mock"
    spec = info["spectrum"]
    assert spec["x_unit"] == "nm"
    assert 440.0 <= spec["x_min"] < spec["x_max"] <= 460.0
    assert spec["z_em"] == pytest.approx(2.0)
    assert spec["resol"] == pytest.approx(45000.0)
    assert spec["n_components"] == 0
    assert info["log_entries"] == 0
    assert info["saved_with_version"] == __version__


def test_open_sets_the_properties_and_logs_them(mock_archive):
    result = call(build_server(), "session_open", path=str(mock_archive),
                  z_em=2.1, resol=40000)
    info = result.structured_content
    assert info["spectrum"]["z_em"] == pytest.approx(2.1)
    assert info["spectrum"]["resol"] == pytest.approx(40000.0)
    assert info["log_entries"] == 1


def test_open_j0103():
    result = call(build_server(), "session_open", path=str(J0103_SESSION))
    spec = result.structured_content["spectrum"]
    assert spec["n_components"] == 764
    assert spec["has_continuum"]


def test_same_file_twice_gets_two_identifiers(mock_archive):
    first, second = calls(build_server(),
                          ("session_open", {"path": str(mock_archive)}),
                          ("session_open", {"path": str(mock_archive)}))
    assert first.structured_content["session_id"] == "mock"
    assert second.structured_content["session_id"] == "mock-2"


def test_open_missing_file(tmp_path):
    missing = tmp_path / "missing.fits"
    assert str(missing) in error_text(
        call(build_server(), "session_open", path=str(missing)))


def test_open_unreadable_file_reports_the_logged_error(tmp_path):
    bad = tmp_path / "bad.fits"
    bad.write_text("not a spectrum")
    text = error_text(call(build_server(), "session_open", path=str(bad)))
    assert "session_open failed" in text
    assert "no reason was logged" not in text


def test_unphysical_properties_are_accepted(mock_archive):
    # Characterization (F46): edit.set_properties does not check the
    # values, and the server passes them on as they are
    result = call(build_server(), "session_open", path=str(mock_archive),
                  z_em=-5, resol=-1)
    spec = result.structured_content["spectrum"]
    assert spec["z_em"] == -5 and spec["resol"] == -1
    assert result.structured_content["warnings"] == []


def test_unknown_session(mock_archive, tmp_path):
    _, result = calls(build_server(),
                      ("session_open", {"path": str(mock_archive)}),
                      ("session_save", {"session_id": "nope",
                                        "path": str(tmp_path / "x")}))
    text = error_text(result)
    assert "nope" in text and "mock" in text


def test_save_round_trip(mock_archive, tmp_path):
    target = tmp_path / "saved"
    _, result = calls(build_server(),
                      ("session_open", {"path": str(mock_archive), "z_em": 2.1}),
                      ("session_save", {"session_id": "mock",
                                        "path": str(target)}))
    saved = result.structured_content
    assert saved["path"] == str(tmp_path / "saved.acs2")
    assert saved["log_entries"] == 1
    # The archive reopens as in the GUI (same loader), with its log
    session = load_session_from_file(saved["path"], "saved", None, "auto")
    assert session.spec._data.z_em == pytest.approx(2.1)
    entries = session.log_manager.to_dict()["entries"]
    assert entries == [{"recipe_name": "set_properties",
                        "params": {"z_em": "2.1"}}]
    assert session.saved_with_version == __version__


def test_save_does_not_overwrite(mock_archive, tmp_path):
    target = tmp_path / "saved.acs2"
    target.write_text("keep me")
    steps = [("session_open", {"path": str(mock_archive)}),
             ("session_save", {"session_id": "mock", "path": str(target)})]
    _, refused = calls(build_server(), *steps)
    assert "overwrite" in error_text(refused)
    assert target.read_text() == "keep me"
    steps[1][1]["overwrite"] = True
    _, done = calls(build_server(), *steps)
    assert not done.is_error


def test_save_to_a_missing_folder(mock_archive, tmp_path):
    _, result = calls(build_server(),
                      ("session_open", {"path": str(mock_archive)}),
                      ("session_save", {"session_id": "mock",
                                        "path": str(tmp_path / "no" / "x")}))
    assert "No such folder" in error_text(result)


# --- run_recipe: the recipe convention (0 and a logged error) -------------

def test_run_recipe_failure_carries_the_logged_error():
    def recipe():
        logging.error("kappa must be positive")
        return 0
    with pytest.raises(ToolError, match="kappa must be positive"):
        run_recipe("continuum_estimate_auto", recipe)


def test_run_recipe_failure_without_a_logged_reason():
    with pytest.raises(ToolError, match="no reason was logged"):
        run_recipe("r", lambda: 0)


def test_run_recipe_exception_keeps_its_message():
    def recipe():
        raise ValueError("bad column")
    with pytest.raises(ToolError, match="ValueError: bad column"):
        run_recipe("r", recipe)


def test_run_recipe_returns_the_warnings_once():
    def recipe():
        logging.info("not returned")
        for _ in range(2):
            logging.warning("z_em is 0.0")
        return "result"
    assert run_recipe("r", recipe) == ("result", ["z_em is 0.0"])


def test_run_recipe_ignores_other_threads():
    started, release = threading.Event(), threading.Event()

    def other():
        started.wait()
        logging.warning("from another call")
        release.set()

    def recipe():
        started.set()
        release.wait(5)
        return "result"

    thread = threading.Thread(target=other)
    thread.start()
    result = run_recipe("r", recipe)
    thread.join()
    assert result == ("result", [])
