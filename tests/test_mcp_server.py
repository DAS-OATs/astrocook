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
from astrocook.mcp_server.outputs import identification_summary  # noqa: E402
from astrocook.mcp_server.server import SCHEMAS, build_server  # noqa: E402
from astrocook.mocks import MockComponent, MockTruth, generate  # noqa: E402

from conftest import J0103_SESSION, J0103_Z_EM  # noqa: E402
from full_chain import mock_session  # noqa: E402

TOOLS = {"session_open", "session_save", "continuum_estimate_auto",
         "absorbers_identify_lines"}


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


# --- Continuum and identification ------------------------------------------

#: CIV systems of the full-chain mock (tests/full_chain.py)
MOCK_CIV_Z = (1.76, 1.83, 1.90)
C_KMS = 299792.458


@pytest.fixture(scope="module")
def forest_archive(tmp_path_factory):
    """The full-chain mock (Lyman forest, CIV, SiIV; z_em = 2.5) as .acs2."""
    path = tmp_path_factory.mktemp("forest") / "forest.acs2"
    mock_session().save(str(path))
    return path


def test_recipe_tool_defaults_are_the_gui_defaults():
    async def go():
        async with Client(build_server()) as client:
            return (await client.list_tools()).tools
    tools = {t.name: t for t in asyncio.run(go())}
    for name, (menu, recipe) in {
            "continuum_estimate_auto": ("continuum", "estimate_auto"),
            "absorbers_identify_lines": ("absorbers", "identify_lines")}.items():
        schema = {p["name"]: p for p in SCHEMAS[menu][recipe]["params"]}
        props = tools[name].input_schema["properties"]
        for param, prop in props.items():
            if param in schema:
                assert prop["default"] == schema[param]["default"], (name, param)
                assert prop["description"] == schema[param]["doc"]
        assert tools[name].description.startswith(SCHEMAS[menu][recipe]["brief"])


def test_continuum(forest_archive):
    _, result = calls(build_server(),
                      ("session_open", {"path": str(forest_archive)}),
                      ("continuum_estimate_auto", {"session_id": "forest"}))
    assert not result.is_error
    info = result.structured_content
    assert info["spectrum"]["has_continuum"]
    assert info["log_entries"] == 1
    cont = info["continuum"]
    assert 0 < cont["absorbed_fraction_red"] < cont["absorbed_fraction_forest"] < 1
    assert cont["snr_median"] == pytest.approx(30, rel=0.2)  # mock S/N 30


def test_continuum_is_logged_and_replays(forest_archive, tmp_path):
    target = tmp_path / "cont.acs2"
    calls(build_server(),
          ("session_open", {"path": str(forest_archive)}),
          ("continuum_estimate_auto", {"session_id": "forest", "kappa": 3.0}),
          ("session_save", {"session_id": "forest", "path": str(target)}))
    saved = load_session_from_file(str(target), "cont", None, "auto")
    (entry,) = saved.log_manager.to_dict()["entries"]
    assert entry["recipe_name"] == "estimate_auto"
    assert entry["params"]["kappa"] == "3.0"
    # The logged parameters give the same continuum when replayed
    original = load_session_from_file(str(forest_archive), "o", None, "auto")
    replayed = original.continuum.estimate_auto(**entry["params"])
    assert replayed.spec.cont.value == pytest.approx(saved.spec.cont.value)


def test_continuum_without_z_em_warns(forest_archive):
    _, result = calls(build_server(),
                      ("session_open", {"path": str(forest_archive), "z_em": 0}),
                      ("continuum_estimate_auto", {"session_id": "forest"}))
    info = result.structured_content
    assert any("z_em" in w for w in info["warnings"])
    assert info["continuum"]["absorbed_fraction_forest"] is None


def test_identify_lines(forest_archive):
    *_, result = calls(
        build_server(),
        ("session_open", {"path": str(forest_archive)}),
        ("continuum_estimate_auto", {"session_id": "forest"}),
        ("absorbers_identify_lines", {"session_id": "forest",
                                      "max_listed": 1000}))
    assert not result.is_error
    info = result.structured_content
    assert info["log_entries"] == 2
    ident = info["identification"]
    assert not ident["truncated"]
    assert ident["n_systems"] == len(ident["systems"])
    assert ident["n_regions"] == len(ident["regions"])
    assert sum(ident["systems_by_series"].values()) == ident["n_systems"]
    assert 0 < ident["absorbed_fraction_identified"] <= 1
    assert info["spectrum"]["n_components"] == 0  # no component created
    civ = [s["z"] for s in ident["systems"] if s["series"] == "CIV"]
    for z in MOCK_CIV_Z:
        assert min(abs(c - z) for c in civ) / (1 + z) * C_KMS < 20, z
    scores = [max(c["score"] for c in r["candidates"]) for r in ident["regions"]]
    assert scores == sorted(scores, reverse=True)
    series = [s["series"] for s in ident["systems"]]
    assert series.index("Ly_a") > max(i for i, x in enumerate(series)
                                      if x == "CIV")


def test_identify_lines_doublets_counted_twice(forest_archive):
    # Characterization (F49): a doublet gives one candidate system per
    # line, at redshifts ~1 km/s apart
    *_, result = calls(
        build_server(),
        ("session_open", {"path": str(forest_archive)}),
        ("continuum_estimate_auto", {"session_id": "forest"}),
        ("absorbers_identify_lines", {"session_id": "forest"}))
    civ = sorted(s["z"] for s in result.structured_content["identification"]
                 ["systems"] if s["series"] == "CIV")
    close = [b - a for a, b in zip(civ, civ[1:]) if (b - a) / (1 + b) * C_KMS < 2]
    assert close


def test_identify_lines_lists_at_most_max_listed(forest_archive):
    *_, result = calls(
        build_server(),
        ("session_open", {"path": str(forest_archive)}),
        ("absorbers_identify_lines", {"session_id": "forest",
                                      "max_listed": 2}))
    ident = result.structured_content["identification"]
    assert ident["truncated"]
    assert len(ident["regions"]) == len(ident["systems"]) == 2
    assert ident["n_systems"] > 2


def test_identify_lines_without_z_em(forest_archive):
    _, result = calls(build_server(),
                      ("session_open", {"path": str(forest_archive), "z_em": 0}),
                      ("absorbers_identify_lines", {"session_id": "forest"}))
    assert "z_em must be set" in error_text(result)


def test_identification_summary_after_a_round_trip(forest_archive, tmp_path):
    # The maps are JSON strings after identify_lines, dicts after an .acs2
    # round trip (F27)
    session = load_session_from_file(str(forest_archive), "f", None, "auto")
    identified = session.continuum.estimate_auto().absorbers.identify_lines()
    path = tmp_path / "ids.acs2"
    identified.save(str(path))
    reloaded = load_session_from_file(str(path), "r", None, "auto")
    before = identification_summary(identified, 1000)
    after = identification_summary(reloaded, 1000)
    assert after.n_systems == before.n_systems
    assert after.n_regions == before.n_regions


def test_recipe_on_unknown_session():
    result = call(build_server(), "continuum_estimate_auto", session_id="nope")
    assert "No session 'nope'" in error_text(result)


@pytest.mark.slow
def test_j0103_continuum_and_identification():
    *_, result = calls(
        build_server(),
        ("session_open", {"path": str(J0103_SESSION), "z_em": J0103_Z_EM}),
        ("continuum_estimate_auto", {"session_id": J0103_SESSION.stem}),
        ("absorbers_identify_lines", {"session_id": J0103_SESSION.stem}))
    ident = result.structured_content["identification"]
    assert {"CIV", "MgII", "Ly_a"} <= set(ident["systems_by_series"])
    assert ident["truncated"]
    assert 0 < ident["absorbed_fraction_identified"] < 1
