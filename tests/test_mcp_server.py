"""MCP server (pilot, ROADMAP.md section 4.1): tools called through the SDK
client, connected to the server in process."""
import asyncio
import logging
import sys
import threading
from importlib.metadata import entry_points
from pathlib import Path

import pytest

pytest.importorskip("mcp")

from mcp import Client, StdioServerParameters  # noqa: E402
from mcp.server.mcpserver.exceptions import ToolError  # noqa: E402

from astrocook import __version__  # noqa: E402
from astrocook.core.session import SessionV2, load_session_from_file  # noqa: E402
from astrocook.mcp_server.errors import run_recipe  # noqa: E402
from astrocook.core.structures import ComponentDataV2  # noqa: E402
from astrocook.mcp_server.outputs import (  # noqa: E402
    component_flags, identification_summary)
from astrocook.mcp_server import server as server_module  # noqa: E402
from astrocook.mcp_server.server import SCHEMAS, build_server  # noqa: E402
from astrocook.mocks import MockComponent, MockTruth, generate  # noqa: E402

from conftest import J0103_SESSION, J0103_Z_EM  # noqa: E402
from full_chain import mock_session, mock_truth  # noqa: E402

TOOLS = {"session_open", "session_save", "continuum_estimate_auto",
         "absorbers_identify_lines", "absorbers_doublets_auto"}


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
            "absorbers_identify_lines": ("absorbers", "identify_lines"),
            "absorbers_doublets_auto": ("absorbers", "doublets_auto")}.items():
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


def test_identify_lines_one_candidate_per_doublet(forest_archive):
    # F49: a doublet gave one candidate system per line, at redshifts ~1
    # km/s apart; now only the stronger line gives the candidates
    *_, result = calls(
        build_server(),
        ("session_open", {"path": str(forest_archive)}),
        ("continuum_estimate_auto", {"session_id": "forest"}),
        ("absorbers_identify_lines", {"session_id": "forest"}))
    civ = sorted(s["z"] for s in result.structured_content["identification"]
                 ["systems"] if s["series"] == "CIV")
    close = [b - a for a, b in zip(civ, civ[1:]) if (b - a) / (1 + b) * C_KMS < 2]
    assert civ and not close


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


# --- Doublet fit -------------------------------------------------------------

#: Largest pull (fit - true) / error accepted for the true CIV components
MAX_PULL = 5.0


def doublets(archive, *more):
    """Open, continuum, doublets_auto (and more steps) on a new server."""
    return calls(build_server(),
                 ("session_open", {"path": str(archive)}),
                 ("continuum_estimate_auto", {"session_id": archive.stem}),
                 ("absorbers_doublets_auto", {"session_id": archive.stem}),
                 *more)


def test_doublets_auto_recovers_the_mock_civ(forest_archive):
    result = doublets(forest_archive)[-1]
    assert not result.is_error
    info = result.structured_content
    fit = info["fit"]
    assert info["log_entries"] == 2  # continuum, doublets
    assert fit["n_removed"] == 0
    assert fit["n_added"] == fit["n_components"] == len(fit["components"])
    assert info["spectrum"]["n_components"] == fit["n_components"]
    zs = [c["z"] for c in fit["components"]]
    assert zs == sorted(zs)
    for true in (c for c in mock_truth().components if c.series == "CIV"):
        comp = min(fit["components"], key=lambda c: abs(c["z"] - true.z))
        assert comp["series"] == "CIV"
        assert abs(comp["z"] - true.z) / (1 + true.z) * C_KMS < 2, true
        for par, value in (("logN", true.logN), ("b", true.b)):
            pull = (comp[par] - value) / comp["d" + par]
            assert abs(pull) < MAX_PULL, (true, par, pull)
        assert not comp["flags"], true
        assert comp["chi2r"] == pytest.approx(1, abs=0.5)


def test_doublets_auto_adds_spurious_components(forest_archive):
    # Characterization (F49): the second candidate of a doublet adds weak
    # components near the real ones (1 and 19 km/s here), with large
    # errors (with kappa = 3 in the continuum one ends at the lower bound
    # of b)
    comps = doublets(forest_archive)[-1].structured_content["fit"]["components"]
    true_z = [c.z for c in mock_truth().components if c.series == "CIV"]
    matched = {min(comps, key=lambda c: abs(c["z"] - z))["uuid"] for z in true_z}
    spurious = [c for c in comps if c["uuid"] not in matched]
    assert len(comps) == 6 and len(matched) == 4  # 4 true components
    for c in spurious:
        assert c["dlogN"] > 0.5


def test_doublets_auto_replaces_the_red_side_components(forest_archive):
    first, second = doublets(
        forest_archive,
        ("absorbers_doublets_auto", {"session_id": "forest",
                                     "max_listed": 1}))[-2:]
    n = first.structured_content["fit"]["n_components"]
    fit = second.structured_content["fit"]
    assert fit["n_removed"] == n
    assert fit["n_components"] == fit["n_added"]
    assert len(fit["components"]) == 1 and fit["truncated"]


def test_doublets_auto_replaces_the_mask_of_the_continuum(forest_archive):
    # Characterization (F50): components_auto makes a new absorbed mask
    # with the default parameters of find_absorbed, so the mask of an
    # estimate_auto with other parameters is lost, and clears the region
    # left out by region_limit
    session = load_session_from_file(str(forest_archive), "f", None, "auto")
    cont = session.continuum.estimate_auto(kappa="3.0")
    fitted = cont.absorbers.doublets_auto()
    mask_cont = cont.spec.get_column("abs_mask").value
    mask_fit = fitted.spec.get_column("abs_mask").value
    assert (mask_cont != mask_fit).any()
    assert fitted.spec.cont.value == pytest.approx(cont.spec.cont.value)
    # The forest (left out by region_limit) is no longer flagged as absorbed
    x = fitted.spec.x.value
    forest = x < 121.567 * (1 + fitted.spec._data.z_em)
    assert mask_cont[forest].any() and not mask_fit[forest].any()


@pytest.mark.parametrize("kwargs, flags", [
    ({}, []),
    ({"logN": 10.0}, ["logN_at_lower_bound"]),
    ({"logN": 23.0}, ["logN_at_upper_bound"]),
    ({"b": 1.0}, ["b_at_lower_bound"]),
    ({"b": 50.0}, ["b_at_max_b"]),
    ({"dlogN": float("nan")}, ["nan_errors"]),
])
def test_component_flags(kwargs, flags):
    values = dict(id=1, z=2.0, dz=1e-5, logN=13.5, dlogN=0.1, b=10.0, db=1.0)
    values.update(kwargs)
    assert component_flags(ComponentDataV2(**values), max_b=50.0) == flags


@pytest.mark.slow
def test_j0103_doublets():
    stem = J0103_SESSION.stem
    *_, result = calls(
        build_server(),
        ("session_open", {"path": str(J0103_SESSION), "z_em": J0103_Z_EM,
                          "resol": 45000}),
        ("continuum_estimate_auto", {"session_id": stem}),
        ("absorbers_doublets_auto", {"session_id": stem}))
    fit = result.structured_content["fit"]
    # The archive already holds components: those on the red side go
    assert fit["n_removed"] > 0 and fit["n_added"] > 0
    assert fit["truncated"]
    assert {c["series"] for c in fit["components"]} <= {"CIV", "SiIV", "MgII"}


# --- The server as a separate process (stdio), as a client starts it ---------

def test_console_script_points_to_the_server():
    (ep,) = [ep for ep in entry_points(group="console_scripts")
             if ep.name == "astrocook-mcp"]
    assert ep.value == "astrocook.mcp_server.server:main"
    assert ep.load() is server_module.main


def stdio_calls(params, *steps):
    """List the tools and run the steps on a server in a subprocess."""
    async def go():
        async with Client(params) as client:
            tools = (await client.list_tools()).tools
            return tools, [await client.call_tool(t, a) for t, a in steps]
    return asyncio.run(go())


def test_stdio(mock_archive, tmp_path):
    # Anything written to stdout by Astrocook would break the protocol;
    # the save logs a warning (astropy), which must go to stderr
    params = StdioServerParameters(
        command=sys.executable, args=["-m", "astrocook.mcp_server.server"])
    tools, (opened, saved, missing) = stdio_calls(
        params,
        ("session_open", {"path": str(mock_archive), "z_em": 2.1}),
        ("session_save", {"session_id": "mock",
                          "path": str(tmp_path / "out.acs2")}),
        ("session_open", {"path": str(tmp_path / "missing.fits")}))
    assert {t.name for t in tools} == TOOLS
    assert opened.structured_content["spectrum"]["z_em"] == pytest.approx(2.1)
    assert saved.structured_content["log_entries"] == 1
    assert (tmp_path / "out.acs2").is_file()
    assert "No such file" in error_text(missing)


def test_stdio_console_script(mock_archive):
    script = Path(sys.executable).with_name("astrocook-mcp")
    if not script.exists():
        pytest.skip("astrocook-mcp not installed next to the interpreter")
    tools, (opened,) = stdio_calls(
        StdioServerParameters(command=str(script)),
        ("session_open", {"path": str(mock_archive)}))
    assert {t.name for t in tools} == TOOLS
    assert not opened.is_error
