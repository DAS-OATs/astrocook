"""The tools that read a session (phase 6-ter, PR 3): components_list,
spectrum_read, spectrum_plot. They change nothing and log nothing."""
import asyncio
import base64

import numpy as np
import pytest

pytest.importorskip("mcp")

from mcp import Client  # noqa: E402

from astrocook.session.session import Session, load_session_from_file  # noqa: E402
from astrocook.mcp_server.server import build_server  # noqa: E402
from astrocook.mocks import MockComponent, MockTruth, generate  # noqa: E402

C_KMS = 299792.458


def calls(server, *steps):
    async def go():
        async with Client(server) as client:
            return [await client.call_tool(tool, args) for tool, args in steps]
    return asyncio.run(go())


def error_text(result):
    assert result.is_error
    return result.content[0].text


def ok(result):
    assert not result.is_error, result.content[0].text
    return result.structured_content


@pytest.fixture(scope="module")
def mock_archive(tmp_path_factory):
    """A small mock (one CIV doublet, z_em = 2.0, R = 45000) as .acs2."""
    truth = MockTruth(components=[MockComponent("CIV", 1.9, 13.8, 10.0)],
                      lam_min_nm=440.0, lam_max_nm=460.0, dv_kms=2.5,
                      resol=45000.0, snr=50.0, z_em=2.0, seed=1)
    path = tmp_path_factory.mktemp("mock") / "mock.acs2"
    Session("mock", None, spec=generate(truth)).save(str(path))
    return path


def open_mock(archive):
    return ("session_open", {"path": str(archive), "z_em": 2.0, "resol": 45000})


def add(series, z, logN=13.0):
    return ("absorbers_add_component", {"session_id": "mock", "series": series, "z": z,
                                        "logN": logN})


# Two CIV ~10 km/s apart (one group), a third far away (its own group)
Z1, Z2, Z3 = 1.9, 1.9 + 10 / C_KMS * 2.9, 1.95


def test_components_list(mock_archive):
    *_, listed, civ, page, history = calls(
        build_server(), open_mock(mock_archive), add("CIV", Z1), add("CIV", Z2),
        add("CIV", Z3),
        ("components_list", {"session_id": "mock"}),
        ("components_list", {"session_id": "mock", "z_min": 1.92}),
        ("components_list", {"session_id": "mock", "offset": 1, "max_listed": 1}),
        ("session_history", {"session_id": "mock"}))
    info = ok(listed)
    assert info["n_components"] == info["n_matched"] == 3 and not info["truncated"]
    comps = info["components"]
    # add_component fits the component near the z given
    assert [c["z"] for c in comps] == pytest.approx([Z1, Z2, Z3], abs=3e-4)
    assert info["n_groups"] == 2
    assert comps[0]["group"] == comps[1]["group"] != comps[2]["group"]
    dv = (comps[1]["z"] - comps[0]["z"]) / (1 + comps[0]["z"]) * C_KMS
    assert comps[0]["dv_kms"] == 0 and comps[1]["dv_kms"] == pytest.approx(dv)
    assert comps[2]["dv_kms"] == 0
    assert all(set(c) >= {"uuid", "logN", "dlogN", "b", "db", "chi2r", "flags",
                          "constraints"} for c in comps)
    assert [c["uuid"] for c in ok(civ)["components"]] == [comps[2]["uuid"]]
    p = ok(page)
    assert [c["uuid"] for c in p["components"]] == [comps[1]["uuid"]] and p["truncated"]
    # Nothing logged by the reading tools
    assert len(ok(history)["states"]) == 5


def test_components_list_filters_by_series_and_flags(mock_archive):
    # doublets_auto replaces the components of the red side: MgII after it
    # (a SiIV at z = 1.9, outside the spectrum, was fitted on the whole
    # spectrum; it is refused since F111, phase 6-ter PR 7-quater)
    *_, fitted, _, by_series, flagged = calls(
        build_server(), open_mock(mock_archive),
        ("absorbers_doublets_auto", {"session_id": "mock"}), add("MgII", 0.6),
        ("components_list", {"session_id": "mock", "series": "SiIV, MgII"}),
        ("components_list", {"session_id": "mock", "flagged_only": True}))
    assert [c["series"] for c in ok(by_series)["components"]] == ["MgII"]
    assert all(c["flags"] for c in ok(flagged)["components"])
    n_flagged = sum(bool(c["flags"]) for c in ok(fitted)["fit"]["components"])
    assert ok(flagged)["n_matched"] >= n_flagged


def test_components_list_shows_the_constraints(mock_archive):
    async def go():
        async with Client(build_server()) as client:
            for tool, args in (open_mock(mock_archive), add("CIV", Z1), add("CIV", Z2)):
                ok(await client.call_tool(tool, args))
            first, second = ok(await client.call_tool(
                "components_list", {"session_id": "mock"}))["components"]
            assert first["constraints"] == second["constraints"] == {}
            for args in ({"uuid": first["uuid"], "param": "b", "is_free": False},
                         {"uuid": second["uuid"], "param": "z",
                          "target_uuid": first["uuid"],
                          "expression": f"p['{first['uuid']}'].z"}):
                ok(await client.call_tool("absorbers_update_constraint",
                                          {"session_id": "mock", **args}))
            return ok(await client.call_tool(
                "components_list", {"session_id": "mock"}))["components"]
    first, second = asyncio.run(go())
    assert first["constraints"] == {"b": "fixed"}
    assert second["constraints"]["z"].startswith("linked: ")
    assert first["uuid"] in second["constraints"]["z"]


def test_spectrum_read(mock_archive):
    _, whole, part, cols, bad_col, bad_range, empty = calls(
        build_server(), ("session_open", {"path": str(mock_archive)}),
        ("spectrum_read", {"session_id": "mock", "max_points": 100}),
        ("spectrum_read", {"session_id": "mock", "x_min": 449, "x_max": 450,
                           "max_points": 5000}),
        ("spectrum_read", {"session_id": "mock", "x_min": 449, "x_max": 450,
                           "columns": ["y"]}),
        ("spectrum_read", {"session_id": "mock", "columns": ["nope"]}),
        ("spectrum_read", {"session_id": "mock", "x_min": 450, "x_max": 449}),
        ("spectrum_read", {"session_id": "mock", "x_min": 900, "x_max": 901}))
    session = load_session_from_file(str(mock_archive), "m", None, "auto")
    x, y = session.spec.x.value, session.spec.y.value
    w = ok(whole)
    assert w["n_pixels"] == len(x) and w["binned_by"] == int(np.ceil(len(x) / 100))
    # The names of the GUI by default (phase 6-quater, PR 8, F120)
    assert len(w["columns"]["λ"]) <= 100
    assert set(w["columns"]) == {"λ", "F", "dF", "cont"}  # no model yet
    # The mean of the first block
    assert w["columns"]["F"][0] == pytest.approx(np.mean(y[:w["binned_by"]]))
    p = ok(part)
    sel = (x >= 449) & (x <= 450)
    assert p["binned_by"] == 1 and p["columns"]["F"] == pytest.approx(list(y[sel]))
    assert p["x_unit"] == "nm" and p["units"]["λ"] == "nm"
    # The axis as λ when it was not asked
    assert set(ok(cols)["columns"]) == {"λ", "y"}
    assert "No column nope" in error_text(bad_col) and "The columns: λ" in error_text(bad_col)
    assert "Empty range" in error_text(bad_range)
    assert "No pixel between 900" in error_text(empty)


def test_spectrum_read_accepts_the_names_of_the_gui(mock_archive):
    # F120 (trial 3 of PR 7): the expression recipes say F and dF
    _, gui, internal = calls(
        build_server(), ("session_open", {"path": str(mock_archive)}),
        ("spectrum_read", {"session_id": "mock", "x_min": 449, "x_max": 450,
                           "columns": ["F", "dF", "cont"]}),
        ("spectrum_read", {"session_id": "mock", "x_min": 449, "x_max": 450,
                           "columns": ["x", "y", "dy", "cont"]}))
    g, i = ok(gui), ok(internal)
    # the values under the names asked (the axis as λ if not asked)
    assert set(g["columns"]) == {"λ", "F", "dF", "cont"}
    assert set(i["columns"]) == {"x", "y", "dy", "cont"}
    assert g["columns"]["λ"] == i["columns"]["x"]
    assert g["columns"]["F"] == i["columns"]["y"]
    assert g["columns"]["dF"] == i["columns"]["dy"]
    assert g["units"]["F"] == i["units"]["y"]


def test_spectrum_plot(mock_archive):
    _, _, plot, normalized = calls(
        build_server(), open_mock(mock_archive), add("CIV", Z1),
        ("spectrum_plot", {"session_id": "mock", "x_min": 445, "x_max": 452}),
        ("spectrum_plot", {"session_id": "mock", "normalized": True}))
    for result in (plot, normalized):
        assert not result.is_error, result.content[0].text
        image, text = result.content
        assert image.type == "image" and image.mime_type == "image/png"
        assert base64.b64decode(image.data)[:8] == b"\x89PNG\r\n\x1a\n"
    assert normalized.content[1].text.endswith("(normalized)")


def test_what_a_recipe_reports_reaches_the_agent_in_details(mock_archive):
    # F118, F119 (trial 3 of PR 7): the same for every recipe, from
    # RecipeResult.extra, no recipe named in the server
    _, _, ident = calls(
        build_server(), open_mock(mock_archive),
        ("continuum_estimate_auto", {"session_id": "mock"}),
        ("absorbers_identify_lines", {"session_id": "mock", "score_threshold": 0.5}),
        )
    assert isinstance(ok(ident)["details"]["dropped_candidates"], list)

    async def optimize():
        async with Client(build_server()) as client:
            await client.call_tool(*open_mock(mock_archive))
            await client.call_tool("continuum_estimate_auto", {"session_id": "mock"})
            await client.call_tool(*add("CIV", Z1))
            listed = ok(await client.call_tool("components_list", {"session_id": "mock"}))
            return await client.call_tool(
                "absorbers_optimize_system",
                {"session_id": "mock", "uuid": listed["components"][0]["uuid"]})
    optimization = ok(asyncio.run(optimize()))["details"]["optimization"]
    assert optimization["stopped"] and isinstance(optimization["tried"], list)
