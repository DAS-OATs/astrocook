"""The gaps of the collection of recipes in the MCP server (phase 6-ter,
PR 5b, decision 6): declared by the agent, or seen by the server."""
import asyncio
import json

import pytest

pytest.importorskip("mcp")

from mcp import Client  # noqa: E402

from astrocook import gaps  # noqa: E402
from astrocook.session.session import Session, load_session_from_file  # noqa: E402
from astrocook.mcp_server.server import build_server  # noqa: E402
from astrocook.mocks import MockComponent, MockTruth, generate  # noqa: E402


def calls(*steps, server=None):
    server = server or build_server()

    async def go():
        async with Client(server) as client:
            return [await client.call_tool(tool, args) for tool, args in steps]
    return asyncio.run(go())


def ok(result):
    assert not result.is_error, result.content[0].text
    return result.structured_content


def lines():
    path = gaps.gaps_path()
    return [json.loads(x) for x in path.read_text().splitlines()] if path.is_file() else []


@pytest.fixture(scope="module")
def mock_archive(tmp_path_factory):
    truth = MockTruth(components=[MockComponent("CIV", 1.9, 13.8, 10.0)],
                      lam_min_nm=440.0, lam_max_nm=460.0, dv_kms=2.5,
                      resol=45000.0, snr=50.0, z_em=2.0, seed=1)
    path = tmp_path_factory.mktemp("mock") / "mock.acs2"
    Session("mock", None, spec=generate(truth)).save(str(path))
    return path


def open_mock(archive, z_em=2.0):
    return ("session_open", {"path": str(archive), "z_em": z_em, "resol": 45000})


def test_a_declared_gap_is_in_the_log_and_the_file(mock_archive, tmp_path):
    saved = tmp_path / "gap.acs2"
    results = calls(
        open_mock(mock_archive),
        ("session_note", {"session_id": "mock", "kind": "gap",
                          "text": "Needed z_em; measured the Ly-a peak with my code."}),
        ("session_note", {"session_id": "mock", "text": "A plain note."}),
        ("session_history", {"session_id": "mock"}),
        ("session_save", {"session_id": "mock", "path": str(saved)}))
    assert ok(results[1])["note"]["kind"] == "gap" and ok(results[1])["unsaved"]
    (line,) = lines()
    assert line["kind"] == "declared" and line["session"] == "mock"
    assert line["text"] == "Needed z_em; measured the Ly-a peak with my code."
    assert [n["kind"] for s in ok(results[3])["states"] for n in s["notes"]] == ["gap", "note"]
    log = load_session_from_file(str(saved), "g", None, "auto").log_manager
    assert [n.kind for n in log.notes] == ["gap", "note"]
    assert [e.recipe_name for e in log.entries] == ["set_properties"]  # the replay is the same


def test_a_kind_of_note_is_one_of_two(mock_archive):
    _, bad = calls(open_mock(mock_archive),
                   ("session_note", {"session_id": "mock", "text": "x", "kind": "gap_seen"}))
    assert bad.is_error and lines() == []


def test_a_prerequisite_that_no_recipe_estimates_is_a_gap(mock_archive, tmp_path):
    results = calls(
        open_mock(mock_archive, z_em=0),
        ("session_save", {"session_id": "mock", "path": str(tmp_path / "s.acs2")}),
        ("continuum_estimate_auto", {"session_id": "mock"}),
        ("continuum_estimate_auto", {"session_id": "mock"}),
        ("session_list", {}),
        ("session_history", {"session_id": "mock"}),
        ("session_close", {"session_id": "mock"}))
    assert results[2].is_error and results[3].is_error
    # Once per session, and no note of the agent
    (line,) = lines()
    assert line["kind"] == "prerequisite"
    assert line["text"] == ("continuum_estimate_auto needs the emission redshift (z_em), "
                            "which no recipe estimates from the data.")
    # Kept with the session, but nothing is lost by closing it
    notes = [n for s in ok(results[5])["states"] for n in s["notes"]]
    assert [n["kind"] for n in notes] == ["gap_seen"]
    assert not ok(results[4])["sessions"][0]["unsaved"]
    assert ok(results[6])["closed"]


def test_a_unsaved_session_stays_unsaved_after_a_gap(mock_archive):
    results = calls(
        open_mock(mock_archive, z_em=0),
        ("flux_smooth", {"session_id": "mock", "sigma_kms": 40}),
        ("continuum_estimate_auto", {"session_id": "mock"}),
        ("session_list", {}))
    assert ok(results[3])["sessions"][0]["unsaved"]
    assert len(lines()) == 1


def test_mistakes_of_the_agent_are_not_gaps(mock_archive):
    results = calls(
        open_mock(mock_archive),
        ("flux_smooth", {"session_id": "mock", "sigma_kms": "fifty"}),
        ("flux_smooth", {"session_id": "mock", "sigma_kms": 40, "nope": 1}),
        ("flux_smooth", {"session_id": "nope", "sigma_kms": 40}))
    assert all(r.is_error for r in results[1:])
    assert lines() == []


def test_an_unexpected_error_on_a_valid_request_is_a_gap(mock_archive, monkeypatch):
    from astrocook.mcp_server import tools
    run = tools.run_registry

    def boom(session, recipe, **kwargs):
        if recipe != "flux.smooth":
            return run(session, recipe, **kwargs)
        raise TypeError("/home/someone/secret/path.fits is not what we expected")
    monkeypatch.setattr("astrocook.mcp_server.tools.run_registry", boom)
    results = calls(
        open_mock(mock_archive),
        ("flux_smooth", {"session_id": "mock", "sigma_kms": 40}),
        ("flux_smooth", {"session_id": "mock", "sigma_kms": 20}))
    assert results[1].is_error and results[2].is_error
    (line,) = lines()  # once per session
    assert line["kind"] == "failed"
    assert line["text"] == "flux_smooth failed on a valid request: TypeError."
    assert "secret" not in json.dumps(line)  # the message is not recorded


def test_a_recipe_that_fails_with_a_reason_is_not_a_gap(mock_archive):
    # RecipeFailed: the data do not allow it, and the agent has the reason
    result, = calls(open_mock(mock_archive),
                    ("flux_rebin", {"session_id": "mock", "xstart": 460.0, "xend": 440.0}))[1:]
    assert result.is_error and lines() == []


def test_an_undone_recipe_is_a_gap(mock_archive):
    results = calls(
        open_mock(mock_archive),
        ("flux_smooth", {"session_id": "mock", "sigma_kms": 40,
                         "reason": "to see the lines"}),
        ("session_undo", {"session_id": "mock", "reason": "it blurred the doublet"}),
        ("flux_smooth", {"session_id": "mock", "sigma_kms": 20}),
        ("session_undo", {"session_id": "mock"}),
        ("session_redo", {"session_id": "mock"}))
    first, second = lines()
    assert first["kind"] == "undone"
    assert first["text"] == ("flux_smooth undone. It was run because: to see the lines "
                             "It was undone because: it blurred the doublet")
    assert second["text"] == "flux_smooth undone."
    assert ok(results[5])["moved"] == ["smooth"]  # redo is no gap
    assert len(lines()) == 2


@pytest.fixture(scope="module")
def physical_archive(tmp_path_factory):
    """The small mock in physical units: a spectrum without a continuum."""
    import astropy.units as au

    from astrocook.core.structures import DataColumn, SpectrumData
    from astrocook.core.spectrum import Spectrum
    truth = MockTruth(components=[MockComponent("CIV", 1.9, 13.8, 10.0)],
                      lam_min_nm=440.0, lam_max_nm=460.0, dv_kms=2.5,
                      resol=45000.0, snr=50.0, z_em=2.0, seed=1)
    mock = generate(truth, truth_columns=False)
    d, k = mock._data, 1e-16 * (mock.x.value / 450.0) ** -1.5
    flux = au.erg / au.cm**2 / au.s / au.nm
    spec = Spectrum(SpectrumData(
        x=d.x, xmin=d.xmin, xmax=d.xmax, y=DataColumn(mock.y.value * k, flux),
        dy=DataColumn(mock.dy.value * k, flux), aux_cols=dict(d.aux_cols), z_em=d.z_em))
    path = tmp_path_factory.mktemp("physical") / "physical.acs2"
    Session("physical", None, spec=spec).save(str(path))
    return path


def test_a_prerequisite_that_a_recipe_provides_is_not_a_gap(physical_archive):
    # The continuum: the agent runs the recipe, or asks for it
    _, result = calls(("session_open", {"path": str(physical_archive)}),
                      ("absorbers_doublets_auto", {"session_id": "physical"}))
    assert result.is_error and "no continuum" in result.content[0].text
    assert lines() == []


def test_provided_steps_are_not_the_recipes_of_the_agent(physical_archive):
    results = calls(
        ("session_open", {"path": str(physical_archive)}),
        ("absorbers_doublets_auto", {"session_id": "physical", "provide_missing": True}),
        ("session_undo", {"session_id": "physical", "steps": 2}))
    assert ok(results[2])["moved"] == ["doublets_auto", "estimate_auto"]
    (line,) = lines()
    assert line["text"] == "absorbers_doublets_auto undone."
    # A copy made after the steps keeps the knowledge
    results = calls(
        ("session_open", {"path": str(physical_archive)}),
        ("absorbers_doublets_auto", {"session_id": "physical", "provide_missing": True}),
        ("session_duplicate", {"session_id": "physical"}),
        ("session_undo", {"session_id": "physical-1", "steps": 2}))
    assert len(lines()) == 2 and lines()[1]["text"] == "absorbers_doublets_auto undone."


def test_gaps_summary(mock_archive):
    calls(open_mock(mock_archive),
          ("session_note", {"session_id": "mock", "kind": "gap", "text": "Needed A."}),
          ("flux_smooth", {"session_id": "mock", "sigma_kms": 40}),
          ("session_undo", {"session_id": "mock"}))
    summary, only, bad = calls(
        ("gaps_summary", {}), ("gaps_summary", {"kind": "undone"}),
        ("gaps_summary", {"since": "yesterday"}))
    data = ok(summary)
    assert data["enabled"] and data["skipped_lines"] == 0
    assert sorted(g["kind"] for g in data["groups"]) == ["declared", "undone"]
    assert [g["kind"] for g in ok(only)["groups"]] == ["undone"]
    assert bad.is_error and "not a date" in bad.content[0].text


def test_gaps_summary_with_the_record_off(monkeypatch):
    monkeypatch.setenv(gaps.ENV_VAR, "off")
    result, = calls(("gaps_summary", {}))
    assert ok(result) == {"enabled": False, "groups": [], "n_groups": 0, "skipped_lines": 0}
