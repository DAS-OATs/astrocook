"""The tools of the sessions (phase 6-ter, PR 2): list, describe, history,
undo and redo, duplicate, close.

Decision 1 of the phase: the server keeps every state of a session, as the
GUI does, so that what a recipe removed can be had back without running
anything again.
"""
import asyncio

import pytest

pytest.importorskip("mcp")

from mcp import Client  # noqa: E402

from astrocook.session.session import Session, load_session_from_file  # noqa: E402
from astrocook.mcp_server.server import build_server  # noqa: E402
from astrocook.mocks import MockComponent, MockTruth, generate  # noqa: E402


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


def test_list_and_describe(mock_archive):
    _, _, listed, described = calls(
        build_server(), open_mock(mock_archive),
        ("continuum_estimate_auto", {"session_id": "mock"}),
        ("session_list", {}),
        ("session_describe", {"session_id": "mock"}))
    (item,) = ok(listed)["sessions"]
    assert item["session_id"] == "mock" and item["origin"] == "opened"
    assert item["path"] == str(mock_archive)
    # As opened, set_properties, estimate_auto
    assert item["n_states"] == 3 and item["current_state"] == 2
    assert item["unsaved"]
    info = ok(described)
    names = [c["name"] for c in info["columns"]]
    assert {"x", "y", "dy", "cont", "abs_mask"} <= set(names)
    assert info["components_by_series"] == {}
    assert info["log_before"] == []
    assert info["spectrum"]["has_continuum"]


def test_the_log_of_an_archive_comes_before_the_history(mock_archive, tmp_path):
    saved = tmp_path / "with_log.acs2"
    calls(build_server(), open_mock(mock_archive),
          ("session_save", {"session_id": "mock", "path": str(saved)}))
    _, history, described = calls(
        build_server(), ("session_open", {"path": str(saved)}),
        ("session_history", {"session_id": "with_log"}),
        ("session_describe", {"session_id": "with_log"}))
    assert [e["recipe"] for e in ok(described)["log_before"]] == ["set_properties"]
    (state,) = ok(history)["states"]
    assert state["recipe"] is None and state["current"]
    assert not ok(described)["unsaved"]


def test_the_states_after_an_archive_log(mock_archive, tmp_path):
    # Each state made in the server is matched with its own entry; the log
    # of the archive is read only, before them (F109: not continued)
    saved = tmp_path / "with_log.acs2"
    calls(build_server(), open_mock(mock_archive),
          ("continuum_estimate_auto", {"session_id": "mock", "kappa": 3}),
          ("session_save", {"session_id": "mock", "path": str(saved)}))
    _, _, _, history, undone, described = calls(
        build_server(), ("session_open", {"path": str(saved)}),
        ("flux_smooth", {"session_id": "with_log", "sigma_kms": 40}),
        ("flux_smooth", {"session_id": "with_log", "sigma_kms": 30}),
        ("session_history", {"session_id": "with_log"}),
        ("session_undo", {"session_id": "with_log"}),
        ("session_describe", {"session_id": "with_log"}))
    states = ok(history)["states"]
    assert [s["recipe"] for s in states] == [None, "smooth", "smooth"]
    assert [s["params"].get("sigma_kms") for s in states] == [None, 40.0, 30.0]
    assert ok(undone)["moved"] == ["smooth"] and ok(undone)["log_entries"] == 1
    assert [e["recipe"] for e in ok(described)["log_before"]] == [
        "set_properties", "estimate_auto"]


def test_history_undo_redo(mock_archive):
    _, _, history, undone, described, redone, again = calls(
        build_server(), open_mock(mock_archive),
        ("continuum_estimate_auto", {"session_id": "mock", "kappa": 3}),
        ("session_history", {"session_id": "mock"}),
        ("session_undo", {"session_id": "mock"}),
        ("session_describe", {"session_id": "mock"}),
        ("session_redo", {"session_id": "mock"}),
        ("session_redo", {"session_id": "mock"}))
    states = ok(history)["states"]
    assert [s["recipe"] for s in states] == [None, "set_properties", "estimate_auto"]
    assert states[2]["params"]["kappa"] == 3.0 and states[2]["current"]
    # The mock has a continuum already (that of the generator)
    assert "column 'cont' changed" in states[2]["summary"]
    assert ok(history)["can_undo"] and not ok(history)["can_redo"]
    u = ok(undone)
    assert u["moved"] == ["estimate_auto"] and u["current_state"] == 1
    assert "abs_mask" not in [c["name"] for c in ok(described)["columns"]]
    # The log follows the state: what is saved now replays without the
    # recipe undone
    assert u["log_entries"] == 1
    r = ok(redone)
    assert r["moved"] == ["estimate_auto"] and r["current_state"] == 2
    assert r["log_entries"] == 2
    assert "Nothing to redo in 'mock'" in error_text(again)


def test_undo_several_steps_and_a_new_recipe_drops_the_redo(mock_archive):
    _, _, _, undone, _, history = calls(
        build_server(), open_mock(mock_archive),
        ("continuum_estimate_auto", {"session_id": "mock"}),
        ("flux_smooth", {"session_id": "mock"}),
        ("session_undo", {"session_id": "mock", "steps": 5}),
        ("flux_rebin", {"session_id": "mock", "dx": 5.0}),
        ("session_history", {"session_id": "mock"}))
    # Stops at the session as opened
    assert ok(undone)["moved"] == ["smooth", "estimate_auto", "set_properties"]
    assert ok(undone)["current_state"] == 0
    assert [s["recipe"] for s in ok(history)["states"]] == [None, "rebin"]


def test_what_a_recipe_removed_is_had_back_with_undo(mock_archive, tmp_path):
    # Decision 1: the result says how; the state before is restored, and
    # what is saved then has the components and the log of that state
    _, first, second, undone, _ = calls(
        build_server(), open_mock(mock_archive),
        ("absorbers_doublets_auto", {"session_id": "mock", "provide_missing": True}),
        ("absorbers_doublets_auto", {"session_id": "mock"}),
        ("session_undo", {"session_id": "mock"}),
        ("session_save", {"session_id": "mock", "path": str(tmp_path / "u.acs2")}))
    added = {c["uuid"] for c in ok(first)["fit"]["components"]}
    assert ok(second)["changes"]["n_components_removed"] == len(added) > 0
    assert "session_undo" in ok(second)["undo_hint"]
    assert ok(first)["undo_hint"] is None
    assert ok(undone)["spectrum"]["n_components"] == len(added)
    saved = load_session_from_file(str(tmp_path / "u.acs2"), "u", None, "auto")
    assert {c.uuid for c in saved.systs.components} == added
    names = [e["recipe_name"] for e in saved.log_manager.to_dict()["entries"]]
    assert names == ["set_properties", "doublets_auto"]


def test_duplicate_is_independent(mock_archive):
    _, _, dup, _, listed, h_orig, h_dup = calls(
        build_server(), open_mock(mock_archive),
        ("continuum_estimate_auto", {"session_id": "mock"}),
        ("session_duplicate", {"session_id": "mock"}),
        ("flux_smooth", {"session_id": "mock-1"}),
        ("session_list", {}),
        ("session_history", {"session_id": "mock"}),
        ("session_history", {"session_id": "mock-1"}))
    d = ok(dup)
    assert d["session_id"] == "mock-1" and d["origin"] == "duplicated"
    assert d["parent"] == "mock" and d["unsaved"]
    assert [s["session_id"] for s in ok(listed)["sessions"]] == ["mock", "mock-1"]
    assert [s["recipe"] for s in ok(h_orig)["states"]] == [None, "set_properties",
                                                            "estimate_auto"]
    assert [s["recipe"] for s in ok(h_dup)["states"]] == [None, "set_properties",
                                                           "estimate_auto", "smooth"]


def test_a_branch_lists_its_origin(mock_archive):
    _, _, listed = calls(
        build_server(), open_mock(mock_archive),
        ("edit_split", {"session_id": "mock", "expression": "x > 450"}),
        ("session_list", {}))
    branch = ok(listed)["sessions"][1]
    assert branch["origin"] == "branched" and branch["parent"] == "mock"
    assert branch["path"] is None and branch["unsaved"]


def test_close_refuses_to_lose_changes(mock_archive, tmp_path):
    _, refused, saved, closed, listed, gone = calls(
        build_server(), open_mock(mock_archive),
        ("session_close", {"session_id": "mock"}),
        ("session_save", {"session_id": "mock", "path": str(tmp_path / "c.acs2")}),
        ("session_close", {"session_id": "mock"}),
        ("session_list", {}),
        ("continuum_estimate_auto", {"session_id": "mock"}))
    assert "has changes not saved" in error_text(refused)
    assert ok(closed) == {"session_id": "mock", "closed": True, "open_sessions": []}
    assert ok(listed)["sessions"] == []
    assert "No session 'mock'" in error_text(gone)


def test_close_discarding_changes_and_a_session_opened_unchanged(mock_archive):
    _, _, discarded, unchanged = calls(
        build_server(), open_mock(mock_archive),
        ("session_open", {"path": str(mock_archive)}),
        ("session_close", {"session_id": "mock", "discard_changes": True}),
        ("session_close", {"session_id": "mock-2"}))
    assert ok(discarded)["open_sessions"] == ["mock-2"]
    assert ok(unchanged)["open_sessions"] == []


def test_undo_back_to_the_saved_state_is_not_unsaved(mock_archive, tmp_path):
    _, _, _, after, _, back = calls(
        build_server(), open_mock(mock_archive),
        ("session_save", {"session_id": "mock", "path": str(tmp_path / "s.acs2")}),
        ("flux_smooth", {"session_id": "mock"}),
        ("session_list", {}),
        ("session_undo", {"session_id": "mock"}),
        ("session_list", {}))
    assert ok(after)["sessions"][0]["unsaved"]
    assert not ok(back)["sessions"][0]["unsaved"]


# --- Reasons and notes (phase 6-ter, PR 5a) ----------------------------------

def test_a_reason_goes_to_the_log(mock_archive, tmp_path):
    saved = tmp_path / "reasons.acs2"
    results = calls(
        build_server(), open_mock(mock_archive),
        ("flux_smooth", {"session_id": "mock", "sigma_kms": 40,
                         "reason": "Too noisy to see the doublet."}),
        ("edit_split", {"session_id": "mock", "expression": "x > 450",
                        "reason": "Only the red half."}),
        ("session_history", {"session_id": "mock"}),
        ("session_history", {"session_id": "mock-1"}),
        ("session_save", {"session_id": "mock", "path": str(saved)}))
    history, branch = ok(results[3])["states"], ok(results[4])["states"]
    assert [s["reason"] for s in history] == [None, None, "Too noisy to see the doublet."]
    assert branch[-1]["recipe"] == "split" and branch[-1]["reason"] == "Only the red half."
    entries = load_session_from_file(str(saved), "r", None, "auto").log_manager.entries
    assert [e.reason for e in entries] == [None, "Too noisy to see the doublet."]
    # Reopened: in the log that came with the archive
    (described,) = calls(build_server(), ("session_open", {"path": str(saved)}),
                         ("session_describe", {"session_id": "reasons"}))[1:]
    assert [e["reason"] for e in ok(described)["log_before"]] == [
        None, "Too noisy to see the doublet."]


def test_a_reason_is_optional_and_checked(mock_archive):
    _, plain, wrong = calls(
        build_server(), open_mock(mock_archive),
        ("flux_smooth", {"session_id": "mock", "reason": None}),
        ("flux_smooth", {"session_id": "mock", "reason": 3}))
    ok(plain)
    assert wrong.is_error


def test_a_note_makes_no_state(mock_archive, tmp_path):
    saved = tmp_path / "noted.acs2"
    results = calls(
        build_server(), open_mock(mock_archive),
        ("session_note", {"session_id": "mock", "text": "z_em from the header."}),
        ("flux_smooth", {"session_id": "mock", "sigma_kms": 40}),
        ("session_note", {"session_id": "mock", "text": "Smoothed: no line visible."}),
        ("session_history", {"session_id": "mock"}),
        ("session_save", {"session_id": "mock", "path": str(saved)}),
        ("session_note", {"session_id": "mock", "text": "After saving."}),
        ("session_list", {}),
        ("session_close", {"session_id": "mock"}))
    first, second = ok(results[1]), ok(results[3])
    assert first["state"] == 1 and first["note"] == {"kind": "note",
                                                     "text": "z_em from the header."}
    assert second["state"] == 2 and second["unsaved"]
    states = ok(results[4])["states"]
    assert len(states) == 3
    assert [[n["text"] for n in s["notes"]] for s in states] == [
        [], ["z_em from the header."], ["Smoothed: no line visible."]]
    # A note is a change to save: the session is not closed without asking
    assert ok(results[7])["sessions"][0]["unsaved"]
    assert "changes not saved" in error_text(results[8])
    log = load_session_from_file(str(saved), "n", None, "auto").log_manager
    assert [n.text for n in log.notes] == ["z_em from the header.",
                                           "Smoothed: no line visible."]
    assert [e.recipe_name for e in log.entries] == ["set_properties", "smooth"]


def test_a_note_is_undone_with_its_state(mock_archive):
    results = calls(
        build_server(), open_mock(mock_archive),
        ("flux_smooth", {"session_id": "mock", "sigma_kms": 40}),
        ("session_note", {"session_id": "mock", "text": "On the smoothed flux."}),
        ("session_undo", {"session_id": "mock"}),
        ("session_note", {"session_id": "mock", "text": "Before smoothing."}),
        ("session_redo", {"session_id": "mock"}),
        ("session_history", {"session_id": "mock"}),
        ("session_undo", {"session_id": "mock"}),
        ("flux_smooth", {"session_id": "mock", "sigma_kms": 20}),
        ("session_history", {"session_id": "mock"}))
    # The note does not drop the redo; the states keep their notes
    assert ok(results[5])["moved"] == ["smooth"]
    # (the undo is a gap seen by the server, PR 5b: not a note of the agent)
    assert [[n["text"] for n in s["notes"] if n["kind"] == "note"]
            for s in ok(results[6])["states"]] == [
        [], ["Before smoothing."], ["On the smoothed flux."]]
    # A new recipe after an undo drops the state undone, and its notes
    states = ok(results[9])["states"]
    assert [s["params"].get("sigma_kms") for s in states] == [None, None, 20.0]
    assert [[n["text"] for n in s["notes"] if n["kind"] == "note"] for s in states] == [
        [], ["Before smoothing."], []]


def test_a_note_needs_a_text(mock_archive):
    _, empty, blank, missing = calls(
        build_server(), open_mock(mock_archive),
        ("session_note", {"session_id": "mock", "text": ""}),
        ("session_note", {"session_id": "mock", "text": "   "}),
        ("session_note", {"session_id": "nope", "text": "x"}))
    assert empty.is_error
    assert "needs a text" in error_text(blank)
    assert "No session 'nope'" in error_text(missing)


def test_notes_of_an_archive(mock_archive, tmp_path):
    # The notes that came with an archive: before an entry of its log, or
    # at the state as opened
    saved = tmp_path / "arch.acs2"
    calls(build_server(), open_mock(mock_archive),
          ("session_note", {"session_id": "mock", "text": "Before smoothing."}),
          ("flux_smooth", {"session_id": "mock", "sigma_kms": 40}),
          ("session_note", {"session_id": "mock", "text": "At the end."}),
          ("session_save", {"session_id": "mock", "path": str(saved)}))
    _, described, history = calls(
        build_server(), ("session_open", {"path": str(saved)}),
        ("session_describe", {"session_id": "arch"}),
        ("session_history", {"session_id": "arch"}))
    before = ok(described)["log_before"]
    assert [(e["recipe"], [n["text"] for n in e["notes_before"]]) for e in before] == [
        ("set_properties", []), ("smooth", ["Before smoothing."])]
    (state,) = ok(history)["states"]
    assert [n["text"] for n in state["notes"]] == ["At the end."]
    assert not ok(described)["unsaved"]


def test_methods_text(mock_archive, tmp_path):
    """session_methods (phase 6-ter, PR 5c): the log as the text of a methods
    section, from the package; it changes nothing."""
    saved = tmp_path / "methods.acs2"
    _, _, _, text, plain, history, _, reopened = calls(
        build_server(), open_mock(mock_archive),
        ("flux_smooth", {"session_id": "mock", "sigma_kms": 40,
                         "reason": "Too noisy to see the doublet."}),
        ("session_note", {"session_id": "mock", "text": "Checked by eye"}),
        ("session_methods", {"session_id": "mock"}),
        ("session_methods", {"session_id": "mock", "menus": False}),
        ("session_history", {"session_id": "mock"}),
        ("session_save", {"session_id": "mock", "path": str(saved)}),
        ("session_open", {"path": str(saved)}))
    out = ok(text)
    assert out["session_id"] == "mock"
    assert out["text"].startswith("The spectrum of mock (")
    smooth = out["paragraphs"][2]
    assert smooth.startswith("Too noisy to see the doublet. The flux, its error and the "
                             "columns like the flux were smoothed with a Gaussian kernel "
                             "of standard deviation 40 km/s (Flux > Smooth spectrum).")
    assert smooth.endswith("Checked by eye.")
    assert out["references"] == ["Cupani G. et al. 2020, Proc. SPIE 11452, 114521U"]
    assert " > " not in ok(plain)["text"]
    assert len(ok(history)["states"]) == 3  # nothing added by the tool
    # An archive opened: its log, without what each recipe changed
    reopened_id = ok(reopened)["session_id"]
    (again,) = calls(build_server(), ("session_open", {"path": str(saved)}),
                     ("session_methods", {"session_id": reopened_id}))[1:]
    assert "For the first two steps, read from the log of the archive" in ok(again)["text"]
