"""The GUI launched by the MCP server (phase 6-ter, PR 5-bis; decision 7 of
the phase, option b): ``gui_open`` saves a session to an archive and launches
the GUI on it; the server saves the archive again after every change, and the
GUI follows it.

The first test is a characterization test written before ``gui_open``:
``session_save`` writes the archive only.
"""
import asyncio
import gc
from pathlib import Path

import astropy.units as au
import pytest

pytest.importorskip("mcp")

from mcp import Client  # noqa: E402

from astrocook import follow  # noqa: E402
from astrocook.core.session import SessionV2, load_session_from_file  # noqa: E402
from astrocook.core.spectrum import SpectrumV2  # noqa: E402
from astrocook.core.structures import DataColumnV2, SpectrumDataV2  # noqa: E402
from astrocook.mcp_server.sessions import SessionStore  # noqa: E402
from astrocook.mcp_server.server import build_server  # noqa: E402
from astrocook.mocks import MockComponent, MockTruth, generate  # noqa: E402


def calls(server, *steps):
    async def go():
        async with Client(server) as client:
            return [await client.call_tool(tool, args) for tool, args in steps]
    return asyncio.run(go())


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
    SessionV2("mock", None, spec=generate(truth)).save(str(path))
    return path


def open_mock(archive):
    return ("session_open", {"path": str(archive), "z_em": 2.0, "resol": 45000})


# --- Characterization -------------------------------------------------------------------------

def test_session_save_writes_the_archive_only(mock_archive, tmp_path):
    out = tmp_path / "out"
    out.mkdir()
    _, saved = calls(build_server(), open_mock(mock_archive),
                     ("session_save", {"session_id": "mock", "path": str(out / "a.acs2")}))
    ok(saved)
    assert [p.name for p in out.iterdir()] == ["a.acs2"]


# --- gui_open and the archive followed -------------------------------------------------------

@pytest.fixture
def launched(monkeypatch):
    """The GUI available, and its launches recorded instead of run."""
    commands = []

    def launch(path):
        commands.append(str(path))
        return ["astrocook", str(path)]
    monkeypatch.setattr(follow, "gui_unavailable_reason", lambda: None)
    monkeypatch.setattr(follow, "launch_gui", launch)
    return commands


def _state_of(path):
    """The session in a followed archive (its spectrum and its log)."""
    return load_session_from_file(str(path), "followed", None, "auto")


def test_gui_open_writes_the_archive_and_its_marker(mock_archive, launched):
    _, opened = calls(build_server(), open_mock(mock_archive),
                      ("gui_open", {"session_id": "mock"}))
    info = ok(opened)
    path = info["path"]
    assert launched == [path] and info["command"] == ["astrocook", path]
    assert path.endswith("mock.acs2")
    marker = follow.read_marker(path)
    assert marker.session_id == "mock" and marker.writer_alive()
    # As opened, with set_properties in its log
    assert [e.recipe_name for e in _state_of(path).log_manager.entries] == ["set_properties"]


def test_the_archive_follows_every_change(mock_archive, launched):
    server = build_server()
    _, opened, _, after_recipe = calls(
        server, open_mock(mock_archive), ("gui_open", {"session_id": "mock"}),
        ("continuum_estimate_auto", {"session_id": "mock"}),
        ("session_list", {}))
    path = ok(opened)["path"]
    session = _state_of(path)
    assert "abs_mask" in session.spec.t.colnames
    assert [e.recipe_name for e in session.log_manager.entries] == [
        "set_properties", "estimate_auto"]
    # The archive is not the user's save: the session stays unsaved
    assert ok(after_recipe)["sessions"][0]["unsaved"]

    calls(server, ("session_undo", {"session_id": "mock"}))
    session = _state_of(path)
    assert "abs_mask" not in session.spec.t.colnames
    assert [e.recipe_name for e in session.log_manager.entries] == ["set_properties"]

    calls(server, ("session_redo", {"session_id": "mock"}),
          ("session_note", {"session_id": "mock", "text": "Looked at it."}))
    session = _state_of(path)
    assert "abs_mask" in session.spec.t.colnames
    # With the gap of the recipe undone, seen by the server (PR 5b)
    assert [n.text for n in session.log_manager.notes if n.kind == "note"] == ["Looked at it."]


def test_the_steps_provided_are_followed(tmp_path, launched):
    truth = MockTruth(components=[MockComponent("CIV", 1.9, 13.8, 10.0)],
                      lam_min_nm=440.0, lam_max_nm=460.0, dv_kms=2.5,
                      resol=45000.0, snr=50.0, z_em=2.0, seed=1)
    # In physical units, without a continuum (as physical_archive of
    # test_mcp_server.py)
    mock = generate(truth, truth_columns=False)
    d, k = mock._data, 1e-16 * (mock.x.value / 450.0) ** -1.5
    flux = au.erg / au.cm**2 / au.s / au.nm
    spec = SpectrumV2(SpectrumDataV2(
        x=d.x, xmin=d.xmin, xmax=d.xmax, y=DataColumnV2(mock.y.value * k, flux),
        dy=DataColumnV2(mock.dy.value * k, flux), aux_cols=dict(d.aux_cols), z_em=d.z_em))
    bare = tmp_path / "bare.acs2"
    SessionV2("bare", None, spec=spec).save(str(bare))
    _, opened, fitted = calls(
        build_server(), ("session_open", {"path": str(bare)}),
        ("gui_open", {"session_id": "bare"}),
        ("absorbers_doublets_auto", {"session_id": "bare", "provide_missing": True}))
    assert ok(fitted)["provided"] == ["continuum.estimate_auto"]
    names = [e.recipe_name for e in _state_of(ok(opened)["path"]).log_manager.entries]
    assert names == ["estimate_auto", "doublets_auto"]


def test_a_new_session_of_a_branching_recipe_is_not_followed(mock_archive, launched):
    server = build_server()
    _, opened, split = calls(
        server, open_mock(mock_archive), ("gui_open", {"session_id": "mock"}),
        ("edit_split", {"session_id": "mock", "expression": "x > 450"}))
    (new,) = ok(split)["new_sessions"]
    path = ok(opened)["path"]
    assert [p.name for p in Path(path).parent.glob("*.acs2")] == ["mock.acs2"]
    assert new["session_id"] == "mock-1"


def test_close_stops_following(mock_archive, launched):
    server = build_server()
    _, opened, closed = calls(
        server, open_mock(mock_archive), ("gui_open", {"session_id": "mock"}),
        ("session_close", {"session_id": "mock", "discard_changes": True}))
    ok(closed)
    path = ok(opened)["path"]
    assert not Path(path).exists() and follow.read_marker(path) is None


def test_gui_open_again_keeps_the_same_archive(mock_archive, launched):
    server = build_server()
    _, first, second = calls(server, open_mock(mock_archive),
                             ("gui_open", {"session_id": "mock"}),
                             ("gui_open", {"session_id": "mock"}))
    assert ok(first)["path"] == ok(second)["path"]
    assert len(launched) == 2


def test_gui_open_where_there_is_no_display(mock_archive, monkeypatch):
    monkeypatch.setattr(follow, "gui_unavailable_reason", lambda: "No display here.")
    monkeypatch.setattr(follow, "launch_gui", lambda path: pytest.fail("launched"))
    _, opened = calls(build_server(), open_mock(mock_archive),
                      ("gui_open", {"session_id": "mock"}))
    assert opened.is_error
    assert "The GUI is not available here: No display here." in opened.content[0].text


def test_gui_open_of_an_unknown_session(launched):
    (opened,) = calls(build_server(), ("gui_open", {"session_id": "nope"}))
    assert opened.is_error and "No session 'nope'" in opened.content[0].text
    assert launched == []


def test_a_failed_save_of_the_archive_is_a_warning(mock_archive, launched, monkeypatch):
    server = build_server()
    calls(server, open_mock(mock_archive), ("gui_open", {"session_id": "mock"}))

    def fail(history, path):
        raise OSError("disk full")
    monkeypatch.setattr(follow, "write_archive", fail)
    (estimated,) = calls(server, ("continuum_estimate_auto", {"session_id": "mock"}))
    info = ok(estimated)
    assert any("not updated (disk full)" in w for w in info["warnings"])
    assert info["spectrum"]["has_continuum"]  # the recipe is kept


def test_the_folder_is_removed_with_the_store(mock_archive):
    store = SessionStore()
    session_id = store.add(load_session_from_file(str(mock_archive), "mock", None, "auto"),
                           str(mock_archive))
    folder = Path(store.follow(session_id)).parent
    assert folder.is_dir()
    del store
    gc.collect()
    assert not folder.exists()
