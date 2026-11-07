"""The GUI launched by the MCP server (phase 6-ter, PR 5-bis; decision 7 of
the phase, option b): ``gui_open`` saves a session to an archive and launches
the GUI on it; the server saves the archive again after every change, and the
GUI follows it.

The first test is a characterization test written before ``gui_open``:
``session_save`` writes the archive only.
"""
import asyncio

import pytest

pytest.importorskip("mcp")

from mcp import Client  # noqa: E402

from astrocook.core.session import SessionV2  # noqa: E402
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
