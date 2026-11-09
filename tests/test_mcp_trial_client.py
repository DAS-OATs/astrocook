"""The minimal MCP client of the trials with fresh agents (phase 6-ter, PR 7).

``benchmarks/mcp_trials/client.py`` keeps a server alive over stdio and
answers ``tools`` and ``call`` on a Unix socket; the trials depend on it
showing the agent what a client shows a model, and on its log.
"""
import json
import subprocess
import sys
import time
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
CLIENT = REPO / "benchmarks" / "mcp_trials" / "client.py"
SERVER = Path(sys.executable).with_name("astrocook-mcp")

pytestmark = pytest.mark.skipif(sys.platform == "win32", reason="Unix socket")


def _client(*args):
    return subprocess.run([sys.executable, str(CLIENT), *args], check=True,
                          capture_output=True, text=True, timeout=60).stdout


def test_client_lists_the_tools_and_calls_one(tmp_path):
    sys.path.insert(0, str(REPO / "tests"))
    from mcp_stdio_check import make_mock
    archive = make_mock(tmp_path)
    # A short path: the length of a Unix socket path is limited
    sock = Path("/tmp") / f"ac-trial-{tmp_path.name[-20:]}.sock"
    log = tmp_path / "trial.jsonl"
    env_cmd = [sys.executable, str(CLIENT), "serve", "--socket", str(sock),
               "--log", str(log), "--images", str(tmp_path / "img"), "--",
               str(SERVER)]
    serve = subprocess.Popen(env_cmd, stdout=subprocess.DEVNULL,
                             stderr=subprocess.DEVNULL)
    try:
        deadline = time.time() + 60
        while not sock.exists():
            assert serve.poll() is None, "the client exited"
            assert time.time() < deadline, "the client did not start"
            time.sleep(0.2)
        listing = _client("tools", "--socket", str(sock))
        assert listing.startswith("# Instructions of the server")
        assert '"name": "session_open"' in listing
        opened = _client("call", "--socket", str(sock), "session_open",
                         json.dumps({"path": str(archive)}))
        session_id = json.loads(opened)["session_id"]
        _client("call", "--socket", str(sock), "continuum_estimate_auto",
                json.dumps({"session_id": session_id}))
        failed = _client("call", "--socket", str(sock), "session_describe",
                         json.dumps({"session_id": "nope"}))
        assert failed.startswith("ERROR")
    finally:
        serve.terminate()
        serve.wait(timeout=30)
        sock.unlink(missing_ok=True)
    events = [json.loads(line) for line in log.read_text().splitlines()]
    assert events[0]["event"] == "tools" and events[0]["n"] > 50
    calls = [e for e in events if e["event"] == "call"]
    assert [c["tool"] for c in calls] == ["session_open", "continuum_estimate_auto",
                                         "session_describe"]
    assert [c["failed"] for c in calls] == [False, False, True]
    # The progress of a recipe reaches the log, not the agent
    assert any(e["event"] == "progress" and e["tool"] == "continuum_estimate_auto"
               for e in events)
