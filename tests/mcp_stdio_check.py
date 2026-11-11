"""A client that checks an MCP server of Astrocook over stdio, from outside.

Not a test of the suite (no ``test_`` in the name): a script that starts the
server as a client does (the command after ``--``), and runs through the
tools that the headless image has to offer (phase 6-ter, PR 6): open a mock
spectrum, run a recipe, save the session, plot it, and ask for the GUI,
which must answer that it is not available.

The CI job ``docker`` runs it inside the image (the script is mounted, the
package and the ``mcp`` client are in the image)::

    docker run --rm -v "$PWD/tests/mcp_stdio_check.py:/check.py:ro" \\
        --entrypoint python astrocook-mcp /check.py -- astrocook-mcp

and ``tests/test_mcp_headless.py`` runs it with this interpreter, with
PySide6 made unavailable. With no command it starts ``astrocook-mcp``. The
exit status is 0 if every step answers as expected, else 1 with the reason.
"""
import argparse
import asyncio
import base64
import sys
import tempfile
from pathlib import Path

from mcp import Client, StdioServerParameters

from astrocook.session.session import SessionV2
from astrocook.mocks import MockComponent, MockTruth, generate


def make_mock(folder: Path) -> Path:
    """A small mock (one CIV doublet, z_em = 2.0, R = 45000) as ``.acs2``."""
    truth = MockTruth(components=[MockComponent("CIV", 1.9, 13.8, 10.0)],
                      lam_min_nm=440.0, lam_max_nm=460.0, dv_kms=2.5,
                      resol=45000.0, snr=50.0, z_em=2.0, seed=1)
    path = folder / "mock.acs2"
    SessionV2("mock", None, spec=generate(truth)).save(str(path))
    return path


def check(condition, message):
    if not condition:
        raise AssertionError(message)


def text(result):
    """The text of a result (the first text content: an image comes first)."""
    return next((c.text for c in result.content if c.type == "text"), "")


async def run(command, args, folder: Path):
    archive = make_mock(folder)
    saved = folder / "saved.acs2"
    async with Client(StdioServerParameters(command=command, args=args)) as client:
        tools = {t.name for t in (await client.list_tools()).tools}
        for name in ("session_open", "session_save", "continuum_estimate_auto",
                     "spectrum_plot", "gui_open"):
            check(name in tools, f"the tool {name} is not offered")
        print(f"ok: {len(tools)} tools")

        async def call(tool, **arguments):
            return await client.call_tool(tool, arguments)

        opened = await call("session_open", path=str(archive), z_em=2.0, resol=45000)
        check(not opened.is_error, f"session_open: {text(opened)}")
        session_id = opened.structured_content["session_id"]
        print(f"ok: session_open ({session_id})")

        ran = await call("continuum_estimate_auto", session_id=session_id,
                         reason="check of the image")
        check(not ran.is_error, f"continuum_estimate_auto: {text(ran)}")
        print("ok: continuum_estimate_auto")

        result = await call("session_save", session_id=session_id, path=str(saved))
        check(not result.is_error, f"session_save: {text(result)}")
        check(saved.is_file(), "session_save wrote no file")
        print("ok: session_save")

        plot = await call("spectrum_plot", session_id=session_id,
                          x_min=445, x_max=452)
        check(not plot.is_error, f"spectrum_plot: {text(plot)}")
        image = plot.content[0]
        check(image.type == "image" and image.mime_type == "image/png",
              "spectrum_plot gave no PNG")
        check(base64.b64decode(image.data)[:8] == b"\x89PNG\r\n\x1a\n",
              "the image is not a PNG")
        print("ok: spectrum_plot")

        gui = await call("gui_open", session_id=session_id)
        check(gui.is_error, "gui_open did not refuse")
        check("The GUI is not available" in gui.content[0].text,
              f"gui_open: {text(gui)}")
        print(f"ok: gui_open refused ({text(gui).split(' Save ')[0]})")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("server", nargs="*", default=["astrocook-mcp"],
                        help="the command of the server, after -- (default: astrocook-mcp)")
    server = (argv if argv is not None else sys.argv[1:])
    server = parser.parse_args(server).server or ["astrocook-mcp"]
    with tempfile.TemporaryDirectory() as tmp:
        try:
            asyncio.run(run(server[0], server[1:], Path(tmp)))
        except AssertionError as e:
            print(f"FAILED: {e}", file=sys.stderr)
            return 1
    print("all steps passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
