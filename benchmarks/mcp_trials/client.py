"""A minimal MCP client for the trials of the server with fresh agents.

Phase 6-ter, PR 7 (``ROADMAP.md``, section 4.3). The agent of a trial sees
only what a client shows a model: the instructions of the server, the name,
description and input schema of each tool, and the content of each result
(the text, and the images, saved as PNG files that the agent opens). The
progress notifications are not shown (a client does not give them to the
model) but are written to the log of the trial, with the time of each, to
measure how long a call stays silent.

A server process lives as long as the trial, so that its sessions do too:
``serve`` starts it (the command after ``--``: ``astrocook-mcp``, or
``docker run -i --rm -v <folder>:/data astrocook-mcp``) and listens on a
Unix socket; ``tools`` and ``call`` talk to it::

    python client.py serve --socket s --log trial.jsonl --images img -- astrocook-mcp
    python client.py tools --socket s
    python client.py call --socket s session_open '{"path": "mock.acs2"}'

The log has one JSON line per event: the list of the tools, each call (its
arguments, the time it took, whether it failed, the size of its text) and
each progress notification.
"""
import argparse
import asyncio
import base64
import json
import socket
import sys
import time
from pathlib import Path

from mcp import Client, StdioServerParameters


def _log(path: Path, **event) -> None:
    event["t"] = round(time.time(), 3)
    with path.open("a") as f:
        f.write(json.dumps(event) + "\n")


def _format_tools(instructions: str, tools) -> str:
    """The instructions and the tools as a client gives them to a model."""
    out = ["# Instructions of the server", "", instructions or "(none)", "",
           f"# Tools ({len(tools)})", ""]
    for t in tools:
        out.append(json.dumps({"name": t.name, "description": t.description,
                               "input_schema": t.input_schema}))
        out.append("")
    return "\n".join(out)


def _format_result(result, images: Path, counter: list) -> str:
    parts = []
    for c in result.content:
        if c.type == "text":
            parts.append(c.text)
        elif c.type == "image":
            counter[0] += 1
            path = images / f"image_{counter[0]:03d}.png"
            path.write_bytes(base64.b64decode(c.data))
            parts.append(f"[image: {path.resolve()}]")
        else:
            parts.append(f"[{c.type} content]")
    text = "\n".join(parts)
    return ("ERROR (the tool failed):\n" + text) if result.is_error else text


async def serve(args) -> None:
    log, images = Path(args.log), Path(args.images)
    images.mkdir(parents=True, exist_ok=True)
    counter = [0]
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    params = StdioServerParameters(command=command[0], args=command[1:])
    async with Client(params, read_timeout_seconds=args.timeout) as client:
        tools = (await client.list_tools()).tools
        listing = _format_tools(client.instructions, tools)
        _log(log, event="tools", n=len(tools), chars=len(listing))
        lock = asyncio.Lock()

        async def handle(reader, writer):
            request = json.loads((await reader.readline()).decode())
            if request["op"] == "tools":
                reply = listing
                _log(log, event="list_tools")
            else:
                name, arguments = request["name"], request["arguments"]
                start = time.time()
                last = [start]

                async def progress(value, total, message):
                    now = time.time()
                    _log(log, event="progress", tool=name, value=value,
                         total=total, message=message,
                         since_last=round(now - last[0], 2))
                    last[0] = now

                async with lock:
                    try:
                        result = await client.call_tool(
                            name, arguments, progress_callback=progress)
                        reply = _format_result(result, images, counter)
                        failed = bool(result.is_error)
                    except Exception as exc:  # the protocol, not the tool
                        reply = f"ERROR (the call failed): {exc}"
                        failed = True
                end = time.time()
                _log(log, event="call", tool=name, arguments=arguments,
                     seconds=round(end - start, 2), failed=failed,
                     chars=len(reply),
                     silent_at_end=round(end - last[0], 2))
            writer.write(reply.encode())
            await writer.drain()
            writer.close()

        server = await asyncio.start_unix_server(handle, path=args.socket)
        print(f"serving {len(tools)} tools on {args.socket}", flush=True)
        async with server:
            await server.serve_forever()


def request(sock_path: str, payload: dict) -> str:
    with socket.socket(socket.AF_UNIX, socket.SOCK_STREAM) as s:
        s.connect(sock_path)
        s.sendall((json.dumps(payload) + "\n").encode())
        s.shutdown(socket.SHUT_WR)
        chunks = []
        while chunk := s.recv(65536):
            chunks.append(chunk)
    return b"".join(chunks).decode()


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="op", required=True)
    p = sub.add_parser("serve")
    p.add_argument("--socket", required=True)
    p.add_argument("--log", required=True)
    p.add_argument("--images", required=True)
    p.add_argument("--timeout", type=float, default=3600.0)
    p.add_argument("command", nargs=argparse.REMAINDER)
    p = sub.add_parser("tools")
    p.add_argument("--socket", required=True)
    p = sub.add_parser("call")
    p.add_argument("--socket", required=True)
    p.add_argument("name")
    p.add_argument("arguments", nargs="?", default="{}")
    args = parser.parse_args(argv)
    if args.op == "serve":
        asyncio.run(serve(args))
    elif args.op == "tools":
        print(request(args.socket, {"op": "tools"}))
    else:
        try:
            arguments = json.loads(args.arguments)
        except json.JSONDecodeError as exc:
            sys.exit(f"the arguments are not valid JSON: {exc}")
        print(request(args.socket, {"op": "call", "name": args.name,
                                    "arguments": arguments}))


if __name__ == "__main__":
    main()
