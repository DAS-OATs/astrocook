# Trials of the MCP server with fresh agents

Phase 6-ter, PR 7 (`ROADMAP.md`, section 4.3; the pilot's trials in
`ROADMAP_HISTORY.md`, "4.2 MCP pilot: report"). An agent with no knowledge of
the code, the tests or the true parameters uses the server through a minimal
MCP client, seeing only what a client shows a model; it gets a request in
plain words, answers it, then writes a critique of the server.

- `client.py`: the client. `serve` starts the server (the command after
  `--`) and keeps it alive on a Unix socket; `tools` prints the instructions
  and the tools as a client gives them to a model; `call` runs a tool and
  prints its content (images saved as PNG files). Its log (`--log`, one JSON
  line per event) holds each call with its time, and the progress
  notifications with the time since the one before, which the agent does not
  see.
- `make_data.py`: the input of each trial (`mock`, `j0103`, `edit`). The
  agent must not read it: it holds the true parameters of the mocks.
- `prompts.md`: the rules given to every agent and the request of each trial.

## Running a trial

```bash
T=/some/folder/trial1; mkdir -p $T/data $T/gaps
python benchmarks/mcp_trials/make_data.py mock $T/data
cp benchmarks/mcp_trials/client.py /elsewhere/client.py   # out of the agent's folder
```

Then the server, in the Docker image (`docker build -t astrocook-mcp .`) or
in a virtual environment, with the file of the gaps out of the user's file:

```bash
python /elsewhere/client.py serve --socket $T/server.sock --log trial1.jsonl \
    --images $T/images -- docker run -i --rm -v $T/data:/data -v $T/gaps:/gaps \
    -e ASTROCOOK_GAPS_FILE=/gaps/gaps.jsonl astrocook-mcp
```

and a wrapper `$T/mcp` (`python /elsewhere/client.py "$op" --socket
$T/server.sock "$@"`) that the agent calls as `./mcp tools` and `./mcp call
<tool> '<json>'`. The agent gets the rules and the request of `prompts.md`;
afterwards its transcript is checked for reads of the code (none allowed).

Behind a proxy that intercepts TLS (as in the cloud sessions of Claude Code)
the image builds only from a local base with the proxy's certificate, e.g. a
`python:3.13-slim` tagged locally from `FROM python:3.13-slim` + `COPY
ca.crt` + `ENV PIP_CERT=...`, with `docker build --network host --build-arg
HTTPS_PROXY=...`; the `Dockerfile` itself needs no change.
