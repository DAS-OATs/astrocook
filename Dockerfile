# Headless image of the MCP server of Astrocook (phase 6-ter, PR 6): the
# package and the [mcp] extra without PySide6, so no Qt libraries and no GUI
# (the tool gui_open says that the GUI is not available). Not published: the
# users build it (docs: "Using Astrocook from an AI Agent").
#
#   docker build -t astrocook-mcp .
#   docker run -i --rm -v /path/to/data:/data astrocook-mcp
#
# The server speaks MCP over stdio, so the container is run with -i and no
# terminal; a client (Claude Code, Claude Desktop) starts it by itself.
FROM python:3.13-slim

ENV PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    MPLBACKEND=Agg \
    ASTROCOOK_GAPS_FILE=/data/.astrocook-gaps.jsonl

# The dependencies of the package, listed here and not read from
# pyproject.toml, so that PySide6 is left out without changing the
# dependencies of the package. tests/test_docker_image.py checks that this
# list is that of pyproject.toml (the dependencies and the [mcp] extra), but
# for what only the GUI needs (PySide6, qtawesome) or nothing does
# (scienceplots): do not edit one without the other. The block between the
# markers is read by that test.
# deps:begin
RUN pip install \
    "astropy>=6.1" \
    "numpy>=1.26" \
    "scipy>=1.14" \
    "matplotlib>=3.9" \
    "numexpr" \
    "mcp>=2.3,<3"
# deps:end

# Layers in order of change: the package last
WORKDIR /src
COPY pyproject.toml README.md MANIFEST.in ./
COPY astrocook ./astrocook
RUN pip install --no-deps . && rm -rf /src

# Not root. /data is where the files of the user are mounted (-v), and where
# the server writes the file of the gaps (decision 6 of the phase) so that it
# is not lost with the container (--rm)
RUN useradd --create-home astrocook && mkdir /data && chown astrocook /data
USER astrocook
WORKDIR /data
VOLUME /data

ENTRYPOINT ["astrocook-mcp"]
