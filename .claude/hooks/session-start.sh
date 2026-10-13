#!/bin/bash
# SessionStart hook for Claude Code cloud sessions: installs the system
# libraries needed by Qt (offscreen) and the package with dev extras
# (test dependencies and the Ruff linter).
set -euo pipefail

if [ "${CLAUDE_CODE_REMOTE:-}" != "true" ]; then
  exit 0
fi

cd "$CLAUDE_PROJECT_DIR"

# Qt runtime libraries (PySide6 fails to import without libEGL).
QT_LIBS="libegl1 libgl1 libxkbcommon0 libfontconfig1 libdbus-1-3"
if ! dpkg -s $QT_LIBS >/dev/null 2>&1; then
  apt-get update -qq || true
  DEBIAN_FRONTEND=noninteractive apt-get install -y -qq $QT_LIBS >/dev/null
fi

# Virtual environment with the package in editable mode.
if [ ! -x .venv/bin/python ]; then
  python3 -m venv .venv
fi
.venv/bin/pip install -q --upgrade pip
.venv/bin/pip install -q -e ".[dev]"

{
  echo "export VIRTUAL_ENV=\"$CLAUDE_PROJECT_DIR/.venv\""
  echo "export PATH=\"$CLAUDE_PROJECT_DIR/.venv/bin:\$PATH\""
  echo "export QT_QPA_PLATFORM=offscreen"
  echo "export MPLBACKEND=Agg"
} >> "$CLAUDE_ENV_FILE"
