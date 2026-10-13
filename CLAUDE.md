# CLAUDE.md

Context for Claude Code sessions on Astrocook. The full plan (analysis,
phases, scientific assessment, v3 vision) is in `ROADMAP.md`: read it before
starting any non-trivial work.

## Language
Everything committed to the repository (code, comments, docs, commit
messages, PR descriptions) is in English. Conversations with Guido are in
Italian.

## Project
Astrocook: a Python package for the analysis of quasar spectra (continuum,
identification of absorption systems, Voigt fitting). PySide6 GUI +
scriptable API. Author and maintainer: Guido Cupani (INAF–OATs).
Current version: 2.0.1. Python ≥3.12.

## Architecture (v2)
- `astrocook/core/` – immutable dataclasses (`structures.py`), `SpectrumV2`,
  pure functions in `spectrum_operations.py`, session, system list.
  Methods return new objects: do not mutate state in place.
- `astrocook/fitting/` – Voigt fitting (`voigt_fitter.py`), constraints
  (`voigt_model.py`), Bayesian engine (`bayesian_fitter.py`).
- `astrocook/recipes/` – user-facing operations, with structured parameters.
- `astrocook/gui/` – PySide6 interface.
- `astrocook/io/` – loaders and V1 → V2 migration.
- `astrocook/legacy/` – V1 code, still imported by several V2 modules.
  To be removed progressively (phase 3); do not use it for new code.

## How to work
- Launch the GUI: `python run_astrocook.py [file ...]`
  (headless: `QT_QPA_PLATFORM=offscreen`).
- Tests: `pytest` (config in `pyproject.toml`). Coverage is still minimal:
  before changing production code, add characterization tests for the
  behaviour involved.
- Work in small PRs; the test suite must be green before and after.
- Do not change numerical behaviour without Guido's agreement: if a
  behaviour looks like a bug, report it instead of silently fixing it.
- User documentation (Sphinx + MyST) lives in `docs/`.

## Direction
- v2: consolidation (tests and CI, cleanup, decoupling from legacy,
  refactoring, error handling).
- Astrocook as a tool for agents: MCP server over the recipes, API with
  structured outputs, reliable natural-language log.
- v3 (horizon): a collection of complete models with alternative
  identifications and credibility, generated automatically in under an hour.

## Current status
Analysis and roadmap done (October 2026). Next step: phase 1 (safety net),
as a pilot. Open GitHub issues are scheduled in `ROADMAP.md`, section 3.1.
Update this section and the table in `ROADMAP.md` at the end of
each phase.
