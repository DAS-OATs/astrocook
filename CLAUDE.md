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
- `astrocook/mocks/` – synthetic spectra with known parameters, independent
  of the fitting code (do not make it depend on `fitting/`).
- `astrocook/legacy/` – V1 code, still imported by several V2 modules.
  To be removed progressively (phase 3); do not use it for new code.

## How to work
- Environment: in cloud sessions `.claude/hooks/session-start.sh` installs
  the Qt system libraries and the package (with `[dev]` extras) into
  `.venv/`, and sets `QT_QPA_PLATFORM=offscreen` and `MPLBACKEND=Agg`.
  Locally: `python -m venv .venv && .venv/bin/pip install -e ".[dev]"`.
- Lint: `ruff check .` (config in `pyproject.toml`; pyflakes rules only for
  now, `legacy/` excluded). New code must not add warnings.
- CI (GitHub Actions, `.github/workflows/tests.yml`): `pytest` on Python
  3.12 and 3.13 (coverage report on 3.13), and `ruff check --select F821`
  (undefined names only). The slow tests (`pytest --run-slow -m slow`) run
  in a separate job, not on pull requests: after each merge to master,
  weekly (Mondays) and on demand from the Actions page.
- Launch the GUI: `python run_astrocook.py [file ...]`
  (headless: `QT_QPA_PLATFORM=offscreen`).
- Tests: `pytest` (config in `pyproject.toml`); coverage with
  `pytest --cov=astrocook`. Coverage is still low (42% of the V2 code):
  before changing production code, add characterization tests for the
  behaviour involved.
- Numerical snapshots: `tests/snapshots/*.npz`, compared by the `snapshot`
  fixture (`tests/conftest.py`). Regenerate only for an intended change
  (`pytest --update-snapshots`) and say so in the PR.
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
Phase 0 done. Phase 1 (safety net) in progress; progress per item in
`ROADMAP.md`, phase 1. Done: snapshots (pure functions, continuum,
identification, Voigt model and fit, SpectrumV2 methods, full chain on
J0103-1305 and on a mock), tests for the "possibly fixed" issues, V1/V2
parity tests, CI, mock generator (`astrocook/mocks/`; sky residuals and
tellurics still missing). Slow tests (the J0103-1305 full chain, ~15
minutes) are marked `slow` and run only with `pytest --run-slow` (and in
a separate CI job, see above).
Coverage measured (V2 without legacy: 38%, 42% after the first targeted
tests; per area and main gaps in `ROADMAP.md`, phase 1, item 6).
Targeted recipe-level tests in progress: `edit.py` and `flux.py` done,
`absorbers.py` and system list operations next.
Next, in this order: GUI smoke tests;
`.acs2` round trips; chi-squared error-scaling check; mutmut. Log replay
and V1 `.acs` archives wait for sessions from Guido.

Questionable behaviour found along the way goes in the findings register
(`ROADMAP.md`, section 3.2), reviewed with Guido at the end of phase 1;
anything that loses user data or crashes is reported at once. Open GitHub
issues are scheduled in `ROADMAP.md`, section 3.1. Wait for each PR to be
merged before opening the next. Update this section and the table in
`ROADMAP.md` at the end of each phase.
