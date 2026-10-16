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
Current version: 2.0.2. Python ≥3.12.

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
  `pytest --cov=astrocook`. Coverage is still low (56% of the V2 code):
  before changing production code, add characterization tests for the
  behaviour involved.
- Numerical snapshots: `tests/snapshots/*.npz`, compared by the `snapshot`
  fixture (`tests/conftest.py`). Regenerate only for an intended change
  (`pytest --update-snapshots`) and say so in the PR.
- Work in small PRs; the test suite must be green before and after.
- Do not change numerical behaviour without Guido's agreement: if a
  behaviour looks like a bug, report it instead of silently fixing it.
- User documentation (Sphinx + MyST) lives in `docs/`.
- Heuristics (thresholds, penalties, bounds, magic numbers) often encode a
  purpose that no test records (e.g. the penalty on b > 50 km/s kept the
  forest fits fast, F32). Before removing or changing one, ask Guido what it
  is for; when the purpose is known, write it as a comment next to the code.
- For PRs that change fitting or the automatic pipelines, run the slow tests
  (`pytest --run-slow -m slow`, here or from the Actions page on the branch)
  before asking for the merge, and report the result in the PR.

## Direction
- v2: consolidation (tests and CI, cleanup, decoupling from legacy,
  refactoring, error handling).
- Astrocook as a tool for agents: MCP server over the recipes, API with
  structured outputs, reliable natural-language log.
- v3 (horizon): a collection of complete models with alternative
  identifications and credibility, generated automatically in under an hour.

## Current status
Phases 0, 1 and 1-bis done (October 2026; details per item in `ROADMAP.md`).
The safety net: numerical snapshots (pure functions, continuum,
identification, Voigt model and fit, SpectrumV2 methods, full chain on
J0103-1305 and on a mock, with intent checks), V1/V2 parity tests, mock
generator, recipe-level and GUI smoke tests, `.acs2` round trips,
error-scaling check, mutation testing (`mutmut`, config in
`pyproject.toml`, run by hand), timing benchmark
(`benchmarks/bench_fitting.py`, run by hand). Slow tests
(`pytest --run-slow`) run in a separate CI job. Still waiting for data from
Guido: sessions for log replay, V1 `.acs` archives, sample files for the
loaders, X-Shooter VIS/NIR data (F15).

Phase 1-bis (findings register, `ROADMAP.md` section 3.2) fixed F1, F2, F5,
F10 (minimal), F11, F16, F19, F20, F23, F24, F25, F26, F27, F29, F32, F33
and issues #480, #442; F3 accepted. Deferred to GitHub issues, to be
scheduled with Guido: #526 (resolution as a single column, exact
convolution: F8, F22, rest of F10, F9, F14), #465 (continuum: spikes, F30,
F17), #517
(components at the bounds, F21), #518 (editable z window, F4). F15 waits
for data.

Release 2.0.2 (October 2026) includes phase 1-bis and #534 (errors scaled
with the chi2 of the final model). Next: to be decided with Guido (the
deferred issues or phase 2). Work in small PRs, one at a time; wait for the merge before the
next.

Questionable behaviour found along the way goes in the findings register
(`ROADMAP.md`, section 3.2), reviewed with Guido;
anything that loses user data or crashes is reported at once. Open GitHub
issues are scheduled in `ROADMAP.md`, section 3.1. Wait for each PR to be
merged before opening the next. Update this section and the table in
`ROADMAP.md` at the end of each phase.
