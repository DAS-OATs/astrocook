# CLAUDE.md

Context for Claude Code sessions on Astrocook. The full plan (analysis,
phases, scientific assessment, v3 vision) is in `ROADMAP.md`: read it before
starting any non-trivial work.

## Language
Everything committed to the repository (code, comments, docs, commit
messages, PR descriptions) is in English. Conversations with Guido are in
Italian.

## Working with Guido
- Do not subscribe to PR activity (no automatic watching of the PRs).
- Short progress messages during the work are fine; at the end of each
  task, a discursive summary of what was done, why, what was found and
  what is left (not just a terse list). Do not wait for long test runs to
  finish before sending it: report their result when it comes.

## Project
Astrocook: a Python package for the analysis of quasar spectra (continuum,
identification of absorption systems, Voigt fitting). PySide6 GUI +
scriptable API. Author and maintainer: Guido Cupani (INAF–OATs).
Current version: 2.1.0. Python ≥3.12.

## Architecture (v2)
- `astrocook/core/` – immutable dataclasses (`structures.py`), `SpectrumV2`,
  pure functions in `spectrum_operations.py`, session, system list.
  Methods return new objects: do not mutate state in place.
- `astrocook/fitting/` – Voigt fitting (`voigt_fitter.py`), constraints
  (`voigt_model.py`), Bayesian engine (`bayesian_fitter.py`).
- `astrocook/recipes/` – user-facing operations, with structured parameters.
- `astrocook/gui/` – PySide6 interface.
- `astrocook/io/` – loaders and V1 → V2 migration.
- `astrocook/mcp_server/` – MCP server for agents (pilot, `ROADMAP.md`
  section 4.1): a thin layer over the recipes, optional `[mcp]` extra.
- `astrocook/mocks/` – synthetic spectra with known parameters, independent
  of the fitting code (do not make it depend on `fitting/`).
- V1 code: removed from the package in phase 3. V1 `.acs` archives are read
  by `astrocook/io/v1_archive.py`; the V1 routines used by the parity tests
  are verbatim copies in `tests/reference/v1/` (test-only, do not edit).

## How to work
- Environment: in cloud sessions `.claude/hooks/session-start.sh` installs
  the Qt system libraries and the package (with `[dev]` extras) into
  `.venv/`, and sets `QT_QPA_PLATFORM=offscreen` and `MPLBACKEND=Agg`.
  Locally: `python -m venv .venv && .venv/bin/pip install -e ".[dev]"`.
- Lint: `ruff check .` (config in `pyproject.toml`; pyflakes rules only for
  now, `tests/reference/v1/` excluded). New code must not add warnings.
- CI (GitHub Actions, `.github/workflows/tests.yml`): `pytest` on Python
  3.12 and 3.13 (coverage report on 3.13), and `ruff check .`. The slow
  tests (`pytest --run-slow -m slow`) run in a separate job, not on pull
  requests: after each merge to master, weekly (Mondays) and on demand
  from the Actions page.
- Launch the GUI: `astrocook [file ...]`, `python -m astrocook` or
  `python run_astrocook.py` (launcher in `astrocook/_run_astrocook.py`)
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
- User documentation (Sphinx + MyST) lives in `docs/`, published by
  `.github/workflows/docs.yml` (`build_site.sh`) on `gh-pages`, one version
  per release: the root holds the docs of the latest release, rebuilt from
  its tag when the release is published; `/dev/` holds the docs of master,
  rebuilt after each merge; `/v1/` the V1 docs. Before publishing a
  release, update the "preferred" entry of `switcher.json` on master (the
  release build checks it). `deploy_docs.sh stable|dev` does the same by
  hand.
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
Phases 0, 1, 1-bis, 2, 2-bis, 3 and 3-bis done (October 2026; details per item in `ROADMAP.md`).
The safety net: numerical snapshots (pure functions, continuum,
identification, Voigt model and fit, SpectrumV2 methods, full chain on
J0103-1305 and on a mock, with intent checks), V1/V2 parity tests, mock
generator, recipe-level and GUI smoke tests, `.acs2` round trips,
error-scaling check, mutation testing (`mutmut`, config in
`pyproject.toml`, run by hand), timing benchmark
(`benchmarks/bench_fitting.py`, run by hand). Slow tests
(`pytest --run-slow`) run in a separate CI job. Still waiting for data from
Guido: sessions for log replay, sample files for the loaders, X-Shooter
VIS/NIR data (F15). V1 `.acs` archives: received in phase 3 (two sessions,
cut down and anonymised in `tests/data/v1_*.acs`).

Phase 1-bis (findings register, `ROADMAP.md` section 3.2) fixed F1, F2, F5,
F10 (minimal), F11, F16, F19, F20, F23, F24, F25, F26, F27, F29, F32, F33
and issues #480, #442; F3 accepted. Deferred to GitHub issues, scheduled
for phase 4: #526 (resolution as a single column, exact
convolution: F8, F22, rest of F10, F9, F14), #465 (continuum: spikes, F30,
F17), #517
(components at the bounds, F21), #518 (editable z window, F4). F15 waits
for data.

Release 2.0.2 (October 2026) includes phase 1-bis and #534 (errors scaled
with the chi2 of the final model). Sequence agreed with Guido (October
2026, revised after phase 3; `ROADMAP.md`, section 3):
1. phase 2 – cleanup and packaging (done, release 2.0.3);
2. phase 3 – decoupling from `legacy/` (done);
3. phase 3-bis – small fixes before the release (done): #575 (F41),
   F6, F13, #559 (F36, brought forward from phase 4), F45 (version in the
   `.acs2` archive), #560 (`norm` of `rebin` no longer offered);
4. release 2.1.0, on GitHub with the macOS app, not on PyPI (done);
5. MCP pilot (3–4 recipes, separate module, installed from GitHub;
   `ROADMAP.md`, section 4);
6. phase 4 – scientific correctness: first a validation on mocks
   (baseline of continuum, identification, doublet fit), then #526,
   #517, #465, #518, #568, #573, each measured on it;
7. release 2.2.0;
8. phase 5 – splitting large files;
9. phase 6 – error handling;
10. phase 6-bis – recipe API (typed parameters, structured results,
    exceptions), from the report of the MCP pilot;
11. full MCP server, with a headless Docker image, its tools generated
    from the recipe registry;
12. phase 7 – robust components: emission model and
    absorption identification, algorithmic choices by Guido, measured
    on the benchmark of phase 4.

LICENSE (#549) and the first publication on PyPI are an independent step:
they go into the first release after INAF answers (2.1.1 if before phase
4, otherwise 2.2.0 or the current version then).

Phase 2 done, released as 2.0.3 (packaging, no numerical changes).
LICENSE (GPLv3 or later) waits for INAF: PR #549 stays a draft. Docker:
no renewal for now; the V1 image (`containers/`, broken) is removed in
phase 3.

Phase 2-bis done (October 2026): F34 (window icon), Ruff warnings (full
`ruff check .` in CI), Sphinx warnings, optional `[bayes]` extra (emcee,
dynesty, corner); found F36 and F37 (#559 and #560, fixed in phase 3-bis).

Phase 3 done (October 2026, PRs #565-#584): `astrocook/legacy/` removed;
V1 `.acs` archives read by `io/v1_archive.py` (saving in `.acs` removed);
the V1 routines of the parity tests in `tests/reference/v1/`; logging left
to the application; the data files used only by V1 removed (telluric data
kept as assets of the pre-release `data-telluric-2026`); docs published per
version. Fixed F42-F44; found F38-F41 (#568 and #573 in phase 4, #575 in
phase 3-bis); #566 (centralized messages) in phase 6.

Phase 3-bis done (October 2026, PRs #587-#592): branching recipes logged
only in the new session (#575, F41); the core no longer imports the GUI
(F6, `core/debug_hooks.py`); link colours of the System Inspector
recomputed and stable (F13); `dbtur` stored after a fit (#559, F36: 0 when
`btur` is frozen, the default); `.acs2` archives record the Astrocook
version (F45, `SessionV2.saved_with_version`); `norm` of `rebin` no longer
offered (#560, F37, Guido's decision). Note on F8 for #526 (`rebin` keeps
the old `resol` meta key).

Release 2.1.0 (October 2026): phases 3 and 3-bis, on GitHub with the
macOS app, not on PyPI (steps of a release in `ROADMAP.md`, "Releases").

MCP pilot in progress (`ROADMAP.md`, section 4.1: scope, tools, PR order;
`astrocook/mcp_server/`, `[mcp]` extra, `astrocook-mcp`). Work in
small PRs, one at a time; wait for the merge before the next. While #549 is open, each new PR goes on its own
`claude/...` branch (agreed with Guido; say which branch).

Questionable behaviour found along the way goes in the findings register
(`ROADMAP.md`, section 3.2), reviewed with Guido;
anything that loses user data or crashes is reported at once. Open GitHub
issues are scheduled in `ROADMAP.md`, section 3.1. Wait for each PR to be
merged before opening the next. Update this section and the table in
`ROADMAP.md` at the end of each phase.
