# CLAUDE.md

Context for Claude Code sessions on Astrocook. The plan is in `ROADMAP.md`
(about 1000 lines, 13 November 2026): do not read it whole. Read the section
of the phase in progress (phase 7, once 2.7.0 is out), "Rules and decisions
in force", and the open rows of the findings register (section 3.2) and of
the issues (section 3.1).
`ROADMAP_HISTORY.md` is the archive of the concluded phases: search it with
`grep` for a detail (a finding, an issue, a PR), never read it whole.

## Language
Everything committed to the repository (code, comments, docs, commit
messages, PR descriptions) is in English. Conversations with Guido are in
Italian.

## Working with Guido
- Do not subscribe to PR activity (no automatic watching of the PRs).
- After opening a PR, do not check its CI and do not schedule reminders
  for it: Guido checks the CI himself (changed in November 2026). Say in
  the summary which checks were run locally.
- Short progress messages during the work are fine; at the end of each
  task, a discursive summary of what was done, why, what was found and
  what is left (not just a terse list). Do not wait for long test runs to
  finish before sending it: report their result when it comes.

## Choice of the model (agreed with Guido, November 2026)
Guido gives each PR (except the smallest) a conversation of its own, and
chooses its model. At the start of each phase, when its plan is written,
say for each PR whether Sonnet is enough, with an estimate of its time,
and write both in the plan of the phase in `ROADMAP.md`.
Calibration (phase 6-quater): PRs 1-7 took one working day in all, about
1 h for a mechanical PR (Sonnet) and 2-4 h for a review (Opus); the
estimates of the plan had been 4-8 times longer. Later in the phase a
mechanical PR estimated at 1 h took about 2 with tests, documentation and
screenshots (PRs 11, 12): estimate 2 h for those.
- Sonnet is enough for mechanical work with a criterion of success that
  a tool checks: lint and type-checker warnings, renamings with aliases,
  moving modules, wiring an existing piece of data to another place (e.g.
  the labels of the menus from the registry), short texts whose content
  is already decided. Tell it the rule "no change of behaviour" in the
  prompt (e.g. for the type checker: `assert` or `cast`, no new branches).
- Opus for the work of judgement: the design of a mechanism used by
  several PRs, the review of the interface of the recipes (names, types,
  defaults, units), anything that can change the numerical results or
  touches the science, the plans of the phases.

## Project
Astrocook: a Python package for the analysis of quasar spectra (continuum,
identification of absorption systems, Voigt fitting). PySide6 GUI +
scriptable API. Author and maintainer: Guido Cupani (INAF–OATs).
Current version: 2.6.0 (2.7.0 is PR 14 of phase 6-quater). Python ≥3.12.

## Architecture (v2)
- `astrocook/core/` – immutable dataclasses (`structures.py`), `Spectrum`,
  pure functions in `spectrum_operations.py`, system list. The bottom layer:
  it does not import `recipes/` or `session/`.
  Methods return new objects: do not mutate state in place. `Spectrum`
  (`spectrum.py`) holds the data, columns and conversions; its other
  methods are in mixins by theme (`spectrum_continuum.py`,
  `spectrum_lines.py`, `spectrum_flux.py`, `spectrum_combine.py`,
  `spectrum_expressions.py`, phase 5): a new method goes in the mixin of
  its theme, and builds its result with `self.__class__(...)`.
- `astrocook/session/` – `Session` (`session.py`) and `SessionHistory`
  (`history.py`), moved out of `core/` in phase 6-quater, PR 5. The old
  paths `astrocook.core.session` and `astrocook.core.session_manager` stay
  as lazy aliases (module `__getattr__`, no warning): new code imports from
  `astrocook.session`.
- Class names: no `V2` suffix since phase 6-quater, PR 6 (`Spectrum`,
  `Session`, `DataColumn`, `RecipeAbsorbers`, `MainWindow`, ...). The old
  names (`SpectrumV2`, ...) are served by a module `__getattr__` made by
  `astrocook/_renamed.py` (`renamed(__name__, {old: new})`): the same class,
  with a `DeprecationWarning`, for scripts only (archives and logs hold no
  class names). A class renamed later gets its entry in the module that
  defines it. The `"AstrocookV2"` of `QSettings` stays (the key of the
  user's saved settings).
- `astrocook/fitting/` – Voigt fitting (`voigt_fitter.py`), constraints
  (`voigt_model.py`), Bayesian engine (`bayesian_fitter.py`).
- `astrocook/recipes/` – user-facing operations, with structured parameters.
- `astrocook/gui/` – PySide6 interface. Split in phase 5, as the spectrum:
  `MainWindow` (`main_window.py`) inherits its mixins by theme
  (`main_menus.py`, `main_panels.py`, `continuum_controls.py`,
  `main_sessions.py`, `main_recipes.py`, `main_scripts.py`), the main plot
  `SpectrumPlotWidget` (`pyside_plot.py`) those of the continuum editor and
  of the mouse (`continuum_editor.py`, `plot_mouse.py`), and its
  `plot_spectrum` calls one method per layer; the System Inspector keeps its
  table and velocity plot in `inspector_table.py` and `velocity_plot.py`. A
  new method goes in the module of its theme; `tests/test_gui_structure.py`
  checks the interface of the split classes.
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
- Lint: `ruff check .` (config in `pyproject.toml`; the pyflakes rules
  and flake8-bugbear, `B`, since phase 6-quater; `tests/reference/v1/`
  excluded). A `zip()` says `strict=True` where its inputs have the same
  length by construction, `itertools.pairwise` for consecutive pairs. New
  code must not add warnings.
- Types: `mypy` (no arguments; config in `[tool.mypy]` of `pyproject.toml`,
  no strict option). It checks the whole package (`recipes/`,
  `mcp_server/`, `core/`, `fitting/` and `gui/`) with zero errors since
  phase 6-quater, PR 4-bis. New code must not add errors. A mixin of the
  spectrum inherits `HostBase` (`core/_spectrum_host.py`) and one of the
  main window `HostBase` of `gui/_window_host.py` (the plot mouse one,
  `gui/_plot_host.py`): what it uses of the class that mixes it in is
  declared there (types only, `object` at run time), and the methods of
  the spectrum that return a new one are annotated `-> Self`; a new
  attribute or method used across mixins is declared there too. To narrow
  an `Optional` use `spec_of(session)` or `column_of(spec, name)`
  (`recipes/_spectrum.py`), `assert` or `cast` (in `gui/` always `cast`:
  an event handler must raise what it raised), no new branches; a Qt enum
  is written with its class (`Qt.CursorShape.WaitCursor`); a
  `# type: ignore[code]` says why, and a real bug found goes in the
  findings register instead of being fixed.
- CI (GitHub Actions, `.github/workflows/tests.yml`): `pytest` on Python
  3.12 and 3.13 (coverage report on 3.13), `ruff check .` and `mypy`
  (the `types` job, on the full `[dev]` install). The slow
  tests (`pytest --run-slow -m slow`) run in a separate job, not on pull
  requests: after each merge to master, weekly (Mondays) and on demand
  from the Actions page.
- Launch the GUI: `astrocook [file ...]`, `python -m astrocook` or
  `python run_astrocook.py` (launcher in `astrocook/_run_astrocook.py`)
  (headless: `QT_QPA_PLATFORM=offscreen`).
- Tests: `pytest` (config in `pyproject.toml`), in parallel with
  `pytest -n auto` (pytest-xdist, ~3 minutes instead of ~7; the CI does
  so); coverage with `pytest --cov=astrocook` (86%, the GUI 85%, the core
  92%, after phase 6-bis): before changing production code, add
  characterization tests for the behaviour involved.
- Which tests to run locally (agreed with Guido, November 2026): the test
  files of the change and those of the code it touches (e.g. the recipes
  changed: their review, inputs, errors, registry and MCP schema tests),
  plus `ruff check .` and `mypy`; the full suite is left to the CI of the
  PR (in the cloud sessions it takes about 25 minutes on 4 cores, not 3).
  Run it in full only for a change across the package (a rename, a moved
  module, a shared mechanism such as `typed_recipe` or the registry). Say
  in the PR which tests were run.
- Warnings are errors in the tests (`filterwarnings = ["error"]` in
  `pyproject.toml`, phase 6): a test that expects a warning says so
  (`pytest.warns`, or `warnings.catch_warnings` with a local filter); a
  warning that a new release of a dependency brings is fixed in the code,
  or ignored by its own entry in `filterwarnings`, with its reason; in the
  package, a cosmetic warning of a library is silenced where it is raised,
  never with a global filter.
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
- Screenshots of the tutorials: `docs/screenshots/make_screenshots.py`
  rebuilds them from the GUI (README in the same folder); rerun the scenes
  touched by a change of the GUI, look at every image, commit those that
  changed.
- Document everything found and used (sources, references, measurements,
  decisions) when it happens, not afterwards: in `docs/` for users (e.g.
  `docs/atomic_data.md`), in `ROADMAP.md` and the findings register for the
  development.
- Heuristics (thresholds, penalties, bounds, magic numbers) often encode a
  purpose that no test records (e.g. the penalty on b > 50 km/s kept the
  forest fits fast, F32). Before removing or changing one, ask Guido what it
  is for; when the purpose is known, write it as a comment next to the code.
- Slow tests and benchmark (agreed with Guido, October 2026): before asking
  for the merge, run the slow tests (`pytest --run-slow -m slow`, here or
  from the Actions page on the branch) and, in phase 4, the full mock
  benchmark (`validate_mocks.py --compare`) only for PRs that can change
  the numerical results of fitting, identification or the automatic
  pipelines; report the results in the PR. For changes that by
  construction leave the results unchanged (e.g. what a pipeline stores
  besides its results, refactoring covered by the snapshots, docs), the
  fast suite is enough: say in the PR why the numbers cannot change, and
  leave the slow tests to the CI job after the merge. When the benchmark
  cannot measure a fix (e.g. F55: no absorption beyond z_em in the mocks),
  a targeted test replaces it; the same when no slow test and no pipeline
  goes through the code changed (check with `grep` first, e.g. PR 10-bis
  of 6-quater: nothing calls `fit_bayesian`): say so in the PR and skip
  them.

## Direction
- v2: consolidation (tests and CI, cleanup, decoupling from legacy,
  refactoring, error handling).
- Astrocook as a tool for agents: MCP server over the recipes, API with
  structured outputs, reliable natural-language log.
- v3 (horizon): a collection of complete models with alternative
  identifications and credibility, generated automatically in under an hour.
- Papers (`ROADMAP.md`, section 6; the full plan in the private
  `astrocook-data`, `papers/plan.md`): a reference paper (MNRAS) with the
  validation on mocks and the automatic re-analysis of HE0940–1050, after
  phase 7. The real spectrum stays blind until then: develop only on the
  mocks and on its synthetic twin.

## Current status
Phases 0, 1, 1-bis, 2, 2-bis, 3, 3-bis, 4, 5, 5-bis, 6, 6-bis, 6-ter and
6-quater are done (October-November 2026; short summaries and the result of
each in `ROADMAP.md`, the detail in `ROADMAP_HISTORY.md`). Release 2.6.0
(10 November 2026) is on GitHub with the macOS app, not on PyPI. The
results on the mock benchmark are those of 2.2.0
(`benchmarks/results/validate_mocks_baseline.json`) except `metals_auto`,
which the benchmark does not run (phase 6-ter, PR 7-quater): phases that
must not change them are checked with `validate_mocks.py --compare`.

Phase 6-quater (review of the recipe parameters, PRs #752-#768, closed in
PR 13 on 13 November 2026) left: aliases for renamed parameters and
classes, mypy at zero errors, one set of column names, the units of the
flux (`arb_flux`), the menus from the registry, the shorter MCP
descriptions (`SHARED_PARAMS`); the numbers of the benchmark unchanged by
construction. What still holds is in "Rules and decisions in force"
(`ROADMAP.md`). Next: PR 14, the release 2.7.0 (the version, the CHANGELOG,
`switcher.json`, the slow tests and the full benchmark on the final master),
then phase 7 (robust components, measured on the benchmark of phase 4; the
real analysis of the MCP server goes there) and phase 7-bis (small fixes: the
register of 3.2 says which findings go to each).

Open constraints:
- #549 (LICENSE and first PyPI publication) waits for INAF; while it is open
  each new PR goes on its own `claude/...` branch (say which).
- Waiting for data from Guido: sessions for the log replay, sample files for
  the loaders, X-Shooter VIS/NIR data (F15).
- One PR at a time: wait for the merge before opening the next.
- Questionable behaviour found along the way goes in the findings register
  (`ROADMAP.md`, section 3.2), reviewed with Guido; anything that loses user
  data or crashes is reported at once. Update this section and the progress
  table of `ROADMAP.md` at the end of each phase.
