# Astrocook – Development Roadmap

Working document, based on an analysis of the v2.0.1 code (October 2026).
It collects the assessment of the current state, the consolidation plan for
v2 and the vision for a future v3. Update it at the end of each phase.
The detail of the concluded parts (phases 0 to 6-bis, closed findings and issues, the MCP pilot) is in `ROADMAP_HISTORY.md`, moved verbatim: search it with `grep`, do not read it whole.

---

## 1. Current state (v2.0.1)

The state of the code when this roadmap was written (October 2026): ~53k lines of Python, of which ~26k active V2 code and ~22.5k in `astrocook/legacy/` (removed in phase 3); I/O, flux operations, continuum, identification and Voigt fitting, a PySide6 GUI and a scriptable API. Its description of the layout, strengths and weaknesses is superseded by phases 1 to 6-bis (the architecture is in `CLAUDE.md`, the weaknesses became the findings register, section 3.2). Full text: `ROADMAP_HISTORY.md`, "1. Current state (v2.0.1)".

## 2. Scientific assessment

### Voigt fitting – solid
- Exact profile via the Faddeeva function (`wofz`); thermal and turbulent b
  added in quadrature.
- Gaussian convolution with constant R, per-pixel R or per-component R.
- `least_squares` (TRF) with physical bounds and `x_scale='jac'`;
  constraints between components; iterative addition of components with an
  AIC criterion.
- Bayesian engine (emcee; dynesty with evidence).

**To check / missing:**
1. **No oversampling:** the model is evaluated at pixel centres. For lines
   that are narrow compared to the pixel, this can bias b and N.
   *First check to do.*
2. **Uncertainties not scaled by the reduced χ²:** the comment in
   `fitting/voigt_fitter.py` (~line 936) says they are, the code does not.
   To be clarified whether this is intended. *First check to do.*
   Checked (phase 1, item 9): with a correct error spectrum the unscaled
   errors are right; with a wrong one they are wrong by the same factor.
   V1 scaled them (lmfit default). Decision pending (F29).
3. Gaussian instrumental profile only (no tabulated LSF, useful for
   ESPRESSO/UVES).
4. No zero-level offset and no free local continuum in the fitting window.

### Continuum – weak point
Current algorithm (`find_absorbed_regions` + `fit_continuum`): running mean
with iterative sigma clipping (different scales inside/outside the Lyα
forest), linear interpolation over absorbed regions, global "fudge" factor,
Gaussian smoothing.

Limitations: in the high-z forest the level starts too low, and a single
global factor cannot correct an error that grows with z; `/8` window
heuristic inherited from V1; QSO template is a composite only, PCA not
implemented; no continuum uncertainty and no propagation to N, b.

Directions (increasing effort):
1. Local continuum (low-order polynomial) fitted together with the lines.
2. Spline with adaptive knots and uncertainty (e.g. Gaussian process).
3. Prediction of the forest continuum from the red side (PCA/ML on mocks).
4. Correction based on the redshift evolution of the mean absorption.

### System identification – sound but limited
Current algorithm (`SpectrumV2.identify_lines`): absorbed regions from the
mask, merged when closer than a velocity threshold
(`merge_regions_by_velocity`); for each multiplet, candidates from the
kinematic overlap of the apparent optical depth profiles
(`_find_kinematic_doublet_candidates`), rated with a score based on the
relation F₁ = F₂^γ (`rate_doublet_candidate`). The V1 approach (product of
significances on a common z grid, `_abs_like`) is no longer used.

Limitations: only the first two lines of each multiplet, one series at a
time; no global conflict resolution; F₁ = F₂^γ holds for the intrinsic
profile, not after convolution; yes/no outcome, not probabilistic; Lyα
interlopers not modelled.

Directions:
1. Probabilistic identification as a comparison between hypotheses.
2. Identify → fit → residuals → identify loop.
3. Detectors trained on mocks for specific species (CIV, MgII, DLAs).

### Validation
Before new algorithms: quantitative validation with simulated spectra
(known parameters, various S/N and resolutions) and comparison with VPFIT or
other tools. Possible basis for a paper on Astrocook 2.

---


## 3. v2 consolidation

**Principle:** build the safety net first, then touch the code. Each step is
a small PR, with the test suite green before and after.

### Rules and decisions in force
Rules and decisions taken in the concluded phases that still hold, copied verbatim from their phase (or shortened where the text was a chronicle), each with its origin. The rest of those phases is in `ROADMAP_HISTORY.md`.

#### Mock benchmark and baseline (phase 4; PRs #602-#604, #631)
`benchmarks/validate_mocks.py` runs the chain of the pilot on a grid of mocks with known parameters. The baseline file is `benchmarks/results/validate_mocks_baseline.json` (that of 2.2.0; the 2.1.0 one is kept as `validate_mocks_baseline_2.1.0.json`). `validate_mocks.py --compare` prints the key numbers of a new run next to it. When a PR must be measured on it (only the PRs that can change the results of fitting, identification or the automatic pipelines), and when the slow tests are enough: `CLAUDE.md`, "Slow tests and benchmark" (agreed with Guido, October 2026).

**Baseline** (version 2.1.0, October 2026; PRs #602, #603 and the PR that
added this paragraph). The full grid of `benchmarks/validate_mocks.py`
(135 mocks, 26 minutes with 4 workers, deterministic: two runs gave the
same numbers). Its summary is in
`benchmarks/results/validate_mocks_baseline_2.1.0.json` (until 2.2.0
in `validate_mocks_baseline.json`, now the baseline of 2.2.0, see
"Result of the phase"); `--compare` prints
the key numbers of a new run next to it, and a slow test runs the reduced
grid (16 mocks) against loose thresholds, to catch regressions.

   - identification: a candidate is correct when its series and z match a
     true system within 20 km/s; precision (purity) and recall
     (completeness) per species, in bins of logN and S/N, recall counted
     on the systems above a detection limit (proposed: the weaker line of
     the doublet detected at 3 sigma in the noiseless mock, from its
     equivalent width and the S/N; to be confirmed with the benchmark); duplicates (F49) and
     candidates beyond z_em (F55) counted apart; the fraction of the
     absorbed pixels explained;
   - fit: true and fitted components matched in z; bias and pulls of z,
     logN and b, spurious components, components at the bounds and with
     absurd errors (#517, F52), reduced chi2.
   Grid: S/N 10, 30, 100; R 8900 (X-Shooter VIS), 45000 (UVES),
   140000 (ESPRESSO HR); z_em 2.5, 3.5, 4.5; 5 random seeds each (135
   mocks, a few minutes: the chain takes ~1 s on a mock).

*Still to decide, when the item comes:* the decisions listed in each
issue (the proximity margin of F55: 3000 km/s, Guido). The detection limit of the recall
stays at 3 sigma (Guido, with the baseline below: the identification
should improve down to it).

**Result of the phase** (master after #631, October 2026; same grid,
24 minutes with 4 workers), against the baseline above:

| | 2.1.0 | End of phase 4 |
|---|---|---|
| Continuum, forest (5th percentile) | 0.936 (0.820) | 0.943 (0.830) |
| Continuum, red side | 1.002 | 1.002 |
| CIV purity / completeness / duplicates | 0.986 / 0.849 / 505 | 0.993 / 0.874 / 168 |
| SiIV purity / completeness / duplicates | 0.967 / 0.686 / 267 | 0.989 / 0.690 / 69 |
| MgII purity / completeness / duplicates | 0.980 / 0.906 / 345 | 0.988 / 0.902 / 96 |
| Fit: matched / spurious / flagged | 1534 / 428 / 413 | 1170 / 11 / 22 |
| Pull width z / logN / b | 1.78 / 2.29 / 2.07 | 4.14 / 5.38 / 4.34 |
| Pull width, single components (one fitted near) | 1.05 / 1.27 / 1.04 | 1.10 / 1.03 / 1.15 |
| System column: fitted / spread / pull width | 828 / 0.10 dex / 1.24 | 840 / 0.077 dex / 2.27 |
| System column spread at R = 8900 | 0.39 dex | 0.095 dex |
| Median reduced chi2 | 1.07 | 1.25 |

Completeness against the significance of the weaker line: 0.10 at 3-5
sigma (was 0.05), 0.44 at 5-10 (0.30), 0.91 above 10 (0.91).

*Left for the next phases*: F52 (absurd errors without NaN, partly
covered by the flags), F56 (forest continuum at high z and low S/N:
the emission model of phase 7), F57 (completeness down to 3 sigma), F58
(decomposition), F59 residual (+0.05 dex at R = 8900), the components
that reach `max_logN` = 19 in the forest fits (phase 7). With the release the
baseline file became that of 2.2.0 (this run; Guido's decision): the
later phases are measured against it. That of 2.1.0 is kept in
`benchmarks/results/validate_mocks_baseline_2.1.0.json`.

#### One default per recipe parameter (phase 4; #573, F40, #605)
A recipe parameter has one default, the same in the GUI, the scripts and the MCP server (the schema of the recipe, now read through the registry); a new parameter follows the rule (for example the z window of #518).

#### GUI benchmark (phase 5; `benchmarks/bench_gui.py`, #634; agreed with Guido, October 2026)
Run by hand, like `bench_fitting.py`, before and after each change of the GUI that can affect its speed, on the same machine; the PR states the numbers.

**GUI benchmark** (`benchmarks/bench_gui.py`, run by hand like
`bench_fitting.py`): the GUI offscreen on the small J0103-1305 archive of
the tests (2156 pixels, 764 components), four cases (a fifth,
`select_row`, added with the fifth fix below), each ended after the
pending Qt events (deferred draws, the result processed with a `QTimer`)
and the median of 5 runs: the first opening of the System Inspector in a
new main window, a second opening, a full redraw of the main plot, and
*Refit Selected* on the first CIV component (z = 3.0085), split into the
recipe in the worker thread and the rest (GUI). `--profile DIR` writes a
`cProfile` file per case. Before the phase (2.2.0, cloud container,
Python 3.13):

| Case | Time (s) |
|---|---|
| `inspector_first` | 3.96 |
| `inspector_again` | 2.25 |
| `main_plot` | 0.21 |
| `refit` | 12.96 |
| `refit_recipe` | 11.19 |
| `refit_gui` | 1.78 |

Runs on the same machine differ by ~10%. In the profile of `refit`, the
table is again the largest part of the GUI time (the link-colour delegate
and `SystemTableModel.data`). The recipe itself (~11 s for one CIV
component, `refit_all` with `group_depth` 2 and `max_nfev` 200) is the
fit, not the GUI: outside this phase.

- GUI benchmark (`bench_gui.py`, the best of two runs, alternating, on one
  machine; the end of phase 5, #663, against master after #688): first
  opening of the System Inspector 0.56 → 0.50 s, next openings
  0.34 → 0.28 s, redraw of the main plot 0.127 → 0.099 s (the legend
  removed, F86), selection of a row 0.174 → 0.144 s, the GUI part of a
  refit 0.61 → 0.57 s; the fit itself unchanged (10.7 s). No regression
  from the new features.

#### Parallel tests (phase 5)
**Parallel tests** (agreed with Guido, October 2026, before the split of
the GUI): `pytest-xdist` in the `test` extra, `pytest -n auto` in the CI.
Measured on 4 cores (as the GitHub runners): the fast suite 440 → 161 s,
the same tests passing and the same coverage (67%). The slow job stays
serial (one test of ~10 minutes dominates it, and each process would
rebuild the costly fixtures). During the split of the GUI the targeted
tests run locally, the full suite in the CI.

#### Screenshots of the tutorials (phase 5)
**Screenshots of the tutorials** (agreed with Guido, October 2026): the
images of the tutorials are rebuilt by a script,
`docs/screenshots/make_screenshots.py`, which runs the real GUI offscreen at
twice the logical resolution, one *scene* per image (spectra, recipes, view,
menus and dialogs, crop), on two spectra provided by Guido ("My QSO", "My
other QSO", kept private in the repository `DAS-OATs/astrocook-data`,
which can also hold the other data that cannot be published). The images no
longer need to be taken by hand on the Mac after each change of the GUI;
`tests/test_screenshots.py` runs two scenes in CI, so a change of the GUI
that breaks the script shows up in its pull request. Writing the scenes
found F63 (data loss in *Set session properties*) and F64 (invisible
checkboxes on Linux and Windows, fixed).

#### Texts and colours of the GUI (phase 5-bis; PRs #664-#688, each change agreed with Guido on screenshots)
One rule for the texts: Title Case in menus and buttons, sentence case in options; "Show ..." instead of "Toggle ..."; whole words instead of abbreviations where there is room ("Set Resolution...", not "Set R..."); the session log in a fixed-width font. The colours of the interface are in one module, `gui/theme.py` (a palette per role, light and dark, dark mode following the system or chosen in *View > Appearance*): none is written in a stylesheet, and a new role is added to both palettes (F89, F90 were the roles missing from the dark one). The plots stay white; *View > Plot Colours* offers the standard palette or a colour-blind safe one. The screenshots of the tutorials are rebuilt at each PR that changes the GUI.

#### Rule of the log messages (phase 6; PR 1, #691; applied by #693-#701)
**Rule of the log messages** (proposed in PR 1, agreed with Guido;
applied package by package from PR 4):
- *Levels.*
  - DEBUG: for the development (internal steps, intermediate values);
    never needed by the user.
  - INFO: what the user wants to read in the log window: what an
    operation did, with the numbers that describe its result ("Fitted 12
    components in 3 groups"); a few lines per recipe, not one per
    iteration.
  - WARNING: the operation gave its result, but the result needs a
    check, or something was done instead of what was asked (a value
    replaced by a default, a column dropped, a step skipped); the message
    says what was done instead.
  - ERROR: the operation did not give its result; the message says why
    and, when it can, what to do.
  - Progress (PR 5a): the lines that the progress dialog of the long
    recipes reads stay INFO, one per step or group ("Fitting group 3
    of 120: CIV at z = 2.3406."): the dialog shows the last INFO as its
    status; `>> PROGRESS: n` sets its bar.
- *One log per error.* An unexpected exception is logged once, where it
  is caught at a boundary (rows B), with the traceback; the code below
  does not log and re-raise (rows R), it puts the context in the
  exception message (`raise ValueError(...) from e`).
- *Texts.* Sentences in English, starting with a capital letter; no
  prefixes with the name of the function or of the step (`"rebin: "`,
  `"Auto-continuum: "`, `"Split: "`): the logger of the module tells
  where a message comes from; no arrows, banners or indentation inside
  the message (`"  -> Estimated R ~ ..."`, `"--- Executing ---"`); values
  with their units, parameters with the names of the dialogs; f-strings
  are fine.
- *The same situation, the same words.* A message used by more than one
  module, or a situation reported in several places (a missing
  continuum, `z_em` not set, a session not found), is defined once in
  `recipes/messages.py` (#566).
- *Library warnings.* A cosmetic warning of a library is silenced where
  it is raised (`warnings.catch_warnings`), never with a global filter.

The texts taken from V1 (`recipes/messages.py`) change too where they do
not follow the rule: the log of a session records the recipes and their
parameters, which the replay reads; the messages are not stored in it.

#### Broad except clauses (phase 6; PR 1, #691, and #694-#701)
**Inventory of the broad except clauses** (PR 1, master after release
2.3.0, `a5b506a`): 189 clauses in `astrocook/` that catch every
exception (12 bare `except:`, 177 `except Exception`, none with
`BaseException`), found by walking the syntax tree.
`tests/test_broad_excepts.py` holds their number per file: a new broad
clause fails the test, and a PR that removes or narrows one lowers the
count there and marks its row here as done. Ruff's `BLE001` was
considered and not used: it does not report a clause that re-raises or
logs with `exc_info=True`, so it misses 59 of the 189 (the entry points
of the recipes among them).

Action proposed for each clause:
- **N** – narrow: catch only the exceptions the code expects; the
  behaviour does not change.
- **B** – boundary: a recipe entry point, a worker thread, a Qt slot, a
  tool of the MCP server, where an unexpected error must not stop the
  application: stays broad, logged once with the traceback
  (`exc_info=True`) and, in the GUI, shown to the user; the `return 0` of
  the recipes becomes an exception in phase 6-bis.
- **C** – cleanup: removing a matplotlib artist, disconnecting a callback,
  closing a dialog, where a failure is expected and harmless: narrowed and
  kept silent, with a comment.
- **R** – remove: the clause adds nothing (`raise e`, or a log line
  before a re-raise, so the same error is logged twice, the second time
  by the caller); the context goes into the exception message.
- **D** – discuss: after a failure the code goes on with a result that
  looks valid (unchanged data, an empty list, a 0, a step skipped) and
  says so only in the log, or not at all (F95).

The inventory of the 189 clauses, with the decisions on each row (rows D of core, fitting and io, of the GUI and the hook, of the recipes), is in `ROADMAP_HISTORY.md`, phase 6. At the end of the phase 47 clauses were left, all boundaries.

#### Conventions of the recipe API (phase 6-bis; PRs #705-#723)
**Result of the phase** (November 2026, PRs #705-#723). The recipe API,
without a change of the numerical results except where Guido decided it
(F102, F103):
- Registry: one read-only description of each recipe (`RecipeSpec`,
  `ParamSpec`), the only reader of the schemas (F79, half of F12), read by
  the GUI, the scripts, the log and the MCP server; per parameter its
  type, bounds, choices, kind (session, column), whether it is offered to
  agents; per recipe what it requires and provides. Every recipe and
  parameter described for a reader and a model (F48).
- Typed parameters: one conversion at the boundary (`typed_recipe`),
  numbers and booleans or the strings of the dialogs and of old logs;
  values without a physical meaning refused (F46); the log records the
  values used, defaults included (#573). Batteries of edge inputs for all
  the recipes, compared before and after each PR (no valid case changed);
  mutation testing of the conversion and of the registry (Guido's
  question on the tests adapted with the code).
- Exceptions: `RecipeError` (`InvalidParameter`, `RecipeFailed`,
  `MissingPrerequisite`) that names the recipe, the steps and the reason;
  no recipe returns 0 any more; the failures logged once, where they are
  caught.
- Structured result: `run` gives a `RecipeResult` (the session, what
  changed by a generic comparison, the warnings, the steps, the values
  used, the reasons of what was removed or discarded, the candidates of
  the identification with region and score, the flag `not_fitted`, the
  steps provided, and the empty fields for phase 7); the GUI says it in the
  status bar and the log.
- Prerequisites: the resolver, the policy of the scripts (`provide`,
  `raise`), the GUI asks for the continuum as for z_em and the resolution,
  the MCP server says what is missing; a step provided is a state and a
  log entry of its own.
- Renamed parameters (phase 6-quater, PR 2): a parameter renamed keeps its
  old names as `aliases` in the registry, silent and permanent; the logs of
  the archives, the scripts and the MCP tools that use an old name still
  run; the log, the methods text, the Log Scripter, the dialogs and the
  schemas show the new name only; the old and the new name together are an
  error. A renaming adds the old name to `aliases` and a row to the table
  below (and to `RENAMED` of `tests/test_recipe_aliases.py`).

  | Recipe | Old name | New name | PR |
  |---|---|---|---|
  | `file.export_ascii` | `path` | `folder` | 6-quater, PR 7 |

- *Prerequisites of the recipes* (proposed in phase 6, PR 5c, after an
  analysis asked by Guido, November 2026; to be refined with him before
  the work starts). What a recipe needs and does not find, and who fills
  the gap, today in three uncoordinated ways: the GUI asks for `z_em` and
  the resolution (`_check_and_handle_requirements`, two lists written by
  hand, out of step with the recipes: e.g. `refit_all` with a region needs
  `z_em` and is not listed), then resumes the recipe; some recipes run
  another one with its defaults (the continuum by `estimate_auto` in
  `populate_from_identification`, `fit_component`, `optimize_system`,
  `refit_all`; the mask by `find_absorbed` in `identify_lines` and
  `components_auto`), in the GUI, the scripts and the MCP server alike;
  others run degraded with a WARNING (`z_em` not set: the forest treated as
  the rest of the spectrum). Guido's goal: a flexible system that
  interprets what the user needs, not one that only says "run that other
  recipe first"; the same principle in a different form for each way of
  use, because who decides how to fill the gap is different:
  - *MCP*: the agent decides. That it notices that the continuum needs
    `z_em`, and tries to estimate it by itself even without a recipe for
    it, is a feature, as long as it is shown clearly to the user; when a
    recipe fills the gap one day (as `estimate_auto` does for the
    continuum), the agent should see it and choose to run it, with a
    notice. The server makes the prerequisites explicit and machine
    readable instead of filling them silently: a tool that cannot run
    answers with what is missing and the recipes that provide it (e.g.
    the continuum is missing; `continuum_estimate_auto` provides it and
    needs `z_em`; `z_em` has no provider, it is set with `set_properties`);
    an option to fill them with the defaults is possible, the step then
    reported in the result.
  - *GUI*: the user decides at the time: the dialog of `z_em` generalized
    to any prerequisite, with the choices taken from the registry (for
    the continuum: estimate with the defaults, open Auto-estimate
    Continuum, cancel), the recipe resumed afterwards, as now.
  - *Scripts*: nobody can decide at the time, so the author of the script
    decides beforehand: a policy declared per session or per call
    (e.g. `on_missing="provide"`: the providers run with their defaults,
    with a WARNING; `"raise"`: a `MissingPrerequisite` exception that names
    what is missing and the call that provides it), possibly per
    prerequisite; the default to be chosen so that the numbers do not
    change in phases 6 and 6-bis (the continuum provided, as now).
  - *Common parts*: the registry declares what each recipe *requires*
    (`z_em`, resolution, continuum, mask) and what it *provides*
    (`estimate_auto` the continuum, `set_properties` `z_em` and the
    resolution, which it can also estimate from "N pixels"), once, for
    the three uses; one resolver in the core (pure: the state of the
    session and a recipe in, what is missing and who provides it out);
    in every use, a step added to fill a gap is a log entry of its own,
    with its reason ("Auto-estimate Continuum, run before Refit All: the
    spectrum had no continuum"), for the reader of the log and for the
    replay (today the log records only the outer recipe).
  - Done in phase 6 (PR 5c), independent of this: the failure of the
    automatic continuum handled, a fit without a resolution said at
    WARNING, `refit_all` saying that it estimates the continuum. The
    lists of the GUI move to the registry here (Guido).

#### Lessons for the full server (section 4.2, MCP pilot report)
**Lessons for the full server** (section 4; not changed in the pilot):
- tools to inspect (list components, describe a session) and to act
  (delete, refit, add a component), and to duplicate, list and close
  sessions;
- destructive recipes report what they remove, and offer a dry run or a
  "keep existing" option, or the server keeps the previous state;
- flags: the tolerance of "at a bound" is too tight (logN = 10.013 is not
  flagged), and the flags above;
- progress notifications for long recipes (the recipes already log
  `>> PROGRESS: n`);
- separate limits per list, filters, a metals-first option.

#### Decisions of phases 2 and 3 still in force (agreed with Guido, October 2026)
- Docker containers (`containers/`, broken since the V1 launchers moved):
  decided with Guido (October 2026): no renewal for now; `containers/` is
  removed in phase 3; a headless image comes with the full MCP server
  (section 4); an image with the GUI only if users ask for it.
- saving in the V1 `.acs` format removed (Guido's decision); V1 archives are
  read by a native reader, `io/v1_archive.py`, checked against the former
  reader on Guido's archives and on two small anonymised test archives
  (`tests/data/v1_*.acs`); F42, F43, F44 fixed along the way;
- logging is left to the application (Guido's decision): importing the
  package no longer configures it; new tutorial *Using Astrocook from
  Python*;
Data files used only by the V1 code, removed from `astrocook/data/`
(October 2026, agreed with Guido): `defaults.json`, `espr_spec_form.dat`,
`table2_inoue.dat` and the telluric data. The telluric data, useful for a
future telluric correction in V2 (e.g. masking of telluric lines), are kept
as assets of the pre-release
[`data-telluric-2026`](https://github.com/DAS-OATs/astrocook/releases/tag/data-telluric-2026)
(not a software release): `sky_telluric.fits` (ESO SkyCalc sky model,
300-2500 nm, 18 columns; V1 used `lam`, `flux_ael`, `trans_ma`) and
`telluric.dat` (360-808 nm, value and T/F flag, origin to be documented),
with their SHA-256. `telluric.fits` was byte-identical to
`sky_telluric.fits`. When V2 needs them, download them on first use and
cache them (e.g. `astropy.utils.data.download_file`), checking the
SHA-256; all files are also in the git history.

### Phase 0 – Environment
Done (October 2026). SessionStart hook for cloud sessions, `pytest-qt` and `pytest-cov` in the `[test]` extra, Ruff with the pyflakes rules only; the four tests out of date with the API fixed, the suite green. #476 and #481 (ASCII export without `import re`, found by Ruff) fixed with a regression test.
Full text: `ROADMAP_HISTORY.md`, "Phase 0 – Environment".

### Phase 1 – Safety net (pilot)
Done (October 2026). The safety net before touching the code: numerical snapshots (pure functions, continuum, identification, Voigt model and fit, `SpectrumV2` methods, the full chain on J0103-1305 and on a mock, with intent checks), V1/V2 parity tests, the mock generator, recipe-level and GUI smoke tests, `.acs2` round trips, the error-scaling check, mutation testing (`mutmut`, by hand), timing benchmark of the fit, and the findings register (section 3.2, F1-F33). The slow tests (full chain on J0103-1305) run in a separate CI job. Still waiting for data from Guido: sessions for the log replay, sample files for the loaders.
Full text: `ROADMAP_HISTORY.md`, "Phase 1 – Safety net (pilot)".

### Phase 2 – Cleanup and packaging
Done (October 2026), PRs #543-#554, release 2.0.3, no numerical change. Obsolete files removed (~54 MB of backup data in every package), `pyproject.toml` cleaned (placeholders, unused dependencies), the `astrocook` command, `extractall` with `filter="data"` (F35), F18 and F31 closed. LICENSE (GPLv3 or later) waits for INAF: PR #549 stays a draft (see "Releases").
Full text: `ROADMAP_HISTORY.md`, "Phase 2 – Cleanup and packaging".

### Phase 2-bis – Leftovers of phase 2
Done (October 2026). F34 (window icon); the ~80 pyflakes warnings cleaned up, each checked by hand, and the full `ruff check .` in CI; the Sphinx warnings; the optional `[bayes]` extra (emcee, dynesty, corner). Found F36 and F37 (#559 and #560, fixed in phase 3-bis).
Full text: `ROADMAP_HISTORY.md`, "Phase 2-bis – Leftovers of phase 2".

### Phase 3 – Decoupling from legacy
Done (October 2026), PRs #565-#584. `astrocook/legacy/` (~22,500 lines), `containers/` and the `[legacy]` extra removed (the wheel from 79.7 to 2.6 MB); V1 `.acs` archives read by `io/v1_archive.py` (saving removed); the V1 routines of the parity tests kept as verbatim copies in `tests/reference/v1/`; logging left to the application; the documentation published per version (#569, #570). Fixed F42-F44; found F38-F41. The decisions still in force are in "Rules and decisions in force".
Full text: `ROADMAP_HISTORY.md`, "Phase 3 – Decoupling from legacy".

### Phase 3-bis – Small fixes before 2.1.0
Done (October 2026), PRs #587-#592, then release 2.1.0 (GitHub and macOS app, not PyPI). Branching recipes logged only in the new session (#575, F41); the core no longer imports the GUI (F6, `core/debug_hooks.py`); stable link colours (F13); `dbtur` stored after a fit (#559, F36); the version in `.acs2` archives (F45); `norm` of `rebin` no longer offered (#560, F37, Guido's decision). Note on F8 for #526.
Full text: `ROADMAP_HISTORY.md`, "Phase 3-bis – Small fixes before 2.1.0".

### Phase 4 – Scientific correctness
Done (October 2026), PRs #601-#631, then release 2.2.0. First a validation on mocks (extended mocks #602, `benchmarks/validate_mocks.py` #603, baseline of 2.1.0), then each fix measured on it: #573 (#605), #518 (#606), #568 (#607-#609, the single table of atomic data), F49 (#611), F55 (#612), F50 (#613), #465 (#614), F14 (#618), F61 (#619), F62 (#624), #526 (#615-#617, #625: resolution as a column, exact convolution), #623 (#626), #517 (#627-#629), F63 (#631).
Measured result against 2.1.0 (same grid): duplicates CIV 505 → 168, SiIV 267 → 69, MgII 345 → 96; purity 0.99 for the three doublets; spurious components 428 → 11; flagged 413 → 22; system column at R = 8900 0.39 → 0.095 dex; forest continuum 0.936 → 0.943. The pulls of the components are wider (1.78 → 4.14 for z) not because of the fitter (isolated components still ~1.1) but because the systems are fitted with fewer components (F58, phase 7). The baseline of the benchmark became that of 2.2.0 (table in "Rules and decisions in force").
Full text: `ROADMAP_HISTORY.md`, "Phase 4 – Scientific correctness".

### Phase 5 – Splitting large files
Done (October 2026), PRs #634-#662, release 2.3.0 with 5-bis. Pure refactoring and performance: the snapshots and the mock benchmark unchanged. The GUI benchmark (`bench_gui.py`, #634); the System Inspector ~7× faster to open (5.07 → 0.74 s) and the selection of a row ~3× (0.59 → 0.21 s); `spectrum.py` 2826 → 774 lines (five mixins), `main_window.py` 4831 → 541, `pyside_plot.py` 2884 → 1266, `system_inspector.py` 2239 → 1085; tests in parallel (440 → 161 s); coverage 67% → 81%; the screenshots of the tutorials rebuilt by a script. Fixed F7, F66, F74, F82 (data loss in the ASCII export), F83, F63, F64 and the calibration (#642); found F66-F84.
Full text: `ROADMAP_HISTORY.md`, "Phase 5 – Splitting large files".

### Phase 5-bis – GUI look and ergonomics
Done (October-November 2026), PRs #664-#688, release 2.3.0. Dark mode and *View > Appearance* (`gui/theme.py`), the toolbar of the plots redesigned, a status bar, one rule for the texts, the main plot without title and legend, the System Inspector refined, the *Identifications* window as a table of doublets, recent files and window geometry, *View > Plot Colours*. Fixed F67, F69-F71, F73, F75-F78, F80, F81, F85-F92 and #630; F93 left for phase 7. Snapshots and mock benchmark unchanged; numbers change only in the calibration (F69, F70) and in the errors of linked parameters (#630). GUI benchmark: no regression (redraw of the main plot 0.127 → 0.099 s). Coverage 82%.
Full text: `ROADMAP_HISTORY.md`, "Phase 5-bis – GUI look and ergonomics".

### Phase 6 – Error handling
Done (November 2026), PRs #691-#701, release 2.4.0. The broad `except` clauses 189 → 47, all boundaries (`tests/test_broad_excepts.py` holds the count per file): the others narrowed or removed, the fallbacks that went on with a result that looked valid (F95, 24 rows, each decided with Guido) now stop with their reason; one logger per module and one rule for the messages (#566; both in "Rules and decisions in force"); a hook shows the unexpected errors of the GUI; the fast suite from 123 warnings to 0, warnings are errors in the tests. Snapshots and mock benchmark unchanged. Fixed F94 (data loss when a save fails), F95, F96, F98, F99; F97 and F100 right after the phase.
Full text: `ROADMAP_HISTORY.md`, "Phase 6 – Error handling".

### Phase 6-bis – Recipe API
Done (November 2026), PRs #705-#723, release 2.5.0. One read-only registry of the recipes (`recipes/registry.py`), typed parameters converted at the boundary, exceptions (`recipes/errors.py`) instead of `return 0`, the structured result (`recipes/result.py`, `run`), the prerequisites (`recipes/prerequisites.py`), all the texts reviewed, the MCP pilot on the new API and the generator of the schemas of the tools. The conventions are in "Rules and decisions in force". Numbers unchanged except F102 and F103 (Guido's decisions). Found and fixed F49 (rest), F50 and F55 (parameters), F102-F105; F12 and F79 closed. Tests 1308 → 3288, coverage 86%.
Full text: `ROADMAP_HISTORY.md`, "Phase 6-bis – Recipe API".

### Phase 6-ter – Full MCP server
Named after 2.5.0 (November 2026, Guido): its own phase, because it is a
new product built over the recipes, not a change of them. Scope in
section 4 ("Full server"): all the recipes as tools, generated from the
registry (`mcp_server/schemas.py`, phase 6-bis PR 8, is the generator of
their parameters), structured outputs from `RecipeResult`, the
prerequisites said by the tools (`provide_missing`), the
natural-language log, a headless Docker image, trials with fresh agents
as in the pilot (section 4.2). Plan of the phase to be agreed with Guido in
its PR 0. The order matters little with the next phase: a parameter
renamed later reaches the tools by the registry.

**Plan of the phase** (agreed with Guido in PR 0, November 2026: the
proposal of each decision below accepted; the options kept for the
record):
- *Starting point* (master after 2.5.0).
  - The registry holds 48 recipes (file 2, edit 15, flux 6, continuum 5,
    absorbers 20); 7 are branching (a new session), 8 hidden in the GUI
    (`stitch`, `populate_from_identification`, `update_component`,
    `fit_component`, `update_constraint`, `update_constraints_column`,
    `add_linked_system`, `fit_bayesian`: called by other recipes or by
    the System Inspector); 6 parameters are not offered to agents. Each
    recipe declares what it requires (`z_em`, `resolution`,
    `continuum`) and provides (`set_properties`: `z_em`, `resolution`;
    `estimate_auto`, `fit_continuum`, `update_from_knots`: `continuum`).
  - The server of the pilot (`astrocook/mcp_server/`, ~900 lines) has
    five tools written by hand: `session_open`, `session_save`,
    `continuum_estimate_auto`, `absorbers_identify_lines`,
    `absorbers_doublets_auto`. Since 6-bis PR 8 they run through `run`
    (the registry, typed parameters, `RecipeResult`), say what a session
    lacks (`MissingPrerequisite`) and fill it only with
    `provide_missing`; their outputs (`outputs.py`) are summaries
    written per tool. `schemas.py` builds the JSON schema of the
    parameters of any recipe from the registry, tested on all of them.
  - Every session of the server is a `SessionHistory`, as in the GUI:
    each recipe a state of its own and a log entry, the previous states
    kept. No tool exposes them yet (no undo, no list of the states).
  - The log entry (`LogEntryV2`) is the name of the recipe and its
    parameters, with the values used (#573); what the recipe did is in
    the `RecipeResult` (changes, warnings, steps provided, reasons of
    the removals), not in the log.
  - Progress: four recipes of the absorbers log `>> PROGRESS: n`, read
    by the GUI (`gui/recipe_progress.py`); the server does not send
    them.
  - PySide6 and matplotlib are required dependencies of the package;
    `core`, `recipes`, `fitting`, `io` and `mcp_server` import neither
    (checked: no import of `PySide6`, `astrocook.gui` or `matplotlib`).
    No Dockerfile in the repository (the V1 image removed in phase 3).
- *Constraints.* No change of the recipes' numbers (the snapshots and
  the mock benchmark unchanged in the whole phase; the slow tests left
  to the CI job after the merge, the PRs touch the server, not the
  pipelines). The pilot's tools keep working until their generated
  replacements are in, with the same names (`<menu>_<recipe>`, agreed in
  the pilot), so that the docs page *Using Astrocook from an AI Agent*
  and the users' configurations stay valid. The server stays a thin
  layer: what an agent needs and a recipe can give goes into the recipe
  or `RecipeResult` (for the GUI and the scripts too), not into the
  server.
- *Principle: a server open to new recipes* (Guido, PR 0). Phase 7 and
  the phases after it will add recipes, to cover more of the analysis
  and give the agents more tools for practical problems: that is also
  why the code was put on solid ground. A recipe added to the registry
  reaches the server as a tool with no server code, as it reaches the
  GUI's menus and dialogs:
  - the tool, its schema and its description come from the registry
    only; nothing in the server names a recipe, except the list of the
    excluded ones, each with its reason;
  - the output of a tool is built from `RecipeResult` only, in the same
    way for every recipe: what is specific to a recipe (the score of a
    candidate, the regions of the forest) reaches the agent through
    `RecipeResult.extra`, filled by the recipe; the summaries of the
    pilot written per tool (continuum, identification, fit) become
    summaries chosen by what changed (the continuum column, the
    components, the identification metadata), not by the name of the
    recipe;
  - a new prerequisite (something a recipe requires and another
    provides) is declared in the registry and in the resolver
    (`recipes/prerequisites.py`), and the tools say it without a change
    of the server;
  - the instructions of the server (what the client shows the model
    first) list the menus and what they are for from the registry, not
    by hand;
  - a test adds a fake recipe to the registry and checks that the
    server offers it, with its schema, runs it and returns its result;
    another that every recipe of the registry is a tool or excluded with
    a reason (a new recipe that is neither fails the test).
  The number of tools grows with the recipes (some 55 at the end of this
  phase, with those of the server). Clients show every tool to the
  model, so many tools cost context and choice; the trials (PR 7)
  measure it. If needed, a fallback that stays generated from the
  registry: the tools grouped (a tool per menu with the recipe as a
  parameter), or a few generic tools (`recipes_list`, `recipe_describe`,
  `recipe_run`) next to the most used ones.
  Where to add recipes: decision 6 below records the gaps the agents
  meet (Guido, PR 0; section 4 had it for later).
- *Decisions* (1-5 agreed with Guido in PR 0, as proposed; 6 and 7 asked
  by Guido in PR 0, their design to be agreed in their PRs).
  1. *Destructive recipes* (the lesson of the pilot: `doublets_auto`
     removed 190 components, the agent protected the work only by
     opening a second copy).
     - (a) *Dry run*: every tool of a recipe takes `dry_run`; the recipe
       runs on the current state, the result (what it would add, remove,
       change) is returned and the state is not kept. Clear for the
       agent, but it costs a second run to apply (minutes for the
       pipelines), unless the server keeps the result of the dry run and
       a tool applies it.
     - (b) *History*: the server exposes the states it already keeps:
       `session_history` (the states, each with its recipe and what it
       changed), `session_undo`, `session_redo`, and `session_duplicate`
       to branch. Nothing is lost, no recipe runs twice; the agent
       must know that it can go back (said in the instructions of the
       server and in the result of a recipe that removed components).
     - (c) both: (b) always, (a) only for the recipes that remove
       components.
     Agreed: (b), with the result of every recipe listing what it
     removed (already in `RecipeResult.changes`) and saying how to undo
     it; (a) later if the trials ask for it.
  2. *Docker image*.
     - (a) A `Dockerfile` in the repository and a CI job that builds it
       and runs the stdio test of the server inside it; the users build
       it (`docker build`), the docs say how to register it in a client
       (`docker run -i --rm -v <data>:/data ...`: the files of the user
       reach the container only through a mounted folder).
     - (b) As (a), and the image published on the GitHub Container
       Registry (`ghcr.io/das-oats/astrocook`) by the release workflow,
       one tag per release. The users only `docker pull`; but it is a
       distribution of the software, like the macOS app, while the
       LICENSE waits for INAF (#549).
     Either way, the image: PySide6 is a required dependency (~0.7 GB installed
     with its Qt libraries), useless without a display. Options: keep it
     (simplest, a large image), or install the package without it in
     the image (`pip install --no-deps` and the other dependencies
     listed by hand, checked by the CI job), without changing the
     dependencies of the package. Making the GUI an optional extra
     (`astrocook[gui]`) would change the installation for every user: not
     proposed in this phase. Agreed: (a) in this phase, without
     PySide6 in the image; (b) with the first release after #549.
  3. *Which recipes become tools.* Guido (6-bis, PR 6): as many as
     possible, the basic ones too. Agreed: all 48, also the 8 hidden
     in the GUI (`fit_component`, `update_component`, ... are exactly the
     editing tools the agents of the pilot missed), except the ones a
     test lists as excluded with a reason. By the same rule, the
     candidates to exclude are offered: `update_from_knots` (section 4
     called it interactive, but the knots are a parameter an agent can
     give), `import_ascii` and `export_ascii` (they read and write files,
     as `session_open` and `session_save`), `fit_bayesian` only when the
     `[bayes]` extra is installed (its description says that it runs for
     minutes or more). The list of exclusions starts empty.
  4. *Agent decisions in the log* (section 4, item 3; Guido's trial: z_em
     estimated from the Ly-a peak, R from the header, nowhere in the
     log).
     - (a) a `reason` parameter, optional, on every tool, stored with
       the log entry of the recipe;
     - (b) a tool `session_note(session_id, text)` that adds an entry
       without a state (a note, skipped by the replay);
     - (c) both.
     Agreed: (c); (a) is where the agent states why it ran a recipe,
     (b) where it records what it did outside the server. Either needs
     a field of `LogEntryV2` (`reason`, and a kind for the notes), saved
     in `.acs2` and shown in the GUI's log: old archives open as before
     (the field absent).
  5. *Natural-language log*: a function of the core (not of the server)
     that writes the log of a session as text, from its entries, the
     reasons and notes of decision 4 and the changes of each state
     (compared as `RecipeResult` does); offered as a tool
     (`session_methods`) and in the GUI (*File > Export Log as Text*?).
     The text: one paragraph per step, with the references of the
     algorithms (e.g. the atomic data of *Atomic Data*). How it reads is
     Guido's to judge, on a draft in its PR. Agreed.
  6. *Gaps of the collection of recipes* (Guido, PR 0: know where to
     direct the work). A gap is something the analysis needed and no
     recipe provided. Two sources:
     - *declared by the agent*: a kind of note of decision 4
       (`session_note(..., kind="gap")`: what was needed, what the agent
       did instead, e.g. "z_em from the Ly-a peak, computed with Python
       on the flux"; Guido's trial had several); the instructions of the
       server ask the agent to declare a gap whenever it works outside
       the tools;
     - *seen by the server*: a prerequisite that no recipe provides (as
       `z_em` today), a tool that failed on a valid request, a recipe
       whose result the agent undid (decision 1: the next step after an
       undo says why, if the agent gave a reason).
     The gaps are saved with the session (its log, decision 4) and also
     appended to a local file of the user (one line each, JSON: date,
     version, kind, text, session), never sent anywhere; a command (and a
     tool) summarizes them, grouped, as the input of the plan of a phase
     (e.g. phase 7, phase 6-quater). Not stored: the data of the user
     beyond the name of the session and the text the agent wrote.
  7. *The GUI launched by the server* (Guido, PR 0: the user sees the
     progress of the analysis as it happens). The GUI is a separate
     process; options:
     - (a) *snapshot*: a tool `gui_open(session_id)` saves the current
       state to an archive and launches `astrocook <archive>`; the user
       sees the state at that moment and asks again for a newer one;
     - (b) *follow*: as (a), and the GUI follows the archive: the server
       saves it after every recipe and the GUI reloads it (a file
       watcher, `QFileSystemWatcher`; the GUI has none today), keeping
       the view (range, zoom, the open System Inspector) where it can.
       The followed session is read-only in the GUI (a banner says so;
       a copy can be made to edit it), so that the user and the agent do
       not change the same session at once;
     - (c) *shared session*: the GUI talks to the server (a local
       socket), the user's edits go back to the agent; the "review
       bench" of section 4, item 4, in full, but a protocol between two
       processes and the conflicts of two editors: not in this phase.
     Proposal: (b). Not available in the Docker image (no display): the
     tool says so. On macOS the app is launched as the trial did
     (`open -a Astrocook`), elsewhere the `astrocook` command.
- *PRs*, one at a time, each with characterization tests first:
  0. this plan (docs only);
  1. *generated tools*: one tool per recipe from the registry (schema
     from `schemas.py`, description from `brief` and `details`),
     calling `run`; the output from `RecipeResult` (summary, changes,
     warnings, notes, steps provided, the new session of a branching
     recipe under its own identifier); a test that every recipe is a
     tool or listed as excluded with a reason, and the test of the fake
     recipe (the principle above); the five tools of the pilot replaced,
     their names kept, their extra summaries (continuum, identification,
     fit) chosen by what changed, what `RecipeResult` lacks moved into
     `extra` by the recipes; the instructions of the server from the
     registry;
  2. *sessions*: `session_list`, `session_describe` (the spectrum, the
     columns, z_em, R, the continuum, the systems, the log),
     `session_duplicate`, `session_close`, and decision 1
     (`session_history`, `session_undo`, `session_redo`, or `dry_run`);
  3. *inspection*: `components_list` (filters by series, z, flags;
     with the fit group and the system, the pilot's question),
     `spectrum_read` (the columns in a range, rebinned to a number of
     points that fits a model's context), `spectrum_plot` (a PNG of a
     range with the model, as MCP image content; matplotlib, already a
     dependency);
  4. *progress*: the `>> PROGRESS` lines of the recipes sent as MCP
     progress notifications, and a message for the recipes that run
     long without them;
  5. *decisions in the log* (decision 4), the *gaps* (decision 6) and
     the *natural-language log* (decision 5), perhaps in three PRs;
  5-bis. *the GUI from the server* (decision 7): the tool, the file
     watcher and the read-only followed session in the GUI (with the GUI
     benchmark, `bench_gui.py`, and the screenshots of the scenes it
     touches);
  6. *Docker* (decision 2): the image, its CI job, the docs;
  7. *trials* with fresh agents, as in section 4.2 (the mock archive
     without its true continuum; J0103-1305; a request that needs the
     editing tools), and a report in section 4, with the measure of the
     number of tools (the principle above). Then *a real analysis*
     (Guido, PR 0): Guido uses the server on a spectrum and a scientific
     question of his own, as in the pilot but with the full server and
     the GUI following it; the criterion of success: a result he would
     sign after reviewing it in the GUI, and a log he could replay and
     describe. What fails becomes a finding (for phase 7 or the server),
     and the gaps recorded (decision 6) the first list of recipes to
     add;
  8. closing: the docs page *Using Astrocook from an AI Agent*
     rewritten for the full server, the result of the phase (docs only);
     then the release 2.6.0 (PR 9).
  Proposed (code review of 2.5.0; see "Code quality" after phase 7): a type
  checker in CI on `recipes/` and `mcp_server/`, ruff with `B`, `SessionV2`
  out of `core/`. Decided at the closing (PR 8, Guido, November 2026): none
  of the three is done in this phase; all go to the start of 6-quater,
  `B` as its first PR (small and mechanical), the type checker and
  `SessionV2` with the renaming of the `V2` suffixes (the same aliases).
- *Done.*
  - PR 1 (generated tools): `mcp_server/tools.py` makes a tool of each
    recipe of the registry (47 of 48), with the schema of
    `schemas.py`, the description from the registry (with what the
    recipe requires and provides, and whether it opens new sessions), and
    a generic output (`RecipeInfo`: what changed, listed, the components
    removed with their reason, `new_sessions`, `value`, `not_fitted`,
    and the summaries of the pilot chosen by what changed); the values
    checked by the registry (`typed_recipe`), unknown parameters
    refused; the instructions of the server list the menus from the
    registry; the messages name the tools (`(tool edit_set_properties)`).
    The recipes that read another session (kind `session`) find it by its
    identifier: the store is their context, as the GUI is
    (`session._gui.session_histories`). Branching recipes open new
    sessions, as in the GUI (`SessionStore.branch`: a copy of the history
    and log of the parent, then the recipe logged in the new session):
    `identify_lines` now opens a new session (the pilot added the
    identification to the same one). Tests: every recipe is a tool or
    excluded with a reason; a recipe added to the registry in the test
    (two, one returning a value) is a tool, checked, run and logged with
    no code of the server. Excluded: `export_ascii` (F106), until PR 1-bis.
  - PR 1-bis (F106, Guido's decision): `export_ascii` declares `prefix`
    and takes `overwrite` (False by default), and is a tool; no recipe is
    excluded. A test checks that every parameter of a method is in the
    registry, or listed with its reason (eight, all with a default:
    `alias_map` filled by the GUI, deprecated ones, `norm` of `rebin`,
    `series` of `update_components`; `template` of `fit_continuum`, not
    offered since 6-bis, to review in 6-quater).
  - PR 2 (sessions, decision 1): `session_list`, `session_describe`,
    `session_history`, `session_undo`, `session_redo`,
    `session_duplicate`, `session_close`, in `mcp_server/session_tools.py`
    (with `session_open` and `session_save`, moved from `server.py`). The
    store keeps per session where it comes from (opened, branched by a
    recipe, duplicated, and its parent), where the entries of this server
    begin in its log (an archive brings its own: `log_before` of
    `session_describe`), and the state last saved or opened: `unsaved`,
    and `session_close` refuses to lose changes unless
    `discard_changes=true` (the server's own guard against data loss).
    Undo and redo move the log with the state (`SessionHistory`, as in
    the GUI): a session saved after an undo replays without the recipe
    undone. The result of a recipe that removed components or columns
    says that `session_undo` restores them (`undo_hint`); the
    instructions of the server say that every state is kept.
  - PR 3 (inspection, `mcp_server/inspect_tools.py`): `components_list`
    (by z; errors, chi2 and flags of the last fit, the fit group, the
    velocity from the first component of the group with the same series,
    the constraints: 'fixed' or 'linked: <expression>'; filters by series,
    z and flags, pages with `offset`), `spectrum_read` (the columns in a
    range of x, averaged in blocks of pixels beyond `max_points`, for
    reading only: the session is not rebinned) and `spectrum_plot` (a PNG
    as MCP image content: flux, error, continuum, model, a tick per line
    of each component). The fit group is that of Refit All before it
    splits the large groups: the pilot inferred it from equal chi2. These
    tools change nothing and log nothing. Guido's question on the PR (why
    these tools, could they be recipes, could existing functions replace
    them): they are needed because the state of a session lives only in
    the server (the pilot's agent read the original FITS with its own
    Python, not the state after the recipes, and only with access to the
    disk); `export_ascii` writes every pixel to disk, without groups or
    constraints; the GUI is Qt and shows the user, not the agent; no plot
    without the GUI existed. Not recipes: they change nothing and would
    only add noise to the log, and the GUI has these views already. But
    what they compute belongs in the package, not in the server (the
    constraint of the phase): moved, before the merge, to
    `SystemListV2.fit_groups` (the grouping of Refit All, which now calls
    it: one criterion, not two copies), `SpectrumV2.columns_in_range` and
    the new `astrocook/plotting.py` (`plot_spectrum`, matplotlib's
    `Figure` without pyplot or Qt, for the headless image; the palettes of
    the plots moved there from `gui/plot_style.py`, which imports them);
    the tools are thin wrappers. Documented for scripts (*Using Astrocook
    from Python*, "Looking at a Session"). Seen along the way: the
    continuum of the J0103-1305 archive of the tests is stored without a
    unit (dimensionless, values in flux units), while a continuum
    estimated now has the unit of the flux; the tools report the units as
    stored.
    Refinements of the decisions, found in the code: `fit_bayesian` is
    offered always (offering it only with the `[bayes]` extra would name a
    recipe in the server; without the extra it fails with the message that
    names the extra to install).
  - PR 4 (progress): the reader of the `>> PROGRESS: n` lines moved from
    the GUI to the package (`recipes/progress.py`: `parse_progress`,
    `watch_progress(callback)`; the dialog of the GUI uses it, same
    behaviour); the tools get the `Context` of the call (an annotated
    parameter, not in the schema) and `mcp_server/progress.py` sends
    `ctx.report_progress(n, 100, last line of the log)`: for every recipe,
    none named. The sync tool runs in a worker thread of the SDK, so the
    notification goes to its loop with `anyio.from_thread.run`; the log
    records are read from the thread of the recipe only (`astrocook`
    logger lowered to INFO for the block and restored: the root is at
    WARNING and the INFO lines would not arrive). The value of a
    notification must increase: repeated lines (the pipelines start
    several steps at 0) are not sent, and the lines of a recipe stop at 99
    so that the end of the tool (100, "done in n s") always has its message.
    A recipe without lines sends a start and an end (decision of the PR,
    with Guido). Without a progress token the SDK ignores the call; if
    the client went away, one INFO line and no more notifications. stdout
    stays clean (tested over stdio). Left out, to add if the trials (PR 7)
    need it: a periodic heartbeat for recipes that are silent and long
    (the elapsed time as the value, no total; some clients restart the
    timeout of a call on each notification, not checked), and the Stop of
    the GUI as a cancellation. Numbers unchanged: no recipe, core or
    fitting code touched, so no slow tests.
  - PR 5a (decision 4, option c: the agent's decisions in the log):
    `LogEntryV2.reason` (why the recipe was run) and the notes,
    `LogNoteV2` (text, kind, `after`: the number of entries before it),
    in `HistoryLogV2.notes`, saved in `.acs2` (the keys `reason` of an
    entry and `notes` of the log, written only when present: a log
    without them is saved as before). Refinement of the decision, agreed
    with Guido on the plan of the PR: the notes are not entries of a kind
    of their own but a list apart, so that the entries stay 1-to-1 with
    the states (`SessionStore.entry_of`, undo and redo unchanged) and an
    archive with notes replays in the versions before (a note among the
    entries would be a recipe `note` that 2.5.0 cannot run). A note is
    undone with its state and dropped with it when a new recipe replaces
    it; a note written after an undo keeps the redo (it is not a state).
    Checked with 2.5.0 (its own code, from the tag): an archive with
    reasons and notes opens, the Log Scripter shows the recipes, the
    entries replay with the same numbers; saved again by 2.5.0 it loses
    the reasons and notes, without a word (expected: it does not know
    them). In the package, for the GUI and the scripts: `add_entry(...,
    reason=)`, `add_note`, `notes_at`, `timeline` (the entries and notes in
    the order written, for the natural-language log of PR 5c). The GUI
    shows them as comments in the Log Scripter (`gui/log_script.py`:
    `# reason: ...` before the recipe, `# note: ...` where it was
    written); *Run All* skips them as steps and keeps them in the log it
    rebuilds (it rebuilt the log from the text: a replay would have lost
    them). The server: `reason` on every recipe tool (a parameter of the
    server, as `provide_missing`), also in the sessions opened by a
    branching recipe (not on the steps provided first); the tool
    `session_note`; `session_history` (reason and notes per state) and
    `session_describe` (those of the log of an archive); a note makes the
    session unsaved, so that `session_close` does not lose it; the
    instructions of the server ask for reasons and notes. Numbers
    unchanged: no recipe, fitting or numerical code touched (the
    snapshots unchanged), so no slow tests.
  - PR 5b (decision 6: the gaps of the collection of recipes), plan agreed
    with Guido: `astrocook/gaps.py` (package, no import of `mcp` or the GUI):
    `record_gap`, `read_gaps`, `summarize`, and the command `astrocook-gaps`
    (`--since`, `--kind`, `--json`, `--path`, `--file`). The file, one JSON
    line each with only `date`, `version`, `kind`, `text`, `session` (the
    name of the session), never sent anywhere: macOS
    `~/Library/Application Support/Astrocook/gaps.jsonl`, Linux
    `$XDG_DATA_HOME/astrocook/gaps.jsonl` (else `~/.local/share/...`), Windows
    `%APPDATA%\Astrocook\gaps.jsonl`; `ASTROCOOK_GAPS_FILE` moves it or, with
    `off`, turns the record off (the tests redirect it to a temporary file:
    `conftest.py`). A write that fails never breaks a tool (one warning).
    Four kinds: `declared` (`session_note(..., kind="gap")`, a note of kind
    `gap` in the log) and, seen by the server (notes of kind `gap_seen` in the
    log, which do not make the session unsaved: nothing of the agent is lost
    by closing it): `prerequisite` (a prerequisite with no `automatic`
    recipe, z_em and the resolution; once per session and tool), `failed`
    and `undone` (`session_undo` takes an optional `reason`; the steps
    provided before a recipe are not counted, the store remembers their log
    entries). *Failed on a valid request*, told by the exception that
    `run_recipe` chains: not a gap the unknown parameters and sessions, the
    `InvalidParameter`, the `MissingPrerequisite` that a recipe provides, and
    the `RecipeFailed` (the data do not allow it, the agent has the reason:
    if it works around, it declares it); a gap the exception that is none of
    these (a bug, once per session and tool), recorded with the type of the
    error and not its message, which may hold a path. `gaps_summary` (tool)
    and `format_summary` group by kind and text with the numbers masked; the
    declared ones are free text, joined by the reader. Seen: the `z_em` or
    `resol` of `session_open` make a state of their own, so undoing it is
    recorded as a recipe undone (`edit_set_properties`); left as it is.
    Numbers unchanged: no recipe, fitting or numerical code touched, so no
    slow tests.
  - PR 5c (decision 5: the natural-language log), plan and draft agreed with
    Guido (the draft on a real session of J0103-1305 through the server:
    continuum, doublets, refit; tone, the reason first as written, the
    median and largest changes, the menus in parentheses):
    `astrocook/methods_text.py` (package level, as `gaps.py` and
    `plotting.py`: it reads the registry and the comparison of the recipes,
    which `core/` should not, see `SessionV2` out of `core/` below):
    `methods_text(history, menus=True)` returns `MethodsText` (paragraphs,
    references, `.text`): an opening paragraph (the spectrum, the version,
    the number of steps), one per recipe (the reason as written; what the
    recipe did; the parameters it does not quote that differ from their
    defaults; what it changed, from `compare_sessions` between its state and
    the one before: components added and removed by series with their range
    of z, the changed ones one by one up to three, else the median and
    largest changes of log N, b and velocity, the columns by name, the
    pixels; the flags of the fit on those components; the notes and declared
    gaps written there, as written; not the gaps seen by the server), and a
    closing one (the fit, if a recipe that cites SciPy ran and components
    are left; the sources of the atomic data of the transitions of the
    components that fall in the spectrum, from the `*_ref` columns of the
    table), then the references, sorted. The text of each recipe is in the
    registry (new keys `methods`, a template of the parameters or a function
    of them, and `references`, keys of `REFERENCES`), written for all 48; a
    test requires it of every recipe (a new recipe too); a recipe without it,
    or unknown (a log of a later version), gets a generic sentence. Nothing
    names a recipe in the module. The text is made on demand, never stored:
    a better wording reaches the old sessions. The states map to the entries
    by the two current indices (the GUI and the server copy the states of a
    branch); the entries before the first state (an archive opened) are
    described by their parameters only, and the text says so. References
    (Guido): Astrocook as Cupani et al. 2020; its algorithms (continuum,
    identification, AIC) as Cupani et al. 2026 (Proc. SPIE 14155, 141551X,
    doi:10.1117/12.3101756; volume and number added after the merge), or the papers after
    phase 7; SciPy for the fit, emcee and dynesty for the Bayesian fit; no
    Astropy or NumPy. The tool `session_methods(session_id, menus)` (a thin
    wrapper, changes nothing; the instructions of the server mention it) and
    *File > Export Log as Text...* (a dialog: the text read only, *Name the
    menus of the recipes*, *Copy*, *Save...*). Found: F107 (the labels of
    the menus are not the `brief` of the registry), F108. Numbers unchanged:
    no recipe, fitting or numerical code touched (new keys in the schemas,
    and their text functions), so no slow tests.
  - After PR 5c (Guido): F108 fixed, with the full reference of Cupani et
    al. 2026 and F107 deferred to phase 6-quater. *Set properties* with the
    object left at `'_current_'` no longer writes it again: it wrote
    `OBJECT = ''` on a spectrum without one (the change "object None → ''"
    of `session_open` with z_em or R), copied `OBJECT` over the ESO target
    name, and added "Object='...'" to the history of the spectrum.
    Characterization tests first. No numbers touched (the header and the
    history only), so no slow tests.
  - PR 5-bis (decision 7, option b: the GUI follows a session of the
    server), plan agreed with Guido: the marker next to the archive, and the
    branching recipes allowed on a followed session. `astrocook/follow.py`
    (package, no Qt or `mcp`, so that any program can be followed):
    `write_archive` (as `session_save`: the current state, the first state as
    the original data, the log; atomic, `_write_replacing`, F94), the marker
    `<archive>.follow.json` (`write_marker`, `read_marker`: the identifier
    of the session, the pid of the writer, the version; `writer_alive`, a
    signal 0 on POSIX, true on Windows), `remove`, `gui_unavailable_reason`
    (no PySide6, as in the Docker image of decision 2; Linux without
    `DISPLAY` or `WAYLAND_DISPLAY`), `gui_command` (macOS: `open -a
    Astrocook` if the application is in `/Applications` or
    `~/Applications`; else `astrocook` on the path, else `python -m
    astrocook`) and `launch_gui` (detached, `start_new_session`, the three
    streams on the null device: stdout is the channel of the server).
    Why a marker and not a flag of the command line: `open -a` passes the
    file as a FileOpen event, and to an application already open no
    argument arrives. The server (`mcp_server/gui_tools.py`): `gui_open`
    saves the session to `<temporary folder>/<session_id>.acs2` with its
    marker and launches the GUI; the store saves it again after every change
    of state or log, in one place (`refresh_gui`: after a recipe and its
    steps provided, undo and redo, a note, a gap seen), and a failed save is
    a warning of the next output of the session, never a failed recipe;
    `session_close` and the end of the server (`weakref.finalize` of the
    folder) remove marker and archive; the sessions opened by a branching
    recipe are not followed. The archive is not the user's save (the session
    stays unsaved). The GUI (`gui/main_follow.py`, a mixin of the main
    window): an archive opened with a marker (from the command line, the
    FileOpen of macOS, the dialog, drag and drop: `_load_session` and the
    launcher keep the path until the session is added) is followed: a
    `QFileSystemWatcher` on the file and on its folder (a replaced file is a
    new inode on Linux: the file is watched again at each event); the events
    of the archive only (not the temporary files of the writer), a debounce
    of 300 ms, a reload only if mtime or size changed, in a worker thread,
    one more try after a failure (the state shown is kept). The new state
    replaces the one shown, with the log of the archive, not autoscaled, so
    that a zoom is kept; the System Inspector keeps its row by the uuid of
    the component (saved in the archive). A banner (a `QToolBar` at the top,
    hidden otherwise; the cards float below it) says *Followed session*,
    read only, with **Make Editable Copy** (a duplicate with a copy of the
    log: *Duplicate* gives a new empty log). Refused on a followed session,
    with the offer of the copy: the recipes that are not branching (in
    `_on_recipe_requested` and before their dialog), the continuum editor,
    Undo and Redo (disabled), the scripts (*Run All* too), renaming.
    Allowed: saving, *Export to ASCII*, the branching recipes (their
    sessions are new and not followed). When the marker disappears (or its
    writer is dead at the opening) the banner says *No longer followed*,
    still read only. Opened again (a second `gui_open` with the application
    open, on macOS), the archive already followed is shown, not opened
    twice. Found: F109 (the original data of an archive are not read back:
    an archive opened and saved again has its state as opened as original
    data, and its log replays twice; data loss, reported to Guido at once).
    `benchmarks/bench_gui.py` did not run since phase 6-bis (its first
    refit asked for z_em, R and the continuum, and the benchmark waited
    forever): it now sets z_em = 4.072 and R = 45000 and answers *Estimate*.
    GUI benchmark, before and after on the same machine (s, medians of 5):
    inspector_first 0.425 / 0.445, inspector_again 0.276 / 0.267, main_plot
    0.087 / 0.085, refit_gui 0.584 / 0.518, select_row 0.119 / 0.118: the
    same within the noise (a session not followed does nothing new). The
    screenshots: a new scene, `agents_followed_session` (the banner), shown
    in *Using Astrocook from an AI Agent*; two scenes of the main window,
    made with and without the PR, are identical pixel by pixel (the banner
    is hidden for a session not followed). Saving the archive costs
    about 0.1 s per change on the J0103-1305 archive of the tests (only for
    a followed session). Numbers unchanged: no recipe, core or fitting code
    touched, so no slow tests.
  - After PR 5-bis (Guido): F109 fixed, by the rule of the sessions (each
    saved archive is a node: A, the `_orig` members; B; the log from A to
    B). `HistoryLogV2.previous` (the log of the archive a session was opened
    from: read only, never replayed, never in `to_dict`) and
    `HistoryLogV2.opened_from(log)` (an empty log with it, the entries undone
    and their notes left out, its own previous dropped: the logs of the
    earlier nodes are in their archives). Used where a session is opened to
    work on: the GUI (`_log_for_new_session`: the dialog, the command line,
    drag and drop, the start), the reload of a followed session, the store of
    the server (`add`); the reader (`load_session_from_file`) still gives the
    log of the archive, for scripts. Carried by the copies of a log (branches,
    the editable copy, *Run All*, which reruns the log of the session only);
    not by *Duplicate* (a new empty log, as before). Shown: in the Log
    Scripter, at the top, as `##` comments, which *Run All* skips and does
    not read as reasons or notes; `log_before` of `session_describe` and the
    notes of the first state in `session_history`; the text of the methods
    (the steps of the archive, by their parameters, as before). The original
    data have their metadata too (`_meta_orig.json`: constraints, identifiers
    of the components, identifications; no log), not read back. `first_entry`
    of the store is now always 0. Characterization tests first
    (`test_session_nodes.py`: the archive saved again had B as original data
    and the old log; *Run All* doubled the flux twice), then the same tests on
    the rule. Docs: *Getting Started* ("Saving and Closing", each saved session
    is a step), *Using Astrocook from an AI Agent* (the log of the agent read
    only; an analysis by an agent in one session, or the agent keeps its own
    log). No numbers of the recipes touched; the core and io touched, so the
    slow tests were run.
  - PR 6 (decision 2, option a: the headless Docker image), plan agreed with
    Guido: the `Dockerfile` (python:3.13-slim; the package with
    `pip install --no-deps .` and the other dependencies listed by hand,
    those of `pyproject.toml` and the `[mcp]` extra without PySide6, qtawesome
    and scienceplots, F110; `astrocook-mcp` over stdio as the entrypoint, not
    root, `/data` for the files of the user, the file of the gaps written there
    so that `--rm` does not lose it) and a `.dockerignore`; the dependencies of
    `pyproject.toml` are not changed. Tests of the suite, with no Docker:
    `tests/test_docker_image.py` reads the Dockerfile (the list is that of
    `pyproject.toml` but for `GUI_ONLY`, each with its reason; a stale entry of
    `GUI_ONLY` fails; nothing outside `gui/` and the launcher of the GUI imports
    PySide6, qtawesome or scienceplots; no PySide6 in the Dockerfile's code; not
    root; stdio entrypoint; the files it copies exist and are not ignored) and
    `tests/test_mcp_headless.py` (`build_server` in a subprocess: no PySide6,
    shiboken6, qtawesome, `astrocook.gui` or pyplot in `sys.modules`; and the
    server over stdio in a process where PySide6 cannot be imported). The
    stdio check, `tests/mcp_stdio_check.py`, is a script, not a test: a client
    that starts the server (the command after `--`), opens a mock, runs
    `continuum_estimate_auto`, saves, plots and asks for `gui_open`, which must
    answer that the GUI is not available (for the lack of PySide6 here, of a
    display elsewhere). The CI job `docker` builds the image, checks that
    PySide6, shiboken6 and qtawesome are not in it, and runs the script inside
    it (mounted; on a pull request only if what goes into the image changed,
    by `git diff` against the base: no third-party action). The image is not
    published (GHCR with the first release after #549). Docs: *Using Astrocook
    from an AI Agent*, section 2.1 (build, `docker run -i --rm -v <data>:/data`
    for Claude Code and the Claude desktop app, the paths inside the container,
    `gui_open` unavailable). Checked here: the contents of the image
    reproduced in a clean venv (the same list, `--no-deps` for the package, no
    PySide6): 64 tools and every step of the script pass; the Docker daemon is
    not available in the cloud session, so `docker build` itself is run only by
    the CI job. Numbers unchanged: no recipe, core or fitting code touched, so
    no slow tests and no benchmark.
  - After PR 6 (Guido): F110 fixed: `scienceplots` removed from the
    dependencies of `pyproject.toml`, from `Astrocook.spec` (the macOS build
    collected it: not rebuilt here, the next release build is the check) and its
    commented line from `gui/plot_canvas.py`; a test checks that nothing
    imports it. The first change of the dependencies of the package in this
    phase. No numbers touched.
  - PR 7 (trials with fresh agents, plan agreed with Guido: the client and
    the scripts committed, J0103-1305 raw): the minimal client
    (`benchmarks/mcp_trials/client.py`, a server kept alive on a Unix socket,
    `tools` and `call`; its log with the times of the calls and of the
    progress notifications; tested in `tests/test_mcp_trial_client.py`), the
    inputs (`make_data.py`: the full-chain mock without its true continuum and
    noiseless flux, J0103-1305 raw, the mock of the editing trial) and the
    prompts (`prompts.md`). Trials 1 (the mock, in the Docker image, built in
    the cloud session for the first time) and 2 (J0103-1305) run, trial 3 (the
    editing tools and undo) after the fixes of the server, in the next week of
    Guido's usage limit; report in section 4.3: F111-F116 found; heartbeat
    and cancellation needed, the tools not grouped; the real analysis moved to
    phase 7 (Guido). Numbers unchanged: nothing of the package touched, so no
    slow tests.
  - PR 7, trial 3 (the editing tools and undo, after PR 7-bis; local server,
    the mock of `make_data.py edit`): the agent reached the true model of the
    blend (two CIV and a SiIV with b tied) through undo, redo and the
    constraints, each step logged with its reason; report in section 4.3:
    F117-F120 found, F112 and F115 seen again; the GUI's Stop and
    `_is_stop_requested` found, to reuse for the cancellation. Numbers
    unchanged: nothing of the package touched, so no slow tests.
  - PR 7-bis (the short PR of the server decided in section 4.3, F113): a
    heartbeat (`progress.py`: after 20 s of silence of a recipe a progress
    notification "Still running <tool>: N s; last: <line>", its value a part
    of what is left before the next percent, so that it always increases;
    sent from its own thread, hence `anyio>=4.11` in the `[mcp]` extra and
    the image) and two locks instead of one (`sessions.py`): `recipe_lock`
    for the recipes and the tools that change the states (open, undo, redo,
    duplicate, close, `gui_open`), `state_lock`, held briefly, for the rest:
    the read-only tools, `session_save` and `session_note` answer while a
    recipe runs, on the state before it (said by `running` in
    `session_list` and `session_describe`, and by a warning of
    `session_save`). Tested with a recipe of the test that waits. Left out:
    the cancellation (with F111). Numbers unchanged: only the server
    touched, so no slow tests.
  - PR 7-ter (the outputs of trial 3, Guido's decisions of section 4.3;
    characterization tests first): F117 (`SystemListV2.update_constraint`
    checks every new link, also with `is_free` left out; the example of
    `update_constraint` is `p['<uuid>'].z`), F118 and F119 (a generic
    `details` field of the output of every recipe, the rest of
    `RecipeResult.extra` as JSON, no recipe named in the server;
    `optimize_hierarchy` takes an optional `trace` that the recipe reports as
    `optimization`; `identify_lines` keeps the candidates it scored and
    dropped in the metadata of the spectrum, `rejected_candidates_json`, and
    reports them as `dropped_candidates`), F120 (`GUI_COLUMN_NAMES`, accepted
    by `spectrum_read`), and the tolerance of `compare_sessions`
    (`SAME_PARAMETER_RTOL = 1e-8`, relative: a parameter listed as changed
    only beyond it; the methods text reads the same comparison, so a refit
    that changed nothing no longer adds "changed" components to it). The
    fits and the identification are not touched (only what they report; the
    snapshots unchanged), so no slow tests and no benchmark.
  - PR 7-quater (F111 and the cancellation, F113 closed; characterization
    tests first; Guido's choices on the diagnosis): the cause of F111
    (register): `metals_auto` refitted the session it started from and lost
    the ions found; the fallback on the whole spectrum was met by the ions
    outside the spectrum. Fixed: the ions found kept; `NoPixelsToFit` in the
    fitter instead of the whole spectrum; in `metals_auto` the ions outside
    the spectrum skipped with a WARNING, a progress line per trial, the AIC
    change of each ion tried (`tried`), the cancellation checked at each
    trial; the texts of the recipe (details, `aic_penalty`, methods) say
    what it does; F121 (Guido's choice): an ion is not kept at the bounds
    of the fit nor with an AIC change above `-aic_penalty`. The
    cancellation as in section 4.3: one signal for the GUI and the server,
    `session_cancel`, the asynchronous tools of the recipes. Numbers
    changed: `metals_auto` (the benchmark does not run it: a targeted test,
    the SiIV of the mock of trial 1 found), and the fits of components with
    no pixel (now not fitted, before at their first values on the whole
    spectrum); slow tests and benchmark run (results in the PR).
    Measured afterwards on J0103-1305 raw (a script outside the repository:
    `make_data.py j0103`, z_em = 4.07, R = 45000, continuum, `doublets_auto`
    with 65 components, `metals_auto` on the six default ions, with `run`):
    320 trials built, 81 not tried for their lines outside the spectrum, 239
    fitted in 6.4 min (26.7 min with the old criterion) and the final
    `refit_all` in 121 s (507.5 s the recipe); 82 ions kept (192 before, 67
    of them at b = 50 km/s; none now), 157 not kept, with the reasons of
    `details.tried`: b at its bound 108, logN below 12 29, AIC change above
    -5 7, b below 2 km/s 7, logN at its bound 6 (F121, same numbers).
  Found along the way: the startup hook of the cloud sessions created
  `.venv` without installing the package once (November 2026, at the
  start of this phase); installing again by hand worked. To watch: if
  it happens again, the hook needs a check that the import works.

**Result of the phase** (November 2026, PRs #726-#747 and the closing PR 8;
the release 2.6.0 follows). The full MCP server, with the principle of
the plan kept: a thin layer over the package, no tool names a recipe.
- Tools: 65, generated from the registry (the 48 recipes, none excluded:
  F106 made `export_ascii` one) and 17 of the server: the sessions
  (`session_open`, `save`, `list`, `describe`, `history`, `undo`, `redo`,
  `duplicate`, `close`, `note`, `methods`, `cancel`), the inspection
  (`components_list`, `spectrum_read`, `spectrum_plot`), `gui_open` and
  `gaps_summary`. A recipe added to the registry is a tool with no code of
  the server (a test adds fake ones). The pilot's five tools keep their names.
- Outputs from `RecipeResult` only (what changed, the components removed
  with their reason, `details` from `extra`, new sessions), the same for
  every recipe. What the tools compute is in the package, for the GUI and
  the scripts too: `SystemListV2.fit_groups`, `SpectrumV2.columns_in_range`,
  `plotting.py`, `recipes/progress.py`, `gaps.py`, `methods_text.py`,
  `follow.py`.
- Nothing is lost: every state kept (`session_undo`, `redo`, `duplicate`),
  `session_close` refuses to lose changes. The agent's decisions in the log
  (`reason`, notes, gaps, recorded also in a local file never sent anywhere),
  the log as the text of a methods section (also *File > Export Log as Text*).
- Long recipes: progress notifications, a heartbeat after 20 s of silence,
  read-only tools that answer during a recipe, and a cancellation shared
  by the GUI's Stop (`RecipeCancelled`).
- The GUI follows a session of the server, read only, with the editable copy
  (PR 5-bis); the headless Docker image, with its CI job (PR 6).
- Trials with fresh agents (section 4.3): the chain works at the first
  try; with undo, redo and the constraints an agent reached the true model
  of a blend, step by step in the log. 65 tools cost about 35 thousand tokens
  of context; not grouped (both agents advised against).
- Found: F106-F121. Fixed: F106, F108-F111, F113, F117-F120, F121 in part.
  Open: F107 (6-quater), F112, F114, F115, F116, and what goes to phase 7:
  the criteria of F118 and F119, the rest of F121.
- Also decided or fixed along the way: F109 (an archive saved again repeated
  its log: data loss, reported at once), `Set properties` no longer writes
  the object, `scienceplots` removed from the dependencies (F110),
  `compare_sessions` with a tolerance, `benchmarks/bench_gui.py` repaired.
- Numbers: the recipes unchanged but for two changes, both of PR 7-quater
  (`metals_auto`: keeps the ions found, skips the ions outside the spectrum
  and those at the bounds of the fit or with a weak AIC change; fits of
  components with no pixel: not fitted instead of fitted on the whole
  spectrum); the benchmark of phase 4 does not run `metals_auto`, a
  targeted test measures it.
- Tests and coverage (master after #747, `pytest -n auto --cov=astrocook`,
  slow tests skipped): 3601 tests collected, none failing, 5 skipped (slow),
  2 expected failures (F17, F28); coverage 87% (core 93%, recipes 94%,
  fitting 88%, GUI 86%, `mcp_server` 98%, `io` 51%: the loaders wait for
  sample files). `ruff check .` clean.
- Left: the real analysis of Guido (moved to phase 7: the first test of its
  components); the descriptions of the parameters shared by every tool shorter
  and a paragraph of typical workflows in the instructions of the server (the
  advice of the trials, not done: for the recipes' review in 6-quater); a
  measure of the fit for the agents (residuals, chi2, AIC and BIC on a fixed
  set of pixels, the gaps of trial 3); the chi2 of the groups in the result
  of every fit; a redo that cancels the gap `undone` of the state undone;
  the publication of the image on GHCR (after #549); the check that the
  macOS app builds without `scienceplots` (F110, the build of the release).
- Code quality (decided at the closing): ruff `B` (a PR of its own), the type
  checker and `SessionV2` out of `core/` moved to the start of 6-quater.
- Before release 2.6.0: the version and the CHANGELOG; the numerical changes
  above in its notes; the slow tests and the mock benchmark on the final
  master; `switcher.json` ("preferred") on master before publishing; the macOS
  app build (F110); the Docker image not published.
- Release 2.6.0 (10 November 2026, PR 9, #751; published by Guido on GitHub
  with the macOS app, not on PyPI). Checked on the final master (`e4c85a3`):
  the mock benchmark (`validate_mocks.py --compare`, full grid, 135 mocks,
  1176 s) identical to the baseline in all 19 rows (it does not run
  `metals_auto`: the targeted test measures it); the slow tests 5 of 5 (on
  the branch, before the merge; the CI job repeats them on master); the fast
  suite green and `ruff check .` clean. The macOS app built with
  `build_app.sh` (without `scienceplots`, F110) opens from the Finder and
  from the command line (Guido). The CHANGELOG of 2.6.0 had been missing the
  entries of PR 7-bis, 7-ter and 7-quater: added in PR 9. Found along the
  way: the CI test `test_the_read_only_tools_answer_while_a_recipe_runs`
  was flaky (it assumed the order of two tasks started together, and the
  recipe's own wait of 30 s ran while the reads of the test ran); fixed in
  #750, reproduced with 8 busy processes on one core.
- Seen after the release, not a defect: the banner of a followed session
  ("Followed session. Updated by an AI agent after every change. Read only.")
  is right when a server (even a forgotten one, here a 2.5.0) is still
  running and writes the archive; it does not say which process. Idea for
  6-quater: pid and version of the marker in the banner.

### Phase 6-quater – Review of the recipe parameters
Named after 2.5.0 (November 2026, Guido): the parameters and interfaces
of the recipes reviewed recipe by recipe, adapted to each procedure (as
the editing recipes already are; Guido's question in phase 6-bis, PR 6),
for the recipes that phase 7 does not rewrite (edit, flux, file, the
manual fits); those of the continuum and the identification go with
phase 7, with the collections of elements to fill some parameters. To be
decided in its plan: how the session logs and the scripts that use a
renamed parameter still replay (e.g. the old name kept as an alias in the
registry). Also F107 (the labels of the menus of the GUI, not the `brief` of
the registry; Guido, phase 6-ter). Plan to be agreed with Guido. Proposed (code review of 2.5.0,
to be agreed; see "Code quality" after phase 7): the `V2` suffixes of the
class names, with the same aliases; the complex functions of the
recipes reviewed.
Brought from 6-ter (release 2.6.0, November 2026): ruff `B` first (a PR of
its own), then the type checker and `SessionV2` out of `core/` with the
`V2` renaming; the descriptions of the parameters shared by every tool
shorter, and a paragraph of typical workflows in the instructions of the
server (the advice of the trials of section 4.3); F107 (the labels of the
menus from the registry); the names of the columns in one set (F120); the
banner of a followed session saying which server writes the archive (pid
and version of the marker). The plan itself is to be written at the start
of the phase, in a conversation of its own.

**Plan of the phase** (agreed with Guido at the start of the phase, 11
November 2026, in the conversation of PR 1; each PR in a conversation of
its own, except the smallest; the model of each PR proposed by the rule
"Choice of the model" in `CLAUDE.md`):
- *Starting point* (master after 2.6.0).
  - Recipes in the scope of the phase: file 2 (5 parameters), edit 15
    (42), flux 6 (14), the manual fits of `absorbers` (about 12, about 60
    parameters); continuum (5) and identification and the pipelines of
    `absorbers` go with phase 7. Examples of what the review is for:
    `rebin(norm: str = 'False')`, `equalize(reference_session: str =
    'None')`, booleans and None written as strings.
  - Ruff `B`: 33 warnings, as in the review of 2.5.0 (21 `B905`, 5
    `B007`, 3 `B008`, 2 `B023`, 2 `B904`).
  - Mypy (`--ignore-missing-imports --follow-imports=silent`, no other
    option): 149 errors in `recipes/` and `mcp_server/`, 336 in `core/`
    and `fitting/`.
  - `V2` suffixes: 21 classes; the most used `DataColumnV2` (~470
    occurrences), `SessionV2` (~370, in 67 files), `SpectrumV2` (~320).
- *PRs*, in this order (estimates in working time of the session, before
  Guido's review):

| PR | Content | Model | Estimate |
|---|---|---|---|
| 1 | Ruff `B` on, its 33 warnings fixed (with the plan of the phase and the rule of the model in `CLAUDE.md`) | Sonnet | 1-2 h |
| 2 | Aliases of renamed names: one mechanism in the registry for a renamed parameter (the old name accepted by the session logs, the scripts and the MCP tools), reused for the `V2` classes | Opus | half a day |
| 3 | Type checker, step 1: mypy in CI on `recipes/` and `mcp_server/` | Sonnet | half a day (done, 11 November 2026: see below) |
| 4 | Type checker, step 2: `core/` and `fitting/` | Sonnet | half a day to a day (done, 11 November 2026: see below) |
| 4-bis | Type checker, step 3: `gui/` (361 errors, 199 of them `attr-defined` from the mixins of `MainWindowV2` and `SpectrumPlotWidget`: the same cure as PR 4, a host class for the type checker) | Sonnet | half a day (done, 11 November 2026: see below) |
| 5 | `SessionV2` (and `SessionHistory`) out of `core/`, e.g. into `astrocook/session/`, the old import kept as an alias | Sonnet | 2-3 h (done, 11 November 2026: see below) |
| 6 | The `V2` suffixes dropped, with the aliases of PR 2 (code, tests, docs) | Sonnet | half a day (done, 11 November 2026: see below) |
| 7 | Review of the recipes of file and flux (8 recipes) | Opus | half a day (done, 11 November 2026: see below) |
| 8 | Review of edit, columns and expressions (`set_properties`, `x_convert`, `y_convert`, `apply_expression`, `mask_expression`, `smooth_column`), with F120 (one set of column names) | Opus | half a day to a day (done, 12 November 2026: see below) |
| 8-bis | The units of the flux and the CSV files (agreed with Guido on PR 8, 12 November 2026): a loader assigns a physical unit to the flux only when it reads it (`BUNIT`, `TUNIT` or a known descriptor of the format), otherwise an Astrocook unit of arbitrary flux (`def_unit`, registered so that the archives and the FITS files read it back); a way to declare the unit of the flux without converting it (`y_convert` converts); `export_ascii` writes the names of the GUI in ASCII (`lam`, `lammin`, `lammax`, `F`, `dF`, ...) and the units, read back by the CSV loader (`x`, `y`, `dy` and the old files still read), and `import_ascii` skips `lam` as it skips `x`; F129 (found on PR 8: the CSV files hold no units, the loader assumes nm and `erg/(Å s cm2)`, so a spectrum exported after `y_convert` or with the axis in km/s comes back with the wrong unit, and `cont` is lost). Today `astrocook_fits`, `simple_csv`, `ascii_telluric_model`, `xshooter_fits` and the V1 archives assign `erg/(Å s cm2)` always, PypeIt, ESPRESSO S2D and `ascii_resvel_header` dimensionless; the iso-magnitude grid is drawn for an uncalibrated spectrum labelled `erg/(Å s cm2)`. Values unchanged (only units), so no number of the continuum, the fits or the identification moves; sample files of the formats needed for the specific loaders | Opus | 3-4 h (done, 12 November 2026: see below) |
| 9 | Review of edit, several spectra (`delete`, `split`, `extract_preset`, `trim_common`, `stitch`, `equalize_and_stitch` with its complexity of 31, `coadd`, `import_systems`, `import_systems_ascii`) | Opus | 2-3 h (done, 12 November 2026, about 4 h with F132-F134: see below) |
| 10 | Review of the manual fits (`add_component`, `update_component(s)`, `delete_component`, `detect_anchor`, `fit_component`, `update_constraint(s)` with F117, `add_linked_system`, `optimize_system`, `fit_bayesian`); possibly split in two | Opus | 3-5 h (two PRs if split) |
| 11 | F107: the labels of the menus of the GUI from the registry (and the screenshots) | Sonnet | 1 h |
| 12 | MCP server: the descriptions of the shared parameters shorter, a paragraph of typical workflows in the instructions, the banner of a followed session with the pid and the version of the marker | Sonnet | 1 h |
| 13 | Closing of the phase (added on PR 8, Guido): "Result of the phase" in this section, the progress table, "Rules and decisions in force" (what of the phase still holds: the aliases, the one set of names of the columns, the units of the flux), `CLAUDE.md`; the concluded detail moved to `ROADMAP_HISTORY.md` | Sonnet | 1 h |
| 14 | Release 2.7.0 (added on PR 8, Guido): the version, the CHANGELOG (the summary of the phase and the numerical results: unchanged on the benchmark), the "preferred" entry of `switcher.json`; the release and the macOS app published by Guido | Sonnet | 1 h |

- *Result of PR 3* (type checker, step 1; 11 November 2026): `mypy`
  (`[tool.mypy]` in `pyproject.toml`, `files` = `recipes/` and
  `mcp_server/`, `ignore_missing_imports`, `follow_imports = "silent"`, no
  strict option) gives zero errors; mypy is in the `[dev]` extra and runs in
  a `types` job of the CI on the full install (the `[mcp]` extra must be
  there, or the server would be checked as `Any`). The 149 messages of the
  starting point (143 errors and 6 `annotation-unchecked` notes, which are
  not errors and stay): 100 `union-attr` and 19 `arg-type` (almost all
  from `SessionV2.spec` typed `Optional`: `spec_of(session)` in the new
  `recipes/_spectrum.py`, an `assert`; `column_of` for `get_column`);
  21 `assignment` (12 `x: T = None` made `Optional[T]` (the recipe
  methods only use their signature for `bind`, never the annotations),
  3 containers typed `dict[str, Any]`, the rest narrowed); 3 `var-annotated`
  (annotations); 3 `attr-defined` (`assert` on the constraints of the
  Bayesian fitter, a `type: ignore` on the marker `typed_recipe`); 2
  `return-value` (the annotation of `coadd` and of `fit_bayesian`, which
  return a session); 1 `operator` (`str(reason)`, the same text). No change
  of behaviour; no bug found, so nothing in the findings register. Three
  `type: ignore` wait for the annotations of `core/` (PR 4):
  `SessionV2.with_new_spectrum` (a session without spectrum can be
  copied), `SystemListV2.delete_component` (`str = None`) and
  `TableAdapterV2.__getitem__`. Only the fast suite was run (3620 passed):
  the corrections are of types, so the numbers cannot change.
- *Result of PR 4* (type checker, step 2; 11 November 2026): `mypy` gives
  zero errors on `recipes/`, `mcp_server/`, `core/` and `fitting/` (`files`
  of `[tool.mypy]`; options unchanged). The 336 errors of the starting
  point (re-measured: 336, in 15 files, 262 of them in the five mixins of
  the spectrum and in `voigt_fitter.py`): 212 `attr-defined`, 33
  `call-arg`, 25 `return-value`, 20 `assignment`, 15 `arg-type`, 9
  `union-attr`, 7 `var-annotated`, 4 `operator`, 4 `misc`, 3 `index`, 2
  `call-overload`, 1 `type-var`, 1 `syntax`. By cause:
  - Root cause 1, the mixins (`spectrum_*.py`): they used `self._data`,
    `x`, `y`, `t`, `history`... that `SpectrumV2` defines, so 150 of the
    `attr-defined` and all 33 `call-arg` (`self.__class__(data=...)`) were
    about `FluxMixin` and the others. New `core/_spectrum_host.py`: a
    `SpectrumHost` class with the declarations, seen by the type checker
    only (`HostBase` is `object` at run time, so nothing changes), and the
    methods that return a new spectrum annotated `-> Self` (25
    `return-value`). One mixin method used by another (`scale_flux`,
    `update_column`...) is declared there too.
  - Root cause 2, `voigt_fitter.py`: attributes started as `None` and
    filled by `prepare_fit_context` (`_x_tau`, `_fit_mask`...) are
    annotated `Any`; `SystemListV2.constraint_model` (created by the
    next line of `__init__`) is annotated as its class, with one
    `type: ignore[assignment]` on the `None`: this removed about 60
    `attr-defined`, `index` and `call-overload` in the fitters.
  - Wrong annotations corrected: `rebin_spectrum` (it returns 6 values,
    the last optional), `_find_kinematic_doublet_candidates` (3-tuples),
    `parsed_constraints` (`Dict[Tuple[str, str], ParameterConstraintV2]`),
    `load_session_from_file` and `SessionV2.open_new` (they return `0` on
    failure: `Union[SessionV2, int]`), `_get_ne_contexts` (a tuple written
    as `(A, B)`), `x: T = None` made `Optional[T]` (`SpectrumV2.__init__`,
    `extra_vars`, `delete_component`, `fitting_context`,
    `set_active_components`, which takes any collection, `with_new_spectrum`),
    `TableAdapterV2.__getitem__` (two `overload`s). The annotations are
    not read at run time (`dataclasses.fields` only reads names).
  - The rest: containers annotated (`aux_map`, `reliable_ids`, `tried`,
    `_clicks`...), a `cast` or an `assert` where a column or a result is
    known to be there (`get_column`, `y_col_new`, `_active_uuids` after a
    fit). Three `type: ignore` in the package: `current_state` of
    `SessionHistory` (returns `None` when there are no states, never
    reached: the history always holds its initial state), plus the one
    on `constraint_model` above and the marker of `typed_recipe` (PR 3).
    The three of PR 3 that waited for this PR are removed.
  No change of behaviour, and no bug in the findings register. One slip
  of PR 3 corrected: the comment of a `type: ignore` in
  `mcp_server/sessions.py` had swallowed the call
  `.with_new_system_list(state.systs)` of the copy of a session; the call
  is back as it was (it is redundant, since `with_new_spectrum` already
  carries the list, so the result was the same). Only the fast suite was
  run: the corrections are of types, so the numbers cannot change.
  `gui/` gives 361 errors in 23 files (199 `attr-defined`, 40
  `union-attr`, 36 `has-type`, 27 `assignment`, 20 `arg-type`, 17
  `call-overload`; 124 in `main_sessions.py`, 61 in `main_recipes.py`, 32
  in `main_scripts.py`): the same mixin pattern, so a PR of its own (4-bis)
  with the same cure.
- *Result of PR 4-bis* (type checker, step 3; 11 November 2026): `mypy`
  gives zero errors on the whole package (`files` adds `astrocook/gui`).
  The 361 errors of the starting point (re-measured, in 23 files), by
  cause:
  - The mixins, 196 errors (`attr-defined` and the `has-type` that came
    with them): as in PR 4, a host class for the type checker.
    `gui/_window_host.py` (`WindowHost`, a `QMainWindow` with the
    attributes and the methods that the seven mixins of `MainWindowV2` use
    of the window or of each other) and `gui/_plot_host.py` (`PlotHost`,
    for `PlotMouseMixin`); `HostBase` is `object` at run time, so nothing
    changes. The methods that other mixins define are declared with a
    loose signature `(*args, **kwargs) -> Any`: the declarations of
    attributes carry the types.
  - Short names of the Qt enums (`Qt.WaitCursor`, `QDialog.Rejected`...,
    35 `attr-defined`): written with their enum
    (`Qt.CursorShape.WaitCursor`), the same objects, found with the
    introspection of PySide6.
  - `Optional` narrowed with `cast`, not `assert` (the handlers of the GUI
    must raise the same exception as before): `active_history`,
    `SessionV2.spec` (about 45 `union-attr` and `arg-type`),
    `QApplication.instance()`, `QSettings.value` (typed `object`).
  - Containers annotated (15 `var-annotated`) and `x: T = None` made
    `Optional[T]` (7 `assignment`).
  - Five `type: ignore`, each with its reason: the two dialogs that keep
    their main layout in `self.layout` (it hides the method of QWidget;
    `layout: Any` on the class), `toolitems` of the toolbar (a `None` marks
    a separator, as in matplotlib), `RecipeDialog.schema = None` (the
    dialog stops at `if not self.schema`) and the decimation of
    `_PlotArrays` (norm and model may be `None`, which the function
    returns as it is).
  No change of behaviour; no bug found, so nothing in the findings
  register. Checks run: `mypy`, `ruff check .`, the fast suite (the GUI
  tests included); slow tests not run, the corrections are of types.
- *How a review PR (7-10) works*: the conversation opens with a table of
  the parameters of its recipes, as they are and as proposed (name,
  type, default, unit, description); Guido agrees it, then the change is
  made, with characterization tests first. No PR of the phase should
  change the numbers of the benchmark: a change of a numerical default
  is said before, and checked with `validate_mocks.py --compare`.
- *Total*: 12-13 PRs, about 6-9 working days of the sessions (phase
  6-bis: 19 PRs in about 3 days; 6-ter: about 15 in 5); 2-3 weeks of
  calendar, the bottleneck being one PR at a time and Guido's review of
  the tables of PRs 7-10.
  *Recalibrated on PR 8* (Guido, 12 November 2026): PRs 1-7 took one
  working day in all, so the estimates above were 4-8 times too long (the
  estimates of PRs 1-8 are left as written, as a record). The PRs left,
  with 8-bis, 13 (closing) and 14 (release 2.7.0) added on PR 8: about
  12-16 h of the sessions, two working days; the calendar is set by
  Guido's review of the tables of 8-bis, 9 and 10.
- *Open decisions*, to be taken in the PR concerned:
  - ~~an old name (parameter or class): a silent alias, or a
    `DeprecationWarning` with a version of removal (the tests treat
    warnings as errors) (PR 2)~~: decided in PR 2 (Guido, 11 November
    2026), below;
  - `refit_all` and `clean_negligible`: in the manual fits (PR 10) or in
    phase 7, which rewrites `refit_all`;
  - F101 (the `app.run_recipe` of the script console fails with a
    TypeError): in PR 10, or with the decision on the console;
  - PRs 3-6 before the reviews (the order proposed), or the renaming
    last, so that the reviews do not meet a change of every file.
- *Result of PR 5* (`SessionV2` out of `core/`; 11 November 2026): new
  package `astrocook/session/`, with `session.py` (`SessionV2`,
  `load_session_from_file`, the rest of the old `core/session.py`, moved with
  `git mv`) and `history.py` (`SessionHistory`, the old
  `core/session_manager.py`); `astrocook.session` re-exports the three
  names. The reason: a session imports the recipes, so it was a module of
  `core/` that depended on the layer above it (and `core/__init__.py`
  imported it at once). Aliases: `core/session.py`, `core/session_manager.py`
  and `core/__init__.py` (for `SessionV2`) hold only a module `__getattr__`
  (PEP 562) that imports the new module **when a name is asked** and returns
  the same object (`from astrocook.core.session import SessionV2` is
  `astrocook.session.session.SessionV2`). Lazy, so `core/` imports nothing
  of `astrocook.session` at import time and no cycle can form; and **without
  `DeprecationWarning`**, unlike the renamed classes of PR 6: the path is
  used by the scripts and by the tests, which have warnings as errors, and
  the plan of PR 2 (the helper `astrocook/_renamed.py`) was not written, so
  there was no mechanism to reuse: PR 6 decided to keep them apart (see
  its result). The names as strings: no pickle, and neither the
  archives nor the logs hold the path of a module (the logs hold recipe
  names and parameters); `SessionV2.__module__` is now
  `astrocook.session.session`. The imports of the package, the tests and the
  docs (`docs/tutorials/scripting.md`, `generate_sphinx.sh`, the docstring
  references) use the new path; `astrocook/session` is in the `files` of
  mypy; `tests/test_broad_excepts.py` follows the file. New
  `tests/test_session_package.py`: the old and the new path give the same
  objects, an unknown name is an `AttributeError`, no module of `core/`
  imports `astrocook.session` at module level. No change of behaviour. Fast
  suite (3623 passed), `ruff check .` and `mypy` clean; no slow test and no
  benchmark: only modules and imports moved, so the numbers cannot change.
- *Result of PR 6* (the `V2` suffixes; 11 November 2026): the 21 classes
  lost the suffix, in the package, the tests, the benchmarks, the docs and
  `CLAUDE.md` (a word-boundary `sed` over the files that held the names,
  then `ruff`, `mypy` and the suite). The names, old → new: `SessionV2`,
  `SessionMetadataV2`, `SpectrumV2`, `TableAdapterV2`, `SpectrumDataV2`,
  `DataColumnV2`, `ComponentDataV2`, `ParameterConstraintV2`,
  `SystemListV2`, `SystemListDataV2`, `LogEntryV2`, `LogNoteV2`,
  `HistoryLogV2`, `VoigtFitterV2`, `VoigtModelConstraintV2`,
  `MainWindowV2`, `RecipeFileV2`, `RecipeEditV2`, `RecipeFluxV2`,
  `RecipeContinuumV2`, `RecipeAbsorbersV2`, each without `V2`. Conflicts
  looked for and none found: no `Session`, `Spectrum`, `SystemList`,
  `MainWindow`, `VoigtFitter`... existed as a name in the package, the
  tests or the imported libraries (Qt has `QMainWindow`; the only
  `Spectrum` is the V1 class in `tests/reference/v1/`, a separate module,
  imported there as `SpectrumV1`). Left as they are, on purpose: the
  `"AstrocookV2"` of the `QSettings` (the key of the user's saved
  settings: renaming it would lose them), the `V2` of prose ("V2
  archive", the format `.acs2`), `ROADMAP_HISTORY.md` and the past entries
  of `CHANGELOG.md`, the V1 copies in `tests/reference/v1/`.
  - *The aliases* (one mechanism, `astrocook/_renamed.py`, written here as
    the decision of PR 2 asked): a module that defined or exported a class
    ends with `__getattr__ = renamed(__name__, {"SpectrumV2": "Spectrum"})`
    (PEP 562). The old name gives **the same class** (`isinstance`
    unchanged) with a `DeprecationWarning` that names the module, the old
    and the new name, pointing at the importer's line; the new name is
    looked up when the old one is asked, so the helper works while the
    module is still importing; a name not in the mapping is the usual
    `AttributeError`. Modules with a mapping: those that defined a class
    (`core/structures.py`, `core/spectrum.py`, `core/system_list.py`,
    `fitting/voigt_model.py`, `fitting/voigt_fitter.py`, the five modules
    of `recipes/`, `gui/main_window.py`, `session/session.py`) and the
    packages that exported one (`astrocook`, `astrocook.session`,
    `astrocook.core`). A class imported from a module that only imports it
    (e.g. `from astrocook.recipes.edit import SpectrumV2`) is not served:
    it was never the home of the name. No version of removal (3.0 at the
    earliest).
  - *The mechanism of PR 5 is not folded into it* (the question PR 5
    left): the two answer different things. A moved module (PR 5: the
    *path* `astrocook.core.session`, the name unchanged) stays silent
    (the scripts and the tests use it, and a path costs nothing); a
    renamed class (PR 6) warns, as decided in PR 2. They compose: the
    forwarders of `core/session.py` serve what the new module serves, so
    `from astrocook.core.session import SessionV2` warns (once, from the
    helper). `core/__init__.py` keeps its lazy `Session` (no import of
    `astrocook.session` at import time) and passes the old names to the
    helper.
  - *Names as strings*: the `.acs2` archives (the three of `tests/data/`
    opened member by member) and the session logs hold no class name (the
    logs hold recipe names and parameters; there is no pickle), so nothing
    to translate and no change of the format. The registry, `__name__`
    and the messages: only the docstrings, comments and log lines
    (`logger.debug("Initialized VoigtModelConstraint ...")`) changed. No
    `mock.patch` of a class by name, and no docs autogenerated from the
    names (`docs/conf.py` only had a comment). A stale comment of
    `structures.py` named `ConstraintDataV2`, which never existed: now
    `ParameterConstraint`. `tests/test_gui_structure.py` follows the new
    names.
  - New `tests/test_renamed_classes.py`: for each of the 21 classes (and
    each module or package that exported it) the old name is the same
    object as the new one and warns with `pytest.warns`, the new name does
    not warn (`simplefilter("error")`), the warning points at the
    importer, `isinstance` with the old name, an unknown name is an
    `AttributeError`, and the helper on a throw-away module.
  - *One snapshot regenerated*: `tests/snapshots/v1_archive_v1_metals_constraints_systems.npz`
    stored `str(ParameterConstraintV2(...))`, the `repr` of the dataclass,
    in its `constraint_values`; the `repr` carries the class name, so the
    stored strings now read `ParameterConstraint(...)`. Checked: of the
    keys of the file only that one differs, and by the `V2` and nothing
    else.
  No change of behaviour. `ruff check .` and `mypy` clean, the fast suite
  green (after the snapshot above); no slow test and no benchmark: only names changed, so the numbers
  cannot change.
- *Result of PR 7* (review of file and flux; 11 November 2026; the table
  of the parameters agreed with Guido in the conversation of the PR): of
  the 20 parameters of the 8 recipes (and `norm` of `rebin`, out of the
  registry), one renamed, two defaults changed, the others kept; no result
  of a valid call changes, so no benchmark and no slow test.
  - *Renamed*: `folder` of `file.export_ascii` (it was `path`, which an
    agent could take for the path of a file, as in `session_save`); `path`
    is its alias (PR 2), the first of the package (the table of the
    renamed parameters above, `RENAMED` in `tests/test_recipe_aliases.py`).
    The export is not logged, so no log holds the old name; scripts and
    MCP calls with `path` give the same files.
  - *Defaults*: `reference_session` of `flux.equalize` and `target_session`
    of `flux.resample`, the string `'None'` → `""`: required, as
    `sibling_name` of `detect_anchor` (the MCP schema lists them in
    `required`; the dialog starts with an empty field and the buttons of
    the sessions). Not `Optional[str]`: without a session the recipe cannot
    run, so None is not a value. `'None'` and `''` still fail as before
    ("No reference session is selected"); a log always holds a name (a
    recipe that fails is not logged). `norm` of `rebin` (out of the
    registry, #560) is annotated `Union[bool, str] = False` and warns that
    it is ignored for any true spelling (`'true'`, `1`: only `'True'`
    before); `'False'` and `'True'` of the V1 logs as before.
  - *Refused before the recipe runs* (`InvalidParameter`; each gave an
    empty or hidden result without a word, characterized before the
    change): `rebin` with `xstart` not below `xend` (a spectrum without
    pixels, F124); `equalize` with `order` below -1 (`min: -1`; -2, -3...
    were a spline, as -1); `calculate_running_std` with `output_col` empty
    or one of the columns of every spectrum (`CORE_COLUMNS` of
    `recipes/flux.py`: the result went to a hidden auxiliary column, F122);
    `calibrate_from_magnitudes` with magnitudes that are not
    `filter=magnitude` pairs or name an unknown filter (a `RecipeFailed` of
    the core before; checked by `_check_magnitudes`, the same reading as
    the core). A log with one of these values would no longer replay: none
    can hold a useful one.
  - *Kept a string*: `magnitudes` (a structured type would need a new type
    of the registry, an `object` in JSON, a widget and a parser, for one
    recipe; one band is the common case, and the string is that of the
    dialog, the log and the methods text). Its description lists the
    filters (from `STANDARD_FILTERS`) and says AB magnitudes.
  - *Descriptions*: `output_col` and `window_pix` of
    `calculate_running_std` (the rule of the window: an even one is
    increased by 1, below 3 is 3), `xend` and `xunit` of `rebin` (as
    `coadd`), the `details` of `export_ascii` (the series and the redshift
    always written; the items not found skipped) and of `import_ascii` (the
    column `x` not read); the docstrings aligned. `export_ascii` names the
    items not found in the session in a WARNING (e.g. `cont` without a
    continuum, the default, or a misspelt `systems`), the outcome
    unchanged. The methods text of `calculate_running_std` gives the window
    used (21 for 20), not the one passed.
  - *Left to later PRs*: the names of the columns (`F`, `dF` of `targets`
    against `y` of `input_col`, F120) to PR 8; the meaning of `-1` in
    `coadd` (F123) and the sessions passed as objects in scripts to PR 9
    (Guido: `equalize` and `resample` read another session only from the
    GUI or the MCP server, `stitch` also takes a `Session`; one rule for the
    parameters of kind `session` there); the *Help* button (F125) to PR 11.
    F30 (the edges of `running_std`, used by `calculate_running_std`)
    stays open, numbers unchanged.
  - *Docs*: the tutorial of Flux (the rebinning is a weighted mean, not a
    conservation of the flux; the sessions are clicked in a list, not a
    dropdown; the limits above), *Getting started* (export and import),
    *Using Astrocook from an AI Agent* (`folder`, the sessions required),
    the CHANGELOG. No screenshot: the dialogs changed (equalize, resample,
    the descriptions) are in no scene.
  - New `tests/test_recipe_review_file_flux.py` (the behaviour kept and
    the behaviour changed, each test saying what it was), new cases in
    `test_recipe_inputs.py`; `test_mcp_schemas.py` follows `folder`.
- *Result of PR 9* (review of edit on several spectra; 12 November 2026;
  the table of the 27 parameters of the 9 recipes and the decisions agreed
  with Guido in the conversation of the PR): one parameter renamed, two
  defaults removed, one type changed; the numbers of a valid call with the
  sessions in one unit do not move (fingerprints taken on master before the
  change, `test_valid_results_unchanged`), those with units mixed are fixed
  (F132-F134, agreed); no recipe of the continuum, the fits or the
  identification is touched, so no benchmark and no slow test.
  - *One rule for the sessions* (left by PR 7; `recipes/_sessions.py`): a
    parameter of kind `session` (one) or `sessions` (several,
    comma-separated; a new kind of the registry, declared to the agents as
    `x-astrocook-kind`, not used by the dialogs) takes names of open
    sessions, the GUI's list or the store of the MCP server; the name of
    the session the recipe runs on is that session also without a GUI; in
    a script also the sessions themselves (`Session` objects), and
    `resolve_params` writes their names in the log. A name not found fails
    naming it, with "In a script, pass the session itself." without a GUI
    (F135, the crash, and F136). Used by `import_systems`, `trim_common`,
    `stitch`, `equalize_and_stitch`, `coadd`, and by `equalize`,
    `resample` (flux) and `detect_anchor` (absorbers), which had the same
    lookup. `stitch.other_sessions` from `list` to `str` of kind
    `sessions` (a list of names is joined; the old logs, a list or its
    string, replay with the same result); `trim_common.others_names`
    renamed `other_sessions` (alias), as in `stitch` and
    `equalize_and_stitch`.
  - *Units* (Guido's rule, F132-F134): the sessions are combined in the
    units of a reference (this session for `stitch`, the bluest arm for
    `equalize_and_stitch`, the first named for `coadd`), the others
    converted to its axis and flux before they are scaled (nothing done
    when the units are the same, hence the same bits); `arb_flux` with a
    physical flux only when scaled (it takes the unit of the reference), a
    normalized flux never with one that is not; an axis in velocity against
    a wavelength refused (`x_convert` first). F138 (`flux.equalize`, the
    same in `equalize_to_reference`) registered, not changed: outside the
    recipes of the PR.
  - *F123* (option A): `equalize_order` of `coadd` `Optional[int]`, None
    no scaling, 0-5 (the maximum measured: `polyfit` is poorly conditioned
    from order 6-7 on overlaps of 20-50 nm, the reason written next to the
    parameter); -1 read as None, silently and for good (the old logs).
  - *Refused before the recipe runs* (`InvalidParameter`, each
    characterized before): `coadd` with a list without this session
    (F137), `xstart` not below `xend` (0 pixels, as F124), `xunit` neither
    a velocity nor a wavelength, `equalize_order` below -1 or above 5;
    `equalize_and_stitch` with cuts not increasing or an arm without pixels
    between its cuts, factors not positive (0 zeroed the flux), a range
    with its start not below its end (all "Split expression resulted in an
    empty spectrum" before, or nothing); `delete` with nothing, a column
    of every spectrum or a column not in the spectrum (a WARNING and a
    state recorded as done); `split` with an empty expression; 
    `extract_preset` with an unknown region (skipped with a WARNING among
    valid ones); `trim_common` with an unknown transition or not one
    transition per session.
  - *Defaults*: `expression` of `split` removed (it was
    `(λ > 400) & (λ < 500)`, an arbitrary cut for an agent that omitted
    it); `source_session` of `import_systems` `'None'` → `""` (required,
    as `reference_session` in PR 7; `'None'` fails as before). In the
    registry a default `""` is "required": the method has none, the dialog
    starts empty, the MCP schema lists it in `required`.
  - *Complexity*: `equalize_and_stitch` 23 (mccabe; 31 radon) and `stitch`
    12 → at most 10 for each function of `recipes/edit.py`, with helpers
    (`_numbers`, `_ranges`, `_blue_to_red`, `_arm_part`, `_equalize_arm`,
    `_in_units_of`), covered by the fingerprints.
  - *Kept*: `append`, `refit` of the imports (`refit_all` with 20 km/s and
    `max_nfev` 100, fixed), `sort` of `stitch`, `alias_map` of `split`
    (GUI only), `'auto'` of `equalize_ranges` and `manual_factors`,
    `window_kms`, `kappa`, `dx`.
  - *Docs*: the tutorial of Edit (the sessions by name or as objects,
    the units of the combination, the limits above), *Using Astrocook from
    an AI Agent* (the tools that read several sessions), the CHANGELOG. No
    screenshot: no dialog of a scene changed (the coadd dialog of the
    scenes is filled by the GUI with this session first).
  - New `tests/test_recipe_review_edit_sessions.py` (the behaviour kept and
    the behaviour changed, each test saying what it was); the batteries of
    `test_recipe_inputs_edit.py` and `test_recipes_edit_flux.py` follow
    (each case changed says what it was).
- *Result of PR 8-bis* (the units of the flux and the CSV files, with F129
  and F130; 12 November 2026; the analysis and the decisions agreed with
  Guido in the conversation of the PR): no value of the flux changes, only
  units, except F130 (the continuum converted with the flux) and the old
  CSV files read by name (F129); F131 changes a unit, not a value; the continuum, the fits and the
  identification do not read the unit, so no benchmark and no slow test.
  - *`arb_flux`*: `core/units.py` (`ARBITRARY_FLUX`, enabled in astropy
    when `core/structures.py` is imported; `parse_flux_unit`, `known_unit`,
    `FLUX_AUX_COLUMNS` moved there from `spectrum_flux.py`). A unit of its
    own, not `represents=dimensionless` (tried: it is `== dimensionless`,
    so it cannot be told from a normalized flux, and the FITS writer drops
    it without a warning). Not a FITS unit: `_write_fits_table` writes it
    in `TUNITn` by hand (astropy dropped it with a warning), the archives
    are read with `unit_parse_strict='silent'` and `known_unit`; Astrocook
    2.6.0 reads it with a `UnitsWarning`.
  - *Loaders*: physical unit only when read (`_flux_unit`): `astrocook_fits`
    (`TUNIT` of `y`), `xshooter_fits` (`TUNIT` of `FLUX`, arbitrary with
    `FLUXCAL = 'UNCALIBRATED'`), `generic_spectrum` (`TUNIT` of the flux),
    `espresso_s2d_fits` (`BUNIT` / Å); arbitrary: PypeIt (the descriptor
    `FLUXED` / `OPT_FLAM` waits for a sample file, Guido), `simple_csv`
    without units, `ascii_telluric_model`, `ascii_resvel_header`, the V1
    archives (no `TUNIT`; Guido: arbitrary; the snapshots of
    `test_v1_archive` regenerated, the strings of the units only). The
    transmission of `pypeit_tellmodel` stays dimensionless, its object
    model is arbitrary. Archives: a flux without a unit is read as
    arbitrary, a flux with a unit keeps it (Guido asked whether the
    `erg/(Å s cm2)` of the archives should become arbitrary too, as most
    were given by default: not done, an archive cannot tell a default from
    a calibration, e.g. `calibrate_from_magnitudes` or an X-Shooter
    product; a few lines in `_v2_table_to_spectrum_data` if wanted).
  - *Declared*: `yunit` of `set_properties` (`_current_`; Guido's choice
    over a recipe of its own): a flux density per unit wavelength or
    `arb_flux`, given to `y`, `dy`, `cont`, `model`, `cont_pl` and to the
    columns in the unit of the flux, without converting
    (`Spectrum.with_flux_unit`). `y_convert` refuses an arbitrary flux and
    names `set_properties` and the calibration; `calibrate_from_magnitudes`
    reads it as `erg/(Å s cm2)`, as a dimensionless one (it failed).
  - *GUI*: the axis and the tooltip say "arbitrary units"; the check box of
    the iso-magnitude grid is disabled, with a tooltip, unless the flux is
    per unit wavelength (it drew the grid for the uncalibrated spectra
    labelled `erg/(Å s cm2)`). The screenshot of *Set Properties*
    (`editing_set_properties`) rebuilt, with the new field; the spectra
    of the scenes are in `erg/(Å s cm2)` and keep it, so no other scene
    changes.
  - *CSV (F129)*: `export_ascii` writes ECSV with `,` (units, z_em in the
    metadata), the names `lam`, `lammin`, `lammax`, `F`, `dF`; the CSV
    loader reads it by name with all its columns; a file with a header
    `x`, `y`, `dy` (the old export) by name too, with the flux arbitrary;
    other files by position as before. `import_ascii` translates the names
    (`internal_column`), skips the axis, keeps the units of the file and
    refuses `F`/`dF` in another unit than the flux.
  - *Found*: F130 (fixed), F131 (`calibrate` labelled per Å a flux
    calibrated in its own unit, 2.5 mag brighter than asked for in nm;
    fixed, Guido: the unit of a physical flux is kept).
  - New `tests/test_flux_units.py` (characterization first, then the
    change, each test saying what it was).
- *Result of PR 8* (review of edit on columns and expressions, with F120;
  12 November 2026; the table of the parameters and the set of names
  agreed with Guido in the conversation of the PR): of the 17 parameters
  of the 6 recipes (and `input_col`, `output_col` of
  `calculate_running_std`, `targets` of `export_ascii`, left by PR 7) none
  renamed, two defaults removed, one changed; no result of a valid call
  changes, so no benchmark and no slow test.
  - *F120, one set of names*: those of the GUI (`λ`, `λmin`, `λmax`, `F`,
    `dF`, the other columns by their own name), wherever a column is
    named: the values of the parameters of kind `column`, the expressions,
    the defaults, the descriptions, the messages, the MCP tools, the docs.
    Chosen because the GUI, the tutorials and the expressions that the GUI
    writes in the logs (the regions, `(λ > ...) & (λ < ...)`) already used
    them; `λ` was chosen for its look (Guido), so `lam`, `lammin`,
    `lammax` are its ASCII spelling (not `lambda`, a keyword of Python:
    `df.lambda` fails in pandas; Guido). The names of the data (`x`, ...,
    `dy`) stay those of the archives, of the CSV files of `export_ascii`
    and of the attributes `Spectrum.x`, `.y`, `.dy`, and are accepted
    everywhere as silent and permanent aliases: the logs, the saved
    expressions and the scripts are not translated (the values are not
    parameters, so `ParamSpec.aliases` does not apply: the translation is
    where a recipe reads a column). The names of the parameters that name
    an axis keep `x` and `y` (`xunit`, `xstart`, `x_min` of the tools,
    `x_convert`). One module, `core/column_names.py` (`internal_column`,
    `internal_expression`, `display_column`, `is_core_column`,
    `is_axis_column`), replaces the five copies of the mapping (edit,
    `export_ascii`, `columns_in_range`, the recipe dialog, the data
    inspector); `GUI_COLUMN_NAMES` is still importable from
    `core/spectrum.py`. Found there: `calculate_running_std` did not find
    `F` (the button of its dialog) and wrote a column `F` that the
    expressions read as `y`. MCP: `spectrum_read` reads `F`, `dF`, `cont`,
    `model` by default and gives the axis as `λ`, or under the name asked
    for it (`columns_in_range` returns `axis`); `session_describe` lists
    the columns by these names; a paragraph of the instructions says it.
  - *Defaults*: `target_col` and `expression` of `apply_expression`
    removed (required, `""`; they were `F` and `F / dF`, which replaced
    the flux with its S/N from the GUI or the MCP server; the method had
    no default already); `input_col` of `calculate_running_std` `'y'` →
    `'F'` (the same column).
  - *Refused before the recipe runs* (`InvalidParameter`, characterized
    before the change): a target of `apply_expression` that is not an
    identifier (`''`, `'a b'` made columns no expression could name); the
    axis as the target of `apply_expression`, `mask_expression`,
    `smooth_column` (F127); a column to mask or smooth that does not exist
    (the message lists the columns; it blamed the expression); `name` of
    `set_properties` empty (ignored without a word); `xunit` neither a
    wavelength nor a velocity, `yunit` not convertible from the unit of the
    flux (both failed in the conversion, the same outcome); `output_col` of
    `calculate_running_std` with a name of the GUI of a core column.
  - *Kept*: `object` of `set_properties` (it hides the builtin in the
    method, but it is clear to a user and an agent: not worth an alias),
    the `_current_` of `z_em`, `resol` and `z_ref`, the defaults of
    `mask_expression` (`F`, `dF <= 0`: a sensible cleaning) and of
    `smooth_column` (`dF`, 100 km/s), `alias_map` (out of the registry, the
    GUI only, as `norm` of `rebin`).
  - *Descriptions*: the zero of the velocity axis of `x_convert` (the
    Lyman-alpha emission at 121.567 nm (1 + z_ref)), `yunit` per unit
    wavelength, `resol` in pixels with the axis in wavelength, the
    functions of NumExpr; the message of a function NumExpr does not have
    (F128).
  - *Left*: F84 (`x_convert` does not keep `z_ref`; a round trip with
    z_em changed in between gives shifted wavelengths, checked again here)
    to phase 7 (Guido: for order), documented as a warning in the
    tutorial of Edit; F126 (the unit of a new column) and F128 (the
    functions of the expressions) registered, not changed.
  - *Docs*: the tutorial of Edit (the names of the columns, the functions,
    the axis not a target, a new section on the conversion of the axes
    with the warning of F84), of Flux (`F`), *Using Astrocook from an AI
    Agent* (`spectrum_read`), the CHANGELOG. No screenshot changed: the
    scenes give their values (`flux_running_std` now passes `F`, the same
    image).
  - New `tests/test_recipe_review_edit_columns.py` (the behaviour kept and
    the behaviour changed, each test saying what it was).
- *Decisions of PR 2* (old names; proposed in the conversation of PR 2,
  agreed with Guido, 11 November 2026):
  - *A renamed parameter*: its old names are `aliases` of the parameter in
    the registry (a key of the parameter, `ParamSpec.aliases`). An alias is
    **silent and permanent**: no warning, no version of removal, because the
    archives are nodes not changed once saved (F109) and their logs must
    replay; an alias costs one line. A DEBUG line says the translation.
  - *One translation*, `RecipeSpec.canonical` (the old names replaced by
    the new one, in their place), called where the parameters enter a
    recipe (`typed_recipe`, before the signature is bound: scripts, Run
    All, the GUI, the MCP tools), where the log records them
    (`resolve_params`: every log written has the new names), by the MCP
    server before it refuses the unknown parameters (an old name taken,
    never listed), and when the log of an archive is read
    (`canonical_params`, in `load_session_from_file`): the methods text, the
    Log Scripter and the log of the MCP server show the new names; the
    archive on disk is not changed (and is never saved again, F109). The
    dialogs and the schemas of the tools are built from the parameters,
    and show the new name only.
  - *The old and the new name together*: `InvalidParameter`, even with the
    same value; in the log of an archive (only a log written by hand) the
    entry is kept as read, and the replay says what is wrong.
  - *Checks*: `build_spec` refuses an alias that is the name or the alias
    of another parameter; the tests check that no alias is in the
    signature of the method nor a parameter of the server. The tests rename
    `sigma_kms` of `flux.smooth` from `sigma_old` (fixture
    `renamed_parameter`, `tests/conftest.py`), no recipe of the package:
    the script, Run All, the archive, the methods text, the Log Scripter,
    the dialog and the MCP tool (`tests/test_recipe_aliases.py`,
    `test_gui_smoke.py`, `test_mcp_tools.py`).
  - *The classes without `V2` (PR 6)*: the same principle in a parallel
    form, since the classes are not in the registry: a mapping old name →
    new name declared in the module that exported the old name, read by a
    module `__getattr__` (PEP 562) made by one small helper (e.g.
    `astrocook/_renamed.py`, `__getattr__ = renamed(__name__, {"SpectrumV2":
    "Spectrum"})`), which gives the same class (`isinstance` unchanged) with
    a **`DeprecationWarning`** (Guido): the old names are used by the
    scripts only, never by the archives (to be checked in PR 6); no version
    of removal fixed (3.0 at the earliest). The tests use the new names, one
    test checks each alias with `pytest.warns`. Not written in PR 2.

### Phase 7 – Robust components: emission model and identification
Agreed with Guido at the start of the MCP pilot (October 2026). The goal of section 5, a complete fit of a spectrum in a single
operation, rests on its components: today the continuum works "well
enough" but is not robust, gives no uncertainty and fails in some cases
(DLAs, among others); the identification is fairly accurate, but its
precision and recall are not measured and it cannot account for *all*
the absorption. Both are redesigned before the full loop of section 5 is
built:
- **Emission model:** continuum (and emission lines) with an
  uncertainty, robust to strong absorbers (DLAs, BALs) and to the forest
  (directions in section 2, "Continuum");
- **Absorption identification:** every absorbed region accounted for,
  identified or explicitly unidentified, with alternatives and a
  credibility (directions in section 2, "System identification").
After #465 the automatic continuum is an iterative clipping with three measured patches (masked running mean in km/s, rejection of upward spikes with two thresholds, local S/N); Guido would like a simpler approach with a statistical basis, to replace it, measured on the same benchmark (directions: spikes recognised by their width, narrower than the resolution element; the continuum from the upper envelope of the flux or from the mean absorption expected at each redshift; the 0.055 sigma bias of one-sided clipping, F60). These are algorithmic choices, possibly complex, and they are Guido's:
for each component a design document first (options, costs, what the
literature does), then a prototype measured on the benchmark of phase
4 (validation on mocks). A new algorithm sits next to the old one, selectable, until the
benchmark shows it better; then it becomes the default, with the change
stated in the release notes. #465 (phase 4) is a robustness fix of the
current continuum, not a replacement. This phase is the first step of the
path of section 5, and a possible basis for a paper.
Proposed (code review of 2.5.0, to be agreed): a ceiling on the
complexity of the code it rewrites (see "Code quality" below).

### Code quality (review of 2.5.0)
A review of the code of 2.5.0 (November 2026, after the question of how
the package compares with similar codes), with measurements, and where
each proposal could go. **Proposals**; the three assigned to 6-ter
were moved to the start of 6-quater at its closing (Guido); the others are
not in a plan yet.

*Measurements* (2.5.0; the package has not changed since on master):
- Coverage 86% (`pytest -n 4 --cov=astrocook`, slow tests skipped).
  Lowest: `io/loaders.py` 17% (waits for sample files),
  `gui/corner_plot.py` 0%, `gui/main_window.py` 61%,
  `gui/system_inspector.py` 69%.
- Cyclomatic complexity (`radon cc`): average 5.1 over 1000 functions,
  methods and classes; 31 above 20, the highest
  `SystemTableModel.data` (68, `gui/inspector_table.py`),
  `RecipeAbsorbersV2.refit_all` (58, 293 lines),
  `PlotMouseMixin.show_data_tooltip` (58),
  `MainPanelsMixin._update_ui_state` (56),
  `MainRecipesMixin._process_recipe_result` (54),
  `VelocityPlotWidget._update_plot` (52),
  `MainRecipesMixin._on_recipe_requested` (48),
  `LinesMixin.identify_lines` (43, 307 lines),
  `VoigtFitterV2.prepare_fit_context` (42). 25 functions are longer
  than 150 lines. Phase 5 split the files, not the functions.
- Types: 45% of the functions annotate their return type; no type
  checker in CI.
- Ruff runs the pyflakes rules only. With flake8-bugbear (`B`): 33
  warnings, 21 of them `zip()` without `strict`; the two of `B023`
  (a function defined in a loop, `voigt_fitter.py` and `io/adapter.py`)
  checked and harmless (the function is called in the same iteration).
- Layers: `core/session.py` imports `recipes` and `io`, so the core is
  not the bottom layer (no cycle today; `SessionV2` is imported from
  `astrocook.core` in 72 places).
- No bare `except`, no `print`, no TODO left; `core`, `fitting`,
  `recipes`, `io` and `mcp_server` import neither Qt nor matplotlib.

*Proposals and where they could go:*

| Proposal | Phase | Why there |
|---|---|---|
| A type checker (mypy or pyright) in CI, gradually: `recipes/` and `mcp_server/` first, then `core/` and `fitting/`; the GUI last, if ever | 6-quater (moved from 6-ter, PR 8) | the server builds the tools from the typed parameters: checked types serve the agents first |
| `SessionV2` (and `SessionHistory`) out of `core/`, e.g. into `astrocook/session/`, the old import kept as an alias | 6-quater, PR 5 (done, 11 November 2026) | the headless server works on the sessions; the core becomes the bottom layer |
| Ruff with `B` (flake8-bugbear), a PR of its own; the style rules (`E`, `W`) left out | first PR of 6-quater (moved from 6-ter, PR 8) | small and mechanical; finds bugs, not style |
| The `V2` suffixes of the class names (`SpectrumV2`, `RecipeAbsorbersV2`, ...) dropped, the old names kept as aliases | 6-quater | changes the interface of the scripts, as a renamed parameter: the same mechanism of aliases |
| The complex functions of the recipes that phase 7 does not rewrite (e.g. `RecipeEditV2.equalize_and_stitch`, 31) | 6-quater | the recipes are reviewed one by one there |
| `identify_lines`, `refit_all` and the rest of the continuum and identification: a ceiling on the complexity of the new code (e.g. no function above 20) among the criteria of the new algorithms | 7 | phase 7 rewrites them: cleaning them before would be work thrown away |
| The complex functions of the GUI (the list above) and the mixins of `MainWindowV2`, which share their state through `self` | a phase after 7 (name to Guido) | phase 7 changes the GUI of the identification (alternatives, credibility) |
| Coverage of `io/loaders.py` | when the sample files arrive | waits for data from Guido |
| A test of `gui/corner_plot.py` | any time | one small PR |

### Releases
So far, releases have been distributed as a macOS app through the GitHub
releases page (not the Mac App Store); Astrocook is not on PyPI yet.
Sequence agreed with Guido (October 2026): phase 2 → 2.0.3 → phase 2-bis
→ phase 3 → phase 3-bis → 2.1.0 → MCP pilot → phase 4 → 2.2.0 → phase 5 →
phase 6 → full MCP server (section 4). Revised after phase 3: 2.1.0 is
published on GitHub (with the macOS app) without waiting for PyPI.
Revised at the start of the MCP pilot (validation on mocks in phase 4,
phases 6-bis and 7 agreed): phase 4 (validation first) →
2.2.0 → phase 5 → phase 6 → phase 6-bis → full MCP server → phase 7.
Revised at the start of phase 5 (October 2026): phase 5 → phase 5-bis →
2.3.0 → phase 6 → phase 6-bis → full MCP server → phase 7.
Revised at the end of phase 6 (November 2026, Guido): phase 6 → 2.4.0 →
phase 6-bis → full MCP server → phase 7 (the fix of a data loss, F94, and
the changed behaviour on failures, published before 6-bis changes the
interface of the recipes).
Revised at the end of phase 6-bis (November 2026, Guido): phase 6-bis →
2.5.0 → full MCP server → phase 7 (the new interface of the recipes,
published stable before the full server is built on it).
Named after 2.5.0 (November 2026, Guido): each step a phase with its own
name, in the order in which they are done: phase 6-ter – full MCP server;
phase 6-quater – review of the recipe parameters (the recipes that phase 7
does not rewrite); phase 7 – robust components.
**LICENSE and PyPI** (agreed with Guido, October 2026): an independent
step, not tied to a phase. The LICENSE (#549) and the first publication on
PyPI go into the first release after INAF answers: 2.1.1 if the answer
comes before phase 4, otherwise 2.2.0 or whatever the current version is
then. Until then Astrocook installs with pip from GitHub.
Documentation (agreed with Guido, October 2026): the root of the website
shows the docs of the latest release, `/dev/` those of master (see
`.github/workflows/docs.yml`). Before publishing a release, update the
"preferred" entry of `switcher.json` on master to the new version.
Steps of a release (as for 2.0.2 and 2.0.3): a release PR on master with
the version in `pyproject.toml`, `astrocook/__init__.py` and `README.md`,
the "preferred" entry of `switcher.json`, and the `CHANGELOG.md` entry
(the "Unreleased" section renamed, with a summary of the numerical
changes); after the merge, the macOS app (`build_app.sh`, on a Mac) and
the GitHub release with tag `v<version>` and the app attached, which
publishes the docs of the tag at the root of the website.
- **2.0.2 (corrective), released after phase 1-bis (October 2026):** the
  fixes of phase 1-bis (F1, F2, F10 minimal, F11, F16, F19, F20, F23-F27,
  F29, F32, #534) and issues #480, #442, together with the earlier fixes
  (ASCII export #476/#481, warning on saving to V1 `.acs`). #465 and the
  rest of F10 were deferred (#465, #526). Release notes in `CHANGELOG.md`
  state which numerical results change.
- **2.0.3 (packaging), after phase 2 (October 2026):** the `astrocook`
  command, installation with pip from GitHub, images in the package,
  macOS build fixes, safe unpacking of archives (F35), automatic docs
  deployment, cleanup; no numerical result changes. No license
  yet (#549).
- **2.1.0 (first consolidated release), after phases 3 and 3-bis (October
  2026):** V2 no
  longer depends on `legacy/` (removed), native reader of the V1 archives,
  `.acs` saving removed, smaller wheel, the fixes of phase 3-bis. Published
  on GitHub with the macOS app, as before; not on PyPI (no LICENSE yet,
  #549).
- **2.1.1 (only if INAF answers before phase 4):** LICENSE and first
  publication on PyPI, no other changes.
- **2.2.0, after phase 4 (October 2026):** the scientific-correctness
  fixes (#526, #517, #465, #518, #568, #573, F14, F49, F50, F55, F61-F63,
  F65); release notes state which numerical results change; new docs
  pages *Science* and *Accuracy* (the results on the mock benchmark);
  baseline of the benchmark moved to 2.2.0. Published on GitHub with the
  macOS app; not on PyPI (INAF has not answered yet, #549).
- **2.3.0, after phases 5 and 5-bis (November 2026):** the faster System Inspector, the split of the large files, the
  GUI review of phase 5-bis, F7; the results of fitting, identification
  and the pipelines unchanged (the *Accuracy* page stays as for 2.2.0);
  the release notes state the three fixes that change numbers elsewhere:
  F69 and F70 (calibration) and #630 (errors of linked parameters).
  A stable point before phases 6 and 6-bis, which change the interface
  of the recipes. Published on GitHub with the macOS app; not on PyPI
  (INAF has not answered yet, #549).
- **2.4.0, after phase 6 (November 2026):** error handling: the broad
  except clauses 189 → 47, all boundaries; the fallbacks that went on with
  a result that looked valid stop with their reason (F95); one logger per
  module and one rule for the messages (#566); the hook for unexpected
  errors in the GUI; F94 (data loss when a save failed), F96-F100. The
  results of fitting, identification and the pipelines unchanged (the
  *Accuracy* page as for 2.2.0); the release notes state the changed
  behaviour on failures. A minor version, not 2.3.1, because the
  behaviour changes (Guido agreed). A stable point before phase 6-bis.
  Published on GitHub with the macOS app; not on PyPI (INAF has not
  answered yet, #549).
- **2.5.0, after phase 6-bis (November 2026):** the recipe API: one
  registry of the recipes, typed parameters, exceptions instead of `0`,
  the structured result (`run`), the prerequisites (estimated, asked or
  said), the texts reviewed, the MCP pilot on the new API; F12, F46, F48,
  F49, F50, F55 (parameters), F79, F102-F105. Numbers: the *Accuracy*
  page unchanged; changed on purpose F103 (spectra without a continuum)
  and F102 (`lya_auto`). A minor version, as 2.4.0 (Guido): the interface
  of the recipes for scripts changes (exceptions), the logs replay
  unchanged. Published on GitHub with the macOS app; not on PyPI (INAF
  has not answered yet, #549).

### 3.1 User-reported issues
GitHub issues are the single source of truth; this section only says *when*
each open issue is addressed. Rule of thumb: bugs before refactoring,
features after, except urgent features in stable areas.

| When | Issues | Notes |
|---|---|---|
| With the natural-language log work | [#424](https://github.com/DAS-OATs/astrocook/issues/424) | logging of manual continuum edits |
| Phase 6-bis (item 3); after the full MCP server (item 2, Guido, November 2026) | #573, items 2 and 3 (closed with item 1 in phase 4) | a user configuration file for the recipe defaults and fit settings (item 2); the session log records the values actually used, not "default" (item 3, PR 2 of phase 6-bis); to be opened as issues when scheduled |
| Phase 6 | [#566](https://github.com/DAS-OATs/astrocook/issues/566) | centralize the log messages (started in phase 3 with `recipes/messages.py`), with the move to per-module loggers; done in phase 6: the loggers per module in PR 3 (#693; `tests/test_logging.py` checks that no module logs on the root logger), the shared messages in `recipes/messages.py` with the PRs per package (#694-#700); the issue to be closed (Guido) |
| After phase 5 | [#471](https://github.com/DAS-OATs/astrocook/issues/471) | combining order-separated spectra; touches `SpectrumV2` (coadd/stitch) |
| After phase 5 | [#478](https://github.com/DAS-OATs/astrocook/issues/478), [#479](https://github.com/DAS-OATs/astrocook/issues/479) | GUI features; wait for the GUI split |
| After phase 5 | [#446](https://github.com/DAS-OATs/astrocook/issues/446), [#437](https://github.com/DAS-OATs/astrocook/issues/437) | flux-calibration filters; independent, not urgent. With them (Guido, October 2026, after F70): a more sophisticated algorithm for the wavelength dependence of the calibration curve (today linear between the effective wavelengths of the filters and constant beyond, solved to reach every magnitude); to be opened as an issue when scheduled |
| When needed (not scheduled; Guido, October 2026) | — | transitions of the V1 table not in the single table of #568 (H2 with non-unique names, CO, Paschen, Brackett, SiI, Balmer beyond H_e, stellar and blended emission markers): added when an analysis needs them, each with its source (the V1 file is in `tests/reference/v1/data/atom_par.dat`) |
| To be decided (Guido) | [#455](https://github.com/DAS-OATs/astrocook/issues/455) (pending question) | standard multiplets (`STANDARD_MULTIPLETS`) for the 10 metal species with 3+ transitions in the table and no entry; unchanged by #568, which touched only the table of transitions; the choice of the multiplets used by the identification belongs to phase 7 |

The closed issues (phases 0 to 6-bis) are in `ROADMAP_HISTORY.md`, "3.1 User-reported issues", one row each.

Each bug fix comes in a small PR with a test reproducing the bug.
Only issues on GitHub are considered; requests received elsewhere should be
opened as issues first.

### 3.2 Findings register
Questionable behaviours and technical debt found during phase 1. Nothing
here is fixed silently: at the end of the phase Guido decides, for each
item, whether to open an issue, fix it within a later phase, or accept it.
Exception: anything that loses user data or crashes is reported at once.

| ID | Area | Finding | Impact | Recorded by | Status |
|---|---|---|---|---|---|
| F15 | `equalize_and_stitch` | The scale factor is a median of the ratio of the two arms in the overlap: with narrow absorption (e.g. tellurics) seen at two different resolutions, it is biased by ~1-1.5% in a mock with equivalent-width-conserving tellurics around 1000 nm (no bias without tellurics). Possibly related to #477; to be checked on real X-Shooter VIS/NIR data | numerical | — (mock in the PR that added this row) | open |
| F17 | continuum (forest) | On the J0103-1305 window, with the same parameters and matched smoothing, V1 `clip_flux` and V2 `find_absorbed` + `fit_continuum` agree on the red side (masks 99.98%, continuum within 0.4%), but in the Lyman-alpha forest V2 marks more pixels as unabsorbed (57.6% vs 54.8% of the window, masks agree on 92.7% of the forest pixels) and its continuum is lower by 1.2% (median; 5th percentile -4.6%). Not explained by the fudge (F3). Relevant to #465 and to the known low forest continuum | numerical | `test_continuum_forest` (strict xfail) | partly addressed by #465 (measured on the mock benchmark (2.1.0 + F49/F50 → #465): forest continuum 0.936 → 0.943 (5th percentile 0.820 → 0.830), red side unchanged (1.002)); V1 parity no longer sought for the continuum (V2 is measured against the truth on mocks instead) |
| F28 | `.acs2` Bayesian samples | The posterior samples stored by `fit_bayesian` in the system list metadata (`bayesian_results`, used by the corner plot) are dropped on save without a warning: the FITS header cannot hold them. The fitted values and errors are kept | data loss (recomputable) | `test_bayesian_results_are_saved` (strict xfail) | open; since phase 6 (PR 7) a WARNING says at the save that the samples are not saved (astropy skipped them with a warning of its own) |
| F30 | `running_mean`, `running_std` | Both use `np.convolve(mode='same')`, which pads with zeros: the first and last h values are biased towards 0 (4/7 of the true mean at the edge for h = 3). V1 `running_mean` repeats the first and last full-window value and ignores NaNs. Used by `find_absorbed_regions` on each interval, with half-windows of hundreds of pixels (5000 km/s / 8 on J0103-1305). Replacing it with the V1 function raises the agreement of the forest mask with V1 from 92.7% to 95.7% and the 5th percentile of the continuum ratio from -4.7% to -3.0% (median unchanged, -1.2%): a partial cause of F17 | numerical, V1/V2 difference | `test_running_mean_and_std_pad_the_edges_with_zeros` (characterization) | fixed in phase 4 for the continuum (#465): `find_absorbed_regions` no longer uses `running_mean`, but a masked mean over a window of the smoothing length in km/s on the pixel grid, exact at the edges (`_masked_window_mean`). The edge bias was worth 0.81 → 0.93 within 2000 km/s of the edges of the forest interval on a mock at S/N 30. The plain normalized convolution of the earlier branch (window in pixels of the clipped array) failed at S/N ≥ 300 (forest continuum up to 1.66 times the true one: the window grows as pixels are clipped, over wavelengths where the continuum varies by more than its noise); the window in km/s does not. `running_mean` and `running_std` themselves still pad with zeros (used by `calculate_running_std`): open, minor |
| F51 | `identify_lines` (forest) | Every absorbed region of the forest without a metal candidate is assigned to Ly-a (score 0, no test), even when Ly-a is not among the requested `multiplets`: `multiplets='CIV'` on the full-chain mock gives 233 Ly-a and 6 CIV candidates. In the MCP pilot the Ly-a candidates crowded out the metals in the output (833 on J0103-1305). Found in the MCP pilot (report) | behaviour | — | open: phase 4 / 7 |
| F52 | `VoigtFitterV2` (uncertainties) | Components that end at the bounds of a fit get uncertainties far beyond any meaning, stored as results without a warning: in the `doublets_auto` fit of J0103-1305 (MCP pilot), four CIV components at z = 3.927 with logN = 10.0 and b = 1.0 km/s have dz up to 131, dlogN up to 3800 and db up to 9e4 km/s; 21 of the 101 new components have dlogN > 1 or db > b. Related to F20, F21. Found in the MCP pilot (report) | numerical, error handling | — | open: with [#517](https://github.com/DAS-OATs/astrocook/issues/517); since its PR 1 the components at a bound carry a flag (`ComponentDataV2.flags`), so their errors can be told apart |
| F53 | continuum (forest at high z) | On an ESPRESSO spectrum of J0959+1312 (z_em ~ 4.06), `estimate_auto` flags 100% of the pixels between 525 and 600 nm as absorbed, so the continuum there is a straight interpolation between the edges of the interval (flat, above the 99th percentile of the flux, ignoring the wing of the Ly-a emission); with `template=True` it follows the template, still with no free pixel, and lies below the flux peaks at 565-590 nm (6% of the pixels above C + 2 sigma). Reported by the agent in Guido's MCP trial; not reproduced here (the spectrum is not in the repository) | numerical | — | open: phase 7 (emission model); with #465 |
| F54 | identification (tellurics) | No telluric masking: the O2 A band (759-770 nm) gives spurious doublets with good scores, CIV at z = 3.927 in J0103-1305 (4 components at the bounds, F52) and CIV z = 3.928 with MgII z = 1.733 in J0959+1312 (fitted together, reduced chi2 12). Found in the MCP pilot (trials) | behaviour | — | open: phase 7 (identification); telluric data in the pre-release `data-telluric-2026` |
| F56 | continuum (forest) | Measured on the mock benchmark (phase 4 baseline, 135 mocks): `estimate_auto` gives a forest continuum (Ly-b to Ly-a) low by 6.4% in the median, 13% at S/N 10, 12% at z_em = 4.5, with 5% of the forest pixels more than 18% low (34% at z_em = 4.5); the red side is within 0.2%. Quantifies F17 and the known low forest continuum (F53 on real data) | numerical | `benchmarks/validate_mocks.py`, baseline | improved by #465 (measured on the mock benchmark (2.1.0 + F49/F50 → #465): forest continuum 0.936 → 0.943 (5th percentile 0.820 → 0.830), red side unchanged (1.002); S/N 10: 0.875; S/N 5: 0.80, the limit of clipping). Also found: spikes falling on absorbed pixels and raised to the continuum level become anchors and raise the forest continuum slightly, i.e. the bias comes from the shortage of anchors at the continuum level. Open: phase 7 (emission model) |
| F57 | identification (completeness) | On the mock benchmark the metal doublets are found only well above the detection limit: none of the systems whose weaker line is below 3 sigma, 5% of those at 3-5 sigma, 30% at 5-10 sigma, 91% above 10; completeness over the systems above 3 sigma 0.85 for CIV, 0.69 for SiIV, 0.91 for MgII. Purity is high (0.97-0.99) | behaviour | `benchmarks/validate_mocks.py`, baseline | improved by #465 (wider window on the red side, 400 km/s instead of ~100: weak lines are flagged better): completeness 3-5 sigma 0.03 → 0.10, 5-10 sigma 0.35 → 0.44, CIV 0.857 → 0.874. Why: `find_absorbed_regions` is both the continuum step and the line detector (`identify_lines` takes its regions from `abs_mask`); a weak line pulls down the running mean of a narrow window by a good fraction of its own depth, so its pixels do not reach 2 sigma; with a window 4 times wider the line weighs 4 times less, and purity is unchanged (clipping still at 2 sigma). Direction for phase 7: decouple the detection from the continuum, e.g. search the lines on the flux normalised by the final, smooth continuum, with a filter matched to the expected line width instead of a per-pixel threshold. Open: phase 7 |
| F58 | `doublets_auto` (decomposition) | On the mock benchmark the systems with several components are decomposed differently from the truth: where one fitted component stands for several true ones the pulls of the matched component are 4-6.5 wide, and spurious components appear next to real ones (F49); the total column of each system is right (median difference -0.002 dex, pull width 1.0-1.1 at R >= 45000). The pipeline adds one component per candidate and never adds components from the residuals | behaviour | `benchmarks/validate_mocks.py`, baseline | open: phase 4 or 7 (Guido, October 2026): since F49 the pipeline starts from one component per peak, so a system of several components is fitted with fewer (true components matched 1534 → 1157, pull width of the system column 1.24 → 2.53); to try: `optimize_system` (adds components while the AIC improves, penalty 5, up to 5 more) after the fit of `doublets_auto`, measured on the benchmark (quality and time); then the fit driven by the residuals of phase 7 |
| F59 | fit at low resolution | At R = 8900 the total column of the metal systems is off by 0.39 dex (spread) with a pull width of 2.5, against 0.03-0.07 dex and ~1 at R = 45000-140000. Possible causes: unresolved saturated components (physics), and the convolution kernel computed at the median wavelength (rest of F10, up to +53%/-26% in width at the edges of a window) | numerical | `benchmarks/validate_mocks.py`, baseline | partly fixed by #526 (PR 3, exact convolution): at R 8900 the spread of the system column goes from 0.175 dex (master after #465; 0.39 in the 2.1.0 baseline) to 0.102, the median reduced chi2 from 2.45 to 1.61 and the flagged components from 30 to 16; R 45000 and 140000 unchanged (0.046 and 0.015 dex). After F61: 0.095 dex at R 8900, 0.041 at R 45000. What is left at R 8900 (median reduced chi2 1.6) is probably the physics of unresolved saturated components; to be looked at in phase 7 |
| F68 | `identify_lines` (`Ly_ab`) | The Ly-a/Ly-b pairs of the forest are never identified: `STANDARD_MULTIPLETS['Ly_ab']` lists the red line first (`['Ly_a', 'Ly_b']`), while `rate_doublet_candidate` expects the blue line first, as for the metal doublets (`CIV_1548`, `CIV_1550`): the two windows in velocity never overlap, and the score is always 0 (checked on clean mocks, logN 13-14.5). `Ly_ab` is in the default list of the identification recipes, yet the Ly-b line of every pair becomes a spurious Ly-a system (z = 2.29 for a pair at z = 2.90); with `bypass_scoring` the pair is recognised. Found in phase 5 while writing the characterization tests of the line methods | identification | `test_ly_ab_pairs_are_never_rated` | open; phase 7 (Guido's decision: simple to fix, not urgent; it changes the identifications, measured on the mock benchmark) |
| F72 | `rebin`, `resample_on_grid` | The boolean masks (e.g. `abs_mask`, the absorbed regions) are dropped with a warning: `get_resampled_column` interpolates only numbers (`rebin` catches its `TypeError`, `resample_on_grid` checks the type first). After a rebin or a resampling the mask of the user is lost and the pipelines make a new one. A mask could follow the new grid (e.g. a pixel masked if any old pixel in it was). Found in phase 5 with F71 | behaviour | `test_rebin_drops_the_masks`, `test_resample_drops_the_masks` | open; phase 7 (agreed with Guido: it changes the results and needs a choice of method, how a mask follows a new grid, on the inputs of the continuum; only the recipes *Rebin*, *Resample* and *Co-add* use these methods, not the pipelines nor the mock benchmark; meanwhile, *Find absorbed regions* after a rebin) |
| F84 | x axis in velocity | The spectrum does not keep the reference redshift of a conversion of the x axis to velocity (`x_convert`; only the history line "X-axis converted to km/s at z_ref=..." records it): every conversion back to wavelength (rebin, smooth, continuum, combination of spectra, the ticks of the systems since F66) uses z_em, so they are wrong if the axis was converted with another z_ref, or if z_em changed after the conversion. Found in phase 5 while fixing F66 | design | — | open (phase 6-bis or 7: keep z_ref with the spectrum) |
| F93 | `identify_lines` (identifications) | The identification is by doublets only: two absorbed regions are a candidate when they overlap once one is shifted by the separation of the doublet in velocity (the ratio of the wavelengths), and the pair is accepted when the R^2 of F_1 = F_2^gamma (gamma from f and lambda of the two lines) is above 0.5 (metals only redward of Ly_a). Three things found while rebuilding the Identifications window: (1) the pairs are not stored: each region keeps only its candidates ([transition, score]), and the window finds the pairs again from the overlap of the redshift ranges and the shared score; (2) the multiplets with more than two lines (FeII, SiII, NI, ...) are never searched, only those with exactly two members; (3) the regions of the forest that no doublet identified are labelled Ly_a with score 0, an assignment and not an identification (lya_auto calls the recipe with "Ly_a" alone, not a doublet, so all its regions are so labelled). Found in phase 5-bis | behaviour | `test_identification_viewer.py` (the window) | to be discussed in phase 7 (identification, Guido's algorithmic choices): store the pairs with their redshift and score; the window shows the doublets as the recipe finds them and the assumed Ly_a as such |
| F101 | script console (`gui/scripting.py`) | `AppController.run_recipe`, the `app.run_recipe(...)` of the scripts of the console, calls `MainWindowV2._on_recipe_requested` with a keyword `session_override` that the method does not take: every call fails with a TypeError. Found in phase 6-bis (PR 3a) while adapting the callers of the recipes | bug | — | open: with the decision on the console (kept in the code, phase 5-bis, Guido) |
| F102 | `lya_auto`, `forest_auto` (cleanup of the narrow lines) | After removing the lines narrower than `min_b`, `_remove_narrow` calls `refit_all` without `region_limit`: every group of the session is refitted, also the components outside the forest (e.g. a CIV doublet on the red side fitted before), whose values change slightly. Found in phase 6-bis (PR 4b) while recording the reasons of the removals. `forest_auto` is not limited to a region at all: it identifies on the whole spectrum and its first refit, before the cleanup, also refits every component of the session | behaviour, numerical | `test_lya_auto_does_not_refit_the_red_side` | `lya_auto` fixed in phase 6-bis (PR 4c, Guido's decision): the refits of the cleanup limited to `lya_forest`, as the pipeline; documented in the tutorial (*Refit All* fits everything again). The mock benchmark does not run `lya_auto`: a targeted test replaces it. `forest_auto`: open, Guido's decision (limiting it means choosing its region, e.g. `all_ly_forest`, for the identification and both refits) |
| F107 | GUI menus and the registry | The labels of the menus of the GUI are written in `gui/main_menus.py` and differ from the `brief` of the recipes in the registry, which the dialogs show as their title (*Auto-detect Doublets...* opens *Auto-fit metal doublets (Red Side)*; *Set Properties...*, *Set session properties*). The methods text (PR 5c) names a recipe by its menu and `brief`, so the name it gives may differ from the label of the menu. Found in phase 6-ter (PR 5c) | cosmetic | — | open: deferred to phase 6-quater (Guido, after PR 5c): the menus from the registry, or a `menu` key in it |
| F111 | `metals_auto` (and the fitter's fallback) | `VoigtFitterV2` fits the whole spectrum when no valid pixel is around the components (`voigt_fitter.py`, warning "No valid pixels around the components to fit: the whole spectrum is used"). `metals_auto` met it in both trials of PR 7: on the full-chain mock (31 thousand pixels) it kept no ion, missing the SiIV at z = 1.90 (logN 13.3, b 12, in the forest; the agent found it by eye and it fits with reduced chi2 1.1), and it does not say which ions it tried nor their AIC; on J0103-1305 (215 thousand pixels, 65 components, six ions) it ran more than 2 h 20 min without finishing, with no progress line. The cause of the missed SiIV is to be confirmed (why no valid pixel: the window of the trial ion, its mask, or the names of the transitions) | numerical, performance | — | fixed (PR 7-quater). Cause found on the mock of trial 1 (a script on the chain of the trial): the trial fit found the SiIV (logN 13.28, b 11.3, as the agent), and the recipe returned `refit_all` of the session it started from, so `metals_auto` never added an ion, whatever it found (the test of the recipe recorded "no ion kept" on a spectrum where all the trial ions were out of range). The fallback on the whole spectrum was not the cause: it was met by the ions whose lines fall outside the spectrum (7 of the 15 trials, NV, OVI and MgII at z = 1.76 and 1.83: no pixel at all), which it fitted on all the pixels (31 thousand; 215 thousand on J0103-1305, hence the hours) and kept at their first values (logN 13.00), so that fixing the result alone would have added them. The fallback that Guido recalls is another one, in `prepare_fit_context`: a window narrower than 30 km/s replaced by one of at least 150 km/s (now said in a comment). Fixed (Guido's choices): the session with the ions found is refitted; the fitter raises `NoPixelsToFit` (a `ValueError`, with the components and the range of the spectrum) instead of fitting the whole spectrum, so that `add_component` of an ion outside the spectrum fails with that message and `refit_all` lists such a group as not fitted; `metals_auto` skips those ions with a WARNING that lists them, reports its progress (one line per trial), checks the cancellation at each trial, and reports each ion tried with its AIC change and the reason when not kept (`tried` in `details`). Its acceptance criterion: F121 |
| F112 | MCP output of the pipelines (`IdentificationSummary`) | The result of `absorbers_doublets_auto` says `identification.n_systems: 0` and `systems_by_series: {}` while its notes say "3 candidate systems found" and three components are fitted (65 on J0103-1305): the summary counts the candidate systems of the metadata of the spectrum, which the pipeline does not leave there. Both agents of PR 7 noted the contradiction | output | — | open: server or recipe (what the pipeline leaves in `RecipeResult.extra`) |
| F113 | MCP server (concurrency, long recipes) | One lock (`SessionStore.lock`) serializes every tool, the read-only ones too (`session_history`, `components_list`, `spectrum_plot`): while a recipe runs no other call is answered. With no heartbeat (PR 4 left it out) and no cancellation, a client that gives up after its timeout leaves the recipe running and the session out of reach: in trial 2 of PR 7 the agent could not save its work, which was lost. | design, data loss (the work of a session) | — | fixed: a heartbeat after 20 s of silence, and the read-only tools, `session_save` and `session_note` outside the lock of the recipes (PR 7-bis); the cancellation (PR 7-quater): the client's `notifications/cancelled`, or `session_cancel`, stops the recipe at its next check, with an error and no new state |
| F114 | Docker image | `absorbers_fit_bayesian` is a tool of the image, but the image does not install the `[bayes]` extra (emcee, dynesty): the tool fails with an ImportError (recorded as a gap `failed`). Found in trial 1 of PR 7 | packaging | `tests/test_docker_image.py` (the list of the image) | open: the extra in the image (size to measure), or the description of the tool saying it |
| F115 | methods text (`methods_text`) | Seen by the agent of trial 1 of PR 7: the declared gaps (meant for the developers) are in the text for the paper, as written (decision of PR 5c); the first sentence names the spectrum by the name of the session ("The spectrum of spectrum"); the reasons open each paragraph as written, which reads oddly in a paper; `session_methods` gives the text only inside its JSON, so the agent wrote it to a file with its own Python (declared as a gap) | output, cosmetic | — | open: Guido's to judge (the text is his, PR 5c) |
| F116 | flags of the components | No flag for components unresolved (b well below the instrumental b: 1.4-2.8 km/s at R = 45000 in J0103-1305) or saturated (MgII at logN 15.06 ± 0.24), which both agents of PR 7 judged unreliable by hand; the flags of #517 cover the bounds and the window only | behaviour | — | open: with #517 and F52 |
| F117 | `update_constraint` (description and checks) | The description of `expression` gives `'p["z"]'` as its example; the expressions that work name the target, `p['<uuid>'].z` (the docstrings of the method and of `SystemListV2` say so). With the example, the constraint is stored with a warning ("the linked z could not be computed") and the next fit ignores it: in trial 3 of PR 7 the SiIV kept its own z after a refit, and `components_list` showed `linked: p["z"]`. `target_uuid` is stored next to the expression but not used to evaluate it. An expression that cannot be evaluated is kept, not refused | doc, behaviour | — | open; decided (Guido, on PR 7, November 2026): the example corrected (`p['<uuid>'].z`) and an expression that cannot be evaluated refused (`InvalidParameter`, nothing stored), not stored with a warning; `target_uuid` not used to expand a bare name. fixed in PR 7-ter (the check runs also with `is_free` left out, where the link was stored unchecked and kept frozen by the fit with an ERROR line; the example corrected) |
| F118 | `optimize_system` (stop criterion, output) | On the mock of trial 3 (PR 7) the optimizer added nothing to a single CIV that hides two components 9.4 km/s apart, where a second component added by hand lowers the chi2 on fixed pixels from 55.4 to 37.1 (dAIC -12): it tries a component only where a residual pixel exceeds `threshold_sigma` (2.5), and the deepest was 2.1 sigma (the reason is an INFO line of the log, not in the result: the agent saw "1 component changed" and nothing else). A systematic misfit spread over many pixels is not seen. Related to F58 (the pipelines never add components from the residuals) | behaviour, output | — | open: the criterion to phase 7 (a heuristic, Guido's); decided (Guido, on PR 7): the reason of the stop (e.g. the deepest residual against `threshold_sigma`) and the AIC of each component tried in `RecipeResult.extra`, numbers unchanged, fixed in PR 7-ter (`details.optimization`: `start_aic`, `tried`, `stopped`), the criterion open: phase 7 |
| F119 | `doublets_auto` (scoring of the candidates) | The SiIV doublet of the z = 2.2 system of the mock of trial 3 (logN 12.9, b 9, S/N 50; depths 0.55 and 0.33) is found as a kinematic candidate but scored R^2 = 0.483, below the threshold of 0.5, and dropped; 0.485 with the true continuum, so the continuum (on the wing of the NV emission) is not the cause. The candidates dropped are not in the result: the agent found the SiIV by eye | behaviour, output | — | open: the scoring to phase 7 (identification); decided (Guido, on PR 7): the candidates dropped, with their score, in `RecipeResult.extra`, fixed in PR 7-ter (`details.dropped_candidates`: series, z, score, threshold, x range; not the forest doublets, which no score accepts), the scoring open: phase 7 |
| F120 | names of the columns (expression recipes, MCP reads) | `edit_apply_expression` and `edit_mask_expression` document the columns by the names of the GUI (`λ`, `F`, `dF`, `cont`), `file_export_ascii` too (`F, dF, cont`); `spectrum_read` and `spectrum_plot` want `x`, `y`, `dy`. In trial 3 of PR 7 a read with `F`, `dF` failed (the error listed the valid names) | cosmetic, output | — | open; decided (Guido, on PR 7): `spectrum_read` and `spectrum_plot` also accept the names of the GUI (`λ`, `F`, `dF`, ...) as aliases, fixed in PR 7-ter (`SpectrumV2.columns_in_range`, values under the name asked; `spectrum_plot` takes no columns); the one set of names in phase 6-quater. Fixed in phase 6-quater, PR 8 (Guido, 12 November 2026): the names of the GUI everywhere a column is named (`core/column_names.py`), `lam` (`lammin`, `lammax`) its ASCII spelling, the names of the data (`x`, ..., `dy`) silent and permanent aliases; two cases that failed found there: `calculate_running_std` did not find `F` (the button of its dialog), and wrote a column `F` that no expression could read |
| F121 | `metals_auto` (acceptance of a trial ion) | `metals_auto` kept a trial ion if its fit gave logN >= 12 and 2 <= b <= 50 km/s; `aic_penalty` was accepted and not used, while its description, the `details` of the registry and the methods text said that an ion is kept if it decreases the AIC by it. Once the ions found were kept (F111), the criterion kept 192 ions of 239 on J0103-1305, 67 of them with b = 50 km/s (the bound), 9 with logN = 19 (the bound), 30 with an AIC change above -5; on the mock of trial 1 a spurious SiIV at z = 1.83 (logN 12.80, b = 50.0) | behaviour | — | partly fixed (Guido's choice, PR 7-quater): an ion is not kept with b or logN at the bounds of the fit (within 0.1 per cent) nor with an AIC change above `-aic_penalty` (against the same model without it, on the pixels of its fit), and each ion tried is reported with the reason (`tried`). Open (phase 6-quater or 7): on the mock of trial 1 a spurious CII at z = 1.83 in the forest is still kept (logN 12.58, b 18, AIC change -33); a trial ion that absorbs a Lyman-alpha line improves the AIC anyway, and the model without it keeps the neighbours refitted with it (a bias towards keeping): the comparison with the alternative identification (Lyman-alpha) is for phase 7; measured on the mocks (the benchmark does not run `metals_auto`). Measured on J0103-1305 raw with the code of master after #747 (November 2026; `edit.set_properties` z_em = 4.07, R = 45000, `continuum.estimate_auto`, `doublets_auto`: 65 components, then `metals_auto` on the six default ions, through `recipes.result.run`): 320 trials built (64 redshifts), 81 not tried for their lines outside the spectrum (as before), 239 fitted, **82 kept** (CII 12, CIV 5, NV 17, OVI 18, SiIV 30; against 192 with the old criterion, 67 of them at b = 50 km/s; none now), 157 not kept: b at its bound 108, logN below 12 29, AIC change above -5 7, b below 2 km/s 7, logN at its bound (19) 6. Time: the trials 386 s (6.4 min) and the final `refit_all` 121 s, 507.5 s the recipe, against 26.7 min for the 239 trials of the old criterion; the whole chain 8.9 min; 147 components at the end. The run interrupted at 70 per cent with the rules of F121 kept 60 of 182 tried (33 per cent), the complete one 82 of 239 (34 per cent) |
| F122 | `calculate_running_std` (`output_col`) | With `output_col` empty or the name of a column of every spectrum (`x`, `xmin`, `xmax`, `y`, `dy`) the result went to an auxiliary column of that name, which the table, the dialogs and the reads do not show (the column of the spectrum unchanged): the result was lost without a word. Found in the review of PR 7 (phase 6-quater) | behaviour | `tests/test_recipe_review_file_flux.py` | fixed (Guido, PR 7 of 6-quater): `InvalidParameter` before the recipe runs |
| F123 | `equalize` and `coadd` (`order = -1`) | `-1` means two things: a spline in `flux.equalize` (`order`), no scaling in `edit.coadd` (`equalize_order`). Found in the review of PR 7 (phase 6-quater) | cosmetic, interface | `test_coadd_minus_one_is_none` | fixed in PR 9 of 6-quater (Guido, option A): `equalize_order` of `coadd` is `Optional[int]`, None no scaling (as `kappa` None no clipping), 0-5; -1, of the old logs, is read as None (silent and permanent), below -1 refused; `equalize` keeps -1 (a spline) as its lowest value (PR 7) |
| F124 | `rebin` (`xstart`, `xend`) | With `xstart` not below `xend` the recipe gave a spectrum without pixels, without a word. Found in the review of PR 7 (phase 6-quater) | behaviour | `tests/test_recipe_review_file_flux.py` | fixed (Guido, PR 7 of 6-quater): `InvalidParameter`; the same check for the range of `coadd` in PR 9 (done, `test_coadd_refused`) |
| F125 | GUI, the *Help* button of the dialogs | The `url` of the recipes in the registry are relative names of pages of the V1 documentation (`file_cb.html#export`, `edit_cb.html#rebin`, ...; some marked "Placeholder URL"), which `webbrowser.open` cannot open as pages of the documentation; some recipes have none. All the recipes. Found in the review of PR 7 (phase 6-quater) | cosmetic | — | open: PR 11 of 6-quater (F107, the menus from the registry): the pages of `docs/` |
| F126 | `apply_expression` (units) | A new column made by an expression is dimensionless whatever the expression: `F_corr = F * 2` loses the unit of the flux (`erg/(nm s cm2)`); a column replaced keeps its unit, also when the expression changes it (`F = F / cont`). Found in the review of PR 8 (phase 6-quater) | behaviour | — | open (Guido, PR 8: registered, not changed; a unit from the expression would need the units of NumExpr's operations); documented in the tutorial of Edit |
| F127 | expression recipes (the axis) | `apply_expression`, `mask_expression` and `smooth_column` wrote the columns of the axis: `λ = λ * 0 + 1` replaced the wavelengths (without `λmin` and `λmax`), a mask put NaN in them, the smoothing smoothed them. Found in the review of PR 8 (phase 6-quater) | behaviour | `test_the_axis_is_not_written` | fixed in PR 8 (Guido's choice): refused (`InvalidParameter`), the axis is changed by `x_convert`; a correction of the wavelengths (`λ * (1 + v/c)`, air to vacuum) would need the three columns together, a recipe of its own if wanted |
| F128 | expressions (functions) | `_NE_GLOBAL_DICT` (`core/spectrum_expressions.py`) offers `median`, `mean`, `std` (NumPy, NaN-aware) and `min`, `max` (`nanmin`, `nanmax`) to NumExpr, which cannot call Python functions: `median(F)`, `mean(F)`, `std(F)` failed with "'VariableNode' object is not callable", while `min(F)`, `max(F)` are NumExpr's own reductions (a NaN gives NaN), as `sum` and `prod`. The comment in the code says otherwise. Found in the review of PR 8 (phase 6-quater) | behaviour | `test_function_that_numexpr_does_not_have` | open: the message says which functions there are (PR 8); the dictionary is left as it is (removing it changes nothing, offering the NaN-aware ones would need evaluating them outside NumExpr: Guido's choice) |
| F129 | CSV files (`export_ascii`, the CSV loader) | The CSV files of `export_ascii` held no units, and the loader of CSV files took the first three columns as wavelength in nm, flux and error in `erg/(Å s cm2)`: a spectrum exported after `y_convert` or with the axis in km/s came back with the wrong unit, `cont` was lost, and the file of the whole spectrum (`targets='spec'`: `x,xmin,xmax,y,dy,...`) came back with `xmin` as the flux. `import_ascii` read a column `F` as a new column `F`. Found on PR 8 (phase 6-quater); the last two in the characterization of PR 8-bis | behaviour, data | `tests/test_flux_units.py` | fixed in PR 8-bis of 6-quater (Guido): ECSV with the units and z_em, the names of the GUI in ASCII (`lam`, `F`, `dF`, ...), read back by name; the old files by name (`x`, `y`, `dy`), the flux in arbitrary units; `import_ascii` translates the names and keeps the units |
| F130 | `y_convert` (columns without a unit) | `convert_y_axis` converted only the auxiliary columns whose unit is equivalent to that of the flux: in an archive whose `cont` and `model` have no unit and the flux has one (J0103-1305 of `tests/data`), `y_convert` from Å to nm multiplied the flux by 10 and left the continuum, so the normalized flux changed by 10. Found in the analysis of PR 8-bis (phase 6-quater) | behaviour, numbers | `test_y_convert_converts_the_continuum` | fixed in PR 8-bis of 6-quater (Guido): `cont`, `model`, `cont_pl` without a unit are in the unit of the flux, both when an archive is read and in the conversion |
| F131 | `calibrate` (unit of the result) | `Spectrum.calibrate` multiplies the values by a curve computed in the unit of the flux, and labels the result `erg/(Å s cm2)`: a flux in `erg/(nm s cm2)` is calibrated in erg/nm and labelled per Å, 10 times brighter than it says (SDSS_g = 15.0 for 17.5 asked). `test_spectrum_calibrate` (`tests/test_photometry.py`) recorded it ("the unit is relabelled"). A flux without a unit or in `arb_flux` is read as `erg/(Å s cm2)`, and its result is right. Found in PR 8-bis (phase 6-quater) | behaviour, units | `test_calibrated_spectrum_has_the_magnitude_asked_for` | fixed in PR 8-bis of 6-quater (Guido): the unit of a physical flux is kept, `erg/(Å s cm2)` for an arbitrary one or one without a unit; the values unchanged |
| F132 | `equalize_and_stitch` (units of the flux) | The arms were equalized on the values of their flux, then joined by `stitch`, which converts to the unit of the bluest: an arm in another physical unit (per nm against per Å) came out off by the ratio of the units (10 times too faint). With `manual_factors` the result was right (the factor applies to the values in the unit of the arm). Found in the analysis of PR 9 (phase 6-quater) | behaviour, units | `test_equalize_and_stitch_converts_the_flux_before_scaling` | fixed in PR 9 of 6-quater (Guido): the arms converted to the units of the bluest before the scaling; same results with arms in one unit (the fingerprints of `test_valid_results_unchanged`) |
| F133 | `coadd` (units of the flux) | The values of the flux were co-added as they were and labelled with the unit of the first session: with `equalize_order = -1` a session per nm was 10 times too bright against one per Å, with the scaling right by chance (the factor absorbs the ratio); a flux in arbitrary units or normalized was co-added with a physical one without a word. Found in the analysis of PR 9 (phase 6-quater) | behaviour, units | `test_coadd_converts_the_flux` | fixed in PR 9 of 6-quater (Guido): converted to the units of the first; an arbitrary flux with a physical one only when scaled (it takes the unit of the reference), a normalized one with one that is not refused |
| F134 | `coadd`, `equalize_and_stitch` (units of the axis) | The overlap for the scaling was looked for in the values of the axis (`compute_flux_scaling` takes values): a session in Å against one in nm had none, so it was not scaled, without a word; `equalize_and_stitch` also applied the cuts and the ranges (nm) in the unit of each arm (an arm in Å kept whole, its pixels twice in the overlap) and ordered the arms by the values of their axis (an arm in Å after one in nm, a part lost). Found in the analysis of PR 9 (phase 6-quater) | behaviour, units | `test_an_arm_in_angstrom_is_equalized`, `test_equalize_and_stitch_orders_the_arms_in_one_unit` | fixed in PR 9 of 6-quater (Guido): the axis converted to the unit of the reference, the arms ordered in nm, the cuts and ranges converted to the unit of the axis; an axis in velocity against one in wavelength refused (`x_convert` first) |
| F135 | recipes that read other sessions, without a GUI | In a script (a session of `load_session_from_file` with `gui=None`) `coadd`, `equalize_and_stitch`, `import_systems` and `trim_common` failed with an `AttributeError` (`None.session_histories`): a crash, reported to Guido at once (PR 9); `stitch` worked only with sessions given as objects, and its log held their `repr`. Found in the analysis of PR 9 (phase 6-quater) | crash | `test_without_a_gui_a_name_fails_with_a_hint`, `test_in_a_script_the_sessions_themselves` | fixed in PR 9 of 6-quater: one rule for the parameters of kind `session` and `sessions` (`recipes/_sessions.py`): names, or in a script the sessions themselves; the name of the session the recipe runs on is that session; the log records names |
| F136 | recipes that read other sessions (a name not found) | Each did its own: `coadd` skipped a name not found without a word (and co-added the others), `stitch` with a WARNING, `trim_common` with a WARNING and then "no range in common" (or "N sessions and M transitions"); `trim_common` also gave "no range in common" for an unknown transition of its session and ignored one of another session (an infinite range). Found in the analysis of PR 9 (phase 6-quater) | behaviour | `test_a_name_not_found_fails_naming_it`, `test_refused` | fixed in PR 9 of 6-quater: `RecipeFailed` naming it; an unknown transition `InvalidParameter` |
| F137 | `coadd` (this session) | The sessions co-added are those named: without this session in the list, the result held the others only, in a new state of this session (its name, its systems), while the description said "with this one"; the first named, not this session, is the reference (units, scaling). Found in the analysis of PR 9 (phase 6-quater) | behaviour, interface | `test_coadd_refused` | fixed in PR 9 of 6-quater (Guido): a list without this session refused; the first named stays the reference (described) |
| F138 | `flux.equalize` (units) | The same as F132 and F134 in `Spectrum.equalize_to_reference`, used by `flux.equalize`: a session per nm equalized to one per Å keeps its unit with the values of the other (10 times too faint); a session with the axis in Å is not scaled to one in nm (no overlap of the values), without a word. Found in PR 9 (phase 6-quater), not changed there (it changes the numbers of a recipe outside the PR; `equalize_and_stitch` and `coadd` convert before calling the core) | behaviour, units | — | open: to Guido (the same conversion in `equalize`, or in `equalize_to_reference`) |

Closed findings (fixed or accepted: all the others, F1 to F106 and F108 to F110) are in `ROADMAP_HISTORY.md`, "3.2 Findings register", one row each, to be found by number (`grep "^| F52 "`).

### 3.3 Phase 1-bis: findings and urgent issues
Review of the findings register with Guido (October 2026): the proposals agreed in general, each item a small PR with its strict xfail turned into a passing test and its numerical effect stated. Done (October 2026, release 2.0.2 with #534): F1, F2, F5, F10 (minimal), F11, F16, F19, F20, F23-F27, F29, F32, F33, #480, #442; F3 accepted; deferred to phase 4 as issues: #526 (F8, F22, rest of F10, F9, F14), #465 (F30, F17), #517 (F21), #518 (F4); F15 waits for data. Full text: `ROADMAP_HISTORY.md`, "3.3 Phase 1-bis: findings and urgent issues".

### Rough estimates
The first estimate (October 2026) was 10–20 work sessions, 2–5 calendar weeks, with human review and test data as the bottleneck; superseded by the phases above. Full text: `ROADMAP_HISTORY.md`, "Rough estimates".

### Decisions reserved to Guido
- License: GPLv3 or later (decided); copyright holder pending INAF (#549).
- Questionable behaviours found by characterization tests (bug or intended?).
- Approach to the continuum spike sensitivity (#465).
- Scientific validation of results.
- Releases (see "Releases" above): timing of 2.1.0 and 2.2.0, macOS app
  signing and notarization. PyPI publication: decided, with the first
  release after INAF answers on the LICENSE.

---

## 4. Direction: Astrocook for agents

Do not embed an LLM inside Astrocook; make Astrocook easy to use for any
LLM.

1. **MCP server:** each recipe as a tool with a name, a description and
   typed parameters (recipes already have structured parameters).
2. **Clean, documented Python API** with structured outputs (JSON:
   parameters with errors, χ², quality flags).
3. **Reliable natural-language log**, ready to paste into a text (e.g. the
   methods section of a paper). Also record the fit *decisions* (components
   added/rejected and the criterion used), not only the results.
4. **GUI as a review bench:** the agent does the repetitive work, the human
   inspects, corrects and approves.

When (agreed with Guido, October 2026):
- **Pilot, right after 2.1.0** (before phase 4): MCP server with 3–4
  recipes (open a spectrum, continuum, identification, fit), as a separate
  module that does not touch the existing code. It installs with `pip`
  from GitHub, without waiting for PyPI; it shows which structured outputs
  and error messages an agent needs, so that phases 4 and 6 can take them
  into account. Scope in section 4.1; done (October 2026), report in
  section 4.2.
- **Full server: phase 6-ter, after phase 6-bis and 2.5.0:** all the recipes, structured outputs,
  natural-language log, and a headless Docker image to distribute it. An
  agent working on its own needs correct and measured results (phase 4),
  clean error handling (phase 6) and a stable recipe API (phase 6-bis).
  The tools are generated from the recipe registry (the schemas the GUI
  already builds its menus and dialogs from), not written by hand: a
  recipe added or changed later reaches the server with no server code,
  and a test checks that every recipe is a tool or is explicitly
  excluded (interactive recipes such as `update_from_knots`). This is
  why the server can come before phase 7: its new and improved recipes
  reach the agents as they reach the GUI. Guido (phase 6-bis, PR 6): the
  server offers as many recipes as possible, also the basic ones that
  today only other recipes call, as blocks for new ones; one day it can
  say which recipes would be worth implementing, from how the agents
  combine the blocks (e.g. the sequences of recipes they run often, and
  what they compute by themselves without a recipe).

### 4.1 MCP pilot: scope
Agreed with Guido (October 2026), done in PRs #595-#600. A thin layer over the recipes, in a subpackage `astrocook/mcp_server/` (optional `[mcp]` extra, `astrocook-mcp`, the official `mcp` SDK, installed from GitHub): five tools (`session_open`, `session_save`, `continuum_estimate_auto`, `absorbers_identify_lines`, `absorbers_doublets_auto`), defaults read from the recipe schemas, a docs page. The goal was to learn what a full server needs. Full text: `ROADMAP_HISTORY.md`, "4.1 MCP pilot: scope".

### 4.2 MCP pilot: report
Two trials with fresh agents (Claude, no knowledge of the code) and Guido's own trial with Claude Code (October 2026, PR #600). Outcome: an agent can carry out an analysis from a request in plain words; found F46-F55 (phases 4, 6-bis, 7); the lessons for the full server are in "Rules and decisions in force", the answers to the questions of 4.1 and the trials in the archive. Full text: `ROADMAP_HISTORY.md`, "4.2 MCP pilot: report".

### 4.3 Full server: trials with fresh agents (PR 7)
November 2026, after PRs #726-#742 (64 tools). As in the pilot (section 4.2):
fresh agents (Claude, no knowledge of the code, the tests or the true
parameters, told not to read them) used the server through a minimal MCP
client (`benchmarks/mcp_trials/client.py`: the instructions, the tools with
their schemas and the content of the results, the images as PNG files; the
progress notifications only in the log of the trial, with their times), each
with a request in plain words (`benchmarks/mcp_trials/prompts.md`), then a
critique. The gaps went to a temporary file (`ASTROCOOK_GAPS_FILE`) and the
agents declared theirs with `session_note(kind="gap")`. The transcripts were
checked for reads of the code: none (one disclosed slip, below). Trials 1
and 2 were run first; trial 3 (the editing tools and undo, on the mock of
`make_data.py edit`) after the heartbeat and the reads during a recipe
(PR 7-bis), before the cancellation.

**Trial 1: the full-chain mock, in the Docker image.** The mock of
`tests/full_chain.py` without the true continuum and the noiseless flux
(`generate(..., truth_columns=False)`; in the pilot the agent fitted on the
true continuum). The server ran in the image of PR 6, built in the cloud
session (the daemon starts there; the base image taken from a mirror and
given the certificate of the proxy, a local workaround only, the
`Dockerfile` unchanged). 35 calls in 4.3 minutes, none over 2.1 s, one
failed (`absorbers_fit_bayesian`: emcee is not in the image, F114); some
120 thousand tokens of context for the agent. Against the truth: CIV at
z = 1.76 and 1.90 right (logN within 0.012, b within 0.2 km/s); CIV at
z = 1.83, two true components 12.7 km/s apart (13.4/8 and 13.9/12), fitted as
one (14.03, b 13.4; the total column right, 14.02), and `optimize_system`
added none (F58); SiIV at z = 1.90 in the forest (13.3/12) missed by
`metals_auto` (F111) and found by the agent by eye on normalized plots,
added and fitted by hand (13.28/11.3). The two spurious components of the
pilot (F49) are gone. Saved session, methods text, notes and reasons as
asked; the only work outside the server, the methods text written to a file
(declared).

**Trial 2: J0103-1305, local server.** The spectrum of the tests raw
(380-780 nm, 215 588 pixels; no continuum, model, telluric model or
components; z_em and R not set). The agent set z_em = 4.07 from the Ly-a peak
by eye and R = 45000 as an assumption (both declared, as a gap: no header to
read, no z_em estimator), estimated the continuum (0.7 s) and ran
`doublets_auto` (25 s, steady progress; 168 s in the pilot): 65 components
in 21 groups, judged one by one with sound reasons (the CIV at z = 3.886 and
3.927 rejected as the O2 A band, F54; the components with b below the
resolution element and a saturated MgII at logN 15 marked unreliable, F116).
Then `metals_auto` on six ions: one notification at the start and none
after; the client gave up after one hour, the server went on at 100% CPU and
had not finished after 2 h 20 min, when it was stopped. Meanwhile the server
answered no other call, read-only ones included (F113): the agent could not
save, and the work was lost, against the request ("I don't want to lose any
work"). Its report lists the steps to redo it. Cause, from the code:
`metals_auto` fits each trial ion with the fitter, which falls back on the
whole spectrum when no valid pixel is around the components (F111). The
agent, waiting, listed the folder of the transcripts of the agents with
`cat` (outside its folder) and so saw the transcript of trial 1 (prompt and
use of the server, no code); it disclosed it in its report; the trial stands.
The client of this trial did not pass `ASTROCOOK_GAPS_FILE` to the server
(the SDK passes only a few variables of the environment): the one declared
gap went to the default file of the cloud session; fixed in the client and
its test.

**Trial 3: the editing tools, local server** (after PR 7-bis). The mock of
`make_data.py edit`: CIV at z = 2.2 in two components 9.4 km/s apart (13.6/9
and 13.3/7), its SiIV (12.9/9), a single CIV at z = 2.25; z_em = 2.6, the
red side only, S/N 50, R = 45000. 55 calls in 3.9 minutes, none over 1.6 s
(8 s in all), one failed (`spectrum_read` with the columns `F`, `dF`: the
names of the GUI, which the expression recipes document, while the reads
want `x`, `y`, `dy`; the error listed the valid names); 129 thousand
characters of results. The agent estimated the continuum, ran
`doublets_auto` (one CIV, 13.79/10.7, and the CIV at 2.25; the SiIV missed,
F119), added the SiIV by hand, judged the single CIV inadequate from its own
residuals (a flat core, narrow wings), tried `optimize_system`, which added
nothing (F118), added a second CIV by hand and kept it (its chi2 on 58 fixed
pixels from 55.4 to 37.1: dAIC -12, dBIC -6), then: a third component (the
optimizer, threshold 2 sigma: none), z(SiIV) tied to the main CIV (kept,
dAIC -2), a SiIV in the second component (logN 10.9 ± 0.5, b at its bound:
undone, dAIC +5), b(SiIV) tied to b(CIV) (the free values 8.67 and 8.98:
kept, dAIC -0.6, dBIC -3.3, and it reasoned on the thermal ratio). Final,
against the truth: CIV 13.622 ± 0.013 / 8.82 ± 0.13 at z = 2.1999996 (true
13.6/9 at 2.2), CIV 13.260 ± 0.032 / 6.88 ± 0.39 at +10.0 km/s (13.3/7 at
+9.4), SiIV 12.890 ± 0.007 (12.9), b tied: every value within 1.5 sigma or
0.04 dex. The model chosen is the true one, and every step of the choice is
in the log with its reason (four undos, three with a reason). The first tie,
made with `update_constraint` and the expression of its description
(`p["z"]`), was stored with a warning, "the linked z could not be computed",
and ignored by the refit (F117); the agent saw it, undid it, and found the
working syntax (`p['<uuid>'].z`) in the constraint that `add_linked_system`
writes. Undo and redo carried the work: the agent went back to compare
models and to drop the alternatives that made things worse, as asked. What it
did outside the server, declared as gaps: the residuals and the chi2, AIC and
BIC of a model on a fixed set of pixels (no tool gives them; the results give
the chi2 of a group only when components are added, on a window that changes
with the model), read with `spectrum_read` and computed with its own Python;
to read the residuals of the model before it undid and redid, which the
server recorded as a recipe undone (a gap `undone` that is not one). It also
asked for a plot of the residuals and of several regions in velocity around
a redshift. Seen again: F112 (`identification.n_systems: 0` with two
candidates found) and F115 (the methods text names the spectrum by the
session, `spectrum`, and carries the gaps and reasons as written; it leaves
out the alternatives undone, which are part of what the astronomer asked).
Also: `components_changed` lists every component a refit touches, also the
ones that moved by 1e-10 (three "changed" after a refit that changed
nothing), and the result of `fit_component` and `optimize_system` gives no
chi2 (the fit summary only when components are added). No slow call: nothing
to say on the heartbeat or the cancellation.

Decisions on the findings of trial 3 (Guido, on PR 7, November 2026): a
short PR of the outputs (PR 7-ter, before the PR of F111, Guido: it is the cheaper and independent one), with F117
(the example corrected, an expression that cannot be evaluated refused),
F118 and F119 (what `optimize_system` tried and why it stopped, and the
candidates `doublets_auto` dropped with their score, in `RecipeResult.extra`;
the criteria to phase 7), F120 (the reads accept the names of the GUI as
aliases) and a tolerance on the changes: `compare_sessions` lists a
parameter as changed only beyond a small tolerance (what `RecipeResult`
reports, not the fits). Not taken now: a redo cancelling the gap `undone` of
the state undone, and the chi2 of the groups in the result of every fit
(also without components added): left open, for the trials or 6-quater.

**The number of tools (64).** The listing a client gives the model is
121 thousand characters (about 35 thousand tokens, estimated), more than half
of the context of an agent at its first call. Both agents found the right
tools at the first try, from the menus of the instructions and the names
`<menu>_<recipe>`, and both advise against grouping them (explicit tools
with precise schemas are easier to call right). What they asked instead: the
parameters repeated on every tool (`session_id`, `reason`,
`provide_missing`, `max_listed`) described more briefly; the tools of other
tasks (stitch, co-add, resample, flux calibration) and the unavailable ones
marked; a paragraph of the usual workflow and the expected time of the slow
recipes in the instructions. The agent of trial 3 the same: the menus made 64
tools workable, some 25 irrelevant to its task cost reading but no wrong
choice; the hard part was the overlap of the absorbers' tools
(`add_component` or `add_linked_system`, `update_constraint` or a linked
system, `fit_component`, `refit_all` or `optimize_system`): the descriptions
should say what happens next (whether a recipe refits) and give examples
that work; a paragraph of typical workflows (fit a system, test a component,
tie ions) would help more than fewer tools.

**Decisions** (the three left open in PR 4 and by the principle of the
phase):
- *Heartbeat for the silent long recipes*: yes. An hour of silence cannot be
  told from a hang, neither by the agent nor by the user.
- *Stop as a cancellation*: yes, with the read-only tools answering while a
  recipe runs (F113): a client that gives up does not stop the recipe, and
  the session is out of reach until it ends. The heartbeat and the read-only
  tools are a short PR of the server (before trial 3); the cancellation needs
  the recipes to check a signal (`metals_auto` logs no progress), so it goes
  with the fix of F111.
  Its design (proposed after PR 7-bis, agreed with Guido, November 2026):
  the server catches the client's `notifications/cancelled`, and a
  `session_cancel` tool serves the clients that do not send it; both raise a
  signal that the recipes check where they already report their progress
  (`recipes/progress.py`, shared with the GUI, which could use it for a Stop
  button); `metals_auto` checks it at each trial ion, and also reports its
  progress and the AIC of the ions tried; a recipe cancelled leaves no new
  state, as an error.
  Found after the design, in the code (trial 3): the GUI already has a Stop
  (`gui/recipe_progress.py`, `_on_stop_recipe`), and the pipelines of the
  absorbers already check it, `RecipeAbsorbersV2._is_stop_requested` (a
  `_stop_flag` of the session, or the worker of the GUI), at their steps and
  in their loops; when stopped, they return the session they started from.
  The signal of the design is to replace or generalize this one, not to add
  a second; to settle in the PR of F111.
  Settled in PR 7-quater (Guido): replaced. The signal is a
  `threading.Event` of the thread that runs the recipe
  (`recipes/progress.py`, `cancellable`, `check_cancelled`, and `progress`,
  which logs the progress line and checks), not of the session, which stays
  immutable; the steps of a pipeline, each on a new session, see it. The
  GUI's Stop sets the event of its worker (no more `_stop_flag` in the
  session nor `_is_stop_requested`); a recipe cancelled raises
  `RecipeCancelled`, a `RecipeError`: the worker emits `stopped` and the
  history is not changed, the server answers with an error and adds no
  state. The tools of the recipes are asynchronous: the SDK cancels the
  scope of a call on `notifications/cancelled`, which does not interrupt a
  thread, so a task next to the thread sets the event (`run_cancellable`).
- *Tools grouped by menu*: no. The cost of the 64 tools is the size of the
  listing, not the choice: to shorten the descriptions of the shared
  parameters, generated from the registry.

**The real analysis** (PR 7 of the plan, Guido's own): moved to phase 7
(Guido, November 2026): with the recipes of today it would mostly measure
limits already known (sections 3.2 and above); it becomes the first test of
the new components of phase 7, with the same criterion (a result he would
sign after reviewing it in the GUI, and a log he could replay).

**Outcome.** The server works as planned for the parts it controls (the
registry, the outputs, history, notes and gaps, the methods, the image), and
an agent chains it at the first try; what failed is a recipe (F111) and the
way the server waits for it (F113). Trial 3 adds the editing tools: with
undo, redo and the constraints an agent reached the true model of a blend
and defended it step by step, and what it missed was a measure (residuals,
fit statistics on fixed pixels) and two descriptions that work (F117,
F120); the automatic steps fell short again (F118, F119). The trials confirm the direction of
phase 7: the agents met the limits of the recipes, not of the interface
(F54, F58, F111, F116, no upper limits, no z_em estimator, no telluric
model).

## 5. v3 vision (horizon, after v2 consolidation)

**Goal:** given a spectrum, produce in less than one hour **a collection of
complete models**, with uncertainties and a reliability map, which also
report the possible alternative identifications of the same feature, each
with a parameter quantifying its credibility.

### A collector of recipes, not an automatic fitter (principle)
Agreed with Guido in phase 6-ter, PR 0 (November 2026), on his question:
with a whole spectrum fitted in a single pass, does Astrocook risk turning
from a collector of recipes into an automatic fitter? The risk is a
monolithic module that takes a spectrum and returns a model: opaque (why
this continuum, these components?), and absorbing the work while the
recipes wither. Three rules keep the single pass a composition of
recipes:
1. every automatic step is a recipe, also offered on its own (in the
   GUI, the scripts and the MCP server), as the pipelines of today
   (`doublets_auto`, `forest_auto`) are made of recipes;
2. every automatic decision goes into the log with its reason (which
   alternative, by which criterion), so that the single pass opens step
   by step in the GUI, replays, and can be corrected at any step;
3. the automatic result is a starting point to review, not a final
   answer: a collection of alternative models with their credibility
   (the goal above), from which the astronomer chooses.
The workflow of the final pass is not one for all spectra (Guido): it
can differ from spectrum to spectrum, optimized for the object and for
the observational parameters as the final spectrum reflects them
(resolution, S/N, wavelength coverage, z_em and the extent of the forest,
DLAs and BALs, tellurics). Whether it is chosen by rules from the
spectrum or by an agent combining the recipes, the choice itself is
logged with its reasons, so the result stays reproducible: the workflow
of a spectrum is its log, and replays as such. A sequence of recipes
that agents or users run again and again on similar spectra becomes a
recipe of its own; the gaps the agents meet (phase 6-ter, decision 6)
say which recipes are missing.

### Iterative loop (deterministic and reproducible)
1. Initial continuum with uncertainty.
2. Detection of absorbed regions and decomposition into independent blocks
   (regions linked by multiplets = one block).
3. Identification as hypothesis testing: metals, then the Lyman series,
   Lyα as the fallback hypothesis.
4. Parallel Voigt fitting per block, with a free local continuum and
   addition of components driven by residuals and an information criterion.
5. Refinement of the continuum using the model.
6. Repeat until convergence.
7. Output: collection of models, components with uncertainties, reliability
   map (local χ², rejected alternatives and their margin, warnings).

### Model credibility
- Akaike/Schwarz weights (AIC/BIC) as a quick estimate; Bayesian evidence
  (dynesty) for important cases.
- Alternatives defined **per block**, not over the whole spectrum, to avoid
  combinatorial explosion.

### Known difficulties
χ² ≈ 1 is necessary but not sufficient: continuum/absorption degeneracy in
the forest; wrong identifications with a perfect χ²; overfitting; unreliable
error spectrum (correlated pixels); pathological cases (DLAs, BALs,
tellurics, artefacts).

### Validation
Mocks with known ground truth: completeness, purity, bias on N/b/z as a
function of S/N and resolution, misidentification rate.

### Path
First the robust components (phase 7: emission model and identification,
measured on the benchmark of phase 4), then the loop.
Pilot on metals outside the forest, then extension to the Lyα forest. The
LLM agent stays at the margins (strategy in stuck cases, guidance during
review); each of its decisions goes into the log. Revised in phase 6-ter,
PR 0 (Guido): the agent may also choose the workflow of a spectrum (the
principle above); the steps it runs stay recipes, logged with their
reasons.

### Implications for the v2 work
- The mock generator from phase 1 is the first core of the v3 test bench.
- When refactoring, do not introduce the assumption "one model per session".
- The log must record the fit decisions.

---

## 6. Papers

Agreed with Guido (October-November 2026). The plan of the papers (which
papers, journals, order, the structure and method of the reference paper)
is in the private repository `DAS-OATs/astrocook-data`,
`papers/plan.md`; here only what drives the development.

**Reference paper (A, after phase 7).** Validation on mocks and an
automatic re-analysis of the UVES Deep Spectrum of HE0940–1050, assisted
by an agent through MCP, compared with the published VPFIT analysis
(D'Odorico et al. 2016) and with a *synthetic twin* of the spectrum (a
mock with the published components, on the grid, with the error
spectrum, resolution and continuum of the real one).

**Rule: the real spectrum of HE0940–1050 stays blind** until the version
of the paper is frozen (after phase 7). Astrocook is developed and
measured only on the mocks and on the twin, never on the real spectrum,
so that it is not tuned on it; the analysis plan is written before the
first run on it.

**What the development owes to the paper:**

| What | When |
|---|---|
| Twin: a mock from an arbitrary component list, on a given pixel grid, with a given error spectrum, resolution and continuum (`astrocook/mocks` has a uniform velocity grid, a power-law continuum and an error from the S/N today) | small, when the material of the published analysis arrives |
| Import of a VPFIT model (fort.26 or the tables of a paper) into a `SystemListV2`, and its evaluation on a spectrum without fitting (chi2, residuals, AIC/BIC per system) | new, small |
| Full MCP server, with the log of the decisions of the agent | phase 6-ter |
| Emission model with an uncertainty; identification with alternatives and credibility | phase 7 |
| Mock benchmark rerun with the frozen version | after phase 7 |

The twin also gives phase 7 a realistic high-S/N development case, which
the mocks of phase 4 do not cover.

**Software paper (B, JOSS)**: as soon as the LICENSE arrives (#549); the
recipe API is stable after phase 6-bis.

---

## 7. Progress

One line per concluded phase; the long rows are in `ROADMAP_HISTORY.md`, "7. Progress".

| Phase | Status |
|---|---|
| Initial analysis and roadmap | done (October 2026) |
| 0 – Environment | done (October 2026) |
| 1 – Safety net | done (October 2026; log replay, V1 `.acs` and loader tests wait for data) |
| 1-bis – Findings and urgent issues (section 3.3) | done (October 2026): see 3.3; F3 accepted; F15 waits for data |
| 2 – Cleanup and packaging | done (October 2026), release 2.0.3; LICENSE pending (#549) |
| 2-bis – Leftovers of phase 2 (F34, Ruff, Sphinx, `[bayes]`) | done (October 2026) |
| 3 – Decoupling from legacy | done (October 2026), PRs #565-#584 |
| 3-bis – Small fixes before 2.1.0 | done (October 2026), PRs #587-#592; release 2.1.0 |
| MCP pilot (sections 4.1, 4.2) | done (October 2026), PRs #595-#600 |
| 4 – Scientific correctness | done (October 2026), PRs #601-#631; release 2.2.0 |
| 5 – Splitting large files (and GUI performance) | done (October 2026), PRs #634-#662 |
| 5-bis – GUI look and ergonomics | done (November 2026), PRs #664-#688; release 2.3.0 with phase 5 |
| 6 – Error handling | done (November 2026), PRs #691-#701; release 2.4.0 |
| 6-bis – Recipe API | done (November 2026), PRs #705-#723; release 2.5.0 |
| 6-ter – Full MCP server | done (November 2026), PRs #726-#747 and the closing PR 8; release 2.6.0 published on GitHub on 10 November 2026 (PR 9, #751: version, CHANGELOG, `switcher.json`); plan of the phase agreed (PR 0, #726); the tools generated from the registry, F106 found (PR 1, #729); F106 fixed (PR 1-bis, #730); the tools of the sessions, undo and redo (PR 2, #731); the tools of inspection (PR 3, #732); progress notifications (PR 4, #733); the agent's decisions in the log (PR 5a); the gaps of the recipes (PR 5b); the log as the text of a methods section (PR 5c); the GUI following a session of the server, F109 found (PR 5-bis); the headless Docker image, F110 found (PR 6, not published); trials with fresh agents, F111-F116 found (PR 7, report in section 4.3); heartbeat and reads during a recipe, F113 in part (PR 7-bis); trial 3, the editing tools, F117-F120 found (PR 7); the outputs of trial 3 (PR 7-ter); F111 and the cancellation, F113 closed, F121 found (PR 7-quater) |
| 6-quater – Review of the recipe parameters | in progress (November 2026): plan of the phase agreed (in PR 1, ruff `B`); PR 2, the aliases of the renamed parameters; PR 3-4-bis, mypy on the whole package; PR 5, `SessionV2` and `SessionHistory` into `astrocook/session/`; PR 6, the `V2` suffixes dropped (#759); PR 7, the review of file and flux, F122-F125 found (#761); PR 8, the review of edit on columns and expressions, F120 fixed (one set of names of the columns), F126-F128 found, F127 fixed; PR 8-bis, the units of the flux (`arb_flux`) and the CSV files in ECSV, F129-F131 fixed; PR 9, the review of edit on several spectra, one rule for the sessions, F123 and F132-F137 fixed, F138 found |
| 7 – Robust components | not started |
