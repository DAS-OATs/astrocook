# Astrocook – Development Roadmap

Working document, based on an analysis of the v2.0.1 code (October 2026).
It collects the assessment of the current state, the consolidation plan for
v2 and the vision for a future v3. Update it at the end of each phase.
The detail of the concluded parts (phases 0 to 6-quater, closed findings and issues, the MCP pilot, the trials of the full server) is in `ROADMAP_HISTORY.md`, moved verbatim: search it with `grep`, do not read it whole.

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
  | `edit.trim_common` | `others_names` | `other_sessions` | 6-quater, PR 9 |

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

#### Decisions of phases 6-ter and 6-quater still in force (November 2026)
The rest of these phases (the plans, the tables reviewed with Guido, the result of each PR) is in `ROADMAP_HISTORY.md`.
- *A server open to new recipes* (6-ter, Guido): a recipe added to the registry reaches the MCP server as a tool with no server code. The tool, its schema and its description come from the registry; nothing in the server names a recipe except the list of the excluded ones, each with its reason (the list is empty); the output of a tool is built from `RecipeResult` only; a test adds a fake recipe and a test checks that every recipe is a tool or excluded. What an agent needs and a recipe can give goes into the recipe or the `RecipeResult`, not into the server. Nothing is lost: every state is kept (`session_undo`, `session_redo`, `session_duplicate`) and `session_close` refuses to lose changes. The tools are not grouped by menu (the trials of 4.3: the cost of the 65 tools in context is accepted, both agents advised against grouping); a recipe silent for more than 20 s sends a heartbeat, and the GUI's Stop and `session_cancel` are one cancellation (`RecipeCancelled`).
- *A renamed parameter keeps its old name as an alias* (6-quater, PR 2): `ParamSpec.aliases`, silent and permanent (no warning, no version of removal: the archives are not changed once saved, F109, and their logs must replay; an alias costs one line). One translation, `RecipeSpec.canonical`, where the parameters enter a recipe, where the log records them (every log written has the new names) and when an archive's log is read. The old and the new name together are an `InvalidParameter`. `build_spec` refuses an alias that is the name or alias of another parameter.
- *A renamed class keeps its old name* (PR 6): `astrocook/_renamed.py`, `__getattr__ = renamed(__name__, {old: new})` in the module that defines it, the same class with a `DeprecationWarning` (for scripts only: archives and logs hold no class names); no version of removal before 3.0. No `V2` suffix in new names. `SessionV2` and the others are served the same way; `astrocook.core.session` and `session_manager` are lazy aliases of `astrocook.session`.
- *One set of names of the columns* (PR 8, F120): those of the GUI (`λ`, `λmin`, `λmax`, `F`, `dF`, the others by their own name; `core/column_names.py`), wherever a column is named (parameters, expressions, MCP tools); `lam`, `lammin`, `lammax` their ASCII spelling (the ASCII export); the names of the data (`x`, `y`, `dy`, ...) are silent and permanent aliases.
- *The units of the flux* (PR 8-bis, F129-F131): a loader assigns a physical unit to the flux only when the file says it (`BUNIT`, `TUNIT`, a known descriptor); otherwise `arb_flux` (`core/units.py`), an Astrocook unit that the archives and the FITS files read back. The unit is declared without converting by `yunit` of `set_properties`; `y_convert` converts and refuses an arbitrary flux. `export_ascii` writes ECSV with the units and `z_em`; a column without a unit (`cont`, `model`) is in the unit of the flux. A recipe that combines spectra converts to the units of its own spectrum or of the reference before comparing or scaling (`recipes/_combine.py`, `core/spectrum_combine.py`).
- *Parameters that name other sessions* (PR 9): kinds `session` and `sessions` of the registry, one rule in `recipes/_sessions.py`: names (in a script also the sessions themselves), a name not found is a `RecipeFailed` naming it, the log records names.
- *Refuse before the recipe runs* (PRs 7-10): an unknown component, a series not in `ATOM_DATA`, an unknown column or a value out of range is an `InvalidParameter` (the checks of the components in `recipes/_components.py`), not a result that looks valid; no default for a parameter that has no sensible one. A review PR does not change the numbers of a valid call: characterization tests first, fingerprints of the valid calls taken on master before the change.
- *Types* (PRs 3-4-bis): mypy gives zero errors on the whole package and new code adds none (rules in `CLAUDE.md`); ruff `B` is on (`zip(..., strict=True)` where the lengths agree by construction).
- *The menus come from the registry* (PR 11, F107, F125): the `menu` key of a recipe is the only place where the label of an entry is written (`"Edit > Set &Properties"`); the order and separators stay in `RECIPE_LAYOUT`; the `brief` is a sentence, used where a message says what a recipe did (log, status bar, MCP descriptions), not a label; every recipe has a `url` (a section of its tutorial) that a test checks. The tooltips (`RECIPE_TOOLTIPS`) are a third text (F151, phase 7-bis).
- *The descriptions of the MCP tools* (PR 12): `SHARED_PARAMS` in `mcp_server/schemas.py` is the only source of the descriptions of `session_id`, `provide_missing`, `max_listed` and `reason`, read by the tools of the recipes and of the sessions; `tests/test_mcp_schemas.py` caps them at 160 characters each, 24 000 in all over the 48 recipe tools and 300 for any description of any tool. A new description is shorter or says something of its recipe. The workflows of the instructions (`INSTRUCTIONS_WORKFLOWS`) are read by a test from the registry: a renamed tool or parameter fails it.

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
Done (November 2026), PRs #726-#747 and the closing PRs, release 2.6.0 (10 November 2026, PR 9, #751; on GitHub with the macOS app, not on PyPI). The full MCP server, a thin layer over the package: 65 tools, 48 generated from the registry of the recipes (none excluded, F106 made `export_ascii` one) and 17 of the server (the sessions with undo, redo and duplicate, the inspection, notes and gaps, the log as text, the cancellation, `gui_open`); the outputs built from `RecipeResult` only; progress notifications, a heartbeat and a cancellation shared with the GUI's Stop; the GUI following a session of the server (read only, with an editable copy); the headless Docker image with its CI job; trials with fresh agents (section 4.3). Found F106-F121. Numbers: the recipes unchanged but for `metals_auto` (PR 7-quater: it keeps the ions found, skips those outside the spectrum and at the bounds of the fit; the benchmark does not run it, a targeted test measures it); the mock benchmark of 2.6.0 identical to the baseline in all 19 rows; 3601 tests collected.
Left from this phase (not in a plan yet): a measure of the fit for the agents (residuals, chi2, AIC and BIC on a fixed set of pixels); the chi2 of the groups in the result of every fit; a redo that cancels the gap `undone` of the state undone; the publication of the image on GHCR (after #549); F112, F114, F115, F116 (register, 3.2).
Full text (plan, decisions, result of each PR and of the phase, trials): `ROADMAP_HISTORY.md`, "Phase 6-ter – Full MCP server" and "4.3 Full server: trials with fresh agents".

### Phase 6-quater – Review of the recipe parameters
Done (November 2026), PRs #752-#768 and the closing PR 13; released as 2.7.0 (PR 14, 13 November 2026). Named after 2.5.0 (Guido): the parameters and interfaces of the recipes reviewed recipe by recipe, adapted to each procedure, for the recipes that phase 7 does not rewrite (file, flux, edit, the manual fits; those of the continuum and the identification go with phase 7), with what the code review of 2.5.0 proposed (ruff `B`, the type checker, `Session` out of `core/`, the `V2` suffixes) and what 6-ter left (F107, the descriptions shared by the tools, the banner of a followed session). The plan was agreed with Guido on 11 November 2026 (PR 1); each PR had a conversation of its own, with the model chosen by the rule of `CLAUDE.md`. Full plan, tables reviewed with Guido, decisions and "Result of PR n" of each PR: `ROADMAP_HISTORY.md`, "Phase 6-quater".

**Plan of the phase** (estimates in working time of the session, before Guido's review; the estimates of PRs 1-8 were 4-8 times too long, recalibrated on PR 8):

| PR | Content | Model | Estimate | Done |
|---|---|---|---|---|
| 1 | Ruff `B` on, its 33 warnings fixed | Sonnet | 1-2 h | 11 Nov |
| 2 | Aliases of renamed names: one mechanism in the registry (the old name accepted by the logs, the scripts and the MCP tools), reused for the classes | Opus | half a day | 11 Nov |
| 3, 4, 4-bis | Type checker: mypy in CI on `recipes/` and `mcp_server/`, then `core/` and `fitting/`, then `gui/` | Sonnet | half a day each | 11 Nov |
| 5 | `Session` and `SessionHistory` out of `core/`, into `astrocook/session/` | Sonnet | 2-3 h | 11 Nov |
| 6 | The `V2` suffixes dropped (21 classes), with the aliases of PR 2 | Sonnet | half a day | 11 Nov |
| 7 | Review of file and flux (8 recipes) | Opus | half a day | 11 Nov |
| 8 | Review of edit on columns and expressions, F120 (one set of column names) | Opus | half a day to a day | 12 Nov |
| 8-bis | The units of the flux and the CSV files (F129-F131) | Opus | 3-4 h | 12 Nov |
| 9 | Review of edit on several spectra (9 recipes) | Opus | 2-3 h | 12 Nov, about 4 h |
| 10, 10-bis | Review of the manual fits (components and constraints), then of `optimize_system` and `fit_bayesian` with the samplers | Opus | 3-5 h (two PRs) | 12 Nov, about 3 h and 2 h |
| 11 | F107: the labels of the menus from the registry; F125: the pages of the docs for the *Help* button | Sonnet | 1 h | 12 Nov, about 2 h |
| 12 | MCP server: the shared descriptions shorter, the workflows in the instructions, the banner of a followed session | Sonnet | 1 h | 13 Nov, about 2 h |
| 13 | Closing of the phase: the review of the findings without a phase with Guido, "Result of the phase", the progress table, the rules in force, `CLAUDE.md`, the concluded detail to `ROADMAP_HISTORY.md` (documents only) | Sonnet | 1 h | 13 Nov |
| 14 | Release 2.7.0: the version, the CHANGELOG (the summary of the phase and the numerical results), the "preferred" entry of `switcher.json`; the slow tests and the full mock benchmark (`validate_mocks.py --compare`) on the final master, as for 2.6.0 (about 20 minutes, to be run in the background); the release and the macOS app published by Guido | Sonnet | 1 h of mechanical work, about 2 h with the CHANGELOG, the checks and the benchmark (the calibration of the phase: the mechanical PRs estimated at 1 h took about 2 with tests and documentation) | after PR 13 |

**Result of the phase** (November 2026, PRs #752-#768; the detail of each PR in `ROADMAP_HISTORY.md`):
- *What was done.* PR 1: ruff `B` on, 33 warnings fixed. PR 2: the old name of a renamed parameter an alias in the registry, silent and permanent (`ParamSpec.aliases`, `RecipeSpec.canonical`). PRs 3, 4, 4-bis: mypy on the whole package, from 840 errors (143 + 336 + 361) to zero, in CI (job `types`). PR 5: `astrocook/session/`. PR 6: the 21 classes without `V2`, the old names served by `astrocook/_renamed.py` with a `DeprecationWarning`. PR 7: file and flux (one parameter renamed: `folder` of `export_ascii`). PR 8: edit on columns and expressions, one set of column names (`core/column_names.py`, F120). PR 8-bis: the units of the flux (`arb_flux`, `core/units.py`) and `export_ascii` in ECSV with units. PR 9: edit on several spectra, one rule for the parameters that name other sessions (`recipes/_sessions.py`), the sessions combined in the units of a reference. PRs 10 and 10-bis: the manual fits, the checks of components, series and links in `recipes/_components.py`, the samplers of `fit_bayesian`. PR 11: the labels of the menus from the registry (`menu` key) and the `url` of every recipe. PR 12: the descriptions of the parameters shared by the tools -34% (`SHARED_PARAMS`), a paragraph of workflows in the instructions, pid and version in the banner of a followed session. PR 13: this closing.
- *Numbers.* 31 findings found (F122-F152), 27 closed (24 fixed, F147 documented, F126 and F128 accepted) and 4 open (F148 and F149 for phase 7, F151 and F152 for phase 7-bis); closed here also F107, F117, F120 and F125, brought from 6-ter. Ruff `B`: 33 warnings to zero. Mypy: zero errors on 93 source files. Tests: 3601 collected at 2.6.0, 4060 now, none failing (the two expected failures, F17 and F28, are strict xfails).
- *Numerical result: unchanged on the benchmark, by construction.* No PR touched the continuum, the identification, the fitting or the automatic pipelines, so the mock benchmark of phase 4 cannot move; the full run on the final master (PR 14) is identical to the baseline in all 19 rows. PR by PR: PRs 1, 3, 4, 4-bis (lint and types: `assert`, `cast`, no new branches), 5 and 6 (modules and names moved; one snapshot of `v1_archive` regenerated, whose `repr` carried the old class name), 11 and 12 (texts, labels, addresses, a banner): nothing numerical. PR 2: names only. PRs 7, 8, 9, 10: the values of a valid call unchanged, checked by the fingerprints and snapshots of `test_valid_results_unchanged` taken on master before each change (the tolerances in the history); what changed is that invalid values are refused before the recipe runs. Changed on purpose (all agreed with Guido, none read by the benchmark): the units of the flux (PR 8-bis: a loader assigns a physical unit only when the file says it, otherwise `arb_flux`; the values do not change, and `cont` converts with the flux, F130); the spectra combined in different units (PR 9: `coadd`, `equalize_and_stitch` and `equalize` convert before scaling, F132-F134 and F138; the same results with the sessions in one unit; `equalize_order=-1` of old logs is read as None); the samplers of PR 10-bis (nested sampling no longer crashes and its posterior is weighted, F139 and F140; MCMC unchanged up to 16 free parameters, `max(32, 2 * ndim)` walkers, F141; no pipeline and no slow test call `fit_bayesian`, a targeted test measures it).
- *Found along the way, reported at once:* F135 (the recipes that read other sessions crashed in a script without a GUI) and F139 (the crash of nested sampling), both fixed.
- *Left.* Phase 7 (the continuum and the identification; F148 and F149 for `optimize_system`) and phase 7-bis (small fixes); the review with Guido of the findings without a phase (PR 13) is in the register, 3.2.
- *Calibration.* PRs 1-7 took one working day in all; a mechanical PR about 1 h estimated and 2 h with tests and documentation, a review PR 2-4 h. The 14 PRs took about three working days of the sessions, against the 6-9 days of the first plan.

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

Findings assigned to phase 7 (register, 3.2): F51, F52 and F116 (with the flags and the uncertainties of the components), F53, F54, F56, F57, F58, F59, F68, F72, F84, F93, F102 (`forest_auto`), F112 (the output of the pipelines), F118, F119, F121, F148, F149 (assigned in the review of PR 13 of 6-quater, Guido, 13 November 2026).

### Phase 7-bis – Small fixes (proposed)
Agreed with Guido in the review of the findings of PR 13 of 6-quater (13 November 2026); not started, its plan to be written at its start, and the order with phase 7 to decide (it does not wait for the science). No algorithm of the continuum or the identification. Assigned:
- F101: the script console is kept and `AppController.run_recipe` corrected (it passes `session_override`, which `_on_recipe_requested` does not accept: every call fails);
- F152: the model redrawn after `update_component(s)` and `update_constraint` without a fit, if it costs little, or the warning "please refit" turned on again (the GUI behaves as the server today: the warning exists but `_suppress_refit_warning` is True by default);
- F151: the tooltips of the entries from `details` or a `tooltip` key of the registry;
- F115: the methods text (Guido's text: the gaps out of it, the name of the spectrum, the reasons, `session_methods` as a file);
- F114: the `[bayes]` extra in the Docker image (size to measure) or the description of the tool saying it;
- F28: the Bayesian samples saved in the archive (an HDU or an `npz`);
- F30: the edges of `running_mean` and `running_std` (they pad with zeros; the numbers of `calculate_running_std` change at the edges: Guido's agreement first);
- F15: the scale factor of `equalize_and_stitch` with narrow absorption, when the X-Shooter VIS/NIR data arrive (without them, closed as not verifiable);
- the comment of `_NE_GLOBAL_DICT` (`core/spectrum_expressions.py`), which says that `median`, `mean`, `std` work in NumExpr (F128: they do not; no change of behaviour).

### Code quality (review of 2.5.0)
A review of the code of 2.5.0 (November 2026, after the question of how
the package compares with similar codes), with measurements, and where
each proposal could go. **Proposals**; the three assigned to 6-ter
were moved to the start of 6-quater at its closing (Guido) and done there
(the table below); the others are not in a plan yet.

*Measurements* (2.5.0, before phases 6-ter and 6-quater changed the package):
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
| A type checker (mypy or pyright) in CI, gradually: `recipes/` and `mcp_server/` first, then `core/` and `fitting/`; the GUI last, if ever | 6-quater (moved from 6-ter, PR 8; done, the whole package, 11 November 2026) | the server builds the tools from the typed parameters: checked types serve the agents first |
| `SessionV2` (and `SessionHistory`) out of `core/`, e.g. into `astrocook/session/`, the old import kept as an alias | 6-quater, PR 5 (done, 11 November 2026) | the headless server works on the sessions; the core becomes the bottom layer |
| Ruff with `B` (flake8-bugbear), a PR of its own; the style rules (`E`, `W`) left out | first PR of 6-quater (moved from 6-ter, PR 8; done) | small and mechanical; finds bugs, not style |
| The `V2` suffixes of the class names (`SpectrumV2`, `RecipeAbsorbersV2`, ...) dropped, the old names kept as aliases | 6-quater, PR 6 (done, 11 November 2026) | changes the interface of the scripts, as a renamed parameter: the same mechanism of aliases |
| The complex functions of the recipes that phase 7 does not rewrite (e.g. `RecipeEditV2.equalize_and_stitch`, 31) | 6-quater, PR 9 (done for `equalize_and_stitch` and `stitch`, 12 November 2026) | the recipes are reviewed one by one there |
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
- **2.6.0, after phase 6-ter (10 November 2026):** the full MCP server (65 tools, 48 generated from the registry, sessions with undo and redo, the log as text, cancellation, the GUI following a session, a headless Docker image); F106-F121. Numbers: unchanged on the benchmark (identical in all 19 rows) except `metals_auto` (not run by it: a targeted test). Published on GitHub with the macOS app; not on PyPI (#549).
- **2.7.0 (13 November 2026, PR 14), after phase 6-quater:** the review of the recipe parameters: the aliases of the renamed parameters (`folder` of `export_ascii`, ...), the classes without `V2` (the old names still served, with a `DeprecationWarning`), mypy on the whole package in CI, one set of names of the columns, the units of the flux (`arb_flux`) and the CSV files in ECSV, the recipes refusing invalid values before they run, the labels of the menus from the registry, the shorter descriptions of the MCP tools; F122-F152, F107, F120, F125. Numbers: unchanged on the benchmark (the full run on the final master, 135 mocks, 1884 s, identical to the baseline in all 19 rows; slow tests 5 of 5); changed on purpose, and said in the notes, the units of the flux, the spectra combined in different units and the samplers of `fit_bayesian`. A minor version (the interface of the scripts changes: names, refusals). Published by Guido on GitHub with the macOS app, not on PyPI.

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

The rows here are the open ones, each with its phase (reviewed with Guido at the closing of phase 6-quater, 13 November 2026: phase 7 for the continuum, the identification and the components; phase 7-bis for the small fixes).

| ID | Area | Finding | Impact | Recorded by | Status |
|---|---|---|---|---|---|
| F15 | `equalize_and_stitch` | The scale factor is a median of the ratio of the two arms in the overlap: with narrow absorption (e.g. tellurics) seen at two different resolutions, it is biased by ~1-1.5% in a mock with equivalent-width-conserving tellurics around 1000 nm (no bias without tellurics). Possibly related to #477; to be checked on real X-Shooter VIS/NIR data | numerical | — (mock in the PR that added this row) | open: phase 7-bis, when the X-Shooter VIS/NIR data arrive (Guido, 13 November 2026); without them by then, closed as not verifiable, the bias documented |
| F28 | `.acs2` Bayesian samples | The posterior samples stored by `fit_bayesian` in the system list metadata (`bayesian_results`, used by the corner plot) are dropped on save without a warning: the FITS header cannot hold them. The fitted values and errors are kept | data loss (recomputable) | `test_bayesian_results_are_saved` (strict xfail) | open; since phase 6 (PR 7) a WARNING says at the save that the samples are not saved (astropy skipped them with a warning of its own); phase 7-bis (Guido, 13 November 2026): the samples in the archive |
| F30 | `running_mean`, `running_std` | Both use `np.convolve(mode='same')`, which pads with zeros: the first and last h values are biased towards 0 (4/7 of the true mean at the edge for h = 3). V1 `running_mean` repeats the first and last full-window value and ignores NaNs. Used by `find_absorbed_regions` on each interval, with half-windows of hundreds of pixels (5000 km/s / 8 on J0103-1305). Replacing it with the V1 function raises the agreement of the forest mask with V1 from 92.7% to 95.7% and the 5th percentile of the continuum ratio from -4.7% to -3.0% (median unchanged, -1.2%): a partial cause of F17 | numerical, V1/V2 difference | `test_running_mean_and_std_pad_the_edges_with_zeros` (characterization) | fixed in phase 4 for the continuum (#465): `find_absorbed_regions` no longer uses `running_mean`, but a masked mean over a window of the smoothing length in km/s on the pixel grid, exact at the edges (`_masked_window_mean`). The edge bias was worth 0.81 → 0.93 within 2000 km/s of the edges of the forest interval on a mock at S/N 30. The plain normalized convolution of the earlier branch (window in pixels of the clipped array) failed at S/N ≥ 300 (forest continuum up to 1.66 times the true one: the window grows as pixels are clipped, over wavelengths where the continuum varies by more than its noise); the window in km/s does not. `running_mean` and `running_std` themselves still pad with zeros (used by `calculate_running_std`): open, minor. The rest, the zero padding of the two functions: phase 7-bis (Guido, 13 November 2026); it changes the numbers of `calculate_running_std` at the edges, with his agreement first |
| F51 | `identify_lines` (forest) | Every absorbed region of the forest without a metal candidate is assigned to Ly-a (score 0, no test), even when Ly-a is not among the requested `multiplets`: `multiplets='CIV'` on the full-chain mock gives 233 Ly-a and 6 CIV candidates. In the MCP pilot the Ly-a candidates crowded out the metals in the output (833 on J0103-1305). Found in the MCP pilot (report) | behaviour | — | open: phase 7 (Guido, 13 November 2026) |
| F52 | `VoigtFitterV2` (uncertainties) | Components that end at the bounds of a fit get uncertainties far beyond any meaning, stored as results without a warning: in the `doublets_auto` fit of J0103-1305 (MCP pilot), four CIV components at z = 3.927 with logN = 10.0 and b = 1.0 km/s have dz up to 131, dlogN up to 3800 and db up to 9e4 km/s; 21 of the 101 new components have dlogN > 1 or db > b. Related to F20, F21. Found in the MCP pilot (report) | numerical, error handling | — | open: with [#517](https://github.com/DAS-OATs/astrocook/issues/517); since its PR 1 the components at a bound carry a flag (`ComponentDataV2.flags`), so their errors can be told apart; phase 7 (Guido, 13 November 2026), with F116 |
| F53 | continuum (forest at high z) | On an ESPRESSO spectrum of J0959+1312 (z_em ~ 4.06), `estimate_auto` flags 100% of the pixels between 525 and 600 nm as absorbed, so the continuum there is a straight interpolation between the edges of the interval (flat, above the 99th percentile of the flux, ignoring the wing of the Ly-a emission); with `template=True` it follows the template, still with no free pixel, and lies below the flux peaks at 565-590 nm (6% of the pixels above C + 2 sigma). Reported by the agent in Guido's MCP trial; not reproduced here (the spectrum is not in the repository) | numerical | — | open: phase 7 (emission model); with #465 |
| F54 | identification (tellurics) | No telluric masking: the O2 A band (759-770 nm) gives spurious doublets with good scores, CIV at z = 3.927 in J0103-1305 (4 components at the bounds, F52) and CIV z = 3.928 with MgII z = 1.733 in J0959+1312 (fitted together, reduced chi2 12). Found in the MCP pilot (trials) | behaviour | — | open: phase 7 (identification); telluric data in the pre-release `data-telluric-2026` |
| F56 | continuum (forest) | Measured on the mock benchmark (phase 4 baseline, 135 mocks): `estimate_auto` gives a forest continuum (Ly-b to Ly-a) low by 6.4% in the median, 13% at S/N 10, 12% at z_em = 4.5, with 5% of the forest pixels more than 18% low (34% at z_em = 4.5); the red side is within 0.2%. Quantifies F17 and the known low forest continuum (F53 on real data) | numerical | `benchmarks/validate_mocks.py`, baseline | improved by #465 (measured on the mock benchmark (2.1.0 + F49/F50 → #465): forest continuum 0.936 → 0.943 (5th percentile 0.820 → 0.830), red side unchanged (1.002); S/N 10: 0.875; S/N 5: 0.80, the limit of clipping). Also found: spikes falling on absorbed pixels and raised to the continuum level become anchors and raise the forest continuum slightly, i.e. the bias comes from the shortage of anchors at the continuum level. Open: phase 7 (emission model) |
| F57 | identification (completeness) | On the mock benchmark the metal doublets are found only well above the detection limit: none of the systems whose weaker line is below 3 sigma, 5% of those at 3-5 sigma, 30% at 5-10 sigma, 91% above 10; completeness over the systems above 3 sigma 0.85 for CIV, 0.69 for SiIV, 0.91 for MgII. Purity is high (0.97-0.99) | behaviour | `benchmarks/validate_mocks.py`, baseline | improved by #465 (wider window on the red side, 400 km/s instead of ~100: weak lines are flagged better): completeness 3-5 sigma 0.03 → 0.10, 5-10 sigma 0.35 → 0.44, CIV 0.857 → 0.874. Why: `find_absorbed_regions` is both the continuum step and the line detector (`identify_lines` takes its regions from `abs_mask`); a weak line pulls down the running mean of a narrow window by a good fraction of its own depth, so its pixels do not reach 2 sigma; with a window 4 times wider the line weighs 4 times less, and purity is unchanged (clipping still at 2 sigma). Direction for phase 7: decouple the detection from the continuum, e.g. search the lines on the flux normalised by the final, smooth continuum, with a filter matched to the expected line width instead of a per-pixel threshold. Open: phase 7 |
| F58 | `doublets_auto` (decomposition) | On the mock benchmark the systems with several components are decomposed differently from the truth: where one fitted component stands for several true ones the pulls of the matched component are 4-6.5 wide, and spurious components appear next to real ones (F49); the total column of each system is right (median difference -0.002 dex, pull width 1.0-1.1 at R >= 45000). The pipeline adds one component per candidate and never adds components from the residuals | behaviour | `benchmarks/validate_mocks.py`, baseline | open: phase 7 (Guido, October 2026; confirmed 13 November 2026) |
| F59 | fit at low resolution | At R = 8900 the total column of the metal systems is off by 0.39 dex (spread) with a pull width of 2.5, against 0.03-0.07 dex and ~1 at R = 45000-140000. Possible causes: unresolved saturated components (physics), and the convolution kernel computed at the median wavelength (rest of F10, up to +53%/-26% in width at the edges of a window) | numerical | `benchmarks/validate_mocks.py`, baseline | partly fixed by #526 (PR 3, exact convolution): at R 8900 the spread of the system column goes from 0.175 dex (master after #465; 0.39 in the 2.1.0 baseline) to 0.102, the median reduced chi2 from 2.45 to 1.61 and the flagged components from 30 to 16; R 45000 and 140000 unchanged (0.046 and 0.015 dex). After F61: 0.095 dex at R 8900, 0.041 at R 45000. What is left at R 8900 (median reduced chi2 1.6) is probably the physics of unresolved saturated components; to be looked at in phase 7 |
| F68 | `identify_lines` (`Ly_ab`) | The Ly-a/Ly-b pairs of the forest are never identified: `STANDARD_MULTIPLETS['Ly_ab']` lists the red line first (`['Ly_a', 'Ly_b']`), while `rate_doublet_candidate` expects the blue line first, as for the metal doublets (`CIV_1548`, `CIV_1550`): the two windows in velocity never overlap, and the score is always 0 (checked on clean mocks, logN 13-14.5). `Ly_ab` is in the default list of the identification recipes, yet the Ly-b line of every pair becomes a spurious Ly-a system (z = 2.29 for a pair at z = 2.90); with `bypass_scoring` the pair is recognised. Found in phase 5 while writing the characterization tests of the line methods | identification | `test_ly_ab_pairs_are_never_rated` | open; phase 7 (Guido's decision: simple to fix, not urgent; it changes the identifications, measured on the mock benchmark) |
| F72 | `rebin`, `resample_on_grid` | The boolean masks (e.g. `abs_mask`, the absorbed regions) are dropped with a warning: `get_resampled_column` interpolates only numbers (`rebin` catches its `TypeError`, `resample_on_grid` checks the type first). After a rebin or a resampling the mask of the user is lost and the pipelines make a new one. A mask could follow the new grid (e.g. a pixel masked if any old pixel in it was). Found in phase 5 with F71 | behaviour | `test_rebin_drops_the_masks`, `test_resample_drops_the_masks` | open; phase 7 (agreed with Guido: it changes the results and needs a choice of method, how a mask follows a new grid, on the inputs of the continuum; only the recipes *Rebin*, *Resample* and *Co-add* use these methods, not the pipelines nor the mock benchmark; meanwhile, *Find absorbed regions* after a rebin) |
| F84 | x axis in velocity | The spectrum does not keep the reference redshift of a conversion of the x axis to velocity (`x_convert`; only the history line "X-axis converted to km/s at z_ref=..." records it): every conversion back to wavelength (rebin, smooth, continuum, combination of spectra, the ticks of the systems since F66) uses z_em, so they are wrong if the axis was converted with another z_ref, or if z_em changed after the conversion. Found in phase 5 while fixing F66 | design | — | open: phase 7 (Guido, 13 November 2026): keep z_ref with the spectrum |
| F93 | `identify_lines` (identifications) | The identification is by doublets only: two absorbed regions are a candidate when they overlap once one is shifted by the separation of the doublet in velocity (the ratio of the wavelengths), and the pair is accepted when the R^2 of F_1 = F_2^gamma (gamma from f and lambda of the two lines) is above 0.5 (metals only redward of Ly_a). Three things found while rebuilding the Identifications window: (1) the pairs are not stored: each region keeps only its candidates ([transition, score]), and the window finds the pairs again from the overlap of the redshift ranges and the shared score; (2) the multiplets with more than two lines (FeII, SiII, NI, ...) are never searched, only those with exactly two members; (3) the regions of the forest that no doublet identified are labelled Ly_a with score 0, an assignment and not an identification (lya_auto calls the recipe with "Ly_a" alone, not a doublet, so all its regions are so labelled). Found in phase 5-bis | behaviour | `test_identification_viewer.py` (the window) | to be discussed in phase 7 (identification, Guido's algorithmic choices): store the pairs with their redshift and score; the window shows the doublets as the recipe finds them and the assumed Ly_a as such |
| F101 | script console (`gui/scripting.py`) | `AppController.run_recipe`, the `app.run_recipe(...)` of the scripts of the console, calls `MainWindowV2._on_recipe_requested` with a keyword `session_override` that the method does not take: every call fails with a TypeError. Found in phase 6-bis (PR 3a) while adapting the callers of the recipes | bug | — | open: phase 7-bis (Guido, 13 November 2026): the console is kept and the call corrected |
| F102 | `lya_auto`, `forest_auto` (cleanup of the narrow lines) | After removing the lines narrower than `min_b`, `_remove_narrow` calls `refit_all` without `region_limit`: every group of the session is refitted, also the components outside the forest (e.g. a CIV doublet on the red side fitted before), whose values change slightly. Found in phase 6-bis (PR 4b) while recording the reasons of the removals. `forest_auto` is not limited to a region at all: it identifies on the whole spectrum and its first refit, before the cleanup, also refits every component of the session | behaviour, numerical | `test_lya_auto_does_not_refit_the_red_side` | `lya_auto` fixed in phase 6-bis (PR 4c, Guido's decision): the refits of the cleanup limited to `lya_forest`, as the pipeline; documented in the tutorial (*Refit All* fits everything again). The mock benchmark does not run `lya_auto`: a targeted test replaces it. `forest_auto`: open, Guido's decision (limiting it means choosing its region, e.g. `all_ly_forest`, for the identification and both refits); `forest_auto`: phase 7 (Guido, 13 November 2026) |
| F112 | MCP output of the pipelines (`IdentificationSummary`) | The result of `absorbers_doublets_auto` says `identification.n_systems: 0` and `systems_by_series: {}` while its notes say "3 candidate systems found" and three components are fitted (65 on J0103-1305): the summary counts the candidate systems of the metadata of the spectrum, which the pipeline does not leave there. Both agents of PR 7 noted the contradiction | output | — | open: phase 7 (Guido, 13 November 2026), with the output of the pipelines: what they leave in `RecipeResult.extra` |
| F114 | Docker image | `absorbers_fit_bayesian` is a tool of the image, but the image does not install the `[bayes]` extra (emcee, dynesty): the tool fails with an ImportError (recorded as a gap `failed`). Found in trial 1 of PR 7 | packaging | `tests/test_docker_image.py` (the list of the image) | open: phase 7-bis (Guido, 13 November 2026): the extra in the image (size to measure), or the description of the tool saying it |
| F115 | methods text (`methods_text`) | Seen by the agent of trial 1 of PR 7: the declared gaps (meant for the developers) are in the text for the paper, as written (decision of PR 5c); the first sentence names the spectrum by the name of the session ("The spectrum of spectrum"); the reasons open each paragraph as written, which reads oddly in a paper; `session_methods` gives the text only inside its JSON, so the agent wrote it to a file with its own Python (declared as a gap) | output, cosmetic | — | open: phase 7-bis (Guido, 13 November 2026): the text is his, PR 5c |
| F116 | flags of the components | No flag for components unresolved (b well below the instrumental b: 1.4-2.8 km/s at R = 45000 in J0103-1305) or saturated (MgII at logN 15.06 ± 0.24), which both agents of PR 7 judged unreliable by hand; the flags of #517 cover the bounds and the window only | behaviour | — | open: phase 7 (Guido, 13 November 2026), with F52 (the flags of #517 cover the bounds and the window only) |
| F118 | `optimize_system` (stop criterion, output) | On the mock of trial 3 (PR 7) the optimizer added nothing to a single CIV that hides two components 9.4 km/s apart, where a second component added by hand lowers the chi2 on fixed pixels from 55.4 to 37.1 (dAIC -12): it tries a component only where a residual pixel exceeds `threshold_sigma` (2.5), and the deepest was 2.1 sigma (the reason is an INFO line of the log, not in the result: the agent saw "1 component changed" and nothing else). A systematic misfit spread over many pixels is not seen. Related to F58 (the pipelines never add components from the residuals) | behaviour, output | — | open: the criterion to phase 7 (a heuristic, Guido's); decided (Guido, on PR 7): the reason of the stop (e.g. the deepest residual against `threshold_sigma`) and the AIC of each component tried in `RecipeResult.extra`, numbers unchanged, fixed in PR 7-ter (`details.optimization`: `start_aic`, `tried`, `stopped`), the criterion open: phase 7 |
| F119 | `doublets_auto` (scoring of the candidates) | The SiIV doublet of the z = 2.2 system of the mock of trial 3 (logN 12.9, b 9, S/N 50; depths 0.55 and 0.33) is found as a kinematic candidate but scored R^2 = 0.483, below the threshold of 0.5, and dropped; 0.485 with the true continuum, so the continuum (on the wing of the NV emission) is not the cause. The candidates dropped are not in the result: the agent found the SiIV by eye | behaviour, output | — | open: the scoring to phase 7 (identification); decided (Guido, on PR 7): the candidates dropped, with their score, in `RecipeResult.extra`, fixed in PR 7-ter (`details.dropped_candidates`: series, z, score, threshold, x range; not the forest doublets, which no score accepts), the scoring open: phase 7 |
| F121 | `metals_auto` (acceptance of a trial ion) | `metals_auto` kept a trial ion if its fit gave logN >= 12 and 2 <= b <= 50 km/s; `aic_penalty` was accepted and not used, while its description, the `details` of the registry and the methods text said that an ion is kept if it decreases the AIC by it. Once the ions found were kept (F111), the criterion kept 192 ions of 239 on J0103-1305, 67 of them with b = 50 km/s (the bound), 9 with logN = 19 (the bound), 30 with an AIC change above -5; on the mock of trial 1 a spurious SiIV at z = 1.83 (logN 12.80, b = 50.0) | behaviour | — | partly fixed (Guido's choice, PR 7-quater): an ion is not kept with b or logN at the bounds of the fit (within 0.1 per cent) nor with an AIC change above `-aic_penalty` (against the same model without it, on the pixels of its fit), and each ion tried is reported with the reason (`tried`). Open (phase 6-quater or 7): on the mock of trial 1 a spurious CII at z = 1.83 in the forest is still kept (logN 12.58, b 18, AIC change -33); a trial ion that absorbs a Lyman-alpha line improves the AIC anyway, and the model without it keeps the neighbours refitted with it (a bias towards keeping): the comparison with the alternative identification (Lyman-alpha) is for phase 7; measured on the mocks (the benchmark does not run `metals_auto`). Measured on J0103-1305 raw with the code of master after #747 (November 2026; `edit.set_properties` z_em = 4.07, R = 45000, `continuum.estimate_auto`, `doublets_auto`: 65 components, then `metals_auto` on the six default ions, through `recipes.result.run`): 320 trials built (64 redshifts), 81 not tried for their lines outside the spectrum (as before), 239 fitted, **82 kept** (CII 12, CIV 5, NV 17, OVI 18, SiIV 30; against 192 with the old criterion, 67 of them at b = 50 km/s; none now), 157 not kept: b at its bound 108, logN below 12 29, AIC change above -5 7, b below 2 km/s 7, logN at its bound (19) 6. Time: the trials 386 s (6.4 min) and the final `refit_all` 121 s, 507.5 s the recipe, against 26.7 min for the 239 trials of the old criterion; the whole chain 8.9 min; 147 components at the end. The run interrupted at 70 per cent with the rules of F121 kept 60 of 182 tried (33 per cent), the complete one 82 of 239 (34 per cent); the rest: phase 7 (Guido, 13 November 2026) |
| F148 | `optimize_system` (window) | The window of the residuals and of the AIC is centred on the first line of the series of the seed (`lam_target_nm`): for a multiplet (CIV) the other lines are out of it (CIV 1550 is 500 km/s from CIV 1548, with a window of ±100 km/s), so a component is tried and judged on one line only, while the fits use them all. Found in PR 10-bis (phase 6-quater), not measured | numerical | – | open: phase 7 (components; changing it changes the results), registered only (Guido) |
| F149 | `optimize_system` (failed fits) | A fit that does not converge within `max_nfev` (`success` False) gives an infinite AIC: its candidate is rejected without the reason (the `tried` of the report says "AIC changes by inf"), and if the starting fit fails, no candidate can be added (inf - inf is NaN, never below `-aic_penalty`). `max_nfev` was fixed at 1000; since PR 10-bis it is a parameter. Found in PR 10-bis (phase 6-quater) | numerical, interface | `test_a_fit_that_reaches_max_nfev_adds_nothing` (characterization) | open: phase 7 (Guido: registered only); `max_nfev` documented |
| F151 | GUI, the tooltips of the entries | The tooltips of the entries of the recipes (`RECIPE_TOOLTIPS` in `gui/main_menus.py`) are a third text for the same recipe, close to its `details` in the registry. The messages that say what a recipe did or why it cannot run (status bar, followed session, continuum required, `RecipeResult`, prerequisites, MCP descriptions) name it by its `brief`, a descriptive sentence: decided with Guido to leave them so (a menu label is for saying where to click, not what happened). Found in PR 11 of 6-quater | cosmetic | — | open: phase 7-bis (Guido, 13 November 2026): the tooltips from `details` or a `tooltip` key |
| F152 | `update_component(s)`, `update_constraint` (model not redrawn) | These recipes change the parameters (or the links, freedom) of components but run no fit and do not recompute the `model` column: after `update_component(logN=14.5)` the model of the session is identical to the one before (checked on the mock of `tests/full_chain.py`; same for `update_constraint`), until the next fit; the components table shows the new value while the plot and `spectrum_read` show the old model. Found while writing the workflows of the instructions of the MCP server (phase 6-quater, PR 12), which now say that these tools do not fit and must be followed by `fit_component`. Not changed: it may be by design (the GUI edits then fits); to decide with Guido whether the model should be redrawn | behaviour | — | open: phase 7-bis (Guido, 13 November 2026): to improve; the model redrawn from the new parameters without a fit, if it costs little, or the warning turned on again. The GUI behaves as the server: editing a table runs only `update_component`, and a dialog ("please refit the components") says so, but `_suppress_refit_warning` is True by default (`gui/main_window.py`), so the new value shows and the model does not change until the next fit |

Closed findings (fixed, documented or accepted: all the others, F1 to F152 not listed above) are in `ROADMAP_HISTORY.md`, "3.2 Findings register", one row each, to be found by number (`grep "^| F52 "`).

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
- **Full server: phase 6-ter, after phase 6-bis and 2.5.0 (done, release 2.6.0; trials in 4.3):** all the recipes, structured outputs,
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
November 2026, after PRs #726-#742 (64 tools): three trials with fresh agents (Claude, no knowledge of the code) through a minimal MCP client (`benchmarks/mcp_trials/`), on the full-chain mock in the Docker image, on a real spectrum and on the editing tools with undo. The chain works at the first try; with undo, redo and the constraints an agent reached the true model of a blend, step by step in the log; what it missed was a measure (residuals, fit statistics on fixed pixels) and two descriptions that work (F117, F120); the automatic steps fell short (F118, F119). The trials confirm the direction of phase 7: the agents met the limits of the recipes, not of the interface (F54, F58, F111, F116, no upper limits, no z_em estimator, no telluric model). Full text of the trials: `ROADMAP_HISTORY.md`, "4.3 Full server: trials with fresh agents".

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
recipe of its own; the gaps the agents meet (phase 6-ter, decision 6, in
`ROADMAP_HISTORY.md`) say which recipes are missing.

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
| 6-ter – Full MCP server | done (November 2026), PRs #726-#747 and the closing PRs; release 2.6.0 published on GitHub on 10 November 2026 (PR 9, #751): 65 tools, sessions with undo, the GUI following a session, the Docker image; F106-F121 found. The long row: `ROADMAP_HISTORY.md`, "7. Progress" |
| 6-quater – Review of the recipe parameters | done (November 2026), PRs #752-#768 and the closing PR 13; released as 2.7.0 on 13 November 2026 (PR 14): ruff `B`, aliases of renamed names, mypy at zero errors, `Session` in `astrocook/session/`, no `V2` suffixes, the review of file, flux, edit and the manual fits, one set of column names, the units of the flux, the menus from the registry, the shorter MCP descriptions; F122-F152 found, 27 closed. The long row: `ROADMAP_HISTORY.md`, "7. Progress" |
| 7 – Robust components | not started |
| 7-bis – Small fixes | proposed (13 November 2026): F101, F152, F151, F115, F114, F28, F30, F15 and the comment of F128; not started |
