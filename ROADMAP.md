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

#### The reference of the benchmark (phase 7, PR 7c.1; the plan of 7c, Guido, 16 November 2026)
From PR 7c.1 the reference of the mock benchmark, on both grids (`full` and `robustness`), is `benchmarks/results/validate_mocks_robustness_baseline.json`, made by the PR 7a.2 and regenerated after the fix of F68 (the scoring of the doublets in the logarithmic frame of their lags). `validate_mocks.py --compare` prints next to it the numbers of phase 4 and the metrics of phase 7 of a run of either grid; `--save-reference` rewrites the entry of the grid run. The baseline of 2.2.0, `validate_mocks_baseline.json`, is kept as history (`--compare-baseline`), with that of 2.1.0. A PR that changes the results on purpose regenerates the reference of both grids and says so; one that must not change them shows `--compare` identical. The Ly-a/Ly-b pairs are measured apart, by `benchmarks/validate_lyab_pairs.py`: the mocks of `validate_mocks.py` start at 104 nm in the rest frame, above the Ly-b emission, and hold no Ly-b line.

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
- *The menus come from the registry* (PR 11, F107, F125): the `menu` key of a recipe is the only place where the label of an entry is written (`"Edit > Set &Properties"`); the order and separators stay in `RECIPE_LAYOUT`; the `brief` is a sentence, used where a message says what a recipe did (log, status bar, MCP descriptions), not a label; every recipe has a `url` (a section of its tutorial) that a test checks. The tooltips are the `details` (F151, B3).
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

### Phase 7 – Robust components: emission model, identification, fit
Agreed with Guido at the start of the MCP pilot (October 2026); scope and plan agreed on 13 November 2026, after 2.7.0. The goal of section 5, a complete fit of a spectrum in a single operation, rests on its components: today the continuum works "well enough" but is not robust, gives no uncertainty and fails in some cases (DLAs, among others); the identification is fairly accurate, but it cannot account for *all* the absorption; the automatic fits miss components that an agent finds by eye. The phase redesigns them, and takes in every limit of the recipes that the agents met in the trials (sections 4.2 and 4.3: "the agents met the limits of the recipes, not of the interface") and the open directions of section 2:
- **Emission model:** continuum (and emission lines) with an uncertainty, robust to strong absorbers (DLAs, BALs) and to the forest, and an estimate of z_em;
- **Absorption identification:** every absorbed region accounted for, identified or explicitly unidentified, with alternatives and a credibility, tellurics recognised;
- **Components and fit:** flags and uncertainties that say when a component is unreliable, components added where the residuals ask for them, a measure of a model on fixed pixels to compare alternatives, upper limits.

**Scope: phase 7 against 7-bis** (Guido, 13 November 2026). Phase 7 takes what changes *what* Astrocook computes or decides (algorithms, criteria, new capabilities of analysis), each PR measured on the benchmark; phase 7-bis what changes *how* the results are exposed, stored, packaged or shown, and the bugs outside the algorithms (no new algorithm; the numbers of the pipelines unchanged, except a fix stated as such). By this rule F72, F84, F102, the immediate fix of F112 and the reason in the report of F149 moved to 7-bis; F68 (a bug that changes the identification) stays here, as the first PR of 7c. One phase, in four blocks; the reference paper (section 6) waits for its end, or, if it must come out first, takes what is done by then (Guido).

**Method** (unchanged from the first plan of the phase): for each component a design document first (options, costs, what the literature does), in `design/` (a folder of the repository, not of the published docs), read and decided by Guido; then a prototype measured on the benchmark, next to the old algorithm and selectable, until the benchmark shows it better; then it becomes the default, with the change stated in the release notes. After #465 the automatic continuum is an iterative clipping with three measured patches (masked running mean in km/s, rejection of upward spikes with two thresholds, local S/N); Guido would like a simpler approach with a statistical basis to replace it (directions: spikes recognised by their width, narrower than the resolution element; the continuum from the upper envelope of the flux or from the mean absorption expected at each redshift; the 0.055 sigma bias of one-sided clipping, F60). Proposed (code review of 2.5.0, to be agreed with the first design document): a ceiling on the complexity of the new code, no function above 20 (`radon cc`), among the criteria of the new algorithms.

**Data and blindness.**
- The mock benchmark of phase 4 keeps its grid and baseline (`validate_mocks.py --compare` stays valid); the new cases (DLAs, tellurics, emission lines, BALs) are a second grid of their own (7a.1).
- Real spectra (Guido, 13 November 2026; they replace the 3-5 proprietary spectra agreed earlier the same day): the 28 UVES quasars of D'Odorico et al. (2010) (z_em 2.14-4.56, S/N per pixel 20-205, pixels of 1.3-2.5 km/s), in the reductions of SQUAD DR1 (Murphy et al. 2019, MNRAS 482, 3458), not in those of the paper: the difference is stated, and its effect on the comparison measured (7a.4). The reference is the published analysis (its CIV and SiIV components), not a truth: on these spectra the benchmark measures the agreement with it, and with the DLAs of the SQUAD catalogue; the continuum has no reference there (the mocks measure it). The spectra are public, but the blindness of the test part is a discipline, not an access control (Guido, 13 November 2026): the spectra and every per-object result live in the private `DAS-OATs/astrocook-data` (folder `benchmark/`), found as the spectra of the screenshots (`--data`, `$ASTROCOOK_DATA/benchmark`, a folder ignored by git); this public repository holds the code and, at most, aggregate numbers Guido chooses to publish. The benchmark on them runs locally, in a cloud session with `astrocook-data` attached, or in a workflow of the private repository that checks out this one (its logs and artifacts stay private; the CI of this repository never reads private data). *Split* (proposed at the planning of 7a, stratified by z_em and S/N; Guido, 13 November 2026): *development*, 8 spectra, used in the phase: Q0329-385, HE1341-1020, PKS0237-23 (z_em < 2.5), PKS0329-255, Q0453-423, HE1347-2457 (2.5-3.2), SDSS J1249-015, PKS2126-158 (> 3.2), S/N 46-140; *test*, the other 20, blind until the version is frozen, among them HE0940–1050 and the extreme cases (Q1209+0919 at S/N 20, Q2237-0608 at z_em 4.56), which the mocks cover in development.
- HE0940–1050 is one of the test spectra and stays blind (section 6). The synthetic twins (7a.3) are made for the spectra whose published analysis gives a model; those of the development part are the realistic cases of the phase.

**Plan** (13 November 2026; estimates in working time of the session, before Guido's review, calibrated on 6-quater: about 2 h for a mechanical PR with tests and documentation, 2-4 h for a review; the prototypes are new ground and their estimates are rougher). The PRs of 7-bis (B1-B11, next section) are interleaved (Guido, 13 November 2026): at most one PR of 7 and one of 7-bis open at a time, each on its own branch; a 7-bis PR touches in `ROADMAP.md` only its own row of the plan and of the register, and is rebased on master before the merge.

*7a – Measure* (first: the prototypes of 7b-7d have nothing to be measured against without it)

| PR | Content | Model | Estimate | Done |
|---|---|---|---|---|
| 7a.1 | Mocks with the cases the grid lacks, in `astrocook/mocks` (independent of `fitting/`; each new field of `MockTruth` off by default, so the mocks made before are unchanged): (a) BAL troughs, a smooth optical depth in velocity on CIV, SiIV, NV and Ly-a in the ratio of their oscillator strengths; (b) a richer set of QSO emission lines (OVI/Ly-b, Ly-a, NV, SiIV+OIV], CIV, CIII] where covered), with asymmetries and a velocity shift per line (the CIV blueshift), z_em the systemic one (`QSO_EMISSION_LINES` unchanged, the new set under another name); (c) tellurics: the molecular transmission (`trans_ma`) of the SkyCalc model of the pre-release `data-telluric-2026` (Paranal, airmass 1.1, PWV 2.5 mm), cut to about 550-1000 nm and resampled into a small file of the package (a few hundred kB), its source and parameters documented; its optical depth split into an O2 part (the gamma, B and A bands, about 627-634, 686-697 and 759-772 nm) and an H2O part (the rest), by wavelength: an approximation, stated, since the cloud sessions cannot reach the ESO hosts for a second SkyCalc run (Guido, 16 November 2026; proper templates per molecule are a choice of 7b.1); the mock scales each part (airmass, PWV) and shifts it by a small velocity, so that the model of 7b.6 is not perfect by construction, and applies it in the observed frame before the convolution; (d) a proximity zone (the density of the forest reduced near z_em) and associated absorbers up to a few thousand km/s beyond z_em (F55); (e) a truth map per pixel: what absorbs each pixel (forest, metal, strong H I, BAL, telluric), from the truth alone, for the metrics of 7a.2. In `validate_mocks.py`, `--grid robustness`: an axis of cases (`dla`, `subdla`, `bal`, `emission`, `telluric`, `proximity`, `combined`) on the base of forest and metals, by S/N 10, 30, 100, R 8900, 45000, z_em 3.5, 4.5 (the tellurics fall in the window) and 3 seeds, 252 mocks; the recipes get a z_em perturbed by a draw of sigma 500 km/s per mock, reproducible from the seed, the true one in the records (Guido, 13 November 2026). The grid `full` and its baseline unchanged (`--compare` identical) | Opus | half a day to a day | done (PR #777, 16 November 2026, branch `claude/7a1-robustness-mocks`): new fields `bal_troughs` (`BALTrough`), `telluric` (`TelluricModel`), `proximity_kms` of `MockTruth`, `shift_kms` and `asymmetry` of `EmissionLine`, `QSO_EMISSION_LINES_FULL`; `lya_forest(z_em=, proximity_kms=)` (columns divided by 1 + omega, omega = (proximity_kms / dv)^2, same draws), `associated_absorbers` (in seed 2 also one at 3500-6000 km/s, beyond the margin of F55: Guido, 16 November 2026); `truth_map` (depth per category: forest, metal, strong_hi, bal, telluric); the telluric template `astrocook/mocks/data/telluric_paranal.npz` (271 kB, tau of `trans_ma` 550-1000 nm at 0.005 nm, made by `benchmarks/make_telluric_template.py` from `sky_telluric.fits`, found in the git history, commit 9740b65c). Mocks made before unchanged (snapshot `validate_mocks_full_mock`) |
| 7a.2 | Metrics of the benchmark for the new goals, with the truth map of 7a.1: the continuum (estimated over true: median, p05/p95, integrated deviation) on the damped wings of the DLAs (the saturated core excluded), across the BAL troughs, on the peaks of the emission lines (within one FWHM) and in the proximity zone; the pull of the continuum against its uncertainty, active when the column `dcont` exists (the name fixed here, for 7b.2); z_em of an estimator against the truth (velocity difference), when the chain has one (F153, 7b.3); the false candidates on tellurics (lines on pixels of true transmission below 0.9 and no true absorber within the tolerance), counted apart; the absorbed pixels of the truth map: identified right, identified wrong, in a region without identification, not detected (the Ly-a by default of F51 counted apart); a confusion table of the flags of the components (#517 and those of the benchmark) against the wrong ones (pull above 3 sigma, spurious, at a bound). The current algorithms on both grids: `benchmarks/results/validate_mocks_robustness_baseline.json`, the reference of the phase; the baseline of `full` unchanged | Opus | 3-4 h | done (PR #783, 16 November 2026, branch `claude/7a2-benchmark-metrics`): `phase7` in the records and summaries of `validate_mocks.py`, on both grids. Definitions (Guido, 16 November 2026): a pixel of the truth map is absorbed above a depth of 0.01 (`ABSORBED_DEPTH`), detectable above max(0.01, 3 sigma per resolution element) (0.19, 0.063, 0.019 at S/N 10, 30, 100), each pixel in the row of the category that absorbs it most; the continuum (ratio, integrated deviation signed and absolute, pull against `dcont`, the name fixed here) in the forest, on the red side, on the damped wings (the Ly-a of the DLA or sub-DLA alone with T 0.1-0.95), across the BAL troughs, within one FWHM of the emission peaks and in the proximity zone; what the truth leaves in the forest (mean T, fraction of pixels with T > 0.97: Guido, 16 November 2026); `z_em_estimate`, None until 7b.3; the false candidates on tellurics (a line on telluric T < 0.9, no true absorber of the series within 20 km/s); the pixel outcomes (undetected, region without identification, Ly-a by default of F51 right or wrong, right, wrong); the flags of #517 and `absurd_errors` against good, pull > 3 sigma, spurious ("at a bound" is a flag, not a class); the false candidates on associated absorbers counted apart (`n_false_associated`). Reference: `benchmarks/results/validate_mocks_robustness_baseline.json` (`full` and `robustness`; `--save-reference`, `--compare`); `full` 42 min, `--compare` identical to the baseline; `robustness` 105 min (the truth map about 50 s per mock). Results of the current algorithms (medians over the mocks): continuum on the DLA wings 0.83 (`full`) and 0.76 (`robustness`), sub-DLA 0.89 and 0.80, across the BAL 0.22 (deviation -71%; continuum <= 0 in 3 of 72 mocks, F160), proximity zone 0.86, emission peaks 0.996 (p05 0.70); forest T > 0.97 on 13% (`full`) and 3% (`robustness`) of the pixels; of the detectable absorbed pixels 73-88% carry only the Ly-a of F51, about 1% an identification tested, 10-22% are in no region (78% of the BAL pixels, 96% of the telluric ones); false candidates on tellurics 60 CIV, 38 SiIV, 41 MgII, 55 Ly-a (nearly all the false metals of the case `telluric`), on associated absorbers 49 CIV and 39 SiIV (nearly all those of `proximity`); the flags mark 7% of the bad components (pull > 3 or spurious), 71% of the flagged are bad, 56% of the unflagged too |
| 7a.3 | The twins and the import of the published analysis: a mock from an arbitrary component list on a given pixel grid, with a given error spectrum, resolution and continuum; the import of a VPFIT model (fort.26, or the tables of the paper) into a `SystemList`; the twins of the red side of the spectra of the sample, from the published CIV and SiIV components (in `astrocook-data`; the other metals and the forest are not in the lists, so a twin covers the CIV and SiIV regions only) | Opus | half a day | waits for the lists of the CIV and SiIV components of D'Odorico et al. (2010), z, logN and b (Guido, 16 November 2026) |
| 7a.4 | The benchmark on the real spectra: the SQUAD DR1 spectra of the sample downloaded into `astrocook-data` (`benchmark/`, with the split); a script (or a mode of `validate_mocks.py`) that reads `$ASTROCOOK_DATA/benchmark` and compares the chain with the published CIV and SiIV components (completeness, purity, N and b, the CIV mass density, on the redshift path searched by the paper) and with the DLAs of the SQUAD catalogue; the effect of the difference between the SQUAD reductions and those of the paper (Guido, 13 November 2026); the results written to `astrocook-data`; the workflow in the private repository | Opus | 3-4 h | waits for the lists of 7a.3 and the 28 SQUAD DR1 spectra in `astrocook-data` (their host, the gDMCP data platform, is not reachable from the cloud sessions: downloaded by Guido) |
| 7a.5 | F154: the evaluation of a model on a fixed set of pixels, without fitting, as a read-only recipe (working name `evaluate_model`, menu Absorbers; so also an MCP tool): the pixels as ranges in x or a mask column; the components of the model all (default) or a selection by system, series or z range; an option to write the normalised residuals in a column. Output: chi2, number of pixels, k (the free parameters of the constraint model: tied and frozen ones do not count), dof, reduced chi2, AIC, AICc, BIC and SpIC (Webb et al. 2021, MNRAS, doi:10.1093/mnras/staa3551, arXiv:2009.08336; implemented from the paper, it needs the footprint in pixels of each component, from the model) (Guido, 13 November 2026), and the mean of the residuals and their fraction beyond 3 sigma. One function for the private `SystemList._compute_fixed_aic` and `VoigtFitter.chi2_on_fit_pixels`, the chi2 of `metals_auto` (F111) unchanged. The basis of the comparison between alternatives (credibility, 7c) and of the evaluation of an imported model (7a.3) | Opus | 3-4 h | done (PR #785, 16 November 2026, branch `claude/7a5-evaluate-model`): `absorbers.evaluate_model` (menu Absorbers > Evaluate Model, MCP tool `absorbers_evaluate_model`; requires the resolution and the continuum) and `astrocook/fitting/model_evaluation.py`. The pixels: ranges in nm (`x_ranges`), a mask column, or both; the model always the whole list (Guido, 16 November 2026), the selection (`uuids`, `series`, `z_min`, `z_max`) only decides who is counted; k the free parameters of the constraint model (frozen and linked do not count). SpIC from the paper (eq. 3, 5, 6): R_a = sum (1 - I^a)/sigma, I^a the profile of the component alone convolved with the resolution (the paper is ambiguous: its eq. 7, the model as the product of the profiles, holds only before the convolution; Guido will ask D. Milakovic, the work goes on meanwhile; another choice touches only `VoigtFitter.component_profiles`, R_a and the influence, not chi2, k, AIC, AICc, BIC), R_a >= k_a + 2 (section 2.4); `spic` with f = 1/2 (SpIC_H, the paper's all-rounder), `spic_a` (f = 1), `spic_b` (f = 0); a velocity component of the paper = the components joined by links of their z (Guido: distinct from the groups of the fits). Counted by default only the components whose influence on the pixels, sum ((1 - I)/sigma)^2 (about the chi2 change of removing them), is at least `min_influence` = 1 (Guido, 16 November 2026: counting all gave k = 420 on 1286 pixels of a 5 nm window of a mock with 140 components, reduced chi2 1.47 for chi2/n 0.99, now k = 24 and 1.009; a rule on the line centres missed the wings); the `uuids` given always counted. One residual function, `VoigtFitter.normalized_residuals`, for `chi2_on_fit_pixels` (`metals_auto`, chi2 identical to the bit, characterization test) and the evaluation; `SystemList._compute_fixed_aic` called by nobody, left as it is. Time: 0.07 s on a 5 nm window, 0.7 s on the whole of a 58 thousand pixel mock with 140 components (the model 0.05 s, the rest the profiles of the components one by one) |

*7b – Emission model*

| PR | Content | Model | Estimate | Done |
|---|---|---|---|---|
| 7b.1 | Design document: the options for the continuum (upper envelope, splines with adaptive knots and uncertainty, Gaussian processes, PCA prediction of the forest from the red side, the mean transmission expected at each redshift, spikes recognised by their width, the F60 bias), the DLAs and BALs (masked or modelled with their wings), the emission lines and z_em (F153), the telluric model (7b.6: a fit with templates per molecule, as Molecfit without its line-by-line radiative transfer, the templates from SkyCalc runs per PWV made by Guido or from the split by wavelength of 7a.1; Molecfit itself, as an external option), the uncertainty of the continuum and its propagation to N and b (a free local continuum and zero level in the fit window, section 2); F53, F56. The common decisions with 7c (`design/identification.md`, "Common decisions with 7b", accepted by Guido on 18 November 2026): the vocabulary of the regions (the enumeration of 7c.3) adopted, or its differences discussed; 7b writes the regions it owns (the masks of the BAL troughs, of the DLA cores and wings) and the column `telluric` (the model transmission, 1 where there is none) | Opus | half a day, then Guido's reading | |
| 7b.2 | Prototype of the continuum, selectable next to `estimate_auto` (the method a parameter of the recipe, its default the old one), with its uncertainty; measured on both grids and on the development spectra | Opus | 1-2 days | |
| 7b.3 | F153: an estimator of z_em from the emission lines (and the onset of the forest), as a recipe; its accuracy on the mocks | Opus | half a day | |
| 7b.4 | The uncertainty of the continuum in the fit: a local continuum (and zero level) free in the fit window, or the continuum uncertainty propagated to N and b, as decided in 7b.1 | Opus | 1 day | |
| 7b.5 | The new continuum the default, if the benchmark shows it better; the release notes | Sonnet | 1 h | |
| 7b.6 | The telluric model, in-house (suggested by V. D'Odorico; Guido, 13 November 2026): a fit of the telluric transmission with the templates of 7a.1 (a depth per part, H2O apart from the rest, a small velocity shift, the resolution of the spectrum), with the local continuum, as a recipe; measured on the telluric case of the robustness grid and on the development spectra | Opus | 1-2 days | |

*7c – Identification*

| PR | Content | Model | Estimate | Done |
|---|---|---|---|---|
| 7c.1 | F68: the Ly-a/Ly-b pairs. The cause (found at the planning of 7c, 16 November 2026, checked on the mock of `test_ly_ab_pairs_are_never_rated`): not the order of the lines of `Ly_ab`, which the kinematic search and `rate_doublet_candidate` read consistently, but the velocity frame of the scoring, linear (`doppler_optical`, c dlambda/lambda), against the lags of the lines, logarithmic (c ln(lambda_2/lambda_1), `get_multiplet_velocity_lags`): about 4100 km/s apart for Ly-a/Ly-b (-46 843 against -50 935), so the two windows never overlap; the same offset misaligns the metal doublets by 0.4 km/s (CIV), 1 (MgII) and 6 (SiIV), perhaps part of the lower completeness of SiIV (F57) and of the case of F119. Content (Guido, 16 November 2026): one logarithmic frame for the scoring of every doublet, metals included; `STANDARD_MULTIPLETS` unchanged; an accepted pair gives a `Ly_a` system as now (a `Ly_ab` system, fitted on both lines, is a choice of 7c.2); the threshold of the pairs 0.5 as for the metals (if the purity measured is bad, reported, not changed); the row of F68 corrected, the characterization test reversed, targeted tests (a clean pair recognised, the Ly-b no longer a Ly-a system, the score of a clean SiIV with small b before and after). To measure: the Ly-a/Ly-b pairs against the truth (the forest of the mocks has the Ly-b): completeness by logN, purity, the false pairs of the dense forest, the Ly-b pixels that were "Ly-a by default, wrong" in 7a.2; the metals: completeness by sigma class (F57) for CIV, SiIV, MgII, purity, the case of F119 again; bias and pulls of `doublets_auto`, the time. The full benchmark on both grids and the slow tests (`doublets_auto` and `components_auto` have `Ly_ab` among their defaults). The reference: `validate_mocks_robustness_baseline.json` (both grids) regenerated, the reference of the identification and the default of `--compare` from then on; `validate_mocks_baseline.json` of 2.2.0 kept as history; the rule written in "Rules and decisions in force" | Opus | 2-3 h, plus about 2.5 h of benchmark |done (PR #793, 16 November 2026, branch `claude/7c1-ly-ab-pairs`): `rate_doublet_candidate` in the logarithmic frame, c ln(lambda / lambda_1), as the lags and the kinematic search; nothing else changed (threshold 0.5, `STANDARD_MULTIPLETS`, a pair gives a `Ly_a` system). Clean mocks (S/N 100, R 45000), scores before and after: Ly-a/Ly-b 0 and 0.99, CIV 0.996 and 0.998, MgII 0.96 and 0.98, SiIV logN 13 with b 5 0.21 and 0.99, with b 3 0 and 0.97, the system of F119 (logN 12.9, b 9) 0.56 and 0.99. Found here: the mocks of `validate_mocks.py` start at 104 nm in the rest frame, above the Ly-b emission, so no Ly-b line falls in them and the grids measure the metals only; the pairs are measured by a new script, `benchmarks/validate_lyab_pairs.py` (the same mocks from 97.5 nm, 54: S/N 10, 30, 100, R 8900, 45000, z_em 2.5, 3.5, 4.5, 3 seeds); `doublets_auto` does not have `Ly_ab` among its defaults (`identify_lines` and `components_auto` do); the pixel outcomes of 7a.2 count a Ly-b pixel with the Ly-a by default as right (an H I row). *Metals*, `full` (against 2.2.0): purity / completeness CIV 0.993 / 0.874 and 0.993 / 0.870, SiIV 0.989 / 0.690 and 0.987 / 0.857, MgII 0.988 / 0.902 and 0.988 / 0.913; by sigma of the weaker line (3-5, 5-10, > 10; identification alone, the old and the new scoring on the same mocks) SiIV 0.08, 0.27, 0.76 and 0.17, 0.45, 0.93, CIV 0.12, 0.58, 0.96 and 0.12, 0.55, 0.96, MgII -, 0.14, 0.95 and -, 0.29, 0.96: the misalignment of 6 km/s was most of the lower completeness of SiIV (F57); fit: matched 1170 and 1242, spurious 11 and 12, flagged 22 and 22, pull widths z / logN / b 4.14 / 5.38 / 4.34 and 3.89 / 5.19 / 4.14, system column 0.077 dex and 0.074 (pull 2.27 and 2.22). `robustness` (against 7a.2): CIV 0.848 / 0.818 and 0.859 / 0.809, SiIV 0.767 / 0.649 and 0.823 / 0.779, MgII 0.894 / 0.802 and 0.864 / 0.808 (its false candidates 48 and 64, nearly all on tellurics, 41 and 57: the narrow telluric lines pass the score more easily once aligned; for 7c.5); false candidates on tellurics 194 and 172 (SiIV 38 and 11, CIV 60 and 49); matched 1753 and 1831, spurious 261 and 252; continuum identical, pixel outcomes in the third decimal. Time: `identify_lines` unchanged, `doublets_auto` about 14% longer (more SiIV systems; scaled by the time of the mocks, the runs on different machines). *Pairs* (`validate_lyab_pairs.py`, old and new scoring): 0 and 674 accepted, purity 0.84 (a pair within 20 km/s of a true component with a detectable Ly-b; the same match for the candidates moved by 300 km/s: 0.38, 0.66 at z_em 4.5, where the pairs are little better than chance); completeness over the components with a detectable Ly-b 0 and 0.087 (0.27 at z_em 2.5, 0.05 at 4.5; by logN 12-13 0.08, 13-14 0.08, 14-15 0.11, > 15 0.10): the regions of the dense forest are merged, and a pair is found per pair of regions, not per component; Ly-b lines with a Ly-a system on them 4063 and 3523 of 10206; the pixels where the Ly-b absorbs most: 91.5% and 74.0% with the Ly-a by default, 0 and 17.5% labelled Ly-b by a pair, 8.5% in no region. The threshold of 0.5 kept (Guido): the purity is not bad, the low completeness is the region approach (for 7c.2). For 7c.2: an accepted pair gives a `Ly_a` system with the score of the pair at every peak of the Ly-a region, and the Ly-b region, labelled `Ly_b`, loses the Ly-a by default also for the true Ly-a lines blended in it. Real spectrum J0103-1305 (snapshots regenerated): 18 SiIV systems where there were none (5 groups; the one at z = 3.51 on known CIV), one CIV (z = 3.2255) from 0.64 to 0.46, dropped; the full chain 14 SiIV components, the rest within the tolerances. The reference `validate_mocks_robustness_baseline.json` (both grids) regenerated, the reference of `--compare` (the rule in "Rules and decisions in force"); `--compare-baseline` for 2.2.0; the completeness by significance per species in the summaries. Tests: `test_ly_ab_pairs_are_rated` (reversed), `test_clean_doublets_score_high`, `test_doublet_score_frame_is_logarithmic`; slow tests pass after the snapshot of the full chain on J0103-1305 |
| 7c.2 | Design document, `design/identification.md`: identification as a comparison between hypotheses. Outline (Guido, 16 November 2026): (1) the current algorithm and its numbers (7a.2, 7c.1); (2) the literature (automatic fitters: AutoVP, VPFIT/RDGEN, AI-VPFIT with the SpIC; doublet searches with matched filters in SDSS, CIV and MgII; classifiers of DLAs and MgII; telluric masks and models), what to reuse; (3) the hypotheses: the regions grouped by the transitions that can link them, for each group alternatives (metal X at z, Ly-a interloper, unidentified, telluric, BAL or DLA from 7b), each a forward model convolved with the resolution and evaluated by `evaluate_model` (so F1 = F2^gamma after the convolution and the multiplets beyond two lines), the credibility as weights of a criterion (exp(-dIC/2)); the evaluation without a fit (N and b from the AOD) for every alternative, a short fit only for the finalists, the two compared; a prior on the incidence of the ions with z (dN/dz) discussed as an option, with its sources and its effect on the experiment; (4) an experiment in the document: on mocks with a known truth, controlled pairs of alternatives (true against false, metal against Ly-a interloper, doublet against chance coincidence) evaluated with chi2, AIC, AICc, BIC and SpIC (H, a, b; the convolved profile of 7a.5, without waiting for D. Milakovic), which picks the true one more often by S/N and R; a script in `benchmarks/`, no production code; (5) the detection decoupled from the continuum (F57: a matched filter on the normalised flux, completeness down to 3 sigma), the scoring (F119); (6) Ly-a only when asked, the assumed Ly-a stated as such, the interlopers (F51, F93 item 3); (7) tellurics as a hypothesis (F54, F155; the model of 7b.6, applied in 7c.5); (8) the standard multiplets (#455); (9) the pairs and alternatives kept (F93, F119): the requirements of 7c.3; (10) the cost: the number of hypotheses per spectrum, the profiles of `evaluate_model` in a window around their lines (a threshold to agree) and cached for the components that do not change between alternatives (0.7 s on a whole mock of 58 thousand pixels and 140 components, the model 0.05 s, 0.07 s on a 5 nm window: PR 7a.5), a time budget consistent with the hour of section 5; (11) the criterion of "better" for 7c.8: which metrics of 7a.2 must improve and which must not worsen; (12) the interface of the prototype, selectable next to `identify_lines`; the system of an accepted Ly-a/Ly-b pair (`Ly_a` or `Ly_ab`, from 7c.1). Common decisions with 7b, taken to the conversation "7 · plan of the phase" (Guido, 16 November 2026): one vocabulary of the regions (BAL, DLA, telluric, forest, metal, unidentified) and who writes it; whether 7c.4 starts from the current continuum (proposed: yes, measured again on that of 7b.2); the telluric model of 7b.6 as a column of the spectrum that the identification reads (proposed: `telluric`) | Opus | half a day to a day, plus a few hours for the experiment, then Guido's reading | done (PR #797, 18 November 2026, branch `claude/7c2-identification-design`). Decided (Guido, 18 November 2026): SpIC_H the criterion, BIC reported; no prior in the prototype, metals of the forest off by default; detection by matched filters at 3 sigma; R^2 kept as the pre-filter; groups of the identification distinct from those of the fits, at most 6 regions, n + 1 hypotheses; pixel window ±(3 FWHM + 5 b); an accepted Ly-a/Ly-b pair an H I system with its Ly-b (`Ly_ab`, `Ly_a` where the Ly-b is not covered; the higher orders of the series not used by the identification: blended, not always accessible, a matter of the fit of the strong systems); the two roles of the multiplets; the three common decisions with 7b (the vocabulary of the regions, written by 7b for what it owns, labelled by the identification; 7c.4 on the current continuum; the column `telluric`). Added after the reading and accepted (Guido, 18 November 2026; PR merged): the sweep of the anchored systems for the rarer and single lines (section 3.5), all the anchors by default, the short fit as starting values (not the model), the display of the graph with F163 before or with 7c.7, a telluric model made outside Astrocook or a flux already corrected (section 7), the replay of the old logs after 7c.8, the calibration of the credibility (F164, PR #798) as a criterion of "better" with the weights tempered for a chi2 that is not calibrated, the forest as a statistical hypothesis as an option after 7b.1. The document: `design/identification.md`, the experiment `benchmarks/identification_hypotheses.py` (no production code). Proposed: detection by matched filters on the normalised flux, groups of regions linked by the transitions, n + 1 hypotheses per group (each candidate as an anchor, plus all by default; unidentified a free line), AOD screen then a short fit of the finalists, the credibility exp(-dIC/2) with SpIC_H (BIC reported), no prior in the prototype (metals of the forest off by default), the Lyman series modelled whole, two roles of the multiplets (search anchors and model lines, #455), `identify_lines(method="hypotheses")`. Experiment (3000 trials, S/N 10, 30, 100 by R 8900, 45000, isolated lines, true continuum; fraction that picks the truth after the fit, above 3 sigma): CIV doublet against two Ly-a, chi2 0.04 above 10 sigma, AIC 0.88-0.94, BIC 0.99-1.00, SpIC_H 0.97-0.99; with the Ly-b of the Ly-a in the model, BIC and SpIC 0.99-1.00; two Ly-a at the CIV lag (chance coincidence), AIC 0.89-0.98, BIC 0.73-0.95, SpIC_H 0.82-0.97 (the misses have columns that differ by about the doublet ratio); the wrong pairing (alias) 0.96-1.00 for every penalised criterion, a doublet at 3-5 sigma against no absorber 1.00, so the completeness of F57 at 3-5 sigma (0.12-0.17) is lost in the detection, not in the comparison; without a fit (AOD) the doublet against two Ly-a only 0.76-0.94 (SpIC_H); prior odds Ly-a : CIV of 10 and 100 keep 0.82 and 0.16 of the true doublets of the forest (SpIC_H, Ly-a alone modelled), 0.97 and 0.88 with the Ly-b. Time: 12-40 ms per trial (two hypotheses, AOD, fit, evaluation) on short mocks, 1 s on a 85 nm mock: `evaluate_model` and the fit work on the whole spectrum, the windows of section 10 needed. |
| 7c.3 | The data structure of the candidates with their alternatives and credibility, and their place in the archive: one mechanism for the structured data that a FITS header cannot hold, used also by the Bayesian samples (F28, from 7-bis). What it holds: the requirements of `design/identification.md`, section 9 (per group the regions and the edges that made it, per hypothesis its components, stage, statistics, weight, status and reason, the dropped candidates, the parameters of the run). Decided at the planning (Guido, 18 November 2026): (a) in the archive, files of their own in the `.acs2` tar, one per datum, with a name and a schema version: JSON for the identification, `.npz` without pickle for the Bayesian samples; the old archives open unchanged, and whether a 2.7 opens an archive with the new files (ignoring them) is checked and stated; (b) each datum with its owner in memory: the identification with the spectrum (as `abs_ids` and its metadata today; it cites the uuids of the systems), the samples with the system list (where `bayesian_results` is today), one mechanism of save and load for both, so that they follow the history of the session and the undo with no extra code; (c) the current method of `identify_lines` fills the structure in a minimal form (the regions, the pairs as edges, one hypothesis per accepted pair, the dropped candidates): F93 item 1 closed, the mechanism tried on real data, the numbers unchanged (the `*_json` maps of the header left as they are; their replacement with 7c.6 or 7c.7); (d) the vocabulary of the regions of 7c.2 ("Common decisions with 7b", item 1: BAL, DLA, telluric, forest assumed or identified, metal, unidentified) as an enumeration shared in `core/`, which 7b uses for the regions it owns. F28 closed (`test_bayesian_results_are_saved` no longer xfail) | Opus | half a day || done (PR #801, 18 November 2026, branch `claude/7c3-candidates-archive`). The structure: `astrocook/core/identification.py`, frozen dataclasses (`Identification` with `method`, `params`, `regions`, `groups`, `dropped`, `x_unit`; `Region` with `kind`, `transitions`, `significance`; `Group` with its `edges` and `hypotheses`; `Hypothesis` with `components`, `stage`, `status`, `reason`, `weight`, `chi2`, `k`, `n`, `criteria`; `Component` with the errors, `candidate` and `uuid`), in the metadata of the spectrum (`identification`), JSON with a schema name and version (a later version refused with a WARNING at the load, an unknown field ignored); `RegionKind` the vocabulary (`bal`, `dla`, `telluric`, `forest`, `forest_assumed`, `metal`, `unidentified`), with `TRUTH_CATEGORY`, its map to the categories of the truth map of 7a.1. The archive: `astrocook/io/archive_data.py`, a table of the structured data (name, owner, metadata key, format, writer, reader), files `<name>_identification.json` and `<name>_bayesian_samples.npz` in the tar (`_orig` for the original data), left out of the FITS headers; a file that cannot be read is left out with a WARNING. Checked: 2.7.1 opens an archive with the new files, ignoring them (a save from 2.7.1 drops them); the archives without them open as before. The current method fills it: every region (also those no identification labels, `unidentified`, with IDs past those of `abs_ids` after a merge), the accepted pairs as edges (transitions, z, lag, R^2), the groups as the regions linked by them with one hypothesis per pair (stage `score`, status `accepted`, the reason in words), a group of one region per region of the forest assumed Ly_a, the dropped candidates; `populate_from_identification` writes the uuids of the components made from the candidates; `components_auto` merges the structures of its steps as it merges `abs_ids`; `forest_auto` replaces it as it replaces `abs_ids`. The `*_json` maps, `region_identifications` and `RecipeResult.extra` unchanged. Size: 5 kB (8 regions) after `doublets_auto` and 194 kB (221 regions) after `lya_auto` on the mock of the full chain. F28 fixed, F93 item 1 fixed (the window still finds the pairs again: 7c.7), F165 found and fixed (a false WARNING at the second save of an identification). Left: the uuids are those given by `populate_from_identification` (a component removed by a later refit keeps its uuid in the structure); the `_orig` members are written but not read, as the other `_orig` members (no reader of the original data exists); the reading of the structure by the GUI (7c.7) and the MCP outputs (7c.6). About 3 hours. |
| 7c.4 | Prototype of the identification, selectable next to `identify_lines`, by `design/identification.md` (sections 3, 5, 6, 10-12 and the decisions of 18 November 2026), in four PRs (Guido, 18 November 2026): each measured alone, the numbers of the pipelines changed only by 7c.4c and only with `method="hypotheses"`. On the current continuum (`estimate_auto`, the common decision with 7b); the tellurics (also the mask of the bands before 7b.6, section 7) in 7c.5, the prior of the incidence not in the prototype; the development spectra when the data of 7a.4 arrive | | | |
| 7c.4a | The evaluation in a window (section 10): `evaluate_model` and the short fit on a slice of the spectrum, the union of ±(3 FWHM + 5 b) around the lines of a group, the components outside the group fixed and their transmission computed once per group (a cache). Infrastructure: no recipe changes its numbers. To measure: the statistics of the window equal to those on the whole spectrum within a tolerance stated (tests); the time per hypothesis (today 1 s on a mock of 85 nm, the aim 10-50 ms) | Opus | 3-4 h | done (PR #804, 18 November 2026, branch `claude/7c4a-window-evaluation`): `astrocook/fitting/evaluation_window.py`, `EvaluationWindow` (`around`: the pixels within ±(3 FWHM + 5 b) of every line of every hypothesis of a group, `window_mask`; `from_mask`: any pixels) cuts the pixels out of the spectrum as a spectrum of their own (flux, error, resolution, the normalization of the parent as its continuum) and keeps the fixed components (the rest of the list, at their values, links resolved) whose lines reach them (`line_reach`, the reach of a line in the model, TAU_CUT, plus 5 sigma of the kernel; factored out of `voigt_optical_depth`, numbers identical), once per group. `evaluate`: `evaluate_model` on the window, the fixed components in the model and never counted. `fit`: `VoigtFitter.fit`, the fixed components frozen, on the pixels the fitter takes from the data on the whole spectrum (`VoigtFitter.feature_mask`, factored out of `prepare_fit_context` without changes: out to the continuum or a saddle between blends) plus the window, each run widened by 60 pixels (the fitter's model slice is ±50, its smoothing about 8), cut out of the parent per fit; it returns the hypothesis alone with its constraints. No recipe changed, no production caller yet (7c.4c). *Two changes from the row.* (1) The fixed components are carried in the model of the window as frozen components (their optical depth computed once per fit, the background of the fitter), not folded in as their transmission convolved apart and multiplied: on a CIV doublet in a forest of 400 Ly-a lines, the product moved the ΔSpIC_H between the two hypotheses (CIV, two Ly-a) by 0.2 at R 45000 and 1.3 at R 8900 (exp(-tau_1 - tau_2) convolved is not the product of the two convolved); the frozen components give the numbers of the whole spectrum at the same cost. (2) The short fit not limited to the window (Guido's concern, 18 November 2026): the window follows the b of the start, and a fit limited to it was wrong for broad saturated lines started narrow (seed b 10, R 45000: logN 16, b 60 → 15.1, b 80; logN 17, b 40 → 19.0, b 18; logN 15, b 120 → b 82; with seed b 25, right) and moved in crowded regions (25 lines within ±1500 km/s, the fixed list imperfect: 20% missing, logN ±0.3, b ×0.7-1.3; 10-23% of the fits by more than 0.01 dex or 1 km/s, up to 2.3 dex and 80 km/s); with the pixels of the data, identical to the whole spectrum in all these cases (tests). The pixels of the comparison (`around`), the same for all the hypotheses of a group: the union of the windows of ±(3 FWHM + 5 b) of the start and of the pixels that the fit of each hypothesis takes from the data (Guido, 18 November 2026: the most robust; with the window alone a broad line started narrow was compared on its core). They stop where the smoothed flux is back above 0.95 of the continuum (the threshold of `_find_feature_limits`, the same as the fit): the outer wings shallower than about 5-10% are left out (8-14 pixels of a line with logN 16, b 60 at R 45000), accepted (Guido, 18 November 2026: the wings at 0.95 as the fit). *Equality* (tests, `tests/test_evaluation_window.py`, R 45000 and 8900): chi2, AIC, AICc, BIC, SpIC within 1e-5 relative of `evaluate_model` on the whole spectrum with the whole list (measured about 1e-7), the short fit within 1e-8 in z, 1e-4 in logN, 1e-3 km/s in b (identical in the forest and broad-line cases). In the experiment of 7c.2 (`identification_hypotheses.py --window`, 1200 trials, 40 per cell, same seeds): the evaluation without a fit identical (chi2 within 6e-7 on all 2400 hypotheses); after the short fit 37 of 2400 hypotheses differ by more than 1e-4 of the chi2, at most 0.45% (before the change (2): 11%, up to 56%). The fits checked are identical; the difference is in the model of a narrow saturated line (e.g. a false CIV fitted to logN 16.8, b 6.6 at R 8900): the fine grid of the model has the same step but starts elsewhere on a cut spectrum, and such a line changes by up to 0.055 sigma per pixel (5e-4 of the continuum at S/N 100) with the phase of the samples, beyond the 1e-5 that the comment of `B_SAMPLES` (`lsf.py`) states for saturated lines: the whole spectrum is not more right, its grid starts elsewhere (registered as F166, Guido, 18 November 2026). *Time*: per trial (two hypotheses, AOD evaluation, short fit, evaluation), `doublet_vs_lyab` (85 nm mocks, about 30 thousand pixels) median 1.34 s → 0.068 s, p95 6.3 → 0.15 s; the short mocks of the other scenarios unchanged (0.02-0.06 s); the whole run 140 → 21 s on 4 cores. In a forest of 400 lines on 32 thousand pixels a hypothesis (window, short fit, evaluation) 17-22 ms against 130-260 ms on the whole spectrum: within the aim of 10-50 ms. Tests run: `test_evaluation_window.py` and the 46 test files of the fitter, the evaluation and the pipelines that fit (`test_evaluate_model`, `test_voigt_fitter_*`, `test_snapshots_voigt`, `test_fitting`, `test_metals_auto`, `test_recipes_absorbers`, ...), `ruff`, `mypy`; no slow tests and no benchmark: no pipeline calls the new module, and `voigt_optical_depth` and `prepare_fit_context` are unchanged (refactorings covered by the snapshots) |
| 7c.4b | The detection by matched filters (section 5, F57): three widths (the instrumental FWHM combined with b = 5, 15, 40 km/s), a line where S > 3 for any width, its region where S > 1, the joint significance of a doublet at its lag (the filter of the weaker scaled by the ratio of f lambda); a function of its own, not yet used by the pipelines (`abs_mask` stays the continuum's mask). To measure, on both grids and against the truth map of 7a.2: the absorbed pixels detected by sigma class, the false regions, the comparison with `abs_mask` of today | Opus | half a day | done, PR #805 (`claude/7c4b-matched-filter-detection`): `core/line_detection.py` (`matched_filter_detection`, `detect_lines`), not used by the pipelines; `benchmarks/validate_detection.py` (the measurement). Two additions to the row: a joint seed of a doublet also needs each line above 1 sigma, the lines consistent with their ratio (chi2 <= 9) and J at its peak (without them a strong line made a doublet with the noise at every lag: 100% of the forest, 26% of the red side covered on the tiny grid); the threshold of J apart (`multiplet_sigma`), and `z_max` (no doublet beyond z_em, F55). Numbers (full grid, 135 mocks; robustness alike): metal systems with the weaker line in a region, by its S, 2-3 sigma: abs_mask 0.08, lines alone 0.69, with the doublets at 3 sigma 1.00, at 4 sigma 1.00; 1-2 sigma: 0, 0.28, 0.72, 0.61; 3-5 sigma: 0.38, 0.97, 1, 1. Red-side metal pixels at true S 3-5 in a region: 0.03 → 0.98. Price: false regions on the red side (true S <= 1) per 1e4 km/s: abs_mask 0.11, lines alone 4.1, doublets at 3 sigma 13.5, at 4 sigma 4.9 (robustness 0.05, 2.4, 8.2, 3.0); on pure noise the lines alone give 0.8-3 per 1e4 km/s (R 8900-45000), the rest is F168 (the bias of the continuum). Forest: absorbed pixels at true S 5-10 in a region 0.02 → 0.41. Median time 0.1 s per spectrum. The threshold of J: 4 sigma by default (`MULTIPLET_SIGMA`; Guido, 19 November 2026); at 4 a doublet of thin lines needs its stronger line above about 3.6 sigma, already a line alone, so the joint search gives its region to the weaker line rather than finding doublets with both lines below 3 sigma. Tests: `tests/test_line_detection.py` (17), ruff, mypy; no slow tests nor `validate_mocks.py --compare` (no pipeline calls the module) |
| 7c.4c | The prototype: `identify_lines(method="hypotheses")` (default `"doublets"`), with the groups (the R^2 pre-filter, at most 6 regions, split at the weakest edges), n + 1 hypotheses per group, the AOD screen, the short fit of the finalists, the weights with SpIC_H (BIC reported); the Ly-a only when asked, assumed apart from identified (F51, F93 items 2 and 3), an accepted Ly-a/Ly-b pair a `Ly_ab` system (`Ly_a` where the Ly-b is not covered); the species by default all those with an anchor, the purity measured against CIV, SiIV, MgII alone (if it drops, the default the common anchors); it fills the structure of 7c.3; the method reaches `components_auto` (which calls `identify_lines`), and `doublets_auto` takes the systems of the other species. The two new heuristics measured, then decided by Guido (18 November 2026): the threshold of the pre-filter (R^2 0, 0.2, 0.3) and of the screen (ΔIC 10, 20, 40), by completeness, purity and time on a subset of the grid; their purpose written next to the code. The tempering of the weights (ΔIC / max(1, chi2_nu) of the accepted model, section 3.3) an option, off, measured on the mocks; on or not when 7a.3 and 7a.4 give twins and development spectra. The calibration of the credibility (F164) in `validate_mocks.py`: a reliability diagram over all the hypotheses, its expected calibration error and Brier score, by category, S/N and R. To measure: the criteria of section 11 on both grids and with `validate_lyab_pairs.py`, against the current method; the calibration; the time per spectrum | Opus | 1-1.5 days || done, PR #PRNUM (`claude/7c4c-hypotheses-identification`, 20 November 2026): `astrocook/fitting/hypotheses.py`, `identify_lines(method="hypotheses")`, also in `components_auto` and `doublets_auto`; the structure of 7c.3 filled (`Hypothesis.anchors`, `DroppedCandidate.reason`); the calibration in `validate_mocks.py` (`phase7["calibration"]`), the options `--method --prefilter --screen --criterion --temper --max-added`; `benchmarks/identification_thresholds.py` (the sweep). The thresholds (Guido, 19 November 2026): pre-filter R^2 >= 0.3 (0 halves the purity and doubles the time, F169); no screen (every threshold dropped half the true CIV and SiIV, ECE 0.31; an option, off). Added to the row: the hypothesis "no line"; a candidate needs peaks of its lines at its lag; the refinement by the residuals (up to 3 components per finalist, F170); the combination of anchored finalists with disjoint peaks; a cap of 60 free parameters for the short fit (F171, Guido). The species by default the common ones (CIV, SiIV, MgII, Ly_ab): "all" leaves the purity unchanged, takes 0.01-0.03 of the completeness and doubles the time on 18 mocks. Numbers, against `doublets` (completeness / purity): grid full CIV 0.804 / 0.975 (0.870 / 0.993), SiIV 0.795 / 0.949 (0.857 / 0.987), MgII 0.845 / 0.961 (0.913 / 0.988), by significance 5-10 sigma 0.55 (0.48), above 10 0.87 (0.95); the fit after it 1494 true components (1242), 47 spurious (12), pull widths z, logN, b 1.6, 2.1, 1.7 (3.9, 5.2, 4.1), system column pull 1.25 (2.22); detectable pixels outside every region 0.009 (0.100); credibility ECE 0.022, Brier 0.038 (R 8900 0.056, 45000 and 140000 0.011); time 15 s per spectrum to identify (0.17). Grid robustness: CIV 0.678 / 0.872 (0.809 / 0.859), SiIV 0.727 / 0.741 (0.779 / 0.823), MgII 0.732 / 0.920 (0.808 / 0.864); fit 2068 true (1831), 242 spurious (252), pulls 1.7, 1.9, 1.7 (4.5, 4.1, 3.4); tellurics taken for lines 111 (172); ECE 0.060 (R 8900 0.096). Without the refinement: F170. Ly-a/Ly-b pairs (`validate_lyab_pairs.py`, 54 mocks): 117 pairs (674), purity 0.85 (0.84), completeness of the detectable Ly-b 0.012 (0.087), none at z_em 4.5, 40 s per spectrum: not yet working (F169, F171, to diagnose). Not better than `doublets` on the identification, much better on the fit after it: whether to join the identification and the fit is F170, to decide in the planning of 7c before 7c.4d. Merged as a prototype (Guido, 20 November 2026); the numbers of the default method unchanged |
| 7c.4d | The sweep of the anchored systems (section 3.5): for each accepted system, every species of the table with lines in the coverage as a hypothesis "Y at z" (z within the width of the anchor's profile, b free, N from the AOD), accepted when it lowers the criterion of the groups involved; the H I of the metal systems (Ly-a and Ly-b); an accepted species a new edge of the graph. To measure: the regions labelled assumed Ly-a or unidentified that pass to a species, the purity of the rarer species on the mocks, the time | Opus | half a day | |
| 7c.5 | Tellurics in the identification and in the fit: the model of 7b.6 (F155) used to flag or mask the telluric lines (F54) | Opus | 2-3 h | |
| 7c.6 | The output of the pipelines: the candidates, their alternatives and the dropped ones with their score (F119, decided on PR 7 of 6-ter) in `RecipeResult.extra` and in the MCP outputs; F112 (`IdentificationSummary`) rebuilt on them. From 7c.3 (PR #801): the outputs read the structure of `core/identification.py` (`identification` in the metadata of the spectrum), not the `*_json` maps of the header nor `region_identifications`, which the recipes (`recipes/absorbers.py`) and `mcp_server/outputs.py` read today; the uuids of the structure are those given by `populate_from_identification`, so a component removed by a later refit keeps its uuid there: the outputs check them against the system list (a stale uuid reported as such, not as a system) | Sonnet | 2 h | |
| 7c.7 | The GUI of the identification: alternatives and credibility in the Identifications window and the System Inspector; screenshots, agreed with Guido. From 7c.3 (PR #801): the window reads the pairs from the edges of the structure of `core/identification.py` instead of finding them again from the overlap of the redshift ranges (the rest of F93, item 1), and the GUI (`identification_viewer_dialog.py`, `plot_mouse.py`, `main_panels.py`, `main_sessions.py`) reads the structure instead of the `*_json` maps and `region_identifications`; the uuids of the structure can be stale after a refit (as in 7c.6); the maps are then removed from the writer (`spectrum_lines.py`) once nothing reads them, while the archives written before 7c.3, which have the maps and no structure, still open with their identifications (the structure built from the maps at the load, or the maps read as a fallback: to decide at the start of 7c.7) | Opus | half a day | |
| 7c.8 | The new identification the default, if the benchmark shows it better; the release notes | Sonnet | 1 h | |

*7d – Components and fit*

| PR | Content | Model | Estimate | Done |
|---|---|---|---|---|
| 7d.1 | Design document (short): flags and uncertainties of the components (F52, F116: unresolved, saturated, absurd errors), the addition of components from the residuals (F58, F118), the acceptance of a component or an ion (F121), the components at `max_logN` = 19 in the forest fits, the fits at low resolution (F59) | Opus | 3-4 h | |
| 7d.2 | F52, F116: the flags for unresolved and saturated components and for errors without meaning, in the results, the MCP outputs and the GUI | Opus | half a day | |
| 7d.3 | `optimize_system`: the window on every line of the multiplet (F148), a stop criterion that sees a misfit spread over many pixels (F118), the fits that do not converge (F149) | Opus | half a day to a day | |
| 7d.4 | F58: `doublets_auto` adds components from the residuals (with 7d.3) | Opus | 1 day | |
| 7d.5 | F121 and the forest fits: the acceptance of a trial ion of `metals_auto` by its AIC, and the components that reach `max_logN` | Opus | half a day | |
| 7d.6 | F59: the system column at R = 8900 (the cause measured: unresolved saturation or the convolution) | Opus | half a day | |
| 7d.7 | F156: upper limits on the column of a species not detected, as a recipe and in the methods text | Opus | half a day | |
| 7d.8 | A tabulated instrumental profile (LSF) next to the Gaussian one (section 2; UVES and ESPRESSO) | Opus | 1 day | |

*Closing*

| PR | Content | Model | Estimate | Done |
|---|---|---|---|---|
| 7.Z | Closing of the phase: "Result of the phase", the progress table, the rules in force, `CLAUDE.md`, the detail to `ROADMAP_HISTORY.md`; the full benchmark on both grids and the slow tests on the final master; the release | Sonnet | 2 h | |

*Order.* 7a.1, 7a.2, 7a.5, then 7c.1; the design documents 7b.1 and 7c.2 next (Guido reads one while the other is written); then the prototypes of 7b and 7c (7b.6 after 7b.2 and before 7c.5), with 7c.3 before 7c.4; 7d after 7c.4 (its flags and the addition of components are measured on the new identification); 7a.3 and 7a.4 when their material arrives. In all about 15-20 working days of the sessions, plus Guido's reading of three design documents.

*Findings* (register, 3.2): F51, F52, F53, F54, F56, F57, F58, F59, F68, F93, F116, F118, F119 (the scoring), F121, F148, F149 (the criterion), F153-F156, F164 (the calibration of the credibility, with 7c.2); F60 (the 0.055 sigma bias, closed in phase 4, its rest here); #455.

**Orchestration of the conversations** (agreed with Guido, 13 November 2026; the criterion that spends the fewest tokens). The state of the phase lives in this section, not in a conversation. Each PR has a conversation of its own, opened by Guido with the model of the plan and a prompt of a few lines that points to its row here; it reports through its PR description and the update of its row. A planning conversation per block (7a, 7b, ...), not one for the whole phase: it starts from this section, writes the details of the block, opens the sessions of its PRs (title `<PR> · <subject>`, the model of the plan, the tag of the block, a prompt that points to the row), reviews only the PRs of judgement (description, `git diff --stat`, the central files) and closes when the block closes. No conversation follows the sessions it opened or reads the transcripts of another (`CLAUDE.md`, "Choice of the model").

### Phase 7-bis – Small fixes
Agreed with Guido in the review of the findings of PR 13 of 6-quater and with the scope of phase 7 (13 November 2026): what changes how the results are exposed, stored, packaged or shown, and the bugs outside the algorithms; no algorithm of the continuum, the identification or the fit. Interleaved with phase 7 (Guido, 13 November 2026; the rules in "Plan" of phase 7). Some PRs wait for a choice of Guido, written in the row.

| PR | Content | Model | Estimate | Done |
|---|---|---|---|---|
| B1 | F101: `AppController.run_recipe` corrected (the console is kept); the comment of `_NE_GLOBAL_DICT` (`core/spectrum_expressions.py`), which says that `median`, `mean`, `std` work in NumExpr (F128: they do not; no change of behaviour) | Sonnet | 1 h | done (13 November 2026): `run_recipe` runs the recipe through the registry (`recipes.run`), synchronously, on the session given, and returns the new session or the value (Guido's choice); `tests/test_gui_scripting.py` |
| B2 | F152: the `model` column recomputed from the new parameters after `update_component(s)` and `update_constraint`, without a fit (proposed), or the warning "please refit" turned on again; choice: Guido | Sonnet | 2 h | done (16 November 2026; choice of Guido: the model redrawn): `update_component`, `update_components` and `update_constraint` (with a link) recompute the `model` column from the new parameters, without a fit (`_redrawn` in `recipes/absorbers.py`, as `delete_component`); nothing is done without a `model` column (never fitted) or without a continuum; freezing or freeing leaves the model as it is. The snapshots of the manual fits with `update_*` or a link changed in the `model` only (components unchanged); the fits do not change. The instructions of the MCP server and `docs/tutorials/agents.md` say it. The GUI dialog stays off (it would say "refit" to see the model: now not needed). Tests: `tests/test_update_redraws_model.py`; the `resol` of the parsing cases is now 20000 instead of 2 (R=2 exhausted the memory of the convolution, now that the model is computed: not a case of use, see the note of F152) |
| B3 | F151: the tooltips of the entries from `details` or from a `tooltip` key of the registry; choice: Guido | Sonnet | 2 h | done (16 November 2026): from `details` (Guido's choice); `RECIPE_TOOLTIPS` removed; the 10 `details` over 200 characters shortened, a test keeps them ≤ 200 |
| B4 | F114: the `[bayes]` extra in the Docker image (its size measured) or the description of the tool saying it is not there | Sonnet | 1-2 h | done: the extra is in the image (+~1.7 MB, measured with pip: no Docker daemon here; Guido chose it, 16 November 2026) |
| B5 | F30: the edges of `running_std` (and of `running_mean`, unused in the package since #465: corrected or removed); changes the numbers of `calculate_running_std` at the edges, with Guido's agreement first | Sonnet | 1-2 h |  done (16 November 2026): truncated window at the edges for both functions (Guido's choice); `running_mean` corrected, not removed; 2 snapshots regenerated (h=1, 10: only the first/last h pixels differ) |
| B6 | F84: the reference redshift of a conversion of the x axis to velocity kept with the spectrum, and used by the conversions back | Opus | 3-4 h | done (16 November 2026; the design agreed with Guido; PR #779): `SpectrumData.z_ref` (None for an axis in wavelength), set by `x_convert`, read by `Spectrum.axis_z_ref()` (z_em when not known: an axis loaded in velocity, an archive saved before), kept in the `.acs2` (FITS key `Z_REF`, written only for an axis in velocity: the old archives read as before) and in the exported CSV. `x_convert` back to wavelength at the z_ref kept (another one refused, `InvalidParameter`); velocity to velocity re-centred at a z_ref given (before: the values unchanged); `'_current_'` is z_em from a wavelength, the z_ref kept from a velocity. The ticks of the systems, `rebin` and `coadd` use it; two axes in velocity with different z_ref are not combined (`stitch`, `coadd`, `equalize`). The MCP summary of the spectrum says it (`z_ref`). No snapshot changed; nothing in the pipelines converts to velocity. Tests: `tests/test_velocity_axis_z_ref.py`, one in `tests/test_main_plot_systems.py`. Found: F159 |
| B7 | F72: the boolean masks follow a `rebin` or a resampling (the rule, e.g. a new pixel masked if any old pixel in it was: Guido) | Sonnet | 2 h | done (16 November 2026; rule chosen by Guido: *any*, a new pixel is masked if any masked old pixel overlaps it; beyond the data and in gaps False; *Co-add* carries no column but now says so): `follow_mask` in `core/spectrum_operations.py`, used by `rebin` and `resample_on_grid`; no snapshot changed (the masks were dropped before) |
| B8 | F102: the refit of the cleanup of `lya_auto` and `forest_auto` limited to their region; the slow tests (it changes the numbers of these pipelines, not of the benchmark) | Opus | 2-3 h | done (16 November 2026; region chosen by Guido; PR #786): `lya_auto` was already limited (6-bis, PR 4c); `forest_auto` works in `all_ly_forest` (Lyman limit to Ly-alpha emission): the identification on a copy of `abs_mask` restricted to it (the mask of the session given back as it was), both refits with `region_limit`; a CIV of the red side is no longer refitted (`test_forest_auto_does_not_refit_the_red_side`). No snapshot, slow test or benchmark pipeline runs `forest_auto` (off the menus: scripts and MCP); F162 found. Guido's iterative scheme for the whole forest (Ly-alpha forest first, then the higher orders, then Ly-alpha in the Ly-beta/gamma/delta forests, to convergence) is for 7c/7d, not here |
| B9 | The outputs: F112 (the immediate contradiction of `IdentificationSummary`, before 7c.6 rebuilds it), F149 (the reason of a rejected candidate in the report of `optimize_system`, not "AIC changes by inf"), F157 (the chi2 of the groups in the result of every fit) | Sonnet | 2 h | done, 16 November 2026 (PR of the branch `claude/b9-fit-outputs`): F112 (the systems of the summary come from `extra['candidates']`), F149 (the reason and `converged: False` in `tried`, `start_fit` in the report), F157 (`chi2r` in `components_changed`, `fit` also when a fit changes components); only outputs, the numbers of the fits unchanged |
| B10 | F115, the mechanical part: the declared gaps out of the methods text, the name of the spectrum, `session_methods` writing the text to a file; the prose of the paragraphs after phase 7 (Guido's text) | Sonnet | 2 h done, 16 November 2026 (PR of the branch `claude/b10-methods-text`): the declared gaps are out of the text and listed in `MethodsText.gaps` and the `gaps` of the result of `session_methods`; a generic name of the session (`spectrum`, `session`, ...) no longer names the spectrum; `session_methods(path, overwrite)` writes the text to a file (`write_methods`), the file and folder rules in the docs; the prose of the paragraphs is untouched (Guido's text, after phase 7). Tests: methods text, MCP sessions/server/gaps, registry, menus, session log; ruff, mypy. Slow tests and benchmark not needed: only the text produced changes |
| B11 | F158: the residuals under `spectrum_plot`, and a plot of several transitions in velocity around a redshift, as MCP tools (the plotting of `astrocook/plotting.py`) | Sonnet | 2-3 h | done (16 November 2026): `spectrum_plot(residuals=)` and the new tool `spectrum_plot_velocity` (names, defaults and results chosen with Guido); the helpers shared with the Inspector in `plotting.py`; the Inspector untouched (F163); no number changes |

**Result of the phase** (13-16 November 2026; one conversation per PR, opened by the planning conversation, each merged before the next opened):
- *What was done.* B1 (#772): the script console's `app.run_recipe` runs the recipe on the session given through the registry (F101); the comment of `_NE_GLOBAL_DICT` corrected (F128). B2 (#774): the `model` column redrawn after `update_component(s)` and `update_constraint` (F152). B3 (#775): the tooltips of the menus from `details` (F151). B4 (#776): the `[bayes]` extra in the Docker image (F114). B5 (#778): the truncated window at the edges of `running_mean` and `running_std` (F30). B6 (#779): the reference redshift kept with the spectrum (F84). B7 (#780): the boolean masks follow `rebin` and the resampling, rule *any* (F72). B8 (#786): `forest_auto` limited to `all_ly_forest` (F102). B9 (#788): the outputs of the identification and of the fits (F112, F149, F157). B10 (#790): the methods text without the declared gaps, `session_methods` to a file (F115, mechanical part). B11 (#791): the residuals under `spectrum_plot` and the new tool `spectrum_plot_velocity` (F158).
- *Numbers that changed, each with Guido's agreement.* B5 (the edges of `calculate_running_std`; two snapshots regenerated), B7 (the masks are kept, not dropped), B8 (`forest_auto`; the slow tests 5/5). The other PRs change what is exposed or stored, not the numbers (each PR says in its description which tests were run and why the numbers cannot change).
- *Found along the way.* F159 (the recipes that fail on an axis in velocity), F161 (the integer columns interpolated by `rebin`), F162 (`z_em_margin_kms` in the pipelines of the forest), F163 (the two pieces of code of the velocity plot): open, to be reviewed with Guido.
- *Left.* F28 (with 7c.3), F15 (the X-Shooter VIS/NIR data), the prose of F115 (after phase 7, Guido's text).

F15 (`equalize_and_stitch` with narrow absorption) waits for the X-Shooter VIS/NIR data; if they have not arrived by the closing of phase 7, it is closed as not verifiable, the bias documented. F28 (the Bayesian samples in the archive) goes with 7c.3 (one mechanism for the structured data of the archive). Not in the phase (not scheduled yet): the technical issues of 3.1 (#424, #471, #478, #479, #446, #437, item 2 of #573) and the test of `gui/corner_plot.py`.

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
- **2.7.1 (17 November 2026), bug-fix release:** 2.7.0 with the fixes of [#781](https://github.com/DAS-OATs/astrocook/issues/781) (saving after the automatic identification, PR #782) and [#794](https://github.com/DAS-OATs/astrocook/issues/794) (Return in a cell editor of the System Inspector, PR #795), reported by a user. Made from the branch `release/2.7.x` (the tag `v2.7.0` and the two commits cherry-picked), not from master, so that it holds nothing of phase 7 (7c.1 changes the identification); on master only the `CHANGELOG.md` entry and `switcher.json`. Agreed with Guido for this release (17 November 2026); the tag `v2.7.1` goes on `release/2.7.x`.

### 3.1 User-reported issues
GitHub issues are the single source of truth; this section only says *when*
each open issue is addressed. Rule of thumb: bugs before refactoring,
features after, except urgent features in stable areas.

| When | Issues | Notes |
|---|---|---|
| Hotfix (17 November 2026, urgent) | [#794](https://github.com/DAS-OATs/astrocook/issues/794) | `Return` in the editor of a cell of the System Inspector also ran *Refit Selected* (`refit_all`), so the value typed was lost and a frozen parameter came back to its old value: the keys of the table, added in phase 5-bis (2.3.0), had the context `WidgetWithChildrenShortcut`, which includes the editor; now `WidgetShortcut` (the table alone). The fit engine was right: a frozen logN changed by hand stays at its value after `refit_all`. Test `test_inspector_keys_while_editing_a_cell`; done (PR #795), confirmed by the reporter, released in 2.7.1 with the fix of #781 |
| With the natural-language log work | [#424](https://github.com/DAS-OATs/astrocook/issues/424) | logging of manual continuum edits |
| Phase 6-bis (item 3); after the full MCP server (item 2, Guido, November 2026) | #573, items 2 and 3 (closed with item 1 in phase 4) | a user configuration file for the recipe defaults and fit settings (item 2); the session log records the values actually used, not "default" (item 3, PR 2 of phase 6-bis); to be opened as issues when scheduled |
| Phase 6 | [#566](https://github.com/DAS-OATs/astrocook/issues/566) | centralize the log messages (started in phase 3 with `recipes/messages.py`), with the move to per-module loggers; done in phase 6: the loggers per module in PR 3 (#693; `tests/test_logging.py` checks that no module logs on the root logger), the shared messages in `recipes/messages.py` with the PRs per package (#694-#700); the issue to be closed (Guido) |
| After phase 5 | [#471](https://github.com/DAS-OATs/astrocook/issues/471) | combining order-separated spectra; touches `SpectrumV2` (coadd/stitch) |
| After phase 5 | [#478](https://github.com/DAS-OATs/astrocook/issues/478), [#479](https://github.com/DAS-OATs/astrocook/issues/479) | GUI features; wait for the GUI split |
| After phase 5 | [#446](https://github.com/DAS-OATs/astrocook/issues/446), [#437](https://github.com/DAS-OATs/astrocook/issues/437) | flux-calibration filters; independent, not urgent. With them (Guido, October 2026, after F70): a more sophisticated algorithm for the wavelength dependence of the calibration curve (today linear between the effective wavelengths of the filters and constant beyond, solved to reach every magnitude); to be opened as an issue when scheduled |
| When needed (not scheduled; Guido, October 2026) | — | transitions of the V1 table not in the single table of #568 (H2 with non-unique names, CO, Paschen, Brackett, SiI, Balmer beyond H_e, stellar and blended emission markers): added when an analysis needs them, each with its source (the V1 file is in `tests/reference/v1/data/atom_par.dat`) |
| To be decided (Guido) | [#455](https://github.com/DAS-OATs/astrocook/issues/455) (pending question) | standard multiplets (`STANDARD_MULTIPLETS`) for the 10 metal species with 3+ transitions in the table and no entry; unchanged by #568, which touched only the table of transitions; the choice of the multiplets used by the identification belongs to phase 7 (design document of PR 7c.2) |

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
| F28 | `.acs2` Bayesian samples | The posterior samples stored by `fit_bayesian` in the system list metadata (`bayesian_results`, used by the corner plot) are dropped on save without a warning: the FITS header cannot hold them. The fitted values and errors are kept | data loss (recomputable) | `test_bayesian_results_are_saved` | open; since phase 6 (PR 7) a WARNING says at the save that the samples are not saved (astropy skipped them with a warning of its own); phase 7-bis (Guido, 13 November 2026): the samples in the archive; with 7c.3 (Guido, 13 November 2026): one mechanism for the structured data of the archive **Fixed in phase 7, PR 7c.3** (18 November 2026): the samples are a file of their own in the `.acs2` tar, `<name>_bayesian_samples.npz` (no pickle, a schema name and version; the components sampled together share one entry, saved once), by the mechanism of the structured data of the archive (`astrocook/io/archive_data.py`), which also saves the identification; the WARNING at the save removed |
| F30 | `running_mean`, `running_std` | Both use `np.convolve(mode='same')`, which pads with zeros: the first and last h values are biased towards 0 (4/7 of the true mean at the edge for h = 3). V1 `running_mean` repeats the first and last full-window value and ignores NaNs. Used by `find_absorbed_regions` on each interval, with half-windows of hundreds of pixels (5000 km/s / 8 on J0103-1305). Replacing it with the V1 function raises the agreement of the forest mask with V1 from 92.7% to 95.7% and the 5th percentile of the continuum ratio from -4.7% to -3.0% (median unchanged, -1.2%): a partial cause of F17 | numerical, V1/V2 difference | `test_running_mean_and_std_pad_the_edges_with_zeros` (characterization) | fixed in phase 4 for the continuum (#465): `find_absorbed_regions` no longer uses `running_mean`, but a masked mean over a window of the smoothing length in km/s on the pixel grid, exact at the edges (`_masked_window_mean`). The edge bias was worth 0.81 → 0.93 within 2000 km/s of the edges of the forest interval on a mock at S/N 30. The plain normalized convolution of the earlier branch (window in pixels of the clipped array) failed at S/N ≥ 300 (forest continuum up to 1.66 times the true one: the window grows as pixels are clipped, over wavelengths where the continuum varies by more than its noise); the window in km/s does not. `running_mean` and `running_std` themselves padded with zeros (used by `calculate_running_std`). Fixed in phase 7-bis, PR B5 (16 November 2026, Guido's choice): within h pixels of the edges both use the truncated window (the pixels that exist); interior values unchanged; `running_mean` kept and corrected (unused in the package). Only `calculate_running_std` (recipe, menu, `running_std` column) changes, at the edges; no pipeline uses it |
| F51 | `identify_lines` (forest) | Every absorbed region of the forest without a metal candidate is assigned to Ly-a (score 0, no test), even when Ly-a is not among the requested `multiplets`: `multiplets='CIV'` on the full-chain mock gives 233 Ly-a and 6 CIV candidates. In the MCP pilot the Ly-a candidates crowded out the metals in the output (833 on J0103-1305). Found in the MCP pilot (report) | behaviour | — | open: phase 7 (Guido, 13 November 2026); PR 7c.4c: with `method="hypotheses"` the Ly-a of the forest only when `Ly_a` or `Ly_ab` is asked, the assumed apart from the identified (kind `FOREST_ASSUMED`, score 0); the default of `doublets` unchanged, still open until 7c.8 |
| F52 | `VoigtFitterV2` (uncertainties) | Components that end at the bounds of a fit get uncertainties far beyond any meaning, stored as results without a warning: in the `doublets_auto` fit of J0103-1305 (MCP pilot), four CIV components at z = 3.927 with logN = 10.0 and b = 1.0 km/s have dz up to 131, dlogN up to 3800 and db up to 9e4 km/s; 21 of the 101 new components have dlogN > 1 or db > b. Related to F20, F21. Found in the MCP pilot (report) | numerical, error handling | — | open: with [#517](https://github.com/DAS-OATs/astrocook/issues/517); since its PR 1 the components at a bound carry a flag (`ComponentDataV2.flags`), so their errors can be told apart; phase 7 (Guido, 13 November 2026), with F116 |
| F53 | continuum (forest at high z) | On an ESPRESSO spectrum of J0959+1312 (z_em ~ 4.06), `estimate_auto` flags 100% of the pixels between 525 and 600 nm as absorbed, so the continuum there is a straight interpolation between the edges of the interval (flat, above the 99th percentile of the flux, ignoring the wing of the Ly-a emission); with `template=True` it follows the template, still with no free pixel, and lies below the flux peaks at 565-590 nm (6% of the pixels above C + 2 sigma). Reported by the agent in Guido's MCP trial; not reproduced here (the spectrum is not in the repository) | numerical | — | open: phase 7 (emission model); with #465 |
| F54 | identification (tellurics) | No telluric masking: the O2 A band (759-770 nm) gives spurious doublets with good scores, CIV at z = 3.927 in J0103-1305 (4 components at the bounds, F52) and CIV z = 3.928 with MgII z = 1.733 in J0959+1312 (fitted together, reduced chi2 12). Found in the MCP pilot (trials) | behaviour | — | open: phase 7 (identification); telluric data in the pre-release `data-telluric-2026` |
| F56 | continuum (forest) | Measured on the mock benchmark (phase 4 baseline, 135 mocks): `estimate_auto` gives a forest continuum (Ly-b to Ly-a) low by 6.4% in the median, 13% at S/N 10, 12% at z_em = 4.5, with 5% of the forest pixels more than 18% low (34% at z_em = 4.5); the red side is within 0.2%. Quantifies F17 and the known low forest continuum (F53 on real data) | numerical | `benchmarks/validate_mocks.py`, baseline | improved by #465 (measured on the mock benchmark (2.1.0 + F49/F50 → #465): forest continuum 0.936 → 0.943 (5th percentile 0.820 → 0.830), red side unchanged (1.002); S/N 10: 0.875; S/N 5: 0.80, the limit of clipping). Also found: spikes falling on absorbed pixels and raised to the continuum level become anchors and raise the forest continuum slightly, i.e. the bias comes from the shortage of anchors at the continuum level. Open: phase 7 (emission model) |
| F57 | identification (completeness) | On the mock benchmark the metal doublets are found only well above the detection limit: none of the systems whose weaker line is below 3 sigma, 5% of those at 3-5 sigma, 30% at 5-10 sigma, 91% above 10; completeness over the systems above 3 sigma 0.85 for CIV, 0.69 for SiIV, 0.91 for MgII. Purity is high (0.97-0.99) | behaviour | `benchmarks/validate_mocks.py`, baseline | improved by #465 (wider window on the red side, 400 km/s instead of ~100: weak lines are flagged better): completeness 3-5 sigma 0.03 → 0.10, 5-10 sigma 0.35 → 0.44, CIV 0.857 → 0.874. Why: `find_absorbed_regions` is both the continuum step and the line detector (`identify_lines` takes its regions from `abs_mask`); a weak line pulls down the running mean of a narrow window by a good fraction of its own depth, so its pixels do not reach 2 sigma; with a window 4 times wider the line weighs 4 times less, and purity is unchanged (clipping still at 2 sigma). Direction for phase 7: decouple the detection from the continuum, e.g. search the lines on the flux normalised by the final, smooth continuum, with a filter matched to the expected line width instead of a per-pixel threshold. Open: phase 7; PR 7c.1 (F68, the doublets scored in the frame of their lags): SiIV completeness 0.690 → 0.857 on the full grid, above 10 sigma 0.76 → 0.93; PR 7c.4b (the matched filters, not yet in the pipelines): the weaker line of the metal systems at 2-3 sigma in a region 0.08 → 1.00 with the doublets (0.69 lines alone), at 3-5 sigma 0.38 → 1.00, at the price of 4-14 false regions per 1e4 km/s on the red side (F168); still open: phase 7 (7c.4c uses them); PR 7c.4c (the prototype, `method="hypotheses"`): detectable pixels outside every region 0.100 → 0.009 (full grid), completeness at 5-10 sigma 0.48 → 0.55 but above 10 sigma 0.95 → 0.87 (F169, F170) |
| F58 | `doublets_auto` (decomposition) | On the mock benchmark the systems with several components are decomposed differently from the truth: where one fitted component stands for several true ones the pulls of the matched component are 4-6.5 wide, and spurious components appear next to real ones (F49); the total column of each system is right (median difference -0.002 dex, pull width 1.0-1.1 at R >= 45000). The pipeline adds one component per candidate and never adds components from the residuals | behaviour | `benchmarks/validate_mocks.py`, baseline | open: phase 7 (Guido, October 2026; confirmed 13 November 2026) |
| F59 | fit at low resolution | At R = 8900 the total column of the metal systems is off by 0.39 dex (spread) with a pull width of 2.5, against 0.03-0.07 dex and ~1 at R = 45000-140000. Possible causes: unresolved saturated components (physics), and the convolution kernel computed at the median wavelength (rest of F10, up to +53%/-26% in width at the edges of a window) | numerical | `benchmarks/validate_mocks.py`, baseline | partly fixed by #526 (PR 3, exact convolution): at R 8900 the spread of the system column goes from 0.175 dex (master after #465; 0.39 in the 2.1.0 baseline) to 0.102, the median reduced chi2 from 2.45 to 1.61 and the flagged components from 30 to 16; R 45000 and 140000 unchanged (0.046 and 0.015 dex). After F61: 0.095 dex at R 8900, 0.041 at R 45000. What is left at R 8900 (median reduced chi2 1.6) is probably the physics of unresolved saturated components; to be looked at in phase 7 |
| F68 | `identify_lines` (`Ly_ab`) | The Ly-a/Ly-b pairs of the forest were never identified: the two windows in velocity never overlapped, and the score was always 0. Cause (corrected at the planning of 7c, 16 November 2026; first attributed to the order of the lines of `Ly_ab`, which is read consistently): `rate_doublet_candidate` measured velocities in a linear frame (`doppler_optical`) and the lags of the lines in a logarithmic one, about 4100 km/s apart for Ly-a/Ly-b; the metal doublets were misaligned by 0.4 (CIV), 1 (MgII) and 6 km/s (SiIV) (checked on clean mocks, logN 13-14.5). `Ly_ab` is in the default list of `identify_lines` and `components_auto` (not of `doublets_auto`), yet the Ly-b line of every pair became a spurious Ly-a system (z = 2.29 for a pair at z = 2.90); with `bypass_scoring` the pair was recognised. Found in phase 5 while writing the characterization tests of the line methods | identification | `test_ly_ab_pairs_are_rated`, `test_clean_doublets_score_high`, `test_doublet_score_frame_is_logarithmic` | fixed in phase 7, PR 7c.1 (16 November 2026): one logarithmic frame for the scoring of every doublet; the numbers in the row of 7c.1 |
| F72 | `rebin`, `resample_on_grid` | The boolean masks (e.g. `abs_mask`, the absorbed regions) are dropped with a warning: `get_resampled_column` interpolates only numbers (`rebin` catches its `TypeError`, `resample_on_grid` checks the type first). After a rebin or a resampling the mask of the user is lost and the pipelines make a new one. A mask could follow the new grid (e.g. a pixel masked if any old pixel in it was). Found in phase 5 with F71 | behaviour | `test_rebin_drops_the_masks`, `test_resample_drops_the_masks` | fixed in phase 7-bis, PR B7 (16 November 2026; rule *any* chosen by Guido; boolean columns only: `follow_mask`; *Co-add* still carries no column, with a warning; the integer columns, e.g. ids, are still interpolated as numbers by `rebin`: F161); was: open; phase 7 (agreed with Guido: it changes the results and needs a choice of method, how a mask follows a new grid, on the inputs of the continuum; only the recipes *Rebin*, *Resample* and *Co-add* use these methods, not the pipelines nor the mock benchmark; meanwhile, *Find absorbed regions* after a rebin); moved to phase 7-bis, PR B7 (Guido, 13 November 2026: it changes how the data are kept, not an algorithm) |
| F84 | x axis in velocity | The spectrum does not keep the reference redshift of a conversion of the x axis to velocity (`x_convert`; only the history line "X-axis converted to km/s at z_ref=..." records it): every conversion back to wavelength (rebin, smooth, continuum, combination of spectra, the ticks of the systems since F66) uses z_em, so they are wrong if the axis was converted with another z_ref, or if z_em changed after the conversion. Found in phase 5 while fixing F66 | design | `tests/test_velocity_axis_z_ref.py` (characterization first, then the new behaviour) | fixed in phase 7-bis, PR B6 (16 November 2026, the design agreed with Guido): the spectrum keeps z_ref (`SpectrumData.z_ref`, archive key `Z_REF`), the conversions back and the ticks use it, z_em only when it is not known. Not affected, contrary to the text above: the smoothing (and the continuum smoothing) needs only differences of velocity, independent of the zero (a test says it). The recipes that do not work on an axis in velocity: F159 |
| F93 | `identify_lines` (identifications) | The identification is by doublets only: two absorbed regions are a candidate when they overlap once one is shifted by the separation of the doublet in velocity (the ratio of the wavelengths), and the pair is accepted when the R^2 of F_1 = F_2^gamma (gamma from f and lambda of the two lines) is above 0.5 (metals only redward of Ly_a). Three things found while rebuilding the Identifications window: (1) the pairs are not stored: each region keeps only its candidates ([transition, score]), and the window finds the pairs again from the overlap of the redshift ranges and the shared score; (2) the multiplets with more than two lines (FeII, SiII, NI, ...) are never searched, only those with exactly two members; (3) the regions of the forest that no doublet identified are labelled Ly_a with score 0, an assignment and not an identification (lya_auto calls the recipe with "Ly_a" alone, not a doublet, so all its regions are so labelled). Found in phase 5-bis | behaviour | `test_identification_viewer.py` (the window) | to be discussed in phase 7 (identification, Guido's algorithmic choices): store the pairs with their redshift and score; the window shows the doublets as the recipe finds them and the assumed Ly_a as such Item 1 **fixed in phase 7, PR 7c.3** (18 November 2026): the pairs are stored, as the edges of the structure of the identification (`core/identification.py`, in the metadata of the spectrum and in the archive), with their regions, transitions, redshift, lag and score; the window still finds them again (its reading of the structure is 7c.7); items 2 and 3 with 7c.4; PR 7c.4c: with `method="hypotheses"` the assumed Ly-a is a kind of its own (`FOREST_ASSUMED`), and the alternatives of each group are stored with their weights (items 2 and 3 for the prototype; in the window with 7c.7) |
| F101 | script console (`gui/scripting.py`) | `AppController.run_recipe`, the `app.run_recipe(...)` of the scripts of the console, calls `MainWindowV2._on_recipe_requested` with a keyword `session_override` that the method does not take: every call fails with a TypeError. Found in phase 6-bis (PR 3a) while adapting the callers of the recipes | bug | `test_run_recipe_returns_the_new_session` | fixed in phase 7-bis, PR B1 (13 November 2026): `run_recipe` no longer goes through the dispatcher of the window (a worker thread, always `None`, the active session only) but runs the recipe on the session given with `recipes.run`; the window and its history change only when the script calls `create_tab` |
| F102 | `lya_auto`, `forest_auto` (cleanup of the narrow lines) | After removing the lines narrower than `min_b`, `_remove_narrow` calls `refit_all` without `region_limit`: every group of the session is refitted, also the components outside the forest (e.g. a CIV doublet on the red side fitted before), whose values change slightly. Found in phase 6-bis (PR 4b) while recording the reasons of the removals. `forest_auto` is not limited to a region at all: it identifies on the whole spectrum and its first refit, before the cleanup, also refits every component of the session | behaviour, numerical | `test_lya_auto_does_not_refit_the_red_side`, `test_forest_auto_does_not_refit_the_red_side` | `lya_auto` fixed in phase 6-bis (PR 4c, Guido's decision): the refits of the cleanup limited to `lya_forest`, as the pipeline; documented in the tutorial (*Refit All* fits everything again). The mock benchmark does not run `lya_auto`: a targeted test replaces it. `forest_auto`: open, Guido's decision (limiting it means choosing its region, e.g. `all_ly_forest`, for the identification and both refits); `forest_auto`: phase 7 (Guido, 13 November 2026); moved to phase 7-bis, PR B8 (Guido, 13 November 2026: a bug of the scope of the refit, not an algorithm); fixed in phase 7-bis, PR B8 (Guido's decision, 16 November 2026): `forest_auto` works in `all_ly_forest`, from the Lyman limit to the Lyman-alpha emission, for the identification (a copy of the mask restricted to the region, as in `components_auto`; `identify_lines` already looked for single Ly-a lines only there) and for both refits; test `test_forest_auto_does_not_refit_the_red_side` (it replaces the characterization `test_forest_auto_refits_the_red_side`: on its mock the CIV went from logN 13.78168 to 13.78178 and b 10.0773 to 10.0707 km/s) |
| F112 | MCP output of the pipelines (`IdentificationSummary`) | The result of `absorbers_doublets_auto` says `identification.n_systems: 0` and `systems_by_series: {}` while its notes say "3 candidate systems found" and three components are fitted (65 on J0103-1305): the summary counts the candidate systems of the metadata of the spectrum, which the pipeline does not leave there. Both agents of PR 7 noted the contradiction | output | — | open: phase 7 (Guido, 13 November 2026), with the output of the pipelines: what they leave in `RecipeResult.extra`; the immediate contradiction **fixed in phase 7-bis, PR B9** (the candidate systems of the summary are those the recipe reports in `extra['candidates']`, the same as its notes; the summary has no region and score of the systems yet), the output rebuilt in phase 7, PR 7c.6 (Guido, 13 November 2026) |
| F114 | Docker image | `absorbers_fit_bayesian` is a tool of the image, but the image does not install the `[bayes]` extra (emcee, dynesty): the tool fails with an ImportError (recorded as a gap `failed`). Found in trial 1 of PR 7 | packaging | `tests/test_docker_image.py` (the list of the image) | fixed in B4 (phase 7-bis, 16 November 2026): emcee, dynesty and corner added to the image, +~1.7 MB (3 pure-Python wheels, 175 KB; no new dependencies: numpy, scipy and matplotlib were already there); `tests/test_docker_image.py` now reads the `[bayes]` extra too |
| F115 | methods text (`methods_text`) | Seen by the agent of trial 1 of PR 7: the declared gaps (meant for the developers) are in the text for the paper, as written (decision of PR 5c); the first sentence names the spectrum by the name of the session ("The spectrum of spectrum"); the reasons open each paragraph as written, which reads oddly in a paper; `session_methods` gives the text only inside its JSON, so the agent wrote it to a file with its own Python (declared as a gap) | output, cosmetic | — | mechanical part done in B10 (16 November 2026: gaps out of the text, the generic name of the session, `session_methods` with `path`); open: the prose of the paragraphs (the reasons open each paragraph as written), after phase 7 (Guido's text, PR 5c) |
| F116 | flags of the components | No flag for components unresolved (b well below the instrumental b: 1.4-2.8 km/s at R = 45000 in J0103-1305) or saturated (MgII at logN 15.06 ± 0.24), which both agents of PR 7 judged unreliable by hand; the flags of #517 cover the bounds and the window only | behaviour | — | open: phase 7 (Guido, 13 November 2026), with F52 (the flags of #517 cover the bounds and the window only) |
| F118 | `optimize_system` (stop criterion, output) | On the mock of trial 3 (PR 7) the optimizer added nothing to a single CIV that hides two components 9.4 km/s apart, where a second component added by hand lowers the chi2 on fixed pixels from 55.4 to 37.1 (dAIC -12): it tries a component only where a residual pixel exceeds `threshold_sigma` (2.5), and the deepest was 2.1 sigma (the reason is an INFO line of the log, not in the result: the agent saw "1 component changed" and nothing else). A systematic misfit spread over many pixels is not seen. Related to F58 (the pipelines never add components from the residuals) | behaviour, output | — | open: the criterion to phase 7 (a heuristic, Guido's); decided (Guido, on PR 7): the reason of the stop (e.g. the deepest residual against `threshold_sigma`) and the AIC of each component tried in `RecipeResult.extra`, numbers unchanged, fixed in PR 7-ter (`details.optimization`: `start_aic`, `tried`, `stopped`), the criterion open: phase 7 |
| F119 | `doublets_auto` (scoring of the candidates) | The SiIV doublet of the z = 2.2 system of the mock of trial 3 (logN 12.9, b 9, S/N 50; depths 0.55 and 0.33) is found as a kinematic candidate but scored R^2 = 0.483, below the threshold of 0.5, and dropped; 0.485 with the true continuum, so the continuum (on the wing of the NV emission) is not the cause. The candidates dropped are not in the result: the agent found the SiIV by eye | behaviour, output | — | open: the scoring to phase 7 (identification); decided (Guido, on PR 7): the candidates dropped, with their score, in `RecipeResult.extra`, fixed in PR 7-ter (`details.dropped_candidates`: series, z, score, threshold, x range; not the forest doublets, which no score accepts), the scoring open: phase 7; PR 7c.1: the scoring misaligned the SiIV lines by 6 km/s (F68): a clean mock SiIV with the parameters of this system now scores 0.99 (0.56 before); the mock of trial 3 not rerun |
| F121 | `metals_auto` (acceptance of a trial ion) | `metals_auto` kept a trial ion if its fit gave logN >= 12 and 2 <= b <= 50 km/s; `aic_penalty` was accepted and not used, while its description, the `details` of the registry and the methods text said that an ion is kept if it decreases the AIC by it. Once the ions found were kept (F111), the criterion kept 192 ions of 239 on J0103-1305, 67 of them with b = 50 km/s (the bound), 9 with logN = 19 (the bound), 30 with an AIC change above -5; on the mock of trial 1 a spurious SiIV at z = 1.83 (logN 12.80, b = 50.0) | behaviour | — | partly fixed (Guido's choice, PR 7-quater): an ion is not kept with b or logN at the bounds of the fit (within 0.1 per cent) nor with an AIC change above `-aic_penalty` (against the same model without it, on the pixels of its fit), and each ion tried is reported with the reason (`tried`). Open (phase 6-quater or 7): on the mock of trial 1 a spurious CII at z = 1.83 in the forest is still kept (logN 12.58, b 18, AIC change -33); a trial ion that absorbs a Lyman-alpha line improves the AIC anyway, and the model without it keeps the neighbours refitted with it (a bias towards keeping): the comparison with the alternative identification (Lyman-alpha) is for phase 7; measured on the mocks (the benchmark does not run `metals_auto`). Measured on J0103-1305 raw with the code of master after #747 (November 2026; `edit.set_properties` z_em = 4.07, R = 45000, `continuum.estimate_auto`, `doublets_auto`: 65 components, then `metals_auto` on the six default ions, through `recipes.result.run`): 320 trials built (64 redshifts), 81 not tried for their lines outside the spectrum (as before), 239 fitted, **82 kept** (CII 12, CIV 5, NV 17, OVI 18, SiIV 30; against 192 with the old criterion, 67 of them at b = 50 km/s; none now), 157 not kept: b at its bound 108, logN below 12 29, AIC change above -5 7, b below 2 km/s 7, logN at its bound (19) 6. Time: the trials 386 s (6.4 min) and the final `refit_all` 121 s, 507.5 s the recipe, against 26.7 min for the 239 trials of the old criterion; the whole chain 8.9 min; 147 components at the end. The run interrupted at 70 per cent with the rules of F121 kept 60 of 182 tried (33 per cent), the complete one 82 of 239 (34 per cent); the rest: phase 7 (Guido, 13 November 2026) |
| F148 | `optimize_system` (window) | The window of the residuals and of the AIC is centred on the first line of the series of the seed (`lam_target_nm`): for a multiplet (CIV) the other lines are out of it (CIV 1550 is 500 km/s from CIV 1548, with a window of ±100 km/s), so a component is tried and judged on one line only, while the fits use them all. Found in PR 10-bis (phase 6-quater), not measured | numerical | – | open: phase 7 (components; changing it changes the results), registered only (Guido) |
| F149 | `optimize_system` (failed fits) | A fit that does not converge within `max_nfev` (`success` False) gives an infinite AIC: its candidate is rejected without the reason (the `tried` of the report says "AIC changes by inf"), and if the starting fit fails, no candidate can be added (inf - inf is NaN, never below `-aic_penalty`). `max_nfev` was fixed at 1000; since PR 10-bis it is a parameter. Found in PR 10-bis (phase 6-quater) | numerical, interface | `test_a_fit_that_reaches_max_nfev_adds_nothing` (characterization) | open: phase 7 (Guido: registered only); `max_nfev` documented; the reason in the report **done in phase 7-bis, PR B9** (a candidate whose fit does not converge is rejected as before, with `outcome` "rejected: the fit of the candidate did not converge within max_nfev = N", `converged: False` in its entry of `tried`, and `start_fit` in the report when the starting fit fails; test `test_a_candidate_whose_fit_does_not_converge_says_so`), the criterion in phase 7, PR 7d.3 (Guido, 13 November 2026) |
| F151 | GUI, the tooltips of the entries | The tooltips of the entries of the recipes (`RECIPE_TOOLTIPS` in `gui/main_menus.py`) are a third text for the same recipe, close to its `details` in the registry. The messages that say what a recipe did or why it cannot run (status bar, followed session, continuum required, `RecipeResult`, prerequisites, MCP descriptions) name it by its `brief`, a descriptive sentence: decided with Guido to leave them so (a menu label is for saying where to click, not what happened). Found in PR 11 of 6-quater | cosmetic | — | closed in B3 (16 November 2026): the tooltip is the `details` of the registry (Guido's choice over a `tooltip` key); `RECIPE_TOOLTIPS` removed; `details` kept short (≤ 200 characters, a test checks it); `split` and `rebin` now have a tooltip |
| F152 | `update_component(s)`, `update_constraint` (model not redrawn) | These recipes change the parameters (or the links, freedom) of components but run no fit and do not recompute the `model` column: after `update_component(logN=14.5)` the model of the session is identical to the one before (checked on the mock of `tests/full_chain.py`; same for `update_constraint`), until the next fit; the components table shows the new value while the plot and `spectrum_read` show the old model. Found while writing the workflows of the instructions of the MCP server (phase 6-quater, PR 12), which now say that these tools do not fit and must be followed by `fit_component`. Not changed: it may be by design (the GUI edits then fits); to decide with Guido whether the model should be redrawn | behaviour | — | open: phase 7-bis (Guido, 13 November 2026): to improve; the model redrawn from the new parameters without a fit, if it costs little, or the warning turned on again. The GUI behaves as the server: editing a table runs only `update_component`, and a dialog ("please refit the components") says so, but `_suppress_refit_warning` is True by default (`gui/main_window.py`), so the new value shows and the model does not change until the next fit; fixed in 7-bis B2 (16 November 2026): the model redrawn without a fit (Guido's choice). Side note: a component with an absurd `resol` (e.g. R=2) now exhausts the memory in `update_component`, as it already did in a fit: not guarded (the profile kernel scales with FWHM/pixel), to be considered with the validation of `resol` |
| F153 | z_em (no estimator) | No recipe estimates the emission redshift: in trial 2 of PR 7 (6-ter) the agent set z_em = 4.07 on J0103-1305 from the Ly-a peak by eye, and declared it as a gap; every pipeline depends on z_em (the forest, the region beyond the emission, the proximity zone of F55). Registered at the planning of phase 7 from the trials (sections 4.2, 4.3) | capability | — | open: phase 7, PR 7b.3 (13 November 2026) |
| F154 | evaluation of a model (fixed pixels) | No recipe gives the residuals, chi2, AIC and BIC of a model on a fixed set of pixels: the results give the chi2 of a group only when components are added, on a window that changes with the model, so two models cannot be compared. In trial 3 of PR 7 (6-ter) the agent computed them with its own Python (a declared gap) to choose between alternatives; the import of a VPFIT model needs the same (section 6). Registered at the planning of phase 7 | capability | — | fixed (phase 7, PR 7a.5, 16 November 2026): `absorbers.evaluate_model`, see the row of 7a.5 |
| F155 | tellurics (no model) | No telluric model is used by the identification or the fit (F54 is about the spurious doublets of the O2 A band); the agents of the trials asked for one (section 4.3). Telluric data in the pre-release `data-telluric-2026`. Registered at the planning of phase 7 | capability | — | open: phase 7, PR 7c.5 (13 November 2026) |
| F156 | upper limits | No recipe gives an upper limit on the column of a species not detected; the agents of the trials noted it (section 4.3). Registered at the planning of phase 7 | capability | — | open: phase 7, PR 7d.7 (13 November 2026) |
| F157 | results of the fits (chi2) | The result of `fit_component` and `optimize_system` gives no chi2: the reduced chi2 of the group is in the summary of the components added only, not of those changed (`ComponentChanged` in `mcp_server/outputs.py`). Seen in trial 3 of PR 7 (6-ter), left open there | output | `test_a_refitted_component_reports_the_chi2_of_its_group` | fixed in phase 7-bis, PR B9 (16 November 2026): `ComponentChanged.chi2r` (reduced chi2 of the group after the call, null if not fitted) and the `fit` summary also when a fit gives a new chi2 to components already in the session; `Changes` and `COMPONENT_FIELDS` unchanged |
| F158 | plots for the agents | `spectrum_plot` (MCP) draws the flux, the continuum, the model and the ticks, without the residuals; there is no plot of several transitions in velocity around a redshift (the velocity plot of the System Inspector is GUI only). The agent of trial 3 of PR 7 (6-ter) asked for both | output | — | fixed in phase 7-bis, PR B11 (16 November 2026) |
| F159 | x axis in velocity (recipes) | Most recipes fail on an axis in velocity, with a clean `RecipeFailed` (no crash): `rebin` always (the new grid is taken in nm for the auxiliary columns, `x_final.to_value(au.nm)` in `spectrum_combine.py`, even with no auxiliary column), `resample`, `split`, `estimate_auto` and `find_absorbed` (a comparison with the Lyman-alpha in nm), `fit_powerlaw` (no point in its regions in nm); the chi2 of the forest and of the red side in the session info compares the axis with a wavelength. What works: `x_convert`, the smoothing, `calculate_running_std`, `smooth_column`, `stitch`, `coadd`, `equalize` (axes with the same z_ref, B6). With the z_ref kept (B6) they could convert the axis to nm and back. Found in phase 7-bis, PR B6 | design | `test_rebin_of_an_axis_in_velocity_fails` (characterization) | open: to be reviewed with Guido (16 November 2026; the user docs say to convert the axis back first) |
| F160 | `estimate_auto` (BAL troughs) | Across a BAL trough the automatic continuum follows the absorption: median ratio to the true continuum 0.22 on the trough pixels of the grid `robustness`, and the 5th percentile of the ratio <= 0 in 3 of 72 mocks (a continuum at or below zero, so a normalised flux infinite or negative there). Found by the metrics of PR 7a.2 | numerical | `validate_mocks.py --grid robustness` (`phase7.continuum.bal`) | open: phase 7, 7b (16 November 2026) |
| F161 | `rebin` | The integer columns (e.g. the ids of the absorbers, `abs_ids`) are interpolated as numbers: after a rebin they are floats with values between two ids (an id that does not exist), and the same for any column of codes. A column of identifiers has no meaning when averaged; it should follow the new pixel as the masks do (e.g. the id of the old pixel with the largest overlap) or be dropped with a warning. `resample_on_grid` interpolates them too. Found in PR B7 (F72), characterized by `test_rebin_interpolates_the_integer_columns` | behaviour | `test_rebin_interpolates_the_integer_columns` | open (16 November 2026): the rule to choose with Guido; the pipelines do not call `rebin` |
| F162 | `identify_lines`, `lya_auto`, `forest_auto` (`z_em_margin_kms`) | `z_em_margin_kms` (F55) applies only to the metal doublets: the single Ly-a lines are looked for only in the forest map of `identify_lines` (`core/spectrum_lines.py`, from the Lyman limit to the Lyman-alpha emission), and the regions of `lya_auto` (`lya_forest`) and `forest_auto` (`all_ly_forest`) stop there too. So no Ly-a beyond z_em is ever found, and the parameter of these two pipelines has no effect, while its description says "How far beyond z_em a candidate system can lie" (associated absorbers). Found in phase 7-bis, PR B8 (F102) | behaviour, interface | — | open (16 November 2026, Guido: recorded, not changed in B8): the choice between letting single Ly-a lines up to the margin (the forest map and the regions widened by it) and saying in the descriptions that the margin is for the doublets only |
| F163 | velocity plot (System Inspector) | The velocity plot of the System Inspector (`gui/velocity_plot.py`) and the plots of the MCP tools (`astrocook/plotting.py`, B11) draw the same panel (flux, error, model, ticks of the components, residuals in sigma) with two pieces of code. The formulas that must agree (velocity of a line, residuals in sigma, the colour of the flags) are in `plotting.py` and a test compares the Inspector's formula with them (`test_the_inspector_ticks_agree_with_the_shared_formula`, `test_flag_colours_agree_with_the_inspector`), but the Inspector does not call them: it would change pixels and layout of the GUI. Found in PR B11 (F158) | design | `test_the_inspector_ticks_agree_with_the_shared_formula` | open (16 November 2026, Guido: agreed to share only the pure helpers and not to touch the Inspector in B11): the Inspector could draw a panel with the functions of `plotting.py`; to be done with the screenshots of the tutorials |
| F164 | credibility of the alternatives (no calibration metric) | The goal of section 5 is a collection of models whose alternative identifications carry a credibility; nothing in the benchmark checks that a credibility means what it says (of the hypotheses given credibility 0.8, about 80% right on the mocks). Completeness and purity do not measure it, and the reliability reported today is weak: the flags mark 7% of the bad components, and 56% of the unflagged ones are bad (PR 7a.2). Proposed: a calibration metric in `validate_mocks.py` for the credibility of 7c (a reliability diagram over bins of credibility against the truth of the mocks, its expected calibration error and Brier score), as also for the flags of the components (7d); since the credibility rests on chi2-based criteria (AIC, BIC, SpIC of 7a.5), it inherits any miscalibration of the error spectrum (correlated pixels of rebinned reductions), which the metric would expose on the twins and the development spectra. Raised in a design conversation with Guido (16 November 2026) | measurement | — | open: phase 7; agreed with the design document of 7c.2 (Guido, 18 November 2026: a criterion of "better" of 7c.8, measured by 7c.4, `design/identification.md` section 11; the weights tempered for a chi2 that is not calibrated, section 3.3); PR 7c.4c: the calibration in `validate_mocks.py` (`phase7["calibration"]`: reliability diagram, ECE, Brier, by kind, S/N and R): the prototype ECE 0.022, Brier 0.038 on the full grid, 0.060 / 0.079 on robustness (0.096 at R 8900); without the refinement 0.18 (F170) |
| F165 | `.acs2` save (identifications of the regions) | After a round trip `region_identifications` is a dict (F27), and a second save warned "The metadata 'region_identifications' (a dictionary) cannot be saved in the archive", though it is saved, in `_meta.json`: `_convert_spec_data_to_table` took it out of the header after `_fits_meta`, which warns about every dict. Found in PR 7c.3 | output (a false WARNING) | `test_a_second_save_has_no_spurious_warning` | fixed in PR 7c.3 (18 November 2026): taken out before `_fits_meta` |
| F166 | fine grid of the model (`fitting/lsf.py`, `B_SAMPLES`) | The model of a narrow, strongly saturated line depends on where the samples of the fine grid fall: the grid has a fixed step (a third of the narrowest b, half the instrumental sigma) but its origin follows the pixels given (`PixelLSF`, the median of their wavelengths and the margin of the kernel), so the same line on the whole spectrum and on a cut of it is computed at another phase of the samples. A CIV at logN 16.8, b 6.6 km/s at R 8900 (step 2.2 km/s, b / 3) changed by up to 0.055 sigma per pixel at S/N 100 (5e-4 of the continuum, 17 pixels above 1e-3 sigma; chi2 644.4 against 646.1), against the 1e-5 of the continuum that the comment of `B_SAMPLES` states also for saturated lines. Neither of the two is more right; the error of the discretization is larger than stated for these lines (the comment was measured on less saturated ones?). Effects: the statistics of `EvaluationWindow` and of the whole spectrum differ for such lines (37 of 2400 hypotheses of the experiment of 7c.2 by more than 1e-4 of the chi2, at most 0.45%), and a fit or a comparison of hypotheses depends at that level on the extent of the spectrum. Found in PR 7c.4a (Guido: a finding, 18 November 2026) | numerical | — | open: to measure the error of the fine grid against a dense grid by logN, b and R (the saturated and the damped lines), then decide whether the step must follow the saturation (e.g. a finer step where tau0 is large) |
| F167 | `VoigtFitter.fit` (time of the pipelines) | Every call of `fit()` ends by computing the model of the whole spectrum with all the components (`voigt_fitter.py`, "Final Flux Model (Un-sliced)"), after iterations that are already local (the pixels of the fitted lines and 50 more on each side). On a mock of 32 thousand pixels with 400 components, a short fit of one component took 0.16 s, of which about 0.12 s the final model (the optical depth of the 400 components and the convolution on the whole spectrum): its cost grows with the components times the length of the spectrum, and the pipelines pay it once per group fitted (`fit_auto`, `refit_all`, the trials of `components_auto`, `SystemList` in `system_list.py`), also where the model returned is discarded (`refit_all` and `system_list.py` take it as `_`). Proposed: measure the share of this time in the pipelines (`components_auto` on the mocks of the benchmark), then compute the final model only where it is used, or update the model of the spectrum only on the pixels the fitted components reach (with the fixed components, as `EvaluationWindow` of PR 7c.4a). Found in PR 7c.4a (Guido: a finding, 18 November 2026) | performance | — | open: phase 7d (the fit) |
| F168 | matched filters on the continuum of `estimate_auto` (`core/line_detection.py`) | The one-sided bias of the continuum (F60) is summed by the broad filter: with the continuum of `estimate_auto` (GUI defaults) high by about 0.06-0.2% in the median on the red side (about 0.06 sigma per pixel), the median S of the b = 40 km/s filter is +0.2 to +0.7 instead of 0, and the regions (S > 1 with a seed above 3) cover 19-21% of the red-side pixels against 7% with the true continuum on the same mocks (S/N 100 R 45000, S/N 30 R 140000; at S/N 10 R 8900 4.6% in both). On the benchmark grids most false regions of the lines alone (4.1 per 1e4 km/s on the red side, against 0.8-3 on pure noise) come from it. Found in PR 7c.4b | behaviour | `benchmarks/validate_detection.py` | open: to decide in 7b (Guido, 19 November 2026), among the continuum of 7b (unbiased), a local baseline of S subtracted (e.g. a running median of S over a few thousand km/s), or a higher extent threshold for the broadest filter |
| F169 | R^2 pre-filter of the hypotheses on the regions of the matched filters (`fitting/hypotheses.py`, `pair_score`) | The pair score of a candidate (the R^2 of the two lines in their regions, as `rate_doublet_candidate`) is computed on the regions of the matched filters (7c.4b), which are wider than those of the mask, the more so with the bias of F168: a true doublet whose two lines fall in one region, or whose regions also hold other lines, gets R^2 <= 0 and never becomes an edge. With the pre-filter at 0 (every matched pair an edge) the completeness comes back (CIV 0.85 on a subset of 12 mocks) but the purity falls (CIV 0.56, MgII 0.52), the ECE rises to 0.20 and the time to 36 s per spectrum (10 minutes at R 140000). Found in PR 7c.4c | behaviour | `benchmarks/identification_thresholds.py` | open: 7c (a score of the pair on the peaks of the de-blending rather than the regions, or the regions of 7b once F168 is solved) |
| F170 | identification and fit (the refinement of the hypotheses, `_refine`) | The components of a hypothesis are the peaks of the de-blending, fewer than the true ones in a saturated or blended system: then a metal hypothesis cannot fit both lines with one ratio and loses against the free lines of the default labels, which fit any profile. The refinement by the residuals (up to 3 components added to each finalist) makes the comparison fair, but it is modeling inside the identification. Without it (`--max-added 0`) the completeness falls (grid full: CIV 0.62, SiIV 0.55, MgII 0.60, against 0.80, 0.80, 0.85 with it and 0.87, 0.86, 0.91 with `doublets`; robustness: 0.55, 0.48, 0.55 against 0.68, 0.73, 0.73 and 0.81, 0.78, 0.81), the fit after it finds fewer true components (full 846 against 1494 and 1242; robustness 1175, 2068, 1831) and the credibility is not calibrated (ECE 0.18 and 0.20 against 0.02 and 0.06: the hypotheses of weight 0.9-1 right in 69% of the cases against 98%); the refinement costs about a third of the time (full 2047 s against 3069 s; robustness 4875 s against 6707 s), and on robustness it makes part of the false SiIV (purity 0.80 without, 0.74 with). The general question (Guido, 20 November 2026): should the identification and the fit be separate at all? The comparison of identifications needs fitted models, the fit after the prototype is much better (pull widths 1.6-2.1 against 3.9-5.2), the automatic fitters of the literature (7c.2) identify and fit together, and the goal of v3 is fitted models with alternative identifications; the separation comes from the V1 chain. Guido's leaning: a joint recipe (identification and fit), with these requirements, so that it does not take flexibility from the agents of the MCP server (a recipe validated on the statistics of the mocks need not be the best for one spectrum): (1) a pipeline of atomic recipes, each exposed and logged (the detection of 7c.4b, the groups and alternatives, the evaluation and short fit of an alternative, the acceptance), as `doublets_auto` and `components_auto` are now; the current recipes kept; (2) its output all the alternatives with their weights (the structure of 7c.3, F164 calibrating them), and a recipe to accept another alternative of a group; (3) the same yardstick for an agent's own models (`evaluate_model` on fixed pixels); (4) the MCP descriptions present the two ways with their trade-offs, no preferred one (the default of 7c.8 is that of the GUI pipeline). `doublets_auto` (and maybe `lya_auto`) would then be superseded, kept for the replay of the logs. Open to an alternative proposer of the hypotheses by machine learning (a classifier trained on the mocks, as the classifiers of DLAs and MgII of 7c.2), judged by the same comparison. Found in PR 7c.4c | behaviour | `benchmarks/validate_mocks.py --max-added` | open: to decide in the planning of 7c (whether and how to join the identification and the fit; the refinement kept, an explicit option, or replaced by better proposals of the components), before 7c.4d |
| F171 | the short fit of large hypotheses (`fitting/hypotheses.py`, `MAX_FIT_K`) | The Jacobian of the short fit is by finite differences (one model per free parameter), so its time grows as k times the evaluations: in a DLA trough (mock "combined" of the grid robustness, S/N 10, R 45000, z_em 4.5) a group had hypotheses of 94 components (k = 281), and one mock did not end in 100 minutes. Fixed in PR 7c.4c by a cap (60 free parameters, Guido, 19 November 2026): a group with a hypothesis above it is not compared and keeps its default labels with the AOD parameters (6 s on that mock); the regions of a damped line are still split by the de-blending into many peaks. Found in PR 7c.4c | performance | `benchmarks/validate_mocks.py --grid robustness` | open: 7d (an analytic Jacobian, F167; the damped lines as one component of 7b) |

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
automatic re-analysis of the 28 UVES quasars of D'Odorico et al. (2010),
in the reductions of SQUAD DR1, assisted by an agent through MCP,
compared with the published analysis and with *synthetic twins* of the
spectra where it gives a model (a mock with the published components, on
the grid, with the error spectrum, resolution and continuum of the real
one). Guido, 13 November 2026; before, the UVES Deep Spectrum of
HE0940–1050 alone, compared with D'Odorico et al. (2016).

**Rule: the test part of the sample stays blind** (20 spectra, among
them HE0940–1050; the split in phase 7, "Data and blindness") until the
version of the paper is frozen (after phase 7). Astrocook is developed
and measured only on the mocks, the development part and the twins,
never on the test spectra, so that it is not tuned on them; the analysis
plan is written before the first run on them.

**What the development owes to the paper:**

| What | When |
|---|---|
| Twins: a mock from an arbitrary component list, on a given pixel grid, with a given error spectrum, resolution and continuum (`astrocook/mocks` has a uniform velocity grid, a power-law continuum and an error from the S/N today) | small, when the material of the published analysis arrives |
| Import of a VPFIT model (fort.26 or the tables of a paper) into a `SystemListV2`, and its evaluation on a spectrum without fitting (chi2, residuals, AIC/BIC per system) | new, small |
| Full MCP server, with the log of the decisions of the agent | phase 6-ter |
| Emission model with an uncertainty; identification with alternatives and credibility | phase 7 |
| Mock benchmark rerun with the frozen version | after phase 7 |

The twins of the development part also give phase 7 realistic
development cases, which the mocks of phase 4 do not cover.

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
| 7 – Robust components | in progress: plan agreed (13 November 2026), four blocks 7a-7d |
| 7-bis – Small fixes | done (16 November 2026), interleaved with phase 7: PRs B1-B11 (#772, #774-#776, #778-#780, #786, #788, #790, #791); F101, F114, F30, F72, F84, F102, F112 (the immediate fix), F115 (the mechanical part), F128, F149 (the reason), F151, F152, F157 and F158 closed; F159, F161, F162 and F163 found |
