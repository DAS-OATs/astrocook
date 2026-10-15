# Astrocook – Development Roadmap

Working document, based on an analysis of the v2.0.1 code (October 2026).
It collects the assessment of the current state, the consolidation plan for
v2 and the vision for a future v3. Update it at the end of each phase.

---

## 1. Current state (v2.0.1)

### What it is
A Python package for the analysis of quasar spectra, with a PySide6 GUI and
a scriptable API. It covers I/O (FITS, ESO/MIDAS, V1 `.acs` and V2 `.acs2`
archives), flux operations, continuum estimation, line identification and
Voigt fitting (least squares + experimental Bayesian engine with
emcee/dynesty).

Size: ~53k lines of Python, of which ~26k are active V2 code and ~22.5k are
in `astrocook/legacy/`.

### Layout
- `astrocook/core/` – immutable data structures (`structures.py`),
  `SpectrumV2` (`spectrum.py`), pure functions (`spectrum_operations.py`),
  session and system list.
- `astrocook/fitting/` – `voigt_fitter.py`, `voigt_model.py` (constraints),
  `bayesian_fitter.py`.
- `astrocook/recipes/` – user-facing operations (continuum, absorbers, edit,
  file, flux), with structured parameters.
- `astrocook/gui/` – PySide6 interface.
- `astrocook/io/` – loaders and V1 → V2 migration.
- `astrocook/legacy/` – V1 code (wxPython), still partly imported.

### Strengths
1. Immutable data structures (frozen dataclasses): reliable undo, consistent
   history, no side effects.
2. Explicit physical units (astropy) on every column.
3. Layered separation: core/fitting, recipes, GUI.
4. Reproducibility: session log with a scripter; automatic V1 migration.
5. Mature Voigt fitting, with a Bayesian engine and asymmetric uncertainties.
6. Numpydoc docstrings, Sphinx documentation with tutorials, curated
   CHANGELOG, macOS and Docker builds.

### Weaknesses
1. **Legacy still coupled:** 17 V2 modules (including `core/session.py` and
   `core/spectrum.py`) import from `astrocook.legacy` (`Defaults`, `Format`,
   `GUILog`, `SpectrumV1`, `SystListV1`, `trans_parse`, `xem_d`,
   `qso_composite`, …).
2. **Few tests:** ~5 `test_*` functions; the rest of `tests/` are
   verification scripts. No CI. No numerical regression tests.
3. **Very large files:** `gui/main_window.py` (4.6k lines),
   `gui/pyside_plot.py` (2.9k), `core/spectrum.py` (2.8k),
   `gui/system_inspector.py` (2.1k). `SpectrumV2` is a "god class".
4. **Broad error handling:** ~200 `except:` / `except Exception` outside the
   legacy code. Inconsistent APIs (e.g. `load_session_from_file` returns `0`
   on failure).
5. **Packaging:** no LICENSE file (field commented out in `pyproject.toml`);
   template placeholders still present; `astrocook-gui` entry point
   commented out; `numpy>=1.17.3` is unrealistic with Python ≥3.12;
   emcee/dynesty/corner are required instead of optional.
6. **Repository clutter:** `docs_backup/`, `*_old.py` files, `*.bck`,
   `telluric_old.fits`, `atom_par_backup.dat` (some shipped in the package).

---

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
Current algorithm: product of significances on a common z grid
(`compute_identification_signal`, from V1); for doublets, kinematic overlap
plus a score based on the relation F₁ = F₂^γ (`rate_doublet_candidate`).

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

### Phase 0 – Environment
Install dependencies; run the GUI with `QT_QPA_PLATFORM=offscreen`.

Done: SessionStart hook for cloud sessions, `pytest-qt`/`pytest-cov` in the
`[test]` extra. The GUI starts offscreen without errors. Baseline of the
existing suite (Python 3.13): 1 passed, 4 failed; all four failures come
from tests that are out of date with the current API (missing
`ComponentDataV2.id`, `rebin_spectrum` now returning 6 values, free
parameter vector without `btur`). Fixed at the start of phase 1: the suite
is green.

Linter: Ruff, pyflakes rules (`F`) only, `legacy/` and backup directories
excluded. First check: 176 warnings (121 unused imports, 24 unused
variables, 16 f-strings without placeholders, 9 undefined names, 6
redefinitions). Three undefined names are real bugs, to be reported to
Guido before fixing:
- `recipes/file.py` (`export_ascii`): `re` not imported, so the recipe
  always raises `NameError`; probably the cause of #476/#481.
- `core/session.py` (V1 `.acs` save): `dataclasses` not imported; the error
  is swallowed by a broad `except`, so the V2 log is silently replaced by an
  empty one. Adding the import alone does not help: `HistoryLogV2` is a
  dataclass without fields, so `asdict()` returns `{}`. Guido's decision:
  no backward support for V1 is needed, the log may be lost (also to
  encourage saving as `.acs2`), but the user must be told. Done: warning in
  the log, and a confirmation dialog in the GUI when saving as `.acs`.
- `fitting/voigt_fitter.py`: `iterator` used before assignment in the debug
  branch taken when the fit mask is empty, so the fit raises
  `UnboundLocalError`.
The `export_ascii` and `voigt_fitter` bugs were fixed right after phase 0,
with regression tests.
The other undefined names are string annotations (`'SpectrumV2'`) without a
`TYPE_CHECKING` import: harmless.

### Phase 1 – Safety net (pilot)
No changes to production code.
1. **Numerical snapshots** (characterization tests) of `SpectrumV2`
   operations and pure functions: outputs of the current code are stored
   and compared with explicit tolerances.
2. **Replay of session logs** as end-to-end tests (requires real sessions
   provided by Guido).
3. **Simulated spectra** with known parameters (extend
   `tests/simulate_and_fit_cii.py`). Design the generator with
   completeness/purity for v3 in mind.
   First part done: `astrocook/mocks/` (new module, no change to existing
   code), independent of the V2 fitting code (own Voigt profile, exact
   convolution in velocity space, pixel integration by oversampling);
   ground truth in `MockTruth` (JSON). Validated in `tests/test_mocks.py`
   (thin-line equivalent width within 0.1% of the analytic value,
   convolution conserving the equivalent width, instrumental FWHM = c/R,
   noise statistics, oversampling convergence). Second part done: random
   Lyman-alpha forest (f(N) ~ N^-1.5 for 12 < logN < 17, log-normal b
   around 25 km/s with b >= 10 km/s, dN/dz = 1700 ((1+z)/4)^2.5, which
   reproduces the mean transmission of Faucher-Giguere et al. 2008 within
   ~0.03 at z = 3 and 3.8), multi-component metal systems, Gaussian
   emission lines (Lyman-alpha + NV, CIV) and single-pixel spikes. The
   default parameters are indicative, to be calibrated before scientific
   use. Still missing: sky residuals and tellurics.
4. **Non-determinism:** fixed seeds, statistical tolerances, UUIDs and
   timestamps excluded from comparisons.
   First case found: the automatic pipelines depend on the UUIDs and on
   the string hash seed (F19); the full-chain snapshots fix both
   (`tests/full_chain.py`).
5. **File compatibility:** `.acs`/`.acs2` round trips (requires real V1
   archives).
   `.acs2` done (`tests/test_acs2_round_trip.py`): a session with
   everything the format is meant to keep (units, float/bool/int auxiliary
   columns, a resolution column, z_em, R, object name, components with
   uncertainties, btur, chi2 and resolution, a custom multiplet, frozen and
   linked parameters, identification metadata), the small J0103-1305
   archive and (slow) the 15 MB one; a second round trip is identical to
   the first. Spectrum, components and constraints come back unchanged.
   Lost or changed: the V2 session log (F26, data loss, reported at once),
   `region_identifications` (F27), the Bayesian samples (F28), the metadata
   key `resol` (back as `RESOL`, the resolution itself is recovered from
   `resol_fwhm`; F8). The state before the changes is stored in the archive
   (`*_orig.fits`) but not loaded back. The pre-existing
   `tests/test_v2_io_integrity.py` catches its own assertion errors and
   always passes; superseded by the new tests, to be removed in phase 2.
   `.acs` (V1) waits for real V1 archives.
6. **Measuring the net:** coverage; mutation testing (`mutmut`) on critical
   modules.
   Coverage measured (October 2026, `pytest --cov=astrocook`, reported in
   the Python 3.13 CI job; the child processes of the full-chain snapshots
   are measured too, `patch = ["subprocess"]`). Statement coverage of the
   V2 code (without `legacy/`): 38%. By area: `mocks/` 99%, `fitting/` 72%
   (`voigt_fitter.py` 84%, `voigt_model.py` 69%, `bayesian_fitter.py`
   14%), `core/` 58% (`spectrum_operations.py` 80%, `spectrum.py` 55%,
   `system_list.py` 41%, `system_list_migration.py` 12%), `recipes/` 37%
   (`absorbers.py` 45%, `continuum.py` 48%, `file.py` 51%, `edit.py` 23%,
   `flux.py` 22%), `io/` 27% (`adapter.py` 64%, `loaders.py` 8%), `gui/`
   25%. The slow tests add less than 1%. Main gaps, by risk:
   - **loaders** (`io/loaders.py`): only `.acs2` is exercised; ESPRESSO
     S2D, X-Shooter, PypeIt (spec1d, coadd, stack, telluric model), generic
     FITS, ASCII with RESVEL header and format detection have no test.
     Needs small sample files of each format (from Guido).
   - **recipes** never run: in `edit.py` `import_systems`, `stitch`,
     `coadd`, `trim_common`, `extract_preset`; in `absorbers.py`
     `metals_auto`, `clean_negligible`, `add_linked_system`,
     `optimize_system`, `fit_bayesian`, `update_constraint`; in `flux.py`
     `resample`. The SpectrumV2 methods behind several of them are covered,
     the recipe layer (parameter parsing, chaining, error handling) is not.
   - **system list** operations: `optimize_hierarchy`, `update_constraint`,
     `filter_by_range`, `delete_component`.
   - **SpectrumV2** methods: `rebin` (only the pure function is tested),
     `apply_expression`, `mask_expression`, `smooth_column`,
     `update_continuum_from_knots`.
   - **V1 migration** (`system_list_migration.py`, `io/v1_stubs.py`): waits
     for V1 archives.
   The GUI smoke tests (item 7) and the `.acs2` round trips (item 5) come
   next and address part of `gui/` and `io/`. The gaps are kept in mind and
   addressed later (agreed with Guido), except for targeted recipe-level
   tests, done first.
   Targeted tests, first part (`tests/test_recipes_edit_flux.py`): the
   recipes of `edit.py` and `flux.py` (parameter parsing as from the GUI
   dialogs, lookup of other sessions by name, chaining, return value on bad
   input; snapshots for smooth column, co-add and rebin). Coverage:
   `edit.py` 23% -> 82%, `flux.py` 22% -> 84%, `spectrum.py` 55% -> 73%;
   V2 code 38% -> 42%. Found F22 and F23.
   Second part (`tests/test_recipes_absorbers.py`): the single-component
   recipes of `absorbers.py` (add, update, delete, fit, constraints and
   links, linked systems, `optimize_system`, `metals_auto`,
   `clean_negligible`, `fit_bayesian`) and the system list operations
   behind them, on the synthetic CIV blend; fitted values compared with the
   known input within loose tolerances (machine-dependent, see item 1).
   Coverage: `absorbers.py` 45% -> 79%, `system_list.py` 41% -> 76%,
   `bayesian_fitter.py` 14% -> 70%, `voigt_model.py` 69% -> 85%; V2 code
   42% -> 46%. Found F24 (reported at once: loses data). Note: with the
   default threshold (2.5 sigma) `optimize_system` does not add the second
   component of the blend (logN 13.1, b = 6 km/s, 11 km/s away, S/N 50);
   with 1.5 sigma it recovers both within 3 km/s and 0.2 dex.
   Mutation testing, pilot on `core/spectrum_operations.py` (`mutmut`,
   configured in `pyproject.toml`, run by hand: ~7 minutes with 4
   processes; GUI and full-chain tests left out): 1621 mutants, 147 not
   reached by any test, of the others 75.8% killed. Most survivors were in
   degenerate input never tested (bad pixels in `rebin_spectrum`, no
   smoothing, windows as long as the array, exact thresholds, no overlap in
   `compute_flux_scaling`): covered by `tests/test_spectrum_operations_edges.py`,
   77.7% killed. The remaining survivors are mostly in the identification
   helpers (`merge_regions_by_velocity`, `rate_doublet_candidate`,
   `find_absorbed_regions`), checked only through snapshots of the final
   result, and equivalent mutants. The 147 unreached mutants are in
   `compute_identification_signal` and `find_signal_peaks`, which nothing
   calls (dead code, F31). Found F30.
   Fitting modules (`fitting/voigt_model.py`, `fitting/voigt_fitter.py`,
   ~40 minutes): 1859 mutants, 226 not reached (mostly `_smart_guess`,
   commented out in `fit`, and the V1 constraint conversion), 68.8% of the
   others killed. Survivors that matter: the turbulent b could be changed
   (e.g. sqrt(b^2 - btur^2)) unnoticed, and no test fitted a line broader
   than 50 km/s; both covered by `tests/test_voigt_fitter_physics.py`,
   which found F32. Many other survivors are in resolution handling (F8,
   F10: to be tested together with the rewrite of the convolution), the
   feature limits of the fit window and log/error messages.
7. **GUI:** smoke tests with `pytest-qt` in offscreen mode.
   Done (`tests/test_gui_smoke.py`, on the small J0103-1305 archive
   `J0103-1305_metallines_lya_new_3.acs2`): the main window starts empty;
   a session opens from a file (missing or invalid files give a message,
   no session); every recipe dialog of the menus is built (the hidden
   recipes are called by the inspector with their own parameters); recipes
   run through the GUI path (string parameters, worker thread, result
   processing, history), including one that opens a new session and the
   z_em requirement check; errors are reported to the user and the window
   stays open (a recipe returning 0 and one raising, F23); undo and redo;
   system inspector, log and script console open; a session closes.
   Message boxes are replaced by stubs that record what would be shown.
   Not covered: what is drawn, mouse interaction on the plot, editing in
   the inspector, the continuum editor. Coverage: `gui/` 25% -> 45%
   (`main_window.py` 52%, `recipe_dialog.py` 69%, `system_inspector.py`
   48%, `pyside_plot.py` 31%); V2 code 46% -> 56%. The suite now takes
   ~2.5 minutes (was ~1). Found F25.
8. **CI:** GitHub Actions (Python 3.12, 3.13). Done (brought forward): tests on both versions;
   Ruff blocking on undefined names (`F821`) only, until the other pyflakes
   warnings are cleaned up in a dedicated PR.
9. **First scientific checks:** oversampling of the Voigt model; scaling of
   uncertainties by the reduced χ².
   Preliminary result on oversampling (mock generator vs V2 model, CIV,
   R = 45000, 2.5 km/s pixels): EW bias of the V2 model +0.05% at
   b = 12 km/s, +0.3% at 3 km/s, +4% at 1 km/s (max flux difference 0.001,
   0.013, 0.115). See F9.
   Scaling of uncertainties, done (`tests/test_error_scaling.py`): Monte
   Carlo on mocks, one CIV component (logN 13.5, b = 12 km/s, R = 45000,
   2.5 km/s pixels, S/N 20), 100-200 noise realizations. With the right
   error spectrum the errors of `VoigtFitterV2` are right: pull standard
   deviation 0.96 (z), 0.98 (logN), 1.05 (b), reduced chi2 0.95. With the
   error spectrum overestimated by 2 (reduced chi2 0.24) the errors are 2
   times too large; underestimated by 2 (reduced chi2 3.8), 2 times too
   small. Scaling by sqrt(reduced chi2), as V1 did, gives pulls of 1.00,
   1.02, 1.09 in all three cases. No bias of the fit: the mean z offset
   (+0.07 +- 0.03 km/s at S/N 20) goes to 0.0001 km/s without noise; b and
   logN are off by +0.025 km/s and -0.0005 dex at high S/N (sampling, F9).
   See F29.
10. **Issues labelled "possibly fixed"** ([#469](https://github.com/DAS-OATs/astrocook/issues/469), [#477](https://github.com/DAS-OATs/astrocook/issues/477), [#455](https://github.com/DAS-OATs/astrocook/issues/455),
    [#450](https://github.com/DAS-OATs/astrocook/issues/450), [#436](https://github.com/DAS-OATs/astrocook/issues/436)): write a test reproducing each reported scenario.
    If it passes, the test is the evidence for closing the issue; if it
    fails, the issue goes back to the open bugs (see section 3.1).
    Done (tests `tests/test_issue_*.py`):
    - #469 blended lines of different systems: grouped and fitted together;
      the inspector highlights the blended component. The blend threshold
      ignores the resolution (F14).
    - #477 / #436 VIS+NIR stitching: one cut-off for two arms, scaling
      recovered within 1% (fitted range, auto, manual factor); the GUI
      passes the other session and asks for N-1 values. With tellurics in
      the overlap the factor is biased by ~1.5% (F15).
    - #455 multiplet subsets: the subset is modelled and fitted alone;
      NiII added. Still open: 10 metal species with 3+ transitions in
      `ATOM_DATA` have no multiplet entry (question pending in the issue).
    - #450 link colours: work within the selected group; stale colours
      after linking with the group still selected (F13).
    Closing the issues is Guido's decision: the reporters of #469, #477,
    #436 and #455 have been asked to confirm with the current version.
11. **V1/V2 parity tests.** The V1 code (`astrocook/legacy/`) is the
    validated reference behind ~50 published papers. For each operation
    reimplemented in V2 (rebin, coadd, continuum, Voigt model and
    convolution, identification), run V1 and V2 on the same input and
    compare the outputs within a tolerance; every difference is either
    documented as intended or treated as a bug. Scheduled right after the
    snapshots of the remaining `SpectrumV2` methods, before the mock
    generator. A first manual check (October 2026) showed that F1, F3 and
    F10 are V2 regressions: V1 behaves correctly in all three cases.
    Done (`tests/test_parity_v1_v2.py`): rebin flux and coadd flux agree
    to 1e-10 on the same bins (V1 takes the first bin centre as start, V2
    its left edge); equalization factors agree within 0.5%. Documented
    differences: rebinned errors (F1), clipping (F2), smoothing width
    (F16). Voigt model (V1 `lines_voigt` + `psf_gauss`, as in `SystModel`):
    agrees within 5e-4 in flux at R = none, 45000, 20000; a constant 0.05%
    equivalent-width offset remains for CIV (identical atomic data, not
    traced; 2e-4 dex in logN), and V2 uses f = 0.4164 for Ly-alpha (V1
    0.416); with a resolution column V2 fails (F10). Continuum (J0103-1305):
    red side agrees, forest differs (F17). Identification: skipped, as
    agreed. V1 `identify_unknown` is unfinished and needs manual input;
    detection alone (V1 `find_lines`, minima at 5 sigma in the raw flux,
    749 lines) is not comparable with the V2 regions (95% of the V1 lines
    in the forest fall in V2 regions, 11% on the red side, where V1 also
    picks tellurics and noise).

At the end of the phase: measure usage and review time to re-estimate the
rest.

Findings are collected in the register below (section 3.2) and reviewed
with Guido at the end of the phase.

Numerical snapshots (item 1): infrastructure in `tests/conftest.py`
(`snapshot` fixture, reference files in `tests/snapshots/`, updated only with
`pytest --update-snapshots`). Done: pure functions of
`spectrum_operations.py` (rebin, smooth, running mean/std, axis conversions,
flux scaling) on a synthetic spectrum; continuum (absorbed regions, auto
fudge, continuum fit, power law) on the raw J0103-1305 spectrum, 570-700 nm,
z_em = 4.072, recipe defaults; line identification on the same window, with
the continuum and mask saved in the session as frozen inputs; Voigt model and
fit on a synthetic spectrum (resolution from the attribute, from a column,
per component, none) and on J0103-1305 (frozen continuum and systems,
R = 45000 given by the maintainer); the SpectrumV2 methods that combine
spectra (stitch, coadd, resample_on_grid, get_resampled_column,
equalize_to_reference, split, scale_flux) on synthetic exposures with
different grids, noise, flux scales and resolutions; the full chain
(`tests/test_snapshots_full_chain.py`: Set properties, Auto-estimate
continuum, Auto-detect doublets, Auto-fit Ly-alpha forest, recipe
defaults, each stage on the output of the previous one) on J0103-1305,
500-700 nm (Lyman-beta emission included), z_em = 4.072, R = 45000, and on
a mock at z_em = 2.5 (355-460 nm, Lyman-alpha + Lyman-beta forest, CIV and
SiIV systems). The J0103-1305 chain takes ~15 minutes: it is marked `slow`
(`pytest --run-slow`) and runs in a separate CI job (after each merge to
master, weekly and on demand, not on pull requests); the mock chain (~40 s)
runs in the default suite. Both run in a child process with a fixed hash
seed and seeded UUIDs, since the automatic pipelines are not reproducible
otherwise (F19). Only the continuum is compared exactly: the fitted stages
amplify floating-point differences (other CPU, or a harmless change in the
order of operations), so they are compared statistically (component counts
per series, fraction matched within 5 km/s, median and 90th percentile of
the logN and b differences, uncertainties, model pixels), with thresholds
calibrated on the observed variation; a resolving power wrong by 11% on the
mock fails them.

### Phase 2 – Cleanup and packaging
Obsolete and backup files; `pyproject.toml` (placeholders, version
constraints, entry point, optional `[bayes]` extra); **LICENSE (Guido's
decision)**.

### Phase 3 – Decoupling from legacy
One dependency at a time (`Defaults`, `trans_parse`, `xem_d`, … become V2
modules). Legacy code remains only behind the V1 migration loader.
The V1 numerical routines used by the parity tests (phase 1, item 11) are
not deleted: they are moved into a test-only reference module (e.g.
`tests/reference/`), so that V2 can still be checked against the validated
V1 behaviour. Only the code needed by the wx GUI and the rest of V1 is
removed.

### Phase 4 – Splitting large files
`SpectrumV2` into thematic modules (continuum, lines, flux) while keeping
the same public API through delegation; same approach for the GUI.

### Phase 5 – Error handling
Broad `except` clauses, case by case: specific exception, explicit log, or
an item to discuss.

### Releases
So far, releases have been distributed as a macOS app through the GitHub
releases page (not the Mac App Store); Astrocook is not on PyPI yet.
Plan (proposed, each release still Guido's decision):
- **2.0.2 (corrective), right after phase 1:** the V2 regressions (F1, F3,
  F10) and the urgent issues #465 and #480, together with the fixes
  already on master (ASCII export #476/#481, Voigt fitter context, warning
  on saving to V1 `.acs`). Release notes state which numerical results
  change. A 2.0.2 with only the fixes already on master can be brought
  forward if users need the ASCII export earlier.
- **2.1 (first consolidated release), after phases 2-3:** clean packaging
  (LICENSE, entry point, realistic dependency versions); first
  publication on PyPI. Also the macOS app, as before.

### 3.1 User-reported issues
GitHub issues are the single source of truth; this section only says *when*
each open issue is addressed. Rule of thumb: bugs before refactoring,
features after, except urgent features in stable areas.

| When | Issues | Notes |
|---|---|---|
| Phase 1 | [#469](https://github.com/DAS-OATs/astrocook/issues/469), [#477](https://github.com/DAS-OATs/astrocook/issues/477), [#455](https://github.com/DAS-OATs/astrocook/issues/455), [#450](https://github.com/DAS-OATs/astrocook/issues/450), [#436](https://github.com/DAS-OATs/astrocook/issues/436) | "possibly fixed": verify with tests, close or reopen |
| Right after phase 1, before phase 2 | F1 ([#490](https://github.com/DAS-OATs/astrocook/issues/490)), F3, F10 (section 3.2) — done in phase 1-bis (F3 accepted, F10 minimal fix) | V2 regressions (V1 is correct): rebinned errors, auto fudge, convolution with a resolution column. Reference behaviour: V1; each fix still agreed with Guido, with a regression test |
| Right after phase 1, before phase 2 | [#465](https://github.com/DAS-OATs/astrocook/issues/465) | continuum sensitive to spikes (urgent); numerical choice, to be agreed with Guido before fixing |
| Phase 0 (brought forward) | [#476](https://github.com/DAS-OATs/astrocook/issues/476), [#481](https://github.com/DAS-OATs/astrocook/issues/481) | ASCII export; same problem (missing `import re`, found by Ruff), fixed with a regression test |
| Right after phase 1, before phase 2 | [#442](https://github.com/DAS-OATs/astrocook/issues/442) | scatter plot not auto-triggering on region overlap |
| Right after phase 1, before phase 2 | [#480](https://github.com/DAS-OATs/astrocook/issues/480) | import of line lists as ASCII (urgent); touches `io/`, little affected by refactoring |
| With the natural-language log work | [#424](https://github.com/DAS-OATs/astrocook/issues/424) | logging of manual continuum edits |
| After phase 1-bis | [#526](https://github.com/DAS-OATs/astrocook/issues/526) | resolution as a single column and exact convolution (F8, F22, rest of F10, F9); design and decisions in the issue |
| After phase 1-bis | [#517](https://github.com/DAS-OATs/astrocook/issues/517), [#518](https://github.com/DAS-OATs/astrocook/issues/518) | from the findings review: flag and refit components at the parameter bounds, with F21 (`lya_auto` output); editable z window in the manual fits (F4) |
| After phase 3/4 | [#471](https://github.com/DAS-OATs/astrocook/issues/471) | combining order-separated spectra; touches `SpectrumV2` (coadd/stitch) |
| After phase 4 | [#478](https://github.com/DAS-OATs/astrocook/issues/478), [#479](https://github.com/DAS-OATs/astrocook/issues/479) | GUI features; wait for the GUI split |
| After phase 4 | [#446](https://github.com/DAS-OATs/astrocook/issues/446), [#437](https://github.com/DAS-OATs/astrocook/issues/437) | flux-calibration filters; independent, not urgent |

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
| F1 | `rebin_spectrum` | Errors do not decrease with the number of pixels per output bin (σ instead of σ/√N) | numerical | `test_rebin_weighted_average_error`, `test_rebin_errors` (V1 parity; regression) | fixed (phase 1-bis): error of the weighted mean as in V1, sqrt(sum w^2 sigma^2) / sum w; issue [#490](https://github.com/DAS-OATs/astrocook/issues/490) |
| F2 | `rebin_spectrum` | Sigma clipping uses the first-pass bin mean, which includes the outlier: with few pixels per bin (e.g. 3) a single spike can make the whole bin NaN. Clipping runs once, docstring says "iterative". Not a pure V2 regression: V1 clips with `sigma_clip` within each bin, which keeps the bin but cannot remove an isolated spike with 3-5 pixels per bin; neither version removes it. Related to #465 | numerical | snapshot `rebin_spectrum_nm_clip`; `test_rebin_clipping_keeps_bins_with_a_spike`, `test_rebin_clipping_keeps_narrow_absorption` (regression) | fixed (phase 1-bis, agreed with Guido): single-pass clipping against a running median of the input pixels, upward only (a two-sided test removed the cores of narrow lines: equivalent width -20% for a 2.5-pixel line; the previous V2 clipping up to -100%), a bin is never emptied; V1 behaviour dropped |
| F3 | `fit_continuum` | `fudge='auto'` (recipe default) is always exactly 1.0: computed on residuals that are zero by construction, before the smoothing. Related to #465 | numerical | snapshots `continuum_j0103_*`; `test_v2_auto_fudge_is_one` (characterization) | accepted (phase 1-bis, decided with Guido): not a regression in substance. The V1 fudge corrects the running mean V1 builds the continuum from; V2 interpolates the unabsorbed pixels, so the fudge is 1 by construction. Local fudges after the smoothing were tried on J0103-1305: none reproduces V1 (effect ±1% or more vs ±0.2%, no correlation), they would correct the bias of the final smoothing instead. Documented in `fit_continuum`; the smoothing bias goes with #465, F17, F30 |
| F4 | `VoigtFitterV2` | z is bounded to ±`z_window_kms` (default 20 km/s) around the initial guess; a guess further away cannot reach the true value. Intended? | behaviour | `test_fitting.py` (guess inside the window) | open |
| F5 | `convert_y_axis` | Annotation declares 2 return values, the function returns 3 | code quality | snapshot `convert_y_axis` | fixed (phase 1-bis): annotation and docstring list the three returned values |
| F6 | architecture | `core/spectrum_operations.py` imports `gui.debug_utils`: core depends on the GUI | architecture | — | phase 3/4 |
| F7 | GUI | `gui/main_window.py` imports `matplotlib.pyplot` before `pyside_plot.py` selects QtAgg: importing the GUI without a `QApplication` (and no display) fails | architecture | workaround in `tests/test_gui_save_session.py` | phase 4 |
| F8 | resolution handling | The resolution lives in five places (per-component value, `resol` column, `SpectrumDataV2.resol`, meta `resol`/`RESOL`, meta `resol_fwhm`), each read with a different priority by the fitter, the `.acs2` loader, coadd and the GUI. The `.acs2` save does not write the attribute: it survives only through the meta keys (set by "Set properties"). The fitter stamps the resolution used into every fitted component, which then overrides later changes. R vs km/s is guessed from a threshold of 500. The J0103-1305 archive (older version) stores no resolution at all. To be analysed by Guido | behaviour, possible data loss | — (map in the PR that added this row) | design agreed, deferred after phase 1-bis: [#526](https://github.com/DAS-OATs/astrocook/issues/526) |
| F9 | Voigt model | Components narrower than a pixel (e.g. Ly-alpha with b = 0.5 km/s, logN = 17.8 in J0103-1305) make the model depend on the sampling; with R = 45000, 1% of the 570-700 nm window still differs from the model saved in the session. Related to the oversampling check (item 9): the V2 model, sampled at pixel centres, overestimates the equivalent width of CIV by 4% at b = 1 km/s and 0.3% at 3 km/s (2.5 km/s pixels, R = 45000) | numerical | — | with the exact convolution, [#526](https://github.com/DAS-OATs/astrocook/issues/526) (timing to be checked) |
| F10 | `convolve_flux` | With a per-pixel resolution column, pixels are grouped by resolution rounded to 100 (R) or 0.5 km/s, but selected by the raw value: pixels whose resolution is not already rounded are never convolved and the model is 0 there (a flat column at R = 45050 gives a zero model everywhere). Affects rebinned and coadded spectra (coadding three synthetic exposures at R = 42000-47500 gives a non-round column on 90% of the pixels and a zero Voigt model on 90% of them, `test_voigt_model_on_coadded_spectrum_is_not_zero`), RESVEL input, non-round "Set properties" values. In addition, even for a constant R the kernel width is computed at the median wavelength and with the median pixel size, so it is wrong far from the centre of the range (+53% / -26% at the ends of 380-780 nm). Proposed: a single exact convolution in velocity space, after the safety net is complete | numerical, wrong results | `test_model_column_resolution_is_not_zero`, `test_voigt_model_on_coadded_spectrum_is_not_zero`, `test_voigt_model_with_resolution_column` (regression) | selection bug fixed (phase 1-bis): pixels selected by the rounded resolution, NaN filled with the unconvolved flux. Kernel width far from the median wavelength and exact convolution: [#526](https://github.com/DAS-OATs/astrocook/issues/526), after phase 1-bis |
| F11 | grouping | `VoigtModelConstraintV2.set_active_components` ignores `group_depth`: it always calls `get_connected_group` with the default depth 2, so the recipe option (1 = neighbours only) and the `group_depth=0` used by the refit have no effect | behaviour | `test_set_active_components_honours_group_depth`, `test_fit_component_group_depth_0_freezes_neighbours` (regression) | fixed (phase 1-bis, decided with Guido): `group_depth` honoured, default 2 everywhere; the internal fits of `refit_all`/`metals_auto` keep the depth-2 expansion they always had (`FIT_EXPANSION_DEPTH`), since depth 0 there was 2.75x faster but fitted blends worse (J0103-1305: chi2 +14% in the forest); the user can reduce the depth in "Fit component" and in the Bayesian fit |
| F12 | global state | Module-level dictionaries are modified at run time: the GUI rewrites the parameter docs of `EDIT_RECIPES_SCHEMAS["equalize_and_stitch"]`, and `add_component` registers custom series in `STANDARD_MULTIPLETS`. The changes persist for the whole process (and leak between tests) | architecture | `restore_multiplets` fixture in `test_issue_455_multiplet_subset.py` | open |
| F13 | inspector | `SystemTableModel.set_highlighted_uuids` returns early when the selection is unchanged, so link colours are not recomputed when a link is created with the group still selected (model-level reproduction; full GUI not driven) | GUI | `test_link_created_with_group_selected_is_coloured` (strict xfail) | open |
| F14 | grouping | Two components are considered blended when their line centres are closer than max(4 b, 30 km/s), independently of the instrumental resolution and of the line strength. At X-Shooter resolution (R ~ 8900) AlII 1670 and CIV 1550 still overlap at 60 km/s but are not grouped; worse in the NIR or for saturated lines. Possibly why #469 is still seen on X-Shooter data | behaviour | `test_blend_threshold_ignores_resolution` (characterization) | open |
| F15 | `equalize_and_stitch` | The scale factor is a median of the ratio of the two arms in the overlap: with narrow absorption (e.g. tellurics) seen at two different resolutions, it is biased by ~1-1.5% in a mock with equivalent-width-conserving tellurics around 1000 nm (no bias without tellurics). Possibly related to #477; to be checked on real X-Shooter VIS/NIR data | numerical | — (mock in the PR that added this row) | open |
| F16 | smoothing | For the same "std" parameter, V1 `_gauss_convolve` uses the kernel exp(-(dv/std)^2), i.e. sigma = std/sqrt(2), while V2 `smooth_spectrum` uses sigma = std: V2 smooths sqrt(2) wider. Affects the continuum smoothing (default 500 km/s: effective sigma 354 km/s in V1, 500 km/s in V2) and any other use of the parameter. Which convention to keep is Guido's decision | numerical | `test_smoothing_width_differs_by_sqrt2` (characterization) | open |
| F17 | continuum (forest) | On the J0103-1305 window, with the same parameters and matched smoothing, V1 `clip_flux` and V2 `find_absorbed` + `fit_continuum` agree on the red side (masks 99.98%, continuum within 0.4%), but in the Lyman-alpha forest V2 marks more pixels as unabsorbed (57.6% vs 54.8% of the window, masks agree on 92.7% of the forest pixels) and its continuum is lower by 1.2% (median; 5th percentile -4.6%). Not explained by the fudge (F3). Relevant to #465 and to the known low forest continuum | numerical | `test_continuum_forest` (strict xfail) | open |
| F18 | dependencies | GitHub Dependabot reports 36 vulnerabilities on the default branch (2 critical, 8 high, 13 moderate, 13 low), seen on push in October 2026; not yet analysed | security | — (Dependabot alerts) | phase 2 |
| F19 | `refit_all` (reproducibility) | With `group_depth >= 1` the groups are built by taking the next seed from a set of component UUIDs (`next(iter(todo_uuids))`), and the depth-limited friends-of-friends group depends on the seed. UUIDs are random (`uuid4`) and string hashes are randomized per process, so the same input with the same parameters gives a different result at every run. Full chain on J0103-1305, 500-700 nm, two runs: `doublets_auto` gives the same 51 CIV components, but 4 differ by more than 0.1 dex in logN (up to 2.35 dex: logN 10.0 vs 12.35 at z = 3.2836); `lya_auto` gives 477 vs 478 components, ~20 not matched within 5 km/s; the normalised model differs by more than 0.01 on 1.8% of the pixels (up to 0.83). On the mock the differences are ~1e-6 dex. Even with a fixed hash seed and seeded UUIDs the result depends on the machine: on a GitHub runner the J0103-1305 chain gives 477 instead of 478 components after `lya_auto` and a normalised model differing by up to 0.84, with a continuum identical to 1e-15; the mock differs by up to 1e-3 dex in logN. Affects `doublets_auto`, `lya_auto`, `components_auto`, `forest_auto`, "Refit all" from the GUI (schema default `group_depth = 2`) and the refit after importing systems (`recipes/edit.py`) | numerical, reproducibility | `test_refit_all_does_not_depend_on_the_uuids` (regression); `tests/full_chain.py` | fixed (phase 1-bis): seeds taken in order of increasing z; two J0103-1305 full chains with different hash seeds and UUIDs now give identical results. Floating-point differences between machines remain (statistical comparison in the full-chain snapshots) |
| F20 | `VoigtFitterV2` (uncertainties) | Errors are `sqrt(diag(pinv(J^T J)))`: a negative diagonal element (numerically ill-conditioned fit) gives a NaN uncertainty with only a `RuntimeWarning`, no log message and no flag on the component. 1 of 478 components in the J0103-1305 full chain. Related to the chi-squared scaling check (item 9) | numerical | snapshot `full_chain_j0103_forest`; `test_nan_uncertainties_are_reported` (regression) | fixed (phase 1-bis): a warning names the components with NaN errors, numpy RuntimeWarning silenced; the NaN error itself is the flag (no new field) |
| F21 | `lya_auto` (output quality) | On the J0103-1305 full chain (500-700 nm): 39 Ly-alpha components end with b < 10 km/s although `min_b = 10` (the filter runs before the final refit, which moves the remaining components again); components stuck at the parameter bounds are kept (12 at logN = 10, 12 at b = 1 km/s, Ly-alpha and CIV), including two Ly-alpha with logN = 20.2 and 21.1 at b = 1 km/s (z = 3.708, 3.506), unphysical in the forest | behaviour | snapshot `full_chain_j0103_forest` | open |
| F22 | resolution handling | `apply_expression`, `mask_expression`, `smooth_column` and `split` rebuild `SpectrumDataV2` field by field and drop the `resol` attribute (reset to 0); `smooth`, `rebin`, `stitch` and `set_properties` keep it. The Voigt model changes (by 0.05 in flux in a synthetic test at R = 45000) only when the resolution is in the attribute alone: `set_properties`, the `.acs2` loader and the mocks also store it in the metadata, which is kept, so the usual paths are not affected. One more case for F8 | behaviour | `test_some_operations_drop_the_resolution_attribute` (characterization) | deferred after phase 1-bis: [#526](https://github.com/DAS-OATs/astrocook/issues/526) |
| F23 | `flux.rebin` recipe | A non-numeric parameter (e.g. `dx = "abc"`) is logged as a parsing error but the recipe goes on and raises `UnboundLocalError`, instead of returning 0 like the other recipes. In the GUI the worker catches the exception, so it is not a crash | error handling | `test_rebin_bad_input_returns_0` (regression) | fixed (phase 1-bis): the recipe returns 0 on a parsing error |
| F24 | `fit_bayesian` (MCMC) | `run_mcmc` discards a fixed burn-in of 200 steps, and the recipe passes only `nsteps` (GUI parameter, default 500): with `nsteps <= 200` the chain is empty, the medians are NaN and they overwrite z, logN and b of the fitted components, and the model becomes NaN. No warning. Reported at once (loses data; undo restores the previous state) | data loss | `test_fit_bayesian_few_steps_keeps_parameters` (regression) | fixed (phase 1-bis): burn-in defaults to min(200, nsteps // 2), an empty chain is refused, non-finite results never overwrite the components |
| F25 | `delete_component` | The recipe recomputes the model after deleting, which needs a continuum: on a spectrum in physical units without one (and with a median flux outside 0.1-10, otherwise it is taken as normalised), `VoigtFitterV2` raises, the error is logged and the session is returned unchanged. The component is not deleted (inspector, "Clean negligible": in the GUI only "Process completed with warnings"). Seen on the small J0103-1305 archive (764 components, no continuum). Fitting a component instead triggers the continuum estimate first | behaviour | `test_delete_component_without_continuum` (regression) | fixed (phase 1-bis): without a continuum the components are deleted, the model is not recomputed and the stale `model` column is removed, with a warning |
| F26 | `.acs2` session log | The `.acs2` archive does not keep the V2 session log: `_serialize_v2_metadata` writes `dataclasses.asdict(HistoryLogV2)`, which is `{}` because the entries are set in `__init__` rather than as dataclass fields (the same cause found in phase 0 for `.acs`), and the loader does not read a log back anyway. Checked through the GUI: after "Set properties" and "Smooth" the archive has `log_history_json: {}` and the reloaded session has an empty log. The GUI tells the user to save as `.acs2` to keep the log. Reported at once (loses data) | data loss | `test_session_log_is_saved`, `test_loaded_log_is_kept_by_the_gui` (regression) | fixed (phase 1-bis): the log entries up to the current one are saved and restored on load, also in the GUI |
| F27 | `.acs2` identifications | `region_identifications` is saved as a JSON string and reloaded as a dict with integer keys. The identification viewer calls `json.loads` on it and shows "Error parsing metadata" (the plot tooltips handle both forms) | GUI | `test_identification_viewer_after_round_trip` (regression) | fixed (phase 1-bis): the viewer accepts both the JSON string and the dict |
| F28 | `.acs2` Bayesian samples | The posterior samples stored by `fit_bayesian` in the system list metadata (`bayesian_results`, used by the corner plot) are dropped on save without a warning: the FITS header cannot hold them. The fitted values and errors are kept | data loss (recomputable) | `test_bayesian_results_are_saved` (strict xfail) | open |
| F29 | `VoigtFitterV2` (uncertainties) | The errors are sqrt(diag((J^T J)^-1)) with residuals weighted by 1/dy, not scaled by the reduced chi2, although the comment in `fit` says they are. V1 scaled them (lmfit `Model.fit`, `scale_covar=True` by default): a V1/V2 difference. Monte Carlo on mocks: unscaled errors are right when the error spectrum is right, and wrong by the same factor when it is not (e.g. 2 times too small with the error spectrum underestimated by 2, reduced chi2 3.8); scaled errors are right in all cases. Scaling, or reporting the reduced chi2 next to the errors, is Guido's decision | numerical, V1/V2 difference | `tests/test_error_scaling.py` (characterization) | open, decision |
| F30 | `running_mean`, `running_std` | Both use `np.convolve(mode='same')`, which pads with zeros: the first and last h values are biased towards 0 (4/7 of the true mean at the edge for h = 3). V1 `running_mean` repeats the first and last full-window value and ignores NaNs. Used by `find_absorbed_regions` on each interval, with half-windows of hundreds of pixels (5000 km/s / 8 on J0103-1305). Replacing it with the V1 function raises the agreement of the forest mask with V1 from 92.7% to 95.7% and the 5th percentile of the continuum ratio from -4.7% to -3.0% (median unchanged, -1.2%): a partial cause of F17 | numerical, V1/V2 difference | `test_running_mean_and_std_pad_the_edges_with_zeros` (characterization) | open |
| F31 | dead code | `compute_identification_signal` and `find_signal_peaks` (`core/spectrum_operations.py`) are not called anywhere: the V2 identification goes through `SpectrumV2.identify_lines`. Section 2 of this roadmap describes the identification as based on the former. Found by mutation testing (no test reaches them) | code quality | — | phase 2 |
| F32 | `VoigtFitterV2` (broad lines) | `_residual_function` appends a penalty 5 (b - 50) for every free b > 50 km/s to the residuals ("soft penalty to prevent runaways"). Broad lines are biased: on mocks (Ly-alpha, logN 13.8, R = 45000), at S/N 20 a line with b = 60, 80 or 120 km/s is fitted with b = 50 km/s and logN low by 0.04, 0.11 and 0.35 dex; even at S/N 1000, b = 120 comes out as 72. The penalty also enters chi2 and the uncertainties. Affects OVI, broad Ly-alpha (BLAs), blends fitted with one component | numerical, wrong results | `test_fit_broad_line`, `test_b_max_bounds_the_fit`, `test_pipelines_bound_b_and_manual_fits_do_not` (regression) | fixed (phase 1-bis): penalty removed; hard bound b <= 50 km/s (`max_b`) in the automatic pipelines, 300 km/s in the manual fits (without any bound `lya_auto` on J0103-1305 takes 4 times longer and gives spurious broad lines); components at the bound: [#517](https://github.com/DAS-OATs/astrocook/issues/517) |

### 3.3 Phase 1-bis: findings and urgent issues
Review of the findings register with Guido (October 2026): the proposals
below are agreed in general; each item is still a separate small PR, one at
a time, with the related strict xfail turned into a passing test (or a new
regression test), and its numerical effect stated in the PR. Items marked
"decision" get a short proposal first, agreed before the code changes.

Order:
1. F24: `fit_bayesian` with `nsteps <= 200` overwrites parameters with NaN
   (burn-in proportional to `nsteps`, or refuse with a message).
2. F26: `.acs2` does not keep the session log (serialize the log entries,
   read them back on load).
3. F32: penalty on b > 50 km/s biases broad lines (remove it; keep a wide
   upper bound on b; decision on the bound).
4. F25: components cannot be deleted without a continuum.
5. F19: automatic pipelines not reproducible (deterministic order of the
   groups in `refit_all`).
6. F1 (#490), F3: V2 regressions in rebinned errors and auto fudge (V1 is
   the reference); F2 reviewed together with F1.
7. F10 with F8 and F22: rewrite of the convolution (exact, in velocity
   space) and a single source of truth for the resolution (decision on the
   design first).
8. F27, F20, F23, F11, F5: small fixes (identification viewer after a
   round trip, warning and flag on NaN uncertainties, `flux.rebin` error
   path, `group_depth`, annotation).
9. Decisions, one at a time: F29 (scale the uncertainties by
   sqrt(reduced chi2), as V1, and report it; proposed), F16 (V1 smoothing
   convention; proposed), F30 then F17 (V1 running mean at the edges;
   proposed), F4, F9, F14, F21, F15 (needs X-Shooter data).
10. Issues of section 3.1 scheduled "right after phase 1": #465, #480,
    #442.
Then release 2.0.2 (Guido's decision), with at least F24, F26 and F32 in
addition to the fixes planned in "Releases". Technical debt goes to its
phase: F18, F31 (phase 2); F6, F7 (phases 3/4); F12, F13 (phase 4); F28
(save the Bayesian samples as a separate file in the archive) when needed.

### Rough estimates
10–20 work sessions, 2–5 calendar weeks; the bottleneck is human review and
the availability of test data.

### Decisions reserved to Guido
- License.
- Questionable behaviours found by characterization tests (bug or intended?).
- Approach to the continuum spike sensitivity (#465).
- Scientific validation of results.
- Releases (see "Releases" above): timing of 2.0.2 and 2.1, PyPI
  publication, macOS app signing and notarization.

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

Proposed pilot: MCP server with 3–4 recipes, after phase 1.

---

## 5. v3 vision (horizon, after v2 consolidation)

**Goal:** given a spectrum, produce in less than one hour **a collection of
complete models**, with uncertainties and a reliability map, which also
report the possible alternative identifications of the same feature, each
with a parameter quantifying its credibility.

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
Pilot on metals outside the forest, then extension to the Lyα forest. The
LLM agent stays at the margins (strategy in stuck cases, guidance during
review); each of its decisions goes into the log.

### Implications for the v2 work
- The mock generator from phase 1 is the first core of the v3 test bench.
- When refactoring, do not introduce the assumption "one model per session".
- The log must record the fit decisions.

---

## 6. Progress

| Phase | Status |
|---|---|
| Initial analysis and roadmap | done (October 2026) |
| 0 – Environment | done (October 2026) |
| 1 – Safety net | done (October 2026; log replay, V1 `.acs` and loader tests wait for data) |
| 1-bis – Findings and urgent issues (section 3.3) | in progress: F24, F26, F32, F25, F19, F1, F2 done, F3 accepted, F10 minimal fix (rest: #526), F27, F20, F23, F11 done |
| 2 – Cleanup and packaging | not started |
| 3 – Decoupling from legacy | not started |
| 4 – Splitting large files | not started |
| 5 – Error handling | not started |
