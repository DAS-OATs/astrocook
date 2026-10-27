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
   emcee/dynesty/corner are required instead of optional (made optional
   in phase 2-bis: extra `[bayes]`).
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
   calls (dead code, F31, removed in phase 2). Found F30.
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
   warnings are cleaned up in a dedicated PR (done in phase 2-bis: CI now
   runs `ruff check .`).
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
      ignored the resolution (F14, fixed in phase 4).
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
Done (October 2026), PRs #543-#554:
- obsolete files removed after checking that nothing uses them:
  `legacy_home/`, backup data files shipped in every package (~54 MB),
  the Gemfile of the V1 docs (`docs_backup/`, still deployed under
  `/v1/`, stays); the `legacy/*_old.py` modules stay for phase 3;
- `pyproject.toml`: placeholders removed, unused dependencies
  (`statsmodels`, `specutils`) dropped, `numpy>=1.26`, `astrocook`
  command (launcher moved into the package, `python -m astrocook` works
  from any directory), extras `[docs]` and `[build]`;
- GUI images inside the package (`astrocook/assets/`);
- macOS build: relative paths in `Astrocook.spec`, single PyInstaller run
  in `build_app.sh`; build and docs scripts check their tools;
- docs: installation with pip from GitHub; the website is published by
  `.github/workflows/docs.yml` after each merge (`build_site.sh`,
  `deploy_docs.sh` kept as a manual fallback);
- F18 (all Dependabot alerts were on the V1 Gemfile), F31 (dead
  identification code), F35 (archives unpacked without a path check)
  fixed; F34 found (phase 2-bis).

Left open:
- **LICENSE**: GPLv3 or later agreed with Guido; PR #549 waits as a draft
  for INAF (copyright holder) and for the contributors. Needed before the
  PyPI publication (independent of the phases: see "Releases").
- Optional `[bayes]` extra: agreed, moved to phase 2-bis.
- Docker containers (`containers/`, broken since the V1 launchers moved):
  decided with Guido (October 2026): no renewal for now; `containers/` is
  removed in phase 3; a headless image comes with the full MCP server
  (section 4); an image with the GUI only if users ask for it.
- `telluric.fits` and `sky_telluric.fits` (~62 MB each, used only by
  `legacy/`) keep the wheel near PyPI's 100 MB limit: to be decided in
  phase 3. Decided: removed from the package (see phase 3).

### Phase 2-bis – Leftovers of phase 2
Agreed with Guido (October 2026), before phase 3; one small PR each, in
this order:
1. **F34**: main window icon from the existing `logo_3d_HR.png` (file to
   be confirmed by Guido in the PR).
2. **Ruff warnings**: clean up the ~80 pyflakes warnings (unused imports
   and variables, f-strings without placeholders, redefinitions), each
   unused variable checked by hand since it may hide a bug; then CI runs
   `ruff check .` with all the configured rules instead of F821 only.
   Done: 162 warnings with Ruff 0.17 (116 unused imports, 24 unused
   variables, 16 f-strings, 6 redefinitions); found F36 and F37 (left
   unchanged, `# noqa` with a reference; deferred to #559 and #560); the unused `aic_penalty` of
   `metals_auto` is already documented as such in its docstring.
3. **Sphinx warnings**: the 4 warnings of the docs build (duplicate object
   descriptions, `api/astrocook.mocks` not in a toctree).
   Done: 0 warnings. Besides the 4, the `__all__` added in the Ruff PR
   made autodoc describe the re-exports of `astrocook.core` twice
   (ambiguous cross-references); `docs/conf.py` now ignores `__all__`
   (`ignore-module-all`), `HistoryLogV2` documents its attributes with
   `#:` comments, `api/astrocook.mocks` is in the API reference.
4. **`[bayes]` extra**: emcee, dynesty and corner become optional; the
   Bayesian fit raises a clear error naming `pip install
   "astrocook[bayes]"`, the GUI disables its entries with a tooltip; the
   macOS app and `[dev]` keep them.
   Done: `fit_bayesian` raises an `ImportError` naming the extra (instead
   of returning 0) before sampling; the system inspector imports the
   corner plot only when it opens it (before, the GUI did not start
   without corner); `[dev]` and `[build]` include `[bayes]`.

Not in phase 2-bis: LICENSE (#549, draft until INAF answers; needed only
for PyPI, 2.1.0), Docker renewal (to be discussed with Guido), size of
the telluric files (phase 3).

### Phase 3 – Decoupling from legacy
Done (October 2026), PRs #565-#584; no V2 module depends on the V1 code
and `astrocook/legacy/` is removed:
- guard test against new imports from `legacy/` (removed with `legacy/`);
  `msg_param_fail` in `recipes/messages.py` (centralizing the messages:
  #566, phase 6);
- `trans_parse` in `core/atomic_data.py`, `x_convert` replaced by
  `convert_axis_velocity`, `qso_composite` in `core/templates.py`; the
  unused V1 `Defaults` and `GUILog` removed from the session and the GUI;
- the V1 routines of the parity tests are verbatim copies in
  `tests/reference/v1/`, with their own data files;
- saving in the V1 `.acs` format removed (Guido's decision); V1 archives are
  read by a native reader, `io/v1_archive.py`, checked against the former
  reader on Guido's archives and on two small anonymised test archives
  (`tests/data/v1_*.acs`); F42, F43, F44 fixed along the way;
- logging is left to the application (Guido's decision): importing the
  package no longer configures it; new tutorial *Using Astrocook from
  Python*;
- `legacy/` (~22,500 lines), `containers/`, the `[legacy]` extra, `lmfit`
  and `cycler` removed (`tqdm` only for the tests); the data files used only
  by V1 removed (see below): the wheel goes from 79.7 to 2.6 MB.

Outside the plan: the documentation website is published per version
(#569, #570): the root holds the docs of the latest release, `/dev/` those
of master. Found in phase 3: F38-F44 (section 3.2); issues opened: #566,
#568, #573, #575 (section 3.1).

The original plan:
One dependency at a time (`Defaults`, `trans_parse`, `xem_d`, … become V2
modules). Legacy code remains only behind the V1 migration loader.
The V1 numerical routines used by the parity tests (phase 1, item 11) are
not deleted: they are moved into a test-only reference module (e.g.
`tests/reference/`), so that V2 can still be checked against the validated
V1 behaviour. Only the code needed by the wx GUI and the rest of V1 is
removed, together with the V1 Docker image (`containers/`, broken: it
installs the V1 dependencies and launches the removed `ac_gui.py`).

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

### Phase 3-bis – Small fixes before 2.1.0
Done (October 2026), PRs #587-#592. Agreed with Guido (October 2026): a
few small fixes before release 2.1.0, one small PR each, with the related
strict xfail turned into a passing test (or a new regression test). Found
along the way and fixed in the same PRs: branching after an undo deleted
the redo entries of the parent log (#575); the docstring of
`rate_doublet_candidate` had the order of the regions reversed (F6); the
link colours changed at random when a link was added (F13);
`evaluate_uncertainties` did not update `dbtur` either (#559). Found and
recorded: F45 (added to the phase), a note on F8 (`rebin` keeps the old
`resol` meta key; #526). The only numerical change is the `dbtur` column
of the fitted components (#559). Order agreed with Guido:
1. [#575](https://github.com/DAS-OATs/astrocook/issues/575) (F41,
   priority): a recipe that opens a new session is also recorded in the
   log of the parent session. Done: the entry goes only into the log of
   each new session, after the copy of the parent log; logging on the
   parent also deleted its redo entries after an undo, fixed with it.
   Logs saved by earlier versions are left as they are: a spurious entry
   in a parent log cannot be told apart from the right one in a branch;
2. F6: `core/spectrum_operations.py` imports `gui.debug_utils`, so the core
   depends on the GUI (importing it also imports PySide6; matters for the
   headless MCP pilot). Proposed: an optional debug hook in `core`, which
   the GUI connects to its signal; no numerical change. Done:
   `core/debug_hooks.py`; importing the core, recipes, fitting, I/O and
   mocks no longer imports PySide6 (checked by a test); the docstring of
   `rate_doublet_candidate` had the order of the two regions reversed;
3. F13: link colours not recomputed in the System Inspector when a link is
   created with the group still selected (strict xfail already in place).
   Done: the colours are recomputed when the links change; each link also
   keeps its colour when another is added (before, the colours followed
   the order of a set of UUIDs and changed at random);
4. [#559](https://github.com/DAS-OATs/astrocook/issues/559) (F36, brought
   forward from phase 4): `dbtur` not updated after a fit; changes only
   the uncertainty of `btur`, slow tests run before the merge. Done:
   `dbtur` stored like the other errors (0 when `btur` is frozen, as it is
   by default; kept for spectators); `evaluate_uncertainties` (#534) did
   not update it either, fixed with it;
5. F45: the `.acs2` archive records the Astrocook version that saved it
   (found with #575; added to the phase by Guido). Done: `astrocook_version`
   in `_meta.json`, read back as `SessionV2.saved_with_version` (`None` for
   the archives saved up to 2.0.3) and kept by the derived states; not in
   the FITS headers (agreed with Guido);
6. [#560](https://github.com/DAS-OATs/astrocook/issues/560) (F37, added
   after discussing it with Guido): `norm` of the `rebin` recipe ignored.
   V1 rebinned y/cont and dy/cont and set `cont` to 1; implementing it in
   V2 would also need `model` (in flux units) and a choice between `cont`
   and `cont × telluric`. Guido's decision: not needed (it served one
   specific case; the plot has a normalized view, and `apply_expression`
   computes `y/cont`), so `norm` is no longer offered by the dialog and
   `norm='True'` logs a warning. Done. An export of the normalized spectrum
   (e.g. to ASCII) may come later.

### Phase 4 – Scientific correctness
First the validation on mocks (below), so that every fix of this phase is
measured on it; then the issues deferred from phase 1-bis, one PR (or a
short series) each, every numerical choice agreed with Guido and every
changed result stated in the PR (also on the benchmark), with the slow
tests run before the merge:
- [#526](https://github.com/DAS-OATs/astrocook/issues/526): resolution as a
  single per-pixel column and exact convolution (F8, F22, rest of F10, F9,
  F14);
- [#517](https://github.com/DAS-OATs/astrocook/issues/517): flag and refit
  components at the parameter bounds (F21);
- [#465](https://github.com/DAS-OATs/astrocook/issues/465): continuum
  sensitivity to spikes, with F30 and F17;
- [#518](https://github.com/DAS-OATs/astrocook/issues/518): editable z window
  in the manual fits (F4);
- [#559](https://github.com/DAS-OATs/astrocook/issues/559): `dbtur` not
  updated after a Voigt fit (F36, found in phase 2-bis); brought forward
  to phase 3-bis;
- [#560](https://github.com/DAS-OATs/astrocook/issues/560): `norm` parameter
  of the `rebin` recipe ignored (F37, found in phase 2-bis); done in phase
  3-bis (no longer offered);
- [#568](https://github.com/DAS-OATs/astrocook/issues/568): a single table of
  atomic parameters, merging `ATOM_DATA` and `atom_par.dat` (F38, F39, found
  in phase 3);
- [#573](https://github.com/DAS-OATs/astrocook/issues/573): recipe defaults
  that differ between the GUI and the API (F40, found in phase 3).

The order of the issues is to be agreed with Guido at the start of the
phase: see "Plan of the phase" below.

**Validation on mocks** (agreed with Guido at the start of the MCP pilot,
October 2026, as the first part of this phase). The
snapshots record what the recipes return, not whether it is right: no
recipe has a measured accuracy today. The validation measures the current
continuum, identification and doublet fit (the recipes of the MCP pilot)
on mocks with known parameters, with no numerical change; the numbers
become the baseline that every later change is measured against
(the fixes of this phase, phase 7, and the "Validation" item of
section 2).
- Mocks: the generator of phase 1, extended where needed with the cases
  the recipes are known to handle badly (DLAs and sub-DLAs, strong
  emission lines, BALs if feasible), over a grid of S/N, resolution and
  z_em.
- Continuum: ratio of estimated to true continuum as a function of
  wavelength, z and S/N, inside and outside the forest; the failure cases.
- Identification: precision and recall per species, as a function of
  logN, b and S/N; the fraction of the absorbed pixels that no
  identification explains.
- Fit (`doublets_auto`): bias and pulls ((fit - true) / error) of z, logN
  and b; components at the bounds (#517).
- A benchmark script (`benchmarks/`, run by hand, and possibly as a slow
  test with loose thresholds) and a report in this section. Weak points
  go to the findings register.

**Plan of the phase** (agreed with Guido at its start, October 2026, after
the MCP pilot: grid and 20 km/s matching tolerance as below, F49, F50 and
F55 added to the phase, #465 brought forward before #526).

*The validation, in three PRs, with no numerical change:*
1. *Mocks.* What the generator has: power-law continuum with Gaussian
   emission lines, Lyman forest (`Ly_a`/`Ly_ab`; series up to `Ly_8` in
   `STANDARD_MULTIPLETS`), metal populations, Voigt profiles with
   damping (so DLAs are already possible), spikes, a resolution per
   pixel. To add: populations of DLAs and sub-DLAs and of metal systems
   with their Lyman lines, the Lyman-limit edge, and an option to leave
   out the truth columns (`cont`, `mock_flux`): in the MCP trial an
   agent fitted on the true continuum carried by the mock. BALs and
   tellurics (F54; telluric data in the pre-release `data-telluric-2026`)
   only if cheap, otherwise phase 7. `mocks/` stays independent of
   `fitting/`.
2. *Benchmark.* `benchmarks/validate_mocks.py`: a grid of mocks, the chain
   of the pilot (`estimate_auto`, `identify_lines`, `doublets_auto`, with
   the GUI defaults of the schemas), and the metrics against the truth,
   written to a JSON file:
   - continuum: ratio to the true one (median, 5th-95th percentiles) in
     the forest and on the red side, in bins of wavelength;
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
3. *Baseline.* The numbers of the current version in this section, the
   weak points in the findings register, and a slow test that runs a
   reduced grid with loose thresholds, to catch regressions.

*Then the fixes, one PR (or a short series) each, every change measured
on the benchmark and stated in the PR:*
4. [#573](https://github.com/DAS-OATs/astrocook/issues/573) (defaults,
   F40): first a test with the 15 differences as strict xfails, then the
   values decided by Guido. The benchmark runs with the schema defaults,
   so its baseline does not move.
5. [#518](https://github.com/DAS-OATs/astrocook/issues/518) (editable z
   window, F4): small, defaults unchanged; its new parameter follows the
   rule of #573.
6. [#568](https://github.com/DAS-OATs/astrocook/issues/568) (atomic data,
   F38, F39): the wrong values first, then the single table. The mocks
   use the same table as the fitter, so its errors cancel in the
   benchmark: it is checked against standard line lists instead.
7. Small fixes of the identification found in the pilot (added to the
   phase by Guido): F55 (candidates limited to z <= z_em plus a proximity
   margin, to be decided), F49 (one candidate per doublet), F50 (the
   pipelines keep the absorbed mask of the user, or say they replace
   it). Measured on purity and on the spurious components.
8. [#465](https://github.com/DAS-OATs/astrocook/issues/465) (continuum:
   spikes, F30, F17; labelled urgent, brought forward by Guido): the
   forest part measured on the benchmark (which continuum is closer to
   the truth, the open question of F30); the spike part on mocks with
   spikes at low resolution, and on real NIR data when they arrive
   (F15).
9. [#526](https://github.com/DAS-OATs/astrocook/issues/526) (resolution
   and exact convolution, F8, F22, F10, F9, F14): the largest change,
   at the core of every fit; four PRs as planned in the issue, after
   Guido's decisions (a)-(c) there; timing benchmark and slow tests.
10. [#517](https://github.com/DAS-OATs/astrocook/issues/517) (components
    at the bounds, F21, with F52): last, because how many components end
    at a bound depends on the resolution (#526), the continuum (#465) and
    the z window (#518).
Then release 2.2.0.

Why this order: the small, independent fixes first (they make the
benchmark measure the right thing: defaults, z window, atomic data,
spurious candidates), then the large changes of the continuum (#465,
urgent, reported by XQR-30) and of the model (#526), and the bounds last,
tuned on the final model.

*Still to decide, when the item comes:* the decisions listed in each
issue (the proximity margin of F55: 3000 km/s, Guido). The detection limit of the recall
stays at 3 sigma (Guido, with the baseline below: the identification
should improve down to it).

**Baseline** (version 2.1.0, October 2026; PRs #602, #603 and the PR that
added this paragraph). The full grid of `benchmarks/validate_mocks.py`
(135 mocks, 26 minutes with 4 workers, deterministic: two runs gave the
same numbers). Its summary is in
`benchmarks/results/validate_mocks_baseline_2.1.0.json` (until 2.2.0
in `validate_mocks_baseline.json`, now the baseline of 2.2.0, see
"Result of the phase"); `--compare` prints
the key numbers of a new run next to it, and a slow test runs the reduced
grid (16 mocks) against loose thresholds, to catch regressions.

*Continuum* (estimated / true, median over the mocks; in brackets the
median 5th percentile):

| | Forest (Ly-b to Ly-a) | Red side |
|---|---|---|
| all | 0.936 (0.820) | 1.002 (0.996) |
| S/N 10 / 30 / 100 | 0.869 / 0.936 / 0.978 | 1.005 / 1.002 / 1.001 |
| z_em 2.5 / 3.5 / 4.5 | 0.970 / 0.935 / 0.883 | 1.002 |
| with a DLA (27 mocks) | 0.920 (lowest 0.657) | |

Excellent on the red side; low in the forest, more at high z and low
S/N (F56).

*Identification of the metal doublets* (`identify_lines`; matched within
20 km/s; completeness over the systems whose weaker line is detected at
3 sigma):

| | Purity | Completeness | Duplicates (F49) |
|---|---|---|---|
| CIV | 0.99 | 0.85 | 505 |
| SiIV | 0.97 | 0.69 | 267 |
| MgII | 0.98 | 0.91 | 345 |

Completeness against the significance of the weaker line: 0 below 3
sigma, 0.05 at 3-5, 0.30 at 5-10, 0.91 above 10 (F57). No candidate
beyond z_em: the mocks have no absorption there, so F55 is not measured
yet.

*Fit* (`doublets_auto`): 1534 of the 2689 true components matched (many
of the others are weak or unresolved, especially at R = 8900), 428
spurious, 413 flagged (351 with absurd errors, F52; 112 with b at a
bound). Bias: none (median pulls of z, logN and b between -0.14 and
-0.07). Width of the pulls (half the 16-84% interval) by kind of
component:

| True component | n | z | logN | b |
|---|---|---|---|---|
| single, one fitted component near it | 53 | 1.05 | 1.27 | 1.04 |
| single, more fitted components near it | 107 | 1.38 | 0.99 | 1.16 |
| in a system of several, one fitted near it | 214 | 4.93 | 6.50 | 4.35 |
| in a system of several, more fitted near it | 1160 | 1.63 | 1.99 | 1.77 |
| all | 1534 | 1.78 | 2.29 | 2.07 |

The errors of the fitter are right for isolated lines; the width comes
from the decomposition of the systems with several components (F58).
The total column of a system (sum of the fitted components within its
velocity range, against the true sum) does not depend on it: fitted for
828 of the 1080 systems, median difference -0.002 dex, spread 0.10 dex,
pull width 1.24; 1.01 at R = 140000, 1.08 at R = 45000, but 2.46 at
R = 8900, with a spread of 0.39 dex (F59). Median reduced chi2 1.07.

*What should improve each number:* the duplicates and the spurious
components, F49 (this phase); the forest continuum, #465 (this phase)
and the emission model (phase 7); the components at the bounds and with
absurd errors, #517 with F52 (this phase); R = 8900, to be checked with
#526 (this phase); the completeness down to 3 sigma and the
decomposition, the identification and the fit driven by the residuals
of phase 7.

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

Where the changes come from: the duplicates, F49 (#611); purity, F49
and F55 (#612); the spurious and flagged components, F49, the cleanup
of the pipelines (#628) and the bounds of #517 (#627-#629); the system
column at R = 8900, the exact convolution of #526 (#617) with F61
(#619), which were the cause of F59 (the residual median at R = 8900
is +0.05 dex, reduced chi2 1.6); the forest continuum, #465 (#614).
The 14 components with absurd errors left (F52) are flagged by the
fitter only when their errors are NaN; the others are recognised only
by the benchmark.

*The pulls of the components are wider, and this is not a
regression of the fitter*: for isolated components they are still ~1.1
(the errors are right), but after F49 the multi-component systems are
fitted with fewer components (one candidate per doublet instead of
several duplicates, which by chance added components), so each fitted
component absorbs several true ones, and its pulls against the nearest
true component grow (F58; the total column of the systems is better,
spread 0.077 dex against 0.10). The decomposition of the systems is the
task of the fit driven by the residuals of phase 7; the system pull
width of 2.27 comes from the same cause (the errors of a fitted
component do not include the components it misses). The median reduced
chi2 grows from 1.07 to 1.25 for the same reason.

*Left for the next phases*: F52 (absurd errors without NaN, partly
covered by the flags), F56 (forest continuum at high z and low S/N:
the emission model of phase 7), F57 (completeness down to 3 sigma), F58
(decomposition), F59 residual (+0.05 dex at R = 8900), the components
that reach `max_logN` = 19 in the forest fits (phase 7). With the release the
baseline file became that of 2.2.0 (this run; Guido's decision): the
later phases are measured against it. That of 2.1.0 is kept in
`benchmarks/results/validate_mocks_baseline_2.1.0.json`.

### Phase 5 – Splitting large files
`SpectrumV2` into thematic modules (continuum, lines, flux) while keeping
the same public API through delegation; same approach for the GUI.

**GUI performance** (agreed with Guido, October 2026, after he found the
System Inspector slow to open): splitting the GUI touches the same code,
so it is measured here. First a timing benchmark of the GUI, like
`benchmarks/bench_fitting.py` for the fit, on J0103-1305 (764
components): opening the System Inspector, redrawing the main plot,
a refit from the inspector; run before and after each change. Measured
in phase 4 (offscreen, `cProfile`): the inspector opens in ~6 s the first
time and ~3 s afterwards, with three causes, each a local change that
does not touch the results:
- the table (~3.7 s): its columns resize to their contents
  (`ResizeToContents`), so Qt queries every row, not only the visible
  ones, several times (138,000 calls of `SystemTableModel.data` for 764
  rows);
- the main plot (~1.2 s): opening the inspector highlights the selection
  and redraws the whole spectrum, rebuilding the x axis as an astropy
  quantity in a loop (970 times, 0.84 s);
- matplotlib: text layout, the best position of the legends
  (`loc="best"`) and the ticks, ~1-3 s between the main and the velocity
  plot.

Order of the phase (agreed with Guido, October 2026): the GUI benchmark
first; then the three local fixes of the System Inspector above, one PR
each, measured before and after; then the split of `SpectrumV2`, covered
by the snapshots; then the split of the GUI (`main_window.py`,
`pyside_plot.py`, `system_inspector.py`), with F7. The phase is pure
refactoring and performance: the snapshots and the mock benchmark stay
unchanged; the visible changes of the GUI go to phase 5-bis. F12 (global
state modified at run time) is moved to phase 6-bis, where the recipe
registry is redesigned.

**GUI benchmark** (`benchmarks/bench_gui.py`, run by hand like
`bench_fitting.py`): the GUI offscreen on the small J0103-1305 archive of
the tests (2156 pixels, 764 components), four cases, each ended after the
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

**GUI look and ergonomics** (review asked by Guido, October 2026, on
screenshots taken offscreen). Done, ahead of the phase (GUI only, no
numerical change): the session list and the plot controls are floating
cards, as tall as their contents and anchored below the menu bar, that
fade out while the mouse is on the plot; their tabs (a list icon and a
cog) sit just outside the cards' edges. The panels keep floating over the
plot instead of shrinking it, because a narrower plot changes the aspect
of the spectrum (Guido's decision); a card as tall as its contents, the
fading, and the anchoring at the top (the controls do not move under the
mouse when a section opens) reduce what they cover. The sections of the
plot controls are collapsible, with their state remembered; while Plot
Elements is open its checkboxes show the plotted styles and those entries
leave the matplotlib legend, which used to show through the panel. Two
columns were considered and rejected: a wider card covers more of the
plot. The proposals of the review still open are in phase 5-bis.

### Phase 5-bis – GUI look and ergonomics
Agreed with Guido (October 2026): the visible changes of the GUI proposed
by the review of phase 5 ("GUI look and ergonomics"), after the split of
the GUI, so that they land in the new modules and phase 5 stays pure
refactoring; each proposal agreed with Guido one by one, with the
screenshots of the tutorials rebuilt before and after. Proposals:
- a status bar (cursor position, outcome of the last recipe) and
  keyboard navigation in the System Inspector;
- the System Inspector shows empty panels without saying why when the
  session has no continuum;
- typed widgets in the recipe dialogs (spin boxes, column combos, a
  selector of transitions) and short labels with the long text in the
  tooltip: with phase 6-bis;
- one theme module instead of colours written in a dozen stylesheets
  (dark mode on macOS);
- the menus: *Smooth Column* (Edit) and *Smooth Spectrum* (Flux),
  *Import Systems* in Edit instead of Absorbers, the View toggles as
  checkable entries; recent files and the window geometry not saved;
- the main plot: the 1-sigma error is drawn as a line by default although
  a comment in the code says band (to be decided), the HI ticks are a
  saturated red block, the title repeats the window title.

Then release 2.3.0 (see "Releases").

### Phase 6 – Error handling
Broad `except` clauses, case by case: specific exception, explicit log, or
an item to discuss.

### Phase 6-bis – Recipe API
Agreed with Guido at the start of the MCP pilot (October 2026). The
interface of the recipes redesigned once, for the GUI, the scripts and
the agents together, starting from the report of the MCP pilot (section
4.1):
- typed parameters instead of strings, with a single set of defaults
  (F40, #573, if not already done in phase 4);
- a structured result for every recipe (the new session, a summary,
  quality flags, warnings), instead of results left in the session or in
  the metadata (e.g. the identifications as JSON strings);
- exceptions with a message instead of `return 0` and a log line (with
  phase 6);
- names and descriptions reviewed where the pilot showed them unclear;
- the recipe registry (the `*_RECIPES_SCHEMAS`, from which the GUI
  builds its menus and dialogs) as the single description of a recipe
  for the GUI, the scripts and the server, with its result type and
  whether it is offered to agents: adding a recipe means adding it
  there, and nothing else;
- room in the results for what the algorithms do not provide yet (the
  uncertainty of the continuum, the credibility of an identification and
  its alternatives, section 5), so that phase 7 changes the algorithms
  without changing the interface.
- F12 (moved from phase 5, agreed with Guido, October 2026): the
  module-level dictionaries modified at run time (the parameter docs of
  `EDIT_RECIPES_SCHEMAS["equalize_and_stitch"]` rewritten by the GUI, the
  custom series of `add_component` in `STANDARD_MULTIPLETS`), with the
  redesign of the registry.
The full MCP server (section 4) is generated from this API.

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
- **2.3.0, after phases 5 and 5-bis (planned, agreed with Guido, October
  2026):** the faster System Inspector, the split of the large files, the
  GUI review of phase 5-bis, F7; no numerical result changes (the
  *Accuracy* page stays as for 2.2.0). A stable point before phases 6 and
  6-bis, which change the interface of the recipes. With the LICENSE and
  PyPI if INAF has answered by then (#549).

### 3.1 User-reported issues
GitHub issues are the single source of truth; this section only says *when*
each open issue is addressed. Rule of thumb: bugs before refactoring,
features after, except urgent features in stable areas.

| When | Issues | Notes |
|---|---|---|
| Phase 1 | [#469](https://github.com/DAS-OATs/astrocook/issues/469), [#477](https://github.com/DAS-OATs/astrocook/issues/477), [#455](https://github.com/DAS-OATs/astrocook/issues/455), [#450](https://github.com/DAS-OATs/astrocook/issues/450), [#436](https://github.com/DAS-OATs/astrocook/issues/436) | "possibly fixed": verify with tests, close or reopen |
| Right after phase 1, before phase 2 | F1 ([#490](https://github.com/DAS-OATs/astrocook/issues/490)), F3, F10 (section 3.2) — done in phase 1-bis (F3 accepted, F10 minimal fix) | V2 regressions (V1 is correct): rebinned errors, auto fudge, convolution with a resolution column. Reference behaviour: V1; each fix still agreed with Guido, with a regression test |
| Phase 4 | [#465](https://github.com/DAS-OATs/astrocook/issues/465) | continuum sensitive to spikes (urgent): done in phase 4 (spikes rejected as continuum anchors, masked running mean in km/s, local S/N; F30, F60; with F17, F56, F57 improved); spikes on real NIR data to be checked when the data arrive (F15) |
| Phase 0 (brought forward) | [#476](https://github.com/DAS-OATs/astrocook/issues/476), [#481](https://github.com/DAS-OATs/astrocook/issues/481) | ASCII export; same problem (missing `import re`, found by Ruff), fixed with a regression test |
| Right after phase 1, before phase 2 | [#442](https://github.com/DAS-OATs/astrocook/issues/442) | scatter plot not auto-triggering on region overlap: done in phase 1-bis (stitch sets the scatter style when the inputs overlap) |
| Right after phase 1, before phase 2 | [#480](https://github.com/DAS-OATs/astrocook/issues/480) | import of line lists as ASCII (urgent): done in phase 1-bis (Edit > Import Systems from ASCII) |
| With the natural-language log work | [#424](https://github.com/DAS-OATs/astrocook/issues/424) | logging of manual continuum edits |
| Phase 4 | [#526](https://github.com/DAS-OATs/astrocook/issues/526) | resolution as a single column and exact convolution (F8, F22, rest of F10, F9, F14, F47); design and decisions in the issue; PR 1 (single source, migration, no stamping, #615), PR 2 (GUI on the same resolver, #616), PR 3 (exact convolution, `fitting/lsf.py`, #617) and PR 4 (no more writes of the legacy sources; the profiles of the System Inspector with the physics and the convolution of the fitter, `convolve_flux` removed) done |
| Before release 2.0.2 (decided with Guido) | [#534](https://github.com/DAS-OATs/astrocook/issues/534) | which chi2 the fit uncertainties are scaled with (follow-up of F29): done, option (a) |
| Phase 4 | [#517](https://github.com/DAS-OATs/astrocook/issues/517), [#518](https://github.com/DAS-OATs/astrocook/issues/518) | from the findings review: flag and refit components at the parameter bounds, with F21 (`lya_auto` output); editable z window in the manual fits (F4). #517 in three PRs (agreed with Guido): flags on the components (PR 1), cleanup of the pipeline output (PR 2), GUI cues and "Select flagged" (PR 3) |
| Phase 3-bis (brought forward) | [#559](https://github.com/DAS-OATs/astrocook/issues/559) | found while cleaning up the Ruff warnings (phase 2-bis): `dbtur` not updated after a fit (F36); done |
| Phase 3-bis | [#560](https://github.com/DAS-OATs/astrocook/issues/560) | found while cleaning up the Ruff warnings (phase 2-bis): `norm` of the `rebin` recipe ignored (F37); done: no longer offered, a warning if set (Guido's decision) |
| Phase 4 | [#568](https://github.com/DAS-OATs/astrocook/issues/568), [#573](https://github.com/DAS-OATs/astrocook/issues/573) | found in phase 3: a single table of atomic parameters (F38, F39); recipe defaults differing between GUI and API (F40) |
| Phase 3-bis (priority) | [#575](https://github.com/DAS-OATs/astrocook/issues/575) | found in phase 3: a recipe that opens a new session is also recorded in the log of the parent session (F41); done |
| Phase 6-bis | #573, items 2 and 3 (closed with item 1 in phase 4) | a user configuration file for the recipe defaults and fit settings; the session log records the values actually used, not "default"; to be opened as issues when scheduled |
| Phase 6 | [#566](https://github.com/DAS-OATs/astrocook/issues/566) | centralize the log messages (started in phase 3 with `recipes/messages.py`), with the move to per-module loggers |
| After phase 5 | [#471](https://github.com/DAS-OATs/astrocook/issues/471) | combining order-separated spectra; touches `SpectrumV2` (coadd/stitch) |
| After phase 5 | [#478](https://github.com/DAS-OATs/astrocook/issues/478), [#479](https://github.com/DAS-OATs/astrocook/issues/479) | GUI features; wait for the GUI split |
| After phase 5 | [#446](https://github.com/DAS-OATs/astrocook/issues/446), [#437](https://github.com/DAS-OATs/astrocook/issues/437) | flux-calibration filters; independent, not urgent |
| When needed (not scheduled; Guido, October 2026) | — | transitions of the V1 table not in the single table of #568 (H2 with non-unique names, CO, Paschen, Brackett, SiI, Balmer beyond H_e, stellar and blended emission markers): added when an analysis needs them, each with its source (the V1 file is in `tests/reference/v1/data/atom_par.dat`) |
| To be decided (Guido) | [#455](https://github.com/DAS-OATs/astrocook/issues/455) (pending question) | standard multiplets (`STANDARD_MULTIPLETS`) for the 10 metal species with 3+ transitions in the table and no entry; unchanged by #568, which touched only the table of transitions; the choice of the multiplets used by the identification belongs to phase 7 |

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
| F4 | `VoigtFitterV2` | z is bounded to ±`z_window_kms` (default 20 km/s) around the initial guess; a guess further away cannot reach the true value. Intended? | behaviour | `test_fitting.py` (guess inside the window) | fixed (phase 4), [#518](https://github.com/DAS-OATs/astrocook/issues/518): default unchanged (20 km/s when adding, 100 when refitting); the System Inspector has a **Fit ±** field (after logN and b) for its manual fits (`test_wider_z_window_reaches_a_distant_guess`, `test_inspector_fit_window`) |
| F5 | `convert_y_axis` | Annotation declares 2 return values, the function returns 3 | code quality | snapshot `convert_y_axis` | fixed (phase 1-bis): annotation and docstring list the three returned values |
| F6 | architecture | `core/spectrum_operations.py` imports `gui.debug_utils`: core depends on the GUI | architecture | `test_debug_hooks.py` (`test_import_does_not_load_the_gui`) | fixed (phase 3-bis): the core requests the debug plots through a hook (`core/debug_hooks.py`) that the GUI connects to its signal |
| F7 | GUI | `gui/main_window.py` imports `matplotlib.pyplot` before `pyside_plot.py` selects QtAgg: importing the GUI without a `QApplication` (and no display) fails | architecture | workaround in `tests/test_gui_save_session.py` | phase 5 |
| F8 | resolution handling | The resolution lives in five places (per-component value, `resol` column, `SpectrumDataV2.resol`, meta `resol`/`RESOL`, meta `resol_fwhm`), each read with a different priority by the fitter, the `.acs2` loader, coadd and the GUI. The `.acs2` save does not write the attribute: it survives only through the meta keys (set by "Set properties"). The fitter stamps the resolution used into every fitted component, which then overrides later changes. R vs km/s is guessed from a threshold of 500. The J0103-1305 archive (older version) stores no resolution at all. To be analysed by Guido. Found in phase 3-bis: `SpectrumV2.rebin` computes the median of the rebinned `resol` column and stores it in `new_meta['resol']`, but then rebuilds `new_meta` from the original metadata, so the meta key keeps the value from before the rebin while the attribute is updated (45000 vs 45008 on a mock with R = 40000-50000) | behaviour, possible data loss | — (map in the PR that added this row) | fixed in phase 4 (#526, PR 1): the `resol` column is the single source (`SpectrumV2.resolution()`, R per pixel; unit from the column: dimensionless R, velocity FWHM); the attribute and the meta keys are migrated into it when a spectrum is built (`core/resolution.py`) and no longer read; the guess of the unit by value (500) survives only in that migration, for data saved before; the fitter no longer stamps its resolution into the components; components keep an override in R (legacy FWHM converted, stamped values within 1% cleared at load: Guido's decisions (a), (b)). PR 4: nothing writes the legacy sources any more (`set_properties`, `rebin`, `coadd`, `split`, the RESVEL loader, the V1 adapter, the mocks write the column only); the migration clears the attribute; the meta keys are left as they came, as they may be keywords of the original file |
| F9 | Voigt model | Components narrower than a pixel (e.g. Ly-alpha with b = 0.5 km/s, logN = 17.8 in J0103-1305) make the model depend on the sampling; with R = 45000, 1% of the 570-700 nm window still differs from the model saved in the session. Related to the oversampling check (item 9): the V2 model, sampled at pixel centres, overestimates the equivalent width of CIV by 4% at b = 1 km/s and 0.3% at 3 km/s (2.5 km/s pixels, R = 45000) | numerical | — | fixed in phase 4 (#526, PR 3): the model is computed on a fine grid (step a third of the narrowest b in the window, at least 0.5 km/s, and half the instrumental sigma) and averaged over each pixel; against the exact flux of mocks, 1.5e-5 to 3e-5 with F61 corrected (`test_lsf.py`). Timing (`bench_fitting.py`): single fit 0.008 → 0.010 s, group 1.44 → 1.65 s, chain on a mock 19.2 → 14.3 s, on J0103-1305 373 → 262 s (each line evaluated only where tau > 1e-6) |
| F10 | `convolve_flux` | With a per-pixel resolution column, pixels are grouped by resolution rounded to 100 (R) or 0.5 km/s, but selected by the raw value: pixels whose resolution is not already rounded are never convolved and the model is 0 there (a flat column at R = 45050 gives a zero model everywhere). Affects rebinned and coadded spectra (coadding three synthetic exposures at R = 42000-47500 gives a non-round column on 90% of the pixels and a zero Voigt model on 90% of them, `test_voigt_model_on_coadded_spectrum_is_not_zero`), RESVEL input, non-round "Set properties" values. In addition, even for a constant R the kernel width is computed at the median wavelength and with the median pixel size, so it is wrong far from the centre of the range (+53% / -26% at the ends of 380-780 nm). Proposed: a single exact convolution in velocity space, after the safety net is complete | numerical, wrong results | `test_model_column_resolution_is_not_zero`, `test_voigt_model_on_coadded_spectrum_is_not_zero`, `test_voigt_model_with_resolution_column` (regression) | selection bug fixed (phase 1-bis): pixels selected by the rounded resolution, NaN filled with the unconvolved flux. Kernel width far from the median wavelength and exact convolution: fixed in phase 4 (#526, PR 3, `astrocook/fitting/lsf.py`): each pixel takes the model convolved with a Gaussian of its own resolution and averaged over its edges, with `test_lsf.py` against the analytic result (1e-8) and against mocks with a flat and a varying resolution (5.9e-3 at R 45000 and 5.7e-2 at R 8900 before) |
| F11 | grouping | `VoigtModelConstraintV2.set_active_components` ignores `group_depth`: it always calls `get_connected_group` with the default depth 2, so the recipe option (1 = neighbours only) and the `group_depth=0` used by the refit have no effect | behaviour | `test_set_active_components_honours_group_depth`, `test_fit_component_group_depth_0_freezes_neighbours` (regression) | fixed (phase 1-bis, decided with Guido): `group_depth` honoured, default 2 everywhere; the internal fits of `refit_all`/`metals_auto` keep the depth-2 expansion they always had (`FIT_EXPANSION_DEPTH`), since depth 0 there was 2.75x faster but fitted blends worse (J0103-1305: chi2 +14% in the forest); the user can reduce the depth in "Fit component" and in the Bayesian fit |
| F12 | global state | Module-level dictionaries are modified at run time: the GUI rewrites the parameter docs of `EDIT_RECIPES_SCHEMAS["equalize_and_stitch"]`, and `add_component` registers custom series in `STANDARD_MULTIPLETS`. The changes persist for the whole process (and leak between tests) | architecture | `restore_multiplets` fixture in `test_issue_455_multiplet_subset.py` | open; phase 6-bis (moved from phase 5, agreed with Guido) |
| F13 | inspector | `SystemTableModel.set_highlighted_uuids` returns early when the selection is unchanged, so link colours are not recomputed when a link is created with the group still selected (model-level reproduction; full GUI not driven) | GUI | `test_link_created_with_group_selected_is_coloured`, `test_link_removed_with_group_selected_is_not_coloured`, `test_link_keeps_its_colour_when_another_is_added` (regression) | fixed (phase 3-bis): the colours are recomputed in `update_data`; a link keeps its colour when others are added or removed |
| F14 | grouping | Two components are considered blended when their line centres are closer than max(4 b, 30 km/s), independently of the instrumental resolution and of the line strength. At X-Shooter resolution (R ~ 8900) AlII 1670 and CIV 1550 still overlap at 60 km/s but are not grouped; worse in the NIR or for saturated lines. Possibly why #469 is still seen on X-Shooter data | behaviour | `test_blend_threshold_*` in `tests/test_issue_469_blends.py` | fixed in phase 4, after PR 1 of [#526](https://github.com/DAS-OATs/astrocook/issues/526) (the single per-pixel resolution), in parallel with its PR 3: the threshold is max(4 b, 30 km/s, k FWHM), with the instrumental FWHM at the lines from the `resol` column (or the resolution of the component, if set, as in the fitter; the larger of the two lines), k = 1 (Guido: the FWHM is the useful check; 4 b and 30 km/s were heuristics, kept as floors and as the fallback without resolution). The spectrum is passed to the grouping by every fit (recipes, `optimize_hierarchy`, Bayesian fit) and by the System Inspector. Only R < ~10000 (FWHM > 30 km/s) changes anything. With k = 1 the example above (60 km/s at R 8900, FWHM 34 km/s) is still not grouped: it would need k ~ 2; Guido kept k = 1 (#618). Mock benchmark: identical results on all 135 mocks (only the timings differ; the benchmark fits only the sparse metal doublets, and on two mocks at R 8900 no pair changed grouping), so a fit-level test replaces it (`test_fit_includes_component_blended_at_low_resolution`). Measured again on the merged code (master after #617 and #618, 3c4089c, against #617 alone, a7b02ce; both before F61, #619, not measured again after it): identical records on all 135 mocks |
| F15 | `equalize_and_stitch` | The scale factor is a median of the ratio of the two arms in the overlap: with narrow absorption (e.g. tellurics) seen at two different resolutions, it is biased by ~1-1.5% in a mock with equivalent-width-conserving tellurics around 1000 nm (no bias without tellurics). Possibly related to #477; to be checked on real X-Shooter VIS/NIR data | numerical | — (mock in the PR that added this row) | open |
| F16 | smoothing | For the same "std" parameter, V1 `_gauss_convolve` uses the kernel exp(-(dv/std)^2), i.e. sigma = std/sqrt(2), while V2 `smooth_spectrum` uses sigma = std: V2 smooths sqrt(2) wider. Affects the continuum smoothing (default 500 km/s: effective sigma 354 km/s in V1, 500 km/s in V2) and any other use of the parameter. Which convention to keep is Guido's decision | numerical | `test_smoothing_width_differs_by_sqrt2` (characterization); `test_estimate_auto_smoothing_is_a_fwhm` | fixed (phase 1-bis, decided with Guido): the continuum smoothing is given as a FWHM (`smooth_fwhm`, default 800 km/s, within 0.5% of the V1 continuum on J0103-1305; FWHM 500 would change it by up to 10% in the forest); `smooth_std` (sigma) kept as a deprecated alias |
| F17 | continuum (forest) | On the J0103-1305 window, with the same parameters and matched smoothing, V1 `clip_flux` and V2 `find_absorbed` + `fit_continuum` agree on the red side (masks 99.98%, continuum within 0.4%), but in the Lyman-alpha forest V2 marks more pixels as unabsorbed (57.6% vs 54.8% of the window, masks agree on 92.7% of the forest pixels) and its continuum is lower by 1.2% (median; 5th percentile -4.6%). Not explained by the fudge (F3). Relevant to #465 and to the known low forest continuum | numerical | `test_continuum_forest` (strict xfail) | partly addressed by #465 (measured on the mock benchmark (2.1.0 + F49/F50 → #465): forest continuum 0.936 → 0.943 (5th percentile 0.820 → 0.830), red side unchanged (1.002)); V1 parity no longer sought for the continuum (V2 is measured against the truth on mocks instead) |
| F18 | dependencies | GitHub Dependabot reports 36 vulnerabilities on the default branch (2 critical, 8 high, 13 moderate, 13 low), seen on push in October 2026; not yet analysed | security | — (Dependabot alerts) | fixed (phase 2): all 36 alerts were on `docs_backup/Gemfile.lock` (Ruby gems of the frozen V1 Jekyll docs, used only for a local rebuild); `Gemfile` and `Gemfile.lock` removed |
| F19 | `refit_all` (reproducibility) | With `group_depth >= 1` the groups are built by taking the next seed from a set of component UUIDs (`next(iter(todo_uuids))`), and the depth-limited friends-of-friends group depends on the seed. UUIDs are random (`uuid4`) and string hashes are randomized per process, so the same input with the same parameters gives a different result at every run. Full chain on J0103-1305, 500-700 nm, two runs: `doublets_auto` gives the same 51 CIV components, but 4 differ by more than 0.1 dex in logN (up to 2.35 dex: logN 10.0 vs 12.35 at z = 3.2836); `lya_auto` gives 477 vs 478 components, ~20 not matched within 5 km/s; the normalised model differs by more than 0.01 on 1.8% of the pixels (up to 0.83). On the mock the differences are ~1e-6 dex. Even with a fixed hash seed and seeded UUIDs the result depends on the machine: on a GitHub runner the J0103-1305 chain gives 477 instead of 478 components after `lya_auto` and a normalised model differing by up to 0.84, with a continuum identical to 1e-15; the mock differs by up to 1e-3 dex in logN. Affects `doublets_auto`, `lya_auto`, `components_auto`, `forest_auto`, "Refit all" from the GUI (schema default `group_depth = 2`) and the refit after importing systems (`recipes/edit.py`) | numerical, reproducibility | `test_refit_all_does_not_depend_on_the_uuids` (regression); `tests/full_chain.py` | fixed (phase 1-bis): seeds taken in order of increasing z; two J0103-1305 full chains with different hash seeds and UUIDs now give identical results. Floating-point differences between machines remain (statistical comparison in the full-chain snapshots) |
| F20 | `VoigtFitterV2` (uncertainties) | Errors are `sqrt(diag(pinv(J^T J)))`: a negative diagonal element (numerically ill-conditioned fit) gives a NaN uncertainty with only a `RuntimeWarning`, no log message and no flag on the component. 1 of 478 components in the J0103-1305 full chain. Related to the chi-squared scaling check (item 9) | numerical | snapshot `full_chain_j0103_forest`; `test_nan_uncertainties_are_reported` (regression) | fixed (phase 1-bis): a warning names the components with NaN errors, numpy RuntimeWarning silenced; the NaN error itself is the flag (no new field) |
| F21 | `lya_auto` (output quality) | On the J0103-1305 full chain (500-700 nm): 39 Ly-alpha components end with b < 10 km/s although `min_b = 10` (the filter runs before the final refit, which moves the remaining components again); components stuck at the parameter bounds are kept (12 at logN = 10, 12 at b = 1 km/s, Ly-alpha and CIV), including two Ly-alpha with logN = 20.2 and 21.1 at b = 1 km/s (z = 3.708, 3.506), unphysical in the forest | behaviour | snapshot `full_chain_j0103_forest` | fixed in phase 4 (#517, agreed with Guido): flags on the components (PR 1, #627); in `lya_auto` the Ly-a below `min_b` are removed and the rest refitted twice, and those still below after the last refit removed without a refit; the pipelines remove the components at logN 10 (flag `logN_at_lower_bound`); those at b = 1 km/s are kept with their flag (PR 2). J0103-1305 chain: Ly-a with b < 10 km/s 33 -> 0, components at logN 10 3 -> 0, at b = 1 km/s 8 -> 2, 40 forest components fewer; forest chi2 per pixel 31.4 -> 30.4, over-absorbed pixels 3.4% -> 2.6%. GUI cues in PR 3: the cell of the parameter at a bound highlighted, the ID of a component without uncertainties in red, a warning sign on the ID, dashed ticks, "Only flagged" and "Select flagged" |
| F22 | resolution handling | `apply_expression`, `mask_expression`, `smooth_column` and `split` rebuild `SpectrumDataV2` field by field and drop the `resol` attribute (reset to 0); `smooth`, `rebin`, `stitch` and `set_properties` keep it. The Voigt model changes (by 0.05 in flux in a synthetic test at R = 45000) only when the resolution is in the attribute alone: `set_properties`, the `.acs2` loader and the mocks also store it in the metadata, which is kept, so the usual paths are not affected. One more case for F8 | behaviour | `test_some_operations_drop_the_resolution_attribute` (characterization) | fixed in phase 4 (#526, PR 1): the resolution is a column, which these methods keep; `rebin`/`coadd` work in R whatever the unit of the column, and an all-NaN column no longer counts as a source in `coadd` |
| F23 | `flux.rebin` recipe | A non-numeric parameter (e.g. `dx = "abc"`) is logged as a parsing error but the recipe goes on and raises `UnboundLocalError`, instead of returning 0 like the other recipes. In the GUI the worker catches the exception, so it is not a crash | error handling | `test_rebin_bad_input_returns_0` (regression) | fixed (phase 1-bis): the recipe returns 0 on a parsing error |
| F24 | `fit_bayesian` (MCMC) | `run_mcmc` discards a fixed burn-in of 200 steps, and the recipe passes only `nsteps` (GUI parameter, default 500): with `nsteps <= 200` the chain is empty, the medians are NaN and they overwrite z, logN and b of the fitted components, and the model becomes NaN. No warning. Reported at once (loses data; undo restores the previous state) | data loss | `test_fit_bayesian_few_steps_keeps_parameters` (regression) | fixed (phase 1-bis): burn-in defaults to min(200, nsteps // 2), an empty chain is refused, non-finite results never overwrite the components |
| F25 | `delete_component` | The recipe recomputes the model after deleting, which needs a continuum: on a spectrum in physical units without one (and with a median flux outside 0.1-10, otherwise it is taken as normalised), `VoigtFitterV2` raises, the error is logged and the session is returned unchanged. The component is not deleted (inspector, "Clean negligible": in the GUI only "Process completed with warnings"). Seen on the small J0103-1305 archive (764 components, no continuum). Fitting a component instead triggers the continuum estimate first | behaviour | `test_delete_component_without_continuum` (regression) | fixed (phase 1-bis): without a continuum the components are deleted, the model is not recomputed and the stale `model` column is removed, with a warning |
| F26 | `.acs2` session log | The `.acs2` archive does not keep the V2 session log: `_serialize_v2_metadata` writes `dataclasses.asdict(HistoryLogV2)`, which is `{}` because the entries are set in `__init__` rather than as dataclass fields (the same cause found in phase 0 for `.acs`), and the loader does not read a log back anyway. Checked through the GUI: after "Set properties" and "Smooth" the archive has `log_history_json: {}` and the reloaded session has an empty log. The GUI tells the user to save as `.acs2` to keep the log. Reported at once (loses data) | data loss | `test_session_log_is_saved`, `test_loaded_log_is_kept_by_the_gui` (regression) | fixed (phase 1-bis): the log entries up to the current one are saved and restored on load, also in the GUI |
| F27 | `.acs2` identifications | `region_identifications` is saved as a JSON string and reloaded as a dict with integer keys. The identification viewer calls `json.loads` on it and shows "Error parsing metadata" (the plot tooltips handle both forms) | GUI | `test_identification_viewer_after_round_trip` (regression) | fixed (phase 1-bis): the viewer accepts both the JSON string and the dict |
| F28 | `.acs2` Bayesian samples | The posterior samples stored by `fit_bayesian` in the system list metadata (`bayesian_results`, used by the corner plot) are dropped on save without a warning: the FITS header cannot hold them. The fitted values and errors are kept | data loss (recomputable) | `test_bayesian_results_are_saved` (strict xfail) | open |
| F29 | `VoigtFitterV2` (uncertainties) | The errors are sqrt(diag((J^T J)^-1)) with residuals weighted by 1/dy, not scaled by the reduced chi2, although the comment in `fit` says they are. V1 scaled them (lmfit `Model.fit`, `scale_covar=True` by default): a V1/V2 difference. Monte Carlo on mocks: unscaled errors are right when the error spectrum is right, and wrong by the same factor when it is not (e.g. 2 times too small with the error spectrum underestimated by 2, reduced chi2 3.8); scaled errors are right in all cases. Scaling, or reporting the reduced chi2 next to the errors, is Guido's decision | numerical, V1/V2 difference | `tests/test_error_scaling.py` (regression) | fixed (phase 1-bis, decided with Guido): errors scaled by sqrt(reduced chi2), as V1; the reduced chi2 is stored in each fitted component. On real data the scaling uses the chi2 of an intermediate pipeline state and of unrelated pixels: [#534](https://github.com/DAS-OATs/astrocook/issues/534) |
| F30 | `running_mean`, `running_std` | Both use `np.convolve(mode='same')`, which pads with zeros: the first and last h values are biased towards 0 (4/7 of the true mean at the edge for h = 3). V1 `running_mean` repeats the first and last full-window value and ignores NaNs. Used by `find_absorbed_regions` on each interval, with half-windows of hundreds of pixels (5000 km/s / 8 on J0103-1305). Replacing it with the V1 function raises the agreement of the forest mask with V1 from 92.7% to 95.7% and the 5th percentile of the continuum ratio from -4.7% to -3.0% (median unchanged, -1.2%): a partial cause of F17 | numerical, V1/V2 difference | `test_running_mean_and_std_pad_the_edges_with_zeros` (characterization) | fixed in phase 4 for the continuum (#465): `find_absorbed_regions` no longer uses `running_mean`, but a masked mean over a window of the smoothing length in km/s on the pixel grid, exact at the edges (`_masked_window_mean`). The edge bias was worth 0.81 → 0.93 within 2000 km/s of the edges of the forest interval on a mock at S/N 30. The plain normalized convolution of the earlier branch (window in pixels of the clipped array) failed at S/N ≥ 300 (forest continuum up to 1.66 times the true one: the window grows as pixels are clipped, over wavelengths where the continuum varies by more than its noise); the window in km/s does not. `running_mean` and `running_std` themselves still pad with zeros (used by `calculate_running_std`): open, minor |
| F31 | dead code | `compute_identification_signal` and `find_signal_peaks` (`core/spectrum_operations.py`) are not called anywhere: the V2 identification goes through `SpectrumV2.identify_lines`. Section 2 of this roadmap describes the identification as based on the former. Found by mutation testing (no test reaches them) | code quality | — | fixed (phase 2): both functions removed, section 2 updated |
| F32 | `VoigtFitterV2` (broad lines) | `_residual_function` appends a penalty 5 (b - 50) for every free b > 50 km/s to the residuals ("soft penalty to prevent runaways"). Broad lines are biased: on mocks (Ly-alpha, logN 13.8, R = 45000), at S/N 20 a line with b = 60, 80 or 120 km/s is fitted with b = 50 km/s and logN low by 0.04, 0.11 and 0.35 dex; even at S/N 1000, b = 120 comes out as 72. The penalty also enters chi2 and the uncertainties. Affects OVI, broad Ly-alpha (BLAs), blends fitted with one component | numerical, wrong results | `test_fit_broad_line`, `test_b_max_bounds_the_fit`, `test_pipelines_bound_b_and_manual_fits_do_not` (regression) | fixed (phase 1-bis): penalty removed; hard bound b <= 50 km/s (`max_b`) in the automatic pipelines, 300 km/s in the manual fits (without any bound `lya_auto` on J0103-1305 takes 4 times longer and gives spurious broad lines); components at the bound: [#517](https://github.com/DAS-OATs/astrocook/issues/517) |
| F33 | tests (`.acs2` I/O) | `tests/test_v2_io_integrity.py` wraps its checks in `try/except` and only prints the error, so it can never fail. It currently prints "INTEGRITY CHECK FAILED": it checks that "CIV" is in `region_identifications` as a string, while after loading it is a dict (same mismatch as F27). Found in phase 1-bis while running the slow tests | test quality | `test_integer_aux_column_round_trip` | fixed (phase 1-bis, decided with Guido): test removed; its only check not covered by `tests/test_acs2_round_trip.py` (integer `abs_ids`) added there |
| F34 | GUI | `gui/main_window.py` sets the main window icon from `assets/logo_3d.png`, which does not exist (the available logos are `logo_3d_HR.png` and `logo_3d_LR.png`): the icon is silently missing. Found in phase 2 while moving `assets/` into the package | GUI | — | fixed (phase 2-bis): icon from `logo_3d_HR.png`, smoke test checks it is not null |
| F35 | security | `.acs2` (`core/session.py`) and `.acs` (`io/v1_stubs.py`) archives were unpacked with a plain `tarfile.extractall`: a crafted archive could write files outside the temporary folder (path traversal). V1 had a check (`legacy/session.py`, `safe_extract`) that V2 lost. Found in phase 2 from a Python 3.14 deprecation warning in the slow tests | security | `tests/test_safe_extract.py` | fixed (phase 2): `extractall(..., filter="data")`, which refuses members outside the folder |
| F36 | `VoigtFitterV2` (errors) | After a fit, `dz`, `dlogN` and `db` are updated, but `dbtur` is computed (`dbtur_new`) and never stored: it keeps its value from before the fit, even when `btur` is free or frozen by the user (where the other errors become 0). Found in phase 2-bis while cleaning up the Ruff warnings | numerical, wrong uncertainties | `test_voigt_fitter_physics.py` (`test_fit_stores_the_btur_error` and three more; regression) | fixed (phase 3-bis), [#559](https://github.com/DAS-OATs/astrocook/issues/559): `dbtur` stored like the other errors, also by `evaluate_uncertainties`; 0 when `btur` is frozen (the default) |
| F37 | `rebin` recipe (`recipes/flux.py`) | The `norm` parameter ("Return normalized spectrum, if continuum exists") is parsed and ignored: the rebinned spectrum is never normalized. Found in phase 2-bis while cleaning up the Ruff warnings | silently ignored parameter | `test_rebin_does_not_offer_norm`, `test_rebin_norm_is_ignored_with_a_warning` | fixed (phase 3-bis, decided with Guido), [#560](https://github.com/DAS-OATs/astrocook/issues/560): not offered by the dialog; `norm='True'` (scripts, logs) logs a warning; the V1 normalization is not implemented |
| F38 | atomic data (`core/atomic_data.py`) | `ATOM_DATA` (V2, 210 transitions) differs from `data/atom_par.dat` (V1, 829) for some values. Many differences look intended (updated Lyman oscillator strengths, masses filled in where V1 had 0), but some look like transcription errors in V2: wavelength of `OIII_832` (832.29 vs 832.93 Å), `PIII_913` (913.40 vs 913.97 Å), `PV_1117` (1117.80 vs 1117.98 Å), `FeI_2501` (2501.19 vs 2501.89 Å), and the damping constant of `ArI_1066` (1.32e7 vs 1.32e8 s^-1). The V1 values agree with the standard line lists. Found in phase 3 (PR 3) | numerical | — | wrong values fixed (phase 4, #568, PR 1): the four wavelengths were mistyped in the conversion from nm to Angstrom (91.39683 nm written 913.39683 A; exactly these four rows of the table show the pattern), the damping constant of ArI 1066 was off by ten; the V1 values, confirmed against Morton (2003) except OIII 832, which Morton does not cover. PR 2: AlII 1670 to the Morton wavelength (1670.7886 A); [NeIII] 3869 and [NeV] 3346 to vacuum (they were in air); TiII_3066, a line no reference lists (Morton, the atom.dat of VPFIT from P. Noterdaeme), kept with f = 0 as in V1 and TiII_3067 added. Decided by Guido: the newer values of the VPFIT atom.dat are adopted line by line, each with its source (f of ZnII 26% higher, of NiII up to 32% lower, of SII 5% lower; CIV 1550 at 1550.77845 A; several gamma), in the PR of the single table (with a source column), together with F39. PR 3 (single table): `data/atomic_data.ecsv` replaces `ATOM_DATA` and `atom_par.dat`, with the source of each wavelength, f and gamma (Morton 2003 for most; `V1`/`V2` where the original source was not recorded); adopted from the atom.dat of VPFIT, with their original sources: f of ZnII 2026, 2062 (Kisielius et al. 2015), NiII 1317-1751 (Boisse & Bergeron 2019), SII 1250-1259 (Kisielius et al. 2014), wavelengths of CIV 1548, 1550 (Petitjean & Aracil 2004); gamma of CI 1560 to Morton (1.27e8, was 1.17e8). Not adopted, as their source is not clear: the other gamma of VPFIT that differ (OI 971, 976, MgI 2026, MgII 2796, 2803, CoII 1466, FeI 3441; placeholders of 2e8 for OII, OIII, FeII 1081), SII 1253 at 1253.811 A, f of NI 1134b (0.0278, Morton 0.0287), SnII 1400. The V1-only transitions (H2 with non-unique names, CO, Paschen, Brackett, SiI, Balmer beyond H_e, stellar and blended emission markers) are not in the table: V2 never used them (not in `xem_d`); adding some is a separate choice |
| F39 | `trans_parse` (GUI) | The V1 series (`trans_parse`, used by the GUI when a series is not in `STANDARD_MULTIPLETS`) keep the order of `atom_par.dat`: `'Ly'` goes from `Ly_lim` to `Ly_a` (also `H`, `D`), and the redshift cursor takes the first transition known to `xem_d` as reference, so with the series `Ly` the cursor z refers to `Ly_lim`, not `Ly_a`. Also, only whole strings are looked up in `STANDARD_MULTIPLETS`: `MgII` gives the doublet, `MgII;Ly_b` the six MgII lines of `atom_par.dat` (1025-2803 Å). Building the series from `ATOM_DATA` instead would change this order and add `OI_988a/b` to `OI`. Kept as in V1 in phase 3 (PR 3) | behaviour | `test_atomic_data.py` (`test_trans_parse_hydrogen_order`), `test_gui_series_lookup.py` | fixed in phase 4 (#568, PR 3): the series are built from the single table, in its order, so `Ly`, `D`, `H` start from the strongest line (Ly_a) and the aliases `Ly-ab`, `Ly-abg` too; `OI` now includes `OI_988a/b`; the lookup of whole strings in `STANDARD_MULTIPLETS` is unchanged |
| F40 | recipe defaults | Each recipe defines its defaults twice, in the schema (GUI dialogs) and in the method signature (scripts, API): 15 optional parameters differ, e.g. `refit_all` `group_depth` 2 (GUI) vs 0 (API), `coadd` `dx`/`xunit` 10 km/s vs 0.01 nm, `optimize_system` `aic_penalty` 5 vs 0, so the same recipe gives different results from the GUI and from a script. Found in phase 3. Also within the GUI: the "Refit" and "Fit component" actions of the System Inspector pass only some parameters, so they use the signature defaults (`refit_all` with `group_depth` 0 and `max_nfev` 100; `fit_component` with `max_nfev` 2000), unlike the "Refit All" dialog (2 and 200); `metals_auto` calls `refit_all` with the signature defaults too (phase 4) | behaviour / numerical | `test_recipe_defaults.py` (15 strict xfails) | issue [#573](https://github.com/DAS-OATs/astrocook/issues/573), phase 4: test added, values to be decided |
| F41 | session log (GUI) | A recipe that opens a new session (`extract_preset`, `split`, `coadd`, ...) is recorded in the log of the parent session too, which does not change: the parent log no longer matches its states, and the extra entry stays in later operations, in the saved `.acs2` and in a replay of the log. Found in phase 3 (PR 6) | reproducibility | `test_gui_smoke.py` (`test_branch_does_not_change_the_parent_log`, `test_branch_after_undo_keeps_the_parent_redo`; regression) | fixed (phase 3-bis), issue [#575](https://github.com/DAS-OATs/astrocook/issues/575): the entry goes only into the log of the new session; logs saved earlier left as they are |
| F42 | opening V1 `.acs` archives | The V1 reader, kept as it was by the native reader of phase 3 (`io/v1_archive.py`), sets `btur` to 0 and drops `dbtur` and `chi2r` from the V1 system table; the flux unit of V1 (erg/(cm2 s Angstrom)) is not kept (dimensionless). In the two sample archives `btur` is 0, `chi2r` is lost. Found in phase 3 (PR 8b) | data | `test_v1_archive.py` (snapshots, `test_v1_columns_are_read`) | fixed (phase 3): the three columns are read (`chi2r` as the reduced `chi2` of V2), the flux in erg/(cm2 s Angstrom) as in the V2 `astrocook_fits` loader (before, `y_convert` failed on V1 sessions) |
| F43 | V1 constraints | V1 links are stored as `lines_voigt_<id>_<par>` (e.g. `lines_voigt_14_b*0.65...`), but the conversion to V2 (`VoigtModelConstraintV2._convert_v1_to_v2`) only translates `p[<id>]`: the linked parameters (124 in the sample archive) are frozen at their V1 values with no target, the expression fails silently at each evaluation, and the System Inspector cannot show the link (reported by Guido). Found in phase 3 (PR 8b) | behaviour | `test_v1_constraints.py` | fixed (phase 3): the V1 reader translates the links into V2 constraints keyed by UUID (target, `p['<uuid>'].<par>` syntax), so that grouping, fit, `.acs2` and System Inspector see them; `logN` constraints are no longer lost (the parameter name was lower-cased). A refit of a V1 session now follows the links |
| F44 | opening archives | A spectrum the V1 reader could not read was replaced with mock data (flux 1 between 100 and 500 nm, no systems) with only a log warning: an `.acs2` archive without `_meta.json`, a damaged `.acs`, or a raw file whose format has no V2 loader. Saving such a session overwrote the data. Found in phase 3 (PR 8b) | data loss | `test_v1_archive.py` | fixed (phase 3, PR 8b): `.acs2` without metadata read with the V2 readers (without constraints, log, identifications), the other cases are errors |
| F45 | `.acs2` archives | The archive does not record the Astrocook version that saved it, so a problem tied to a version cannot be recognised on load: e.g. the spurious log entries of F41 in archives saved with 2.0.2 and 2.0.3, left as they are because they cannot be told apart from correct ones. Proposed: write the version in the archive metadata from now on (small, separate PR). Found in phase 3-bis (#575) | data provenance | `test_acs2_round_trip.py` (`test_archive_records_the_version`, `test_archive_without_version`) | fixed (phase 3-bis): `astrocook_version` in `_meta.json`, read back as `SessionV2.saved_with_version` |
| F46 | `edit.set_properties` | The values are not checked: z_em = -5 and R = -1 are accepted without a warning (also from the GUI dialog). An agent passing a wrong value gets no signal. Found in the MCP pilot (PR 2) | behaviour, error handling | `test_mcp_server.py` (`test_unphysical_properties_are_accepted`, characterization) | open: phase 6-bis (typed, checked parameters) |
| F47 | `set_properties` (GUI) | When R is changed, the GUI (`_process_recipe_result` in `gui/main_window.py`) writes the new value into every component; the recipe does not, so a script or the MCP server setting R on a session with components leaves the old per-component resolution, which the fitter uses (F8). Same kind of logic as the session log, kept in the GUI instead of the recipe. Found in the MCP pilot (PR 2) | behaviour | — | fixed in phase 4 (#526, PR 1): neither the GUI nor the recipe writes the resolution into the components; `set_properties` keeps a column that varies with wavelength unless `overwrite_resol_col` (default now False, Guido's decision (d)); with `'_current_'` it no longer flattens it |
| F48 | `identify_lines` (docs) | The schema and the docstring of the recipe describe `mask_col` as "True=Unabsorbed"; the `abs_mask` column made by `find_absorbed` is True for the absorbed pixels (its own description, and how `identify_lines` uses it). The schema texts are what the GUI dialogs show and what the tools of the full MCP server will show to a model. Found in the MCP pilot (PR 3) | documentation | — | open: phase 6-bis (descriptions reviewed) |
| F49 | `identify_lines` (candidate systems) | A doublet gives one candidate system per line: the two redshifts differ by ~0.5-1 km/s and survive the de-duplication, which rounds z to 5 decimals. On the full-chain mock the 3 CIV systems (4 components) give 6 candidates; on J0103-1305, 19 pairs of CIV candidates and 2 of MgII are closer than 2 km/s. `populate_from_identification` would add a component for each. Also, the candidate systems keep neither their region nor their score. Found in the MCP pilot (PR 3) | behaviour, structured output | `test_identify_lines.py`, `test_mcp_server.py` (`test_identify_lines_one_candidate_per_doublet`, `test_doublets_auto_fits_a_blend_with_one_component`) | fixed in phase 4: both regions of a doublet keep their label, only the region of the stronger line gives the candidate systems. Measured on the mock benchmark (2.1.0 → fix): duplicates CIV/SiIV/MgII 505/267/345 → 165/69/95 (those left are distinct components of a true system), purity 0.986/0.967/0.980 → 0.993/0.989/0.988, completeness unchanged; fit: spurious components 428 → 11, flagged 413 → 38, but true components matched 1534 → 1157, pull widths z/logN/b 1.8/2.3/2.1 → 4.3/5.5/4.5, system column spread 0.100 → 0.093 dex with pull width 1.24 → 2.53. The duplicate candidate, within a resolution element of the first, worked as a two-component start: the fit moved it onto a second true component when there was one, and made a spurious one otherwise. A variant keeping the candidates of the second line farther than c/R from those of the first (1171 matched, 17 spurious) was not adopted. The candidate systems still keep neither their region nor their score (phase 6-bis). Decomposition: F58 |
| F50 | `components_auto` (`doublets_auto`, `lya_auto`, ...) | The pipeline always makes a new absorbed mask with the default parameters of `find_absorbed` (5000 and 400 km/s, kappa 2) and stores it in `abs_mask`, so the mask of an earlier "Auto-estimate continuum" with other parameters is replaced without notice (the continuum is kept); the GUI dialog and the tool expose no mask parameter. Also, the region left out by `region_limit` is cleared in the stored `abs_mask`: after `doublets_auto` on the full-chain mock no forest pixel is flagged as absorbed (36% before), so a later "Fit continuum" on that mask would treat the whole forest as unabsorbed (the automatic pipelines make a new mask each time and are not affected). Found in the MCP pilot (PR 4) | behaviour | `test_mcp_server.py` (`test_doublets_auto_keeps_the_mask_of_the_continuum`, `test_doublets_auto_without_a_mask_stores_the_whole_new_one`) | fixed in phase 4 (Guido's choice): the pipeline still makes a fresh mask with the default parameters (to avoid stale masks between recipes, as its comment says) and restricts it to `region_limit`, but uses it only inside; the session it returns keeps the mask of the user, or, if there was none, the whole new mask. Results unchanged. A mask parameter for the pipelines: phase 6-bis |
| F51 | `identify_lines` (forest) | Every absorbed region of the forest without a metal candidate is assigned to Ly-a (score 0, no test), even when Ly-a is not among the requested `multiplets`: `multiplets='CIV'` on the full-chain mock gives 233 Ly-a and 6 CIV candidates. In the MCP pilot the Ly-a candidates crowded out the metals in the output (833 on J0103-1305). Found in the MCP pilot (report) | behaviour | — | open: phase 4 / 7 |
| F52 | `VoigtFitterV2` (uncertainties) | Components that end at the bounds of a fit get uncertainties far beyond any meaning, stored as results without a warning: in the `doublets_auto` fit of J0103-1305 (MCP pilot), four CIV components at z = 3.927 with logN = 10.0 and b = 1.0 km/s have dz up to 131, dlogN up to 3800 and db up to 9e4 km/s; 21 of the 101 new components have dlogN > 1 or db > b. Related to F20, F21. Found in the MCP pilot (report) | numerical, error handling | — | open: with [#517](https://github.com/DAS-OATs/astrocook/issues/517); since its PR 1 the components at a bound carry a flag (`ComponentDataV2.flags`), so their errors can be told apart |
| F53 | continuum (forest at high z) | On an ESPRESSO spectrum of J0959+1312 (z_em ~ 4.06), `estimate_auto` flags 100% of the pixels between 525 and 600 nm as absorbed, so the continuum there is a straight interpolation between the edges of the interval (flat, above the 99th percentile of the flux, ignoring the wing of the Ly-a emission); with `template=True` it follows the template, still with no free pixel, and lies below the flux peaks at 565-590 nm (6% of the pixels above C + 2 sigma). Reported by the agent in Guido's MCP trial; not reproduced here (the spectrum is not in the repository) | numerical | — | open: phase 7 (emission model); with #465 |
| F54 | identification (tellurics) | No telluric masking: the O2 A band (759-770 nm) gives spurious doublets with good scores, CIV at z = 3.927 in J0103-1305 (4 components at the bounds, F52) and CIV z = 3.928 with MgII z = 1.733 in J0959+1312 (fitted together, reduced chi2 12). Found in the MCP pilot (trials) | behaviour | — | open: phase 7 (identification); telluric data in the pre-release `data-telluric-2026` |
| F55 | `identify_lines` (z range) | The doublet candidates are not limited to z <= z_em (plus a proximity margin): `z_em` is passed to `_find_kinematic_doublet_candidates` but not used. In Guido's MCP trial `doublets_auto` fitted a SiIV at z = 4.33 with z_em ~ 4.06 (16000 km/s beyond, in a telluric band). The schema of `doublets_auto` says "red side (z < z_em)". Found in the MCP pilot (Guido's trial) | behaviour | `test_identify_lines.py` (`test_candidates_beyond_z_em`) | fixed in phase 4: no candidate system beyond z_em by more than 3000 km/s (`Z_EM_MARGIN_KMS` in `core/spectrum.py`, a constant with its purpose in a comment; Guido: it may become a parameter in phase 6-bis); associated absorbers and the uncertainty of z_em (~1000 km/s) stay within it. Not measured by the mock benchmark, which has no absorption beyond z_em |
| F56 | continuum (forest) | Measured on the mock benchmark (phase 4 baseline, 135 mocks): `estimate_auto` gives a forest continuum (Ly-b to Ly-a) low by 6.4% in the median, 13% at S/N 10, 12% at z_em = 4.5, with 5% of the forest pixels more than 18% low (34% at z_em = 4.5); the red side is within 0.2%. Quantifies F17 and the known low forest continuum (F53 on real data) | numerical | `benchmarks/validate_mocks.py`, baseline | improved by #465 (measured on the mock benchmark (2.1.0 + F49/F50 → #465): forest continuum 0.936 → 0.943 (5th percentile 0.820 → 0.830), red side unchanged (1.002); S/N 10: 0.875; S/N 5: 0.80, the limit of clipping). Also found: spikes falling on absorbed pixels and raised to the continuum level become anchors and raise the forest continuum slightly, i.e. the bias comes from the shortage of anchors at the continuum level. Open: phase 7 (emission model) |
| F57 | identification (completeness) | On the mock benchmark the metal doublets are found only well above the detection limit: none of the systems whose weaker line is below 3 sigma, 5% of those at 3-5 sigma, 30% at 5-10 sigma, 91% above 10; completeness over the systems above 3 sigma 0.85 for CIV, 0.69 for SiIV, 0.91 for MgII. Purity is high (0.97-0.99) | behaviour | `benchmarks/validate_mocks.py`, baseline | improved by #465 (wider window on the red side, 400 km/s instead of ~100: weak lines are flagged better): completeness 3-5 sigma 0.03 → 0.10, 5-10 sigma 0.35 → 0.44, CIV 0.857 → 0.874. Why: `find_absorbed_regions` is both the continuum step and the line detector (`identify_lines` takes its regions from `abs_mask`); a weak line pulls down the running mean of a narrow window by a good fraction of its own depth, so its pixels do not reach 2 sigma; with a window 4 times wider the line weighs 4 times less, and purity is unchanged (clipping still at 2 sigma). Direction for phase 7: decouple the detection from the continuum, e.g. search the lines on the flux normalised by the final, smooth continuum, with a filter matched to the expected line width instead of a per-pixel threshold. Open: phase 7 |
| F58 | `doublets_auto` (decomposition) | On the mock benchmark the systems with several components are decomposed differently from the truth: where one fitted component stands for several true ones the pulls of the matched component are 4-6.5 wide, and spurious components appear next to real ones (F49); the total column of each system is right (median difference -0.002 dex, pull width 1.0-1.1 at R >= 45000). The pipeline adds one component per candidate and never adds components from the residuals | behaviour | `benchmarks/validate_mocks.py`, baseline | open: phase 4 or 7 (Guido, October 2026): since F49 the pipeline starts from one component per peak, so a system of several components is fitted with fewer (true components matched 1534 → 1157, pull width of the system column 1.24 → 2.53); to try: `optimize_system` (adds components while the AIC improves, penalty 5, up to 5 more) after the fit of `doublets_auto`, measured on the benchmark (quality and time); then the fit driven by the residuals of phase 7 |
| F59 | fit at low resolution | At R = 8900 the total column of the metal systems is off by 0.39 dex (spread) with a pull width of 2.5, against 0.03-0.07 dex and ~1 at R = 45000-140000. Possible causes: unresolved saturated components (physics), and the convolution kernel computed at the median wavelength (rest of F10, up to +53%/-26% in width at the edges of a window) | numerical | `benchmarks/validate_mocks.py`, baseline | partly fixed by #526 (PR 3, exact convolution): at R 8900 the spread of the system column goes from 0.175 dex (master after #465; 0.39 in the 2.1.0 baseline) to 0.102, the median reduced chi2 from 2.45 to 1.61 and the flagged components from 30 to 16; R 45000 and 140000 unchanged (0.046 and 0.015 dex). After F61: 0.095 dex at R 8900, 0.041 at R 45000. What is left at R 8900 (median reduced chi2 1.6) is probably the physics of unresolved saturated components; to be looked at in phase 7 |
| F60 | continuum (low S/N) | `find_absorbed_regions` excluded from the continuum anchors each pixel with y < 3 dy, meant to exclude the regions with no signal: at S/N ~3 this drops the lower half of the noise, and the continuum came out 26% too high on the mocks (the mean of a half Gaussian, +0.80 sigma: 26.6% expected, 26.5% measured). Also measured: the one-sided clipping at -2 sigma biases the continuum up by 0.055 sigma (truncated Gaussian), 0.6% at S/N 10, as found on the red side; the V1 fudge probably corrected it (F3) | numerical | `test_continuum_robustness.py` (`test_low_snr_red_side`) | fixed in phase 4 (#465): regions with no signal are those where the median S/N over 101 pixels is below 1; S/N 3: 1.265 → 1.02. The 0.055 sigma bias: phase 7 |
| F61 | Voigt model | The damping parameter of the Voigt profile is computed with the observed wavelength of the line centre, a = Gamma lambda_c / (4 pi b) with lambda_c = lambda_0 (1 + z), instead of the rest one, a = Gamma lambda_0 / (4 pi b) (V1 `_voigt_par_convert`, the mocks): the Lorentzian part is (1 + z) times too strong, 4 times at z = 3. The peak optical depth comes out lower (2% on SiIV with b = 1.5 km/s, logN 13.8 at z = 3.1) and the damping wings too strong: a column measured from the wings (DLA, LLS) comes out low by ~log(1 + z), 0.6 dex at z = 3. In the fitter (`_voigt_optical_depth`, `_voigt_optical_depth_vectorized`) and in the profiles of the System Inspector; there since the first V2 version. Found in #526 (PR 3): the residual of the fitter model against the exact flux of the mocks, 1.3e-3 to 1.8e-3 on CIV and Ly-a and 3.6% on a saturated SiIV line, drops to 3e-5 with the rest wavelength | numerical, wrong results | `test_damping.py` (fit of a damped Ly-a on mocks: logN 20.3 recovered within 0.01 dex, 19.75 at z = 2.5 and 19.64 at z = 3.5 before), `test_parity_v1_v2.py` (damped Ly-a, equivalent width 1.6 times that of V1 before), `test_lsf.py` (tolerance 1e-4) | fixed in phase 4 (separate PR after #526 PR 3, agreed with Guido): a with the rest wavelength, in the fitter and in the System Inspector. Benchmark on mocks (no damped systems there): spread of the system column 0.102 → 0.095 dex at R 8900, 0.046 → 0.041 at R 45000, the rest unchanged |
| F62 | automatic fits (forest) | The pixels of a fit are chosen at the start, around the lines of the active components, from their initial logN (`prepare_fit_context`, `_find_feature_limits`): if logN grows during the fit, the damping wings that fall outside cost nothing, and logN can run away to the bound (23). The groups fitted later see the component as a fixed background. On J0103-1305 `lya_auto` gives components with logN 21-22.4 (b 1 or ~50 km/s) that black out the model where the data are at 0.3-0.8 (540-590 nm) and lower it by 2-13% on the red side. Already present before F61 (logN 20.3-21 with b = 1 km/s, model black at 540-550 nm), made worse by F61: the damping wings, (1 + z) too strong, penalised the growth of logN inside the window. The mocks of the benchmark have no saturated forest lines strong enough to show it | numerical, wrong results | `test_snapshots_full_chain.py` (J0103-1305, characterization) | fixed in phase 4 (after F61): upper bound on logN in the automatic pipelines, `max_logN` = 19 (Guido's decision; the Lyman-limit systems stay within reach, damped systems are fitted by hand; `refit_all` and the manual fits keep 23). Tried first and dropped: extending the fit to the pixels reached by the wings at the fitted logN (within 1000 km/s). In the pipelines the absorption around a candidate is not modelled yet, so a wider window rewards the wings: more runaways (up to logN 22.97) and a slower chain. Tried also the bound at 17.5. On the J0103-1305 chain, forest 540-600 nm: chi2 per pixel 98 before F61, 305 after, 31 with the bound at 19 (28 at 17.5); over-absorbed pixels (model below the data by more than 5 sigma) 28%, 68%, 3.4% (2.3%); red side as before F61. The components at the bound are what a richer pipeline must handle (phase 7) |
| F63 | GUI (Set session properties) | The dialog fills each parameter whose default is `_current_` with the attribute of the same name of the spectrum data (`recipe_dialog.py`), and with "0.0" when there is none: the session name and the object have no such attribute (and the resolution no longer has one after #526), so the dialog shows "0.0" in them. Pressing Run without editing renames the session "0.0" and writes "0.0" into the OBJECT keyword of the header; the resolution is kept (a value of 0 is ignored), z_em is shown correctly. Found while scripting the screenshots of the tutorials (October 2026) | data loss (name, object) | — | fixed in phase 4 (Guido's agreement): the dialog shows the current value of the session (`RecipeDialog._current_value`: name, OBJECT, z_em, the resolution, left empty with the median as a hint when it varies), and an empty field keeps it (`'_current_'`). The same "0.0" was in z_ref of *Convert x axis* (velocities relative to z = 0 instead of z_em; the note of #573 that the dialog used z_em was wrong) and z_target of *Trim to common coverage*. `test_recipe_dialog_current.py` |
| F64 | GUI (plot controls) | The stylesheet of the checkboxes of the plot controls set `background-color: transparent`, which with the non-native Qt styles (Fusion: Linux, Windows) also clears the frame of the indicator: unchecked boxes were invisible, checked ones showed only the tick. macOS (native style) was not affected. Found while scripting the screenshots | GUI | `test_unchecked_boxes_of_the_plot_controls_are_visible` (regression) | fixed: the background colour is no longer set |
| F65 | `forest_auto` | The recipe takes `min_b` (schema and docstring: "Minimum b parameter (km/s)") but never uses it: no filter of narrow lines, unlike `lya_auto`. Found in #517 (PR 2) | behaviour | — | fixed in phase 4 (#517, PR 3, Guido's decision): the filter of `lya_auto` (`_remove_narrow`) applied to the series `Lyman`; `test_pipeline_cleanup.py` |

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
phase: F18, F31 (phase 2); F6, F7 (phases 3/5); F12, F13 (phase 5); F28
(save the Bayesian samples as a separate file in the archive) when needed.

### Rough estimates
10–20 work sessions, 2–5 calendar weeks; the bottleneck is human review and
the availability of test data.

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
- **Full server, after phase 6-bis:** all the recipes, structured outputs,
  natural-language log, and a headless Docker image to distribute it. An
  agent working on its own needs correct and measured results (phase 4),
  clean error handling (phase 6) and a stable recipe API (phase 6-bis).
  The tools are generated from the recipe registry (the schemas the GUI
  already builds its menus and dialogs from), not written by hand: a
  recipe added or changed later reaches the server with no server code,
  and a test checks that every recipe is a tool or is explicitly
  excluded (interactive recipes such as `update_from_knots`). This is
  why the server can come before phase 7: its new and improved recipes
  reach the agents as they reach the GUI.

### 4.1 MCP pilot: scope
Proposed in October 2026, to be agreed with Guido before the code.

**What MCP is.** The Model Context Protocol is an open standard for
connecting an LLM application (Claude Code, Claude Desktop, other agents:
the *client*) to external programs (the *servers*). A server declares a
list of *tools*, each with a name, a description and a JSON schema of its
parameters and of its output; the client shows them to the model, which
decides when to call them; the client runs the call and gives the result
back to the model. For a local program like Astrocook the client launches
the server as a subprocess and talks to it over stdin/stdout (`stdio`
transport): no network, no account, the data stay on the machine. The
model never sees Astrocook's code, only the tool descriptions and what the
tools return: those are the whole interface, and what the pilot is about.

**Tools (4 recipes and a save).** Each tool calls one existing recipe of
`SessionV2` and returns a structured summary, not the spectrum. The
sessions live in the server process, under an identifier the tools return
and take. A tool is named as its recipe, prefixed with its menu
(`<menu>_<recipe>`, agreed with Guido): a single vocabulary in the tools,
the log, the archives and the docs, and the full server can generate its
tools from the recipe schemas (`*_RECIPES_SCHEMAS`: names, descriptions,
parameters). An underscore rather than a dot, since some clients accept
only letters, digits, `_` and `-`. A recipe name unclear to a model is a
finding on the recipe, fixed at the source (phase 6-bis), not hidden by
the server. Opening and saving are methods of `SessionV2`, not recipes:
`session_open` and `session_save`.

| Tool | Recipe | Main parameters | Structured output |
|---|---|---|---|
| `session_open` | `SessionV2.open_new` (+ `edit.set_properties` if `z_em` or `resol` are given) | path, format (`auto`), `z_em`, `resol` | session id; number of pixels; wavelength range and unit; `z_em`; resolution; continuum and components present; warnings |
| `continuum_estimate_auto` | `continuum.estimate_auto` | session id; `kappa`, `smooth_fwhm` (others at the recipe default) | fraction of pixels flagged as absorbed; median S/N; warnings (e.g. `z_em` = 0) |
| `absorbers_identify_lines` | `absorbers.identify_lines` | session id; multiplets, `score_threshold` | list of candidates: region, wavelength range, series, z, score; no component is created |
| `absorbers_doublets_auto` | `absorbers.doublets_auto` | session id; multiplets, `score_threshold`, `max_b` | table of the components: series, z, logN, b with errors, reduced chi2; flags: NaN errors (F20), at a bound (F21, #517); warnings |
| `session_save` | `SessionV2.save` | session id; path (`.acs2`) | path; number of log entries |

The fit is `doublets_auto` (identification and fit on the red side in one
recipe, the most tested path: full-chain snapshots) rather than a fit of
the candidates chosen by the agent, because that path goes through
`populate_from_identification`, which takes the identifications as JSON
strings stored in the spectrum metadata: a first finding on structured
outputs. `lya_auto` and the forest are left out of the pilot (slow, and
the forest is where the continuum is weakest, #465). `save_session` is
not a recipe, but without it the human cannot open the agent's work in
the GUI (the review bench of item 4).

**A thin layer.** The server translates, it does not redesign: the
recipes are good starting points, but not robust and not validated
scientifically, and a server that corrected them on its own would
diverge from the GUI and the scripts (as the double defaults of F40, on a
larger scale). Their weak points are exposed and recorded, and fixed in
the recipes themselves (phases 4, 6-bis and 7, section 3).

**Module.** A subpackage `astrocook/mcp_server/` in this repository, with
an optional extra `[mcp]` and an entry point `astrocook-mcp`: same
version and same tests as the rest, nothing to install separately. The
existing code is not changed: what the pilot needs and the recipes do not
offer goes into the findings register, not into a quick fix.
- `server.py`: the server and the tool definitions (names, descriptions,
  typed parameters);
- `sessions.py`: the sessions in memory, one `SessionHistory` per
  identifier, so that the log is kept as the GUI keeps it (the log entries
  are written by the GUI, not by the recipes: a headless caller must add
  them itself);
- `outputs.py`: the output models (pydantic, from which the SDK derives
  the output schemas);
- `errors.py`: runs a recipe while collecting its log records; a recipe
  returning 0 becomes an MCP error carrying the logged message, the
  warnings go into the output.
Tests (`tests/test_mcp_server.py`) call the server in process through the
SDK client, on the mock and on J0103-1305; skipped if `mcp` is not
installed, the `[dev]` extra includes it so that CI runs them.

**SDK.** The official Python SDK (`mcp` on PyPI, MIT license, version 2.3
in October 2026, Python >= 3.10), with its high-level `MCPServer` class
(called `FastMCP` in the 1.x series): a decorated function is a tool, the
parameter types give the input schema and the return type the output
schema. Pinned to `mcp>=2.3,<3`, since 2.0 changed the API. Two
behaviours to keep in mind, checked on 2.3: an unexpected exception
reaches the model only as "Error executing tool <name>", without the
message, so the expected failures must be raised as `ToolError`; and
anything written to stdout breaks the `stdio` transport (no `print` in
the non-GUI code today; the logging goes to stderr).

**Installation.** From GitHub, as Astrocook itself until PyPI:

    pip install "astrocook[mcp] @ git+https://github.com/DAS-OATs/astrocook"

then the server is registered in the client, e.g. for Claude Code
`claude mcp add astrocook -- astrocook-mcp` (Claude Desktop: the same
command in its JSON configuration). PySide6 is installed too, as a
dependency of Astrocook, though the server does not use the GUI; a
lighter headless installation is for the Docker image of the full server.
The docs get a short page "Astrocook from an agent".

**What we want to learn.** The pilot ends with a short report (in this
section), written after using the server from an agent on the mock and on
J0103-1305, with answers to:
1. *Outputs (for phase 6-bis and the full server):* which fields the agent
   needs to decide (quality flags, chi2, components at the bounds, F21,
   #517) and which the recipes do not return today (identifications in
   the metadata, results only in the session); whether a common output
   format for the recipes is worth adding to the API.
2. *Errors (for phases 6 and 6-bis):* how many failures reach the agent with a
   useful message; which recipes return 0 or log only a warning (e.g.
   `z_em` = 0 before the continuum) where an exception or a structured
   warning would be better; a list of cases for #566 and for the broad
   `except` clauses.
3. *Parameters (for #573):* recipe parameters are strings
   (`kappa: str = '2.0'`) and have two sets of defaults (F40); which
   types and defaults the tools need.
4. *Log:* whether the log of a session built by an agent is readable and
   replayable in the GUI, and what is missing for a natural-language log
   (item 3).
5. *Usability:* whether the descriptions are enough for the model to
   chain the tools correctly, and how long the calls take (`doublets_auto`
   on J0103-1305).
The answers become findings, issues, or items of phases 4 and 6, agreed
with Guido as usual; none changes the recipes within the pilot.

**PRs, in order (one at a time, each merged before the next):**
1. This scope in `ROADMAP.md` (no code).
2. Skeleton: `astrocook/mcp_server/`, `[mcp]` extra, `astrocook-mcp`
   entry point, the sessions, the error capture, `session_open` and
   `session_save`, with tests (the `.acs2` saved by the server reopens in
   the GUI with its log).
3. `continuum_estimate_auto` and `absorbers_identify_lines`, with tests.
4. `absorbers_doublets_auto`, with tests on the mock (components recovered within
   their errors) and on J0103-1305 (slow test).
5. Docs page (installation, registration in the clients, an example
   session) and a test of the server launched as a subprocess over stdio.
6. Report: a session with an agent on the mock and on J0103-1305, the
   answers above in this section, the new findings; `CLAUDE.md` and the
   progress table updated. The report is the starting point of the
   recipe API (phase 6-bis). Then phase 4.

### 4.2 MCP pilot: report
October 2026, after PRs #595-#599. Two trials with fresh agents, then
Guido's own with Claude Code on his computer (below).

**Trials.** Two fresh agents (Claude, with no knowledge of the code, the
tests or the true parameters, told not to read them) used the server
through a minimal MCP client, seeing only what a client shows a model:
the server instructions, the tool descriptions and schemas, the results.
Each got a request in plain words and answered it, then wrote a critique
of the server. Transcripts and reports are summarized here.
- *Full-chain mock* ("find the metal systems on the red side, fit them,
  say which components are reliable"). 10 calls, none failed, none over
  1 s. The agent found the 3 CIV systems and judged 6 components: the 2
  spurious ones of F49 "not reliable" (0.7 km/s from the main component,
  below the resolution element, dlogN 0.6 and 1.1), 3 "reliable", and
  one real component (z = 1.83000, logN 13.4, b 8) only "tentative". It
  explored by itself (lower `score_threshold`, wider `z_window_kms`).
  The mock archive also carries the true continuum (`cont` from the
  generator), which the agent used in its first session: future trials
  should remove it.
- *J0103-1305* ("estimate the continuum, tell me the metal systems, fit
  the doublets; I do not want to lose my work"). 7 calls, none failed;
  `absorbers_doublets_auto` took 168 s. The agent worked on a second copy
  of the session to keep the original, set z_em and R (the archive had
  neither), noticed and reported that the doublet fit removed 190 of the
  764 components (the count is all the output gives), and singled out the
  bad fits: CIV at z = 3.927 (4 components at logN = 10 and b = 1 km/s,
  dz up to 131, dlogN up to 3800), z = 3.885 and 3.28 (reduced chi2 up to
  135), several MgII with b of 1.8-2.2 km/s. It suggested telluric
  absorption (O2 A band at 760 nm) for the worst ones, which no output
  could confirm: to be checked by Guido.

**Answers to the questions of section 4.1.**
1. *Outputs.* The summaries were enough to chain the tools and to judge
   most components; the agents asked for:
   - quality beyond the bounds: today a component with dlogN of 1 dex or
     db > b has no flag, so `n_flagged = 0` reads as "all good" (21 of
     the 101 new components on J0103-1305); flags for large relative
     errors, high chi2, components closer than the resolution element,
     and a verdict per system;
   - structure: the fit group of a component (now inferred from equal
     chi2), its system, the transitions fitted, the velocity offsets in
     km/s, the link between the systems of `identify_lines` and those of
     `doublets_auto`, the score of each candidate system;
   - provenance: what a session held when opened (an existing continuum,
     components with an empty log, z_em and R missing), and what a recipe
     replaced or removed (the 190 components, the continuum, the mask:
     F50), listed, not only counted;
   - lists that fit a model's context: 833 Ly-a candidates crowd out the
     metals in `absorbers_identify_lines` even when Ly-a is not asked for
     (F51); filters by series or range.
   A common result format for the recipes (summary, what changed, flags,
   warnings) is worth adding to the API (phase 6-bis): the server built
   it here by comparing states, which a recipe can do better.
2. *Errors.* No call failed in the trials; in the tests, every failure
   reached the model with Astrocook's message (missing or unreadable
   file, z_em not set, unknown session), thanks to the capture of the
   logged errors (`errors.py`). The recipes also fail silently in ways no
   message reports: unphysical values accepted (F46), absurd
   uncertainties stored as results (F52), a mask rebuilt with defaults
   (F50), a mask built when missing (INFO only). The only warning seen
   was noise (astropy's `VerifyWarning` on the `resol_fwhm` FITS keyword,
   at every save). For phase 6: structured warnings with a code, the
   cosmetic library warnings filtered.
3. *Parameters.* Typed parameters and the schema defaults worked; the
   descriptions of the schemas are too short for a model (`z_window_kms`:
   half or full width? `fudge`: what?, the "Fit window"), sometimes wrong
   (F48), and some parameters must not be offered (column names, GUI
   plots): the schema needs a flag for them (phase 6-bis).
4. *Log.* The log of a session built by the agents is the one the GUI
   would write (`set_properties`, `estimate_auto`, `identify_lines`,
   `doublets_auto`, with their parameters as strings) and replays: the
   entries of the mock trial, run again on the original archive, give
   the same 6 components. For a natural-language log (section 4, item 3)
   the entries lack what the recipe *did* (components added and removed,
   flags): the same information as the structured outputs.
5. *Usability.* Both agents chained the tools correctly at the first
   attempt, from the descriptions alone, and explored parameters by
   themselves. They could not act on their judgement: no tool to remove
   or refit a component, compare models, list components, duplicate or
   close a session. `absorbers_doublets_auto` is the one dangerous tool
   (it removes components); the agent protected the user's work only by
   opening a second copy. Calls of minutes need progress notifications.

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

**Guido's trial.** Claude Code on Guido's computer, with the server
registered as in the docs, on an ESPRESSO spectrum of J0959+1312 (a raw
FITS, z_em not in the file), in a free conversation in Italian: fit the
continuum, show it, open it in the GUI, why is the forest wrong, try the
template, fit the absorbers, look for the Ly-a of the CIV systems, fit the
forest. Guido's verdict: far beyond expectations, "a different approach
to the analysis". What it showed:
- *The agent works around the server.* With shell commands and Python on
  the FITS file, outside the server, it estimated z_em from the Ly-a
  emission peak (4.06, said to be approximate), inferred R from the
  header (VLT-U1, so ESPRESSO HR, R = 140000, asking for confirmation),
  measured the continuum against the flux envelope, measured the Ly-a and
  Ly-b of the CIV systems on the normalized flux (all saturated, only
  lower limits on N(HI)), made an interactive plot of spectrum and
  continua, opened the archives in the GUI (`open -a Astrocook`), and
  offered to write its own Voigt fit of the forest. Powerful, but none of
  this is in the session log or reproducible from it: the server needs
  tools for what the agent did by hand (reading the spectrum in a range,
  a preview, opening the GUI, the forest fit, editing components), and
  the log needs the agent's decisions and their reasons (z_em estimated
  from the Ly-a peak, R from the header; section 4, item 3).
- *It diagnosed the recipes correctly.* The automatic continuum is a
  straight interpolation over the forest because no pixel there is left
  unabsorbed (F53); the template gives it a shape but sits below the flux
  peaks at 565-590 nm; the doublet fit gives spurious systems in the O2 A
  band (CIV z = 3.928 and MgII z = 1.733 fitted together at 763 nm, as
  CIV z = 3.927 in J0103-1305: F54) and a SiIV at z = 4.33 > z_em (F55).
- *It misread one cause.* `absorbers_identify_lines` after the doublet
  fit found no forest region, which the agent explained by blended
  absorption; the cause is F50 (the fit cleared the forest in the stored
  mask). Outputs that say what a recipe changed (answer 1 above) would
  have shown it.
- It kept the files apart (a new archive for each step, the original
  untouched) and stated its assumptions and what it could not do.
- The docs did not say how to install Claude Code (fixed in this PR).

**Outcome.** The pilot reached its goal: an agent can carry out an
analysis with Astrocook from a request in plain words, and the trials
show what the full server and the recipes need (answers above). The
findings F46-F55 go to phases 4, 6-bis and 7; the lessons for the
server to section 4 and phase 6-bis.

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
First the robust components (phase 7: emission model and identification,
measured on the benchmark of phase 4), then the loop.
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
| 1-bis – Findings and urgent issues (section 3.3) | done (October 2026): F1, F2, F5, F10 (minimal), F11, F16, F19, F20, F23, F24, F25, F26, F27, F29, F32, F33, #480, #442 fixed; F3 accepted; deferred to issues: #526 (F8, F22, F10, F9, F14), #534, #465 (F30, F17), #517 (F21), #518 (F4); F15 waits for data |
| 2 – Cleanup and packaging | done (October 2026), release 2.0.3; LICENSE pending (#549) |
| 2-bis – Leftovers of phase 2 (F34, Ruff, Sphinx, `[bayes]`) | done (October 2026): F34 fixed; Ruff and Sphinx warnings cleaned up, full `ruff check .` in CI; `[bayes]` extra; found F36, F37 (#559, #560, phase 4) |
| 3 – Decoupling from legacy | done (October 2026): `legacy/` removed, native reader of the V1 archives, `.acs` saving removed; F42, F43, F44 fixed; found F38-F41 (#568, #573, #575); docs published per version |
| 3-bis – Small fixes before 2.1.0 (#575, F6, F13, #559, F45, #560) | done (October 2026): F41 (#575), F6, F13, F36 (#559), F45, F37 (#560) fixed; note on F8 |
| MCP pilot (sections 4.1, 4.2) | done (October 2026): scope (#595); skeleton with `session_open`, `session_save` (#596); `continuum_estimate_auto`, `absorbers_identify_lines` (#597); `absorbers_doublets_auto` (#598); docs page and stdio test (#599); report (section 4.2, #600): two agent trials and Guido's; found F46-F55 |
| 4 – Scientific correctness (validation on mocks; #526, #517, #465, #518, #568, #573) | done (October 2026): plan agreed (#601); validation: mocks extended (#602), benchmark (#603), baseline (#604); found F56-F59; #573 done (#605: one default per parameter); #518 done (#606: fit window in the System Inspector); #568: wrong atomic values (#607), AlII, [NeIII], [NeV], TiII (#608), single table with sources and newer values, F39 (#609); F49: one candidate per doublet (#611); F55: candidates within 3000 km/s of z_em (#612); F50: the pipelines keep the mask of the user (#613); #465: robust automatic continuum (spikes, masked running mean, local S/N) (#614); #526 PR 1: resolution as a single column (#615); PR 2: the GUI on the same resolver (#616); PR 3: exact convolution with the LSF of each pixel (#617); F61: damping parameter (1 + z) too large (#619); F62: `max_logN` = 19 in the automatic pipelines (#624); PR 4: legacy sources no longer written, System Inspector profiles as the fitter (#625); #623: refit of sessions with pixel edges without a unit (regression of #617, #626); F14: blend threshold of the grouping from the instrumental FWHM (#618); #517 with F52 partly: flags of the fit (#627), cleanup of the pipelines (#628), GUI cues and F65 (#629); F63 (#631); result against the baseline in "Result of the phase" |
| 5 – Splitting large files (and GUI performance) | in progress (October 2026): order of the phase agreed, F12 moved to 6-bis; GUI benchmark (`benchmarks/bench_gui.py`) |
| 5-bis – GUI look and ergonomics | not started |
| 6 – Error handling | not started |
| 6-bis – Recipe API | not started |
| 7 – Robust components | not started |
