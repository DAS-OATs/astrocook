# Astrocook – Development roadmap, history

Archive of the concluded parts of `ROADMAP.md`, moved verbatim: the detail of phases 0 to 6-quater, the closed findings and issues, the MCP pilot, the trials of the full server, the first state of v2. Do not read it whole (about 4600 lines): search it with `grep` for a finding (`F52`), an issue (`#526`), a PR or a section title.

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


---

## 3. v2 consolidation (concluded phases)

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

Fixes of the System Inspector, each measured with the benchmark before and
after on the same machine:
1. *The table* (first fix): the model answered every role Qt asks for,
   and the columns were fitted with `ResizeToContents`, so Qt measured
   every cell (nine roles each) three times while the inspector opened
   (138,000 calls of `SystemTableModel.data`). Now the roles the model
   does not use are answered at once, and the columns are fitted once
   (`resizeColumnsToContents`) when the rows shown change (new session,
   group view, sorting, filters): 55,700 calls, of which 13,500 are the
   sorting. Same widths (characterization test) and the same images of the
   tutorials; the columns can now also be resized by hand. Before → after:
   `inspector_first` 3.74 → 1.34 s, `inspector_again` 2.42 → 0.89 s,
   `refit_gui` 1.73 → 1.33 s; `main_plot` and `refit_recipe` unchanged.
2. *The main plot* (second fix): the ticks of the systems were converted
   to the x unit one transition at a time (an astropy quantity each, 954
   conversions on J0103-1305, ~65% of a redraw); now the observed
   wavelengths are collected and converted at once. Same positions, bit for
   bit (snapshots `main_plot_system_ticks_nm`, `_kms`). Before → after:
   `main_plot` 0.19 → 0.13 s; `inspector_first` 1.42 → 1.33 s (within the
   noise). The 1.2 s measured in phase 4 for this redraw is not seen on
   this machine (0.24 s inside the opening of the inspector, before the
   fix). Found F66 (ticks in velocity always at z_ref = 0). What is left
   of a redraw is matplotlib (ticks, text layout): the third fix. In the
   opening of the inspector the largest part is now the column fit of the
   table (~0.7 s, every row measured once).
3. *Matplotlib* (third fix): the best position of the legends
   (`loc="best"`), suspected in phase 4, costs ~0.01 s here and is left as
   it is. The cost was in the velocity plot of the inspector: on opening,
   it was rebuilt twice (the selection, then the first resize of the
   widget) and its canvas drawn three times (twice by `_update_plot`, once
   by the resize), ~0.1 s per draw. `_update_plot` now asks for
   `draw_idle`, so the draws in a row are coalesced into one: same pixels
   once the events are processed (characterization test), same images of
   the tutorials. Before → after: `inspector_first` 1.43 → 1.19 s;
   `inspector_again` (one rebuild and one draw already) and the other
   cases unchanged within the noise. What is left of an opening is mostly
   the table: its column fit (~0.7 s on the first opening, ~0.8 s of a
   second one, which refills the table).
4. *The column fit of the table* (fourth fix, agreed with Guido). The
   cause was not the fit itself but the cost of each call of
   `SystemTableModel.data`, ~29 us: with PySide6 every access to a Qt
   enum through its namespace (`Qt.DisplayRole`, `Qt.AlignCenter`...)
   costs ~4 us, and the model made several per call. The enums of the
   table (model, delegate, proxy) are now read once, at module level:
   ~2 us per call, and the column fit of Qt 0.52 → 0.16 s. A fit of our
   own (the text widths with the font metrics, Qt's measure only for the
   widest cells: the same widths, 0.06 s) was tried and left out: 0.1 s
   for ~40 lines of code, and slower than Qt's above 1000 rows, where Qt
   measures a sample. The two other ways considered (only the visible
   rows; widths from the text alone) were not needed. Same values of the
   model (characterization tests), same images of the tutorials. Before
   → after: `inspector_first` 1.22 → 0.54 s, `inspector_again` 0.93 →
   0.37 s, `refit_gui` 1.20 → 0.64 s. The same cost of the Qt enums may
   be found in other code called per pixel or per cell (to check with
   the split of the GUI).

5. *Selecting a row of the table* (fifth fix, agreed with Guido after his
   trial of the four fixes: the inspector opened faster, but browsing the
   table was still slow; his benchmark on a MacBook Air after #638:
   `inspector_first` 0.45 s, `inspector_again` 0.33 s, `main_plot`
   0.07 s, `refit_gui` 0.54 s). Measured with a new case of the benchmark,
   `select_row`, with the window at the size of a laptop screen (the cost
   of a draw grows with the size of the plots): ~0.3 s per selection, two
   full matplotlib draws (the main plot, rebuilt to highlight the
   selection, and the velocity plot). The main plot now redraws only what
   changes: the guide lines, the markers and the legend are drawn over the
   pixels of the rest of the axes, kept from the first selection after
   each full draw (blitting, as the cursor of z already did). A full draw
   is still made when the "Selected" entry of the legend appears or
   disappears, while the continuum is edited, and whenever anything else
   changes (any full draw discards the kept pixels). Same pixels as a full
   rebuild after a sequence of selections (characterization test, checked
   to fail if the markers are not drawn), same images of the tutorials.
   Before → after: `select_row` 0.32 → 0.22 s; the other cases unchanged
   within the noise. What is left of a selection is the velocity plot
   (rebuilt for the new component, ~0.15 s). Seen along the way: after a
   refit the main plot is drawn twice (already before this fix).
   Guido's benchmark on his MacBook Air with this fix: `select_row`
   0.135 s offscreen, 0.148 s with the real screen (`QT_QPA_PLATFORM=cocoa`:
   the windows are shown and the times include the drawing on the screen;
   `inspector_first` 0.67 s, `inspector_again` 0.41 s, `refit_gui`
   0.60 s). What is left (agreed with Guido, not urgent: with the split of
   the GUI): the draw of the velocity plot (~half of a selection: its
   constrained layout measures the ticks and the text several times per
   draw; a layout computed once per number of panels may move the margins
   by a few pixels), its rebuild (the panels updated instead of recreated
   when their number does not change), and the updates of the plots only
   for the last row reached while an arrow key is held.

Result of the four fixes on the benchmark (same machine, against the
numbers before the phase): `inspector_first` 3.96 → 0.54 s,
`inspector_again` 2.25 → 0.37 s, `main_plot` 0.21 → 0.14 s, `refit_gui`
1.78 → 0.64 s.

**Split of `SpectrumV2`** (plan agreed with Guido, October 2026):
`core/spectrum.py` (2826 lines) keeps the core (data, columns, properties,
conversions of the axes, `TableAdapterV2`); the other methods go into
thematic mixin classes that `SpectrumV2` inherits, so that its interface
stays the same (names, signatures, docstrings, the references of the
docs). Five modules, one PR each: `spectrum_continuum.py`
(`ContinuumMixin`), `spectrum_lines.py` (detection and identification of
lines, doublets), `spectrum_flux.py` (smoothing, running std, scaling,
calibration), `spectrum_combine.py` (`stitch`, `coadd`, `rebin`,
resampling), `spectrum_expressions.py` (expressions, masks, `split`). The
methods are moved verbatim (reviewed with `git diff --color-moved`); the
only change is `SpectrumV2(...)` → `self.__class__(...)` for the results
(no subclass exists; `coadd` already did so). Each mixin lists in its
docstring what it uses of `SpectrumV2`. `tests/test_spectrum_structure.py`
checks that the interface is the one of 2.2.0 and that no name is defined
twice (with multiple inheritance the first class would win silently).
Before each move, characterization tests for the paths of the moved
methods that no test reached. Done: continuum (`find_absorbed`,
`fit_continuum`, `fit_powerlaw`, `update_continuum_from_knots`; the last
had no test; found F67, a crash with a model and no continuum; #640);
lines (`detect_absorptions`, `identify_lines`, `get_velocity_bounds`,
`detect_doublet_z`, with `Z_EM_MARGIN_KMS`; three of them had no test;
found F68, the Ly-a/Ly-b pairs never identified; phase 7, Guido's
decision). Before the flux methods: *Flux > Calibrate* failed at every
call (`np.trapz`, removed from recent NumPy, in `core/photometry.py`; no
test reached it), fixed with SciPy's `trapezoid` in a PR of its own, with
the first tests of the photometry; found F69 (the "flux" columns of the
flux methods) and F70 (a two-band calibration ~0.1 mag from its targets;
to discuss, Guido) (#642). Flux (`scale_flux`, `calculate_running_std`,
`smooth`, `smooth_column`, `calibrate`; the calibration had no test until
#642, the others were reached only through the recipes; #643). Combination
(`stitch`, `rebin`, `equalize_to_reference`, `coadd`, `resample_on_grid`,
`get_resampled_column`; found F71, an error message that fails, and F72,
the masks dropped by a rebin; #644). Expressions, masks and regions
(`apply_expression`, `mask_expression`, `prepare_expression_context`,
`calc_intersection_expression`, `split`, `mask_region`,
`get_region_bounds`, with `_NE_GLOBAL_DICT` and `REGION_PRESETS`). The
split is complete: `spectrum.py` 2826 → 774 lines (the core: data,
columns, properties, conversions of the axes, `TableAdapterV2`), and the
five mixins `spectrum_continuum.py` (341), `spectrum_lines.py` (486),
`spectrum_flux.py` (338), `spectrum_combine.py` (553),
`spectrum_expressions.py` (485). Before the split the tests reached 78%
of `spectrum.py`; every moved method now has characterization tests
(`tests/test_spectrum_*.py`). Found along the way: F67-F72, and the
calibration that failed with recent NumPy (#642).

**Parallel tests** (agreed with Guido, October 2026, before the split of
the GUI): `pytest-xdist` in the `test` extra, `pytest -n auto` in the CI.
Measured on 4 cores (as the GitHub runners): the fast suite 440 → 161 s,
the same tests passing and the same coverage (67%). The slow job stays
serial (one test of ~10 minutes dominates it, and each process would
rebuild the costly fixtures). During the split of the GUI the targeted
tests run locally, the full suite in the CI.

**Split of the GUI** (plan agreed with Guido, October 2026): F7 first;
then `system_inspector.py` (the table, the velocity plot, then the two
deferred speed-ups of the velocity plot), `pyside_plot.py` (helpers,
canvas, toolbar and data inspector in their modules; the continuum editor
and the mouse handlers in mixins; `plot_spectrum`, ~1000 lines in one
method, broken into steps, a real refactoring checked pixel by pixel on a
set of view settings), `main_window.py` (accessory widgets and dialogs,
menus, panels, sessions, recipes and scripts), then F66 and the closing of
the phase. As for `SpectrumV2`: characterization tests before each move,
verbatim moves into modules and mixins; checked by pixel comparisons, all
the images of the tutorials and `bench_gui.py`. Done: F7 (#647); the
table of the System Inspector in `gui/inspector_table.py` (columns, cues
of the flags, model, delegate and proxy, 410 lines; `system_inspector.py`
2239 → 1884; editing a cell, the tooltip of a link by expression and the
filters and sorting of the proxy had no test; #648); the velocity plot in
`gui/velocity_plot.py` (`VelocityPlotWidget` and `calc_voigt_profile`, 820
lines; `system_inspector.py` 1884 → 1085; the tooltip of the ticks, the
profile under the cursor and the menu of the right click had no test;
#649, with a test of the GUI that waited a fixed time after a recipe and
failed on a loaded machine). Speed-ups of the velocity plot (the two
deferred with Guido): the rest wavelengths converted once per rebuild
(the ticks of every panel converted ~1000 of them with astropy, half of
a rebuild), and the constrained layout of matplotlib kept by (panels,
residuals, size of the figure, x limits) and reused, instead of laying
out at every draw (half of a draw, also at each move of the mouse, which
redraws the cursor profile). The same pixels as a rebuild laid out from
scratch (test checked to fail if the stored layout is not applied); no
margin moves, unlike what was expected. `select_row` 0.23-0.26 → 0.18 s,
`inspector_again` 0.43-0.55 → 0.36 s. The panels are still rebuilt at
each selection (~0.06 s): updating them in place would mean rewriting
`_update_plot`, not worth it for now (#650). The parts of the main
plot out of `pyside_plot.py`: `plot_style.py` (style, colour cycle,
decimation), `plot_canvas.py` (`MatplotlibCanvas`), `plot_toolbar.py`
(`AstrocookToolbar`, also used by the velocity and corner plots),
`data_inspector.py` (`DataInspectorDialog`); `pyside_plot.py` 2884 →
2255 lines; the decimation, the cursor, the tooltip after a pause, the
toolbar and the data inspector had no test; found F73 (the selection of
a region with the mouse, dead code).
The continuum editor of the main plot (knots spaced in velocity, the
draft curve patched locally around a moved knot, the stride, reset, save
and cancel) in a mixin, `continuum_editor.py` (`ContinuumEditorMixin`),
moved verbatim; `pyside_plot.py` 2255 → 1932 lines. The editor had no
test: characterization tests in `tests/test_continuum_editor.py`, down to
an edit saved from the window; `tests/test_gui_structure.py` checks that
the classes split into mixins keep their interface and define no name
twice. Found F74 (a missing method in the save path of the editor).
The mouse on the main plot (the redshift cursor, the tooltip after a
pause, the menu of the right click, the knots moved, added and removed
with the mouse, the data inspector opened from the menu) in a mixin,
`plot_mouse.py` (`PlotMouseMixin`), moved verbatim; `pyside_plot.py` 1932
→ 1330 lines. The tooltip, the menu and the knots moved with the mouse
had no test: characterization tests in `tests/test_main_plot_mouse.py`
(checked against eight mutations of the handlers, all caught).
`plot_spectrum` (~1000 lines) split into its parts, a real refactoring:
the state of the plot (caches, zoom kept, limit callbacks, legend corner,
view settings in a `_PlotView`), the arrays (`_PlotArrays`: raw,
normalized, SNR, cut to the view; the 17 cache attributes replaced by a
dictionary), one method per layer (other sessions, error, flux,
continuum and model, auxiliary column, systems, cursor), legend, axes
(units, coordinates under the mouse, rest-frame axis, emission lines)
and limits; `pyside_plot.py` 1330 → 1263 lines, `plot_spectrum` ~100.
Checked pixel by pixel and artist by artist (data, colours, styles,
labels, limits, legend, ticks, the text of the coordinates) on a
sequence of 60 view settings (normalized, SNR with two error columns,
log axes, three x units, three kinds of auxiliary column, highlights,
cursor, iso-magnitude grid, zoom, scatter style, other sessions, no
session) and of 11 on a spectrum of 2 million pixels (decimation): all
identical, and the same messages in the log. Characterization tests in
`tests/test_main_plot_layers.py`. GUI benchmark unchanged (`main_plot`
0.15 → 0.14 s). Found F75 (the default error column of the SNR view
looked up as a column) and F76 (decimation of long spectra only at the
first plot).
Split of `main_window.py` (4831 lines), in five steps: the accessory
widgets and dialogs, then mixins for the menus, the panels, the sessions,
and the recipes and scripts; `tests/test_gui_structure.py` checks the
interface of `MainWindowV2` too. First step: the floating cards (fader,
collapsible sections, wrapped label) in `floating_cards.py`, the progress
of a recipe (dialog and log handlers) in `recipe_progress.py`, the About
dialog in `about_dialog.py`, moved verbatim; `main_window.py` 4831 →
4470 lines. The progress dialog, the log handlers and the About dialog
had no test: `tests/test_gui_dialogs.py`.
Second step: the menus (the menu bar and its entries, undo and redo,
the windows of the View menu, export and import of ASCII files, the
documentation) in a mixin, `main_menus.py` (`MainMenusMixin`), moved
verbatim; `main_window.py` 4470 → 3984 lines. Characterization tests in
`tests/test_main_menus.py` (the entries of each menu and their order,
the shortcuts, the actions enabled with a session, the recipe dialog of
each entry, undo and redo; checked against five mutations, all caught).
Found F77 (three placeholder dialogs, dead code).
Third step: the panels (the floating cards and their tabs, opened and
closed with an animation; the plot controls, each redrawing the plot;
axis limits, redshift cursor, view; the welcome view; the state of the
window with and without a session, `_update_ui_state`) in a mixin,
`main_panels.py` (`MainPanelsMixin`), moved verbatim with the constants
of the sidebars and cards; `main_window.py` 3984 → 2920 lines. The slots
of the continuum editor's buttons stay for the next step.
Characterization tests in `tests/test_main_panels.py` (checked against
six mutations: five caught; the sixth, the end position of the tab in
the opening animation, is reset by the layout at the end). Found F78
(the "Logarithmic F" box is not in the panel).
Fourth step: the controls of the continuum editor (Start/Save, Cancel,
Reset, the spacing of the knots on a log slider or typed) in a mixin,
`continuum_controls.py` (`ContinuumControlsMixin`), moved verbatim;
`main_window.py` 2920 → 2645 lines. Characterization tests in
`tests/test_continuum_controls.py` (checked against five mutations: four
caught; the fifth, the lines of the cursor not removed, is hidden by the
redraw of the editor). F74 fixed in the same PR (a separate commit): the
view is refreshed with `_update_view_for_session`, as the switch of
session does; `update_gui_session_state`, first proposed, would have
added the state twice (checked with a test, Guido agreed).
Fifth step: the sessions (opening them from the file dialog, a path or
drag and drop; the session list, switching, renaming and its context
menu; the view of a session; a recipe result added to the history or as
a new branch; duplicating, stitching and co-adding; the session
inspector; saving and closing) in a mixin, `main_sessions.py`
(`MainSessionsMixin`), moved verbatim; `main_window.py` 2645 → 1543
lines. Characterization tests in `tests/test_main_sessions.py` (checked
against six mutations, all caught). Found F79 (the parameters of
`equalize_and_stitch` documented in the shared schema).
Sixth and last step: the recipes (requirements and the recipe resumed
after Set Properties, the recipe dialog, the checks before a recipe
runs, the progress dialog and the worker, the result, errors and Stop)
in `main_recipes.py` (`MainRecipesMixin`) and the scripts (log scripter,
script run) in `main_scripts.py` (`MainScriptsMixin`), moved verbatim;
`main_window.py` 1543 → 541 lines (4831 before the phase: `__init__`,
styles, the custom message box, the identification viewer, the debug
plots and the links to the System Inspector). Characterization tests in
`tests/test_main_recipes.py` (checked against six mutations: five
caught; the sixth only matters for a linear recipe returning several
sessions, which none does). Found F80 (dead code), F81 (opening a recipe
dialog switches the zoom on), F82 and F83 (export of ASCII files: an
existing file overwritten without asking, or the question failing),
reported to Guido at once.
F82 and F83 fixed right after (agreed with Guido): the files of an
export of ASCII files are decided by the recipe for both the export and
the GUI (`export_ascii_outputs`), and the question before overwriting
works; the target `spec` writes the whole spectrum, as the user guide
says.
F66 fixed (option chosen by Guido among two): the ticks of the systems
on a velocity axis are referred to z_em, like the other conversions; the
spectrum does not keep the z_ref of a conversion (F84, for phase 6-bis
or 7).

**Result of the phase** (October 2026). Pure refactoring and performance:
the snapshots and the mock benchmark unchanged (no PR touched the
numerical code; the fixes of the phase, F7, F66, F74, F82, F83 and the
calibration of #642, touch the GUI and the export).
- GUI benchmark (`bench_gui.py`, the best of two runs, on one machine,
  against the code at the start of the phase, the merge of #634, measured
  with the same script): first opening of the System Inspector
  5.07 → 0.74 s, next openings 3.22 → 0.50 s, selection of a row (laptop
  size) 0.59 → 0.21 s, redraw of the main plot 0.24 → 0.17 s, the GUI
  part of a refit 2.33 → 0.85 s; the fit itself unchanged. (That machine
  is slower than the one of the earlier measurements: 3.96 s for the
  first opening before the phase.)
- Large files: `core/spectrum.py` 2826 → 774 lines (five mixins of 340-550
  lines); `gui/main_window.py` 4831 → 541 (six mixins and three modules);
  `gui/pyside_plot.py` 2884 → 1266 (`plot_spectrum` ~1000 → ~100 lines,
  two mixins and four modules); `gui/system_inspector.py` 2239 → 1085
  (the table and the velocity plot in their modules). The largest files
  are now the recipes (`absorbers.py` 2116, `edit.py` 1373 lines), whose
  split belongs to phase 6-bis (the recipe API), and `system_list.py`
  (1284) and `spectrum_operations.py` (1142), pure functions.
- Tests: 73 test files; coverage 67% → 81% (`pytest --cov=astrocook`; the
  GUI 81%, the core 90%), because each moved part got its
  characterization tests first, checked with mutations; the interface of
  the split classes is checked by `test_spectrum_structure.py` and
  `test_gui_structure.py`. The suite runs in parallel (~4 min in the CI).
- Documentation: the API reference lists the methods of `SpectrumV2` on
  its page again (they had moved to pages of the mixins with the split,
  and 19 references of the recipes no longer resolved; found at the
  closing of the phase), and the build has no warnings.
- Findings: F66-F84 found in the phase. Fixed: F7, F66, F74, F82, F83
  (and the calibration, #642). For phase 5-bis: F67, F69, F70, F71
  (agreed), and proposed F75 (the SNR view and its error column), F78
  (Logarithmic F not in the panel), F81 (the zoom switched on by the
  recipe dialogs); for Guido's decision: F73, F77, F80 (dead code), F76
  (decimation of long spectra). For phase 6-bis: F79 (a schema rewritten
  at run time), F84 (the z_ref of a velocity axis), with F12. For phase
  7: F68, F72.

F66 (ticks in velocity) is taken up at the end of the phase (agreed with
Guido).

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

Also in this phase (agreed with Guido, October 2026), small fixes to be in
2.3.0: F67 (the continuum from manual knots fails on a spectrum with a
model and no continuum), F69 and F70 (the calibration: the columns it
scales, and the target magnitudes missed by ~0.1 mag with two or more
bands), F71 (an error message that fails) and
[#630](https://github.com/DAS-OATs/astrocook/issues/630) (the errors of
linked parameters are written as 0 in the ASCII list of the systems: a
problem when the list is used for statistics; what error a linked
parameter should carry, e.g. that of its target, scaled through the
expression of the link, to be agreed with Guido).

Found in phase 5 and proposed for this phase (each to be agreed with
Guido): F75 (the default error column of the SNR view looked up as a
column, an error in the log at each plot), F78 (the "Logarithmic F" box
not in the panel), F81 (the recipe dialogs switch the zoom of the main
plot on); and Guido's decision on the dead code (F73, F77, F80) and on
the decimation of long spectra (F76).

Order of the phase (agreed with Guido, October 2026): all the findings
above are in the phase, F75, F78, F81 and F73, F76, F77, F80 included.
Decisions: F69, the columns changed by the flux methods in two lists by
name (first proposed by unit, as `y_convert`, then changed with Guido:
a flux without unit is common, and then the unit does not tell the
columns apart); #630, the error of a linked parameter
propagated from its target through the expression of the link (numerical
derivative); F73, the selection of a region restored (to be discussed
when its turn comes); F76, the decimation kept and redone; F78, the box
added to the panel (to be discussed). One PR at a time:
1. F71 + F67 (core);
2. F69 (columns of the flux methods);
3. F70 (calibration with two or more bands);
4. #630 (errors of linked parameters in the ASCII list);
5. F81 + F75 (zoom after a recipe dialog, SNR view);
6. F78 ("Logarithmic F" in the panel);
6-bis. F85, the System Inspector: the boxes above the table touch (no spacing in their row) and "Group View (Selected Component Only)" does not say what it shows (the group of the selection: linked and blended components); "Only group" with a tooltip, explicit spacing (agreed with Guido);
7. F77 + F80 (dead code);
8. F73 (selection of a region restored);
9. F76 (decimation of long spectra redone);
9-bis. F86, the legend of the main plot removed (only in the panel), with swatches for the other sessions and the auxiliary column (agreed with Guido);
10. the proposals of the review, each agreed first: one theme module
    (agreed with Guido: PR A, the colours of the interface in
    `gui/theme.py` with no visible change, and the icon of Select from
    FontAwesome like those of Matplotlib; PR B, a real dark mode, the plots
    staying white, with F87; then the toolbar of the plots redesigned in a
    PR of its own, keeping all its functions, Guido's request); the
    main plot; the menus (possibly two PRs: order, then recent files and
    window geometry); the System Inspector without a continuum; status
    bar and keyboard navigation of the System Inspector;
10-bis. the visual review (Guido's request, October 2026): a screenshot of
    every window and dialog, light and dark, with observations and a
    proposed place for each (agreed with Guido: when each turn comes, decide
    whether to postpone it). First the two bugs it found, F89 and F90, with
    the two roles of the dark palette missing (links: "Claudia Cupani" in
    *About* dark blue on dark; the placeholder of the script console
    invisible); then a pass on the texts: one rule for capitals (Title Case
    in menus and buttons, sentence case in options), "Show ..." instead of
    "Toggle ...", whole words instead of "Trans.", "Aux. Column",
    "Error Col.", "Red. χ²" where there is room, "Set R..." → "Set
    Resolution...", the session log in a fixed-width font, the median step
    of *Session Info* said plainly. Then, with the proposals above: the
    main plot (the labels of the emission lines over the flux and with
    internal names, e.g. Ly_b, NV_1238: typographic names, above the data;
    the unit of the y axis in the usual order); the System Inspector (the
    error of z in scientific notation, e.g. "1.943529 ± 4.54e-06": in km/s
    or with the digits of z; the column btur, almost always a frozen 0,
    hidden by default; the row of fields above the plot in two groups,
    view and new component, with the keyboard navigation); the *Data
    Inspector* (every column in scientific notation, also λ: a format per
    column); the script console (no titles, the output dark also in the
    light appearance); the empty window (no direct action: a button to
    open, the recent files, drop a file, with the menus); the
    *Identifications* window as a table (species, z, region). For phase
    6-bis: the long labels of the recipe dialogs (e.g. `forest_auto`
    wider than 800 px; short labels with the explanation in the tooltip,
    possibly earlier), their descriptions written for developers, the
    typed widgets, the selectable buttons of the ASCII export;
11. closing of the phase.
None of them changes the results of fitting, identification or the
automatic pipelines: the fast suite and targeted tests (F70: the
magnitudes reached) are enough; the GUI benchmark before and after F76
and the theme module.

**Result of the phase** (November 2026, PRs #664-#688). The GUI changed
in look and ergonomics, each change agreed with Guido on screenshots
before the PR; the screenshots of the tutorials rebuilt at each PR that
changed the GUI.
- Numbers: the snapshots and the mock benchmark unchanged (no snapshot
  regenerated in the phase). Three fixes change numbers outside them, as
  agreed: F69 (the columns scaled, calibrated and smoothed by the flux
  methods), F70 (the calibration with two or more bands reaches the
  target magnitudes, ~0.1 mag before) and #630 (the errors of linked
  parameters propagated through the link instead of 0; the fitted values
  and the errors of the free parameters unchanged).
- The interface: dark mode following the system, or chosen in
  *View > Appearance*, with the colours in one module (`gui/theme.py`);
  the toolbar of the plots redesigned (one family of icons, keys) and
  floating over the plot as the side cards, with its tab, remembered
  shown or hidden also in the System Inspector; a status bar with the
  outcome of the last recipe and the coordinates of the mouse, in both
  windows; the texts in one rule (Title Case in menus and buttons,
  whole words); the main plot without title and legend (swatches in the
  panels instead, F86), the emission lines in spectroscopic notation, the
  ticks of the systems in two softer rows, the redshift of the cursor on
  the plot (F92); the System Inspector with the error of z in km/s, btur
  hidden while unused (a menu of the columns), the boxes in groups, keys
  for refit and delete, a message instead of empty panels without a
  continuum; the *Identifications* window a table of the doublets as the
  recipe finds them; the *Data Inspector* with a format per column; the
  welcome screen, *File > Open Recent*, the window geometry remembered,
  several files loaded in the background with their progress, the File
  menu in groups with the imports; *View > Plot Colours*, the standard
  palette with a thinner model or a colour-blind safe one.
- Restored or redone: the selection of a region (F73, a Select tool), the
  decimation of long spectra (F76), "Logarithmic F" in the panel (F78);
  dead code removed (F77, F80).
- GUI benchmark (`bench_gui.py`, the best of two runs, alternating, on one
  machine; the end of phase 5, #663, against master after #688): first
  opening of the System Inspector 0.56 → 0.50 s, next openings
  0.34 → 0.28 s, redraw of the main plot 0.127 → 0.099 s (the legend
  removed, F86), selection of a row 0.174 → 0.144 s, the GUI part of a
  refit 0.61 → 0.57 s; the fit itself unchanged (10.7 s). No regression
  from the new features.
- Tests: 77 test files (73 at the end of phase 5); coverage 81% → 82%
  (the GUI 81% → 84%, the core 90%). The documentation builds without
  warnings; its review at the closing fixed two names of menu entries
  left from before the phase (*Identify Absorption Lines*, *Export to
  ASCII...*).
- Findings: fixed in the phase F67, F69, F70, F71, F73, F75, F76, F77,
  F78, F80, F81 (from phase 5) and F85-F92 (found in the phase; F90 also
  touched the identification: the identifications of old regions were
  lost when a pipeline ran again); #630. For phase 7: F93 (the
  identification by doublets only, the pairs not stored). For phase
  6-bis: the long labels and the typed widgets of the recipe dialogs, the
  buttons of the ASCII export. Kept as it is: the script console (Guido:
  in the code until it is clear whether it is needed with the MCP
  server, or can be given another use).

Then release 2.3.0 (see "Releases"; released in November 2026).

### Phase 6 – Error handling
Broad `except` clauses, case by case: specific exception, explicit log, or
an item to discuss.

**Plan of the phase** (agreed with Guido, November 2026):
- *Scope.* The broad `except` clauses (inventory below); #566 (loggers
  per module, the messages shared by several recipes in
  `recipes/messages.py`); the log messages rationalized, levels and
  texts, by the rule below (Guido's question at the start of the phase:
  done here and not in 6-bis, because the move to loggers per module
  touches every logging call anyway, and 6-bis then finds the warnings
  that matter for a result already identified); the cosmetic warnings of
  the libraries filtered (the `VerifyWarning` of astropy on the
  `resol_fwhm` keyword at every save, section 4.2). Left to phase 6-bis:
  exceptions instead of `return 0` in the recipes, the structured
  warnings with a code in the result of a recipe, the log of the values
  actually used (#573).
- *Order*, one PR at a time: (1) the inventory, a test that keeps it
  up to date, the rule of the messages (no change of behaviour, #691);
  (2) F94, found with the inventory (data loss: brought forward, agreed
  with Guido, #692); (3) the loggers per module, with characterization
  tests of what reaches the log window of the GUI and the MCP server first
  (F96 fixed with it);
  then one PR per package, each with its `except` clauses and its
  messages: (4) core, fitting, io, in two PRs (agreed with Guido): (4a)
  the rows N and R and the messages, no change of behaviour; (4b) the
  rows D, as decided with Guido (below); (5) recipes, in the same way:
  (5a) the rows N, R and B and the messages, (5b) the rows D; (6) GUI, in
  two PRs (agreed with Guido): (6a) the rows N and C and the messages,
  no change of behaviour, (6b) the rows B and D and a hook that shows an
  unexpected error to the user; (7) the warnings of the libraries. The rows marked D are reviewed with Guido before their
  PR.
  The rule of the messages agreed with Guido as proposed (November 2026).
- *Numbers.* The phase does not change the results: the snapshots and
  the fast suite cover each PR; the slow tests run in the CI job after
  the merge. A fallback that changes a result (rows D) is a finding,
  fixed only with Guido's agreement and with a regression test.

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

| Package | N | B | C | R | D | Total |
|---|---|---|---|---|---|---|
| `core` | 7 | 2 | – | 8 | 8 | 25 |
| `fitting` | 4 | – | – | – | 2 | 6 |
| `io` | 4 | – | – | 12 | 2 (1 done) | 18 |
| `recipes` | 7 | 28 | – | 1 | 10 | 46 |
| `gui` | 53 | 17 | 17 | – | 2 | 89 |
| `mcp_server` | – | 1 | – | – | – | 1 |
| launcher (`_run_astrocook.py`) | 3 | 1 | – | – | – | 4 |
| **Total** | **78** | **49** | **17** | **21** | **24** | **189** |

Done: the rows N and R of core, fitting and io (35, PR 4a, #694), their
rows D (11, PR 4b), row D of `save_archive_v2` (F94, PR 2): core,
fitting and io keep only the two boundaries of `load_session_from_file`;
the rows N, R and B of the recipes (36, PR 5a, #696: the B rows logged
once, with the traceback, by the rule), their rows D (10, PR 5b): of the
recipes, 4 clauses of `absorbers.py` and 2 of `edit.py` removed or
narrowed, 4 kept as boundaries (the cleaning of `lya_auto` and
`forest_auto`, `delete_component`, `detect_anchor`); the rows N and C of
the GUI (70, PR 6a): narrowed to the exceptions of `gui/errors.py`
(`DATA_ERRORS`, the errors of the data shown; `QT_ERRORS`, those of the
Qt and matplotlib objects already gone or not ready) or to the specific
ones (a missing file, a series that cannot be parsed); where a missing
column was caught as the `AttributeError` of `None.value`, the check is
explicit; with them, the fallback of the matplotlib style, which applied
`'fast'` again in the `except` (it would fail again); the rows B and D of
the GUI and the rows of the launcher (23, PR 6b): the boundaries of the
slots that only logged (opening the Log Scripter and the Identifications
window, switching and renaming a session, reading the parameters of a
dialog) removed, the hook shows their errors; the export without its own
clause (`export_ascii` is a boundary: its False is now shown in a dialog);
the save of the knots catches only their `ValueError`; kept as
boundaries, logged with the traceback: the update after a recipe, the
opening and the save of a file, the workers, the mouse and the cursor of
the plot (an error there would otherwise show a dialog at each motion),
the debug plot, the script console, the Dock icon of macOS (pyobjc) and
the files of the command line. 47 clauses left, all boundaries (B):
the entry points of the recipes (28), the workers (3), the tool of the
MCP server, the GUI and the launcher; phase 6-bis turns those of the
recipes into exceptions.

**Warnings of the libraries** (PR 7): the fast suite gave 123 warnings,
now none. The `VerifyWarning` on `resol_fwhm` of the MCP pilot was gone
already (#526 no longer writes it). What was left, each silenced where it
is raised or removed at its source:
- the `AstropyUserWarning` on the keywords of the FITS structure (`NAXIS`,
  `TTYPEn`, `TFORMn`...) at every save of a session from a V1 archive (its
  metadata keep the first lines of the V1 header): left out of the
  metadata written (`_fits_meta`, with the rule of astropy:
  `REMOVE_KEYWORDS`, `is_column_keyword`), the header written is the same;
- the `VerifyWarning` on the long keywords (`series_map_json`...) written
  as HIERARCH cards, the standard way: silenced in `_write_fits_table`;
- the warning of astropy that skipped `bayesian_results` (F28): replaced by
  a WARNING of Astrocook that says what is lost;
- two warnings of Astrocook itself: `invalidateFilter()` of Qt, deprecated
  since 6.10 (the minimum of PySide6), replaced by
  `beginFilterChange`/`endFilterChange`; the `All-NaN slice` of
  `compute_auto_fudge`, a case already handled (silenced there).
  The suite then turns the warnings into errors (`filterwarnings =
  ["error"]` in `pyproject.toml`, Guido's choice among three: all of them,
  only those of the code of Astrocook, only the deprecations), so that a
  new one is seen at once: a test that expects a warning says so
  (`pytest.warns`, or `warnings.catch_warnings` with a local filter); a
  warning brought by a new release of a dependency is fixed in the code or
  ignored by an entry of its own, with its reason. The lines are those of `a5b506a`; the
function identifies the clause when the lines move.

**Rows D of core, fitting and io** (reviewed with Guido, November 2026;
applied in PR 4b, each with a test that reproduces the failure, in
`tests/test_error_handling_core.py`; with them the two fallbacks of
`fit_continuum` and `update_continuum_from_knots` found in PR 4a, F95):
- *A result that looks valid: the error reaches the recipe* (which logs
  it and gives no new session): `smooth_spectrum` (it returned the
  unsmoothed array), `merge_regions_by_velocity` (the regions unmerged),
  `_find_kinematic_doublet_candidates` (no candidates for the whole
  multiplet), `_renormalize_model` (it can fail only by a bug);
  `update_continuum_from_knots` raises a `ValueError` with a clear
  message (two knots at the same wavelength, the likely cause), shown by
  the GUI, instead of returning the spectrum unchanged.
- *A reasonable fallback, narrowed*: `rate_doublet_candidate` keeps the
  score 0 ("not a doublet") for a candidate that cannot be rated, only
  for the expected errors, at DEBUG (a WARNING with traceback for each
  candidate before); the others reach the recipe.
- *Data loss: fail at once*: the fallback constraint model of
  `SystemListV2` removed (with an empty model a save wrote the archive
  without links and frozen parameters, with `None` the save crashed);
  an unreadable V1 system list makes the opening of the `.acs` fail with
  its message (it opened the session without systems, and a save in
  `.acs2` then lost them), as for the V2 archives (Guido's choice
  between failing and opening with a warning shown to the user).
- *Links*: `update_constraint` refuses a link that does not compile or
  names a component that does not exist, with an error shown to the
  user; an old archive with such a link still opens, with an ERROR that
  says the parameter is kept frozen at its value (`SyntaxError` only);
  a link that fails to evaluate during a fit makes the fit fail with
  its message instead of keeping the old value silently. As done: each
  link is checked when it is compiled, by evaluating it with the
  current values (`SyntaxError`, and the errors of an expression,
  `LINK_ERRORS` in `fitting/voigt_model.py`), so a link that names a
  missing component is caught there too, not only one that does not
  compile; `SystemListV2.update_constraint` evaluates a new link in the
  same way and raises a `ValueError` (only when an expression is given:
  freezing a column does not re-check the old links).
- *fit_continuum, update_continuum_from_knots* (F95, found in PR 4a):
  fewer than two unabsorbed pixels, no knots, knots and fluxes of
  different lengths, a single knot or two knots at the same wavelength
  raise a `ValueError` with its message, instead of returning the
  spectrum unchanged; the continuum editor of the GUI shows it in a
  dialog (the knots edited are lost, as before: keeping the editor
  open on a failure is for the PR of the GUI).

**Rows D of the GUI and the hook** (reviewed with Guido, November 2026;
applied in PR 6b): the two `except Exception: pass` around the entries of
the session log in `_process_recipe_result` removed (`add_entry` only
appends to a list: an error there is a bug, which reaches the boundary of
`_process_recipe_result` and its dialog; the entry written before the new
state is applied, so that the log and the state stay together; a failure
dropped the entry silently while the result was applied, and the replay
diverged); the refresh of the Log Scripter in a clause of its own (C,
`RuntimeError`); a log manager that is not a `HistoryLogV2` an ERROR
instead of a silent skip. `V1LogArtifact` (a wrapper of V1 logs that was
never created) removed in PR 6a (Guido). The hook: a `sys.excepthook`
installed by the launcher (PySide6 calls it for the exceptions not
handled in the slots) logs the error with its traceback and shows a
dialog with the message and a *Details* button with the traceback.

**Rows D of the recipes** (reviewed with Guido, November 2026; applied in
PR 5b, each with a test):
- `components_auto`: `region_limit` checked at the start (an unknown
  region makes the recipe fail with its message; it cleaned nothing and
  then masked nothing, so the pipeline worked on the whole spectrum and
  added duplicates); the `except` of the cleaning and of the merge of the
  `abs_ids` removed (an error there is a bug). The same check of
  `region_limit` in `refit_all` (it fitted the whole spectrum, with a
  WARNING), found in PR 5a.
- `lya_auto`, `forest_auto` (the cleaning of the lines narrower than
  `min_b`, at the end of long pipelines): option (b), Guido's choice:
  the result is kept, with a WARNING that says how many narrow lines it
  still contains, the traceback at DEBUG; the clause becomes a boundary
  (B) of the step.
- `metals_auto`: only the numerical errors of a trial fit (`ValueError`,
  `ArithmeticError`) reject the candidate; a WARNING at the end counts
  the failed trial fits, when there are any; the other errors reach the
  recipe.
- `refit_all`: only the numerical errors of a group; one WARNING at the
  end with the number of the groups that failed and their redshifts
  (a flag `not_fitted` on their components: phase 6-bis, Guido).
- `delete_component` returns 0 on a failure (it returned the session
  unchanged, which the GUI took as a success); `detect_anchor` becomes
  a boundary (B): its 0.0 is the failure of the recipe convention.
- import of systems with a refit: option (b), Guido's choice: the
  import stays, with a WARNING that the systems keep their imported
  values when the refit fails; the broad `except` removed.
- `equalize_and_stitch`: ranges or factors that cannot be read, or not
  one per session, make the recipe fail with the message of the
  parameters (they were replaced by `auto`, with a WARNING).
- Found in PR 5a, agreed for 5b: `import_systems` from a session that
  is not open, `trim_common` without a common range and `stitch`
  without sessions, or without valid ones, return 0 with an ERROR (they
  returned the session unchanged, as `delete_component` did).
- As done in PR 5b: the trial fits of `metals_auto` and the groups of
  `refit_all` catch `FIT_ERRORS` (`ValueError`, `ArithmeticError`);
  `metals_auto` has no boundary of its own, so another error reaches its
  caller (the worker of the GUI); an unreadable `min_b` of `lya_auto`
  and `forest_auto` makes them fail with the message of the parameters
  (it skipped the cleaning with a WARNING); the regions are checked by
  `_region_is_valid`, so `lya_auto` with `z_em` not set now fails (the
  region of the forest cannot be computed: it ran on the whole
  spectrum).

<details>
<summary>The 189 clauses</summary>

| File | Line | Function | Catches | Action | What it does |
|---|---|---|---|---|---|
| `_run_astrocook.py` | 103 | `main` | `Exception` | N: done | Dock icon (macOS, `AppKit`) |
| `_run_astrocook.py` | 106 | `main` | `Exception` | N: done | window icon |
| `_run_astrocook.py` | 153 | `main` | `Exception` | N: done | splash screen |
| `_run_astrocook.py` | 207 | `main` | `Exception` | B: done | files given on the command line |
| `core/session.py` | 74 | `load_session_from_file` | `Exception` | N: done | format detection, falls back to `generic_spectrum` |
| `core/session.py` | 82 | `load_session_from_file` | `Exception` | B | V2 loader, returns 0 (removed in 6-bis, PR 3: an exception) |
| `core/session.py` | 172 | `load_session_from_file` | `Exception` | B | archive or spectrum loading, returns 0 (removed in 6-bis, PR 3: an exception) |
| `core/session.py` | 296 | `SessionV2.open_new` | `Exception` | R: done | wraps `load_session_from_file`, which already catches and logs |
| `core/spectrum.py` | 483 | `SpectrumV2._renormalize_model` | `Exception` | D: done | model left against the old continuum, with a WARNING |
| `core/spectrum_combine.py` | 205 | `CombineMixin.rebin` | `Exception` | N: done | aux column dropped by `rebin`, with a WARNING |
| `core/spectrum_continuum.py` | 312 | `ContinuumMixin.update_continuum_from_knots` | `Exception` | D: done | continuum from knots: returns the spectrum unchanged, the edit is lost (ERROR only) |
| `core/spectrum_expressions.py` | 155 | `prepare_expression_context.replacer` | `Exception` | N: done | re-raised as `ValueError`: add `from e` |
| `core/spectrum_expressions.py` | 194 | `ExpressionsMixin.apply_expression` | `Exception` | R: done | log and re-raise: the expression goes in the exception message |
| `core/spectrum_expressions.py` | 283 | `ExpressionsMixin.mask_expression` | `Exception` | R: done | log and re-raise: the expression goes in the exception message |
| `core/spectrum_expressions.py` | 420 | `ExpressionsMixin.split` | `Exception` | R: done | log and re-raise: the expression goes in the exception message |
| `core/spectrum_flux.py` | 200 | `FluxMixin.smooth` | `Exception` | N: done | aux column left unsmoothed, with a WARNING |
| `core/spectrum_flux.py` | 321 | `FluxMixin.calibrate` | `Exception` | R: done | log and re-raise |
| `core/spectrum_lines.py` | 162 | `LinesMixin.identify_lines` | `Exception` | R: done | log and re-raise; the recipe logs the traceback again |
| `core/spectrum_lines.py` | 177 | `LinesMixin.identify_lines` | `Exception` | R: done | log and re-raise; the recipe logs the traceback again |
| `core/spectrum_operations.py` | 150 | `smooth_spectrum` | `Exception` | D: done | `smooth_spectrum` returns the unsmoothed array (ERROR only) |
| `core/spectrum_operations.py` | 795 | `fit_powerlaw_to_regions` | `Exception` | R: done | log and re-raise |
| `core/spectrum_operations.py` | 851 | `merge_regions_by_velocity` | `Exception` | D: done | regions returned unmerged (ERROR only) |
| `core/spectrum_operations.py` | 1022 | `_find_kinematic_doublet_candidates` | `Exception` | D: done | no doublet candidates for the multiplet (ERROR only) |
| `core/spectrum_operations.py` | 1140 | `rate_doublet_candidate` | `Exception` | D: done | candidate rated 0 (WARNING only) |
| `core/system_list.py` | 249 | `SystemListV2._ensure_clean_data` | bare | N: done | `literal_eval` of a key (bare `except`) |
| `core/system_list.py` | 298 | `SystemListV2._initialize_constraint_model` | `Exception` | D: done | empty fallback constraint model: links and frozen parameters lost (ERROR only) |
| `core/system_list.py` | 306 | `SystemListV2._initialize_constraint_model` | `Exception` | D: done | empty fallback constraint model: links and frozen parameters lost (ERROR only) |
| `core/system_list.py` | 431 | `SystemListV2.fitting_context` | `Exception` | N: done | restore of the active components after a fit |
| `core/system_list.py` | 436 | `SystemListV2.fitting_context` | `Exception` | N: done | restore of the active components after a fit |
| `fitting/bayesian_fitter.py` | 108 | `BayesianVoigtFitter.ln_likelihood` | `Exception` | N: done | likelihood −inf for a model that cannot be evaluated: catch the numerical errors only |
| `fitting/voigt_fitter.py` | 414 | `VoigtFitterV2.evaluate_uncertainties` | `Exception` | N: done | uncertainties set to NaN when the covariance fails |
| `fitting/voigt_fitter.py` | 940 | `VoigtFitterV2.fit` | `Exception` | N: done | uncertainties set to NaN when the covariance fails |
| `fitting/voigt_model.py` | 160 | `VoigtModelConstraintV2._sanitize_data_structures` | bare | N: done | parsing of a key (bare `except`) |
| `fitting/voigt_model.py` | 260 | `VoigtModelConstraintV2._compile_links` | `Exception` | D: done | a link that does not compile is skipped (ERROR only) |
| `fitting/voigt_model.py` | 351 | `VoigtModelConstraintV2.map_p_free_to_full` | `Exception` | D: done | a link that fails to evaluate keeps the old value, silently, at every model evaluation |
| `gui/about_dialog.py` | 41 | `AboutDialog.__init__` | `Exception` | N: done | logo or icon |
| `gui/continuum_controls.py` | 97 | `ContinuumControlsMixin._toggle_continuum_editor` | bare | C: done | removal of a matplotlib artist |
| `gui/continuum_controls.py` | 201 | `ContinuumControlsMixin._toggle_continuum_editor` | `Exception` | B: done | saving the edited continuum: ERROR in the log only, the user should see it |
| `gui/continuum_editor.py` | 210 | `ContinuumEditorMixin._remove_editor_artists` | bare | C: done | removal of a matplotlib artist |
| `gui/continuum_editor.py` | 215 | `ContinuumEditorMixin._remove_editor_artists` | bare | C: done | removal of a matplotlib artist |
| `gui/continuum_editor.py` | 344 | `ContinuumEditorMixin.refresh_edit_layer` | `Exception` | N: done | refresh of the editor layer (DEBUG) |
| `gui/data_inspector.py` | 223 | `DataInspectorDialog._populate_range` | bare | N: done | value of a column (bare `except`) |
| `gui/identification_viewer_dialog.py` | 185 | `IdentificationViewerDialog._populate_results` | `Exception` | N: done | `abs_ids` column |
| `gui/identification_viewer_dialog.py` | 211 | `IdentificationViewerDialog._populate_results` | `Exception` | N: done | JSON of the identifications, shown in the dialog |
| `gui/identification_viewer_dialog.py` | 267 | `IdentificationViewerDialog._get_ref_transition` | `Exception` | N: done | `trans_parse` of a series |
| `gui/inspector_table.py` | 217 | `SystemTableModel.data` | bare | N: done | cell value as a number (bare `except`, hot path) |
| `gui/log_scripter_dialog.py` | 193 | `LogScripterDialog._populate_v1_script` | `Exception` | N: done | V1 log JSON |
| `gui/log_scripter_dialog.py` | 290 | `LogScripterDialog._on_save_log` | `Exception` | N: done | saving the log: `OSError`, shown in a dialog |
| `gui/log_scripter_dialog.py` | 313 | `LogScripterDialog._on_insert_recipe` | `Exception` | N: done | recipe menu and template of the log scripter |
| `gui/log_scripter_dialog.py` | 342 | `LogScripterDialog._insert_template` | `Exception` | N: done | recipe menu and template of the log scripter |
| `gui/main_panels.py` | 189 | `MainPanelsMixin._setup_empty_view` | `Exception` | N: done | logo or icon |
| `gui/main_panels.py` | 754 | `MainPanelsMixin._update_limit_boxes_from_plot` | `Exception` | N: done | limit boxes from the plot |
| `gui/main_panels.py` | 784 | `MainPanelsMixin._force_restack_floating_widgets` | `Exception` | N: done | re-stacking of the floating widgets |
| `gui/main_recipes.py` | 160 | `MainRecipesMixin._launch_recipe_dialog` | `Exception` | C: done | reset of the plot toolbar |
| `gui/main_recipes.py` | 273 | `MainRecipesMixin._on_recipe_requested` | `Exception` | B: done | export from the dialog: ERROR in the log only |
| `gui/main_recipes.py` | 497 | `MainRecipesMixin._process_recipe_result` | `Exception` | D: done | the entry of a recipe in the session log is lost silently |
| `gui/main_recipes.py` | 542 | `MainRecipesMixin._process_recipe_result` | `Exception` | D: done | the entry of a recipe in the session log is lost silently |
| `gui/main_recipes.py` | 627 | `MainRecipesMixin._process_recipe_result` | `Exception` | B: done | processing of a recipe result, shown in a dialog |
| `gui/main_recipes.py` | 711 | `MainRecipesMixin._cleanup_progress_ui` | `Exception` | C: done | closing of the progress dialog and of its log bridge |
| `gui/main_recipes.py` | 810 | `MainRecipesMixin.launch_split_from_current_view` | `Exception` | N: done | plot limits for Split |
| `gui/main_recipes.py` | 844 | `_safely_close_progress_dialog._destroy` | bare | C: done | closing of the progress dialog and of its log bridge |
| `gui/main_recipes.py` | 853 | `_safely_close_progress_dialog._destroy` | `Exception` | C: done | closing of the progress dialog and of its log bridge |
| `gui/main_scripts.py` | 57 | `MainScriptsMixin._launch_log_scripter` | `Exception` | B: done | opening of the log scripter |
| `gui/main_scripts.py` | 71 | `MainScriptsMixin._on_log_scripter_closed` | `Exception` | C: done | re-raising of the side buttons |
| `gui/main_scripts.py` | 98 | `MainScriptsMixin.run_script` | `Exception` | C: done | reset of the plot toolbar |
| `gui/main_sessions.py` | 99 | `MainSessionsMixin.update_gui_session_state` | `Exception` | N: done | copy of the session log for a branch: aborts with a dialog |
| `gui/main_sessions.py` | 201 | `MainSessionsMixin.update_gui_session_state` | `Exception` | N: done | refresh of the log scripter |
| `gui/main_sessions.py` | 314 | `MainSessionsMixin._update_view_for_session` | `Exception` | N: done | aux column combobox |
| `gui/main_sessions.py` | 364 | `MainSessionsMixin._update_view_for_session` | `Exception` | N: done | selection in the session list |
| `gui/main_sessions.py` | 400 | `MainSessionsMixin._on_session_switched` | `Exception` | B: done | Qt slot (session switched, renamed) |
| `gui/main_sessions.py` | 823 | `MainSessionsMixin._on_session_info` | bare | N: done | S/N and chi2 of the info panel |
| `gui/main_sessions.py` | 898 | `MainSessionsMixin._on_session_info` | `Exception` | N: done | S/N and chi2 of the info panel |
| `gui/main_sessions.py` | 1021 | `MainSessionsMixin._on_session_name_changed` | `Exception` | B: done | Qt slot (session switched, renamed) |
| `gui/main_sessions.py` | 1069 | `MainSessionsMixin._save_session` | `Exception` | B: done | saving or opening a session, shown in a dialog |
| `gui/main_sessions.py` | 1186 | `MainSessionsMixin.open_session_from_path` | `Exception` | B: done | saving or opening a session, shown in a dialog |
| `gui/main_window.py` | 124 | `MainWindowV2.__init__` | `Exception` | N: done | logo or icon |
| `gui/main_window.py` | 249 | `MainWindowV2._apply_styles` | `Exception` | N: done | palette query, fallback colours |
| `gui/main_window.py` | 369 | `MainWindowV2._show_custom_message` | `Exception` | N: done | icon of a message box |
| `gui/main_window.py` | 463 | `MainWindowV2._launch_identification_viewer` | `Exception` | C: done | closing of a previous dialog |
| `gui/main_window.py` | 480 | `MainWindowV2._launch_identification_viewer` | `Exception` | B: done | opening of the identifications window |
| `gui/main_window.py` | 506 | `MainWindowV2._flux_limits_in` | `Exception` | N: done | x axis in nm for the flux limits |
| `gui/main_window.py` | 600 | `MainWindowV2._on_debug_plot_requested` | `Exception` | B: done | debug plot |
| `gui/plot_canvas.py` | 48 | `MatplotlibCanvas.__init__` | `Exception` | N: done | matplotlib style (the message still names `scienceplots`) |
| `gui/plot_canvas.py` | 149 | `MatplotlibCanvas._capture_background` | `Exception` | N: done | background copy for blitting |
| `gui/plot_canvas.py` | 155 | `MatplotlibCanvas._capture_background` | `Exception` | C: done | disconnection of a matplotlib callback |
| `gui/plot_canvas.py` | 171 | `MatplotlibCanvas.draw_animated` | `Exception` | N: done | blitting |
| `gui/plot_canvas.py` | 181 | `MatplotlibCanvas.draw_animated` | `Exception` | N: done | blitting |
| `gui/plot_canvas.py` | 268 | `MatplotlibCanvas.on_motion` | `Exception` | B: done | mouse motion slot |
| `gui/plot_canvas.py` | 303 | `MatplotlibCanvas.cleanup_selection` | `Exception` | C: done | removal of a matplotlib artist |
| `gui/plot_mouse.py` | 108 | `PlotMouseMixin.get_cursor_line_positions_at_z` | `Exception` | N: done | `trans_parse` of a series |
| `gui/plot_mouse.py` | 144 | `PlotMouseMixin.calculate_z_from_x` | `Exception` | N: done | `trans_parse` of a series |
| `gui/plot_mouse.py` | 162 | `PlotMouseMixin.calculate_z_from_x` | `Exception` | N: done | redshift from the cursor position |
| `gui/plot_mouse.py` | 272 | `PlotMouseMixin.show_data_tooltip` | `Exception` | N: done | data tooltip |
| `gui/plot_mouse.py` | 298 | `PlotMouseMixin.show_data_tooltip` | `Exception` | N: done | data tooltip |
| `gui/plot_mouse.py` | 662 | `PlotMouseMixin._open_data_inspector` | `Exception` | C: done | closing of a previous dialog |
| `gui/plot_style.py` | 189 | `z_convert_inverse` | `Exception` | N: done | inverse redshift conversion of the axis |
| `gui/plot_style.py` | 200 | `get_color_cycle` | `Exception` | N: done | colour map, fallback |
| `gui/plot_style.py` | 206 | `get_color_cycle` | `Exception` | N: done | colour map, fallback |
| `gui/plot_style.py` | 237 | `decimate_y_min_max` | `Exception` | N: done | decimation for plotting |
| `gui/plot_style.py` | 263 | `decimate_x_for_min_max` | `Exception` | N: done | decimation for plotting |
| `gui/pyside_plot.py` | 416 | `SpectrumPlotWidget._connect_background_capture` | `Exception` | C: done | disconnection of a matplotlib callback |
| `gui/pyside_plot.py` | 431 | `SpectrumPlotWidget._was_zoomed` | `Exception` | N: done | zoom state |
| `gui/pyside_plot.py` | 440 | `SpectrumPlotWidget._reconnect_limit_callbacks` | `Exception` | C: done | disconnection of a matplotlib callback |
| `gui/pyside_plot.py` | 443 | `SpectrumPlotWidget._reconnect_limit_callbacks` | `Exception` | C: done | disconnection of a matplotlib callback |
| `gui/pyside_plot.py` | 449 | `SpectrumPlotWidget._reconnect_limit_callbacks` | `Exception` | N: done | limit callbacks |
| `gui/pyside_plot.py` | 531 | `SpectrumPlotWidget._plot_arrays` | `Exception` | N: done | error column of the S/N view, fallback to `dy` |
| `gui/pyside_plot.py` | 570 | `SpectrumPlotWidget._slice_to_view` | `Exception` | N: done | zoom slice |
| `gui/pyside_plot.py` | 635 | `SpectrumPlotWidget._plot_other_sessions` | bare | N: done | S/N of the other sessions (bare `except`) |
| `gui/pyside_plot.py` | 807 | `SpectrumPlotWidget._plot_aux_column` | `Exception` | N: done | aux column |
| `gui/pyside_plot.py` | 853 | `SpectrumPlotWidget._plot_systems` | `Exception` | N: done | series of a system that `trans_parse` cannot read: not drawn, silently |
| `gui/pyside_plot.py` | 878 | `SpectrumPlotWidget._plot_systems` | `Exception` | N: done | ticks of the systems in velocity |
| `gui/pyside_plot.py` | 951 | `SpectrumPlotWidget._plot_cursor` | `Exception` | B: done | redshift cursor |
| `gui/pyside_plot.py` | 1103 | `SpectrumPlotWidget._plot_emission_lines` | `Exception` | N: done | emission line labels |
| `gui/pyside_plot.py` | 1155 | `SpectrumPlotWidget._draw_isomag_grid` | bare | C: done | removal of a matplotlib artist |
| `gui/pyside_plot.py` | 1221 | `SpectrumPlotWidget._draw_isomag_grid` | `Exception` | N: done | iso-magnitude grid |
| `gui/qt_workers.py` | 114 | `RecipeWorker.run` | `Exception` | B: done | worker thread (recipe, script, session loading) |
| `gui/qt_workers.py` | 249 | `ScriptWorker.run` | `Exception` | B: done | worker thread (recipe, script, session loading) |
| `gui/qt_workers.py` | 293 | `SessionLoadWorker.run` | `Exception` | B: done | worker thread (recipe, script, session loading) |
| `gui/recipe_dialog.py` | 45 | `RecipeDialog.__init__` | `Exception` | N: done | schema of the recipe |
| `gui/recipe_dialog.py` | 111 | `RecipeDialog._get_available_columns` | `Exception` | N: done | columns offered in the dialog |
| `gui/recipe_dialog.py` | 551 | `RecipeDialog._check_oversample` | `Exception` | N: done | oversampling warning |
| `gui/recipe_dialog.py` | 597 | `RecipeDialog.accept` | `Exception` | B: done | reading the parameters of the dialog, shown in a dialog |
| `gui/recipe_dialog.py` | 618 | `RecipeDialog._show_help` | `Exception` | N: done | opening of the help page |
| `gui/script_dialog.py` | 63 | `ScriptDialog.run_script` | `Exception` | B: done | script console (`exec`) |
| `io/adapter.py` | 106 | `load_and_migrate_structure` | `Exception` | D: done | unreadable V1 system list: the session opens with no systems (ERROR only) |
| `io/adapter.py` | 376 | `save_archive_v2` | `Exception` | D: done | original data not written, the save goes on and fails later (F94): removed, a failure of the original data stops the save (phase 6, F94 PR) |
| `io/adapter.py` | 390 | `save_archive_v2` | `Exception` | R: done | log and re-raise; the GUI shows the error |
| `io/loaders.py` | 121 | `detect_file_format` | `Exception` | N: done | format detection from the first lines |
| `io/loaders.py` | 221 | `load_astrocook_fits` | `Exception` | R: done | log and re-raise; `load_session_from_file` logs again |
| `io/loaders.py` | 286 | `load_pypeit_coadd` | `Exception` | R: done | log and re-raise; `load_session_from_file` logs again |
| `io/loaders.py` | 310 | `load_pypeit_spec1d_multi` | bare | N: done | order number from the HDU name, 0 if missing (bare `except`) |
| `io/loaders.py` | 385 | `load_pypeit_spec1d_multi` | `Exception` | R: done | log and re-raise; `load_session_from_file` logs again |
| `io/loaders.py` | 440 | `load_pypeit_stack` | `Exception` | R: done | log and re-raise; `load_session_from_file` logs again |
| `io/loaders.py` | 504 | `load_pypeit_tellmodel` | `Exception` | R: done | log and re-raise; `load_session_from_file` logs again |
| `io/loaders.py` | 544 | `load_simple_csv` | `Exception` | R: done | log and re-raise; `load_session_from_file` logs again |
| `io/loaders.py` | 599 | `load_ascii_telluric_model` | `Exception` | R: done | log and re-raise; `load_session_from_file` logs again |
| `io/loaders.py` | 642 | `load_ascii_with_resvel_header` | `Exception` | N: done | non-numeric header line skipped |
| `io/loaders.py` | 699 | `load_ascii_with_resvel_header` | `Exception` | R: done | log and re-raise; `load_session_from_file` logs again |
| `io/loaders.py` | 813 | `load_espresso_s2d_fits` | `Exception` | R: done | log and re-raise; `load_session_from_file` logs again |
| `io/loaders.py` | 906 | `load_xshooter_fits` | `Exception` | R: done | log and re-raise; `load_session_from_file` logs again |
| `io/loaders.py` | 1006 | `load_generic_spectrum` | `Exception` | R: done | log and re-raise; `load_session_from_file` logs again |
| `io/v1_archive.py` | 67 | `V1ArchiveManager.unpack` | `Exception` | N: done | `.acs` archive that cannot be unpacked |
| `mcp_server/errors.py` | 85 | `run_recipe` | `Exception` | B | a tool call: the exception becomes a `ToolError` |
| `recipes/absorbers.py` | 507 | `RecipeAbsorbersV2.components_auto` | `Exception` | B: done | `components_auto`, mask step |
| `recipes/absorbers.py` | 530 | `RecipeAbsorbersV2.components_auto` | `Exception` | D: done | `components_auto`: old systems of the region not removed, the pipeline goes on (WARNING only) |
| `recipes/absorbers.py` | 586 | `RecipeAbsorbersV2.components_auto` | `Exception` | D: done | `components_auto`: `abs_ids` not merged (WARNING only) |
| `recipes/absorbers.py` | 734 | `RecipeAbsorbersV2.lya_auto` | `Exception` | D: done | `lya_auto`: cleanup of the narrow components skipped (WARNING only) |
| `recipes/absorbers.py` | 817 | `RecipeAbsorbersV2.forest_auto` | `Exception` | D: done | `forest_auto`: cleanup of the narrow components skipped (WARNING only) |
| `recipes/absorbers.py` | 916 | `RecipeAbsorbersV2.metals_auto` | `Exception` | D: done | `metals_auto`: a failed trial fit logged at DEBUG only |
| `recipes/absorbers.py` | 990 | `RecipeAbsorbersV2.identify_lines` | `Exception` | B: done | recipe entry point, returns 0 (removed in 6-bis, PR 3: an exception); parameter parsing to `msg_param_fail` |
| `recipes/absorbers.py` | 1029 | `RecipeAbsorbersV2.populate_from_identification` | `Exception` | N: done | JSON of the identification maps, returns 0 |
| `recipes/absorbers.py` | 1134 | `RecipeAbsorbersV2.add_component` | `Exception` | B: done | recipe entry point, returns 0 (removed in 6-bis, PR 3: an exception); parameter parsing to `msg_param_fail` |
| `recipes/absorbers.py` | 1179 | `RecipeAbsorbersV2.update_component` | `Exception` | B: done | recipe entry point, returns 0 (removed in 6-bis, PR 3: an exception); parameter parsing to `msg_param_fail` |
| `recipes/absorbers.py` | 1230 | `RecipeAbsorbersV2.update_components` | `Exception` | B: done | recipe entry point, returns 0 (removed in 6-bis, PR 3: an exception); parameter parsing to `msg_param_fail` |
| `recipes/absorbers.py` | 1283 | `RecipeAbsorbersV2.delete_component` | `Exception` | D: done | `delete_component` returns the unchanged session: a failure looks like a success |
| `recipes/absorbers.py` | 1327 | `RecipeAbsorbersV2.detect_anchor` | `Exception` | D: done | `detect_anchor` returns 0.0 on failure |
| `recipes/absorbers.py` | 1394 | `RecipeAbsorbersV2.fit_component` | `Exception` | B: done | recipe entry point, returns 0 (removed in 6-bis, PR 3: an exception); parameter parsing to `msg_param_fail` |
| `recipes/absorbers.py` | 1449 | `RecipeAbsorbersV2.update_constraint` | `Exception` | N: done | auto-update of a linked value |
| `recipes/absorbers.py` | 1454 | `RecipeAbsorbersV2.update_constraint` | `Exception` | B: done | recipe entry point, returns 0 (removed in 6-bis, PR 3: an exception); parameter parsing to `msg_param_fail` |
| `recipes/absorbers.py` | 1658 | `RecipeAbsorbersV2.optimize_system` | `Exception` | B: done | recipe entry point, returns 0 (removed in 6-bis, PR 3: an exception); parameter parsing to `msg_param_fail` |
| `recipes/absorbers.py` | 1753 | `RecipeAbsorbersV2.fit_bayesian` | `Exception` | B: done | recipe entry point, returns 0 (removed in 6-bis, PR 3: an exception); parameter parsing to `msg_param_fail` |
| `recipes/absorbers.py` | 1841 | `RecipeAbsorbersV2.refit_all` | `Exception` | B: done | `refit_all`, consolidated fit, returns 0 |
| `recipes/absorbers.py` | 2006 | `RecipeAbsorbersV2.refit_all` | `Exception` | D: done | `refit_all`: a group whose fit fails is left unfitted (WARNING only) |
| `recipes/absorbers.py` | 2043 | `RecipeAbsorbersV2.refit_all` | `Exception` | B: done | recipe entry point, returns 0 (removed in 6-bis, PR 3: an exception); parameter parsing to `msg_param_fail` |
| `recipes/absorbers.py` | 2070 | `RecipeAbsorbersV2._final_uncertainties` | `Exception` | N: done | final uncertainties of a group |
| `recipes/continuum.py` | 181 | `RecipeContinuumV2.estimate_auto` | `Exception` | B: done | recipe entry point, returns 0 (removed in 6-bis, PR 3: an exception) |
| `recipes/continuum.py` | 236 | `RecipeContinuumV2.find_absorbed` | `Exception` | B: done | recipe entry point, returns 0 (removed in 6-bis, PR 3: an exception) |
| `recipes/continuum.py` | 297 | `RecipeContinuumV2.fit_continuum` | `Exception` | B: done | recipe entry point, returns 0 (removed in 6-bis, PR 3: an exception) |
| `recipes/continuum.py` | 338 | `RecipeContinuumV2.fit_powerlaw` | `Exception` | B: done | recipe entry point, returns 0 (removed in 6-bis, PR 3: an exception) |
| `recipes/edit.py` | 393 | `RecipeEditV2.set_properties` | `Exception` | R: done | `except Exception as e: raise e` |
| `recipes/edit.py` | 515 | `RecipeEditV2.apply_expression` | `Exception` | B: done | recipe entry point, returns 0 (removed in 6-bis, PR 3: an exception) |
| `recipes/edit.py` | 561 | `RecipeEditV2.mask_expression` | `Exception` | B: done | recipe entry point, returns 0 (removed in 6-bis, PR 3: an exception) |
| `recipes/edit.py` | 610 | `RecipeEditV2.smooth_column` | `Exception` | B: done | recipe entry point, returns 0 (removed in 6-bis, PR 3: an exception) |
| `recipes/edit.py` | 689 | `RecipeEditV2._merge_imported_systems` | `Exception` | N: done | imported systems filtered by range |
| `recipes/edit.py` | 725 | `RecipeEditV2._merge_imported_systems` | `Exception` | D: done | import of systems: the automatic refit fails, the unfitted systems are returned (ERROR only) |
| `recipes/edit.py` | 768 | `RecipeEditV2.import_systems_ascii` | `Exception` | N: done | ASCII table of systems that cannot be read |
| `recipes/edit.py` | 952 | `RecipeEditV2.extract_preset` | `Exception` | N: done | `extract_preset`: a region that cannot be extracted is skipped |
| `recipes/edit.py` | 1035 | `RecipeEditV2.trim_common` | `Exception` | B: done | recipe entry point, returns 0 (removed in 6-bis, PR 3: an exception) |
| `recipes/edit.py` | 1107 | `RecipeEditV2.stitch` | `Exception` | N: done | `stitch`: cleanup of the NaNs |
| `recipes/edit.py` | 1113 | `RecipeEditV2.stitch` | `Exception` | B: done | recipe entry point, returns 0 (removed in 6-bis, PR 3: an exception) |
| `recipes/edit.py` | 1191 | `RecipeEditV2.equalize_and_stitch` | `Exception` | D: done | `equalize_and_stitch`: invalid ranges replaced by `auto` (WARNING only) |
| `recipes/edit.py` | 1372 | `RecipeEditV2.coadd` | `Exception` | B: done | recipe entry point, returns 0 (removed in 6-bis, PR 3: an exception) |
| `recipes/file.py` | 105 | `RecipeFileV2.export_ascii` | `Exception` | B: done | `export_ascii`, returns False |
| `recipes/file.py` | 207 | `RecipeFileV2.import_ascii` | `Exception` | B: done | `import_ascii`, returns 0 |
| `recipes/flux.py` | 143 | `RecipeFluxV2.calculate_running_std` | `Exception` | B: done | recipe entry point, returns 0 (removed in 6-bis, PR 3: an exception) |
| `recipes/flux.py` | 177 | `RecipeFluxV2.smooth` | `Exception` | B: done | recipe entry point, returns 0 (removed in 6-bis, PR 3: an exception) |
| `recipes/flux.py` | 290 | `RecipeFluxV2.equalize` | `Exception` | B: done | recipe entry point, returns 0 (removed in 6-bis, PR 3: an exception) |
| `recipes/flux.py` | 341 | `RecipeFluxV2.resample` | `Exception` | B: done | recipe entry point, returns 0 (removed in 6-bis, PR 3: an exception) |
| `recipes/flux.py` | 369 | `RecipeFluxV2.calibrate_from_magnitudes` | `Exception` | B: done | recipe entry point, returns 0 (removed in 6-bis, PR 3: an exception) |

</details>

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

**Result of the phase** (November 2026, PRs #691-#701). Error handling,
without a change of the numerical results:
- Numbers: the snapshots and the mock benchmark unchanged (no snapshot
  regenerated in the phase); the slow tests run before the merge of each
  PR that touched the pipelines (#695, #697, #701). Only failures behave
  differently: they stop the operation with their reason, where a result
  that looked valid went on before.
- The broad `except` clauses: 189 → 47, all boundaries declared as such
  (the entry points of the recipes, 28, which phase 6-bis turns into
  exceptions; the workers of the GUI, 3; the tool of the MCP server; the
  update of the GUI after a recipe, the opening and the save of a file,
  the mouse and the cursor of the plot, the script console; the Dock icon
  and the files of the command line of the launcher).
  `tests/test_broad_excepts.py` holds the count per file, so that no new
  one is added. The others narrowed to the exceptions expected (in the
  GUI, `DATA_ERRORS` and `QT_ERRORS` of `gui/errors.py`), or removed: the
  log-and-raise clauses (an error is logged once, where it is caught, with
  its context), the fallbacks that went on with a result that looked
  valid (F95, the 24 rows D, each decided with Guido and reproduced by a
  test), the clauses of the slots that only logged.
- Logging: one logger per module (`getLogger(__name__)`, #566), checked by
  a test; the messages of the whole package by one rule (the levels by
  what the user needs, sentences without prefixes, one log per error,
  progress lines at INFO), the shared ones in `recipes/messages.py`
  (#566); a script can set the level of Astrocook, or of a part of it,
  alone (scripting tutorial).
- Unexpected errors: a hook (`install_excepthook`, from the launcher)
  logs them with the traceback and shows them in a dialog with the
  details, instead of the terminal only (*Getting Started*, "When
  Something Goes Wrong").
- Warnings: the fast suite 123 → 0, each silenced where it is raised or
  removed at its source; warnings are errors in the tests
  (`filterwarnings = ["error"]`, Guido's choice), with one targeted entry
  (the verbatim V1 copy).
- Found and fixed: F94 (a failed save replaced the archive with a partial
  one: data loss), F95 (fallbacks after a failure), F96 (a script could not
  configure logging after a recipe), F98 (the functions of the links),
  F99 (fits without a resolution said at DEBUG only; the failure of the
  automatic continuum); with them a bug of the fallback of the matplotlib
  style. Found and left: F97 (a short CSV read as a generic spectrum),
  fixed right after the phase (Guido), with F100. `V1LogArtifact` (never created) removed.
- For phase 6-bis: the scheme of the prerequisites of the recipes (one
  principle, a form for each way of use: MCP, GUI, scripts; "Prerequisites
  of the recipes", to be refined with Guido), the lists of the GUI moved
  to the registry, the `return 0` of the recipes turned into exceptions,
  the flag `not_fitted` of the groups that `refit_all` could not fit.
- Tests 1211 → 1293; coverage 84% (core 91%, GUI 84%, recipes 86%;
  `io` 49%: the loaders of formats without sample files, which wait for
  data from Guido).

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
The full MCP server (section 4) is generated from this API.

**Plan of the phase** (agreed with Guido at its start, November 2026):
- *Starting point.* The registry is five dictionaries
  (`*_RECIPES_SCHEMAS`, 48 recipes: file 2, flux 6, continuum 5, edit 15,
  absorbers 20); a parameter has a Python type, a default and a short
  `doc`, and the signatures of the methods repeat the defaults as strings
  (`smooth_len_lya: str = '5000.0'`), parsed inside each recipe. What
  describes a recipe is also kept elsewhere: `BRANCHING_RECIPES`
  (`core/utils.py`), `_RECIPES_REQUIRING_Z_EM` and
  `_RECIPES_REQUIRING_RESOL` (`gui/main_recipes.py`, written by hand), the
  tools of the MCP server (written by hand). A recipe that succeeds
  returns a new `SessionV2` (sessions are immutable): for the 7 branching
  recipes the GUI opens it as a new session, for the others it becomes
  the next state of the same session (undo and redo); on failure it
  returns 0 (84 `return 0` in the recipes, 28 of them the entry points
  declared as boundaries in phase 6), which every caller checks (the
  workers of the GUI, `scripting.run_recipe`, the MCP server). The file
  recipes return other values (the saves a boolean); a recipe stopped by
  the user returns the session unchanged. The log (`LogEntryV2.params`)
  stores the parameters as strings, and is what a replay runs.
- *Constraints.* The snapshots and the mock benchmark unchanged in the
  whole phase: no PR changes the numbers, by construction; the slow tests
  before the merge only for the PRs that touch the pipelines (the
  prerequisites). The logs of the existing `.acs2` archives still replay:
  the conversion from string to type stays one step at the boundary,
  driven by the registry. The scripts of the users and the tutorials
  (`s = s.continuum.estimate_auto(...)`) do not break.
- *Decisions.*
  1. What a recipe returns to a script (Guido): the methods keep returning the new session on success
     (for a branching recipe, the new session; for the others, the next
     state), failure becomes an exception instead of 0; one entry point
     from the registry, `run(session, "continuum.estimate_auto", ...)`,
     returns a `RecipeResult` (the session, and whether it is a new one,
     with the rest of the structured result) and is used by the GUI, the
     MCP server and the scripts that want the full result. The
     alternative, the methods returning a `RecipeResult`, is cleaner but
     would break the scripts.
  2. Prerequisites in the scripts: the default policy `provide` (the
     providers run with their defaults, with a WARNING), so that the
     numbers do not change (Guido).
  3. The user configuration file of the recipe defaults (#573, item 2):
     after the full MCP server, not in this phase (Guido).
  4. F84 (the z_ref of a conversion to velocity): not part of the recipe
     API; phase 7, or a PR of its own (Guido).
  5. The scheme of the prerequisites refined with Guido while PRs 1-3 are
     done, which do not depend on it (Guido).
  Items deferred from this phase that do not belong to phase 7 (such as
  decision 3) can go to a phase 6-ter, to be defined (Guido). (Named
  after 2.5.0: the full MCP server is phase 6-ter, the review of the
  parameters phase 6-quater.)
- *PRs*, one at a time, each with characterization tests first:
  0. this plan (docs only);
  1. *typed registry*: `RecipeSpec` and `ParamSpec` (dataclasses) built
     from the present schemas, with the flags of a recipe (branching,
     interactive, offered to agents) and of a parameter (hidden in the
     GUI); `BRANCHING_RECIPES` and the two lists of the GUI move there; a
     test that the registry and the methods agree (every method has its
     spec, the defaults are the same); the schemas no longer rewritten at
     run time (F79, the half of F12 about the schemas). Done: the
     registry (`recipes/registry.py`: `RECIPES`, `get_recipe`,
     `find_recipe`, frozen `RecipeSpec` and `ParamSpec`) built at import
     from the schemas, which only it reads (a test checks it); each entry
     declares `branching` and `requires` (`z_em`, `resolution`: the lists
     of the GUI as they were, `refit_all` still missing, PR 5); unknown
     keys rejected. The flags "interactive" and "offered to agents" are
     left to PR 6, where the MCP server is their first reader;
  2. *typed parameters*: signatures with real types and one conversion
     from string to type at the boundary (log, GUI, MCP); the values
     checked by the registry (bounds, allowed values; F46); the log
     records the values actually used (#573, item 3); replay of old
     `.acs2` archives tested. In three PRs, as PR 3 (the bodies of the
     recipes change with their types): 2a the conversion and the recipes
     of continuum, flux and file; 2b edit, with F46 (the bounds of
     `set_properties`); 2c absorbers, then the scripting tutorial (the
     parameters typed or as strings). Done in 2a:
     `recipes/params.py`, the decorator `typed_recipe` (the parameters
     passed are converted to the types of the registry before the recipe
     runs; a value that cannot be, or is out of its bounds, stops the
     recipe with `msg_param_fail` and the reason, and 0 until PR 3; the
     parameters that the registry does not declare pass unchanged) and
     `resolve_params` (the values used, defaults included, written by the
     GUI, the script worker and the MCP server into the log); new keys of a
     parameter: `optional` (None is a value, `'None'` or an empty field in
     a dialog), `keep` (special strings, e.g. `'_current_'`), `min`, `max`,
     `allow_nan`; `rebin` with `xstart`, `xend` and `kappa` optional and
     `filling` a number (NaN allowed). Booleans are read without regard
     to case (`'true'` was False for the recipes that compared with
     `'True'`, True for the others), and 1 and 0 are booleans (Guido); an
     empty field is None only where the default is None (an empty `kappa`
     of *Rebin* is an error, not "no clipping": Guido); a boolean given as
     a string is no longer truthy (`update_from_knots` renormalized the
     model with `renorm_model='False'`; a float `fudge` failed with an
     AttributeError). No bounds added yet: proposed to Guido with 2b.
     Coverage (Guido's question): the suite saw none of the changes on
     wrong inputs, so a battery compares each typed recipe with 2.4.0 on
     358 edge inputs (numbers as strings and typed, empty, None, NaN,
     infinities, booleans written in other ways, unknown sessions and
     units): identical on 284, the 74 changes all on invalid or edge
     inputs (refused instead of read as False or 1, or of a crash; NaN
     refused; `'true'`, `'1'` and 1 now True; an empty or None start and
     end of *Rebin* now the data range). `tests/test_recipe_inputs.py`
     holds the table, each change with what 2.4.0 did; 2b and 2c add
     their recipes, with the comparison made before the code changes.
     Done in 2b: the 15 recipes of *Edit* typed (`set_properties` with
     `z_em` a number or `'_current_'` and `resol` a string, R or 'Npx',
     checked by the recipe; *Co-add* as *Rebin*: optional start, end and
     kappa); the bounds decided by Guido before the code changed: F46 (a
     redshift above -1 and finite; a resolution positive and finite),
     and widths and thresholds positive and finite everywhere, the
     recipes of 2a too (smoothing lengths, FWHM and sigma, kappa, dx, the
     window of *Running std* and of *Trim*), no infinity as a value (new
     key `gt`, a bound excluded). The battery of Edit,
     `tests/test_recipe_inputs_edit.py`, run on master before the change:
     324 cases, 82 changed (values out of the bounds refused; None, that
     crashed the recipe, refused with its reason; numbers and booleans
     that were not read now are, e.g. `stitch_wavelengths=398`;
     `append='true'` now adds the systems, it replaced them); the battery
     of 2a gains 38 changes from the bounds (112 of 358 from 2.4.0).
     Done in 2c: the 20 recipes of *Absorbers* typed, so all the 48;
     `TYPED_RECIPES` (the recipes typed so far) removed, a test checks
     that every recipe of the registry is decorated; the scripting
     tutorial (the parameters typed or as strings, the values refused).
     The bounds decided by Guido before the code changed: the upper
     bounds of the fits within those of the fitter (max_logN in (10, 23],
     max_b in (1, 300]); logN, b, btur of a component within them (10-23,
     1-300, 0-100), its z above -1, its resol positive; score_threshold in
     [0, 1]; 0 allowed where it means "none" (merge_dv, min_b,
     aic_penalty, min_dv, max_components, patience), counts at least 1;
     group_depth at least 0, not 1 as proposed: 0 is the "iterative
     single" strategy of `refit_all`, which `metals_auto` uses (found
     while applying the bounds). The battery of Absorbers
     (`tests/test_recipe_inputs_absorbers.py`, a small mock with a CIV
     doublet and two Ly-a lines), run on master before the change: 1071
     cases, 446 changed, all on invalid or edge inputs; among them
     `add_linked_system` failed on any z given as a string,
     `identify_lines` and `clean_negligible` on a boolean given as a bool,
     the pipelines went on with a value that could not be read (the
     components returned unfitted), `update_component(series=None)`
     renamed the series `'None'`. `fit_bayesian` is not seeded: its valid
     cases check only that the value is accepted. Found and fixed in the
     conversion: an optional string parameter refused None;
  3. *exceptions* (3a continuum, flux and file; 3b edit; 3c absorbers):
     `RecipeError`, with `InvalidParameter`, `MissingPrerequisite` and
     `RecipeFailed`; the `return 0` replaced and the callers adapted (the
     workers of the GUI, the scripting, the MCP server); the count of
     `tests/test_broad_excepts.py` goes down. Decisions of Guido at the
     start: the errors of the core that mean "this data cannot be used"
     (`ValueError`, `KeyError`, `OSError`) become `RecipeFailed` at the
     boundary of the recipe, the decorator `typed_recipe`, and the broad
     `except` clauses of the recipes go (the others, usually bugs, go on
     with a note naming the recipe); a recipe raises and does not log, the
     error is logged once where it is caught (the GUI, the worker of the
     scripts, the MCP server). Every error names the recipe and the reason
     (Guido's question): the decorator records each recipe the error goes
     through, so that the message of a step of a pipeline names both
     ("Populate components from identification → Auto-estimate continuum
     (single-click): ..."). Done in 3a: `recipes/errors.py`
     (`RecipeError`, `InvalidParameter`, `RecipeFailed`,
     `MissingPrerequisite`), the decorator as the boundary, the callers
     (the worker of the GUI shows the message, the worker of the scripts
     logs it without a traceback, the MCP server passes it to the agent,
     the export of the GUI and the continuum editor show it), the recipes
     of continuum, flux and file (their 11 broad excepts removed;
     `fit_powerlaw` without z_em raises `MissingPrerequisite`;
     `export_ascii` raises instead of returning False, which now means
     "nothing to export"), `_with_continuum` lets the failure of the
     estimate go on (the fits of absorbers still turn it into 0 and a
     log line until 3c, the message already names the step). The batteries
     record the outcome of each failure (`tests/recipe_outcomes.py`:
     REFUSED, FAILS, RETURNS_0 until 3b and 3c, RAISES) and check that the
     message names the recipe; no valid case changed. A trade-off of the
     decision: a `KeyError` that is a bug reaches the user as an error of
     the recipe (its cause attached). Found: F101 (the script console).
     Done in 3b: the recipes of edit (33 `return 0` and its 6 broad
     excepts): the parameters read by the recipes themselves raise
     `InvalidParameter` naming them (`resol`, `xunit`, `yunit`,
     `stitch_wavelengths`, `manual_factors`, `equalize_ranges`), the
     other failures `RecipeFailed` with the texts of before; an expression
     that cannot be evaluated (apply, mask, split) fails with the
     expression in the message (a `SyntaxError` reached the GUI as a bug);
     `extract_preset` fails with the reason of each region when none can
     be extracted. The battery of edit: no valid case changed, the 50
     RETURNS_0 now exceptions, the last RAISES (`split('λ >')`) a
     RecipeFailed. Done in 3c, absorbers (38 `return 0`, 13 of its 15
     broad excepts; the 2 of the cleaning of the narrow lines stay, phase
     6, Guido), with Guido's decisions: nothing to do is not a failure
     (`refit_all`, `metals_auto`, `clean_negligible` without components:
     the session back, a WARNING; the pipelines that find no candidate go
     on with it); a final fit that fails fails the pipeline (its 0 made
     the pipeline give back the components unfitted, found by the battery
     of 2c); without components or a selection, RecipeFailed with the
     reason. The checks of "no components" looked at the list object,
     always true, so the recipes ran on an empty list: they look at the
     components now. `_remove_narrow`: nothing left to refit is checked
     first; a refit that fails goes to the cleaning boundary (a WARNING,
     the result kept) instead of removing the lines without a refit,
     silently. `detect_anchor` without the GUI failed with an
     AttributeError hidden by its broad except. No recipe returns 0 any
     more: the workers of the GUI and the MCP server treat a result of 0
     or None as a bug; `RETURNS_0` removed from the batteries; the 26 rows
     of the recipes in the inventory of phase 6 marked as removed. The
     battery of absorbers: no valid case changed;
  4. *structured result* (4a, 4b): `RecipeResult` (the session, a summary,
     what changed, listed and not only counted, flags, warnings, the steps
     added to fill a prerequisite, and the empty fields for phase 7: the
     uncertainty of the continuum, the credibility of an identification
     and its alternatives); with it, the candidate systems of
     `identify_lines` keep their region and score (the rest of F49), and
     the flag `not_fitted` of the groups that `refit_all` could not fit.
     Guido's choice at the start: what changed by a generic comparison of the session before and after (one
     place for every recipe), the warnings captured from the log as text
     (codes can come later without changing the form). Done in 4a:
     `recipes/result.py`, `RecipeResult` (session, branching, changes,
     the values used, warnings, steps, value, extra, and the empty fields
     of phase 7), `compare_sessions` (components added, removed and changed
     with their parameters, columns added, removed and changed, the
     properties: name, object, pixels, z_em, resolution, units; on another
     grid the pixels change, not the columns), `run` (the recipes run as
     steps collected by the decorator while `run` watches); the worker of
     the GUI runs the recipes through `run` and the window keeps the last
     result (the status bar unchanged); the scripting tutorial. Found: the
     description of *Smooth spectrum* said the error is not smoothed (it
     is), fixed. Mutation testing of `params.py`, `errors.py`,
     `registry.py` (Guido's question on the tests adapted with the code;
     run by hand, not in CI): 154 of 359 mutants survived, mostly in the
     code that builds the registry at import, which no test checked: a
     test now builds it again and compares, another reads every key of an
     entry, and the messages of values that are neither strings nor
     numbers are checked; 18 left, equivalent (`bool(get(k, None))` and
     `bool(get(k, False))`, `copy` and `deepcopy` of a default copied again
     when read) or out of reach of mutmut (the decorator applied at
     import) or the message of a duplicate recipe name;
     Done in 4b: the candidate systems of `identify_lines` keep their
     region and score (meta `region_map_json`, `score_map_json`, beside the
     maps of the identification, and `extra['candidates']`,
     `extra['regions']` of the result, with the x range of each region);
     `refit_all` flags `not_fitted` the components of a group whose fit
     failed (they keep their values; `extra['not_fitted']`; a cue in the
     System Inspector), as agreed in phase 6; the result says why each
     component was removed (`extra['removed']`: deleted, too weak by a
     criterion of *Clean negligible*, narrower than `min_b`, logN at the
     lower bound of the fit, replaced by the components of a region) and
     which components a pipeline made and dropped itself
     (`extra['discarded']`). Found: F102. Done in 4c: F102 for `lya_auto`
     (Guido); what the user should know of a result (Guido: everything the
     user does well to know is reported): `RecipeResult.notes` and
     `counts` (the components removed and why, discarded, not fitted, the
     candidates), counted in the status bar of the GUI, in full in its
     tooltip and in the log;
  5. *prerequisites*: *requires* and *provides* in the registry, the
     resolver in the core, the policy of the scripts, the dialog of the
     GUI generalized to any prerequisite, a log entry for each step
     added (possibly two PRs: core and scripts, then the GUI); the slow
     tests before the merge. Guido's choices at the start of the PR: the
     policy of the scripts set for a block (`with policy("raise")`) or a
     call (`run(..., on_missing="raise")`), the signatures of the recipes
     unchanged; with `provide`, z_em and the resolution, which no recipe
     can provide by itself, left to the recipe as before (its WARNING);
     the GUI asks for the continuum (estimate with the defaults, open
     Auto-estimate Continuum, cancel); the step added is a log entry of
     its own, before the recipe; the identification and the pipelines
     need the continuum (F103). Done in 5a: `recipes/prerequisites.py`
     (`PREREQUISITES`, `missing`, `policy`, `fill`, `ProvidedStep`), the
     keys `provides` and the requirement `continuum` in the registry
     (13 recipes need it: those that estimated it themselves, those that
     fit, the identification and the pipelines; every recipe that fits
     needs the resolution too, F99), the prerequisites filled by
     `typed_recipe` before the recipe (only the outermost provides: the
     steps find them), `_with_continuum` removed; the estimate said at
     WARNING (it was INFO); `RecipeResult.provided`, a state and a log
     entry of their own in the GUI and in the script console; a region
     of `refit_all` without z_em is a `MissingPrerequisite`. Checked
     with a spectrum in physical units: the continuum estimated first
     gives the same numbers as before for the recipes that estimated it
     themselves (`add_component`, `fit_component`, `optimize_system`,
     `refit_all`, `fit_bayesian`, `lya_auto`, `forest_auto`), different
     ones, as agreed, for `identify_lines`, `components_auto`,
     `doublets_auto` and `metals_auto` (F103). Done in 5b: the GUI asks
     through the resolver (`missing`): z_em and the resolution as before
     (Set Properties), the continuum with three choices (Estimate, with the
     defaults of Auto-estimate Continuum; Open Auto-estimate Continuum…;
     Cancel); the pending recipe resumed after the recipe that provides what
     it lacked (`_pending_recipe`, with `after`), dropped when that recipe
     fails or its dialog is cancelled;
  6. *descriptions and names* reviewed (F48, the texts too short for a
     model), the flag of the parameters not offered to agents (column
     names, plots), the parameters of the mask of the pipelines (F50) and
     of the margin on z_em (F55), with defaults equal to the present
     behaviour. In two PRs: 6a the texts and the flags, 6b the new
     parameters. Guido's choices: every recipe is offered to agents, also
     the basic ones called only by other recipes, which are blocks of new
     recipes (no flag of the recipes; "interactive" not needed:
     `update_from_knots` takes its knots as lists); not offered only the
     parameters of debug and the questions of the GUI; all the texts
     reviewed in one PR, no name changed. Done in 6a: the key `agents` of
     a parameter (`debug_rating`, `bypass_scoring`, the four
     `renorm_model`); the details of every recipe and a sentence for every
     parameter (none described by its name any more), errors fixed (F48,
     `mask_unabs`, `import_ascii` said it reads the system list, the
     meaning of `group_depth`, of `order` of *Equalize*); `stitch` declares
     its parameters (`other_sessions`, `sort`), so that it can be offered;
     a test that every recipe and parameter is described. Found after
     the merge, while rebuilding the screenshots: with the longer texts the
     dialogs grew wider than the window (the descriptions are their
     labels); fixed in PR 6a-bis (labels and tooltips that wrap, a test on
     the width of every dialog, the five screenshots of dialogs rebuilt);
     done in 6b: the parameters `mask` (`'new'`, `'session'`; F50) of
     `components_auto`, `doublets_auto`, `lya_auto` and `z_em_margin_kms`
     (F55) of `identify_lines` and of the pipelines, at the end of their
     signatures, with the defaults of before (results unchanged). Guido
     asked whether the parameters and the interfaces of the recipes could
     be adapted to each procedure, as for the editing recipes: in 6-bis the
     kinds of the parameters in the registry (PR 6c: fixed choices, a
     session, a column), so that the GUI shows lists and the agents the
     values allowed; the review of the parameters recipe by recipe later,
     in phase 6-quater for the recipes that phase 7 does not rewrite (named
     after 2.5.0, Guido), with phase 7 for the others; also (Guido) the collections of elements to
     fill some parameters, as the buttons of the columns and sessions of
     the editing recipes, in phase 7. Done in 6c: the keys `choices` (the
     values allowed: `region_limit` with `all_ly_forest`, which was
     accepted and not said, `mask`, `engine`, the `param` of the
     constraints; a list in the GUI, an unknown value refused, the case not
     checked) and `kind` (`'session'`, `'column'`) of a parameter; the
     panel of the sessions of a dialog shown for every recipe with a
     parameter of kind session (`detect_anchor` too), its buttons fixed
     (F105);
  7. F12, the other half: the custom series of `add_component` kept with
     the session instead of `STANDARD_MULTIPLETS`. Done: no table is needed,
     since a custom series (`'FeII_2586,FeII_2600'`) carries its
     transitions in its name (`atomic_data.series_lines`, `is_multiplet`,
     used by the 19 readers of the series of a component, in core,
     fitting, recipes and GUI); the two writers (`add_component`, the
     loaders of `io/adapter.py`) removed. The readers that go through the
     whole table (the siblings of a transition, the pairs of doublets of
     the *Identifications* window) no longer see a custom series there, as
     they should not;
  8. the MCP pilot moved onto the new API (still a thin layer: the test
     that the phase serves the full server), then the closing of the
     phase. Done: the recipe tools run through `run` (typed parameters,
     the warnings, `summary` and `notes` of the result in the output); the
     prerequisites said, not filled (a tool that cannot run says what is
     missing and who provides it; `provide_missing=true` provides it with
     the defaults, reported in `provided` and logged as a step of its
     own); the new parameters (`mask` as a list, `z_em_margin_kms`);
     `mcp_server/schemas.py`, the JSON schema of the parameters of any
     recipe from the registry (types, defaults, bounds, choices, kinds,
     the parameters offered to agents), tested on all 48 recipes: the
     generator of the full server; the tutorial for agents (it said that
     `doublets_auto` uses the mask of the session: F50, it made a new
     one).

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
- Numbers: the snapshots and the mock benchmark unchanged; changed on
  purpose, with Guido's decision: F102 (`lya_auto` no longer refits the
  red side) and F103 (the identification and the pipelines on a spectrum
  without a continuum: it is now estimated first); the slow tests run
  before the merge of the PRs that touched the pipelines.
- Found and fixed: F49 (the rest), F50 and F55 (their parameters), F102,
  F103, F104 (a pending recipe left after a failed Set Properties), F105
  (the buttons of the sessions of three dialogs did nothing), the dialogs
  wider than the window with the new texts; `stitch` without its
  parameters in the registry; several recipes that failed on values a
  script or the GUI could give (#710-#712). Left: F101 (the script
  console), F102 for `forest_auto` (its region, Guido's decision).
- Decided for later: the review of the parameters recipe by recipe (phase
  6-quater, for the recipes that phase 7 does not rewrite; phase 7 for the
  others), the collections of
  elements to fill some parameters (phase 7), the user configuration of
  the defaults (after the full MCP server).
- Tests 1308 → 3288 (mostly the cases of the batteries of inputs, one
  test each); coverage 84% → 86% (core 92%, GUI 85%, recipes 95%, the MCP
  server 99%; `io` 50%, the loaders without sample files).

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
| 10 | Review of the manual fits (`add_component`, `update_component(s)`, `delete_component`, `detect_anchor`, `fit_component`, `update_constraint(s)` with F117, `add_linked_system`, `optimize_system`, `fit_bayesian`); possibly split in two | Opus | 3-5 h (two PRs if split) (split, Guido, 12 November 2026: PR 10 the components and the constraints, done, 12 November 2026, about 3 h: see below) |
| 10-bis | Review of `optimize_system` and `fit_bayesian` (split from PR 10, Guido, 12 November 2026): the unknown components refused as in PR 10 (F142), `uuids` of `fit_bayesian` required, the GUI no longer passing `"500"`, `"250"`, `"2"` as strings, the description of `z_window_kms` of `optimize_system` (the window analysed and the freedom of z in its fits); the samplers: F139 (the crash of nested sampling), F140 (the posterior of nested sampling without its weights), F141 (MCMC with more than 16 free parameters: its 32 walkers are a heuristic, to ask Guido before changing them), the minimum of `nlive`. The numbers of nested sampling change (slow tests, a measure). F28 and F114 not touched | Opus | 2 h (done, 12 November 2026: see below) |
| 11 | F107: the labels of the menus of the GUI from the registry (and the screenshots), with F125 (the pages of `docs/` for the *Help* button) | Sonnet | 1 h (done, 12 November 2026, about 2 h with the tests, the documentation and the screenshots: see below) |
| 12 | MCP server: the descriptions of the shared parameters shorter, a paragraph of typical workflows in the instructions, the banner of a followed session with the pid and the version of the marker | Sonnet | 1 h (done, 13 November 2026, about 2 h with the tests, the screenshot and the documentation: see below) |
| 13 | Closing of the phase (added on PR 8, Guido): "Result of the phase" in this section, the progress table, "Rules and decisions in force" (what of the phase still holds: the aliases, the one set of names of the columns, the units of the flux), `CLAUDE.md`; the concluded detail moved to `ROADMAP_HISTORY.md`; before it, the review with Guido of the findings without a phase (added on PR 10-bis, Guido: F15, F17, F28, F30, F52, F101, F102, F112, F114, F115, F116, F126, F128; a phase for each, or closed: F52 and F116 wait for #517, closed in phase 4) | Sonnet | 1 h |
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
  - ~~`refit_all` and `clean_negligible`: in the manual fits (PR 10) or in
    phase 7, which rewrites `refit_all`~~: phase 7 (Guido, PR 10);
  - ~~F101 (the `app.run_recipe` of the script console fails with a
    TypeError): in PR 10, or with the decision on the console~~: with the
    decision on the console (Guido, PR 10);
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
- *Result of PR 12* (MCP server: shared parameters, workflows, banner; 13
  November 2026; the table of the shared descriptions and the tests agreed
  with Guido in the conversation of the PR). (1) **The descriptions of the
  four parameters of every recipe tool** (`session_id`, `provide_missing`,
  `max_listed`, `reason`) are shorter: 34 018 characters in the 65 tools
  (the `session_id` of the sessions' tools included) became 22 522 (-34%),
  of the 51 078 of all the descriptions of parameters; no constraint lost
  (the default of `provide_missing`, that the counts are complete, that a
  recipe that only returns a value logs nothing). One source:
  `SHARED_PARAMS` in `mcp_server/schemas.py`, read by `SERVER_PARAMS`
  (`tools.py`) and by `SESSION_ID` (`session_tools.py`), which had a second
  text. Tests (`test_mcp_schemas.py`): one text in every tool, at most 160
  characters per shared description, 24 000 in all over the 48 recipe tools
  (21 360 now), and at most 300 for any description of any tool (the longest
  is 257); the registry's own descriptions (`z_window_kms`, `max_b`, ...) were
  not shortened: each says something of its recipe, and appears in 1-10 tools.
  (2) **A paragraph of workflows** in the instructions of the server
  (`INSTRUCTIONS_WORKFLOWS`, 1 860 characters): which tools of *Absorbers*
  fit (`add_component`, `add_linked_system`, `fit_component`, then
  `refit_all`, `optimize_system`) and which do not (`update_component(s)`,
  `update_constraint`, `update_constraints_column`), and three workflows with
  examples (adapt a system, try a component, link ions); the examples were
  run on the mock, and a test reads them (tools, parameters, required ones)
  from the registry. The advice of the trials of section 4.3 (undo to drop
  alternatives, `session_duplicate` to try) is in it. Found, and not
  changed: F152 (the update recipes do not redraw the model). (3) **The
  banner of a followed session** says the pid and the version of Astrocook
  of the process that writes the archive. The marker already held them (PR
  5-bis of 6-ter): the GUI did not show them; `FollowMarker.pid` is now
  optional, because `read_marker` discarded a marker without it, so one of an
  earlier version is read and its banner is as before. The scene
  `agents_followed_session` of the screenshots writes a fixed marker (pid
  4121, 2.6.0) so that the image does not change at every run. Numbers and
  behaviour of the recipes unchanged: only texts, the instructions and a
  banner (no recipe, fit or numerical code touched), so neither the slow
  tests nor the benchmark apply. Tests run here: `test_follow`,
  `test_gui_follow`, `test_mcp_*` (but the slow ones), `ruff check .` and
  `mypy`.
- *Result of PR 11* (the labels of the menus from the registry, F107, and
  the *Help* button, F125; 12 November 2026; the table of the entries and of
  the pages agreed with Guido in the conversation of the PR). A new key
  `menu` of the recipes in the registry, `"Menu > Label with one &"`
  (`"Edit > Set &Properties"`): the **only** place where the label of an
  entry is written (Guido's question: does a second text not break the
  single source of truth? Not if `main_menus.py` writes none). The GUI
  builds the actions of the recipes from it (`MainMenusMixin._recipe_action`:
  the text with its `...` and mnemonic, the tooltip, the attribute
  `<recipe>_action`, an error if the registry puts the recipe in another menu
  than the one the layout builds), the order and the separators stay in
  `RECIPE_LAYOUT` (names of recipes only); the dialog is titled by the
  label (`RecipeSpec.title`), and the methods text names the recipe by it
  (`RecipeSpec.menu_path`), instead of deriving `Category > brief`. The
  `brief` stays what it was (log, status bar, MCP descriptions): a sentence,
  not a label. The 27 recipes of the menus have a `menu`; the others
  (`x_convert`, `trim_common`, `components_auto`, ... run from the plots, the
  System Inspector or scripts; `forest_auto` and `metals_auto`, whose entries
  are commented out of the menu, kept as actions in `RECIPES_OFF_MENU`: Guido,
  after the PR, took their `menu` away until they come back) have none, and the methods text names no menu
  for them (it named `Edit > Convert X axis`, an entry that does not
  exist; and `Edit > Import Systems from Session` for an entry of File).
  The labels did not change (the documentation, the tests and the
  screenshots quote them); `Auto-detect Doublets...` now opens *Recipe:
  Auto-detect Doublets* (it was *Auto-fit metal doublets (Red Side)*).
  F125: the `url` of every recipe (all have one, 11 had none or a V1 name)
  is a page of the documentation relative to its root, with the section
  (`tutorials/editing_data.html#masking-data`); `RecipeSpec.doc_url` joins it
  to `registry.DOCS_URL`, which *Help > Online Documentation* uses too; the
  *Help* button opens it. `test_every_recipe_points_to_a_section_of_the_docs`
  checks the page and the section (the id Sphinx gives to the title) of
  every recipe. Four recipes have no section of their own and point to the
  nearest (`trim_common`, `forest_auto`, `metals_auto`, `components_auto`);
  the sections exist on the site only from the next release (the root of
  the documentation is the latest release: until then the button opens the
  page, at its top). Renamed, with the same meaning, the attributes of the
  window `auto_cont_action`, `find_abs_action`, `fit_cont_action` and
  `calibrate_action` (now `estimate_auto_action`, `find_absorbed_action`,
  `fit_continuum_action`, `calibrate_from_magnitudes_action`: one name rule,
  `<recipe>_action`). Added to the tutorial *Getting Started* a section on
  the dialog of a recipe and its *Help* button. Numbers unchanged: only
  texts, labels and addresses (no recipe, fit or numerical code touched), so
  neither the slow tests nor the benchmark apply. Found: F151 (the tooltips only: the messages that describe what a recipe did keep the `brief`, decided with Guido). Tests run
  here: `test_recipe_registry`, `test_main_menus`, `test_methods_text`,
  `test_gui_structure`, `test_recipe_dialog_current`, `test_gui_smoke`,
  `test_gui_dialogs`, `ruff check .` and `mypy`.
- *Result of PR 10-bis* (review of `optimize_system` and `fit_bayesian`;
  12 November 2026; the table of their 15 parameters and the decisions A-H
  agreed with Guido in the conversation of the PR): no parameter renamed;
  two parameters added (`seed` of `fit_bayesian`, `max_nfev` of
  `optimize_system`, whose default is the value that was fixed); `uuids` of
  `fit_bayesian` required. The numbers of `optimize_system` and of the MCMC
  do not move (seven valid calls taken on master before the change, the
  MCMC with the global generator of NumPy seeded: snapshots
  `optimize_bayesian_*.npz` of `test_valid_results_unchanged`, the
  tolerances of PR 10 for `optimize_system`; the MCMC to the bit here, but
  on the CI, Python 3.12, the chain of 500 steps took another move and b
  came out 0.0013 km/s away, 0.01 of its error: the MCMC is compared within
  0.2 of its errors, the errors within 20%, the model within 1e-3); those of nested sampling are new (it crashed, F139)
  and weighted (F140). No fit of the benchmark, of the automatic pipelines
  or of the slow tests goes through these recipes (none calls them), so
  neither the slow tests nor the benchmark can measure the change: the
  targeted test `test_nested_posterior_is_weighted` replaces them (the
  nested and the MCMC agree within 2 sigma in the values and within a
  factor 2 in the errors; before, the errors of nested were 300 times
  those of the MCMC).
  - *The samplers* (`fitting/bayesian_fitter.py`): F139 (the remaining
    evidence read from the iteration), F140 (weighted medians and standard
    deviations, equal-weight samples for the corner plot), F141 (decision
    A: `max(32, 2 * ndim)` walkers; Guido asked whether fewer walkers, e.g.
    6 for 3 parameters, would be faster: they would, about 5 times, but 2
    per parameter is the floor of emcee, not a choice, the ensemble mixes
    worse, a chain of 300 steps after the burn-in gives 1800 correlated
    samples instead of 9600, and the numbers of every fit would change;
    the MCMC takes 3.4 s on a CIV component, the slow sampler is the nested
    one); a `RuntimeError` of a sampler is a `RecipeFailed`; the warnings of
    dynesty to the log, once each (decision F; with the warnings as errors
    of the tests they failed every nested sampling).
  - *Refused* (`InvalidParameter`, before the sampling): an unknown
    component (F142, `check_uuids`), fewer than 10 live points per free
    parameter (decision B, F150; `NLIVE_PER_PARAMETER`, `min` 10 in the
    registry).
  - *`seed`* (decision D): `RandomState(seed)` for the walkers and the
    moves of emcee, `default_rng(seed)` for dynesty; left out, the global
    generator as before (the same bits).
  - *`optimize_system`*: `max_nfev` (decision E), the descriptions of
    `z_window_kms` (also the freedom of z in its fits) and of
    `max_components` (the candidates tried, F147); F148 and F149
    registered, not changed (decision G).
  - *The GUI*: `_fit_bayesian` passes the components and the engine only.
  - *Docs*: *Absorption Systems* (the prior, the walkers, nested sampling
    and `seed` in scripts), *Using Astrocook from an AI Agent* (`seed`),
    the CHANGELOG. No screenshot: no dialog changed.
  - *Found* (in the register): F147-F150. F28 and F114 not touched.
  - *Asked by Guido*: the findings without a phase (13) to be reviewed
    before the closing of the phase (row 13 of the plan).
  - New `tests/test_recipe_review_optimize_bayesian.py`; the expected
    messages of `fit_bayesian` without components in `test_recipe_errors.py`
    and `test_recipes_absorbers.py` follow (the session is checked first,
    `uuids` required).
- *Result of PR 10* (review of the manual fits, components and
  constraints; 12 November 2026; the table of the parameters of the 9
  recipes and the decisions A-F agreed with Guido in the conversation of
  the PR; `optimize_system` and `fit_bayesian` split to PR 10-bis): no
  parameter renamed; five defaults removed; the numbers of a valid call do
  not move (the values of 24 calls taken on master before the change,
  snapshots `manual_fits_*.npz` of `test_valid_results_unchanged`, the
  values and the model within 1e-4 relative, the errors within 1e-2, and the
  anchor of `detect_anchor` on the mock within 1e-5; on the runners of the
  CI the fits differ from here by up to 2.5e-6 in b and 3e-4 in the errors
  (a b at its bound of 300 km/s), with the same code: an md5 of the values,
  as in PR 9, failed on Python 3.12, a tolerance of 1e-6 on 3.13; here
  3.12 and 3.13 give the same bits), except two cases agreed: the links made without `target_uuid` (F146)
  and `detect_anchor` on a flux not normalized (F145). No fit of the
  benchmark or of the automatic pipelines goes through the changed paths
  (the pipelines do not call these recipes; `add_linked_system` and the GUI
  pass `target_uuid`), so neither the slow tests nor the benchmark were run.
  - *The checks* (`recipes/_components.py`): `check_uuids` (an unknown
    component, F142), `check_series` (a series not in `ATOM_DATA`, F143; a
    list of series to `add_component` names `add_linked_system`),
    `link_target` (F146). `InvalidParameter` before the recipe runs,
    naming the parameter, each characterized before (the probes of the
    conversation; each test says what it was).
  - *Nothing to change refused* (decision F): `update_component(s)` with
    no value, `update_constraint` without `is_free` and expression,
    `delete_component` without components (a state that changed nothing;
    `update_components` gave the same session). The batteries of
    `test_recipe_inputs_absorbers.py` that passed one value None now pass
    another parameter with it, so that they still test the conversion.
  - *F117 closed*: the expression that cannot be evaluated is an
    `InvalidParameter` (`RecipeFailed` from the core before).
  - *F146* (decision C): `target_uuid` is the first component the
    expression names when not given (it decides the group of the fit,
    `system_list.py`), refused if the expression does not name it; a link to
    the component itself and an expression with `is_free=True` refused; a
    link stored with `is_free=False` (the fit treats it so: only the stored
    flag changes).
  - *F145* (decision D): `Spectrum.detect_doublet_z` on `y / norm` (the
    same anchor on the mock, whose continuum is within 0.5% of 1, to the
    bit), the axis converted to nm (it was read as nm), a `ValueError`
    without a common range (0.0) or without a continuum on a flux not
    normalized.
  - *Defaults* (F144): `z` of `add_component` and of `add_linked_system`,
    `series_list`, the transitions of `detect_anchor`, `is_free` of
    `update_constraints_column`: required. A number or a bool is required
    in the registry by a default None that is not `optional`; the dialog
    shows it empty (`recipe_dialog.py`); `add_component` keeps `z=None` in
    its signature (it follows `series`, which has a default) and refuses it
    in its body. `update_components.series` in the registry,
    `Optional[str] = None` (the `'None'` of the logs is None);
    `NOT_IN_REGISTRY` loses its row.
  - *`resol`* (decision E): 0 removes the resolution of the component
    (`min 0`, `gt 0` before); the docstrings said "R or FWHM", it is R
    since #526.
  - *Strings*: `add_component` and `add_linked_system` pass `z_window_kms`
    to `fit_component` as a number, and so does the System Inspector (its
    field of the fit window); the log records the converted value either
    way. The strings of `_fit_bayesian` and `_refit_list` stay (PR 10-bis,
    and `refit_all` in phase 7).
  - *Docs*: *Absorption Systems* (the resolution of a component, "Linking
    in Scripts and for Agents", "Doublet Anchor Across Two Sessions"),
    *Using Astrocook from an AI Agent*, the CHANGELOG. No screenshot: no
    scene shows the dialogs changed.
  - *Found* (in the register): F139-F141 (the samplers, PR 10-bis; F139 a
    crash, reported at once), F142-F146 (fixed here).
  - *Logs*: a session log with a step that the recipes now refuse (e.g.
    `update_component` of a component not in the session, which changed
    nothing) stops its replay at that step, with the message; the logs of
    the GUI cannot hold one (it names the components of the table).
  - New `tests/test_recipe_review_absorbers_manual.py`.
- *Result of PR 9* (review of edit on several spectra; 12 November 2026;
  the table of the 27 parameters of the 9 recipes and the decisions agreed
  with Guido in the conversation of the PR): one parameter renamed, two
  defaults removed, one type changed; the numbers of a valid call with the
  sessions in one unit do not move (fingerprints taken on master before the
  change, `test_valid_results_unchanged`), those with units mixed are fixed
  (F132-F134 and F138, agreed); no recipe of the continuum, the fits or the
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
    a wavelength refused (`x_convert` first). The rule is in
    `recipes/_combine.py`. F138 (`flux.equalize` and the core of
    `equalize_to_reference` and `coadd`, the same faults), found in the
    PR and fixed before the merge (Guido): the core compares in the units
    of the spectrum scaled (a script gets the right numbers too),
    `flux.equalize` keeps the unit of its flux (an arbitrary one takes
    that of the reference); fingerprints of `equalize` (orders 0, 1, -1)
    and of the core `coadd` taken before.
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
    and `in_units_of` of `recipes/_combine.py`), covered by the
    fingerprints.
  - *Kept*: `append`, `refit` of the imports (`refit_all` with 20 km/s and
    `max_nfev` 100, fixed), `sort` of `stitch`, `alias_map` of `split`
    (GUI only), `'auto'` of `equalize_ranges` and `manual_factors`,
    `window_kms`, `kappa`, `dx`.
  - *Docs*: the tutorial of Edit (the sessions by name or as objects,
    the units of the combination, the limits above), of Flux (the units
    of *Equalize*), *Using Astrocook from
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
| Phase 4 | [#526](https://github.com/DAS-OATs/astrocook/issues/526) | resolution as a single column and exact convolution (F8, F22, rest of F10, F9, F14, F47); design and decisions in the issue; PR 1 (single source, migration, no stamping, #615), PR 2 (GUI on the same resolver, #616), PR 3 (exact convolution, `fitting/lsf.py`, #617) and PR 4 (no more writes of the legacy sources; the profiles of the System Inspector with the physics and the convolution of the fitter, `convolve_flux` removed) done |
| Before release 2.0.2 (decided with Guido) | [#534](https://github.com/DAS-OATs/astrocook/issues/534) | which chi2 the fit uncertainties are scaled with (follow-up of F29): done, option (a) |
| Phase 4 | [#517](https://github.com/DAS-OATs/astrocook/issues/517), [#518](https://github.com/DAS-OATs/astrocook/issues/518) | from the findings review: flag and refit components at the parameter bounds, with F21 (`lya_auto` output); editable z window in the manual fits (F4). #517 in three PRs (agreed with Guido): flags on the components (PR 1), cleanup of the pipeline output (PR 2), GUI cues and "Select flagged" (PR 3) |
| Phase 3-bis (brought forward) | [#559](https://github.com/DAS-OATs/astrocook/issues/559) | found while cleaning up the Ruff warnings (phase 2-bis): `dbtur` not updated after a fit (F36); done |
| Phase 3-bis | [#560](https://github.com/DAS-OATs/astrocook/issues/560) | found while cleaning up the Ruff warnings (phase 2-bis): `norm` of the `rebin` recipe ignored (F37); done: no longer offered, a warning if set (Guido's decision) |
| Phase 4 | [#568](https://github.com/DAS-OATs/astrocook/issues/568), [#573](https://github.com/DAS-OATs/astrocook/issues/573) | found in phase 3: a single table of atomic parameters (F38, F39); recipe defaults differing between GUI and API (F40) |
| Phase 3-bis (priority) | [#575](https://github.com/DAS-OATs/astrocook/issues/575) | found in phase 3: a recipe that opens a new session is also recorded in the log of the parent session (F41); done |
| Phase 5-bis | [#630](https://github.com/DAS-OATs/astrocook/issues/630) | errors of linked parameters written as 0 in the ASCII list of the systems (reported by a user, October 2026); to be fixed before 2.3.0 (Guido's decision); done in PR 4 of the phase: the error of a linked parameter is propagated from the free parameters through the expression of the link (numerical gradient and the covariance of the fit, `VoigtFitterV2._linked_errors`), in the fits and in `evaluate_uncertainties`; 0 when its targets are frozen or outside the fitted group |

### 3.2 Findings register
| ID | Area | Finding | Impact | Recorded by | Status |
|---|---|---|---|---|---|
| F1 | `rebin_spectrum` | Errors do not decrease with the number of pixels per output bin (σ instead of σ/√N) | numerical | `test_rebin_weighted_average_error`, `test_rebin_errors` (V1 parity; regression) | fixed (phase 1-bis): error of the weighted mean as in V1, sqrt(sum w^2 sigma^2) / sum w; issue [#490](https://github.com/DAS-OATs/astrocook/issues/490) |
| F2 | `rebin_spectrum` | Sigma clipping uses the first-pass bin mean, which includes the outlier: with few pixels per bin (e.g. 3) a single spike can make the whole bin NaN. Clipping runs once, docstring says "iterative". Not a pure V2 regression: V1 clips with `sigma_clip` within each bin, which keeps the bin but cannot remove an isolated spike with 3-5 pixels per bin; neither version removes it. Related to #465 | numerical | snapshot `rebin_spectrum_nm_clip`; `test_rebin_clipping_keeps_bins_with_a_spike`, `test_rebin_clipping_keeps_narrow_absorption` (regression) | fixed (phase 1-bis, agreed with Guido): single-pass clipping against a running median of the input pixels, upward only (a two-sided test removed the cores of narrow lines: equivalent width -20% for a 2.5-pixel line; the previous V2 clipping up to -100%), a bin is never emptied; V1 behaviour dropped |
| F3 | `fit_continuum` | `fudge='auto'` (recipe default) is always exactly 1.0: computed on residuals that are zero by construction, before the smoothing. Related to #465 | numerical | snapshots `continuum_j0103_*`; `test_v2_auto_fudge_is_one` (characterization) | accepted (phase 1-bis, decided with Guido): not a regression in substance. The V1 fudge corrects the running mean V1 builds the continuum from; V2 interpolates the unabsorbed pixels, so the fudge is 1 by construction. Local fudges after the smoothing were tried on J0103-1305: none reproduces V1 (effect ±1% or more vs ±0.2%, no correlation), they would correct the bias of the final smoothing instead. Documented in `fit_continuum`; the smoothing bias goes with #465, F17, F30 |
| F4 | `VoigtFitterV2` | z is bounded to ±`z_window_kms` (default 20 km/s) around the initial guess; a guess further away cannot reach the true value. Intended? | behaviour | `test_fitting.py` (guess inside the window) | fixed (phase 4), [#518](https://github.com/DAS-OATs/astrocook/issues/518): default unchanged (20 km/s when adding, 100 when refitting); the System Inspector has a **Fit ±** field (after logN and b) for its manual fits (`test_wider_z_window_reaches_a_distant_guess`, `test_inspector_fit_window`) |
| F5 | `convert_y_axis` | Annotation declares 2 return values, the function returns 3 | code quality | snapshot `convert_y_axis` | fixed (phase 1-bis): annotation and docstring list the three returned values |
| F6 | architecture | `core/spectrum_operations.py` imports `gui.debug_utils`: core depends on the GUI | architecture | `test_debug_hooks.py` (`test_import_does_not_load_the_gui`) | fixed (phase 3-bis): the core requests the debug plots through a hook (`core/debug_hooks.py`) that the GUI connects to its signal |
| F7 | GUI | `gui/main_window.py` imports `matplotlib.pyplot` before `pyside_plot.py` selects QtAgg: importing the GUI without a `QApplication` (and no display) fails | architecture | `test_gui_imports_after_pyplot` (regression) | fixed (phase 5): reproduced as "pyplot imported first, without a QApplication (a script, a test), then the GUI": `matplotlib.use('QtAgg')` at the import of `pyside_plot.py` failed on the headless backend already chosen. The GUI modules no longer choose the backend of pyplot nor import it (colormaps, rcParams and style from `matplotlib`; the canvases are QtAgg anyway); the launcher sets QtAgg for pyplot (the debug plots) once the `QApplication` exists. The workarounds of the tests (the GUI imported inside the fixtures) removed; the images of the tutorials unchanged |
| F8 | resolution handling | The resolution lives in five places (per-component value, `resol` column, `SpectrumDataV2.resol`, meta `resol`/`RESOL`, meta `resol_fwhm`), each read with a different priority by the fitter, the `.acs2` loader, coadd and the GUI. The `.acs2` save does not write the attribute: it survives only through the meta keys (set by "Set properties"). The fitter stamps the resolution used into every fitted component, which then overrides later changes. R vs km/s is guessed from a threshold of 500. The J0103-1305 archive (older version) stores no resolution at all. To be analysed by Guido. Found in phase 3-bis: `SpectrumV2.rebin` computes the median of the rebinned `resol` column and stores it in `new_meta['resol']`, but then rebuilds `new_meta` from the original metadata, so the meta key keeps the value from before the rebin while the attribute is updated (45000 vs 45008 on a mock with R = 40000-50000) | behaviour, possible data loss | — (map in the PR that added this row) | fixed in phase 4 (#526, PR 1): the `resol` column is the single source (`SpectrumV2.resolution()`, R per pixel; unit from the column: dimensionless R, velocity FWHM); the attribute and the meta keys are migrated into it when a spectrum is built (`core/resolution.py`) and no longer read; the guess of the unit by value (500) survives only in that migration, for data saved before; the fitter no longer stamps its resolution into the components; components keep an override in R (legacy FWHM converted, stamped values within 1% cleared at load: Guido's decisions (a), (b)). PR 4: nothing writes the legacy sources any more (`set_properties`, `rebin`, `coadd`, `split`, the RESVEL loader, the V1 adapter, the mocks write the column only); the migration clears the attribute; the meta keys are left as they came, as they may be keywords of the original file |
| F9 | Voigt model | Components narrower than a pixel (e.g. Ly-alpha with b = 0.5 km/s, logN = 17.8 in J0103-1305) make the model depend on the sampling; with R = 45000, 1% of the 570-700 nm window still differs from the model saved in the session. Related to the oversampling check (item 9): the V2 model, sampled at pixel centres, overestimates the equivalent width of CIV by 4% at b = 1 km/s and 0.3% at 3 km/s (2.5 km/s pixels, R = 45000) | numerical | — | fixed in phase 4 (#526, PR 3): the model is computed on a fine grid (step a third of the narrowest b in the window, at least 0.5 km/s, and half the instrumental sigma) and averaged over each pixel; against the exact flux of mocks, 1.5e-5 to 3e-5 with F61 corrected (`test_lsf.py`). Timing (`bench_fitting.py`): single fit 0.008 → 0.010 s, group 1.44 → 1.65 s, chain on a mock 19.2 → 14.3 s, on J0103-1305 373 → 262 s (each line evaluated only where tau > 1e-6) |
| F10 | `convolve_flux` | With a per-pixel resolution column, pixels are grouped by resolution rounded to 100 (R) or 0.5 km/s, but selected by the raw value: pixels whose resolution is not already rounded are never convolved and the model is 0 there (a flat column at R = 45050 gives a zero model everywhere). Affects rebinned and coadded spectra (coadding three synthetic exposures at R = 42000-47500 gives a non-round column on 90% of the pixels and a zero Voigt model on 90% of them, `test_voigt_model_on_coadded_spectrum_is_not_zero`), RESVEL input, non-round "Set properties" values. In addition, even for a constant R the kernel width is computed at the median wavelength and with the median pixel size, so it is wrong far from the centre of the range (+53% / -26% at the ends of 380-780 nm). Proposed: a single exact convolution in velocity space, after the safety net is complete | numerical, wrong results | `test_model_column_resolution_is_not_zero`, `test_voigt_model_on_coadded_spectrum_is_not_zero`, `test_voigt_model_with_resolution_column` (regression) | selection bug fixed (phase 1-bis): pixels selected by the rounded resolution, NaN filled with the unconvolved flux. Kernel width far from the median wavelength and exact convolution: fixed in phase 4 (#526, PR 3, `astrocook/fitting/lsf.py`): each pixel takes the model convolved with a Gaussian of its own resolution and averaged over its edges, with `test_lsf.py` against the analytic result (1e-8) and against mocks with a flat and a varying resolution (5.9e-3 at R 45000 and 5.7e-2 at R 8900 before) |
| F11 | grouping | `VoigtModelConstraintV2.set_active_components` ignores `group_depth`: it always calls `get_connected_group` with the default depth 2, so the recipe option (1 = neighbours only) and the `group_depth=0` used by the refit have no effect | behaviour | `test_set_active_components_honours_group_depth`, `test_fit_component_group_depth_0_freezes_neighbours` (regression) | fixed (phase 1-bis, decided with Guido): `group_depth` honoured, default 2 everywhere; the internal fits of `refit_all`/`metals_auto` keep the depth-2 expansion they always had (`FIT_EXPANSION_DEPTH`), since depth 0 there was 2.75x faster but fitted blends worse (J0103-1305: chi2 +14% in the forest); the user can reduce the depth in "Fit component" and in the Bayesian fit |
| F12 | global state | Module-level dictionaries are modified at run time: the GUI rewrites the parameter docs of `EDIT_RECIPES_SCHEMAS["equalize_and_stitch"]`, and `add_component` registers custom series in `STANDARD_MULTIPLETS`. The changes persist for the whole process (and leak between tests) | architecture | `restore_multiplets` fixture in `test_issue_455_multiplet_subset.py`; `test_recipe_registry.py` | half fixed in phase 6-bis, PR 1: the recipes are described by a read-only registry (`recipes/registry.py`, frozen `RecipeSpec`), the only reader of the schemas, and the GUI no longer rewrites them (F79); fixed in PR 7: a custom series carries its transitions in its name, so it needs no table at all (not even one kept with the session): `atomic_data.series_lines` and `is_multiplet` resolve it, and the readers of the series of a component use them; `add_component` and the loaders no longer add it to `STANDARD_MULTIPLETS`; the fixtures that restored the table removed, the tests check that it does not change |
| F13 | inspector | `SystemTableModel.set_highlighted_uuids` returns early when the selection is unchanged, so link colours are not recomputed when a link is created with the group still selected (model-level reproduction; full GUI not driven) | GUI | `test_link_created_with_group_selected_is_coloured`, `test_link_removed_with_group_selected_is_not_coloured`, `test_link_keeps_its_colour_when_another_is_added` (regression) | fixed (phase 3-bis): the colours are recomputed in `update_data`; a link keeps its colour when others are added or removed |
| F14 | grouping | Two components are considered blended when their line centres are closer than max(4 b, 30 km/s), independently of the instrumental resolution and of the line strength. At X-Shooter resolution (R ~ 8900) AlII 1670 and CIV 1550 still overlap at 60 km/s but are not grouped; worse in the NIR or for saturated lines. Possibly why #469 is still seen on X-Shooter data | behaviour | `test_blend_threshold_*` in `tests/test_issue_469_blends.py` | fixed in phase 4, after PR 1 of [#526](https://github.com/DAS-OATs/astrocook/issues/526) (the single per-pixel resolution), in parallel with its PR 3: the threshold is max(4 b, 30 km/s, k FWHM), with the instrumental FWHM at the lines from the `resol` column (or the resolution of the component, if set, as in the fitter; the larger of the two lines), k = 1 (Guido: the FWHM is the useful check; 4 b and 30 km/s were heuristics, kept as floors and as the fallback without resolution). The spectrum is passed to the grouping by every fit (recipes, `optimize_hierarchy`, Bayesian fit) and by the System Inspector. Only R < ~10000 (FWHM > 30 km/s) changes anything. With k = 1 the example above (60 km/s at R 8900, FWHM 34 km/s) is still not grouped: it would need k ~ 2; Guido kept k = 1 (#618). Mock benchmark: identical results on all 135 mocks (only the timings differ; the benchmark fits only the sparse metal doublets, and on two mocks at R 8900 no pair changed grouping), so a fit-level test replaces it (`test_fit_includes_component_blended_at_low_resolution`). Measured again on the merged code (master after #617 and #618, 3c4089c, against #617 alone, a7b02ce; both before F61, #619, not measured again after it): identical records on all 135 mocks |
| F16 | smoothing | For the same "std" parameter, V1 `_gauss_convolve` uses the kernel exp(-(dv/std)^2), i.e. sigma = std/sqrt(2), while V2 `smooth_spectrum` uses sigma = std: V2 smooths sqrt(2) wider. Affects the continuum smoothing (default 500 km/s: effective sigma 354 km/s in V1, 500 km/s in V2) and any other use of the parameter. Which convention to keep is Guido's decision | numerical | `test_smoothing_width_differs_by_sqrt2` (characterization); `test_estimate_auto_smoothing_is_a_fwhm` | fixed (phase 1-bis, decided with Guido): the continuum smoothing is given as a FWHM (`smooth_fwhm`, default 800 km/s, within 0.5% of the V1 continuum on J0103-1305; FWHM 500 would change it by up to 10% in the forest); `smooth_std` (sigma) kept as a deprecated alias |
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
| F29 | `VoigtFitterV2` (uncertainties) | The errors are sqrt(diag((J^T J)^-1)) with residuals weighted by 1/dy, not scaled by the reduced chi2, although the comment in `fit` says they are. V1 scaled them (lmfit `Model.fit`, `scale_covar=True` by default): a V1/V2 difference. Monte Carlo on mocks: unscaled errors are right when the error spectrum is right, and wrong by the same factor when it is not (e.g. 2 times too small with the error spectrum underestimated by 2, reduced chi2 3.8); scaled errors are right in all cases. Scaling, or reporting the reduced chi2 next to the errors, is Guido's decision | numerical, V1/V2 difference | `tests/test_error_scaling.py` (regression) | fixed (phase 1-bis, decided with Guido): errors scaled by sqrt(reduced chi2), as V1; the reduced chi2 is stored in each fitted component. On real data the scaling uses the chi2 of an intermediate pipeline state and of unrelated pixels: [#534](https://github.com/DAS-OATs/astrocook/issues/534) |
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
| F46 | `edit.set_properties` | The values are not checked: z_em = -5 and R = -1 are accepted without a warning (also from the GUI dialog). An agent passing a wrong value gets no signal. Found in the MCP pilot (PR 2) | behaviour, error handling | `test_mcp_server.py` (`test_unphysical_properties_are_refused`), `test_recipe_inputs_edit.py` | fixed in phase 6-bis, PR 2b (bounds decided by Guido): z_em, the z_ref of *Convert X axis* and the z_target of *Trim* above -1 and finite; a resolution R or N pixels positive and finite (0, negative and NaN were "no change" silently, '0px' an infinite R); the value is refused with an error that names it |
| F47 | `set_properties` (GUI) | When R is changed, the GUI (`_process_recipe_result` in `gui/main_window.py`) writes the new value into every component; the recipe does not, so a script or the MCP server setting R on a session with components leaves the old per-component resolution, which the fitter uses (F8). Same kind of logic as the session log, kept in the GUI instead of the recipe. Found in the MCP pilot (PR 2) | behaviour | — | fixed in phase 4 (#526, PR 1): neither the GUI nor the recipe writes the resolution into the components; `set_properties` keeps a column that varies with wavelength unless `overwrite_resol_col` (default now False, Guido's decision (d)); with `'_current_'` it no longer flattens it |
| F48 | `identify_lines` (docs) | The schema and the docstring of the recipe describe `mask_col` as "True=Unabsorbed"; the `abs_mask` column made by `find_absorbed` is True for the absorbed pixels (its own description, and how `identify_lines` uses it). The schema texts are what the GUI dialogs show and what the tools of the full MCP server will show to a model. Found in the MCP pilot (PR 3) | documentation | `test_every_recipe_and_parameter_is_described` | fixed in phase 6-bis (PR 6a): with all the descriptions of the registry reviewed; also `find_absorbed` said it writes a column `mask_unabs` (it writes `abs_mask`, True for the absorbed pixels) |
| F49 | `identify_lines` (candidate systems) | A doublet gives one candidate system per line: the two redshifts differ by ~0.5-1 km/s and survive the de-duplication, which rounds z to 5 decimals. On the full-chain mock the 3 CIV systems (4 components) give 6 candidates; on J0103-1305, 19 pairs of CIV candidates and 2 of MgII are closer than 2 km/s. `populate_from_identification` would add a component for each. Also, the candidate systems keep neither their region nor their score. Found in the MCP pilot (PR 3) | behaviour, structured output | `test_identify_lines.py`, `test_mcp_server.py` (`test_identify_lines_one_candidate_per_doublet`, `test_doublets_auto_fits_a_blend_with_one_component`) | fixed in phase 4: both regions of a doublet keep their label, only the region of the stronger line gives the candidate systems. Measured on the mock benchmark (2.1.0 → fix): duplicates CIV/SiIV/MgII 505/267/345 → 165/69/95 (those left are distinct components of a true system), purity 0.986/0.967/0.980 → 0.993/0.989/0.988, completeness unchanged; fit: spurious components 428 → 11, flagged 413 → 38, but true components matched 1534 → 1157, pull widths z/logN/b 1.8/2.3/2.1 → 4.3/5.5/4.5, system column spread 0.100 → 0.093 dex with pull width 1.24 → 2.53. The duplicate candidate, within a resolution element of the first, worked as a two-component start: the fit moved it onto a second true component when there was one, and made a spurious one otherwise. A variant keeping the candidates of the second line farther than c/R from those of the first (1171 matched, 17 spurious) was not adopted. Decomposition: F58; the region and the score of the candidate systems, fixed in phase 6-bis (PR 4b) |
| F50 | `components_auto` (`doublets_auto`, `lya_auto`, ...) | The pipeline always makes a new absorbed mask with the default parameters of `find_absorbed` (5000 and 400 km/s, kappa 2) and stores it in `abs_mask`, so the mask of an earlier "Auto-estimate continuum" with other parameters is replaced without notice (the continuum is kept); the GUI dialog and the tool expose no mask parameter. Also, the region left out by `region_limit` is cleared in the stored `abs_mask`: after `doublets_auto` on the full-chain mock no forest pixel is flagged as absorbed (36% before), so a later "Fit continuum" on that mask would treat the whole forest as unabsorbed (the automatic pipelines make a new mask each time and are not affected). Found in the MCP pilot (PR 4) | behaviour | `test_mcp_server.py` (`test_doublets_auto_keeps_the_mask_of_the_continuum`, `test_doublets_auto_without_a_mask_stores_the_whole_new_one`) | fixed in phase 4 (Guido's choice): the pipeline still makes a fresh mask with the default parameters (to avoid stale masks between recipes, as its comment says) and restricts it to `region_limit`, but uses it only inside; the session it returns keeps the mask of the user, or, if there was none, the whole new mask. Results unchanged. A mask parameter for the pipelines in phase 6-bis (PR 6b, Guido's choice): `mask` of `components_auto`, `doublets_auto`, `lya_auto`, `'new'` (the default, as before) or `'session'` (the `abs_mask` of the session; a `MissingPrerequisite` if there is none) |
| F55 | `identify_lines` (z range) | The doublet candidates are not limited to z <= z_em (plus a proximity margin): `z_em` is passed to `_find_kinematic_doublet_candidates` but not used. In Guido's MCP trial `doublets_auto` fitted a SiIV at z = 4.33 with z_em ~ 4.06 (16000 km/s beyond, in a telluric band). The schema of `doublets_auto` says "red side (z < z_em)". Found in the MCP pilot (Guido's trial) | behaviour | `test_identify_lines.py` (`test_candidates_beyond_z_em`) | fixed in phase 4: no candidate system beyond z_em by more than 3000 km/s (`Z_EM_MARGIN_KMS`, in `core/spectrum_lines.py` since the split of phase 5, a constant with its purpose in a comment; a parameter since phase 6-bis, PR 6b: `z_em_margin_kms` of `identify_lines` and of the pipelines, default 3000); associated absorbers and the uncertainty of z_em (~1000 km/s) stay within it. Not measured by the mock benchmark, which has no absorption beyond z_em |
| F60 | continuum (low S/N) | `find_absorbed_regions` excluded from the continuum anchors each pixel with y < 3 dy, meant to exclude the regions with no signal: at S/N ~3 this drops the lower half of the noise, and the continuum came out 26% too high on the mocks (the mean of a half Gaussian, +0.80 sigma: 26.6% expected, 26.5% measured). Also measured: the one-sided clipping at -2 sigma biases the continuum up by 0.055 sigma (truncated Gaussian), 0.6% at S/N 10, as found on the red side; the V1 fudge probably corrected it (F3) | numerical | `test_continuum_robustness.py` (`test_low_snr_red_side`) | fixed in phase 4 (#465): regions with no signal are those where the median S/N over 101 pixels is below 1; S/N 3: 1.265 → 1.02. The 0.055 sigma bias: phase 7 |
| F61 | Voigt model | The damping parameter of the Voigt profile is computed with the observed wavelength of the line centre, a = Gamma lambda_c / (4 pi b) with lambda_c = lambda_0 (1 + z), instead of the rest one, a = Gamma lambda_0 / (4 pi b) (V1 `_voigt_par_convert`, the mocks): the Lorentzian part is (1 + z) times too strong, 4 times at z = 3. The peak optical depth comes out lower (2% on SiIV with b = 1.5 km/s, logN 13.8 at z = 3.1) and the damping wings too strong: a column measured from the wings (DLA, LLS) comes out low by ~log(1 + z), 0.6 dex at z = 3. In the fitter (`_voigt_optical_depth`, `_voigt_optical_depth_vectorized`) and in the profiles of the System Inspector; there since the first V2 version. Found in #526 (PR 3): the residual of the fitter model against the exact flux of the mocks, 1.3e-3 to 1.8e-3 on CIV and Ly-a and 3.6% on a saturated SiIV line, drops to 3e-5 with the rest wavelength | numerical, wrong results | `test_damping.py` (fit of a damped Ly-a on mocks: logN 20.3 recovered within 0.01 dex, 19.75 at z = 2.5 and 19.64 at z = 3.5 before), `test_parity_v1_v2.py` (damped Ly-a, equivalent width 1.6 times that of V1 before), `test_lsf.py` (tolerance 1e-4) | fixed in phase 4 (separate PR after #526 PR 3, agreed with Guido): a with the rest wavelength, in the fitter and in the System Inspector. Benchmark on mocks (no damped systems there): spread of the system column 0.102 → 0.095 dex at R 8900, 0.046 → 0.041 at R 45000, the rest unchanged |
| F62 | automatic fits (forest) | The pixels of a fit are chosen at the start, around the lines of the active components, from their initial logN (`prepare_fit_context`, `_find_feature_limits`): if logN grows during the fit, the damping wings that fall outside cost nothing, and logN can run away to the bound (23). The groups fitted later see the component as a fixed background. On J0103-1305 `lya_auto` gives components with logN 21-22.4 (b 1 or ~50 km/s) that black out the model where the data are at 0.3-0.8 (540-590 nm) and lower it by 2-13% on the red side. Already present before F61 (logN 20.3-21 with b = 1 km/s, model black at 540-550 nm), made worse by F61: the damping wings, (1 + z) too strong, penalised the growth of logN inside the window. The mocks of the benchmark have no saturated forest lines strong enough to show it | numerical, wrong results | `test_snapshots_full_chain.py` (J0103-1305, characterization) | fixed in phase 4 (after F61): upper bound on logN in the automatic pipelines, `max_logN` = 19 (Guido's decision; the Lyman-limit systems stay within reach, damped systems are fitted by hand; `refit_all` and the manual fits keep 23). Tried first and dropped: extending the fit to the pixels reached by the wings at the fitted logN (within 1000 km/s). In the pipelines the absorption around a candidate is not modelled yet, so a wider window rewards the wings: more runaways (up to logN 22.97) and a slower chain. Tried also the bound at 17.5. On the J0103-1305 chain, forest 540-600 nm: chi2 per pixel 98 before F61, 305 after, 31 with the bound at 19 (28 at 17.5); over-absorbed pixels (model below the data by more than 5 sigma) 28%, 68%, 3.4% (2.3%); red side as before F61. The components at the bound are what a richer pipeline must handle (phase 7) |
| F63 | GUI (Set session properties) | The dialog fills each parameter whose default is `_current_` with the attribute of the same name of the spectrum data (`recipe_dialog.py`), and with "0.0" when there is none: the session name and the object have no such attribute (and the resolution no longer has one after #526), so the dialog shows "0.0" in them. Pressing Run without editing renames the session "0.0" and writes "0.0" into the OBJECT keyword of the header; the resolution is kept (a value of 0 is ignored), z_em is shown correctly. Found while scripting the screenshots of the tutorials (October 2026) | data loss (name, object) | — | fixed in phase 4 (Guido's agreement): the dialog shows the current value of the session (`RecipeDialog._current_value`: name, OBJECT, z_em, the resolution, left empty with the median as a hint when it varies), and an empty field keeps it (`'_current_'`). The same "0.0" was in z_ref of *Convert x axis* (velocities relative to z = 0 instead of z_em; the note of #573 that the dialog used z_em was wrong) and z_target of *Trim to common coverage*. `test_recipe_dialog_current.py` |
| F64 | GUI (plot controls) | The stylesheet of the checkboxes of the plot controls set `background-color: transparent`, which with the non-native Qt styles (Fusion: Linux, Windows) also clears the frame of the indicator: unchecked boxes were invisible, checked ones showed only the tick. macOS (native style) was not affected. Found while scripting the screenshots | GUI | `test_unchecked_boxes_of_the_plot_controls_are_visible` (regression) | fixed: the background colour is no longer set |
| F65 | `forest_auto` | The recipe takes `min_b` (schema and docstring: "Minimum b parameter (km/s)") but never uses it: no filter of narrow lines, unlike `lya_auto`. Found in #517 (PR 2) | behaviour | — | fixed in phase 4 (#517, PR 3, Guido's decision): the filter of `lya_auto` (`_remove_narrow`) applied to the series `Lyman`; `test_pipeline_cleanup.py` |
| F66 | main plot (velocity) | The ticks of the absorption systems are converted to the x unit with `getattr(spec, '_zem', 0.0)`, but `SpectrumV2` has `_z_em`, not `_zem`: the reference redshift of the ticks is always 0. After *Convert X axis* to km/s with z_ref ≠ 0 (the default is z_em) the ticks fall outside the plot and are not drawn; and a correct fix needs the z_ref of the conversion, which the spectrum does not keep (only its history line). Found in phase 5 while writing the snapshot of the ticks | GUI | snapshot `main_plot_system_ticks_kms`, `test_system_ticks_in_velocity_with_z_em` | fixed at the end of phase 5 (option chosen by Guido): the ticks are referred to z_em, as the conversion of the axis by default and all the other conversions (rebin, smooth, continuum, combination); the spectrum does not keep the z_ref of the conversion, F84 |
| F67 | `update_continuum_from_knots` | With `renorm_model=True` (the default, also of the recipe *Update from knots*, used when a continuum edited by hand is saved) the method renormalizes the model with the old continuum without checking that there is one: on a spectrum with a `model` column and no `cont` it raises `AttributeError` ('NoneType' object has no attribute 'value'). Rare (a model without a continuum, e.g. after removing the column). Found in phase 5 while writing the characterization tests of the continuum methods | crash | `test_knots_with_a_model_and_no_continuum_keep_the_model` | fixed in phase 5-bis (PR 1): without an old continuum the model is kept as it is, as in `fit_continuum` |
| F69 | flux methods (columns) | The methods that change the flux pick the "flux" columns by three different lists: `scale_flux` scales `cont` and `model`; `smooth` smooths `cont`, `model` and `telluric_model`; `calibrate` scales `cont`, `model`, `cont_pl` and `telluric_model`, and gives them all the unit of the flux. The telluric model is a dimensionless transmission (`SpectrumV2.norm` multiplies the continuum by it): after a calibration it carries the calibration curve and a flux unit, so the normalization is wrong; `cont_pl` is calibrated but not scaled. `y_convert` instead converts the auxiliary columns by their unit. Found in phase 5 while writing the characterization tests of the flux methods | behaviour | `test_spectrum_flux.py` (`test_scale_flux_by_a_number`, `test_smooth`, `test_calibrate_leaves_the_transmission`, `test_flux_columns_with_a_flux_without_unit`) | fixed in phase 5-bis (PR 2, agreed with Guido): two lists by name in `spectrum_flux.py`, the flux columns (`cont`, `model`, `cont_pl`: scaled, calibrated and smoothed) and the transmissions (`telluric_model`: only smoothed). Not by unit, as first proposed: many loaders give the flux no unit, and then the unit does not tell a continuum from a transmission or a resolution. Changes: `calibrate` no longer touches the transmission, `scale_flux` and `smooth` also change `cont_pl`. The loader of molecfit models (`ascii_telluric_model`) stores the transmission also as `model`, which is now scaled as a flux: left as it is (a telluric model, not a spectrum to calibrate; Guido's decision) |
| F70 | `generate_calibration_curve` | With two or more bands the scale factors are interpolated between the effective wavelengths of the filters, so the calibrated spectrum does not reach the target magnitudes: g = 17.5 and r = 17.0 give 17.40 and 17.08 on a flat spectrum (~0.1 mag); with one band the result is exact. Found in phase 5 while fixing the calibration (`np.trapz`) | numerical | `test_calibration_with_more_bands_is_exact`, `test_calibration_that_cannot_be_exact` | fixed in phase 5-bis (PR 3): the curve keeps its shape (linear between the effective wavelengths, constant beyond), and its values there are solved together instead of iterated: the mean f_nu in each band is linear in them, so one linear system gives the curve that reaches every target (to 1e-9 mag on flat and power-law spectra, two and three bands). When the solution is not positive (e.g. strongly overlapping filters with very different targets, or the same filter twice) a warning and the old factors of the bands alone |
| F71 | `get_resampled_column` | The error for a column that does not exist builds its message with `self.name`, which `SpectrumV2` does not have: the caller gets `AttributeError: 'SpectrumV2' object has no attribute 'name'` instead of the `ValueError` naming the column. Found in phase 5 while writing the characterization tests of the methods that combine spectra | error handling | `test_resampled_column_errors` | fixed in phase 5-bis (PR 1): the `ValueError` names the column |
| F73 | main plot (region selection) | The selection of a region with the mouse is dead code: `AstrocookToolbar.toggle_select_mode` (→ `SpectrumPlotWidget.toggle_region_selector`, then the span drawn by the canvas and `on_release`) is called by nothing, the toolbar has no 'Select' action (it is not in its `toolitems`), and the last step calls `main_window.launch_split_from_region`, which does not exist. Splitting by the current view is in the context menu (`launch_split_from_current_view`). Found in phase 5 while writing the characterization tests of the plot parts | code quality | `test_select_tool`, `test_region_menu`, `test_a_click_is_not_a_region`, `test_select_needs_a_wavelength_axis` | fixed in phase 5-bis (PR 8, agreed with Guido): restored, with a 'Select' button in the toolbar of the main plot after Zoom, exclusive with Pan and Zoom; a drag selects a region, which stays while a menu offers Split region, Mask region..., Rebin region... (the recipe dialogs filled in) and Zoom to region; disabled on a velocity axis. Found with it: the span did not follow the mouse (the start of the drag was stored on the widget, read on the canvas), and the old switch-off of Pan/Zoom had the bug of F81. Later (Guido): the statistics of a region (mean flux, rms, S/N, pixels); adding the region to the mask of the continuum (`abs_mask`), which needs a new recipe, in a later phase |
| F74 | continuum editor (save) | When the edited continuum is saved, `MainWindowV2._toggle_continuum_editor` refreshes the GUI by re-emitting `clicked` on the current index of the session list; when that index is not valid it calls `self.update_gui_from_state(new_session)`, which does not exist: the new state is already in the history, but the `AttributeError` is caught and logged as "Failed to save continuum", and the plot is not refreshed. Probably not reachable (an open session has a current index). Found in phase 5 while writing the characterization tests of the continuum editor | code quality | `test_save_with_no_current_session_in_the_list` | fixed in phase 5 with the split of `main_window.py`: the view is refreshed with `_update_view_for_session` (`update_gui_session_state` would have added the state twice) |
| F75 | main plot (SNR view) | The default error column of the SNR view is shown as "dF" and mapped to "dy" for the label and the check of the cache, but the arrays are computed with the text of the box: `spec.get_column('dF')` returns None, so each plot in the SNR view logs an error ("Could not get SNR error column 'dF'") before falling back to dy (the right values), and the caches of the plot are emptied at each plot ("dF" stored, "dy" compared). Found in phase 5 while checking the split of `plot_spectrum` | code quality, performance | `test_signal_to_noise`, `test_other_selected_sessions` | fixed in phase 5-bis (PR 5): the arrays use the mapped name of the view (`_PlotView.snr_col`), so no error in the log and the SNR arrays stay cached. Found with it: the other selected sessions ("Ghost") were left in flux in the SNR view with the default "dF" (the lookup failed silently); now divided by their error |
| F76 | main plot (decimation) | Spectra longer than 2 million pixels (`DECIMATION_THRESHOLD`; the first comment said 20k) are decimated (min-max, 2 points every 10) only in the full view of the first plot after the axes are cleared: Home (`force_autoscale`) excludes it, and at any later plot the limits are kept (`was_zoomed`). The decimated arrays stay in the cache and are then drawn as steps with the label "Spectrum", and an auxiliary column (not decimated, other length) fails with a warning, until the cache is emptied. Found in phase 5 while checking the split of `plot_spectrum` | code quality, performance | `test_decimation_of_long_spectra` | fixed in phase 5-bis (PR 9, agreed with Guido): the arrays in view are decimated at every plot when they are more than `DECIMATION_THRESHOLD`, after the cut to the view (full view, Home, a wide zoom; not a narrow zoom); the cache keeps the whole arrays, the decimated ones are cached by view; the auxiliary column cut and decimated as the flux. Measured (offscreen, synthetic spectrum, legend in a fixed corner): a redraw of 3 million pixels 3.9 s whole, 1.0 s decimated; of 300 000 pixels 0.77 and 0.57 s. Threshold lowered from 2 to 1 million (Guido). The measure found F86 |
| F77 | main window (placeholders) | `_launch_x_convert_dialog`, `_launch_y_convert_dialog` and `_launch_rebin_dialog` only print "Launching ... Dialog..." on the standard output: nothing calls them (the conversions of the axes are in the plot controls, rebin is a recipe of the Flux menu). Found in phase 5 while writing the characterization tests of the menus | code quality | `test_placeholder_dialogs_are_gone` | fixed in phase 5-bis (PR 7): removed (Guido's decision) |
| F78 | plot controls (log F) | The "Logarithmic F" box (`log_y_checkbox`) is created and connected (the plot reads it), but it is never added to the View Controls section (only "Logarithmic λ" is): the log scale of the flux cannot be chosen from the GUI. Found in phase 5 while writing the characterization tests of the panels | GUI | `test_log_f_is_in_the_panel`, `test_logarithmic_flux`, `test_custom_limits_with_log_f` | fixed in phase 5-bis (PR 6, agreed with Guido): the box below "Logarithmic λ", with a tooltip; pixels with F ≤ 0 clipped to the bottom (as Matplotlib does); in the normalized view the default limits are (0.01, 1.3) instead of (−0.3, 1.3), which a log axis cannot show (the axis was inverted, with a warning); a zoom with a lower limit ≤ 0 is not kept (normalized: from 0.01, otherwise autoscale); a non-positive F Min typed in *Axis Limits* is not applied, with a message. The iso-mag grid works with it (straight lines on a calibrated spectrum) |
| F79 | sessions (Equalize and Stitch) | `_on_equalized_stitch_dialog_requested` writes the number of values expected ("requires N values") into the documentation of three parameters of `equalize_and_stitch` in `EDIT_RECIPES_SCHEMAS`, the schema shared by the whole application: the text stays until the next call, and the dialog opened from the Edit menu (or anything reading the schema in the same process) shows the count of the last stitch. Found in phase 5 while writing the characterization tests of the sessions | code quality | `test_equalize_and_stitch`, `test_gui_dialog_asks_for_one_cut_with_two_arms` | fixed in phase 6-bis, PR 1: the counts go to the dialog alone (`param_docs` of `RecipeDialog`), the registry is read only |
| F80 | main window (dead code) | `_show_refit_warning_with_checkbox` (replaced by `_show_custom_message` with `checkbox_text`) and `_SCRIPT_LINE_REGEX` are not used. Found in phase 5 while writing the characterization tests of the recipes | code quality | `test_dead_code_is_gone` | fixed in phase 5-bis (PR 7): removed (Guido's decision), with the imports used only by them |
| F81 | recipe dialog, script run (toolbar) | Before opening a recipe dialog (and before running a script) the code calls `toolbar.pan(False)` and `toolbar.zoom(False)` to switch the tools off, but matplotlib's `pan` and `zoom` toggle the tool and ignore the argument: with no tool on, pan is switched on and then replaced by zoom, so after any recipe dialog the main plot is in zoom mode. Found in phase 5 while writing the characterization tests of the recipes | GUI | `test_recipe_dialog` | fixed in phase 5-bis (PR 5): `AstrocookToolbar.deactivate_tools` switches off only the active tool, before a recipe dialog and before a script run |
| F82 | export of ASCII files | The export decides whether a spectrum file is written (for the suggested name and the check of existing files) from the internal names of the columns (x, y, dy, cont, ...), while the export itself also takes the names shown in the GUI (λ, F, dF): with targets named only as in the GUI (e.g. "F, dF", or the default "F, cont, systems" on a session without continuum) the spectrum file is not expected, the suggested name is the bare prefix and an existing `<prefix>_spec.csv` is overwritten without asking. **Loses user data** (a file on disk). Found in phase 5 while writing the characterization tests of the recipes | data loss | `test_export_with_the_names_of_the_gui`, `test_export_ascii_outputs` | fixed in phase 5 (agreed with Guido): the recipe tells which files it writes (`RecipeFileV2.export_ascii_outputs`, the same rule as `export_ascii`) and the GUI asks it; the target `spec` (all the columns of the spectrum, as documented) wrote nothing and now writes the spectrum file. Not #630 (errors of linked parameters in the exported systems, phase 5-bis) |
| F83 | export of ASCII files | When the files to export already exist, the question before overwriting raises `UnboundLocalError`: `_on_recipe_requested` imports `QMessageBox` again further down (`from PySide6.QtWidgets import QMessageBox, QApplication`, in the checks of `set_properties` and `coadd`), so the name is local to the whole method and not yet bound in the export branch. Nothing is written and no message is shown (the error goes to the console). Found in phase 5 while writing the characterization tests of the recipes | crash | `test_export_over_existing_files` | fixed in phase 5 (agreed with Guido): the two local imports removed (the names are imported by the module) |
| F85 | System Inspector (boxes above the table) | The row of the boxes above the table ("Group View (Selected Component Only)", "Only flagged", "Select flagged") had no spacing of its own: a nested layout takes that of its parent, 0 here, so the boxes touched, and overlapped on macOS (reported by Guido). The label said "Selected Component Only", the tutorial "the components physically linked to your selection": the box shows the group of the selection (`get_connected_group`: linked, blended or close in velocity), the components fitted together. Reported by Guido in phase 5-bis | GUI | `test_inspector_boxes_above_the_table` | fixed in phase 5-bis (PR 6-bis, agreed with Guido): "Only group" with a tooltip, "Inspect ... (its group)" in the context menu of the main plot, spacing 12 px |
| F86 | main plot (legend) | When both side panels are closed, and at the first plot, the legend of the main plot is placed with `loc='best'` (`_legend_location`): Matplotlib then counts the points of every line at each candidate position, which dominates the redraw of a long spectrum: 300 000 pixels (an ordinary ESPRESSO spectrum, below the threshold of the decimation) 6.3 s instead of 0.77 s with a fixed corner; 3 million pixels 60 s whole, 8.4 s decimated. Found in phase 5-bis while measuring F76 | performance | `test_no_legend`, `test_no_legend_on_the_plot`, `test_swatch_of_the_auxiliary_column`, `test_swatch_of_a_shaded_area`, `test_other_selected_sessions` | fixed in phase 5-bis (PR 9-bis, Guido's decision): no legend on the plot; the swatches of the panels replace it: the Plot Elements section as before, a swatch next to **Aux. Column**, and one next to the name of each other session shown in the list of the sessions (`SessionListModel`); the redraws that only moved the legend (the card of the plot controls opened or closed, its section toggled) removed. Measured as for F76: a redraw of 300 000 pixels 0.65 s (6.3 s with the legend at 'best'), of 3 million pixels decimated 0.83 s (8.4 s) |
| F87 | dialog *Session info* (stylesheet) | The stylesheet of the dialog (`MainSessionsMixin._on_session_info`) is a plain string, not an f-string, but is written as one: it contains `{palette.color(...) if 'palette' in locals() else '#FFFFFF'}`, `{text_color}` and doubled braces, so it is not valid CSS and its fields do not get the style of the other dialogs (rounded border, accent when focused). Found in phase 5-bis while gathering the colours into the theme module | GUI | `test_session_info_stylesheet` | fixed in phase 5-bis (theme PR B): an f-string with the palette and `FOCUS_BORDER` |
| F88 | plot controls (drop-down menus) | The three rows with a drop-down menu (Aux. Column in Plot Elements; λ Unit and Error Col. in View Controls) were three form layouts, each with the width of its own label (and "Error Col." indented with two spaces), so the menus started at three different places and had different widths. Reported by Guido in phase 5-bis | GUI | `test_the_drop_down_menus_line_up` | fixed in phase 5-bis (PR of View > Appearance): one kind of row (`_combo_row`), the labels as wide as the widest, measured with the font of the stylesheet; the swatch of the auxiliary column moved into its menu (the icon of the column chosen) |
| F89 | dialog *About* (tagline) | The tagline "A thousand recipes to cook a spectrum" is in 16pt bold in a dialog 400 px wide: wider than the dialog (offscreen with Fusion 400 px against 360 available), it is clipped on both sides ("thousand recipes to cook a spectru"). Found in phase 5-bis by the visual review | GUI | `test_about_tagline_is_not_clipped` | fixed in phase 5-bis (PR of the review bugs): the tagline wraps when it does not fit |
| F90 | `components_auto` (identifications of the regions) | The identifications of the absorbed regions (`region_identifications` in the metadata, {region ID: candidates}, written by `identify_lines`) are keyed by the IDs of the regions in `abs_ids`. When a pipeline runs on a session that already has regions (e.g. `lya_auto` after `doublets_auto`), `components_auto` shifts the IDs of its new regions after the old ones in `abs_ids`, but replaced the metadata with those of its own run, keyed by the unshifted IDs: the identifications of the old regions (the metals of `doublets_auto`) were lost, and those of the new regions were attached to the wrong IDs. Seen in the *Identifications* window ("222 regions" but only Ly_a at 17 redshifts, some beyond z_em) and in the tooltips of the absorbed regions in the main plot (`plot_mouse`), which showed the identification of another region. The fits are not affected (the systems are built from the identifications of their own run, keyed by system), nor the identification summary of the MCP server (one run). Found in phase 5-bis by the visual review | behaviour | `test_merge_region_ids_shifts_the_identifications`, `test_merge_region_ids_without_old_regions`, the intent of `test_full_chain_mock` (one identification per region, those of the doublets kept after `lya_auto`) | fixed in phase 5-bis (PR of the review bugs): `merge_region_ids` shifts the keys of the new identifications as their regions, keeps the old ones, and drops those of regions no longer in the mask. Sessions saved before keep the wrong identifications until the pipelines are run again |
| F91 | main plot (redshift cursor) | The lines of the redshift cursor are blitted over the pixels of the plot kept after a full draw (`SpectrumPlotCanvas.background`), captured once and renewed only when the plot is rebuilt or its limits change: when the plot was resized (e.g. *View > Show Toolbar*, the window resized) the cursor restored the old pixels, with the old size, over the new plot. Reported by Guido in phase 5-bis | GUI | `test_cursor_after_a_resize` | fixed in phase 5-bis (PR of the review bugs): the pixels are captured again at the redraw after each resize |
| F92 | main plot (redshift cursor) | While the redshift cursor follows the mouse over the plot, the card of the plot controls fades (`CardFader`, opacity 0.35 while the mouse is away), and the redshift written in its box is hard to read. Reported by Guido in phase 5-bis | GUI | `test_cursor_label` | fixed in phase 5-bis (PR of the text pass; Guido chose the first of three options: the redshift on the plot, the card kept opaque, the redshift next to λ and F in the toolbar): the redshift written at the top of the first line of the cursor in view, blitted with the lines |
| F94 | `save_archive_v2` (saving a session) | The `.acs2` archive is written directly over the target file: `tarfile.open(final_path, 'w:gz')` truncates it before its members are added, so a save that fails while the archive is being written leaves a partial archive (or an empty file) in place of the previous one, and the previous session is lost. One way to get there: if writing the original data (`session_v2_initial`) fails, the broad `except` logs an ERROR and goes on with the name of a file that does not exist, and `tar.add` then fails; reproduced in phase 6 (PR 1) by making the conversion of the original spectrum raise: the save raises `FileNotFoundError`, and the target file, an existing archive, is replaced by one with the current spectrum, systems and metadata only. Data loss, with a low probability (a failure while writing). Found in phase 6 while making the inventory of the broad except clauses | data loss, error handling | — | `test_save_archive.py` | fixed in phase 6 (PR after the inventory, agreed with Guido): the archive is written to a temporary file in the same folder and replaces the target only when complete (`os.replace`; symbolic links followed, permissions of the file replaced kept), and a failure of the original data stops the save before the target is touched |
| F95 | fallbacks after a failure (24 places) | After a failure, the code goes on with a result that looks valid and says so only in the log, or not at all: e.g. a link of a constraint that fails to evaluate keeps its old value silently at every evaluation of the model (`VoigtModelConstraintV2.map_p_free_to_full`); `smooth_spectrum` returns the unsmoothed array; `update_continuum_from_knots` returns the spectrum unchanged, so the edit of the continuum is lost; `delete_component` returns the unchanged session, which the GUI takes as a success; the entry of a recipe in the session log is dropped silently if `add_entry` fails; in the pipelines, a cleanup or a fit that fails is skipped with a WARNING (`lya_auto`, `forest_auto`, `refit_all`) or a DEBUG line (`metals_auto`). Found in phase 6 (PR 1); in PR 4a, without a broad `except`, also `fit_continuum` (an ERROR and the spectrum returned unchanged when fewer than two pixels are unabsorbed) and `update_continuum_from_knots` (a WARNING and the spectrum unchanged without knots) | error handling, behaviour | — | fixed in phase 6: the rows D of the inventory, each decided with Guido and reproduced by a test, those of core, fitting and io (with `fit_continuum` and the knots) in PR 4b (#695), those of the recipes (with the four sessions returned unchanged and the regions of `refit_all`) in PR 5b (#697), those of the GUI (the entries of the session log) in PR 6b (#700) |
| F96 | logging (scripts) | The package logged with the module-level functions of `logging` (`logging.warning(...)`), which call `logging.basicConfig()` when the root logger has no handler: the first message of a recipe in a script that had not configured logging yet configured it (level WARNING, format `WARNING:root:...`), and the script's own `basicConfig(level=logging.INFO)` called afterwards did nothing, so its level and format were ignored. Found in phase 6 (PR 3) while writing the characterization tests of the loggers per module | behaviour | `test_script_can_configure_logging_after_a_recipe` (regression) | fixed in phase 6 (PR 3): the loggers per module (#566) do not configure the root logger; a script that configures nothing still sees the warnings (Python's last-resort handler), without the `WARNING:root:` prefix |
| F97 | `detect_file_format` (loaders) | The detection of an ASCII format reads the first 10 lines with `next(f)`, which raises `StopIteration` for a file with fewer lines: the exception was caught by a broad `except`, and a `.csv` file of fewer than 10 lines is read as a generic spectrum instead of a CSV table (by the generic loader, which may then fail on it). Found in phase 6 (PR 4a) while narrowing the clause | behaviour | `test_short_csv_is_read_as_csv` | fixed after phase 6 (Guido): at most 10 lines read (`itertools.islice`); the generic loader reads FITS only, so such a file did not open at all ("No SIMPLE card found"), not only as a generic spectrum; `test_short_csv_is_read_as_csv` (regression) |
| F98 | links (functions) | The functions offered to the expression of a link were not the same everywhere: the fit (`VoigtModelConstraintV2`) offered `np`, `log` (log10) and `ln`, the automatic update of the linked value in the `update_constraint` recipe `np`, `sqrt` and `log`. A link with `sqrt(...)` was accepted and set its value once, and was then skipped silently by every fit (the parameter stayed at that value); one with `ln(...)` was accepted with a WARNING. Since PR 4b the links are checked with the functions of the fit, so `sqrt` is refused; in PR 5a the update uses the same functions (`LINK_FUNCTIONS`). The System Inspector writes the thermal link of b with a number (the square root of the mass ratio), so only the links typed by hand were affected. Found in phase 6 (PR 5a) | behaviour | — | fixed in PR 5a (one set of functions); no functions added (Guido, November 2026: a link is the same value, or the same temperature, which the System Inspector writes as a number; the turbulent part has its own parameter, `btur`) |
| F99 | fits without a resolution, automatic continuum | Found in phase 6 (PR 5c) while analysing the prerequisites of the recipes (Guido's question): (1) a fit on a spectrum without a resolution uses unconvolved profiles, so b includes the instrumental broadening, and said so only at DEBUG (`VoigtFitterV2`); the GUI asks for the resolution first, a script or an agent got no notice; (2) when the automatic continuum, run by `populate_from_identification`, `fit_component`, `optimize_system` and `refit_all` on a spectrum without one, failed, its 0 went on as the session and the recipe failed with "'int' object has no attribute 'systs'"; (3) `refit_all` estimated the continuum without saying so | error handling, behaviour | `tests/test_error_handling_recipes.py` | fixed in PR 5c: a WARNING (`msg_no_resolution`) at the start of the fit recipes (`refit_all`, `fit_component`, `optimize_system`, `fit_bayesian`, `metals_auto`), once per call (the GUI and the MCP server drop the repetitions); `_with_continuum` stops the recipe with an ERROR when the estimate fails, and says that it estimates; the way of filling the prerequisites redesigned in 6-bis (phase 6-bis, "Prerequisites of the recipes") |
| F100 | `load_simple_csv` (loaders) | The loader of the `.csv` files reads them with `Table.read(format='ascii.csv')`, which takes the first line as the names of the columns, and then asks for the columns `col1`, `col2`, `col3`: a CSV file of numbers only (x, y, dy, no header) does not open (`KeyError: 'col1'`), whatever its length; only a file with the header `col1,col2,col3` (or one whose first line is taken as such) does. The units are fixed (nm, erg/cm2/s/Å). Found after phase 6 while fixing F97 | behaviour | `tests/test_loader_csv.py` | fixed after phase 6 (Guido): the first three columns, by position, with or without a header (any names, comment lines with `#`); fewer than three columns a `ValueError` that says so; the units stay fixed (*Getting Started*) |
| F103 | `identify_lines`, the pipelines, `metals_auto` (spectra without a continuum) | On a spectrum without a continuum (physical flux), `identify_lines` computes the optical depth with a continuum of 1 (`spectrum_lines.py`, `spec.cont` or ones), i.e. on the raw flux, without a word; the pipelines (`components_auto`, `doublets_auto`, ...) identify so, and the continuum is estimated only later, by `populate_from_identification`, for the fit; `metals_auto` tries its metals with fits on the raw flux before the `refit_all` that estimates it. On a mock in physical units `doublets_auto` found no component (the CIV with the continuum first). Found in phase 6-bis (PR 5a) while checking that the continuum provided first gives the same numbers | behaviour, numerical | `test_provided_first_as_the_replay_will` | fixed in phase 6-bis (PR 5a, Guido's decision): these recipes require the continuum, provided first (or a `MissingPrerequisite` with the policy `raise`, a question of the GUI in 5b). The benchmark of the mocks is not affected: the mocks are normalized (continuum 1) |
| F104 | GUI, the recipe resumed after a prerequisite | A recipe kept pending while Set Properties provided z_em or the resolution stayed pending when Set Properties failed, and ran at the next Set Properties that succeeded, perhaps much later, unexpected. Found in phase 6-bis (PR 5b) while generalizing the pending recipe | behaviour | `test_a_failed_estimate_drops_the_pending_recipe` | fixed in phase 6-bis (PR 5b, #717): the pending recipe is dropped when the recipe that provides what it lacked fails or its dialog is cancelled |
| F105 | recipe dialogs, the buttons of the sessions | In the dialogs of *Import systems*, *Equalize* and *Resample*, the buttons of the open sessions wrote the name into a field 'expression' that these dialogs do not have: a click did nothing, and the name had to be typed. Found in phase 6-bis (PR 6c) while showing that panel by the kind of the parameter | bug | `test_the_session_buttons_fill_the_session` | fixed in phase 6-bis (PR 6c): the name goes in the parameter of kind session |
| F106 | `export_ascii` (*File > Export to ASCII*) | The recipe takes a `prefix` (the name of the files) that the registry does not declare: the GUI gives it from its file dialog, so the recipe cannot be run from the registry alone (a script, the MCP server). It writes `<prefix>_spec.csv` and `<prefix>_systems.csv` with `overwrite=True`: only the GUI checks before for files that would be overwritten (F82), a script or an agent overwrites them without a word. Found in phase 6-ter (PR 1) while generating the tools from the registry | behaviour (data loss outside the GUI) | `test_every_recipe_is_a_tool_or_excluded_with_a_reason` | fixed in phase 6-ter (PR 1-bis, Guido's decision): `prefix` in the registry (hidden in the GUI, as `path`) and `overwrite` (False by default: the recipe stops, writes nothing and names the files that exist; the GUI asks as before and passes True); a tool of the MCP server again. Scripts that export twice to the same files must pass `overwrite=True`. Why 6-bis missed it: its test checked that every parameter of the registry is in the signature of the method, not the converse; the converse is now tested (`test_every_parameter_of_the_method_is_in_the_registry`, the parameters left out listed with their reason) |
| F108 | `session_open` (MCP), `set_properties` | `session_open` with `z_em` or `resol` runs `set_properties`, which also writes the OBJECT of the spectrum as `''` when it had none: the result lists "object None → ''", a change that is none. Found in phase 6-ter (PR 5c) on the J0103-1305 session of the draft | cosmetic, metadata | `test_set_properties_current_object_without_one`, `test_set_properties_current_object_keeps_the_header` | fixed (phase 6-ter, after PR 5c): `'_current_'` leaves the object and the header as they are; it also copied `OBJECT` over the ESO target name and added "Object='...'" to the history of the spectrum (`test_set_properties_current_object_*`) |
| F109 | `.acs2` original data (`save_archive_v2`, `load_session_from_file`) | An archive keeps the original data (`*_spec_orig.fits`, `*_systs_orig.fits`, the first state of the history) so that its log can be replayed from them, but they are not read back: an archive opened becomes a history whose first state is its current state, with the log of the archive. Saved again (GUI *Save Session*, the MCP `session_save` of a session opened from an archive, the editable copy of a followed session), the archive gets that state as its original data, and its log still lists the recipes that made it: a replay from the original data runs them a second time, on data already processed (reproduced in phase 6-ter PR 5-bis: `set_properties z_em=4.072` saved, reopened and saved again; the original data of the second archive have z_em = 4.072). *Run All* of the Log Scripter, which replays from the first state of the history, does the same on an archive just opened. Known since phase 1 ("stored but not loaded back", `ROADMAP_HISTORY.md`, step 5 of phase 1), never recorded as a finding with this consequence. Found while writing the follow of the GUI (phase 6-ter PR 5-bis) | data loss (the original data of an archive saved again), reproducibility of the log | `test_session_nodes.py` | fixed (phase 6-ter, after PR 5-bis, #740): reported to Guido at once (PR 5-bis); decided with Guido (November 2026), to fix in its own PR after PR 5-bis: the rule of the sessions is restored, each saved archive is a node: A (the `_orig` copy, completed with the metadata of A: constraints, identifications, in a `_meta_orig.json`), B and the log from A to B. When an archive is opened, B becomes the new A and the log starts empty (GUI, `session_open`, the editable copy); the log of the archive is shown read only while the session is open (Log Scripter, `session_describe`), never replayed nor saved again. The log of the earlier nodes is not kept in the archive (Guido: an analysis by an agent is done in one session, or the agent saves its own log of the operations to resume from; not needed in the archive). Considered and set aside: reading `_orig` back as a fixed origin with the log continued (the log grows without limit, and contradicts the rule); a reference to the previous node with its hash (depends on files being kept) |
| F110 | dependencies (`pyproject.toml`) | `scienceplots` is a required dependency of the package, but no code imports it: the only mention is a commented line in `gui/plot_canvas.py` (`#plt.style.use(['science', 'fast'])`; the code uses matplotlib's own `fast` style). Found in phase 6-ter PR 6: the headless image leaves it out (`tests/test_docker_image.py`, `GUI_ONLY`) and the server works. | Every installation pulls a package it never uses (and its LaTeX-related style files). Harmless; to remove with the review of the dependencies (F18). | Phase 6-ter PR 6 | Fixed (November 2026, after phase 6-ter PR 6, Guido's decision): removed from `pyproject.toml` and from `Astrocook.spec` (the macOS build collected it), and the commented line from `gui/plot_canvas.py`; a test checks that nothing imports it. The image's list and `pyproject.toml` stay aligned. |
| F17 | continuum (forest) | On the J0103-1305 window, with the same parameters and matched smoothing, V1 `clip_flux` and V2 `find_absorbed` + `fit_continuum` agree on the red side (masks 99.98%, continuum within 0.4%), but in the Lyman-alpha forest V2 marks more pixels as unabsorbed (57.6% vs 54.8% of the window, masks agree on 92.7% of the forest pixels) and its continuum is lower by 1.2% (median; 5th percentile -4.6%). Not explained by the fudge (F3). Relevant to #465 and to the known low forest continuum | numerical | `test_continuum_forest` (strict xfail) | partly addressed by #465 (measured on the mock benchmark (2.1.0 + F49/F50 → #465): forest continuum 0.936 → 0.943 (5th percentile 0.820 → 0.830), red side unchanged (1.002)); V1 parity no longer sought for the continuum (V2 is measured against the truth on mocks instead). Closed in PR 13 of 6-quater (Guido, 13 November 2026): accepted; the strict xfail of `test_continuum_forest` stays as the record of the difference with V1; the forest continuum is measured against the truth on the mocks (F56, phase 7) |
| F107 | GUI menus and the registry | The labels of the menus of the GUI are written in `gui/main_menus.py` and differ from the `brief` of the recipes in the registry, which the dialogs show as their title (*Auto-detect Doublets...* opens *Auto-fit metal doublets (Red Side)*; *Set Properties...*, *Set session properties*). The methods text (PR 5c) names a recipe by its menu and `brief`, so the name it gives may differ from the label of the menu. Found in phase 6-ter (PR 5c) | cosmetic | — | fixed in PR 11 of 6-quater (Guido): a `menu` key in the registry, the only place of the label, which the menus, the title of the dialogs and the methods text read; the recipes with no entry name no menu |
| F111 | `metals_auto` (and the fitter's fallback) | `VoigtFitterV2` fits the whole spectrum when no valid pixel is around the components (`voigt_fitter.py`, warning "No valid pixels around the components to fit: the whole spectrum is used"). `metals_auto` met it in both trials of PR 7: on the full-chain mock (31 thousand pixels) it kept no ion, missing the SiIV at z = 1.90 (logN 13.3, b 12, in the forest; the agent found it by eye and it fits with reduced chi2 1.1), and it does not say which ions it tried nor their AIC; on J0103-1305 (215 thousand pixels, 65 components, six ions) it ran more than 2 h 20 min without finishing, with no progress line. The cause of the missed SiIV is to be confirmed (why no valid pixel: the window of the trial ion, its mask, or the names of the transitions) | numerical, performance | — | fixed (PR 7-quater). Cause found on the mock of trial 1 (a script on the chain of the trial): the trial fit found the SiIV (logN 13.28, b 11.3, as the agent), and the recipe returned `refit_all` of the session it started from, so `metals_auto` never added an ion, whatever it found (the test of the recipe recorded "no ion kept" on a spectrum where all the trial ions were out of range). The fallback on the whole spectrum was not the cause: it was met by the ions whose lines fall outside the spectrum (7 of the 15 trials, NV, OVI and MgII at z = 1.76 and 1.83: no pixel at all), which it fitted on all the pixels (31 thousand; 215 thousand on J0103-1305, hence the hours) and kept at their first values (logN 13.00), so that fixing the result alone would have added them. The fallback that Guido recalls is another one, in `prepare_fit_context`: a window narrower than 30 km/s replaced by one of at least 150 km/s (now said in a comment). Fixed (Guido's choices): the session with the ions found is refitted; the fitter raises `NoPixelsToFit` (a `ValueError`, with the components and the range of the spectrum) instead of fitting the whole spectrum, so that `add_component` of an ion outside the spectrum fails with that message and `refit_all` lists such a group as not fitted; `metals_auto` skips those ions with a WARNING that lists them, reports its progress (one line per trial), checks the cancellation at each trial, and reports each ion tried with its AIC change and the reason when not kept (`tried` in `details`). Its acceptance criterion: F121 |
| F113 | MCP server (concurrency, long recipes) | One lock (`SessionStore.lock`) serializes every tool, the read-only ones too (`session_history`, `components_list`, `spectrum_plot`): while a recipe runs no other call is answered. With no heartbeat (PR 4 left it out) and no cancellation, a client that gives up after its timeout leaves the recipe running and the session out of reach: in trial 2 of PR 7 the agent could not save its work, which was lost. | design, data loss (the work of a session) | — | fixed: a heartbeat after 20 s of silence, and the read-only tools, `session_save` and `session_note` outside the lock of the recipes (PR 7-bis); the cancellation (PR 7-quater): the client's `notifications/cancelled`, or `session_cancel`, stops the recipe at its next check, with an error and no new state |
| F117 | `update_constraint` (description and checks) | The description of `expression` gives `'p["z"]'` as its example; the expressions that work name the target, `p['<uuid>'].z` (the docstrings of the method and of `SystemListV2` say so). With the example, the constraint is stored with a warning ("the linked z could not be computed") and the next fit ignores it: in trial 3 of PR 7 the SiIV kept its own z after a refit, and `components_list` showed `linked: p["z"]`. `target_uuid` is stored next to the expression but not used to evaluate it. An expression that cannot be evaluated is kept, not refused | doc, behaviour | — | open; decided (Guido, on PR 7, November 2026): the example corrected (`p['<uuid>'].z`) and an expression that cannot be evaluated refused (`InvalidParameter`, nothing stored), not stored with a warning; `target_uuid` not used to expand a bare name. fixed in PR 7-ter (the check runs also with `is_free` left out, where the link was stored unchecked and kept frozen by the fit with an ERROR line; the example corrected); closed in PR 10 of 6-quater: the refusal is an `InvalidParameter` (it was a `RecipeFailed` from the core), the last item of the decision |
| F120 | names of the columns (expression recipes, MCP reads) | `edit_apply_expression` and `edit_mask_expression` document the columns by the names of the GUI (`λ`, `F`, `dF`, `cont`), `file_export_ascii` too (`F, dF, cont`); `spectrum_read` and `spectrum_plot` want `x`, `y`, `dy`. In trial 3 of PR 7 a read with `F`, `dF` failed (the error listed the valid names) | cosmetic, output | — | open; decided (Guido, on PR 7): `spectrum_read` and `spectrum_plot` also accept the names of the GUI (`λ`, `F`, `dF`, ...) as aliases, fixed in PR 7-ter (`SpectrumV2.columns_in_range`, values under the name asked; `spectrum_plot` takes no columns); the one set of names in phase 6-quater. Fixed in phase 6-quater, PR 8 (Guido, 12 November 2026): the names of the GUI everywhere a column is named (`core/column_names.py`), `lam` (`lammin`, `lammax`) its ASCII spelling, the names of the data (`x`, ..., `dy`) silent and permanent aliases; two cases that failed found there: `calculate_running_std` did not find `F` (the button of its dialog), and wrote a column `F` that no expression could read |
| F122 | `calculate_running_std` (`output_col`) | With `output_col` empty or the name of a column of every spectrum (`x`, `xmin`, `xmax`, `y`, `dy`) the result went to an auxiliary column of that name, which the table, the dialogs and the reads do not show (the column of the spectrum unchanged): the result was lost without a word. Found in the review of PR 7 (phase 6-quater) | behaviour | `tests/test_recipe_review_file_flux.py` | fixed (Guido, PR 7 of 6-quater): `InvalidParameter` before the recipe runs |
| F123 | `equalize` and `coadd` (`order = -1`) | `-1` means two things: a spline in `flux.equalize` (`order`), no scaling in `edit.coadd` (`equalize_order`). Found in the review of PR 7 (phase 6-quater) | cosmetic, interface | `test_coadd_minus_one_is_none` | fixed in PR 9 of 6-quater (Guido, option A): `equalize_order` of `coadd` is `Optional[int]`, None no scaling (as `kappa` None no clipping), 0-5; -1, of the old logs, is read as None (silent and permanent), below -1 refused; `equalize` keeps -1 (a spline) as its lowest value (PR 7) |
| F124 | `rebin` (`xstart`, `xend`) | With `xstart` not below `xend` the recipe gave a spectrum without pixels, without a word. Found in the review of PR 7 (phase 6-quater) | behaviour | `tests/test_recipe_review_file_flux.py` | fixed (Guido, PR 7 of 6-quater): `InvalidParameter`; the same check for the range of `coadd` in PR 9 (done, `test_coadd_refused`) |
| F125 | GUI, the *Help* button of the dialogs | The `url` of the recipes in the registry are relative names of pages of the V1 documentation (`file_cb.html#export`, `edit_cb.html#rebin`, ...; some marked "Placeholder URL"), which `webbrowser.open` cannot open as pages of the documentation; some recipes have none. All the recipes. Found in the review of PR 7 (phase 6-quater) | cosmetic | — | fixed in PR 11 of 6-quater: every recipe has a `url` relative to the root of the documentation, a section of its tutorial, and a test checks that the section exists |
| F126 | `apply_expression` (units) | A new column made by an expression is dimensionless whatever the expression: `F_corr = F * 2` loses the unit of the flux (`erg/(nm s cm2)`); a column replaced keeps its unit, also when the expression changes it (`F = F / cont`). Found in the review of PR 8 (phase 6-quater) | behaviour | — | open (Guido, PR 8: registered, not changed; a unit from the expression would need the units of NumExpr's operations); documented in the tutorial of Edit; closed in PR 13 of 6-quater (Guido, 13 November 2026): accepted |
| F127 | expression recipes (the axis) | `apply_expression`, `mask_expression` and `smooth_column` wrote the columns of the axis: `λ = λ * 0 + 1` replaced the wavelengths (without `λmin` and `λmax`), a mask put NaN in them, the smoothing smoothed them. Found in the review of PR 8 (phase 6-quater) | behaviour | `test_the_axis_is_not_written` | fixed in PR 8 (Guido's choice): refused (`InvalidParameter`), the axis is changed by `x_convert`; a correction of the wavelengths (`λ * (1 + v/c)`, air to vacuum) would need the three columns together, a recipe of its own if wanted |
| F128 | expressions (functions) | `_NE_GLOBAL_DICT` (`core/spectrum_expressions.py`) offers `median`, `mean`, `std` (NumPy, NaN-aware) and `min`, `max` (`nanmin`, `nanmax`) to NumExpr, which cannot call Python functions: `median(F)`, `mean(F)`, `std(F)` failed with "'VariableNode' object is not callable", while `min(F)`, `max(F)` are NumExpr's own reductions (a NaN gives NaN), as `sum` and `prod`. The comment in the code says otherwise. Found in the review of PR 8 (phase 6-quater) | behaviour | `test_function_that_numexpr_does_not_have` | open: the message says which functions there are (PR 8); the dictionary is left as it is (removing it changes nothing, offering the NaN-aware ones would need evaluating them outside NumExpr: Guido's choice); closed in PR 13 of 6-quater (Guido, 13 November 2026): accepted; the comment of `_NE_GLOBAL_DICT` that says otherwise is corrected in phase 7-bis |
| F129 | CSV files (`export_ascii`, the CSV loader) | The CSV files of `export_ascii` held no units, and the loader of CSV files took the first three columns as wavelength in nm, flux and error in `erg/(Å s cm2)`: a spectrum exported after `y_convert` or with the axis in km/s came back with the wrong unit, `cont` was lost, and the file of the whole spectrum (`targets='spec'`: `x,xmin,xmax,y,dy,...`) came back with `xmin` as the flux. `import_ascii` read a column `F` as a new column `F`. Found on PR 8 (phase 6-quater); the last two in the characterization of PR 8-bis | behaviour, data | `tests/test_flux_units.py` | fixed in PR 8-bis of 6-quater (Guido): ECSV with the units and z_em, the names of the GUI in ASCII (`lam`, `F`, `dF`, ...), read back by name; the old files by name (`x`, `y`, `dy`), the flux in arbitrary units; `import_ascii` translates the names and keeps the units |
| F130 | `y_convert` (columns without a unit) | `convert_y_axis` converted only the auxiliary columns whose unit is equivalent to that of the flux: in an archive whose `cont` and `model` have no unit and the flux has one (J0103-1305 of `tests/data`), `y_convert` from Å to nm multiplied the flux by 10 and left the continuum, so the normalized flux changed by 10. Found in the analysis of PR 8-bis (phase 6-quater) | behaviour, numbers | `test_y_convert_converts_the_continuum` | fixed in PR 8-bis of 6-quater (Guido): `cont`, `model`, `cont_pl` without a unit are in the unit of the flux, both when an archive is read and in the conversion |
| F131 | `calibrate` (unit of the result) | `Spectrum.calibrate` multiplies the values by a curve computed in the unit of the flux, and labels the result `erg/(Å s cm2)`: a flux in `erg/(nm s cm2)` is calibrated in erg/nm and labelled per Å, 10 times brighter than it says (SDSS_g = 15.0 for 17.5 asked). `test_spectrum_calibrate` (`tests/test_photometry.py`) recorded it ("the unit is relabelled"). A flux without a unit or in `arb_flux` is read as `erg/(Å s cm2)`, and its result is right. Found in PR 8-bis (phase 6-quater) | behaviour, units | `test_calibrated_spectrum_has_the_magnitude_asked_for` | fixed in PR 8-bis of 6-quater (Guido): the unit of a physical flux is kept, `erg/(Å s cm2)` for an arbitrary one or one without a unit; the values unchanged |
| F132 | `equalize_and_stitch` (units of the flux) | The arms were equalized on the values of their flux, then joined by `stitch`, which converts to the unit of the bluest: an arm in another physical unit (per nm against per Å) came out off by the ratio of the units (10 times too faint). With `manual_factors` the result was right (the factor applies to the values in the unit of the arm). Found in the analysis of PR 9 (phase 6-quater) | behaviour, units | `test_equalize_and_stitch_converts_the_flux_before_scaling` | fixed in PR 9 of 6-quater (Guido): the arms converted to the units of the bluest before the scaling; same results with arms in one unit (the fingerprints of `test_valid_results_unchanged`) |
| F133 | `coadd` (units of the flux) | The values of the flux were co-added as they were and labelled with the unit of the first session: with `equalize_order = -1` a session per nm was 10 times too bright against one per Å, with the scaling right by chance (the factor absorbs the ratio); a flux in arbitrary units or normalized was co-added with a physical one without a word. Found in the analysis of PR 9 (phase 6-quater) | behaviour, units | `test_coadd_converts_the_flux` | fixed in PR 9 of 6-quater (Guido): converted to the units of the first; an arbitrary flux with a physical one only when scaled (it takes the unit of the reference), a normalized one with one that is not refused |
| F134 | `coadd`, `equalize_and_stitch` (units of the axis) | The overlap for the scaling was looked for in the values of the axis (`compute_flux_scaling` takes values): a session in Å against one in nm had none, so it was not scaled, without a word; `equalize_and_stitch` also applied the cuts and the ranges (nm) in the unit of each arm (an arm in Å kept whole, its pixels twice in the overlap) and ordered the arms by the values of their axis (an arm in Å after one in nm, a part lost). Found in the analysis of PR 9 (phase 6-quater) | behaviour, units | `test_an_arm_in_angstrom_is_equalized`, `test_equalize_and_stitch_orders_the_arms_in_one_unit` | fixed in PR 9 of 6-quater (Guido): the axis converted to the unit of the reference, the arms ordered in nm, the cuts and ranges converted to the unit of the axis; an axis in velocity against one in wavelength refused (`x_convert` first) |
| F135 | recipes that read other sessions, without a GUI | In a script (a session of `load_session_from_file` with `gui=None`) `coadd`, `equalize_and_stitch`, `import_systems` and `trim_common` failed with an `AttributeError` (`None.session_histories`): a crash, reported to Guido at once (PR 9); `stitch` worked only with sessions given as objects, and its log held their `repr`. Found in the analysis of PR 9 (phase 6-quater) | crash | `test_without_a_gui_a_name_fails_with_a_hint`, `test_in_a_script_the_sessions_themselves` | fixed in PR 9 of 6-quater: one rule for the parameters of kind `session` and `sessions` (`recipes/_sessions.py`): names, or in a script the sessions themselves; the name of the session the recipe runs on is that session; the log records names |
| F136 | recipes that read other sessions (a name not found) | Each did its own: `coadd` skipped a name not found without a word (and co-added the others), `stitch` with a WARNING, `trim_common` with a WARNING and then "no range in common" (or "N sessions and M transitions"); `trim_common` also gave "no range in common" for an unknown transition of its session and ignored one of another session (an infinite range). Found in the analysis of PR 9 (phase 6-quater) | behaviour | `test_a_name_not_found_fails_naming_it`, `test_refused` | fixed in PR 9 of 6-quater: `RecipeFailed` naming it; an unknown transition `InvalidParameter` |
| F137 | `coadd` (this session) | The sessions co-added are those named: without this session in the list, the result held the others only, in a new state of this session (its name, its systems), while the description said "with this one"; the first named, not this session, is the reference (units, scaling). Found in the analysis of PR 9 (phase 6-quater) | behaviour, interface | `test_coadd_refused` | fixed in PR 9 of 6-quater (Guido): a list without this session refused; the first named stays the reference (described) |
| F138 | `flux.equalize` (units) | The same as F132 and F134 in `Spectrum.equalize_to_reference`, used by `flux.equalize`: a session per nm equalized to one per Å kept its unit with the values of the other (10 times too faint); a session with the axis in Å was not scaled to one in nm (no overlap of the values), without a word. The core of `Spectrum.coadd`, called by a script, had F133 and F134 too. Found in PR 9 (phase 6-quater) | behaviour, units | `test_equalize_compares_the_flux_in_one_unit`, `test_equalize_compares_the_axis_in_one_unit`, `test_the_core_in_a_script` | fixed in PR 9 of 6-quater (Guido, added before the merge): the core compares in the units of this spectrum (`_in_units_of_this` of `core/spectrum_combine.py`: the reference converted for the scaling, the result in the unit of this flux; an axis in velocity against a wavelength a `ValueError`), `coadd` of the core converts the flux to the unit of the reference; `flux.equalize` applies the rule of `recipes/_combine.py` (an arbitrary flux takes the unit of the reference, a normalized one refused); same results with one unit (fingerprints) |
| F139 | `fit_bayesian` (nested sampling) | With dynesty 3.2 the line of the log on the progress (`bayesian_fitter.py`, `run_nested`) formats `dsampler.dlogz` as a number, and it is a dict: a `TypeError` as soon as the sampling lasts more than 0.5 s, so almost always (tried with `nlive` 10 and 250). The GUI hides the action (commented out), the MCP tool and the scripts offer it. Found in the analysis of PR 10 (phase 6-quater) | crash | — | fixed in PR 10-bis of 6-quater: the remaining evidence read from the iteration (`delta_logz`), not from the sampler; the nested sampling runs (`test_nested_sampling_runs_and_logs_its_progress`) |
| F140 | `fit_bayesian` (nested sampling) | The medians and standard deviations are taken from `res.samples`, the dead points of dynesty, without their weights: they are not a sample of the posterior (`samples_equal()`, or `resample_equal` with the weights), so the values written to the components are wrong when the sampler does not fail (F139). Found reading the code in PR 10, not measured | numerical | — | fixed in PR 10-bis of 6-quater (Guido): the medians are the weighted quantiles (`dynesty.utils.quantile`), the standard deviations weighted, the samples stored for the corner plot `samples_equal()`. Measured on a CIV component of the mock of `test_recipe_inputs_absorbers` (3 free parameters): without the weights dlogN 1.68 and db 42 km/s, with them 0.005 and 0.14, the MCMC 0.0049 and 0.14; the medians within 0.002 in logN and 0.3 km/s in b (`test_nested_posterior_is_weighted`) |
| F141 | `fit_bayesian` (MCMC) | The number of walkers is fixed at 32: emcee needs at least twice the free parameters, so a group with more than 16 (6 components with z, logN and b free) fails with a `RuntimeError` ("It is unadvisable to use a red-blue move with fewer walkers than twice the number of dimensions"), which the boundary of the recipes does not convert (it reaches the user as a bug). Found in PR 10 | crash, interface | — | fixed in PR 10-bis of 6-quater (Guido): `max(32, 2 * ndim)` walkers, so the same sampling up to 16 free parameters (the snapshots of PR 10-bis); fewer walkers than `2 * ndim` refused by `run_mcmc` (`ValueError`); a `RuntimeError` of a sampler is a `RecipeFailed` (`test_mcmc_with_more_than_16_free_parameters`) |
| F142 | manual fits (unknown components) | A UUID that names no component was taken: `update_component` logged a WARNING and recorded a state that changed nothing, `update_components`, `delete_component`, `fit_component` and `optimize_system` recorded one silently, `update_constraint` stored a constraint for no component, `fit_bayesian` failed with "cond is not defined on empty arrays". Found in PR 10 (phase 6-quater) | interface | `test_refused_before_the_recipe_runs` | fixed in PR 10 of 6-quater (Guido) for the first five (`check_uuids` of `recipes/_components.py`); `optimize_system` and `fit_bayesian` in PR 10-bis (`uuids` of `fit_bayesian` required) |
| F143 | manual fits (unknown series) | A series not in the atomic data was stored: `update_component(s)` gave it to the component, `add_linked_system` added a component without lines (`'CIV,nope'`), `add_component` added it and failed with "No valid pixels" (also for a list of series, `'CIV,SiIV'`, read as a multiplet of two transitions that do not exist). Found in PR 10 | interface | `test_refused_before_the_recipe_runs` | fixed in PR 10 of 6-quater (Guido): refused (`check_series`); the list of series names `add_linked_system`; a series given twice to `add_linked_system` refused |
| F144 | registry and methods (defaults) | Defaults of the registry that the method did not have: `detect_anchor` (`trans_self`, `trans_sibling`, the OVI doublet), `add_linked_system` (`series_list`, `z`), `update_constraints_column` (`is_free`). The dialogs and the MCP schemas offered them; an agent or a script that left them out got a `TypeError`. `test_recipe_defaults.py` skips the parameters the method requires, so nothing saw it. Found in PR 10 | interface | `test_typed_recipe_defaults_are_those_of_the_signature` | fixed in PR 10 of 6-quater (Guido): required on both sides, also `z` of `add_component` (its 0.0 never fell on the spectrum) |
| F145 | `detect_anchor` | `Spectrum.detect_doublet_z` took the depth as `1 - y`, right only for a normalized flux (with a flux of 1e-16 every pixel was "absorbed" and the anchor was the start of the common range); it read the axis as nm whatever its unit; it returned 0.0, as a redshift, when the two transitions cover no redshift in common. Found in PR 10 | numerical | `test_detect_anchor_reads_the_normalized_flux` and the two after | fixed in PR 10 of 6-quater (Guido): `y / norm`, the axis in nm, a `ValueError` (a `RecipeFailed`) without a common range or a continuum; the same anchor on a normalized flux |
| F146 | `update_constraint` (links) | `target_uuid`, which puts the target in the group of the fit of the linked component, was not checked against the expression, and was left empty by a link made without it (scripts, agents: the target was fitted without it); a link to the component itself was accepted; an expression with `is_free=True` was dropped without a word; with `is_free` left out the link was stored free (the fit froze it). Found in PR 10 | behaviour, interface | `test_link_without_target_takes_that_of_the_expression`, `test_refused_before_the_recipe_runs` | fixed in PR 10 of 6-quater (Guido): the target is the first component the expression names, a different one refused, the other cases refused, the link stored frozen; the links of the GUI unchanged (fingerprints) |
| F147 | `optimize_system` (`max_components`) | The loop of `optimize_hierarchy` runs `max_components` times, and a candidate too close to a component (`min_dv`) or rejected by the AIC takes one turn: `max_components` is the number of candidates tried, not of the components added, as its description said. Found in PR 10-bis (phase 6-quater) | interface (documentation) | – | documented in PR 10-bis (Guido): the descriptions of the registry, the method and the core say so; the behaviour unchanged |
| F150 | `fit_bayesian` (nested sampling, `nlive`) | With few live points dynesty fails or does not end: with 3 free parameters `nlive` 1 an `IndexError` (not converted by the recipes), 3 a `RecipeFailed` ("zero-size array"), 4 and 7 no end (12 and 2.5 minutes, stopped); 10, 20 and 250 end. Its time also varies with the seed, 7-45 s for one CIV component with any `nlive` from 30 to 250 (the warnings of dynesty on the enlargement of its bounds), the values the same within their errors. Found in PR 10-bis (phase 6-quater) | crash, performance | `test_refused_before_the_recipe_runs` | fixed in PR 10-bis (Guido): at least 10 live points per free parameter (`NLIVE_PER_PARAMETER` of `recipes/absorbers.py`), `min` 10 in the registry; the time left as it is (the GUI does not offer nested sampling) |

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

---

## 7. Progress (long rows of the concluded phases)

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
| 5 – Splitting large files (and GUI performance) | done (October 2026): order of the phase agreed, F12 moved to 6-bis; GUI benchmark (`benchmarks/bench_gui.py`, #634); inspector table (first opening 3.7 → 1.3 s, #635); ticks of the main plot (redraw 0.19 → 0.13 s, #636), found F66; draws of the velocity plot coalesced (first opening 1.43 → 1.19 s, #637); Qt enums of the table read once (first opening 1.22 → 0.54 s; from 3.96 s before the phase, #638); highlights of the main plot blitted (selection of a row 0.32 → 0.22 s, #639); split of `SpectrumV2`: continuum (#640), found F67; lines (#641), found F68; calibration fixed (`np.trapz` removed from NumPy, #642), found F69, F70; flux (#643); combination (#644), found F71, F72; expressions and masks: split of `SpectrumV2` complete (2826 → 774 lines, #645); tests in parallel (#646); split of the GUI: F7 fixed (#647); inspector table (#648); velocity plot (#649); its rebuild and layout faster (selection 0.23 → 0.18 s, #650); parts of the main plot in their modules, found F73 (#651); continuum editor in a mixin, found F74 (#652); mouse handlers and tooltip in a mixin (#653); `plot_spectrum` split into its parts, found F75, F76 (#654); split of `main_window.py`: accessory widgets and dialogs (#655), menus, found F77 (#656); panels, found F78 (#657); controls of the continuum editor, F74 fixed (#658); sessions, found F79 (#659); recipes and scripts, found F80-F83: split of `main_window.py` complete (4831 → 541 lines) (#660); F82, F83 fixed (#661); F66 fixed, found F84 (#662); closing: documentation of the API (the methods of `SpectrumV2` back on its page), result of the phase (System Inspector opened ~7× faster, selection ~3×; `spectrum.py` 2826 → 774 lines, `main_window.py` 4831 → 541; coverage 67% → 81%) |
| 5-bis – GUI look and ergonomics | done (November 2026): order of the phase agreed (all the findings proposed, F73, F76, F77, F80 included); F67, F71 fixed (PR 1, #664); F69 (PR 2, #665); F70 (PR 3, #666); #630 (PR 4, #667); F81, F75 (PR 5, #668); F78 (PR 6, #669); F85, the boxes of the System Inspector (PR 6-bis, #670); F77, F80 (PR 7, #671); F73 (PR 8, #672); F76 (PR 9, #673), found F86; F86 (PR 9-bis, #674); theme module, PR A (no visible change; Select icon, #675), found F87; PR B (dark mode, F87, #676); the toolbar of the plots redesigned (Guido's request, #677); View > Appearance and View > Show Toolbar, F88 (#678); the visual review, found F89, F90; the bugs of the review (F89, F90, the links and placeholders of the dark palette) and F91 (reported by Guido), found F92 (#679); the pass on the texts, with F92 (#680); the main plot (names of the emission lines in spectroscopic notation on a white ground, the ticks of the systems in two rows and softer, no title, the unit of the flux in the usual order; agreed with Guido, #681); the System Inspector (the error of z in km/s, btur hidden while not in use with a menu of the columns, the boxes above the plot in three groups on one row, to keep the plot aligned with the table: Guido; #682); the System Inspector without a continuum or systems (a message with a button instead of empty panels, the fits disabled, the transitions typed shown without a component, #683); the status bar (outcome of the last recipe, coordinates moved from the toolbars) and the keys of the System Inspector (#684; with it, the axes of the main plot on whole pixels: the constrained layout moved them by ~1e-13 px between draws); the toolbar floating over the plot and fading as the cards (Guido's request; as a card with its tab, as wide as the cards, #685); the minor points of the review: the Identifications window as a table, the formats of the Data Inspector (found F93: the identification by doublets, for phase 7; #686); the main window and the menus (the welcome screen with Open Session and the recent files, File > Open Recent, the window geometry remembered, Import Systems in Absorbers, several files loaded in the background with their progress, Guido's request; #687); the script console left as it is (Guido: kept in the code until it is clear whether it is needed with the MCP server, or can be given another use); the colours of the plots (Guido's choice after a comparison on the main plot and the System Inspector, with simulations of colour blindness: View > Plot Colours, the standard palette with a thinner model, or a colour-blind safe one from the Okabe-Ito palette, where model, selection and residuals of the velocity plot no longer turn one olive for a deuteranope; with it, the CIV/SiIV system of the screenshots of the System Inspector refined with Optimize System, Guido's request; #688); closing: result of the phase in "Result of the phase" |
| 6 – Error handling | done (November 2026): plan, inventory of the broad except clauses, rule of the log messages (PR 1, #691); F94 (PR 2, #692); loggers per module, F96 (PR 3, #693); rows N and R of core, fitting and io and their messages, found F97 (PR 4a, #694); their rows D (PR 4b, #695); rows N, R and B of the recipes and their messages, found F98 (PR 5a, #696); their rows D (PR 5b, #697); prerequisites of the fits, F99, the scheme of the prerequisites for 6-bis (PR 5c, #698); rows N and C of the GUI and its messages, `V1LogArtifact` removed (PR 6a, #699); rows B and D of the GUI and of the launcher, the hook for unexpected errors (PR 6b, #700); the warnings of the libraries, warnings as errors in the tests (PR 7, #701); closing: result of the phase in "Result of the phase" |
| 6-bis – Recipe API | done (November 2026): plan of the phase (PR 0, #705); the typed registry, F79, half of F12 (PR 1, #706); typed parameters of continuum, flux and file, the log of the values used (PR 2a, #707); Edit, F46 and the bounds (PR 2b, #708); Absorbers, all the recipes typed (PR 2c, #709); exceptions, continuum, flux and file (PR 3a, #710); edit (PR 3b, #711); absorbers, no recipe returns 0 (PR 3c, #712); structured result, `run` (PR 4a, #713); region and score of the candidates, `not_fitted`, the reasons of the removals (PR 4b, #714); F102, the notes of a result in the status bar and the log (PR 4c, #715); prerequisites, core and scripts, F103 (PR 5a, #716); the GUI asks for the continuum, F104 (PR 5b, #717); descriptions reviewed, F48, the parameters not for agents (PR 6a, #718); the dialogs with the longer texts (PR 6a-bis, #719); the mask and the margin on z_em of the pipelines (PR 6b, #720); the kinds of the parameters, F105 (PR 6c, #721); the custom multiplets without the global table, F12 (PR 7, #722); the MCP pilot on the new API, the generator of the schemas of the tools (PR 8); closing: result of the phase in "Result of the phase" |
| 6-ter – Full MCP server | done (November 2026), PRs #726-#747 and the closing PR 8; release 2.6.0 published on GitHub on 10 November 2026 (PR 9, #751: version, CHANGELOG, `switcher.json`); plan of the phase agreed (PR 0, #726); the tools generated from the registry, F106 found (PR 1, #729); F106 fixed (PR 1-bis, #730); the tools of the sessions, undo and redo (PR 2, #731); the tools of inspection (PR 3, #732); progress notifications (PR 4, #733); the agent's decisions in the log (PR 5a); the gaps of the recipes (PR 5b); the log as the text of a methods section (PR 5c); the GUI following a session of the server, F109 found (PR 5-bis); the headless Docker image, F110 found (PR 6, not published); trials with fresh agents, F111-F116 found (PR 7, report in section 4.3); heartbeat and reads during a recipe, F113 in part (PR 7-bis); trial 3, the editing tools, F117-F120 found (PR 7); the outputs of trial 3 (PR 7-ter); F111 and the cancellation, F113 closed, F121 found (PR 7-quater) |
| 6-quater – Review of the recipe parameters | done (November 2026): plan of the phase agreed (in PR 1, ruff `B`); PR 2, the aliases of the renamed parameters; PR 3-4-bis, mypy on the whole package; PR 5, `SessionV2` and `SessionHistory` into `astrocook/session/`; PR 6, the `V2` suffixes dropped (#759); PR 7, the review of file and flux, F122-F125 found (#761); PR 8, the review of edit on columns and expressions, F120 fixed (one set of names of the columns), F126-F128 found, F127 fixed; PR 8-bis, the units of the flux (`arb_flux`) and the CSV files in ECSV, F129-F131 fixed; PR 9, the review of edit on several spectra, one rule for the sessions, F123 and F132-F138 fixed; PR 10, the review of the manual fits (components and constraints), F117 closed, F142-F146 fixed, F139-F141 found (the samplers, PR 10-bis); PR 10-bis, the review of `optimize_system` and `fit_bayesian`, F139-F142 and F150 fixed, F147 documented, F148-F149 registered; PR 11, the labels of the menus from the registry (`menu` key, F107) and the pages of the docs for the *Help* button (F125), F151 found; PR 12, the MCP server: the shared descriptions shorter (-34%), the workflows in the instructions, the pid and version of the marker in the banner, F152 found; PR 13, the closing: the findings without a phase reviewed with Guido (F15, F28, F30, F101, F114, F115, F151, F152 to phase 7-bis; F52, F102, F112, F116 and the stale ones to phase 7; F17, F126, F128 closed), "Result of the phase", the rules in force, `CLAUDE.md`, the detail moved here; release 2.7.0 follows (PR 14) |
| 7 – Robust components | not started |
