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
4. **Non-determinism:** fixed seeds, statistical tolerances, UUIDs and
   timestamps excluded from comparisons.
5. **File compatibility:** `.acs`/`.acs2` round trips (requires real V1
   archives).
6. **Measuring the net:** coverage; mutation testing (`mutmut`) on critical
   modules.
7. **GUI:** smoke tests with `pytest-qt` in offscreen mode.
8. **CI:** GitHub Actions (Python 3.12, 3.13). Done (brought forward): tests on both versions;
   Ruff blocking on undefined names (`F821`) only, until the other pyflakes
   warnings are cleaned up in a dedicated PR.
9. **First scientific checks:** oversampling of the Voigt model; scaling of
   uncertainties by the reduced χ².
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
    (F16). Next: Voigt model, continuum, identification (statistical).

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
different grids, noise, flux scales and resolutions. Next: the V1/V2 parity
tests (item 11), then the full chain.

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

### 3.1 User-reported issues
GitHub issues are the single source of truth; this section only says *when*
each open issue is addressed. Rule of thumb: bugs before refactoring,
features after, except urgent features in stable areas.

| When | Issues | Notes |
|---|---|---|
| Phase 1 | [#469](https://github.com/DAS-OATs/astrocook/issues/469), [#477](https://github.com/DAS-OATs/astrocook/issues/477), [#455](https://github.com/DAS-OATs/astrocook/issues/455), [#450](https://github.com/DAS-OATs/astrocook/issues/450), [#436](https://github.com/DAS-OATs/astrocook/issues/436) | "possibly fixed": verify with tests, close or reopen |
| Right after phase 1, before phase 2 | F1 ([#490](https://github.com/DAS-OATs/astrocook/issues/490)), F3, F10 (section 3.2) | V2 regressions (V1 is correct): rebinned errors, auto fudge, convolution with a resolution column. Reference behaviour: V1; each fix still agreed with Guido, with a regression test |
| Right after phase 1, before phase 2 | [#465](https://github.com/DAS-OATs/astrocook/issues/465) | continuum sensitive to spikes (urgent); numerical choice, to be agreed with Guido before fixing |
| Phase 0 (brought forward) | [#476](https://github.com/DAS-OATs/astrocook/issues/476), [#481](https://github.com/DAS-OATs/astrocook/issues/481) | ASCII export; same problem (missing `import re`, found by Ruff), fixed with a regression test |
| Right after phase 1, before phase 2 | [#442](https://github.com/DAS-OATs/astrocook/issues/442) | scatter plot not auto-triggering on region overlap |
| Right after phase 1, before phase 2 | [#480](https://github.com/DAS-OATs/astrocook/issues/480) | import of line lists as ASCII (urgent); touches `io/`, little affected by refactoring |
| With the natural-language log work | [#424](https://github.com/DAS-OATs/astrocook/issues/424) | logging of manual continuum edits |
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
| F1 | `rebin_spectrum` | Errors do not decrease with the number of pixels per output bin (σ instead of σ/√N) | numerical | `test_rebin_weighted_average_error` (strict xfail) | issue [#490](https://github.com/DAS-OATs/astrocook/issues/490); V2 regression (V1 gives σ/√N), scheduled as a bug (section 3.1) |
| F2 | `rebin_spectrum` | Sigma clipping uses the first-pass bin mean, which includes the outlier: with few pixels per bin (e.g. 3) a single spike can make the whole bin NaN. Clipping runs once, docstring says "iterative". Not a pure V2 regression: V1 clips with `sigma_clip` within each bin, which keeps the bin but cannot remove an isolated spike with 3-5 pixels per bin; neither version removes it. Related to #465 | numerical | snapshot `rebin_spectrum_nm_clip`; `test_rebin_clipping_keeps_bins_with_a_spike` (strict xfail) | open |
| F3 | `fit_continuum` | `fudge='auto'` (recipe default) is always exactly 1.0: computed on residuals that are zero by construction, before the smoothing. Related to #465 | numerical | snapshots `continuum_j0103_*` | V2 regression (V1 uses a running mean, so the fudge is effective and local), scheduled as a bug (section 3.1) |
| F4 | `VoigtFitterV2` | z is bounded to ±`z_window_kms` (default 20 km/s) around the initial guess; a guess further away cannot reach the true value. Intended? | behaviour | `test_fitting.py` (guess inside the window) | open |
| F5 | `convert_y_axis` | Annotation declares 2 return values, the function returns 3 | code quality | snapshot `convert_y_axis` | open |
| F6 | architecture | `core/spectrum_operations.py` imports `gui.debug_utils`: core depends on the GUI | architecture | — | phase 3/4 |
| F7 | GUI | `gui/main_window.py` imports `matplotlib.pyplot` before `pyside_plot.py` selects QtAgg: importing the GUI without a `QApplication` (and no display) fails | architecture | workaround in `tests/test_gui_save_session.py` | phase 4 |
| F8 | resolution handling | The resolution lives in five places (per-component value, `resol` column, `SpectrumDataV2.resol`, meta `resol`/`RESOL`, meta `resol_fwhm`), each read with a different priority by the fitter, the `.acs2` loader, coadd and the GUI. The `.acs2` save does not write the attribute: it survives only through the meta keys (set by "Set properties"). The fitter stamps the resolution used into every fitted component, which then overrides later changes. R vs km/s is guessed from a threshold of 500. The J0103-1305 archive (older version) stores no resolution at all. To be analysed by Guido | behaviour, possible data loss | — (map in the PR that added this row) | open, priority |
| F9 | Voigt model | Components narrower than a pixel (e.g. Ly-alpha with b = 0.5 km/s, logN = 17.8 in J0103-1305) make the model depend on the sampling; with R = 45000, 1% of the 570-700 nm window still differs from the model saved in the session. Related to the oversampling check (item 9) | numerical | — | open |
| F10 | `convolve_flux` | With a per-pixel resolution column, pixels are grouped by resolution rounded to 100 (R) or 0.5 km/s, but selected by the raw value: pixels whose resolution is not already rounded are never convolved and the model is 0 there (a flat column at R = 45050 gives a zero model everywhere). Affects rebinned and coadded spectra (coadding three synthetic exposures at R = 42000-47500 gives a non-round column on 90% of the pixels and a zero Voigt model on 90% of them, `test_voigt_model_on_coadded_spectrum_is_not_zero`), RESVEL input, non-round "Set properties" values. In addition, even for a constant R the kernel width is computed at the median wavelength and with the median pixel size, so it is wrong far from the centre of the range (+53% / -26% at the ends of 380-780 nm). Proposed: a single exact convolution in velocity space, after the safety net is complete | numerical, wrong results | `test_model_column_resolution_is_not_zero` (strict xfail) | V2 regression (V1 takes R at the system centre, no rounding; from reading the code), scheduled as a bug (section 3.1), high priority |
| F11 | grouping | `VoigtModelConstraintV2.set_active_components` ignores `group_depth`: it always calls `get_connected_group` with the default depth 2, so the recipe option (1 = neighbours only) and the `group_depth=0` used by the refit have no effect | behaviour | — | open |
| F12 | global state | Module-level dictionaries are modified at run time: the GUI rewrites the parameter docs of `EDIT_RECIPES_SCHEMAS["equalize_and_stitch"]`, and `add_component` registers custom series in `STANDARD_MULTIPLETS`. The changes persist for the whole process (and leak between tests) | architecture | `restore_multiplets` fixture in `test_issue_455_multiplet_subset.py` | open |
| F13 | inspector | `SystemTableModel.set_highlighted_uuids` returns early when the selection is unchanged, so link colours are not recomputed when a link is created with the group still selected (model-level reproduction; full GUI not driven) | GUI | `test_link_created_with_group_selected_is_coloured` (strict xfail) | open |
| F14 | grouping | Two components are considered blended when their line centres are closer than max(4 b, 30 km/s), independently of the instrumental resolution and of the line strength. At X-Shooter resolution (R ~ 8900) AlII 1670 and CIV 1550 still overlap at 60 km/s but are not grouped; worse in the NIR or for saturated lines. Possibly why #469 is still seen on X-Shooter data | behaviour | `test_blend_threshold_ignores_resolution` (characterization) | open |
| F15 | `equalize_and_stitch` | The scale factor is a median of the ratio of the two arms in the overlap: with narrow absorption (e.g. tellurics) seen at two different resolutions, it is biased by ~1-1.5% in a mock with equivalent-width-conserving tellurics around 1000 nm (no bias without tellurics). Possibly related to #477; to be checked on real X-Shooter VIS/NIR data | numerical | — (mock in the PR that added this row) | open |
| F16 | smoothing | For the same "std" parameter, V1 `_gauss_convolve` uses the kernel exp(-(dv/std)^2), i.e. sigma = std/sqrt(2), while V2 `smooth_spectrum` uses sigma = std: V2 smooths sqrt(2) wider. Affects the continuum smoothing (default 500 km/s: effective sigma 354 km/s in V1, 500 km/s in V2) and any other use of the parameter. Which convention to keep is Guido's decision | numerical | `test_smoothing_width_differs_by_sqrt2` (characterization) | open |

### Rough estimates
10–20 work sessions, 2–5 calendar weeks; the bottleneck is human review and
the availability of test data.

### Decisions reserved to Guido
- License.
- Questionable behaviours found by characterization tests (bug or intended?).
- Approach to the continuum spike sensitivity (#465).
- Scientific validation of results.
- Releases (PyPI, macOS app signing). In particular: whether to publish a
  corrective 2.0.2 once the V2 regressions (F1, F3, F10) are fixed.

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
| 1 – Safety net | in progress |
| Open bugs and urgent features (section 3.1) | not started |
| 2 – Cleanup and packaging | not started |
| 3 – Decoupling from legacy | not started |
| 4 – Splitting large files | not started |
| 5 – Error handling | not started |
