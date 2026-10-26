## [Unreleased]

### Added
- Documentation: the screenshots of the tutorials are now made by a script, `docs/screenshots/make_screenshots.py`, which runs the GUI offscreen (one scene per image, see `docs/screenshots/README.md`), and all of them show the current interface.
- Voigt fit: each fitted component records the warnings of its last fit, in a new `flags` field saved in the session archives (#517): a parameter at a bound of the fit (`z_at_window`, `logN_at_lower_bound`, `logN_at_upper_bound`, `b_at_lower_bound`, `b_at_upper_bound`; the upper bounds of the automatic pipelines are **Max logN** and **Max b**), or uncertainties that could not be computed (`nan_errors`). Editing a parameter of the component clears them. *Refit All* and *Refit Selected* report in the log how many components have warnings. For agents, `absorbers_doublets_auto` of the MCP server returns these flags, renamed accordingly (`b_at_max_b` and `logN_at_max_logN` are now `b_at_upper_bound` and `logN_at_upper_bound`). Earlier versions of Astrocook ignore the field. No change of results.

- MCP server for AI agents (pilot): Claude Code, the Claude desktop app or other MCP clients can open a spectrum, estimate the continuum, identify lines, fit the metal doublets and save the session, through the `astrocook-mcp` command (install with `pip install "astrocook[mcp] @ git+https://github.com/DAS-OATs/astrocook.git"`). The tools return structured summaries and Astrocook's error messages, and record their recipes in the session log, so that the result can be reviewed in the GUI. See the new tutorial *Using Astrocook from an AI Agent*. No existing behaviour changes.

- System Inspector: a **Fit ±** field, after logN and b above the velocity plot, sets how far (km/s) the redshift of each component may move in the fits run from the inspector (adding a component in the velocity plot, Refit Selected, Refit All); left empty, each fit keeps its window (20 km/s when adding, 100 km/s when refitting), so nothing changes unless it is set (#518).

- Documentation: a new page, *Atomic Data*, describes the table of atomic data, the source of each value, the choices made and the series and multiplets that Astrocook knows (#568).

### Fixed
- GUI (Linux, Windows): the unchecked boxes of the plot controls were invisible, and the checked ones showed only the tick; macOS was not affected (F64).
- GUI: the card of the plot controls could cut the last line of the hints of the continuum editor.
- Voigt fit: sessions saved by earlier versions with pixel edges (`xmin`, `xmax`) without a unit could no longer be refitted, with the error "'' (dimensionless) and 'Angstrom' (length) are not convertible" (#623, a regression of the exact convolution of #526, which uses the edges). Edges without a unit now take the unit of the wavelengths when a session is opened; if their values are not around the wavelengths, they are recomputed halfway between the pixels. No change of results for other sessions.
- System Inspector: the profiles drawn for the selected components and under the cursor now use the physics and the convolution of the fit (#526): the same Voigt function, and each pixel convolved with its own resolution and averaged over its edges. Before, they had their own copy of the Voigt function and a kernel taken at the median wavelength, so they could differ from the model of the fit at low resolution or for narrow lines. No change of results.
- Voigt fit and System Inspector: the damping wings of the lines were (1 + z) times too strong, because the damping parameter of the Voigt profile was computed with the observed wavelength of the line instead of the rest one (F61). **Results change** for absorbers whose damping wings matter: the column density of a damped Ly-alpha or Lyman-limit system fitted from its wings was too low by log(1 + z), 0.5-0.7 dex at z = 2.5-4 (on a mock damped Ly-alpha with logN = 20.3: 19.75 at z = 2.5 and 19.64 at z = 3.5 before, 20.29 now). Narrow saturated metal lines change by up to 2% in depth. Astrocook 1 was not affected. Systems fitted with earlier versions should be refitted.
- Voigt fit: the model is now convolved with the resolution of each pixel and averaged over the pixel (#526). It is computed on a grid finer than the pixels, uniform in velocity, then each pixel takes it convolved with a Gaussian of its own resolution and averaged over its edges. Before, the model was sampled at the pixel centres and convolved with one kernel per fitting window, its width taken at the median wavelength (up to +53%/-26% wrong at the edges of a wide window) and a resolution column that varied was rounded to 100. Against the exact flux of mock spectra the model was off by up to 0.6% of the continuum at R = 45000 and 6% at R = 8900; now the difference is below 0.003% (with the correction of the damping parameter, see above). **Results change**, mostly at low resolution and for lines narrower than a pixel: on the validation mocks at R = 8900 the spread of the column density of the systems goes from 0.18 to 0.10 dex, the median reduced chi2 from 2.45 to 1.61, the components with suspicious errors from 30 to 16; at R = 45000 and above the results are the same within 0.002 dex. The fit is not slower: a full automatic chain is 25-30% faster, a single blended group up to 15% slower.
- GUI: the System Inspector read the resolution column always as a FWHM in km/s, so a column of resolving powers (R) gave profiles with the wrong instrumental width; it now reads it by its unit, as the fit does (#526). The session panel, the dialogs and the checks use the same resolution (the median R when it varies with wavelength); *Set Properties* asks before replacing a resolution only when it varies with wavelength. In tooltips, a component without its own resolution shows "spectrum".
- Automatic continuum (*Auto-estimate continuum*, *Find absorbed regions*, and the pipelines that use them): upward spikes (cosmic rays, bad pixels, frequent in NIR spectra) no longer pull the continuum up; with 1% of the pixels in spikes it came out 2 to 20 times too high, as reported for X-Shooter and FIRE spectra (#465). At S/N ~3 the continuum on the red side is no longer 26% too high. **Results change**: the running mean of the clipping is now taken over the smoothing length in km/s, exact at the edges (it was biased low there), and pixels with no signal and spikes are flagged in `abs_mask` as not continuum. On the validation mocks the forest continuum moves closer to the true one (0.936 → 0.943 of it; 5th percentile 0.820 → 0.830), the red side is unchanged, and more weak metal doublets are found (CIV completeness 0.857 → 0.874; systems detected at 5-10 sigma 35% → 44%). The forest continuum is still ~5% low at S/N 30.
- Automatic pipelines (`doublets_auto`, `lya_auto`, `components_auto`): the absorbed mask of the session (`abs_mask`, e.g. from *Auto-estimate continuum* with other parameters) is no longer replaced by the mask the pipeline makes for itself, nor cleared outside the region of the pipeline (the forest after `doublets_auto`), so a later *Fit continuum* uses the mask of the user (F50). The results of the pipelines do not change.
- Atomic data: the wavelengths of OIII 832, PIII 913, PV 1117 and FeI 2501 were wrong by 0.2-0.7 A (mistyped in the conversion from nm to Angstrom, e.g. 913.397 instead of 913.968 A), and the damping constant of ArI 1066 was ten times too small; they now agree with Morton (2003) and with the table of V1 (#568). **Results change** for absorbers fitted or identified with these lines.
- Atomic data: AlII 1670 at 1670.7886 A (was 1670.7874, 0.2 km/s); the forbidden emission lines [NeIII] 3869 and [NeV] 3346 at their vacuum wavelengths, as the rest of the table (they were in air, 1.1 and 1.0 A shorter); `TiII_3066`, a line that no reference lists, has f = 0 as in V1, and the real `TiII_3067` (3067.2379 A) is added (#568).
- Redshift cursor and identification viewer (GUI): with a series of the V1 syntax such as `Ly`, `D` or `H` the redshift refers to the strongest line (Ly-a), not to the Lyman limit or to the weakest line, as the series followed the order of the V1 table; the aliases `Ly-ab` and `Ly-abg` also start from Ly-a (#568, F39).

### Changed
- Automatic pipelines (`components_auto`, `doublets_auto`, `lya_auto`, `forest_auto`, `metals_auto`): a new parameter, **Max logN** (`max_logN`, default `19`), bounds the column density in their fits, as **Max b** bounds the Doppler parameter (F62). In a saturated trough of the forest that the detected lines do not fill, a single component could fill it with its damping wings and grow to logN 20-23, blacking out the model for hundreds of km/s around it. This was already possible before and became worse once the damping wings were corrected (see Fixed). On J0103-1305 the model of the forest is now much closer to the data: chi2 per pixel 31 against 98 before the correction of the wings, and 2-3% of the pixels over-absorbed against 28%. Damped Ly-alpha systems and sub-DLAs are to be fitted by hand: *Refit All* (whose **Max logN** defaults to `23`), the System Inspector and the other manual fits keep the physical limit. **Results change** for the automatic fits of the forest; the metal doublets are unchanged.
- GUI: the session list and the plot controls are now cards as tall as their contents, below the menu bar, that fade out while the mouse is on the plot and come back when it returns; they still float over the plot, so that opening them does not change the aspect of the spectrum. They open with a list icon (sessions) and a cog (plot controls) instead of the arrows. All the sections of the plot controls can be opened and closed with a click on their title, and Astrocook remembers their state. While *Plot Elements* is open, each checkbox shows a sample of how its element is drawn, and those entries leave the legend of the plot, which used to show through the panel. Long session names are shortened in the middle, with the full name on hover. No change of results.
- Voigt fit: components whose lines are closer than the instrumental FWHM are now fitted together (the FWHM from the resolution of the spectrum at the lines, or of the component if set), besides those closer than 4 b or 30 km/s as before (F14). Before, at low resolution (R below ~10000, e.g. X-Shooter) lines of different systems whose profiles overlapped could be fitted one at a time, with the other frozen. **Results change** only for such blends at R below ~10000; the results on the validation mocks are identical (they fit only the sparse metal doublets). The System Inspector highlights the same groups.
- Resolution: the `resol` column of the spectrum is now its only source, R or a FWHM in km/s according to its unit (#526). The resolution stored elsewhere by earlier versions (spectrum attribute, header keywords) is moved into the column when a session is opened, so existing sessions open as before. The Voigt fit no longer writes the resolution it used into the fitted components: a component keeps its own resolution only if set on purpose, and those written by earlier fits are cleared when a session is opened, so that changing the resolution of the spectrum applies to all its components. *Set Properties* keeps a resolution column that varies with wavelength unless asked to replace it, and no longer flattens it when the resolution is left unchanged. No change of results for spectra with a single resolution. Astrocook no longer writes the resolution into the header keywords `resol`, `resol_fwhm` (and `RESOL` in archives): the resolution is saved in the `resol` column only, which earlier versions of Astrocook 2 also read; the keywords of files from earlier versions are still read when there is no column.
- Line identification: no candidate system lies beyond the emission redshift by more than 3000 km/s, a margin for absorbers associated with the quasar and for the uncertainty of z_em (F55). Before, the metal doublets had no limit: a SiIV 16000 km/s beyond z_em, in a telluric band, was identified and fitted. **Results change** only for spectra with such candidates.
- Line identification: a metal doublet gives one candidate system per absorption component, from the stronger line (e.g. CIV 1548); before, the other line gave a second candidate ~1 km/s away, and the automatic pipelines (`doublets_auto`, `components_auto`) added a component for each (F49). **Results change**: on the validation mocks the spurious components drop from 428 to 11 and those with meaningless errors from 413 to 38; the total column of each system stays right, but systems with several close components are now fitted with fewer components (the second candidate sometimes ended on one of them).
- Atomic data: a single table, `astrocook/data/atomic_data.ecsv`, replaces the table in the code and the V1 file `atom_par.dat`, and records the source of each wavelength, oscillator strength and damping constant (#568). Newer values replace those of Morton (2003), as in the atom.dat of VPFIT (decided by Guido): the oscillator strengths of ZnII 2026 and 2062 (Kisielius et al. 2015; 26% higher, **logN(ZnII) 0.10 dex lower**), of NiII 1317-1751 (Boissé & Bergeron 2019; up to 32% lower, **logN(NiII) up to 0.17 dex higher**) and of SII 1250-1259 (Kisielius et al. 2014; 5% lower, logN 0.02 dex higher), and the wavelengths of CIV 1548 and 1550 (Petitjean & Aracil 2004; 0.2 and 0.5 km/s). The damping constant of CI 1560 is the one of Morton (2003), 1.27e8 s^-1 (was 1.17e8). **Results change** for absorbers with these lines. The transitions of the V1 file that V2 never used (H2, CO, Paschen and Brackett series, ...) are not in the table; the V1 series names made only of them (`H2`, `COJ0`, ...) are no longer recognised.
- The defaults of the recipes are now the same in the GUI dialogs and in scripts (#573): before, 15 parameters had a different default in each, so the same recipe gave different results from a dialog, from a script, or from the System Inspector (which passes only some parameters). The values of the dialogs are kept, with three exceptions decided by Guido. **Results change** where a default changes:
  - `refit_all` from a script and from the **Refit** button of the System Inspector: grouping depth 2 and 200 iterations per group (were 0 and 100), as the *Refit All* dialog; slower, and better on blends. The automatic pipelines (`metals_auto`, `forest_auto`, the refit after importing systems) keep the values they always used.
  - `optimize_system` from a script: an AIC improvement of at least 5 to accept a new component (was 0); the dialog now offers up to 5 extra components and a residual threshold of 2.5 sigma, the values of the API.
  - `fit_component`: the dialog now offers 2000 iterations, the value the System Inspector and scripts always used.
  - `rebin`: the dialog now clips upward outliers at 5 sigma by default (was off), as scripts did.
  - `coadd` from a script: a grid of 10 km/s (was 0.01 nm), as the dialog; choose the native pixel size for high-resolution data.
  - `x_convert` from a script: velocities relative to the emission redshift of the session (was z = 0), as the dialog.
  - No change of results: `estimate_auto` takes `fudge='auto'` (the same as `'1.0'`), two parameter names of the dialogs follow the methods (`uuids` of `fit_bayesian`, `file_path` of `import_ascii`), and the maps of `populate_from_identification` default to an empty string, as before treated as missing.

## [2.1.0] - 2026-10-21

First consolidated release (phases 3 and 3-bis of the consolidation): Astrocook no longer depends on the V1 code, which is removed, and V1 `.acs` archives are read by a new native reader. Published on GitHub with the macOS app; not on PyPI yet (install with pip from GitHub). **Some results change**: refitting a V1 session with links between parameters (now kept as constraints), and the `dbtur` column of the fitted components (now updated; 0 when `btur` is frozen). Nothing else changes numerically with respect to 2.0.3.

### Added
- `.acs2` archives record the version of Astrocook that saved them (`astrocook_version` in the metadata file); a session opened from an archive reports it as `saved_with_version` (`None` for archives saved by 2.0.3 or earlier). Older versions of Astrocook ignore it.

### Fixed
- *Data loss*: a session that could not be read (an `.acs2` archive without its metadata file, a damaged `.acs` archive, a file in a format without a loader) was opened with mock data in place of the spectrum, and saving it overwrote the file. An `.acs2` archive without metadata now opens with its spectrum and systems (without constraints, session log and identifications); the other cases report an error.
- V1 `.acs` archives: the turbulent broadening (`btur`) and its error, and the reduced chi2 of each component are now read (they were set to 0 or dropped); the flux keeps its unit, erg/(cm2 s Angstrom), so that unit conversions (Edit > Convert y axis) work on V1 sessions.
- V1 `.acs` archives: the links between parameters (e.g. SiIV z and b tied to CIV) are now kept as V2 constraints. Before, the linked parameters were frozen at their V1 values, did not follow their target in a refit, and the System Inspector showed no target; constraints on logN were lost. **Refitting a V1 session with links gives different (correct) results.**
- Session log (GUI): a recipe that opens a new session (`extract_preset`, `split`, `stitch`, `equalize_and_stitch`, `coadd`, `resample`, `identify_lines`) was also recorded in the log of the session it was run from, whose data do not change; after an undo, it also deleted that log's redo entries. It is now recorded only in the log of the new session. *Note*: `.acs2` archives saved with 2.0.2 or 2.0.3 may still contain such an entry in their log; it is left as it is (it cannot be told apart from a correct one), so check the log before replaying it as a script.
- System Inspector: a link created or removed while its group is selected is now coloured (or uncoloured) at once, and the colours of the other links no longer change when a link is added.
- Voigt fit: the uncertainty of the turbulent broadening (`dbtur`) was never updated and kept its value from before the fit. It is now set like the other errors: computed when `btur` is free, 0 when it is frozen (the default), unchanged for the components outside the fitted group. **The `dbtur` column of the fitted components changes** (usually from empty to 0); no other result changes.

### Changed
- Scripts that use Astrocook show only warnings and errors unless they configure logging (`logging.basicConfig(level=logging.INFO)` shows the progress messages too, see the new tutorial *Using Astrocook from Python*): the package no longer configures logging when imported, as is usual for libraries. The GUI is unchanged.
- Using Astrocook from Python (scripts, notebooks) no longer imports the GUI toolkit (PySide6): the core requested its debug plots through the GUI.
- V1 `.acs` archives are read by a new reader that does not use the V1 code; apart from the fixes above, the result is the same as before.
- Rebin (Flux > Rebin): the option "Return normalized spectrum, if continuum exists" (`norm`) is no longer offered: it was never applied in V2, the spectrum was not normalized. Scripts that pass `norm='True'` still run, with a warning; the normalized flux is shown by the normalized view of the plot, or computed with `apply_expression` (`y / cont`).

### Removed
- The V1 code (`astrocook/legacy/`, ~22,500 lines, with the old wxPython interface): Astrocook no longer depends on it. Dependencies no longer installed: `lmfit`, `tqdm` and `cycler` (the latter still comes with matplotlib); the `[legacy]` extra (wxPython) and the broken V1 Docker image (`containers/`) are gone too.
- Data files used only by the V1 code (~140 MB): `telluric.fits`, `sky_telluric.fits` (identical to each other), `telluric.dat`, `defaults.json`, `espr_spec_form.dat`, `table2_inoue.dat`. The telluric data are kept for a future telluric correction as assets of the GitHub pre-release `data-telluric-2026`.
- Saving sessions in the V1 `.acs` format: sessions are saved as `.acs2` only (`.acs2` is added to a file name without it, and an `.acs` name is refused with a message). V1 `.acs` archives can still be opened; Astrocook 1.9.0 can still write them.

## [2.0.3] - 2026-10-19

Packaging release. **No numerical result changes** with respect to 2.0.2.

### Added
- `astrocook [file ...]` command to launch the GUI; `python -m astrocook` now works from any directory (it needed the root of a source checkout).
- Installation with pip directly from GitHub: `pip install git+https://github.com/DAS-OATs/astrocook.git` (see the installation page; Astrocook is not on PyPI yet).
- Optional extras `[docs]` (to build the documentation) and `[build]` (to build the macOS app).

### Fixed
- *Security*: opening a crafted `.acs2` or `.acs` archive could write files outside the temporary folder; such archives are now refused.
- Installed with pip, the GUI started without its window icon and splash screen: the images now ship with the package.
- The macOS build (`Astrocook.spec`) worked only on one machine (absolute paths); `build_app.sh` ran PyInstaller twice.
- The online documentation is now published automatically after each change.
- The documentation and build scripts stop at once with the command to run when Sphinx or PyInstaller is missing.

### Removed
- Unused dependencies `statsmodels` and `specutils`; the minimum NumPy version is now 1.26.
- Backup data files that were shipped in every package (~54 MB), the old V1 launchers (`legacy_home/`), unused identification code, and the Gemfile of the V1 documentation (the source of all the Dependabot security alerts).

## [2.0.2] - 2026-10-16

Corrective release. **Some numerical results change**: see "Changed behaviour" below.

### Fixed
- *Data loss*: Bayesian fit (MCMC) with 200 steps or fewer overwrote the component parameters with NaN; the burn-in now adapts to the number of steps.
- *Data loss*: the `.acs2` archive did not keep the session log; it is now saved and restored (also in the GUI).
- Voigt fits of broad lines (b > 50 km/s) were biased (e.g. b = 80 fitted as ~50, logN 0.1 dex low): the soft penalty is removed (see below).
- Rebinned and co-added errors did not decrease with the number of pixels per bin (#490).
- Rebin clipping could empty bins and removed the cores of narrow absorption lines; it now clips spikes only (upward, against a running median) and never empties a bin.
- The Voigt model was zero on pixels whose resolution column held non-round or NaN values (e.g. co-added spectra).
- The automatic pipelines (`doublets_auto`, `lya_auto`, ...) gave different results at every run; they are now reproducible.
- Components could not be deleted from a spectrum without a continuum.
- `group_depth` was ignored in the fits ("Fit component", Bayesian fit).
- Rebin with a non-numeric parameter raised an error instead of aborting cleanly.
- The identification viewer failed on sessions reloaded from `.acs2`.
- NaN uncertainties (ill-conditioned fits) are now reported in the log with the components involved.
- ASCII export (#476, #481); warning when saving a session with a V2 log to the V1 `.acs` format.

### Added
- *Edit > Import Systems from ASCII*: read absorption components from a CSV or whitespace-separated table (columns `series`, `z`, optionally `logN`, `b`, ...) (#480).
- Stitching overlapping regions switches the plot to scatter style (#442).
- `max_b` parameter of the automatic pipelines (default 50 km/s) and of *Refit All* (default 300 km/s).

### Changed behaviour
- **Uncertainties**: the errors of the fitted parameters are scaled by sqrt(reduced chi2), as in Astrocook v1; the reduced chi2 is computed on the pixels that constrain the fit and, after the automatic pipelines, on the final model. Errors become larger where the model does not describe the data well (e.g. unresolved blends).
- **Broad lines**: the automatic pipelines bound b to 50 km/s (as the former penalty did, but without biasing the fits); manual fits use only the physical limit (300 km/s). Components that end at 50 km/s may be broader and can be refitted by hand.
- **Continuum smoothing**: the width is now given as a FWHM (`smooth_fwhm`, default 800 km/s, which reproduces the v1 continuum). The former `smooth_std` (a standard deviation, default 500 km/s) is still accepted. The default continuum is smoothed less than in 2.0.1, so continua and the results of the automatic pipelines change slightly.

## [2.0.1] - 2026-04-21

- *Feature*: Added "Equalize and Stitch" to the GUI with dynamic parameters.
- *Fix*: Prevented FITS save crashes from Ly-α identifications and ensured tooltips render safely upon reloading.

## [2.0.0] - 2026-04-07

- Added an experimental Bayesian fitting engine supporting MCMC (`emcee`).
- Introduced interactive Corner Plots to visualize Bayesian posteriors with asymmetric uncertainties and auto-refresh capabilities.
- Fixed a critical group definition bug to ensure physically linked components are always correctly grouped during fitting.
- Improved continuum fitting and smoothing methods, including dynamic re-normalization of the Voigt model.
- Polished the GUI with publication-ready plot formatting, stacked asymmetric error labels, and scientific notation for redshift uncertainties.
- Finalized the transition to the immutable V2 data architecture and removed legacy dead code.

## [2.0.0-beta.3] - 2026-02-18

- Improved robustness (crash recovery) and fitting (weak lines, Lyman forest)
- Updated line lists (`FeII`, `SiII`, etc.) and set turbulent broadening (`btur`) to default to frozen.
- Added batch actions (refit, freeze, delete), group highlighting, and visual feedback for frozen parameters
- Introduced a Data Inspector and a recipe to equalize flux across sessions
- Added Drag-and-Drop support to the main window

## [2.0.0-beta.2] - 2026-01-22

- More flexible way to manually edit continuum
- More robust fitting engine (when parameters hit boundary values)
- More options in adding/fitting components (single doublet component; stop button for refitting)
- Optional refitting of systems upon import into a session

## [2.0.0-beta.1] - 2026-01-02

- Fixed a typo in the `Continuum` menu.

## [2.0.0-beta] - 2025-12-23

Includes

- a fully redesigned architecture, with immutable structures and exposed APIs;
- a new PySide6 QUI
- countless improvements to the V1 algorithms and new algorithms (especially for continuum and absorber fitting)

Also downloadable as a macOS App.