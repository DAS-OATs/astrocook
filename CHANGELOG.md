## [Unreleased]

### Fixed
- *Data loss*: a session that could not be read (an `.acs2` archive without its metadata file, a damaged `.acs` archive, a file in a format without a loader) was opened with mock data in place of the spectrum, and saving it overwrote the file. An `.acs2` archive without metadata now opens with its spectrum and systems (without constraints, session log and identifications); the other cases report an error.
- V1 `.acs` archives: the turbulent broadening (`btur`) and its error, and the reduced chi2 of each component are now read (they were set to 0 or dropped); the flux keeps its unit, erg/(cm2 s Angstrom), so that unit conversions (Edit > Convert y axis) work on V1 sessions.

- V1 `.acs` archives: the links between parameters (e.g. SiIV z and b tied to CIV) are now kept as V2 constraints. Before, the linked parameters were frozen at their V1 values, did not follow their target in a refit, and the System Inspector showed no target; constraints on logN were lost. **Refitting a V1 session with links gives different (correct) results.**

### Changed
- Scripts that use Astrocook show only warnings and errors unless they configure logging (`logging.basicConfig(level=logging.INFO)` shows the progress messages too, see the new tutorial *Using Astrocook from Python*): the package no longer configures logging when imported, as is usual for libraries. The GUI is unchanged.
- V1 `.acs` archives are read by a new reader that does not use the V1 code; apart from the fixes above, the result is the same as before.

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