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