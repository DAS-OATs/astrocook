## [2.0.0-beta.3] - 2026-02-18

- Improved robustness (crash recovery) and fitting (weak lines, Lyman forest)
- Updated line lists (FeII, SiII, etc.) and set turbulent broadening (btur) to default to frozen.
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