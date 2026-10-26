# Screenshots of the tutorials

The images of the tutorials (`docs/_static/*.png`, except the logo and the
icon) are made by `make_screenshots.py`, which runs the real GUI offscreen.
Rebuild them after a change of the interface, look at every image, and
commit the ones that changed.

```bash
python docs/screenshots/make_screenshots.py --download   # all scenes, ~4 min
python docs/screenshots/make_screenshots.py --only continuum_auto flux_smooth
python docs/screenshots/make_screenshots.py --list
```

## Data

Two spectra, anonymised as "My QSO" and "My other QSO" (provided by Guido,
October 2026; rebuilt with only the spectrum columns, x to resol, and a
minimal header: no coordinates, dates, archive identifiers or earlier
analysis), are kept out of the repository as assets of the GitHub
pre-release `data-docs-2026`: `My_QSO.acs2` and `My_other_QSO.acs2`.
`--download` fetches the missing ones into `docs/screenshots/data/`
(ignored by git); `--data` reads them from another folder. The session takes
the name of the file, with spaces for the underscores.

They are raw spectra: the scenes set the properties of "My QSO"
(z_em = 2.730, R = 135000), estimate the continuum and run the pipelines
themselves. The scenes of the System Inspector need the whole chain
(continuum, `doublets_auto`, `lya_auto`, about 90 s): it runs once per
call and its result is reused.

## How the images are made

- Qt runs on the `offscreen` platform with a large virtual screen (the
  default one, 400 x 400 points, squeezes the dialogs) and a scale factor of
  2, so the images have the density of a retina screen. The style is Fusion,
  the same on every machine.
- Windows, dialogs, menus and tooltips are grabbed one by one and composed:
  the offscreen platform draws no window frame and no mouse pointer, so the
  script draws a neutral frame (title bar, rounded corners, shadow) and the
  pointer where a scene puts it.
- Menus that the GUI opens with `exec()` (context menus) are grabbed from
  inside their own event loop by a timer, with the chosen entry
  highlighted, and then closed.
- Each scene is a function named after its image, registered with
  `@scene`; positions are in logical pixels (points) of the main window.

## Adding or changing a scene

Write a function decorated with `@scene` that returns a Pillow image, using
the helpers of `Studio` (open sessions, run recipes, set axis limits, open
panels and sections, dialogs, context menus, combo boxes) and the image
helpers (`framed`, `popup`, `paste`, `pointer`, `crop`, `label`). Run it
with `--only`, look at the result, and adjust.

`tests/test_screenshots.py` runs two quick scenes on a test archive, so that
a change of the GUI that breaks the script shows up in its pull request.
