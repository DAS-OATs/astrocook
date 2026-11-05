# Using Astrocook from Python

Everything you do in the GUI can also be done from your own Python
scripts, for instance to process many spectra in the same way. This
tutorial shows the basics: opening a session, running recipes, saving the
result and seeing Astrocook's messages.

## 1. Opening a Session

Open a spectrum or a session file (`.acs2`, V1 `.acs`, or the spectrum
formats the GUI can open) with `load_session_from_file`:
```python
from astrocook.core.session import load_session_from_file

session = load_session_from_file("my_session.acs2", "my_session", None, "auto")
```
The arguments are the file path, the name of the session, the GUI (`None`
in a script) and the file format (`"auto"` detects it, as the GUI does).
The function returns `0` if the file cannot be opened.

The spectrum and the system list are then available as `session.spec` and
`session.systs`, e.g. `session.spec.x` (wavelengths, with units) or
`len(session.systs.components)`.

## 2. Running Recipes

The recipes of the menus are grouped as in the GUI: `session.edit`,
`session.flux`, `session.continuum`, `session.absorbers` and
`session.file`. Each recipe returns a **new** session and leaves the
original unchanged, so the steps can be chained:
```python
smoothed = session.flux.smooth(sigma_kms=50.0)
```
The parameters are the same as in the recipe dialogs. Give them with their
types (numbers, `True` and `False`, `None` where a parameter can be left
out, e.g. the start of *Rebin*), or as strings, as typed in the dialogs and
written in the session log (`sigma_kms="50.0"`): they are converted before
the recipe runs. A value that cannot be converted, or that has no physical
meaning (a negative smoothing width, a redshift of -1 or below, a column
density outside 10-23), stops the recipe with an exception that says in
which recipe it happened and why; so does a recipe that cannot do what was
asked (a session that is not open, a file that cannot be read). Catch
`RecipeError` to go on after it:
```python
from astrocook.recipes.errors import RecipeError

try:
    smoothed = session.flux.smooth(sigma_kms="fifty")
except RecipeError as e:
    print(e)  # Smooth spectrum: sigma_kms: 'fifty' is not a number.
```
The error of a step of a longer recipe names both, e.g. *Populate
components from identification → Auto-estimate continuum (single-click):
...*. Some recipes of the *Absorbers* menu still report a failure by
logging an error and returning `0` instead of a session, until the next
release: check `result == 0` for them.
The recipes and their parameters are listed in the
[API reference](../api/astrocook.recipes).

## 3. Saving the Result

Save a session as an `.acs2` archive; passing the original session as
`initial_session` stores the original data in the archive too, as the GUI
does:
```python
smoothed.save("my_session_smoothed.acs2", initial_session=session)
```

## 4. Seeing Astrocook's Messages

Astrocook reports its progress with Python's `logging` module and, like
other libraries, leaves its configuration to your script: by default only
warnings and errors are shown. To see the progress messages as well, add
these lines at the start of the script:
```python
import logging
logging.basicConfig(level=logging.INFO)
```
The GUI shows them anyway.

Each module of Astrocook has its own logger, named after the module
(`astrocook.recipes.absorbers`, `astrocook.fitting.voigt_fitter`, ...), so
the level can be chosen for a part of the package alone. For example, to
see the progress of Astrocook without that of the other libraries, or to
hide the messages of the fits:
```python
logging.basicConfig()                                    # warnings and errors
logging.getLogger("astrocook").setLevel(logging.INFO)     # all of Astrocook
logging.getLogger("astrocook.fitting").setLevel(logging.WARNING)
```
