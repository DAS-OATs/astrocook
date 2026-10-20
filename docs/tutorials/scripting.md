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
smoothed = session.flux.smooth(sigma_kms="50.0")
```
The parameters are the same as in the recipe dialogs, and are given as
strings, as typed in the dialogs. If a parameter cannot be understood, the
recipe logs an error and returns `0` instead of a session, so check the
result before going on:
```python
result = session.flux.smooth(sigma_kms="fifty")
if result == 0:
    print("The recipe failed")
```
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
