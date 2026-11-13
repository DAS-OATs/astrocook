# Using Astrocook from Python

Everything you do in the GUI can also be done from your own Python
scripts, for instance to process many spectra in the same way. This
tutorial shows the basics: opening a session, running recipes, saving the
result and seeing Astrocook's messages.

## 1. Opening a Session

Open a spectrum or a session file (`.acs2`, V1 `.acs`, or the spectrum
formats the GUI can open) with `load_session_from_file`:
```python
from astrocook.session.session import load_session_from_file

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
The error of a step of a longer recipe names both, e.g. *Auto-fit metal
doublets (Red Side) → Auto-detect and fit components → Refit all systems:
...*. A recipe with nothing to do (e.g. *Refit All* on a session without
components) is not a failure: it gives the session back and logs a
warning.

To know what a recipe did, run it with `run`, which returns its whole
result: the new session, what changed (the components added, removed and
changed, the columns, the properties), the warnings, the recipes it ran as
steps and the values of the parameters it used:
```python
from astrocook.recipes.result import run

result = run(session, "absorbers.doublets_auto")
print(result.summary)   # Auto-fit metal doublets (Red Side): 3 components added; ...
for c in result.changes.components_added:
    print(c.series, c.z, c.logN, c.b)
print(result.warnings, result.steps)
fitted = result.session
```
`result.notes` says what else is worth knowing, one line each: the
components removed and why (e.g. *1 component removed: logN < min_logN =
14.0 (too weak)*), those the recipe made and discarded, those whose fit
failed (flagged `not_fitted`, with the values from before), the candidate
systems found; `result.extra` has the same as data (by component). The GUI
counts them in the status bar and writes them in the log.
The recipes and their parameters are listed in the
[API reference](../api/astrocook.recipes).

## 3. What a Recipe Needs

Some recipes need something from the session: the emission redshift
(`z_em`), the resolution, a continuum. The [API reference](../api/astrocook.recipes)
lists them for each recipe (`requires`), and
`astrocook.recipes.prerequisites.missing` says what a session lacks for a
recipe, and which recipe provides it:
```python
from astrocook.recipes.prerequisites import missing
from astrocook.recipes.registry import find_recipe

for gap in missing(session, find_recipe("doublets_auto")):
    print(gap.message)
# The spectrum has no continuum: Auto-estimate continuum (single-click)
# (continuum.estimate_auto) provides it with its defaults
```
By default a script goes on: a continuum that is missing is estimated
first with *Auto-estimate Continuum* and its defaults, with a warning; a
recipe that needs `z_em` or the resolution and does not find them goes on
too, with a warning that says what it does without them. To stop instead,
with a `MissingPrerequisite` exception that says what is missing and how to
provide it, choose the policy `"raise"`, for a block of the script or for
one call:
```python
from astrocook.recipes.prerequisites import policy

with policy("raise"):
    fitted = session.absorbers.doublets_auto()

result = run(session, "absorbers.doublets_auto", on_missing="raise")
```
A recipe run first to provide something is in `result.provided`, with the
session it gave, and the GUI writes it in the session log
as a step of its own, before the recipe.

## 4. Looking at a Session

A few functions read a session without changing it, as the System
Inspector and the main plot of the GUI do (the AI agents use them too):
```python
from astrocook.plotting import plot_spectrum

# The flux, error, continuum and model in a range (nm), averaged in
# blocks of pixels if there are more than 200
data = fitted.spec.columns_in_range(["y", "dy", "cont"], 449, 451, max_points=200)
data["columns"]["y"], data["binned_by"]

# The groups of components that Refit All fits together
uuids = [c.uuid for c in fitted.systs.components]
groups = fitted.systs.fit_groups(uuids, max_depth=2, spec=fitted.spec)

# A plot of a range, normalized to the continuum, with a tick at each line
# of each component; in a notebook, the figure shows itself
fig = plot_spectrum(fitted, 445, 452, normalized=True)
fig.savefig("doublet.png")
```

The log of a session as the text of a methods section, as *File > Export
Log as Text* writes it (one paragraph per recipe, with the reasons, what it
changed and the notes, then the references):
```python
from astrocook.session.history import SessionHistory
from astrocook.methods_text import methods_text

text = methods_text(history)        # a SessionHistory: the states and the log
print(text)                          # or text.paragraphs, text.references
methods_text(history, menus=False)  # without the menus, for a paper
```
A session alone (`methods_text(session)`) gives its log, without what each
recipe changed: the states before and after it are kept by the history (the
GUI and the MCP server keep one per session).

## 5. Saving the Result

Save a session as an `.acs2` archive; passing the original session as
`initial_session` stores the original data in the archive too, as the GUI
does:
```python
smoothed.save("my_session_smoothed.acs2", initial_session=session)
```

## 6. Seeing Astrocook's Messages

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
