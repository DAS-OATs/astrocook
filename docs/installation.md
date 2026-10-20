# Installation

Astrocook is available as a standalone app for macOS or as a Python package for advanced users.

## Quick Start ( macOS)

If you are on a Mac and want to use the GUI immediately without managing Python environments:

1. **Download the App**: Go to our **[Latest Release Page](https://github.com/DAS-OATs/astrocook/releases/latest)** and look for `Astrocook.dmg`.
2. **Install**: Open `Astrocook.dmg` and drag the app to your `Applications` folder.
3. **Run**: Double-click **Astrocook** in your Applications.

:::{note}
If you see a security warning, you may need to Right-Click the app and select 'Open' the first time, or explicitly authorize the app on System Settings).
:::

## Installing with pip (Linux, Windows, macOS)

Astrocook requires **Python 3.12 or later**. `pip` installs all the other
dependencies (PySide6, astropy, numpy, scipy, matplotlib, ...).

Astrocook is not on PyPI yet; install it directly from GitHub, preferably in
a virtual environment:
```
python3 -m venv ac2
source ac2/bin/activate  # On Windows: ac2\Scripts\activate
pip install git+https://github.com/DAS-OATs/astrocook.git
```

The Bayesian fit (MCMC and nested sampling, corner plots) needs optional
dependencies (emcee, dynesty, corner), which the macOS app already
includes. To install them with the package:
```
pip install "astrocook[bayes] @ git+https://github.com/DAS-OATs/astrocook.git"
```
Without them, the Bayesian entries of the system inspector are disabled and
the `fit_bayesian` recipe raises an error naming this command.

:::{note}
If the `venv` command fails, you may need to install the module first (e.g., `sudo apt install python3-venv` on Ubuntu/Debian). Check the Python version with `python3 --version`; if it is older than 3.12, use a newer interpreter (e.g., `python3.12 -m venv ac2`).
:::

To update to the latest version, run the same `pip install` command with
`--upgrade`.

## Running the App

With the virtual environment active, launch the GUI from any directory:
```
astrocook [file ...]
```
`python -m astrocook [file ...]` is equivalent. The optional arguments are
spectra or sessions to open at startup.

:::{note}
If the app fails to launch with a Qt/xcb error, ensure your system has the necessary Qt libraries installed (e.g., `sudo apt install libxcb-cursor0`).
:::

## Using Astrocook from Python

The package can also be used from your own scripts, without the GUI. For
example, to open a session and smooth its spectrum:
```python
from astrocook.core.session import load_session_from_file

session = load_session_from_file("my_session.acs2", "my_session", None, "auto")
smoothed = session.flux.smooth(sigma_kms="50.0")
```

Astrocook reports its progress with Python's `logging` module and, like
other libraries, leaves its configuration to your script: by default only
warnings and errors are shown. To see the progress messages as well, add
these lines at the start of the script:
```python
import logging
logging.basicConfig(level=logging.INFO)
```
The GUI shows them anyway.

## Developer Installation

To modify the code, clone the repository and install it in editable mode,
with the test and development tools:
```
git clone https://github.com/DAS-OATs/astrocook.git
cd astrocook
python3 -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
pip install -e ".[dev]"
```

Code changes take effect without reinstalling. After a `git pull` that
changes `pyproject.toml` (dependencies, commands), run
`pip install -e ".[dev]"` again: commands such as `astrocook` are created at
installation time. The tests run with `pytest`.

Optional extras for maintenance tasks: `pip install -e ".[docs]"` to build
this documentation (`generate_html.sh`, `deploy_docs.sh`) and
`pip install -e ".[build]"` to build the macOS app (`build_app.sh`; the DMG
also needs `brew install create-dmg`).
