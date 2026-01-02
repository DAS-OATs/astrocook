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

## Developer Installation (Python)

For Linux/Windows users, or developers who want to modify the code, install from source.

### Prerequisites

Astrocook requires Python 3.9 or later and the following libraries:
- **GUI**: `PySide6`
- **Science**: `astropy`, `numpy`, `scipy`
- **Plotting**: `matplotlib`, `scienceplots`
- **Utils**: `qtawesome`, `numexpr`

### Installing from Source

1. Clone the repository:
```
git clone https://github.com/das-oats/astrocook.git
cd astrocook
```

2. Move to the `develop` branch (where V2 is temporarily hosted):
```
git checkout develop
git pull
```

3. Create a virtual environment (Recommended):
```
python -m venv ac2
source ac2/bin/activate  # On Windows: ac2\Scripts\activate
```

:::{note}
If the `venv` command fails, you may need to install the module first (e.g., `sudo apt install python3-venv` on Ubuntu/Debian).
:::

4. Install dependencies:
```
pip install -e .
```

### Running the App

To launch the GUI from the source code, run this from the root directory:
```
python -m astrocook
```

:::{note}
If the app fails to launch with a Qt/xcb error, ensure your system has the necessary Qt libraries installed (e.g., `sudo apt install libxcb-cursor0`).
:::