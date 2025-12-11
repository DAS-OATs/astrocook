# Installation

Astrocook is available as a standalone app for macOS or as a Python package for advanced users.

## Quick Start ( macOS)

If you are on a Mac and want to use the GUI immediately without managing Python environments:

1. **Download the App**: Go to our **[Latest Release Page](https://github.com/DAS-OATs/astrocook/releases/tag/v2.0.0-alpha)** and look for `Astrocook.dmg`.
2. **Install**: Open `Astrocook.dmg` and drag the app to your `Applications` folder.
3. **Run**: Double-click **Astrocook** in your Applications.

*(Note: If you see a security warning, you may need to Right-Click the app and select 'Open' the first time).*

## Developer Installation (Python)

For Linux/Windows users, or developers who want to modify the code, install from source.

### Prerequisites

Astrocook requires Python 3.9 or later and the following libraries:
- **GUI**: `PySide6`
- **Science**: `astropy`, `numpy`, `scipy`
- **Plotting**: `matplotlib`, `scienceplots`
- **Utils**: `qtawesome`, `numexpr`

### Installing from Source

1. Clone the repository and move to the `develop` branch:
```
git clone https://github.com/das-oats/astrocook.git
git checkout develop
git pull
cd astrocook
```

2. Create a virtual environment (Recommended):
```
python -m venv ac2
source ac2/bin/activate  # On Windows: ac2\Scripts\activate
```

3. Install dependencies:
```
pip install -e .
```

### Running the App

To launch the GUI from the source code, run this from the root directory:
```
python -m astrocook
```