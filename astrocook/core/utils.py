from pathlib import Path
import sys


def resource_path(relative_path):
    if getattr(sys, 'frozen', False):
        base_path = Path(sys._MEIPASS)
    else:
        # The package directory (astrocook/assets/ is installed with it;
        # Astrocook.spec copies it to 'assets' in the frozen app)
        base_path = Path(__file__).resolve().parent.parent
    return str(base_path / relative_path)
