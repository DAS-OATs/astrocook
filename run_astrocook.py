"""Launch the Astrocook GUI from a source checkout.

The launcher lives in ``astrocook/_run_astrocook.py``; installed packages
provide the ``astrocook`` command and ``python -m astrocook``. This script is
also the PyInstaller entry point (``Astrocook.spec``).
"""
from astrocook._run_astrocook import main

if __name__ == '__main__':
    main()
