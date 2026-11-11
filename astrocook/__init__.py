# astrocook/__init__.py
from astrocook._renamed import renamed

# Expose the Session as the main entry point
from .session.session import Session

__all__ = ["Session"]

# Expose formatting/version info if needed
__version__ = "2.6.0"

# Old name (phase 6-quater, PR 6): see astrocook/_renamed.py
__getattr__ = renamed(__name__, {"SessionV2": "Session"})
