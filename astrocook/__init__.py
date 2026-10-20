# astrocook/__init__.py

import logging

# INFO messages on stderr for scripts using the API, as when astrocook.legacy
# was imported with the package (legacy/message.py called basicConfig);
# kept in phase 3 until Guido decides whether the package should configure
# logging. It does nothing if logging is already configured.
logging.basicConfig(level=logging.INFO,
                    format="[%(levelname)s] %(module)s: %(message)s")

# Expose the V2 Session as the main entry point
from .core.session import SessionV2

__all__ = ["SessionV2"]

# Expose formatting/version info if needed
__version__ = "2.0.3"