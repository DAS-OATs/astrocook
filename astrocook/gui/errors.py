"""Exceptions that the GUI expects and handles (phase 6).

The GUI catches only the errors it expects; an unexpected one is a bug, and
reaches the boundary of the action that raised it (``ROADMAP.md``, phase 6,
inventory of the broad except clauses).
"""
import astropy.units as au

#: Errors of the data being drawn or shown (a value missing, of the wrong
#: type or unit, an index out of range): the GUI shows what it can and goes
#: on. ``ValueError`` covers ``np.linalg.LinAlgError``. ``AttributeError`` and
#: ``NameError`` are left out: they mean a bug.
DATA_ERRORS = (ValueError, TypeError, KeyError, IndexError, ArithmeticError,
               au.UnitsError)

#: Errors of Qt and matplotlib objects that are already gone or not ready
#: (an artist removed twice, a connection that does not exist, a dialog or a
#: widget whose C++ object was deleted, axes not drawn yet): harmless when
#: removing, disconnecting, closing or blitting them.
QT_ERRORS = (ValueError, TypeError, NotImplementedError, AttributeError,
             RuntimeError)
