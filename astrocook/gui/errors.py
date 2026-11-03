"""Exceptions that the GUI expects and handles (phase 6).

The GUI catches only the errors it expects; an unexpected one is a bug, and
reaches the boundary of the action that raised it (``ROADMAP.md``, phase 6,
inventory of the broad except clauses); what nothing handles is shown to
the user by the hook of :func:`install_excepthook`.
"""
import logging
import sys
import traceback

import astropy.units as au

logger = logging.getLogger(__name__)

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


# --- Unexpected errors -----------------------------------------------------------------

_showing = False


def report_unexpected_error(exc_type, exc, tb) -> None:
    """Log an unexpected error with its traceback and show it to the user.

    Called by the hook installed by :func:`install_excepthook` for the
    exceptions that nothing handled (PySide6 passes those of the slots and
    of the timers to ``sys.excepthook``): before phase 6 they reached only
    the terminal. The dialog shows the message, with the traceback under
    *Details*; it is shown only in the main thread, and not while another
    one is open (an error in a redraw could otherwise open one per frame).
    """
    logger.error(f"Unexpected error: {exc_type.__name__}: {exc}",
                 exc_info=(exc_type, exc, tb))
    global _showing
    from PySide6.QtCore import QThread
    from PySide6.QtWidgets import QApplication, QMessageBox
    app = QApplication.instance()
    if _showing or app is None or QThread.currentThread() is not app.thread():
        return
    _showing = True
    try:
        box = QMessageBox(QMessageBox.Icon.Critical, "Unexpected Error",
                          f"An unexpected error occurred: {exc_type.__name__}: {exc}",
                          QMessageBox.StandardButton.Ok, app.activeWindow())
        box.setInformativeText("The operation was not completed; Astrocook goes on. "
                               "The details are also in the log.")
        box.setDetailedText("".join(traceback.format_exception(exc_type, exc, tb)))
        box.exec()
    finally:
        _showing = False


def install_excepthook():
    """Install :func:`report_unexpected_error` as ``sys.excepthook``.

    A ``KeyboardInterrupt`` goes to the previous hook. Returns the previous
    hook (the launcher installs it once, at the start of the GUI).
    """
    previous = sys.excepthook

    def hook(exc_type, exc, tb):
        if issubclass(exc_type, KeyboardInterrupt):
            previous(exc_type, exc, tb)
            return
        report_unexpected_error(exc_type, exc, tb)

    sys.excepthook = hook
    return previous
