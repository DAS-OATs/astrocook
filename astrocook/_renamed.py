"""Old names of renamed classes, served by a module ``__getattr__`` (PEP 562).

In phase 6-quater (PR 6) the classes lost the ``V2`` suffix of their names
(``SpectrumV2`` became ``Spectrum``). A module that defined or exported one
declares the old names it must still answer to, and installs this helper::

    __getattr__ = renamed(__name__, {"SpectrumV2": "Spectrum"})

``from astrocook.core.spectrum import SpectrumV2`` still works: it gives the
same class as the new name (so ``isinstance`` is unchanged) and emits a
``DeprecationWarning``. The old names are used by scripts only: the archives
and the session logs hold no class names. No version of removal is fixed
(3.0 at the earliest). Names that are not in the mapping raise the usual
``AttributeError``.
"""
import sys
import warnings
from typing import Any, Callable, Mapping


def renamed(module_name: str, mapping: Mapping[str, str]) -> Callable[[str], Any]:
    """Build the ``__getattr__`` of the module ``module_name``.

    ``mapping`` goes from each old name to the new name, an attribute of
    the same module (looked up when the old name is asked, so that the
    module may be still importing when the helper is called).
    """
    def __getattr__(name: str) -> Any:
        new = mapping.get(name)
        if new is None:
            raise AttributeError(f"module {module_name!r} has no attribute {name!r}")
        warnings.warn(
            f"{module_name}.{name} is deprecated: the class is now called "
            f"{new} (the old name will be removed in a future major release)",
            DeprecationWarning, stacklevel=2)
        return getattr(sys.modules[module_name], new)
    return __getattr__
