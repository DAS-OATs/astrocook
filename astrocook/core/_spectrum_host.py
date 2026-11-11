"""What the mixins of :class:`~astrocook.core.spectrum.Spectrum` use of
the class that mixes them in, declared for the type checker (phase
6-quater, PR 4).

Each mixin (``spectrum_flux.py`` and the others) uses attributes and methods
that ``Spectrum`` defines, or another mixin does. The mixins inherit
``SpectrumHost`` only when the type checker reads them; at run time their
base is ``object``, so this module changes nothing: see
``HostBase`` below. The declarations here are of types only: where one
changes, change the real definition with it.
"""
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Self

import astropy.units as au
import numpy as np

from astrocook.core.structures import SpectrumData

if TYPE_CHECKING:
    from astrocook.core.spectrum import TableAdapter

    class SpectrumHost:
        """The part of ``Spectrum`` that the mixins use."""
        _data: SpectrumData
        history: List[Any]

        def __init__(self, data: SpectrumData,
                     history: Optional[list] = None) -> None: ...

        @property
        def x(self) -> au.Quantity: ...
        @property
        def y(self) -> au.Quantity: ...
        @property
        def dy(self) -> au.Quantity: ...
        @property
        def meta(self) -> Dict[str, Any]: ...
        @property
        def t(self) -> 'TableAdapter': ...

        @property
        def xmin(self) -> au.Quantity: ...
        @property
        def xmax(self) -> au.Quantity: ...
        def resolution(self) -> np.ndarray: ...
        @property
        def cont(self) -> Optional[au.Quantity]: ...

        def get_column(self, name: str) -> Optional[au.Quantity]: ...
        def has_aux_column(self, name: str) -> bool: ...
        def update_column(self, name: str, values: np.ndarray,
                          unit: Optional[au.Unit] = None,
                          description: str = "") -> Self: ...
        # Defined in the mixins, used by another one
        def scale_flux(self, scale_factor: Any) -> Self: ...
        def find_absorbed(self, *args: Any, **kwargs: Any) -> Self: ...
        def _renormalize_model(self, *args: Any, **kwargs: Any) -> Any: ...

    HostBase = SpectrumHost
else:
    HostBase = object
