"""Spectral templates shipped with the package (``astrocook/data``)."""
import functools
from pathlib import Path
from typing import Tuple

import numpy as np

DATA_DIR = Path(__file__).resolve().parent.parent / 'data'


@functools.cache
def qso_composite() -> Tuple[np.ndarray, np.ndarray]:
    """
    QSO composite spectrum, used to trace the Lyman-alpha emission when
    fitting the continuum (``SpectrumV2.fit_continuum(template=True)``).

    Read once from ``qso_composite.dat`` (the same file as V1).

    Returns
    -------
    x : numpy.ndarray
        Rest-frame wavelength, in Angstrom.
    y : numpy.ndarray
        Flux of the composite, in arbitrary units.

    Both arrays are read-only, since they are shared by all the callers.
    """
    x, y = np.loadtxt(DATA_DIR / 'qso_composite.dat', unpack=True)
    x.flags.writeable = False
    y.flags.writeable = False
    return x, y
