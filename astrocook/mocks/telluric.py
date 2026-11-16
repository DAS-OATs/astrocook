"""Telluric absorption for mocks (phase 7, PR 7a.1).

The template is the molecular absorption (``trans_ma``) of an ESO SkyCalc
sky model (SkyCalc 2.0.9; Noll et al. 2012, A&A 543, A92; Jones et al.
2013, A&A 560, A91), for Paranal at airmass 1.1 and PWV 2.5 mm, without a
line spread function: the file ``sky_telluric.fits`` of the pre-release
``data-telluric-2026`` of the repository. ``data/telluric_paranal.npz``
holds its optical depth, tau = -ln(trans_ma), between 550 and 1000 nm on
the native grid of the model (0.005 nm, vacuum wavelengths), in float32,
with tau < 1e-5 set to zero; it is made by
``benchmarks/make_telluric_template.py``, and its metadata say where it
comes from.

The optical depth is split by wavelength into an O2 part, the gamma, B
and A bands (:data:`O2_BANDS_NM`), and an H2O part, the rest. This is an
approximation: the H2O lines within the O2 bands are scaled as O2, and
the weak O2 lines outside them as H2O; a template per molecule needs one
SkyCalc run each, which the cloud sessions cannot make (the ESO hosts are
not reachable; Guido, 16 November 2026; templates per molecule are a
choice of PR 7b.1). The optical depth is scaled linearly with the airmass
and the PWV, which holds for the line optical depth at the resolution of
the template (the lines are resolved), not for the transmission after the
instrumental convolution. Outside the template (below 550 nm, above 1000
nm) there is no telluric absorption.
"""
import functools
import json
from importlib import resources
from typing import Tuple

import numpy as np

from astrocook.mocks.physics import C_KMS

#: The O2 bands (nm, vacuum): gamma, B and A
O2_BANDS_NM = ((627.0, 634.0), (686.0, 697.0), (759.0, 772.0))
#: Airmass and PWV (mm) of the template
TEMPLATE_AIRMASS = 1.1
TEMPLATE_PWV_MM = 2.5


@functools.lru_cache(maxsize=1)
def template() -> Tuple[np.ndarray, np.ndarray, dict]:
    """Wavelengths (nm), optical depth and metadata of the template."""
    path = resources.files("astrocook.mocks").joinpath("data/telluric_paranal.npz")
    with resources.as_file(path) as p, np.load(p, allow_pickle=False) as data:
        tau = data["tau"].astype(float)
        lam = float(data["lam_start_nm"]) + float(data["lam_step_nm"]) * np.arange(len(tau))
        meta = json.loads(str(data["meta"]))
    return lam, tau, meta


def is_o2(lam_nm: np.ndarray) -> np.ndarray:
    """True within the O2 bands of :data:`O2_BANDS_NM`."""
    lam = np.asarray(lam_nm, dtype=float)
    out = np.zeros(lam.shape, dtype=bool)
    for lo, hi in O2_BANDS_NM:
        out |= (lam >= lo) & (lam <= hi)
    return out


def telluric_optical_depth(lam_nm: np.ndarray, model) -> np.ndarray:
    """Telluric optical depth of a :class:`~astrocook.mocks.TelluricModel`
    at observed wavelengths (nm): the template shifted by
    ``model.shift_kms``, its O2 part scaled by the airmass and its H2O part
    by the airmass and the PWV."""
    lam_t, tau_t, _ = template()
    lam_rest = np.asarray(lam_nm, dtype=float) / (1.0 + model.shift_kms / C_KMS)
    tau = np.interp(lam_rest, lam_t, tau_t, left=0.0, right=0.0)
    airmass = model.airmass / TEMPLATE_AIRMASS
    scale = np.where(is_o2(lam_rest), airmass,
                     airmass * model.pwv_mm / TEMPLATE_PWV_MM)
    return tau * scale
