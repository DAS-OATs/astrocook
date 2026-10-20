"""Numerical snapshots of the continuum estimation on a real spectrum.

Input: raw J0103-1305 spectrum in a window covering the Lyman-alpha forest,
the Lyman-alpha emission and the red side (see ``conftest.py``). Parameters
are the defaults of the continuum recipes. Continua are stored every
``STEP`` pixels: they are smooth, and this keeps the snapshot files small.
"""
import numpy as np
import astropy.units as au
import pytest

from astrocook.core.spectrum_operations import (
    compute_auto_fudge, find_absorbed_regions, interpolate_continuum_mask,
)

STEP = 10

# Defaults of recipes/continuum.py
SMOOTH_LEN_LYA = 5000.0 * au.km / au.s
SMOOTH_LEN_OUT = 400.0 * au.km / au.s
KAPPA = 2.0
SMOOTH_STD_KMS = 500.0
POWERLAW_REGIONS = "128.0-129.0, 131.5-132.5, 134.5-136.0"


@pytest.fixture(scope="module")
def absorbed(j0103_window):
    return j0103_window.find_absorbed(SMOOTH_LEN_LYA, SMOOTH_LEN_OUT, KAPPA, False)


def test_find_absorbed(snapshot, absorbed):
    mask = absorbed.get_column("abs_mask").value.astype(bool)
    snapshot.check("continuum_j0103_abs_mask", {"abs_mask": mask})


def test_find_absorbed_regions_matches_method(j0103_window, absorbed):
    # The pure function and the SpectrumV2 method must agree
    spec = j0103_window
    mask = find_absorbed_regions(spec.x, spec.y.value, spec.dy.value,
                                 spec._data.z_em, SMOOTH_LEN_LYA,
                                 SMOOTH_LEN_OUT, KAPPA, False)
    np.testing.assert_array_equal(
        mask, absorbed.get_column("abs_mask").value.astype(bool))


def test_auto_fudge(snapshot, absorbed):
    mask = absorbed.get_column("abs_mask").value.astype(bool)
    y = absorbed.y.value
    cont_raw = interpolate_continuum_mask(absorbed.x.value, y, mask)
    fudge = compute_auto_fudge(y, mask, cont_raw)
    snapshot.check("continuum_j0103_auto_fudge",
                   {"cont_raw": cont_raw[::STEP], "fudge": np.float64(fudge)})


@pytest.mark.parametrize("fudge", ["auto", 1.0])
def test_fit_continuum(snapshot, absorbed, fudge):
    cont = absorbed.fit_continuum(fudge, SMOOTH_STD_KMS).get_column("cont")
    snapshot.check(f"continuum_j0103_fit_fudge_{fudge}",
                   {"cont": cont[::STEP]})


def test_fit_powerlaw(snapshot, j0103_window):
    cont_pl = j0103_window.fit_powerlaw(POWERLAW_REGIONS, 3.0).get_column("cont_pl")
    snapshot.check("continuum_j0103_powerlaw", {"cont_pl": cont_pl[::STEP]})


def test_estimate_auto_smoothing_is_a_fwhm(j0103_window):
    # F16 (findings register, decided with Guido): the final smoothing is
    # given as a FWHM (default 800 km/s); the deprecated smooth_std is a sigma
    from astrocook.core.session import SessionV2
    session = SessionV2("w", None, spec=j0103_window)
    default = session.continuum.estimate_auto().spec.get_column("cont").value
    as_sigma = session.continuum.estimate_auto(
        smooth_std=str(800.0 / 2.354820045)).spec.get_column("cont").value
    np.testing.assert_allclose(default, as_sigma, rtol=1e-12)
    wider = session.continuum.estimate_auto(smooth_fwhm="1600").spec.get_column("cont").value
    assert not np.allclose(default, wider)


def test_fit_continuum_template(snapshot, absorbed):
    # template=True divides the flux by the QSO composite (data file
    # qso_composite.dat, Angstrom, shifted to z_em) before the interpolation
    # and multiplies the continuum back: the window covers the Lyman-alpha
    # emission, where the template matters
    cont = absorbed.fit_continuum(1.0, SMOOTH_STD_KMS, template=True).get_column("cont")
    plain = absorbed.fit_continuum(1.0, SMOOTH_STD_KMS).get_column("cont")
    assert not np.allclose(cont.value, plain.value)
    snapshot.check("continuum_j0103_fit_template", {"cont": cont[::STEP]})
