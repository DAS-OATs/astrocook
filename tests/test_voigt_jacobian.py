"""The analytic Jacobian of the Voigt fit and the final model only where it is
used (phase 7, PR 7c.4d-bis, F171 and F167): the derivatives of the model
against finite differences, the fits with the analytic Jacobian against those
with the finite differences of ``least_squares`` (the Jacobian before), and
``fit(model=False)``. Data from the mock generator."""
import dataclasses

import astropy.units as au
import numpy as np
import pytest
from scipy.optimize._numdiff import approx_derivative

from astrocook.core.spectrum import Spectrum
from astrocook.core.structures import DataColumn, SystemListData
from astrocook.core.system_list import SystemList
from astrocook.fitting.lsf import CHUNK, PixelLSF
from astrocook.fitting.voigt_fitter import VoigtFitter, voigt_optical_depth
from astrocook.mocks import MockComponent, MockTruth, generate

LYA = 1215.6701


def _empty():
    return SystemList(SystemListData(components=[]))


def _spec(resol, comps, lam=(455.0, 470.0), seed=1, snr=50.0):
    return generate(MockTruth(components=comps, lam_min_nm=lam[0], lam_max_nm=lam[1],
                              dv_kms=2.5 if (resol or 0) > 20000 else 10.0, resol=resol,
                              snr=snr, z_em=3.9, seed=seed))


def _with_resol(spec, resol):
    """``spec`` with its ``resol`` column replaced (None: no column)."""
    aux = {k: v for k, v in spec._data.aux_cols.items() if k != "resol"}
    if resol is not None:
        aux["resol"] = DataColumn(resol, au.dimensionless_unscaled)
    return Spectrum(dataclasses.replace(spec._data, aux_cols=aux))


def _doublet_and_lya():
    """A CIV doublet and a Ly-a in 455-470 nm, the start a bit off."""
    comps = [MockComponent("CIV", 2.0, 13.5, 12.0), MockComponent("Ly_a", 2.784, 14.0, 25.0)]
    systs = (_empty().add_component("CIV", 2.00001, logN=13.3, b=10.0)
             .add_component("Ly_a", 2.784, logN=13.8, b=20.0))
    return comps, systs


def _numerical(fitter, p, lower, upper):
    """Central differences with small steps (1e-7 in z, 1e-5 in logN and
    b), within the bounds."""
    is_free = fitter._constraints._param_map['is_free']
    kind = (np.flatnonzero(is_free) % 4)[:len(p)]
    step = np.where(kind == 0, 1e-7, 1e-5)
    return approx_derivative(fitter._residual_function, p, method='3-point',
                             abs_step=step, bounds=(lower, upper))


def _check_jacobian(spec, systs, rtol=1e-5):
    fitter = VoigtFitter(spec, systs)
    lower, upper = fitter.prepare_fit_context()
    p = np.clip(fitter._constraints.p_free_vector, lower, upper)
    analytic = fitter._jacobian(p)
    numerical = _numerical(fitter, p, lower, upper)
    assert analytic.shape == numerical.shape
    scale = np.maximum(np.max(np.abs(numerical), axis=0), 1e-12)
    err = np.max(np.abs(analytic - numerical), axis=0) / scale
    assert np.all(err < rtol), err


# --- The convolution as a matrix ---------------------------------------------

@pytest.mark.parametrize("resol", [45000.0, None, "varying"])
@pytest.mark.parametrize("n_pix", [300, CHUNK + 700])
def test_the_matrix_of_the_lsf_is_its_convolution(resol, n_pix):
    edges = 6000.0 * np.exp(np.arange(n_pix + 1) * 2.5 / 299792.458)
    r = np.linspace(40000.0, 50000.0, n_pix) if resol == "varying" else resol
    lsf = PixelLSF(edges[:-1], edges[1:], r_min=40000.0, r_max=50000.0, b_min=8.0)
    rng = np.random.default_rng(0)
    flux = np.exp(-rng.uniform(0.0, 2.0, len(lsf.v_fine)))
    key = "column" if resol == "varying" else None
    np.testing.assert_allclose(lsf.matrix(r, key=key) @ flux, lsf.apply(flux, r, key=key),
                               rtol=0, atol=1e-13)


def test_the_optical_depth_with_its_derivatives_is_that_of_the_model():
    from astrocook.fitting.voigt_fitter import voigt_optical_depth_derivatives
    grid = np.linspace(4590.0, 4610.0, 4000)
    args = (np.array([LYA, 1548.204]), np.array([0.4164, 0.19]), np.array([6.265e8, 2.643e8]),
            np.array([2.78, 1.97]), np.array([1e14, 1e13]), np.array([20.0, 8.0]))
    tau = voigt_optical_depth_derivatives(grid, *args)[0]
    assert np.array_equal(tau, voigt_optical_depth(grid, *args))


# --- The Jacobian ---------------------------------------------------------------

@pytest.mark.parametrize("resol", [45000.0, 8900.0, None])
def test_the_jacobian_is_that_of_the_model(resol):
    comps, systs = _doublet_and_lya()
    _check_jacobian(_spec(resol, comps), systs)


def test_the_jacobian_with_a_resolution_per_pixel():
    comps, systs = _doublet_and_lya()
    spec = _spec(8900.0, comps)
    x = spec.x.to(au.nm).value
    _check_jacobian(_with_resol(spec, 8900.0 * (x / 462.0) ** 4), systs)


def test_the_jacobian_with_two_resolutions():
    # A component with its own resolution is convolved apart: the model is
    # the product of the two groups
    comps, systs = _doublet_and_lya()
    uuid = systs.components[1].uuid
    systs = SystemList(dataclasses.replace(systs._data, components=[
        c if c.uuid != uuid else dataclasses.replace(c, resol=20000.0)
        for c in systs.components]))
    _check_jacobian(_spec(45000.0, comps), systs)


def test_the_jacobian_with_links_and_turbulence():
    comps, systs = _doublet_and_lya()
    a, b = (c.uuid for c in systs.components)
    systs = (systs.add_component("Ly_a", 2.7843, logN=13.0, b=15.0))
    c = systs.components[2].uuid
    systs = (systs.update_constraint(c, "b", is_free=False, expression=f"p['{b}'].b * 0.8",
                                     target_uuid=b)
             .update_constraint(c, "z", is_free=False, expression=f"p['{b}'].z + 3e-4",
                                target_uuid=b)
             .update_constraint(a, "btur", is_free=True))
    systs = SystemList(dataclasses.replace(systs._data, components=[
        dataclasses.replace(x, btur=5.0) if x.uuid == a else x for x in systs.components]))
    _check_jacobian(_spec(45000.0, comps), systs)


def test_the_jacobian_of_a_damped_line():
    # The damping wings: the derivative in a, through b
    comps = [MockComponent("Ly_a", 2.784, 20.3, 20.0)]
    systs = _empty().add_component("Ly_a", 2.7841, logN=20.0, b=25.0)
    _check_jacobian(_spec(45000.0, comps, snr=20.0), systs)


def test_the_jacobian_with_fixed_components_around():
    rng = np.random.default_rng(3)
    comps = [MockComponent("Ly_a", lam / LYA - 1, rng.uniform(12.8, 14.5), rng.uniform(15, 35))
             for lam in np.sort(rng.uniform(4560.0, 4640.0, 20))]
    systs = _empty()
    for c in comps:
        systs = systs.add_component("Ly_a", c.z, logN=c.logN + 0.1, b=c.b * 0.9)
    fit = [c.uuid for c in systs.components[8:12]]
    with systs.fitting_context(fit, group_depth=0):
        _check_jacobian(_spec(45000.0, comps), systs)


# --- The fits -------------------------------------------------------------------

def _fit_both(spec, systs, **kwargs):
    out = []
    for jac in ("analytic", "2-point"):
        fitted, _, res = VoigtFitter(spec, systs).fit(jac=jac, **kwargs)
        out.append((fitted, res))
    return out


@pytest.mark.parametrize("resol", [45000.0, 8900.0])
def test_the_fit_is_that_with_finite_differences(resol):
    comps, systs = _doublet_and_lya()
    (a, ra), (n, rn) = _fit_both(_spec(resol, comps), systs)
    assert ra.success and rn.success
    # The same minimum: chi2, values within 1e-3 of their errors, errors
    # within 1e-3
    assert 2 * ra.cost == pytest.approx(2 * rn.cost, rel=1e-8)
    for ca, cn in zip(a.components, n.components, strict=True):
        for p in ("z", "logN", "b"):
            assert abs(getattr(ca, p) - getattr(cn, p)) < 1e-3 * getattr(cn, "d" + p)
            assert getattr(ca, "d" + p) == pytest.approx(getattr(cn, "d" + p), rel=1e-3)
        assert ca.chi2 == pytest.approx(cn.chi2, rel=1e-6)


@pytest.mark.parametrize("resol", [45000.0, 8900.0])
def test_a_blended_group_is_fitted_as_with_finite_differences(resol):
    # 12 Ly-a lines in 2 nm: blends whose parameters are degenerate; the
    # optimizer stops (ftol) at the same chi2 within 1e-6, the values
    # within a small fraction of their errors (measured with 40 lines: up
    # to 0.1 sigma at R 8900, from the path of the optimizer)
    rng = np.random.default_rng(2)
    comps = [MockComponent("Ly_a", lam / LYA - 1, rng.uniform(12.8, 14.5), rng.uniform(15, 35))
             for lam in np.sort(rng.uniform(4590.0, 4610.0, 12))]
    systs = _empty()
    for c in comps:
        systs = systs.add_component("Ly_a", c.z + rng.normal(0, 1e-5), logN=c.logN + 0.2,
                                    b=c.b * 1.1)
    (a, ra), (n, rn) = _fit_both(_spec(resol, comps, lam=(455.0, 466.0)), systs,
                                 max_nfev=300, b_max=50.0, logN_max=19.0)
    assert 2 * ra.cost == pytest.approx(2 * rn.cost, rel=1e-6)
    for ca, cn in zip(a.components, n.components, strict=True):
        for p in ("z", "logN", "b"):
            assert abs(getattr(ca, p) - getattr(cn, p)) < 0.1 * getattr(cn, "d" + p)


# --- The final model ---------------------------------------------------------

def test_without_the_model_the_fit_is_the_same():
    comps, systs = _doublet_and_lya()
    spec = _spec(45000.0, comps)
    with_model, model, _ = VoigtFitter(spec, systs).fit()
    without, none, _ = VoigtFitter(spec, systs).fit(model=False)
    assert none is None and model is not None
    assert with_model.components == without.components
    np.testing.assert_array_equal(model, VoigtFitter(spec, with_model).compute_model_flux()[1])
