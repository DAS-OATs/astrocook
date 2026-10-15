"""Numerical snapshots of the SpectrumV2 methods that combine spectra.

Covers stitch, coadd, resample_on_grid, get_resampled_column,
equalize_to_reference, split and scale_flux on synthetic exposures of the
same object: different grids, noise levels and flux scales, with a
``resol`` column and a ``cont`` column.
"""
import numpy as np
import astropy.units as au
import pytest

from astrocook.core.spectrum import SpectrumV2
from astrocook.core.structures import DataColumnV2, SpectrumDataV2

FLUX_UNIT = au.erg / (au.cm**2 * au.s * au.nm)


def _true_flux(x):
    cont = 1.0 + 0.03 * (x - 500.0)
    lines = np.ones_like(x)
    for centre, depth, width in [(501.0, 0.7, 0.02), (503.2, 0.4, 0.05),
                                 (506.5, 0.9, 0.015)]:
        lines -= depth * np.exp(-0.5 * ((x - centre) / width) ** 2)
    return cont * lines, cont


def _exposure(start, stop, step, scale=1.0, noise=0.02, resol=45000.0, seed=0):
    """One exposure: its own grid, flux scale, noise and resolution."""
    rng = np.random.default_rng(seed)
    x = np.arange(start, stop, step)
    flux, cont = _true_flux(x)
    y = scale * flux * (1 + noise * rng.standard_normal(len(x)))
    dy = scale * noise * flux + 1e-3
    return SpectrumV2(SpectrumDataV2(
        x=DataColumnV2(x, au.nm), xmin=DataColumnV2(x - step / 2, au.nm),
        xmax=DataColumnV2(x + step / 2, au.nm),
        y=DataColumnV2(y, FLUX_UNIT), dy=DataColumnV2(dy, FLUX_UNIT),
        aux_cols={
            "cont": DataColumnV2(scale * cont, FLUX_UNIT),
            "resol": DataColumnV2(np.full(len(x), resol), au.dimensionless_unscaled),
        },
        z_em=3.0))


@pytest.fixture(scope="module")
def exposures():
    return [
        _exposure(499.0, 508.0, 0.004, scale=1.0, noise=0.02, resol=45000.0, seed=1),
        _exposure(499.5, 508.5, 0.005, scale=0.8, noise=0.03, resol=42000.0, seed=2),
        _exposure(498.5, 507.5, 0.0045, scale=1.1, noise=0.025, resol=47500.0, seed=3),
    ]


def _columns(spec, names=("x", "xmin", "xmax", "y", "dy")):
    out = {}
    for name in names:
        out[name] = getattr(spec, name) if name in ("x", "xmin", "xmax", "y", "dy") \
            else spec.get_column(name)
    return out


def test_stitch(snapshot, exposures):
    blue, red = exposures[0].split("x <= 503.0"), exposures[1].split("x > 503.0")
    out = blue.stitch([red], sort=True)
    assert np.all(np.diff(out.x.value) > 0)
    snapshot.check("spectrum_stitch",
                   _columns(out, ("x", "y", "dy", "cont", "resol")))


@pytest.mark.parametrize("label, dx, kappa, equalize_order", [
    ("kms_default", 2.0 * au.km / au.s, 5.0, 0),
    ("kms_no_clip_no_equalize", 2.0 * au.km / au.s, None, -1),
    ("nm_equalize_linear", 0.005 * au.nm, 5.0, 1),
])
def test_coadd(snapshot, exposures, label, dx, kappa, equalize_order):
    ref, *others = exposures
    out = ref.coadd(others, xstart=None, xend=None, dx=dx, kappa=kappa,
                    equalize_order=equalize_order)
    resol = out.get_column("resol")
    values = _columns(out)
    values["resol"] = resol
    snapshot.check(f"spectrum_coadd_{label}", values)


def test_resample_on_grid(snapshot, exposures):
    out = exposures[1].resample_on_grid(exposures[0])
    np.testing.assert_array_equal(out.x.value, exposures[0].x.value)
    snapshot.check("spectrum_resample_on_grid",
                   _columns(out, ("x", "y", "dy", "cont", "resol")))


def test_get_resampled_column(snapshot, exposures):
    target = np.linspace(500.0, 507.0, 300)
    snapshot.check("spectrum_get_resampled_column", {
        "y": exposures[2].get_resampled_column("y", target),
        "cont": exposures[2].get_resampled_column("cont", target),
    })


@pytest.mark.parametrize("order", [0, 1])
def test_equalize_to_reference(snapshot, exposures, order):
    out = exposures[1].equalize_to_reference(exposures[0], order=order)
    snapshot.check(f"spectrum_equalize_order{order}",
                   _columns(out, ("y", "dy", "cont")))


def test_split_and_scale_flux(snapshot, exposures):
    part = exposures[0].split("(x > 501.0) & (x < 504.0)")
    scaled = part.scale_flux(2.5)
    snapshot.check("spectrum_split_scale",
                   _columns(scaled, ("x", "y", "dy", "cont", "resol")))


def test_voigt_model_on_coadded_spectrum_is_not_zero(exposures):
    from astrocook.core.structures import ComponentDataV2, SystemListDataV2
    from astrocook.core.system_list import SystemListV2
    from astrocook.fitting.voigt_fitter import VoigtFitterV2

    ref, *others = exposures
    coadd = ref.coadd(others, xstart=None, xend=None, dx=2.0 * au.km / au.s)
    comp = ComponentDataV2(id=1, z=501.0 / 154.8204 - 1, dz=None, logN=13.5,
                           dlogN=None, b=10.0, db=None, series="CIV")
    _, model = VoigtFitterV2(
        coadd, SystemListV2(SystemListDataV2(components=[comp]))).compute_model_flux()
    assert np.mean(model == 0) < 0.01
