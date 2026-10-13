"""Numerical snapshots of the pure functions in core/spectrum_operations.py.

Characterization tests: they record the current output on a fixed synthetic
spectrum (see ``conftest.py`` for how snapshots are stored and updated).
"""
import numpy as np
import astropy.units as au
import pytest

from astrocook.core.spectrum_operations import (
    compute_flux_scaling, convert_x_axis, convert_y_axis, rebin_spectrum,
    running_mean, running_std, smooth_spectrum,
)
from astrocook.core.structures import DataColumnV2, SpectrumDataV2

FLUX_UNIT = au.erg / (au.cm**2 * au.s * au.nm)


def _synthetic_spectrum(seed=0, n=2000, z_em=0.0):
    """Sloped continuum, three Gaussian absorptions, Gaussian noise, two spikes."""
    rng = np.random.default_rng(seed)
    x = np.linspace(400.0, 410.0, n)
    dx = x[1] - x[0]
    cont = 1.0 + 0.02 * (x - 400.0)
    absorption = np.ones_like(x)
    for centre, depth, width in [(402.0, 0.6, 0.02), (405.5, 0.3, 0.05),
                                 (408.0, 0.9, 0.01)]:
        absorption -= depth * np.exp(-0.5 * ((x - centre) / width) ** 2)
    dy = 0.02 * cont
    y = cont * absorption + rng.normal(0.0, 1.0, n) * dy
    y[[500, 1500]] += 1.0  # spikes, to exercise sigma clipping
    return SpectrumDataV2(
        x=DataColumnV2(x, au.nm),
        xmin=DataColumnV2(x - dx / 2, au.nm),
        xmax=DataColumnV2(x + dx / 2, au.nm),
        y=DataColumnV2(y, FLUX_UNIT),
        dy=DataColumnV2(dy, FLUX_UNIT),
        z_em=z_em,
    )


def _rebin(spec, dx, kappa, xstart=None, xend=None):
    x_out, xmin_out, xmax_out, y_out, dy_out, resol_out = rebin_spectrum(
        spec.x.quantity, spec.xmin.quantity, spec.xmax.quantity, spec,
        xstart=xstart, xend=xend, dx=dx, calc_unit=spec.x.quantity.unit,
        kappa=kappa, filling=np.nan)
    return {"x": x_out, "xmin": xmin_out, "xmax": xmax_out,
            "y": y_out, "dy": dy_out, "resol": resol_out}


@pytest.mark.parametrize("label, dx, kappa", [
    ("nm_no_clip", 0.015 * au.nm, None),
    ("nm_clip", 0.015 * au.nm, 3.0),
    ("kms_no_clip", 10.0 * au.km / au.s, None),
])
def test_rebin_spectrum(snapshot, label, dx, kappa):
    spec = _synthetic_spectrum(z_em=2.3)
    snapshot.check(f"rebin_spectrum_{label}", _rebin(spec, dx, kappa))


def test_rebin_spectrum_with_range(snapshot):
    spec = _synthetic_spectrum()
    out = _rebin(spec, 0.02 * au.nm, None,
                 xstart=401.0 * au.nm, xend=403.0 * au.nm)
    snapshot.check("rebin_spectrum_range", out)


@pytest.mark.parametrize("sigma_kms", [5.0, 50.0])
def test_smooth_spectrum(snapshot, sigma_kms):
    spec = _synthetic_spectrum()
    y = smooth_spectrum(spec.x.quantity, spec.y.value, sigma_kms)
    snapshot.check(f"smooth_spectrum_{sigma_kms:g}kms", {"y": y})


@pytest.mark.parametrize("h", [1, 10, 5000])  # 5000: window larger than array
def test_running_mean_and_std(snapshot, h):
    y = _synthetic_spectrum().y.value
    snapshot.check(f"running_mean_std_h{h}",
                   {"mean": running_mean(y, h), "std": running_std(y, h)})


@pytest.mark.parametrize("label, z_ref, xunit", [
    ("kms", 2.3, au.km / au.s),
    ("angstrom", 0.0, au.Angstrom),
])
def test_convert_x_axis(snapshot, label, z_ref, xunit):
    x, xmin, xmax = convert_x_axis(_synthetic_spectrum(), z_ref, xunit)
    snapshot.check(f"convert_x_axis_{label}", {"x": x, "xmin": xmin, "xmax": xmax})


def test_convert_y_axis(snapshot):
    spec = _synthetic_spectrum()
    # One flux-like auxiliary column (converted) and one that is not (copied)
    spec.aux_cols["cont"] = DataColumnV2(np.full(len(spec.x.value), 1.5), FLUX_UNIT)
    spec.aux_cols["resol"] = DataColumnV2(np.full(len(spec.x.value), 4.5e4),
                                          au.dimensionless_unscaled)
    y, dy, aux = convert_y_axis(spec, au.erg / (au.cm**2 * au.s * au.Angstrom))
    snapshot.check("convert_y_axis", {
        "y": y.quantity, "dy": dy.quantity,
        "cont": aux["cont"].quantity, "resol": aux["resol"].quantity})


@pytest.mark.parametrize("order", [0, 2, -1])
def test_compute_flux_scaling(snapshot, order):
    ref = _synthetic_spectrum(seed=1)
    other = _synthetic_spectrum(seed=2)
    # The second spectrum is offset in flux by a smooth factor
    y_in = other.y.value * (0.8 + 0.01 * (other.x.value - 400.0))
    correction = compute_flux_scaling(ref.x.value, ref.y.value,
                                      other.x.value, y_in, order=order)
    snapshot.check(f"compute_flux_scaling_order{order}", {"correction": correction})
