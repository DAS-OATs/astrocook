"""The reference redshift of an x axis in velocity (F84, phase 7-bis, B6).

``x_convert`` to a velocity puts the zero of the axis at the Lyman-alpha
emission at ``z_ref``, 121.567 nm * (1 + z_ref). Characterization of the
behaviour before B6: the spectrum does not keep ``z_ref`` (only the line of
its history does), and every conversion back uses ``z_em``.
"""
import logging

import numpy as np
import astropy.units as au
import pytest

from astrocook.core.spectrum_operations import convert_axis_velocity
from astrocook.recipes.errors import RecipeFailed

from test_acs2_round_trip import _round_trip
from test_recipes_edit_flux import Z_EM, _session, _spec

KMS = au.km / au.s
C_KMS = 299792.458
Z_REF = 3.0  # not Z_EM (4.1)


def _velocity(x_nm, z_ref):
    return np.log(x_nm / (121.567 * (1 + z_ref))) * C_KMS


@pytest.fixture(autouse=True)
def _quiet():
    logging.disable(logging.CRITICAL)
    yield
    logging.disable(logging.NOTSET)


def test_the_forward_conversion_uses_z_ref():
    spec = _spec()
    conv = spec.x_convert(Z_REF, KMS)
    assert conv.x.unit == KMS
    np.testing.assert_allclose(conv.x.value, _velocity(spec.x.value, Z_REF), rtol=1e-12)
    np.testing.assert_allclose(conv.xmin.value, _velocity(spec.xmin.value, Z_REF), rtol=1e-12)
    assert conv.history[-1] == f"X-axis converted to {KMS} at z_ref={Z_REF}"


def test_the_round_trip_at_the_same_z_ref_is_exact():
    spec = _spec()
    back = spec.x_convert(Z_REF, KMS).x_convert(Z_REF, au.nm)
    np.testing.assert_allclose(back.x.value, spec.x.value, rtol=1e-12)


def test_the_recipe_converts_back_at_z_em():
    """Before B6: the recipe's default z_ref ('_current_') is z_em in both
    directions, so an axis converted at another z_ref comes back shifted
    by (1 + z_em) / (1 + z_ref)."""
    session = _session()
    x_nm = session.spec.x.value
    conv = session.edit.x_convert(z_ref=Z_REF, xunit='km/s')
    back = conv.edit.x_convert(xunit='nm')
    np.testing.assert_allclose(back.spec.x.value, x_nm * (1 + Z_EM) / (1 + Z_REF), rtol=1e-12)


def test_a_change_of_z_em_after_the_conversion_moves_the_conversion_back():
    session = _session()
    x_nm = session.spec.x.value
    conv = session.edit.x_convert(xunit='km/s')  # at z_em = 4.1
    moved = conv.edit.set_properties(z_em='3.5')
    back = moved.edit.x_convert(xunit='nm')
    np.testing.assert_allclose(back.spec.x.value, x_nm * 4.5 / 5.1, rtol=1e-12)


def test_a_velocity_to_velocity_conversion_ignores_z_ref():
    """Before B6: km/s to km/s at another z_ref leaves the values as they
    are (no new zero), though the history says the new z_ref."""
    conv = _spec().x_convert(Z_REF, KMS)
    again = conv.x_convert(1.0, KMS)
    np.testing.assert_array_equal(again.x.value, conv.x.value)
    assert again.history[-1].endswith("at z_ref=1.0")
    np.testing.assert_allclose(conv.x_convert(Z_REF, au.m / au.s).x.value,
                               conv.x.value * 1e3, rtol=1e-12)


def test_the_spectrum_keeps_nothing_of_z_ref():
    spec = _spec()
    conv = spec.x_convert(Z_REF, KMS)
    assert conv._data.z_em == spec._data.z_em
    assert conv._data.meta == spec._data.meta


def test_the_archive_keeps_the_unit_but_not_z_ref(tmp_path):
    conv = _session().edit.x_convert(z_ref=Z_REF, xunit='km/s')
    reloaded, _ = _round_trip(conv, tmp_path)
    assert reloaded.spec.x.unit == KMS
    np.testing.assert_array_equal(reloaded.spec.x.value, conv.spec.x.value)
    assert reloaded.spec._data.z_em == Z_EM
    assert not any(k.upper() in ('Z_REF', 'X_ZREF') for k in reloaded.spec._data.meta)


def test_smoothing_does_not_depend_on_the_reference():
    """The kernel in km/s needs only differences of velocity: the same
    result for any zero of the axis, and for the axis in wavelength."""
    spec = _spec()
    smoothed = spec.smooth(20.0).y.value
    for z in (0.0, Z_REF, Z_EM):
        np.testing.assert_allclose(spec.x_convert(z, KMS).smooth(20.0).y.value,
                                   smoothed, rtol=1e-10)


def test_the_reference_does_not_change_velocity_differences():
    x = _spec().x
    d0 = np.diff(convert_axis_velocity(x, 0.0, KMS).value)
    d1 = np.diff(convert_axis_velocity(x, Z_REF, KMS).value)
    np.testing.assert_allclose(d1, d0, rtol=1e-9)


def test_rebin_of_an_axis_in_velocity_fails():
    """Before B6: the rebinned grid is taken in nm for the auxiliary
    columns, so a spectrum in velocity cannot be rebinned (a clean
    RecipeFailed, not a crash)."""
    conv = _session().edit.x_convert(z_ref=Z_REF, xunit='km/s')
    with pytest.raises(RecipeFailed, match="not convertible"):
        conv.flux.rebin(dx=5.0, xunit='km/s')


def test_coadd_joins_axes_in_velocity_with_different_references():
    """Before B6: two spectra in km/s are joined by their values, though
    their zeros differ (here by (1 + 4.1) / (1 + 3.0) in wavelength)."""
    session = _session()
    a = session.edit.x_convert(z_ref=Z_REF, xunit='km/s')
    b = session.edit.x_convert(xunit='km/s')  # at z_em
    out = a.edit.coadd(session_names=[a, b], dx=5.0)
    assert out.spec.x.unit == KMS
    lo = min(a.spec.x.value.min(), b.spec.x.value.min())
    hi = max(a.spec.x.value.max(), b.spec.x.value.max())
    assert out.spec.x.value.min() < lo + 10 and out.spec.x.value.max() > hi - 10
