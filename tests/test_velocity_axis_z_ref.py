"""The reference redshift of an x axis in velocity (F84, phase 7-bis, B6).

``x_convert`` to a velocity puts the zero of the axis at the Lyman-alpha
emission at ``z_ref``, 121.567 nm * (1 + z_ref). Before B6 the spectrum
did not keep ``z_ref`` (only the line of its history did), and every
conversion back used ``z_em``: wrong when the axis had been converted at
another z_ref, or z_em had changed since. Since B6 the spectrum keeps it
(``SpectrumData.z_ref``, ``Spectrum.axis_z_ref``), the archive too (FITS key
``Z_REF``), and the conversions back use it (the design agreed with Guido,
16 November 2026). Each test that changed with B6 says what it was.
"""
import logging

import numpy as np
import astropy.units as au
import pytest

from astrocook.core.spectrum_operations import convert_axis_velocity
from astrocook.recipes.errors import InvalidParameter, RecipeFailed

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


def test_the_recipe_converts_back_at_the_z_ref_of_the_axis():
    """Before B6: at z_em ('_current_' in both directions), so the axis came
    back shifted by (1 + z_em) / (1 + z_ref)."""
    session = _session()
    x_nm = session.spec.x.value
    conv = session.edit.x_convert(z_ref=Z_REF, xunit='km/s')
    assert conv.spec._data.z_ref == Z_REF and conv.spec.axis_z_ref() == Z_REF
    back = conv.edit.x_convert(xunit='nm')
    np.testing.assert_allclose(back.spec.x.value, x_nm, rtol=1e-12)
    assert back.spec._data.z_ref is None
    # The same z_ref, given: accepted
    again = conv.edit.x_convert(z_ref=str(Z_REF), xunit='Angstrom')
    np.testing.assert_allclose(again.spec.x.value, x_nm * 10, rtol=1e-12)


def test_back_at_another_z_ref_is_refused():
    conv = _session().edit.x_convert(z_ref=Z_REF, xunit='km/s')
    with pytest.raises(InvalidParameter, match="converted at z_ref = 3"):
        conv.edit.x_convert(z_ref='2.0', xunit='nm')
    with pytest.raises(ValueError, match="converted back at the same"):
        conv.spec.x_convert(2.0, au.nm)


def test_a_change_of_z_em_after_the_conversion_does_not_move_the_axis():
    """Before B6: the conversion back used the new z_em (4.5 / 5.1 here)."""
    session = _session()
    x_nm = session.spec.x.value
    conv = session.edit.x_convert(xunit='km/s')  # at z_em = 4.1
    assert conv.spec._data.z_ref == Z_EM
    moved = conv.edit.set_properties(z_em='3.5')
    assert moved.spec._data.z_ref == Z_EM
    back = moved.edit.x_convert(xunit='nm')
    np.testing.assert_allclose(back.spec.x.value, x_nm, rtol=1e-12)


def test_a_velocity_to_velocity_conversion_recentres_the_axis():
    """Before B6: km/s to km/s at another z_ref left the values as they
    were (no new zero), though the history said the new z_ref."""
    spec = _spec()
    conv = spec.x_convert(Z_REF, KMS)
    again = conv.x_convert(1.0, KMS)
    np.testing.assert_allclose(again.x.value, _velocity(spec.x.value, 1.0), rtol=1e-12)
    np.testing.assert_allclose(again.xmax.value, _velocity(spec.xmax.value, 1.0), rtol=1e-12)
    assert again._data.z_ref == 1.0
    assert again.history[-1].endswith("at z_ref=1.0")
    np.testing.assert_allclose(again.x_convert(None, au.nm).x.value, spec.x.value, rtol=1e-12)
    # Only the unit changes without a z_ref ('_current_' in the recipe)
    in_m_s = conv.x_convert(None, au.m / au.s)
    np.testing.assert_allclose(in_m_s.x.value, conv.x.value * 1e3, rtol=1e-12)
    assert in_m_s._data.z_ref == Z_REF


def test_the_recipe_recentres_and_keeps_the_unit_by_default():
    conv = _session().edit.x_convert(z_ref=Z_REF, xunit='km/s')
    same = conv.edit.x_convert(xunit='km/s')
    np.testing.assert_array_equal(same.spec.x.value, conv.spec.x.value)
    recentred = conv.edit.x_convert(z_ref='0', xunit='km/s')
    assert recentred.spec._data.z_ref == 0.0
    np.testing.assert_allclose(recentred.spec.x.value - conv.spec.x.value,
                               np.log(1 + Z_REF) * C_KMS, rtol=1e-9)


def test_an_axis_in_velocity_without_a_reference_falls_back_on_z_em(caplog):
    """An axis loaded in velocity, or from an archive saved before B6:
    z_em, as all the conversions before B6, with a warning; a z_ref given
    for the conversion back declares the zero."""
    import dataclasses
    from astrocook.core.spectrum import Spectrum
    spec = _spec()
    conv = spec.x_convert(Z_EM, KMS)
    unknown = Spectrum(dataclasses.replace(conv._data, z_ref=None))
    assert unknown.axis_z_ref() == Z_EM
    logging.disable(logging.NOTSET)
    with caplog.at_level(logging.WARNING, logger="astrocook.core.spectrum"):
        back = unknown.x_convert(None, au.nm)
    assert "not known: z_em (4.1) is used" in caplog.text
    np.testing.assert_allclose(back.x.value, spec.x.value, rtol=1e-12)
    declared = Spectrum(dataclasses.replace(spec.x_convert(Z_REF, KMS)._data, z_ref=None))
    np.testing.assert_allclose(declared.x_convert(Z_REF, au.nm).x.value, spec.x.value,
                               rtol=1e-12)


def test_the_axis_of_a_wavelength_has_no_reference():
    spec = _spec()
    assert spec._data.z_ref is None and spec.axis_z_ref() == Z_EM
    assert spec.x_convert(None, au.AA)._data.z_ref is None
    np.testing.assert_allclose(spec.x_convert(None, KMS).x.value,
                               _velocity(spec.x.value, Z_EM), rtol=1e-12)


def test_the_reference_follows_the_operations_on_the_spectrum():
    conv = _session().edit.x_convert(z_ref=Z_REF, xunit='km/s')
    steps = [
        conv.flux.smooth(sigma_kms=10.0),
        conv.edit.smooth_column(),
        conv.flux.calculate_running_std(),
        conv.edit.y_convert(yunit='erg/(Angstrom s cm^2)'),
        conv.edit.set_properties(z_em='3.5'),
    ]
    for step in steps:
        assert step.spec._data.z_ref == Z_REF, step.spec.history[-1:]
    spec = conv.spec
    assert spec.apply_expression('model', 'model * 1')._data.z_ref == Z_REF


def test_undo_restores_the_reference():
    from astrocook.core.structures import HistoryLog
    from astrocook.session.history import SessionHistory
    session = _session()
    history = SessionHistory(session, HistoryLog())
    history.add_state(session.edit.x_convert(z_ref=Z_REF, xunit='km/s'))
    assert history.current_state.spec._data.z_ref == Z_REF
    history.undo()
    assert history.current_state.spec._data.z_ref is None


def test_the_archive_keeps_the_reference(tmp_path):
    """Before B6: the unit was kept, z_ref was not."""
    conv = _session().edit.x_convert(z_ref=Z_REF, xunit='km/s')
    reloaded, _ = _round_trip(conv, tmp_path)
    assert reloaded.spec.x.unit == KMS
    np.testing.assert_array_equal(reloaded.spec.x.value, conv.spec.x.value)
    assert reloaded.spec._data.z_em == Z_EM
    assert reloaded.spec._data.z_ref == Z_REF
    assert not any(k.upper() == 'Z_REF' for k in reloaded.spec._data.meta)
    back = reloaded.edit.x_convert(xunit='nm')
    np.testing.assert_allclose(back.spec.x.value, _session().spec.x.value, rtol=1e-12)


def test_an_archive_in_wavelength_writes_no_reference(tmp_path):
    from astropy.io import fits
    import tarfile
    reloaded, path = _round_trip(_session(), tmp_path)
    assert reloaded.spec._data.z_ref is None
    with tarfile.open(path) as tar:
        names = [m for m in tar.getmembers() if m.name.endswith('_spec.fits')]
        assert names
        with fits.open(tar.extractfile(names[0])) as hdul:
            assert 'Z_REF' not in hdul[1].header and 'Z_EM' in hdul[1].header


def test_the_exported_csv_keeps_the_reference(tmp_path):
    from astrocook.io.loaders import load_simple_csv
    conv = _session().edit.x_convert(z_ref=Z_REF, xunit='km/s')
    conv.file.export_ascii(targets='spec', folder=str(tmp_path), prefix='p')
    data = load_simple_csv(str(tmp_path / 'p_spec.csv'))
    assert data.x.unit == KMS and data.z_ref == Z_REF
    _session().file.export_ascii(targets='spec', folder=str(tmp_path), prefix='w')
    assert load_simple_csv(str(tmp_path / 'w_spec.csv')).z_ref is None


def test_the_mcp_summary_says_the_reference():
    from astrocook.mcp_server.outputs import spectrum_summary
    assert spectrum_summary(_session()).z_ref is None
    conv = _session().edit.x_convert(z_ref=Z_REF, xunit='km/s')
    assert spectrum_summary(conv).z_ref == Z_REF


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


def test_combining_axes_in_velocity_with_different_references_is_refused():
    """Before B6: two spectra in km/s were joined by their values, though
    their zeros differed (here by (1 + 4.1) / (1 + 3.0) in wavelength)."""
    session = _session()
    a = session.edit.x_convert(z_ref=Z_REF, xunit='km/s')
    b = session.edit.x_convert(xunit='km/s')  # at z_em
    with pytest.raises(RecipeFailed, match="zero at z_ref = 4.1 and that of main at z_ref = 3"):
        a.edit.coadd(session_names=[a, b], dx=5.0)
    with pytest.raises(RecipeFailed, match="zero at z_ref"):
        a.edit.stitch(other_sessions=[b])
    with pytest.raises(RecipeFailed, match="zero at z_ref"):
        a.flux.equalize(reference_session=b)


def test_coadd_of_axes_in_velocity_with_the_same_reference():
    """The zero of the velocity of the coadd is that of the axes (before
    B6: meta['z_em'], 0 for these spectra): an xstart in nm is converted at
    it, and the result keeps it."""
    session = _session()
    a = session.edit.x_convert(z_ref=Z_REF, xunit='km/s')
    out = a.edit.coadd(session_names=[a, a], xstart=614.0, dx=5.0)
    assert out.spec.x.unit == KMS and out.spec._data.z_ref == Z_REF
    assert out.spec.xmin.value[0] == pytest.approx(_velocity(614.0, Z_REF), abs=1e-6)
