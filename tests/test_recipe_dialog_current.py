"""Parameters of the recipe dialogs whose default is '_current_' (F63):
they show the current value of the session, and Run without editing them
changes nothing. Before, the dialog showed "0.0" for the parameters
without an attribute of the same name in the spectrum data, so that Run
renamed the session "0.0", wrote "0.0" into OBJECT and converted the
velocities relative to z = 0."""
import astropy.units as au
import numpy as np
import pytest

from astrocook.core.session import SessionV2
from astrocook.core.structures import DataColumnV2
from astrocook.mocks import MockComponent, MockTruth, generate


def _session(resol_col=None):
    spec = generate(MockTruth(components=[MockComponent("CIV", 2.97, 13.5, 10.0)],
                              lam_min_nm=612.0, lam_max_nm=618.0, dv_kms=2.5,
                              resol=45000.0, snr=30.0, z_em=3.2, seed=1))
    spec = spec.with_properties(z_em=3.2, object_name="QSO J1234")
    if resol_col is not None:
        import dataclasses
        from astrocook.core.spectrum import SpectrumV2
        aux = dict(spec._data.aux_cols)
        aux["resol"] = DataColumnV2(resol_col(len(spec.x)), au.dimensionless_unscaled)
        spec = SpectrumV2(dataclasses.replace(spec._data, aux_cols=aux))
    return SessionV2("my session", None, spec=spec)


def _run(session, category, name):
    """Open the dialog, press Run without editing; the parameters it sends."""
    from astrocook.gui.recipe_dialog import RecipeDialog
    dialog = RecipeDialog(category, name, session)
    sent = []
    dialog.recipe_requested.connect(lambda cat, rec, params, aliases: sent.append(params))
    dialog.accept()
    return dialog, sent[0]


def test_set_properties_unchanged(qapp):
    session = _session()
    dialog, params = _run(session, "edit", "set_properties")
    assert dialog.input_widgets["name"].text() == "my session"
    assert dialog.input_widgets["object"].text() == "QSO J1234"
    assert dialog.input_widgets["resol"].text() == "45000"
    out = session.edit.set_properties(**params)
    assert out.name == "my session"
    assert out.spec.meta["OBJECT"] == "QSO J1234"
    assert out.spec._data.z_em == pytest.approx(3.2)
    np.testing.assert_allclose(out.spec.resolution(), 45000.0)


def test_set_properties_keeps_a_varying_resolution(qapp):
    session = _session(resol_col=lambda n: np.linspace(40000.0, 50000.0, n))
    dialog, params = _run(session, "edit", "set_properties")
    assert dialog.input_widgets["resol"].text() == ""
    assert "varies" in dialog.input_widgets["resol"].placeholderText()
    assert params["resol"] == "_current_"
    out = session.edit.set_properties(**params)
    np.testing.assert_allclose(out.spec.resolution(), session.spec.resolution())


@pytest.mark.parametrize("name, param", [("x_convert", "z_ref"), ("trim_common", "z_target")])
def test_redshifts_default_to_z_em(qapp, name, param):
    dialog, params = _run(_session(), "edit", name)
    assert float(params[param]) == pytest.approx(3.2)
