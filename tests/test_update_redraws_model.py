"""F152 (phase 7-bis, B2): ``update_component(s)`` and ``update_constraint``
redraw the ``model`` column from the new parameters, without a fit.

Before B2 the model stayed as it was until the next fit, while the
components table showed the new values."""
import warnings

import numpy as np
import pytest

from test_recipe_inputs_absorbers import _sessions


@pytest.fixture(scope="module")
def fitted():
    sessions, values = _sessions()
    return sessions["fitted"], values["UUID"]


def _model(session):
    return np.asarray(session.spec.t["model"], dtype=float)


def _quiet(run):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return run()


def _expected(session):
    """The model of the components of ``session``, computed afresh."""
    from astrocook.fitting.voigt_fitter import VoigtFitter
    return VoigtFitter(session.spec, session.systs).compute_model_flux()[1]


def test_update_component_redraws_the_model(fitted):
    s, u = fitted
    new = _quiet(lambda: s.absorbers.update_component(u, logN=14.9))
    assert not np.allclose(_model(new), _model(s))
    np.testing.assert_allclose(_model(new), _expected(new))


def test_update_components_redraws_the_model(fitted):
    s, u = fitted
    new = _quiet(lambda: s.absorbers.update_components([u], b=25.0))
    assert not np.allclose(_model(new), _model(s))
    np.testing.assert_allclose(_model(new), _expected(new))


def test_update_constraint_with_a_link_redraws_the_model(fitted):
    s, u = fitted
    others = [c for c in s.systs.components if c.uuid != u and c.series != "CIV"]
    target = others[0].uuid
    new = _quiet(lambda: s.absorbers.update_constraint(
        u, "b", expression=f"p['{target}'].b * 2"))
    np.testing.assert_allclose(_model(new), _expected(new))


def test_freezing_leaves_the_model_as_it_is(fitted):
    s, u = fitted
    new = _quiet(lambda: s.absorbers.update_constraint(u, "b", is_free=False))
    np.testing.assert_array_equal(_model(new), _model(s))


def test_update_without_a_model_column_adds_none(fitted):
    """A session never fitted has no model column: update_component does not
    create one (it would show a model nobody asked for)."""
    s, u = fitted
    bare = s.with_new_spectrum(s.spec.remove_columns(["model"]))
    new = _quiet(lambda: bare.absorbers.update_component(u, logN=14.9))
    assert not new.spec.has_aux_column("model")


def test_the_parameters_are_those_asked_for(fitted):
    s, u = fitted
    new = _quiet(lambda: s.absorbers.update_component(u, logN=14.9))
    comp = next(c for c in new.systs.components if c.uuid == u)
    assert comp.logN == pytest.approx(14.9)
