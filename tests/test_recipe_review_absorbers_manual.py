"""Review of the manual fits of Absorbers (phase 6-quater, PR 10).

``add_component``, ``update_component``, ``update_components``,
``delete_component``, ``detect_anchor``, ``fit_component``,
``update_constraint``, ``update_constraints_column``,
``add_linked_system``: the behaviour kept and the behaviour changed, each
test saying what it was. The values of the valid results were taken on
master before the change (``tests/snapshots/manual_fits_*.npz``): the
numbers of a valid call do not move, except where a finding says so (F145,
F146). They are compared with a relative tolerance of 1e-6, not to the
bit: the fits give the last digits differently on different processors
(an md5 of the values, as in PR 9, failed on a runner of the CI).
"""
import hashlib
import re
import warnings

import numpy as np
import pytest

from test_recipe_inputs_absorbers import _sessions


def _quiet(run):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return run()


def _values(session):
    """The components (values, errors, resolution: an array, NaN for None,
    and the series), the constraints (strings) and the model of a session,
    for a snapshot, with the UUIDs replaced by the positions of the
    components (a new component has a new random UUID)."""
    comps = session.systs.components if session.systs else []
    index = {c.uuid: str(i) for i, c in enumerate(comps)}

    def plain(text):
        if text is None:
            return "None"
        return re.sub(r"[0-9a-f]{8}-[0-9a-f-]{27}", lambda m: index.get(m.group(), "?"), text)

    numbers = np.array([[np.nan if v is None else float(v)
                         for v in (c.z, c.dz, c.logN, c.dlogN, c.b, c.db, c.btur, c.dbtur,
                                   getattr(c, "resol", None))] for c in comps]).reshape(-1, 9)
    cons = sorted(f"{index.get(u, '?')} {p} {k.is_free} {plain(k.expression)} "
                  f"{plain(k.target_uuid)}"
                  for (u, p), k in (session.systs.parsed_constraints.items()
                                    if session.systs else []))
    values = {"components": numbers, "series": np.array([c.series for c in comps], dtype=str),
              "constraints": np.array(cons, dtype=str)}
    if session.spec.has_aux_column("model"):
        values["model"] = np.asarray(session.spec.t["model"], dtype=float)
    return values


def _check(snapshot, case, session):
    snapshot.check(f"manual_fits_{case}", _values(session), rtol=1e-6, atol=1e-12)


def _fingerprint(session):
    """The components (values, errors, resolution), the constraints and the
    model of a session, with the UUIDs replaced by the positions of the
    components (a new component has a new random UUID), as an md5: to
    compare two sessions of the same run only (across processors the last
    digits of a fit differ: :func:`_values` and the snapshots)."""
    comps = session.systs.components if session.systs else []
    index = {c.uuid: str(i) for i, c in enumerate(comps)}

    def plain(text):
        if text is None:
            return None
        return re.sub(r"[0-9a-f]{8}-[0-9a-f-]{27}", lambda m: index.get(m.group(), "?"), text)

    digest = hashlib.md5()
    for c in comps:
        values = [c.z, c.dz, c.logN, c.dlogN, c.b, c.db, c.btur, c.dbtur,
                  getattr(c, "resol", None)]
        digest.update(c.series.encode())
        digest.update(repr([None if v is None else round(float(v), 9) for v in values]).encode())
    if session.systs:
        cons = sorted((index.get(u, "?"), p, k.is_free, plain(k.expression), plain(k.target_uuid))
                      for (u, p), k in session.systs.parsed_constraints.items())
        digest.update(repr(cons).encode())
    if session.spec.has_aux_column("model"):
        model = np.asarray(session.spec.t["model"], dtype=float)
        digest.update(np.round(np.nan_to_num(model, nan=-9e9), 9).tobytes())
    return digest.hexdigest()


def _ids():
    sessions, values = _sessions()
    fitted = sessions["fitted"]
    u = values["UUID"]
    lya = [c.uuid for c in fitted.systs.components if c.series == "Ly_a"]
    return sessions, u, lya


def _linked(s, u, lya):
    """The GUI's link: the b of the second Ly-a to the first, with its target."""
    return s.absorbers.update_constraint(lya[1], "b", is_free=False,
                                         expression=f"p['{lya[0]}'].b", target_uuid=lya[0])


#: Valid calls, on the sessions of test_recipe_inputs_absorbers (the mock
#: of a CIV doublet and two Ly-a, fitted by components_auto)
KEPT = {
    "add_component": lambda S, u, lya: S["base"].absorbers.add_component("CIV", 1.931),
    "add_component_typed": lambda S, u, lya: S["fitted"].absorbers.add_component(
        series="CIV", z=1.9302, logN=13.0, b=8.0, btur=2.0, z_window_kms=30.0),
    "add_component_custom": lambda S, u, lya: S["base"].absorbers.add_component(
        "CIV_1548,CIV_1550", 1.931),
    "add_component_strings": lambda S, u, lya: S["base"].absorbers.add_component(
        "CIV", "1.931", "13.5", "10.0", "0.0", "20.0"),
    "update_component": lambda S, u, lya: S["fitted"].absorbers.update_component(
        u, z=1.9301, logN=13.7, b=12.0, btur=1.0),
    "update_component_series": lambda S, u, lya: S["fitted"].absorbers.update_component(
        lya[0], series="Ly_b"),
    "update_component_resol": lambda S, u, lya: S["fitted"].absorbers.update_component(
        u, resol=30000.0),
    "update_components": lambda S, u, lya: S["fitted"].absorbers.update_components(
        lya, b=20.0),
    "update_components_series": lambda S, u, lya: S["fitted"].absorbers.update_components(
        lya, series="Ly_b"),
    "delete_component": lambda S, u, lya: S["fitted"].absorbers.delete_component(uuid=u),
    "delete_components": lambda S, u, lya: S["fitted"].absorbers.delete_component(uuids=lya),
    "delete_both": lambda S, u, lya: S["fitted"].absorbers.delete_component(
        uuid=u, uuids=[lya[0]]),
    "fit_component": lambda S, u, lya: S["fitted"].absorbers.fit_component(u),
    "fit_component_options": lambda S, u, lya: S["fitted"].absorbers.fit_component(
        lya[0], max_nfev=50, z_window_kms=50.0, group_depth=0),
    "freeze": lambda S, u, lya: S["fitted"].absorbers.update_constraint(u, "z", is_free=False),
    "free_btur": lambda S, u, lya: S["fitted"].absorbers.update_constraint(
        u, "btur", is_free=True),
    "link_gui": lambda S, u, lya: _linked(S["fitted"], u, lya),
    "link_gui_fitted": lambda S, u, lya: _linked(S["fitted"], u, lya).absorbers.fit_component(
        lya[1]),
    "unlink": lambda S, u, lya: _linked(S["fitted"], u, lya).absorbers.update_constraint(
        lya[1], "b", is_free=True),
    "link_then_delete_target": lambda S, u, lya: _linked(
        S["fitted"], u, lya).absorbers.delete_component(uuid=lya[0]),
    "freeze_column": lambda S, u, lya: S["fitted"].absorbers.update_constraints_column(
        "b", False),
    "free_column": lambda S, u, lya: S["fitted"].absorbers.update_constraints_column(
        "btur", True),
    "add_linked_system": lambda S, u, lya: S["base"].absorbers.add_linked_system(
        "CIV_1548,CIV_1550", 1.931),
    "add_linked_system_one": lambda S, u, lya: S["base"].absorbers.add_linked_system(
        "CIV", 1.931, logN=13.0, b=12.0, z_window_kms=30.0),
}



@pytest.mark.parametrize("case", sorted(KEPT))
def test_valid_results_unchanged(case, snapshot):
    """The checks added by PR 10 leave the results of the valid calls as
    they were (to the bit here, within 1e-6 on another processor)."""
    sessions, u, lya = _ids()
    result = _quiet(lambda: KEPT[case](sessions, u, lya))
    _check(snapshot, case, result)


def test_detect_anchor_on_a_normalized_flux_unchanged():
    # The mock has a continuum at 1: y and y/cont give the same anchor
    sessions, u, lya = _ids()
    fitted = sessions["fitted"]
    z = fitted.absorbers.detect_anchor(fitted, "CIV_1548", "CIV_1550")
    assert z == pytest.approx(DETECT_ANCHOR_Z, abs=0, rel=1e-6)


DETECT_ANCHOR_Z = 1.930000147166357  # master, 12 November 2026


# --- Kept: the old logs ---------------------------------------------------

def test_update_components_series_none_string_is_unchanged(snapshot):
    # 'None' of the old logs is "unchanged", as before (series was not in
    # the registry, a string); with another parameter, the same result
    sessions, u, lya = _ids()
    fitted = sessions["fitted"]
    old = _quiet(lambda: fitted.absorbers.update_components(lya, b=20.0, series="None"))
    _check(snapshot, "update_components", old)


def test_typed_recipe_defaults_are_those_of_the_signature():
    # One default for the GUI and the API: a parameter the method requires
    # is required in the registry ("" or None), and the MCP schema lists it
    # (the registry had defaults the method did not, F144)
    import inspect

    from astrocook.mcp_server.schemas import parameters_schema
    from astrocook.recipes.absorbers import RecipeAbsorbers
    from astrocook.recipes.registry import get_recipe
    for name in RECIPES:
        spec = get_recipe("absorbers", name)
        signature = inspect.signature(getattr(RecipeAbsorbers, name)).parameters
        required = parameters_schema(spec)["required"]
        for p in spec.params:
            if signature[p.name].default is inspect.Parameter.empty:
                assert p.default in ("", None) and not p.optional, (name, p.name)
            if p.agents and p.default in ("", None) and not p.optional:
                assert p.name in required, (name, p.name)
    assert "z" in parameters_schema(get_recipe("absorbers", "add_component"))["required"]


RECIPES = ("add_component", "update_component", "update_components", "delete_component",
           "detect_anchor", "fit_component", "update_constraint",
           "update_constraints_column", "add_linked_system")


# --- Changed: refused before the recipe runs (F142-F146) ------------------

from astrocook.recipes.errors import InvalidParameter, RecipeFailed  # noqa: E402

UNKNOWN = "00000000-0000-0000-0000-000000000000"

REFUSED = {
    # F142: an unknown component; each a state that changed nothing before
    # (update_component with a WARNING; update_constraint stored a
    # constraint for no component)
    "update_component_unknown": (
        lambda S, u, lya: S["fitted"].absorbers.update_component(UNKNOWN, z=1.93),
        "uuid: there is no component"),
    "update_components_unknown": (
        lambda S, u, lya: S["fitted"].absorbers.update_components([u, UNKNOWN], z=1.93),
        "uuids: there is no component '0000"),
    "update_components_empty": (
        lambda S, u, lya: S["fitted"].absorbers.update_components([], z=1.93),
        "uuids: give at least one"),
    "delete_unknown": (
        lambda S, u, lya: S["fitted"].absorbers.delete_component(uuid=UNKNOWN),
        "uuid: there is no component"),
    "delete_unknown_list": (
        lambda S, u, lya: S["fitted"].absorbers.delete_component(uuids=[UNKNOWN]),
        "uuids: there is no component"),
    "delete_nothing": (
        lambda S, u, lya: S["fitted"].absorbers.delete_component(),
        "give the component to delete"),
    "fit_unknown": (
        lambda S, u, lya: S["fitted"].absorbers.fit_component(UNKNOWN),
        "uuid: there is no component"),
    "constraint_unknown": (
        lambda S, u, lya: S["fitted"].absorbers.update_constraint(UNKNOWN, "z", is_free=False),
        "uuid: there is no component"),
    # Nothing to change: a state that changed nothing (update_components
    # returned the same session)
    "update_component_nothing": (
        lambda S, u, lya: S["fitted"].absorbers.update_component(u),
        "give at least one parameter"),
    "update_components_nothing": (
        lambda S, u, lya: S["fitted"].absorbers.update_components([u]),
        "give at least one parameter"),
    "constraint_nothing": (
        lambda S, u, lya: S["fitted"].absorbers.update_constraint(u, "z"),
        "give is_free, or an expression"),
    # F143: an unknown series, stored (a component without lines)
    "update_component_series": (
        lambda S, u, lya: S["fitted"].absorbers.update_component(u, series="nope"),
        "series: 'nope' is not a known"),
    "update_components_series": (
        lambda S, u, lya: S["fitted"].absorbers.update_components([u], series="nope"),
        "series: 'nope' is not a known"),
    "add_linked_unknown": (
        lambda S, u, lya: S["base"].absorbers.add_linked_system("CIV,nope", 1.931),
        "series_list: 'nope' is not a known"),
    "add_linked_twice": (
        lambda S, u, lya: S["base"].absorbers.add_linked_system("CIV,CIV", 1.931),
        "series_list: CIV given twice"),
    # add_component: failed with "No valid pixels" (F111) after adding it
    "add_component_unknown": (
        lambda S, u, lya: S["base"].absorbers.add_component("nope", 1.931),
        "series: 'nope' is not a known"),
    "add_component_series_list": (
        lambda S, u, lya: S["base"].absorbers.add_component("CIV,SiIV", 1.931),
        "use add_linked_system"),
    "add_component_empty_series": (
        lambda S, u, lya: S["base"].absorbers.add_component("", 1.931),
        "series: give a transition"),
    # z required: the default 0.0 failed with "No valid pixels"
    "add_component_no_z": (
        lambda S, u, lya: S["base"].absorbers.add_component("CIV"),
        "z: a value is needed"),
    # detect_anchor: an unknown transition failed with "'nope' not found"
    "anchor_unknown": (
        lambda S, u, lya: S["fitted"].absorbers.detect_anchor(S["fitted"], "nope", "CIV_1550"),
        "trans_self: 'nope' is not a known transition"),
    # F146: the links of update_constraint
    "link_and_free": (
        lambda S, u, lya: S["fitted"].absorbers.update_constraint(
            lya[1], "b", is_free=True, expression=f"p['{lya[0]}'].b"),
        "is_free: a linked parameter is not free"),
    "link_to_itself": (
        lambda S, u, lya: S["fitted"].absorbers.update_constraint(
            u, "z", is_free=False, expression=f"p['{u}'].z + 0.001"),
        "links the component to itself"),
    "link_other_target": (
        lambda S, u, lya: S["fitted"].absorbers.update_constraint(
            lya[1], "b", is_free=False, expression=f"p['{lya[0]}'].b", target_uuid=u),
        "target_uuid: .* is not the component of the expression"),
    "link_unknown_target": (
        lambda S, u, lya: S["fitted"].absorbers.update_constraint(
            lya[1], "b", is_free=False, expression=f"p['{lya[0]}'].b", target_uuid=UNKNOWN),
        "target_uuid: there is no component"),
    "target_without_link": (
        lambda S, u, lya: S["fitted"].absorbers.update_constraint(
            lya[1], "b", is_free=False, target_uuid=lya[0]),
        "target_uuid: there is no expression"),
    # F117: InvalidParameter (a RecipeFailed from the core before)
    "link_not_evaluated": (
        lambda S, u, lya: S["fitted"].absorbers.update_constraint(
            lya[1], "z", expression='p["z"]'),
        "expression: The link of z"),
}


@pytest.mark.parametrize("case", sorted(REFUSED))
def test_refused_before_the_recipe_runs(case, caplog):
    sessions, u, lya = _ids()
    run, message = REFUSED[case]
    before = _fingerprint(sessions["fitted"])
    with pytest.raises(InvalidParameter, match=message):
        _quiet(lambda: run(sessions, u, lya))
    assert _fingerprint(sessions["fitted"]) == before
    assert "WARNING" not in caplog.text


def test_required_parameters_of_the_api():
    # The registry offered defaults that the method did not have: an agent
    # that omitted them got a TypeError (F144); now required on both sides
    sessions, u, lya = _ids()
    fitted = sessions["fitted"]
    for call in (lambda: fitted.absorbers.detect_anchor(fitted),
                 lambda: fitted.absorbers.update_constraints_column("b"),
                 lambda: sessions["base"].absorbers.add_linked_system("CIV")):
        with pytest.raises(TypeError, match="missing a required argument"):
            call()


# --- Changed: what a valid call gives -------------------------------------

def test_resol_zero_removes_the_override():
    # 0 was refused (gt 0): the override of a component could not be removed
    sessions, u, lya = _ids()
    fitted = sessions["fitted"]
    set_ = _quiet(lambda: fitted.absorbers.update_component(u, resol=30000.0))
    assert set_.systs.get_component_by_uuid(u).resol == 30000.0
    back = _quiet(lambda: set_.absorbers.update_component(u, resol=0))
    assert back.systs.get_component_by_uuid(u).resol is None
    assert _fingerprint(back) == _fingerprint(fitted)
    both = _quiet(lambda: set_.absorbers.update_components([u], resol="0"))
    assert both.systs.get_component_by_uuid(u).resol is None


def test_link_without_target_takes_that_of_the_expression(snapshot):
    # F146: without target_uuid the target was not stored, and the fits
    # left it out of the group of the component; now it is the component
    # the expression names, as the GUI gives it
    sessions, u, lya = _ids()
    fitted = sessions["fitted"]
    new = fitted.absorbers.update_constraint(lya[1], "b", expression=f"p['{lya[0]}'].b")
    _check(snapshot, "link_gui", new)
    constraint = new.systs.parsed_constraints[(lya[1], "b")]
    # stored frozen, as the fit treats it (it was stored free, with the
    # expression, when is_free was left out)
    assert (constraint.is_free, constraint.target_uuid) == (False, lya[0])


def test_constant_link_has_no_target():
    sessions, u, lya = _ids()
    new = sessions["fitted"].absorbers.update_constraint(u, "b", expression="12.5")
    constraint = new.systs.parsed_constraints[(u, "b")]
    assert (constraint.is_free, constraint.expression, constraint.target_uuid) == (
        False, "12.5", None)
    assert new.systs.get_component_by_uuid(u).b == 12.5


def _anchor_session(scale=1.0, cont=True, xunit=None):
    """The fitted mock with its flux (and continuum) times ``scale``,
    without the continuum, or with its axis in ``xunit``."""
    import astropy.units as au
    sessions, u, lya = _ids()
    spec = sessions["fitted"].spec
    if xunit is not None:
        spec = spec.x_convert(zem=2.7, xunit=xunit)
    data = spec._data
    aux = dict(data.aux_cols)
    if not cont:
        aux.pop("cont")
    elif scale != 1.0:
        aux["cont"] = type(aux["cont"])(aux["cont"].values * scale, au.erg / au.cm**2 / au.s / au.AA)
    import dataclasses
    data = dataclasses.replace(
        data, y=type(data.y)(data.y.values * scale, au.erg / au.cm**2 / au.s / au.AA
                             if scale != 1.0 else data.y.unit),
        dy=type(data.dy)(data.dy.values * scale, au.erg / au.cm**2 / au.s / au.AA
                         if scale != 1.0 else data.dy.unit),
        aux_cols=aux)
    return sessions["fitted"].with_new_spectrum(type(spec)(data))


def test_detect_anchor_reads_the_normalized_flux():
    # F145: the depths were 1 - y: with a flux of 1e-16 every pixel was
    # "absorbed", and the anchor was the start of the common range
    session = _anchor_session(scale=1e-16)
    z = session.absorbers.detect_anchor(session, "CIV_1548", "CIV_1550")
    assert z == pytest.approx(DETECT_ANCHOR_Z, abs=0, rel=1e-6)


def test_detect_anchor_without_a_continuum_fails():
    session = _anchor_session(scale=1e-16, cont=False)
    with pytest.raises(RecipeFailed, match="not normalized and it has no continuum"):
        session.absorbers.detect_anchor(session, "CIV_1548", "CIV_1550")


def test_detect_anchor_without_a_common_range_fails():
    # F145: 0.0 was returned, as a redshift
    sessions, u, lya = _ids()
    fitted = sessions["fitted"]
    with pytest.raises(RecipeFailed, match="cover no redshift in common"):
        fitted.absorbers.detect_anchor(fitted, "CIV_1548", "Ly_a")
