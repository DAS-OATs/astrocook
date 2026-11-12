"""Review of the manual fits of Absorbers (phase 6-quater, PR 10).

``add_component``, ``update_component``, ``update_components``,
``delete_component``, ``detect_anchor``, ``fit_component``,
``update_constraint``, ``update_constraints_column``,
``add_linked_system``: the behaviour kept and the behaviour changed, each
test saying what it was. The fingerprints of the valid results were taken
on master before the change: the numbers of a valid call do not move,
except where a finding says so (F145, F146).
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


def _fingerprint(session):
    """The components (values, errors, resolution), the constraints and the
    model of a session, with the UUIDs replaced by the positions of the
    components (a new component has a new random UUID)."""
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

FINGERPRINTS = {
    "add_component": "f8d2dc34194a4d87c8a771dc43675bff",
    "add_component_custom": "ec13762c37461f96c5ef4f27b8035cb1",
    "add_component_strings": "f8d2dc34194a4d87c8a771dc43675bff",
    "add_component_typed": "a23ff0e347430aa4573614acd327ead4",
    "add_linked_system": "bff326baa8545f06916ca4d18de48947",
    "add_linked_system_one": "3785c6201f703b97bcaf2bd676b58969",
    "delete_both": "627ce184e915b90c69f9b5b8518d6ec5",
    "delete_component": "f891f55c60baabc7663cb841086575dc",
    "delete_components": "30149d000d79c9ed863649cc5104a313",
    "fit_component": "c3df43f3ba7e464ade09fc7660c68b6b",
    "fit_component_options": "6b464d33f10262a300da88735e855f81",
    "free_btur": "8da507d3ab33dba98fda8819ef072fb4",
    "free_column": "f787ea5af9cc1922715bde2d00e8fffd",
    "freeze": "397f1d494793583eb3e5e1642add23b4",
    "freeze_column": "61a4cfcd06a9cf6f184dacf30526a9e3",
    "link_gui": "a3a2076d8c085eca0e94cdd5f39b6cb9",
    "link_gui_fitted": "cd925dea0e279f7d75d5f1ac3d9dbaaa",
    "link_then_delete_target": "8f1c43a39c64499e75cab6e97a4c2370",
    "unlink": "04f86bf667925c62c59c78d35fa9cfab",
    "update_component": "8f53ac8f875d4cf2d821bec1b15a42e1",
    "update_component_resol": "764bb5a36ee6b46db8605f4dab134be9",
    "update_component_series": "322c8b3408ff8a1a0a27dfbf4e3af7b0",
    "update_components": "671811cbcbee1329f4775a5dc0464b50",
    "update_components_series": "43e31765c602db13a19cda09d1029246",
}


@pytest.mark.parametrize("case", sorted(KEPT))
def test_valid_results_unchanged(case):
    """The checks added by PR 10 leave the results of the valid calls as
    they were, to the bit."""
    sessions, u, lya = _ids()
    result = _quiet(lambda: KEPT[case](sessions, u, lya))
    assert _fingerprint(result) == FINGERPRINTS[case]


def test_detect_anchor_on_a_normalized_flux_unchanged():
    # The mock has a continuum at 1: y and y/cont give the same anchor
    sessions, u, lya = _ids()
    fitted = sessions["fitted"]
    z = fitted.absorbers.detect_anchor(fitted, "CIV_1548", "CIV_1550")
    assert z == pytest.approx(DETECT_ANCHOR_Z, abs=0, rel=0)


DETECT_ANCHOR_Z = 1.930000147166357  # master, 12 November 2026


# --- Kept: the old logs ---------------------------------------------------

def test_update_components_series_none_string_is_unchanged():
    # 'None' of the old logs is "unchanged", as before (series was not in
    # the registry, a string); with another parameter, the same result
    sessions, u, lya = _ids()
    fitted = sessions["fitted"]
    old = _quiet(lambda: fitted.absorbers.update_components(lya, b=20.0, series="None"))
    assert _fingerprint(old) == FINGERPRINTS["update_components"]


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


def test_link_without_target_takes_that_of_the_expression():
    # F146: without target_uuid the target was not stored, and the fits
    # left it out of the group of the component; now it is the component
    # the expression names, as the GUI gives it
    sessions, u, lya = _ids()
    fitted = sessions["fitted"]
    new = fitted.absorbers.update_constraint(lya[1], "b", expression=f"p['{lya[0]}'].b")
    assert _fingerprint(new) == FINGERPRINTS["link_gui"]
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
    assert z == pytest.approx(DETECT_ANCHOR_Z, abs=0, rel=0)


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
