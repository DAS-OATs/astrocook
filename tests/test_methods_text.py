"""The log of a session as the text of a methods section (phase 6-ter, PR 5c:
astrocook.methods_text, decision 5 of the phase)."""
import string
from types import MappingProxyType

import pytest
from astropy.table import Table

from astrocook.core.atomic_data import ATOMIC_DATA_PATH
from astrocook.core.session import SessionV2
from astrocook.core.session_manager import SessionHistory
from astrocook.core.structures import HistoryLogV2
from astrocook.methods_text import (
    REFERENCES, MethodsValues, fmt, methods_file_name, methods_text, series_list,
    series_name)
from astrocook.mocks import MockComponent, MockTruth, generate
from astrocook.recipes import registry
from astrocook.recipes.registry import all_recipes, build_spec
from astrocook.recipes.result import run

Z = 1.9


def _mock():
    """A small mock: one CIV doublet at z = 1.9, z_em = 2.0, R = 45000."""
    truth = MockTruth(components=[MockComponent("CIV", Z, 13.8, 10.0)],
                      lam_min_nm=440.0, lam_max_nm=460.0, dv_kms=2.5,
                      resol=45000.0, snr=50.0, z_em=2.0, seed=1)
    return SessionV2("mock", None, spec=generate(truth))


def step(history, recipe, reason=None, **params):
    """Run a recipe on the current state and keep it, as the GUI and the MCP
    server do: the log entry, then the state."""
    result = run(history.current_state, recipe, **params)
    history.log_manager.add_entry(result.recipe.split('.')[1], dict(result.params),
                                  reason=reason)
    history.add_state(result.session)
    return result


@pytest.fixture(scope="module")
def history():
    """Continuum, a note, a CIV component, a refit; a component added and
    undone, with a note written after it (undone with it)."""
    h = SessionHistory(_mock(), HistoryLogV2())
    step(h, "continuum.estimate_auto", reason="To normalize the flux")
    h.log_manager.add_note("The continuum was checked by eye.")
    h.log_manager.add_note("A gap the server saw.", kind="gap_seen")
    step(h, "absorbers.add_component", reason="The doublet at 440.9 nm.",
         series="CIV", z=Z, logN=13.5, b=10.0)
    step(h, "absorbers.refit_all")
    step(h, "absorbers.add_component", series="SiIV", z=Z)
    h.log_manager.add_note("Written after the step undone.")
    h.undo()
    return h


def paragraphs(h, **kw):
    return [p.text for p in methods_text(h, **kw).paragraphs]


def test_one_paragraph_per_recipe(history):
    out = methods_text(history)
    kinds = [p.kind for p in out.paragraphs]
    assert kinds == ["opening", "step", "step", "step", "closing"]
    assert [p.recipe for p in out.paragraphs if p.kind == "step"] == [
        "continuum.estimate_auto", "absorbers.add_component", "absorbers.refit_all"]
    opening = out.paragraphs[0].text
    assert opening.startswith("The spectrum of mock (")
    assert "pixels, 440-460 nm) was analysed with Astrocook" in opening
    assert "(Cupani et al. 2020) in the three steps described below." in opening
    assert out.text == "\n\n".join(p.text for p in out.paragraphs) + \
        "\n\nReferences\n" + "\n".join(out.references)
    assert str(out) == out.text


def test_reason_text_changes_and_notes(history):
    _, cont, comp, refit, closing = paragraphs(history)
    # The reason first, as written (with a full stop), then the recipe with
    # its menu and its parameters, then what it changed, then the notes
    assert cont.startswith("To normalize the flux. The continuum was estimated "
                           "automatically (Continuum > Auto-estimate continuum).")
    assert "kappa = 2;" in cont and "Cupani et al. 2026" in cont
    # (the mock has its true continuum already)
    assert "This added the mask of the absorbed pixels to the spectrum. This changed " \
        "the continuum." in cont
    assert cont.endswith("The continuum was checked by eye.")
    # The gaps seen by the server are for the development: left out
    assert "A gap the server saw" not in "\n".join(paragraphs(history))
    assert comp.startswith("The doublet at 440.9 nm. A C IV component was added at "
                           "z = 1.9 (starting from log N = 13.5, b = 10 km/s)")
    assert "This added a C IV component at z = 1.9" in comp
    # Three components or fewer changed: each listed
    assert refit.startswith("All the components were refitted group by group "
                            "(Absorbers > Refit all systems).")
    assert refit.endswith("This changed the C IV component at z = 1.90000 by less "
                          "than the digits quoted.")
    assert "SciPy (Virtanen et al. 2020)" in closing
    assert "from Morton (2003), except the rest wavelengths of C IV 1548, 1550 " \
        "(Petitjean & Aracil 2004)" in closing


def test_the_recipes_undone_are_left_out(history):
    text = methods_text(history).text
    assert "Si IV" not in text and "Written after the step undone" not in text


def test_references(history):
    refs = methods_text(history).references
    assert refs == tuple(sorted(refs))
    for key in ("C20", "C26", "SciPy", "M03", "PA04"):
        assert REFERENCES[key][1] in refs
    assert REFERENCES["emcee"][1] not in refs


def test_without_the_menus(history):
    text = methods_text(history, menus=False).text
    assert " > " not in text
    assert "The continuum was estimated automatically. First," in text


def test_a_session_alone_has_only_its_log(history):
    """A session without its history (e.g. an archive just opened): the
    parameters, not what each recipe changed."""
    session = history.current_state
    text = methods_text(session).text
    assert "For the first three steps, read from the log of the archive, only " \
        "the parameters are known" in text
    assert "This added" not in text
    assert "The continuum was estimated automatically" in text


def test_an_archive_then_new_steps(history):
    """The entries of an archive have no states; those made after it do."""
    log = HistoryLogV2()
    log.add_entry("estimate_auto", dict(history.log_manager.entries[0].params))
    session = history.states[1].with_new_spectrum(history.states[1].spec)
    h = SessionHistory(session, log)
    step(h, "absorbers.add_component", series="CIV", z=Z)
    out = paragraphs(h)
    assert "For the first step, read from the log of the archive" in out[0]
    assert "This added" not in out[1]
    assert "This added a C IV component" in out[2]


def test_an_empty_log():
    out = methods_text(SessionHistory(_mock(), HistoryLogV2()))
    assert [p.kind for p in out.paragraphs] == ["opening"]
    assert "was opened with Astrocook" in out.text and "its log has no recipes" in out.text
    assert out.references == (REFERENCES["C20"][1],)


def test_notes_before_the_first_recipe_and_gaps():
    h = SessionHistory(_mock(), HistoryLogV2())
    h.log_manager.add_note("z_em from the Ly-a peak")
    h.log_manager.add_note("Needed z_em: measured by hand.", kind="gap")
    out = methods_text(h).paragraphs
    assert out[1].kind == "note"
    assert out[1].text == "z_em from the Ly-a peak. Needed z_em: measured by hand."


def test_an_unknown_recipe():
    """A recipe of a later version (not in the registry): its name and
    parameters."""
    session = _mock()
    session.log_manager = HistoryLogV2()
    session.log_manager.add_entry("future_recipe", {"a": 1.0, "b": "x"})
    assert "The recipe future_recipe was run (a = 1, b = x)." in methods_text(session).text


@pytest.fixture
def new_recipe(monkeypatch):
    """Two recipes added to the registry: one with its text, one without."""
    flux = dict(registry.RECIPES["flux"])
    params = [{"name": "factor", "type": float, "default": 2.0, "doc": "Factor."},
              {"name": "column", "type": str, "default": "y_scaled", "doc": "Column."},
              {"name": "hidden", "type": bool, "default": False, "gui_hidden": True}]
    flux["scale_test"] = build_spec("flux", "scale_test", {
        "brief": "Scale the flux (test).", "params": params,
        "methods": "The flux was multiplied by {factor}{menu}", "references": ("C26",)})
    flux["bare_test"] = build_spec("flux", "bare_test", {
        "brief": "Do something (test).", "params": params})
    recipes = dict(registry.RECIPES)
    recipes["flux"] = MappingProxyType(flux)
    monkeypatch.setattr(registry, "RECIPES", MappingProxyType(recipes))


def test_a_new_recipe_is_described_with_no_code_here(new_recipe):
    session = _mock()
    session.log_manager = HistoryLogV2()
    session.log_manager.add_entry("scale_test", {"factor": 3.0, "column": "y3",
                                                 "hidden": True})
    session.log_manager.add_entry("bare_test", {"factor": 2.0, "column": "y4"})
    out = methods_text(session)
    scale, bare = (p.text for p in out.paragraphs if p.kind == "step")
    # Its template, its menu, and the parameters it does not quote that
    # differ from their defaults (not the hidden ones)
    assert scale == ("The flux was multiplied by 3 (Flux > Scale the flux). "
                     "Other parameters: column = y3.")
    assert REFERENCES["C26"][1] in out.references
    # Without a text: a generic sentence
    assert bare == "Do something (test) was run (Flux > Do something). Parameters: column = y4."


@pytest.mark.parametrize("spec", list(all_recipes()), ids=lambda s: f"{s.category}.{s.name}")
def test_every_recipe_has_its_text(spec):
    """Every recipe of the registry says what it did (a recipe added later
    too), with its defaults, names its menu, and cites known works."""
    assert spec.methods is not None
    params = spec.defaults()
    p = MethodsValues(params, menu=" (MENU)")
    if callable(spec.methods):
        text = spec.methods(p)
    else:
        text = string.Formatter().vformat(spec.methods, (), p)
    assert text and "{" not in text
    assert "(MENU)" in text
    assert set(spec.references) <= set(REFERENCES)
    # The whole paragraph, from a log entry with the defaults
    session = _mock()
    session.log_manager = HistoryLogV2()
    session.log_manager.add_entry(spec.name, params)
    step_text = [x for x in methods_text(session).paragraphs if x.kind == "step"][0].text
    assert "Other parameters" not in step_text
    assert ("(" + spec.category.capitalize() + " > ") in step_text or spec.gui_hidden


def test_every_source_of_the_atomic_data_has_a_reference():
    table = Table.read(ATOMIC_DATA_PATH, format="ascii.ecsv")
    keys = set(table.meta["references"]) - {"-"}
    assert keys <= set(REFERENCES)


def test_values_and_names():
    assert [fmt(v) for v in (5000.0, 0.5, 2, True, None, [1.0, 2.5], "a")] == \
        ["5000", "0.5", "2", "yes", "none", "1, 2.5", "a"]
    assert series_name("CIV") == "C IV"
    assert series_name("CIV_1548") == "C IV 1548"
    assert series_name("Ly_a") == "Ly-alpha"
    assert series_name("FeII_2586,FeII_2600") == "Fe II 2586, Fe II 2600"
    assert series_list("CIV,SiIV,MgII") == "C IV, Si IV and Mg II"
    assert methods_file_name("J0103-1305 (2)") == "J0103-1305_2_methods.txt"


def _change(at, **fields):
    from astrocook.recipes.result import ComponentChange, ComponentRef
    return ComponentChange(ComponentRef(uuid=f"u{at}", series="CIV", z=at, logN=13.0,
                                        b=10.0), fields)


def test_what_changed_in_the_components():
    from types import SimpleNamespace
    from astrocook.methods_text import _changed_sentences, _flags_sentence
    # Up to three: each component, with the changes the digits show
    one = _change(2.0, logN=(13.0, 13.25), b=(10.0, 12.0), z=(2.0, 2.0001))
    assert _changed_sentences([one], 5) == [
        "This changed the C IV component at z = 2.00000: log N from 13.00 to 13.25, "
        "b from 10.0 to 12.0 km/s and z by 10.0 km/s."]
    # More: the median and the largest changes
    many = [_change(2.0 + i, logN=(13.0, 13.0 + 0.1 * i), b=(10.0, 10.0 + i))
            for i in range(1, 5)]
    assert _changed_sentences(many, 4) == [
        "The parameters of all 4 components changed.",
        "The median changes were 0.25 dex in log N and 2.5 km/s in b; the largest "
        "were 0.4 dex and 4 km/s."]
    assert _changed_sentences(many, 9)[0] == "The parameters of 4 components changed."
    # The flags the fit left
    comps = [SimpleNamespace(uuid="a", flags=("b_at_lower_bound",)),
             SimpleNamespace(uuid="b", flags=("b_at_lower_bound", "z_at_window")),
             SimpleNamespace(uuid="c", flags=("b_at_lower_bound",))]
    after = SimpleNamespace(systs=SimpleNamespace(components=comps))
    assert _flags_sentence(after, ["a", "b"]) == (
        "At the end, 2 components have b at the lower bound of the fit and 1 has z at the "
        "edge of the window of the fit: a parameter at a bound has an unreliable value "
        "and error.")
    assert _flags_sentence(after, []) is None
