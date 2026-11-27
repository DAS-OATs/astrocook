"""The upper limits on the column of a species not detected (phase 7, PR 7d.7;
F156): the computation (profile likelihood and equivalent width), the recipes
``upper_limit`` and ``delete_upper_limit``, where they are stored (the
metadata of the system list) and saved (the archive), their removal when the
data or the components on their pixels change, the MCP output."""
import dataclasses
import functools
import json
import logging
import warnings

import numpy as np
import pytest

from astrocook.core.structures import SystemListData
from astrocook.core.system_list import SystemList
from astrocook.core.upper_limits import (
    META_KEY, SCHEMA_VERSION, UpperLimit, from_file_dict, limits_of, to_file_dict, with_limits)
from astrocook.fitting.upper_limits import (
    COG_LINEAR, components_fingerprint, stale_reason, upper_limit)
from astrocook.mocks import MockComponent, MockTruth, generate
from astrocook.recipes.errors import InvalidParameter, RecipeFailed
from astrocook.recipes.result import run
from astrocook.session.session import Session, load_session_from_file

C_KMS = 299792.458
#: A fitted CIV, a weak SiIV at its redshift (below the noise), a SiIV
#: detected at another redshift; R 45000, S/N 50, the true continuum
TRUTH = MockTruth(components=[MockComponent("CIV", 1.93, 13.8, 10.0),
                              MockComponent("SiIV", 1.93, 11.0, 8.0),
                              MockComponent("SiIV", 1.95, 12.8, 8.0)],
                  lam_min_nm=400.0, lam_max_nm=460.0, dv_kms=2.5, resol=45000.0,
                  snr=50.0, z_em=2.7, seed=1)


@functools.lru_cache(maxsize=1)
def _base():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        s = Session("mock", None, spec=generate(TRUTH, truth_columns=False))
        s = s.continuum.estimate_auto()
        s = s.absorbers.add_component(series="CIV", z=1.93, logN=13.8, b=10.0)
    return s


@pytest.fixture
def base():
    return _base()


def _civ(session):
    return next(c for c in session.systs.components if c.series == "CIV")


# --- The structure ---------------------------------------------------------------

def _limit(**kw):
    d = dict(series="SiIV", z=1.93, b=4.0, logN_lim=11.3, sigma=3.0, method="profile",
             transitions=("SiIV_1393", "SiIV_1402"), ranges_nm=((408.1, 408.5),),
             reach_nm=((408.0, 408.6),), refitted=("a",), n_pix=30)
    d.update(kw)
    return UpperLimit(**d)


def test_the_structure_reads_back():
    lim = _limit(logN_best=10.9, skipped=("SiIV_1402",))
    assert UpperLimit.from_dict(json.loads(json.dumps(lim.to_dict()))) == lim
    assert from_file_dict(json.loads(json.dumps(to_file_dict((lim,))))) == (lim,)
    # A key of a later version is ignored; a later schema refused
    assert UpperLimit.from_dict({**lim.to_dict(), "new_key": 1}) == lim
    with pytest.raises(ValueError, match="version"):
        from_file_dict({**to_file_dict((lim,)), "schema_version": SCHEMA_VERSION + 1})
    with pytest.raises(ValueError, match="not the upper limits"):
        from_file_dict({"schema": "other"})
    assert lim.label() == "< 11.30"


def test_the_limits_live_in_the_metadata_of_the_system_list():
    systs = SystemList(SystemListData())
    assert limits_of(systs) == () and limits_of(None) == ()
    lim = _limit()
    new = with_limits(systs, (lim,))
    assert limits_of(new) == (lim,)
    assert META_KEY not in systs._data.meta  # the old list unchanged


# --- The computation -------------------------------------------------------------

def test_ew_is_the_noise_of_the_equivalent_width(base):
    """The limit of 'ew' from its definition, on the pixels it reports."""
    spec, z = base.spec, 1.93
    lim = upper_limit(spec, base.systs, "SiIV_1393", z, method="ew")
    x = spec.x.to_value("Angstrom")
    dl = (spec.xmax - spec.xmin).to_value("Angstrom")
    norm = spec.norm.value
    lo, hi = lim.ranges_nm[0]
    sel = (x >= lo * 10 * (1 - 1e-12)) & (x <= hi * 10 * (1 + 1e-12))
    sigma_w = np.sqrt(np.sum((spec.dy.value[sel] / norm[sel] * dl[sel]) ** 2)) / (1 + z)
    f, lam = 0.5280, 1393.76
    from astrocook.core.atomic_data import ATOM_DATA
    f, lam = ATOM_DATA["SiIV_1393"]["f"], ATOM_DATA["SiIV_1393"]["wave"]
    expected = np.log10(3 * COG_LINEAR * sigma_w / (f * lam ** 2))
    assert lim.logN_lim == pytest.approx(expected, abs=1e-9)
    assert lim.n_pix == int(sel.sum()) and lim.method == "ew" and lim.refitted == ()
    # b = b_inst at R 45000
    assert lim.b == pytest.approx(C_KMS / (2 * np.sqrt(np.log(2)) * 45000))


def test_profile_and_ew_agree_on_a_clean_region(base):
    """Where nothing absorbs (CIV at z 1.6, 402-403 nm), the chi2 profile
    and the noise of the equivalent width give the same limit to about
    0.1 dex."""
    prof = upper_limit(base.spec, base.systs, "CIV", 1.6)
    ew = upper_limit(base.spec, base.systs, "CIV", 1.6, method="ew")
    assert prof.logN_lim == pytest.approx(ew.logN_lim, abs=0.1)
    assert prof.significance < 3 and prof.refitted == ()
    assert prof.transitions == ("CIV_1548", "CIV_1550") and prof.skipped == ()


def test_a_weak_line_is_below_the_limit(base):
    """The weak SiIV of the mock (log N 11.0) is below the limit; the
    minimum of the chi2 near it raises the limit over the noise alone."""
    prof = upper_limit(base.spec, base.systs, "SiIV", 1.93)
    ew = upper_limit(base.spec, base.systs, "SiIV", 1.93, method="ew")
    assert 11.0 < ew.logN_lim < prof.logN_lim < 11.8
    assert prof.significance < 3


def test_the_limit_grows_with_the_confidence_and_with_b(base):
    one = upper_limit(base.spec, base.systs, "SiIV", 1.93, sigma=1.0)
    three = upper_limit(base.spec, base.systs, "SiIV", 1.93, sigma=3.0)
    wide = upper_limit(base.spec, base.systs, "SiIV", 1.93, b=30.0)
    assert one.logN_lim < three.logN_lim < wide.logN_lim


def test_a_species_detected_says_so(base):
    lim = upper_limit(base.spec, base.systs, "SiIV", 1.95, b=8.0)
    assert lim.significance > 10
    # The minimum near the truth, the limit just above it
    assert lim.logN_best == pytest.approx(12.8, abs=0.1)
    assert lim.logN_best < lim.logN_lim < lim.logN_best + 0.2


def test_a_blend_is_refitted(base):
    """A CIV 20 km/s from the fitted one: the fitted CIV is refitted at
    each step, so the limit is not set by its wing."""
    civ = _civ(base)
    prof = upper_limit(base.spec, base.systs, "CIV", 1.93 + 20 / C_KMS * 2.93)
    ew = upper_limit(base.spec, base.systs, "CIV", 1.93 + 20 / C_KMS * 2.93, method="ew")
    assert prof.refitted == (civ.uuid,)
    # The noise alone sees the wing of the CIV as absorption
    assert ew.significance > 3
    assert np.isfinite(prof.logN_lim) and prof.logN_lim < 13.8


def test_the_lines_outside_the_spectrum_are_left_out(base):
    # MgII at z 0.46: 2796 at 408.3 nm in the spectrum, 2803 at 409.3 too;
    # CIV at z 1.6: 1548 at 402.5, 1550 at 403.2 in; at z 1.97 1550 at 460.6 out
    lim = upper_limit(base.spec, base.systs, "CIV", 1.968, method="ew")
    assert lim.transitions == ("CIV_1548",) and lim.skipped == ("CIV_1550",)
    with pytest.raises(ValueError, match="falls in the spectrum"):
        upper_limit(base.spec, base.systs, "NV", 1.93)


def test_the_model_does_not_change(base):
    """A limit is outside the model: the components and the model column
    are those before."""
    out = base.absorbers.upper_limit(series="SiIV", z=1.93)
    assert out.systs.components == base.systs.components
    np.testing.assert_array_equal(out.spec.get_column("model").value,
                                  base.spec.get_column("model").value)


# --- The recipes -------------------------------------------------------------------

def test_the_recipe_stores_the_limit_and_reports_it(base):
    civ = _civ(base)
    result = run(base, "absorbers.upper_limit", series="SiIV", uuid=civ.uuid)
    (lim,) = limits_of(result.session.systs)
    assert lim.ref_uuid == civ.uuid and lim.z == civ.z and lim.series == "SiIV"
    assert lim.recipe == "absorbers.upper_limit"
    reported = result.extra["upper_limit"]
    assert reported["logN_lim"] == lim.logN_lim and reported["uuid"] == lim.uuid
    assert "components" not in reported and "data" not in reported
    assert result.warnings == ()


def test_the_recipe_refuses(base):
    civ = _civ(base)
    for kw, match in (({"series": "", "z": 1.93}, "series"),
                      ({"series": "XX_1", "z": 1.93}, "series"),
                      ({"series": "SiIV"}, "give one of them"),
                      ({"series": "SiIV", "z": 1.93, "uuid": civ.uuid}, "give one of them"),
                      ({"series": "SiIV", "uuid": "nope"}, "uuid"),
                      ({"series": "SiIV", "z": 1.93, "sigma": 0}, "sigma"),
                      ({"series": "SiIV", "z": 1.93, "method": "fit"}, "method")):
        with pytest.raises(InvalidParameter, match=match):
            base.absorbers.upper_limit(**kw)
    with pytest.raises(RecipeFailed, match="falls in the spectrum"):
        base.absorbers.upper_limit(series="NV", z=1.93)


def test_without_the_resolution_b_is_needed(base):
    from astrocook.core.structures import DataColumn
    import astropy.units as au
    spec = base.spec
    nan = DataColumn(np.full(len(spec.x), np.nan), au.dimensionless_unscaled)
    no_r = spec.__class__(dataclasses.replace(
        spec._data, aux_cols={**spec._data.aux_cols, "resol": nan}, resol=0.0))
    with pytest.raises(ValueError, match="give b"):
        upper_limit(no_r, base.systs, "SiIV", 1.93)
    assert np.isfinite(upper_limit(no_r, base.systs, "SiIV", 1.93, b=5.0).logN_lim)


def test_warnings_of_a_detection(base, caplog):
    result = run(base, "absorbers.upper_limit", series="SiIV", z=1.95, b=8.0)
    assert any("is detected at" in w for w in result.warnings)
    result = run(base, "absorbers.upper_limit", series="SiIV", z=1.95, b=8.0, method="ew")
    assert any("a sensitivity, not a limit" in w for w in result.warnings)


def test_a_limit_at_the_same_redshift_is_replaced_and_deleted(base):
    s = base.absorbers.upper_limit(series="SiIV", z=1.93, method="ew")
    s = s.absorbers.upper_limit(series="SiIV", z=1.93 + 0.5 / C_KMS * 2.93)
    s = s.absorbers.upper_limit(series="SiIV_1393", z=1.93, method="ew")
    limits = limits_of(s.systs)
    assert [(lim.series, lim.method) for lim in limits] == [("SiIV", "profile"),
                                                            ("SiIV_1393", "ew")]
    s = s.absorbers.delete_upper_limit(uuid=limits[0].uuid)
    assert limits_of(s.systs) == (limits[1],)
    with pytest.raises(InvalidParameter, match="no upper limit"):
        s.absorbers.delete_upper_limit(uuid=limits[0].uuid)


def test_methods_text():
    from astrocook.methods_text import MethodsValues
    from astrocook.recipes.absorbers import _methods_delete_upper_limit, _methods_upper_limit
    defaults = {"series": "SiIV", "z": None, "uuid": None, "b": None, "sigma": 3.0,
                "method": "profile"}
    text = _methods_upper_limit(MethodsValues({**defaults, "z": 1.93}, menu=" (MENU)"))
    assert text == (
        "An upper limit on the column of Si IV at z = 1.93 was computed with b equal to the "
        "instrumental b at its strongest line, by the profile likelihood on fixed pixels: a "
        "component with z and b fixed, the components with a line on the pixels refitted at "
        "each log N, the limit where the chi2 rises by 3^2 above its minimum (MENU)")
    text = _methods_upper_limit(MethodsValues(
        {**defaults, "uuid": "u", "b": 8.0, "method": "ew", "sigma": 2.0}))
    assert "the redshift of the selected component" in text and "b = 8 km/s" in text
    assert "equivalent width" in text and "at 2 sigma" in text
    assert _methods_delete_upper_limit(MethodsValues({"uuid": "u"})) == "An upper limit was deleted"


# --- Removed when the state changes (Guido, 27 November 2026) ---------------------

def test_a_change_on_the_pixels_removes_the_limit(base):
    civ = _civ(base)
    near = 1.93 + 20 / C_KMS * 2.93
    s = base.absorbers.upper_limit(series="CIV", z=near)
    s = s.absorbers.upper_limit(series="SiIV", z=1.93, method="ew")
    result = run(s, "absorbers.update_component", uuid=civ.uuid, b=12.0)
    # The CIV limit refitted the CIV: removed; the SiIV, far, kept
    (kept,) = limits_of(result.session.systs)
    assert kept.series == "SiIV"
    assert any("upper limit on CIV" in w and "the components on its pixels changed" in w
               for w in result.warnings)
    (removed,) = result.extra["upper_limits_removed"]
    assert removed["series"] == "CIV" and removed["reason"].startswith("the components")


def test_a_refit_on_the_pixels_removes_the_limit(base):
    civ = _civ(base)
    edited = base.absorbers.update_component(uuid=civ.uuid, b=12.0)
    s = edited.absorbers.upper_limit(series="CIV", z=1.93 + 20 / C_KMS * 2.93)
    assert len(limits_of(s.systs)) == 1
    s = s.absorbers.fit_component(uuid=civ.uuid)
    assert limits_of(s.systs) == ()


def test_a_component_added_on_the_pixels_removes_the_limit(base):
    s = base.absorbers.upper_limit(series="SiIV", z=1.93)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        s = s.absorbers.add_component(series="SiIV", z=1.93, logN=11.5, b=8.0)
    assert limits_of(s.systs) == ()


def test_deleting_the_component_of_the_redshift_removes_the_limit(base):
    civ = _civ(base)
    s = base.absorbers.upper_limit(series="SiIV", uuid=civ.uuid, method="ew")
    assert stale_reason(s.spec, s.systs, limits_of(s.systs)[0]) is None
    s = s.absorbers.delete_component(uuid=civ.uuid)
    assert limits_of(s.systs) == ()


def test_a_new_continuum_removes_the_limit(base):
    s = base.absorbers.upper_limit(series="SiIV", z=1.93, method="ew")
    result = run(s, "continuum.estimate_auto", smooth_len_lya=2000.0)
    assert limits_of(result.session.systs) == ()
    assert any("the data on its pixels changed" in w for w in result.warnings)


def test_a_change_elsewhere_keeps_the_limit(base):
    s = base.absorbers.upper_limit(series="SiIV", z=1.93)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        s = s.absorbers.add_component(series="CIV", z=1.95, logN=12.0, b=10.0)
    assert len(limits_of(s.systs)) == 1
    # A list built anew with the same components keeps them, checked
    same = SystemList(SystemListData(components=list(s.systs.components)))
    assert len(limits_of(s.with_new_system_list(same).systs)) == 1


def test_a_fingerprint_follows_the_constraints(base):
    civ = _civ(base)
    before = components_fingerprint(base.systs, [civ.uuid])
    frozen = base.systs.update_constraint(civ.uuid, "b", is_free=False)
    assert components_fingerprint(frozen, [civ.uuid]) != before
    # The defaults written out are the defaults
    explicit = base.systs.update_constraint(civ.uuid, "b", is_free=True)
    assert components_fingerprint(explicit, [civ.uuid]) == before


# --- The archive ---------------------------------------------------------------------

def test_the_limits_are_saved_and_read_back(base, tmp_path):
    s = base.absorbers.upper_limit(series="CIV", z=1.93 + 20 / C_KMS * 2.93)
    s = s.absorbers.upper_limit(series="SiIV", z=1.93, method="ew")
    path = tmp_path / "ul.acs2"
    s.save(str(path))
    back = load_session_from_file(archive_path=str(path), name="ul", gui_context=None,
                                  format_name="acs2")
    assert limits_of(back.systs) == limits_of(s.systs)
    # Read back, they still describe the state
    for lim in limits_of(back.systs):
        assert stale_reason(back.spec, back.systs, lim) is None


def test_an_archive_without_limits(base, tmp_path):
    import tarfile
    path = tmp_path / "none.acs2"
    base.save(str(path))
    with tarfile.open(path) as tar:
        assert not any("upper_limits" in n for n in tar.getnames())


# --- MCP -------------------------------------------------------------------------------

def test_components_list_gives_the_limits(base, tmp_path):
    pytest.importorskip("mcp")
    from test_mcp_inspect import calls, ok
    from astrocook.mcp_server.server import build_server
    path = tmp_path / "m.acs2"
    s = base.absorbers.upper_limit(series="SiIV", z=1.93, method="ew")
    s.save(str(path))
    server = build_server()
    _, listed, only_civ = calls(
        server, ("session_open", {"path": str(path), "name": "m"}),
        ("components_list", {"session_id": "m"}),
        ("components_list", {"session_id": "m", "series": "CIV"}))
    (info,) = ok(listed)["upper_limits"]
    assert info["series"] == "SiIV" and info["method"] == "ew"
    assert info["logN_lim"] == pytest.approx(limits_of(s.systs)[0].logN_lim)
    assert ok(only_civ)["upper_limits"] == []


def test_the_mcp_tool_computes_a_limit(base, tmp_path):
    pytest.importorskip("mcp")
    from test_mcp_inspect import calls, ok
    from astrocook.mcp_server.server import build_server
    path = tmp_path / "m.acs2"
    base.save(str(path))
    server = build_server()
    _, out = calls(server, ("session_open", {"path": str(path), "name": "m"}),
                   ("absorbers_upper_limit", {"session_id": "m", "series": "SiIV",
                                              "z": 1.93, "method": "ew"}))
    assert ok(out)
    logging.getLogger(__name__).debug(out)
