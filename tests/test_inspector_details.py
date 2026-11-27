"""The Details panel of the System Inspector (phase 7, PR 7d.9): the group
fit, the component with the quantities of the criteria of its flags, the
parameter of the current cell; it follows the selection, its members select
their component, it collapses."""
import dataclasses

import numpy as np
import pytest
from PySide6.QtCore import QItemSelectionModel, QUrl

from astrocook.core.structures import ComponentData, FitRecord, ParameterConstraint
from astrocook.fitting.component_flags import FitPixels
from astrocook.gui.inspector_details import (
    MEMBER_SCHEME, DetailsPanel, component_html, group_html, parameter_html)

from test_inspector_flags import _gui


def _record(members, **kw):
    base = dict(chi2=120.5, n_pix=110, n_free=6, aic=132.5, bic=148.7, chi2_reduced=1.12,
                n_constraining=95, members=tuple(members), converged=True, nfev=17,
                max_nfev=200, recipe="absorbers.refit_all", time="2026-11-27T10:00:00+00:00",
                bounds=((members[0], "z", 2.99980, 3.00020), (members[0], "b", 1.0, 50.0)))
    base.update(kw)
    return FitRecord(**base)


def _pair(**kw):
    a = ComponentData(id=1, z=3.0, dz=1e-5, logN=13.5, dlogN=0.05, b=12.0, db=1.0,
                      series="CIV", chi2=1.12)
    b = ComponentData(id=2, z=3.0001, dz=2e-5, logN=13.0, dlogN=0.1, b=8.0, db=1.5,
                      series="CIV", chi2=1.12)
    rec = _record([a.uuid, b.uuid], **kw)
    return dataclasses.replace(a, fit=rec), dataclasses.replace(b, fit=rec)


def _text(html_):
    from PySide6.QtGui import QTextDocument
    doc = QTextDocument()
    doc.setHtml(html_)
    return doc.toPlainText()


# --- The sections ------------------------------------------------------------------

def test_group(qapp):
    a, b = _pair()
    text = _text(group_html(a, [a, b]))
    for part in ("Reduced χ²", "1.12", "95 constraining pixels", "120.50 on 110 pixels",
                 "132.50 / 148.70", "Free parameters", "6", "yes (17 evaluations)",
                 "refit_all", "Fitted together (2)", "CIV at z = 3.00010 (ID 2)"):
        assert part in text, part
    # The other member is a link that selects it
    assert f"href='{MEMBER_SCHEME}:{b.uuid}'" in group_html(a, [a, b])
    assert f"{MEMBER_SCHEME}:{a.uuid}" not in group_html(a, [a, b])


def test_group_not_converged_and_final(qapp):
    a, b = _pair(converged=False, final=True)
    text = _text(group_html(a, [a, b]))
    assert "no: stopped at 200 evaluations" in text
    assert "computed again after the neighbours were fitted" in text


def test_a_member_fitted_again_since(qapp):
    a, b = _pair()
    b = dataclasses.replace(b, fit=_record([b.uuid]))
    assert "(fitted again since)" in _text(group_html(a, [a, b]))
    assert "a component since removed" in _text(group_html(a, [a]))


def test_group_without_a_record(qapp):
    a, _ = _pair()
    old = dataclasses.replace(a, fit=None)
    assert "were not stored" in _text(group_html(old, [old]))
    assert group_html(dataclasses.replace(old, chi2=None), [old]) == ""


def test_component_and_the_criteria_of_its_flags(qapp):
    comp = ComponentData(id=3, z=3.0, dz=1e-5, logN=15.2, dlogN=0.3, b=1.5, db=0.4,
                         series="CIV", flags=("unresolved", "saturated"))
    wave = np.linspace(6150.0, 6250.0, 4000)
    pixels = FitPixels(wave=wave, resol=np.full(wave.size, 45000.0),
                       telluric=np.where(np.abs(wave - 6193.0) < 1.0, 0.9, 1.0))
    text = _text(component_html(comp, pixels, {}))
    for part in ("15.200 ± 0.300", "1.50 ± 0.40 km/s", "R at CIV_1548", "45000",
                 "inst = 4.00 km/s", "unresolved below 0.5", "τ0 of CIV_1550",
                 "saturated above 3", "Telluric depth", "0.100", "Unresolved:", "Saturated:"):
        assert part in text, part
    # Without the spectrum, the values and the flags only
    assert "criteria" not in _text(component_html(comp, None, {}))


def test_the_error_of_a_frozen_or_linked_parameter_is_not_shown(qapp):
    a, _ = _pair()
    text = _text(component_html(a, None, {"b": ParameterConstraint(is_free=False)}))
    assert "12.00 km/s" in text and "12.00 ± 1.00" not in text


def test_parameter(qapp):
    a, b = _pair()
    text = _text(parameter_html(a, "z", None, [a, b]))
    assert "free" in text and "the window of the fit, 30.0 km/s" in text
    assert "1 … 50 km/s" in _text(parameter_html(a, "b", None, [a, b]))
    link = ParameterConstraint(is_free=False, target_uuid=a.uuid,
                               expression=f"p['{a.uuid}'].b * 0.5")
    text = _text(parameter_html(b, "b", link, [a, b]))
    assert "linked" in text and "CIV at z = 3.00000 (ID 1).b * 0.5" in text
    assert "frozen" in _text(parameter_html(a, "logN", ParameterConstraint(is_free=False),
                                            [a, b]))
    # Not a parameter
    assert parameter_html(a, "series", None, [a, b]) == ""
    assert parameter_html(a, None, None, [a, b]) == ""


def test_sections_say_only_what_they_have(qapp):
    comp = ComponentData(id=1, z=3.0, dz=None, logN=13.0, dlogN=None, b=10.0, db=None,
                         series="CIV")
    parts = DetailsPanel.sections(comp, [comp], {}, None, None)
    assert len(parts) == 1 and "Component" in parts[0]


# --- The widget in the inspector -----------------------------------------------------

@pytest.fixture
def inspector(qapp):
    from astrocook.core.structures import SystemListData
    from astrocook.core.system_list import SystemList
    from astrocook.mocks import MockComponent, MockTruth, generate
    from astrocook.session.session import Session
    a, b = _pair()
    c = ComponentData(id=3, z=3.005, dz=None, logN=13.0, dlogN=None, b=10.0, db=None,
                      series="SiIV")
    systs = SystemList(SystemListData(components=[a, b, c]))
    spec = generate(MockTruth(
        components=[MockComponent("CIV", 3.0, 13.5, 12.0)], lam_min_nm=555.0,
        lam_max_nm=625.0, dv_kms=2.5, resol=45000.0, snr=30.0, z_em=3.5, seed=1))
    insp = _gui().SystemInspector(None)
    insp.current_session = Session("s", None, spec=spec, systs=systs)
    insp.table_model.update_data(systs.components, {})
    insp.details_cb.setChecked(True)
    return insp


def _select(insp, row, column="ID"):
    index = insp.proxy_model.index(row, _gui().COLUMNS.index(column))
    insp.table_view.selectionModel().setCurrentIndex(
        index, QItemSelectionModel.ClearAndSelect | QItemSelectionModel.Rows)


def test_the_panel_follows_the_selection(inspector):
    _select(inspector, 0)
    text = inspector.details.text()
    assert "Group fit" in text and "ID 1" in text and "Parameter" not in text
    _select(inspector, 0, "b")
    assert "Parameter: b" in inspector.details.text()
    _select(inspector, 2)
    text = inspector.details.text()
    assert "Group fit" not in text and "SiIV" in text
    # The criteria of the flags, on the spectrum of the session
    assert "R at SiIV_1393" in text


def test_a_member_clicked_is_selected(inspector):
    _select(inspector, 0)
    b = inspector.table_model.get_component_at(1)
    inspector.details._on_anchor(QUrl(f"{MEMBER_SCHEME}:{b.uuid}"))
    assert inspector._selected_component().uuid == b.uuid


def test_the_panel_collapses(inspector):
    inspector.details_cb.setChecked(False)
    assert not inspector.details.isVisibleTo(inspector)
    _select(inspector, 0)
    # Not updated while collapsed; updated when expanded again
    inspector.details_cb.setChecked(True)
    assert "ID 1" in inspector.details.text()
