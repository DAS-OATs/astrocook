"""The Details panel of the System Inspector (phase 7, PR 7d.9).

Below the table, shown or hidden by the box "Details" above it, it follows the selection of the
table (not the cursor) and says what the tooltip of a cell could not hold,
in three sections, each shown only when it has something to say:

1. the *group*: the statistics of the last fit of the selected component
   (:class:`~astrocook.core.structures.FitRecord`: reduced chi2, chi2,
   AIC, BIC, pixels, free parameters, convergence, the recipe and the
   time), and the components fitted with it, which select them when
   clicked;
2. the *component*: its values with their errors, the quantities of the
   criteria of its flags (b_inst and R at the strongest line, tau_0 of the
   weakest, the telluric depth), so that a flag can be checked, and each
   flag in one line;
3. the *parameter*, when a cell of a parameter is the current one: its
   state, the formula of its link with the names of the components, the
   bounds of its last fit.

The upper limits (7d.7) have a section of their own: for a row of a limit,
the limit alone (:func:`limit_html`); for a component, after the others,
the limits at its redshift (:func:`limits_here_html`). The plot of the
tabulated LSF (7d.8) will add its own: :meth:`DetailsPanel.sections` is
the one place that lists them.

The HTML of each section is built by a function of the data only, so that
the tests read it without a window.
"""
from __future__ import annotations

import html
from datetime import datetime
from typing import Callable, Dict, List, Optional, Sequence

import numpy as np
from PySide6.QtCore import QUrl, Signal
from PySide6.QtWidgets import QSizePolicy, QTextBrowser, QVBoxLayout, QWidget

from ..core.structures import ComponentData, FitRecord
from ..core.upper_limits import UpperLimit
from ..fitting.component_flags import (
    SATURATED_TAU, TELLURIC_DEPTH, UNRESOLVED_RATIO, FitPixels, flag_quantities)
from .inspector_table import C_KMS, dz_kms, flag_cue, pretty_expression

#: The width of the panel (pixels) from which the sections are in two columns
TWO_COLUMNS_WIDTH = 640

#: The text of the panel when no component is selected
NOTHING_SELECTED = "Select a component in the table to see its details."

#: The anchor of a member of the group, followed by its uuid
MEMBER_SCHEME = "component"

_PARAMS = ("z", "logN", "b", "btur")
_UNITS = {"z": "", "logN": "", "b": " km/s", "btur": " km/s"}


def _num(value: Optional[float], fmt: str = "{:.2f}") -> str:
    if value is None or not np.isfinite(value):
        return "n/a"
    return fmt.format(value)


def _row(label: str, value: str) -> str:
    return f"<tr><td style='padding-right: 10px'>{label}</td><td>{value}</td></tr>"


def _table(rows: Sequence[str]) -> str:
    return "<table cellspacing='0' cellpadding='1'>" + "".join(rows) + "</table>"


def _heading(text: str) -> str:
    return f"<p style='margin: 6px 0 2px 0'><b>{text}</b></p>"


def component_name(comp: ComponentData) -> str:
    """How a component is named in the panel: its series, z and ID."""
    return f"{html.escape(comp.series)} at z = {comp.z:.5f} (ID {comp.id})"


def _local_time(iso: str) -> str:
    try:
        return datetime.fromisoformat(iso).astimezone().strftime("%Y-%m-%d %H:%M:%S")
    except (TypeError, ValueError):
        return html.escape(iso or "unknown")


# --- 1. The group ----------------------------------------------------------------

def group_html(comp: ComponentData, components: Sequence[ComponentData]) -> str:
    """The statistics of the last fit of ``comp`` and its members; empty if
    it has none and no reduced chi2."""
    rec: Optional[FitRecord] = comp.fit
    if rec is None:
        if comp.chi2 is None:
            return ""
        # An archive saved before 7d.9 has only the reduced chi2
        return (_heading("Group fit")
                + _table([_row("Reduced χ²", _num(comp.chi2))])
                + "<p style='margin: 2px 0'>The other statistics of the fit were not "
                  "stored: refit the component to see them.</p>")
    if rec.converged is None:
        conv = "unknown"
    elif rec.converged:
        conv = f"yes ({rec.nfev} evaluations)" if rec.nfev is not None else "yes"
    else:
        conv = (f"<b>no</b>: stopped at {rec.max_nfev} evaluations, the values may not "
                "be at the minimum" if rec.max_nfev is not None else "<b>no</b>")
    rows = [
        _row("Reduced χ²", f"{_num(rec.chi2_reduced)} "
             f"({rec.n_constraining} constraining pixels)"),
        _row("χ²", f"{_num(rec.chi2)} on {rec.n_pix} pixels"),
        _row("AIC / BIC", f"{_num(rec.aic)} / {_num(rec.bic)}"),
        _row("Free parameters", str(rec.n_free)),
        _row("Converged", conv),
        # A break allowed after the dot: the column can be narrow
        _row("Recipe", html.escape(rec.recipe or "none (a script)").replace(".", ".\u200b")),
        _row("When", _local_time(rec.time)),
    ]
    out = _heading("Group fit") + _table(rows)
    if rec.final:
        out += ("<p style='margin: 2px 0'>The statistics at the final values, computed "
                "again after the neighbours were fitted (Refit All).</p>")
    by_uuid = {c.uuid: c for c in components}
    members = []
    for uuid in rec.members:
        other = by_uuid.get(uuid)
        if other is None:
            members.append("a component since removed")
        elif uuid == comp.uuid:
            members.append(f"<b>{component_name(other)}</b>")
        else:
            later = "" if other.fit is rec else " (fitted again since)"
            members.append(f"<a href='{MEMBER_SCHEME}:{uuid}'>{component_name(other)}</a>{later}")
    if members:
        out += (f"<p style='margin: 4px 0 2px 0'>Fitted together ({len(members)}):</p>"
                + "<br>".join(members))
    return out


# --- 2. The component ------------------------------------------------------------

def _value(comp: ComponentData, param: str, free: bool) -> str:
    val = getattr(comp, param)
    err = getattr(comp, "d" + param, None)
    if param == "z":
        text = f"{val:.6f}"
        if free and err is not None and np.isfinite(err) and err > 0:
            text += f" ± {err:.2e} ({dz_kms(comp.z, err):.2f} km/s)"
        return text
    fmt = "{:.3f}" if param == "logN" else "{:.2f}"
    text = fmt.format(val)
    if free and err is not None and np.isfinite(err) and err > 0:
        text += " ± " + fmt.format(err)
    elif free and err is not None and not np.isfinite(err):
        text += " ± NaN"
    return text + _UNITS[param]


def component_html(comp: ComponentData, pixels: Optional[FitPixels],
                   constraints: Dict[str, object]) -> str:
    """The values of ``comp``, the quantities of the criteria of its flags
    (on ``pixels``, the whole spectrum) and its flags, one per line."""
    rows = []
    for param in _PARAMS:
        if param == "btur" and not comp.btur and "btur" not in constraints:
            continue
        c = constraints.get(param)
        free = c is None or (getattr(c, "is_free", True)
                             and not getattr(c, "target_uuid", None)
                             and not getattr(c, "expression", None))
        rows.append(_row(param, _value(comp, param, free)))
    rows.append(_row("R", f"its own, {comp.resol:.0f}" if comp.resol else "that of the spectrum"))
    out = _heading(f"Component: {component_name(comp)}") + _table(rows)

    q = None
    if pixels is not None:
        q, _ = flag_quantities(comp, pixels)
    if q is not None:
        b_tot = float(np.hypot(comp.b, comp.btur or 0.0))
        at = html.escape(q.strongest or 'the strongest line')
        if q.b_inst > 0:
            resol = _row(f"R at {at}",
                         f"{_num(q.resol, '{:.0f}')}, b<sub>inst</sub> = {_num(q.b_inst)} km/s "
                         f"(unresolved below {UNRESOLVED_RATIO:g} b<sub>inst</sub> = "
                         f"{_num(UNRESOLVED_RATIO * q.b_inst)} km/s; b = {b_tot:.2f} km/s)")
        else:
            resol = _row(f"R at {at}", "not known: never unresolved")
        crit = [
            resol,
            _row(f"τ<sub>0</sub> of {html.escape(q.weakest or 'the weakest line')}",
                 f"{_num(q.tau0_weakest, '{:.3g}')} (saturated above {SATURATED_TAU:g})"),
        ]
        if q.telluric_depth is not None:
            crit.append(_row("Telluric depth", f"{_num(q.telluric_depth, '{:.3f}')} "
                             f"(on a telluric line above {TELLURIC_DEPTH:g})"))
        out += (_heading("The criteria of the flags")
                + "<p style='margin: 0 0 2px 0'>On the lines in the spectrum (the fit: "
                  "those in its pixels).</p>" + _table(crit))
    _, messages = flag_cue(comp)
    if messages:
        out += _heading("Flags") + "<br>".join(f"⚠ {html.escape(m)}" for m in messages)
    # 7d.8: the plot of the LSF at the strongest line goes here
    return out


# --- 3. The parameter ------------------------------------------------------------

def parameter_html(comp: ComponentData, param: Optional[str], constraint,
                   components: Sequence[ComponentData]) -> str:
    """The state of the parameter ``param`` of ``comp`` (None: no cell of a
    parameter selected), the formula of its link and its bounds."""
    if param not in _PARAMS:
        return ""
    rows = []
    target = getattr(constraint, "target_uuid", None)
    expression = getattr(constraint, "expression", None)
    if target or expression:
        rows.append(_row("State", "linked"))
        if target:
            ref = next((c for c in components if c.uuid == target), None)
            rows.append(_row("To", component_name(ref) if ref else "a component since removed"))
        if expression:
            pretty = pretty_expression(expression, components, name=component_name)
            rows.append(_row("Formula", f"<code>{pretty}</code>"))
    elif constraint is not None and not getattr(constraint, "is_free", True):
        rows.append(_row("State", "frozen: not fitted"))
    else:
        rows.append(_row("State", "free"))
    rec = comp.fit
    bounds = rec.bounds_of(comp.uuid, param) if rec is not None else None
    if bounds is not None:
        lo, hi = bounds
        if param == "z":
            width = (hi - lo) / (1 + comp.z) * C_KMS
            text = f"{lo:.6f} … {hi:.6f} (the window of the fit, {width:.1f} km/s)"
        else:
            text = f"{lo:g} … {hi:g}{_UNITS[param]}"
        rows.append(_row("Bounds of the fit", text))
    _, messages = flag_cue(comp, param)
    out = _heading(f"Parameter: {param}") + _table(rows)
    if messages:
        out += "<br>".join(f"⚠ {html.escape(m)}" for m in messages)
    return out


# --- 4. The upper limits (7d.7) ---------------------------------------------------

#: A limit is shown with a component when it took its redshift, or is within
#: this velocity of it (km/s)
LIMIT_HERE_KMS = 10.0

_METHOD_WORDS = {
    "profile": ("the χ² profile on fixed pixels: a component with z and b fixed, the "
                "components with a line on the pixels refitted at each log N; the limit "
                "where χ² rises by {s:g}² above its minimum"),
    "ew": ("the noise of the equivalent width within one FWHM of each line, on the "
           "linear part of the curve of growth; {s:g} times the error of the column "
           "(what absorbs there is not accounted for)"),
}


def limit_html(lim: UpperLimit, components: Sequence[ComponentData]) -> str:
    """The details of an upper limit: how it was computed, on which lines
    and pixels, how far from a detection, the components refitted."""
    b_text = f"{lim.b:.2f} km/s, fixed"
    if lim.b_inst > 0 and abs(lim.b - lim.b_inst) < 1e-6:
        b_text += " (the instrumental b at the strongest line)"
    elif lim.b_inst > 0:
        b_text += f" (b<sub>inst</sub> = {lim.b_inst:.2f} km/s)"
    detected = lim.significance > lim.sigma
    if lim.method == "profile":
        best = ("below the scan" if lim.logN_best is None
                else f"log N = {lim.logN_best:.2f}")
    else:
        best = ("none" if lim.logN_best is None
                else f"log N = {lim.logN_best:.2f} (the equivalent width measured)")
    sig = f"{lim.significance:.1f} σ (best: {best})"
    if detected:
        sig = f"<b>{sig}</b>: detected, the limit is not meaningful"
    ranges = ", ".join(f"{lo:.2f}–{hi:.2f}" for lo, hi in lim.ranges_nm[:4])
    if len(lim.ranges_nm) > 4:
        ranges += f" and {len(lim.ranges_nm) - 4} more"
    rows = [
        _row("log N", f"&lt; {lim.logN_lim:.2f} ({lim.sigma:g} σ)"),
        _row("z", f"{lim.z:.6f}"),
        _row("b", b_text),
        _row("Method", html.escape(lim.method) + ": "
             + _METHOD_WORDS.get(lim.method, "").format(s=lim.sigma)),
        _row("Lines", html.escape(", ".join(lim.transitions))),
    ]
    if lim.skipped:
        rows.append(_row("Left out", html.escape(", ".join(lim.skipped))
                         + " (outside the spectrum or without valid pixels)"))
    rows += [
        _row("Pixels", f"{lim.n_pix} ({ranges} nm)"),
        _row("Over no absorption", sig),
        _row("Recipe", html.escape(lim.recipe or "none (a script)").replace(".", ".\u200b")),
        _row("When", _local_time(lim.time)),
    ]
    out = _heading(f"Upper limit: {html.escape(lim.series)}") + _table(rows)
    by_uuid = {c.uuid: c for c in components}
    named = []
    if lim.ref_uuid and lim.ref_uuid in by_uuid:
        named.append(("Redshift of", [lim.ref_uuid]))
    if lim.refitted:
        named.append(("Refitted at each step", list(lim.refitted)))
    for label, uuids in named:
        links = [f"<a href='{MEMBER_SCHEME}:{u}'>{component_name(by_uuid[u])}</a>"
                 if u in by_uuid else "a component since removed" for u in uuids]
        out += f"<p style='margin: 4px 0 2px 0'>{label}:</p>" + "<br>".join(links)
    out += ("<p style='margin: 4px 0 2px 0; color: gray'>Not a component: outside the "
            "model. Removed when the data or the components on its pixels change.</p>")
    return out


def limits_here_html(comp: ComponentData, limits: Sequence[UpperLimit]) -> str:
    """The upper limits at the redshift of ``comp``, one line each; empty if
    none."""
    here = [lim for lim in limits
            if lim.ref_uuid == comp.uuid
            or abs(lim.z - comp.z) / (1 + comp.z) * C_KMS < LIMIT_HERE_KMS]
    if not here:
        return ""
    lines = [f"{html.escape(lim.series)}: log N &lt; {lim.logN_lim:.2f} "
             f"({lim.sigma:g} σ, {html.escape(lim.method)})" for lim in here]
    return _heading("Upper limits at this redshift") + "<br>".join(lines)


# --- The widget ------------------------------------------------------------------

class DetailsPanel(QWidget):
    """The Details panel: the text of the sections of the selected
    component. Shown or hidden by the box "Details" of the inspector,
    above the table (no header of its own: it took a row and a wide
    button, Guido, 27 November 2026)."""

    #: The uuid of a member of the group, clicked
    component_clicked = Signal(str)

    def __init__(self, parent: Optional[QWidget] = None):
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        self.browser = QTextBrowser()
        self.browser.setOpenLinks(False)
        self.browser.anchorClicked.connect(self._on_anchor)
        self.browser.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        layout.addWidget(self.browser, 1)
        self._comp: Optional[ComponentData] = None
        self._limit: Optional[UpperLimit] = None
        self._parts: List[str] = []
        self._has_group = False
        self._expanded = True
        self.browser.setHtml(f"<p style='color: gray'>{NOTHING_SELECTED}</p>")

    def set_expanded(self, expanded: bool) -> None:
        """Show or hide the panel (hidden, the splitter gives its height to
        the table)."""
        self._expanded = bool(expanded)
        self.setVisible(self._expanded)

    def is_expanded(self) -> bool:
        return self._expanded

    def _on_anchor(self, url: QUrl) -> None:
        text = url.toString()
        if text.startswith(MEMBER_SCHEME + ":"):
            self.component_clicked.emit(text[len(MEMBER_SCHEME) + 1:])

    @staticmethod
    def sections(comp: ComponentData, components: Sequence[ComponentData],
                 constraints: Dict[str, object], pixels: Optional[FitPixels],
                 param: Optional[str], limits: Sequence[UpperLimit] = ()) -> List[str]:
        """The HTML of the sections that have something to say, in order:
        group, component, parameter, the upper limits at its redshift."""
        parts: List[Callable[[], str]] = [
            lambda: group_html(comp, components),
            lambda: component_html(comp, pixels, constraints),
            lambda: parameter_html(comp, param, constraints.get(param) if param else None,
                                   components),
            lambda: limits_here_html(comp, limits),
        ]
        return [s for s in (p() for p in parts) if s]

    def show_component(self, comp: Optional[ComponentData],
                       components: Sequence[ComponentData] = (),
                       constraints: Optional[Dict[str, object]] = None,
                       pixels: Optional[FitPixels] = None,
                       param: Optional[str] = None,
                       limits: Sequence[UpperLimit] = ()) -> None:
        """Show the details of ``comp`` (None: nothing selected)."""
        self._comp = comp
        self._limit = None
        if comp is None:
            self._parts = []
            self.browser.setHtml(f"<p style='color: gray'>{NOTHING_SELECTED}</p>")
            return
        self._parts = self.sections(comp, components, constraints or {}, pixels, param,
                                    limits)
        self._has_group = comp.fit is not None or comp.chi2 is not None
        self._render()

    def show_limit(self, lim: UpperLimit,
                   components: Sequence[ComponentData] = ()) -> None:
        """Show the details of an upper limit (a row of the table, 7d.7)."""
        self._comp = None
        self._limit = lim
        self._parts = [limit_html(lim, components)]
        self._has_group = False
        self._render()

    def _two_columns(self) -> bool:
        return self.browser.viewport().width() >= TWO_COLUMNS_WIDTH

    def _render(self) -> None:
        """Lay the sections out: in two columns when the panel is wide (the
        group on the left, the component and the parameter on the right),
        else one after the other."""
        parts = self._parts
        self._columns = self._two_columns()
        if self._columns and self._has_group and len(parts) > 1:
            cells = [parts[0], "".join(parts[1:])]
            row = "".join(f"<td valign='top' style='padding-right: 18px' width='50%'>{c}</td>"
                          for c in cells)
            body = f"<table width='100%'><tr>{row}</tr></table>"
        else:
            body = "<hr>".join(parts)
        self.browser.setHtml(f"<div style='font-size: 12px'>{body}</div>")

    def resizeEvent(self, event) -> None:
        super().resizeEvent(event)
        # Laid out again only when the number of columns changes
        if ((self._comp is not None or self._limit is not None) and self._parts
                and self._two_columns() != getattr(self, '_columns', None)):
            self._render()

    def text(self) -> str:
        """The text shown (for the tests)."""
        return self.browser.toPlainText()


__all__ = ["DetailsPanel", "group_html", "component_html", "parameter_html",
           "limit_html", "limits_here_html", "component_name", "MEMBER_SCHEME"]

