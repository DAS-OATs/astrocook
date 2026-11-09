"""The log of a session as the text of a methods section (phase 6-ter, PR 5c).

:func:`methods_text` writes what was done to a session as English text, to
paste in the methods section of a paper (decision 5 of phase 6-ter): an
opening paragraph, one paragraph per recipe of the log, a closing paragraph
on the fit and the atomic data, and the references.

A paragraph of a recipe is made of:

- the reason given for the recipe (``LogEntryV2.reason``, e.g. by an agent
  through the MCP server), as written;
- what the recipe did, from the key ``methods`` of the registry (a template
  of its parameters, or a function of them), with the menu of the recipe in
  parentheses (``menus=False`` leaves it out), and the parameters the text
  does not quote that differ from their defaults;
- what the recipe changed, from the comparison of the state before and
  after it (:func:`~astrocook.recipes.result.compare_sessions`, as
  ``RecipeResult`` does): the components added, removed and changed, the
  columns, the pixels; and the warnings the fit left on the components;
- the notes written at its state (``LogNoteV2``: notes and gaps declared,
  as written; the gaps seen by the MCP server are for the development and
  are left out).

The recipes undone are left out. The text is made from the log each time it
is asked for, never stored: a better wording reaches the sessions saved
before it. Nothing here names a recipe: a recipe added to the registry with
its ``methods`` is described with no code here, one without it with a
generic sentence.
"""
from __future__ import annotations

import os
import re
import string
from collections import Counter
from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

#: The works the texts cite: key -> (the citation in the text, the entry of
#: the list of references). The keys of the atomic data are those of the
#: table (``astrocook/data/atomic_data.ecsv``, docs page *Atomic Data*); a
#: test checks that every key of the table is here.
REFERENCES: Dict[str, Tuple[str, str]] = {
    'C20': ("Cupani et al. 2020",
            "Cupani G. et al. 2020, Proc. SPIE 11452, 114521U"),
    'C26': ("Cupani et al. 2026",
            "Cupani G. et al. 2026, Proc. SPIE 14155, 141551X, doi:10.1117/12.3101756"),
    'SciPy': ("Virtanen et al. 2020",
              "Virtanen P. et al. 2020, Nature Methods 17, 261"),
    'emcee': ("Foreman-Mackey et al. 2013",
              "Foreman-Mackey D. et al. 2013, PASP 125, 306"),
    'dynesty': ("Speagle 2020", "Speagle J. S. 2020, MNRAS 493, 3132"),
    # The sources of the atomic data
    'M03': ("Morton 2003", "Morton D. C. 2003, ApJS 149, 205"),
    'K14': ("Kisielius et al. 2014", "Kisielius R. et al. 2014, ApJ 780, 76"),
    'K15': ("Kisielius et al. 2015", "Kisielius R. et al. 2015, ApJ 804, 76"),
    'BB19': ("Boisse & Bergeron 2019", "Boisse P., Bergeron J. 2019, A&A 622, A140"),
    'PA04': ("Petitjean & Aracil 2004", "Petitjean P., Aracil B. 2004, A&A 422, 523"),
    'VPFIT': ("the atom.dat file of VPFIT; Carswell & Webb 2014",
              "Carswell R. F., Webb J. K. 2014, Astrophysics Source Code "
              "Library, ascl:1408.015"),
    'V1': ("the table of Astrocook 1, original source not recorded", None),
    'V2': ("the table of Astrocook 2.0-2.1, original source not recorded", None),
}

#: The key of the work that describes Astrocook, cited by every text
SOFTWARE_REFERENCE = 'C20'
#: The key whose presence among the references of the recipes run says that
#: Voigt profiles were fitted (the closing paragraph says how)
FIT_REFERENCE = 'SciPy'

#: How the text names the columns of the spectrum; None: not named (they
#: follow from what the text says, e.g. the resolution of Set properties)
COLUMN_NAMES: Dict[str, Optional[str]] = {
    'x': "the wavelengths", 'y': "the flux", 'dy': "the error of the flux",
    'cont': "the continuum", 'abs_mask': "the mask of the absorbed pixels",
    'model': "the model", 'cont_pl': "the power-law continuum (cont_pl)",
    'resol': None, 'abs_ids': None, 'fit_mask': None,
}

#: How the text names the flags the fit leaves on the components
#: (``voigt_fitter.COMPONENT_FLAGS``, and ``not_fitted`` of Refit all)
FLAG_NAMES: Dict[str, str] = {
    'z_at_window': "z at the edge of the window of the fit",
    'logN_at_lower_bound': "log N at the lower bound of the fit",
    'logN_at_upper_bound': "log N at the upper bound of the fit",
    'b_at_lower_bound': "b at the lower bound of the fit",
    'b_at_upper_bound': "b at the upper bound of the fit",
    'nan_errors': "no uncertainties (they could not be computed)",
    'not_fitted': "the values from before the fit (the fit of the group failed)",
}

#: The kinds of notes left out of the text
NOTES_LEFT_OUT = ('gap_seen',)

_C_KMS = 299792.458
_QUANTITIES = (('wave', "rest wavelengths"), ('f', "oscillator strengths"),
               ('gamma', "damping constants"))


# --- The values of the parameters, as the templates read them ---------------

def fmt(value: Any) -> str:
    """A value of a parameter as the text writes it: ``5000.0`` as
    ``5000``, a list as its items separated by commas."""
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, (float, np.floating)):
        value = float(value)
        if value.is_integer() and abs(value) < 1e15:
            return str(int(value))
        return f"{value:g}"
    if isinstance(value, (list, tuple)):
        return ", ".join(fmt(v) for v in value)
    if value is None:
        return "none"
    return str(value)


def given(value: Any) -> bool:
    """The value says something: not None, empty, 'None' or '_current_'."""
    if value is None:
        return False
    if isinstance(value, str):
        return value.strip() not in ("", "None", "none", "_current_", "auto")
    if isinstance(value, (list, tuple)):
        return len(value) > 0
    return True


class MethodsValues:
    """The parameters of a log entry, as the ``methods`` of the registry read
    them: ``p['kappa']`` is the value as text (:func:`fmt`), ``p.raw('kappa')``
    the value as logged, ``p.menu`` the menu of the recipe in parentheses
    (with a space before it) or "". The parameters read are remembered: the
    others that differ from their defaults are listed after the text."""

    def __init__(self, params: Dict[str, Any], menu: str = ""):
        self._params = dict(params)
        self.menu = menu
        self.used: set = set()

    def raw(self, name: str) -> Any:
        self.used.add(name)
        return self._params.get(name)

    def given(self, name: str) -> bool:
        """The parameter has a value that says something (:func:`given`)."""
        return given(self.raw(name))

    def __getitem__(self, name: str) -> str:
        if name == 'menu':
            return self.menu
        return fmt(self.raw(name))

    def __contains__(self, name: str) -> bool:
        return name in self._params


# --- Words shared by the texts of the recipes -------------------------------

#: The regions of the spectrum, as the parameters name them
REGION_WORDS = {
    'red_side': "redward of the Lyman-alpha emission",
    'lya_forest': "in the Lyman-alpha forest",
    'all_ly_forest': "in the Lyman forest (from the Lyman limit to Lyman-alpha)",
}


def region_words(region: Any, before: str = " ") -> str:
    """`` redward of the Lyman-alpha emission`` for ``'red_side'``; "" for
    no region."""
    if not given(region):
        return ""
    words = [REGION_WORDS.get(r.strip(), f"in the region {r.strip()}")
             for r in str(region).split(',') if r.strip()]
    return before + _join(words)


def fit_bounds(p: MethodsValues) -> str:
    """``b <= 50 km/s, log N <= 19 and the redshift free within +-20 km/s``,
    from the parameters ``max_b``, ``max_logN`` and ``z_window_kms`` that the
    recipe has."""
    parts = []
    if 'max_b' in p:
        parts.append(f"b <= {p['max_b']} km/s")
    if 'max_logN' in p:
        parts.append(f"log N <= {p['max_logN']}")
    if 'z_window_kms' in p:
        parts.append(f"the redshift free within +-{p['z_window_kms']} km/s")
    return _join(parts)


def mask_words(p: MethodsValues) -> str:
    """The mask of the absorbed pixels that a recipe used (its parameter
    ``mask``, and ``mask_col`` if it has one)."""
    if 'mask' in p and p.raw('mask') == 'new':
        return ("a new mask of the absorbed pixels (made with the defaults "
                "of Find absorbed regions)")
    column = p['mask_col'] if 'mask_col' in p else "abs_mask"
    return f"the mask of the absorbed pixels of the session ({column})"


def file_name(p: MethodsValues, name: str) -> str:
    """The name of a file, without its folder (which may say who the user
    is, and means nothing to the reader)."""
    value = p.raw(name)
    return os.path.basename(str(value)) if given(value) else "a file"


# --- The text ------------------------------------------------------------

@dataclass(frozen=True)
class MethodsParagraph:
    """A paragraph of the text."""
    text: str
    #: ``opening``, ``step`` (a recipe of the log), ``note`` (notes written
    #: before the first recipe) or ``closing``
    kind: str
    #: The recipe of a step, as ``category.name``
    recipe: Optional[str] = None
    #: The index of the log entry of a step
    entry: Optional[int] = None


@dataclass(frozen=True)
class MethodsText:
    """The text of the methods of a session (:func:`methods_text`)."""
    paragraphs: Tuple[MethodsParagraph, ...]
    #: The list of references, one entry each, sorted
    references: Tuple[str, ...]

    @property
    def text(self) -> str:
        """The paragraphs separated by blank lines, then the references."""
        out = "\n\n".join(p.text for p in self.paragraphs)
        if self.references:
            out += "\n\nReferences\n" + "\n".join(self.references)
        return out

    def __str__(self) -> str:
        return self.text


def _join(items: Sequence[str]) -> str:
    items = list(items)
    if len(items) <= 1:
        return "".join(items)
    return ", ".join(items[:-1]) + " and " + items[-1]


def _sentence(text: str) -> str:
    """The text with a capital letter and a full stop."""
    text = _as_written(text)
    return text[:1].upper() + text[1:]


def _as_written(text: str) -> str:
    """A text written by who ran the analysis (a reason, a note), with a
    full stop if it has none, and otherwise as written."""
    text = text.strip()
    if text and text[-1] not in ".!?":
        text += "."
    return text


_NUMBERS = ("no", "one", "two", "three", "four", "five", "six", "seven", "eight",
            "nine", "ten")


def _number(n: int) -> str:
    return _NUMBERS[n] if 0 <= n < len(_NUMBERS) else str(n)


_ION = re.compile(r"^([A-Z][a-z]?)([IVX]+)(\*{0,2})$")
_LYMAN = {'Ly_a': "Ly-alpha", 'Ly_b': "Ly-beta", 'Ly_ab': "Ly-alpha and Ly-beta",
          'Ly': "H I (Lyman series)"}


def series_list(text: str) -> str:
    """A list of series separated by commas (a parameter, e.g.
    ``'CIV,SiIV,MgII'``) as ``C IV, Si IV and Mg II``."""
    return _join([series_name(s.strip()) for s in str(text).split(',') if s.strip()])


def series_name(series: str) -> str:
    """A series as a paper writes it: ``CIV`` as ``C IV``, ``CIV_1548`` as
    ``C IV 1548``, ``Ly_a`` as ``Ly-alpha``."""
    if ',' in series:
        return ", ".join(series_name(s.strip()) for s in series.split(',') if s.strip())
    if series in _LYMAN:
        return _LYMAN[series]
    ion, _, line = series.partition('_')
    m = _ION.match(ion)
    if m:
        ion = f"{m.group(1)} {m.group(2)}{m.group(3)}"
    if ion == 'Ly' and line:
        return f"Ly-{line}"
    return f"{ion} {line}".strip()


def _menu(spec) -> str:
    """`` (Continuum > Auto-estimate continuum)``, or "" for a recipe not in
    the menus."""
    if spec.gui_hidden:
        return ""
    brief = re.sub(r"\s*\([^)]*\)\s*$", "", spec.brief.rstrip('.').strip())
    return f" ({spec.category.capitalize()} > {brief})"


def _recipe_sentence(spec, params: Dict[str, Any], menus: bool) -> Tuple[str, List[str]]:
    """What the recipe did, and the keys of the works it cites."""
    p = MethodsValues(params, _menu(spec) if menus else "")
    if spec.methods is None:
        text = f"{spec.brief.rstrip('.')} was run{p.menu}"
    elif callable(spec.methods):
        text = spec.methods(p)
    else:
        text = string.Formatter().vformat(spec.methods, (), p)
    text = _sentence(text)
    others = []
    for ps in spec.params:
        if ps.name in p.used or ps.gui_hidden or not ps.agents or ps.name not in params:
            continue
        value = params[ps.name]
        if fmt(value) != fmt(ps.default):
            others.append(f"{ps.name} = {fmt(value)}")
    if others:
        label = "Other parameters" if spec.methods is not None else "Parameters"
        text += f" {label}: {', '.join(others)}."
    return text, list(spec.references)


def _unknown_sentence(name: str, params: Dict[str, Any]) -> str:
    values = ", ".join(f"{k} = {fmt(v)}" for k, v in params.items())
    return _sentence(f"The recipe {name} was run" + (f" ({values})" if values else ""))


# --- What a recipe changed ----------------------------------------------------

def _z_range(zs: Sequence[float]) -> str:
    if len(zs) == 1:
        return f"z = {zs[0]:.3f}"
    lo, hi = min(zs), max(zs)
    if f"{lo:.3f}" == f"{hi:.3f}":
        return f"z = {lo:.3f}"
    return f"{lo:.3f} < z < {hi:.3f}"


def _components_by_series(refs) -> str:
    by_series: Dict[str, List[float]] = {}
    for r in refs:
        by_series.setdefault(r.series, []).append(r.z)
    return _join([f"{len(zs)} {series_name(s)} ({_z_range(zs)})"
                  for s, zs in by_series.items()])


def _one_component(r) -> str:
    return (f"a {series_name(r.series)} component at z = {r.z:.5f} "
            f"(log N = {r.logN:.2f}, b = {r.b:.1f} km/s)")


def _small(v: float) -> str:
    if v == 0:
        return "0"
    if v < 0.001:
        return "< 0.001"
    return f"{v:.2g}"


def _changed_sentences(changed, n_after: int) -> List[str]:
    if not changed:
        return []
    if len(changed) <= 3:
        out = []
        for ch in changed:
            r = ch.component
            parts = []
            for field, (old, new) in ch.changed.items():
                if field == 'flags' or old is None or new is None:
                    continue
                if field == 'z':
                    dv = abs(new - old) / (1 + old) * _C_KMS
                    if dv >= 0.05:
                        parts.append(f"z by {dv:.1f} km/s")
                elif field == 'logN':
                    if f"{old:.2f}" != f"{new:.2f}":
                        parts.append(f"log N from {old:.2f} to {new:.2f}")
                elif field in ('b', 'btur'):
                    label = 'b' if field == 'b' else 'b_tur'
                    if f"{old:.1f}" != f"{new:.1f}":
                        parts.append(f"{label} from {old:.1f} to {new:.1f} km/s")
                elif field == 'series':
                    parts.append(f"the series from {series_name(old)} to {series_name(new)}")
                elif field == 'resol':
                    parts.append(f"the resolving power from {fmt(old)} to {fmt(new)}")
            where = f"the {series_name(r.series)} component at z = {r.z:.5f}"
            if parts:
                out.append(f"This changed {where}: {_join(parts)}.")
            elif any(f != 'flags' for f in ch.changed):
                out.append(f"This changed {where} by less than the digits quoted.")
        return out
    d = {'logN': [], 'b': [], 'v': []}
    for ch in changed:
        c = ch.changed
        d['logN'].append(abs(c['logN'][1] - c['logN'][0]) if 'logN' in c else 0.0)
        d['b'].append(abs(c['b'][1] - c['b'][0]) if 'b' in c else 0.0)
        d['v'].append(abs(c['z'][1] - c['z'][0]) / (1 + c['z'][0]) * _C_KMS
                      if 'z' in c else 0.0)
    n = len(changed)
    which = f"all {n}" if n == n_after else str(n)
    out = [f"The parameters of {which} components changed."]
    units = {'logN': (" dex in log N", " dex"), 'b': (" km/s in b", " km/s"),
             'v': (" km/s in velocity", " km/s")}
    keys = [k for k in ('logN', 'b', 'v') if any(d[k])]
    if keys:
        med = _join([f"{_small(float(np.median(d[k])))}{units[k][0]}" for k in keys])
        top = _join([f"{_small(float(np.max(d[k])))}{units[k][1]}" for k in keys])
        out.append(f"The median changes were {med}; the largest were {top}.")
    return out


def _flags_sentence(after, uuids: Iterable[str]) -> Optional[str]:
    uuids = set(uuids)
    systs = getattr(after, 'systs', None)
    if not uuids or systs is None:
        return None
    counts = Counter(f for c in systs.components if c.uuid in uuids
                     for f in (c.flags or ()))
    if not counts:
        return None
    parts = []
    for flag, n in counts.items():
        what = FLAG_NAMES.get(flag, flag.replace('_', ' '))
        noun = f" component{'s' if n != 1 else ''}" if not parts else ""
        parts.append(f"{n}{noun} {'have' if n != 1 else 'has'} {what}")
    tail = (": a parameter at a bound has an unreliable value and error"
            if any('bound' in f or 'window' in f for f in counts) else "")
    return f"At the end, {_join(parts)}{tail}."


def _columns_sentence(names: Sequence[str], verb: str, components_changed: bool) -> Optional[str]:
    words = []
    for n in names:
        if n == 'model' and components_changed:
            continue
        word = COLUMN_NAMES.get(n, f"the column '{n}'")
        if word:
            words.append(word)
    if not words:
        return None
    where = {"added": " to the spectrum", "removed": " from the spectrum", "changed": ""}[verb]
    return f"This {verb} {_join(words)}{where}."


def changes_text(before, after) -> List[str]:
    """What changed from the state ``before`` to ``after``, in sentences."""
    from astrocook.recipes.result import compare_sessions
    ch = compare_sessions(before, after)
    out = []
    comps = bool(ch.components_added or ch.components_removed or ch.components_changed)
    for refs, verb in ((ch.components_added, "added"), (ch.components_removed, "removed")):
        if len(refs) == 1:
            out.append(f"This {verb} {_one_component(refs[0])}.")
        elif refs:
            out.append(f"This {verb} {len(refs)} components: {_components_by_series(refs)}.")
    n_after = len(after.systs.components) if getattr(after, 'systs', None) else 0
    out += _changed_sentences(ch.components_changed, n_after)
    for names, verb in ((ch.columns_added, "added"), (ch.columns_removed, "removed"),
                        (ch.columns_changed, "changed")):
        s = _columns_sentence(names, verb, comps)
        if s:
            out.append(s)
    if 'pixels' in ch.properties_changed:
        old, new = ch.properties_changed['pixels']
        if old is not None and new is not None:
            out.append(f"The spectrum has now {new} pixels ({old} before).")
    flags = _flags_sentence(after, [r.uuid for r in ch.components_added]
                            + [c.component.uuid for c in ch.components_changed])
    if flags:
        out.append(flags)
    return out


# --- The atomic data ------------------------------------------------------------

def _observed_range_angstrom(spec) -> Optional[Tuple[float, float]]:
    import astropy.units as au
    try:
        x = spec.x.to(au.Angstrom).value
    except au.UnitConversionError:  # an x axis in velocity
        return None
    x = x[np.isfinite(x)]
    return (float(x.min()), float(x.max())) if len(x) else None


def atomic_data_sentence(state) -> Tuple[Optional[str], List[str]]:
    """Where the atomic data of the transitions of the components of
    ``state`` come from, and the keys of the works it cites. Only the
    transitions that fall in the spectrum are counted (all, if the x axis
    is not a wavelength)."""
    from astropy.table import Table
    from astrocook.core.atomic_data import ATOMIC_DATA_PATH, series_lines
    systs = getattr(state, 'systs', None)
    if systs is None or not systs.components:
        return None, []
    table = Table.read(ATOMIC_DATA_PATH, format='ascii.ecsv')
    rows = {str(r['name']): r for r in table}
    xr = _observed_range_angstrom(state.spec) if getattr(state, 'spec', None) else None
    used = []
    for c in systs.components:
        for t in series_lines(c.series):
            r = rows.get(t)
            if r is None or t in used:
                continue
            if xr is not None:
                obs = float(r['wave']) * (1 + c.z)
                if not xr[0] <= obs <= xr[1]:
                    continue
            used.append(t)
    if not used:
        return None, []
    sources = Counter(str(rows[t][f"{q}_ref"]) for t in used for q, _ in _QUANTITIES
                      if str(rows[t][f"{q}_ref"]) != '-')
    if not sources:
        return None, []
    main = sources.most_common(1)[0][0]
    exceptions = []
    keys = [main]
    for q, label in _QUANTITIES:
        by_key: Dict[str, List[str]] = {}
        for t in used:
            k = str(rows[t][f"{q}_ref"])
            if k not in ('-', main):
                by_key.setdefault(k, []).append(t)
        for k, ts in by_key.items():
            exceptions.append(f"the {label} of {_transitions(ts)} ({_cite(k)})")
            keys.append(k)
    text = f"The atomic data are those of the Astrocook table: from {_cite(main, True)}"
    if exceptions:
        text += ", except " + _join(exceptions)
    return text + ".", keys


def _transitions(names: Sequence[str]) -> str:
    """``CIV_1548, CIV_1550`` as ``C IV 1548, 1550``."""
    groups: Dict[str, List[str]] = {}
    for n in names:
        ion, _, line = n.partition('_')
        groups.setdefault(ion, []).append(line)
    out = []
    for ion, lines in groups.items():
        name = series_name(f"{ion}_{lines[0]}")
        out.append(", ".join([name] + lines[1:]))
    return _join(out)


def _cite(key: str, narrative: bool = False) -> str:
    short = REFERENCES.get(key, (key, None))[0]
    m = re.match(r"^(.*) (\d{4})$", short)
    if narrative and m and ';' not in short:
        return f"{m.group(1)} ({m.group(2)})"
    return short


# --- The log and the states -------------------------------------------------------

def _source(source):
    """The states, the index of the current one and the log of ``source``:
    a ``SessionHistory`` or a session (its log only, no states)."""
    from astrocook.core.structures import HistoryLogV2
    if hasattr(source, 'states') and hasattr(source, 'log_manager'):
        states, current, log = list(source.states), source.current_index, source.log_manager
    else:
        states, current, log = [source], 0, getattr(source, 'log_manager', None)
    if not isinstance(log, HistoryLogV2):
        log = None
    return states, current, log


def _spectrum_sentence(name: str, state, n_steps: int) -> str:
    from astrocook import __version__
    spec = getattr(state, 'spec', None)
    what = f"The spectrum of {name}"
    if spec is not None and len(spec.x):
        x = np.asarray(spec.x.value, dtype=float)
        x = x[np.isfinite(x)]
        if len(x):
            what += (f" ({len(spec.x)} pixels, {x.min():.0f}-{x.max():.0f} "
                     f"{spec.x.unit.to_string()})")
    software = f"Astrocook {__version__} ({_cite(SOFTWARE_REFERENCE)})"
    if n_steps == 0:
        return f"{what} was opened with {software}; its log has no recipes."
    steps = f"in the {_number(n_steps)} steps described below" if n_steps != 1 else \
        "in the step described below"
    return f"{what} was analysed with {software} {steps}."


def _name(source, state) -> str:
    spec = getattr(state, 'spec', None)
    obj = spec.meta.get('OBJECT') if spec is not None else None
    if isinstance(obj, str) and obj.strip() and not obj.strip().isdigit():
        return obj.strip()
    return getattr(source, 'display_name', None) or getattr(state, 'name', None) or "the object"


def methods_text(source, *, menus: bool = True) -> MethodsText:
    """The log of a session as the text of a methods section.

    Parameters
    ----------
    source : SessionHistory or SessionV2
        A history (the states and the log, as the GUI and the MCP server
        keep them), or a session alone: then only the log is known, not what
        each recipe changed.
    menus : bool
        Name the menu of each recipe, in parentheses (for who reruns the
        analysis in the GUI; not for a paper).

    Returns
    -------
    MethodsText
        The paragraphs and the references; ``.text`` (or ``str()``) is the
        text.
    """
    from astrocook.recipes.registry import find_recipe
    states, current, log = _source(source)
    timeline = log.timeline() if log is not None else []
    # The state made by the entry i is states[i - offset]; the entries up to
    # offset were made before the first state (e.g. read from an archive)
    offset = (log.current_index - current) if log is not None else 0
    # The log of the archive the session was opened from (F109): its steps
    # come first, made before the first state
    previous = getattr(log, 'previous', None)
    if previous is not None:
        earlier = previous.timeline()
        timeline = earlier + timeline
        offset += sum(1 for item, _ in earlier if hasattr(item, 'recipe_name'))
    first = states[0] if states else None
    entries = [item for item, _ in timeline if hasattr(item, 'recipe_name')]
    paragraphs: List[MethodsParagraph] = []
    keys: List[str] = [SOFTWARE_REFERENCE]

    opening = [_spectrum_sentence(_name(source, first), first, len(entries))]
    if entries:
        opening.append("Parameters not quoted were left at their default values. The "
                       "session log, saved with the session, replays the analysis "
                       "exactly.")
    before_states = sum(1 for i in range(len(entries)) if i - offset < 1)
    if before_states == 1:
        opening.append("For the first step, read from the log of the archive, only "
                       "the parameters are known, not what it changed.")
    elif before_states:
        opening.append(f"For the first {_number(before_states)} steps, read from the "
                       f"log of the archive, only the parameters are known, not what "
                       f"they changed.")
    paragraphs.append(MethodsParagraph(" ".join(opening), "opening"))

    index = -1
    for item, _ in timeline:
        if not hasattr(item, 'recipe_name'):
            if item.kind in NOTES_LEFT_OUT:
                continue
            note = _as_written(item.text)
            if paragraphs[-1].kind in ("step", "note"):
                last = paragraphs[-1]
                paragraphs[-1] = MethodsParagraph(f"{last.text} {note}", last.kind,
                                                  last.recipe, last.entry)
            else:
                paragraphs.append(MethodsParagraph(note, "note"))
            continue
        index += 1
        spec = find_recipe(item.recipe_name)
        parts = [_as_written(item.reason)] if item.reason else []
        if spec is None:
            parts.append(_unknown_sentence(item.recipe_name, item.params))
            recipe = item.recipe_name
        else:
            text, refs = _recipe_sentence(spec, item.params, menus)
            parts.append(text)
            keys += refs
            recipe = f"{spec.category}.{spec.name}"
        s = index - offset
        if 1 <= s < len(states):
            parts += changes_text(states[s - 1], states[s])
        paragraphs.append(MethodsParagraph(" ".join(parts), "step", recipe, index))

    final = states[current] if states else None
    closing = []
    has_components = bool(final is not None and getattr(final, 'systs', None)
                          and final.systs.components)
    if has_components and FIT_REFERENCE in keys:
        closing.append("The Voigt profiles were computed with the Faddeeva function "
                       "and fitted by least squares, both from SciPy "
                       f"({_cite(FIT_REFERENCE)}).")
    if has_components:
        atomic, atomic_keys = atomic_data_sentence(final)
        if atomic:
            closing.append(atomic)
            keys += atomic_keys
    if closing:
        paragraphs.append(MethodsParagraph(" ".join(closing), "closing"))

    references = sorted({REFERENCES[k][1] for k in keys
                         if k in REFERENCES and REFERENCES[k][1]})
    return MethodsText(paragraphs=tuple(paragraphs), references=tuple(references))


def methods_file_name(session_name: str) -> str:
    """A file name for the text of a session: ``<name>_methods.txt``."""
    stem = re.sub(r"[^\w.-]+", "_", session_name or "session").strip("_") or "session"
    return f"{stem}_methods.txt"


__all__ = ['MethodsParagraph', 'MethodsText', 'MethodsValues', 'REFERENCES',
           'atomic_data_sentence', 'changes_text', 'file_name', 'fit_bounds', 'fmt',
           'given', 'mask_words', 'series_list', 'methods_file_name', 'methods_text',
           'region_words', 'series_name']

