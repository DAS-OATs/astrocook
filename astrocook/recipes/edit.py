from astrocook._renamed import renamed
import contextlib
from itertools import pairwise
from astropy import units as au
import logging
import numpy as np
import re
from typing import TYPE_CHECKING, Dict, cast, Optional, Union, List, Tuple

from astrocook.core.atomic_data import ATOM_DATA
from astrocook.core.column_names import (display_columns, internal_column,
                                         internal_expression, is_axis_column,
                                         is_core_column)
from astrocook.core.spectrum_expressions import REGION_PRESETS
from astrocook.core.structures import SystemListData
from astrocook.core.units import (has_no_unit, is_arbitrary_flux,
                                  is_flux_density_per_wavelength)
from astrocook.core.system_list import SystemList
from astrocook.recipes.errors import (InvalidParameter, RecipeError, RecipeFailed,
                                      reason_of)
from astrocook.recipes._combine import in_units_of_first
from astrocook.recipes._sessions import is_session, open_session, open_sessions, split_names
from astrocook.recipes._spectrum import spec_of
from astrocook.recipes.params import DATA_ERRORS, typed_recipe
from astrocook.recipes.messages import msg_session_not_found

if TYPE_CHECKING:
    from astrocook.session.session import Session

logger = logging.getLogger(__name__)

# --- The text of the methods (phase 6-ter, PR 5c: astrocook.methods_text) ---
from astrocook.methods_text import file_name, fmt, region_words
def _methods_set_properties(p):
    parts = []
    if p.given('name'):
        parts.append(f"the name of the session to {p['name']}")
    if p.given('object'):
        parts.append(f"the name of the object to {p['object']}")
    if p.given('z_em'):
        parts.append(f"the emission redshift to z_em = {p['z_em']}")
    if p.given('resol'):
        resol = str(p.raw('resol')).strip()
        if resol.endswith('px'):
            parts.append(f"the resolution to a line spread function of FWHM {resol[:-2]} pixels")
        else:
            try:
                resol = fmt(float(resol))
            except ValueError:
                pass
            parts.append(f"the resolving power to R = {resol}")
    if p.given('yunit'):
        unit = str(p['yunit'])
        parts.append("the flux was declared in arbitrary units" if unit.strip() == 'arb_flux'
                     else f"the unit of the flux was declared as {unit}, the values unchanged")
    p.raw('overwrite_resol_col')
    if not parts:
        return f"The properties of the session were left as they were{p.menu}"
    first = parts[0].removeprefix("the ").replace(" to ", " was set to ", 1)
    text = "The " + _join_and([first] + parts[1:])
    return text + p.menu


def _join_and(items):
    return items[0] if len(items) == 1 else ", ".join(items[:-1]) + " and " + items[-1]


def _methods_x_convert(p):
    unit = p['xunit']
    try:
        velocity = au.Unit(unit).physical_type == 'speed'
    except ValueError:
        velocity = False
    if velocity:
        z = f"z = {p['z_ref']}" if p.given('z_ref') else "z_em"
        return f"The x axis was converted to velocity ({unit}) relative to {z}{p.menu}"
    p.raw('z_ref')
    return f"The x axis was converted to {unit}{p.menu}"


def _methods_import_systems(p):
    how = "added to those" if p.raw('append') else "put in place of those"
    refit = ", and refitted on its spectrum" if p.raw('refit') else ""
    return (f"The absorption systems of the session {p['source_session']} were "
            f"copied and {how} of this session{refit}{p.menu}")


def _methods_import_systems_ascii(p):
    how = "added to those of the session" if p.raw('append') else "put in place of those of the session"
    refit = ", and refitted" if p.raw('refit') else ""
    return (f"Absorption components were read from the table "
            f"{file_name(p, 'file_path')} and {how}{refit}{p.menu}")


def _methods_delete(p):
    items = [t.strip() for t in str(p.raw('targets') or '').split(',') if t.strip()]
    parts = []
    if 'systems' in items:
        parts.append("all the components")
    columns = [t for t in items if t != 'systems']
    if columns:
        parts.append(f"the column{'s' if len(columns) > 1 else ''} {', '.join(columns)}")
    what = _join_and(parts) if parts else "nothing"
    return f"{what[:1].upper()}{what[1:]} {'was' if what == 'nothing' else 'were'} deleted{p.menu}"


def _methods_extract_preset(p):
    regions = [r.strip() for r in str(p.raw('region')).split(',') if r.strip()]
    words = [region_words(r, before="") for r in regions]
    each = "each into a new session" if len(regions) > 1 else "into a new session"
    return f"The part of the spectrum {_join_and(words)} was copied {each}{p.menu}"


def _methods_trim_common(p):
    z = f"z = {p['z_target']}" if p.given('z_target') else "z_em"
    return (f"The spectrum was trimmed to the velocity range of +-{p['window_kms']} km/s "
            f"around {z} that all the sessions cover ({p['trans_self']} in this "
            f"session; {p['trans_others']} in {p['others_names']}){p.menu}")


def _methods_equalize_and_stitch(p):
    if p.given('manual_factors'):
        scale = f"was multiplied by the factors {p['manual_factors']}"
    elif p.given('equalize_ranges'):
        scale = f"was scaled to match in the ranges {p['equalize_ranges']} nm"
    else:
        p.raw('equalize_ranges')
        scale = "was scaled to match in their overlaps"
    return (f"The flux of this session and of {p['other_sessions']} {scale}, and the "
            f"sessions were joined at {p['stitch_wavelengths']} nm{p.menu}")


def _grid(p):
    text = f"a grid with a constant step of {p['dx']} {p['xunit']}"
    if p.given('xstart') or p.given('xend'):
        text += f", from {p['xstart'] if p.given('xstart') else 'the start of the data'}"
        text += f" to {p['xend'] if p.given('xend') else 'the end of the data'} nm"
    return text


def _clipping(p):
    return (f", with a {p['kappa']}-sigma clipping of outliers" if p.given('kappa')
            else ", without clipping")


def _methods_coadd(p):
    order = p.raw('equalize_order')
    if order is None or order < 0:
        scale = ""
    elif order == 0:
        scale = ", each scaled to the first by a constant factor"
    else:
        scale = f", each scaled to the first by a polynomial of order {order}"
    return (f"The sessions {p['session_names']} were co-added with this one on "
            f"{_grid(p)}, weighting each pixel by its inverse variance"
            f"{_clipping(p)}{scale}{p.menu}")


EDIT_RECIPES_SCHEMAS = {
    "set_properties": {
        "methods": _methods_set_properties,
        "brief": "Set session properties.",
        "url": "tutorials/editing_data.html#setting-session-properties",
        "menu": "Edit > Set &Properties",
        "provides": ("z_em", "resolution"),
        "details": "Sets the name of the session, the name of the object, the emission redshift, the resolution of the spectrum and the unit of its flux (a flux density per unit wavelength, declared, not converted: y_convert converts). A field left at '_current_' keeps its value.",
        "params": [
            {"name": "name", "type": str, "default": "_current_", "keep": ("_current_",), "doc": "Name of the session, as shown in the session list (not empty)."},
            {"name": "object", "type": str, "default": "_current_", "keep": ("_current_",), "doc": "Name of the object (the OBJECT keyword of the header)."},
            # Above -1 (1 + z > 0) and finite (Guido, phase 6-bis, F46)
            {"name": "z_em", "type": float, "default": "_current_", "keep": ("_current_",), "gt": -1, "doc": "Emission redshift of the quasar (z_em); needed to locate the Lyman-alpha forest."},
            # A number (R) or 'Npx', checked by the recipe: R > 0, N > 0
            {"name": "resol", "type": str, "default": "_current_", "keep": ("_current_",), "doc": "Resolution of the spectrum: the resolving power R (e.g. '50000'), or the FWHM of the line spread function in pixels (e.g. '3px', with the x axis in wavelength); written in the 'resol' column."},
            # Replace a resolution column that varies with wavelength by a
            # flat one; False (#526) keeps it. The GUI asks the user
            # Declared, not converted (PR 8-bis, Guido: a parameter of
            # set_properties, as z_em and resol); y_convert converts
            {"name": "yunit", "type": str, "default": "_current_", "keep": ("_current_",), "doc": "Unit of the flux, declared without converting the values, e.g. 'erg/(Angstrom s cm^2)', or 'arb_flux' (arbitrary units)."},
            {"name": "overwrite_resol_col", "type": bool, "default": False, "gui_hidden": True, "doc": "Replace an existing 'resol' column that varies along the spectrum (asked by the GUI)."}
        ],
    },
    "x_convert": {
        "methods": _methods_x_convert,
        "brief": "Convert X axis.",
        "url": "tutorials/editing_data.html#converting-the-axes",
        "details": "Converts the x axis to another unit of wavelength (e.g. nm, Angstrom) or to velocity (e.g. km/s); the zero of the velocity is the Lyman-alpha emission at the reference redshift, 121.567 nm * (1 + z_ref).",
        "params": [
            {"name": "z_ref", "type": float, "default": "_current_", "keep": ("_current_",), "gt": -1, "doc": "Redshift of the zero of the velocity axis, the Lyman-alpha emission at 121.567 nm * (1 + z_ref) ('_current_': z_em); not used for a unit of wavelength."},
            {"name": "xunit", "type": str, "default": "km/s", "doc": "Unit of the new x axis: a wavelength (e.g. 'nm', 'Angstrom') or a velocity (e.g. 'km/s')."},
        ],
    },
    "y_convert": {
        "methods": "The flux and its error were converted to {yunit}{menu}",
        "brief": "Convert Y axis units.",
        "url": "tutorials/editing_data.html#converting-the-axes",
        "details": "Converts the flux, its error and the columns in the unit of the flux (e.g. cont, model) to another unit of flux density per unit wavelength.",
        "params": [
            {"name": "yunit", "type": str, "default": "erg/(nm s cm^2)", "doc": "Unit of flux density per unit wavelength, e.g. 'erg/(nm s cm^2)' or 'erg/(Angstrom s cm^2)'; a unit per frequency (Jy) is not converted."},
        ],
    },
    "apply_expression": {
        "methods": "The column {target_col} was computed as {expression}{menu}",
        "brief": "Apply expression to columns.",
        "url": "tutorials/editing_data.html#arithmetic-on-columns",
        "menu": "Edit > Apply &Expression",
        "details": "Evaluates an expression of the columns (NumExpr: + - * / **, comparisons, sqrt, exp, log, log10, abs, where, sum, min, max, ...) and writes it in a column, new or existing. The columns are named as in the GUI: λ (or lam), F, dF, cont, model, ...; e.g. 'F / cont', 'F * 2.0', 'log10(F)'. The axis (λ, λmin, λmax) is not written: x_convert converts it.",
        "params": [
            # Required (PR 8): with the defaults 'F' and 'F / dF' it
            # replaced the flux with its S/N
            {"name": "target_col", "type": str, "default": "", "doc": "Column to write, e.g. 'snr': replaced if it exists, created otherwise; a name of letters, digits and _ (not the axis, λ, λmin, λmax).", "kind": "column"},
            {"name": "expression", "type": str, "default": "", "doc": "Expression of the columns, e.g. 'F / dF' or 'F / cont'."},
        ],
    },
    "mask_expression": {
        "methods": "The column {target_col} was masked where {expression}{menu}",
        "brief": "Mask column by expression.",
        "url": "tutorials/editing_data.html#masking-data",
        "menu": "Edit > &Mask by Expression",
        "details": "Sets a column to NaN where a boolean expression of the columns is true, e.g. '(λ < 300) | (λ > 400)' or 'dF <= 0'. The columns are named as in the GUI (λ or lam, F, dF, cont, ...). The axis (λ, λmin, λmax) is not masked.",
        "params": [
            {"name": "target_col", "type": str, "default": "F", "doc": "Column to mask, e.g. 'F' (not the axis, λ, λmin, λmax).", "kind": "column"},
            {"name": "expression", "type": str, "default": "dF <= 0", "doc": "Boolean expression of the columns: the pixels where it is true are masked."},
        ],
    },
    "smooth_column": {
        "methods": "The column {target_col} was smoothed with a Gaussian kernel of standard deviation {sigma_kms} km/s{menu}",
        "brief": "Smooth a single column.",
        "url": "tutorials/editing_data.html#smoothing-single-columns",
        "menu": "Edit > Smooth &Column",
        "details": "Smooths one column of the spectrum (e.g. 'dF' or 'cont') with a Gaussian kernel; the other columns are unchanged.",
        "params": [
            {"name": "target_col", "type": str, "default": "dF", "doc": "Column to smooth, e.g. 'dF' (not the axis, λ, λmin, λmax).", "kind": "column"},
            {"name": "sigma_kms", "type": float, "default": 100.0, "gt": 0, "doc": "Standard deviation of the Gaussian kernel (km/s)."},
            {"name": "renorm_model", "type": bool, "default": False, "doc": "When smoothing 'cont', rescale the model to the new continuum (asked by the GUI).", "gui_hidden": True, "agents": False},
        ],
    },
    "import_systems": {
        "methods": _methods_import_systems,
        "brief": "Import systems from session.",
        "url": "tutorials/editing_data.html#importing-systems",
        "menu": "File > Import Systems from &Session",
        "details": "Copies the absorption systems (components and constraints) of another open session into this one.",
        "params": [
            {"name": "source_session", "type": str, "default": "", "doc": "Name of the open session to copy the systems from (in a script, also the session itself).", "kind": "session"},
            {"name": "append", "type": bool, "default": True, "doc": "True: add the systems to those of this session; False: replace them."},
            {"name": "refit", "type": bool, "default": False, "doc": "Refit the systems on the spectrum of this session after copying them."}
        ],
    },
    "import_systems_ascii": {
        "methods": _methods_import_systems_ascii,
        "brief": "Import systems from an ASCII file.",
        "url": "tutorials/absorption_systems.html#importing-a-line-list",
        "menu": "File > Import Systems from &ASCII",
        "details": "Read absorption components from an ASCII/CSV table with a header (columns series, z and optionally logN, b, btur, dz, dlogN, db), e.g. the _systems.csv written by Export ASCII.",
        "params": [
            {"name": "file_path", "type": str, "default": "", "doc": "Path of the ASCII/CSV file to read.", "gui_hidden": True},
            {"name": "append", "type": bool, "default": True, "doc": "True: add the components to those of this session; False: replace them."},
            {"name": "refit", "type": bool, "default": False, "doc": "Refit the components on the spectrum after reading them."}
        ],
    },
    "delete": {
        "methods": _methods_delete,
        "brief": "Delete elements.",
        "url": "tutorials/editing_data.html#deleting-elements",
        "menu": "Edit > &Delete Elements",
        "details": "Removes columns from the spectrum, or all the components of the system list ('systems'). The columns of every spectrum (λ, λmin, λmax, F, dF) are not removed.",
        "params": [
            {"name": "targets", "type": str, "default": "", "doc": "Comma-separated items to delete: names of columns of the spectrum (e.g. 'cont, model') and/or 'systems' for the whole system list."}
        ],
    },
    "split": {
        "methods": "The pixels where {expression} were copied into a new session{menu}",
        "brief": "Split spectrum by expression.",
        "url": "tutorials/editing_data.html#manual-split",
        "menu": "Edit > S&plit Out Region",
        "branching": True,
        "details": "Copies the pixels where a boolean expression of the columns is true into a new session, e.g. '(λ > 400) & (λ < 500)'.",
        "params": [
            {"name": "expression", "type": str, "default": "", "doc": "Boolean expression of the columns, e.g. '(λ > 400) & (λ < 500)': the pixels where it is true go to the new session."}
        ],
    },
    "extract_preset": {
        "methods": _methods_extract_preset,
        "brief": "Extract region by preset.",
        "url": "tutorials/editing_data.html#quick-extraction-presets",
        "menu": "Edit > Extract &Preset Region",
        "branching": True,
        "requires": ("z_em",),
        "details": "Copies a standard region into a new session, from the emission redshift: 'lya_forest' (between Ly-beta and Ly-alpha), 'all_ly_forest' (between the Lyman limit and Ly-alpha), 'red_side' (redward of Ly-alpha). Several regions, separated by commas, give one session each.",
        "params": [
            # We use a string for now. Options: 'lya_forest', 'all_ly_forest', 'red_side'
            {"name": "region", "type": str, "default": "lya_forest", "doc": "Region(s) to extract, comma-separated: 'lya_forest', 'all_ly_forest', 'red_side'."}
        ],
    },
    "trim_common": {
        "methods": _methods_trim_common,
        "brief": "Trim to common velocity coverage.",
        "url": "tutorials/editing_data.html#combining-sessions-stitch-co-add",
        "details": "Trims this session to the velocity range around a redshift that all the given sessions cover, each with its own transition (e.g. the arms of an instrument, or the lines of a doublet).",
        "params": [
            {"name": "other_sessions", "type": str, "default": "", "kind": "sessions", "aliases": ["others_names"], "doc": "Comma-separated names of the other sessions (in a script, also a list of the sessions themselves)."},
            {"name": "z_target", "type": float, "default": "_current_", "keep": ("_current_",), "gt": -1, "doc": "Redshift at the centre of the range ('_current_': z_em)."},
            {"name": "trans_self", "type": str, "default": "Ly_a", "doc": "Transition of this session that sets its velocity zero, as named in the atomic data, e.g. 'Ly_a', 'CIV_1548'."},
            {"name": "trans_others", "type": str, "default": "", "doc": "Comma-separated transitions of the other sessions, one per session, in their order."},
            {"name": "window_kms", "type": float, "default": 500.0, "gt": 0, "doc": "Half-width of the velocity range kept (km/s)."}
        ],
    },
    "stitch": {
        "methods": "This session and {other_sessions} were joined into one spectrum, without equalizing their flux{menu}",
        "brief": "Merge multiple sessions.",
        "url": "tutorials/editing_data.html#stitching",
        "details": "Joins this session and other sessions into one spectrum, without equalizing their flux (see Equalize and stitch arms); the step that Equalize and stitch arms ends with.",
        "branching": True,
        "params": [
            {"name": "other_sessions", "type": str, "default": "", "kind": "sessions", "doc": "Comma-separated names of the other sessions (in a script, also a list of the sessions themselves). Their flux and axis are converted to the units of this session; a flux that cannot be converted (arbitrary, or normalized, against a physical one) is refused."},
            {"name": "sort", "type": bool, "default": True, "doc": "Sort the joined spectrum by wavelength."},
        ],
        "gui_hidden": True  # Usually called via script
    },
    "equalize_and_stitch": {
        "methods": _methods_equalize_and_stitch,
        "brief": "Equalize and stitch arms.",
        "url": "tutorials/editing_data.html#equalize-and-stitch-arms",
        "branching": True,
        "details": "Scales the flux of the sessions (the arms of an instrument) to match in their overlaps, then joins them into one spectrum, cutting at the given wavelengths. The arms are taken from blue to red, in any order given; their flux and axis are converted to the units of the bluest, and an arm in arbitrary units of flux takes those of the bluest once scaled.",
        "params": [
            {"name": "other_sessions", "type": str, "default": "", "kind": "sessions", "doc": "Comma-separated names of the other sessions (in a script, also a list of the sessions themselves).", "gui_hidden": True},
            {"name": "stitch_wavelengths", "type": str, "default": "", "doc": "Comma-separated wavelengths (nm), increasing, where each arm ends and the next begins: one fewer than the sessions."},
            {"name": "equalize_ranges", "type": str, "default": "auto", "doc": "Comma-separated wavelength ranges (nm) where the flux of adjacent arms is matched (e.g. '545-555, 795-805'), or 'auto' for their overlaps."},
            {"name": "manual_factors", "type": str, "default": "auto", "doc": "Comma-separated positive factors that multiply the flux of each arm but the bluest instead of the equalization, or 'auto'."}
        ],
    },
    "coadd": {
        "methods": _methods_coadd,
        "brief": "Co-add multiple sessions.",
        "url": "tutorials/editing_data.html#co-adding-spectra",
        "menu": "Edit > Co-&add Spectra",
        "branching": True,
        "details": "Combines several sessions on one new grid with a constant step, weighting each pixel with its inverse variance, in one pass. The first session is the reference: the others are converted to its units and scaled to its flux; an arbitrary flux is co-added with a physical one only when scaled.",
        "params": [
            {"name": "session_names", "type": str, "default": "", "kind": "sessions", "doc": "Comma-separated names of the sessions to co-add, this one among them (in a script, also a list of the sessions themselves); the first is the reference."},
            {"name": "xstart", "type": float, "default": None, "optional": True, "doc": "Start of the new grid (nm); None: the start of the data."},
            {"name": "xend", "type": float, "default": None, "optional": True, "doc": "End of the new grid (nm), above xstart; None: the end of the data."},
            {"name": "dx", "type": float, "default": 10.0, "gt": 0, "doc": "Step of the new grid, in the unit of xunit."},
            {"name": "xunit", "type": str, "default": "km/s", "doc": "Unit of the step: a velocity (km/s, a grid uniform in log wavelength) or a wavelength (nm)."},
            {"name": "kappa", "type": float, "default": 5.0, "optional": True, "gt": 0, "doc": "Threshold (in sigma) of the clipping of outliers; None: no clipping."},
            # max 5: polyfit is poorly conditioned from order 6-7 on overlaps
            # of 20-50 nm (measured in PR 9 of 6-quater); -1 (no scaling
            # before PR 9, F123) is still read as None
            {"name": "equalize_order", "type": int, "default": 0, "optional": True, "min": -1, "max": 5, "doc": "Order of the polynomial (0 to 5) that scales each session to the first before co-adding: 0 a constant factor; None no scaling (-1, of the old logs, is the same)."},
        ],
    }
}

@contextlib.contextmanager
def _expression_errors(expression: str):
    """An expression that cannot be evaluated (a syntax error, a name that
    is not a column, ...) fails the recipe with the expression in the
    message (phase 6-bis, PR 3; a SyntaxError reached the GUI as a bug)."""
    try:
        yield
    except (*DATA_ERRORS, SyntaxError, TypeError) as e:
        reason = reason_of(e)
        if isinstance(e, TypeError) and "is not callable" in str(e):
            # NumExpr names a function it does not have so (PR 8; F128)
            reason = ("it uses a function that the expressions do not have (they "
                      "have those of NumExpr: sqrt, exp, log, log10, abs, sin, cos, "
                      "where, sum, prod, min, max, ...).")
        raise RecipeFailed(f"The expression '{expression}' cannot be evaluated: "
                           f"{reason}") from e


SESSION_VAR_REGEX = re.compile(r'([a-zA-Z0-9_]+)\s*\.\s*([a-zA-Z0-9_]+)')

class RecipeEdit:
    """
    Recipes for general structure editing and arithmetic operations.

    These methods are accessed via the ``session.edit`` attribute and delegate
    heavy logic to the Core structures :class:`~astrocook.core.spectrum.Spectrum`
    and :class:`~astrocook.core.system_list.SystemList`.
    """
    def __init__(self, session_v2: 'Session'):
        self._session = session_v2
        self._tag = 'cb'

    def _translate_ui_names(self, name_str: str) -> str:
        """The names of the columns of every spectrum in ``name_str`` (an
        expression or a list of names) written as in the data (λ, lambda ->
        x, F -> y, dF -> dy, ...; :mod:`astrocook.core.column_names`, F120).
        The names of the data (x, y, dy) are left as they are."""
        if not isinstance(name_str, str):
            return name_str
        return internal_expression(name_str)

    def _column_to_write(self, target_col: str, must_exist: bool) -> str:
        """The name in the data of the column a recipe writes, checked
        (phase 6-quater, PR 8): a name usable in an expression, not the
        axis (F127), and one of the spectrum if ``must_exist``.
        :class:`InvalidParameter` otherwise."""
        name = target_col.strip()
        internal = internal_column(name)
        if not internal.isidentifier():
            # '' and 'a b' made columns that no expression could name
            raise InvalidParameter(f"target_col: {target_col!r} is not a name of a column "
                                   "(letters, digits and _, not starting with a digit)")
        if is_axis_column(name):
            raise InvalidParameter(f"target_col: {name} is the axis, which this recipe "
                                   "does not change (x_convert converts it)")
        colnames = spec_of(self._session).t.colnames
        if must_exist and internal not in colnames:
            raise InvalidParameter(f"target_col: there is no column {name}. The columns: "
                                   f"{', '.join(display_columns(colnames))}.")
        return internal

    def _prepare_expression_contexts(self, expression: str, alias_map: Optional[Dict[str, str]]) -> Tuple[str, Dict[str, np.ndarray]]:
        """
        Internal helper to resolve session aliases and delegate parsing to the Core.

        This prepares the context for multi-session arithmetic (e.g., subtracting
        one spectrum from another) by resolving GUI session names to actual
        Spectrum objects.

        Parameters
        ----------
        expression : str
            The raw expression string provided by the user (e.g. ``"F - s1.F"``).
        alias_map : dict
            A dictionary mapping aliases (e.g. ``"s1"``) to full session names.

        Returns
        -------
        tuple
            A tuple containing:
            1. **final_expression** (*str*): The processed expression compatible with ``numexpr``.
            2. **extra_vars** (*dict*): A dictionary of numpy arrays representing the
               resampled data from other sessions.

        See Also
        --------
        :meth:`astrocook.core.spectrum.Spectrum.prepare_expression_context`
        """
        if not alias_map:
            return expression, {}

        # 1. Resolve Aliases to Objects (GUI Logic)
        resolved_specs = {}
        
        # Access the GUI context
        gui = getattr(self._session, '_gui', None)
        if not gui or not hasattr(gui, 'session_histories'):
            # If running in script without GUI, we can't resolve aliases easily 
            # unless we passed the objects directly. For now, assume GUI context or fail gracefully.
            if alias_map:
                logger.warning("Other sessions can be named in an expression only in the GUI: "
                               "the aliases are not resolved.")
            return expression, {}

        for alias, full_name in alias_map.items():
            found = False
            for hist in gui.session_histories:
                if hist.display_name == full_name:
                    resolved_specs[alias] = hist.current_state.spec
                    found = True
                    break
            if not found:
                logger.warning(msg_session_not_found.format(name=full_name)
                               + f" The alias '{alias}' is not resolved.")

        # 2. Delegate to Core (Business Logic)
        return spec_of(self._session).prepare_expression_context(expression, resolved_specs)
    
    def _parse_resolution_string(self, resol_str: str) -> float:
        """
        Parses user input for resolution.

        Accepts either a raw floating point number (R) or a pixel width string
        (e.g., "3px") which is converted to R based on the current spectral sampling.

        Parameters
        ----------
        resol_str : str
            Input string, either a raw number (R) or a pixel width (e.g., "3px").

        Returns
        -------
        float
            The resolving power R. Returns 0.0 if parsing fails or spectrum is invalid.
        """
        clean_str = resol_str.strip().lower()
        
        # Check for pixel syntax
        if 'px' in clean_str or 'pixel' in clean_str:
            # Extract the number part
            match = re.match(r"([\d\.]+)", clean_str)
            if not match:
                raise ValueError(f"resol: '{resol_str}' is not a number of pixels ('Npx')")
            
            pix_fwhm = float(match.group(1))
            if not (np.isfinite(pix_fwhm) and pix_fwhm > 0):
                # '0px' gave an infinite R (Guido, phase 6-bis, F46)
                raise ValueError(f"resol: the FWHM in pixels must be positive ('{resol_str}')")
            logger.debug(f"Estimating the resolution for a FWHM of {pix_fwhm:.1f} pixels.")
            
            # 1. Get raw x array
            x_vals = spec_of(self._session).x.value
            
            if len(x_vals) > 1:
                # 2. Calculate Median Sampling (dx)
                dx = np.nanmedian(np.diff(x_vals))
                
                # 3. Calculate Median Wavelength (lambda)
                mid_x = np.nanmedian(x_vals)
                
                if dx > 0 and mid_x > 0:
                    # 4. Calculate R = lambda / (N_pix * dx)
                    fwhm_wave = pix_fwhm * dx
                    estimated_R = mid_x / fwhm_wave
                    
                    logger.warning(f"Resolution estimated from the sampling: R = {estimated_R:.0f} "
                                   f"at {mid_x:.1f}, for a FWHM of {pix_fwhm:.1f} pixels "
                                   "(an approximation: the sampling is taken as constant).")
                    return estimated_R
                else:
                    logger.warning("The resolution cannot be estimated (the wavelengths or their "
                                   "spacing are not positive): it is not changed.")
                    return 0.0
            else:
                logger.warning("The resolution cannot be estimated (fewer than two pixels): "
                               "it is not changed.")
                return 0.0
                
        # Fallback: Try standard float conversion
        try:
            resol_f = float(clean_str)
        except ValueError:
            raise ValueError(f"resol: '{resol_str}' is not a number (R) or 'Npx'") from None
        if not (np.isfinite(resol_f) and resol_f > 0):
            # 0, a negative R and NaN were taken as "no change", silently
            # (Guido, phase 6-bis, F46): '_current_' is the way to say it
            raise ValueError(f"resol: R must be positive and finite ('{resol_str}'); "
                             "leave it empty to keep the current resolution.")
        return resol_f
        

    @typed_recipe
    def set_properties(self, name: str = '_current_', object: str = '_current_',
                       z_em: Union[float, str] = '_current_', resol: str = '_current_',
                       overwrite_resol_col: bool = False,
                       yunit: str = '_current_') -> Optional['Session']:
        """
        Set session properties.

        Set core session properties like name, object, redshift, and resolution.
        Delegates to :meth:`astrocook.core.spectrum.Spectrum.with_properties`.

        Parameters
        ----------
        name : str, optional
            The session Name (display label). Defaults to ``'_current_'`` (no change).
        object : str, optional
            The object Name (header target). Defaults to ``'_current_'`` (no change).
        z_em : float or str, optional
            Emission Redshift, above -1 and finite. Defaults to
            ``'_current_'`` (no change).
        resol : str, optional
            Resolving Power R (e.g., ``50000``) OR Pixel FWHM (e.g., ``'3px'``),
            positive and finite. Defaults to ``'_current_'`` (no change).
        overwrite_resol_col : bool, optional
            If ``True``, a ``resol`` column that varies with wavelength is
            replaced by a flat one at the new resolution; if ``False`` it is
            kept (a flat or missing column is always set). Defaults to
            ``False`` (#526: before, a variable column was flattened even
            with ``resol='_current_'``).
        yunit : str, optional
            The unit of the flux, declared without converting the values
            (phase 6-quater, PR 8-bis): a flux density per unit wavelength
            (e.g. ``'erg/(Angstrom s cm^2)'``) or ``'arb_flux'``, the
            arbitrary units of Astrocook. Defaults to ``'_current_'`` (no
            change). :meth:`y_convert` converts a flux that has a unit.

        Returns
        -------
        Session
            A new :class:`~astrocook.session.session.Session` with updated properties.
        """
        if name != '_current_' and not name.strip():
            # An empty name was ignored without a word (PR 8)
            raise InvalidParameter("name: the name of a session cannot be empty "
                                   "('_current_' keeps it)")
        yunit_obj = None
        if str(yunit) != '_current_':
            try:
                yunit_obj = au.Unit(str(yunit))
            except ValueError:
                raise InvalidParameter(f"yunit: {yunit!r} is not a unit") from None
            if not (is_arbitrary_flux(yunit_obj) or is_flux_density_per_wavelength(yunit_obj)):
                raise InvalidParameter(f"yunit: {yunit!r} is not a flux density per unit "
                                       "wavelength (e.g. 'erg/(Angstrom s cm^2)') nor "
                                       "'arb_flux' (arbitrary units)")
        # 1. Parse numerical inputs (z_em, resol)
        try:
            if str(z_em) == '_current_':
                z_em_f = spec_of(self._session)._data.z_em
            else:
                z_em_f = float(z_em)
                
            if str(resol) == '_current_':
                resol_f = 0.0  # leave the resolution as it is (#526)
            else:
                resol_f = self._parse_resolution_string(str(resol))
        except ValueError as e:
            raise InvalidParameter(str(e)) from None

        # 2. Handle text inputs
        # Determine new values, falling back to current if '_current_' is passed
        new_name = name if name != '_current_' else self._session.name
        # The object left at '_current_' is not written again: None leaves the
        # header as it is (F108: it wrote OBJECT = '' when there was none, and
        # copied OBJECT over the ESO target name)
        new_object = object if object != '_current_' else None
        
        # 3. Apply changes via Spectrum.with_properties
        new_spec = spec_of(self._session).with_properties(
            z_em=z_em_f, 
            resol=resol_f,
            object_name=new_object,
            overwrite_resol_col=overwrite_resol_col 
        )

        if yunit_obj is not None:
            new_spec = new_spec.with_flux_unit(yunit_obj)

        # 4. Handle Session Name (Lives on the Session wrapper)
        new_session = self._session.with_new_spectrum(new_spec)
        new_session.name = new_name

        return new_session

    @typed_recipe
    def x_convert(self, z_ref: Union[float, str] = '_current_', xunit: str = 'km/s') -> Optional['Session']:
        """
        Convert X axis.

        Convert the x axis units and zero point to wavelength or velocity units.
        Delegates to :meth:`astrocook.core.spectrum.Spectrum.x_convert`.

        Parameters
        ----------
        z_ref : float or str, optional
            Reference redshift for velocity conversion, above -1 and finite;
            ``'_current_'`` uses the session's z_em. Defaults to
            ``'_current_'``, as the GUI.
        xunit : str, optional
            Unit of wavelength or velocity (e.g., ``'nm'``, ``'km/s'``).
            Defaults to ``'km/s'``.

        Returns
        -------
        Session
            A new :class:`~astrocook.session.session.Session` with the converted X-axis.
        """
        # Handle default logic here
        if z_ref == '_current_':
            z_ref_float = spec_of(self._session)._data.z_em
        else:
            z_ref_float = cast(float, z_ref)
            
        try:
            xunit_obj = au.Unit(xunit)
        except ValueError:
            raise InvalidParameter(f"xunit: {xunit!r} is not a unit") from None
        if not (xunit_obj.is_equivalent(au.nm) or xunit_obj.is_equivalent(au.km / au.s)):
            # It failed in the conversion, with the same outcome (PR 8)
            raise InvalidParameter(f"xunit: {xunit!r} is neither a wavelength (e.g. nm) "
                                   "nor a velocity (e.g. km/s)")

        # 1. Get the new Spectrum instance (the immutable operation)
        new_spec_v2 = spec_of(self._session).x_convert(z_ref=z_ref_float, xunit=xunit_obj)
        
        # 2. Return a NEW Session instance with the updated spectrum
        return self._session.with_new_spectrum(new_spec_v2) 
    
    @typed_recipe
    def y_convert(self, yunit: str = 'erg/(nm s cm^2)') -> Optional['Session']:
        """
        Convert Y axis units.

        Converts the flux and error arrays to the specified unit.
        Delegates to :meth:`astrocook.core.spectrum.Spectrum.y_convert`.

        Parameters
        ----------
        yunit : str, optional
            Target unit of flux density. Defaults to ``'erg/(nm s cm^2)'``.

        Returns
        -------
        Session
            A new :class:`~astrocook.session.session.Session` with the converted Y-axis.
        """
        # C. Parameter Validation and Type Casting
        try:
            yunit_obj = au.Unit(yunit)
        except ValueError:
            raise InvalidParameter(f"yunit: {yunit!r} is not a unit") from None
        y_unit = spec_of(self._session).y.unit
        if is_arbitrary_flux(y_unit) or has_no_unit(y_unit):
            # PR 8-bis: the message says how to give the flux a unit
            raise InvalidParameter("yunit: the flux is in arbitrary units and cannot be "
                                   "converted: declare its unit with set_properties "
                                   "(yunit), or calibrate it with "
                                   "calibrate_from_magnitudes")
        if not yunit_obj.is_equivalent(y_unit):
            # It failed in the conversion, with the same outcome (PR 8)
            raise InvalidParameter(f"yunit: the flux ({y_unit}) cannot be converted to "
                                   f"{yunit!r}: a flux density per unit wavelength, e.g. "
                                   "'erg/(nm s cm^2)'")

        # 1. Execute the immutable V2 operation
        new_spec_v2 = spec_of(self._session).y_convert(yunit=yunit_obj)
        
        # 2. Return a NEW Session instance (the new state)
        return self._session.with_new_spectrum(new_spec_v2)
    
    @typed_recipe
    def apply_expression(self, target_col: str, expression: str, alias_map: Optional[Dict[str, str]] = None) -> 'Session':
        """
        Apply expression to columns.

        Apply a NumPy-style expression. Use column names (λ, F, cont...) as variables.
        Supports multi-session arithmetic via aliases.
        Delegates to :meth:`astrocook.core.spectrum.Spectrum.apply_expression`.

        Parameters
        ----------
        target_col : str
            Target column (to be overwritten or created).
        expression : str
            Expression to evaluate (e.g., ``'F / dF'``, ``'F * 2.0'``).
        alias_map : dict, optional
            Mapping of aliases to session names (e.g. ``{'s1': 'Session 1'}``).

        Returns
        -------
        Session
            A new :class:`~astrocook.session.session.Session` with the computed column.
        """
        # The names of the GUI (or of the data) to those of the data
        target_col = self._column_to_write(target_col, must_exist=False)
        expression = self._translate_ui_names(expression.strip())
        if not expression:
            raise RecipeFailed("The expression is empty.")
        # 1. Prepare multi-session vars
        final_expression, extra_vars = self._prepare_expression_contexts(expression, alias_map)
        logger.debug(f"Original expression: '{expression}'")
        logger.debug(f"Final expression: '{final_expression}'")
        logger.debug(f"Extra vars: {list(extra_vars.keys())}")

        # 2. Call the API with the prepared data
        with _expression_errors(expression):
            new_spec_v2 = spec_of(self._session).apply_expression(
                target_col=target_col,
                expression=final_expression,
                extra_vars=extra_vars
            )
        return self._session.with_new_spectrum(new_spec_v2)
    @typed_recipe
    def mask_expression(self, target_col: str, expression: str, alias_map: Optional[Dict[str, str]] = None) -> 'Session':
        """
        Mask column by expression.

        Mask a target column by setting values to NaN where the expression is True.
        Delegates to :meth:`astrocook.core.spectrum.Spectrum.mask_expression`.

        Parameters
        ----------
        target_col : str
            Target column to mask (e.g., ``'F'``).
        expression : str
            Boolean expression (e.g., ``'(λ < 300) | (λ > 400)'``, ``'dF <= 0'``).
        alias_map : dict, optional
            Mapping of aliases to session names.

        Returns
        -------
        Session
            A new :class:`~astrocook.session.session.Session` with the masked column.
        """
        # The names of the GUI (or of the data) to those of the data
        target_col = self._column_to_write(target_col, must_exist=True)
        expression = self._translate_ui_names(expression.strip())
        if not expression:
            raise RecipeFailed("The expression is empty.")
        # 1. Prepare multi-session vars
        final_expression, extra_vars = self._prepare_expression_contexts(expression, alias_map)
        logger.debug(f"Original expression: '{expression}'")
        logger.debug(f"Final expression: '{final_expression}'")
        logger.debug(f"Extra vars: {list(extra_vars.keys())}")
        
        # 2. Call the API with the prepared data
        with _expression_errors(expression):
            new_spec_v2 = spec_of(self._session).mask_expression(
                target_col=target_col,
                expression=final_expression,
                extra_vars=extra_vars
            )
        return self._session.with_new_spectrum(new_spec_v2)
    @typed_recipe
    def smooth_column(self, target_col: str = 'dF', sigma_kms: float = 100.0, renorm_model: bool = False) -> 'Session':
        """
        Smooth a single column.

        Apply a Gaussian filter to a single data column.
        Delegates to :meth:`astrocook.core.spectrum.Spectrum.smooth_column`.

        Parameters
        ----------
        target_col : str, optional
            Target column to smooth. Defaults to ``'dF'``.
        sigma_kms : float, optional
            Standard deviation for Gaussian kernel (km/s), positive.
            Defaults to ``100.0``.
        renorm_model : bool, optional
            If ``True`` and smoothing ``'cont'``, the ``'model'`` column is
            re-normalized. Defaults to ``False``.

        Returns
        -------
        Session
            A new :class:`~astrocook.session.session.Session` with the smoothed column.
        """
        # The names of the GUI (or of the data) to those of the data
        target_col = self._column_to_write(target_col, must_exist=True)
        # 1. Call the immutable V2 operation
        new_spec_v2 = spec_of(self._session).smooth_column(
            target_col=target_col,
            sigma_kms=sigma_kms,
            renorm_model=renorm_model
        )
        
        # 2. Return a NEW Session instance
        return self._session.with_new_spectrum(new_spec_v2)
        
    @typed_recipe
    def import_systems(self, source_session: Union[str, 'Session'], append: bool = True,
                       refit: bool = False) -> 'Session':
        """
        Import systems from session.

        Copy absorption systems from another open session into this one, filtering
        those that fall outside the current spectral range.
        Delegates filtering to :meth:`astrocook.core.system_list.SystemList.filter_by_range`.

        Parameters
        ----------
        source_session : str
            Name of the open session to import from (required; in a script,
            also the session itself).
        append : bool, optional
            If ``True``, adds to existing systems. If ``False``, replaces them.
            Defaults to ``True``.
        refit : bool, optional
            If ``True``, runs a refit on the imported systems. Defaults to ``False``.

        Returns
        -------
        Session
            A new :class:`~astrocook.session.session.Session` with imported systems.
        """

        # A name not found fails (it returned the session unchanged, which
        # looked like a success; without a GUI, an AttributeError, F135)
        source = open_session(self._session, source_session, 'source_session')
        label = source.name if is_session(source_session) else source_session
        source_systs = source.systs
        if not source_systs or not source_systs.components:
            logger.warning(f"The session '{label}' has no systems: nothing is imported.")
            return self._session

        return self._merge_imported_systems(source_systs, append, refit, label=f"'{label}'")

    def _merge_imported_systems(self, source_systs, do_append: bool, do_refit: bool,
                                label: str) -> 'Session':
        """Filter by spectral range, merge and optionally refit imported systems
        (shared by :meth:`import_systems` and :meth:`import_systems_ascii`)."""
        import numpy as np
        from astrocook.core.structures import SystemListData
        from astrocook.core.system_list import SystemList

        # 2. Prepare current list
        if do_append and self._session.systs:
            base_systs = self._session.systs
        else:
            base_systs = SystemList(SystemListData())

        # 3. Filter Source Components (Delegated to Core)
        if self._session.spec and len(spec_of(self._session).x) > 0:
            try:
                xmin = np.min(spec_of(self._session).x.value)
                xmax = np.max(spec_of(self._session).x.value)
                
                filtered_source_systs = source_systs.filter_by_range(xmin, xmax)
                
                skipped = len(source_systs.components) - len(filtered_source_systs.components)
                if skipped > 0:
                    logger.warning(f"{skipped} systems outside the spectrum are not imported.")
                    
                source_systs = filtered_source_systs 
                
            except (ValueError, TypeError, KeyError) as e:
                logger.warning(f"The systems could not be filtered by the range of the "
                               f"spectrum ({e}): all of them are imported.")

        if not source_systs.components:
            logger.warning("No systems fall within the spectrum: nothing is imported.")
            return self._session

        # 4. Merge
        logger.info(f"Importing {len(source_systs.components)} systems from {label}.")
        new_systs = base_systs.merge(source_systs, copy_uuids=False)
        
        # Create intermediate session with imported (but unfitted) systems
        intermediate_session = self._session.with_new_system_list(new_systs)

        # 5. Optional Refit
        if do_refit:
            logger.info("Refitting the imported systems.")
            # Local import to avoid circular dependency
            from astrocook.recipes.absorbers import RecipeAbsorbers

            # Wrap the intermediate session in the Absorbers recipe
            rec_abs = RecipeAbsorbers(intermediate_session)

            # Run refit_all
            # You can tweak z_window_kms or max_nfev here if needed, 
            # or expose them as advanced params in the schema.
            # max_nfev as before the defaults were aligned with the GUI (#573)
            try:
                return rec_abs.refit_all(z_window_kms=20.0, group_depth=1, max_nfev='100')
            except RecipeError as e:
                # The import stays (Guido's choice, phase 6), with the reason
                logger.warning(f"The imported systems could not be refitted ({e.reason}): "
                               "they keep their imported values.")

        # Return unfitted result
        return intermediate_session
        
    @typed_recipe
    def import_systems_ascii(self, file_path: str, append: bool = True,
                             refit: bool = False) -> 'Session':
        """
        Import systems from an ASCII file (issue #480).

        Reads a table with a header, in any format astropy can guess (CSV,
        whitespace-separated, ...), e.g. the ``_systems.csv`` written by
        :meth:`astrocook.recipes.file.RecipeFile.export_ascii`. Required
        columns: ``series`` and ``z``; optional: ``logN`` (default 13.0), ``b``
        (default 10.0 km/s), ``btur`` (default 0), ``dz``, ``dlogN``, ``db``.
        Other columns are ignored. Components outside the spectral range are
        skipped.

        Parameters
        ----------
        file_path : str
            Path of the file.
        append : bool, optional
            If ``True``, add to the existing systems; otherwise replace them.
            Defaults to ``True``.
        refit : bool, optional
            If ``True``, refit the systems after importing. Defaults to ``False``.

        Returns
        -------
        Session
            A new session with the imported systems (the session as it was
            if the file holds no valid row).

        Raises
        ------
        RecipeFailed
            If the file cannot be read or lacks the required columns.
        """
        import numpy as np
        from astropy.io import ascii as ap_ascii
        from astrocook.core.structures import ComponentData, SystemListData
        from astrocook.core.system_list import SystemList

        try:
            table = ap_ascii.read(file_path, guess=True)
        except (OSError, ValueError) as e:
            # ValueError covers the table errors of astropy and the decoding
            raise RecipeFailed(f"Could not read '{file_path}': {e}") from e
        missing = [c for c in ('series', 'z') if c not in table.colnames]
        if missing:
            raise RecipeFailed(f"'{file_path}' lacks the column(s) {', '.join(missing)}.")

        def value(row, name, default):
            if name not in table.colnames:
                return default
            v = row[name]
            if np.ma.is_masked(v):
                return default
            try:
                v = float(v)
            except (TypeError, ValueError):
                return default
            return default if not np.isfinite(v) else v

        comps = []
        for i, row in enumerate(table):
            try:
                comps.append(ComponentData(
                    id=i + 1, series=str(row['series']).strip(), z=float(row['z']),
                    dz=value(row, 'dz', None), logN=value(row, 'logN', 13.0),
                    dlogN=value(row, 'dlogN', None), b=value(row, 'b', 10.0),
                    db=value(row, 'db', None), btur=value(row, 'btur', 0.0)))
            except (TypeError, ValueError) as e:
                logger.warning(f"Row {i + 1} of '{file_path}' skipped: {e}")
        if not comps:
            logger.warning(f"No valid systems in '{file_path}': nothing is imported.")
            return self._session

        source_systs = SystemList(SystemListData(components=comps))
        return self._merge_imported_systems(source_systs, append, refit,
                                            label=f"'{file_path}'")

    @typed_recipe
    def delete(self, targets: str) -> 'Session':
        """
        Delete elements.

        Remove columns from the spectrum or clear the line list.
        Delegates to :meth:`astrocook.core.spectrum.Spectrum.remove_columns`.

        Parameters
        ----------
        targets : str
            Comma-separated list of items to delete (e.g. ``'cont, model'``),
            columns of the spectrum other than its axis, flux and error.
            Use ``'systems'`` (or ``'lines'``) to clear the component list.

        Returns
        -------
        Session
            A new :class:`~astrocook.session.session.Session` with the elements removed.

        Raises
        ------
        InvalidParameter
            If ``targets`` names nothing, a column of every spectrum (λ,
            λmin, λmax, F, dF) or a column that the spectrum does not have
            (phase 6-quater, PR 9: each gave a WARNING and an unchanged
            state, recorded as done).
        """
        target_list = [t.strip() for t in targets.split(',') if t.strip()]
        if not target_list:
            raise InvalidParameter("targets: name at least one column, or 'systems'")
        clear_systems = any(t in ('lines', 'systems') for t in target_list)
        columns = [t for t in target_list if t not in ('lines', 'systems')]
        spec = spec_of(self._session)
        core = [c for c in columns if is_core_column(c)]
        if core:
            raise InvalidParameter(f"targets: {', '.join(core)} cannot be deleted: λ, λmin, "
                                   "λmax, F and dF are the data of every spectrum")
        aux = {k.strip() for k in spec._data.aux_cols}
        missing = [c for c in columns if c not in aux]
        if missing:
            raise InvalidParameter(f"targets: there is no column {', '.join(missing)}. The "
                                   f"columns that can be deleted: {', '.join(sorted(aux)) or 'none'}.")

        new_systs = self._session.systs
        if clear_systems:
            new_systs = SystemList(SystemListData())
            logger.info("The systems are removed.")
        new_spec = spec.remove_columns(columns) if columns else spec
        return self._session.with_new_spectrum(new_spec).with_new_system_list(new_systs)
    
    @typed_recipe
    def split(self, expression: str, alias_map: Optional[Dict[str, str]] = None) -> 'Session':
        """
        Split spectrum by expression.

        Extract a portion of the spectrum into a new session using a boolean expression.
        Only rows where the expression is True are kept.
        Delegates to :meth:`astrocook.core.spectrum.Spectrum.split`.

        Parameters
        ----------
        expression : str
            Boolean expression to select data (e.g., ``'(λ > 400) & (λ < 500)'``);
            required (no default since phase 6-quater, PR 9).
        alias_map : dict, optional
            Mapping for multi-session variables.

        Returns
        -------
        Session
            A new :class:`~astrocook.session.session.Session` containing the subset of data.
        """
        expression = self._translate_ui_names(expression.strip())
        if not expression:
            raise InvalidParameter("expression: give a boolean expression of the columns, "
                                   "e.g. '(λ > 400) & (λ < 500)'")
            
        # 1. Prepare multi-session vars (re-using our existing helper!)
        final_expression, extra_vars = self._prepare_expression_contexts(expression, alias_map)
        
        # 2. Call the immutable V2 operation
        with _expression_errors(expression):
            new_spec_v2 = spec_of(self._session).split(
                expression=final_expression,
                extra_vars=extra_vars
            )

        new_systs = self._session.systs
        if new_spec_v2 and len(new_spec_v2.x) > 0 and new_systs:
            # Safely get min/max in nanometers, regardless of current plot units
            xmin_nm = np.min(new_spec_v2.x.to(au.nm).value)
            xmax_nm = np.max(new_spec_v2.x.to(au.nm).value)
            
            original_count = len(new_systs.components)
            new_systs = new_systs.filter_by_range(xmin_nm, xmax_nm)
            filtered_count = len(new_systs.components)
            
            if original_count != filtered_count:
                logger.info(f"{original_count - filtered_count} systems outside the new range removed.")
        
        # 4. Return a NEW Session instance (GUI handles branching)
        return self._session.with_new_spectrum(new_spec_v2).with_new_system_list(new_systs)
    
    @typed_recipe
    def extract_preset(self, region: str = 'lya_forest') -> Union['Session', List['Session']]:
        """
        Extract region by preset.

        Extract a specific standard region (like the Ly-alpha forest) into a new
        session, based on the current emission redshift (z_em).
        Delegates to :meth:`astrocook.core.spectrum.Spectrum.get_region_bounds`
        and :meth:`split`.

        Parameters
        ----------
        region : str, optional
            Preset region identifier. Options: ``'lya_forest'``, ``'all_ly_forest'``,
            ``'red_side'``. Can be a comma-separated list. Defaults to ``'lya_forest'``.

        Returns
        -------
        Session or List[Session]
            One or more new :class:`~astrocook.session.session.Session` instances.

        Raises
        ------
        InvalidParameter
            If ``region`` names no region, or one that is not a preset (phase
            6-quater, PR 9: an unknown one was skipped with a WARNING).
        """
        region_list = [r.strip() for r in region.split(',') if r.strip()]
        options = ", ".join(REGION_PRESETS)
        if not region_list:
            raise InvalidParameter(f"region: name at least one region: {options}")
        unknown = [r for r in region_list if r not in REGION_PRESETS]
        if unknown:
            raise InvalidParameter(f"region: {', '.join(unknown)} is not a region; the "
                                   f"regions are {options}")
        
        results = []
        failures = []

        for reg in region_list:
            try:
                # Calculate bounds
                obs_min, obs_max = spec_of(self._session).get_region_bounds(reg)
                
                # Convert to native units for expression
                spec_unit = spec_of(self._session).x.unit
                v_min = obs_min.to(spec_unit).value
                v_max = obs_max.to(spec_unit).value
                
                # Build expression
                expression = f"(λ > {v_min:.4f}) & (λ < {v_max:.4f})"
                logger.info(f"Extracting {reg} ({expression}, z_em = {spec_of(self._session)._data.z_em:.4f}).")
                
                results.append(self.split(expression=expression))
                
            except (ValueError, KeyError, TypeError, RecipeError) as e:
                # The other regions are extracted; with none, the recipe
                # fails with the reasons (phase 6-bis, PR 3)
                reason = e.reason if isinstance(e, RecipeError) else reason_of(e)
                failures.append(f"{reg}: {reason}")
                logger.warning(f"The region {reg} could not be extracted ({reason}): "
                               "it is skipped.")

        if not results:
            raise RecipeFailed("No region could be extracted (" + "; ".join(failures) + ")")

        # Return single object if one, or list if multiple
        if len(results) == 1:
            return results[0]
        
        return results
    
    @typed_recipe
    def trim_common(self, other_sessions: Union[str, List[Union['Session', str]]],
                    z_target: Union[float, str], trans_self: str,
                    trans_others: str, window_kms: float = 500.0) -> 'Session':
        """
        Trim to common velocity coverage.

        Trims the current session to the velocity range shared by a list of other
        sessions (intersection).
        Delegates to :meth:`astrocook.core.spectrum.Spectrum.calc_intersection_expression`.

        Parameters
        ----------
        other_sessions : str or list
            Comma-separated names of other sessions (in a script, also a
            list of the sessions themselves). ``others_names`` until phase
            6-quater, PR 9, still accepted.
        z_target : float
            Target redshift, above -1 and finite.
        trans_self : str
            Transition for current session (e.g. ``'Ly_a'``), as named in the
            atomic data.
        trans_others : str
            Comma-separated transitions for other sessions, one per session.
        window_kms : float, optional
            Symmetric window width in km/s, positive. Defaults to ``500.0``.

        Returns
        -------
        Session
            A new :class:`~astrocook.session.session.Session` trimmed to the common coverage.

        Raises
        ------
        InvalidParameter
            If a transition is not in the atomic data, or the transitions are
            not one per session (phase 6-quater, PR 9: an unknown transition
            of this session gave "no range in common", one of another session
            was ignored).
        """
        names = split_names(other_sessions)
        trans_list_others = [t.strip() for t in trans_others.split(',') if t.strip()]
        unknown = [t for t in [trans_self.strip(), *trans_list_others] if t not in ATOM_DATA]
        if unknown:
            raise InvalidParameter(f"trans_self, trans_others: {', '.join(unknown)} is not a "
                                   "transition of the atomic data (e.g. 'Ly_a', 'CIV_1548')")
        if len(names) != len(trans_list_others):
            raise InvalidParameter(f"trans_others: there are {len(names)} sessions and "
                                   f"{len(trans_list_others)} transitions: give one per session")
        others = open_sessions(self._session, names, 'other_sessions')

        split_expression = spec_of(self._session).calc_intersection_expression(
            others_specs=[spec_of(s) for s in others],
            z_target=float(z_target),
            trans_self=trans_self.strip(),
            trans_others=trans_list_others,
            window_kms=float(window_kms)
        )
        
        if not split_expression:
            # (A WARNING and the session unchanged, phase 6)
            raise RecipeFailed("The sessions have no range in common: nothing is trimmed.")

        logger.info(f"Trimming to the range in common: {split_expression}.")
        return self.split(split_expression)

    @typed_recipe
    def stitch(self, other_sessions: Union[str, List[Union['Session', str]]],
               sort: bool = True) -> 'Session':
        """
        Merge multiple sessions.

        Concatenates the current session with a list of other sessions.
        Delegates to :meth:`astrocook.core.spectrum.Spectrum.stitch`.

        Parameters
        ----------
        other_sessions : str or list
            Comma-separated names of the other sessions, or a list of names
            or (in a script) of :class:`~astrocook.session.session.Session`
            objects. Their flux and axis are converted to the units of this
            session.
        sort : bool, optional
            If ``True``, sorts the final spectrum by wavelength. Defaults to ``True``.

        Returns
        -------
        Session
            A new :class:`~astrocook.session.session.Session` containing the merged data.

        Raises
        ------
        RecipeFailed
            If a name is not that of an open session (phase 6-quater, PR 9:
            skipped with a WARNING), or a flux or an axis cannot be converted
            to the units of this session.
        """
        items = split_names(other_sessions)
        others = open_sessions(self._session, items, 'other_sessions')
        names = [self._session.name] + _labels(items, others)
        specs = in_units_of_first([spec_of(self._session)] + [spec_of(s) for s in others],
                                   names, scaled=False)
        new_spec = specs[0].stitch(specs[1:], sort=sort)

        # Remove NaNs which break the 'step' plot (causing "scatter" look)
        # We use the public split API to keep valid data only.
        try:
            # Check for NaNs in X or Y
            x_val = new_spec.x.value
            y_val = new_spec.y.value
            if not (np.all(np.isfinite(x_val)) and np.all(np.isfinite(y_val))):
                logger.debug("Removing the pixels with NaN values.")
                new_spec = new_spec.split("isfinite(x) & isfinite(y)")
        except (ValueError, KeyError, TypeError) as clean_err:
            logger.warning(f"The pixels with NaN values could not be removed "
                           f"({clean_err}): they are kept.")
        
        return self._session.with_new_spectrum(new_spec)
        
    @typed_recipe
    def equalize_and_stitch(self, other_sessions: Union[str, List[Union['Session', str]]],
                            stitch_wavelengths: str, equalize_ranges: str = 'auto',
                            manual_factors: str = 'auto') -> 'Session':
        """
        Equalize and stitch arms.

        Equalize flux levels and stitch multiple sessions at specific cut-off wavelengths.
        This recipe is designed to seamlessly combine data from multi-arm spectrographs 
        (e.g., X-Shooter's UVB, VIS, and NIR arms).

        The procedure enforces a strict blue-to-red processing order regardless of 
        input order. For $N$ total sessions, exactly $N-1$ parameters must be provided 
        for cuts, ranges, and factors. The arms are converted to the units of
        the axis and of the flux of the bluest (an arm in arbitrary units of
        flux takes those of the bluest once scaled; phase 6-quater, PR 9,
        F132, F134).

        Parameters
        ----------
        other_sessions : str or list
            Comma-separated names of the other sessions to stitch to the current one
            (in a script, also a list of the sessions themselves).
            (This parameter is hidden in the GUI dialog).
        stitch_wavelengths : str
            Comma-separated cut-off wavelengths (in nm), increasing, where the
            arms join. Data bluer than a cut is taken from the bluer arm; data
            redder is taken from the redder arm. Requires $N-1$ values.
        equalize_ranges : str, optional
            Comma-separated overlap ranges formatted as ``'min-max, min-max'`` (in nm) 
            used to calculate the scaling factor between adjacent arms. If set to 
            ``'auto'``, the recipe attempts a global match. Defaults to ``'auto'``.
        manual_factors : str, optional
            Comma-separated fixed positive multipliers to scale adjacent arms. If provided, 
            these override the `equalize_ranges` calculation. Defaults to ``'auto'``.

        Returns
        -------
        Session
            A new :class:`~astrocook.session.session.Session` containing the equalized 
            and stitched spectrum.
        """
        names = split_names(other_sessions)
        if not names:
            raise InvalidParameter("other_sessions: name at least one session")
        n = len(names)
        cuts = _numbers('stitch_wavelengths', stitch_wavelengths, n)
        if any(b <= a for a, b in pairwise(cuts)):
            raise InvalidParameter("stitch_wavelengths: give them increasing, from blue to red")
        # Factors or ranges that cannot be read make the recipe fail: they
        # were replaced by 'auto', with a WARNING (phase 6)
        factors: List[Optional[float]] = [None] * n
        if manual_factors.strip().lower() != 'auto':
            factors = list(_numbers('manual_factors', manual_factors, n, positive=True))
        ranges: List[Optional[Tuple[float, float]]] = [None] * n
        if equalize_ranges.strip().lower() != 'auto' and factors[0] is None:
            ranges = list(_ranges('equalize_ranges', equalize_ranges, n))

        others = open_sessions(self._session, names, 'other_sessions')
        arms = _blue_to_red([(self._session.name, spec_of(self._session))]
                            + list(zip(_labels(names, others), [spec_of(s) for s in others],
                                       strict=True)))
        arm_names = [name for name, _ in arms]
        specs = in_units_of_first([spec for _, spec in arms], arm_names, scaled=True)
        unit = specs[0].x.unit
        bounds = [None, *cuts, None]

        # Each arm scaled to the arms bluer than it, joined, then cut
        reference = specs[0]
        parts = [_arm_part(specs[0], arm_names[0], None, cuts[0], unit)]
        for i, arm in enumerate(specs[1:]):
            equalized = _equalize_arm(arm, arm_names[i + 1], reference, factors[i],
                                      ranges[i], unit)
            parts.append(_arm_part(equalized, arm_names[i + 1], bounds[i + 1],
                                   bounds[i + 2], unit))
            reference = reference.stitch([equalized], sort=True)

        logger.debug("Stitching the parts together.")
        final_spec = parts[0].stitch(parts[1:], sort=True)
        new_session = self._session.with_new_spectrum(final_spec)
        new_session.name = "Stitched_" + "_".join(arm_names)
        return new_session
        
    @typed_recipe
    def coadd(self, session_names: Union[str, List[Union['Session', str]]],
              xstart: Optional[float] = None, xend: Optional[float] = None,
              dx: float = 10.0, xunit: str = 'km/s', kappa: Optional[float] = 5.0,
              equalize_order: Optional[int] = 0) -> 'Session':
        """
        Co-add multiple sessions.

        Stitches multiple sessions together and rebins them onto a single common
        grid in one pass (Inverse Variance Weighting).
        Delegates to :meth:`astrocook.core.spectrum.Spectrum.coadd`.

        Parameters
        ----------
        session_names : str or list
            Comma-separated names of sessions to co-add, this one among them
            (in a script, also a list of the sessions themselves). The first
            is the reference: the others are converted to the units of its
            axis and flux, and scaled to its flux (F133, F134).
        xstart : float or None, optional
            Start wavelength (nm). None for auto. Defaults to None.
        xend : float or None, optional
            End wavelength (nm), above ``xstart``. None for auto. Defaults to None.
        dx : float, optional
            Step size for the final grid, positive. Defaults to ``10.0`` (km/s), as
            the GUI: choose the native pixel size for high-resolution data
            (e.g. ~2.5 km/s for UVES, ~0.8 km/s for ESPRESSO).
        xunit : str, optional
            Unit for the step size: a velocity or a wavelength (e.g.,
            ``'km/s'`` or ``'nm'``). Defaults to ``'km/s'``.
        kappa : float or None, optional
            Sigma clipping threshold, positive. None for off. Defaults to ``5.0``.
        equalize_order : int or None, optional
            Order of the polynomial (0 to 5) that scales each session to the
            first: 0 a constant factor, None no scaling (``-1``, the value of
            no scaling until phase 6-quater, PR 9, F123, is read as None).
            Defaults to ``0``.

        Returns
        -------
        Session
            A new :class:`~astrocook.session.session.Session` containing the co-added spectrum.

        Raises
        ------
        InvalidParameter
            If ``session_names`` names no session or not this one, ``xstart``
            is not below ``xend`` (an empty spectrum) or ``xunit`` is neither
            a velocity nor a wavelength (phase 6-quater, PR 9).
        RecipeFailed
            If a name is not that of an open session (it was skipped without
            a word), or a flux or an axis cannot be converted.
        """
        names = split_names(session_names)
        if not names:
            raise InvalidParameter("session_names: name at least one session")
        dx_unit = _grid_unit(xunit)
        if xstart is not None and xend is not None and not xstart < xend:
            raise InvalidParameter(f"xstart: {xstart} is not below xend ({xend})")
        sessions = open_sessions(self._session, names, 'session_names')
        if not (any(s is self._session for s in sessions) or self._session.name in names):
            raise InvalidParameter(f"session_names: name this session too "
                                   f"('{self._session.name}'): the result is a new state of it")
        eq_order = -1 if equalize_order is None or equalize_order < 0 else equalize_order
        specs = in_units_of_first([spec_of(s) for s in sessions], _labels(names, sessions),
                                   scaled=eq_order >= 0)

        xstart_q = au.Quantity(xstart, au.nm) if xstart is not None else None
        xend_q = au.Quantity(xend, au.nm) if xend is not None else None
        dx_q = au.Quantity(dx, dx_unit)
        logger.info(f"Combining {len(specs)} sessions on a grid of step {dx_q}.")
        new_spec = specs[0].coadd(others=specs[1:], xstart=xstart_q, xend=xend_q, dx=dx_q,
                                  kappa=kappa, equalize_order=eq_order)
        return self._session.with_new_spectrum(new_spec)


# --- Helpers of the recipes on several spectra (phase 6-quater, PR 9) ---

def _labels(items: list, sessions: list) -> List[str]:
    """The names of the sessions in the messages: those given (a session
    given as itself, by its own name)."""
    return [item if isinstance(item, str) else s.name
            for item, s in zip(items, sessions, strict=True)]


def _numbers(param: str, text: str, n: int, positive: bool = False) -> List[float]:
    """The ``n`` comma-separated numbers of ``text`` (finite, and positive
    if ``positive``); InvalidParameter otherwise."""
    try:
        values = [float(v.strip()) for v in text.split(',')]
    except ValueError:
        raise InvalidParameter(f"{param}: give numbers separated by commas"
                               + ("" if param == 'stitch_wavelengths' else ", or 'auto'")) from None
    if len(values) != n:
        raise InvalidParameter(f"{param}: give {n}, one per other session ({len(values)} given)")
    if not all(np.isfinite(values)) or (positive and not all(v > 0 for v in values)):
        raise InvalidParameter(f"{param}: give finite numbers" + (" above 0" if positive else ""))
    return values


def _ranges(param: str, text: str, n: int) -> List[Tuple[float, float]]:
    """The ``n`` ranges ``'min-max'`` of ``text`` (nm), each with min below
    max; InvalidParameter otherwise."""
    ranges = []
    try:
        for r in text.split(','):
            vmin, vmax = map(float, r.split('-'))
            ranges.append((vmin, vmax))
    except ValueError:
        raise InvalidParameter(f"{param}: write them 'min-max, min-max' (in nm), "
                               "or 'auto'") from None
    if len(ranges) != n:
        raise InvalidParameter(f"{param}: give {n}, one per other session ({len(ranges)} given)")
    if any(not vmin < vmax for vmin, vmax in ranges):
        raise InvalidParameter(f"{param}: the start of each range must be below its end")
    return ranges


def _grid_unit(xunit: str) -> au.UnitBase:
    """The unit of the step of a grid: a velocity or a wavelength."""
    try:
        unit = au.Unit(xunit)
    except (ValueError, TypeError):
        raise InvalidParameter(f"xunit: {xunit!r} is not a unit") from None
    if not (unit.is_equivalent(au.km / au.s) or unit.is_equivalent(au.nm)):
        raise InvalidParameter(f"xunit: {xunit!r} is neither a velocity nor a wavelength")
    return unit


def _blue_to_red(arms: list) -> list:
    """The arms ``(name, spectrum)`` from blue to red, by their bluest
    pixel in nm (it was compared in the unit of each arm, F134); an arm
    whose axis is not a wavelength cannot be cut at wavelengths."""
    for name, spec in arms:
        if not au.Unit(spec.x.unit).is_equivalent(au.nm):
            raise RecipeFailed(f"The axis of {name} is in {spec.x.unit}: the arms are cut at "
                               "wavelengths, so convert it to a wavelength with x_convert.")
    return sorted(arms, key=lambda arm: np.nanmin(arm[1].x.to_value(au.nm)))


def _in_axis_unit(value_nm: float, unit) -> float:
    """A wavelength in nm, in the unit of the axis (as it is for nm)."""
    return value_nm if au.Unit(unit) == au.nm else (value_nm * au.nm).to_value(unit)


def _arm_part(arm, name: str, lo: Optional[float], hi: Optional[float], unit):
    """The part of ``arm`` above the cut ``lo`` and up to the cut ``hi`` (nm;
    None: no cut); InvalidParameter if it holds no pixel (it failed with
    "Split expression resulted in an empty spectrum")."""
    terms = []
    if lo is not None:
        terms.append(f"x > {_in_axis_unit(lo, unit)}")
    if hi is not None:
        terms.append(f"x <= {_in_axis_unit(hi, unit)}")
    where = (f"from {lo} to {hi}" if lo is not None and hi is not None
             else f"from {lo}" if lo is not None else f"up to {hi}")
    x = arm.x.value
    inside = np.ones(len(x), dtype=bool)
    if lo is not None:
        inside &= x > _in_axis_unit(lo, unit)
    if hi is not None:
        inside &= x <= _in_axis_unit(hi, unit)
    if not inside.any():
        raise InvalidParameter(f"stitch_wavelengths: {name} has no pixel {where} nm; each "
                               "arm must cover the range between its cuts")
    logger.info(f"{name}: kept {where}.")
    expression = terms[0] if len(terms) == 1 else f"({terms[0]}) & ({terms[1]})"
    return arm.split(expression)


def _equalize_arm(arm, name: str, reference, factor: Optional[float],
                  rng: Optional[Tuple[float, float]], unit):
    """``arm`` scaled to ``reference``: by ``factor`` if given, else to match
    in the range ``rng`` (nm), else in their whole overlap."""
    if factor is not None:
        logger.info(f"{name}: scaled by {factor:.4f} (factor given).")
        return arm.scale_flux(factor)
    if rng is None:
        logger.debug(f"{name}: equalizing on the whole overlap.")
        return arm.equalize_to_reference(reference, order=0)
    vmin, vmax = rng
    logger.debug(f"{name}: equalizing between {vmin} and {vmax} nm.")
    where = f"(x >= {_in_axis_unit(vmin, unit)}) & (x <= {_in_axis_unit(vmax, unit)})"
    ref_overlap = reference.split(where)
    other_overlap = arm.split(where)
    if len(ref_overlap.x) > 5 and len(other_overlap.x) > 5:
        from astrocook.core.spectrum_operations import compute_flux_scaling
        scale_model = compute_flux_scaling(
            ref_overlap.x.value, ref_overlap.y.value,
            other_overlap.x.value, other_overlap.y.value, order=0
        )
        scalar_factor = np.nanmedian(scale_model)
        logger.info(f"{name}: scaled by {scalar_factor:.4f} (between {vmin} and {vmax} nm).")
        return arm.scale_flux(scalar_factor)
    logger.warning(f"{name}: not enough pixels between {vmin} and {vmax} nm: "
                   "equalized on the whole overlap instead.")
    return arm.equalize_to_reference(reference, order=0)


# Old names (phase 6-quater, PR 6): see astrocook/_renamed.py
__getattr__ = renamed(__name__, {"RecipeEditV2": "RecipeEdit"})
