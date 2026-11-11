import astropy.units as au
import json
import logging
import numpy as np
from typing import TYPE_CHECKING, Any, Optional

from astrocook.core.atomic_data import STANDARD_MULTIPLETS, is_multiplet, series_lines
from astrocook.fitting.voigt_fitter import NoPixelsToFit, VoigtFitterV2

from astrocook.fitting.voigt_model import LINK_ERRORS, LINK_FUNCTIONS
from astrocook.recipes.errors import InvalidParameter, MissingPrerequisite, RecipeFailed
from astrocook.recipes._spectrum import spec_of
from astrocook.recipes.params import typed_recipe
from astrocook.recipes.progress import check_cancelled, progress
from astrocook.recipes.result import report, report_removed
from astrocook.recipes.messages import (msg_no_resolution, msg_session_not_found,
                                       msg_z_em_unset)

if TYPE_CHECKING:
    from astrocook.core.session import SessionV2

logger = logging.getLogger(__name__)

# Depth of the friends-of-friends expansion used when fitting the groups
# built by refit_all and metals_auto. These calls asked for depth 0, but
# set_active_components ignored the depth and always used 2 (findings
# register F11). Depth 2 is kept on purpose: with 0 the neighbours stay
# frozen at their initial values, and on the J0103-1305 full chain the fits
# got 2.75 times faster but worse (chi2 +14% in the forest, CIV logN off by
# up to 0.5 dex). Decided with Guido.
FIT_EXPANSION_DEPTH = 2

# --- The text of the methods (phase 6-ter, PR 5c: astrocook.methods_text) ---
from astrocook.methods_text import (fit_bounds, mask_words, region_words, series_list,
                                    series_name)
def _beyond(p):
    return f"up to {p['z_em_margin_kms']} km/s beyond z_em"


def _methods_components_auto(p):
    return (f"The candidate systems of {series_list(p.raw('multiplets'))} were "
            f"searched for{region_words(p.raw('region_limit'))}{p.menu}. The absorbed "
            f"regions (merged within {p['merge_dv']} km/s, at least "
            f"{p['min_pix_region']} pixels wide) were found on {mask_words(p)}, and "
            f"each candidate was scored by the correlation of its lines (accepted for "
            f"R^2 >= {p['score_threshold']}, {_beyond(p)}; Cupani et al. 2026). The "
            f"accepted candidates were fitted with Voigt profiles, with {fit_bounds(p)}")


def _methods_doublets_auto(p):
    return (f"The doublets of {series_list(p.raw('multiplets'))} were searched "
            f"for redward of the Lyman-alpha emission{p.menu}. The absorbed regions "
            f"were found on {mask_words(p)}, and each candidate was scored by the "
            f"correlation of its two lines (accepted for R^2 >= {p['score_threshold']}, "
            f"{_beyond(p)}; Cupani et al. 2026). The accepted candidates were fitted "
            f"with Voigt profiles, with {fit_bounds(p)}")


def _methods_lya_auto(p):
    return (f"Single Lyman-alpha lines were searched for in the forest, between "
            f"Lyman-beta and Lyman-alpha{p.menu}. The absorbed regions were found on "
            f"{mask_words(p)}, and the lines accepted for a score >= "
            f"{p['score_threshold']} ({_beyond(p)}; Cupani et al. 2026) were fitted "
            f"with Voigt profiles, with {fit_bounds(p)}. Lines narrower than "
            f"b = {p['min_b']} km/s were removed as interlopers, and the others refitted")


def _methods_forest_auto(p):
    return (f"Lyman-alpha lines were searched for and fitted each as a Lyman series, "
            f"from Lyman-alpha to the Lyman limit, so that the higher-order lines "
            f"constrain the column density{p.menu}. The lines were accepted for a "
            f"score >= {p['score_threshold']} ({_beyond(p)}; Cupani et al. 2026) and "
            f"fitted with Voigt profiles, with {fit_bounds(p)}; systems narrower than "
            f"b = {p['min_b']} km/s were removed")


def _methods_metals_auto(p):
    return (f"For each system, the ions {series_list(p.raw('ions'))} were tried "
            f"at its redshift, and each was kept if its fit gave log N >= 12 and "
            f"2 <= b <= 50 km/s, neither at the bounds of the fit, and decreased the "
            f"Akaike information criterion by at least {p['aic_penalty']} (Cupani et "
            f"al. 2026){p.menu}; the systems were then refitted, with {fit_bounds(p)}")


def _methods_identify_lines(p):
    for name in ('bypass_scoring', 'debug_rating'):
        p.raw(name)
    return (f"The candidate systems of {series_list(p.raw('multiplets'))} were "
            f"identified in the absorbed regions of the mask {p['mask_col']} (merged "
            f"within {p['merge_dv']} km/s, at least {p['min_pix_region']} pixels wide), "
            f"each scored by the correlation of its lines (accepted for R^2 >= "
            f"{p['score_threshold']}, {_beyond(p)}; Cupani et al. 2026), and written "
            f"in a new session{p.menu}")


def _methods_add_component(p):
    btur = f", b_tur = {p['btur']} km/s" if p.raw('btur') else ""
    return (f"A {series_name(str(p.raw('series')))} component was added at "
            f"z = {p['z']} (starting from log N = {p['logN']}, b = {p['b']} km/s"
            f"{btur}) and fitted with its neighbours, with {fit_bounds(p)}{p.menu}")


_SET_WORDS = (('z', "z = {}"), ('logN', "log N = {}"), ('b', "b = {} km/s"),
              ('btur', "b_tur = {} km/s"), ('series', "the series {}"),
              ('resol', "the resolving power R = {}"))


def _set_parts(p):
    parts = []
    for name, words in _SET_WORDS:
        if name in p and p.given(name):
            value = (series_name(str(p.raw(name))) if name == 'series' else p[name])
            parts.append(words.format(value))
    return ", ".join(parts) or "nothing"


def _methods_update_component(p):
    p.raw('uuid')
    return f"The parameters of a component were set by hand: {_set_parts(p)}{p.menu}"


def _methods_update_components(p):
    n = len(p.raw('uuids') or [])
    return (f"The parameters of {n} component{'s' if n != 1 else ''} were set by "
            f"hand: {_set_parts(p)}{p.menu}")


def _methods_delete_component(p):
    uuids = list(p.raw('uuids') or [])
    if p.given('uuid'):
        uuids.append(p.raw('uuid'))
    n = len(set(uuids))
    if n == 1:
        return f"A component was deleted by hand{p.menu}"
    return f"{n} components were deleted by hand{p.menu}"


def _methods_fit_component(p):
    p.raw('uuid')
    return (f"A component was fitted with its neighbours (group depth "
            f"{p['group_depth']}), with {fit_bounds(p)} and at most {p['max_nfev']} "
            f"evaluations of the model{p.menu}")


def _methods_update_constraint(p):
    p.raw('uuid'), p.raw('target_uuid')
    param = p['param'].replace('logN', 'log N').replace('btur', 'b_tur')
    if p.given('expression'):
        return (f"The parameter {param} of a component was linked to another "
                f"component by the expression {p['expression']}{p.menu}")
    state = "freed" if p.raw('is_free') else "frozen"
    return f"The parameter {param} of a component was {state}{p.menu}"


def _methods_update_constraints_column(p):
    param = p['param'].replace('logN', 'log N').replace('btur', 'b_tur')
    state = "freed" if p.raw('is_free') else "frozen"
    return f"The parameter {param} of all the components was {state}{p.menu}"


def _methods_add_linked_system(p):
    return (f"Components of {series_list(p.raw('series_list'))} were added at "
            f"z = {p['z']}, with their redshifts linked (starting from log N = "
            f"{p['logN']}, b = {p['b']} km/s), and fitted, with {fit_bounds(p)}{p.menu}")


def _methods_optimize_system(p):
    p.raw('uuid')
    return (f"Components were added to a system one at a time where the residuals "
            f"exceeded {p['threshold_sigma']} sigma, at least {p['min_dv']} km/s from "
            f"the others and within +-{p['z_window_kms']} km/s of the system, each kept "
            f"only if it decreased the Akaike information criterion by at least "
            f"{p['aic_penalty']} (at most {p['max_components']} components, stopping "
            f"after {p['patience']} tries without an improvement; group depth "
            f"{p['group_depth']}; Cupani et al. 2026){p.menu}")


def _methods_fit_bayesian(p):
    n = len(p.raw('uuids') or [])
    group = f"a group of {n} component{'s' if n != 1 else ''}" if n else "a group of components"
    if p.raw('engine') == 'nested':
        p.raw('nsteps')
        how = (f"nested sampling (dynesty; Speagle 2020), with {p['nlive']} live "
               f"points")
    else:
        p.raw('nlive')
        how = (f"MCMC (emcee; Foreman-Mackey et al. 2013), with {p['nsteps']} steps "
               f"per walker")
    return (f"The posterior of the parameters of {group} (group depth "
            f"{p['group_depth']}) was sampled with {how}{p.menu}")


def _methods_refit_all(p):
    selected = p.raw('selected_uuids')
    if selected:
        which = (f"The {len(selected)} selected component{'s' if len(selected) != 1 else ''}"
                 f" and their groups were refitted")
    else:
        which = "All the components were refitted"
    region = region_words(p.raw('region_limit'), before=", ")
    return (f"{which}{region} group by group{p.menu}. A group joins the components "
            f"blended with each other, up to {p['group_depth']} steps away, with at "
            f"most {p['max_group_size']} components and {p['max_nfev']} evaluations "
            f"of the model. The fit allowed {fit_bounds(p)}")


def _methods_clean_negligible(p):
    broad = (f"b > {p['max_b']} km/s and log N < 13" if p.raw('combined_check')
             else f"b > {p['max_b']} km/s")
    return (f"The components with log N < {p['min_logN']}, b < {p['min_b']} km/s, or "
            f"{broad} were removed as negligible or unphysical{p.menu}")


ABSORBERS_RECIPES_SCHEMAS = {
    # --- 1. THE MASTER PIPELINE ---
    "components_auto": {
        "methods": _methods_components_auto,
        "references": ("C26", "SciPy"),
        "brief": "Auto-detect and fit components.",
        "requires": ("z_em", "resolution", "continuum"),
        "details": "Finds the absorbed regions, identifies the candidate systems of the given multiplets (Identify lines), adds them as components and fits them all (Refit all); the pipeline behind the other automatic recipes.",
        "params": [
            {"name": "multiplets", "type": str, "default": "CIV,SiIV,MgII,Ly_ab", "doc": "Comma-separated multiplets to search for, e.g. 'CIV,SiIV,MgII,Ly_ab'."},
            {"name": "score_threshold", "type": float, "default": 0.5, "min": 0, "max": 1, "doc": "Minimum score (0-1, the R^2 of the match of the lines of a multiplet) to accept a candidate."},
            {"name": "merge_dv", "type": float, "default": 10.0, "min": 0, "doc": "Absorbed regions closer than this (km/s) are merged."},
            {"name": "min_pix_region", "type": int, "default": 3, "min": 1, "doc": "Minimum width of an absorbed region (pixels)."},
            {"name": "mask_col", "type": str, "default": "abs_mask", "doc": "Column of the mask of the absorbed pixels.", "kind": "column"},
            {"name": "z_window_kms", "type": float, "default": 20.0, "gt": 0, "doc": "Half-width of the window in which the redshift of a component can move in the fit (km/s)."},
            {"name": "max_b", "type": float, "default": 50.0, "gt": 1, "max": 300, "doc": "Upper bound on b in the fits (km/s); components at the bound may be broader and can be refitted."},
            {"name": "max_logN", "type": float, "default": 19.0, "gt": 10, "max": 23, "doc": "Upper bound on logN in the fits; damped systems (logN > 19) are to be fitted by hand."},
            # New generic region limiter
            {"name": "region_limit", "type": str, "default": "None", "doc": "Search and fit only in a region: 'red_side', 'lya_forest', 'all_ly_forest', or 'None' for the whole spectrum.", "choices": ["None", "lya_forest", "all_ly_forest", "red_side"]},
            {"name": "mask", "type": str, "default": "new", "doc": "Mask of the absorbed pixels used to find the regions: 'new' (made again with the defaults of Find absorbed regions) or 'session' (the 'abs_mask' already in the session, e.g. from Auto-estimate continuum with other parameters).", "choices": ["new", "session"]},
            {"name": "z_em_margin_kms", "type": float, "default": 3000.0, "min": 0, "doc": "How far beyond z_em a candidate system can lie (km/s): the associated absorbers and the uncertainty of z_em."},
        ],
        "url": "absorbers_cb.html#components_auto"
    },
    # --- 2. WRAPPERS ---
    "doublets_auto": {
        "methods": _methods_doublets_auto,
        "references": ("C26", "SciPy"),
        "brief": "Auto-fit metal doublets (Red Side).",
        "requires": ("z_em", "resolution", "continuum"),
        "details": "Finds and fits metal doublets (CIV, SiIV, MgII by default) on the red side of the Lyman-alpha emission, where the forest does not confuse them; candidates are kept up to z_em_margin_kms beyond z_em.",
        "params": [
            {"name": "multiplets", "type": str, "default": "CIV,SiIV,MgII", "doc": "Comma-separated doublets to search for."},
            {"name": "score_threshold", "type": float, "default": 0.5, "min": 0, "max": 1, "doc": "Minimum score (0-1, the R^2 of the match of the two lines) to accept a candidate."},
            {"name": "z_window_kms", "type": float, "default": 20.0, "gt": 0, "doc": "Half-width of the window in which the redshift of a component can move in the fit (km/s)."},
            {"name": "max_b", "type": float, "default": 50.0, "gt": 1, "max": 300, "doc": "Upper bound on b in the fits (km/s); components at the bound may be broader and can be refitted."},
            {"name": "max_logN", "type": float, "default": 19.0, "gt": 10, "max": 23, "doc": "Upper bound on logN in the fits; damped systems (logN > 19) are to be fitted by hand."},
            {"name": "mask", "type": str, "default": "new", "doc": "Mask of the absorbed pixels used to find the regions: 'new' (made again with the defaults of Find absorbed regions) or 'session' (the 'abs_mask' already in the session, e.g. from Auto-estimate continuum with other parameters).", "choices": ["new", "session"]},
            {"name": "z_em_margin_kms", "type": float, "default": 3000.0, "min": 0, "doc": "How far beyond z_em a candidate system can lie (km/s): the associated absorbers and the uncertainty of z_em."},
        ],
        "url": "absorbers_cb.html#doublets_auto"
    },
    "lya_auto": {
        "methods": _methods_lya_auto,
        "references": ("C26", "SciPy"),
        "brief": "Auto-fit Ly-alpha forest.",
        "requires": ("z_em", "resolution", "continuum"),
        "details": "Finds and fits single Lyman-alpha lines in the forest (between Ly-beta and Ly-alpha); lines narrower than min_b (metal lines or noise) are removed and the rest refitted.",
        "params": [
            {"name": "score_threshold", "type": float, "default": 0.1, "min": 0, "max": 1, "doc": "Minimum score (0-1) to accept a line; low, as single lines have no partner to match."},
            {"name": "min_b", "type": float, "default": 10.0, "min": 0, "doc": "Lines narrower than this (km/s) are removed as interlopers."},
            {"name": "z_window_kms", "type": float, "default": 100.0, "gt": 0, "doc": "Half-width of the window in which the redshift of a component can move in the fit (km/s)."},
            {"name": "max_b", "type": float, "default": 50.0, "gt": 1, "max": 300, "doc": "Upper bound on b in the fits (km/s); components at the bound may be broader and can be refitted."},
            {"name": "max_logN", "type": float, "default": 19.0, "gt": 10, "max": 23, "doc": "Upper bound on logN in the fits; damped systems (logN > 19) are to be fitted by hand."},
            {"name": "mask", "type": str, "default": "new", "doc": "Mask of the absorbed pixels used to find the regions: 'new' (made again with the defaults of Find absorbed regions) or 'session' (the 'abs_mask' already in the session, e.g. from Auto-estimate continuum with other parameters).", "choices": ["new", "session"]},
            {"name": "z_em_margin_kms", "type": float, "default": 3000.0, "min": 0, "doc": "How far beyond z_em a candidate system can lie (km/s): the associated absorbers and the uncertainty of z_em."},
        ],
        "url": "absorbers_cb.html#lya_auto"
    },
    "forest_auto": {
        "methods": _methods_forest_auto,
        "references": ("C26", "SciPy"),
        "brief": "Auto-fit full Lyman series.",
        "requires": ("z_em", "resolution", "continuum"),
        "details": "Finds Lyman-alpha lines and fits each as a Lyman series (Ly-alpha to the Lyman limit), so that the higher-order lines constrain the column density; systems narrower than min_b are removed.",
        "params": [
            {"name": "score_threshold", "type": float, "default": 0.1, "min": 0, "max": 1, "doc": "Minimum score (0-1) to accept a line."},
            {"name": "min_b", "type": float, "default": 10.0, "min": 0, "doc": "Lyman systems narrower than this (km/s) are removed."},
            {"name": "z_window_kms", "type": float, "default": 20.0, "gt": 0, "doc": "Half-width of the window in which the redshift of a component can move in the fit (km/s)."},
            {"name": "max_b", "type": float, "default": 50.0, "gt": 1, "max": 300, "doc": "Upper bound on b in the fits (km/s); components at the bound may be broader and can be refitted."},
            {"name": "max_logN", "type": float, "default": 19.0, "gt": 10, "max": 23, "doc": "Upper bound on logN in the fits; damped systems (logN > 19) are to be fitted by hand."},
            {"name": "z_em_margin_kms", "type": float, "default": 3000.0, "min": 0, "doc": "How far beyond z_em a candidate system can lie (km/s): the associated absorbers and the uncertainty of z_em."},
        ],
        "url": "absorbers_cb.html#forest_auto"
    },
    "metals_auto": {
        "methods": _methods_metals_auto,
        "references": ("C26", "SciPy"),
        "brief": "Check for associated metals.",
        "requires": ("z_em", "resolution", "continuum"),
        "details": "For each system already in the list, tries to add other metal ions at its redshift and keeps each whose fit gives logN >= 12 and 2 <= b <= 50 km/s, not at the bounds of the fit, and decreases the AIC by the AIC penalty; skips the ions whose lines fall outside the spectrum; then refits all.",
        "params": [
            {"name": "ions", "type": str, "default": "CIV,SiIV,NV,OVI,CII,MgII", "doc": "Comma-separated ions to try, e.g. 'CIV,SiIV,NV,OVI,CII,MgII'."},
            {"name": "aic_penalty", "type": float, "default": 5.0, "min": 0, "doc": "Decrease of the AIC needed to keep a new ion (against the same model without it, on the pixels of its fit)."},
            {"name": "z_window_kms", "type": float, "default": 20.0, "gt": 0, "doc": "Half-width of the window in which the redshift of a component can move in the fit (km/s)."},
            {"name": "max_b", "type": float, "default": 50.0, "gt": 1, "max": 300, "doc": "Upper bound on b in the fits (km/s); components at the bound may be broader and can be refitted."},
            {"name": "max_logN", "type": float, "default": 19.0, "gt": 10, "max": 23, "doc": "Upper bound on logN in the fits; damped systems (logN > 19) are to be fitted by hand."}
        ],
        "url": "absorbers_cb.html#metals_auto"
    },
    "identify_lines": {
        "methods": _methods_identify_lines,
        "references": ("C26",),
        "brief": "Identify likely absorption systems.",
        "requires": ("z_em", "resolution", "continuum"),
        "branching": True,
        "details": "Finds the absorbed regions from the mask, scores how well each multiplet matches them (the R^2 of the correlation of its lines) and keeps the best candidate systems; writes them in a new session, with the region and the score of each.",
        "params": [
            {"name": "multiplets", "type": str, "default": "CIV,SiIV,MgII,Ly_ab", "doc": "Comma-separated multiplets to search for, e.g. 'CIV,SiIV,MgII,Ly_ab'."},
            {"name": "mask_col", "type": str, "default": "abs_mask", "doc": "Column of the mask (True for the absorbed pixels), e.g. 'abs_mask' of Find absorbed regions.", "kind": "column"},
            {"name": "min_pix_region", "type": int, "default": 3, "min": 1, "doc": "Minimum width of an absorbed region (pixels)."},
            {"name": "merge_dv", "type": float, "default": 10.0, "min": 0, "doc": "Absorbed regions closer than this (km/s) are merged."},
            {"name": "score_threshold", "type": float, "default": 0.5, "min": 0, "max": 1, "doc": "Minimum score (0-1, the R^2 of the match of the lines) to accept a candidate."},
            {"name": "bypass_scoring", "type": bool, "default": False, "doc": "Accept every candidate found by the positions of its lines, without the score (for debugging).", "agents": False},
            {"name": "debug_rating", "type": bool, "default": False, "doc": "Show the plots of the scores (for debugging).", "agents": False},
            {"name": "z_em_margin_kms", "type": float, "default": 3000.0, "min": 0, "doc": "How far beyond z_em a candidate system can lie (km/s): the associated absorbers and the uncertainty of z_em."},
        ],
        "url": "absorbers_cb.html#identify_lines"
    },
    "populate_from_identification": {
        "methods": "Components were added at the candidate systems of the identification{menu}",
        "brief": "Populate components from identification.",
        "requires": ("continuum",),
        "details": "Populates the system list using the regions and candidates found by 'identify_lines'.",
        "params": [
            {"name": "region_id_col", "type": str, "default": "abs_ids", "doc": "Column of the region of each pixel, written by Identify lines.", "kind": "column"},
            {"name": "series_map_json", "type": str, "default": "", "doc": "JSON object: region -> series of its candidate."},
            {"name": "z_map_json", "type": str, "default": "", "doc": "JSON object: region -> redshift of its candidate."}
        ],
        "url": "absorbers_cb.html#populate_from_identification",
        "gui_hidden": True
    },
    "add_component": {
        "methods": _methods_add_component,
        "references": ("SciPy",),
        "brief": "Add a single component.",
        "requires": ("resolution", "continuum"),
        "details": "Adds a Voigt component of a transition or multiplet at a redshift and fits it with its neighbours.",
        "params": [
            {"name": "series", "type": str, "default": "Ly_a", "doc": "Transition or multiplet, e.g. 'Ly_a', 'CIV'."},
            {"name": "z", "type": float, "default": 0.0, "gt": -1, "doc": "Redshift."},
            {"name": "logN", "type": float, "default": 13.5, "min": 10, "max": 23, "doc": "Starting column density (log cm^-2)."},
            {"name": "b", "type": float, "default": 10.0, "min": 1, "max": 300, "doc": "Starting Doppler parameter (km/s)."},
            {"name": "btur", "type": float, "default": 0.0, "min": 0, "max": 100, "doc": "Turbulent broadening (km/s)."},
            {"name": "z_window_kms", "type": float, "default": 20.0, "gt": 0, "doc": "Half-width of the window in which the redshift can move in the fit (km/s)."}
        ],
        "url": "absorbers_cb.html#add_component" 
    },
    "update_component": {
        "methods": _methods_update_component,
        "brief": "Update component.",
        "details": "Changes the parameters of a component; a parameter left at None keeps its value.",
        "params": [
            {"name": "uuid", "type": str, "default": "", "doc": "Identifier of the component.", "gui_hidden": True},
            {"name": "z", "type": float, "default": None, "optional": True, "gt": -1, "doc": "Redshift."},
            {"name": "logN", "type": float, "default": None, "optional": True, "min": 10, "max": 23, "doc": "Column density (log cm^-2)."},
            {"name": "b", "type": float, "default": None, "optional": True, "min": 1, "max": 300, "doc": "Doppler parameter (km/s)."},
            {"name": "btur", "type": float, "default": None, "optional": True, "min": 0, "max": 100, "doc": "Turbulent broadening (km/s)."},
            {"name": "series", "type": str, "default": None, "optional": True, "doc": "Transition or multiplet, e.g. 'CIV'."},
            {"name": "resol", "type": float, "default": None, "optional": True, "gt": 0, "doc": "Resolving power R of this component (overrides the spectrum)."}
        ],
        "url": "absorbers_cb.html#update_component",
        "gui_hidden": True
    },
    "update_components": {
        "methods": _methods_update_components,
        "brief": "Batch update components.",
        "details": "Sets the same parameters on several components; a parameter left at None keeps its value.",
        "params": [
            {"name": "uuids", "type": list, "default": None, "gui_hidden": True, "doc": "Identifiers of the components."},
            {"name": "z", "type": float, "default": None, "optional": True, "gt": -1, "doc": "Redshift."},
            {"name": "logN", "type": float, "default": None, "optional": True, "min": 10, "max": 23, "doc": "Column density (log cm^-2)."},
            {"name": "b", "type": float, "default": None, "optional": True, "min": 1, "max": 300, "doc": "Doppler parameter (km/s)."},
            {"name": "btur", "type": float, "default": None, "optional": True, "min": 0, "max": 100, "doc": "Turbulent broadening (km/s)."},
            {"name": "resol", "type": float, "default": None, "optional": True, "gt": 0, "doc": "Resolving power R of these components (overrides the spectrum)."}
        ]
    },
    "delete_component": {
        "methods": _methods_delete_component,
        "brief": "Delete component(s).",
        "details": "Deletes one or more components and recomputes the model.",
        "params": [
            {"name": "uuid", "type": str, "default": None, "optional": True, "gui_hidden": True, "doc": "Identifier of the component to delete."},
            {"name": "uuids", "type": list, "default": None, "optional": True, "gui_hidden": True, "doc": "Identifiers of the components to delete."}
        ]
    },
    "detect_anchor": {
        "methods": "The redshift at which {trans_self} in this session and {trans_sibling} in the session {sibling_name} overlap best was found{menu}",
        "brief": "Detect doublet anchor redshift.",
        "details": "Finds the redshift at which a transition of this session and one of a sibling session (e.g. the two lines of a doublet, in two sessions) overlap best.",
        "params": [
            {"name": "sibling_name", "type": str, "default": "", "doc": "Name of the sibling session.", "kind": "session"},
            {"name": "trans_self", "type": str, "default": "OVI_1031", "doc": "Transition in this session, e.g. 'OVI_1031'."},
            {"name": "trans_sibling", "type": str, "default": "OVI_1037", "doc": "Transition in the sibling session, e.g. 'OVI_1037'."}
        ],
        "url": "absorbers_cb.html#detect_anchor"
    },
    "fit_component": {
        "methods": _methods_fit_component,
        "references": ("SciPy",),
        "brief": "Fit component.",
        "requires": ("resolution", "continuum"),
        "details": "Fits a component with its neighbours (the components blended with it) using the Voigt fitter.",
        "params": [
            {"name": "uuid", "type": str, "default": "", "doc": "Identifier of the component.", "gui_hidden": True},
            {"name": "max_nfev", "type": int, "default": 2000, "min": 1, "doc": "Maximum number of evaluations of the model."},
            {"name": "z_window_kms", "type": float, "default": 20.0, "gt": 0, "doc": "Half-width of the window in which the redshift can move in the fit (km/s)."},
            {"name": "group_depth", "type": int, "default": 2, "min": 0, "doc": "How far the group extends: 0 the component alone, 1 the components blended with it, 2 also those blended with them."} 
        ],
        "url": "absorbers_cb.html#fit_component",
        "gui_hidden": True
    },
    "update_constraint": {
        "methods": _methods_update_constraint,
        "brief": "Update constraint.",
        "details": "Freezes or frees a parameter of a component, or links it to a parameter of another component with an expression.",
        "params": [
            {"name": "uuid", "type": str, "default": "", "doc": "Identifier of the component.", "gui_hidden": True},
            {"name": "param", "type": str, "default": "", "doc": "Parameter: 'z', 'logN', 'b' or 'btur'.", "choices": ["z", "logN", "b", "btur"]},
            {"name": "is_free", "type": bool, "default": None, "optional": True, "doc": "True: free to vary in the fits; False: frozen."},
            {"name": "expression", "type": str, "default": None, "optional": True, "doc": "Expression that links the parameter to another component, e.g. \"p['<uuid of the other component>'].z\" (None: no link); one that cannot be evaluated is refused."},
            {"name": "target_uuid", "type": str, "default": None, "optional": True, "doc": "Identifier of the component it is linked to."} 
        ],
        "url": "absorbers_cb.html#update_constraint",
        "gui_hidden": True
    },
    "update_constraints_column": {
        "methods": _methods_update_constraints_column,
        "brief": "Bulk update constraints.",
        "details": "Freezes or frees a parameter of all the components.",
        "params": [
            {"name": "param", "type": str, "default": "", "doc": "Parameter: 'z', 'logN', 'b' or 'btur'.", "choices": ["z", "logN", "b", "btur"]},
            {"name": "is_free", "type": bool, "default": True, "doc": "True: free to vary in the fits; False: frozen."}
        ],
        "url": "absorbers_cb.html#update_constraints_column",
        "gui_hidden": True
    },
    "add_linked_system": {
        "methods": _methods_add_linked_system,
        "references": ("SciPy",),
        "brief": "Add linked system.",
        "requires": ("resolution", "continuum"),
        "details": "Adds components of several series at the same redshift, with the redshifts linked, and fits them.",
        "params": [
            {"name": "series_list", "type": str, "default": "CIV,SiIV", "doc": "Comma-separated series, e.g. 'CIV,SiIV'."},
            {"name": "z", "type": float, "default": 0.0, "gt": -1, "doc": "Redshift."},
            {"name": "logN", "type": float, "default": 13.5, "min": 10, "max": 23, "doc": "Starting column density (log cm^-2)."},
            {"name": "b", "type": float, "default": 10.0, "min": 1, "max": 300, "doc": "Starting Doppler parameter (km/s)."},
            {"name": "z_window_kms", "type": float, "default": 20.0, "gt": 0, "doc": "Half-width of the window in which the redshift can move in the fit (km/s)."}
        ],
        "url": "absorbers_cb.html#add_linked_system",
        "gui_hidden": True
    },
    "optimize_system": {
        "methods": _methods_optimize_system,
        "references": ("C26", "SciPy"),
        "brief": "Optimize system (Iterative Fit).",
        "requires": ("resolution", "continuum"),
        "details": "Adds components to a system one at a time where the residuals are significant, keeping each only if it improves the fit by the AIC penalty.",
        "params": [
            {"name": "uuid", "type": str, "default": "", "doc": "Identifier of a component of the system.", "gui_hidden": True},
            {"name": "max_components", "type": int, "default": 5, "min": 0, "doc": "Maximum number of components to add."},
            {"name": "threshold_sigma", "type": float, "default": 2.5, "gt": 0, "doc": "Residual (in sigma) above which a new component is tried."},
            {"name": "aic_penalty", "type": float, "default": 5.0, "min": 0, "doc": "Decrease of the AIC needed to keep a new component."},
            {"name": "z_window_kms", "type": float, "default": 100.0, "gt": 0, "doc": "Half-width of the velocity window around the system that is analysed (km/s)."},
            {"name": "min_dv", "type": float, "default": 10.0, "min": 0, "doc": "Minimum separation of a new component from the others (km/s); lower to resolve blends."},
            {"name": "group_depth", "type": int, "default": 2, "min": 0, "doc": "How far the group extends: 0 the component alone, 1 the components blended with it, 2 also those blended with them."},
            {"name": "patience", "type": int, "default": 2, "min": 0, "doc": "Number of tries without an improvement before stopping."}
        ],
        "url": "absorbers_cb.html#optimize_system"
    },
    "fit_bayesian": {
        "methods": _methods_fit_bayesian,
        "references": ("emcee", "dynesty"),
        "brief": "Fit Bayesian.",
        "requires": ("resolution", "continuum"),
        "details": "Samples the posterior of the parameters of a group of components, with MCMC (emcee) or nested sampling (dynesty); needs the [bayes] extra.",
        "params": [
            {"name": "uuids", "type": list, "default": None, "optional": True, "gui_hidden": True, "doc": "Identifiers of the components of the group."},
            {"name": "engine", "type": str, "default": "mcmc", "doc": "Sampler: 'mcmc' (emcee) or 'nested' (dynesty).", "choices": ["mcmc", "nested"]},
            {"name": "nsteps", "type": int, "default": 500, "min": 1, "doc": "Number of steps of each walker (MCMC)."},
            {"name": "nlive", "type": int, "default": 250, "min": 1, "doc": "Number of live points (nested sampling)."},
            {"name": "group_depth", "type": int, "default": 2, "min": 0, "doc": "How far the group extends: 0 the components alone, 1 the components blended with them, 2 also those blended with these."}
        ],
        "gui_hidden": True
    },
    "refit_all": {
        "methods": _methods_refit_all,
        "references": ("SciPy",),
        "brief": "Refit all systems.",
        "requires": ("resolution", "continuum"),
        "details": "Fits all the components, group by group (the components blended with each other), each group on its own.",
        "params": [
            {"name": "max_nfev", "type": int, "default": 200, "min": 1, "doc": "Maximum number of evaluations of the model, per group."},
            {"name": "z_window_kms", "type": float, "default": 100.0, "gt": 0, "doc": "Half-width of the window in which the redshift of a component can move in the fit (km/s)."},
            {"name": "group_depth", "type": int, "default": 2, "min": 0, "doc": "How far a group extends: 0 each component alone, 1 the components blended with it, 2 also those blended with them."},
            {"name": "max_group_size", "type": int, "default": 20, "min": 1, "doc": "Maximum number of components fitted together: larger groups are split."},
            {"name": "max_b", "type": float, "default": 300.0, "gt": 1, "max": 300, "doc": "Upper bound on b in the fits (km/s)."},
            {"name": "max_logN", "type": float, "default": 23.0, "gt": 10, "max": 23, "doc": "Upper bound on logN in the fits (log cm^-2)."},
            {"name": "region_limit", "type": str, "default": "None", "doc": "Fit only the components in a region: 'red_side', 'lya_forest', 'all_ly_forest', or 'None' for all.", "choices": ["None", "lya_forest", "all_ly_forest", "red_side"]},
            {"name": "selected_uuids", "type": list, "default": None, "optional": True, "gui_hidden": True, "doc": "Fit only these components, with their groups (None: all)."}
        ],
        "url": "absorbers_cb.html#refit_all"
    },
    "clean_negligible": {
        "methods": _methods_clean_negligible,
        "brief": "Remove unphysical/negligible components.",
        "details": "Removes the components that are too weak, too narrow, or too broad (and weak) to be real absorbers.",
        "params": [
            {"name": "min_logN", "type": float, "default": 11.5, "doc": "Components with a lower column density (log cm^-2) are removed."},
            {"name": "max_b", "type": float, "default": 80.0, "gt": 0, "doc": "Components broader than this (km/s) are removed (see combined_check)."},
            {"name": "min_b", "type": float, "default": 3.0, "min": 0, "doc": "Components narrower than this (km/s) are removed."},
            {"name": "combined_check", "type": bool, "default": True, "doc": "True: remove a component broader than max_b only if it is also weak (logN < 13); False: remove all of them."}
        ],
        "url": "absorbers_cb.html#clean_negligible"
    },
}

#: Rounds of "remove the Ly-a narrower than min_b, refit the rest" in
#: lya_auto (#517, F21, agreed with Guido): the refit can narrow other
#: lines below min_b, so the filter runs again after it, with one more
#: refit; what is still below min_b after the last refit is removed
#: without a refit. Before, the filter ran once, before the last refit,
#: and 39 Ly-a with b < 10 km/s were left on J0103-1305
MIN_B_ROUNDS = 2


def _region_ident_dict(value) -> dict:
    """The identifications of the regions ({region ID: [[transition,
    score], ...]}), stored as a JSON string or, after an .acs2 round trip, as
    a dict (F27)."""
    if not value:
        return {}
    raw = json.loads(value) if isinstance(value, str) else value
    return {int(float(k)): v for k, v in raw.items()}


def _identification_extra(spec) -> dict:
    """The identification of :meth:`RecipeAbsorbersV2.identify_lines`, for
    the structured result (phase 6-bis, PR 4b; F49: the candidate systems
    keep their region and their score): the candidate systems and the
    regions, with their range in x and their identifications."""
    meta = spec.meta
    series = _region_ident_dict(meta.get('series_map_json'))
    maps = {k: _region_ident_dict(meta.get(f"{k}_map_json"))
            for k in ('z', 'logN', 'b', 'region', 'score')}
    candidates = [{'series': name, 'z': maps['z'].get(sid),
                   'logN_guess': maps['logN'].get(sid), 'b_guess': maps['b'].get(sid),
                   'region': maps['region'].get(sid), 'score': maps['score'].get(sid)}
                  for sid, name in sorted(series.items())]
    regions = []
    idents = _region_ident_dict(meta.get('region_identifications') or {})
    if spec.has_aux_column('abs_ids'):
        ids = np.asarray(spec.get_column('abs_ids').value)
        x = np.asarray(spec.x.value)
        for rid in sorted(int(r) for r in np.unique(ids) if r > 0):
            inside = x[ids == rid]
            regions.append({'region': rid, 'x_min': float(inside.min()),
                            'x_max': float(inside.max()), 'x_unit': str(spec.x.unit),
                            'identifications': [
                                {'transition': t, 'score': float(sc)}
                                for t, sc in sorted(idents.get(rid, []),
                                                    key=lambda ts: -ts[1])]})
    dropped = meta.get('rejected_candidates_json') or '[]'
    dropped = json.loads(dropped) if isinstance(dropped, str) else list(dropped)
    return {'candidates': candidates, 'regions': regions,
            'dropped_candidates': dropped}


def merge_region_ids(old_ids, old_ident, new_ids, new_ident):
    """Merge the regions of a new identification with those of earlier
    steps of a pipeline.

    The old regions are kept; the new ones, with their IDs shifted past the
    largest old ID, fill the pixels that no old region covers. The
    identifications follow the same IDs: the old ones are kept, the new ones
    get the shifted IDs, and only the regions left in the merged column are
    kept (F90).

    Returns
    -------
    tuple
        The merged region IDs (array) and identifications (dict).
    """
    old_ids = np.asarray(old_ids)
    new_ids = np.asarray(new_ids)
    max_old = np.max(old_ids) if len(old_ids) > 0 else 0
    shifted_new = np.zeros_like(new_ids)
    mask_new = new_ids > 0
    shifted_new[mask_new] = new_ids[mask_new] + max_old
    final_ids = old_ids.copy()
    fill_mask = (shifted_new > 0) & (final_ids == 0)
    final_ids[fill_mask] = shifted_new[fill_mask]

    merged = dict(_region_ident_dict(old_ident))
    for rid, cands in _region_ident_dict(new_ident).items():
        merged[int(rid + max_old)] = cands
    present = set(np.unique(final_ids[final_ids > 0]).astype(int).tolist())
    merged = {k: v for k, v in merged.items() if k in present}
    return final_ids, merged


def _without_components(session, uuids):
    """The session without the components ``uuids`` (and the constraints
    on them), with the model recomputed from the remaining ones."""
    from astrocook.fitting.voigt_fitter import VoigtFitterV2
    systs = session.systs.delete_component(uuids=list(uuids))
    session = session.with_new_system_list(systs)
    if systs.components:
        _, model_flux = VoigtFitterV2(session.spec, systs).compute_model_flux()
        session = session.with_new_spectrum(session.spec.update_model(model_flux))
    return session


def _remove_narrow(session, series: str, min_b: float, refit_kwargs: dict):
    """Remove the components of ``series`` narrower than ``min_b`` and refit
    the rest, MIN_B_ROUNDS times; remove without a refit those still
    narrower after the last refit (#517, F21; F65 for ``forest_auto``)."""
    def narrow(s):
        return [c.uuid for c in s.systs.components
                if c.series == series and c.b < min_b] if s.systs else []

    for _ in range(MIN_B_ROUNDS):
        uuids = narrow(session)
        if not uuids:
            break
        logger.info(f"{len(uuids)} lines narrower than {min_b} km/s removed; refitting the others.")
        report_removed(uuids, f"narrower than min_b = {min_b} km/s")
        if len(uuids) == len(session.systs.components):
            # Nothing left to refit (its 0 said so before phase 6-bis, PR 3)
            return _without_components(session, uuids)
        trimmed = session.with_new_system_list(session.systs.delete_component(uuids=uuids))
        # A failure of the refit is a RecipeError, which the cleaning
        # boundary of the pipeline turns into a WARNING (phase 6); before
        # phase 6-bis (PR 3) its 0 removed the lines without a refit, silently
        session = RecipeAbsorbersV2(trimmed).refit_all(**refit_kwargs)
    uuids = narrow(session)
    if uuids:
        logger.info(f"{len(uuids)} lines still narrower than {min_b} km/s after the last refit removed.")
        report_removed(uuids, f"narrower than min_b = {min_b} km/s after the last refit")
        session = _without_components(session, uuids)
    return session


#: The errors of a fit that does not converge or cannot be evaluated
#: (ValueError covers np.linalg.LinAlgError and a start outside the bounds)
FIT_ERRORS = (ValueError, ArithmeticError)


#: A fitted value within this fraction of its bound is at the bound
AT_BOUND_RTOL = 1e-3


def _metals_rejection(comp, delta_aic: float, aic_penalty: float, b_max: float,
                      logN_max: float) -> Optional[str]:
    """Why ``metals_auto`` does not keep a trial ion, or None if it does.

    The thresholds on logN and b are those the recipe always used. The
    values at the bounds of the fit and the AIC (F121, Guido's choice,
    phase 6-ter PR 7-quater): once the ions found were kept (F111), the
    criterion kept 192 ions of 239 on J0103-1305, 67 with b at its bound
    of 50 km/s; a value at its bound is not a measure (the fit wanted
    more), and ``aic_penalty`` was documented and not used.
    """
    if comp.logN < 12.0:
        return "logN below 12"
    if comp.b < 2.0:
        return "b below 2 km/s"
    if comp.b > 50.0:
        return "b above 50 km/s"
    if comp.b >= b_max * (1 - AT_BOUND_RTOL):
        return f"b at its bound ({b_max:g} km/s)"
    if comp.logN >= logN_max * (1 - AT_BOUND_RTOL):
        return f"logN at its bound ({logN_max:g})"
    if delta_aic > -aic_penalty:
        return f"AIC change above -{aic_penalty:g}"
    return None


def _fit_cancelled() -> bool:
    """The stop check of the fitter (called at each evaluation of the
    residuals): raises RecipeCancelled when the recipe was cancelled, so
    that a long fit stops too (phase 6-ter, PR 7-quater)."""
    check_cancelled()
    return False


def _warn_if_unconvolved(spec) -> None:
    """WARNING when a fit runs without a resolution (it was at DEBUG only,
    in the fitter): the profiles are not convolved, so b includes the
    instrumental broadening."""
    from astrocook.core.resolution import resolution_summary
    if not resolution_summary(spec)[0] > 0:
        logger.warning(msg_no_resolution)


def _warn_cleaning_failed(session, series: str, min_b: float, error: Exception):
    """The cleaning of the narrow lines failed: say what the result keeps."""
    n = (sum(1 for c in session.systs.components
             if c.series == series and c.b < min_b) if session.systs else 0)
    logger.warning(f"The lines narrower than {min_b} km/s could not be removed "
                   f"({type(error).__name__}: {error}): {n} of them are in the result.")
    logger.debug("Traceback of the cleaning:", exc_info=error)


def _check_region(spec, region_limit: str) -> None:
    """RecipeFailed if ``region_limit`` cannot be used.

    An unknown region (or ``z_em`` not set) used to be skipped with a
    WARNING, and the pipeline or the refit worked on the whole spectrum
    (phase 6); it was an ERROR and 0 until phase 6-bis (PR 3).
    """
    if region_limit == "None":
        return
    if spec._data.z_em == 0.0:
        # A prerequisite of the region, not of the recipe (PR 5)
        raise MissingPrerequisite(f"The region '{region_limit}' needs the emission "
                                  f"redshift: {msg_z_em_unset}")
    try:
        spec.get_region_bounds(region_limit)
    except ValueError as e:
        raise RecipeFailed(f"The region '{region_limit}' cannot be used: {e}") from e


def _drop_negligible(session):
    """The automatic pipelines remove the components that end at the lower
    bound of logN (10, flag ``logN_at_lower_bound``, #517): they absorb
    nothing measurable, and were 12 of the forest of J0103-1305 (F21)."""
    if not session or not session.systs:
        return session
    uuids = [c.uuid for c in session.systs.components if 'logN_at_lower_bound' in c.flags]
    if not uuids:
        return session
    logger.info(f"Removed {len(uuids)} component{'s' if len(uuids) > 1 else ''} "
                 "at the lower bound of logN (negligible).")
    report_removed(uuids, "logN at the lower bound of the fit (negligible)")
    return _without_components(session, uuids)


#: The flag of the components whose group refit_all could not fit (phase
#: 6-bis, PR 4b): they keep their values from before the fit
NOT_FITTED = 'not_fitted'


def _flag_not_fitted(systs, uuids):
    """The system list with the flag ``not_fitted`` on the components
    ``uuids``."""
    import dataclasses
    from astrocook.core.system_list import SystemListV2
    targets = set(uuids)
    comps = [dataclasses.replace(c, flags=tuple(c.flags or ()) + (NOT_FITTED,))
             if c.uuid in targets and NOT_FITTED not in (c.flags or ()) else c
             for c in systs.components]
    return SystemListV2(dataclasses.replace(systs._data, components=comps))


def _log_flags(components) -> None:
    """Log how many of the fitted ``components`` have flags (#517): a
    parameter at a bound of the fit, or uncertainties not computed."""
    from collections import Counter
    counts = Counter(f for c in components for f in c.flags)
    n = sum(1 for c in components if c.flags)
    if not n:
        return
    listed = ", ".join(f"{flag.replace('_', ' ')}: {k}" for flag, k in sorted(counts.items()))
    logger.warning(f"{n} component{'s' if n > 1 else ''} with fit warnings ({listed}). "
                    "A parameter at a bound has an unreliable value and error: refit "
                    "with wider bounds (e.g. Refit Selected) or check the component.")


class RecipeAbsorbersV2:
    """
    Recipes for identification, fitting, and management of absorption lines.

    These methods are accessed via the ``session.absorbers`` attribute and delegate
    logic to :class:`~astrocook.core.system_list.SystemListV2`,
    :class:`~astrocook.core.spectrum.SpectrumV2`, and the fitting engine.
    """
    def __init__(self, session_v2: 'SessionV2'):
        self._session = session_v2
        self._tag = 'abs'

    @typed_recipe
    def components_auto(self, multiplets: str = "CIV,SiIV,MgII,Ly_ab",
                        score_threshold: float = 0.5,
                        merge_dv: float = 10.0,
                        min_pix_region: int = 3,
                        mask_col: str = "abs_mask",
                        z_window_kms: float = 20.0,
                        region_limit: str = "None",
                        max_b: float = 50.0, max_logN: float = 19.0,
                        mask: str = "new",
                        z_em_margin_kms: float = 3000.0) -> 'SessionV2':
        """
        Auto-detect and fit components.

        The Master Pipeline: Identifies candidates (doublets or singles), populates
        them as systems, and refits the spectrum. Supports masking by region.
        Delegates to :meth:`identify_lines`, :meth:`populate_from_identification`,
        and :meth:`refit_all`.

        Parameters
        ----------
        multiplets : str, optional
            Comma-separated list of multiplets to search for. Defaults to ``"CIV,SiIV,MgII,Ly_ab"``.
        score_threshold : float, optional
            Min R^2 score (0-1) to accept a candidate. Defaults to ``0.5``.
        merge_dv : float, optional
            Max velocity (km/s) to merge adjacent regions. Defaults to ``10.0``.
        min_pix_region : int, optional
            Minimum width in pixels to count as a region. Defaults to ``3``.
        mask_col : str, optional
            Name of the mask column (True=Absorbed). Defaults to ``"abs_mask"``.
        z_window_kms : float, optional
            Velocity window for fitting (km/s). Defaults to ``20.0``.
        region_limit : str, optional
            Limit search to ``'red_side'``, ``'lya_forest'``, or ``'None'``.
            Defaults to ``"None"``.
        max_b : float, optional
            Upper bound on b (km/s) in the fits. It keeps the components from
            broadening into their neighbours, which makes the fits of the
            forest slow and gives spurious broad lines; components that end
            at the bound may be broader and can be refitted with
            :meth:`refit_all` or :meth:`fit_component`. Defaults to ``50.0``.
        max_logN : float, optional
            Upper bound on logN in the fits. In a saturated trough that the
            candidates do not fill, a component can otherwise fill it with
            its damping wings and run away to logN 20-23, blacking out the
            model around it (findings register F62). With 19, the Lyman-limit
            systems stay within reach; damped systems are to be fitted by
            hand. Components that end at the bound can be refitted with
            :meth:`refit_all`. Defaults to ``19.0``.
        mask : str, optional
            The mask of the absorbed pixels used to find the regions:
            ``'new'`` (made again with the defaults of
            :meth:`~astrocook.recipes.continuum.RecipeContinuumV2.find_absorbed`)
            or ``'session'`` (the ``abs_mask`` already in the session, e.g.
            from ``estimate_auto`` with other parameters). Defaults to
            ``'new'`` (F50, phase 6-bis, PR 6b).
        z_em_margin_kms : float, optional
            How far beyond ``z_em`` a candidate system can lie (km/s).
            Defaults to ``3000.0`` (F55, phase 6-bis, PR 6b).

        Returns
        -------
        SessionV2
            A new :class:`~astrocook.core.session.SessionV2` with detected systems.
        """
        progress(logger, 0)
        logger.info(f"Searching for components (region: {region_limit}).")
        _check_region(self._session.spec, region_limit)

        # --- 1. PREPARE REGION MASK ---
        # We create a temporary session with a restricted mask to force identify_lines
        # to ignore specific regions.
        
        session_current = self._session

        # The mask of the user (e.g. from "Auto-estimate continuum" with other
        # parameters) is given back with the result: the fresh mask below,
        # restricted by region_limit, is used only inside the pipeline (F50)
        user_mask = spec_of(session_current)._data.aux_cols.get(mask_col)
        
        if mask == "session":
            # The mask of the session, e.g. of Auto-estimate continuum with
            # other parameters (F50, phase 6-bis, PR 6b)
            if user_mask is None:
                raise MissingPrerequisite(
                    f"The session has no mask '{mask_col}' to use (mask='session'): "
                    "make it with Find absorbed regions or Auto-estimate continuum, "
                    "or use mask='new'.")
            logger.info(f"Using the absorbed regions of the session (column '{mask_col}').")
        else:
            # --- 0. FORCE MASK REGENERATION ---
            # Always generate a fresh mask to ensure it matches current z_em/flux state.
            # This prevents "stale mask" bugs when running recipes sequentially.
            logger.info(f"Finding the absorbed regions (column '{mask_col}').")
            from astrocook.recipes.continuum import RecipeContinuumV2
            rec_cont = RecipeContinuumV2(session_current)
            session_current = rec_cont.find_absorbed(
                smooth_len_lya="5000.0", 
                smooth_len_out="400.0",
                kappa="2.0",
                template="False"
            )
        kept_mask = (user_mask if user_mask is not None
                     else spec_of(session_current)._data.aux_cols[mask_col])

        def with_kept_mask(session):
            spec = session.spec.update_column(
                mask_col, kept_mask.values, unit=kept_mask.unit,
                description=kept_mask.description)
            return session.with_new_spectrum(spec)
        
        # --- 0.5 SMART CLEAN (Refactored) ---
        if region_limit != "None" and session_current.systs:
            spec = spec_of(session_current)
            # Get bounds as Quantities
            obs_min, obs_max = spec.get_region_bounds(region_limit)

            # Delegate logic to SystemListV2
            new_systs = session_current.systs.remove_in_region(obs_min, obs_max)
            kept = {c.uuid for c in new_systs.components}
            report_removed([c.uuid for c in session_current.systs.components
                            if c.uuid not in kept],
                           f"in the region of the pipeline ({region_limit}): replaced by "
                           "its components")
            session_current = session_current.with_new_system_list(new_systs)
        
        # --- 1. PREPARE REGION MASK ---
        target_mask_col = mask_col
        if region_limit != "None":
            # (Checked at the start: it was skipped with a WARNING, and the
            # pipeline worked on the whole spectrum, phase 6)
            obs_min, obs_max = spec_of(session_current).get_region_bounds(region_limit)
            new_spec = spec_of(session_current).mask_region(mask_col, obs_min, obs_max)
            session_current = session_current.with_new_spectrum(new_spec)
        
        # --- 2. IDENTIFY (using the masked session) ---
        logger.info("Identifying the doublets (step 1 of 3).")
        
        # Capture the OLD abs_ids before they are overwritten
        old_abs_ids = None
        old_ident = spec_of(session_current).meta.get('region_identifications')
        if spec_of(session_current).has_aux_column('abs_ids'):
            old_col = spec_of(session_current).get_column('abs_ids')
            assert old_col is not None  # has_aux_column('abs_ids') just said so
            old_abs_ids = old_col.value.copy()
            
        # Create a temp recipe wrapper for the masked session
        rec_temp = RecipeAbsorbersV2(session_current)
        
        session_step1 = rec_temp.identify_lines(
            multiplets=multiplets, 
            mask_col=target_mask_col, # Use the restricted mask
            min_pix_region=min_pix_region,
            merge_dv=merge_dv, 
            score_threshold=score_threshold,
            bypass_scoring="False", 
            debug_rating="False",
            z_em_margin_kms=z_em_margin_kms
        )
        
        progress(logger, 30)

        # --- 2.5 MERGE IDENTIFICATION IDs ---
        if old_abs_ids is not None and session_step1.spec.has_aux_column('abs_ids'):
            # The regions of the earlier steps are kept, the new ones fill
            # the gaps with shifted IDs; the identifications of the
            # regions are merged in the same way (F90: they kept the
            # unshifted IDs of the last step only, and pointed to the
            # wrong regions)
            final_ids, merged_ident = merge_region_ids(
                old_abs_ids, old_ident,
                session_step1.spec.get_column('abs_ids').value,
                session_step1.spec.meta.get('region_identifications'))
            final_spec = session_step1.spec.update_column('abs_ids', final_ids)
            final_spec._data.meta['region_identifications'] = json.dumps(
                {str(k): v for k, v in merged_ident.items()})
            session_step1 = session_step1.with_new_spectrum(final_spec)

        spec_meta = session_step1.spec.meta
        series_map_json = spec_meta.get('series_map_json')
        z_map_json = spec_meta.get('z_map_json')
        
        if not series_map_json:
            logger.warning("No candidates found: no components are added.")
            return self._session

        # --- 3. POPULATE ---
        logger.info("Adding the components (step 2 of 3).")
        rec_step2 = RecipeAbsorbersV2(session_step1)
        session_step2 = rec_step2.populate_from_identification(
            region_id_col='abs_ids',
            series_map_json=series_map_json,
            z_map_json=z_map_json
        )

        progress(logger, 50)
        
        # --- 4. REFIT ALL ---
        logger.info("Fitting the components (step 3 of 3).")
        rec_step3 = RecipeAbsorbersV2(session_step2)
        final_session = rec_step3.refit_all(
            max_nfev='200',
            z_window_kms=z_window_kms,
            group_depth='1',
            max_group_size='25',
            region_limit=region_limit,
            max_b=max_b,
            max_logN=max_logN
        )
        
        # With no components refit_all gives the session back with a
        # WARNING; a fit that fails is a RecipeError, the failure of the
        # pipeline (Guido, phase 6-bis, PR 3: its 0 made the pipeline give
        # back the components unfitted, as a result)

        # Check for stop request on the valid final session

        check_cancelled()
        
        final_session = _drop_negligible(final_session)
        logger.debug("Component pipeline complete.")
        return with_kept_mask(final_session)

    @typed_recipe
    def doublets_auto(self, multiplets: str = "CIV,SiIV,MgII",
                      score_threshold: float = 0.5,
                      z_window_kms: float = 20.0,
                      max_b: float = 50.0, max_logN: float = 19.0,
                      mask: str = "new",
                      z_em_margin_kms: float = 3000.0) -> 'SessionV2':
        """
        Auto-fit metal doublets (Red Side).

        Wrapper for :meth:`components_auto` targeted at the red side (z < z_em).
        Finds CIV, SiIV, MgII.

        Parameters
        ----------
        multiplets : str, optional
            Multiplets to search for. Defaults to ``"CIV,SiIV,MgII"``.
        score_threshold : float, optional
            Min R^2 score (0-1). Defaults to ``0.5``.
        z_window_kms : float, optional
            Fit window (km/s). Defaults to ``20.0``.
        max_b : float, optional
            Upper bound on b (km/s) in the fits, see :meth:`components_auto`.
            Defaults to ``50.0``.
        max_logN : float, optional
            Upper bound on logN in the fits, see :meth:`components_auto`.
            Defaults to ``19.0``.
        mask : str, optional
            The mask of the absorbed pixels used to find the regions:
            ``'new'`` (made again with the defaults of
            :meth:`~astrocook.recipes.continuum.RecipeContinuumV2.find_absorbed`)
            or ``'session'`` (the ``abs_mask`` already in the session, e.g.
            from ``estimate_auto`` with other parameters). Defaults to
            ``'new'`` (F50, phase 6-bis, PR 6b).
        z_em_margin_kms : float, optional
            How far beyond ``z_em`` a candidate system can lie (km/s).
            Defaults to ``3000.0`` (F55, phase 6-bis, PR 6b).

        Returns
        -------
        SessionV2
            A new :class:`~astrocook.core.session.SessionV2` with detected doublets.
        """
        return self.components_auto(
            multiplets=multiplets,
            score_threshold=score_threshold,
            z_window_kms=z_window_kms,
            region_limit="red_side",  # <--- FORCE RED SIDE
            max_b=max_b,
            max_logN=max_logN,
            mask=mask,
            z_em_margin_kms=z_em_margin_kms
        )

    @typed_recipe
    def lya_auto(self, score_threshold: float = 0.1, min_b: float = 10.0,
                 z_window_kms: float = 100.0, max_b: float = 50.0,
                 max_logN: float = 19.0, mask: str = "new",
                 z_em_margin_kms: float = 3000.0) -> 'SessionV2':
        """
        Auto-fit Ly-alpha forest.

        Wrapper for :meth:`components_auto` targeted at the forest. Includes
        filtering of interlopers based on the b-parameter.
        Delegates to :meth:`astrocook.core.system_list.SystemListV2.filter_by_criteria`.

        Parameters
        ----------
        score_threshold : float, optional
            Lower threshold for single lines (R^2). Defaults to ``0.1``.
        min_b : float, optional
            Minimum b parameter (km/s) to keep. Defaults to ``10.0``.
        z_window_kms : float, optional
            Fit window (km/s). Defaults to ``100.0``.
        max_b : float, optional
            Upper bound on b (km/s) in the fits, see :meth:`components_auto`.
            Defaults to ``50.0``.
        max_logN : float, optional
            Upper bound on logN in the fits, see :meth:`components_auto`.
            Defaults to ``19.0``.
        mask : str, optional
            The mask of the absorbed pixels used to find the regions:
            ``'new'`` (made again with the defaults of
            :meth:`~astrocook.recipes.continuum.RecipeContinuumV2.find_absorbed`)
            or ``'session'`` (the ``abs_mask`` already in the session, e.g.
            from ``estimate_auto`` with other parameters). Defaults to
            ``'new'`` (F50, phase 6-bis, PR 6b).
        z_em_margin_kms : float, optional
            How far beyond ``z_em`` a candidate system can lie (km/s).
            Defaults to ``3000.0`` (F55, phase 6-bis, PR 6b).

        Returns
        -------
        SessionV2
            A new :class:`~astrocook.core.session.SessionV2` with the forest model.
        """
        progress(logger, 0)
        logger.info("Fitting the Lyman-alpha forest.")
        
        # 1. Run Pipeline (FORCE region='lya_forest')
        final_session = self.components_auto(
            multiplets="Ly_a",
            score_threshold=score_threshold,
            z_window_kms=z_window_kms,
            mask_col='abs_mask',
            region_limit='lya_forest', # <--- FORCE FOREST
            max_b=max_b,
            max_logN=max_logN,
            mask=mask,
            z_em_margin_kms=z_em_margin_kms
        )
        check_cancelled()
        
        # 2. Cleanup Phase (Interlopers): the Ly-a narrower than min_b are
        # removed and the rest refitted, MIN_B_ROUNDS times; those still
        # narrower after the last refit are removed without a refit (#517)
        progress(logger, 0)
        logger.info(f"Removing the lines narrower than {min_b} km/s (interlopers).")
        
        min_b_val = min_b
            # A boundary of the step (B, Guido's choice, phase 6): the result of
        # the pipeline is kept if the cleaning fails, with a WARNING
        try:
            # The refit is limited to the forest, as the pipeline: before
            # phase 6-bis (PR 4c, F102, Guido's decision) it refitted every
            # component of the session, e.g. the metals of the red side
            final_session = _remove_narrow(final_session, 'Ly_a', min_b_val, dict(
                max_nfev='200', z_window_kms=z_window_kms, group_depth='1',
                max_group_size='25', max_b=max_b, max_logN=max_logN,
                region_limit='lya_forest'))
        except Exception as e:
            _warn_cleaning_failed(final_session, 'Ly_a', min_b_val, e)
        final_session = _drop_negligible(final_session)

        return final_session

    @typed_recipe
    def forest_auto(self, score_threshold: float = 0.1, min_b: float = 10.0,
                    z_window_kms: float = 20.0, max_b: float = 50.0,
                    max_logN: float = 19.0,
                    z_em_margin_kms: float = 3000.0) -> 'SessionV2':
        """
        Auto-fit full Lyman series.

        Pipeline that models the entire Lyman series (Ly-alpha to Ly-limit)
        simultaneously. Uses 'Ly_a' detection but maps it to 'Lyman' series.

        Parameters
        ----------
        score_threshold : float, optional
            Threshold for detection. Defaults to ``0.1``.
        min_b : float, optional
            Minimum b parameter (km/s). Defaults to ``10.0``.
        z_window_kms : float, optional
            Fit window (km/s). Defaults to ``20.0``.
        max_b : float, optional
            Upper bound on b (km/s) in the fits, see :meth:`components_auto`.
            Defaults to ``50.0``.
        max_logN : float, optional
            Upper bound on logN in the fits, see :meth:`components_auto`.
            Defaults to ``19.0``.
        z_em_margin_kms : float, optional
            How far beyond ``z_em`` a candidate system can lie (km/s).
            Defaults to ``3000.0`` (F55, phase 6-bis, PR 6b).

        Returns
        -------
        SessionV2
            A new :class:`~astrocook.core.session.SessionV2` with Lyman series systems.
        """
        logger.info("Fitting the Lyman series in the forest.")
        
        # 1. IDENTIFY 'Ly_a' lines
        # (We can use components_auto logic, but we need to intercept the map before populate)
        # So we manually call steps 1 & 2.
        
        # Force forest region
        session_id_only = self.identify_lines(
            multiplets="Ly_a", 
            score_threshold=score_threshold,
            z_em_margin_kms=z_em_margin_kms,
            # We can't easily pass region mask to identify_lines directly via public API without
            # duplicating the logic in components_auto. 
            # Ideally, identify_lines should support region_limit, but for now, 
            # let's trust the user or just run on whole spectrum (since forest is dominant).
            # OR better: Use components_auto but hacking the 'multiplets' is hard.
        )
        
        # 2. MODIFY MAPS: Swap 'Ly_a' -> 'Lyman'
        import json
        spec_meta = session_id_only.spec.meta
        series_json = spec_meta.get('series_map_json')
        
        if series_json:
            series_map = json.loads(series_json)
            for k, v in series_map.items():
                if v == 'Ly_a': series_map[k] = 'Lyman'
            new_series_json = json.dumps(series_map)
        else:
            new_series_json = None

        # 3. POPULATE & FIT
        rec_step2 = RecipeAbsorbersV2(session_id_only)
        session_step2 = rec_step2.populate_from_identification(
            series_map_json=new_series_json,
            z_map_json=spec_meta.get('z_map_json')
        )
        
        # 4. Refit All
        rec_step3 = RecipeAbsorbersV2(session_step2)
        # max_nfev as before the defaults were aligned with the GUI (#573)
        refit_kwargs = dict(z_window_kms=z_window_kms, group_depth='1',
                            max_nfev='100', max_b=max_b, max_logN=max_logN)
        fitted = rec_step3.refit_all(**refit_kwargs)
        # The Lyman systems narrower than min_b are removed as in lya_auto
        # (F65: before, min_b was taken and not used)
        min_b_val = min_b
        # A boundary of the step, as in lya_auto (phase 6)
        try:
            fitted = _remove_narrow(fitted, 'Lyman', min_b_val, refit_kwargs)
        except Exception as e:
            _warn_cleaning_failed(fitted, 'Lyman', min_b_val, e)
        return _drop_negligible(fitted)

    @typed_recipe
    def metals_auto(self, ions: str = "CIV,SiIV,NV,OVI,CII,MgII",
                    aic_penalty: float = 5.0, z_window_kms: float = 20.0,
                    max_b: float = 50.0, max_logN: float = 19.0) -> 'SessionV2':
        """
        Check for associated metals.

        Scans existing systems and attempts to add other common metal ions
        (e.g., NV, SiIV) at the same redshift, keeping them if fit quality holds.
        Delegates to :meth:`astrocook.core.system_list.SystemListV2.add_component`.

        Parameters
        ----------
        ions : str, optional
            Comma-separated list of ions to check. Defaults to ``"CIV,SiIV,NV,OVI,CII,MgII"``.
        aic_penalty : float, optional
            Decrease of the AIC needed to keep a new ion, against the same
            model without it on the pixels of its fit. An ion is also kept
            only if its fit gives logN >= 12 and 2 <= b <= 50 km/s, with
            neither b nor logN at the bounds of the fit (F121). Each ion
            tried is reported (``tried`` in the extra of the result), with
            the reason when it is not kept. Defaults to ``5.0``.
        z_window_kms : float, optional
            Fit window (km/s) for local island fitting. Defaults to ``20.0``.
        max_b : float, optional
            Upper bound on b (km/s) in the fits, see :meth:`components_auto`.
            Defaults to ``50.0``.
        max_logN : float, optional
            Upper bound on logN in the fits, see :meth:`components_auto`.
            Defaults to ``19.0``.

        Returns
        -------
        SessionV2
            A new :class:`~astrocook.core.session.SessionV2` with added metals.
        """
        ion_list = [i.strip() for i in ions.split(',')]
        z_win = z_window_kms
        b_max = max_b
        logN_max = max_logN
        
        if not self._session.systs or not self._session.systs.components:
            # Nothing to do, not a failure (Guido, phase 6-bis, PR 3)
            logger.warning("There are no systems to search metals at: nothing is added.")
            return self._session

        logger.info(f"Searching {len(ion_list)} ions at the redshifts of the systems.")
        _warn_if_unconvolved(self._session.spec)
        
        current_session = self._session
        # Optimization: Group components by Z
        from collections import defaultdict
        systems_by_z = defaultdict(list)
        for c in current_session.systs.components:
            z_round = round(c.z, 4)
            systems_by_z[z_round].append(c)
        trials = [(z_val, ion) for z_val, comps in systems_by_z.items()
                  for ion in ion_list if ion not in {c.series for c in comps}]

        spec = spec_of(current_session)
        running_systs = current_session.systs
        count_added = 0
        failed_trials = 0
        out_of_spectrum = []
        tried = []
        progress(logger, 0)

        for i, (z_val, ion) in enumerate(trials):
            check_cancelled()
            logger.debug(f"  Checking {ion} at z={z_val}...")

            # A. Add
            trial_systs = running_systs.add_component(ion, z_val, logN=13.0, b=10.0)
            new_uuid = trial_systs.components[-1].uuid

            # B. Fit (Local Island)
            group = trial_systs.get_connected_group([new_uuid], spec=spec)

            try:
                with trial_systs.fitting_context(list(group), group_depth=FIT_EXPANSION_DEPTH,
                                                 spec=spec):
                    fitter = VoigtFitterV2(spec, trial_systs, check_stop=_fit_cancelled)
                    fitted_systs, model_flux, res = fitter.fit(max_nfev=100, z_window_kms=z_win,
                                                               verbose=0, b_max=b_max,
                                                               logN_max=logN_max)
            except NoPixelsToFit:
                # Its lines are outside the spectrum: it cannot be tried, and
                # it was fitted on the whole spectrum, for hours on a long
                # one, and kept at its first values (F111, Guido's choice)
                out_of_spectrum.append(f"{ion} at z = {z_val}")
                continue
            except FIT_ERRORS as e:
                # A trial that cannot be fitted is a candidate rejected;
                # the other errors reach the recipe (phase 6)
                failed_trials += 1
                logger.debug(f"The trial fit of {ion} at z = {z_val} failed: {e}")
                continue
            finally:
                progress(logger, int(100 * (i + 1) / len(trials)))

            # C. Evaluate
            fitted_comp = fitted_systs.get_component_by_uuid(new_uuid)
            if not fitted_comp:
                continue
            # The AIC with the ion against the same model without it, on
            # the pixels of its fit (3 more free parameters)
            chi2_with, n_pix = fitter.chi2_on_fit_pixels(model_flux)
            _, flux_without = VoigtFitterV2(
                spec, fitted_systs.delete_component(uuid=new_uuid)).compute_model_flux()
            chi2_without, _ = fitter.chi2_on_fit_pixels(flux_without)
            delta_aic = chi2_with - chi2_without + 2 * 3
            reason = _metals_rejection(fitted_comp, delta_aic, aic_penalty, b_max, logN_max)
            keep = reason is None
            tried.append({"ion": ion, "z": float(z_val), "logN": round(float(fitted_comp.logN), 3),
                          "b": round(float(fitted_comp.b), 2), "pixels": n_pix,
                          "delta_aic": round(float(delta_aic), 2), "kept": keep,
                          "reason": reason or ""})
            logger.info(f"Tried {ion} at z = {z_val}: logN = {fitted_comp.logN:.2f}, "
                        f"b = {fitted_comp.b:.1f} km/s, AIC change {delta_aic:+.1f} on "
                        f"{n_pix} pixels: {'kept' if keep else 'not kept, ' + str(reason)}.")
            if keep:
                logger.info(f"Found {ion} at z = {z_val} (logN = {fitted_comp.logN:.2f}).")
                running_systs = fitted_systs
                count_added += 1

        report("tried", tried)
        if out_of_spectrum:
            listed = ", ".join(out_of_spectrum[:10])
            more = f" and {len(out_of_spectrum) - 10} more" if len(out_of_spectrum) > 10 else ""
            logger.warning(f"{len(out_of_spectrum)} trial ion{'s' if len(out_of_spectrum) > 1 else ''} "
                           f"not tried, their lines outside the spectrum: {listed}{more}.")
        if failed_trials:
            logger.warning(f"{failed_trials} trial fit{'s' if failed_trials > 1 else ''} "
                           "failed: those candidates were not added.")
        if count_added > 0:
            logger.info(f"Added {count_added} metal components.")
            # The values the pipeline always used, before the defaults of
            # refit_all were aligned with the GUI (#573). The session with
            # the ions found: it refitted the session it started from, and
            # lost them (F111)
            found = current_session.with_new_system_list(running_systs)
            return RecipeAbsorbersV2(found).refit_all(max_b=max_b, max_logN=max_logN,
                                                      group_depth='0', max_nfev='100')
        else:
            logger.info("No metals found.")
            return current_session

    @typed_recipe
    def identify_lines(self, 
                       multiplets: str = "CIV,SiIV,MgII,Ly_ab",
                       mask_col: str = 'abs_mask', 
                       min_pix_region: int = 3,
                       merge_dv: float = 10.0,
                       score_threshold: float = 0.5,
                       bypass_scoring: bool = False,
                       debug_rating: bool = False,
                       z_em_margin_kms: float = 3000.0) -> 'SessionV2':
        """
        Identify likely absorption systems.

        Finds absorption regions, computes correlation signals for multiplets, and
        identifies the most likely candidates using a scoring system.
        Delegates to :meth:`astrocook.core.spectrum.SpectrumV2.identify_lines`.

        Parameters
        ----------
        multiplets : str, optional
            Multiplets to search for (comma-separated). Defaults to ``"CIV,SiIV,MgII,Ly_ab"``.
        mask_col : str, optional
            Column of the mask (True for the absorbed pixels, used to find the absorbed regions). Defaults to ``"abs_mask"``.
        min_pix_region : int, optional
            Minimum width in pixels to count as a region. Defaults to ``3``.
        merge_dv : float, optional
            Max velocity (km/s) to merge adjacent regions. Defaults to ``10.0``.
        score_threshold : float, optional
            Min R^2 score (0 to 1) to accept a doublet candidate. Defaults to ``0.5``.
        bypass_scoring : bool, optional
            If ``True``, accept all kinematic candidates ignoring score.
            Defaults to ``False``.
        debug_rating : bool, optional
            If ``True``, show debug plots. Defaults to ``False``.
        z_em_margin_kms : float, optional
            How far beyond ``z_em`` a candidate system can lie (km/s).
            Defaults to ``3000.0`` (F55, phase 6-bis, PR 6b).

        Returns
        -------
        SessionV2
            A new :class:`~astrocook.core.session.SessionV2` with identification maps
            in metadata.
        """
        multiplet_list = [m.strip() for m in multiplets.split(',') if m.strip()]
        
        logger.debug("identify_lines: Calling spec.identify_lines...")
        new_spec_v2 = spec_of(self._session).identify_lines(
            multiplet_list=multiplet_list,
            mask_col=mask_col,
            min_pix_region=min_pix_region,
            merge_dv=merge_dv,
            score_threshold=score_threshold,
            bypass_scoring=bypass_scoring,
            debug_rating=debug_rating,
            z_em_margin_kms=z_em_margin_kms
        )
        for key, value in _identification_extra(new_spec_v2).items():
            report(key, value)
        return self._session.with_new_spectrum(new_spec_v2)
    @typed_recipe
    def populate_from_identification(self, region_id_col: str = 'abs_ids',
                       series_map_json: str = '', z_map_json: str = '') -> 'SessionV2':
        """
        Populate components from identification.

        Populates the system list using the regions and candidates found by
        :meth:`identify_lines`. Delegates to :meth:`astrocook.core.system_list.SystemListV2.add_component`
        and :meth:`astrocook.core.spectrum.SpectrumV2.update_model`.

        Parameters
        ----------
        region_id_col : str, optional
            Column containing region IDs. Defaults to ``"abs_ids"``.
        series_map_json : str, optional
            JSON string mapping Region ID -> Series Name.
        z_map_json : str, optional
            JSON string mapping Region ID -> Redshift.

        Returns
        -------
        SessionV2
            A new :class:`~astrocook.core.session.SessionV2` with populated systems
            and updated model.
        """
        from astrocook.fitting.voigt_fitter import VoigtFitterV2

        if not series_map_json or not z_map_json:
            raise RecipeFailed("There is no identification to build the components from: "
                         "run Identify Lines first.")
        
        try:
            series_map = json.loads(series_map_json)
            z_map = json.loads(z_map_json)
            series_map = {int(k): v for k, v in series_map.items()}
            z_map = {int(k): v for k, v in z_map.items()}
        except (ValueError, TypeError, AttributeError) as e:
            # ValueError covers json.JSONDecodeError
            raise RecipeFailed(f"The identification cannot be read: {e}") from e

        running_systs = self._session.systs
        spec = spec_of(self._session)

        # 1. Trigger Continuum if needed
        # The continuum, if missing, was provided before the recipe (PR 5)
        cont_sess = self._session
        spec = spec_of(cont_sess)
        
        logger.info(f"Adding {len(series_map)} components from the identification.")

        # 2. Iterate and Build
        for sys_id in sorted(series_map.keys()):
            series = series_map[sys_id]
            z_val = z_map.get(sys_id, 0.0)
            
            # A. Add Single Component (e.g. 'CIV')
            running_systs = running_systs.add_component(series, z_val, logN=13.5, b=15.0)

        # 3. Cleanup and Model Generation (Refactored)
        if running_systs.constraint_model:
            running_systs.constraint_model.set_active_components(None)

        final_fitter = VoigtFitterV2(spec, running_systs)
        _, model_flux = final_fitter.compute_model_flux()
        
        # Delegate to SpectrumV2
        new_spec = spec.update_model(model_flux)
        
        # Remove temporary JSON maps to prevent FITS header overflow
        # These keys are no longer needed once components are populated.
        if hasattr(new_spec, 'remove_meta_keys'):
            new_spec = new_spec.remove_meta_keys([
                'series_map_json', 
                'z_map_json', 
                'logN_map_json', 
                'b_map_json'
            ])

        return self._session.with_new_spectrum(new_spec).with_new_system_list(running_systs)
        
    @typed_recipe
    def add_component(self, series: str = 'Ly_a', z: float = 0.0,
                      logN: float = 13.5, b: float = 10.0, btur: float = 0.0,
                      z_window_kms: float = 20.0) -> 'SessionV2':
        """
        Add a single component.

        Adds a new manual Voigt component at a specific redshift and auto-fits it.
        Delegates to :meth:`astrocook.core.system_list.SystemListV2.add_component`.

        Parameters
        ----------
        series : str, optional
            Transition or Multiplet name (e.g. ``'CIV'``). Defaults to ``"Ly_a"``.
        z : float, optional
            Redshift. Defaults to ``0.0``.
        logN : float, optional
            Column Density (log cm^-2). Defaults to ``13.5``.
        b : float, optional
            Doppler parameter (km/s). Defaults to ``10.0``.
        btur : float, optional
            Turbulent broadening (km/s). Defaults to ``0.0``.
        z_window_kms : float, optional
            Max shift for initial fit (km/s). Defaults to ``20.0``.

        Returns
        -------
        SessionV2
            A new :class:`~astrocook.core.session.SessionV2` with the added component.
        """
        z_f, logN_f, b_f, btur_f, z_win_f = z, logN, b, btur, z_window_kms

        # A custom multiplet ("FeII_2586,FeII_2600") carries its transitions in
        # its name (atomic_data.series_lines): it was added to
        # STANDARD_MULTIPLETS here, for the whole process (F12, phase 6-bis PR 7)
        if is_multiplet(series) and series not in STANDARD_MULTIPLETS:
            logger.info(f"Custom multiplet {series} ({', '.join(series_lines(series))}).")

        if self._session.systs is None:
            from astrocook.core.structures import SystemListDataV2
            from astrocook.core.system_list import SystemListV2
            empty_systs = SystemListV2(SystemListDataV2())
            new_systs = empty_systs.add_component(series, z_f, logN_f, b_f, btur_f)
        else:
            new_systs = self._session.systs.add_component(series, z_f, logN_f, b_f, btur_f)
            
        added_comp = new_systs.components[-1]
        temp_session = self._session.with_new_system_list(new_systs)
        temp_recipe = RecipeAbsorbersV2(temp_session)
            
        logger.info(f"Fitting the new {added_comp.series} component.")
        final_session = temp_recipe.fit_component(added_comp.uuid, z_window_kms=str(z_win_f))
        return final_session
    @typed_recipe
    def update_component(self, uuid: str, z: Optional[float] = None, logN: Optional[float] = None,
                         b: Optional[float] = None, btur: Optional[float] = None,
                         series: Optional[str] = None, resol: Optional[float] = None) -> 'SessionV2':
        """
        Update component parameters.

        Modifies the physical parameters of an existing component.
        Delegates to :meth:`astrocook.core.system_list.SystemListV2.update_component`.

        Parameters
        ----------
        uuid : str
            Component UUID.
        z : float or None, optional
            Redshift. Defaults to None (unchanged).
        logN : float or None, optional
            Column Density (log). Defaults to None (unchanged).
        b : float or None, optional
            Doppler parameter (km/s). Defaults to None (unchanged).
        btur : float or None, optional
            Turbulent broadening (km/s). Defaults to None (unchanged).
        series : str, optional
            Series name. Defaults to ``'None'``.
        resol : float or None, optional
            Resolution (R or FWHM). Defaults to None (unchanged).

        Returns
        -------
        SessionV2
            A new :class:`~astrocook.core.session.SessionV2` with updated component.
        """
        changes: dict[str, Any] = {}
        # None: unchanged ('None' in the dialogs and the log)
        if z is not None: changes['z'] = z
        if logN is not None: changes['logN'] = logN
        if b is not None: changes['b'] = b
        if btur is not None: changes['btur'] = btur
        if series is not None: changes['series'] = series
        if resol is not None: changes['resol'] = resol
            
        new_systs = self._session.systs.update_component(uuid, **changes)
        return self._session.with_new_system_list(new_systs)
    @typed_recipe
    def update_components(self, uuids: list, z: Optional[float] = None,
                          logN: Optional[float] = None, b: Optional[float] = None,
                          btur: Optional[float] = None, series: str = 'None',
                          resol: Optional[float] = None) -> 'SessionV2':
        """
        Batch update component parameters.

        Modifies the physical parameters of multiple existing components simultaneously.
        Delegates to :meth:`astrocook.core.system_list.SystemListV2.update_components`.

        Parameters
        ----------
        uuids : list
            List of component UUIDs to update.
        z : float or None, optional
            Redshift. Defaults to None (unchanged).
        logN : float or None, optional
            Column Density (log). Defaults to None (unchanged).
        b : float or None, optional
            Doppler parameter (km/s). Defaults to None (unchanged).
        btur : float or None, optional
            Turbulent broadening (km/s). Defaults to None (unchanged).
        series : str, optional
            Series name. Defaults to ``'None'``.
        resol : float or None, optional
            Resolution (R or FWHM). Defaults to None (unchanged).

        Returns
        -------
        SessionV2
            A new :class:`~astrocook.core.session.SessionV2` with updated components.
        """
        # Parse values
        changes: dict[str, Any] = {}
        # None: unchanged ('None' in the dialogs and the log); series is
        # not in the registry and stays a string, 'None' for unchanged
        if z is not None: changes['z'] = z
        if logN is not None: changes['logN'] = logN
        if b is not None: changes['b'] = b
        if btur is not None: changes['btur'] = btur
        if series not in (None, 'None'): changes['series'] = str(series)
        if resol is not None: changes['resol'] = resol
            
        if not changes: return self._session

        # Delegate to the new bulk method
        new_systs = self._session.systs.update_components(uuids, **changes)
        return self._session.with_new_system_list(new_systs)
            
    @typed_recipe
    def delete_component(self, uuid: Optional[str] = None, uuids: Optional[list] = None) -> 'SessionV2':
        """
        Delete component(s).

        Removes specified components and updates the model. Without a
        continuum the model cannot be computed: the components are deleted
        anyway and the ``model`` column, if any, is removed.
        Delegates to :meth:`astrocook.core.system_list.SystemListV2.delete_component`
        and :meth:`astrocook.core.spectrum.SpectrumV2.update_model`.

        Parameters
        ----------
        uuid : str, optional
            UUID of a single component to delete.
        uuids : list, optional
            List of UUIDs to delete.

        Returns
        -------
        SessionV2
            A new :class:`~astrocook.core.session.SessionV2` with components removed.
        """
        from astrocook.fitting.voigt_fitter import VoigtFitterV2
        
        # 1. Delegate to SystemListV2 (which we updated to handle lists)
        # We pass both arguments; SystemListV2 handles the set logic.
        new_systs = self._session.systs.delete_component(uuid=uuid, uuids=uuids)  # type: ignore[arg-type]  # the core annotates them `str`/`list` with a None default (core/: PR 4)
        kept = {c.uuid for c in new_systs.components}
        report_removed([c.uuid for c in self._session.systs.components if c.uuid not in kept],
                       "deleted")
        spec = spec_of(self._session)

        # 2. Without a continuum the model cannot be computed: delete the
        # components anyway and drop the model column, which would still
        # include them.
        if spec.norm is None:
            logger.warning("There is no continuum: the model is not recomputed, and "
                           "its column is removed.")
            if spec.has_aux_column('model'):
                spec = spec.remove_columns(['model'])
            return self._session.with_new_spectrum(spec).with_new_system_list(new_systs)

        # 3. Recompute Model (Only once!)
        fitter = VoigtFitterV2(spec, new_systs)
        _, model_flux = fitter.compute_model_flux()
            
        # 4. Update Spectrum
        new_spec = spec.update_model(model_flux)
            
        return self._session.with_new_spectrum(new_spec).with_new_system_list(new_systs)
            
    @typed_recipe
    def detect_anchor(self, sibling_name: str, trans_self: str, trans_sibling: str) -> float:
        """
        Detect doublet anchor redshift.

        Finds the redshift of maximum overlap between this session and a sibling
        session (useful for doublet identification across orders/exposures).
        Delegates to :meth:`astrocook.core.spectrum.SpectrumV2.detect_doublet_z`.

        Parameters
        ----------
        sibling_name : str
            Name of the sibling session.
        trans_self : str
            Transition in this session (e.g., ``'OVI_1031'``).
        trans_sibling : str
            Transition in sibling session (e.g., ``'OVI_1037'``).

        Returns
        -------
        float
            The detected redshift.
        """
        sibling_session = None
        # Without the GUI (a script) there are no other sessions: it failed
        # with an AttributeError, hidden by a broad except until phase 6-bis
        gui = getattr(self._session, '_gui', None)
        for hist in getattr(gui, 'session_histories', None) or []:
            if hist.display_name == sibling_name:
                sibling_session = hist.current_state
                break

        if not sibling_session:
            raise RecipeFailed(msg_session_not_found.format(name=sibling_name))
            
        z_found = spec_of(self._session).detect_doublet_z(
            sibling_session.spec, trans_self, trans_sibling
        )
        logger.info(f"Anchor redshift: z = {z_found:.5f}.")
        return z_found
            
    @typed_recipe
    def fit_component(self, uuid: str, max_nfev: int = 2000, z_window_kms: float = 20.0,
                      group_depth: int = 2) -> 'SessionV2':
        """
        Fit component.

        Fits the selected component (and its connected neighbors) using the Voigt engine.
        Delegates to :meth:`astrocook.core.system_list.SystemListV2.fitting_context`
        and :class:`astrocook.fitting.voigt_fitter.VoigtFitterV2`.

        Parameters
        ----------
        uuid : str
            Target Component UUID.
        max_nfev : int, optional
            Max iterations for the optimizer. Defaults to ``2000``.
        z_window_kms : float, optional
            Max velocity shift allowed during fit (km/s). Defaults to ``20.0``.
        group_depth : int, optional
            Grouping depth (1=Neighbors, 2=FoF). Defaults to ``2``.

        Returns
        -------
        SessionV2
            A new :class:`~astrocook.core.session.SessionV2` with updated parameters
            and model.
        """
        max_nfev_i, z_win_f, depth_i = max_nfev, z_window_kms, group_depth

        if not self._session.systs or not self._session.systs.components:
            raise RecipeFailed("There are no components to fit.")

        current_session = self._session

        _warn_if_unconvolved(current_session.spec)
        with current_session.systs.fitting_context([uuid], group_depth=depth_i,
                                                   spec=current_session.spec):
            fitter = VoigtFitterV2(spec_of(current_session), current_session.systs,
                                   check_stop=_fit_cancelled)
            new_systs, model_flux, res = fitter.fit(max_nfev=max_nfev_i, z_window_kms=z_win_f)

        if new_systs.constraint_model:
            new_systs.constraint_model.set_active_components(None)

        # Delegate to SpectrumV2
        new_spec = spec_of(current_session).update_model(model_flux)
            
        return current_session.with_new_spectrum(new_spec).with_new_system_list(new_systs)
            
    @typed_recipe
    def update_constraint(self, uuid: str, param: str, is_free: Optional[bool] = None, 
                          expression: Optional[str] = None, target_uuid: Optional[str] = None) -> 'SessionV2':
        """
        Update constraint.

        Freeze, unfreeze, or link component parameters.
        Delegates to :meth:`astrocook.core.system_list.SystemListV2.update_constraint`.

        Parameters
        ----------
        uuid : str
            Component UUID.
        param : str
            Parameter name (``z``, ``logN``, ``b``).
        is_free : bool, optional
            Is the parameter free to vary?
        expression : str, optional
            Math expression for linking (e.g., ``"p['uuid'].z"``).
        target_uuid : str, optional
            UUID of the source component (for linking tracking).

        Returns
        -------
        SessionV2
            A new :class:`~astrocook.core.session.SessionV2` with updated constraints.
        """
        new_systs = self._session.systs.update_constraint(
            uuid, param, is_free, expression, target_uuid
        )
            
        if expression and not is_free:
            class CompProxy:
                def __init__(self, c): self.c = c
                @property
                def z(self): return self.c.z
                @property
                def logN(self): return self.c.logN
                @property
                def b(self): return self.c.b
                @property
                def btur(self): return self.c.btur

            comp_map = {c.uuid: CompProxy(c) for c in new_systs.components}
            # The functions of the fit (F98: 'sqrt' was offered here and
            # not by the fit, which then ignored such a link)
            eval_globals = {'p': comp_map, **LINK_FUNCTIONS}
                
            try:
                new_val = float(eval(expression, eval_globals))
                logger.info(f"Linked: {param} set to {new_val:.4f}.")
                new_systs = new_systs.update_component(uuid, **{param: new_val})
            except LINK_ERRORS as e:
                # The link is checked by SystemListV2.update_constraint:
                # this should not happen (phase 6)
                logger.warning(f"The linked {param} could not be computed ({e}): it "
                               "keeps its value until the next fit.")

        return self._session.with_new_system_list(new_systs)
            
    @typed_recipe
    def update_constraints_column(self, param: str, is_free: bool) -> 'SessionV2':
        """
        Bulk update constraints.

        Updates the constraint status (free/frozen) for a specific parameter
        across ALL components in the session.

        Parameters
        ----------
        param : str
            Parameter name (``z``, ``logN``, ``b``, ``btur``).
        is_free : bool
            Set to ``True`` (free) or ``False`` (frozen).

        Returns
        -------
        SessionV2
            A new :class:`~astrocook.core.session.SessionV2` with updated constraints.
        """
        if not self._session.systs or not self._session.systs.components:
            raise RecipeFailed("There are no components.")

        current_systs = self._session.systs
        # Iterate over all components and apply the update_constraint logic
        # Since update_constraint returns a NEW system list, we chain the updates.
        
        # Optimization: We can modify the map directly if we access the underlying dict,
        # but let's use the API loop to be safe, modifying the wrapper in place.
        
        running_systs = current_systs
        
        for c in current_systs.components:
            # We explicitly pass None for expression/target to clear links if setting to free,
            # or preserve them if just freezing? 
            # Usually "Freeze All" implies breaking links or just fixing values.
            # The `update_constraint` logic handles "is_free=True/False" by updating the flag.
            
            # Note: If is_free=True, update_constraint clears links.
            # If is_free=False, it keeps links if they exist, or just freezes value.
            running_systs = running_systs.update_constraint(
                c.uuid, param, is_free=is_free
            )
            
        return self._session.with_new_system_list(running_systs)

    @typed_recipe
    def add_linked_system(self, series_list: str, z: float, logN: float = 13.5, 
                          b: float = 10.0, z_window_kms: float = 20.0) -> 'SessionV2':
        """
        Add linked system.

        Add multiple components (e.g. CIV, SiIV) and link their redshifts to the
        first component in the list.

        Parameters
        ----------
        series_list : str
            Comma-separated list of series (e.g., ``"CIV,SiIV"``).
        z : float
            Redshift.
        logN : float, optional
            Column Density. Defaults to ``13.5``.
        b : float, optional
            Doppler parameter. Defaults to ``10.0``.
        z_window_kms : float, optional
            Fit window. Defaults to ``20.0``.

        Returns
        -------
        SessionV2
            A new :class:`~astrocook.core.session.SessionV2` with linked components
            fitted.
        """
        series_names = [s.strip() for s in series_list.split(',') if s.strip()]
        if not series_names:
            raise InvalidParameter("series_list: give at least one series")
        
        current_session = self._session
        
        # 1. Add Primary (The Anchor)
        primary_series = series_names[0]
        logger.info(f"Adding {primary_series} at z = {z:.5f}.")
        
        # Use existing add_component recipe logic to handle creation
        # (This handles empty list creation etc.)
        # But we need the UUID, so we call the SystemList API directly
        
        if current_session.systs is None:
            from astrocook.core.structures import SystemListDataV2
            from astrocook.core.system_list import SystemListV2
            new_systs = SystemListV2(SystemListDataV2())
        else:
            new_systs = current_session.systs
            
        # Add Primary
        new_systs = new_systs.add_component(primary_series, z, logN, b)
        primary_uuid = new_systs.components[-1].uuid
        
        # 2. Add Secondaries (Linked)
        secondary_uuids = []
        
        for s_name in series_names[1:]:
            logger.info(f"Adding {s_name}, linked to it.")
            new_systs = new_systs.add_component(s_name, z, logN, b)
            sec_uuid = new_systs.components[-1].uuid
            secondary_uuids.append(sec_uuid)
            
            # Link Z to Primary
            expr = f"p['{primary_uuid}'].z"
            new_systs = new_systs.update_constraint(
                sec_uuid, 'z', is_free=False, expression=expr, target_uuid=primary_uuid
            )
            
        # 3. Create Temp Session & Fit
        # We fit the Primary (which pulls in the linked secondaries via grouping)
        temp_session = current_session.with_new_system_list(new_systs)
        temp_recipe = RecipeAbsorbersV2(temp_session)
        
        logger.info("Fitting the linked components.")
        # Fit centered on Primary
        final_session = temp_recipe.fit_component(primary_uuid, z_window_kms=str(z_window_kms))
        
        return final_session

    @typed_recipe
    def optimize_system(self, uuid: str, max_components: int = 5,
                        threshold_sigma: float = 2.5, aic_penalty: float = 5.0,
                        z_window_kms: float = 100.0, min_dv: float = 10.0,
                        group_depth: int = 2, patience: int = 2) -> 'SessionV2':
        """
        Optimize system (Iterative Fit).

        Iteratively adds components to a system to minimize residuals using the
        Akaike Information Criterion (AIC).
        Delegates to :meth:`astrocook.core.system_list.SystemListV2.optimize_hierarchy`.

        Parameters
        ----------
        uuid : str
            Target Component UUID (seed).
        max_components : int, optional
            Max extra components to try. Defaults to ``5``.
        threshold_sigma : float, optional
            Residual threshold (sigma) to trigger new component. Defaults to ``2.5``.
        aic_penalty : float, optional
            Min AIC improvement to accept new component (the AIC already
            counts 2 per free parameter; 5 asks for substantial support).
            Defaults to ``5.0``.
        z_window_kms : float, optional
            Velocity window (km/s) to analyze. Defaults to ``100.0``.
        min_dv : float, optional
            Min separation (km/s). Defaults to ``10.0``.
        group_depth : int, optional
            Grouping depth. Defaults to ``2``.
        patience : int, optional
            Trials allowed without AIC improvement. Defaults to ``2``.

        Returns
        -------
        SessionV2
            A new :class:`~astrocook.core.session.SessionV2` with optimized systems.
        """
        max_c, thresh, aic_p = max_components, threshold_sigma, aic_penalty
        z_win, min_dv_f, depth, pat = z_window_kms, min_dv, group_depth, patience

        if not self._session.systs or not self._session.systs.components:
            raise RecipeFailed("There are no components to optimize.")

        current_session = self._session

        _warn_if_unconvolved(current_session.spec)

        trace: dict[str, Any] = {}
        new_systs = current_session.systs.optimize_hierarchy(
            spec=spec_of(current_session),
            uuid_seed=uuid,
            max_components=max_c,
            threshold_sigma=thresh,
            aic_penalty=aic_p,
            z_window_kms=z_win,
            min_dv=min_dv_f,
            group_depth=depth,
            patience=pat,
            trace=trace
        )
        # What the search tried and why it stopped, for the result (F118)
        report('optimization', trace)
            
        if new_systs.constraint_model:
            new_systs.constraint_model.set_active_components(None)

        from astrocook.fitting.voigt_fitter import VoigtFitterV2
        fitter = VoigtFitterV2(spec_of(current_session), new_systs)
        _, model_flux = fitter.compute_model_flux()
            
        # Delegate to SpectrumV2
        new_spec = spec_of(current_session).update_model(model_flux)
            
        return current_session.with_new_spectrum(new_spec).with_new_system_list(new_systs)
                        
    @typed_recipe
    def fit_bayesian(self, uuids: Optional[list[str]] = None, engine: str = 'mcmc',
                     nsteps: int = 500, nlive: int = 250,
                     group_depth: int = 2) -> 'SessionV2':
        """
        Run Bayesian sampling (MCMC or Nested) on a group of components.

        Parameters
        ----------
        uuids : list of str, optional
            A list of component UUIDs to fit. If not provided, a single component 
            must be specified via other means (not supported here).
        engine : {'mcmc', 'nested'}, optional
            The sampling engine to use (default: 'mcmc').
        nsteps : int, optional
            Number of MCMC steps (default: 500).
        nlive : int, optional
            Number of live points for Nested Sampling (default: 250).
        group_depth : int, optional
            Depth of connected group search (FoF) to define the fitting group (default: 2).

        Raises
        ------
        ImportError
            If the sampler of ``engine`` is not installed (optional
            dependencies: ``pip install "astrocook[bayes]"``).
        """
        if not uuids:
            raise RecipeFailed("No components are selected.")
        if isinstance(uuids, str): uuids = [uuids] # Fallback for single uuid string

        from ..fitting.bayesian_fitter import BayesianVoigtFitter, require_bayes_engine
        # Raised, not turned into a return value of 0 like the other
        # failures: the user must know what to install
        require_bayes_engine(str(engine))

        nsteps_i, nlive_i, depth_i = nsteps, nlive, group_depth
            
        if not self._session.systs or not self._session.systs.components:
            raise RecipeFailed("There are no components to fit.")
        current_session = self._session
        _warn_if_unconvolved(current_session.spec)

        with current_session.systs.fitting_context(uuids, group_depth=depth_i,
                                                   spec=current_session.spec):
            b_fitter = BayesianVoigtFitter(current_session.spec, current_session.systs)
                
            if engine.lower() == 'nested':
                logger.debug(f"Running nested sampling (nlive={nlive_i}).")
                res = b_fitter.run_nested(nlive=nlive_i)
            else:
                logger.debug(f"Running MCMC (nsteps={nsteps_i}).")
                res = b_fitter.run_mcmc(nsteps=nsteps_i)

            # Never overwrite the components with non-finite values
            if not np.all(np.isfinite(res['medians'])):
                raise RecipeFailed("The sampler returned non-finite parameters: the "
                             "components are not changed.")

            # 1. Update component parameters with medians
            # We do this INSIDE the context so the fitter's constraints match the sampling results
            new_systs = b_fitter.apply_results(res)
                
            # 2. Store samples in metadata for later corner plot
            new_meta = dict(new_systs._data.meta)
            if 'bayesian_results' not in new_meta:
                new_meta['bayesian_results'] = {}
                
            # We store results for ALL active components in the context
            # so the corner plot is accessible from any of them.
            constraints = b_fitter.fitter._constraints
            assert constraints is not None  # set by the fitter from the system list
            group_results = {
                'samples': res['samples'],
                'labels': constraints.get_active_param_labels(),
                'engine': engine
            }
                
            active_uuids = constraints._active_uuids
            for active_uuid in active_uuids:
                new_meta['bayesian_results'][active_uuid] = group_results
                
            from dataclasses import replace
            final_systs_data = replace(new_systs._data, meta=new_meta)
            from ..core.system_list import SystemListV2
            final_systs = SystemListV2(final_systs_data)

        # 3. Update Model flux
        from ..fitting.voigt_fitter import VoigtFitterV2
        fitter = VoigtFitterV2(spec_of(current_session), final_systs)
        _, model_flux = fitter.compute_model_flux()
        new_spec = spec_of(current_session).update_model(model_flux)

        return current_session.with_new_spectrum(new_spec).with_new_system_list(final_systs)

    @typed_recipe
    def refit_all(self, max_nfev: int = 200, z_window_kms: float = 100.0,
                  group_depth: int = 2, max_group_size: int = 20,
                  region_limit: str = "None",
                  selected_uuids: Optional[list] = None,
                  max_b: float = 300.0, max_logN: float = 23.0) -> 'SessionV2':
        """
        Refit all systems.

        Iteratively fits all disjoint groups of components (islands) in the spectrum.
        Uses chunking for large groups.
        Delegates to :meth:`astrocook.core.system_list.SystemListV2.get_connected_group`
        and :class:`astrocook.fitting.voigt_fitter.VoigtFitterV2`.

        Parameters
        ----------
        max_nfev : int, optional
            Max iterations per group. Defaults to ``200``.
        z_window_kms : float, optional
            Fit window around lines (km/s). Defaults to ``100.0``.
        group_depth : int, optional
            Grouping depth (0=Iterative Single, 1=Neighbors, 2=FoF).
            Defaults to ``2``, as the GUI dialog (#573, F11).
        max_group_size : int, optional
            Max components per fit group (chunking). Defaults to ``20``.
        region_limit : str, optional
            Limit fitting to ``'red_side'`` or ``'lya_forest'``. Defaults to ``"None"``.
        max_b : float, optional
            Upper bound on b (km/s). Defaults to ``300.0``, the physical
            limit of the constraint model; the automatic pipelines pass
            ``"50.0"`` (see :meth:`components_auto`).
        max_logN : float, optional
            Upper bound on logN. Defaults to ``23.0``, the physical limit
            of the constraint model; the automatic pipelines pass ``"19.0"``
            (see :meth:`components_auto`).

        Returns
        -------
        SessionV2
            A new :class:`~astrocook.core.session.SessionV2` with all systems refitted.
        """
        max_nfev_i, z_win_f, depth_i = max_nfev, z_window_kms, group_depth
        max_size_i, b_max_f, logN_max_f = max_group_size, max_b, max_logN

        if not self._session.systs or not self._session.systs.components:
            # Nothing to do, not a failure (Guido, phase 6-bis, PR 3): the
            # pipelines that find no candidates go on with it
            logger.warning("There are no components to fit.")
            return self._session

        # The continuum, if missing, was provided before the recipe (PR 5)
        current_session = self._session
        _warn_if_unconvolved(current_session.spec)

        # --- START NEW CONSOLIDATED LOGIC FOR BETA.4 ---
        # If selected_uuids is provided (e.g., from 'Refit Selected' in UI),
        # we fit the entire connected group once instead of looping.
        if selected_uuids:
            logger.info(f"Refitting {len(selected_uuids)} component{'s' if len(selected_uuids) > 1 else ''}.")
                
            # Use our new spatial-index-powered grouping logic
            combined_active_set = current_session.systs.get_connected_group(
                selected_uuids, max_depth=depth_i, spec=current_session.spec
            )

            # Execute a single fit pass for the entire 'Super-Group'
            with current_session.systs.fitting_context(list(combined_active_set), group_depth=FIT_EXPANSION_DEPTH,
                                                       spec=current_session.spec):
                fitter = VoigtFitterV2(spec_of(current_session), current_session.systs, 
                                      check_stop=_fit_cancelled)
                new_systs, model_flux, res = fitter.fit(
                    max_nfev=max_nfev_i, z_window_kms=z_win_f, b_max=b_max_f,
                    logN_max=logN_max_f
                )
                    
            # Update model and return immediately
            new_spec = spec_of(current_session).update_model(model_flux)
            _log_flags([c for c in new_systs.components
                        if c.uuid in combined_active_set])
            return current_session.with_new_spectrum(new_spec).with_new_system_list(new_systs)
        # --- END NEW CONSOLIDATED LOGIC ---

        # --- 0. DETERMINE FILTER BOUNDS ---
        obs_min_limit = 0.0 * au.nm
        obs_max_limit = 1e6 * au.nm
        use_region_filter = False

        if region_limit != "None":
            # (It fitted the whole spectrum, with a WARNING, phase 6)
            _check_region(current_session.spec, region_limit)
            obs_min_limit, obs_max_limit = spec_of(current_session).get_region_bounds(region_limit)
            use_region_filter = True
            logger.info(f"Refit limited to {region_limit} ({obs_min_limit:.1f} to {obs_max_limit:.1f}).")

        # --- 1. IDENTIFY VALID COMPONENTS (Filtering Step) ---
        from astrocook.core.atomic_data import ATOM_DATA
            
        uuids_in_region = set()
            
        # Prepare selection set
        target_set = set(selected_uuids) if selected_uuids else None

        for c in current_session.systs.components:
            # Check selection first
            if target_set is not None and c.uuid not in target_set:
                continue

            is_valid = True
            if use_region_filter:
                # Resolve primary rest wavelength
                series = c.series
                lam_rest = None
                if series in ATOM_DATA:
                    lam_rest = ATOM_DATA[series]['wave']
                elif is_multiplet(series):
                    prim = series_lines(series)[0]
                    if prim in ATOM_DATA:
                        lam_rest = ATOM_DATA[prim]['wave']
                    
                if lam_rest:
                    lam_obs_ang = lam_rest * (1 + c.z)
                    lam_obs_q = lam_obs_ang * au.Angstrom
                    if not (obs_min_limit <= lam_obs_q <= obs_max_limit):
                        is_valid = False
                
            if is_valid:
                uuids_in_region.add(c.uuid)

        if not uuids_in_region:
            logger.info(f"There are no components to fit in {region_limit}.")
            return current_session

        # --- 2. BUILD WORK UNITS (Grouping Step) ---
        work_units = []

        # 1. Define Work Units
        if depth_i == 0:
            logger.debug("Refit All: Strategy = Iterative Single-Component (Robust).")
            # Filter the list, then sort by Z
            valid_comps = [c for c in current_session.systs.components if c.uuid in uuids_in_region]
            sorted_comps = sorted(valid_comps, key=lambda c: c.z)
            work_units = [[c.uuid] for c in sorted_comps]
        else:
            logger.debug(f"Refit All: Strategy = Connected Groups (depth={depth_i}).")
            # The groups of the components in the region; their neighbours
            # outside it can be pulled in (moved to SystemListV2.fit_groups
            # in phase 6-ter, PR 3: the MCP server shows the same groups)
            work_units = current_session.systs.fit_groups(
                uuids_in_region, max_depth=depth_i, spec=current_session.spec,
                max_size=max_size_i)

        total_units = len(work_units)
        logger.info(f"Fitting {total_units} groups of components.")
        progress(logger, 0)

        # 2. Sequential Local Fit Loop
        reference_systs = current_session.systs
        updated_components_map = {} 

        import dataclasses
        from astrocook.core.system_list import SystemListV2

        failed_units = []
        not_fitted: list[Any] = []  # the components of the groups that could not be fitted
        for i, target_uuids in enumerate(work_units):

            progress(logger, int(100 * (i + 1) / total_units))
                
            if target_uuids:
                peek_comp = next((c for c in reference_systs.components if c.uuid == target_uuids[0]), None)
                if peek_comp:
                    series_str = peek_comp.series
                    z_str = f"{peek_comp.z:.4f}"
                    logger.info(f"Fitting group {i+1} of {total_units}: {series_str} at z = {z_str}.")
                else:
                    logger.info(f"Fitting group {i+1} of {total_units}.")

            if not target_uuids: continue
                
            # A. Define Local Context Window
            target_comps = [c for c in reference_systs.components if c.uuid in target_uuids]
            if not target_comps: continue
                
            min_z = min(c.z for c in target_comps)
            max_z = max(c.z for c in target_comps)
                
            c_kms = 299792.458
            dz_buffer = (1 + max_z) * (1000.0 / c_kms)
            z_start = min_z - dz_buffer
            z_end = max_z + dz_buffer
                
            local_comps = [c for c in reference_systs.components if z_start <= c.z <= z_end]
            local_systs_data = dataclasses.replace(reference_systs._data, components=local_comps)
            local_systs = SystemListV2(local_systs_data)
                
            try:
                fitter = VoigtFitterV2(spec_of(current_session), local_systs)
                    
                with local_systs.fitting_context(target_uuids, group_depth=FIT_EXPANSION_DEPTH,
                                                 spec=current_session.spec):
                    new_local_systs, _, res = fitter.fit(max_nfev=max_nfev_i, z_window_kms=z_win_f, verbose=0,
                                                                 b_max=b_max_f, logN_max=logN_max_f)
                    
                for fitted_comp in new_local_systs.components:
                    if fitted_comp.uuid in target_uuids:
                        updated_components_map[fitted_comp.uuid] = fitted_comp
                            
            except FIT_ERRORS as e:
                # The group keeps its values; reported once below, the
                # other errors reach the recipe (phase 6)
                failed_units.append(min(c.z for c in target_comps))
                not_fitted.extend(c.uuid for c in target_comps)
                logger.debug(f"The fit of group {i + 1} failed: {e}")

        if failed_units:
            listed = ", ".join(f"{z:.5f}" for z in failed_units[:10])
            more = f" and {len(failed_units) - 10} more" if len(failed_units) > 10 else ""
            logger.warning(f"{len(failed_units)} of {total_units} groups could not be "
                           f"fitted (z = {listed}{more}): their components keep "
                           "their previous values.")

        # 3. Batch Merge
        progress(logger, 100)
        logger.debug(f"Refit All: Merging {len(updated_components_map)} fitted components...")
            
        final_comps = []
        for c in reference_systs.components:
            if c.uuid in updated_components_map:
                final_comps.append(updated_components_map[c.uuid])
            else:
                final_comps.append(c)
            
        new_systs_data = dataclasses.replace(reference_systs._data, components=final_comps)
        final_systs = SystemListV2(new_systs_data)

        # 3b. Uncertainties on the final model (issue #534): each unit was
        # fitted while its neighbours still had their initial values, so
        # the errors are recomputed at the final parameters (one Jacobian
        # per unit, no refit)
        final_systs = self._final_uncertainties(current_session.spec, final_systs,
                                                work_units, z_win_f)
        if not_fitted:
            # Flagged (phase 6-bis PR 4b, agreed with Guido in phase 6): the
            # GUI shows them with the other flags, and a later fit clears it
            final_systs = _flag_not_fitted(final_systs, not_fitted)
            report('not_fitted', list(not_fitted))

        # 4. Final Global Model
        if final_systs.constraint_model:
            final_systs.constraint_model.set_active_components(None)
            
        final_fitter = VoigtFitterV2(spec_of(current_session), final_systs)
        _, model_flux = final_fitter.compute_model_flux()

        # Delegate to SpectrumV2
        new_spec = spec_of(current_session).update_model(model_flux)
        _log_flags(final_systs.components)
            
        return current_session.with_new_spectrum(new_spec).with_new_system_list(final_systs)

    def _final_uncertainties(self, spec, systs, work_units, z_window_kms):
        """Recompute the uncertainties of each work unit on the final model."""
        import dataclasses
        from astrocook.core.system_list import SystemListV2
        from astrocook.fitting.voigt_fitter import VoigtFitterV2
        c_kms = 299792.458
        updated = {}
        for target_uuids in work_units:
            targets = [c for c in systs.components if c.uuid in target_uuids]
            if not targets:
                continue
            dz_buffer = (1 + max(c.z for c in targets)) * (1000.0 / c_kms)
            z_lo = min(c.z for c in targets) - dz_buffer
            z_hi = max(c.z for c in targets) + dz_buffer
            local = SystemListV2(dataclasses.replace(
                systs._data, components=[c for c in systs.components if z_lo <= c.z <= z_hi]))
            try:
                with local.fitting_context(list(target_uuids), group_depth=FIT_EXPANSION_DEPTH,
                                           spec=spec):
                    evaluated = VoigtFitterV2(spec, local).evaluate_uncertainties(z_window_kms)
                for c in evaluated.components:
                    if c.uuid in target_uuids:
                        updated[c.uuid] = c
            except (ValueError, ArithmeticError, KeyError, TypeError) as e:
                # The errors of a fit (ValueError covers LinAlgError)
                logger.warning(f"The uncertainties of a group could not be recomputed at "
                               f"the final parameters ({e}): it keeps those of its fit.")
        comps = [updated.get(c.uuid, c) for c in systs.components]
        return SystemListV2(dataclasses.replace(systs._data, components=comps))

    @typed_recipe
    def clean_negligible(self, min_logN: float = 11.5, max_b: float = 80.0,
                         min_b: float = 3.0, combined_check: bool = True) -> 'SessionV2':
        r"""
        Clean negligible components.

        Removes unphysical or negligible components from the system list. This is useful 
        for cleaning up results after automated fitting procedures, which may produce 
        artifacts such as extremely broad, shallow lines (ghosts) or narrow noise spikes.

        The cleaning criteria are:
        1. **Narrow Spikes**: Components with ``b < min_b``.
        2. **Weak Lines**: Components with ``logN < min_logN``.
        3. **Broad/Ghost Lines**: Components with ``b > max_b``.
        
        If ``combined_check`` is ``True``, broad lines are only removed if they 
        are *also* weak (``logN < 13.0``). This protects strong, broad features 
        like DLA wings or OVI absorbers.

        Parameters
        ----------
        min_logN : float, optional
            Minimum column density (log cm\ :sup:`-2`) to keep. Defaults to ``11.5``.
        max_b : float, optional
            Maximum Doppler parameter (km/s). Components broader than this may be removed. 
            Defaults to ``80.0``.
        min_b : float, optional
            Minimum Doppler parameter (km/s). Components narrower than this are always removed. 
            Defaults to ``3.0``.
        combined_check : bool, optional
            If ``True``, broad components (``b > max_b``) are only removed if they are 
            also weak (``logN < 13.0``). If ``False``, all broad components are removed. 
            Defaults to ``True``.

        Returns
        -------
        SessionV2
            A new :class:`~astrocook.core.session.SessionV2` with the flagged components removed.
        """
        lim_N, lim_b_max, lim_b_min = min_logN, max_b, min_b
        use_combined = combined_check

        if not self._session.systs or not self._session.systs.components:
            # Nothing to do, not a failure (Guido, phase 6-bis, PR 3)
            logger.warning("There are no components to clean.")
            return self._session

        to_remove = []
        for c in self._session.systs.components:
            # 1. Unphysically Narrow (Noise spikes) -> Always Remove
            if c.b < lim_b_min:
                logger.debug(f"{c.series} at z = {c.z:.4f} removed: b = {c.b:.1f} km/s, too narrow.")
                to_remove.append(c.uuid)
                report_removed([c.uuid], f"b < min_b = {lim_b_min} km/s (too narrow)")
                continue

            # 2. Too Weak (Negligible)
            if c.logN < lim_N:
                logger.debug(f"{c.series} at z = {c.z:.4f} removed: logN = {c.logN:.2f}, too weak.")
                to_remove.append(c.uuid)
                report_removed([c.uuid], f"logN < min_logN = {lim_N} (too weak)")
                continue

            # 3. Too Broad (Continuum drifts / Ghosts)
            if c.b > lim_b_max:
                if use_combined:
                    # Only remove if ALSO weak (e.g. logN < 13.0)
                    # This protects strong, broad features (like DLA wings or broad OVI)
                    if c.logN < 13.0: 
                        logger.debug(f"{c.series} at z = {c.z:.4f} removed: b = {c.b:.1f} km/s and "
                                     f"logN = {c.logN:.2f}, broad and weak.")
                        to_remove.append(c.uuid)
                        report_removed([c.uuid], f"b > max_b = {lim_b_max} km/s and "
                                                 "logN < 13 (broad and weak)")
                else:
                    # Unconditional removal
                    logger.debug(f"{c.series} at z = {c.z:.4f} removed: b = {c.b:.1f} km/s, too broad.")
                    to_remove.append(c.uuid)
                    report_removed([c.uuid], f"b > max_b = {lim_b_max} km/s (too broad)")

        if to_remove:
            logger.info(f"{len(to_remove)} components removed (too narrow, weak or broad).")
            # Use existing delete method which handles model updates
            return self.delete_component(uuids=to_remove)
        else:
            logger.info("No components to remove.")
            return self._session