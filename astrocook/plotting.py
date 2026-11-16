"""Plots of a session without the GUI (phase 6-ter, PR 3): for the MCP
server, which shows them to an agent as images, and for scripts (e.g. in a
notebook).

Matplotlib's :class:`~matplotlib.figure.Figure` is used directly, without
pyplot: no backend and no display are needed, so the plots work in the
headless Docker image of the server (no PySide6). The colours of the plots
of the GUI are defined here (moved from ``gui/plot_style.py``, which
imports them), so that the two look alike.
"""
import io
from typing import Optional

import astropy.units as au
import numpy as np

from astrocook.core.atomic_data import ATOM_DATA, is_multiplet, series_lines
from astrocook.core.system_list import get_all_w_rest

#: The speed of light (km/s)
C_KMS = 299792.458

#: The residuals panel shows the residuals in sigma between these limits (as
#: the residuals of the velocity plot of the System Inspector)
RESIDUALS_YLIM = 5.0

#: The error of the flux (band, line or bars), in every plot
ERROR_COLOR = '#aaaaaa'

#: The ticks of the HI components on the main plot: a softer red than pure
#: red, drawn translucent (the standard palette, ``PLOT_PALETTES``)
HI_TICK_COLOR = '#d62728'

#: The palettes of the plots, chosen in View > Plot Colours (phase 5-bis):
#: "standard", the colours of 2.2.0 with a thinner model (lw 1.6 -> 1.1: the
#: thick model covered the data in the forest), and "colorblind", from the
#: Okabe-Ito palette, whose blue model, dark spectrum and orange selection
#: stay apart for every type of colour blindness (in the standard palette,
#: model, selection and residuals of the velocity plot become one olive for
#: a deuteranope). The selection (``theme.SELECTION``, an amber close to
#: the orange of Okabe-Ito) is the same in both: it is shared with the
#: widgets. Compared on the "My QSO" session with simulations of colour
#: blindness (Machado et al. 2009) before the choice (Guido).
PLOT_PALETTES = {
    "standard": {
        "flux": "#1f77b4", "model": "#2ca02c", "model_lw": 1.1, "model_alpha": 0.9,
        "cursor": "#e377c2", "residuals": "#ff0000", "hi_ticks": HI_TICK_COLOR,
        "flag_bound": "#9467bd", "flag_nan": "#d62728",
    },
    "colorblind": {
        "flux": "#333333", "model": "#0072b2", "model_lw": 1.1, "model_alpha": 0.9,
        "cursor": "#56b4e9", "residuals": "#d55e00", "hi_ticks": "#d55e00",
        "flag_bound": "#cc79a7", "flag_nan": "#d55e00",
    },
}


#: The continuum, as in the main plot of the GUI
CONTINUUM_COLOR = "black"


def plot_spectrum(session, x_min: Optional[float] = None, x_max: Optional[float] = None,
                  normalized: bool = False, palette: str = "standard",
                  residuals: bool = False):
    """A figure of the spectrum of a session in a range of x.

    The flux and its error, the continuum and the model (if any), and a
    tick at each line of each component, labelled with its series (the HI
    series in the colour of the HI ticks of the GUI).

    Parameters
    ----------
    session : Session
    x_min, x_max : float, optional
        The range, in the unit of x; None: the start (end) of the spectrum.
    normalized : bool, optional
        Divide by the continuum, if there is one.
    palette : str, optional
        ``'standard'`` or ``'colorblind'`` (:data:`PLOT_PALETTES`).
    residuals : bool, optional
        Add a panel below with the residuals, ``(flux - model) / error``
        (in sigma, so the same with or without ``normalized``; F158).

    Returns
    -------
    matplotlib.figure.Figure

    Raises
    ------
    ValueError
        If the range is empty or outside the spectrum, or the residuals are
        asked for and the spectrum has no model.
    """
    from matplotlib.figure import Figure

    colours = PLOT_PALETTES.get(palette, PLOT_PALETTES["standard"])
    spec = session.spec
    data = spec.columns_in_range(["x", "y", "dy"], x_min, x_max)
    x = data["columns"]["x"]
    lo, hi = x.min(), x.max()
    sel = (spec.x.value >= lo) & (spec.x.value <= hi)
    cont = spec.cont.value[sel] if spec.cont is not None else None
    norm = cont if (normalized and cont is not None) else 1.0
    if residuals and not spec.has_aux_column("model"):
        raise ValueError("The spectrum has no model: the residuals need one "
                         "(fit or add components first)")
    if residuals:
        fig = Figure(figsize=(10, 4.8), dpi=100)
        ax, ax_res = fig.subplots(2, 1, sharex=True, gridspec_kw={
            "height_ratios": [3, 1]})
    else:
        fig = Figure(figsize=(10, 3.6), dpi=100)
        ax, ax_res = fig.add_subplot(), None
    ax.step(x, data["columns"]["y"] / norm, where="mid", lw=0.7, color=colours["flux"],
            label="flux")
    ax.step(x, data["columns"]["dy"] / norm, where="mid", lw=0.6, color=ERROR_COLOR,
            label="error")
    if cont is not None:
        ax.plot(x, cont / norm, lw=0.8, color=CONTINUUM_COLOR, label="continuum")
    if spec.has_aux_column("model"):
        ax.plot(x, spec.get_column("model").value[sel] / norm, lw=colours["model_lw"],
                alpha=colours["model_alpha"], color=colours["model"], label="model")
    y = data["columns"]["y"] / norm
    top = np.nanpercentile(y, 99) if np.isfinite(y).any() else 1.0
    comps = session.systs.components if session.systs else []
    for c in comps:
        colour = colours["hi_ticks"] if c.series.startswith("Ly") else colours["model"]
        for w in get_all_w_rest(c.series):
            xo = (w * (1 + c.z) * au.AA).to_value(spec.x.unit, equivalencies=au.spectral())
            if lo <= xo <= hi:
                ax.axvline(xo, ymin=0.9, ymax=0.97, color=colour, lw=0.8)
                ax.text(xo, top * 1.12, c.series, fontsize=6, rotation=90,
                        ha="center", va="bottom", color=colour)
    ax.set_xlim(lo, hi)
    ax.set_ylabel("normalized flux" if np.ndim(norm) else f"flux ({spec.y.unit})")
    ax.legend(fontsize=7, loc="lower left")
    if ax_res is not None:
        _draw_residuals(ax_res, x, residuals_sigma(
            data["columns"]["y"], data["columns"]["dy"],
            spec.get_column("model").value[sel]), colours)
        ax.tick_params(labelbottom=False)
        ax_res.set_xlabel(f"x ({spec.x.unit})")
    else:
        ax.set_xlabel(f"x ({spec.x.unit})")
    fig.tight_layout(h_pad=0.4)
    return fig


def residuals_sigma(y, dy, model):
    """The residuals ``(y - model) / dy``, in sigma (NaN where ``dy`` is 0).

    The same as the residuals of the velocity plot of the System Inspector;
    dividing flux, model and error by the continuum leaves them unchanged.
    """
    y, dy, model = (np.asarray(a, dtype=float) for a in (y, dy, model))
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.divide(y - model, dy, out=np.full_like(y, np.nan), where=dy != 0)


def _draw_residuals(ax, x, resid, colours):
    """The panel of the residuals: the band of 1 sigma, the residuals and the
    zero line, between -5 and +5 sigma."""
    ax.fill_between(x, -1.0, 1.0, step="mid", color=ERROR_COLOR, alpha=0.15)
    ax.step(x, resid, where="mid", color=colours["residuals"], lw=0.8)
    ax.axhline(0, color=colours["model"], lw=1.0, alpha=0.5)
    ax.set_ylim(-RESIDUALS_YLIM, RESIDUALS_YLIM)
    ax.set_ylabel("(F-model)/dF", fontsize=8)


# --- Transitions in velocity ---------------------------------------------------

#: The most transitions of one plot (panels): more would be unreadable
MAX_PANELS = 12


def line_velocity(w_rest, z, w_rest_ref, z_ref):
    """The velocity (km/s) of the line of rest wavelength ``w_rest`` at
    redshift ``z`` from the line ``w_rest_ref`` at redshift ``z_ref``.

    The one formula of the ticks of the velocity plots, the System
    Inspector's among them (for the same transition it is
    ``c (z - z_ref) / (1 + z_ref)``).
    """
    lam_ref = w_rest_ref * (1 + z_ref)
    return C_KMS * (w_rest * (1 + z) - lam_ref) / lam_ref


def component_lines(series):
    """The transitions (names with a wavelength in the atomic data) of the
    series of a component."""
    names = series_lines(series) if is_multiplet(series) else [series]
    return [n for n in names if n in ATOM_DATA]


def tick_velocities(components, transition, z, v_min, v_max):
    """The lines of the components that fall in the velocity window of a
    panel: ``(v, component)`` for each, ``v`` from ``transition`` at ``z``."""
    w_ref = ATOM_DATA[transition]["wave"]
    ticks = []
    for c in components:
        for name in component_lines(c.series):
            v = line_velocity(ATOM_DATA[name]["wave"], c.z, w_ref, z)
            if v_min <= v <= v_max:
                ticks.append((v, c))
    return ticks


def velocity_transitions(session, z, transitions=None, half_width_kms=300.0,
                         max_panels=MAX_PANELS):
    """The transitions of a plot in velocity, and those left out.

    Parameters
    ----------
    session : Session
    z : float
        The redshift at v = 0.
    transitions : list of str, optional
        Names of transitions (``'Ly_a'``, ``'CIV_1548'``) or of multiplets
        (``'CIV'``, expanded); None: the transitions of the components of the
        session that have a line within the window around ``z``.
    half_width_kms : float, optional
        Half of the window.
    max_panels : int, optional

    Returns
    -------
    plotted : list of str
        In the order given (by rest wavelength if chosen here).
    skipped : dict of str to str
        The transitions with no pixel in the window, and why.

    Raises
    ------
    ValueError
        An unknown transition, nothing to plot or too many transitions.
    """
    spec = session.spec
    unit = spec.x.unit
    if transitions is None:
        found = {}
        for c in (session.systs.components if session.systs else []):
            for name in component_lines(c.series):
                if abs(line_velocity(ATOM_DATA[name]["wave"], c.z,
                                     ATOM_DATA[name]["wave"], z)) <= half_width_kms:
                    found[name] = ATOM_DATA[name]["wave"]
        names = sorted(found, key=found.get)
        if not names:
            raise ValueError(
                f"No component has a line within {half_width_kms:g} km/s of z = {z:g}: "
                "give the transitions")
    else:
        names, unknown = [], []
        for t in transitions:
            for name in (series_lines(t) if is_multiplet(t) else [t]):
                (names if name in ATOM_DATA else unknown).append(name)
        if unknown:
            raise ValueError(f"Unknown transitions: {', '.join(unknown)}")
        names = list(dict.fromkeys(names))
        if not names:
            raise ValueError("No transitions given")
    plotted, skipped = [], {}
    x = spec.x.value
    for name in names:
        lam_0 = (ATOM_DATA[name]["wave"] * au.AA).to_value(unit, equivalencies=au.spectral())
        v = C_KMS * (x - lam_0 * (1 + z)) / (lam_0 * (1 + z))
        if np.any(np.abs(v) <= half_width_kms):
            plotted.append(name)
        else:
            skipped[name] = "outside the spectrum"
    if not plotted:
        raise ValueError("No transition has pixels in the window: "
                         + ", ".join(f"{n} ({r})" for n, r in skipped.items()))
    if len(plotted) > max_panels:
        raise ValueError(f"{len(plotted)} transitions: at most {max_panels} per plot; "
                         "give fewer")
    return plotted, skipped


def plot_velocity(session, z, transitions=None, half_width_kms=300.0,
                  normalized=True, residuals=False, palette="standard"):
    """A figure of several transitions in velocity around a redshift.

    One panel per transition, with ``v = 0`` at ``z``, laid out in 1, 2 or 3
    columns: the flux and its error, the model (if any), the continuum
    (if not normalized), the levels 0 and 1 (if normalized) and, at the bottom
    of each panel, a tick for each line of each component in the window (HI in
    the colour of the HI ticks). The panels are those of
    :func:`velocity_transitions`.

    Parameters
    ----------
    session : Session
    z : float
    transitions : list of str, optional
        As in :func:`velocity_transitions`; None: those of the components.
    half_width_kms : float, optional
    normalized : bool, optional
        Divide by the continuum, if there is one.
    residuals : bool, optional
        Under each panel, the residuals in sigma (needs a model).
    palette : str, optional

    Returns
    -------
    matplotlib.figure.Figure

    Raises
    ------
    ValueError
        As :func:`velocity_transitions`, or the residuals are asked for
        and the spectrum has no model.
    """
    from matplotlib.figure import Figure

    spec = session.spec
    if residuals and not spec.has_aux_column("model"):
        raise ValueError("The spectrum has no model: the residuals need one "
                         "(fit or add components first)")
    names, _ = velocity_transitions(session, z, transitions, half_width_kms)
    colours = PLOT_PALETTES.get(palette, PLOT_PALETTES["standard"])
    unit = spec.x.unit
    x_all = spec.x.value
    y_all, dy_all = spec.y.value, spec.dy.value
    cont = spec.cont.value if spec.cont is not None else None
    model = spec.get_column("model").value if spec.has_aux_column("model") else None
    norm = cont if (normalized and cont is not None) else np.ones_like(y_all)
    with np.errstate(divide="ignore", invalid="ignore"):
        y_n, dy_n = y_all / norm, dy_all / norm
        model_n = model / norm if model is not None else None
        cont_n = cont / norm if cont is not None else None
    comps = session.systs.components if session.systs else []

    n = len(names)
    ncols = 1 if n <= 3 else (2 if n <= 8 else 3)
    nrows = -(-n // ncols)
    fig = Figure(figsize=(4.6 * ncols + 0.6, (2.5 if residuals else 1.9) * nrows + 0.5), dpi=100)
    grid = fig.add_gridspec(nrows, ncols, left=0.09 if ncols == 1 else 0.07,
                            right=0.99, top=0.98, bottom=0.12 / max(nrows, 1) + 0.04,
                            hspace=0.12 if residuals else 0.1, wspace=0.18)
    for i, name in enumerate(names):
        cell = grid[i // ncols, i % ncols]
        if residuals:
            sub =cell.subgridspec(2, 1, height_ratios=[3, 1], hspace=0.05)
            ax = fig.add_subplot(sub[0])
            ax_res = fig.add_subplot(sub[1], sharex=ax)
        else:
            ax, ax_res = fig.add_subplot(cell), None
        lam_0 = (ATOM_DATA[name]["wave"] * au.AA).to_value(unit, equivalencies=au.spectral())
        lam_obs = lam_0 * (1 + z)
        v = C_KMS * (x_all - lam_obs) / lam_obs
        m = np.abs(v) <= half_width_kms * 1.05
        if m.any():
            ax.step(v[m], y_n[m], where="mid", lw=0.7, color=colours["flux"])
            ax.step(v[m], dy_n[m], where="mid", lw=0.6, color=ERROR_COLOR)
            if cont_n is not None and not normalized:
                ax.plot(v[m], cont_n[m], lw=0.8, color=CONTINUUM_COLOR)
            if model_n is not None:
                ax.plot(v[m], model_n[m], lw=colours["model_lw"], alpha=colours["model_alpha"],
                        color=colours["model"])
            if cont is not None and normalized:
                ax.axhline(1.0, color="gray", ls=":", lw=0.8)
            ax.axhline(0.0, color="gray", ls=":", lw=0.8)
            top = np.nanpercentile(np.maximum(y_n[m], cont_n[m] if cont_n is not None else y_n[m]), 99)
            top = top if np.isfinite(top) and top > 0 else 1.0
            ax.set_ylim(-0.1, 1.4 if (normalized and cont is not None) else 1.3 * top)
        ax.axvline(0, color="gray", ls=":", lw=0.8)
        trans = ax.get_xaxis_transform()
        for dv, c in tick_velocities(comps, name, z, -half_width_kms, half_width_kms):
            hi = c.series.startswith("Ly") or c.series.startswith("H")
            ax.plot([dv, dv], [0.02, 0.1], transform=trans,
                    color=colours["hi_ticks"] if hi else "gray", lw=1.5 if hi else 1.0, alpha=0.85)
        ax.text(0.98, 0.85, name, transform=ax.transAxes, ha="right", fontweight="bold",
                fontsize=8, bbox=dict(facecolor="white", alpha=0.7, edgecolor="none"))
        ax.set_xlim(-half_width_kms, half_width_kms)
        ax.tick_params(labelsize=7)
        last_row = i // ncols == nrows - 1 or i + ncols >= n
        if ax_res is not None:
            if m.any():
                resid = residuals_sigma(y_all[m], dy_all[m], model[m])
                _draw_residuals(ax_res, v[m], resid, colours)
            ax_res.set_ylabel("σ", fontsize=7)
            ax_res.tick_params(labelsize=7, labelbottom=last_row)
            ax.tick_params(labelbottom=False)
        else:
            ax.tick_params(labelbottom=last_row)
    fig.supxlabel("Relative velocity (Δv, km/s)", fontsize=8)
    fig.supylabel("Normalized flux" if (normalized and cont is not None)
                  else f"Flux ({spec.y.unit})", fontsize=8)
    return fig


def figure_png(fig) -> bytes:
    """A figure as PNG bytes."""
    buf = io.BytesIO()
    fig.savefig(buf, format="png")
    return buf.getvalue()
