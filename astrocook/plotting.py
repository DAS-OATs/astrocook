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

from astrocook.core.system_list import get_all_w_rest

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
                  normalized: bool = False, palette: str = "standard"):
    """A figure of the spectrum of a session in a range of x.

    The flux and its error, the continuum and the model (if any), and a
    tick at each line of each component, labelled with its series (the HI
    series in the colour of the HI ticks of the GUI).

    Parameters
    ----------
    session : SessionV2
    x_min, x_max : float, optional
        The range, in the unit of x; None: the start (end) of the spectrum.
    normalized : bool, optional
        Divide by the continuum, if there is one.
    palette : str, optional
        ``'standard'`` or ``'colorblind'`` (:data:`PLOT_PALETTES`).

    Returns
    -------
    matplotlib.figure.Figure

    Raises
    ------
    ValueError
        If the range is empty or outside the spectrum.
    """
    from matplotlib.figure import Figure

    colours = PLOT_PALETTES.get(palette, PLOT_PALETTES["standard"])
    spec = session.spec
    data = spec.columns_in_range(["y", "dy"], x_min, x_max)
    x = data["columns"]["x"]
    lo, hi = x.min(), x.max()
    sel = (spec.x.value >= lo) & (spec.x.value <= hi)
    cont = spec.cont.value[sel] if spec.cont is not None else None
    norm = cont if (normalized and cont is not None) else 1.0
    fig = Figure(figsize=(10, 3.6), dpi=100)
    ax = fig.add_subplot()
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
    ax.set_xlabel(f"x ({spec.x.unit})")
    ax.set_ylabel("normalized flux" if np.ndim(norm) else f"flux ({spec.y.unit})")
    ax.legend(fontsize=7, loc="lower left")
    fig.tight_layout()
    return fig


def figure_png(fig) -> bytes:
    """A figure as PNG bytes."""
    buf = io.BytesIO()
    fig.savefig(buf, format="png")
    return buf.getvalue()
