"""A model brought onto the grid of a spectrum (phase 7-ter, PR T2).

An external model (a telluric transmission, and in T3-T4 a continuum or an
absorption model) is given on a grid of its own: the pixels of the model,
with their edges. :func:`model_on_grid` gives its values on the pixels of
the spectrum:

- a *copy* where the model is on the grid of the spectrum: over the
  pixels where the two overlap, every centre of the model within
  :data:`SAME_GRID_TOLERANCE` of a pixel of the spectrum (the axes written
  in a text file lose a few decimals); the model may cover a part of the
  spectrum or more than it;
- otherwise an *average over the pixels*: for each pixel of the spectrum
  (``xmin``-``xmax``) the mean of the model over that interval, each pixel
  of the model weighted by its overlap. Where the model does not cover a
  pixel of the spectrum (beyond its ends, in a gap, over pixels of the
  model whose value is not finite) it takes the value of the role beyond
  the model (``fill``): with a finite one (1 for a transmission) a pixel
  covered in part is the mean of the model and of ``fill`` weighted by
  their widths; with NaN (a continuum) only the pixels covered in full
  (:data:`COVERED_FRACTION`) have a value.

The model is a step function (one value per pixel of the model), so that
the average conserves its integral: the equivalent width of a telluric
line is the same on the two grids. Pure functions on arrays in one unit;
the units are converted by the caller.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

#: Two grids are the same when every centre of the model is within this
#: fraction of the width of the pixel of the spectrum (T2: the axis of a
#: model written in a text file with fewer decimals than the spectrum, e.g.
#: four decimals in nm at 500 nm on a grid of 1 km/s, ~0.0017 nm, is off by
#: up to 3% of a pixel; a model shifted by less than this is copied).
#: Confirmed by Guido (20 November 2026)
SAME_GRID_TOLERANCE = 0.1
#: A pixel of the spectrum is covered by the model when the model covers
#: this fraction of its width (the edges of a model written in a text file
#: with six decimals in nm are off by up to 3e-4 of a pixel of 1 km/s at
#: 500 nm; those converted between units differ in the last digits).
#: Confirmed by Guido (20 November 2026)
COVERED_FRACTION = 1.0 - 1e-3


@dataclass(frozen=True)
class ModelOnGrid:
    """The values of a model on the pixels of a spectrum."""
    #: One value per pixel of the spectrum (``fill`` where not covered)
    values: np.ndarray
    #: True for the pixels of the spectrum the model covers
    covered: np.ndarray
    #: ``'copy'`` or ``'average'``
    method: str


def auto_edges(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """The edges of contiguous pixels from their centres: each pixel
    extends half the distance to the previous one on both sides (the first
    as the second), as the loaders do for a file without edges."""
    x = np.asarray(x, dtype=float)
    if len(x) > 1:
        diff = np.diff(x, prepend=x[0] - (x[1] - x[0]))
    else:
        diff = np.array([1.0])
    return x - diff / 2.0, x + diff / 2.0


def same_grid(x_model: np.ndarray, x: np.ndarray, xmin: np.ndarray,
              xmax: np.ndarray) -> bool:
    """Whether the centres of the model are those of the spectrum, within
    :data:`SAME_GRID_TOLERANCE` of the width of each pixel."""
    x_model = np.asarray(x_model, dtype=float)
    if x_model.shape != np.shape(x):
        return False
    width = np.abs(np.asarray(xmax, dtype=float) - np.asarray(xmin, dtype=float))
    with np.errstate(invalid='ignore'):
        return bool(np.all(np.abs(x_model - np.asarray(x, dtype=float))
                           <= SAME_GRID_TOLERANCE * width))


def _cumulative(m_xmin, m_xmax, values):
    """The knots of the integral of the model and of its coverage, as
    functions of the axis (piecewise linear between the knots)."""
    order = np.argsort(0.5 * (m_xmin + m_xmax), kind='stable')
    lo, hi, v = m_xmin[order], m_xmax[order], values[order]
    lo, hi = np.minimum(lo, hi), np.maximum(lo, hi)
    # Overlapping pixels (edges computed from rounded centres) meet halfway
    overlap = hi[:-1] > lo[1:]
    if np.any(overlap):
        mid = 0.5 * (hi[:-1] + lo[1:])
        hi[:-1] = np.where(overlap, mid, hi[:-1])
        lo[1:] = np.where(overlap, mid, lo[1:])
        hi = np.maximum(hi, lo)
    good = np.isfinite(v) & np.isfinite(lo) & np.isfinite(hi)
    width = np.where(good, hi - lo, 0.0)
    area = np.where(good, v, 0.0) * width
    knots = np.empty(2 * len(lo))
    knots[0::2], knots[1::2] = lo, hi
    cum_area = np.concatenate([[0.0], np.cumsum(area)])
    cum_width = np.concatenate([[0.0], np.cumsum(width)])
    area_knots = np.empty_like(knots)
    width_knots = np.empty_like(knots)
    area_knots[0::2], area_knots[1::2] = cum_area[:-1], cum_area[1:]
    width_knots[0::2], width_knots[1::2] = cum_width[:-1], cum_width[1:]
    return knots, area_knots, width_knots


def average_over_pixels(m_xmin: np.ndarray, m_xmax: np.ndarray, values: np.ndarray,
                        xmin: np.ndarray, xmax: np.ndarray, fill: float) -> ModelOnGrid:
    """The mean of the model (pixels ``m_xmin``-``m_xmax``, ``values``) over
    each pixel of the spectrum (``xmin``-``xmax``), with ``fill`` where the
    model does not cover it (module docstring). ``covered``: the pixels
    covered in full."""
    m_xmin, m_xmax = np.asarray(m_xmin, dtype=float), np.asarray(m_xmax, dtype=float)
    values = np.asarray(values, dtype=float)
    xmin, xmax = np.asarray(xmin, dtype=float), np.asarray(xmax, dtype=float)
    out = np.full(len(xmin), fill, dtype=float)
    if len(values) == 0:
        return ModelOnGrid(out, np.zeros(len(xmin), dtype=bool), 'average')
    knots, area_knots, width_knots = _cumulative(m_xmin, m_xmax, values)
    a, b = np.minimum(xmin, xmax), np.maximum(xmin, xmax)
    area = np.interp(b, knots, area_knots) - np.interp(a, knots, area_knots)
    covered_width = np.interp(b, knots, width_knots) - np.interp(a, knots, width_knots)
    with np.errstate(invalid='ignore'):
        covered = (b > a) & (covered_width >= COVERED_FRACTION * (b - a))
    if np.isfinite(fill):
        width = b - a
        some = (width > 0) & (covered_width > 0)
        out[some] = (area[some] + (width[some] - covered_width[some]) * fill) / width[some]
        out[covered] = area[covered] / covered_width[covered]
    else:
        out[covered] = area[covered] / covered_width[covered]
    return ModelOnGrid(out, covered, 'average')


def _aligned_run(m_x: np.ndarray, x: np.ndarray, xmin: np.ndarray, xmax: np.ndarray):
    """The pixels of the spectrum and of the model that coincide, when the
    model is on the grid of the spectrum over their overlap (the same
    grid, a part of it, or a wider one): the slices of the two, or None."""
    if len(m_x) == 0 or len(x) == 0:
        return None
    if not (np.all(np.diff(m_x) > 0) and np.all(np.diff(x) > 0)):
        return None
    # The pixel of the other grid nearest to the first of each
    if m_x[0] >= x[0]:
        k = int(np.argmin(np.abs(x - m_x[0])))
        j = 0
    else:
        k = 0
        j = int(np.argmin(np.abs(m_x - x[0])))
    n = min(len(x) - k, len(m_x) - j)
    if n <= 0:
        return None
    s_spec, s_model = slice(k, k + n), slice(j, j + n)
    if not same_grid(m_x[s_model], x[s_spec], xmin[s_spec], xmax[s_spec]):
        return None
    return s_spec, s_model


def model_on_grid(m_x: np.ndarray, m_xmin: np.ndarray, m_xmax: np.ndarray,
                  values: np.ndarray, x: np.ndarray, xmin: np.ndarray,
                  xmax: np.ndarray, fill: float) -> ModelOnGrid:
    """The model on the pixels of the spectrum: a copy where the model is on
    the grid of the spectrum (the same grid, or a part of it, or a wider
    one: :func:`same_grid` over their overlap; ``fill`` beyond), otherwise
    the average over the pixels (:func:`average_over_pixels`). In a copy, a
    value of the model that is not finite gets ``fill``, as in the
    average."""
    m_x = np.asarray(m_x, dtype=float)
    x = np.asarray(x, dtype=float)
    xmin, xmax = np.asarray(xmin, dtype=float), np.asarray(xmax, dtype=float)
    values = np.asarray(values, dtype=float)
    run = _aligned_run(m_x, x, xmin, xmax)
    if run is not None:
        s_spec, s_model = run
        out = np.full(len(x), fill, dtype=float)
        covered = np.zeros(len(x), dtype=bool)
        part = values[s_model]
        covered[s_spec] = np.isfinite(part)
        out[s_spec] = np.where(covered[s_spec], part, fill)
        return ModelOnGrid(out, covered, 'copy')
    return average_over_pixels(m_xmin, m_xmax, values, xmin, xmax, fill)
