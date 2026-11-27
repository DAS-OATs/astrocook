"""The statistics of a group fit (phase 7, PR 7d.9).

:class:`~astrocook.core.structures.FitRecord` holds the numbers of one fit
of a group of components: chi2, AIC, BIC, pixels, free parameters,
convergence, evaluations, the recipe that ran it and when. The fitter makes
one at each fit (:func:`fit_record`) and gives it to each component it
fitted; the Details panel of the System Inspector and ``components_list``
of the MCP server show it.

The recipe is read from :func:`running_recipe`, which the boundary of the
recipes (``recipes.params.typed_recipe``) sets for the outermost recipe:
the fitter does not know who called it, and ``fitting/`` does not import
``recipes/``.
"""
from __future__ import annotations

import contextlib
import contextvars
import dataclasses
from datetime import datetime, timezone
from typing import Iterable, Iterator, Optional, Sequence, Tuple

import numpy as np

from astrocook.core.structures import ComponentData, FitRecord

_RECIPE: contextvars.ContextVar[Optional[str]] = contextvars.ContextVar(
    'astrocook_fit_recipe', default=None)


@contextlib.contextmanager
def running_recipe(name: str) -> Iterator[None]:
    """Name the recipe that runs the fits made inside, ``category.name``;
    a recipe run by another keeps the name of the outer one (the one the
    user ran)."""
    token = _RECIPE.set(_RECIPE.get() or name)
    try:
        yield
    finally:
        _RECIPE.reset(token)


def current_recipe() -> Optional[str]:
    """The recipe that is running (:func:`running_recipe`), or None."""
    return _RECIPE.get()


def information_criteria(chi2: float, n: int, k: int
                         ) -> Tuple[float, Optional[float], Optional[float]]:
    """(AIC, AICc, BIC): AIC = chi2 + 2k, AICc = AIC + 2k(k+1)/(n-k-1) (None if n - k - 1 <= 0)
    and BIC = chi2 + k ln n (None without pixels); chi2 stands for
    -2 ln L, whose constant does not matter between models on the same
    pixels (as in ``metals_auto``)."""
    aic = chi2 + 2 * k
    aicc = aic + 2 * k * (k + 1) / (n - k - 1) if n - k - 1 > 0 else None
    bic = chi2 + k * np.log(n) if n > 0 else None
    return (float(aic), None if aicc is None else float(aicc),
            None if bic is None else float(bic))


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec='seconds')


def fit_record(residuals: np.ndarray, n_free: int, chi2_reduced: Optional[float],
               n_constraining: int, members: Iterable[str], *,
               converged: Optional[bool] = None, nfev: Optional[int] = None,
               max_nfev: Optional[int] = None, message: str = "",
               bounds: Sequence[Tuple[str, str, float, float]] = (),
               final: bool = False) -> FitRecord:
    """The record of a fit from its normalized residuals on the pixels of
    the fit (those the optimizer saw, the non-finite already set to 0)."""
    r = np.asarray(residuals, dtype=float)
    r = r[np.isfinite(r)]
    chi2 = float(np.sum(r ** 2))
    n_pix = int(r.size)
    aic, _, bic = information_criteria(chi2, n_pix, int(n_free))
    return FitRecord(
        chi2=chi2, n_pix=n_pix, n_free=int(n_free), aic=aic, bic=bic,
        chi2_reduced=None if chi2_reduced is None else float(chi2_reduced),
        n_constraining=int(n_constraining), members=tuple(members),
        converged=converged, nfev=None if nfev is None else int(nfev),
        max_nfev=None if max_nfev is None else int(max_nfev), message=str(message),
        recipe=current_recipe(), time=_now(), final=final,
        bounds=tuple((str(u), str(p), float(lo), float(hi)) for u, p, lo, hi in bounds))


def previous_fit(components: Sequence[ComponentData]) -> Optional[FitRecord]:
    """The record most of ``components`` hold (that of their last fit), for
    the statistics recomputed at the final values; None if none holds one."""
    counts: dict = {}
    for c in components:
        if c.fit is not None:
            counts.setdefault(c.fit.id, [c.fit, 0])[1] += 1
    if not counts:
        return None
    return max(counts.values(), key=lambda v: v[1])[0]


def final_record(previous: Optional[FitRecord], residuals: np.ndarray, n_free: int,
                 chi2_reduced: Optional[float], n_constraining: int,
                 members: Iterable[str]) -> FitRecord:
    """The record of the statistics computed again at the final values
    (``evaluate_uncertainties``): the numbers of the model now, the
    convergence and the evaluations of the fit before (``previous``)."""
    new = fit_record(residuals, n_free, chi2_reduced, n_constraining, members, final=True)
    if previous is None:
        return new
    return dataclasses.replace(new, converged=previous.converged, nfev=previous.nfev,
                               max_nfev=previous.max_nfev, message=previous.message,
                               recipe=previous.recipe or new.recipe, bounds=previous.bounds)


def records_of(components: Iterable[ComponentData]) -> dict:
    """The records held by ``components``, by id, in their order."""
    out: dict = {}
    for c in components:
        if c.fit is not None:
            out.setdefault(c.fit.id, c.fit)
    return out
