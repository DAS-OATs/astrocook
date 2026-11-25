"""The evaluation of alternative models in a window of the spectrum (phase 7,
PR 7c.4a; ``design/identification.md``, sections 3.3 and 10).

The identification compares the hypotheses of a group of absorbed regions on
the same pixels: the union of the windows of +-(3 FWHM + 5 b) around every
line of every hypothesis of the group, and of the pixels that a fit of each
hypothesis takes from the data (out to the continuum or to a saddle between
blended lines), so that a broad line started with a small b is compared on
its whole profile, not on its core (Guido, 18 November 2026). :func:`evaluate_model` and the short
fit of :class:`VoigtFitter` compute the model and the profiles on the whole
spectrum, so their cost follows the length of the spectrum, not the pixels of
the comparison (about 1 s per hypothesis on a mock of 85 nm, PR 7c.2).

An :class:`EvaluationWindow` cuts those pixels out of the spectrum once per
group, as a spectrum of their own, and selects once the components that stay
fixed (the rest of the spectrum, outside the group) whose lines reach them:
both are shared by every hypothesis of the group. The evaluation and the
short fit of a hypothesis then run on the window, with its components and
those fixed ones (frozen, not counted), through the same functions as on the
whole spectrum; within a fit the optical depth of the fixed components is
computed once (the background of :meth:`VoigtFitter.prepare_fit_context`).

The short fit is not limited to the window: it takes its pixels from the data
as on the whole spectrum (:meth:`VoigtFitter.feature_mask`, out to the
continuum or to a saddle between blended lines), cut out of the parent per
fit. Limited to the window, which follows the b of the start, a broad or
saturated line started narrow was fitted on its core alone (logN 16, b 60
came out as 15.1, b 80 at R 45000) and a line in a crowded region, with the
fixed list imperfect, moved by more than 0.01 dex or 1 km/s in 10-23% of the
fits (Guido's concern, PR 7c.4a); with the pixels of the data the fit is that
on the whole spectrum.

The statistics agree with those on the whole spectrum to about 1e-5 of the
chi2 (the fine grid of the model follows the narrowest line in the pixels,
which on the whole spectrum can be another one). Folding the fixed
components into the continuum instead, as their transmission convolved apart
and multiplied, was measured and dropped (PR 7c.4a): it was not faster, and
where a fixed line is blended with a line of the group it moved the
difference of the criteria between two hypotheses by 0.2 at R 45000 and 1.3
at R 8900 (exp(-tau_1 - tau_2) convolved is not the product of the two
convolved).
"""
import dataclasses
import functools
import logging
from typing import Any, Dict, List, Optional, Sequence

import astropy.units as au
import numpy as np
from scipy import ndimage

from astrocook.core.atomic_data import ATOM_DATA
from astrocook.core.spectrum import Spectrum
from astrocook.core.structures import ComponentData, DataColumn, SystemListData
from astrocook.core.system_list import SystemList
from astrocook.fitting.lsf import C_KMS, FWHM_TO_SIGMA, KERNEL_SIGMAS
from astrocook.fitting.model_evaluation import MIN_INFLUENCE, ModelEvaluation, evaluate_model
from astrocook.fitting.voigt_fitter import VoigtFitter, _transitions, line_reach

logger = logging.getLogger(__name__)

#: The half-width of the window around a line, in units of the instrumental
#: FWHM and of the b of its component: +-(3 FWHM + 5 b) (decided with Guido,
#: 18 November 2026, ``design/identification.md`` section 10): the wings of
#: the line convolved with the resolution are below the noise beyond it, and
#: the pixels of the continuum on both sides are part of the comparison (a
#: hypothesis that puts absorption where there is none pays for it there)
WINDOW_FWHM = 3.0
WINDOW_B = 5.0

#: Pixels added on both sides of the pixels a fit takes from the data, in the
#: spectrum cut for the fit: the fitter computes its model on its pixels and
#: 50 more on each side (``prepare_fit_context``) and smooths the flux over
#: about 8 when it looks for the limits of a line; with 60 the fit on the cut
#: is that on the whole spectrum
FIT_PAD = 60


def _lines(systs: SystemList) -> List[tuple]:
    """(observed wavelength in Angstrom, b in km/s) of every line of every
    component."""
    out = []
    for comp in systs.components:
        b = float(np.hypot(comp.b, comp.btur or 0.0))
        for t in _transitions(comp.series):
            out.append((ATOM_DATA[t]['wave'] * (1.0 + comp.z), b))
    return out


def window_mask(spec: Spectrum, hypotheses: Sequence[SystemList],
                fwhm_factor: float = WINDOW_FWHM, b_factor: float = WINDOW_B) -> np.ndarray:
    """The pixels within +-(``fwhm_factor`` FWHM + ``b_factor`` b) of a line
    of a component of any of ``hypotheses``: the FWHM the instrumental one
    at the line (c / R of the nearest pixel; 0 where the resolution is not
    known), b that of its component (with btur)."""
    x_ang = spec.x.to_value(au.Angstrom)
    resol = spec.resolution()
    mask = np.zeros(len(x_ang), dtype=bool)
    for systs in hypotheses:
        for lam, b in _lines(systs):
            i = int(np.clip(np.searchsorted(x_ang, lam), 0, len(x_ang) - 1))
            fwhm = C_KMS / resol[i] if np.isfinite(resol[i]) and resol[i] > 0 else 0.0
            half = fwhm_factor * fwhm + b_factor * b
            lo, hi = lam * np.exp(-half / C_KMS), lam * np.exp(half / C_KMS)
            mask[np.searchsorted(x_ang, lo, side='left'):
                 np.searchsorted(x_ang, hi, side='right')] = True
    return mask


def _take(spec: Spectrum, idx: np.ndarray, norm: np.ndarray, norm_unit: Any) -> Spectrum:
    """The pixels ``idx`` of ``spec`` as a spectrum of their own, with the
    flux, its error, the resolution and ``norm`` (the normalization of the
    parent: continuum times telluric model) as its continuum; the other
    columns left out (a fit and an evaluation read only these)."""
    data = spec._data

    def col(c: DataColumn) -> DataColumn:
        return dataclasses.replace(c, values=np.asarray(c.values)[idx])

    aux = {'cont': DataColumn(norm, norm_unit, "normalization of the parent spectrum")}
    if 'resol' in data.aux_cols:
        aux['resol'] = col(data.aux_cols['resol'])
    sub = dataclasses.replace(data, x=col(data.x), xmin=col(data.xmin), xmax=col(data.xmax),
                              y=col(data.y), dy=col(data.dy), aux_cols=aux)
    return Spectrum(sub)


def _static(systs: SystemList) -> List[ComponentData]:
    """The components of ``systs`` at the values of their parameters as the
    model computes them (the links resolved), with no constraints: fixed
    components, not fitted."""
    if not systs.components:
        return []
    cm = systs.constraint_model
    assert cm is not None  # a list with components has its constraint model
    p_full = cm.map_p_free_to_full(cm.p_free_vector)
    return [dataclasses.replace(c, z=float(p_full[4 * i]), logN=float(p_full[4 * i + 1]),
                                b=float(p_full[4 * i + 2]), btur=float(p_full[4 * i + 3]))
            for i, c in enumerate(systs.components)]


def _reaching(sub: Spectrum, comps: List[ComponentData]) -> List[ComponentData]:
    """The components with a line whose reach (its optical depth above the
    cut of the model, :func:`line_reach`, plus 5 sigma of the kernel of the
    lowest resolution) falls on a pixel of ``sub``: the others do not change
    the model there."""
    rows = [(k, ATOM_DATA[t], c) for k, c in enumerate(comps) for t in _transitions(c.series)]
    if not rows:
        return []
    xmin = sub.xmin.to_value(au.Angstrom)
    xmax = sub.xmax.to_value(au.Angstrom)
    r = sub.resolution()
    resols = [float(v) for v in r[np.isfinite(r) & (r > 0)]] + [
        float(c.resol) for c in comps if c.resol is not None and c.resol > 0]
    kernel = KERNEL_SIGMAS * C_KMS / min(resols) * FWHM_TO_SIGMA if resols else 0.0
    owner = np.array([k for k, _, _ in rows])
    lam_c, half = line_reach(
        np.array([a['wave'] for _, a, _ in rows]), np.array([a['f'] for _, a, _ in rows]),
        np.array([a['gamma'] for _, a, _ in rows]), np.array([c.z for _, _, c in rows]),
        np.array([10.0 ** c.logN for _, _, c in rows]),
        np.array([np.hypot(c.b, c.btur or 0.0) for _, _, c in rows]))
    pad = lam_c * kernel / C_KMS
    lo, hi = lam_c - half - pad, lam_c + half + pad
    # A line reaches the window when a pixel has xmax >= lo and xmin <= hi
    i = np.searchsorted(xmax, lo, side='left')
    hits = (i < len(xmin)) & (xmin[np.minimum(i, len(xmin) - 1)] <= hi)
    keep = set(owner[hits].tolist())
    return [c for k, c in enumerate(comps) if k in keep]


@dataclasses.dataclass(frozen=True)
class EvaluationWindow:
    """The pixels of a group cut out of the spectrum, with the fixed
    components that reach them (computed once, shared by the hypotheses of
    the group).

    Build it with :meth:`around` (the window of the lines of the hypotheses)
    or :meth:`from_mask` (any pixels).
    """
    #: The spectrum of the window (its pixels, with the normalization of
    #: the parent as its continuum)
    spectrum: Spectrum
    #: The index of the pixels of the window in the parent spectrum
    index: np.ndarray
    #: The fixed components that reach the window, at their values, without
    #: constraints (the others are left out)
    fixed: List[ComponentData]
    #: The number of pixels of the parent spectrum
    n_parent: int
    #: The parent spectrum and all the fixed components (for the pixels of
    #: the fits, which are taken from the data as on the whole spectrum)
    parent: Spectrum = dataclasses.field(repr=False, compare=False)
    all_fixed: List[ComponentData] = dataclasses.field(repr=False, compare=False)

    @functools.cached_property
    def _parent_fitter(self) -> VoigtFitter:
        """A fitter of the parent spectrum, for :meth:`VoigtFitter.feature_mask`
        (made at the first fit, once per window)."""
        return VoigtFitter(self.parent, SystemList(SystemListData(components=[])))

    @classmethod
    def from_mask(cls, spec: Spectrum, mask: np.ndarray,
                  fixed: Optional[SystemList] = None) -> 'EvaluationWindow':
        """The window of the pixels ``mask`` of ``spec``, with the components
        of ``fixed`` (at their current values) that reach it.

        Raises
        ------
        ValueError
            If the spectrum has no continuum, or no pixel is in ``mask``.
        """
        norm_q = spec.norm
        if norm_q is None:
            raise ValueError("Normalization vector is missing (Flux is physical and "
                             "no continuum found).")
        index = np.flatnonzero(np.asarray(mask, dtype=bool))
        if not index.size:
            raise ValueError("no pixel in the window")
        sub = _take(spec, index, np.asarray(norm_q.value, dtype=float)[index], norm_q.unit)
        static = _static(fixed) if fixed is not None else []
        return cls(spectrum=sub, index=index, fixed=_reaching(sub, static), n_parent=len(spec.x),
                   parent=spec, all_fixed=static)

    @classmethod
    def around(cls, spec: Spectrum, hypotheses: Sequence[SystemList],
               fixed: Optional[SystemList] = None,
               fwhm_factor: float = WINDOW_FWHM, b_factor: float = WINDOW_B,
               z_window_kms: float = 20.0) -> 'EvaluationWindow':
        """The window of the comparison of ``hypotheses``, the same pixels
        for all of them: :func:`window_mask` (+-(3 FWHM + 5 b) around their
        lines) and the pixels that a fit of each of them takes from the data
        (:meth:`VoigtFitter.feature_mask`, out to the continuum or to a
        saddle between blended lines, whatever the b of the start; with
        ``z_window_kms`` as in the fit), with ``fixed``. Without the second,
        a broad line started narrow was compared on its core alone (Guido,
        18 November 2026)."""
        fitter = VoigtFitter(spec, SystemList(SystemListData(components=[])))
        mask = window_mask(spec, hypotheses, fwhm_factor, b_factor)
        for systs in hypotheses:
            if systs.components:
                mask |= fitter.feature_mask(_static(systs), z_window_kms)
        win = cls.from_mask(spec, mask, fixed)
        # The fitter of the parent, for the pixels of the fits: made once
        # (the storage of the cached property)
        win.__dict__['_parent_fitter'] = fitter
        return win

    @property
    def mask(self) -> np.ndarray:
        """The pixels of the window in the parent spectrum."""
        out = np.zeros(self.n_parent, dtype=bool)
        out[self.index] = True
        return out

    def to_parent(self, values: np.ndarray, fill: float = np.nan) -> np.ndarray:
        """An array on the pixels of the window, on those of the parent
        (``fill`` elsewhere)."""
        out = np.full(self.n_parent, fill, dtype=float)
        out[self.index] = values
        return out

    def _with_fixed(self, systs: SystemList,
                    fixed_comps: Optional[List[ComponentData]] = None) -> SystemList:
        """``systs`` with the fixed components (those of the window, or
        ``fixed_comps``) before its own (those with the uuid of one of its
        own left out), its constraints kept."""
        own = {c.uuid for c in systs.components}
        fixed = [c for c in (self.fixed if fixed_comps is None else fixed_comps)
                 if c.uuid not in own]
        return SystemList(dataclasses.replace(
            systs._data, components=fixed + list(systs.components)))

    def evaluate(self, systs: SystemList, uuids: Optional[Sequence[str]] = None,
                 min_influence: float = MIN_INFLUENCE) -> ModelEvaluation:
        """:func:`evaluate_model` on every pixel of the window, the model the
        components of ``systs`` (a hypothesis) and the fixed ones; counted
        the components of ``uuids`` (None: all those of ``systs``, never the
        fixed ones). The residuals on the pixels of the parent spectrum."""
        counted = [c.uuid for c in systs.components] if uuids is None else list(uuids)
        ev = evaluate_model(self.spectrum, self._with_fixed(systs),
                            np.ones(len(self.index), dtype=bool),
                            uuids=counted, min_influence=min_influence)
        return dataclasses.replace(ev, residuals=self.to_parent(ev.residuals))

    def fit_pixels(self, systs: SystemList, z_window_kms: float = 20.0) -> np.ndarray:
        """The pixels of the parent on which :meth:`fit` fits ``systs``: the
        window, and the pixels that a fit of its components takes from the
        data on the whole spectrum (:meth:`VoigtFitter.feature_mask`, out to
        the continuum or to a saddle between blended lines, whatever the b
        of the start), each run widened by ``FIT_PAD`` pixels."""
        feature = self._parent_fitter.feature_mask(_static(systs), z_window_kms)
        if FIT_PAD > 0:
            feature = ndimage.binary_dilation(feature, iterations=FIT_PAD)
        return feature | self.mask

    def fit(self, systs: SystemList, **kwargs: Any) -> SystemList:
        """The fit of the components of ``systs`` (a hypothesis) by
        :meth:`VoigtFitter.fit` with ``kwargs`` (e.g. the ``max_nfev``,
        ``b_max`` and ``logN_max`` of a short fit), the fixed components that
        reach its pixels frozen. Its pixels are those of :meth:`fit_pixels`:
        the fit is that on the whole spectrum, not limited to the window (a
        broad or saturated line started with a small b, or a line in a
        crowded region, takes its pixels from the data, PR 7c.4a). Returns
        ``systs`` with its components fitted (a list without components as it
        is).

        Raises
        ------
        NoPixelsToFit
            If no valid pixel is around the components.
        """
        if not systs.components:
            return systs
        index = np.flatnonzero(self.fit_pixels(systs, kwargs.get('z_window_kms', 20.0)))
        norm_q = self.parent.norm
        assert norm_q is not None  # checked when the window was made
        sub = _take(self.parent, index, np.asarray(norm_q.value, dtype=float)[index],
                    norm_q.unit)
        own = [c.uuid for c in systs.components]
        model = self._with_fixed(systs, _reaching(sub, self.all_fixed))
        with model.fitting_context(own, group_depth=0):
            fitted, _, _ = VoigtFitter(sub, model).fit(**kwargs, model=False)
        keep = set(own)
        return SystemList(dataclasses.replace(
            fitted._data, components=[c for c in fitted.components if c.uuid in keep]))

    def as_dict(self) -> Dict[str, Any]:
        """A summary, for the logs."""
        return {"n_pixels": int(self.index.size), "n_fixed": len(self.fixed)}
