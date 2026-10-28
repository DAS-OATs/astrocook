"""The interface of SpectrumV2 across the split of ``spectrum.py``.

Phase 5 moves the methods of SpectrumV2 into thematic mixin classes
(``spectrum_continuum.py``...), which SpectrumV2 inherits. The interface
must stay the same: the same methods and properties as before the split,
and no name defined twice (with multiple inheritance, the first class
would win without a word).
"""
from astrocook.core.spectrum import SpectrumV2

#: Methods and properties of SpectrumV2 before the split (2.2.0)
INTERFACE = {
    '_get_ne_contexts', '_renormalize_model', '_safe', '_xunit', '_xunit_old',
    '_yunit', '_z_em', 'apply_expression', 'calc_intersection_expression',
    'calculate_running_std', 'calibrate', 'coadd', 'cont', 'detect_absorptions',
    'detect_doublet_z', 'dy', 'equalize_to_reference', 'find_absorbed',
    'fit_continuum', 'fit_powerlaw', 'get_column', 'get_region_bounds',
    'get_resampled_column', 'get_velocity_bounds', 'has_aux_column',
    'identify_lines', 'mask_expression', 'mask_region', 'merge_column', 'meta',
    'model', 'norm', 'prepare_expression_context', 'rebin', 'remove_columns',
    'remove_meta_keys', 'resample_on_grid', 'resolution',
    'sanitize_legacy_tellurics', 'scale_flux', 'smooth', 'smooth_column',
    'split', 'stitch', 't', 'update_column', 'update_continuum_from_knots',
    'update_model', 'with_properties', 'x', 'x_convert', 'xmax', 'xmin', 'y',
    'y_convert',
}


def _own_names(cls):
    return {n for n in vars(cls) if not n.startswith('__')}


def test_the_interface_is_unchanged():
    assert {n for n in dir(SpectrumV2) if not n.startswith('__')} == INTERFACE


def test_no_name_is_defined_twice():
    classes = [c for c in SpectrumV2.__mro__ if c is not object]
    assert len(classes) > 1  # SpectrumV2 and its mixins
    seen = {}
    for cls in classes:
        for name in _own_names(cls):
            assert name not in seen, f"{name} in {seen.get(name)} and {cls.__name__}"
            seen[name] = cls.__name__
