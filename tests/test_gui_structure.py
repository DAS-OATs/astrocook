"""The interface of the GUI classes split into mixins (phase 5, split of the
GUI): the same methods as before the split, and no name defined twice
(with multiple inheritance the first class would win without a word).
Only the classes of Astrocook are looked at, not the Qt bases.
"""
import pytest

from astrocook.gui.pyside_plot import SpectrumPlotWidget

#: Methods of SpectrumPlotWidget before the split (staticMetaObject is
#: added by PySide to each class)
SPECTRUM_PLOT_WIDGET = {
    '_apply_local_patch', '_blit_highlights', '_draw_isomag_grid',
    '_generate_velocity_knots', '_get_sanitized_knots',
    '_get_sanitized_knots_with_anchors', '_highlight_guide', '_on_full_draw',
    '_open_data_inspector', '_remove_editor_artists', '_update_tooltip_display',
    'begin_stride_interaction', 'calculate_z_from_x', 'cancel_continuum_edit',
    'end_stride_interaction', 'get_cursor_line_positions_at_z',
    'on_background_draw', 'on_mouse_move', 'on_press', 'on_release',
    'plot_spectrum', 'refresh_edit_layer', 'reset_continuum_to_original',
    'set_highlights', 'show_data_tooltip', 'start_continuum_edit',
    'stop_continuum_edit', 'toggle_isomag_grid', 'toggle_region_selector',
    'update_continuum_stride', 'update_plot',
}


def _astrocook_classes(cls):
    return [c for c in cls.__mro__ if c.__module__.startswith("astrocook.")]


def _own_names(cls):
    return {n for n in vars(cls) if not n.startswith('__') and n != 'staticMetaObject'}


@pytest.mark.parametrize("cls, interface", [(SpectrumPlotWidget, SPECTRUM_PLOT_WIDGET)])
def test_the_interface_is_unchanged(cls, interface):
    names = set().union(*(_own_names(c) for c in _astrocook_classes(cls)))
    assert names == interface


@pytest.mark.parametrize("cls", [SpectrumPlotWidget])
def test_no_name_is_defined_twice(cls):
    seen = {}
    for c in _astrocook_classes(cls):
        for name in _own_names(c):
            assert name not in seen, f"{name} in {seen.get(name)} and {c.__name__}"
            seen[name] = c.__name__
