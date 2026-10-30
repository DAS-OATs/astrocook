"""The interface of the GUI classes split into mixins (phase 5, split of the
GUI): the same methods as before the split, and no name defined twice
(with multiple inheritance the first class would win without a word).
Only the classes of Astrocook are looked at, not the Qt bases.
"""
import pytest

from astrocook.gui.main_window import MainWindowV2
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
#: ... and the parts of plot_spectrum, split in phase 5
SPECTRUM_PLOT_WIDGET |= {
    '_add_legend', '_connect_background_capture', '_format_axes',
    '_invalidate_plot_caches', '_legend_location', '_plot_arrays',
    '_plot_aux_column', '_plot_continuum_and_model', '_plot_cursor',
    '_plot_emission_lines', '_plot_error', '_plot_flux', '_plot_other_sessions',
    '_plot_systems', '_raw_plot_arrays', '_read_view',
    '_reconnect_limit_callbacks', '_set_limits', '_slice_to_view', '_was_zoomed',
}
#: ... and the selection of a region, restored in phase 5-bis (F73)
SPECTRUM_PLOT_WIDGET |= {'_region_menu', 'set_region_selection'}
#: ... and the decimation by view, redone in phase 5-bis (F76)
SPECTRUM_PLOT_WIDGET |= {'_decimate'}


#: Methods and class attributes of MainWindowV2 before the split (without
#: the dead code removed in phase 5-bis, F77 and F80)
MAIN_WINDOW = {
    '_RECIPES_REQUIRING_RESOL', '_RECIPES_REQUIRING_Z_EM', '_add_session_internal',
    '_apply_session_info_changes', '_apply_styles',
    '_ask_to_renormalize_model', '_blank_swatch', '_check_and_handle_requirements',
    '_cleanup_progress_ui', '_close_session', '_confirm_discard_continuum_edits',
    '_create_menubar', '_force_restack_floating_widgets', '_get_log_stride',
    '_import_systems_ascii', '_launch_identification_viewer', '_launch_log_scripter',
    '_launch_recipe_dialog', '_log_for_new_session', '_make_card_section',
    '_on_about_requested', '_on_cancel_continuum_clicked',
    '_on_close_multiple_sessions', '_on_close_session_requested', '_on_coadd_requested',
    '_on_cursor_checkbox_toggled', '_on_debug_plot_requested', '_on_duplicate_session',
    '_on_equalized_stitch_dialog_requested', '_on_export_requested',
    '_on_identification_viewer_closed', '_on_import_requested',
    '_on_log_scripter_closed', '_on_open_spectrum', '_on_recipe_dialog_finished',
    '_on_recipe_error', '_on_recipe_finished', '_on_recipe_requested',
    '_on_reset_continuum_clicked', '_on_save_session', '_on_save_session_context',
    '_on_script_error', '_on_script_finished', '_on_script_progress',
    '_on_session_info', '_on_session_list_context_menu', '_on_session_name_changed',
    '_on_session_switched', '_on_set_custom_limits', '_on_stitch_requested',
    '_on_stop_recipe', '_on_stride_change', '_on_stride_slider_pressed',
    '_on_stride_slider_released', '_on_stride_text_changed', '_on_view_identifications',
    '_on_view_log', '_on_view_script_console', '_on_view_system_inspector',
    '_on_view_toggle', '_open_docs', '_populate_snr_combo', '_process_recipe_result',
    '_process_script_result', '_redo_last_action', '_replot_if_session',
    '_reposition_floating_widgets', '_reset_continuum_ui',
    '_safely_close_progress_dialog', '_save_session', '_session_list_height',
    '_set_pending_action', '_set_plot_controls_faded', '_set_sessions_faded',
    '_set_slider_from_stride', '_setup_collapse_buttons', '_setup_empty_view',
    '_setup_left_sidebar', '_setup_plot_view', '_setup_right_sidebar',
    '_show_custom_message', '_swatch',
    '_swatch_checkboxes', '_toggle_continuum_editor', '_toggle_left_sidebar',
    '_toggle_right_sidebar', '_trigger_replot', '_undo_last_action',
    '_update_cursor_and_replot', '_update_limit_boxes_from_plot',
    '_update_sidebar_button_icon', '_update_ui_state', '_update_undo_redo_actions',
    '_update_view_for_session', '_update_window_title', 'add_session', 'dragEnterEvent',
    'dropEvent', 'launch_split_from_current_view', 'open_inspector_on_component',
    'open_session_from_path', 'resizeEvent', 'run_script', 'update_gui_session_state',
    'update_plot_swatches', 'update_session_highlights',
}
#: ... and the split of a selected region (phase 5-bis, F73)
MAIN_WINDOW |= {'launch_split_from_region'}

def _astrocook_classes(cls):
    return [c for c in cls.__mro__ if c.__module__.startswith("astrocook.")]


def _own_names(cls):
    return {n for n in vars(cls) if not n.startswith('__') and n != 'staticMetaObject'}


@pytest.mark.parametrize("cls, interface", [(SpectrumPlotWidget, SPECTRUM_PLOT_WIDGET),
                                            (MainWindowV2, MAIN_WINDOW)])
def test_the_interface_is_unchanged(cls, interface):
    names = set().union(*(_own_names(c) for c in _astrocook_classes(cls)))
    assert names == interface


@pytest.mark.parametrize("cls", [SpectrumPlotWidget, MainWindowV2])
def test_no_name_is_defined_twice(cls):
    seen = {}
    for c in _astrocook_classes(cls):
        for name in _own_names(c):
            assert name not in seen, f"{name} in {seen.get(name)} and {c.__name__}"
            seen[name] = c.__name__
