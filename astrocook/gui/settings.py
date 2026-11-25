"""The user's settings of the GUI (phase 7-quater, PR U2).

One place for ``QSettings("Astrocook", "AstrocookV2")`` (the application name
is the key of the settings already saved by the users: it stays), the keys
and defaults of the *Settings...* dialog, and small readers that never
raise on a value of the wrong type.

Rule of the PR: a setting is about the interface only. None changes a number,
none applies to scripts, the MCP server or the log (a recipe default stays in
the registry; the atomic-data lists and the thresholds of identification or
fit are never settings).
"""
import os
from datetime import datetime, timedelta
from typing import Any, Optional

from PySide6.QtCore import QSettings

ORGANIZATION = "Astrocook"
APPLICATION = "AstrocookV2"  # the key of the saved settings (not the version)

#: The format of ``app_settings``: the native one of the platform (INI file
#: on Linux, registry on Windows, plist on macOS), as ``QSettings(org, app)``
#: always used. The tests put it to INI in a folder of their own (conftest).
SETTINGS_FORMAT = QSettings.Format.NativeFormat


def app_settings() -> QSettings:
    """The settings of the application."""
    return QSettings(SETTINGS_FORMAT, QSettings.Scope.UserScope, ORGANIZATION, APPLICATION)


# --- Keys (the first five were already in use: unchanged) ---
APPEARANCE_KEY = "appearance"
PLOT_COLOURS_KEY = "plot_colours"
RECENT_FILES_KEY = "recent_files"
LAST_DIR_KEY = "last_working_dir"
# General
CHECK_UPDATES_KEY = "check_updates_at_start"
CHECK_FREQUENCY_KEY = "check_updates_frequency"
LAST_CHECK_KEY = "last_update_check"
MAX_RECENT_KEY = "max_recent_files"
CONFIRM_CLOSE_KEY = "confirm_close_unsaved"
DIR_MODE_KEY = "working_dir_mode"
FIXED_DIR_KEY = "working_dir_fixed"
# Appearance
FONT_SIZE_KEY = "font_size"
LINE_WIDTH_KEY = "plot_line_width"
INITIAL_X_UNIT_KEY = "initial_x_unit"
DECIMALS_KEY = "table_decimals"
# Performance
DOWNSAMPLE_KEY = "downsample_threshold"

#: The values that were hard-coded before the dialog: the defaults keep the
#: behaviour of the previous version (and its numbers) exactly
DEFAULTS: dict[str, Any] = {
    CHECK_UPDATES_KEY: False,         # no network access unless asked for
    CHECK_FREQUENCY_KEY: "daily",
    MAX_RECENT_KEY: 10,
    CONFIRM_CLOSE_KEY: True,
    DIR_MODE_KEY: "last",
    FIXED_DIR_KEY: "",
    APPEARANCE_KEY: "system",
    PLOT_COLOURS_KEY: "standard",
    FONT_SIZE_KEY: 0,                 # 0: the font of the system
    LINE_WIDTH_KEY: 1.0,              # a factor on the widths of the plot
    INITIAL_X_UNIT_KEY: "nm",
    DECIMALS_KEY: 5,                  # of the wavelengths of the Data Inspector
    DOWNSAMPLE_KEY: 1000000,          # pixels in view (pyside_plot.DECIMATION_THRESHOLD)
}

CHECK_FREQUENCIES = {"start": "Every start", "daily": "Once a day", "weekly": "Once a week"}
FREQUENCY_DELTA = {"start": timedelta(0), "daily": timedelta(days=1),
                   "weekly": timedelta(days=7)}
X_UNITS = ("nm", "Angstrom", "micron")
DIR_MODES = {"last": "The last one used", "fixed": "A fixed folder"}

#: Limits of the numeric settings (also the ranges of the dialog)
RECENT_RANGE = (0, 30)
FONT_RANGE = (8, 24)
LINE_WIDTH_RANGE = (0.5, 3.0)
DECIMALS_RANGE = (2, 8)
DOWNSAMPLE_RANGE = (10_000, 20_000_000)


def _as_bool(value: Any, default: bool) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, str):  # INI files give back "true" / "false"
        return value.lower() == "true" if value.lower() in ("true", "false") else default
    return default


def get_bool(key: str) -> bool:
    return _as_bool(app_settings().value(key, DEFAULTS[key]), DEFAULTS[key])


def get_choice(key: str, allowed) -> str:
    value: Any = app_settings().value(key, DEFAULTS[key])
    return value if value in allowed else DEFAULTS[key]


def get_number(key: str, limits: tuple, kind=int):
    """The number of ``key`` inside ``limits`` (the default when it is not one)."""
    value: Any = app_settings().value(key, DEFAULTS[key])
    try:
        number = kind(float(value))
    except (TypeError, ValueError):
        return DEFAULTS[key]
    return max(limits[0], min(limits[1], number))


def max_recent_files() -> int:
    return get_number(MAX_RECENT_KEY, RECENT_RANGE)


def confirm_close_unsaved() -> bool:
    return get_bool(CONFIRM_CLOSE_KEY)


def initial_x_unit() -> str:
    return get_choice(INITIAL_X_UNIT_KEY, X_UNITS)


def table_decimals() -> int:
    return get_number(DECIMALS_KEY, DECIMALS_RANGE)


def line_width_factor() -> float:
    return get_number(LINE_WIDTH_KEY, LINE_WIDTH_RANGE, float)


def downsample_threshold(default: Optional[int] = None) -> int:
    """The pixels in view above which the plot is downsampled. ``default``
    (the constant of ``pyside_plot``) is what is used while the user has not
    chosen one."""
    if default is not None and app_settings().value(DOWNSAMPLE_KEY) is None:
        return default
    return get_number(DOWNSAMPLE_KEY, DOWNSAMPLE_RANGE)


def starting_folder() -> str:
    """The folder where the file dialogs start: the last one used, or the
    fixed one of the settings (when it still exists)."""
    settings = app_settings()
    if get_choice(DIR_MODE_KEY, DIR_MODES) == "fixed":
        fixed = str(settings.value(FIXED_DIR_KEY, "") or "")
        if fixed and os.path.isdir(fixed):
            return fixed
    return str(settings.value(LAST_DIR_KEY, os.getcwd()))


def update_check_due(now: Optional[datetime] = None) -> bool:
    """Whether the automatic check for updates is on and its frequency says
    it is time (``LAST_CHECK_KEY`` holds the ISO time of the last one)."""
    if not get_bool(CHECK_UPDATES_KEY):
        return False
    frequency = get_choice(CHECK_FREQUENCY_KEY, CHECK_FREQUENCIES)
    last = app_settings().value(LAST_CHECK_KEY, "")
    try:
        last_time = datetime.fromisoformat(str(last))
    except ValueError:
        return True
    return (now or datetime.now()) - last_time >= FREQUENCY_DELTA[frequency]


def record_update_check(now: Optional[datetime] = None) -> None:
    app_settings().setValue(LAST_CHECK_KEY, (now or datetime.now()).isoformat())
