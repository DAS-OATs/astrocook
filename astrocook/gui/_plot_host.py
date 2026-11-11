"""What ``PlotMouseMixin`` uses of the widget that mixes it in, declared
for the type checker (phase 6-quater, PR 4-bis; see ``_window_host.py``)."""
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from PySide6.QtWidgets import QWidget

    from astrocook.gui.main_window import MainWindowV2

    class PlotHost(QWidget):
        """The part of ``SpectrumPlotWidget`` that the mouse mixin uses."""
        main_window: MainWindowV2
        canvas: Any
        toolbar: Any
        _is_selecting_region: bool

    PlotHostBase = PlotHost
else:
    PlotHostBase = object
