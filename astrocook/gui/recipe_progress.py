"""The progress of a running recipe: the dialog with its bar, its last
status and its Stop button, the log handler that feeds it from any thread,
and the handler that collects the warnings to show at the end.

Moved out of ``main_window.py`` in phase 5 (split of the GUI).
"""
import logging

from PySide6.QtCore import QObject, Qt, Signal
from PySide6.QtWidgets import (QApplication, QDialog, QHBoxLayout, QLabel, QProgressBar,
                               QPushButton, QVBoxLayout)


class LogSignalBridge(QObject):
    """Thread-safe bridge to emit log messages as Qt signals."""
    new_log_message = Signal(str)

class QtLogHandler(logging.Handler):
    """Python Logging Handler that pipes messages to a Qt Signal."""
    def __init__(self, bridge):
        super().__init__()
        self.bridge = bridge
    
    def emit(self, record):
        # Format the message
        msg = self.format(record)
        # Emit signal (thread-safe way to talk to GUI)
        self.bridge.new_log_message.emit(msg)

class RecipeProgressDialog(QDialog):
    """
    A minimal progress dialog for long-running recipes.
    Shows a determinate bar and the last status message.
    """
    stop_requested = Signal()

    def __init__(self, title, parent=None):
        super().__init__(parent)
        self.setWindowTitle(f"Running {title}")
        self.setModal(True)
        # Resize to be compact but wide enough for text
        self.resize(400, 120)
        
        # Remove [X] to prevent unsafe closing
        self.setWindowFlags(self.windowFlags() & ~Qt.WindowCloseButtonHint)

        layout = QVBoxLayout(self)
        layout.setSpacing(15)
        layout.setContentsMargins(20, 20, 20, 20)
        
        # Header
        self.header_label = QLabel("Recipe Running...")
        self.header_label.setStyleSheet("font-size: 14px;")
        layout.addWidget(self.header_label)

        # Progress Bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setValue(0)
        self.progress_bar.setTextVisible(True)
        # Optional: Make it look a bit thicker/modern
        self.progress_bar.setStyleSheet("QProgressBar { min-height: 20px; }")
        layout.addWidget(self.progress_bar)

        # Status Label (Last log line)
        self.status_label = QLabel("Initializing...")
        self.status_label.setStyleSheet("font-size: 11px;")
        self.status_label.setWordWrap(True)
        layout.addWidget(self.status_label)

        # Add a stretch to keep things tight at the top
        layout.addStretch()

        # Use a horizontal layout to manage alignment
        btn_layout = QHBoxLayout()
        btn_layout.setContentsMargins(0, 0, 0, 0)

        # 1. Add a "spring" to push the button to the right
        btn_layout.addStretch()

        # [NEW] Stop Button
        self.stop_btn = QPushButton("Stop")
        self.stop_btn.setMaximumWidth(100)
        self.stop_btn.setToolTip("Stop execution gracefully")
        self.stop_btn.clicked.connect(self._on_stop_clicked)
        
        # 3. Add button to the horizontal row
        btn_layout.addWidget(self.stop_btn)

        # 4. Add the row to the main vertical layout
        layout.addLayout(btn_layout)

    def update_log(self, message):
        """
        Parses logs. 
        - '>> PROGRESS: 50' -> Updates bar.
        - Other text -> Updates status label.
        """
        if ">> PROGRESS:" in message:
            try:
                parts = message.split(">> PROGRESS:")
                val = int(parts[1].strip())
                self.progress_bar.setValue(val)
            except ValueError:
                pass
            return

        clean_msg = message.strip()
        if clean_msg:
            # Update the text label
            self.status_label.setText(clean_msg)
            # Force UI update immediately so text doesn't lag behind heavy calculations
            QApplication.processEvents()
    
    def _on_stop_clicked(self):
        self.stop_btn.setEnabled(False)
        self.stop_btn.setText("Stopping...")
        self.status_label.setText("Attempting to stop gracefully...")
        # Emit signal to parent
        self.stop_requested.emit()
        
class WarningCaptureHandler(logging.Handler):
    """Captures warning messages to a list for GUI display."""
    def __init__(self):
        super().__init__()
        self.warnings = []
        
    def emit(self, record):
        # Only capture meaningful warnings (skip debug/info)
        if record.levelno >= logging.WARNING:
            self.warnings.append(self.format(record))
