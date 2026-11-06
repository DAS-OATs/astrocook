"""The dialog of *File > Export Log as Text* (phase 6-ter, PR 5c).

It shows the log of the active session as the text of a methods section
(:func:`astrocook.methods_text.methods_text`), to copy or to save as a text
file. The text is made by the package, as for the scripts and the MCP
server; the dialog only shows it.
"""
import logging
from pathlib import Path
from typing import Optional

from PySide6.QtWidgets import (
    QApplication, QCheckBox, QDialog, QDialogButtonBox, QFileDialog, QHBoxLayout,
    QLabel, QMessageBox, QPlainTextEdit, QPushButton, QVBoxLayout)

from astrocook.methods_text import methods_file_name, methods_text

logger = logging.getLogger(__name__)


class MethodsTextDialog(QDialog):
    """The text of the methods of a session history, read only."""

    def __init__(self, history, parent=None, start_dir: Optional[str] = None):
        super().__init__(parent)
        self.history = history
        self.start_dir = start_dir or ""
        name = getattr(history, 'display_name', None) or "session"
        self.setWindowTitle(f"Log as Text: {name}")
        self.resize(760, 560)

        layout = QVBoxLayout(self)
        layout.addWidget(QLabel(
            "The session log as the text of a methods section: one paragraph per "
            "recipe, with its reason, what it did and changed, and the notes. "
            "The recipes undone are left out."))
        self.menus_box = QCheckBox("Name the menus of the recipes")
        self.menus_box.setToolTip("For who reruns the analysis in the GUI; uncheck for a paper.")
        self.menus_box.setChecked(True)
        self.menus_box.toggled.connect(self.refresh)
        layout.addWidget(self.menus_box)

        self.text_edit = QPlainTextEdit(self)
        self.text_edit.setReadOnly(True)
        self.text_edit.setLineWrapMode(QPlainTextEdit.LineWrapMode.WidgetWidth)
        layout.addWidget(self.text_edit)

        buttons = QHBoxLayout()
        self.copy_button = QPushButton("Copy")
        self.copy_button.clicked.connect(self.copy)
        buttons.addWidget(self.copy_button)
        self.save_button = QPushButton("Save...")
        self.save_button.clicked.connect(self.save)
        buttons.addWidget(self.save_button)
        buttons.addStretch()
        close_box = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        close_box.rejected.connect(self.reject)
        buttons.addWidget(close_box)
        layout.addLayout(buttons)
        self.refresh()

    @property
    def text(self) -> str:
        return self.text_edit.toPlainText()

    def refresh(self) -> None:
        """Write the text again (e.g. after the menus are switched)."""
        self.text_edit.setPlainText(
            methods_text(self.history, menus=self.menus_box.isChecked()).text)

    def copy(self) -> None:
        QApplication.clipboard().setText(self.text)

    def save(self, path: Optional[str] = None) -> Optional[str]:
        """Save the text to ``path``, or to the file chosen in a dialog; the
        path written, or None."""
        if path is None:
            name = getattr(self.history, 'display_name', None) or "session"
            default = str(Path(self.start_dir) / methods_file_name(name))
            path, _ = QFileDialog.getSaveFileName(
                self, "Save Log as Text", default, "Text Files (*.txt);;All Files (*)")
            if not path:
                return None
        try:
            Path(path).write_text(self.text + "\n", encoding="utf-8")
        except OSError as e:
            logger.error(f"The text could not be saved to {path}: {e}")
            QMessageBox.critical(self, "Save Log as Text", f"The text could not be saved:\n{e}")
            return None
        logger.info(f"Log as text saved to {path}.")
        return path
