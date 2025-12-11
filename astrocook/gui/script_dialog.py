from PySide6.QtCore import Qt
from PySide6.QtGui import QFont, QSyntaxHighlighter, QTextCharFormat, QColor
from PySide6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QTextEdit, 
                               QPushButton, QLabel, QSplitter)
from .scripting import AppController

class ScriptDialog(QDialog):
    def __init__(self, main_window):
        super().__init__(main_window)
        self.resize(800, 600)
        self.setWindowTitle("Astrocook Script Console")
        self.controller = AppController(main_window)
        
        layout = QVBoxLayout(self)
        
        # Editor
        self.editor = QTextEdit()
        self.editor.setFont(QFont("Monospace", 12))
        self.editor.setPlaceholderText("# Access loaded data via 'app'\n"
                                       "# s = app.get_active_session()\n"
                                       "# s.spec.y_convert(...)")
        
        # Output Log
        self.log_view = QTextEdit()
        self.log_view.setReadOnly(True)
        self.log_view.setStyleSheet("background-color: #222; color: #EEE;")
        
        splitter = QSplitter()
        splitter.setOrientation(Qt.Vertical) # Vertical
        splitter.addWidget(self.editor)
        splitter.addWidget(self.log_view)
        
        layout.addWidget(splitter)
        
        # Controls
        btn_layout = QHBoxLayout()
        self.btn_run = QPushButton("Run Script")
        self.btn_run.clicked.connect(self.run_script)
        btn_layout.addStretch()
        btn_layout.addWidget(self.btn_run)
        layout.addLayout(btn_layout)
        
        # Connect Controller Log
        self.controller.log_message.connect(self.append_log)

    def append_log(self, text):
        self.log_view.append(text)

    def run_script(self):
        script = self.editor.toPlainText()
        
        # Define the execution context
        # 'app' is the magic variable exposed to the user
        local_scope = {'app': self.controller}
        
        try:
            self.append_log("--- Executing ---")
            exec(script, {}, local_scope)
            self.append_log("--- Done ---")
        except Exception as e:
            self.append_log(f"Error: {e}")