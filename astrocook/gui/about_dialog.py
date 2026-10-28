"""The About dialog: version, credits, and the links to report an issue and
to cite Astrocook.

Moved out of ``main_window.py`` in phase 5 (split of the GUI).
"""
import os

from PySide6.QtCore import QUrl, Qt
from PySide6.QtGui import QDesktopServices, QPixmap
from PySide6.QtWidgets import QDialog, QDialogButtonBox, QLabel, QPushButton, QVBoxLayout

from astrocook.core.utils import resource_path


class AboutDialog(QDialog):
    """
    Dialog to display application information, credits, and citation.
    """
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("About Astrocook")
        self.setFixedWidth(400)
        
        layout = QVBoxLayout(self)
        layout.setSpacing(5)
        layout.setContentsMargins(20, 20, 20, 20)

        layout.addSpacing(20)

        # --- Logo ---
        try:
            # Reusing the resource path logic from MainWindow
            icon_path = resource_path(os.path.join("assets", "logo_3d_HR.png"))
            if os.path.exists(icon_path):
                logo_lbl = QLabel()
                pixmap = QPixmap(icon_path).scaled(160, 160, Qt.KeepAspectRatio, Qt.SmoothTransformation)
                logo_lbl.setPixmap(pixmap)
                logo_lbl.setAlignment(Qt.AlignCenter)
                layout.addWidget(logo_lbl)
        except Exception: 
            pass

        # --- Title & Version ---
        title_lbl = QLabel("Astrocook V2")
        title_lbl.setStyleSheet("font-size: 20pt; font-weight: bold; color: #296bff;")
        title_lbl.setAlignment(Qt.AlignCenter)
        #layout.addWidget(title_lbl)

        copy_lbl = QLabel("© 2017-2026 Guido Cupani")
        copy_lbl.setStyleSheet("font-size: 14pt; color: #296bff;")
        copy_lbl.setAlignment(Qt.AlignCenter)
        layout.addWidget(copy_lbl)

        # Try to import version, else fallback
        try:
            from astrocook import __version__ as ver
        except ImportError:
            ver = "unknown"

        ver_lbl = QLabel(f"Version: {ver}")
        ver_lbl.setStyleSheet("font-weight: bold; color: #f5a100;")
        ver_lbl.setAlignment(Qt.AlignCenter)
        layout.addWidget(ver_lbl)

        layout.addSpacing(20)

        tagline_lbl = QLabel("A thousand recipes to cook a spectrum")
        tagline_lbl.setStyleSheet("font-size: 16pt; font-weight: bold;")
        tagline_lbl.setAlignment(Qt.AlignCenter)
        layout.addWidget(tagline_lbl)

        layout.addSpacing(20)

        # --- Description / Credits ---
        credits_html = (
            "<p align='center'>Developed by <b>Guido Cupani</b><br>at INAF–Osservatorio Astronomico di Trieste.</p>"
            "<p align='center'><b>Many thanks</b> to (in alphabetical order):<br>Giorgio Calderone, Stefano Cristiani, Simona Di Stefano,<br>Valentina D'Odorico, Francesco Guarneri, Elena Marcuzzo,<br> Christian Piscitelli, Stefano Alberto Russo, Andrea Trost.</p>"
            "<p align='center'>Logo and icon by <a href='https://claudiacupani.com'><b>Claudia Cupani</b></a>.</p>"
            "<p align='center'>Partly funded by the PRIN PNRR project<br>“Next generation computing & data technologies<br>to probe the cosmic metal content”.</p>"
        )
        credits_lbl = QLabel(credits_html)
        credits_lbl.setWordWrap(True)
        credits_lbl.setAlignment(Qt.AlignCenter)
        credits_lbl.setOpenExternalLinks(True)
        layout.addWidget(credits_lbl)

        layout.addSpacing(30)

        # --- Buttons ---
        # Create the standard Dialog Button Box
        btn_box = QDialogButtonBox(QDialogButtonBox.Ok)
        btn_box.accepted.connect(self.accept)

        # 1. Report Issue Button (GitHub)
        issue_btn = QPushButton("Report Issue")
        issue_btn.setToolTip("Open the GitHub Issues page to report bugs or request features")
        issue_btn.clicked.connect(self._open_github_issues)
        btn_box.addButton(issue_btn, QDialogButtonBox.ActionRole)
        
        # 2. Cite Button (ADS)
        cite_btn = QPushButton("Cite on ADS")
        cite_btn.setToolTip("Open the NASA ADS citation export page for Astrocook")
        cite_btn.clicked.connect(self._open_ads_citation)
        
        # Add "Cite" to the button box with ActionRole 
        # (This usually places it to the left of the OK button)
        btn_box.addButton(cite_btn, QDialogButtonBox.ActionRole)

        layout.addWidget(btn_box)

    def _open_ads_citation(self):
        """Opens the ADS export citation page in the default browser."""
        # URL for the 2020 Astronomy and Computing paper
        url = "https://ui.adsabs.harvard.edu/abs/2020SPIE11452E..1UC/exportcitation"
        QDesktopServices.openUrl(QUrl(url))

    def _open_github_issues(self):
        """Opens the GitHub Issues page in the default browser."""
        # Update this URL if the repo is under a different organization
        url = "https://github.com/DAS-OATs/astrocook/issues"
        QDesktopServices.openUrl(QUrl(url))