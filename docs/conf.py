# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

# -- Path setup --------------------------------------------------------------
sys.path.insert(0, os.path.abspath('..'))

# -- Project information -----------------------------------------------------
project = 'Astrocook'
copyright = '2025, Guido Cupani'
author = 'Guido Cupani'
release = '2.0.0'

# -- General configuration ---------------------------------------------------
extensions = [
    'sphinx.ext.autodoc',      # Core library for html generation from docstrings
    'sphinx.ext.autosummary',  # Create neat summary tables
    'sphinx.ext.viewcode',     # Add links to highlighted source code
    'sphinx.ext.napoleon',     # Support for NumPy and Google style docstrings
    'myst_parser',             # Markdown support
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']



# -- Options for HTML output -------------------------------------------------
html_theme = 'pydata_sphinx_theme'
html_static_path = ['_static']

# -- MyST Parser Configuration -----------------------------------------------
# This allows you to use colon fences for directives (like ::: note)
myst_enable_extensions = [
    "colon_fence",
    "dollarmath",
]

# -- PyData Theme Configuration ----------------------------------------------
html_theme_options = {
    "github_url": "https://github.com/das-oats/astrocook",
    "show_nav_level": 2,
    "navigation_depth": 4,
}