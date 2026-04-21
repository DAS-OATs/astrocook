# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

from astrocook import __version__ as version

myst_substitutions = {
    "version": version
}

# -- Path setup --------------------------------------------------------------
sys.path.insert(0, os.path.abspath('..'))

# -- Project information -----------------------------------------------------
project = 'Astrocook'
copyright = '2026, Guido Cupani'
author = 'Guido Cupani'
release = version

# -- General configuration ---------------------------------------------------
extensions = [
    'sphinx.ext.autodoc',      # Core library for html generation from docstrings
    'sphinx.ext.autosummary',  # Create neat summary tables
    'sphinx.ext.viewcode',     # Add links to highlighted source code
    'sphinx.ext.napoleon',     # Support for NumPy and Google style docstrings
    'myst_parser',             # Markdown support
    'sphinx_design',           # Design components
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']



# -- Options for HTML output -------------------------------------------------
html_theme = 'pydata_sphinx_theme'
html_static_path = ['_static']

# Add the custom CSS file
html_css_files = [
    'css/custom.css',
]

# -- MyST Parser Configuration -----------------------------------------------
# This allows you to use colon fences for directives (like ::: note)
myst_enable_extensions = [
    "colon_fence",    # Allows ::: for directives
    "dollarmath",     # Allows $$ math $$
    "deflist",        # Definition lists
    "html_image",     # Better image handling
    "substitution",
]

# -- PyData Theme Configuration ----------------------------------------------
# The path is relative to the configuration directory (docs/)
html_logo = "_static/icon_3d_HR.png"

# Optional: PyData Theme specific logo configuration
html_theme_options = {
    "github_url": "https://github.com/das-oats/astrocook",
    "show_nav_level": 2,
    "navigation_depth": 4,

    "logo": {
        "text": "Astrocook",  # Text to show next to the logo
        # You can specify different images for light/dark modes if you have them
        # "image_light": "_static/logo_3d.png",
        # "image_dark": "_static/logo_3d_dark.png",
    },
    
    # --- Add these lines ---
    "navbar_align": "content", # Align navbar items to the content
    "switcher": {
        # This must be the absolute URL to the JSON file you just uploaded
        "json_url": "https://das-oats.github.io/astrocook/switcher.json",
        
        # This tells the theme which version THIS build represents.
        "version_match": "2.0.1",
    },
    "navbar_end": ["theme-switcher", "version-switcher", "navbar-icon-links"],
}

# -- Options for autodoc ----------------------------------------------------
# This controls what appears in the API docs
autodoc_default_options = {
    'members': True,           # Document all members (functions/classes)
    'undoc-members': True,     # Even those without docstrings
    'show-inheritance': True,  # Show parent classes
    'special-members': '__init__', # Document the __init__ method
    # 'private-members': True, # Uncomment if you also want _functions
}

# Ensure order follows source code, not alphabetical
autodoc_member_order = 'bysource'