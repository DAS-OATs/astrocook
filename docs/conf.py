# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

from astrocook import __version__

# Version of these docs: "dev" for the docs of master (published under
# /dev/), the package version for a release (published at the root); set
# by build_site.sh. It is shown in the pages, and the version switcher and
# its warning banner compare it with switcher.json: with the version of the
# last release, the theme would take the dev docs for the release docs and
# show no banner.
version = os.environ.get("ASTROCOOK_DOCS_VERSION", __version__)

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
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', 'screenshots']



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
    # At most four links in the header, the others under "More": with
    # more, the header wrapped onto two rows on a laptop screen (one row
    # from ~1280 px with four, ~1400 px with five, none with six)
    "header_links_before_dropdown": 4,
    "switcher": {
        # One list for all the versions, at the root of the site; its
        # "preferred" entry is the latest release (switcher.json)
        "json_url": "https://das-oats.github.io/astrocook/switcher.json",
        
        # This tells the theme which version THIS build represents.
        "version_match": version,
    },
    # Banner on the pages of any version other than the latest release
    # (the dev docs), with a link to it
    "show_version_warning_banner": True,
    "navbar_end": ["theme-switcher", "version-switcher", "navbar-icon-links"],
}

# -- Options for autodoc ----------------------------------------------------
# This controls what appears in the API docs
autodoc_default_options = {
    'members': True,           # Document all members (functions/classes)
    'undoc-members': True,     # Even those without docstrings
    'show-inheritance': True,  # Show parent classes
    'special-members': '__init__', # Document the __init__ method
    # Document only what a module defines: the re-exports listed in the
    # __all__ of a package (astrocook.core, astrocook.mocks) would be
    # described twice
    'ignore-module-all': True,
    # The methods of SpectrumV2 live in mixin classes, one per theme
    # (phase 5): its page lists them as before, and the references to them
    # resolve; the members of the bases outside Astrocook are left out
    'inherited-members': True,
    # 'private-members': True, # Uncomment if you also want _functions
}

# Ensure order follows source code, not alphabetical
autodoc_member_order = 'bysource'