import os
import sys

# -- Project information -------------------------------------------------------
project = 'LBCgo'
author = 'Chris Howk'
copyright = '2024, Chris Howk'
release = '0.1.6'

# -- General configuration -----------------------------------------------------
extensions = [
    'sphinx.ext.napoleon',   # NumPy/Google docstring parsing
    'autoapi.extension',     # Static API doc generation (no import needed)
    'sphinx.ext.viewcode',   # Source code links in API docs
]

# sphinx-autoapi: point at the package source, generate docs statically
autoapi_dirs = ['../LBCgo']
autoapi_type = 'python'
# Exclude files with dots in their names (not valid Python modules)
autoapi_ignore = ['*/whatidid.astrom.py', '*.bak', '*/examples/*']
autoapi_options = [
    'members',
    'undoc-members',
    'show-inheritance',
    'show-module-summary',
]
# Don't show private (_foo) or special (__foo__) members
autoapi_member_order = 'bysource'
autoapi_python_class_content = 'both'

# Napoleon settings for NumPy docstrings
napoleon_numpy_docstring = True
napoleon_google_docstring = False
napoleon_include_private_with_doc = False

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- HTML output ---------------------------------------------------------------
html_theme = 'furo'
html_static_path = ['_static']

# Furo options
html_theme_options = {
    "sidebar_hide_name": False,
}
