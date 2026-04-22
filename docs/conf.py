# Configuration file for the Sphinx documentation builder.
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

# -- Path setup ---------------------------------------------------------------
# Allow Sphinx to find the pycmplot package (needed for autodoc)
sys.path.insert(0, os.path.abspath(".."))

# -- Project information -------------------------------------------------------
project = "pycmplot"
copyright = "2026, Kevin Esoh"
author = "Kevin Esoh"
release = "0.2.6"  # update to match PyPI version

# -- General configuration -----------------------------------------------------
extensions = [
    "sphinx.ext.autodoc",       # auto-generate docs from docstrings
    "sphinx.ext.autosummary",   # summary tables for modules/classes
    "sphinx.ext.napoleon",      # NumPy / Google docstring styles
    "sphinx.ext.viewcode",      # [source] links in API docs
    "sphinx.ext.intersphinx",   # cross-links to numpy, pandas, matplotlib docs
    "numpydoc",                 # richer NumPy-style rendering
    "nbsphinx",                 # embed Jupyter notebooks
    "sphinx_copybutton",        # copy-button on code blocks
    "myst_parser",              # allow Markdown (.md) pages alongside .rst
]

# Napoleon settings (NumPy docstring style)
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = False
napoleon_use_param = True
napoleon_use_rtype = True

# numpydoc settings
numpydoc_show_class_members = False

# Autosummary: auto-generate stub files
autosummary_generate = True
autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "show-inheritance": True,
}

# Intersphinx: link to external package docs
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "pandas": ("https://pandas.pydata.org/docs", None),
    "matplotlib": ("https://matplotlib.org/stable", None),
}

# nbsphinx: do not re-execute notebooks during docs build
nbsphinx_execute = "never"

# Source file suffixes
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "**.ipynb_checkpoints"]

# -- Options for HTML output ---------------------------------------------------
html_theme = "sphinx_rtd_theme"
html_theme_options = {
    "logo_only": False,
    "prev_next_buttons_location": "bottom",
    "style_nav_header_background": "#2980B9",
    "collapse_navigation": False,
    "sticky_navigation": True,
    "navigation_depth": 4,
    "includehidden": True,
    "titles_only": False,
}

html_static_path = ["_static"]
html_css_files = ["custom.css"]

# Optional: path to logo image (add docs/_static/logo.png if you have one)
# html_logo = "_static/logo.png"

html_show_sourcelink = True
html_show_sphinx = True
html_show_copyright = True
