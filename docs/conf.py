# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "SNF-simulations"
copyright = "2026, Zuzanna Leliwa, Abigail Power, Liz Kneale, Martin Dyer"
author = "Zuzanna Leliwa, Abigail Power, Liz Kneale, Martin Dyer"
from importlib.metadata import version as get_version
release = get_version("snf_simulations")
release = ".".join(release.split('.')[:3])  # Only use major.minor.patch for docs

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "pydata_sphinx_theme"
html_static_path = ["_static"]

# -- Extensions configuration ------------------------------------------------

extensions = [
#    "myst_parser",  # Markdown support
    "myst_nb",  # MyST Notebooks support (also includes myst_parser)
    "autodoc2",  # Automatic API documentation generation
    "sphinx.ext.napoleon",  # Parse Google and NumPy style docstrings
    "sphinx.ext.viewcode",  # Add links to source code
]

# -- Options for MyST Parser -------------------------------------------------
# https://myst-parser.readthedocs.io/en/latest/index.html

myst_enable_extensions = [
    "fieldlist",
    "dollarmath",
    "colon_fence",
]

# -- Options for autodoc2 ---------------------------------------------------
# https://sphinx-autodoc2.readthedocs.io/en/latest/index.html

autodoc2_packages = [
    "../src/snf_simulations",
]
autodoc2_render_plugin = "myst"

autodoc2_sort_names = True
autodoc2_class_docstring = "both"
autodoc2_hidden_objects = {"inherited", "dunder"}

# Below is needed for autodoc2_docstrings_parser,
# see https://github.com/sphinx-extensions2/sphinx-autodoc2/issues/33#issuecomment-2684177928
import os
import sys
sys.path.insert(0, os.path.abspath("."))
autodoc2_docstring_parser_regexes = [
    (r".*", "autodoc2_docstrings_parser"),
]

# -- Options for Napoleon -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/extensions/napoleon.html

napoleon_google_docstring = True
napoleon_use_ivar = True
