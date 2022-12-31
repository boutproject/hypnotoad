# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
import hypnotoad

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "hypnotoad"
copyright = "2022, J.T. Omotani, B.D. Dudson and the hypnotoad team"
author = "J.T. Omotani, B.D. Dudson and the hypnotoad team"
release = hypnotoad.__version__

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

sys.path.insert(0, os.path.abspath(".."))

# Are we running on readthedocs?
on_rtd = os.environ.get("READTHEDOCS") == "True"

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
]

# templates_path = ['_templates']
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

if on_rtd:
    html_theme = "default"
else:
    html_theme = "sphinx_rtd_theme"
# html_static_path = ['_static']
