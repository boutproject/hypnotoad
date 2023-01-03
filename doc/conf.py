# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
from pathlib import Path
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
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.napoleon",
]

# templates_path = ['_templates']
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

autosectionlabel_prefix_document = True


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

if on_rtd:
    html_theme = "default"
else:
    html_theme = "sphinx_rtd_theme"

# These folders are copied to the documentation's HTML output
html_static_path = ["_static"]

# These paths are either relative to html_static_path
# or fully qualified paths (eg. https://...)
html_css_files = [
    "custom.css",
]

# -- Create temporary files with programatically generated content that can be included
# in the .rst files
tempdir = Path("_temp")
if not tempdir.exists():
    tempdir.mkdir()


def reformat_Texttable(tt_string):
    """
    Adds a 'blank line' between lines of multi-line entries in the table, so that when
    sphinx renders the ReStructuredText, the entry actually covers multiple lines.
    """
    tt_lines = tt_string.split("\n")
    blank_line = tt_lines[0].replace("+", "|").replace("-", " ")
    result = ""
    is_first = True
    for line in tt_lines:
        if line[0] == "+":
            result = result + line + "\n"
            is_first = True
        else:
            if not is_first:
                result = result + blank_line + "\n"
            result = result + line + "\n"
            is_first = False
    return result


def create_options_rst(filename, title, options_factory):
    with open(tempdir.joinpath(filename), "w") as f:
        f.write(":orphan:\n\n")
        f.write(
            ".. This file is auto-generated. To update the content, edit the "
            "OptionsFactory constructors in the hypnotoad source code. To change the "
            "formatting, edit the ``create_options_rst()`` function in ``conf.py``.\n\n"
        )
        f.write(f"{title}\n{'='*len(title)}\n\n")
        tt = options_factory.get_help_table(as_Texttable=True)
        tt.set_cols_width([62, 80, 20])
        f.write(reformat_Texttable(tt.draw()))


create_options_rst(
    "options.rst",
    "Tokamak options",
    hypnotoad.tokamak.TokamakEquilibrium.user_options_factory,
)
create_options_rst(
    "nonorthogonal-options.rst",
    "Nonorthogonal options",
    hypnotoad.tokamak.TokamakEquilibrium.nonorthogonal_options_factory,
)
create_options_rst(
    "circular-options.rst",
    "Circular options",
    hypnotoad.circular.CircularEquilibrium.user_options_factory,
)
create_options_rst(
    "torpex-options.rst",
    "TORPEX options",
    hypnotoad.torpex.TORPEXMagneticField.user_options_factory,
)
