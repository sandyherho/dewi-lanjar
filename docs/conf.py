"""Sphinx configuration for dewi-lanjar."""

import os
import sys

sys.path.insert(0, os.path.abspath("../src"))

project = "dewi-lanjar"
author = "Sandy H. S. Herho, Rima Rachmayani, Rusmawan Suwarman, Dasapta E. Irawan"
copyright = "2026, Applied Geology Research Group, Bandung Institute of Technology"
release = version = "0.0.1"

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.mathjax",
]
myst_enable_extensions = ["dollarmath", "amsmath", "colon_fence", "deflist"]
myst_heading_anchors = 3
source_suffix = {".md": "markdown", ".rst": "restructuredtext"}
master_doc = "index"
exclude_patterns = ["_build"]

autodoc_member_order = "bysource"
autodoc_typehints = "description"
napoleon_numpy_docstring = True

html_theme = "furo"
html_title = "dewi-lanjar 0.0.1"
html_static_path = ["_static"]
html_theme_options = {
    "source_repository": "https://github.com/sandyherho/dewi-lanjar",
    "source_branch": "main",
    "source_directory": "docs/",
}
