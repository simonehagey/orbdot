import os
import sys

sys.path.insert(0, os.path.abspath("../.."))

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx_togglebutton",
    "sphinxcontrib.bibtex",
    "sphinx_copybutton",
    "sphinx.ext.autosummary",
    "sphinx.ext.viewcode",
]

templates_path = ["_templates"]
source_suffix = ".rst"
master_doc = "index"
pygments_style = "sphinx"

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_logo = "_static/orbdot_logo.png"

html_theme_options = {
    "logo_only": False,
    "collapse_navigation": False,
    "sticky_navigation": True,
    "navigation_depth": -1,
    "prev_next_buttons_location": "both",
    "style_external_links": True,
    "style_nav_header_background": "white",
}

html_context = {
    "display_github": True,
    "github_user": "simonehagey",
    "github_repo": "orbdot",
    "github_version": "main",
    "conf_py_path": "/docs/source/",
}

napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_rtype = False

bibtex_bibfiles = ["references.bib"]
bibtex_reference_style = "author_year"
autosummary_generate = True
