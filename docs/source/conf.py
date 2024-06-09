import os
import sys


project = 'orbdot'
copyright = '2024, Simone R. Hagey'
author = 'Simone R. Hagey'
version = '0.1.0'
release = '0.1.0'

import os
import sys

sys.path.insert(0, os.path.abspath('../..'))
# sys.path.insert(0, os.path.abspath('../../orbdot/'))

# 'sphinx.ext.coverage', 'numpydoc',
extensions = ['sphinx.ext.autodoc', 'sphinx.ext.napoleon',
              'sphinx_togglebutton', 'sphinxcontrib.bibtex', 'sphinx_copybutton']

templates_path = ['_templates']
source_suffix = '.rst'
master_doc = 'index'
pygments_style = 'sphinx'

html_theme = "sphinx_rtd_theme"
html_static_path = ['_static']
html_logo = "_static/orbdot_logo.png"

html_theme_options = {
    'logo_only': False,
    'display_version': True,
    'collapse_navigation': False,
    'sticky_navigation': True,
    'navigation_depth': -1,
    'prev_next_buttons_location': 'both',
    'style_external_links': True,
    'style_nav_header_background': 'white',
}

napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True
napoleon_use_admonition_for_notes = True
napoleon_use_admonition_for_references = True

bibtex_bibfiles = ["references.bib"]