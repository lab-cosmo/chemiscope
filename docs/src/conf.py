import sphinx_bootstrap_theme
import os
import sys

# -- Project information -----------------------------------------------------

project = 'chemiscope'
copyright = '2020, Guillaume Fraux'
author = 'Guillaume Fraux'

# The full version, including alpha/beta/rc tags
release = '0.0'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
]

sys.path.append(os.path.join(os.path.dirname(__file__), "..", "..", "utils"))

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

html_theme = 'bootstrap'
html_theme_path = sphinx_bootstrap_theme.get_html_theme_path()

html_theme_options = {
    'navbar_links': [
        ("See it in action!", "https://chemiscope.org/", True),
    ],
    'navbar_pagenav': False,
    'bootswatch_theme': "flatly",
    'source_link_position': None,
}

html_sidebars = {
    '**': ['sidebar-toc.html', 'searchbox.html']
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = [os.path.join('..', '_static')]

templates_path = [os.path.join('..', '_templates')]


def setup(app):
    app.add_stylesheet("css/chemiscope.css")
