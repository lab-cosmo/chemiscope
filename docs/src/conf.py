import os

import sphinx_bootstrap_theme

import chemiscope
from chemiscope.sphinx import ChemiscopeScraper

ROOT = os.path.abspath(os.path.join("..", ".."))

# -- Project information -----------------------------------------------------

project = "chemiscope"
copyright = "2021, Guillaume Fraux"
author = "Guillaume Fraux"

# The full version, including alpha/beta/rc tags
release = chemiscope.__version__

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx_argparse_cli",
    "sphinx_gallery.gen_gallery",
    "chemiscope.sphinx",
]

intersphinx_mapping = {
    "sphinx_gallery": ("https://sphinx-gallery.github.io/stable/", None),
}

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


sphinx_gallery_conf = {
    "examples_dirs": os.path.join(ROOT, "python", "examples"),
    "gallery_dirs": "examples",
    "filename_pattern": ".*",
    "within_subsection_order": "FileNameSortKey",
    "image_scrapers": ("matplotlib", ChemiscopeScraper()),
}

# -- Options for HTML output -------------------------------------------------

html_theme = "bootstrap"
html_theme_path = sphinx_bootstrap_theme.get_html_theme_path()

html_theme_options = {
    "navbar_links": [
        ("Documentation", "index.html", True),
    ],
    "navbar_pagenav": False,
    "bootswatch_theme": "flatly",
    "source_link_position": None,
}

html_sidebars = {"**": ["sidebar-toc.html", "searchbox.html"]}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = [os.path.join("..", "_static")]
html_favicon = "../_static/chemiscope-icon.png"

templates_path = [os.path.join("..", "_templates")]


def setup(app):
    app.add_css_file("css/chemiscope.css")
