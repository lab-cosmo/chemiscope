import os
import re
import warnings

import chemiscope
from chemiscope.sphinx import ChemiscopeScraper
from sphinx_gallery.sorting import ExplicitOrder

# Filter out warning to hide it in the documentation. It comes from the import of MACE
warnings.filterwarnings("ignore", message="Can't initialize NVML")

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
html_extra_path = []

examples_dirs = os.path.join(ROOT, "python", "examples")

def get_sorted_examples(directory):
    files = [f for f in os.listdir(directory) if f.endswith('.py')]
    
    def sort_key(filename):
        match = re.search(r'(\d+)', filename)
        if match:
            # numbered files
            return (0, int(match.group(1)))
    
        # unnumbered files
        return (1, filename)

    return sorted(files, key=sort_key)

sorted_examples = get_sorted_examples(examples_dirs)

sphinx_gallery_conf = {
    "examples_dirs": examples_dirs,
    "gallery_dirs": "examples",
    "filename_pattern": ".*",
    "within_subsection_order": ExplicitOrder(sorted_examples),
    "image_scrapers": ("matplotlib", ChemiscopeScraper()),
}

# -- Options for HTML output -------------------------------------------------

html_theme = "furo"

html_theme_options = {
    "light_css_variables": {
        "color-background-secondary": "#f4f8ff",
    },
    "dark_css_variables": {
        "color-background-secondary": "#303d4f",
    },
}

html_static_path = [os.path.join("..", "_static")]
html_favicon = "../_static/chemiscope-icon.png"
html_logo = "../_static/chemiscope-icon-no-name.svg"

templates_path = [os.path.join("..", "_templates")]


def setup(app):
    app.add_css_file("css/chemiscope.css")
