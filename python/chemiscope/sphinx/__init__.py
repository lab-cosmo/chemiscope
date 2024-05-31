from .directive import Chemiscope
from .nodes import (
    chemiscope,
    visit_chemiscope_html,
    depart_chemiscope_html,
    visit_chemiscope_latex,
    depart_chemiscope_latex,
)
from .utils import copy_chemiscope_min_js, copy_additional_files
from .scraper import ChemiscopeScraper  # noqa: F401


def setup(app):
    app.connect("builder-inited", copy_chemiscope_min_js)
    app.connect("build-finished", copy_additional_files)
    app.add_directive("chemiscope", Chemiscope)
    app.add_node(
        chemiscope,
        override=True,
        html=(visit_chemiscope_html, depart_chemiscope_html),
        latex=(visit_chemiscope_latex, depart_chemiscope_latex),
    )
