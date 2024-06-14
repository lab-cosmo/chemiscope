import os

from .directive import ChemiscopeDirective
from .nodes import (
    chemiscope,
    visit_chemiscope_html,
    depart_chemiscope_html,
    visit_chemiscope_latex,
    depart_chemiscope_latex,
)
from .scraper import ChemiscopeScraper  # noqa: F401


def setup(app):
    app.add_directive("chemiscope", ChemiscopeDirective)
    app.add_node(
        chemiscope,
        override=True,
        html=(visit_chemiscope_html, depart_chemiscope_html),
        latex=(visit_chemiscope_latex, depart_chemiscope_latex),
    )

    def add_static_path(app):
        static_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "static"))
        # put our path first, so the user can override our files if they want to
        app.config.html_static_path.insert(0, static_path)

    app.connect("builder-inited", add_static_path)
