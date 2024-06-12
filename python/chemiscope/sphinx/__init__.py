import os

from .directive import ChemiscopeDirective
from .nodes import (
    chemiscope,
    visit_chemiscope_html,
    depart_chemiscope_html,
    visit_chemiscope_latex,
    depart_chemiscope_latex,
)


def setup(app):
    app.add_directive("chemiscope", ChemiscopeDirective)
    app.add_node(
        chemiscope,
        override=True,
        html=(visit_chemiscope_html, depart_chemiscope_html),
        latex=(visit_chemiscope_latex, depart_chemiscope_latex),
    )
    static_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "static"))
    app.config.html_static_path.append(static_path)
