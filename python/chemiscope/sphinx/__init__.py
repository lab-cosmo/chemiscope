from .directive import ChemiscopeDirective
from .nodes import (
    chemiscope,
    visit_chemiscope_html,
    depart_chemiscope_html,
    visit_chemiscope_latex,
    depart_chemiscope_latex,
)
from .utils import copy_static_folder


def setup(app):
    app.connect("build-finished", copy_static_folder)
    app.add_directive("chemiscope", ChemiscopeDirective)
    app.add_node(
        chemiscope,
        override=True,
        html=(visit_chemiscope_html, depart_chemiscope_html),
        latex=(visit_chemiscope_latex, depart_chemiscope_latex),
    )
