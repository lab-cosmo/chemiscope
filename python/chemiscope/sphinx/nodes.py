import os
import uuid

from docutils import nodes


class chemiscope(nodes.Element):
    pass


def visit_chemiscope_latex(self, node):
    pass


def depart_chemiscope_latex(self, node):
    pass


def visit_chemiscope_html(self, node):
    self.body.append(generate_html_content(node["filepath"], node["mode"]))


def depart_chemiscope_html(self, node):
    pass


def generate_html_content(filepath, mode="default"):
    """
    Generate HTML content for displaying a chemiscope visualization.

    Parameters:
    - filepath (str): The path to the generated dataset
    - mode (str): The display mode for the visualization

    Returns:
    - str: The chemiscope widget wrapped in the html
    """
    # Generate a unique id for the chemiscope div
    div_id = f"sphinx-gallery-{uuid.uuid4()}"

    # Read the template html file
    current_folder_path = os.path.dirname(__file__)
    template_path = os.path.join(
        current_folder_path, "static", "html", "chemiscope-sphinx-gallery.html"
    )
    with open(template_path, "r") as file:
        html_template = file.read()

    # Replace html placeholders with actual values
    return (
        html_template.replace("{{div_id}}", div_id)
        .replace("{{filepath}}", filepath)
        .replace("{{mode}}", mode)
    )
