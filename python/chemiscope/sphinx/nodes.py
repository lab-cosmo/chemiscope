import os
import re
import uuid

from docutils import nodes


class chemiscope(nodes.Element):
    pass


def visit_chemiscope_latex(self, node):
    self.body.append(
        "\\\\(The chemiscope visualisation cannot be rendered to LaTeX)\\\\\n"
    )


def depart_chemiscope_latex(self, node):
    pass


def visit_chemiscope_html(self, node):

    # determines the path of the static resources
    # relative to the output document
    html_file_path = self.builder.get_outfilename(self.builder.current_docname)
    html_file_dir = os.path.dirname(html_file_path)
    # _static is the hardcoded output folder for static files
    static_dir = os.path.join(self.builder.outdir, "_static")
    static_path = os.path.relpath(static_dir, html_file_dir)

    self.body.append(
        generate_html_content(
            node["filepath"], node["mode"], node["include_headers"], static_path
        )
    )


def depart_chemiscope_html(self, node):
    pass


def generate_html_content(
    filepath, mode="default", include_headers=True, static_path="_static"
):
    """
    Generate HTML content for displaying a chemiscope visualization.

    Parameters:
    - filepath (str): The path to the generated dataset
    - mode (str): The display mode for the visualization

    Returns:
    - str: The chemiscope widget wrapped in the html
    """
    # Generate a unique id for the chemiscope div
    div_id = f"chsp-{uuid.uuid4()}"

    # Read the template html file
    current_folder_path = os.path.dirname(__file__)
    template_path = os.path.join(
        current_folder_path, "static", "chemiscope-sphinx.html"
    )
    with open(template_path, "r") as file:
        html_template = file.read()

    # Choose whether to discard headers
    if not include_headers:
        html_template = re.sub(
            r"<!-- begin headers -->.*?<!-- end headers -->",
            "",
            html_template,
            flags=re.DOTALL,
        )

    # Replace html placeholders  with actual values
    return (
        html_template.replace("{{div_id}}", div_id)
        .replace("{{filepath}}", filepath)
        .replace("{{mode}}", mode)
        .replace("{{static_path}}", static_path)
    )
