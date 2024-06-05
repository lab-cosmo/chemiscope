from docutils.parsers.rst import Directive

from .nodes import chemiscope
from .utils import copy_file_to_build_dir


class Chemiscope(Directive):
    """Directive to handle chemiscope visualizations in documentation"""

    node_class = chemiscope
    has_content = True
    required_arguments = 0
    optional_arguments = 0
    final_argument_whitespace = False
    option_spec = {"filename": str, "mode": str}

    def run(self):
        # Create the chemiscope node
        node = chemiscope()
        node["filename"] = self.options.get("filename")
        node["mode"] = self.options.get("mode")
        self.state.nested_parse(self.content, self.content_offset, node)

        # Copy the dataset file to the build directory
        try:
            app = self.state.document.settings.env.app
            copy_file_to_build_dir(app, node["filename"])
        except Exception as e:
            print(f"Error copying files: {e}")

        return [node]
