import os

from docutils.parsers.rst import Directive

from .nodes import chemiscope
from .utils import copy_file


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
        filePath = self.options.get("filename")
        node["filename"] = self.get_rel_dest_path(filePath)
        node["mode"] = self.options.get("mode")
        self.state.nested_parse(self.content, self.content_offset, node)

        # Copy the dataset file to the build directory
        try:
            app = self.state.document.settings.env.app
            dst_path = os.path.join(app.outdir, node["filename"])
            copy_file(filePath, dst_path)
        except Exception as e:
            print(f"Error copying files: {e}")

        return [node]

    def get_rel_dest_path(self, path):
        """Get the last two folders and the filename from a given path"""
        path = path.rstrip(os.sep)
        path, data_dir = os.path.split(path)
        path, examples_dir = os.path.split(path)
        result = os.path.join(examples_dir, data_dir)
        return result
