import os

from docutils.parsers.rst import Directive

from .nodes import chemiscope
from .utils import copy_file


class Chemiscope(Directive):
    """Directive to handle chemiscope visualizations in documentation

    e.g.::

        .. chemiscope::
            :filepath: [chemiscope_path]/docs/src/examples/datasets/fig_base_001.json.gz
            :mode: default

    The resulting html is the chemiscope widget wrapped in the
    chemiscope-sphinx.html template.
    """

    node_class = chemiscope
    has_content = True
    required_arguments = 0
    optional_arguments = 0
    final_argument_whitespace = False
    option_spec = {"filepath": str, "mode": str}

    def run(self):
        try:
            # Path to the saved dataset in the .rst files folder
            rst_file_path = self.options.get("filepath")

            # Copy dataset to the docs/build/html/_datasets folder
            filename = os.path.basename(rst_file_path)
            build_file_path, rel_file_path = self.get_build_file_path(filename)
            copy_file(rst_file_path, build_file_path)

            # Create the chemiscope node
            node = self.create_node(rel_file_path)
            return [node]
        except Exception as e:
            print(f"Error during run: {e}")
            return []

    def get_build_file_path(self, filename):
        """
        Construct the path to the build directory

        Parameters:
        - filename (str): The name of the file

        Returns:
        - tuple: A tuple containing the build file path and the relative file path
        """
        # Get the destination folder
        outdir = self.state.document.settings.env.app.outdir
        target_dir = os.path.join(outdir, "_datasets")
        os.makedirs(target_dir, exist_ok=True)

        # Get destination paths
        build_file_path = os.path.join(target_dir, filename)
        rel_file_path = os.path.relpath(build_file_path, outdir)
        return build_file_path, rel_file_path

    def create_node(self, rel_file_path):
        """
        Create a chemiscope node with the specified file path and mode

        Parameters:
        - rel_file_path (str): The dataset path relative to the build directory

        Returns:
        - chemiscope: The created chemiscope node
        """
        node = chemiscope()
        node["filepath"] = rel_file_path
        node["mode"] = self.options.get("mode")
        self.state.nested_parse(self.content, self.content_offset, node)
        return node
