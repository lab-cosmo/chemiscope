import hashlib
import os

from docutils.parsers.rst import Directive

from ._nodes import chemiscope
from ._utils import copy_external_structures, copy_file


class ChemiscopeDirective(Directive):
    """Custom RST directive to include a chemiscope visualizer in a
    sphinx documentation file. It has two options: `filepath` - the path to a chemiscope
    JSON or gzipped JSON file, relative to the path of the RST file, and `mode`, which
    can be `default`, `structure` or `map` depending on the desired type of
    visualization. `warning_timeout` specifies the time in milliseconds to display
    a warning (set it to 0 to make them persistent, and to -1 to disable).

    e.g.::

        .. chemiscope:: datasets/polarizability.json.gz
            :mode: map
            :warning_timeout: 1000

    The resulting html is the chemiscope widget wrapped in the chemiscope-sphinx.html
    template.
    """

    node_class = chemiscope
    has_content = True
    required_arguments = 1
    optional_arguments = 0
    final_argument_whitespace = False
    option_spec = {"mode": str, "warning_timeout": float}
    pages_with_headers = []

    def run(self):
        # Get the source path from the document
        source = self.state.document["source"]
        source_path = os.path.dirname(source) + "/"

        # Path to the saved dataset in the .rst files folder
        dataset_rel_path = self.arguments[0].strip()
        dataset_path = source_path + dataset_rel_path

        # Ensure unique file name to avoid clashes if the same file name is used in
        # different directives.
        with open(dataset_path, "rb") as fd:
            shasum = hashlib.sha256(fd.read()).hexdigest()
        filename = f"{shasum}-{os.path.basename(dataset_path)}"

        # Copy dataset to the docs/build/html/_datasets folder
        build_file_path, rel_file_path = self.get_build_file_path(filename)
        copy_file(dataset_path, build_file_path)
        # moves files referenced by external structures in a folder with
        # a name based on the dataset file
        copy_external_structures(dataset_path, f"{build_file_path}-ext/")

        # Create the chemiscope node
        node = self.create_node(
            rel_file_path,
            include_headers=source not in ChemiscopeDirective.pages_with_headers,
        )

        # indicates that the page has already been added headers
        if source not in ChemiscopeDirective.pages_with_headers:
            ChemiscopeDirective.pages_with_headers.append(source)
        return [node]

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

        # Get path of output dataset relative to output HTML
        env = self.state.document.settings.env
        builder = env.app.builder
        html_file_path = builder.get_outfilename(env.docname)
        html_file_dir = os.path.dirname(html_file_path)

        # Relative path *for the output files*
        rel_file_path = os.path.relpath(build_file_path, html_file_dir)
        return build_file_path, rel_file_path

    def create_node(self, rel_file_path, include_headers=True):
        """
        Create a chemiscope node with the specified file path and mode

        Parameters:
        - rel_file_path (str): The dataset path relative to the build directory

        Returns:
        - chemiscope: The created chemiscope node
        """

        node = chemiscope()
        node["filepath"] = rel_file_path
        node["mode"] = self.options.get("mode", "default")
        node["warning_timeout"] = self.options.get("warning_timeout", 4000)
        node["include_headers"] = include_headers

        self.state.nested_parse(self.content, self.content_offset, node)

        return node
