import os
import warnings

from chemiscope.jupyter import ChemiscopeWidget, MapWidget, StructureWidget

from ._file_path_iterator import FilePathIterator
from ._utils import copy_external_structures


warnings.filterwarnings(
    "ignore",
    message="chemiscope.show only works in a jupyter notebook or a sphinx build",
)


class ChemiscopeScraper:
    """Custom scraper for Chemiscope visualizations"""

    def __init__(self):
        # Create an iterator to generate the file name
        self.iterator = FilePathIterator()

        # Process the widgets once by storing their `id()` here
        self.seen = set()

    def __repr__(self):
        return "ChemiscopeScraper"

    def __call__(self, _block, block_vars, _gallery_conf):
        """
        Extracts Chemiscope widget data, saves its .json dataset in the .rst directory
        and generates a .rst directive for embedding. Also copies external structures
        if present.

        Triggered on the each output of the script.

        Parameters:
        - _block: Unused parameter, a .rst code block
        - block_vars (dict): Variables from the executed code block
        - _gallery_conf: Unused parameter, sphinx-gallery config

        Returns:
        - str: .rst directive for embedding the chemiscope visualization
        """
        # Retrieve the chemiscope widget from block variables
        widget = block_vars.get("example_globals", {}).get("___")
        mode = self.get_widget_mode(widget)

        if (id(widget)) not in self.seen and mode is not None:
            self.seen.add(id(widget))

            # Get the target directory to save the dataset next to the .rst files
            target_file = block_vars.get("target_file")  # Python file
            rst_dataset_dir = os.path.join(os.path.dirname(target_file), "_datasets")
            os.makedirs(rst_dataset_dir, exist_ok=True)

            # Save the related dataset to the directory next to the related .rst file
            infix = os.path.splitext(os.path.basename(target_file))[0]
            dataset_filename = self.get_dataset_filename(infix)
            dataset_file_path = os.path.join(rst_dataset_dir, dataset_filename)
            widget.save(dataset_file_path)

            # also saves in-place to facilitate copying external structures
            src_file = block_vars.get("src_file")
            src_dataset = os.path.join(os.path.dirname(src_file), dataset_filename)
            widget.save(src_dataset)

            copy_external_structures(src_dataset, rst_dataset_dir)

            rel_file_path = os.path.join("_datasets", dataset_filename)

            # Return the .rst directive
            return self.generate_rst(rel_file_path, mode)
        else:
            # chemiscope.show is not called
            return ""

    def generate_rst(self, filepath, mode):
        """Generate a .rst directive for embedding chemiscope widget"""
        return f""".. chemiscope:: {filepath}
            :mode: {mode}
            :warning_timeout: {2000}
        """

    def get_widget_mode(self, widget):
        """Determine the mode of the chemiscope widget"""
        if isinstance(widget, ChemiscopeWidget):
            return "default"
        elif isinstance(widget, StructureWidget):
            return "structure"
        elif isinstance(widget, MapWidget):
            return "map"
        else:
            return None

    def get_dataset_filename(self, infix):
        """Generate a dataset filename using the file path infix"""
        self.iterator.set_infix(infix)
        return self.iterator.next()
