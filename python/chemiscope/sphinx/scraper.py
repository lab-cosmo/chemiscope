import os

from chemiscope.jupyter import ChemiscopeWidget, MapWidget, StructureWidget

from .file_path_iterator import FilePathIterator
from .utils import get_raw_filename


class ChemiscopeScraper:
    """Custom scraper for Chemiscope visualizations"""

    def __init__(self):
        # Create an iterator to generate the file name
        self.iterator = FilePathIterator()

    def __repr__(self):
        return "ChemiscopeScraper"

    def __call__(self, _block, block_vars, _gallery_conf):
        # Retrieve the chemiscope widget from block variables
        widget = block_vars.get("example_globals", {}).get("___")
        mode = self.get_widget_mode(widget)

        if mode is not None:
            # Get the target directory to save the dataset next to the .rst files
            src_file = block_vars.get("target_file")  # Python file
            rst_dataset_dir = os.path.join(os.path.dirname(src_file), "datasets")
            os.makedirs(rst_dataset_dir, exist_ok=True)

            # Save the related dataset to the directory next to the related .rst file
            infix = get_raw_filename(src_file)
            dataset_filename = self.get_dataset_filename(infix)
            dataset_file_path = os.path.join(rst_dataset_dir, dataset_filename)
            widget.save(dataset_file_path)

            # Return the .rst directive
            return self.generate_rst(dataset_file_path, mode)
        else:
            # Return an empty string if the widget type is not recognized
            return ""

    def generate_rst(self, filepath, mode):
        """Generate a .rst directive for embedding chemiscope widget"""
        return f""".. chemiscope::
            :filepath: {filepath}
            :mode: {mode}
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
