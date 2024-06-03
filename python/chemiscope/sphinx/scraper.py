import os

from chemiscope.jupyter import ChemiscopeWidget, MapWidget, StructureWidget

from .file_path_iterator import FilePathIterator


class ChemiscopeScraper:
    """Custom scraper for Chemiscope visualizations"""

    def __init__(self, examples_dir):
        if not examples_dir:
            raise ValueError("examples_dir must be provided.")

        # Get a target directory with the source files
        target_dir = os.path.join(examples_dir, "data")
        os.makedirs(target_dir, exist_ok=True)

        # Create an iterator to generate the file name
        self.iterator = FilePathIterator(target_dir)

    def __repr__(self):
        return "ChemiscopeScraper"

    def __call__(self, _block, block_vars, _gallery_conf):
        # Retrieve the chemiscope widget from block variables
        widget = block_vars.get("example_globals", {}).get("___")
        mode = self.get_widget_mode(widget)

        if mode is not None:
            dataset_file_path = self.iterator.next()
            widget.save(dataset_file_path)

            return f""".. chemiscope::
                :filename: {os.path.basename(dataset_file_path)}
                :mode: {mode}
            """
        else:
            return ""

    def get_widget_mode(self, widget):
        if type(widget) is ChemiscopeWidget:
            return "default"
        elif type(widget) is StructureWidget:
            return "structure"
        elif type(widget) is MapWidget:
            return "map"
        else:
            return None
