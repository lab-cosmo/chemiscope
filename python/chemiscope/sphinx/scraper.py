import os
from PIL import Image
from chemiscope.jupyter import ChemiscopeWidget, MapWidget, StructureWidget


class ChemiscopeScraper:
    """Custom scraper for Chemiscope visualizations"""

    def __repr__(self):
        return "ChemiscopeScraper"

    def __call__(self, block, block_vars, gallery_conf):
        variables = block_vars.get("example_globals", {})
        iterator = block_vars.get("image_path_iterator")
        widget = variables.get("___")

        if widget:
            path_png = iterator.next()
            path_json = path_png.replace(".png", ".json.gz")
            widget.save(path_json)

            if type(widget) is ChemiscopeWidget:
                mode = "default"
            elif type(widget) is StructureWidget:
                mode = "structure"
            elif type(widget) is MapWidget:
                mode = "map"
            else:
                raise TypeError("Scraped widget is not a chemiscope widget")

            self.save_empty_png(path_png)

            return f""".. chemiscope::
                :filename: {os.path.basename(path_json)}
                :mode: {mode}
            """
        else:
            return ""

    def save_empty_png(self, path):
        img = Image.new("RGBA", (1, 1), color=(0, 0, 0, 0))
        img.save(path)
        return path
