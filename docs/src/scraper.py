import os
import shutil

from sphinx_gallery.scrapers import matplotlib_scraper

class SphinxGalleryScraper:
    def __repr__(self):
        return 'SphinxGalleryScraper'

    def __call__(self, block, block_vars, gallery_conf):
        variables = block_vars.get("example_globals", {})
        widget = variables.get("___")
        if widget:
            # Copy js scripts
            target_dir = gallery_conf.get("gallery_dirs")
            if not target_dir:
                raise ValueError("gallery_conf does not contain 'gallery_dirs' key.")
            source_files = [
                "./_static/js/chemiscope.min.js",
                "./_static/js/chemischope-sphinx-gallery.js"
            ]
            self.copy_chemiscope_files(target_dir, source_files)

            # Use the custom html content
            widget._repr_html_ = lambda: self.generate_html_content(widget.value)

        # Use sphinx-gallery standard scraper
        return matplotlib_scraper(block, block_vars, gallery_conf)

    def copy_chemiscope_files(self, target_dir, files):
        # Ensure target_dir exists
        if not os.path.exists(target_dir):
            os.makedirs(target_dir)

        # Copy each file to the target directory
        for file in files:
            source_path = os.path.abspath(file)
            target_path = os.path.join(target_dir, os.path.basename(file))
            shutil.copyfile(source_path, target_path)
            print(f"Copied {source_path} to {target_path}")

    def generate_html_content(self, dataset):
        div_id = "sphinx-gallery-" + str(id(self))
        return f"""
            <!-- Load all dependencies -->
            <!-- JQuery -->
            <script
                src="https://cdnjs.cloudflare.com/ajax/libs/jquery/3.4.1/jquery.min.js"
                integrity="sha256-CSXorXvZcTkaix6Yvo6HppcZGetbYMGWSFlBw8HfCJo=" crossorigin="anonymous"
            ></script>

            <!-- JS scripts -->
            <script src="../../../_static/js/chemiscope.min.js"></script>
            <script async="async" src="../../../_static/js/chemischope-sphinx-gallery.js"></script>

            <!-- Font-awesome -->
            <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.11.2/css/all.min.css" />
            <link rel="icon" type="image/png" href="../../../_static/chemiscope-icon.png" sizes="32x32" />

            <!-- CSS styles -->
            <link rel="stylesheet" href="../../../_static/chemischope-sphinx-gallery.css" type="text/css" />

            <!-- Warning Display -->
            <div id="{div_id}-warning-display"
                style="display: none;
                flex-wrap: wrap;
                border: 1px solid #f0ad4e;
                background-color: #fcf8e3;
                color: #8a6d3b;
                padding: 15px;
                border-radius: 4px;
                margin-bottom: 20px;">
                    <p style="flex: 1"></p>
                    <button
                        type="button"
                        style="position: relative; float: right; font-size: 3em; line-height: 1; color: inherit; text-shadow: none; background-color: transparent; border: 0; -webkit-appearance: none;"
                        aria-label="Close"
                        onclick="hideElement('{div_id}-warning-display')"
                    >&times;</button>`
            </div>
            <div style="padding: 0.5rem 1.25rem">&nbsp;</div>

            <!-- Loading Spinner -->
            <div id="{div_id}-loading">
                <i class="fa fa-spinner fa-spin"></i>
            </div>

            <!-- Chemiscope Visualization -->
            <div id="{div_id}"></div>

            <!-- Load Visualization Script -->
            <script>
                loadChemiscopeSphinxGallery("{div_id}", {dataset});
            </script>
        """
