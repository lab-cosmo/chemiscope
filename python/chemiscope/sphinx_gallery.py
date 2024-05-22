import os
import shutil

from sphinx_gallery.scrapers import matplotlib_scraper

class ChemiscopeScraper:
    def __repr__(self):
        return 'ChemiscopeScraper'

    def __call__(self, block, block_vars, gallery_conf):
        variables = block_vars.get("example_globals", {})
        widget = variables.get("___")
        if widget:
            # Save the dataset to json
            path_png = block_vars.get("image_path_iterator").next()
            path_json = path_png.replace('.png', '.json.gz')
            widget.save(path_json)

            # Use the custom html content
            widget._repr_html_ = lambda: self.generate_html_content(os.path.basename(path_json))

        # Use sphinx-gallery standard scrapper
        return matplotlib_scraper(block, block_vars, gallery_conf)

    def generate_html_content(self, filename):
        div_id = "sphinx-gallery-" + str(id(self))
        return f"""
            <!-- Load all dependencies -->
            <!-- JQuery -->
            <script
                src="https://cdnjs.cloudflare.com/ajax/libs/jquery/3.4.1/jquery.min.js"
                integrity="sha256-CSXorXvZcTkaix6Yvo6HppcZGetbYMGWSFlBw8HfCJo=" crossorigin="anonymous"
            ></script>
            <script src="https://cdnjs.cloudflare.com/ajax/libs/pako/2.0.4/pako.min.js"></script>

            <!-- JS scripts -->
            <script src="../../../_static/js/chemiscope.min.js"></script>
            <script src="../../../_static/js/chemischope-sphinx-gallery.js"></script>

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
                loadChemiscopeSphinxGallery("{div_id}", "{filename}");
            </script>
        """

def copy_additional_files(app, exception):
    if exception:
        return

    # TODO it might be not necessary to have the gallery_dirs configuration
    gallery_dirs = app.config.sphinx_gallery_conf.get("gallery_dirs")
    if gallery_dirs is None:
        print("'gallery_dirs' configuration is not defined.")
        return

    # Copy files from source to build directory
    src_dir = os.path.join(app.srcdir, gallery_dirs)
    build_dir = os.path.join(app.outdir, gallery_dirs)
    try:
        copy_files(src_dir, build_dir)
    except Exception as e:
        print(f"Error copying files: {e}")

def copy_files(src_dir, build_dir):
    for root, _, files in os.walk(src_dir):
        for file in files:
            if file.endswith('.json.gz'):
                src_file = os.path.join(root, file)
                dst_file = os.path.join(build_dir, file)
                copy_file(src_file, dst_file)

def copy_file(src_file, dst_file):
    dst_dir = os.path.dirname(dst_file)
    os.makedirs(dst_dir, exist_ok=True)
    shutil.copyfile(src_file, dst_file)
    print(f"Copied {src_file} to {dst_file}")
