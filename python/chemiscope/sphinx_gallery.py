import os
import shutil

from sphinx_gallery.scrapers import matplotlib_scraper

class ChemiscopeScraper:
    def __repr__(self):
        return "ChemiscopeScraper"

    def __call__(self, block, block_vars, gallery_conf):
        variables = block_vars.get("example_globals", {})
        widget = variables.get("___")
        if widget:
            # Save the dataset to json
            path_png = block_vars.get("image_path_iterator").next()
            path_json = path_png.replace(".png", ".json.gz")
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
            <script src="static/js/chemiscope.min.js"></script>
            <script src="static/js/chemischope-sphinx-gallery.js"></script>

            <!-- CSS styles -->
            <link rel="stylesheet" href="static/css/chemischope-sphinx-gallery.css" type="text/css" />

            <!-- Warning Display -->
            <div id="{div_id}-warning-display" class="chemiscope-sphinx-gallery-warning">
                <p></p>
                <button
                    type="button"
                    aria-label="Close"
                    onclick="hideElement('{div_id}-warning-display')"
                >&times;</button>
            </div>
            <div style="padding: 0.5rem 1.25rem" />

            <!-- Loading Spinner -->
            <div id="{div_id}-loading" class="chemiscope-sphinx-gallery-spinner">
                <img src="static/loading-icon.svg" alt="Loading icon" />
            </div>

            <!-- Chemiscope Visualization -->
            <div id="{div_id}" />

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
        copy_static_files(build_dir)
    except Exception as e:
        print(f"Error copying files: {e}")

def copy_files(src_dir, build_dir):
    for root, _, files in os.walk(src_dir):
        for file in files:
            if file.endswith(".json.gz"):
                src_file = os.path.join(root, file)
                dst_file = os.path.join(build_dir, file)
                copy_file(src_file, dst_file)

def copy_file(src_file, dst_file):
    dst_dir = os.path.dirname(dst_file)
    os.makedirs(dst_dir, exist_ok=True)
    shutil.copyfile(src_file, dst_file)

def copy_static_files(build_dir):
    current_file_dir = os.path.dirname(__file__)
    static_dir = os.path.join(current_file_dir, "static")
    files_to_copy = [
        "js/chemischope-sphinx-gallery.js", "js/chemiscope.min.js",
        "loading-icon.svg", "css/chemischope-sphinx-gallery.css"
    ]

    # Copy each static file to the static directory in the build directory
    for file_name in files_to_copy:
        src_file = os.path.join(static_dir, file_name)
        dst_file = os.path.join(build_dir, 'static', file_name)
        try:
            copy_file(src_file, dst_file)
        except Exception as e:
            print(f"Error copying {src_file} to {dst_file}: {e}")
