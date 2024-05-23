import os
import shutil
from PIL import Image

from sphinx_gallery.scrapers import matplotlib_scraper

class ChemiscopeScraper:
    def __repr__(self):
        return "ChemiscopeScraper"

    def __call__(self, block, block_vars, gallery_conf):
        variables = block_vars.get("example_globals", {})
        iterator = block_vars.get("image_path_iterator")
        widget = variables.get("___")

        if widget:
            # Save the dataset to json
            path_png = iterator.next()
            path_json = path_png.replace(".png", ".json.gz")
            widget.save(path_json)

            # Work around to use the sphynx-gallery generated name for the non-image extension (.json.gz)
            self.save_empty_png(path_png)

            # Use the custom html content
            widget._repr_html_ = lambda: self.generate_html_content(os.path.basename(path_json))

        # Use sphinx-gallery standard scrapper
        return matplotlib_scraper(block, block_vars, gallery_conf)

    def generate_html_content(self, filename):
        # Generate the unique id
        div_id = f"sphinx-gallery-{id(self)}"

        # Read the template html file
        current_folder_path = os.path.dirname(__file__)
        template_path = os.path.join(current_folder_path, "static", "html", "chemischope-sphinx-gallery.html")
        with open(template_path, "r") as file:
            html_template = file.read()

        # Replace html placeholders with actual values
        return html_template.replace("{{div_id}}", div_id).replace("{{filename}}", filename)

    def save_empty_png(self, path):
        img = Image.new('RGBA', (1, 1), color=(0, 0, 0, 0))
        img.save(path)
        return path

def copy_additional_files(app, exception):
    if exception:
        return

    # TODO it might be not necessary to have the gallery_dirs configuration
    gallery_dirs = app.config.sphinx_gallery_conf.get("gallery_dirs")
    if gallery_dirs is None:
        print("'gallery_dirs' configuration is not defined.")
        return

    # Copy files from source to build directory
    try:
        copy_files(app.srcdir, app.outdir, gallery_dirs)
        copy_static_files(app.outdir)
    except Exception as e:
        print(f"Error copying files: {e}")

def copy_files(src_dir, build_dir, gallery_dirs):
    source_gallery_dir = os.path.join(src_dir, gallery_dirs)
    build_gallery_dir = os.path.join(build_dir, gallery_dirs)
    for root, _, files in os.walk(source_gallery_dir):
        for file in files:
            if file.endswith(".json.gz"):
                src_file = os.path.join(root, file)
                dst_file = os.path.join(build_gallery_dir, file)
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
