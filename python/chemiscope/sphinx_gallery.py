import os
import shutil
from PIL import Image
import warnings
import uuid
import sysconfig
from chemiscope.jupyter import ChemiscopeWidget, StructureWidget, MapWidget

from sphinx_gallery.scrapers import matplotlib_scraper

warnings.filterwarnings(
    "ignore",
    category=UserWarning,
    message="chemiscope.show only works inside a jupyter notebook",
)


def setup(app):
    setup_image_scrapers(app)
    app.connect("builder-inited", copy_chemiscope_min_js)
    app.connect("build-finished", copy_additional_files)


def setup_image_scrapers(app):
    # Check if the 'image_scrapper' was already added to conf
    sphinx_gallery_conf = app.config.sphinx_gallery_conf
    if "image_scrapers" not in sphinx_gallery_conf:
        sphinx_gallery_conf["image_scrapers"] = ()

    # Add custom scrapper
    sphinx_gallery_conf["image_scrapers"] += (ChemiscopeScraper(),)
    app.config.sphinx_gallery_conf.update(sphinx_gallery_conf)


class ChemiscopeScraper:
    """Custom scraper for Chemiscope visualizations"""

    def __repr__(self):
        return "ChemiscopeScraper"

    def __call__(self, block, block_vars, gallery_conf):
        """Process code blocks to generate chemiscope widget"""
        variables = block_vars.get("example_globals", {})
        iterator = block_vars.get("image_path_iterator")
        widget = variables.get(
            "___"
        )

        if widget:
            # Save the dataset to json
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

            # Work around to use the sphynx-gallery generated name for the non-image extension (.json.gz)
            self.save_empty_png(path_png)

            # Use the custom html content
            widget._repr_html_ = lambda: self.generate_html_content(
                os.path.basename(path_json), mode
            )

        # Use sphinx-gallery standard scrapper
        return matplotlib_scraper(block, block_vars, gallery_conf)

    def generate_html_content(self, filename, mode="default"):
        # Generate a unique id for the chemiscope div
        div_id = f"sphinx-gallery-{uuid.uuid4()}"

        # Read the template html file
        current_folder_path = os.path.dirname(__file__)
        template_path = os.path.join(
            current_folder_path, "static", "html", "chemischope-sphinx-gallery.html"
        )
        with open(template_path, "r") as file:
            html_template = file.read()
            
        chemiscope_src_file = get_chemiscope_src_file()

        # Replace html placeholders with actual values
        return (
            html_template.replace("{{div_id}}", div_id)
            .replace("{{filename}}", filename)
            .replace("{{mode}}", mode)
        )

    def save_empty_png(self, path):
        img = Image.new("RGBA", (1, 1), color=(0, 0, 0, 0))
        img.save(path)
        return path


def copy_additional_files(app, exception):
    """Copy additional files after the build is finished"""
    if exception:
        return

    gallery_dirs = app.config.sphinx_gallery_conf.get("gallery_dirs")
    if gallery_dirs is None:
        print("'gallery_dirs' configuration is not defined.")
        return

    src_gallery_dir = os.path.join(app.srcdir, gallery_dirs)
    build_gallery_dir = os.path.join(app.outdir, gallery_dirs)

    # Copy files from source to build directory
    try:
        copy_files_from_folder(src_gallery_dir, build_gallery_dir, ".json.gz")
        copy_static_files(build_gallery_dir)
    except Exception as e:
        print(f"Error copying files: {e}")


def copy_files_from_folder(src_dir, dest_dir, file_extension):
    """Copy files from source to build directory"""
    for root, _, files in os.walk(src_dir):
        for file in files:
            if file.endswith(file_extension):
                src_file = os.path.join(root, file)
                dst_file = os.path.join(dest_dir, file)
                copy_file(src_file, dst_file)


def copy_file(src_file, dst_file):
    dst_dir = os.path.dirname(dst_file)
    os.makedirs(dst_dir, exist_ok=True)
    shutil.copyfile(src_file, dst_file)


def copy_static_files(build_gallery_dir):
    current_file_dir = os.path.dirname(__file__)
    static_dir = os.path.join(current_file_dir, "static")
    files_to_copy = [
        "js/chemischope-sphinx-gallery.js",
        "js/chemiscope.min.js",
        "loading-icon.svg",
        "css/chemischope-sphinx-gallery.css",
    ]
    # Copy each static file to the static directory in the build directory
    for file_name in files_to_copy:
        src_file = os.path.join(static_dir, file_name)
        dst_file = os.path.join(build_gallery_dir, "static", file_name)
        try:
            copy_file(src_file, dst_file)
        except Exception as e:
            print(f"Error copying {src_file} to {dst_file}: {e}")


def copy_chemiscope_min_js(_app):
    src_file = get_chemiscope_src_file()
    current_file_dir = os.path.dirname(__file__)
    dst_file = os.path.join(current_file_dir, "static", "js", "chemiscope.min.js")
    try:
        copy_file(src_file, dst_file)
    except Exception as e:
        print(f"Error copying {src_file} to {dst_file}: {e}")


def get_install_prefix():
    return sysconfig.get_paths()["data"]


def get_chemiscope_src_file():
    prefix = get_install_prefix()
    chemiscope_dir = os.path.join(prefix, "share", "chemiscope")
    src_file = os.path.join(chemiscope_dir, "chemiscope.min.js")
    if not os.path.exists(src_file):
        print(f"Chemiscope file not found at: {src_file}")
        return
    return src_file
