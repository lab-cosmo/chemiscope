import os
import shutil
import sysconfig


def copy_additional_files(app, exception):
    if exception:
        return
    gallery_dirs = app.config.sphinx_gallery_conf.get("gallery_dirs")
    examples_dirs = app.config.sphinx_gallery_conf.get("examples_dirs")
    src_gallery_dir = os.path.join(examples_dirs, "data")
    build_gallery_dir = os.path.join(app.outdir, gallery_dirs)
    try:
        copy_files_from_folder(src_gallery_dir, build_gallery_dir, ".json.gz")
        copy_static_files(build_gallery_dir)
    except Exception as e:
        print(f"Error copying files: {e}")


def copy_files_from_folder(src_dir, dest_dir, file_extension):
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


def get_chemiscope_src_file():
    prefix = get_install_prefix()
    chemiscope_dir = os.path.join(prefix, "share", "chemiscope")
    src_file = os.path.join(chemiscope_dir, "chemiscope.min.js")
    if not os.path.exists(src_file):
        print(f"Chemiscope file not found at: {src_file}")
        return
    return src_file


def get_install_prefix():
    return sysconfig.get_paths()["data"]
