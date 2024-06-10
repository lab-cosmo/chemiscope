import os
import shutil
import sysconfig


def copy_files_from_folder(src_dir, dest_dir, file_extension):
    """Copy files with a specific extension from source to destination"""
    for root, _, files in os.walk(src_dir):
        for file in files:
            if file.endswith(file_extension):
                src_file = os.path.join(root, file)
                dst_file = os.path.join(dest_dir, file)
                copy_file(src_file, dst_file)


def copy_file(src_file, dst_file):
    """Copy one file from source to destination"""
    dst_dir = os.path.dirname(dst_file)
    os.makedirs(dst_dir, exist_ok=True)
    shutil.copyfile(src_file, dst_file)


def copy_static_folder(app, exception):
    """Copy static files required for the chemiscope widget"""
    if exception:
        return
    current_file_dir = os.path.dirname(__file__)
    static_dir = os.path.join(current_file_dir, "static")
    files_to_copy = [
        "js/chemiscope-sphinx-gallery.js",
        "js/chemiscope.min.js",
        "loading-icon.svg",
        "css/chemiscope-sphinx-gallery.css",
    ]
    for file_name in files_to_copy:
        src_file = os.path.join(static_dir, file_name)
        dst_file = os.path.join(app.outdir, "static", file_name)
        try:
            copy_file(src_file, dst_file)
        except Exception as e:
            print(f"Error copying {src_file} to {dst_file}: {e}")


def copy_chemiscope_min_js(_app):
    """Copy the chemiscope.min.js  to the static directory"""
    src_file = get_chemiscope_src_file()
    current_file_dir = os.path.dirname(__file__)
    dst_file = os.path.join(current_file_dir, "static", "js", "chemiscope.min.js")
    try:
        copy_file(src_file, dst_file)
    except Exception as e:
        print(f"Error copying {src_file} to {dst_file}: {e}")


def get_chemiscope_src_file():
    """Find the path of get the path of the chemiscope.min.js"""
    prefix = get_install_prefix()
    chemiscope_dir = os.path.join(prefix, "share", "chemiscope")
    src_file = os.path.join(chemiscope_dir, "chemiscope.min.js")
    if not os.path.exists(src_file):
        print(f"Chemiscope file not found at: {src_file}")
        return
    return src_file


def get_install_prefix():
    """Get the installation prefix directory"""
    return sysconfig.get_paths()["data"]


def get_raw_filename(file_path):
    """
    Extract the filename from a file path without extension

    Parameters:
        file_path (str): The full path of the file

    Returns:
        str: The filename without extension
    """
    return os.path.splitext(os.path.basename(file_path))[0]
