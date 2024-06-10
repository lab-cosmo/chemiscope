import os
import shutil


def copy_file(src_file, dst_file):
    """
    Copy a single file from source to destination

    Parameters:
    - src_file (str): The path to the source file
    - dst_file (str): The path to the destination file
    """
    dst_dir = os.path.dirname(dst_file)
    os.makedirs(dst_dir, exist_ok=True)
    shutil.copyfile(src_file, dst_file)


def copy_static_folder(app, exception):
    """
    Copy static files required for the chemiscope widget

    Parameters:
    - app (Sphinx.application.Sphinx): The sphinx application object
    - exception (Exception): The exception that occurred during Sphinx build, if any
    """
    if exception:
        return
    current_file_dir = os.path.dirname(__file__)
    static_dir = os.path.join(current_file_dir, "static")
    build_static_dir = os.path.join(app.outdir, "static")

    # Copy the entire static directory to the destination directory
    try:
        shutil.copytree(static_dir, build_static_dir, dirs_exist_ok=True)
    except Exception as e:
        print(f"Error copying {static_dir} to {build_static_dir}: {e}")


def get_raw_filename(file_path):
    """
    Extract the filename from a file path without extension

    Parameters:
    - file_path (str): The full path of the file

    Returns:
    - str: The filename without extension
    """
    return os.path.splitext(os.path.basename(file_path))[0]
