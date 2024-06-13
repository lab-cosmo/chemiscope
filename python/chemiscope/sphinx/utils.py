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
