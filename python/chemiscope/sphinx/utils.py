import gzip
import json
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


def copy_external_structures(src_json: str, destination: str) -> None:
    """
    Copy external structure JSON files referenced by ``src_json``
    into the same ``_datasets`` tree as ``destination``, preserving relative paths.
    """

    try:
        if src_json.endswith(".gz"):
            fd = gzip.open(src_json, "rt", encoding="utf8")
        else:
            fd = open(src_json, "r", encoding="utf8")
        dataset = json.load(fd)
        fd.close()
    except Exception:
        # if anything goes wrong reading the dataset, just skip
        print("Warning: could not read dataset JSON to copy external structures")
        return

    if not isinstance(dataset, dict):
        return

    structures = dataset.get("structures", [])
    if not isinstance(structures, list):
        return

    # return if there are no external structures
    if "data" not in structures[0]:
        return

    src_base = os.path.dirname(src_json)

    for struct in structures:
        if not isinstance(struct, dict):
            continue

        data = struct["data"]

        # Expecting data to be just the path of the structure.json
        if isinstance(data, str):
            rel_path = data
        else:
            raise ValueError(
                "unsupported external structure data format: "
                f"expected string but got {data.__class__.__name__}"
            )

        rel_path = os.path.normpath(rel_path)

        # Source file: relative to *original* JSON location
        src_struct = os.path.join(src_base, rel_path)

        if not os.path.exists(src_struct):
            raise FileNotFoundError(f"External structure file {src_struct} missing")

        # Destination: under the target destination folder
        dst_struct = os.path.join(destination, rel_path)
        copy_file(src_struct, dst_struct)
