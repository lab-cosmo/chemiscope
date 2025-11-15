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


def copy_external_structures(src_json: str, dst_json: str) -> None:
    """Copy external structure JSON files referenced by `src_json`
    into the same `_datasets` tree as `dst_json`, preserving relative paths.
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
    dst_base = os.path.dirname(dst_json)

    for struct in structures:
        if not isinstance(struct, dict):
            continue

        data = struct.get("data")

        rel_path: str | None = None

        # Case 1: data is a simple string: "frames/frame_0001.json"
        if isinstance(data, str):
            rel_path = data

        # Case 2: data is a dict with a filename/path/file key
        elif isinstance(data, dict):
            candidate = data.get("filename") or data.get("path") or data.get("file")
            if isinstance(candidate, str):
                rel_path = candidate

        if rel_path is None:
            continue

        # Ignore absolute paths; only support relative ones
        if os.path.isabs(rel_path):
            continue

        # Normalise to avoid "../" nonsense
        rel_path = os.path.normpath(rel_path)

        # Source file: relative to *original* JSON location
        src_struct = os.path.join(src_base, rel_path)

        if not os.path.exists(src_struct):
            # Optional: log a warning if you have access to a logger
            # logger.warning("Missing external structure file %s", src_struct)
            continue

        # Destination: under the same `_datasets` dir as dst_json
        dst_struct = os.path.join(dst_base, rel_path)
        copy_file(src_struct, dst_struct)
