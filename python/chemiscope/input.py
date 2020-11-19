# -*- coding: utf-8 -*-
"""
Generate JSON input files for the default chemiscope visualizer.
"""
import warnings
import json
import gzip
import os

import numpy as np

from .adapters import frames_to_json, atom_properties, structure_properties


def _typetransform(data):
    assert isinstance(data, list) and len(data) > 0
    if isinstance(data[0], str):
        return list(map(str, data))
    elif isinstance(data[0], bytes):
        return list(map(lambda u: u.decode("utf8"), data))
    else:
        try:
            return [float(value) for value in data]
        except ValueError:
            raise Exception("unsupported type in property")


def _linearize(name, property):
    """
    Transform 2D arrays in multiple 1D arrays, converting types to fit json as
    needed.
    """
    data = {}
    if isinstance(property["values"], list):
        data[name] = {
            "target": property["target"],
            "values": _typetransform(property["values"]),
        }
    elif isinstance(property["values"], np.ndarray):
        if len(property["values"].shape) == 1:
            data[name] = {
                "target": property["target"],
                "values": _typetransform(list(property["values"])),
            }
        elif len(property["values"].shape) == 2:
            if property["values"].shape[1] == 1:
                data[name] = {
                    "target": property["target"],
                    "values": _typetransform(list(property["values"])),
                }
            else:
                for i in range(property["values"].shape[1]):
                    data[f"{name}[{i + 1}]"] = {
                        "target": property["target"],
                        "values": _typetransform(list(property["values"][:, i])),
                    }
        else:
            raise Exception("unsupported ndarray property")
    else:
        raise Exception(f"unknown type for property {name}")

    return data


def _generate_environments(frames, cutoff):
    environments = []
    for frame_id, frame in enumerate(frames):
        for center in range(len(frame)):
            environments.append(
                {"structure": frame_id, "center": center, "cutoff": cutoff}
            )
    return environments


def create_input(frames, meta=None, extra=None, cutoff=None):
    """
    Create a dictionary that can be saved to JSON using the format used by
    the default chemiscope visualizer.

    :param list frames: list of atomic structures. For now, only `ase.Atoms`_
                        objects are supported
    :param dict meta: optional metadata of the dataset, see below
    :param dict extra: optional dictionary of additional properties, see below
    :param float cutoff: optional. If present, will be used to generate
                         atom-centered environments

    The dataset metadata should be given in the ``meta`` dictionary, the
    possible keys are:

    .. code-block:: python

        meta = {
            'name': '...',         # str, dataset name
            'description': '...',  # str, dataset description
            'authors': [           # list of str, dataset authors, OPTIONAL
                '...',
            ],
            'references': [        # list of str, references for this dataset,
                '...',             # OPTIONAL
            ],
        }

    The returned dictionary will contain all the properties defined on the
    `ase.Atoms`_ objects. Values in ``ase.Atoms.arrays`` are mapped to
    ``target = "atom"`` properties; while values in ``ase.Atoms.info`` are
    mapped to ``target = "structure"`` properties. The only exception is
    ``ase.Atoms.arrays["numbers"]``, which is always ignored. If you want to
    have the atomic numbers as a property, you should add it to ``extra``
    manually.

    Additional properties can be added with the ``extra`` parameter. This
    parameter should be a dictionary containing one entry for each property.
    Each entry contains a ``target`` attribute (``'atom'`` or ``'structure'``)
    and a set of values. ``values`` can be a Python list of float or string; a
    1D numpy array of numeric values; or a 2D numpy array of numeric values. In
    the later case, multiple properties will be generated along the second axis.
    For example, passing

    .. code-block:: python

        extra = {
            'cheese': {
                'target': 'atom',
                'values': np.zeros((300, 4)),
                # optional: property unit
                'unit': 'random / fs',
                # optional: property description
                'description': 'a random property for example',
            }
        }

    will generate four properties named ``cheese[1]``, ``cheese[2]``,
    ``cheese[3]``,  and ``cheese[4]``, each containing 300 values.

    .. _`ase.Atoms`: https://wiki.fysik.dtu.dk/ase/ase/atoms.html
    """

    data = {"meta": {}}

    if meta is not None:
        if "name" in meta:
            data["meta"]["name"] = str(meta["name"])

        if "description" in meta:
            data["meta"]["description"] = str(meta["description"])

        if "authors" in meta:
            data["meta"]["authors"] = list(map(str, meta["authors"]))

        if "references" in meta:
            data["meta"]["references"] = list(map(str, meta["references"]))

        for key in meta.keys():
            if key not in ["name", "description", "authors", "references"]:
                warnings.warn("ignoring unexpected metadata: {}".format(key))

    if "name" not in data["meta"] or not data["meta"]["name"]:
        data["meta"]["name"] = "<unknown>"

    properties = {}
    if extra is not None:
        for name, value in extra.items():
            properties.update(_linearize(name, value))

    # Read properties coming from the frames
    for name, value in atom_properties(frames):
        properties.update(_linearize(name, value))

    for name, value in structure_properties(frames):
        properties.update(_linearize(name, value))

    data["properties"] = properties
    data["structures"] = frames_to_json(frames)

    if cutoff is not None:
        data["environments"] = _generate_environments(frames, cutoff)

    return data


def write_input(path, frames, meta=None, extra=None, cutoff=None):
    """
    Create the input JSON file used by the default chemiscope visualizer, and
    save it to the given ``path``.

    :param str path: name of the file to use to save the json data. If it
                     ends with '.gz', a gzip compressed file will be written
    :param list frames: list of atomic structures. For now, only `ase.Atoms`_
                        objects are supported
    :param dict meta: optional metadata of the dataset
    :param dict extra: optional dictionary of additional properties
    :param float cutoff: optional. If present, will be used to generate
                         atom-centered environments

    This function uses :py:func:`create_input` to generate the input data, see
    the documentation of this function for more information.
    """

    if not (path.endswith(".json") or path.endswith(".json.gz")):
        raise Exception("path should end with .json or .json.gz")

    data = create_input(frames, meta, extra, cutoff)

    if "name" not in data["meta"] or data["meta"]["name"] == "<unknown>":
        data["meta"]["name"] = os.path.basename(path).split(".")[0]

    if path.endswith(".gz"):
        with gzip.open(path, "w", 9) as file:
            file.write(json.dumps(data).encode("utf8"))
    else:
        with open(path, "w") as file:
            json.dump(data, file)
