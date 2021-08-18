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


def _expand_properties(properties, n_structures, n_atoms):
    """
    Convert a shortened entries of properties into the expanded form.
    Entries in already expanded form are not changed.

    :param dict properties: properties to handle
    :param int n_structures: number of structures in the dataset
    :param int n_atoms: total number of atoms in the whole dataset

    For example this property dict:
    .. code-block:: python

        properties = {
            'aple': {
                'target': 'atom',
                'values': np.zeros((300, 4)),
                'unit': 'random / fs',
            }
            'orange' : np.zeros((100, 42)),
            'banana' : np.zeros((300, 17)),
        }

    will be converted to
    .. code-block:: python

        properties = {
            'aple': {
                'target': 'atom',
                'values': np.zeros((300, 4)),
                'unit': 'random / fs',
            }
            'orange': {
                'target': 'structure'
                'values': np.zeros((100, 42)),
            }
            'banana': {
                'target': 'atom',
                'values': np.zeros((300, 17)),
           }
        }

    assuming that number of structures in the dataset is 100 and
    total number of atoms in the dataset is 300.
    """

    for key, value in properties.items():
        if not isinstance(value, dict):
            if (not isinstance(value, list)) and (not isinstance(value, np.ndarray)):
                raise ValueError(
                    "Property values should be either list or numpy array, "
                    + f"got {type(value)} instead"
                )
            if n_structures == n_atoms:
                raise ValueError(
                    "Unable to guess the property target when the number of "
                    + "structures is equal to the number of atoms. We have "
                    + f"n_structures = n_atoms = {n_atoms} for the '{key}' property"
                )
            if (len(value) != n_structures) and (len(value) != n_atoms):
                raise ValueError(
                    "The length of property values is different from the "
                    + "number of structures and the number of atoms, we can not "
                    + f"guess the target. Got n_atoms = {n_atoms}, n_structures = "
                    + f"{n_structures}, the length of property values is "
                    + f"{len(value)}, for the '{key}' property"
                )
            property = {"values": value}
            if len(value) == n_structures:
                property["target"] = "structure"
            if len(value) == n_atoms:
                property["target"] = "atom"
            properties[key] = property
    return properties


def create_input(frames, meta=None, properties=None, cutoff=None, composition=False):
    """
    Create a dictionary that can be saved to JSON using the format used by
    the default chemiscope visualizer.

    :param list frames: list of atomic structures. For now, only `ase.Atoms`_
                        objects are supported
    :param dict meta: optional metadata of the dataset, see below
    :param dict properties: optional dictionary of additional properties, see below
    :param float cutoff: optional. If present, will be used to generate
                         atom-centered environments
    :param bool composition: optional. False by default. If True, will add to
                            the structure and atom properties information
                            about chemical composition

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
    have the atomic numbers as a property, you should add it to ``properties``
    manually.

    Additional properties can be added with the ``properties`` parameter. This
    parameter should be a dictionary containing one entry for each property.
    Each entry contains a ``target`` attribute (``'atom'`` or ``'structure'``)
    and a set of values. ``values`` can be a Python list of float or string; a
    1D numpy array of numeric values; or a 2D numpy array of numeric values. In
    the later case, multiple properties will be generated along the second axis.
    For example, passing

    .. code-block:: python

        properties = {
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

    It is also possible to pass shortened representation of the properties, for
    instance:

    .. code-block:: python

        properties = {
            'cheese':  np.zeros((300, 4)),
            }
        }

    In this case, the type of property (structure or atom) would be deduced
    by comparing the numbers atoms and structures in the dataset to the
    length of provided list/np.ndarray.

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
                warnings.warn(f"ignoring unexpected metadata: {key}")

    if "name" not in data["meta"] or not data["meta"]["name"]:
        data["meta"]["name"] = "<unknown>"

    data["structures"] = frames_to_json(frames)
    n_structures = len(data["structures"])
    n_atoms = sum(s["size"] for s in data["structures"])

    data["properties"] = {}
    if properties is not None:
        properties = _expand_properties(properties, n_structures, n_atoms)
        for name, value in properties.items():
            _validate_property(name, value)
            data["properties"].update(_linearize(name, value, n_structures, n_atoms))

    # Read properties coming from the frames
    for name, value in atom_properties(frames, composition).items():
        _validate_property(name, value)
        for name, value in _linearize(name, value, n_structures, n_atoms).items():
            if name in data["properties"]:
                warnings.warn(
                    f"ignoring the '{name}' atom property coming from the "
                    "structures since it is already part of the properties"
                )
                continue

            data["properties"][name] = value

    for name, value in structure_properties(frames, composition).items():
        _validate_property(name, value)
        for name, value in _linearize(name, value, n_structures, n_atoms).items():
            if name in data["properties"]:
                warnings.warn(
                    f"ignoring the '{name}' structure property coming from the "
                    "structures since it is already part of the properties"
                )
                continue

            data["properties"][name] = value

    if cutoff is not None:
        data["environments"] = _generate_environments(frames, cutoff)

    return data


def write_input(
    path, frames, meta=None, properties=None, cutoff=None, composition=False
):
    """
    Create the input JSON file used by the default chemiscope visualizer, and
    save it to the given ``path``.

    :param str path: name of the file to use to save the json data. If it ends
                     with '.gz', a gzip compressed file will be written
    :param list frames: list of atomic structures. For now, only `ase.Atoms`_
                        objects are supported
    :param dict meta: optional metadata of the dataset
    :param dict properties: optional dictionary of additional properties
    :param float cutoff: optional. If present, will be used to generate
                         atom-centered environments
    :param bool composition: optional. False by default. If True, will add to
                                the structure and atom properties information
                                about chemical composition

    This function uses :py:func:`create_input` to generate the input data, see
    the documentation of this function for more information.

    Here is a quick example of generating a chemiscope input reading the
    structures from a file that `ase <ase-io_>`_ can read, and performing PCA
    using `sklearn`_ on a descriptor computed with another package.

    .. code-block:: python

        import ase
        from ase import io
        import numpy as np
        import sklearn
        from sklearn import decomposition
        from chemiscope import write_input

        frames = ase.io.read('trajectory.xyz', ':')

        # example property 1: list containing the energy of each structure,
        # from calculations performed beforehand
        energies = [ ... ]


        # example property 2: PCA projection computed using sklearn.
        # X contains a multi-dimensional descriptor of the structure
        X = np.array( ... )
        pca = sklearn.decomposition.PCA(n_components=3).fit_transform(X)

        properties = {
            "PCA": {
                "target": "atom",
                "values": pca,
            },
            "energies": {
                "target": "structure",
                "values": energies,
                "units": "kcal/mol",
            },
        }

        write_input("chemiscope.json.gz", frames=frames, properties=properties)

    .. _ase-io: https://wiki.fysik.dtu.dk/ase/ase/io/io.html
    .. _sklearn: https://scikit-learn.org/
    """

    if not (path.endswith(".json") or path.endswith(".json.gz")):
        raise Exception("path should end with .json or .json.gz")

    data = create_input(frames, meta, properties, cutoff, composition)

    if "name" not in data["meta"] or data["meta"]["name"] == "<unknown>":
        data["meta"]["name"] = os.path.basename(path).split(".")[0]

    if path.endswith(".gz"):
        with gzip.open(path, "w", 9) as file:
            file.write(json.dumps(data).encode("utf8"))
    else:
        with open(path, "w") as file:
            json.dump(data, file)


def _typetransform(data, name):
    """
    Transform the given data to either a list of string or a list of floats.

    :param data: list of unknown type to be converted
    :param name: name of the property related to this data, to be used in
                 error messages
    """
    assert isinstance(data, list) and len(data) > 0
    if isinstance(data[0], str):
        return list(map(str, data))
    elif isinstance(data[0], bytes):
        return list(map(lambda u: u.decode("utf8"), data))
    else:
        try:
            return [float(value) for value in data]
        except Exception:
            raise Exception(
                f"unsupported type in property '{name}' values: "
                + "should be string or number"
            )


def _linearize(name, property, n_structures, n_atoms):
    """
    Transform a property dict (containing "value", "target", "units",
    "description") with potential multi-dimensional "values" key to data that
    chemiscope can load.

    Multi-dimensional "value" generate multiple properties named "XXX [1]", "XXX
    [2]", "XXX [3]", etc. Data in "values" are converted to either string or
    float, to ensure it is compatible with JSON.

    :param name: name of the property related to this data, to be used in error
                 messages
    :param property: dictionary containing the property data and metadata
    :param n_structures: total number of structures, to validate the array sizes
    :param n_atoms: total number of atoms, to validate the array sizes
    """
    data = {}
    if isinstance(property["values"], list):
        data[name] = {
            "target": property["target"],
            "values": _typetransform(property["values"], name),
        }
    elif isinstance(property["values"], np.ndarray):
        if len(property["values"].shape) == 1:
            data[name] = {
                "target": property["target"],
                "values": _typetransform(list(property["values"]), name),
            }
        elif len(property["values"].shape) == 2:
            if property["values"].shape[1] == 1:
                data[name] = {
                    "target": property["target"],
                    "values": _typetransform(list(property["values"]), name),
                }
            else:
                for i in range(property["values"].shape[1]):
                    data[f"{name}[{i + 1}]"] = {
                        "target": property["target"],
                        "values": _typetransform(list(property["values"][:, i]), name),
                    }
        else:
            raise Exception("unsupported ndarray property")
    else:
        raise Exception(
            f"unknown type ({type(property['values'])}) for property '{name}'"
        )

    # get property metadata
    if "units" in property:
        for item in data.values():
            item["units"] = str(property["units"])

    if "description" in property:
        for i, item in enumerate(data.values()):
            # add [component XX] to the description if values was a ndarray
            extra = f" [component {i + 1}]" if len(data) > 1 else ""
            item["description"] = str(property["description"]) + extra

    # Validate the properties size
    for prop in data.values():
        if prop["target"] == "atom" and len(prop["values"]) != n_atoms:
            raise Exception(
                f"wrong size for the property '{name}' with target=='atom': "
                + f"expected {n_atoms} values, got {len(prop['values'])}"
            )

        if prop["target"] == "structure" and len(prop["values"]) != n_structures:
            raise Exception(
                f"wrong size for the property '{name}' with target=='structure': "
                + f"expected {n_structures} values, got {len(prop['values'])}"
            )

    return data


def _generate_environments(frames, cutoff):
    if not isinstance(cutoff, float):
        raise Exception(
            f"cutoff must be a float, got '{cutoff}' of type {type(cutoff)}"
        )

    environments = []
    for frame_id, frame in enumerate(frames):
        for center in range(len(frame)):
            environments.append(
                {"structure": frame_id, "center": center, "cutoff": cutoff}
            )
    return environments


def _validate_property(name, property):
    if name == "":
        raise Exception("the name of a property can not be the empty string")
    elif not isinstance(name, str):
        raise Exception(
            "the name of a property name must be a string, "
            + f"got '{name}' of type {type(name)}"
        )

    if "target" not in property:
        raise Exception(f"missing 'target' for the '{name}' property")
    elif property["target"] not in ["atom", "structure"]:
        raise Exception(
            f"the target must be 'atom' or 'structure' for the '{name}' property"
        )

    if "values" not in property:
        raise Exception(f"missing 'values' for the '{name}' property")

    for key in property.keys():
        if key not in ["target", "values", "description", "units"]:
            warnings.warn(f"ignoring unexpected property key: {key}")
