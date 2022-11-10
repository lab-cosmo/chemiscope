# -*- coding: utf-8 -*-
"""
Generate JSON input files for the default chemiscope visualizer.
"""
import gzip
import json
import os
import warnings

import numpy as np

from .structures import (
    _list_atom_properties,
    _list_structure_properties,
    frames_to_json,
)
from .version import __version__

if tuple(map(int, __version__.split("."))) >= (0, 6, 0):
    raise Exception(
        "this is a reminder to remove the warning about automatic properties extraction"
    )


def create_input(
    frames=None,
    meta=None,
    properties=None,
    environments=None,
    settings=None,
):
    """
    Create a dictionary that can be saved to JSON using the format used by
    the default chemiscope visualizer.

    :param list frames: list of atomic structures. For now, only `ase.Atoms`_
                        objects are supported
    :param dict meta: optional metadata of the dataset, see below
    :param dict properties: optional dictionary of properties, see below
    :param list environments: optional list of (structure id, atom id, cutoff)
        specifying which atoms have properties attached and how far out
        atom-centered environments should be drawn by default. Functions like
        :py:func:`all_atomic_environments` or :py:func:`librascal_atomic_environments`
        can be used to generate the list of environments in simple cases.
    :param dict settings: optional dictionary of settings to use when displaying
        the data. Possible entries for the ``settings`` dictionary are documented
        in the chemiscope input file reference.

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

    Properties can be added with the ``properties`` parameter. This parameter
    should be a dictionary containing one entry for each property. Properties
    can be extracted from structures with :py:func:`extract_properties` or
    :py:func:`composition_properties`, or manually defined by the user.

    Each entry in the ``properties`` dictionary contains a ``target`` attribute
    (``'atom'`` or ``'structure'``) and a set of values. ``values`` can be a
    Python list of float or string; a 1D numpy array of numeric values; or a 2D
    numpy array of numeric values. In the later case, multiple properties will
    be generated along the second axis. For example, passing

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

    data = {
        "meta": _normalize_metadata(meta if meta is not None else {}),
    }

    if settings is not None:
        # dump/load as json to catch possible json incompatibility in settings
        # early
        if not isinstance(settings, dict):
            raise ValueError(
                f"expected 'settings' to be a dict, got {type(settings)} instead"
            )

        data["settings"] = json.loads(json.dumps(settings))

    data["structures"] = []
    n_structures = None

    if frames is not None:
        data["structures"] = frames_to_json(frames)
        n_structures = len(data["structures"])
        n_atoms = sum(s["size"] for s in data["structures"])
    else:
        n_atoms = 0

        # if frames are not given, we create a dataset with only properties.
        # In that case, all properties should be structure properties
        for name, value in properties.items():
            if not isinstance(value, dict):
                if n_structures is None:
                    n_structures = len(value)
                else:
                    if len(value) != n_structures:
                        raise ValueError(
                            f"wrong size for property '{name}': expected "
                            f"{n_structures} elements, but got an array with "
                            f"{len(value)} entries"
                        )
            else:
                if value["target"] != "structure":
                    raise ValueError(
                        f"property '{name}' has a non-structure target, "
                        "which is not allowed if frames are not provided"
                    )

                n_structures = len(value["values"])

    if environments is not None:
        if "structures" not in data:
            raise ValueError("can not have environments without structures")

        data["environments"] = _normalize_environments(environments, data["structures"])
        n_atoms = len(data["environments"])

    data["properties"] = {}
    if properties is not None:
        properties = _expand_properties(properties, n_structures, n_atoms)
        for name, value in properties.items():
            data["properties"].update(_linearize(name, value, n_structures, n_atoms))

    # Check to tell the user they might have forgotten some properties coming
    # from the frames (that chemiscope used to automatically extract). This code
    # should be removed in version 0.6 of chemiscope.
    if frames is not None:
        found_one_from_frame = False
        atom_properties = _list_atom_properties(frames)
        for name in atom_properties:
            if name in data["properties"]:
                found_one_from_frame = True

        structure_properties = _list_structure_properties(frames)
        for name in structure_properties:
            if name in data["properties"]:
                found_one_from_frame = True

        if not found_one_from_frame:
            properties_list = ""

            if len(structure_properties) != 0:
                properties_list += "[" + ", ".join(structure_properties) + "]"

            if len(atom_properties) != 0:
                if len(properties_list) != 0:
                    properties_list += " and "
                properties_list += "[" + ", ".join(atom_properties) + "]"

            if len(properties_list) != 0:
                warnings.warn(
                    "chemiscope behavior changed to no longer include properties "
                    "from the structure objects. Use `chemiscope.extract_properties` "
                    f"to also visualize these properties ({properties_list})"
                )

    return data


def write_input(
    path,
    frames,
    meta=None,
    properties=None,
    environments=None,
    settings=None,
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
    :param list environments: optional list of (structure id, atom id, cutoff)
        specifying which atoms have properties attached and how far out
        atom-centered environments should be drawn by default.
    :param dict settings: optional dictionary of settings to use when displaying
        the data. Possible entries for the ``settings`` dictionary are documented
        in the chemiscope input file reference.

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
        import chemiscope

        frames = ase.io.read('trajectory.xyz', ':')

        # example property 1: list containing the energy of each structure,
        # from calculations performed beforehand
        energies = np.loadtxt('energies.txt')

        # example property 2: PCA projection computed using sklearn.
        # X contains a multi-dimensional descriptor of the structure
        X = np.array( ... )
        pca = sklearn.decomposition.PCA(n_components=3).fit_transform(X)

        properties = {
            "PCA": {
                "target": "atom",
                "values": pca,
                "description": "PCA of per-atom representation of the structures",
            },
            "energies": {
                "target": "structure",
                "values": energies,
                "units": "kcal/mol",
            },
        }

        # additional properties coming from the trajectory
        frame_properties = chemiscope.extract_properties(
            frames,
            only=["temperature", "classification"]
        )

        # merge all properties together
        properties.extend(frame_properties)

        chemiscope.write_input(
            path="chemiscope.json.gz",
            frames=frames,
            properties=properties,
            # This is required to display properties with `target: "atom"`
            environments=chemiscope.all_atomic_environments(frames),
        )

    .. _ase-io: https://wiki.fysik.dtu.dk/ase/ase/io/io.html
    .. _sklearn: https://scikit-learn.org/
    """

    if not (path.endswith(".json") or path.endswith(".json.gz")):
        raise Exception("path should end with .json or .json.gz")

    data = create_input(
        frames=frames,
        meta=meta,
        properties=properties,
        environments=environments,
        settings=settings,
    )

    if "name" not in data["meta"] or data["meta"]["name"] == "<unknown>":
        data["meta"]["name"] = os.path.basename(path).split(".")[0]

    if path.endswith(".gz"):
        with gzip.open(path, "w", 9) as file:
            file.write(json.dumps(data).encode("utf8"))
    else:
        with open(path, "w") as file:
            json.dump(data, file, indent=2)


def _normalize_environments(environments, structures):
    cleaned = []
    for environment in environments:
        if len(environment) != 3:
            raise ValueError(
                f"expected environments to contain three values, got {environment}"
            )

        structure, center, cutoff = environment
        structure = int(structure)
        center = int(center)
        cutoff = float(cutoff)

        if structure >= len(structures):
            raise ValueError(
                f"invalid structure index in environments: got {structure}, "
                f"but we have {len(structures)} structures"
            )

        if center >= structures[structure]["size"]:
            raise ValueError(
                f"invalid center index in environments: got {center} in structure "
                f"which only contains {structures[structure]['size']} atoms"
            )

        if cutoff <= 0:
            raise ValueError("negative cutoff in environments is not valid")

        cleaned.append(
            {
                "structure": structure,
                "center": center,
                "cutoff": float(cutoff),
            }
        )

    return cleaned


def _normalize_metadata(meta):
    cleaned = {}
    if "name" in meta and str(meta["name"]) != "":
        cleaned["name"] = str(meta["name"])
    else:
        cleaned["name"] = "<unknown>"

    if "description" in meta:
        cleaned["description"] = str(meta["description"])

    if "authors" in meta:
        cleaned["authors"] = list(map(str, meta["authors"]))

    if "references" in meta:
        cleaned["references"] = list(map(str, meta["references"]))

    for key in meta.keys():
        if key not in ["name", "description", "authors", "references"]:
            warnings.warn(f"ignoring unexpected metadata: {key}")

    return cleaned


def _expand_properties(short_properties, n_structures, n_atoms):
    """
    Convert a shortened entries of properties into the expanded form.
    Entries in already expanded form are not changed.

    :param dict short_properties: properties to handle
    :param int n_structures: number of structures in the dataset
    :param int n_atoms: total number of atoms in the whole dataset

    For example this property dict:
    .. code-block:: python

        properties = {
            'apple': {
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
    properties = {}
    for key, value in short_properties.items():
        if isinstance(value, dict):
            properties[key] = value
        else:
            if (not isinstance(value, list)) and (not isinstance(value, np.ndarray)):
                raise ValueError(
                    "Property values should be either list or numpy array, "
                    f"got {type(value)} instead"
                )
            if n_structures == n_atoms:
                warnings.warn(
                    f"The target of the property '{key}' is ambiguous because "
                    "there is the same number of atoms and structures. "
                    "We will assume target=structure"
                )

            dict_property = {"values": value}

            # heuristically determines the type of target
            if len(value) == n_structures:
                dict_property["target"] = "structure"
            elif len(value) == n_atoms:
                dict_property["target"] = "atom"
            else:
                raise ValueError(
                    "The length of property values is different from the "
                    "number of structures and the number of atoms, we can not "
                    f"guess the target. Got n_atoms = {n_atoms}, n_structures = "
                    f"{n_structures}, the length of property values is "
                    f"{len(value)}, for the '{key}' property"
                )

            properties[key] = dict_property
    return properties


def _linearize(name, property, n_structures, n_centers):
    """
    Transform a single property dict (containing "value", "target", "units",
    "description") with potential multi-dimensional "values" key to data that
    chemiscope can load.

    Multi-dimensional "value" generate multiple properties named "XXX [1]", "XXX
    [2]", "XXX [3]", etc. Data in "values" are converted to either string or
    float, to ensure it is compatible with JSON.

    :param name: name of the property related to this data, to be used in error
                 messages
    :param property: dictionary containing the property data and metadata
    :param n_structures: total number of structures, to validate the array sizes
    :param n_centers: total number of atoms, to validate the array sizes
    """
    _validate_property(name, property)

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
        if prop["target"] == "atom" and len(prop["values"]) != n_centers:
            raise Exception(
                f"wrong size for the property '{name}' with target=='atom': "
                f"expected {n_centers} values, got {len(prop['values'])}"
            )

        if prop["target"] == "structure" and len(prop["values"]) != n_structures:
            raise Exception(
                f"wrong size for the property '{name}' with target=='structure': "
                f"expected {n_structures} values, got {len(prop['values'])}"
            )

    return data


def _validate_property(name, property):
    if name == "":
        raise Exception("the name of a property can not be the empty string")
    elif not isinstance(name, str):
        raise Exception(
            "the name of a property name must be a string, "
            f"got '{name}' of type {type(name)}"
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
                "should be string or number"
            )
