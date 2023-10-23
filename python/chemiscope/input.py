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


def create_input(
    frames=None,
    meta=None,
    properties=None,
    environments=None,
    settings=None,
    shapes=None,
    parameters=None,
):
    """
    Create a dictionary that can be saved to JSON using the format used by the default
    chemiscope visualizer.

    :param list frames: list of atomic structures. For now, only `ase.Atoms`_ objects
        are supported

    :param dict meta: optional metadata of the dataset, see below

    :param dict properties: optional dictionary of properties, see below

    :param list environments: optional list of ``(structure id, atom id, cutoff)``
        specifying which atoms have properties attached and how far out atom-centered
        environments should be drawn by default. Functions like
        :py:func:`all_atomic_environments` or :py:func:`librascal_atomic_environments`
        can be used to generate the list of environments in simple cases.

    :param dict shapes: optional dictionary of shapes to have available for display,
        see below.

    :param dict settings: optional dictionary of settings to use when displaying the
        data. Possible entries for the ``settings`` dictionary are documented in the
        chemiscope input file reference.

    :param dict parameters: optional dictionary of parameters for multidimensional
        properties, see below


    Dataset metadata
    ----------------

    The dataset metadata should be given in the ``meta`` dictionary, the possible keys
    are:

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

    Dataset properties
    ------------------

    Properties can be added with the ``properties`` parameter. This parameter should be
    a dictionary containing one entry for each property. Properties can be extracted
    from structures with :py:func:`extract_properties` or
    :py:func:`composition_properties`, or manually defined by the user.

    Each entry in the ``properties`` dictionary contains a ``target`` attribute
    (``'atom'`` or ``'structure'``) and a set of values. ``values`` can be a Python list
    of float or string; a 1D numpy array of numeric values; or a 2D numpy array of
    numeric values. In the later case, multiple properties will be generated along the
    second axis. For example, passing

    .. code-block:: python

        properties = {
            'cheese': {
                'target': 'atom',
                'values': np.zeros((300, 4)),
                # optional: property unit
                'units': 'random / fs',
                # optional: property description
                'description': 'a random property for example',
            }
        }

    will generate four properties named ``cheese[1]``, ``cheese[2]``, ``cheese[3]``, and
    ``cheese[4]``, each containing 300 values.

    It is also possible to pass shortened representation of the properties, for
    instance:

    .. code-block:: python

        properties = {
            'cheese':  np.zeros((300, 4)),
        }

    In this case, the type of property (structure or atom) would be deduced by comparing
    the numbers atoms and structures in the dataset to the length of provided
    list/np.ndarray.

    Multi-dimensional properties
    ----------------------------

    One can give 2D properties to be displayed as curves in the info panel by setting a
    ``parameters`` in the property, and giving the corresponding ``parameters`` values
    to this function. The previous example becomes:

    .. code-block:: python

        properties = {
            'cheese': {
                'target': 'atom',
                'values': np.zeros((300, 4)),
                # optional: property unit
                'units': 'random / fs',
                # optional: property description
                'description': 'a random property for example',
                'parameters': ['origin'],
            }
        }

    This input describes a 2D property ``cheese`` with 300 samples and 4 values taken by
    the ``origin`` parameter. We also need to provide the ``parameters`` values to this
    function:

    .. code-block:: python

        parameters = {
            'origin' : {
                # an array of numbers containing the values of the parameter
                # the size should correspond to the second dimension
                # of the corresponding multidimensional property
                'values': [0, 1, 2, 3],
                # optional free-form description of the parameter as a string
                'name': 'a short description of this parameter',
                # optional units of the values in the values array
                'units': 'eV',
            }
        }

    Custom shapes
    -------------

    The ``shapes`` option should have the format ``{"<name>": shape_definition }``,
    where each shape is defined as a dictionary containing the kind of shape, and its
    parameters

    .. code-block:: python

        shapes = {
            "shape name": {
                "kind" : "sphere",
                "parameters" : shape_parameters
            }
        }

    Each parameters block defines `global`, `structure` and `atom` - level parameters.

    .. code-block:: python

        parameters = {
            "global" : global_parameters,
            "structure" : [ structure_1, structure_2, .... ],
            "atom"" : [ atom_1, atom_2, .... ]
        }

    Each of these can contain some or all of the parameters associated with each shape,
    and the parameters for each shape are obtained by combining the parameters from the
    most general to the most specific, i.e., if there is a duplicate key in the
    `global` and `atom` fields, the value within the `atom` field will supersede the
    `global` field for that atom. The parameters for atom `k` that is part of structure
    `j` are obtained as

    .. code-block:: python

        global_parameters.update(structure_j).update(atom_k)

    If given, the `structure` parameters list should contain one entry per structure,
    and the `atom` parameters list should be a flat list corresponding to the atoms of
    each consecutive structure. All shapes accept a few general parameters, and some
    specific ones

    .. code-block:: python

        # general parameters
        {
            # centering (defaults to origin for structure, atom position for atom)
            "position" : [float, float, float],
            # scaling of the size of the shape
            "scale" : float,
            # optional, given as quaternion in (x, y, z, w) format
            "orientation" [float, float, float, float],
            "color" : string | hex code # e.g. 0xFF0000
        }

        # "kind" : "sphere"
        {
            "radius": float,
        }

        # "kind" : "ellipsoid"
        {
            "semiaxes": [float, float, float],
        }

        # "kind" : "arrow"
        {   # "orientation" is redundant and hence ignored
            "vector" : [float, float, float],  # orientation and shape of the arrow
            "baseRadius" : float,
            "headRadius" : float,
            # the tip of the arrow is at the end of the segment.
            # It will extend past the base point if the arrow is not long enough
            "headLength" : float,
        }

        # "kind" : "custom"
        {
            "vertices": [ # list of vertices
                [float, float, float],
                ...
            ],
            # mesh triangulation (optional); computed via convex triangulation
            # where omitted
            "simplices": [
                [int, int, int],  # indices refer to the list of vertices
                ...
            ],
        }


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

    if shapes is not None:
        data["shapes"] = _validate_shapes(data["structures"], shapes)

    data["properties"] = {}
    if properties is not None:
        properties = _expand_properties(properties, n_structures, n_atoms)
        for name, value in properties.items():
            multidimensional = False
            if "parameters" in value:
                multidimensional = True
            data["properties"].update(
                _linearize(name, value, n_structures, n_atoms, multidimensional)
            )

    if parameters is not None:
        if not isinstance(parameters, dict):
            raise ValueError(
                f"expecting parameters to be a of type 'dict' not '{type(parameters)}'"
            )

        data["parameters"] = {}
        for key in parameters:
            param = {}
            if not isinstance(parameters[key], dict):
                raise ValueError(
                    f"expecting parameter {key} to be of type 'dict' not "
                    f"'{type(parameters[key])}'"
                )

            if not isinstance(np.asarray(parameters[key]["values"]), np.ndarray):
                raise ValueError(
                    f"expecting the 'values' of parameter {key} to be of type "
                    f"'numpy.ndarray' not '{type(parameters[key]['values'])}'"
                )
            if len(parameters[key]["values"].shape) != 1:
                raise ValueError(
                    f"the 'values' of parameter {key} should be a 1D numpy array"
                )
            param["values"] = _typetransform(
                list(parameters[key]["values"]), f"parameter {key}"
            )

            if "name" in parameters[key]:
                if not isinstance(parameters[key]["name"], str):
                    raise ValueError(
                        f"optional 'name': {parameters[key]['name']} of "
                        f"parameter {key} should a string"
                    )
                param["name"] = parameters[key]["name"]

            if "units" in parameters[key]:
                if not isinstance(parameters[key]["units"], str):
                    raise ValueError(
                        f"optional 'units': {parameters[key]['units']} of "
                        f"parameter {key} should a string"
                    )
                param["units"] = parameters[key]["units"]
            data["parameters"][key] = param

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

    return data


def write_input(
    path,
    frames,
    meta=None,
    properties=None,
    environments=None,
    shapes=None,
    settings=None,
    parameters=None,
):
    """
    Create the input JSON file used by the default chemiscope visualizer, and save it to
    the given ``path``.

    :param str path: name of the file to use to save the json data. If it ends with
        '.gz', a gzip compressed file will be written

    :param list frames: list of atomic structures. For now, only `ase.Atoms`_ objects
        are supported

    :param dict meta: optional metadata of the dataset

    :param dict properties: optional dictionary of additional properties

    :param list environments: optional list of ``(structure id, atom id, cutoff)``
        specifying which atoms have properties attached and how far out atom-centered
        environments should be drawn by default.

    :param dict shapes: optional dictionary of shapes to have available for display.
        See :py:func:`create_input` for more information on how to define shapes.

    :param dict settings: optional dictionary of settings to use when displaying the
        data. Possible entries for the ``settings`` dictionary are documented in the
        chemiscope input file reference.

    :param dict parameters: optional dictionary of parameters of multidimensional
        properties

    This function uses :py:func:`create_input` to generate the input data, see the
    documentation of this function for more information.

    Here is a quick example of generating a chemiscope input reading the structures from
    a file that `ase <ase-io_>`_ can read, and performing PCA using `sklearn`_ on a
    descriptor computed with another package.

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

        # additional multidimensional properties to be plotted
        dos = np.loadtxt(...) # load the 2D data
        dos_energy_grid = np.loadtxt(...)
        multidimensional_properties = {
            "DOS": {
                "target": "structure",
                "values": dos,
                "parameters": ["energy"],
            }
        }

        multidimensional_parameters = {
            "energy": {
                "values": dos_energy_grid,
                "units": "eV",
            }
        }

        # merge all properties together
        properties.extend(frame_properties)
        properties.extend(multidimensional_properties)

        chemiscope.write_input(
            path="chemiscope.json.gz",
            frames=frames,
            properties=properties,
            # This is required to display properties with `target: "atom"`
            environments=chemiscope.all_atomic_environments(frames),
            # this is necessary to plot the multidimensional data
            parameters=multidimensional_parameters,
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
        shapes=shapes,
        settings=settings,
        parameters=parameters,
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
                'units': 'random / fs',
            }
            'orange' : np.zeros((100, 42)),
            'banana' : np.zeros((300, 17)),
        }

    will be converted to
    .. code-block:: python

        properties = {
            'apple': {
                'target': 'atom',
                'values': np.zeros((300, 4)),
                'units': 'random / fs',
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


def _linearize(name, property, n_structures, n_centers, multidimensional=False):
    """
    Transform a single property dict (containing "value", "target", "units",
    "description", "parameters") with potential multi-dimensional "values" key
    to data that chemiscope can load.

    Multi-dimensional "value" generate multiple properties named "XXX [1]", "XXX
    [2]", "XXX [3]", etc. Data in "values" are converted to either string or
    float, to ensure it is compatible with JSON.

    :param name: name of the property related to this data, to be used in error
                 messages
    :param property: dictionary containing the property data and metadata
    :param n_structures: total number of structures, to validate the array sizes
    :param n_centers: total number of atoms, to validate the array sizes
    :param param: describes if the property is multisimensional to be plotted
    """
    _validate_property(name, property)

    data = {}
    values = property["values"]
    if not multidimensional:
        if isinstance(values, list) and isinstance(values[0], list):
            # convert to ndarray so we can use the parser below
            values = np.array(values)
        if isinstance(values, list):
            # deal with simple list
            data[name] = {
                "target": property["target"],
                "values": _typetransform(values, name),
            }
        elif isinstance(values, np.ndarray):
            if len(values.shape) == 1:
                data[name] = {
                    "target": property["target"],
                    "values": _typetransform(list(values), name),
                }
            elif len(values.shape) == 2:
                if values.shape[1] == 1:
                    data[name] = {
                        "target": property["target"],
                        "values": _typetransform(list(values), name),
                    }
                else:
                    for i in range(values.shape[1]):
                        data[f"{name}[{i + 1}]"] = {
                            "target": property["target"],
                            "values": _typetransform(list(values[:, i]), name),
                        }
            else:
                raise Exception("unsupported ndarray property")
        else:
            raise Exception(
                f"unknown type ({type(property['values'])}) for property '{name}'"
            )
    else:
        assert isinstance(property["parameters"][0], str)

        ndvalues = []
        for i in range(len(values)):
            ndvalues += [_typetransform(list(values[i]), name)]
        data[name] = {
            "target": property["target"],
            "values": ndvalues,
            "parameters": property["parameters"],
        }

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

    if "parameters" in property:
        if len(property["parameters"]) != 1 or not isinstance(
            property["parameters"][0], str
        ):
            raise Exception(
                f"the parameter of property {name} should be a list containing "
                "a single string"
            )

    for key in property.keys():
        if key not in ["target", "values", "description", "units", "parameters"]:
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


def _validate_shapes(structures, shapes):
    if not isinstance(shapes, dict):
        raise TypeError(f"`shapes` must be a dictionary, got {type(shapes)} instead")

    # validate type and number of element for each entries in the shapes
    for key, shapes_for_key in shapes.items():
        if not isinstance(key, str):
            raise TypeError(
                f"the `shapes` dictionary keys must be strings, got {type(key)}"
            )

        if not isinstance(shapes_for_key, dict):
            raise TypeError(
                "Each entry in `shapes` must be a dictionary, "
                f"got {type(shapes_for_key)} instead for '{key}'"
            )

        for shape_key in shapes_for_key:
            if shape_key not in ["kind", "parameters"]:
                raise ValueError(
                    f"Invalid entry `{shape_key}` in the specification "
                    f"for shape `{key}`"
                )

        base_shape = shapes_for_key["parameters"].get("global", {})
        structure_parameters = shapes_for_key["parameters"].get("structure", None)
        atom_parameters = shapes_for_key["parameters"].get("atom", None)
        atom_counter = 0

        for structure_i in range(len(structures)):
            if (
                structure_parameters is not None
                and len(structure_parameters) <= structure_i
            ):
                raise TypeError(
                    f"structure_parameters must be a list with {(len(structures))} "
                    f"elements, got {len(structure_parameters)}"
                )
            for _ in range(structures[structure_i]["size"]):
                if atom_parameters is not None and len(atom_parameters) <= atom_counter:
                    raise TypeError(
                        f"atom_parameters must be a list coinciding to the atomic "
                        f"environments, got {len(atom_parameters)} elements"
                    )

                shape = {
                    "kind": shapes_for_key["kind"],
                    "parameters": {
                        "global": base_shape,
                        "structure": structure_parameters[structure_i]
                        if structure_parameters is not None
                        else {},
                        "atom": atom_parameters[atom_counter]
                        if atom_parameters is not None
                        else {},
                    },
                }
                _check_valid_shape(shape)
                atom_counter += 1

    for key, shape in shapes.items():
        if (
            shape["kind"] == "custom"
            and "vertices" in shape["parameters"]
            and "simplices" not in shape["parameters"]
        ):
            try:
                import scipy.spatial
            except ImportError as e:
                raise RuntimeError(
                    "Missing simplices in custom shape, and scipy is not " "installed"
                ) from e

            convex_hull = scipy.spatial.ConvexHull(shape["parameters"]["vertices"])
            shape["parameters"]["simplices"] = [
                s.tolist() for s in convex_hull.simplices
            ]
    return shapes


def _check_valid_shape(shape):
    if not isinstance(shape, dict):
        raise TypeError(
            f"individual shapes must be dictionaries, got {type(shape)} instead"
        )

    always_okay = ["orientation", "scale", "position", "color"]
    parameters = {}
    if "parameters" in shape:
        parameters.update(shape["parameters"]["global"])
    if "structure" in shape["parameters"]:
        parameters.update(shape["parameters"]["structure"])
    if "atom" in shape["parameters"]:
        parameters.update(shape["parameters"]["atom"])

    if len(parameters) == 0:
        raise ValueError(f"no parameters provided for {shape['kind']} shape")
    if shape["kind"] == "sphere":
        for parameter in parameters:
            if parameter not in ["radius", *always_okay]:
                raise ValueError(
                    f"unknown shape parameter '{parameter}' for 'sphere' shape kind"
                )

        if not isinstance(parameters["radius"], float):
            raise TypeError(
                "sphere shape 'radius' must be a float, got "
                f"{type(parameters['radius'])}"
            )

    elif shape["kind"] == "ellipsoid":
        for parameter in parameters.keys():
            if parameter not in ["semiaxes", *always_okay]:
                raise ValueError(
                    f"unknown shape parameter '{parameter}' for 'ellipsoid' shape kind"
                )

        semiaxes_array = np.asarray(parameters["semiaxes"]).astype(
            np.float64, casting="safe", subok=False, copy=False
        )

        if not semiaxes_array.shape == (3,):
            raise ValueError(
                "'semiaxes' must be an array with 3 values for 'ellipsoid' shape kind"
            )

    elif shape["kind"] == "custom":
        for parameter in parameters.keys():
            if parameter not in ["vertices", "simplices", *always_okay]:
                raise ValueError(
                    f"unknown shape parameter '{parameter}' for 'custom' shape kind"
                )

        vertices_array = np.asarray(parameters["vertices"]).astype(
            np.float64, casting="safe", subok=False, copy=False
        )

        if len(vertices_array.shape) != 2 or vertices_array.shape[1] != 3:
            raise ValueError(
                "'vertices' must be an Nx3 array values for 'custom' shape kind"
            )

        if "simplices" in parameter:
            simplices_array = np.asarray(parameters["vertices"]).astype(
                int, casting="safe", subok=False, copy=False
            )

            if len(simplices_array.shape) != 2 or simplices_array.shape[1] != 3:
                raise ValueError(
                    "'simplices' must be an Nx3 array values for 'custom' shape kind"
                )
    elif shape["kind"] == "arrow":
        if not isinstance(parameters["baseRadius"], float):
            raise TypeError(
                "sphere shape 'baseRadius' must be a float, "
                f"got {type(parameters['baseRadius'])}"
            )
        if not isinstance(parameters["headRadius"], float):
            raise TypeError(
                "sphere shape 'headRadius' must be a float, "
                f"got {type(parameters['headRadius'])}"
            )
        if not isinstance(parameters["headLength"], float):
            raise TypeError(
                "sphere shape 'headLength' must be a float, "
                f"got {type(parameters['headLength'])}"
            )

        vector_array = np.asarray(parameters["vector"]).astype(
            np.float64, casting="safe", subok=False, copy=False
        )

        if not vector_array.shape == (3,):
            raise ValueError(
                "'vector' must be an array with 3 values for 'arrow' shape kind"
            )
    else:
        raise ValueError(f"unknown shape kind '{shape['kind']}'")

    if "orientation" in parameters:
        orientation_array = np.asarray(parameters["orientation"]).astype(
            np.float64, casting="safe", subok=False, copy=False
        )

        if not orientation_array.shape == (4,):
            raise ValueError(
                "semiaxes must be an array with 4 values for 'ellipsoid' shape kind"
            )
