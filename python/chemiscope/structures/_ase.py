import warnings

import numpy as np

from ._shapes import arrow_from_vector, ellipsoid_from_tensor


try:
    import ase

    HAVE_ASE = True
except ImportError:
    HAVE_ASE = False


def _ase_valid_structures(structures):
    try:
        structures_list = list(structures)
    except TypeError:
        return [], False

    if HAVE_ASE and isinstance(structures_list[0], ase.Atoms):
        for structure in structures_list:
            assert isinstance(structure, ase.Atoms)
        return structures_list, True
    elif HAVE_ASE and isinstance(structures_list[0], ase.Atom):
        # deal with the user passing a single structure
        return [structures], True
    else:
        return structures_list, False


def _ase_normalize_value(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    elif isinstance(value, (np.int32, np.int64)):
        return int(value)
    elif isinstance(value, (np.float32, np.float64)):
        return float(value)
    else:
        return value


def _ase_get_atom_properties(structures):
    IGNORED_ASE_ARRAYS = ["positions", "numbers", "center_atoms_mask"]
    # extract the set of common properties between all structures
    all_properties = {}
    extra = set()

    arrays = structures[0].arrays.copy()
    # workaround for ase >= 3.23
    if structures[0].calc is not None and "forces" in structures[0].calc.results:
        arrays["forces"] = structures[0].calc.results["forces"]

    for name in arrays.keys():
        if name in IGNORED_ASE_ARRAYS:
            continue
        all_properties[name] = [arrays[name]]

    for structure in structures[1:]:
        arrays = structure.arrays.copy()
        # workaround for ase >= 3.23
        if structure.calc is not None and "forces" in structure.calc.results:
            arrays["forces"] = structure.calc.results["forces"]

        for name in arrays.keys():
            if name in IGNORED_ASE_ARRAYS:
                continue

            if name in all_properties:
                all_properties[name].append(arrays[name])
            else:
                extra.add(name)

        for name in list(all_properties.keys()):
            if name not in arrays.keys():
                all_properties.pop(name, None)
                extra.add(name)

    if len(extra) != 0:
        warnings.warn(
            "the following atomic properties are only defined for a subset "
            f"of structures: {list(sorted(extra))}; they will be ignored",
            stacklevel=2,
        )

    return all_properties


def _ase_get_structure_properties(structures):
    # extract the set of common properties between all structures
    all_properties = {}
    extra = set()

    info = structures[0].info.copy()
    # workaround for ase >= 3.23
    if structures[0].calc is not None and "energy" in structures[0].calc.results:
        info["energy"] = structures[0].calc.results["energy"]

    for name in info.keys():
        all_properties[name] = [_ase_normalize_value(info[name])]

    for structure in structures[1:]:
        info = structure.info.copy()
        # workaround for ase >= 3.23
        if structure.calc is not None and "energy" in structure.calc.results:
            info["energy"] = structure.calc.results["energy"]

        for name in info.keys():
            if name in all_properties:
                all_properties[name].append(_ase_normalize_value(info[name]))
            else:
                extra.add(name)

        for name in list(all_properties.keys()):
            if name not in info.keys():
                all_properties.pop(name, None)
                extra.add(name)

    if len(extra) != 0:
        warnings.warn(
            "the following structure properties are only defined for a subset "
            f"of structures: {list(sorted(extra))}; they will be ignored",
            stacklevel=2,
        )

    # ensures that if a property is a mix of strings and numbers, everything
    # is converted to string (as these should be categorical properties)
    for name in all_properties.keys():
        property = all_properties[name]
        if any(isinstance(x, str) for x in property):
            all_properties[name] = [str(x) for x in property]

    return all_properties


def _ase_atom_properties(structures, only=None, atoms_mask=None):
    all_properties = _ase_get_atom_properties(structures)
    if only is None:
        selected = all_properties
    else:
        selected = {}
        for name in only:
            if name in all_properties.keys():
                selected[name] = all_properties[name]

    # create property in the format expected by create_input
    properties = {
        name: {"target": "atom", "values": np.concatenate(value)}
        for name, value in selected.items()
    }

    _remove_invalid_properties(properties, "ASE")

    if atoms_mask is not None:
        # only include values for requested atoms
        for property in properties.values():
            property["values"] = property["values"][atoms_mask]

    return properties


def _ase_structure_properties(structures, only=None):
    all_properties = _ase_get_structure_properties(structures)
    if only is None:
        selected = all_properties
    else:
        selected = {}
        for name in only:
            if name in all_properties.keys():
                selected[name] = all_properties[name]

    # create property in the format expected by create_input
    properties = {
        name: {"target": "structure", "values": value}
        for name, value in selected.items()
    }

    _remove_invalid_properties(properties, "ASE")

    return properties


def _ase_extract_properties(structures, only=None, environments=None):
    """implementation of ``extract_properties`` for ASE"""
    properties = _ase_structure_properties(structures, only)

    if environments is not None:
        atoms_mask = [[False] * len(f) for f in structures]
        for structure, center, _ in environments:
            atoms_mask[structure][center] = True

        atoms_mask = np.concatenate(atoms_mask)
    else:
        atoms_mask = None

    atom_properties = _ase_atom_properties(structures, only, atoms_mask)

    for name, values in atom_properties.items():
        if name in properties:
            warnings.warn(
                f"a property named '{name}' is defined for both atoms and structures, "
                "the atom one will be ignored",
                stacklevel=2,
            )
        else:
            properties[name] = values

    return properties


def _ase_all_atomic_environments(structures, cutoff):
    "Extract all atomic environments out of a set of ASE Atoms objects"
    environments = []
    for structure_i, structure in enumerate(structures):
        for atom_i in range(len(structure)):
            environments.append((structure_i, atom_i, cutoff))
    return environments


def _ase_to_json(structure):
    """Implementation of structures_to_json for ase.Atoms"""
    data = {}
    data["size"] = len(structure)
    data["names"] = list(structure.symbols)
    data["x"] = structure.positions[:, 0].tolist()
    data["y"] = structure.positions[:, 1].tolist()
    data["z"] = structure.positions[:, 2].tolist()

    if (structure.cell.lengths() != [0.0, 0.0, 0.0]).all():
        data["cell"] = list(np.concatenate(structure.cell))

    return data


def _remove_invalid_properties(properties, origin):
    """
    Remove invalid properties from the ``properties`` dictionary. ``origin`` is
    used in error messages as the property origin
    """
    to_remove = set()

    for name, property in properties.items():
        values = property["values"]

        for value in values:
            # check type consistency
            if not _is_convertible_to_property(value):
                warnings.warn(
                    f"value '{value}' of type '{type(value)}' for the '{name}' "
                    f"property from {origin} is not convertible to float, array or "
                    "string, this property will be ignored.",
                    stacklevel=2,
                )
                to_remove.add(name)
                break

        if (
            name in to_remove
            or not hasattr(values[0], "__len__")
            or isinstance(values[0], str)
        ):
            continue

        # check length consistency
        lengths = [len(v) for v in property["values"]]
        if len(set(lengths)) != 1:
            to_remove.add(name)
            warnings.warn(
                f"values of the property '{name}' have inconsistent length across "
                "different structures, it will be ignored",
                stacklevel=2,
            )

    for name in to_remove:
        del properties[name]


def _is_convertible_to_property(value):
    """
    Check whether a value is convertible to a chemiscope property, i.e. if it is
    a string or something convertible to float.
    """
    if isinstance(value, (bytes, str)) or (
        (isinstance(value, list) or isinstance(value, np.ndarray))
        and isinstance(value[0], str)
    ):
        # string types
        return True
    else:
        # everything convertible to float
        try:
            float(value)
            return True
        except Exception:
            try:
                np.array(value, dtype=float)
                return True
            except Exception:
                return False


def _extract_key_from_ase(structure, key, target=None):
    """
    Extract a property from either the atomic array of info fields of an ``ase.Atoms``
    structure. Defaults to atoms if no target is specified, and also returns the actual
    target it picked the key from.
    """

    arrays = structure.arrays.copy()
    info = structure.info.copy()
    # ase >= 3.23 workaround
    if structure.calc is not None and "forces" in structure.calc.results:
        arrays["forces"] = structure.calc.results["forces"]
    if structure.calc is not None and "energy" in structure.calc.results:
        info["energy"] = structure.calc.results["energy"]

    if target is None:
        if key in arrays:  # defaults to atom target
            target = "atom"
        elif key in info:
            target = "structure"
        else:
            raise IndexError(
                f"Key {key} not found in neither `Atoms.arrays` or `Atoms.info`"
            )

    if target == "atom":
        try:
            values = arrays[key]
        except IndexError:
            raise IndexError(f"Key {key} not found in `Atoms.arrays`")
    if target == "structure":
        try:
            values = info[key]
        except IndexError:
            raise IndexError(f"Key {key} not found in `Atoms.info`")

    return values, target


def ase_vectors_to_arrows(
    structures=None,
    key=None,
    target=None,
    *,
    scale=1.0,
    radius=0.1,
    head_radius_scale=1.75,
    head_length_scale=2.0,
    frames=None,
):
    """
    Extract a vectorial atom property from a list of ``ase.Atoms``, and returns
    a list of arrow shapes. Besides the specific parameters it also accepts the same
    parameters as :py:func:`arrow_from_vector`, which are used to define the style of
    the arrows.

    :param structures: list of ASE Atoms
    :param key: name of the ASE atom property. Should contain three components
       corresponding to x,y,z
    :param target: whether the properties should be associated with the entire
       structure, or each atom (``structure`` or ``atom``). defaults to autodetection
    """
    if frames is not None:
        warnings.warn(
            "`frames` parameter is deprecated, please use `structures` instead",
            stacklevel=2,
        )

        if structures is not None:
            raise ValueError("cannot use both `structures` and `frames` parameters")

        structures = frames

    if key is None:
        # this is only required while we keep support for the deprecated `frames`
        # parameter. We can remove the default None key when we remove `frames`
        raise ValueError("`key` parameter must be provided")

    vectors = []

    # set shape parameters globally to reduce size of the arrow list
    globs = {}
    globs["baseRadius"] = radius
    globs["headRadius"] = globs["baseRadius"] * head_radius_scale
    globs["headLength"] = globs["baseRadius"] * head_length_scale

    for structure in structures:
        values, target = _extract_key_from_ase(structure, key, target)
        if target == "structure":
            values = values.reshape(1, -1)
        if len(values.shape) != 2 or values.shape[1] != 3:
            raise ValueError(
                f"Property array {key} has not the shape of a list of 3-vectors"
            )

        # makes a list of arrows to visualize the property
        vectors = vectors + [
            arrow_from_vector(
                v,
                scale=scale,
                radius=None,
                head_radius_scale=None,
                head_length_scale=None,
            )
            for v in values
        ]

    if target == "atom":
        return {"kind": "arrow", "parameters": {"global": globs, "atom": vectors}}
    else:
        return {"kind": "arrow", "parameters": {"global": globs, "structure": vectors}}


def ase_tensors_to_ellipsoids(
    structures=None,
    key=None,
    target=None,
    *,
    scale=1.0,
    force_positive=False,
    frames=None,
):
    """
    Extract a 2-tensor atom property from a list of ``ase.Atoms``, and returns a list of
    ellipsoids shapes. Besides the specific parameters it also accepts the same
    parameters as :py:func:`ellipsoid_from_tensor`, which are used to draw the shapes

    :param structures: list of ASE Atoms
    :param key: name of the ASE atom property. Should contain nine components
       corresponding to xx,xy,xz,yx,yy,yz,zx,zy,zz or six components corresponding to
       xx,yy,zz,xy,xz,yz
    :param target: whether the properties should be associated with the entire
       structure, or each atom (``structure`` or ``atom``). defaults to autodetection

    :param scale: see :py:func:`ellipsoid_from_tensor`
    :param force_positive: see :py:func:`ellipsoid_from_tensor`
    """

    if frames is not None:
        warnings.warn(
            "`frames` parameter is deprecated, please use `structures` instead",
            stacklevel=2,
        )

        if structures is not None:
            raise ValueError("cannot use both `structures` and `frames` parameters")

        structures = frames

    if key is None:
        # this is only required while we keep support for the deprecated `frames`
        # parameter. We can remove the default None key when we remove `frames`
        raise ValueError("`key` parameter must be provided")

    tensors = []

    for structure in structures:
        values, target = _extract_key_from_ase(structure, key, target)
        if target == "structure":
            values = values.reshape(1, -1)
        if len(values.shape) != 2 or (values.shape[1] != 6 and values.shape[1] != 9):
            raise ValueError(
                f"Property array {key} has not the shape of a list of 6 or 9-vectors"
            )

        # makes a list of ellipsoids to visualize the property
        tensors = tensors + [
            ellipsoid_from_tensor(v, scale=scale, force_positive=force_positive)
            for v in values
        ]

    if target == "atom":
        return dict(kind="ellipsoid", parameters={"global": {}, "atom": tensors})
    else:
        return dict(kind="ellipsoid", parameters={"global": {}, "structure": tensors})


# Required parameters from different kinds of shapes
SHAPE_PARAMETERS = {
    "ellipsoid": "semiaxes",
    "sphere": "radius",
}


def _extract_lammps_shapes(structure, key):
    if key in structure.info:
        if structure.info[key] not in SHAPE_PARAMETERS:
            raise KeyError(
                "The currently-supported shape in `extract_lammps_shapes_from_ase` are "
                f"{list(SHAPE_PARAMETERS.keys())}, received '{structure.info[key]}'"
            )

        shape = _get_shape_params(key, structure.info[key], structure.info)
        if "orientation" in structure.arrays:
            return {
                key: [
                    {**shape, "orientation": list(o)}
                    for o in structure.arrays["orientation"]
                ]
            }
        else:
            return {key: [shape for _ in structure]}

    elif key in structure.arrays:
        shapes = []
        for atom_i, shape_key in enumerate(structure.arrays[key]):
            if shape_key not in SHAPE_PARAMETERS:
                raise KeyError(
                    "The currently-supported shape types are {}, received {}.".format(
                        ", ".join(SHAPE_PARAMETERS.keys()), shape_key
                    )
                )
            shape = _get_shape_params_atom(
                key,
                shape_key,
                structure.arrays,
                atom_i,
            )

            if "orientation" in structure.arrays:
                shape["orientation"] = list(structure.arrays["orientation"][atom_i])

            shapes.append(shape)

        return {key: shapes}


def _get_shape_params(prefix, shape_kind, dictionary):
    shape = {"kind": shape_kind}
    parameter = SHAPE_PARAMETERS[shape_kind]
    try:
        shape[parameter] = dictionary[f"{prefix}_{parameter}"]
    except KeyError:
        raise KeyError(
            f"Missing required parameter '{prefix}_{parameter}' for "
            f"'{shape_kind}' shape"
        )

    return shape


def _get_shape_params_atom(prefix, shape_kind, dictionary, atom_i):
    """Extract shape parameters for a single atom"""

    shape = {"kind": shape_kind}
    parameter = SHAPE_PARAMETERS[shape_kind]
    try:
        shape[parameter] = dictionary[f"{prefix}_{parameter}"][atom_i]
    except KeyError:
        raise KeyError(
            f"Missing required parameter '{prefix}_{parameter}' for "
            f"'{shape_kind}' shape"
        )

    return shape
