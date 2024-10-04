import warnings
from inspect import signature

import numpy as np

from ._shapes import arrow_from_vector, ellipsoid_from_tensor

try:
    import ase

    HAVE_ASE = True
except ImportError:
    HAVE_ASE = False


def _ase_valid_structures(frames):
    try:
        frames_list = list(frames)
    except TypeError:
        return [], False

    if HAVE_ASE and isinstance(frames_list[0], ase.Atoms):
        for frame in frames_list:
            assert isinstance(frame, ase.Atoms)
        return frames, True
    elif HAVE_ASE and isinstance(frames_list[0], ase.Atom):
        # deal with the user passing a single frame
        return [frames], True
    else:
        return frames, False


def _ase_list_atom_properties(frames):
    IGNORED_ASE_ARRAYS = ["positions", "numbers", "center_atoms_mask"]
    # extract the set of common properties between all frames
    all_names = set()
    extra = set()

    for name in frames[0].arrays.keys():
        if name in IGNORED_ASE_ARRAYS:
            continue
        all_names.add(name)

    for frame in frames[1:]:
        for name in frame.arrays.keys():
            if name in IGNORED_ASE_ARRAYS:
                continue

            if name not in all_names:
                extra.add(name)

        remove = []
        for name in all_names:
            if name not in frame.arrays.keys():
                remove.append(name)

        for name in remove:
            all_names.remove(name)
            extra.add(name)

    if len(extra) != 0:
        warnings.warn(
            "the following atomic properties are only defined for a subset "
            f"of frames: {list(sorted(extra))}; they will be ignored",
            stacklevel=2,
        )

    return all_names


def _ase_list_structure_properties(frames):
    # extract the set of common properties between all frames
    all_names = set()
    extra = set()

    for name in frames[0].info.keys():
        all_names.add(name)

    for frame in frames[1:]:
        for name in frame.info.keys():
            if name not in all_names:
                extra.add(name)

        remove = []
        for name in all_names:
            if name not in frame.info.keys():
                remove.append(name)

        for name in remove:
            all_names.remove(name)
            extra.add(name)

    if len(extra) != 0:
        warnings.warn(
            "the following structure properties are only defined for a subset "
            f"of frames: {list(sorted(extra))}; they will be ignored",
            stacklevel=2,
        )

    return all_names


def _ase_atom_properties(frames, only, atoms_mask=None):
    all_names = _ase_list_atom_properties(frames)
    if only is not None:
        all_names = [name for name in all_names if name in only]

    # create property in the format expected by create_input
    properties = {
        name: {"target": "atom", "values": value}
        for name, value in frames[0].arrays.items()
        if name in all_names
    }

    for frame in frames[1:]:
        for name, value in frame.arrays.items():
            if name not in all_names:
                continue
            properties[name]["values"] = np.concatenate(
                [properties[name]["values"], value]
            )

    _remove_invalid_properties(properties, "ASE")

    if atoms_mask is not None:
        # only include values for requested atoms
        for property in properties.values():
            property["values"] = property["values"][atoms_mask]

    return properties


def _ase_structure_properties(frames, only=None):
    all_names = _ase_list_structure_properties(frames)
    if only is not None:
        all_names = [name for name in all_names if name in only]

    # create property in the format expected by create_input
    properties = {name: {"target": "structure", "values": []} for name in all_names}

    for frame in frames:
        for name, value in frame.info.items():
            if name in all_names:
                properties[name]["values"].append(value)

    _remove_invalid_properties(properties, "ASE")

    return properties


def _ase_extract_properties(frames, only=None, environments=None):
    """implementation of ``extract_properties`` for ASE"""

    properties = _ase_structure_properties(frames, only)

    if environments is not None:
        atoms_mask = [[False] * len(f) for f in frames]
        for structure, center, _ in environments:
            atoms_mask[structure][center] = True

        atoms_mask = np.concatenate(atoms_mask)
    else:
        atoms_mask = None

    atom_properties = _ase_atom_properties(frames, only, atoms_mask)

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


def _ase_all_atomic_environments(frames, cutoff):
    "Extract all atomic environments out of a set of ASE Atoms objects"
    environments = []
    for structure_i, frame in enumerate(frames):
        for atom_i in range(len(frame)):
            environments.append((structure_i, atom_i, cutoff))
    return environments


def _ase_to_json(frame):
    """Implementation of frame_to_json for ase.Atoms"""
    data = {}
    data["size"] = len(frame)
    data["names"] = list(frame.symbols)
    data["x"] = [float(value) for value in frame.positions[:, 0]]
    data["y"] = [float(value) for value in frame.positions[:, 1]]
    data["z"] = [float(value) for value in frame.positions[:, 2]]

    if (frame.cell.lengths() != [0.0, 0.0, 0.0]).all():
        data["cell"] = list(np.concatenate(frame.cell))

    return data


def _remove_invalid_properties(properties, origin):
    """
    Remove invalid properties from the ``properties`` dictionary. ``origin`` is
    used in error messages as the property origin
    """
    to_remove = []
    for name, property in properties.items():
        for value in property["values"]:
            if not _is_convertible_to_property(value):
                warnings.warn(
                    f"value '{value}' of type '{type(value)}' for the '{name}' "
                    f"property from {origin} is not convertible to float, array or "
                    "string, this property will be ignored.",
                    stacklevel=2,
                )
                to_remove.append(name)
                break

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


def _extract_key_from_ase(frame, key, target=None):
    """
    Extract a property from either the atomic array of info fields
    of an ase.Atoms frame. Defaults to atoms if no target is specified,
    and also returns the actual target it picked the key from.
    """

    if target is None:
        if key in frame.arrays:  # defaults to atom target
            target = "atom"
        elif key in frame.info:
            target = "structure"
        else:
            raise IndexError(
                f"Key {key} not found in neither `Atoms.arrays` or `Atoms.info`"
            )

    if target == "atom":
        try:
            values = frame.arrays[key]
        except IndexError:
            raise IndexError(f"Key {key} not found in `Atoms.arrays`")
    if target == "structure":
        try:
            values = frame.info[key]
        except IndexError:
            raise IndexError(f"Key {key} not found in `Atoms.info`")

    return values, target


def ase_vectors_to_arrows(frames, key="forces", target=None, **kwargs):
    """
    Extract a vectorial atom property from a list of ase.Atoms
    objects, and returns a list of arrow shapes. Besides the specific
    parameters it also accepts the same parameters as
    :py:func:`arrow_from_vector`, which are used to define the
    style of the arrows.

    :param frames: list of ASE Atoms objects
    :param key: name of the ASE atom property. Should contain
       three components corresponding to x,y,z
    :param target: whether the properties should be associated with
       the entire structure, or each atom (`structure` or `atom`).
       defaults to autodetection
    """

    vectors = []

    # set shape parameters globally to reduce size of the arrow list
    globs = {}
    defaults = signature(arrow_from_vector).parameters
    globs["baseRadius"] = (
        kwargs.pop("radius") if ("radius" in kwargs) else defaults["radius"].default
    )
    globs["headRadius"] = globs["baseRadius"] * (
        kwargs.pop("head_radius_scale")
        if ("head_radius_scale" in kwargs)
        else defaults["head_radius_scale"].default
    )
    globs["headLength"] = globs["baseRadius"] * (
        kwargs.pop("head_length_scale")
        if ("head_length_scale" in kwargs)
        else defaults["head_length_scale"].default
    )

    for f in frames:
        values, target = _extract_key_from_ase(f, key, target)
        if target == "structure":
            values = values.reshape(1, -1)
        if len(values.shape) != 2 or values.shape[1] != 3:
            raise ValueError(
                f"Property array {key} has not the shape of a list of 3-vectors"
            )

        # makes a list of arrows to visualize the property
        vectors = vectors + [
            arrow_from_vector(
                v, radius=None, head_radius_scale=None, head_length_scale=None, **kwargs
            )
            for v in values
        ]

    if target == "atom":
        return {"kind": "arrow", "parameters": {"global": globs, "atom": vectors}}
    else:
        return {"kind": "arrow", "parameters": {"global": globs, "structure": vectors}}


def ase_tensors_to_ellipsoids(frames, key, target=None, **kwargs):
    """
    Extract a 2-tensor atom property from a list of ase.Atoms
    objects, and returns a list of ellipsoids shapes. Besides the specific
    parameters it also accepts the same parameters as
    `ellipsoid_from_tensor`, which are used to draw the shapes

    :param frames: list of ASE Atoms objects
    :param key: name of the ASE atom property. Should contain
       nine components corresponding to xx,xy,xz,yx,yy,yz,zx,zy,zz or
       six components corresponding to xx,yy,zz,xy,xz,yz
    :param target: whether the properties should be associated with
       the entire structure, or each atom (`structure` or `atom`).
       defaults to autodetection
    """

    tensors = []

    for f in frames:
        values, target = _extract_key_from_ase(f, key, target)
        if target == "structure":
            values = values.reshape(1, -1)
        if len(values.shape) != 2 or (values.shape[1] != 6 and values.shape[1] != 9):
            raise ValueError(
                f"Property array {key} has not the shape of a list of 6 or 9-vectors"
            )

        # makes a list of ellipsoids to visualize the property
        tensors = tensors + [ellipsoid_from_tensor(v, **kwargs) for v in values]

    if target == "atom":
        return dict(kind="ellipsoid", parameters={"global": {}, "atom": tensors})
    else:
        return dict(kind="ellipsoid", parameters={"global": {}, "structure": tensors})


# Required parameters from different kinds of shapes
SHAPE_PARAMETERS = {
    "ellipsoid": "semiaxes",
    "sphere": "radius",
}


def _extract_lammps_shapes(frame, key):
    if key in frame.info:
        if frame.info[key] not in SHAPE_PARAMETERS:
            raise KeyError(
                "The currently-supported shape in `extract_lammps_shapes_from_ase` are "
                f"{list(SHAPE_PARAMETERS.keys())}, received '{frame.info[key]}'"
            )

        shape = _get_shape_params(key, frame.info[key], frame.info)
        if "orientation" in frame.arrays:
            return {
                key: [
                    {**shape, "orientation": list(o)}
                    for o in frame.arrays["orientation"]
                ]
            }
        else:
            return {key: [shape for _ in frame]}

    elif key in frame.arrays:
        shapes = []
        for atom_i, shape_key in enumerate(frame.arrays[key]):
            if shape_key not in SHAPE_PARAMETERS:
                raise KeyError(
                    "The currently-supported shape types are {}, received {}.".format(
                        ", ".join(SHAPE_PARAMETERS.keys()), shape_key
                    )
                )
            shape = _get_shape_params_atom(
                key,
                shape_key,
                frame.arrays,
                atom_i,
            )

            if "orientation" in frame.arrays:
                shape["orientation"] = list(frame.arrays["orientation"][atom_i])

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


def ase_merge_pi_frames(bead_frames):
    """
    Takes a list of lists of ase.Atoms objects, corresponding to
    the trajectories of the beads of a path integral calculation.
    Returns a consolidated trajectory in which all beads are in the
    same frame, together with properties and shapes that help visualize
    the path integral. The various pieces are returned as a dictionary
    and can be passed directly to `create_input` and related functions.

    :param bead_frames: list of lists of ASE Atoms objects. the outer
       list corresponds to the bead index
    """

    nbeads, nframes, natoms = (
        len(bead_frames),
        len(bead_frames[0]),
        len(bead_frames[0][0]),
    )

    # creates frames with all beads
    full_frames = []
    for i in range(nframes):
        struc = bead_frames[0][i].copy()
        for k in range(1, nbeads):
            struc += bead_frames[k][i]
        full_frames.append(struc)

    atom_id = [
        j
        for i in range(nframes)  # frame
        for k in range(nbeads)  # bead
        for j in range(natoms)  # atom
    ]
    bead_id = [
        k
        for i in range(nframes)  # frame
        for k in range(nbeads)  # bead
        for j in range(natoms)  # atom
    ]

    paths_shapes = dict(
        kind="cylinder",
        parameters={
            "global": {"radius": 0.2},
            "atom": [
                {
                    "vector": full_frames[i]
                    .get_distance(
                        k * natoms + j,
                        ((k + 1) % nbeads) * natoms + j,
                        mic=True,
                        vector=True,
                    )
                    .tolist(),
                }
                for i in range(nframes)  # frame
                for k in range(nbeads)  # bead
                for j in range(natoms)  # atom
            ],
        },
    )
    return {
        "frames": full_frames,
        "shapes": {"paths": paths_shapes},
        "properties": {"atom_id": atom_id, "bead_id": bead_id},
        "settings": {
            "structure": [
                {
                    "atoms": True,
                    "bonds": False,
                    "shape": "paths",
                    "unitCell": bool(full_frames[0].pbc[0]),
                    "environments": {"activated": False},
                    "color": {"property": "atom_id", "palette": "hsv (periodic)"},
                }
            ]
        },
        "environments": _ase_all_atomic_environments(full_frames, 4.0),
    }
