# -*- coding: utf-8 -*-
import warnings
from collections import Counter

import numpy as np

try:
    import ase

    HAVE_ASE = True
except ImportError:
    HAVE_ASE = False


def _guess_adapter(frames):
    """
    Guess which adapter to use for the given frames. This function return the
    frames as a list and a string describing which adapter should be used.
    """
    frames_list = list(frames)
    if HAVE_ASE and isinstance(frames_list[0], ase.Atoms):
        for frame in frames_list:
            assert isinstance(frame, ase.Atoms)
        return frames, "ASE"
    elif HAVE_ASE and isinstance(frames_list[0], ase.Atom):
        # deal with the user passing a single frame
        return [frames], "ASE"
    else:
        raise Exception(f"unknown frame type: '{frames_list[0].__class__.__name__}'")


def frames_to_json(frames):
    """
    Convert the given ``frames`` to the JSON structure used by chemiscope.

    This function is a shim calling specialized implementations for all the
    supported frame types. Currently only `ase.Atoms` frames are supported.

    :param frames: iterable over structures (typically a list of frames)
    """
    frames, adapter = _guess_adapter(frames)

    if adapter == "ASE":
        return [_ase_to_json(frame) for frame in frames]
    else:
        raise Exception("reached unreachable code")


def atom_properties(frames, composition, atoms_mask=None):
    """
    Extract "atom" properties from the given ``frames``, and give them as a
    dictionary compatible with :py:func:`create_input`.

    This function is a shim calling specialized implementations for all the
    supported frame types. Currently only `ase.Atoms` frames are supported.

    :param frames: iterable over structures (typically a list of frames)
    :param composition: whether to also add properties containing information
                        about the chemical composition of the system
    :param atoms_mask: optional list of booleans containing which atoms should
                       be include in the output
    """
    frames, adapter = _guess_adapter(frames)

    if adapter == "ASE":
        return _ase_atom_properties(frames, composition, atoms_mask)
    else:
        raise Exception("reached unreachable code")


def structure_properties(frames, composition):
    """
    Extract "structure" properties from the given ``frames``, and give them as a
    dictionary compatible with :py:func:`create_input`.

    This function is a shim calling specialized implementations for all the
    supported frame types. Currently only `ase.Atoms` frames are supported.

    :param frames: iterable over structures (typically a list of frames)
    :param composition: whether to also add properties containing information
                        about the chemical composition of the system
    """
    frames, adapter = _guess_adapter(frames)

    if adapter == "ASE":
        return _ase_structure_properties(frames, composition)
    else:
        raise Exception("reached unreachable code")


def all_atomic_environments(frames, cutoff=3.5):
    """
    Generate a list of environments containing all the atoms in the given
    ``frames``. The optional spherical ``cutoff`` radius is used to display the
    environments in chemiscope.

    :param frames: iterable over structures (typically a list of frames)
    :param float cutoff: spherical cutoff radius used when displaying the
                         environments
    """
    frames, adapter = _guess_adapter(frames)

    if adapter == "ASE":
        environments = []
        for structure_i, frame in enumerate(frames):
            for atom_i in range(len(frame)):
                environments.append((structure_i, atom_i, cutoff))
        return environments
    else:
        raise Exception("reached unreachable code")


def librascal_atomic_environments(frames, cutoff=3.5):
    """
    Generate the list of environments for the given ``frames``, matching the
    behavior used by librascal when computing descriptors for only a subset of
    the atomic centers. The optional spherical ``cutoff`` radius is used to
    display the environments in chemiscope.

    Only ``ase.Atoms`` are supported for the ``frames`` since that's what
    librascal uses.

    :param frames: iterable over ``ase.Atoms``
    :param float cutoff: spherical cutoff radius used when displaying the
                         environments
    """
    frames, adapter = _guess_adapter(frames)

    if adapter != "ASE":
        raise Exception("librascal_atomic_environments only supports ASE frames")

    environments = []
    for structure_i, frame in enumerate(frames):
        if "center_atoms_mask" in frame.arrays:
            atoms_iter = np.where(frame.arrays["center_atoms_mask"])[0]
        else:
            # use all atoms
            atoms_iter = range(len(frame))

        for atom_i in atoms_iter:
            environments.append((structure_i, atom_i, cutoff))

    return environments


def _add_structure_chemical_composition(frames):
    """Adds information about the chemical composition and count for every chemical species
    to the info dictionary for every frame as a structure property"""
    all_elements = set()
    for frame in frames:
        all_elements.update(frame.symbols)
    all_elements = set(all_elements)
    for frame in frames:
        frame.info["composition"] = str(frame.symbols)
        dict_composition = dict(Counter(frame.symbols))
        for element in all_elements:
            if element in dict_composition:
                frame.info[f"n_{element}"] = dict_composition[element]
            else:
                frame.info[f"n_{element}"] = 0
    return frames


def _add_atom_chemical_composition(frames):
    """Adds information about the atoms chemical species to the frame.arrays
    for every frame as an atom property"""
    for frame in frames:
        dict_elements = {"element": list(frame.symbols)}
        frame.arrays.update(dict_elements)
    return frames


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


def _ase_atom_properties(frames, composition, atoms_mask):
    """Implementation of atom_properties for ase.Atoms"""
    IGNORED_ASE_ARRAYS = ["positions", "numbers", "center_atoms_mask"]
    # extract the set of common properties between all frames
    all_names = set()
    extra = set()

    if composition:
        _add_atom_chemical_composition(frames)

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
            "the following atomic properties properties are only defined "
            + f"for a subset of frames: {list(sorted(extra))}; they will be ignored"
        )

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


def _ase_structure_properties(frames, composition):
    """Implementation of structure_properties for ase.Atoms"""
    # extract the set of common properties between all frames
    all_names = set()
    extra = set()

    if composition:
        _add_structure_chemical_composition(frames)

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
            "the following structure properties properties are only defined "
            + f"for a subset of frames: {list(sorted(extra))}; they will be ignored"
        )

    # create property in the format expected by create_input
    properties = {name: {"target": "structure", "values": []} for name in all_names}

    for frame in frames:
        for name, value in frame.info.items():
            if name in all_names:
                properties[name]["values"].append(value)

    _remove_invalid_properties(properties, "ASE")

    return properties


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
                    + f"property from {origin} is not convertible to float or "
                    + "string, this property will be ignored."
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
    if isinstance(value, (bytes, str)):
        # string types
        return True
    else:
        # everything convertible to float
        try:
            float(value)
            return True
        except Exception:
            return False
