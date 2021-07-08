# -*- coding: utf-8 -*-
import numpy as np
import warnings

from collections import Counter

try:
    import ase

    HAVE_ASE = True
except ImportError:
    HAVE_ASE = False


def frames_to_json(frames):
    """
    Convert the given ``frames`` to the JSON structure used by chemiscope.

    This function is a shim calling specialized implementations for all the
    supported frame types. Currently only `ase.Atoms` frames are supported.

    :param frames: iterable over structures (typically a list of frames)
    """
    frames_list = list(frames)

    if HAVE_ASE and isinstance(frames_list[0], ase.Atoms):
        return [_ase_to_json(frame) for frame in frames_list]
    elif HAVE_ASE and isinstance(frames_list[0], ase.Atom):
        # deal with the user passing a single frame
        return frames_to_json([frames])
    else:
        raise Exception(f"unknown frame type: '{frames_list[0].__class__.__name__}'")


def atom_properties(frames, composition):
    """
    Extract "atom" properties from the given ``frames``, and give them as a
    dictionary compatible with :py:func:`create_input`.

    This function is a shim calling specialized implementations for all the
    supported frame types. Currently only `ase.Atoms` frames are supported.

    :param frames: iterable over structures (typically a list of frames)
    """
    frames_list = list(frames)

    if HAVE_ASE and isinstance(frames_list[0], ase.Atoms):
        return _ase_atom_properties(frames_list, composition)
    elif HAVE_ASE and isinstance(frames_list[0], ase.Atom):
        # deal with the user passing a single frame
        return atom_properties([frames], composition)
    else:
        raise Exception(f"unknown frame type: '{frames_list[0].__class__.__name__}'")


def structure_properties(frames, composition):
    """
    Extract "structure" properties from the given ``frames``, and give them as a
    dictionary compatible with :py:func:`create_input`.

    This function is a shim calling specialized implementations for all the
    supported frame types. Currently only `ase.Atoms` frames are supported.

    :param frames: iterable over structures (typically a list of frames)
    """
    frames_list = list(frames)

    if HAVE_ASE and isinstance(frames_list[0], ase.Atoms):
        return _ase_structure_properties(frames_list, composition)
    elif HAVE_ASE and isinstance(frames_list[0], ase.Atom):
        # deal with the user passing a single frame
        return structure_properties([frames], composition)
    else:
        raise Exception(f"unknown frame type: '{frames_list[0].__class__.__name__}'")


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


def _ase_atom_properties(frames, composition):
    """Implementation of atom_properties for ase.Atoms"""
    IGNORED_ASE_ARRAYS = ["positions", "numbers"]
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
