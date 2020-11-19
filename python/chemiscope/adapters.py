# -*- coding: utf-8 -*-
import numpy as np
import warnings

try:
    import ase

    HAVE_ASE = True
except ImportError:
    HAVE_ASE = False


def frames_to_json(frames):
    frames = list(frames)

    if HAVE_ASE and isinstance(frames[0], ase.Atoms):
        return [_ase_to_json(frame) for frame in frames]
    elif HAVE_ASE and isinstance(frames[0], ase.Atom):
        raise Exception(
            "expected ase.Atoms, got ase.Atom. "
            + "Try passing a list of frames instead of a single frame "
            + "(`[frame]` instead of `frame`)"
        )
    else:
        raise Exception(f"unknown frame type: '{frames[0].__class__.__name__}'")


def atom_properties(frames):
    frames = list(frames)

    if HAVE_ASE and isinstance(frames[0], ase.Atoms):
        return _ase_atom_properties(frames)
    elif HAVE_ASE and isinstance(frames[0], ase.Atom):
        raise Exception(
            "expected ase.Atoms, got ase.Atom. "
            + "Try passing a list of frames instead of a single frame "
            + "(`[frame]` instead of `frame`)"
        )
    else:
        raise Exception(f"unknown frame type: '{frames[0].__class__.__name__}'")


def structure_properties(frames):
    frames = list(frames)

    if HAVE_ASE and isinstance(frames[0], ase.Atoms):
        return _ase_structure_properties(frames)
    elif HAVE_ASE and isinstance(frames[0], ase.Atom):
        raise Exception(
            "expected ase.Atoms, got ase.Atom. "
            + "Try passing a list of frames instead of a single frame "
            + "(`[frame]` instead of `frame`)"
        )
    else:
        raise Exception(f"unknown frame type: '{frames[0].__class__.__name__}'")


def _ase_to_json(frame):
    data = {}
    data["size"] = len(frame)
    data["names"] = list(frame.symbols)
    data["x"] = [float(value) for value in frame.positions[:, 0]]
    data["y"] = [float(value) for value in frame.positions[:, 1]]
    data["z"] = [float(value) for value in frame.positions[:, 2]]

    if (frame.cell.lengths() != [0.0, 0.0, 0.0]).all():
        data["cell"] = list(np.concatenate(frame.cell))

    return data


IGNORED_ASE_ARRAYS = ["positions", "numbers"]


def _ase_atom_properties(frames):
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
    return properties


def _ase_structure_properties(frames):
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
            "the following structure properties properties are only defined "
            + f"for a subset of frames: {list(sorted(extra))}; they will be ignored"
        )

    # create property in the format expected by create_input
    properties = {name: {"target": "structure", "values": []} for name in all_names}

    for frame in frames:
        for name, value in frame.info.items():
            if name in all_names:
                properties[name]["values"].append(value)

    return properties
