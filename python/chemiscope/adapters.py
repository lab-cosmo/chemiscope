# -*- coding: utf-8 -*-
import numpy as np

try:
    import ase

    HAVE_ASE = True
except ImportError:
    HAVE_ASE = False


def frames_to_json(frames):
    frames = list(frames)

    if HAVE_ASE and isinstance(frames[0], ase.Atoms):
        return [_ase_to_json(frame) for frame in frames]
    else:
        raise Exception(f"unknown frame type: '{frames[0].__class__.__name__}'")


def atom_properties(frames):
    frames = list(frames)

    if HAVE_ASE and isinstance(frames[0], ase.Atoms):
        return [_ase_atom_properties(frame) for frame in frames]
    else:
        raise Exception(f"unknown frame type: '{frames[0].__class__.__name__}'")


def structure_properties(frames):
    frames = list(frames)

    if HAVE_ASE and isinstance(frames[0], ase.Atoms):
        return [_ase_structure_properties(frame) for frame in frames]
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
    properties = {
        name: {"target": "atom", "values": value}
        for name, value in frames[0].arrays.items()
        if name not in IGNORED_ASE_ARRAYS
    }

    for frame in frames[1:]:
        for name, value in frame.arrays.items():
            if name in IGNORED_ASE_ARRAYS:
                continue
            properties[name]["values"] = np.concatenate(
                [properties[name]["values"], value]
            )
    return properties


def _ase_structure_properties(frames):
    properties = {
        name: {"target": "structure", "values": []} for name in frames[0].info.keys()
    }

    for frame in frames:
        for name, value in frame.info.items():
            properties[name]["values"].append(value)

    return properties
