import numpy as np

try:
    import MDAnalysis as mda

    HAVE_MDA = True
except ImportError:
    HAVE_MDA = False


def _mda_valid_structures(frames: mda.Universe):
    if HAVE_MDA and isinstance(frames, mda.Universe):
        return frames, True
    else:
        return [], False


def _mda_to_json(ag):
    data = {}
    data["size"] = len(ag)
    data["names"] = [atom.name for atom in ag]
    data["x"] = [float(value) for value in ag.positions[:, 0]]
    data["y"] = [float(value) for value in ag.positions[:, 1]]
    data["z"] = [float(value) for value in ag.positions[:, 2]]
    if ag.dimensions is not None:
        data["cell"] = list(
            np.concatenate(
                mda.lib.mdamath.triclinic_vectors(ag.dimensions), dtype=np.float64
                # should be np.float64 otherwise not serializable
            )
        )

    return data


def _mda_list_atom_properties(frames) -> list:
    # mda cannot have atom properties or structure properties, so skipping.
    return []


def _mda_list_structure_properties(frames) -> list:
    # mda cannot have atom properties or structure properties, so skipping.
    return []
