import warnings

import numpy as np


try:
    import MDAnalysis as mda

    HAVE_MDA = True
except ImportError:
    HAVE_MDA = False


def _mda_valid_structures(frames):
    if HAVE_MDA and isinstance(frames, mda.AtomGroup):
        return frames, True
    else:
        return [], False


def _mda_to_json(ag):
    data = {}
    data["size"] = len(ag)
    data["elements"] = [
        atom.element if hasattr(atom, "element") else atom.type
        for atom in ag
        # `element` is better, but not always available, e.g. xyz file
    ]
    if hasattr(ag, "names"):
        data["names"] = [atom.name for atom in ag]
    else:
        data["names"] = data["elements"]
    data["x"] = [float(value) for value in ag.positions[:, 0]]
    data["y"] = [float(value) for value in ag.positions[:, 1]]
    data["z"] = [float(value) for value in ag.positions[:, 2]]
    if ag.dimensions is not None:
        data["cell"] = list(
            np.concatenate(
                mda.lib.mdamath.triclinic_vectors(ag.dimensions),
                dtype=np.float64,
                # should be np.float64 otherwise not serializable
            )
        )
    data["hetatom"] = [True for _ in ag]
    if hasattr(ag, "chainIDs") and ag.chainIDs is not None:
        data["chains"] = [atom.chainID for atom in ag]
    if hasattr(ag, "resnames") and ag.resnames is not None:
        data["residues"] = [atom.resname for atom in ag]
        # atom selection requires the `resname`
        for idx in ag.select_atoms("protein or nucleic").indices:
            data["hetatom"][idx] = False
    if hasattr(ag, "resids") and ag.resids is not None:
        data["resids"] = [int(atom.resid) for atom in ag]

    return data


def _mda_get_atom_properties(ag) -> dict:
    POSSIBLE_PROPERTIES = ["velocities", "forces"]
    all_properties = {}
    extra = set()
    for property in POSSIBLE_PROPERTIES:
        try:
            all_properties[property] = []
            all_properties[property].append([getattr(ag, property)])
        except mda.exceptions.NoDataError:
            del all_properties[property]

    if len(ag.universe.trajectory) > 1:
        for _ in ag.universe.trajectory[1:]:
            for property in all_properties:
                try:
                    all_properties[property].append([getattr(ag, property)])
                except mda.exceptions.NoDataError:
                    extra.add(property)

    if len(extra) != 0:
        warnings.warn(
            "the following per-atom properties are only defined for a subset "
            f"of frames: {list(sorted(extra))}; they will be ignored",
            stacklevel=2,
        )

    return all_properties


def _mda_atom_properties(ag, only=None, atoms_mask=None):
    all_properties = _mda_get_atom_properties(ag)
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

    if atoms_mask is not None:
        # only include values for requested atoms
        for property in properties.values():
            property["values"] = property["values"][atoms_mask]

    return properties


def _mda_extract_properties(ag, only=None, environments=None):
    """implementation of ``extract_properties`` for MDAnalysis"""
    properties = {}

    if environments is not None:
        atoms_mask = [[False] * len(f) for f in ag.universe.trajectory]
        for structure, center, _ in environments:
            atoms_mask[structure][center] = True

        atoms_mask = np.concatenate(atoms_mask)
    else:
        atoms_mask = None

    atom_properties = _mda_atom_properties(ag, only, atoms_mask)

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
