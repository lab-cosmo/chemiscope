import warnings

import numpy as np

from ._chemfiles import STANDARD_RESIDUES


try:
    import MDAnalysis as mda

    HAVE_MDA = True
except ImportError:
    HAVE_MDA = False

BIO_PROPERTIES = ["hetatom", "resids", "resnames", "chains"]


def _mda_valid_structures(structures):
    if HAVE_MDA and isinstance(structures, mda.AtomGroup):
        return structures, True
    else:
        return [], False


def _mda_to_json(ag):
    data = {}
    data["size"] = len(ag)
    data["elements"] = (
        ag.elements.tolist() if hasattr(ag, "elements") else ag.types.tolist()
    )
    # `element` is better, but not always available, e.g. xyz file
    if hasattr(ag, "names"):
        data["names"] = ag.names.tolist()
    else:
        data["names"] = data["elements"]
    x, y, z = ag.positions.T
    data["x"] = x.tolist()
    data["y"] = y.tolist()
    data["z"] = z.tolist()
    hetatom = np.full(ag.n_atoms, True)

    if ag.dimensions is not None:
        data["cell"] = list(
            np.concatenate(
                mda.lib.mdamath.triclinic_vectors(ag.dimensions),
                dtype=np.float64,
                # should be np.float64 otherwise not serializable
            )
        )
    if hasattr(ag, "chainIDs") and ag.chainIDs is not None:
        data["chains"] = ag.chainIDs.tolist()
    elif hasattr(ag, "segids") and ag.segids is not None:
        # segids are sometimes abused to store chain ids in PDBs, so we use them here
        data["chains"] = ag.segids.tolist()

    if hasattr(ag, "resnames") and ag.resnames is not None:
        data["resnames"] = [
            resname if resname is not None and resname != "" else "UNK"
            for resname in ag.resnames
        ]
        # atom selection requires the `resname`. `select_atoms` numbers its result
        # in the universe, not in this group, so a group that is not the whole
        # universe has to be matched by index instead of sliced with them
        selected = ag.select_atoms("protein or nucleic")
        hetatom[np.isin(ag.indices, selected.indices)] = False
        # MDAnalysis does not know all force-field names of nucleotides
        hetatom[np.isin(data["resnames"], list(STANDARD_RESIDUES))] = False

    if hasattr(ag, "resids") and ag.resids is not None:
        data["resids"] = ag.resids.view(dtype=int).tolist()

    if hasattr(ag, "bonds") and ag.bonds is not None:
        # bond indices are numbered in the universe as well, so they have to be
        # mapped into this group; a bond with one atom outside the group is dropped
        local = np.full(ag.universe.atoms.n_atoms, -1)
        local[ag.indices] = np.arange(ag.n_atoms)
        pairs = local[ag.bonds.indices] if len(ag.bonds) > 0 else np.empty((0, 2), int)
        pairs = pairs[(pairs >= 0).all(axis=1)]
        # a PDB only has to CONECT its non-standard residues, so MDAnalysis
        # often ends up knowing the bonds of a ligand and nothing else. The
        # viewer works bonds out from the geometry when it is given none, but
        # trusts the list when there is one, and would then draw the rest of
        # the structure as loose atoms. So only pass the bonds on when every
        # residue that could have one does.
        bonded = np.zeros(ag.n_atoms, dtype=bool)
        bonded[pairs.ravel()] = True
        residue = getattr(ag, "resindices", np.zeros(ag.n_atoms, dtype=int))
        _, residue = np.unique(residue, return_inverse=True)
        sizes = np.bincount(residue)
        if len(pairs) > 0 and np.all((sizes == 1) | (np.bincount(residue, bonded) > 0)):
            data["bonds"] = np.hstack((pairs, np.full((len(pairs), 1), 1))).tolist()

    data["hetatom"] = hetatom.tolist()

    # remove bio-related properties if any of them are missing
    existing_properties = []
    for prop in BIO_PROPERTIES:
        if prop in data:
            existing_properties.append(prop)

    if len(existing_properties) != len(BIO_PROPERTIES):
        for prop in existing_properties:
            del data[prop]

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
            f"of structures: {list(sorted(extra))}; they will be ignored",
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


def _mda_all_atomic_environments(ag, cutoff):
    "Extract all atomic environments out of a set of MDAnalysis AtomGroup objects"
    environments = []
    for structure_i, structure in enumerate(ag.universe.trajectory):
        for atom_i in range(len(structure)):
            environments.append((structure_i, atom_i, cutoff))
    return environments
