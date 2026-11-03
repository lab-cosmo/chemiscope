import numpy as np

from ._chemfiles import STANDARD_RESIDUES


try:
    import openmm.app
    import openmm.unit

    HAVE_OPENMM = True
except ImportError:
    HAVE_OPENMM = False


def _is_openmm(structure):
    # duck-typed: openmm.app.Modeller, PDBFile, PDBxFile, ...
    return isinstance(
        getattr(structure, "topology", None), openmm.app.Topology
    ) and hasattr(structure, "positions")


def _openmm_valid_structures(structures):
    if not HAVE_OPENMM:
        return structures, False
    elif _is_openmm(structures):
        # deal with the user passing a single structure
        return [structures], True
    elif (
        isinstance(structures, list)
        and len(structures) > 0
        and _is_openmm(structures[0])
    ):
        for structure in structures:
            assert _is_openmm(structure)
        return structures, True
    else:
        return structures, False


def _openmm_to_json(structure):
    """
    Implementation of ``structures_to_json`` for objects carrying an OpenMM
    ``topology`` and ``positions``.
    """
    topology = structure.topology
    atoms = list(topology.atoms())
    residues = [atom.residue for atom in atoms]
    angstrom = openmm.unit.angstrom
    x, y, z = np.asarray(structure.positions.value_in_unit(angstrom), dtype=float).T

    data = {
        "size": len(atoms),
        "names": [atom.name for atom in atoms],
        "elements": [
            "X" if atom.element is None else atom.element.symbol for atom in atoms
        ],
        "x": x.tolist(),
        "y": y.tolist(),
        "z": z.tolist(),
        "resnames": [residue.name or "UNK" for residue in residues],
        "resids": [int(residue.id) for residue in residues],
        "chains": [residue.chain.id for residue in residues],
        "hetatom": [residue.name not in STANDARD_RESIDUES for residue in residues],
    }

    box = topology.getPeriodicBoxVectors()
    if box is not None:
        cell = np.asarray(box.value_in_unit(angstrom), dtype=float)
        data["cell"] = cell.flatten().tolist()

    bonds = [[bond.atom1.index, bond.atom2.index, 1] for bond in topology.bonds()]
    if bonds:
        data["bonds"] = bonds

    return data


def _openmm_all_atomic_environments(structures, cutoff):
    "Extract all atomic environments out of a set of OpenMM structures"
    return [
        (structure_i, atom_i, cutoff)
        for structure_i, structure in enumerate(structures)
        for atom_i in range(structure.topology.getNumAtoms())
    ]
