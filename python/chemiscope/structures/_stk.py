from typing import Dict, List, Tuple, Union


try:
    from stk import Molecule

    HAVE_STK = True
except ImportError:

    class Molecule:
        pass

    HAVE_STK = False


def _stk_valid_structures(
    structures: Union[Molecule, List[Molecule]],
) -> Tuple[List[Molecule], bool]:
    if HAVE_STK and isinstance(structures, Molecule):
        # deal with the user passing a single structure
        return [structures], True
    elif (
        HAVE_STK
        and isinstance(structures, list)
        and isinstance(structures[0], Molecule)
    ):
        for structure in structures:
            assert isinstance(structure, Molecule)
        return structures, True
    else:
        return structures, False


def _stk_to_json(molecule: Molecule) -> Dict[str, Union[int, list]]:
    """Implementation of structures_to_json for stk's ``Molecule``.

    The current implementation assumes no periodic information, which is safe
    for the majority of stk molecules. If necessary, we can add cell information
    in the future.

    """

    # stk specific bond order map.
    bond_map = {
        # stk int to chemiscope int.
        1: 1,
        2: 2,
        3: 3,
        # Represents dative bond types in stk.
        9: 1,
    }

    pos_mat = molecule.get_position_matrix()
    data = {}
    data["size"] = molecule.get_num_atoms()
    data["names"] = [atom.__class__.__name__ for atom in molecule.get_atoms()]
    data["x"] = [float(pos_mat[atom.get_id()][0]) for atom in molecule.get_atoms()]
    data["y"] = [float(pos_mat[atom.get_id()][1]) for atom in molecule.get_atoms()]
    data["z"] = [float(pos_mat[atom.get_id()][2]) for atom in molecule.get_atoms()]
    try:
        data["bonds"] = [
            (
                bond.get_atom1().get_id(),
                bond.get_atom2().get_id(),
                bond_map[bond.get_order()],
            )
            for bond in molecule.get_bonds()
        ]
    except KeyError as e:
        raise ValueError(f"Unknown bond order {e} (1, 2, 3, 9 (in stk) allowed)") from e

    return data


def _stk_all_atomic_environments(
    structures: List[Molecule],
    cutoff: float,
) -> List[Tuple[int, int, float]]:
    "Extract all atomic environments out of a set of stk Molecule objects"
    environments = []
    for structure_i, structure in enumerate(structures):
        for atom in structure.get_atoms():
            environments.append((structure_i, atom.get_id(), cutoff))

    return environments
