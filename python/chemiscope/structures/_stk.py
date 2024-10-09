from typing import Dict, List, Tuple, Union

try:
    from stk import Molecule

    HAVE_STK = True
except ImportError:

    class Molecule:
        pass

    HAVE_STK = False


def _stk_valid_structures(
    frames: Union[Molecule, List[Molecule]]
) -> Tuple[List[Molecule], bool]:
    if HAVE_STK and isinstance(frames, Molecule):
        # deal with the user passing a single frame
        return [frames], True
    elif HAVE_STK and isinstance(frames[0], Molecule):
        for frame in frames:
            assert isinstance(frame, Molecule)
        return frames, True
    else:
        return frames, False


def _stk_to_json(molecule: Molecule) -> Dict[str, Union[int, list]]:
    """Implementation of frame_to_json for stk's ``Molecule``.

    The current implementation assumes no periodic information, which is safe
    for the majority of stk molecules. If necessary, we can add cell information
    in the future.

    """
    pos_mat = molecule.get_position_matrix()
    data = {}
    data["size"] = molecule.get_num_atoms()
    data["names"] = [atom.__class__.__name__ for atom in molecule.get_atoms()]
    data["x"] = [float(pos_mat[atom.get_id()][0]) for atom in molecule.get_atoms()]
    data["y"] = [float(pos_mat[atom.get_id()][1]) for atom in molecule.get_atoms()]
    data["z"] = [float(pos_mat[atom.get_id()][2]) for atom in molecule.get_atoms()]

    return data


def _stk_all_atomic_environments(
    frames: List[Molecule],
    cutoff: float,
) -> List[Tuple[int, int, float]]:
    "Extract all atomic environments out of a set of stk Molecule objects"
    environments = []
    for structure_i, frame in enumerate(frames):
        for atom in frame.get_atoms():
            environments.append((structure_i, atom.get_id(), cutoff))

    return environments


def convert_stk_bonds_as_shapes(
    frames: List[Molecule],
    bond_color: str,
    bond_radius: float,
) -> Dict[str, Dict]:
    """Convert connections between atom ids in each structure to shapes.

    Parameters:

        frames:
            List of Molecule objects, which each are structures in
            chemiscope.

        bond_color:
            How to color the bonds added.

        bond_radius:
            Radius of bonds to add.


    """

    shape_dict: Dict[str, dict] = {}
    max_length = 0
    for molecule in frames:
        bonds_to_add = tuple(
            (bond.get_atom1().get_id(), bond.get_atom2().get_id())
            for bond in molecule.get_bonds()
        )

        for bid, bond_info in enumerate(bonds_to_add):
            bname = f"bond_{bid}"

            # Compute the bond vector.
            position_matrix = molecule.get_position_matrix()
            bond_geometry = {
                "vector": (
                    position_matrix[bond_info[1]] - position_matrix[bond_info[0]]
                ).tolist(),
                "position": (position_matrix[bond_info[0]]).tolist(),
            }

            # Add the bond name to the dictionary to be iterated through.
            if bname not in shape_dict:
                if bname == "bond_0":
                    shape_dict[bname] = {
                        "kind": "cylinder",
                        "parameters": {
                            "global": {"radius": bond_radius, "color": bond_color},
                            "structure": [],
                        },
                    }

                else:
                    num_to_add = len(shape_dict["bond_0"]["parameters"]["structure"])
                    shape_dict[bname] = {
                        "kind": "cylinder",
                        "parameters": {
                            "global": {"radius": bond_radius, "color": bond_color},
                            # Add zero placements for previously non-existant
                            # bond shapes up to the length of bond_0 -1,
                            # because that should already be at the current
                            # length that the new one should be.
                            "structure": [
                                {"vector": [0, 0, 0], "position": [0, 0, 0]}
                                for i in range(num_to_add - 1)
                            ],
                        },
                    }

            # Add vector to the shape dictionary.
            shape_dict[bname]["parameters"]["structure"].append(bond_geometry)
            max_length = max(
                (max_length, len(shape_dict[bname]["parameters"]["structure"]))
            )

        # Fill in bond shapes that are not the same length as the max length.
        for bname in shape_dict.keys():
            missing = max_length - len(shape_dict[bname]["parameters"]["structure"])
            if missing == 0:
                continue
            for _ in range(missing):
                fake_bond = {"vector": [0, 0, 0], "position": [0, 0, 0]}
                shape_dict[bname]["parameters"]["structure"].append(fake_bond)

    return shape_dict


def _stk_list_atom_properties(frames: List[Molecule]) -> list:
    # stk cannot have atom properties or structure properties, so skipping.
    return []


def _stk_list_structure_properties(frames: List[Molecule]) -> list:
    # stk cannot have atom properties or structure properties, so skipping.
    return []
