"""Module for handling bonding in chemiscope."""

try:
    import ase

    HAVE_ASE = True
except ImportError:
    HAVE_ASE = False


def convert_bonds_as_shapes(
    frames: list[ase.Atoms],
    bonds_per_frames: dict[int, list[tuple[int, int]]],
) -> dict[str, dict]:
    """Convert connections between atom ids in each structure to shapes.

    Parameters:

        frames:
            List of ase.Atoms object, which each are an molecule in chemiscope.

        bonds_per_frames:
            The bonds you want to add between atoms with the given ids. The key
            to the dictionary should align with the id of the molecule in
            frames.

    """
    shape_dict: dict[str, dict] = {}
    max_length = 0
    for fid, molecule in enumerate(frames):
        bonds_to_add = bonds_per_frames[fid]

        for bid, bond_info in enumerate(bonds_to_add):
            bname = f"bond_{bid}"

            # Compute the bond vector.
            position_matrix = molecule.get_positions()
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
                            "global": {"radius": 0.12, "color": "#fc5500"},
                            "structure": [],
                        },
                    }

                else:
                    num_to_add = len(shape_dict["bond_0"]["parameters"]["structure"])
                    shape_dict[bname] = {
                        "kind": "cylinder",
                        "parameters": {
                            "global": {"radius": 0.12, "color": "#fc5500"},
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
            for i in range(missing):
                fake_bond = {"vector": [0, 0, 0], "position": [0, 0, 0]}
                shape_dict[bname]["parameters"]["structure"].append(fake_bond)

    return shape_dict
