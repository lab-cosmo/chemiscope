import unittest

import MDAnalysis as mda
import numpy as np
from MDAnalysis.coordinates.memory import MemoryReader

import chemiscope


POSITIONS = np.array(
    [[[0, 0, 0], [1, 1, 1], [2, 2, 5]], [[0.1, 0, 0], [1, 1.1, 1], [2, 2, 5.1]]]
)
VELOCITIES = np.array(
    [[[0, 0, 1], [1, 1, 1], [2, 2, 5]], [[0.1, 0, 2], [1, 1.1, 1], [2, 2, 5.1]]]
)
FORCES = np.array(
    [[[0, 0, 3], [1, 1, 1], [2, 2, 5]], [[0.1, 0, 4], [1, 1.1, 1], [2, 2, 5.1]]]
)


def create_test_universe(has_velocities=True, has_forces=True):
    u = mda.Universe.empty(
        n_atoms=3, trajectory=True, velocities=has_velocities, forces=has_forces
    )
    u.add_TopologyAttr("type", ["C", "O", "O"])
    properties = {
        "velocities": VELOCITIES,
        "forces": FORCES,
    }
    if not has_forces:
        del properties["forces"]

    if not has_velocities:
        del properties["velocities"]

    u.trajectory = MemoryReader(
        coordinate_array=POSITIONS,
        order="fac",
        n_atoms=3,
        n_frames=2,
        **properties,
    )

    return u


class TestStructures(unittest.TestCase):
    """Conversion of structure data to chemiscope JSON"""

    def test_structures(self):
        u = create_test_universe()
        data = chemiscope.create_input(u.atoms)
        self.assertTrue(np.allclose(len(data["structures"]), 2))
        self.assertTrue(np.allclose(data["structures"][0]["size"], 3))
        self.assertEqual(data["structures"][0]["names"], ["C", "O", "O"])
        self.assertTrue(np.allclose(data["structures"][0]["x"], [0, 1, 2]))
        self.assertTrue(np.allclose(data["structures"][0]["y"], [0, 1, 2]))
        self.assertTrue(np.allclose(data["structures"][0]["z"], [0, 1, 5]))
        self.assertTrue(np.allclose(data["structures"][1]["x"], [0.1, 1, 2]))
        self.assertTrue(np.allclose(data["structures"][1]["y"], [0, 1.1, 2]))
        self.assertTrue(np.allclose(data["structures"][1]["z"], [0, 1, 5.1]))
        self.assertEqual(data["structures"][0].get("cell"), None)
        structure = u.copy()
        structure.dimensions = [23, 22, 11, 90, 90, 90]
        data = chemiscope.create_input(structure.atoms)
        self.assertEqual(len(data["structures"]), 2)
        self.assertEqual(data["structures"][0]["cell"], [23, 0, 0, 0, 22, 0, 0, 0, 11])

        structure = u.copy()
        structure.dimensions = [23, 22, 11, 120, 90, 70]
        data = chemiscope.create_input(structure.atoms)
        self.assertEqual(len(data["structures"]), 2)

        cell = [
            23.0,
            0.0,
            0.0,
            7.5244431531647145,
            20.673237657289985,
            0.0,
            0.0,
            -5.852977748617515,
            9.313573507209156,
        ]
        self.assertTrue(np.allclose(data["structures"][0]["cell"], cell))

    def test_nucleotide_residues(self):
        # MDAnalysis' `nucleic` selection misses some force-field names
        names = ["G5", "C3", "RA5", "DA", "NAG", "HOH"]
        u = mda.Universe.empty(
            n_atoms=6, n_residues=6, atom_resindex=range(6), trajectory=True
        )
        u.add_TopologyAttr("type", ["C"] * 6)
        u.add_TopologyAttr("resname", names)
        u.add_TopologyAttr("resid", range(1, 7))
        u.add_TopologyAttr("chainID", ["A"] * 6)

        structure = chemiscope.create_input(u.atoms)["structures"][0]
        self.assertEqual(structure["resnames"], names)
        self.assertEqual(structure["hetatom"], [False, False, False, False, True, True])

    def test_bonds(self):
        # four atoms per residue, bonded in a chain inside each residue
        u = mda.Universe.empty(
            n_atoms=8, n_residues=2, atom_resindex=[0] * 4 + [1] * 4, trajectory=True
        )
        u.add_TopologyAttr("type", ["C"] * 8)
        u.add_TopologyAttr("bonds", [(0, 1), (1, 2), (2, 3), (4, 5), (5, 6), (6, 7)])

        structure = chemiscope.create_input(u.atoms)["structures"][0]
        self.assertEqual(
            structure["bonds"],
            [[0, 1, 1], [1, 2, 1], [2, 3, 1], [4, 5, 1], [5, 6, 1], [6, 7, 1]],
        )

        # a selection is numbered in the universe, the bonds have to follow it
        structure = chemiscope.create_input(u.atoms[4:])["structures"][0]
        self.assertEqual(structure["bonds"], [[0, 1, 1], [1, 2, 1], [2, 3, 1]])

    def test_partial_bonds_are_left_to_the_viewer(self):
        # a PDB with CONECT records for one residue only ends up like this
        u = mda.Universe.empty(
            n_atoms=8, n_residues=2, atom_resindex=[0] * 4 + [1] * 4, trajectory=True
        )
        u.add_TopologyAttr("type", ["C"] * 8)
        u.add_TopologyAttr("bonds", [(0, 1), (1, 2), (2, 3)])

        structure = chemiscope.create_input(u.atoms)["structures"][0]
        self.assertNotIn("bonds", structure)

        # a lone atom, an ion or a water oxygen, does not make them partial
        u = mda.Universe.empty(
            n_atoms=5, n_residues=2, atom_resindex=[0] * 4 + [1], trajectory=True
        )
        u.add_TopologyAttr("type", ["C"] * 4 + ["Na"])
        u.add_TopologyAttr("bonds", [(0, 1), (1, 2), (2, 3)])

        structure = chemiscope.create_input(u.atoms)["structures"][0]
        self.assertEqual(structure["bonds"], [[0, 1, 1], [1, 2, 1], [2, 3, 1]])


class TestExtractProperties(unittest.TestCase):
    """Properties extraction"""

    def test_properties(self):
        u = create_test_universe()
        properties = chemiscope.extract_properties(u.atoms)
        self.assertTrue(np.allclose(properties["velocities"]["values"], VELOCITIES))
        self.assertTrue(np.allclose(properties["forces"]["values"], FORCES))

    def test_missing_properties(self):
        u = create_test_universe(has_forces=False)
        properties = chemiscope.extract_properties(u.atoms)
        self.assertFalse("forces" in properties)


if __name__ == "__main__":
    unittest.main()
