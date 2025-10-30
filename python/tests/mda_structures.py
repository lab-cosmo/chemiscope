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
        frame = u.copy()
        frame.dimensions = [23, 22, 11, 90, 90, 90]
        data = chemiscope.create_input(frame.atoms)
        self.assertEqual(len(data["structures"]), 2)
        self.assertEqual(data["structures"][0]["cell"], [23, 0, 0, 0, 22, 0, 0, 0, 11])

        frame = u.copy()
        frame.dimensions = [23, 22, 11, 120, 90, 70]
        data = chemiscope.create_input(frame.atoms)
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
