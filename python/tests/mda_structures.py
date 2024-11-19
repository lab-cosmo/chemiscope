import unittest

import MDAnalysis as mda
import numpy as np

import chemiscope

BASE_FRAME = mda.Universe.empty(n_atoms=3, trajectory=True)
BASE_FRAME.add_TopologyAttr("type", ["C", "O", "O"])
BASE_FRAME.atoms.positions = np.array([[0, 0, 0], [1, 1, 1], [2, 2, 5]])


class TestStructures(unittest.TestCase):
    """Conversion of structure data to chemiscope JSON"""

    def test_structures(self):
        data = chemiscope.create_input(BASE_FRAME)
        self.assertEqual(len(data["structures"]), 1)
        self.assertEqual(data["structures"][0]["size"], 3)
        self.assertEqual(data["structures"][0]["names"], ["C", "O", "O"])
        self.assertEqual(data["structures"][0]["x"], [0, 1, 2])
        self.assertEqual(data["structures"][0]["y"], [0, 1, 2])
        self.assertEqual(data["structures"][0]["z"], [0, 1, 5])
        self.assertEqual(data["structures"][0].get("cell"), None)

        frame = BASE_FRAME.copy()
        frame.dimensions = [23, 22, 11, 90, 90, 90]
        data = chemiscope.create_input(frame)
        self.assertEqual(len(data["structures"]), 1)
        self.assertEqual(data["structures"][0]["cell"], [23, 0, 0, 0, 22, 0, 0, 0, 11])

        frame = BASE_FRAME.copy()
        frame.dimensions = [23, 22, 11, 120, 90, 70]
        data = chemiscope.create_input(frame)
        self.assertEqual(len(data["structures"]), 1)

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

    def test_exception(self):
        with self.assertRaises(RuntimeError):
            chemiscope.extract_properties(BASE_FRAME)


if __name__ == "__main__":
    unittest.main()
