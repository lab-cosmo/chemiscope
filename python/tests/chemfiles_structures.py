import copy
import unittest

import numpy as np
from chemfiles import Atom, Frame, MemoryTrajectory, UnitCell

import chemiscope


def get_frame():
    frame = Frame()
    frame.add_atom(Atom("C"), [0, 0, 0])
    frame.add_atom(Atom("O"), [1, 1, 1])
    frame.add_atom(Atom("O"), [2, 2, 5])

    return frame


BASE_FRAME = get_frame()


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

        frame = copy.copy(BASE_FRAME)
        frame.cell = UnitCell([23, 22, 11])
        data = chemiscope.create_input([frame])
        self.assertEqual(len(data["structures"]), 1)
        self.assertTrue(
            np.allclose(
                data["structures"][0]["cell"],
                [23, 0, 0, 0, 22, 0, 0, 0, 11],
            )
        )

        frame = copy.copy(BASE_FRAME)
        frame.cell = UnitCell([23, 22, 11], [120, 90, 70])
        data = chemiscope.create_input([frame])
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

        # check passing a trajectory directly
        data = """3

C 0 0 0
O 1 1 1
O 2 2 5
3

C 10 10 10
O 11 11 11
O 12 12 15
"""
        t = MemoryTrajectory(format="XYZ", mode="r", data=data)

        data = chemiscope.create_input(t)
        self.assertEqual(len(data["structures"]), 2)
        self.assertEqual(data["structures"][0]["size"], 3)
        self.assertEqual(data["structures"][0]["names"], ["C", "O", "O"])

        self.assertEqual(data["structures"][1]["size"], 3)
        self.assertEqual(data["structures"][1]["names"], ["C", "O", "O"])


class TestExtractProperties(unittest.TestCase):
    """Properties extraction"""

    def test_structure_properties(self):
        frame_1 = copy.copy(BASE_FRAME)
        frame_1["number"] = 4
        frame_1["string"] = "test"
        frame_1["bool"] = True
        frame_1["vec"] = [1, 2, 3]
        frame_1["only"] = 33

        frame_2 = copy.copy(BASE_FRAME)
        frame_2["number"] = 6
        frame_2["string"] = "here"
        frame_2["bool"] = False
        frame_2["vec"] = [4, 5, 6]
        frame_2["there"] = -12

        with self.assertWarns(UserWarning) as cm:
            properties = chemiscope.extract_properties([frame_1, frame_2])

        self.assertEqual(
            cm.warning.args,
            (
                "the following structure properties are only defined for a subset "
                "of structures: ['only', 'there']; they will be ignored",
            ),
        )

        self.assertEqual(len(properties.keys()), 4)
        self.assertEqual(properties["number"]["target"], "structure")
        self.assertEqual(properties["number"]["values"].tolist(), [4, 6])

        self.assertEqual(properties["string"]["target"], "structure")
        self.assertEqual(properties["string"]["values"].tolist(), ["test", "here"])

        self.assertEqual(properties["bool"]["target"], "structure")
        self.assertEqual(properties["bool"]["values"].tolist(), [True, False])

        self.assertEqual(properties["vec"]["target"], "structure")
        self.assertIsInstance(properties["vec"]["values"], np.ndarray)
        self.assertEqual(properties["vec"]["values"].tolist(), [[1, 2, 3], [4, 5, 6]])

        # make sure the data is valid for create_input
        _ = chemiscope.create_input([frame_1, frame_2], properties=properties)

    def test_atom_properties(self):
        frame_1 = copy.copy(BASE_FRAME)
        for index, atom in enumerate(frame_1.atoms):
            atom["number"] = index + 4
            atom["string"] = "ABC"[index]
            atom["bool"] = index % 2 == 0
            atom["vec"] = [index, index + 1, index + 2]
            atom["only"] = index + 4

        frame_2 = copy.copy(BASE_FRAME)
        for index, atom in enumerate(frame_2.atoms):
            atom["number"] = index * 2 - 2
            atom["string"] = "DEF"[index]
            atom["bool"] = index % 2 == 1
            atom["vec"] = [index + 3, index + 4, index + 5]
            atom["there"] = index * 11 + 33

        with self.assertWarns(UserWarning) as cm:
            properties = chemiscope.extract_properties([frame_1, frame_2])

        self.assertEqual(
            cm.warning.args,
            (
                "the following atomic properties are only defined for a subset "
                "of structures: ['only', 'there']; they will be ignored",
            ),
        )

        self.assertEqual(len(properties.keys()), 4)
        self.assertEqual(properties["number"]["target"], "atom")
        self.assertEqual(properties["number"]["values"].tolist(), [4, 5, 6, -2, 0, 2])

        self.assertEqual(properties["string"]["target"], "atom")
        self.assertEqual(
            properties["string"]["values"].tolist(), ["A", "B", "C", "D", "E", "F"]
        )

        self.assertEqual(properties["bool"]["target"], "atom")
        self.assertEqual(
            properties["bool"]["values"].tolist(),
            [True, False, True, False, True, False],
        )

        self.assertEqual(properties["vec"]["target"], "atom")
        self.assertIsInstance(properties["vec"]["values"], np.ndarray)
        self.assertEqual(
            properties["vec"]["values"].tolist(),
            [[0, 1, 2], [1, 2, 3], [2, 3, 4], [3, 4, 5], [4, 5, 6], [5, 6, 7]],
        )

        # make sure the data is valid for create_input
        _ = chemiscope.create_input([frame_1, frame_2], properties=properties)


if __name__ == "__main__":
    unittest.main()
