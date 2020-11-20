import unittest
import ase

from chemiscope import create_input

BASE_FRAME = ase.Atoms("CO2", positions=[[0, 0, 0], [1, 1, 1], [2, 2, 5]])


class TestASEAdapter(unittest.TestCase):
    def test_structures(self):
        data = create_input(BASE_FRAME)
        self.assertEqual(len(data["structures"]), 1)
        self.assertEqual(data["structures"][0]["size"], 3)
        self.assertEqual(data["structures"][0]["names"], ["C", "O", "O"])
        self.assertEqual(data["structures"][0]["x"], [0, 1, 2])
        self.assertEqual(data["structures"][0]["y"], [0, 1, 2])
        self.assertEqual(data["structures"][0]["z"], [0, 1, 5])
        self.assertEqual(data["structures"][0].get("cell"), None)

        frame = BASE_FRAME.copy()
        frame.cell = [23, 22, 11]
        data = create_input([frame])
        self.assertEqual(len(data["structures"]), 1)
        self.assertEqual(data["structures"][0]["cell"], [23, 0, 0, 0, 22, 0, 0, 0, 11])

        frame = BASE_FRAME.copy()
        frame.cell = [23, 22, 11, 120, 90, 70]
        data = create_input([frame])
        self.assertEqual(len(data["structures"]), 1)
        self.assertEqual(
            data["structures"][0]["cell"],
            [
                23.0,
                0.0,
                0.0,
                7.5244431531647145,
                20.673237657289985,
                0.0,
                0.0,
                -5.852977748617515,
                9.313573507209156,
            ],
        )

    def test_arrays_numbers_postions_ignored(self):
        data = create_input(BASE_FRAME)
        self.assertEqual(len(BASE_FRAME.arrays.keys()), 2)
        self.assertEqual(len(data["properties"].keys()), 0)

    def test_arrays_as_atom_properties(self):
        frame = BASE_FRAME.copy()
        frame.arrays["bar"] = [4, 5, 6]
        data = create_input([frame])
        self.assertEqual(len(data["properties"].keys()), 1)
        self.assertEqual(data["properties"]["bar"]["target"], "atom")
        self.assertEqual(data["properties"]["bar"]["values"], [4, 5, 6])
        self.assertEqual(data["properties"]["bar"].get("units"), None)
        self.assertEqual(data["properties"]["bar"].get("description"), None)

    def test_info_as_frame_properties(self):
        frame = BASE_FRAME.copy()
        frame.info["bar"] = 4

        frame2 = BASE_FRAME.copy()
        frame2.info["bar"] = 6
        data = create_input([frame, frame2])
        self.assertEqual(len(data["properties"].keys()), 1)
        self.assertEqual(data["properties"]["bar"]["target"], "structure")
        self.assertEqual(data["properties"]["bar"]["values"], [4, 6])
        self.assertEqual(data["properties"]["bar"].get("units"), None)
        self.assertEqual(data["properties"]["bar"].get("description"), None)

    def test_different_arrays(self):
        frame = BASE_FRAME.copy()
        frame.arrays["bar"] = [4, 5, 6]
        frame.arrays["foo"] = [4, 5, 6]

        frame2 = BASE_FRAME.copy()
        frame2.arrays["bar"] = [-1, 2, 3]
        frame2.arrays["baz"] = [33, 44, 55]

        with self.assertWarns(UserWarning) as cm:
            data = create_input([frame, frame2])
        self.assertEqual(
            cm.warning.args,
            (
                "the following atomic properties properties are only defined "
                + "for a subset of frames: ['baz', 'foo']; they will be ignored",
            ),
        )

        self.assertEqual(len(data["properties"].keys()), 1)
        self.assertEqual(data["properties"]["bar"]["target"], "atom")
        self.assertEqual(data["properties"]["bar"]["values"], [4, 5, 6, -1, 2, 3])
        self.assertEqual(data["properties"]["bar"].get("units"), None)
        self.assertEqual(data["properties"]["bar"].get("description"), None)

    def test_different_info(self):
        frame = BASE_FRAME.copy()
        frame.info["bar"] = 4
        frame.info["foo"] = False

        frame2 = BASE_FRAME.copy()
        frame2.info["bar"] = -1
        frame2.info["baz"] = "test"

        with self.assertWarns(UserWarning) as cm:
            data = create_input([frame, frame2])
        self.assertEqual(
            cm.warning.args,
            (
                "the following structure properties properties are only defined "
                + "for a subset of frames: ['baz', 'foo']; they will be ignored",
            ),
        )

        self.assertEqual(len(data["properties"].keys()), 1)
        self.assertEqual(data["properties"]["bar"]["target"], "structure")
        self.assertEqual(data["properties"]["bar"]["values"], [4, -1])
        self.assertEqual(data["properties"]["bar"].get("units"), None)
        self.assertEqual(data["properties"]["bar"].get("description"), None)

    def test_wrong_property_type(self):
        frame = BASE_FRAME.copy()
        frame.arrays["bar"] = [{"f": 3}, {"f": 3}, {"f": 3}]

        with self.assertWarns(UserWarning) as cm:
            data = create_input(frame)
        self.assertEqual(
            cm.warning.args,
            (
                "value '{'f': 3}' of type '<class 'dict'>' for the 'bar' property "
                + "from ASE is not convertible to float or string, this property "
                + "will be ignored.",
            ),
        )

        self.assertEqual(len(data["properties"].keys()), 0)


if __name__ == "__main__":
    unittest.main()
