import unittest

import ase
import numpy as np

import chemiscope

BASE_FRAME = ase.Atoms("CO2", positions=[[0, 0, 0], [1, 1, 1], [2, 2, 5]])


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
        frame.cell = [23, 22, 11]
        data = chemiscope.create_input([frame])
        self.assertEqual(len(data["structures"]), 1)
        self.assertEqual(data["structures"][0]["cell"], [23, 0, 0, 0, 22, 0, 0, 0, 11])

        frame = BASE_FRAME.copy()
        frame.cell = [23, 22, 11, 120, 90, 70]
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


class TestExtractProperties(unittest.TestCase):
    """Properties extraction"""

    def test_arrays_numbers_postions_ignored(self):
        properties = chemiscope.extract_properties(BASE_FRAME)
        self.assertEqual(len(properties.keys()), 0)

    def test_arrays_as_atom_properties(self):
        frame = BASE_FRAME.copy()
        frame.arrays["bar"] = [4, 5, 6]
        properties = chemiscope.extract_properties(frame)

        self.assertEqual(len(properties.keys()), 1)
        self.assertEqual(properties["bar"]["target"], "atom")
        self.assertEqual(properties["bar"]["values"], [4, 5, 6])
        self.assertEqual(properties["bar"].get("units"), None)
        self.assertEqual(properties["bar"].get("description"), None)

    def test_info_as_frame_properties(self):
        frame_1 = BASE_FRAME.copy()
        frame_1.info["bar"] = 4

        frame_2 = BASE_FRAME.copy()
        frame_2.info["bar"] = 6

        properties = chemiscope.extract_properties([frame_1, frame_2])
        self.assertEqual(len(properties.keys()), 1)
        self.assertEqual(properties["bar"]["target"], "structure")
        self.assertEqual(properties["bar"]["values"], [4, 6])
        self.assertEqual(properties["bar"].get("units"), None)
        self.assertEqual(properties["bar"].get("description"), None)

    def test_different_arrays(self):
        frame_1 = BASE_FRAME.copy()
        frame_1.arrays["bar"] = [4, 5, 6]
        frame_1.arrays["foo"] = [4, 5, 6]

        frame_2 = BASE_FRAME.copy()
        frame_2.arrays["bar"] = [-1, 2, 3]
        frame_2.arrays["baz"] = [33, 44, 55]

        with self.assertWarns(UserWarning) as cm:
            properties = chemiscope.extract_properties([frame_1, frame_2])

        self.assertEqual(
            cm.warning.args,
            (
                "the following atomic properties properties are only defined "
                "for a subset of frames: ['baz', 'foo']; they will be ignored",
            ),
        )

        self.assertEqual(len(properties.keys()), 1)
        self.assertEqual(properties["bar"]["target"], "atom")
        self.assertEqual(properties["bar"]["values"].tolist(), [4, 5, 6, -1, 2, 3])
        self.assertEqual(properties["bar"].get("units"), None)
        self.assertEqual(properties["bar"].get("description"), None)

    def test_different_info(self):
        frame_1 = BASE_FRAME.copy()
        frame_1.info["bar"] = 4
        frame_1.info["foo"] = False

        frame_2 = BASE_FRAME.copy()
        frame_2.info["bar"] = -1
        frame_2.info["baz"] = "test"

        with self.assertWarns(UserWarning) as cm:
            properties = chemiscope.extract_properties([frame_1, frame_2])

        self.assertEqual(
            str(cm.warning),
            "the following structure properties properties are only defined "
            "for a subset of frames: ['baz', 'foo']; they will be ignored",
        )

        self.assertEqual(len(properties.keys()), 1)
        self.assertEqual(properties["bar"]["target"], "structure")
        self.assertEqual(properties["bar"]["values"], [4, -1])
        self.assertEqual(properties["bar"].get("units"), None)
        self.assertEqual(properties["bar"].get("description"), None)

    def test_wrong_property_type(self):
        frame = BASE_FRAME.copy()
        frame.arrays["bar"] = [{"f": 3}, {"f": 3}, {"f": 3}]

        with self.assertWarns(UserWarning) as cm:
            properties = chemiscope.extract_properties(frame)

        self.assertEqual(
            str(cm.warning),
            "value '{'f': 3}' of type '<class 'dict'>' for the 'bar' property "
            "from ASE is not convertible to float, array or string, this property "
            "will be ignored.",
        )

        self.assertEqual(len(properties.keys()), 0)


class TestCompositionProperties(unittest.TestCase):
    """Composition properties"""

    def test_composition(self):
        properties = chemiscope.composition_properties([BASE_FRAME, BASE_FRAME])
        self.assertEqual(len(properties.keys()), 5)

        self.assertEqual(properties["composition"]["target"], "structure")
        self.assertEqual(properties["composition"]["values"], ["CO2", "CO2"])

        self.assertEqual(properties["n_C"]["target"], "structure")
        self.assertEqual(properties["n_C"]["values"], [1, 1])

        self.assertEqual(properties["n_O"]["target"], "structure")
        self.assertEqual(properties["n_O"]["values"], [2, 2])

        self.assertEqual(properties["symbol"]["target"], "atom")
        self.assertEqual(properties["symbol"]["values"], ["C", "O", "O", "C", "O", "O"])

        self.assertEqual(properties["number"]["target"], "atom")
        self.assertEqual(properties["number"]["values"], [6, 8, 8, 6, 8, 8])


class TestEnvironments(unittest.TestCase):
    """Generate the list of environments"""

    def test_librascal_environments(self):
        frames = [ase.Atoms("CO2"), ase.Atoms("NH3")]
        for frame in frames:
            frame.arrays["atomic number"] = frame.numbers

        # center_atoms_mask is used by librascal to specify which atoms to consider
        frames[1].arrays["center_atoms_mask"] = [True, False, False, False]

        environments = chemiscope.librascal_atomic_environments(frames)
        self.assertEqual(len(environments), 4)

        properties = chemiscope.extract_properties(frames, environments=environments)
        atomic_number = properties["atomic number"]
        self.assertEqual(atomic_number["target"], "atom")
        self.assertEqual(len(atomic_number["values"]), 4)
        self.assertEqual(atomic_number["values"][0], 6)  # C in CO2
        self.assertEqual(atomic_number["values"][1], 8)  # O1 in CO2
        self.assertEqual(atomic_number["values"][2], 8)  # O2 in CO2
        self.assertEqual(atomic_number["values"][3], 7)  # N in NH3


if __name__ == "__main__":
    unittest.main()
