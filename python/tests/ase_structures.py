import unittest

import ase
import numpy as np

import chemiscope


BASE_ATOMS = ase.Atoms("CO2", positions=[[0, 0, 0], [1, 1, 1], [2, 2, 5]])


class TestStructures(unittest.TestCase):
    """Conversion of structure data to chemiscope JSON"""

    def test_structures(self):
        data = chemiscope.create_input(BASE_ATOMS)
        self.assertEqual(len(data["structures"]), 1)
        self.assertEqual(data["structures"][0]["size"], 3)
        self.assertEqual(data["structures"][0]["names"], ["C", "O", "O"])
        self.assertEqual(data["structures"][0]["x"], [0, 1, 2])
        self.assertEqual(data["structures"][0]["y"], [0, 1, 2])
        self.assertEqual(data["structures"][0]["z"], [0, 1, 5])
        self.assertEqual(data["structures"][0].get("cell"), None)

        atoms = BASE_ATOMS.copy()
        atoms.cell = [23, 22, 11]
        data = chemiscope.create_input([atoms])
        self.assertEqual(len(data["structures"]), 1)
        self.assertEqual(data["structures"][0]["cell"], [23, 0, 0, 0, 22, 0, 0, 0, 11])

        atoms = BASE_ATOMS.copy()
        atoms.cell = [23, 22, 11, 120, 90, 70]
        data = chemiscope.create_input([atoms])
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
        properties = chemiscope.extract_properties(BASE_ATOMS)
        self.assertEqual(len(properties.keys()), 0)

    def test_arrays_as_atom_properties(self):
        atoms = BASE_ATOMS.copy()
        atoms.arrays["bar"] = [4, 5, 6]
        properties = chemiscope.extract_properties(atoms)

        self.assertEqual(len(properties.keys()), 1)
        self.assertEqual(properties["bar"]["target"], "atom")
        self.assertEqual(properties["bar"]["values"].tolist(), [4, 5, 6])
        self.assertEqual(properties["bar"].get("units"), None)
        self.assertEqual(properties["bar"].get("description"), None)

    def test_info_as_structure_properties(self):
        atoms_1 = BASE_ATOMS.copy()
        atoms_1.info["bar"] = 4

        atoms_2 = BASE_ATOMS.copy()
        atoms_2.info["bar"] = 6

        properties = chemiscope.extract_properties([atoms_1, atoms_2])
        self.assertEqual(len(properties.keys()), 1)
        self.assertEqual(properties["bar"]["target"], "structure")
        self.assertEqual(properties["bar"]["values"], [4, 6])
        self.assertEqual(properties["bar"].get("units"), None)
        self.assertEqual(properties["bar"].get("description"), None)

    def test_different_arrays(self):
        atoms_1 = BASE_ATOMS.copy()
        atoms_1.arrays["bar"] = [4, 5, 6]
        atoms_1.arrays["foo"] = [4, 5, 6]

        atoms_2 = BASE_ATOMS.copy()
        atoms_2.arrays["bar"] = [-1, 2, 3]
        atoms_2.arrays["baz"] = [33, 44, 55]

        with self.assertWarns(UserWarning) as cm:
            properties = chemiscope.extract_properties([atoms_1, atoms_2])

        self.assertEqual(
            cm.warning.args,
            (
                "the following atomic properties are only defined for a subset "
                "of structures: ['baz', 'foo']; they will be ignored",
            ),
        )

        self.assertEqual(len(properties.keys()), 1)
        self.assertEqual(properties["bar"]["target"], "atom")
        self.assertEqual(properties["bar"]["values"].tolist(), [4, 5, 6, -1, 2, 3])
        self.assertEqual(properties["bar"].get("units"), None)
        self.assertEqual(properties["bar"].get("description"), None)

    def test_different_info(self):
        atoms_1 = BASE_ATOMS.copy()
        atoms_1.info["bar"] = 4
        atoms_1.info["foo"] = False

        atoms_2 = BASE_ATOMS.copy()
        atoms_2.info["bar"] = -1
        atoms_2.info["baz"] = "test"

        with self.assertWarns(UserWarning) as cm:
            properties = chemiscope.extract_properties([atoms_1, atoms_2])

        self.assertEqual(
            str(cm.warning),
            "the following structure properties are only defined for a subset "
            "of structures: ['baz', 'foo']; they will be ignored",
        )

        self.assertEqual(len(properties.keys()), 1)
        self.assertEqual(properties["bar"]["target"], "structure")
        self.assertEqual(properties["bar"]["values"], [4, -1])
        self.assertEqual(properties["bar"].get("units"), None)
        self.assertEqual(properties["bar"].get("description"), None)

    def test_wrong_property_type(self):
        atoms = BASE_ATOMS.copy()
        atoms.arrays["bar"] = [{"f": 3}, {"f": 3}, {"f": 3}]

        with self.assertWarns(UserWarning) as cm:
            properties = chemiscope.extract_properties(atoms)

        self.assertEqual(
            str(cm.warning),
            "value '{'f': 3}' of type '<class 'dict'>' for the 'bar' property "
            "from ASE is not convertible to float, array or string, this property "
            "will be ignored.",
        )

        self.assertEqual(len(properties.keys()), 0)

    def test_different_lengths(self):
        atoms_1 = BASE_ATOMS.copy()
        atoms_1.info["valid_prop"] = [1, 2]
        atoms_1.info["valid_string"] = "hello"
        atoms_1.info["invalid_prop"] = [1, 2]
        atoms_1.info["non_iter_prop"] = 3

        atoms_2 = BASE_ATOMS.copy()
        atoms_2.info["valid_prop"] = [3, 4]
        # this will be turned into a string because it is a string in the other
        # structure
        atoms_2.info["valid_string"] = 1345
        atoms_2.info["invalid_prop"] = [1]
        atoms_2.info["non_iter_prop"] = 4

        with self.assertWarns(UserWarning) as cm:
            properties = chemiscope.extract_properties([atoms_1, atoms_2])

        self.assertEqual(
            str(cm.warning),
            "values of the property 'invalid_prop' have inconsistent length across "
            "different structures, it will be ignored",
        )

        self.assertEqual(
            set(properties.keys()), {"valid_prop", "valid_string", "non_iter_prop"}
        )
        self.assertEqual(properties["valid_prop"]["values"], [[1, 2], [3, 4]])
        self.assertEqual(properties["valid_string"]["values"], ["hello", "1345"])
        self.assertEqual(properties["non_iter_prop"]["values"], [3, 4])


if __name__ == "__main__":
    unittest.main()
