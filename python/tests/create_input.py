import unittest
import numpy as np
import copy
import ase

from chemiscope import create_input
from chemiscope import all_atomic_environments, librascal_atomic_environments

TEST_FRAMES = [ase.Atoms("CO2")]


class TestCreateInputMeta(unittest.TestCase):
    def test_meta(self):
        meta = {}
        data = create_input(frames=TEST_FRAMES, meta=meta)
        self.assertEqual(data["meta"]["name"], "<unknown>")
        self.assertEqual(len(data["meta"].keys()), 1)

        meta = {"name": ""}
        data = create_input(frames=TEST_FRAMES, meta=meta)
        self.assertEqual(data["meta"]["name"], "<unknown>")
        self.assertEqual(len(data["meta"].keys()), 1)

        meta = {"name": "foo"}
        data = create_input(frames=TEST_FRAMES, meta=meta)
        self.assertEqual(data["meta"]["name"], "foo")
        self.assertEqual(len(data["meta"].keys()), 1)

        meta = {"name": "foo", "description": "bar"}
        data = create_input(frames=TEST_FRAMES, meta=meta)
        self.assertEqual(data["meta"]["name"], "foo")
        self.assertEqual(data["meta"]["description"], "bar")
        self.assertEqual(len(data["meta"].keys()), 2)

        meta = {"name": "foo", "references": ["bar"]}
        data = create_input(frames=TEST_FRAMES, meta=meta)
        self.assertEqual(data["meta"]["name"], "foo")
        self.assertEqual(len(data["meta"]["references"]), 1)
        self.assertEqual(data["meta"]["references"][0], "bar")
        self.assertEqual(len(data["meta"].keys()), 2)

        meta = {"name": "foo", "authors": ["bar"]}
        data = create_input(frames=TEST_FRAMES, meta=meta)
        self.assertEqual(data["meta"]["name"], "foo")
        self.assertEqual(len(data["meta"]["authors"]), 1)
        self.assertEqual(data["meta"]["authors"][0], "bar")
        self.assertEqual(len(data["meta"].keys()), 2)

    def test_meta_unknown_keys_warning(self):
        meta = {"name": "foo", "what_is_this": "I don't know"}
        with self.assertWarns(UserWarning) as cm:
            data = create_input(frames=TEST_FRAMES, meta=meta)

        self.assertEqual(data["meta"]["name"], "foo")
        self.assertEqual(len(data["meta"].keys()), 1)

        self.assertEqual(str(cm.warning), "ignoring unexpected metadata: what_is_this")

    def test_meta_conversions(self):
        meta = {"name": 33}
        data = create_input(frames=TEST_FRAMES, meta=meta)
        self.assertEqual(data["meta"]["name"], "33")
        self.assertEqual(len(data["meta"].keys()), 1)

        meta = {"name": ["foo", "bar"], "description": False}
        data = create_input(frames=TEST_FRAMES, meta=meta)
        self.assertEqual(data["meta"]["name"], "['foo', 'bar']")
        self.assertEqual(data["meta"]["description"], "False")
        self.assertEqual(len(data["meta"].keys()), 2)

        meta = {"name": "foo", "references": (3, False)}
        data = create_input(frames=TEST_FRAMES, meta=meta)
        self.assertEqual(data["meta"]["name"], "foo")
        self.assertEqual(len(data["meta"]["references"]), 2)
        self.assertEqual(data["meta"]["references"][0], "3")
        self.assertEqual(data["meta"]["references"][1], "False")
        self.assertEqual(len(data["meta"].keys()), 2)

        meta = {"name": "foo", "authors": (3, False)}
        data = create_input(frames=TEST_FRAMES, meta=meta)
        self.assertEqual(data["meta"]["name"], "foo")
        self.assertEqual(len(data["meta"]["authors"]), 2)
        self.assertEqual(data["meta"]["authors"][0], "3")
        self.assertEqual(data["meta"]["authors"][1], "False")
        self.assertEqual(len(data["meta"].keys()), 2)


class TestCreateInputProperties(unittest.TestCase):
    def test_properties(self):
        # values are numbers
        properties = {"name": {"target": "atom", "values": [2, 3, 4]}}
        data = create_input(frames=TEST_FRAMES, properties=properties)
        self.assertEqual(data["properties"]["name"]["target"], "atom")
        self.assertEqual(data["properties"]["name"]["values"], [2, 3, 4])
        self.assertEqual(len(data["properties"]["name"].keys()), 2)

        # values are strings
        properties = {"name": {"target": "atom", "values": ["2", "3", "4"]}}
        data = create_input(frames=TEST_FRAMES, properties=properties)
        self.assertEqual(data["properties"]["name"]["target"], "atom")
        self.assertEqual(data["properties"]["name"]["values"], ["2", "3", "4"])
        self.assertEqual(len(data["properties"]["name"].keys()), 2)

        properties = {
            "name": {
                "target": "atom",
                "values": [2, 3, 4],
                "description": "foo",
            },
        }
        data = create_input(frames=TEST_FRAMES, properties=properties)
        self.assertEqual(data["properties"]["name"]["target"], "atom")
        self.assertEqual(data["properties"]["name"]["description"], "foo")
        self.assertEqual(data["properties"]["name"]["values"], [2, 3, 4])
        self.assertEqual(len(data["properties"]["name"].keys()), 3)

        properties = {
            "name": {
                "target": "atom",
                "values": [2, 3, 4],
                "units": "foo",
            },
        }
        data = create_input(frames=TEST_FRAMES, properties=properties)
        self.assertEqual(data["properties"]["name"]["target"], "atom")
        self.assertEqual(data["properties"]["name"]["units"], "foo")
        self.assertEqual(data["properties"]["name"]["values"], [2, 3, 4])
        self.assertEqual(len(data["properties"]["name"].keys()), 3)

    def test_ndarray_properties(self):
        # shape N
        properties = {"name": {"target": "atom", "values": np.array([2, 3, 4])}}
        data = create_input(frames=TEST_FRAMES, properties=properties)
        self.assertEqual(data["properties"]["name"]["target"], "atom")
        self.assertEqual(data["properties"]["name"]["values"], [2, 3, 4])
        self.assertEqual(len(data["properties"].keys()), 1)

        # shape N
        properties = {"name": {"target": "atom", "values": np.array(["2", "3", "4"])}}
        data = create_input(frames=TEST_FRAMES, properties=properties)
        self.assertEqual(data["properties"]["name"]["target"], "atom")
        self.assertEqual(data["properties"]["name"]["values"], ["2", "3", "4"])
        self.assertEqual(len(data["properties"].keys()), 1)

        # shape N x 1
        properties = {"name": {"target": "atom", "values": np.array([[2], [3], [4]])}}
        data = create_input(frames=TEST_FRAMES, properties=properties)
        self.assertEqual(data["properties"]["name"]["target"], "atom")
        self.assertEqual(data["properties"]["name"]["values"], [2, 3, 4])
        self.assertEqual(len(data["properties"].keys()), 1)

        # shape N x 3
        properties = {
            "name": {
                "target": "atom",
                "values": np.array([[1, 2, 4], [1, 2, 4], [1, 2, 4]]),
            }
        }
        data = create_input(frames=TEST_FRAMES, properties=properties)
        self.assertEqual(data["properties"]["name[1]"]["target"], "atom")
        self.assertEqual(data["properties"]["name[1]"]["values"], [1, 1, 1])
        self.assertEqual(data["properties"]["name[2]"]["target"], "atom")
        self.assertEqual(data["properties"]["name[2]"]["values"], [2, 2, 2])
        self.assertEqual(data["properties"]["name[3]"]["target"], "atom")
        self.assertEqual(data["properties"]["name[3]"]["values"], [4, 4, 4])
        self.assertEqual(len(data["properties"].keys()), 3)

    def test_shortened_properties(self):
        # atom property
        properties = {"name": [2, 3, 4]}
        data = create_input(frames=TEST_FRAMES, properties=properties)
        self.assertEqual(data["properties"]["name"]["target"], "atom")
        self.assertEqual(data["properties"]["name"]["values"], [2, 3, 4])
        self.assertEqual(len(data["properties"]["name"].keys()), 2)

        # frame property
        properties = {"name": [2]}
        data = create_input(frames=TEST_FRAMES, properties=properties)
        self.assertEqual(data["properties"]["name"]["target"], "structure")
        self.assertEqual(data["properties"]["name"]["values"], [2])
        self.assertEqual(len(data["properties"]["name"].keys()), 2)

        # ndarray frame property
        properties = {"name": np.array([[2, 4]])}
        data = create_input(frames=TEST_FRAMES, properties=properties)
        self.assertEqual(data["properties"]["name[1]"]["target"], "structure")
        self.assertEqual(data["properties"]["name[1]"]["values"], [2])
        self.assertEqual(len(data["properties"]["name[1]"].keys()), 2)

        self.assertEqual(data["properties"]["name[2]"]["target"], "structure")
        self.assertEqual(data["properties"]["name[2]"]["values"], [4])
        self.assertEqual(len(data["properties"]["name[2]"].keys()), 2)

        # the initial properties object must not be changed
        self.assertEqual(type(properties["name"]), np.ndarray)

    def test_shortened_properties_errors(self):
        properties = {"name": ["2", "3"]}
        with self.assertRaises(ValueError) as cm:
            create_input(frames=TEST_FRAMES, properties=properties)
        self.assertEqual(
            str(cm.exception),
            "The length of property values is different from the number of "
            + "structures and the number of atoms, we can not guess the target. "
            + "Got n_atoms = 3, n_structures = 1, the length of property values "
            + "is 2, for the 'name' property",
        )

        properties = {"name": ase.Atoms("CO2")}
        with self.assertRaises(ValueError) as cm:
            create_input(frames=TEST_FRAMES, properties=properties)
        self.assertEqual(
            str(cm.exception),
            "Property values should be either list or numpy array, got "
            + "<class 'ase.atoms.Atoms'> instead",
        )

        properties = {"name": ["2", "3"]}
        frames_single_atoms = [ase.Atoms("C"), ase.Atoms("H")]
        with self.assertWarns(UserWarning) as cm:
            data = create_input(frames=frames_single_atoms, properties=properties)

        self.assertEqual(data["properties"]["name"]["target"], "structure")

        self.assertEqual(
            cm.warning.args[0],
            "The target of the property 'name' is ambiguous because there is the same "
            + "number of atoms and structures. Will assume target=structure. ",
        )

    def test_invalid_name(self):
        properties = {"": {"target": "atom", "values": [2, 3, 4]}}
        with self.assertRaises(Exception) as cm:
            create_input(frames=TEST_FRAMES, properties=properties)
        self.assertEqual(
            str(cm.exception), "the name of a property can not be the empty string"
        )

        properties = {False: {"target": "atom", "values": [2, 3, 4]}}
        with self.assertRaises(Exception) as cm:
            create_input(frames=TEST_FRAMES, properties=properties)
        self.assertEqual(
            str(cm.exception),
            "the name of a property name must be a string, "
            + "got 'False' of type <class 'bool'>",
        )

    def test_invalid_target(self):
        properties = {"name": {"values": [2, 3, 4]}}
        with self.assertRaises(Exception) as cm:
            create_input(frames=TEST_FRAMES, properties=properties)
        self.assertEqual(str(cm.exception), "missing 'target' for the 'name' property")

        properties = {"name": {"target": "atoms", "values": [2, 3, 4]}}
        with self.assertRaises(Exception) as cm:
            create_input(frames=TEST_FRAMES, properties=properties)
        self.assertEqual(
            str(cm.exception),
            "the target must be 'atom' or 'structure' for the 'name' property",
        )

    def test_invalid_types_metadata(self):
        properties = {"name": {"target": "atom", "values": [2, 3, 4], "units": False}}
        data = create_input(frames=TEST_FRAMES, properties=properties)
        self.assertEqual(data["properties"]["name"]["units"], "False")

        properties = {
            "name": {"target": "atom", "values": [2, 3, 4], "description": False}
        }
        data = create_input(frames=TEST_FRAMES, properties=properties)
        self.assertEqual(data["properties"]["name"]["description"], "False")

    def test_property_unknown_keys_warning(self):
        properties = {"name": {"target": "atom", "values": [2, 3, 4], "what": False}}
        with self.assertWarns(UserWarning) as cm:
            create_input(frames=TEST_FRAMES, properties=properties)
        self.assertEqual(str(cm.warning), "ignoring unexpected property key: what")

    def test_invalid_values_types(self):
        properties = {"name": {"target": "atom", "values": 3}}
        with self.assertRaises(Exception) as cm:
            create_input(frames=TEST_FRAMES, properties=properties)
        self.assertEqual(
            str(cm.exception), "unknown type (<class 'int'>) for property 'name'"
        )

        properties = {"name": {"target": "atom", "values": {"test": "bad"}}}
        with self.assertRaises(Exception) as cm:
            create_input(frames=TEST_FRAMES, properties=properties)
        self.assertEqual(
            str(cm.exception), "unknown type (<class 'dict'>) for property 'name'"
        )

        properties = {"name": {"target": "atom", "values": [{}, {}, {}]}}
        with self.assertRaises(Exception) as cm:
            create_input(frames=TEST_FRAMES, properties=properties)
        self.assertEqual(
            str(cm.exception),
            "unsupported type in property 'name' values: should be string or number",
        )

    def test_wrong_number_of_values(self):
        properties = {"name": {"target": "atom", "values": [2, 3]}}
        environments = [(0, 0, 3), (0, 1, 3), (0, 2, 3)]
        with self.assertRaises(Exception) as cm:
            create_input(
                frames=TEST_FRAMES, properties=properties, environments=environments
            )
        self.assertEqual(
            str(cm.exception),
            "wrong size for the property 'name' with target=='atom': "
            + "expected 3 values, got 2",
        )

        properties = {"name": {"target": "structure", "values": [2, 3, 5]}}
        with self.assertRaises(Exception) as cm:
            create_input(frames=TEST_FRAMES, properties=properties)
        self.assertEqual(
            str(cm.exception),
            "wrong size for the property 'name' with target=='structure': "
            + "expected 1 values, got 3",
        )

    def test_properties_also_in_frame(self):
        properties = {"name": {"target": "atom", "values": [2, 3, 4]}}
        frames = copy.deepcopy(TEST_FRAMES)
        for frame in frames:
            frame.info["name"] = "test"

        with self.assertWarns(UserWarning) as cm:
            data = create_input(frames=frames, properties=properties)

        self.assertEqual(
            str(cm.warning),
            "ignoring the 'name' structure property coming from the "
            + "structures since it is already part of the properties",
        )
        self.assertEqual(data["properties"]["name"]["target"], "atom")

        properties = {"name": {"target": "structure", "values": [2]}}
        frames = copy.deepcopy(TEST_FRAMES)
        for frame in frames:
            frame.arrays["name"] = [0] * len(frame)

        with self.assertWarns(UserWarning) as cm:
            data = create_input(frames=frames, properties=properties)

        self.assertEqual(
            str(cm.warning),
            "ignoring the 'name' atom property coming from the "
            + "structures since it is already part of the properties",
        )
        self.assertEqual(data["properties"]["name"]["target"], "structure")

    def test_property_only(self):
        properties = {"name": [2, 3, 4]}
        data = create_input(properties=properties)
        self.assertEqual(data["properties"]["name"]["target"], "structure")

        properties = {
            "first": [2, 3, 4],
            "second": np.array([[1, 2], [1, 2], [1, 2]]),
        }
        data = create_input(properties=properties)
        self.assertEqual(data["properties"]["second[1]"]["target"], "structure")
        self.assertEqual(data["properties"]["second[2]"]["target"], "structure")

        # error: different size
        properties = {"first": [2, 3, 4], "second": [2, 3, 4, 5]}
        with self.assertRaises(Exception) as cm:
            data = create_input(properties=properties)

        self.assertEqual(
            str(cm.exception),
            "wrong size for property 'second': expected 3 elements, but got an array with 4 entries",
        )

        # error: target is not "structure"
        properties = {"name": {"target": "atom", "values": [2]}}
        with self.assertRaises(Exception) as cm:
            data = create_input(properties=properties)

        self.assertEqual(
            str(cm.exception),
            "property 'name' has a non-structure target, which is not allowed if frames are not provided",
        )


class TestCreateInputEnvironments(unittest.TestCase):
    def test_manual_environments_list(self):
        environments = [
            (0, 0, 3.5),
            (1, 1, 2.5),
            (1, 2, 3),
        ]
        data = create_input(frames=TEST_FRAMES + TEST_FRAMES, environments=environments)
        self.assertEqual(len(data["environments"]), 3)

        for i, env in enumerate(data["environments"]):
            self.assertEqual(env["structure"], environments[i][0])
            self.assertEqual(env["center"], environments[i][1])
            self.assertEqual(env["cutoff"], environments[i][2])

    def test_all_environments(self):
        environments = all_atomic_environments(TEST_FRAMES, cutoff=6)
        for i, (structure, center, cutoff) in enumerate(environments):
            self.assertEqual(structure, 0)
            self.assertEqual(center, i)
            self.assertEqual(cutoff, 6)

    def test_librascal_environments(self):
        frames = [ase.Atoms("CO2"), ase.Atoms("NH3")]
        for frame in frames:
            frame.arrays["atomic number"] = frame.numbers

        # center_atoms_mask is used by librascal to specify which atoms to consider
        frames[1].arrays["center_atoms_mask"] = [True, False, False, False]

        environments = librascal_atomic_environments(frames)
        data = create_input(frames=frames, environments=environments)

        self.assertEqual(len(data["environments"]), 4)

        atomic_number = data["properties"]["atomic number"]
        self.assertEqual(atomic_number["target"], "atom")
        self.assertEqual(len(atomic_number["values"]), 4)
        self.assertEqual(atomic_number["values"][0], 6)  # C in CO2
        self.assertEqual(atomic_number["values"][1], 8)  # O1 in CO2
        self.assertEqual(atomic_number["values"][2], 8)  # O2 in CO2
        self.assertEqual(atomic_number["values"][3], 7)  # N in NH3


if __name__ == "__main__":
    unittest.main()
