import unittest

import stk

import chemiscope


BASE_STRUCTURE = stk.BuildingBlock("N#CC")


class TestStructures(unittest.TestCase):
    """Conversion of structure data to chemiscope JSON"""

    def test_structures(self):
        data = chemiscope.create_input(BASE_STRUCTURE)
        self.assertEqual(len(data["structures"]), 1)
        self.assertEqual(data["structures"][0]["size"], 6)
        self.assertEqual(
            data["structures"][0]["names"],
            ["N", "C", "C", "H", "H", "H"],
        )
        self.assertEqual(len(data["structures"][0]["x"]), 6)
        self.assertEqual(len(data["structures"][0]["y"]), 6)
        self.assertEqual(len(data["structures"][0]["z"]), 6)
        self.assertEqual(data["structures"][0].get("cell"), None)

        # Not testing cell because stk implementation does not have that yet.

        structure = BASE_STRUCTURE.clone()
        data = chemiscope.create_input([structure])
        self.assertEqual(len(data["structures"]), 1)
        self.assertEqual(data["structures"][0]["size"], 6)
        self.assertEqual(
            data["structures"][0]["names"],
            ["N", "C", "C", "H", "H", "H"],
        )
        self.assertEqual(len(data["structures"][0]["x"]), 6)
        self.assertEqual(len(data["structures"][0]["y"]), 6)
        self.assertEqual(len(data["structures"][0]["z"]), 6)
        self.assertEqual(data["structures"][0].get("cell"), None)


class TestExtractProperties(unittest.TestCase):
    """Properties extraction"""

    def test_exception(self):
        with self.assertRaises(RuntimeError):
            chemiscope.extract_properties(BASE_STRUCTURE)


if __name__ == "__main__":
    unittest.main()
