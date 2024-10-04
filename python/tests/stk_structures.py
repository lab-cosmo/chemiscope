import unittest

import stk

import chemiscope

BASE_FRAME = stk.BuildingBlock("N#CC")


class TestStructures(unittest.TestCase):
    """Conversion of structure data to chemiscope JSON"""

    def test_structures(self):
        data = chemiscope.create_input(BASE_FRAME)
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

        frame = BASE_FRAME.clone()
        data = chemiscope.create_input([frame])
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
            chemiscope.extract_properties(BASE_FRAME)


if __name__ == "__main__":
    unittest.main()
