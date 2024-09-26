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
        self.assertEqual(
            data["structures"][0]["x"],
            [
                1.6991195138834223,
                0.7737143493209756,
                -0.41192204250544034,
                -0.7778845126633998,
                -1.1777543806588109,
                -0.10527292738297804,
            ],
        )
        self.assertEqual(
            data["structures"][0]["y"],
            [
                -1.2265369887154756,
                -0.5721898035707434,
                0.28832060028277334,
                0.6076276888433211,
                -0.27163665176706653,
                1.1744151549238042,
            ],
        )
        self.assertEqual(
            data["structures"][0]["z"],
            [
                -0.19321573000005213,
                -0.10192268845612924,
                0.03435599430880268,
                -0.9630155400427929,
                0.6165952621860082,
                0.6072027020039786,
            ],
        )
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
        self.assertEqual(
            data["structures"][0]["x"],
            [
                1.6991195138834223,
                0.7737143493209756,
                -0.41192204250544034,
                -0.7778845126633998,
                -1.1777543806588109,
                -0.10527292738297804,
            ],
        )
        self.assertEqual(
            data["structures"][0]["y"],
            [
                -1.2265369887154756,
                -0.5721898035707434,
                0.28832060028277334,
                0.6076276888433211,
                -0.27163665176706653,
                1.1744151549238042,
            ],
        )
        self.assertEqual(
            data["structures"][0]["z"],
            [
                -0.19321573000005213,
                -0.10192268845612924,
                0.03435599430880268,
                -0.9630155400427929,
                0.6165952621860082,
                0.6072027020039786,
            ],
        )
        self.assertEqual(data["structures"][0].get("cell"), None)


class TestExtractProperties(unittest.TestCase):
    """Properties extraction"""

    def test_exception(self):
        with self.assertRaises(RuntimeError):
            chemiscope.extract_properties(BASE_FRAME)


if __name__ == "__main__":
    unittest.main()
