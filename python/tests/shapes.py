import unittest

import ase
import numpy as np

import chemiscope

CUBE_VERTICES = [
    [0, 0, 0],
    [1, 0, 0],
    [0, 1, 0],
    [1, 1, 0],
    [0, 0, 1],
    [1, 0, 1],
    [0, 1, 1],
    [1, 1, 1],
]

CUBE_SIMPLICES = [
    [0, 1, 2],
    [1, 2, 3],
    [4, 5, 6],
    [5, 6, 7],
    [0, 1, 4],
    [0, 2, 4],
    [1, 3, 5],
    [2, 3, 6],
    [1, 4, 5],
    [2, 4, 6],
    [3, 5, 7],
    [3, 6, 7],
]

SHAPE_DEFAULTS = {
    "ellipsoid": ("semiaxes", [1, 2, 1]),
    "sphere": ("radius", 1.0),
}


class TestShapesFromASE(unittest.TestCase):
    """Conversion of shape data in ASE to chemiscope JSON"""

    def setUp(self):
        self.frame = ase.Atoms(
            numbers=np.zeros(3), positions=[[0, 0, 0], [1, 1, 1], [2, 2, 5]]
        )
        self.frame.arrays["orientation"] = [
            [1, 0, 0, 0] for _ in range(len(self.frame))
        ]

    def test_no_shape(self):
        data = chemiscope.create_input(frames=[self.frame])
        self.assertNotIn("shape", data["structures"][0])

    def test_bad_shape(self):
        frame = self.frame.copy()
        frame.info["shape"] = "invalid"

        with self.assertRaises(KeyError) as cm:
            chemiscope.extract_lammps_shapes_from_ase([frame])

        self.assertEqual(
            cm.exception.args[0],
            "The currently-supported shape in `extract_lammps_shapes_from_ase` are "
            "['ellipsoid', 'sphere'], received 'invalid'",
        )

    def test_shape_by_frame(self):
        for key in SHAPE_DEFAULTS:
            with self.subTest(shape=key):
                frame = self.frame.copy()
                frame.info["shape"] = key

                parameter, value = SHAPE_DEFAULTS[key]
                frame.info[f"shape_{parameter}"] = value
                shapes = chemiscope.extract_lammps_shapes_from_ase([frame])

                data = chemiscope.create_input(frames=[frame], shapes=shapes)
                self.assertIn("shape", data["structures"][0])

                self.assertTrue(
                    all(
                        [
                            parameter in ss
                            for s in data["structures"]
                            for ss in s["shape"]["shape"]
                        ]
                    )
                )

                self.assertTrue(
                    all(
                        [
                            "orientation" in ss
                            for s in data["structures"]
                            for ss in s["shape"]["shape"]
                        ]
                    )
                )

                self.assertTrue(
                    all(
                        [
                            len(s["shape"]["shape"]) == len(frame)
                            for s in data["structures"]
                        ]
                    )
                )

                frame.info.pop(f"shape_{parameter}")
                with self.assertRaises(KeyError) as cm:
                    chemiscope.create_input(
                        frames=[frame],
                        shapes=chemiscope.extract_lammps_shapes_from_ase([frame]),
                    )

                self.assertEquals(
                    cm.exception.args[0],
                    f"Missing required parameter 'shape_{parameter}' for '{key}' shape",
                )

    def test_by_index(self):
        for key in SHAPE_DEFAULTS:
            with self.subTest(shape=key):
                frame = self.frame.copy()
                frame.arrays["shape"] = [key] * len(frame)

                parameter, value = SHAPE_DEFAULTS[key]
                frame.arrays[f"shape_{parameter}"] = [value] * len(frame)
                shapes = chemiscope.extract_lammps_shapes_from_ase([frame])

                data = chemiscope.create_input(frames=[frame], shapes=shapes)
                self.assertIn("shape", data["structures"][0])

                self.assertTrue(
                    all(
                        [
                            parameter in ss
                            for s in data["structures"]
                            for ss in s["shape"]["shape"]
                        ]
                    )
                )

                self.assertTrue(
                    all(
                        [
                            "orientation" in ss
                            for s in data["structures"]
                            for ss in s["shape"]["shape"]
                        ]
                    )
                )

                self.assertTrue(
                    all(
                        [
                            len(s["shape"]["shape"]) == len(frame)
                            for s in data["structures"]
                        ]
                    )
                )

                frame.arrays.pop(f"shape_{parameter}")
                with self.assertRaises(KeyError) as cm:
                    chemiscope.create_input(
                        frames=[frame],
                        shapes=chemiscope.extract_lammps_shapes_from_ase([frame]),
                    )

                self.assertEquals(
                    cm.exception.args[0],
                    f"Missing required parameter 'shape_{parameter}' for '{key}' shape",
                )

    def test_no_shapes(self):
        with self.assertRaises(ValueError) as cm:
            chemiscope.extract_lammps_shapes_from_ase([self.frame])

        self.assertEqual(
            cm.exception.args[0],
            "1 frame(s) do not contain shape information",
        )


class TestShapesValidation(unittest.TestCase):
    def test_custom_shapes(self):
        pass


if __name__ == "__main__":
    unittest.main()
