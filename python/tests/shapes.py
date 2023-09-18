import unittest

import ase

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


class TestShapes(unittest.TestCase):
    def test_custom_shapes(self):
        frame = ase.Atoms(
            numbers=[1, 1, 1], positions=[[0, 0, 0], [1, 1, 1], [2, 2, 5]]
        )

        # TODO re-enable after having fixed custom shapes
        """ 
            "cubes": [
                [
                    {"kind": "custom", "vertices": CUBE_VERTICES},
                    {
                        "kind": "custom",
                        "vertices": CUBE_VERTICES,
                        "orientation": [1, 0, 0, 0],
                    },
                    {
                        "kind": "custom",
                        "vertices": CUBE_VERTICES,
                        "simplices": CUBE_SIMPLICES,
                    },
                ],
            ],"""
        shapes = {
            "spheres_structure": {
                "kind": "sphere",
                "parameters": {
                    "global": {"radius": 0.1},
                    "structure": [{"position": [1, 2, 3]}],
                },
            },
            "ellipsoids_atoms": {
                "kind": "ellipsoid",
                "parameters": {
                    "global": {"semiaxes": [0.3, 0.2, 0.1]},
                    "atom": [
                        {},
                        {"semiaxes": [0.1, 0.2, 0.3]},
                        {"orientation": [0.2, 0.3, 0.4, 1]},
                    ],
                },
            },
        }

        data = chemiscope.create_input(frames=[frame], shapes=shapes)

        result = data["shapes"]
        print(result)
        self.assertEqual(list(result.keys()), ["spheres_structure", "ellipsoids_atoms"])

        for key, values in result.items():
            self.assertEqual(shapes[key], values)


class TestShapesFromASE(unittest.TestCase):
    """Conversion of shape data in ASE to chemiscope JSON"""

    def setUp(self):
        self.frame = ase.Atoms(
            numbers=[1, 1, 1], positions=[[0, 0, 0], [1, 1, 1], [2, 2, 5]]
        )
        self.frame.arrays["orientation"] = [
            [1, 0, 0, 0] for _ in range(len(self.frame))
        ]

    def test_no_shape(self):
        data = chemiscope.create_input(frames=[self.frame])
        self.assertNotIn("shapes", data["structures"][0])

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
        for shape_kind in SHAPE_DEFAULTS:
            with self.subTest(shape=shape_kind):
                frame = self.frame.copy()
                shape_name = "shape_name"
                frame.info[shape_name] = shape_kind

                parameter, value = SHAPE_DEFAULTS[shape_kind]
                frame.info[f"{shape_name}_{parameter}"] = value
                shapes = chemiscope.extract_lammps_shapes_from_ase(
                    [frame], key=shape_name
                )

                data = chemiscope.create_input(frames=[frame], shapes=shapes)
                self.assertIn("shapes", data["structures"][0])

                self.assertTrue(
                    all(
                        [
                            parameter in ss
                            for s in data["structures"]
                            for ss in s["shapes"][shape_name]
                        ]
                    )
                )

                self.assertTrue(
                    all(
                        [
                            "orientation" in ss
                            for s in data["structures"]
                            for ss in s["shapes"][shape_name]
                        ]
                    )
                )

                self.assertTrue(
                    all(
                        [
                            len(s["shapes"][shape_name]) == len(frame)
                            for s in data["structures"]
                        ]
                    )
                )

                frame.info.pop(f"{shape_name}_{parameter}")
                with self.assertRaises(KeyError) as cm:
                    chemiscope.create_input(
                        frames=[frame],
                        shapes=chemiscope.extract_lammps_shapes_from_ase(
                            [frame], key=shape_name
                        ),
                    )

                self.assertEquals(
                    cm.exception.args[0],
                    f"Missing required parameter '{shape_name}_{parameter}' for "
                    f"'{shape_kind}' shape",
                )

    def test_by_index(self):
        for shape_kind in SHAPE_DEFAULTS:
            with self.subTest(shape=shape_kind):
                frame = self.frame.copy()
                shape_name = "shape_name"
                frame.arrays[shape_name] = [shape_kind] * len(frame)

                parameter, value = SHAPE_DEFAULTS[shape_kind]
                frame.arrays[f"{shape_name}_{parameter}"] = [value] * len(frame)
                shapes = chemiscope.extract_lammps_shapes_from_ase(
                    [frame], key=shape_name
                )

                data = chemiscope.create_input(frames=[frame], shapes=shapes)
                self.assertIn("shapes", data["structures"][0])

                self.assertTrue(
                    all(
                        [
                            parameter in ss
                            for s in data["structures"]
                            for ss in s["shapes"][shape_name]
                        ]
                    )
                )

                self.assertTrue(
                    all(
                        [
                            "orientation" in ss
                            for s in data["structures"]
                            for ss in s["shapes"][shape_name]
                        ]
                    )
                )

                self.assertTrue(
                    all(
                        [
                            len(s["shapes"][shape_name]) == len(frame)
                            for s in data["structures"]
                        ]
                    )
                )

                frame.arrays.pop(f"{shape_name}_{parameter}")
                with self.assertRaises(KeyError) as cm:
                    chemiscope.create_input(
                        frames=[frame],
                        shapes=chemiscope.extract_lammps_shapes_from_ase(
                            [frame], key=shape_name
                        ),
                    )

                self.assertEquals(
                    cm.exception.args[0],
                    f"Missing required parameter '{shape_name}_{parameter}' for "
                    f"'{shape_kind}' shape",
                )

    def test_no_shapes(self):
        with self.assertRaises(ValueError) as cm:
            chemiscope.extract_lammps_shapes_from_ase([self.frame])

        self.assertEqual(
            cm.exception.args[0],
            "1 frame(s) do not contain shape information",
        )


if __name__ == "__main__":
    unittest.main()
