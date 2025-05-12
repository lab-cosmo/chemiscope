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

DEFAULT_FRAME = ase.Atoms(
    numbers=[1, 1, 1], positions=[[0, 0, 0], [1, 1, 1], [2, 2, 5]]
)
SHAPE_DEFAULTS = {
    "spheres_structure": {
        "kind": "sphere",
        "parameters": {
            "global": {"radius": 0.1},
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
    "cubes1_structure": {
        "kind": "custom",
        "parameters": {"global": {"vertices": CUBE_VERTICES}},
    },
    "cubes2_structure": {
        "kind": "custom",
        "parameters": {
            "global": {"vertices": CUBE_VERTICES, "scale": 1.0},
            "atom": [{"orientation": [1, 0, 0, 0]}, {"scale": 0.5}, {}],
        },
    },
    "cubes3_structure": {
        "kind": "custom",
        "parameters": {
            "global": {"vertices": CUBE_VERTICES, "simplices": CUBE_SIMPLICES},
        },
    },
}
INVALID_SHAPES = [
    [
        {
            "kind": "sphere",
            "parameters": {
                "global": {"diameter": 0.1},
            },
        },
        ValueError,
        "unknown shape parameter 'diameter' for 'sphere' shape kind",
    ],
    [
        {
            "kind": "sphere",
            "parameters": {
                "global": {"radius": "A"},
            },
        },
        TypeError,
        "sphere shape 'radius' must be a float, got str",
    ],
    [
        {
            "kind": "ellipsoid",
            "parameters": {
                "global": {"radius": [0.3, 0.2, 0.1]},
            },
        },
        ValueError,
        "unknown shape parameter 'radius' for 'ellipsoid' shape kind",
    ],
    [
        {
            "kind": "ellipsoid",
            "parameters": {
                "global": {"semiaxes": ["A"]},
            },
        },
        TypeError,
        "'semiaxes' must be an array with 3 values for 'ellipsoid' shape kind",
    ],
    [
        {"kind": "custom", "parameters": {"global": {"radius": 1.0}}},
        ValueError,
        "unknown shape parameter 'radius' for 'custom' shape kind",
    ],
    [
        {
            "kind": "custom",
            "parameters": {"global": {"vertices": np.array(CUBE_VERTICES)[:, :2]}},
        },
        ValueError,
        "'vertices' must be an Nx3 array values for 'custom' shape kind",
    ],
    [
        {
            "kind": "custom",
            "parameters": {"global": {"simplices": np.array(CUBE_SIMPLICES)[:, :2]}},
        },
        KeyError,
        "'vertices'",
    ],
    [
        {
            "kind": "custom",
            "parameters": {
                "global": {
                    "vertices": np.array(CUBE_VERTICES)[:, :2],
                    "simplices": np.array(CUBE_SIMPLICES)[:, :2],
                }
            },
        },
        ValueError,
        "'simplices' must be an Nx3 array values for 'custom' shape kind",
    ],
    [
        {"kind": "garbage", "parameters": {"radius": None}},
        ValueError,
        "unknown shape kind 'garbage'",
    ],
]


class TestShapes(unittest.TestCase):
    def test_custom_shapes(self):
        data = chemiscope.create_input(frames=[DEFAULT_FRAME], shapes=SHAPE_DEFAULTS)

        result = data["shapes"]
        self.assertEqual(
            list(result.keys()),
            [
                "spheres_structure",
                "ellipsoids_atoms",
                "cubes1_structure",
                "cubes2_structure",
                "cubes3_structure",
            ],
        )

        for key, values in result.items():
            self.assertEqual(SHAPE_DEFAULTS[key], values)

    def test_shape_errors(self):
        for shape, errortype, message in INVALID_SHAPES:
            with self.assertRaises(errortype) as cm:
                chemiscope.create_input(frames=[DEFAULT_FRAME], shapes={"shape": shape})
                self.assertEqual(cm.message, message)


if __name__ == "__main__":
    unittest.main()
