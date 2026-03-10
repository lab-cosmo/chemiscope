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

DEFAULT_STRUCTURE = ase.Atoms(
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
    # cylinders: global scalar broadcast for radii
    "cylinders_global": {
        "kind": "cylinders",
        "parameters": {
            "global": {
                "vectors": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                "radii": 0.05,
            },
        },
    },
    # cylinders: per-atom with bases, per-element radii and colors
    "cylinders_atoms": {
        "kind": "cylinders",
        "parameters": {
            "global": {
                "vectors": [[1, 0, 0], [0, 1, 0]],
                "bases": [[0, 0, 0], [1, 0, 0]],
                "radii": [0.05, 0.1],
                "colors": [0xFF0000, 0x00FF00],
            },
            "atom": [{}, {"scale": 2.0}, {}],
        },
    },
    # cylinders: structure-level with scale
    "cylinders_structure": {
        "kind": "cylinders",
        "parameters": {
            "global": {
                "vectors": [[1, 0, 0]],
                "radii": 0.1,
            },
            "structure": [{"position": [0, 0, 0], "scale": 1.5}],
        },
    },
    # spheres: global scalar broadcast for radii
    "spheres_multi_global": {
        "kind": "spheres",
        "parameters": {
            "global": {
                "centers": [[0, 0, 0], [1, 1, 1], [2, 0, 0]],
                "radii": 0.2,
            },
        },
    },
    # spheres: per-atom with per-element radii and colors
    "spheres_multi_atoms": {
        "kind": "spheres",
        "parameters": {
            "global": {
                "centers": [[0, 0, 0], [1, 0, 0]],
                "radii": [0.1, 0.2],
                "colors": [0xFF0000, 0x0000FF],
            },
            "atom": [{}, {"scale": 0.5}, {}],
        },
    },
    # spheres: structure-level with scale
    "spheres_multi_structure": {
        "kind": "spheres",
        "parameters": {
            "global": {
                "centers": [[0, 0, 0], [1, 1, 1]],
            },
            "structure": [{"position": [5, 5, 5], "scale": 2.0}],
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
        "sphere shape 'radius' must be a numeric value, got <class 'str'>",
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
        "ellipsoid shape 'semiaxes' must be an array with 3 values",
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
        "custom shape 'vertices' must be an Nx3 array",
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
        "custom shape 'simplices' must be an Nx3 array",
    ],
    [
        {"kind": "garbage", "parameters": {"radius": None}},
        ValueError,
        "unknown shape kind 'garbage'",
    ],
    # cylinders: unknown parameter
    [
        {
            "kind": "cylinders",
            "parameters": {
                "global": {"vectors": [[1, 0, 0]], "vector": [1, 0, 0]},
            },
        },
        ValueError,
        "unknown shape parameter 'vector' for 'cylinders' shape kind",
    ],
    # cylinders: vectors not Nx3
    [
        {
            "kind": "cylinders",
            "parameters": {
                "global": {"vectors": [[1, 0], [0, 1]]},
            },
        },
        ValueError,
        "cylinders shape 'vectors' must be an Nx3 array",
    ],
    # cylinders: bases wrong length
    [
        {
            "kind": "cylinders",
            "parameters": {
                "global": {
                    "vectors": [[1, 0, 0], [0, 1, 0]],
                    "bases": [[0, 0, 0]],
                },
            },
        },
        ValueError,
        "cylinders shape 'bases' must be an 2x3 array",
    ],
    # cylinders: radii wrong length
    [
        {
            "kind": "cylinders",
            "parameters": {
                "global": {
                    "vectors": [[1, 0, 0], [0, 1, 0]],
                    "radii": [0.1],
                },
            },
        },
        ValueError,
        "cylinders shape 'radii' must be a scalar or array of length 2",
    ],
    # spheres: unknown parameter
    [
        {
            "kind": "spheres",
            "parameters": {
                "global": {"centers": [[0, 0, 0]], "radius": 0.1},
            },
        },
        ValueError,
        "unknown shape parameter 'radius' for 'spheres' shape kind",
    ],
    # spheres: centers not Nx3
    [
        {
            "kind": "spheres",
            "parameters": {
                "global": {"centers": [[1, 0], [0, 1]]},
            },
        },
        ValueError,
        "spheres shape 'centers' must be an Nx3 array",
    ],
    # spheres: radii wrong length
    [
        {
            "kind": "spheres",
            "parameters": {
                "global": {
                    "centers": [[0, 0, 0], [1, 1, 1]],
                    "radii": [0.1, 0.2, 0.3],
                },
            },
        },
        ValueError,
        "spheres shape 'radii' must be a scalar or array of length 2",
    ],
]


class TestShapes(unittest.TestCase):
    def test_custom_shapes(self):
        data = chemiscope.create_input(
            structures=[DEFAULT_STRUCTURE], shapes=SHAPE_DEFAULTS
        )

        result = data["shapes"]
        self.assertEqual(
            list(result.keys()),
            list(SHAPE_DEFAULTS.keys()),
        )

        for key, values in result.items():
            self.assertEqual(SHAPE_DEFAULTS[key], values)

    def test_shape_errors(self):
        for shape, errortype, message in INVALID_SHAPES:
            with self.assertRaises(errortype) as cm:
                chemiscope.create_input(
                    structures=[DEFAULT_STRUCTURE], shapes={"shape": shape}
                )
                self.assertEqual(cm.message, message)


if __name__ == "__main__":
    unittest.main()
