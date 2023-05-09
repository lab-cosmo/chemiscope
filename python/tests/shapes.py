import unittest

import ase
import numpy as np

import chemiscope

from chemiscope.structures._shapes import SHAPE_PARAMS
from chemiscope.structures._ase import _ase_extract_shapes

BASE_FRAME = ase.Atoms(numbers=np.zeros(3), positions=[[0, 0, 0], [1, 1, 1], [2, 2, 5]])
BASE_FRAME.arrays["orientation"] = [[1, 0, 0, 0] for _ in BASE_FRAME]
CUBE = [
    [0, 0, 0],
    [1, 0, 0],
    [0, 1, 0],
    [1, 1, 0],
    [0, 0, 1],
    [1, 0, 1],
    [0, 1, 1],
    [1, 1, 1],
]
CUBE_SPLX = [
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
    "ellipsoid": {"semiaxes": [1, 2, 1]},
    "custom": {"vertices": CUBE, "simplices": CUBE_SPLX},
    "sphere": {"radius": 1.0},
}


class TestStructures(unittest.TestCase):
    """Conversion of structure data to chemiscope JSON"""

    def test_no_shape(self):
        frame = BASE_FRAME.copy()
        data = chemiscope.create_input(frames=[frame])
        self.assertTrue("shape" not in data["structures"][0])

    def test_bad_shape(self):
        frame = BASE_FRAME.copy()
        frame.info["shape"] = "invalid"

        with self.assertRaises(KeyError) as cm:
            _ = _ase_extract_shapes([frame])
            self.assertEquals(
                cm.message,
                "The currently-supported shape types are ellipsoid, sphere, custom, received invalid.",
            )

    def test_by_frame(self):
        for key in SHAPE_DEFAULTS:
            with self.subTest(shape=key):
                frame = BASE_FRAME.copy()
                frame.info["shape"] = key
                for k, v in SHAPE_DEFAULTS[key].items():
                    frame.info[f"shape_{k}"] = v
                shapes = _ase_extract_shapes([frame])
                print(shapes)
                data = chemiscope.create_input(frames=[frame], shapes=shapes)
                self.assertTrue("shape" in data["structures"][0])
                self.assertTrue(
                    all(
                        [
                            k in ss
                            for s in data["structures"]
                            for ss in s["shape"]["shape"]
                            for k in SHAPE_PARAMS[key]["required"]
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
                for k in SHAPE_PARAMS[key]["optional"]:
                    frame.info.pop(f"shape_{k}", None)
                    with self.assertWarns(Warning) as cm:
                        _ = chemiscope.create_input(
                            frames=[frame], shapes=_ase_extract_shapes([frame])
                        )
                        self.assertEquals(
                            cm.message,
                            "Missing optional parameter {} for {}.".format(k, key),
                        )
                    frame.info[f"shape_{k}"] = SHAPE_DEFAULTS[key][k]

                for k in SHAPE_PARAMS[key]["computable"]:
                    frame.arrays.pop(f"shape_{k}", None)
                    data = chemiscope.create_input(
                        frames=[frame], shapes=_ase_extract_shapes([frame])
                    )
                    self.assertTrue(
                        all(
                            [
                                k in ss
                                for s in data["structures"]
                                for ss in s["shape"]["shape"]
                            ]
                        )
                    )
                    frame.arrays[f"shape_{k}"] = [SHAPE_DEFAULTS[key][k] for _ in frame]

                for k in SHAPE_PARAMS[key]["required"]:
                    frame.info.pop(f"shape_{k}")
                    with self.assertRaises(KeyError) as cm:
                        _ = chemiscope.create_input(
                            frames=[frame], shapes=_ase_extract_shapes([frame])
                        )
                        self.assertEquals(
                            cm.message,
                            "Missing required parameter {} for {}".format(k, key),
                        )
                    frame.info[f"shape_{k}"] = SHAPE_DEFAULTS[key][k]

    def test_by_index(self):
        for key in SHAPE_DEFAULTS:
            with self.subTest(shape=key):
                frame = BASE_FRAME.copy()
                frame.arrays["shape"] = [key for _ in frame]
                for k, v in SHAPE_DEFAULTS[key].items():
                    frame.arrays[f"shape_{k}"] = [v for _ in frame]
                shapes = _ase_extract_shapes([frame])
                data = chemiscope.create_input(frames=[frame], shapes=shapes)
                self.assertTrue("shape" in data["structures"][0])
                self.assertTrue(
                    all(
                        [
                            k in ss
                            for s in data["structures"]
                            for ss in s["shape"]["shape"]
                            for k in SHAPE_PARAMS[key]["required"]
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
                for k in SHAPE_PARAMS[key]["optional"]:
                    frame.arrays.pop(f"shape_{k}", None)
                    with self.assertWarns(Warning) as cm:
                        _ = chemiscope.create_input(
                            frames=[frame], shapes=_ase_extract_shapes([frame])
                        )
                        self.assertEquals(
                            cm.message,
                            "Missing optional parameter {} for {}.".format(k, key),
                        )
                    frame.arrays[f"shape_{k}"] = [SHAPE_DEFAULTS[key][k] for _ in frame]

                for k in SHAPE_PARAMS[key]["computable"]:
                    frame.arrays.pop(f"shape_{k}", None)
                    data = chemiscope.create_input(
                        frames=[frame], shapes=_ase_extract_shapes([frame])
                    )
                    self.assertTrue(
                        all(
                            [
                                k in ss
                                for s in data["structures"]
                                for ss in s["shape"]["shape"]
                            ]
                        )
                    )
                    frame.arrays[f"shape_{k}"] = [SHAPE_DEFAULTS[key][k] for _ in frame]

                for k in SHAPE_PARAMS[key]["required"]:
                    frame.arrays.pop(f"shape_{k}")
                    with self.assertRaises(KeyError) as cm:
                        _ = chemiscope.create_input(
                            frames=[frame], shapes=_ase_extract_shapes([frame])
                        )
                        self.assertEquals(
                            cm.message,
                            "Missing required parameter {} for {}".format(k, key),
                        )
                    frame.arrays[f"shape_{k}"] = [SHAPE_DEFAULTS[key][k] for _ in frame]

    def test_no_shapes(self):
        frame = BASE_FRAME.copy()
        with self.assertRaises(ValueError) as cm:
            _ = _ase_extract_shapes([frame])
            self.assertEqual(cm.message, "1 frame(s) do not contain shape information.")


if __name__ == "__main__":
    unittest.main()
