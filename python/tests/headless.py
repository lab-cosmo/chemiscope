import gzip
import json
import os
import shutil
import struct
import tempfile
import unittest

from chemiscope.headless import headless


class TestHeadless(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Create a simple dataset for testing
        cls.structures = [
            {
                "size": 2,
                "names": ["H", "H"],
                "x": [0.0, 1.0],
                "y": [0.0, 0.0],
                "z": [0.0, 0.0],
            },
            {
                "size": 3,
                "names": ["H", "H", "O"],
                "x": [0.0, 1.0, 0.5],
                "y": [0.0, 0.0, 1.0],
                "z": [0.0, 0.0, 0.0],
            },
        ]
        cls.properties = {
            "index": {"target": "structure", "values": [0, 1]},
            "energy": {"target": "structure", "values": [-13.4, -15.2]},
            "dummy": {"target": "atom", "values": [1, 2, 3, 4, 5]},
        }

    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp()
        self.widget = headless(
            structures=self.structures, properties=self.properties, mode="default"
        )

    def tearDown(self):
        self.widget.close()
        shutil.rmtree(self.tmp_dir)

    def _is_valid_png(self, path):
        if not os.path.exists(path):
            return False
        with open(path, "rb") as f:
            header = f.read(8)
        # Check PNG signature
        return header == b"\x89PNG\r\n\x1a\n"

    def test_save_structure_image(self):
        path = os.path.join(self.tmp_dir, "structure.png")
        self.widget.save_structure_image(path)
        self.assertTrue(self._is_valid_png(path))

    def test_save_map_image(self):
        path = os.path.join(self.tmp_dir, "map.png")
        self.widget.save_map_image(path)
        self.assertTrue(self._is_valid_png(path))

    def test_get_structure_sequence(self):
        indices = [0, 1]
        settings = [{"structure": [{"bonds": False}]}, {"structure": [{"bonds": True}]}]

        images = self.widget.get_structure_sequence(indices, settings=settings)
        self.assertEqual(len(images), 2)

        # Check that we got valid PNG data
        for img_data in images:
            self.assertTrue(img_data.startswith(b"\x89PNG\r\n\x1a\n"))

    def test_save_structure_sequence(self):
        indices = [0, 1]
        paths = [
            os.path.join(self.tmp_dir, "seq_0.png"),
            os.path.join(self.tmp_dir, "seq_1.png"),
        ]

        self.widget.save_structure_sequence(indices, paths)

        for path in paths:
            self.assertTrue(self._is_valid_png(path))

    def test_save_json(self):
        path = os.path.join(self.tmp_dir, "data.json")
        self.widget.save(path)

        with open(path, "r") as f:
            data = json.load(f)

        self.assertIn("meta", data)
        self.assertIn("structures", data)
        self.assertIn("properties", data)

    def test_save_gz_json(self):
        path = os.path.join(self.tmp_dir, "data.json.gz")
        self.widget.save(path)

        self.assertTrue(os.path.exists(path))
        self.assertGreater(os.path.getsize(path), 0)

        with gzip.open(path, "rt") as f:
            data = json.load(f)

        self.assertIn("meta", data)

    def test_missing_visualizer_errors(self):
        # Test structure-only mode
        widget = headless(structures=self.structures, mode="structure")
        try:
            with self.assertRaisesRegex(
                RuntimeError,
                "Cannot retrieve map image: this widget is a structure-only viewer.",
            ):
                widget.get_map_image()

            with self.assertRaisesRegex(
                RuntimeError,
                "Cannot retrieve map image: this widget is a structure-only viewer.",
            ):
                widget.save_map_image("dummy.png")
        finally:
            widget.close()

        # Test map-only mode
        widget = headless(
            structures=self.structures, properties=self.properties, mode="map"
        )
        try:
            with self.assertRaisesRegex(
                RuntimeError,
                "Cannot retrieve structure image: this widget is a map-only viewer.",
            ):
                widget.get_structure_image()

            with self.assertRaisesRegex(
                RuntimeError,
                "Cannot retrieve structure image: this widget is a map-only viewer.",
            ):
                widget.save_structure_image("dummy.png")

            with self.assertRaisesRegex(
                RuntimeError,
                "Cannot retrieve structure image: this widget is a map-only viewer.",
            ):
                widget.get_structure_sequence([0])

            with self.assertRaisesRegex(
                RuntimeError,
                "Cannot retrieve structure image: this widget is a map-only viewer.",
            ):
                widget.save_structure_sequence([0], ["dummy.png"])

        finally:
            widget.close()


class TestHeadlessDimensions(unittest.TestCase):
    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp()
        self.structures = [
            {
                "size": 1,
                "names": ["H"],
                "x": [0.0],
                "y": [0.0],
                "z": [0.0],
            }
        ]

    def tearDown(self):
        shutil.rmtree(self.tmp_dir)

    def _get_dims(self, path):
        with open(path, "rb") as f:
            f.seek(16)
            width = struct.unpack(">I", f.read(4))[0]
            height = struct.unpack(">I", f.read(4))[0]
        return width, height

    def test_default_vs_explicit(self):
        # Default (800x600 container)
        widget = headless(structures=self.structures, mode="structure")
        path1 = os.path.join(self.tmp_dir, "default.png")
        widget.save_structure_image(path1)
        widget.close()
        w1, h1 = self._get_dims(path1)

        # Explicit smaller (400x300 container)
        widget = headless(
            structures=self.structures, mode="structure", width=400, height=300
        )
        path2 = os.path.join(self.tmp_dir, "explicit.png")
        widget.save_structure_image(path2)
        widget.close()
        w2, h2 = self._get_dims(path2)

        # Ensure explicit is smaller than default
        self.assertLess(w2, w1 * 0.6)
        self.assertLess(h2, h1 * 0.6)

        # Ensure aspect ratio is roughly preserved (4:3)
        # Allow 20% tolerance due to scrollbars/padding/DPI effects
        ratio1 = w1 / h1
        ratio2 = w2 / h2
        self.assertAlmostEqual(ratio1, 4 / 3, delta=0.2)
        self.assertAlmostEqual(ratio2, 4 / 3, delta=0.2)

    def test_single_dimension(self):
        # width=400 -> height=300 (4:3)
        widget = headless(structures=self.structures, mode="structure", width=400)
        path = os.path.join(self.tmp_dir, "width.png")
        widget.save_structure_image(path)
        widget.close()
        w, h = self._get_dims(path)

        # Check aspect ratio
        self.assertAlmostEqual(w / h, 4 / 3, delta=0.2)


if __name__ == "__main__":
    unittest.main()
