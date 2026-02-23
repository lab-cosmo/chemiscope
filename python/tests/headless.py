import gzip
import json
import os
import shutil
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
                "Cannot save map image: this widget is a structure-only viewer.",
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
                "Cannot save structure image: this widget is a map-only viewer.",
            ):
                widget.save_structure_image("dummy.png")

            with self.assertRaisesRegex(
                RuntimeError,
                "Cannot retrieve structure image: this widget is a map-only viewer.",
            ):
                widget.get_structure_sequence([0])

            with self.assertRaisesRegex(
                RuntimeError,
                "Cannot save structure image: this widget is a map-only viewer.",
            ):
                widget.save_structure_sequence([0], ["dummy.png"])

        finally:
            widget.close()


if __name__ == "__main__":
    unittest.main()
