import unittest
import tempfile
import json
import ase
import os

from chemiscope import write_input

TEST_FRAMES = [ase.Atoms("CO2")]


def is_gz_file(filepath):
    # check for the magic number at the beggining of a file
    with open(filepath, "rb") as test_f:
        return test_f.read(2) == b"\x1f\x8b"


class TestWriteInput(unittest.TestCase):
    def test_file_path_as_dataset_name(self):
        with tempfile.TemporaryDirectory() as dirname:
            path = os.path.join(dirname, "test.json")
            write_input(path, TEST_FRAMES)

            with open(path) as fd:
                data = json.load(fd)

            self.assertEqual(data["meta"]["name"], "test")

    def test_create_gz_file(self):
        with tempfile.TemporaryDirectory() as dirname:
            path = os.path.join(dirname, "test.json")
            write_input(path, TEST_FRAMES)
            self.assertFalse(is_gz_file(path))

            path = os.path.join(dirname, "test.json.gz")
            write_input(path, TEST_FRAMES)
            self.assertTrue(is_gz_file(path))

    def test_wrong_path(self):
        with tempfile.TemporaryDirectory() as dirname:
            path = os.path.join(dirname, "test.tmp")
            with self.assertRaises(Exception) as cm:
                write_input(path, TEST_FRAMES)

            self.assertEqual(
                str(cm.exception), "path should end with .json or .json.gz"
            )


if __name__ == "__main__":
    unittest.main()
