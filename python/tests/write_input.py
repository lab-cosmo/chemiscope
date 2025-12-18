import json
import os
import tempfile
import unittest

import ase
import stk

from chemiscope import write_input


TEST_STRUCTURES = [ase.Atoms("CO2")]
TEST_STRUCTURES_STK = [stk.BuildingBlock("NCCN")]


def is_gz_file(filepath):
    # check for the magic number at the beginning of a file
    with open(filepath, "rb") as test_f:
        return test_f.read(2) == b"\x1f\x8b"


class TestWriteInput(unittest.TestCase):
    def test_file_path_as_dataset_name(self):
        for structure in (TEST_STRUCTURES, TEST_STRUCTURES_STK):
            with tempfile.TemporaryDirectory() as dirname:
                path = os.path.join(dirname, "test.json")
                write_input(path, structure)

                with open(path) as fd:
                    data = json.load(fd)

                self.assertEqual(data["meta"]["name"], "test")

    def test_create_gz_file(self):
        for structure in (TEST_STRUCTURES, TEST_STRUCTURES_STK):
            with tempfile.TemporaryDirectory() as dirname:
                path = os.path.join(dirname, "test.json")
                write_input(path, structure)
                self.assertFalse(is_gz_file(path))

                path = os.path.join(dirname, "test.json.gz")
                write_input(path, structure)
                self.assertTrue(is_gz_file(path))

    def test_wrong_path(self):
        for structure in (TEST_STRUCTURES, TEST_STRUCTURES_STK):
            with tempfile.TemporaryDirectory() as dirname:
                path = os.path.join(dirname, "test.tmp")
                with self.assertRaises(Exception) as cm:
                    write_input(path, structure)
                self.assertEqual(
                    str(cm.exception), "path should end with .json or .json.gz"
                )


if __name__ == "__main__":
    unittest.main()
