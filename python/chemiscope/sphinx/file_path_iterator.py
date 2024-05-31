import os

from sphinx.errors import ExtensionError


class FilePathIterator:
    def __init__(self, base_dir, prefix="fig", extension=".json.gz"):
        """
        Initialize the iterator

        Parameters:
        - base_dir (str): The directory where files will be saved
        - prefix (str): The prefix for image filenames
        - extension (str): The file extension for files
        """
        self.base_dir = base_dir
        self.prefix = prefix
        self.extension = extension
        self.counter = 0
        self._stop = 1000000

        # Ensure the base directory exists
        os.makedirs(base_dir, exist_ok=True)

    def __iter__(self):
        """Iterate over paths"""
        for _i in range(self._stop):
            yield self.next()
        else:
            raise ExtensionError(f"Generated over {self._stop} files")

    def __next__(self):
        """Generates the file name"""
        self.counter += 1
        filename = f"{self.prefix}_{self.counter:03d}{self.extension}"
        return os.path.join(self.base_dir, filename)

    def next(self):
        """Return the next file path, with numbering starting at 1"""
        return self.__next__()
