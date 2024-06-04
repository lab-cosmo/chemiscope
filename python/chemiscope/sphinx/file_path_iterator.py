import os
import re

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

        # Continue the counter
        self._set_counter_based_on_existing_files()

    def _set_counter_based_on_existing_files(self):
        """
        Set the counter to one more than the highest existing file number.
        If the directory is used by the multiple sphinx configurations, we don't want to
        erase the files starting the generation of filename from the zero again
        """
        pattern = re.compile(
            rf"{re.escape(self.prefix)}_(\d+){re.escape(self.extension)}"
        )
        existing_files = [f for f in os.listdir(self.base_dir) if pattern.match(f)]
        if existing_files:
            highest_number = max(int(pattern.match(f).group(1)) for f in existing_files)
            self.counter = highest_number

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
