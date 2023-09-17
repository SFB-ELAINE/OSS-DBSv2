import os


class OutputDirectory:
    """Represnts output directory for result files and generates output
    directory if directory does not exist.

    Attributes
    ----------
    directory : str
        Output directory.
    """

    def __init__(self, directory) -> None:
        self._directory = directory

    def directory(self):
        """Return output directory and generates directory if given directory
        does not exist.
        """
        file_dir = self._directory if self._directory else "result"
        os.makedirs(file_dir, exist_ok=True)
        return file_dir
