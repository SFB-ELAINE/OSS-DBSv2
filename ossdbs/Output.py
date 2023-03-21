
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
        self.__directory = directory

    def output_directory(self):
        """Return output directory and generates directory if given directory
        does not exist.
        """
        file_dir = self.__directory if self.__directory else 'result'
        os.makedirs(file_dir, exist_ok=True)
        return file_dir
