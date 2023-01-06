from dataclasses import dataclass
import numpy as np


@dataclass
class Voxels:
    """Represents the collection of Voxels in space.

    Attributes
    ----------
    data : np.ndarray
    start : tuple
    end : tuple

    """
    data: np.ndarray = np.empty(0)
    start: tuple = (0, 0, 0)
    end: tuple = (0, 0, 0)
