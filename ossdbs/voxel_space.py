from dataclasses import dataclass
import numpy as np


@dataclass
class VoxelSpace:
    data: np.ndarray = np.empty(0)
    start: tuple = (0, 0, 0)
    end: tuple = (0, 0, 0)
