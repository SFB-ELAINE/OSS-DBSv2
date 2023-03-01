
from typing import List
import numpy as np


class Region:

    def __init__(self, start: tuple, shape: tuple, voxel_size: tuple) -> None:
        self.start = start
        self.shape = shape
        self.voxel_size = voxel_size

    @property
    def end(self):
        return tuple(np.array(self.voxel_size) * self.shape + self.start)

    def points(self, offset: tuple, voxel_size: tuple) -> List:
        offset = np.array(offset) % self.start
        print(offset)
        shape = self.shape * np.divide(self.voxel_size, voxel_size)
        start = np.floor(np.divide(self.start, self.voxel_size)) + offset

        x_values = start[0] + np.arange(shape[0]) * voxel_size[0]
        y_values = start[1] + np.arange(shape[1]) * voxel_size[1]
        z_values = start[2] + np.arange(shape[2]) * voxel_size[2]

        return [(x, y, z)
                for x in x_values
                for y in y_values
                for z in z_values]
