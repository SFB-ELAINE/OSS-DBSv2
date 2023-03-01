
from dataclasses import dataclass
from typing import List
import numpy as np


@dataclass
class BoundingBox:

    def __init__(self,
                 start: tuple = (0, 0, 0),
                 end: tuple = (0, 0, 0)) -> None:
        self.start = start
        self.end = end

    @property
    def shape(self):
        return tuple(np.round(np.subtract(self.end, self.start)).astype(int))

    def intersection(self, bbox: 'BoundingBox') -> 'BoundingBox':
        x_s = max(self.start[0], bbox.start[0])
        y_s = max(self.start[1], bbox.start[1])
        z_s = max(self.start[2], bbox.start[2])
        x_e = min(self.end[0], bbox.end[0])
        y_e = min(self.end[0], bbox.end[1])
        z_e = min(self.end[1], bbox.end[2])
        return BoundingBox(start=(x_s, y_s, z_s), end=(x_e, y_e, z_e))

    def points(self, offset: tuple, voxel_size: tuple) -> List:
        start, end = self.start, self.end
        start_index = np.floor(np.subtract(start, offset) / voxel_size)
        end_index = np.ceil(np.subtract(end, offset) / voxel_size)
        shape = end_index - start_index
        x_s, y_s, z_s = start_index * voxel_size + offset
        x_values = x_s + np.arange(shape[0]) * voxel_size[0]
        y_values = y_s + np.arange(shape[1]) * voxel_size[1]
        z_values = z_s + np.arange(shape[2]) * voxel_size[2]
        return [(x, y, z)
                for x in x_values
                for y in y_values
                for z in z_values]

    def __eq__(self, other: 'BoundingBox') -> bool:
        return self.start == other.start and self.end == other.end
