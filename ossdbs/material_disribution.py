
from ossdbs.materials import Material
import numpy as np


class MaterialDistribution:

    def __init__(self,
                 data: np.ndarray,
                 offset: tuple,
                 voxel_size: tuple
                 ) -> None:
        self.__data = data
        self.__offset = offset
        self.__voxel_size = voxel_size

    def where(self, material: Material) -> tuple:
        x_position, y_position, z_position = np.where(self.__data == material)
        x_start = x_position * self.__voxel_size[0] + self.__offset[0]
        y_start = y_position * self.__voxel_size[1] + self.__offset[1]
        z_start = z_position * self.__voxel_size[2] + self.__offset[2]
        starts = np.array([x_start, y_start, z_start]).T
        ends = starts + self.__voxel_size
        return starts, ends
