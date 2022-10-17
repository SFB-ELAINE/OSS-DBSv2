import nibabel
import numpy as np
from typing import Any


class Nifti1Image():

    def __init__(self, file_path) -> None:
        self.__image = self.__load_image(file_path)

    @staticmethod
    def __load_image(file_path):
        try: 
            return nibabel.load(file_path)
        except FileNotFoundError:
            raise IOError('File Not Found.')

    def data_file(self) -> np.memmap:
        return self.__image.get_fdata()

    def dimensions(self) -> tuple[int]:
        return self.__image.header.get_zooms()

    def bounding_box(self) -> np.ndarray:
        starts = np.array([self.__image.header['qoffset_x'],
                            self.__image.header['qoffset_y'],
                            self.__image.header['qoffset_z']])
        ends = starts + np.multiply(self.shape(), self.dimensions())[:3]
        return np.array([starts, ends]) * self.__scaling()

    def header(self) -> nibabel.nifti1.Nifti1Header:
        return self.__image.header

    def shape(self) -> tuple[int]:
        return self.__image.header.get_data_shape()  

    def values_at(self, positions: np.ndarray) -> list[Any]:
        positions_reshaped = positions.reshape((len(positions), 1, 1, 1, 3)) 
        voxel_positions = self.__voxel_positions()
        
        lower_limits = voxel_positions * self.dimensions()
        upper_limits = lower_limits + self.dimensions()

        is_above = np.all(lower_limits <= positions_reshaped, axis=4)
        is_below = np.all(positions_reshaped < upper_limits, axis=4)
        return [self.data_file()[inside][0]
                for inside in np.logical_and(is_above, is_below)]

    def __scaling(self) -> float:
        xyz_unit = self.__image.header.get_xyzt_units()[0]
        return {'unknown': 1.0, 
                'meter': 1000.0, 
                'mm': 1.0, 
                'micron': 0.001}[xyz_unit]

    def __voxel_positions(self) -> np.ndarray:
        x_max, y_max, z_max = self.shape()
        coordinates = [(x,y,z) for x in range(x_max) 
                               for y in range(y_max) 
                               for z in range(z_max)]
        positions = np.empty((x_max, y_max, z_max, 3))

        for coordinate in coordinates:
             positions[coordinate] = coordinate

        return positions