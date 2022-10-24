import nibabel
import numpy as np


class Nifti1Image:

    __VOXEL_DIMENSION = 3

    def __init__(self, file_path) -> None:
        self.__image = self.__load_image(file_path)

    def data_file(self) -> np.memmap:
        return self.__image.get_fdata()

    def bounding_box(self) -> np.ndarray:
        starts = np.array([self.__image.header['qoffset_x'],
                            self.__image.header['qoffset_y'],
                            self.__image.header['qoffset_z']])
        ends = starts + np.multiply(self.__xyz_shape(), self.__xyz_dimension())
        return np.array([starts, ends]) * self.__scaling()

    def header(self) -> nibabel.nifti1.Nifti1Header:
        return self.__image.header

    def values_at(self, positions: list) -> list:
        shape = (len(positions), 1, 1, 1, self.__VOXEL_DIMENSION)
        positions_reshaped = np.array(positions).reshape(shape) 
        lower_limits = self.__voxel_positions()
        upper_limits = lower_limits + self.__xyz_dimension()

        is_above = np.all(lower_limits <= positions_reshaped, axis=4)
        is_below = np.all(positions_reshaped < upper_limits, axis=4)

        return [self.data_file()[inside][0]
                for inside in np.logical_and(is_above, is_below)]

    @staticmethod
    def __load_image(file_path):
        try: 
            return nibabel.load(file_path)
        except FileNotFoundError:
            raise IOError('File Not Found.')

    def __scaling(self) -> float:
        xyz_unit = self.__image.header.get_xyzt_units()[0]
        return {'unknown': 1.0, 
                'meter': 1000.0, 
                'mm': 1.0, 
                'micron': 0.001}[xyz_unit]

    def __voxel_positions(self) -> np.ndarray:
        x_max, y_max, z_max = self.__xyz_shape()
        shape = (x_max, y_max, z_max,self.__VOXEL_DIMENSION)
        coordinates = np.array([(x,y,z) for x in range(x_max)
                                        for y in range(y_max) 
                                        for z in range(z_max)]).reshape(shape)
        return coordinates * self.__xyz_dimension()

    def __xyz_dimension(self) -> tuple:
        return self.__image.header.get_zooms()[:self.__VOXEL_DIMENSION]

    def __xyz_shape(self):
        return self.__image.header.get_data_shape()[:self.__VOXEL_DIMENSION]
