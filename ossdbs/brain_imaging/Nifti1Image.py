import nibabel
import numpy as np


class Nifti1Image:

    __VOXEL_DIMENSION = 3

    def __init__(self, file_path: str) -> None:
        self.__image = self.__load_image(file_path)

    def data_map(self) -> np.memmap:
        return self.__image.get_fdata()

    def bounding_box(self) -> np.ndarray:
        starts = np.array([self.__image.header['qoffset_x'],
                           self.__image.header['qoffset_y'],
                           self.__image.header['qoffset_z']])
        ends = starts + np.multiply(self._xyz_shape(), self._xyz_dimension())
        return tuple(np.array([starts, ends]) * self.__scaling())

    def header(self) -> nibabel.nifti1.Nifti1Header:
        return self.__image.header

    @staticmethod
    def __load_image(file_path: str):
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

    def _xyz_dimension(self) -> tuple:
        return self.__image.header.get_zooms()[:self.__VOXEL_DIMENSION]

    def _xyz_shape(self):
        return self.__image.header.get_data_shape()[:self.__VOXEL_DIMENSION]
