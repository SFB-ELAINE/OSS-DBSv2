import nibabel
import numpy as np


class Nifti1Image:

    __VOXEL_DIMENSION = 3

    def __init__(self, file_path: str) -> None:
        self._image = self.__load_image(file_path)
        self._shift = (0.0, 0.0, 0.0)

    def data_map(self) -> np.memmap:
        return self._image.get_fdata()

    def bounding_box(self) -> tuple:
        start = self.offset()
        shape = np.array(self.xyz_shape(), dtype=np.float64)
        ends = start + shape * self.voxel_size()
        return tuple(start), tuple(ends)

    def header(self) -> nibabel.nifti1.Nifti1Header:
        return self._image.header

    def offset(self) -> tuple:
        offset = np.array([self._image.header['qoffset_x'],
                           self._image.header['qoffset_y'],
                           self._image.header['qoffset_z']
                           ], dtype=np.float64)
        return offset * self.__scaling() + self._shift

    def set_offset(self, shift: tuple) -> None:
        self._shift = shift

    def voxel_size(self) -> tuple:
        x, y, z = self._image.header.get_zooms()[:self.__VOXEL_DIMENSION]
        return tuple(np.array((x, y, z), dtype=np.float64) * self.__scaling())

    def xyz_shape(self) -> tuple:
        return self._image.header.get_data_shape()[:self.__VOXEL_DIMENSION]

    @staticmethod
    def __load_image(file_path: str) -> nibabel.nifti1.Nifti1Image:
        try:
            return nibabel.load(file_path)
        except FileNotFoundError:
            raise IOError('File Not Found.')

    def __scaling(self) -> float:
        xyz_unit = self._image.header.get_xyzt_units()[0]
        return {'unknown': 1e-3,
                'meter': 1.0,
                'mm': 1e-3,
                'micron': 1e-6}[xyz_unit]
