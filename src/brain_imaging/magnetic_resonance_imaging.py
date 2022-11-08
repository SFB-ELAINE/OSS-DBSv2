from src.brain_imaging.Nifti1Image import Nifti1Image
import numpy as np


class MagneticResonanceImage(Nifti1Image):
    pass


class DefaultMagneticResonanceImage(Nifti1Image):

    def __init__(self, file_path: str) -> None:
        self.__x = 50
        self.__y = 50
        self.__z = 50
        box = np.full((self.__x, self.__y, self.__z), 1)
        self.__array = np.empty((2*self.__x, 2*self.__y, 2*self.__z))
        self.__array[:self.__x, :self.__y, :self.__z] = box * 0
        self.__array[self.__x:, :self.__y, :self.__z] = box * 1
        self.__array[:self.__x, self.__y:, :self.__z] = box * 2
        self.__array[self.__x:, self.__y:, :self.__z] = box * 3
        self.__array[:self.__x, :self.__y, self.__z:] = box * 3
        self.__array[self.__x:, :self.__y, self.__z:] = box * 2
        self.__array[:self.__x, self.__y:, self.__z:] = box * 1
        self.__array[self.__x:, self.__y:, self.__z:] = box * 0

    def data_map(self) -> np.array:
        return self.__array

    def bounding_box(self) -> np.ndarray:
        starts = np.array([0, 0, 0])
        ends = starts + np.multiply(self._xyz_shape(), self._xyz_dimension())
        return tuple(np.array([starts, ends]))

    def _xyz_dimension(self) -> tuple:
        return (0.1, 0.1, 0.1)

    def _xyz_shape(self):
        return (2 * self.__x, 2 * self.__y, 2 * self.__z)
