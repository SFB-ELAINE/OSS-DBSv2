from ossdbs.brain_imaging.Nifti1Image import Nifti1Image
import numpy as np


class DiffusionTensorImage(Nifti1Image):

    def diffusion(self) -> np.ndarray:
        x_max, y_max, z_max = self._xyz_shape()
        diffusion = self.data_map().reshape((x_max * y_max * z_max, 6))
        return np.array([(xx, xy, xz, xy, yy, yz, xz, yz, zz)
                         for xx, xy, xz, yy, yz, zz in diffusion]
                        ).reshape(x_max, y_max, z_max, 3, 3).astype(float)


class DefaultDiffusionTensorImage(Nifti1Image):

    def __init__(self) -> None:
        self.__x = 50
        self.__y = 50
        self.__z = 50

    def diffusion(self) -> np.ndarray:
        return np.array([np.eye(3)])

    def bounding_box(self) -> np.ndarray:
        starts = np.array([0, 0, 0])
        ends = starts + np.multiply(self._xyz_shape(), self._xyz_dimension())
        return tuple(np.array([starts, ends]))

    def _xyz_dimension(self) -> tuple:
        return (0.1, 0.1, 0.1)

    def _xyz_shape(self):
        return (2 * self.__x, 2 * self.__y, 2 * self.__z)
