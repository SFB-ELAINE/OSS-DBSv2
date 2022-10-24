from src.Nifti1Image import Nifti1Image
import numpy as np


class DiffusionTensorImage(Nifti1Image):

    def diffusion(self) -> np.ndarray:
        x_max, y_max, z_max = self._xyz_shape()
        diffusion = self.data_map().reshape((x_max * y_max * z_max, 6))
        return np.array([(xx, xy, xz, xy, yy, yz, xz, yz, zz)
                         for xx, xy, xz, yy, yz, zz in diffusion]
                        ).reshape(x_max, y_max, z_max, 3, 3).astype(float)
