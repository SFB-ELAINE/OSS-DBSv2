from src.brain_imaging.Nifti1Image import Nifti1Image
import numpy as np


class DiffusionTensorImage(Nifti1Image):

    def __init__(self, file_path: str) -> None:
        super().__init__(file_path)
        if not self._image.get_fdata().ndim == 4:
            raise IOError('Data Shape not four dimensional.')
        if not self._image.get_fdata().shape[-1] == 6:
            raise IOError('Data Shape is not (x,y,z,6).')

    def diffusion(self) -> np.ndarray:
        x_max, y_max, z_max = self.xyz_shape()
        diffusion = self.data_map().reshape((x_max * y_max * z_max, 6))
        return np.array([(xx, xy, xz, xy, yy, yz, xz, yz, zz)
                         for xx, xy, xz, yy, yz, zz in diffusion]
                        ).reshape(x_max, y_max, z_max, 3, 3).astype(float)
