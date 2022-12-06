from src.brain_imaging.Nifti1Image import Nifti1Image
from src.brainsubstance import Material
from src.voxels import Voxels
import numpy as np


class MagneticResonanceImage(Nifti1Image):

    def __init__(self, file_path: str, material_coding: dict = None) -> None:
        super().__init__(file_path)

        if not self._image.get_fdata().ndim == 3:
            raise IOError('Data Shape not three dimensional.')

        if not material_coding:
            material_coding = {Material.GRAY_MATTER: 1,
                               Material.WHITE_MATTER: 2,
                               Material.CSF: 3,
                               Material.UNKNOWN: 0}

        self.__coding = material_coding

    def data_map(self) -> np.memmap:
        data = self._image.get_fdata()
        new_data = np.empty(data.shape)
        for material in self.__coding.keys():
            new_data[data == self.__coding[material]] = material
        return new_data

    def material_distribution(self, material: Material) -> Voxels:
        start, end = self.bounding_box()
        is_material = self.data_map() == material
        return Voxels(data=is_material, start=tuple(start), end=tuple(end))
