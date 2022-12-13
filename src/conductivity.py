
from src.brain_imaging.magnetic_resonance_imaging import MagneticResonanceImage
from src.brainsubstance import Material
from src.dielectric_model.dielectric_model_var1 import ModelCreator
from src.voxels import Voxels
import numpy as np


class Conductivity:

    __MATERIALS = [Material.CSF, Material.GRAY_MATTER, Material.WHITE_MATTER]

    def __init__(self, mri: MagneticResonanceImage) -> None:
        self.__mri = mri

    def complex_conductivity(self, frequency: float) -> Voxels:
        omega = 2 * np.pi * frequency
        default = ModelCreator.create(Material.GRAY_MATTER).conductivity(omega)
        data = np.full(self.__mri.xyz_shape(), default)

        for material in self.__MATERIALS:
            position = self.__mri.material_distribution(material=material)
            conductivity = ModelCreator.create(material).conductivity(omega)
            data[position.data] = conductivity

        start, end = self.__mri.bounding_box()
        return Voxels(data=data, start=start, end=end)
