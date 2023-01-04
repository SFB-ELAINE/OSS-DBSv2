
from ossdbs.brainsubstance import Material
from ossdbs.brain_imaging.mri import MagneticResonanceImage
from ossdbs.dielectric_model import ColeColeFourModelCreator
from ossdbs.voxels import Voxels
import numpy as np


class Conductivity:

    __MATERIALS = [Material.CSF, Material.GRAY_MATTER, Material.WHITE_MATTER]

    def __init__(self, mri: MagneticResonanceImage) -> None:
        self.__mri = mri
        self.__complex = False

    def conductivity(self, frequency: float) -> Voxels:
        omega = 2 * np.pi * frequency
        colecole_model = ColeColeFourModelCreator.create(Material.GRAY_MATTER)
        default = colecole_model.conductivity(omega)
        data = np.full(self.__mri.xyz_shape(), default)

        for material in self.__MATERIALS:
            position = self.__mri.material_distribution(material=material)
            colecole_model = ColeColeFourModelCreator.create(material)
            data[position.data] = colecole_model.conductivity(omega)

        start, end = self.__mri.bounding_box()

        if not self.__complex:
            return Voxels(data=np.real(data), start=start, end=end)

        return Voxels(data=data, start=start, end=end)

    def set_complex(self, value: bool) -> None:
        self.__complex = value
